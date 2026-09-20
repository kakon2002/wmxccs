"""`python -m wmxccs`: what an operator sees before the server takes over.

The entry point had no tests. It is thin - it parses three arguments, prints a banner and
hands off to uvicorn - but the banner is the only thing that tells an operator whether a
model was found, and an installation with an empty seed directory serves a working API whose
/harmonize answers 501. That is correct rather than broken, and it must not be a surprise.

THE BUG THESE WERE WRITTEN FOR. `uvicorn.run` never returns, and Python line-buffers stdout
only when stdout is a terminal. Redirected to a file or a pipe it block-buffers, so the
banner sat unflushed for the lifetime of the process: started with stdout to a file, the file
was EMPTY while uvicorn's own logging - which goes to stderr - came through normally. An
operator redirecting output to check what had loaded saw nothing, and the one reader who most
needs the banner is the one who sent it somewhere other than a screen.

WHY ONE OF THESE STARTS A REAL PROCESS. The property is "survives a stdout that is not a
terminal", and pytest's own capture replaces stdout with something that is not a terminal
either - so a test that patched `sys.stdout` would be asserting against a stand-in for the
thing that breaks. It was written that way first and passed with `-s` while failing under
capture, which is the clearest possible sign that it was measuring the harness. The
subprocess redirects a real file descriptor, which is what an operator does.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
import types

import pytest

from wmxccs.__main__ import main

BANNER_TIMEOUT_SECONDS = 60
UNLIKELY_PORT = 8765


@pytest.fixture
def stubbed_uvicorn(monkeypatch):
    """`main` with the server stubbed out, so it returns instead of blocking."""
    stub = types.ModuleType("uvicorn")
    stub.called_with = None

    def run(app, **kwargs):
        stub.called_with = (app, kwargs)

    stub.run = run
    monkeypatch.setitem(sys.modules, "uvicorn", stub)
    return stub


# --- the banner survives a stdout that is not a terminal ----------------------------------------


def test_the_banner_reaches_a_redirected_stdout_before_the_server_is_stopped(tmp_path):
    """The real thing: a real process, a real file descriptor, no monkeypatching.

    Without the flush this file is empty for as long as the process lives, because
    uvicorn.run never returns and there is nothing to force the buffer out.
    """
    log = tmp_path / "banner.txt"
    with log.open("wb") as handle:
        process = subprocess.Popen(
            [sys.executable, "-m", "wmxccs", "--port", str(UNLIKELY_PORT)],
            stdout=handle,
            stderr=subprocess.DEVNULL,
            cwd=os.fspath(tmp_path.parent),
        )
        try:
            deadline = time.monotonic() + BANNER_TIMEOUT_SECONDS
            text = ""
            while time.monotonic() < deadline:
                text = log.read_text(encoding="utf-8", errors="replace")
                if "serving on" in text:
                    break
                if process.poll() is not None:
                    break
                time.sleep(0.25)
        finally:
            process.terminate()
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:  # pragma: no cover - defensive
                process.kill()
                process.wait(timeout=30)

    text = log.read_text(encoding="utf-8", errors="replace")
    if "uvicorn is not installed" in text:
        pytest.skip("the serving extra is not installed in this environment")
    assert "wmxccs " in text, (
        "nothing reached the redirected stdout: the banner is buffered and uvicorn.run never"
        f" returns to flush it. File held {text!r}"
    )
    assert "serving on" in text


# --- what the banner says ------------------------------------------------------------------------


def test_the_banner_says_what_is_being_served(stubbed_uvicorn, capsys):
    assert main([]) == 0
    printed = capsys.readouterr().out
    assert "wmxccs " in printed
    assert "model loaded:" in printed or "NO MODEL LOADED" in printed
    assert f"serving on http://127.0.0.1:8000" in printed
    assert "docs at /docs" in printed


def test_the_banner_carries_the_model_version_an_answer_can_be_reproduced_from(
    stubbed_uvicorn, capsys
):
    main([])
    printed = capsys.readouterr().out
    if "NO MODEL LOADED" in printed:
        pytest.skip("this deployment found no seed data, so there is no version to print")
    assert "model version:" in printed
    assert "SCOPE:" in printed
    assert "NOT interlaboratory reproducibility" in printed


def test_the_default_host_is_loopback_and_not_every_interface(stubbed_uvicorn, capsys):
    """There is no authentication on this service, so the default bind is the only guard.

    Recorded in LIMITATIONS 7G as a deliberate omission rather than an oversight: the
    mechanism depends on a decision about who may call this. Until that decision is made,
    binding every interface by default would publish an unauthenticated service to the
    network, so it has to be asked for explicitly.
    """
    main([])
    assert stubbed_uvicorn.called_with[1]["host"] == "127.0.0.1"
    assert "127.0.0.1" in capsys.readouterr().out


def test_a_host_the_operator_names_is_honoured(stubbed_uvicorn):
    """The other direction. A hard-coded loopback would pass the test above and be a lie."""
    main(["--host", "0.0.0.0", "--port", "9999"])
    assert stubbed_uvicorn.called_with[1]["host"] == "0.0.0.0"
    assert stubbed_uvicorn.called_with[1]["port"] == 9999


def test_a_missing_uvicorn_is_reported_as_a_missing_extra_rather_than_a_traceback(monkeypatch):
    """The library does not need a server; only this entry point does.

    A ModuleNotFoundError traceback would read as a broken installation rather than as a
    missing optional dependency, and the fix is one pip command.
    """
    real_import = __import__

    def no_uvicorn(name, *args, **kwargs):
        if name == "uvicorn":
            raise ModuleNotFoundError("No module named 'uvicorn'")
        return real_import(name, *args, **kwargs)

    monkeypatch.delitem(sys.modules, "uvicorn", raising=False)
    monkeypatch.setattr("builtins.__import__", no_uvicorn)
    assert main([]) == 2
