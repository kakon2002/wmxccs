"""Applying the catalogue, judging each result, and refusing to round up.

THE REPOSITORY IS NEVER WRITTEN.

An earlier version of this engine mutated `src/wmxccs` in place and restored
it in a finally block. That is careful, and it is still a foot-gun: a hard kill
between the write and the restore leaves the working tree holding deliberate
defects, and nothing in the tree says so. So the package is copied to a
temporary directory, the mutation is applied to the COPY, and the copy is put
first on PYTHONPATH where it shadows the editable install. Killing the process
at any moment loses a temporary directory and nothing else.

That design brings its own silent failure, and it is the worst one available
here: if the shadow does not take, every run measures UNMUTATED code, every
mutation survives, and the report reads as a catastrophic loss of coverage
rather than as a broken harness. So the runner asks where `import wmxccs`
actually resolves and refuses to continue unless the answer is inside the
shadow.

The engine is separated from the suite it runs: `sweep` takes a callable, so
the decision logic can be tested in milliseconds without starting pytest. That
matters here more than usual, because this tool's whole job is to notice when
something silently stops working, and a tool with no tests of its own would be
the least defensible thing in the repository.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

from . import Anchor, Expect, Mutation, Verdict, repeated_labels
from .catalogue import MUTATIONS

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC = PROJECT_ROOT / "src" / "wmxccs"

# pytest's own: 0 all passed, 1 some test failed, 2 interrupted, 3 internal
# error, 4 usage error, 5 nothing collected. Only 1 means a test NOTICED.
PYTEST_TESTS_FAILED = 1

SuiteRun = Callable[[Path], tuple[int, list[str]]]

USAGE = """Run the mutation catalogue against the test suite.

    python -m tools.mutation                 the full sweep, roughly three quarters of an hour
    python -m tools.mutation --check         anchors only, about a second
    python -m tools.mutation --list          what is in the catalogue
    python -m tools.mutation "[S] " conform  only mutations whose label contains one of these

The repository is never written: mutations are applied to a copy that shadows
the installed package. Exit status is zero only when every mutation selected
identified exactly one site, changed the file, and did what it was expected
to do."""

KNOWN_FLAGS = ("--check", "--list", "--help")


class ShadowFailed(RuntimeError):
    """The mutated copy is not what `import wmxccs` resolves to.

    Refused rather than reported, because the alternative is a sweep measuring
    unmutated code in which every mutation survives - a result that looks like
    a total collapse of the test suite and is really a broken harness.
    """


@dataclass(frozen=True)
class Outcome:
    """What happened to one mutation, and whether that is acceptable."""

    mutation: Mutation
    anchor: Anchor
    verdict: Verdict | None = None  # None when the mutation was never applied
    killer: str = ""  # the first FAILED line, so a false kill can be spotted
    exit_code: int | None = None  # pytest's own, so a crash can be told from a failure

    @property
    def ok(self) -> bool:
        return self.anchor is Anchor.OK and self.verdict is self.mutation.expected_verdict

    @property
    def problem(self) -> str | None:
        """Why this outcome fails the run, or None if it is fine."""
        if self.anchor is Anchor.STALE:
            return (
                f"anchor matches nothing in {self.mutation.file}: the guard moved or went away."
                " Re-anchor it, or delete it and say what now covers its intent"
            )
        if self.anchor is Anchor.AMBIGUOUS:
            return (
                f"anchor matches more than once in {self.mutation.file}, so which site it would hit"
                " is arbitrary. Lengthen it until it identifies one"
            )
        if self.anchor is Anchor.INERT:
            return f"applied to {self.mutation.file} and changed nothing, so it tests nothing"
        if self.verdict is Verdict.UNRUNNABLE:
            return (
                f"the suite could not run with this applied (pytest exit {self.exit_code}), so nothing was"
                " tested. A mutation that breaks collection reports a kill it did not earn; make it a change"
                " to behaviour rather than one that stops the code importing"
            )
        if self.verdict is Verdict.SURVIVED and self.mutation.expect is Expect.KILLED:
            return "NO TEST CAUGHT THIS: the behaviour has nothing behind it"
        if self.verdict is Verdict.KILLED and self.mutation.expect is Expect.SURVIVES:
            return (
                "a test now catches this, but it is recorded as an expected survivor."
                f" The exemption has rotted and should be deleted. It said: {self.mutation.reason}"
            )
        return None


# --- the shadow ---------------------------------------------------------------------


def shadow_of(src: Path, into: Path) -> Path:
    """A fresh copy of the package under `into`, ready to go first on PYTHONPATH."""
    root = into / "shadow"
    package = root / src.name
    if package.exists():
        shutil.rmtree(package)
    shutil.copytree(src, package, ignore=shutil.ignore_patterns("__pycache__"))
    return root


def _env_for(shadow: Path) -> dict[str, str]:
    env = dict(os.environ)
    existing = env.get("PYTHONPATH")
    env["PYTHONPATH"] = str(shadow) + (os.pathsep + existing if existing else "")
    # Writing .pyc files into the shadow would be harmless but pointless, and a
    # stale one in the real tree would not be.
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return env


def where_package_resolves(shadow: Path, package_name: str) -> str:
    result = subprocess.run(
        [sys.executable, "-c", f"import {package_name} as p; print(p.__file__)"],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        env=_env_for(shadow),
    )
    return (result.stdout or "").strip()


def shadowing_worked(resolved: str, package: Path) -> bool:
    """Whether the import really landed inside the shadow.

    Split out so the judgement can be tested on its own. The failure that
    matters here - the import succeeding but finding the INSTALLED package
    instead of the copy - is not something a test can honestly arrange, because
    PYTHONPATH winning is the very thing being relied upon. What a test can do
    is check that the right conclusion is drawn from each possible answer.
    """
    if not resolved:
        return False
    return Path(resolved).resolve().is_relative_to(package.resolve())


def verified_shadow(src: Path, into: Path) -> Path:
    """A shadow that `import <package>` demonstrably resolves to, or ShadowFailed."""
    shadow = shadow_of(src, into)
    package = shadow / src.name
    resolved = where_package_resolves(shadow, src.name)
    if not shadowing_worked(resolved, package):
        raise ShadowFailed(
            f"`import {src.name}` resolves to {resolved or '(nothing)'}, which is not inside"
            f" {package.resolve()}. Every result would be measured against unmutated code, so every"
            " mutation would look like a survivor and the report would read as a collapse of the test"
            " suite rather than a broken tool"
        )
    return shadow


# --- running the suite ----------------------------------------------------------------


def pytest_command() -> list[str]:
    """The whole suite, stopping at the first failure.

    The whole suite, deliberately. An earlier harness passed a hand-written list
    of test files and warned in its own comment that a missing one would score a
    mutation SURVIVED for the wrong reason. That list is exactly the kind of
    thing that goes stale quietly, so there is no list: `-x` keeps a killed
    mutation cheap, because it stops at the first test that notices.

    Nothing is excluded either. While the harness mutated the repository in
    place, the catalogue's own anchor tests had to be ignored - they assert that
    every anchor matches, which is false mid-run when a neighbouring line has
    moved. Mutating a copy removes the conflict rather than working around it.
    """
    return [sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider"]


def run_suite(shadow: Path) -> tuple[int, list[str]]:
    result = subprocess.run(
        pytest_command(), cwd=PROJECT_ROOT, capture_output=True, text=True, env=_env_for(shadow)
    )
    return result.returncode, (result.stdout or "").strip().splitlines()


def _first_failure(lines: Sequence[str]) -> str:
    return next((line for line in lines if line.startswith("FAILED")), lines[-1] if lines else "")


# Wider than the longest status word, so a label can never abut it. AMBIGUOUS is
# nine characters and a nine-wide field ran it straight into the label.
_STATUS_WIDTH = 10


def _line(status: str, label: str) -> str:
    return f"  {status:<{_STATUS_WIDTH}}{label}"


def _detail(text: str) -> str:
    return " " * (2 + _STATUS_WIDTH) + str(text)


# --- the sweep ------------------------------------------------------------------------


def sweep(
    mutations: Sequence[Mutation],
    *,
    src: Path = SRC,
    run: SuiteRun = run_suite,
    log: Callable[[str], None] = print,
    verify: bool = True,
) -> list[Outcome]:
    """Apply each mutation to a shadowed copy, never to the repository.

    `src` is read and never written. The copy is mutated, the suite is run
    against it, and the copy is put back to pristine between mutations. A crash
    at any point loses a temporary directory; the working tree is untouched
    throughout, which is the property the in-place version could not offer
    however carefully it restored.
    """
    outcomes: list[Outcome] = []
    with tempfile.TemporaryDirectory(prefix="wmx-mutation-") as temporary:
        shadow = verified_shadow(src, Path(temporary)) if verify else shadow_of(src, Path(temporary))
        package = shadow / src.name
        pristine = {name: (package / name).read_text(encoding="utf-8") for name in {m.file for m in mutations}}
        for mutation in mutations:
            original = pristine[mutation.file]
            anchor = mutation.anchor_in(original)
            if anchor is not Anchor.OK:
                outcome = Outcome(mutation=mutation, anchor=anchor)
                outcomes.append(outcome)
                log(_line(anchor.value.upper(), mutation.label))
                log(_detail(outcome.problem))
                continue
            target = package / mutation.file
            target.write_text(mutation.apply_to(original), encoding="utf-8")
            try:
                code, lines = run(shadow)
            finally:
                target.write_text(original, encoding="utf-8")
            # Exactly one exit code means "a test noticed": pytest's 1. Zero is a
            # survivor; anything above one means the suite could not run, which
            # is not a kill however non-zero it looks.
            if code == 0:
                verdict = Verdict.SURVIVED
            elif code == PYTEST_TESTS_FAILED:
                verdict = Verdict.KILLED
            else:
                verdict = Verdict.UNRUNNABLE
            outcome = Outcome(
                mutation=mutation,
                anchor=anchor,
                verdict=verdict,
                killer=_first_failure(lines) if verdict is Verdict.KILLED else "",
                exit_code=code,
            )
            outcomes.append(outcome)
            if outcome.ok and verdict is Verdict.KILLED:
                log(_line("killed", mutation.label))
                log(_detail(outcome.killer[:120]))
            elif outcome.ok:
                log(_line("survived", f"{mutation.label}   (expected, documented)"))
                log(_detail(mutation.reason))
            else:
                log(_line("FAILED", mutation.label))
                log(_detail(outcome.problem))
    return outcomes


def select(mutations: Sequence[Mutation], wanted: Sequence[str]) -> list[Mutation]:
    """Mutations whose label contains any of the given substrings; all of them if none given."""
    if not wanted:
        return list(mutations)
    return [m for m in mutations if any(text in m.label for text in wanted)]


def report(outcomes: Iterable[Outcome], *, selected: int, total: int, log: Callable[[str], None] = print) -> int:
    """Print the summary and return the exit code. Non-zero if anything is not as expected."""
    outcomes = list(outcomes)
    failures = [outcome for outcome in outcomes if not outcome.ok]
    killed = sum(1 for outcome in outcomes if outcome.verdict is Verdict.KILLED)
    log("")
    if selected != total:
        log(f"FILTERED RUN: {selected} of {total} mutations selected. This is not a full sweep.")
    log(f"{killed} killed, {len(outcomes) - killed} not killed, of {len(outcomes)} run")
    if not failures:
        log("every mutation did what it was expected to do")
        return 0
    log(f"\n{len(failures)} PROBLEM(S):")
    for outcome in failures:
        log(f"  {outcome.mutation.label}")
        log(f"    {outcome.problem}")
    return 1


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    # An unrecognised option must be refused, never dropped. Selecting on the
    # remaining arguments and ignoring the rest would quietly start a full sweep
    # in answer to a typo - this tool's own failure mode, turned on itself.
    unknown = [arg for arg in argv if arg.startswith("--") and arg not in KNOWN_FLAGS]
    if unknown:
        print(f"unrecognised option(s): {' '.join(unknown)}\n")
        print(USAGE)
        return 1
    if "--help" in argv:
        print(USAGE)
        return 0
    check_only = "--check" in argv
    listing = "--list" in argv
    wanted = [arg for arg in argv if not arg.startswith("--")]

    duplicates = repeated_labels(MUTATIONS)
    if duplicates:
        print(f"the catalogue uses a label more than once, so a filtered run is ambiguous: {duplicates}")
        return 1

    chosen = select(MUTATIONS, wanted)
    if not chosen:
        print(f"no mutation label contains any of {wanted}")
        return 1

    if listing:
        for mutation in chosen:
            marker = "  " if mutation.expect is Expect.KILLED else "~ "
            print(f"{marker}{mutation.label}  [{mutation.file}]")
        print(f"\n{len(chosen)} of {len(MUTATIONS)} mutations")
        return 0

    if check_only:
        # Anchors only: no suite, no copy, about a second.
        problems = 0
        for mutation in chosen:
            anchor = mutation.anchor_in((SRC / mutation.file).read_text(encoding="utf-8"))
            if anchor is not Anchor.OK:
                problems += 1
                print(_line(anchor.value.upper(), mutation.label))
                print(_detail(Outcome(mutation=mutation, anchor=anchor).problem))
        print(f"\n{len(chosen) - problems} of {len(chosen)} anchors identify exactly one site")
        return 1 if problems else 0

    # The baseline runs through the shadow too, so a broken shadow or a copy
    # that cannot import fails in half a minute rather than after the first
    # mutation has been blamed for it.
    try:
        with tempfile.TemporaryDirectory(prefix="wmx-baseline-") as temporary:
            shadow = verified_shadow(SRC, Path(temporary))
            code, lines = run_suite(shadow)
    except ShadowFailed as failure:
        print(f"the shadow copy is not in effect, so nothing could be measured:\n  {failure}")
        return 1
    if code != 0:
        print("the suite is not green before mutating, so nothing here would mean anything; stopping")
        print("\n".join(lines[-8:]))
        return 1
    print(f"baseline: {lines[-1] if lines else '(no output)'}  (through the shadow)\n")

    outcomes = sweep(chosen)
    return report(outcomes, selected=len(chosen), total=len(MUTATIONS))
