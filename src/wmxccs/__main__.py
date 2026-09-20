"""Run the API. `python -m wmxccs`, or `wmxccs-serve` once installed.

Deliberately thin. Everything it does is start uvicorn on `wmxccs.api:app`, which is
already built with a model fitted from the seed corpus - so there is no separate
"load the model" step to forget, and no state this module holds.

It prints what it is serving BEFORE it starts, because an operator needs to know whether
a model was found. An installation whose data directory is empty serves a working API
whose /harmonize answers 501, and that is correct rather than broken - but it should not
be a surprise.

AND IT FLUSHES THAT BANNER, which is not decoration. `uvicorn.run` never returns, and
Python only line-buffers stdout when stdout is a terminal - so to a pipe, a log file or a
systemd journal the banner stayed in the buffer for the lifetime of the process. The one
reader who most needs it is the one who redirected it somewhere.
"""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m wmxccs",
        description="Serve the cross-platform CCS harmonization API.",
    )
    parser.add_argument("--host", default="127.0.0.1", help="default 127.0.0.1, not 0.0.0.0")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument(
        "--reload", action="store_true", help="development only; reloads on source changes"
    )
    args = parser.parse_args(argv)

    try:
        import uvicorn
    except ModuleNotFoundError:
        print(
            "uvicorn is not installed, so there is no server to start. Install the serving extra:\n"
            "    pip install -e .[serve]\n"
            "The library itself does not need it; only this entry point does.",
            file=sys.stderr,
        )
        return 2

    from . import __version__
    from .api import app

    model = app.state.model
    print(f"wmxccs {__version__}")
    if model is None:
        print(
            "  NO MODEL LOADED. /health and /confidence/rules answer normally;"
            " /harmonize will answer 501.\n"
            "  A model is fitted from the CSV files in data/seed - check that directory exists"
            " and holds seed data."
        )
    else:
        print(
            f"  model loaded: {len(model.applied)} applicable correction(s),"
            f" {model.maturity.matched_ion_count} matched ions,"
            f" maturity {model.maturity.data_maturity.value}"
        )
        print(f"  model version: {model.fingerprint.short}  (corpus/parameters, sha256)")
        print(f"  SCOPE: {model.scope.caveat()}")
    print(f"  serving on http://{args.host}:{args.port}  (docs at /docs)")

    # FLUSHED BEFORE THE SERVER BLOCKS, and this line is the whole reason the banner is
    # readable anywhere but a terminal. Python line-buffers stdout only when it is a tty;
    # to a pipe or a file it block-buffers, and uvicorn.run() below never returns - so the
    # banner sat in an 8 KB buffer until the process was killed. Measured on 20 September
    # 2026 by starting the server with stdout redirected to a file: uvicorn's own logging
    # goes to stderr and appeared, the banner did not appear at all, and an operator
    # checking whether a model had loaded saw an empty file. Which is exactly the surprise
    # this module's docstring says it exists to prevent.
    sys.stdout.flush()

    uvicorn.run("wmxccs.api:app", host=args.host, port=args.port, reload=args.reload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
