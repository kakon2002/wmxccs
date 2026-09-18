"""Entry point for `python -m tools.mutation`. See --help, or runner.USAGE."""

import sys

from .runner import main

if __name__ == "__main__":
    sys.exit(main())
