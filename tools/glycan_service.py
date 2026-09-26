"""The composition root: the glycan service wired to the CCS core's evidence.

THE ONLY PLACE THE TWO PACKAGES MEET, together with the adapter it uses. `src/wmxglycan` and
`src/wmxccs` never import each other in either direction, and
`tests/test_glycan_boundary.py` enforces that. This file is in neither package, imports both,
and hands one to the other through the `CCSEvidenceLookup` Protocol the glycan layer declares.

Run it:

    .venv/Scripts/python -m uvicorn tools.glycan_service:app --port 8010
    .venv/Scripts/python tools/glycan_service.py            # prints what is wired, serves nothing

WHY THE DEFAULT `wmxglycan.api:app` IS NOT THIS ONE

`wmxglycan.api:app` is built with no evidence source, so it reports NOT_CONSULTED for every
cross section. That is the honest state of a deployment nobody configured, and it is
distinguishable from an absence - which is the entire reason the evidence enum has five states
rather than two. Wiring the adapter inside `wmxglycan.api` would have meant that module
importing `wmxccs`, which is the one thing the two-package structure exists to prevent.

So there are two apps and the difference between them is a fact about the deployment rather
than about the code: one has been pointed at a CCS corpus and one has not, and every response
says which.

WHAT THIS WIRING CAN AND CANNOT FIND, and it is counted rather than claimed

`SeedCorpusEvidence.reachable_states()` reports it. Against the corpus as it stands:
HELD_NOT_RELEASABLE is reachable (89 Struwe 2015 N-glycan rows, 8 of them for Hex5HexNAc2
[M-H]-), NONE_IN_CORPUS_SEARCHED is reachable, LOOKUP_FAILED is reachable, and
MEASURED_REFERENCE is NOT - because all 24 cleared glycan records are milk oligosaccharides on
a lactose core carrying no composition, so no N-glycan composition matches one. That is derived
by counting, so the statement changes by itself the day such a record arrives.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
# Both, and both are needed: `src` so `wmxglycan` imports, and the repository root so
# `tools.glycan_ccs_evidence` does when this file is run directly as a script rather than
# imported as `tools.glycan_service`. Running it directly failed on exactly that.
for entry in (REPO / "src", REPO):  # pragma: no cover - a convenience for `python tools/...`
    if str(entry) not in sys.path:
        sys.path.insert(0, str(entry))

from wmxglycan.api import DEFAULT_DATABASE, create_app  # noqa: E402  (after the path fix)
from wmxglycan.store import Store  # noqa: E402

from tools.glycan_ccs_evidence import SeedCorpusEvidence  # noqa: E402


def build(database: Path | str = DEFAULT_DATABASE, seed: Path | None = None):
    """The wired application. `database` is a file, so a restart keeps every frozen prediction."""
    return create_app(store=Store(database), evidence=SeedCorpusEvidence(seed))


# The served application, wired. `uvicorn tools.glycan_service:app` from the repository root.
app = build()


def main() -> int:  # pragma: no cover - a hand-run report
    evidence = SeedCorpusEvidence()
    store = Store(DEFAULT_DATABASE)
    print("wmxglycan service, wired to the CCS seed corpus")
    print(f"  database                 {DEFAULT_DATABASE}")
    print(f"  predictions frozen       {len(store.predictions(limit=1_000_000))}")
    print(f"  glycan rows consulted    {evidence.rows_consulted}")
    print(f"  cleared records          {evidence.cleared_consulted}")
    print()
    print("  what the CCS evidence lookup can produce against this corpus:")
    for state, why in evidence.reachable_states().items():
        print(f"    {state:26} {why}")
    print()
    print("  serve it with:")
    print("    .venv/Scripts/python -m uvicorn tools.glycan_service:app --port 8010")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
