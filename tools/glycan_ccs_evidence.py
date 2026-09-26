"""The adapter between the CCS core and the glycan layer. The ONLY thing that imports both.

`src/wmxccs` and `src/wmxglycan` never import each other, in either direction, and
tests/test_glycan_boundary.py enforces that. This module is not in either package: it is
outside both, it imports both, and it translates. That is the whole design - the wall stays
up and something outside it carries the traffic.

WHAT IT CAN AND CANNOT FIND, MEASURED RATHER THAN ASSUMED

`reachable_states()` computes which evidence states this adapter can actually produce
against the corpus it is pointed at, and the answer today is uncomfortable and worth
publishing:

  HELD_NOT_RELEASABLE      REACHABLE. The 89 Struwe 2015 records are N-glycans,
                           Hex3-Hex10HexNAc2. Hex5HexNAc2 [M-H]- alone has 8 of them.
  NONE_IN_CORPUS_SEARCHED  REACHABLE, for any ion the seed files do not hold.
  MEASURED_REFERENCE       NOT REACHABLE. Not because the code path is missing, but because
                           the two seed files are exactly complementary: struwe2016 has 24
                           CLEARED records with structure strings and NO composition, and
                           every one is a milk oligosaccharide on a lactose core, which no
                           N-glycan composition matches; struwe2015 has 89 records WITH
                           compositions and all of them are held.
  LOOKUP_FAILED            REACHABLE, and deliberately so: a lookup that raises must not
                           degrade into an absence.

The unreachability of MEASURED_REFERENCE is held to the same standard as AI_ONLY in the
ranker: it is DERIVED from the corpus by counting, so the statement changes on its own the
day a cleared record with a composition arrives, and a test fails at that moment rather than
the claim quietly becoming false.

WHY A HELD VALUE IS NEVER SHOWN

The 89 records are licence-clean - `open_attribution` - and still may not be shown, because
what blocks them is not the licence. Their drift gas is `UNSTATED` and their
`uncertainty_type` is `unknown` with no uncertainty value, and hard constraint 7 refuses a
number without an uncertainty type. Reading Hofmann 2014 for the gas would clear one blocker
and leave the other exactly where it is, because the uncertainty is absent from the source.
So `any_blocker_resolves_alone` is False, and it is a field rather than a sentence because
"one blocker" invites the reading "resolve this and you have data".
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO / "src") not in sys.path:  # pragma: no cover - a convenience for `python tools/...`
    sys.path.insert(0, str(REPO / "src"))

from wmxccs.loader import load_measurements_file  # noqa: E402  (after the path fix, deliberately)
from wmxglycan.ccs_evidence import (  # noqa: E402
    CCSEvidence,
    CCSEvidenceState,
    CCSReference,
    EvidenceLevel,
    HeldValues,
)
from wmxglycan.composition import Composition, parse_composition  # noqa: E402

SEED = REPO / "data" / "seed"
CORPUS_NAME = "wmxccs seed corpus (data/seed)"

# What blocks the Struwe 2015 values, named once. Both, always, because naming one invites the
# reading that clearing it is enough.
STRUWE_2015_BLOCKERS = (
    "the drift gas is UNSTATED in the source, so the value is not comparable across platforms"
    " (hard constraint 4)",
    "uncertainty_type is `unknown` and no uncertainty value is given; the uncertainty is absent"
    " from the source rather than unread, and a number without a type is refused (hard"
    " constraint 7)",
)
STRUWE_2015_RELEASE = (
    "Nothing available today. The gas could be resolved by reading Hofmann 2014, but the"
    " uncertainty cannot be resolved by reading anything: it is not in the source. Releasing"
    " these would mean either inventing an uncertainty or dropping the requirement for one, and"
    " the owner has ruled out loosening the gate"
)


@dataclass(frozen=True)
class _Row:
    composition: str
    adduct: str
    charge: int
    doi: str
    source: str
    conformer: str
    conformers_total: str


class SeedCorpusEvidence:
    """CCS evidence for a glycan ion, read from the CCS core's seed corpus.

    Satisfies `wmxglycan.ccs_evidence.CCSEvidenceLookup` without either package knowing about
    the other.
    """

    def __init__(self, seed: Path | None = None) -> None:
        self.seed = SEED if seed is None else seed
        self._rows: list[_Row] = []
        self._cleared: list[object] = []
        for path in sorted(self.seed.glob("*.csv")):
            report = load_measurements_file(path)
            self._cleared.extend(report.cleared)
            with path.open(encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle):
                    if (row.get("analyte_kind") or "").strip() != "glycan":
                        continue
                    composition = (row.get("analyte_composition") or "").strip()
                    if not composition:
                        continue
                    try:
                        charge = int((row.get("charge") or "").strip())
                    except ValueError:
                        continue
                    self._rows.append(
                        _Row(
                            composition=parse_composition(composition).canonical,
                            adduct=(row.get("adduct") or "").strip(),
                            charge=charge,
                            doi=(row.get("doi") or "").strip(),
                            source=(row.get("source") or "").strip(),
                            conformer=(row.get("conformer") or "").strip(),
                            conformers_total=(row.get("conformers_total") or "").strip(),
                        )
                    )

    # --- what the ranker asks for -----------------------------------------------------------------

    def evidence_for(self, composition: Composition, adduct: str, charge: int) -> CCSEvidence:
        wanted = composition.canonical
        cleared = self._cleared_for(wanted, adduct, charge)
        if cleared is not None:
            return cleared

        held = [
            row
            for row in self._rows
            if row.composition == wanted and row.adduct == adduct and row.charge == charge
        ]
        if not held:
            return CCSEvidence(
                state=CCSEvidenceState.NONE_IN_CORPUS_SEARCHED,
                level=EvidenceLevel.ION,
                composition=wanted,
                adduct=adduct,
                charge=charge,
                corpus_searched=CORPUS_NAME,
                records_consulted=len(self._rows),
            )

        conformers = sum(1 for row in held if row.conformer)
        note = STRUWE_2015_RELEASE
        if conformers:
            note += (
                f". {conformers} of these {len(held)} rows are numbered conformers of one ion"
                " rather than independent measurements, so the record count is not a count of"
                " separate determinations"
            )
        return CCSEvidence(
            state=CCSEvidenceState.HELD_NOT_RELEASABLE,
            level=EvidenceLevel.ION,
            composition=wanted,
            adduct=adduct,
            charge=charge,
            held=HeldValues(
                records=len(held),
                blockers=STRUWE_2015_BLOCKERS,
                any_blocker_resolves_alone=False,
                doi=held[0].doi or None,
                what_would_release_them=note,
            ),
        )

    def _cleared_for(self, wanted: str, adduct: str, charge: int) -> CCSEvidence | None:
        """A cleared record for this ion, if one exists and is fit to be shown.

        Returns None rather than a state, so the caller distinguishes "no cleared record" from
        "a cleared record that cannot be rendered", which are different problems.
        """
        for record in self._cleared:
            analyte = getattr(record, "analyte", None)
            composition = getattr(analyte, "composition", None)
            if composition is None:
                continue
            text = composition if isinstance(composition, str) else str(composition)
            if parse_composition(text).canonical != wanted:
                continue
            if str(record.adduct) != adduct or int(record.charge) != charge:
                continue
            reference = CCSReference(
                ccs=float(record.ccs),
                uncertainty=float(record.ccs_uncertainty),
                uncertainty_type=str(record.uncertainty_type),
                adduct=str(record.adduct),
                charge=int(record.charge),
                polarity=str(record.polarity),
                ims_type=str(record.ims_type),
                drift_gas=str(record.drift_gas),
                calibrant=record.calibrant,
                source=str(record.source),
                doi=str(record.doi or ""),
                source_locator=str(record.source_locator or ""),
                reuse_status=record.reuse_status,
                measured_on_canonical_key=None,
            )
            return CCSEvidence(
                state=CCSEvidenceState.MEASURED_REFERENCE,
                # COMPOSITION level, not STRUCTURE: the record gives a composition and an ion, so
                # the value belongs to all the isomers equally and to none of them in particular.
                level=EvidenceLevel.COMPOSITION,
                composition=wanted,
                adduct=adduct,
                charge=charge,
                reference=reference,
                discriminates_between_candidates=False,
            )
        return None

    # --- what this adapter can actually produce, counted rather than claimed ------------------------

    def reachable_states(self) -> dict[str, str]:
        """Which states this corpus can produce, derived by counting. See the module docstring."""
        with_composition = sum(
            1
            for record in self._cleared
            if getattr(getattr(record, "analyte", None), "composition", None) is not None
        )
        return {
            CCSEvidenceState.HELD_NOT_RELEASABLE.value: (
                f"REACHABLE: {len(self._rows)} glycan row(s) in the seed corpus carry a"
                " composition and are held"
            ),
            CCSEvidenceState.NONE_IN_CORPUS_SEARCHED.value: "REACHABLE for any ion not in the corpus",
            CCSEvidenceState.MEASURED_REFERENCE.value: (
                f"REACHABLE: {with_composition} cleared record(s) carry a composition"
                if with_composition
                else (
                    f"NOT REACHABLE: 0 of {len(self._cleared)} cleared record(s) carry a"
                    " composition. The cleared glycan values are milk oligosaccharides on a"
                    " lactose core identified only by a structure string, and no N-glycan"
                    " composition matches one. This is counted, not assumed, so it changes by"
                    " itself the day such a record arrives"
                )
            ),
            CCSEvidenceState.LOOKUP_FAILED.value: (
                "REACHABLE, and deliberately: a lookup that raises is reported as a failure and"
                " never as an absence"
            ),
        }

    @property
    def rows_consulted(self) -> int:
        return len(self._rows)

    @property
    def cleared_consulted(self) -> int:
        return len(self._cleared)


def main() -> int:  # pragma: no cover - a hand-run report
    evidence = SeedCorpusEvidence()
    print(f"{CORPUS_NAME}: {evidence.rows_consulted} glycan rows with a composition,"
          f" {evidence.cleared_consulted} cleared records")
    print("\nwhat this adapter can produce:")
    for state, why in evidence.reachable_states().items():
        print(f"  {state:26} {why}")
    print("\nthe four reported compositions, as [M-H]- singly charged:")
    for text in ("Hex3HexNAc4Fuc1", "Hex4HexNAc4Fuc1", "Hex5HexNAc4Fuc1", "Hex5HexNAc2"):
        found = evidence.evidence_for(parse_composition(text), "[M-H]-", -1)
        print(f"  {text:18} {found.state.value:24} {found.summary()[:110]}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
