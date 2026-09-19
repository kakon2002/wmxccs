"""One ion, end to end, with real numbers from the steroid corpus.

Run: python tools/demo_end_to_end.py

WHAT THIS IS FOR
----------------
Every stage of this platform is tested in isolation. This runs ONE ion through all of
them in order and prints what each stage did to it, so that somebody who has not read the
code can see where a harmonized number comes from and what travels with it.

It is a demonstration, not a test: it asserts nothing and proves nothing. The assertions
are in tests/. What this does is make the pipeline legible.

Nothing here is invented. Every number printed is read from data/seed/steroid_jasms2022.csv,
which was converted from the published supporting information of
DOI 10.1021/jasms.2c00196 by tools/ingest_steroid.py.
"""

from __future__ import annotations

import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from wmxccs.harmonization import fit_harmonization, harmonize, measure_coverage  # noqa: E402
from wmxccs.loader import load_measurements_file  # noqa: E402
from wmxccs.matching import build_matched_ions  # noqa: E402
from wmxccs.scope import Claim, ScopeExceededError, assert_may_be_quoted_as  # noqa: E402
from wmxccs.statistics import compare_platforms  # noqa: E402

SEED = REPO / "data" / "seed" / "steroid_jasms2022.csv"
RULE = "=" * 92


def head(number: int, title: str) -> None:
    print(f"\n{RULE}\n{number}. {title}\n{RULE}")


def main() -> int:
    # --- 1. INGESTION ---------------------------------------------------------------------
    head(1, "INGESTION: a published table becomes records, or is refused with a reason")
    report = load_measurements_file(SEED)
    print(f"   file                {SEED.relative_to(REPO)}")
    print(f"   rows read           {report.rows_read}")
    print(f"   records built       {report.records_built}   (rows failed: {report.rows_failed})")
    print(f"   cleared to train    {report.records_cleared}")
    print(f"   held back           {report.records_refused}")
    for reason, count in report.gate_counts.items():
        print(f"      {count} x {reason}")
    print("   The licence gate ran FIRST. Only cleared records go any further.")

    # --- 2. IDENTITY ----------------------------------------------------------------------
    head(2, "IDENTITY: what makes two measurements measurements OF THE SAME ION")
    # One compound measured on all four platforms, chosen for the demonstration by being
    # the first such compound in the file rather than by being convenient.
    matching = build_matched_ions(report.cleared)
    ion = next(i for i in matching.matched if len(i.platforms) == 4)
    compound = ion.key.analyte[2].removeprefix("steroid_jasms2022:")
    print(f"   compound            {compound}")
    print(f"   adduct              {ion.key.adduct}   charge {ion.key.charge:+d}")
    print(f"   gas the value refers to  {ion.key.drift_gas}")
    print(f"   matched-ion key     {ion.key}")
    print("   The key carries NO platform, which is why values from four platforms land on it.")
    print(f"   identity is the dataset-scoped id {ion.key.analyte[2]!r},")
    print("   which pairs inside this source and matches nothing outside it.")

    # --- 3. MATCHED-ION CONSTRUCTION -------------------------------------------------------
    head(3, "MATCHED IONS: the same ion on four platforms, held apart, never merged")
    for measurement in sorted(ion.measurements, key=lambda m: m.platform):
        record = measurement.record
        spread = (
            f"+/- {record.ccs_uncertainty:g} ({record.uncertainty_type.value})"
            if record.ccs_uncertainty
            else "no spread reported"
        )
        print(
            f"   {measurement.platform:20s} {record.ccs:8.3f} A^2   {spread:28s}"
            f" calibrant: {record.calibrant or '(none: primary)'}"
        )
    print(f"   cross-platform matched ions in the corpus: {len(matching.matched)}")
    print(f"   of which span all four platforms:          "
          f"{sum(1 for i in matching.matched if len(i.platforms) == 4)}")

    # --- 4. STATISTICS ---------------------------------------------------------------------
    head(4, "STATISTICS: association and agreement, per platform pair AND per calibration group")
    comparison = compare_platforms(matching)
    print(f"   platform pairs {len(comparison.pairs)}   strata {sum(len(p.strata) for p in comparison.pairs)}"
          f"   paired points {comparison.n_points}")
    print("   NO POOLED FIGURE for any pair: the three adducts are three calibration groups.")
    for pair in comparison.pairs:
        for stratum in pair.strata:
            if stratum.points[0].ion.key.adduct != ion.key.adduct:
                continue
            print(
                f"   {stratum.other_platform:20s} vs {stratum.reference_platform:20s}"
                f" n={stratum.n:3d}  median {stratum.agreement.median_difference_percent:+6.3f}%"
                f"  Lin CCC {stratum.agreement.lins_ccc:+.4f}"
            )

    # --- 5. THE CORRECTION -----------------------------------------------------------------
    head(5, "THE CORRECTION: three ways, and which one the data supports")
    model = fit_harmonization(comparison)
    twims = next(m.record for m in ion.measurements if m.platform == "TWIMS")
    result = harmonize(model, twims)
    correction = result.correction
    print(f"   correcting the {twims.ims_type.value} value of this ion, referred to"
          f" {correction.reference_platform}")
    print(f"   stratum             n={correction.n}, fitted over this adduct only")
    print(f"   Deming slope        {correction.stratum.agreement.deming_slope:.5f}"
          f"   intercept {correction.stratum.agreement.deming_intercept:+.3f}")
    print(f"   robust slope        {correction.robust.slope:.5f}"
          f"   interval {correction.robust.slope_interval}")
    print(f"   median offset       {correction.median_offset_percent:+.3f}%")
    print(f"   slope differs from 1?  "
          f"{'yes' if correction.robust.slope_distinguishable_from_unity else 'NO'}"
          f"  ->  headline basis: {result.basis.value}")
    if not correction.robust.slope_distinguishable_from_unity:
        print("      The slope cannot be told from 1, so a slope-derived correction and a constant")
        print("      offset describe the same data. The offset is the one that does not extrapolate.")

    # --- 6. THE ANSWER ---------------------------------------------------------------------
    head(6, "THE ANSWER: the original untouched, the correction beside it, and its interval")
    band = result.interval
    print(f"   ORIGINAL            {result.original_ccs:8.3f} A^2   <- returned unchanged, always")
    print(f"   harmonized          {result.harmonized_ccs:8.3f} A^2   ({result.basis.value})")
    print(f"      slope-derived    {result.slope_derived_ccs:8.3f}")
    print(f"      median-derived   {result.median_derived_ccs:8.3f}")
    print(f"      robust-derived   {result.robust_derived_ccs:8.3f}")
    print(f"   shift               {100 * (result.harmonized_ccs - result.original_ccs) / result.original_ccs:+.3f}%")
    print(f"   interval            {band.low:.3f} to {band.high:.3f} A^2"
          f"   (width {band.width:.3f}, {100 * band.width / band.at_value:.2f}% of CCS)")
    print(f"   coverage            nominal {band.nominal_coverage:.0%},"
          f" GUARANTEED {band.guaranteed_coverage:.0%}   [jackknife+ proves 1-2*alpha]")
    print(f"   informative?        {band.interval_is_informative}"
          f"   tail ratio {band.tail_ratio:.2f}")
    print(f"   grade               {result.confidence.grade.value}")
    for demotion in result.confidence.demotions:
        print(f"      {demotion.rule}: {demotion.detail[:96]}")
    for note in result.confidence.not_checked:
        print(f"      NOT CHECKED - {note[:96]}")

    # --- 7. WHAT IT MAY BE CALLED -----------------------------------------------------------
    head(7, "SCOPE: what this number is, and what it is not")
    print(f"   {result.scope.caveat()}")
    print()
    for claim in Claim:
        try:
            assert_may_be_quoted_as(result.scope, claim)
            print(f"   MAY be quoted as    {claim.value}")
        except ScopeExceededError:
            print(f"   may NOT be quoted as {claim.value}")
    print()
    check = measure_coverage(correction.stratum)
    print(f"   leave-one-out coverage on this stratum: {check.covered}/{check.n} = {check.rate:.3f}")
    print(f"   arithmetically pinned at {check.pinned_rate:.3f} -> identical: {check.is_pinned}")
    print(f"   is that evidence of calibration?  {check.is_evidence_of_calibration}")
    print("   It is not. The informative figures are the interval width and the tail ratio above.")
    print()
    print(f"   maturity            {model.maturity.data_maturity.value},"
          f" {model.maturity.matched_ion_count} matched ions")
    print("   Provisional, and it cannot be otherwise while the corpus is one study.")
    print()
    print(f"   model version       {model.fingerprint.short}")
    print(f"      corpus     {model.fingerprint.corpus}")
    print(f"      parameters {model.fingerprint.parameters}")
    print("   Two digests: the records behind the fit, and everything that decides an answer.")
    print("   Quote them with any number from here, or it cannot be reproduced - the model is")
    print("   refitted at every startup and nothing else records that the data or code moved.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
