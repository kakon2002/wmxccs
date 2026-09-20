"""One ion, end to end, with real numbers from the steroid corpus.

Run: .venv/Scripts/python tools/demo_end_to_end.py   (.venv/bin/python on macOS and Linux)

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
import textwrap

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from wmxccs.harmonization import (  # noqa: E402
    fit_harmonization,
    harmonize,
    measure_coverage,
    outward,
    to_places,
)
from wmxccs.loader import load_measurements_file  # noqa: E402
from wmxccs.matching import build_matched_ions  # noqa: E402
from wmxccs.readiness import TARGET_MATCHED_IONS  # noqa: E402
from wmxccs.scope import Claim, ScopeExceededError, assert_may_be_quoted_as  # noqa: E402
from wmxccs.statistics import compare_platforms  # noqa: E402

SEED = REPO / "data" / "seed" / "steroid_jasms2022.csv"
RULE_WIDTH = 92
RULE = "=" * RULE_WIDTH


def head(number: int, title: str) -> None:
    print(f"\n{RULE}\n{number}. {title}\n{RULE}")


def wrap(text: str, indent: int = 3) -> None:
    """Print `text` inside the rule, wrapped on word boundaries.

    The grade reasons used to be cut at 96 characters with no ellipsis, which ended the
    first one at "though within 10% of i" and the second at "without a cave". A demonstration
    that truncates its own explanation mid-word is worse than one that omits it: the reader
    cannot tell whether the sentence was unimportant or whether the tool is broken.
    """
    lead = " " * indent
    for line in textwrap.wrap(text, width=RULE_WIDTH - indent, subsequent_indent="  "):
        print(f"{lead}{line}")


def explain_the_grade(model, result) -> None:
    """Why `weak` is the ordinary answer on this corpus and not a sign of breakage.

    216 of the 417 corrected records grade weak. A reader who meets one and assumes
    something failed will be wrong, and preventing that reading is this screen's job.

    Every figure here is DERIVED from the model in front of it rather than typed. Nothing
    tests this file, so a typed number would be a literal with nothing behind it - which is
    the failure class LIMITATIONS 4.5 exists to record.
    """
    applied = [correction.n for correction in model.corrections if correction.is_applied]
    print("   why weak, and why it is the ordinary answer here rather than a failure:")
    wrap(f"Every stratum this model applies holds between {min(applied)} and {max(applied)}"
         f" matched ions, short of the {TARGET_MATCHED_IONS} the scheme asks for before a"
         " platform-pair figure is quoted without a caveat. That rule therefore fires on EVERY"
         " corrected record in this corpus, so nothing here grades better than `qualified` and"
         " `supported` is out of reach by construction rather than by any fault of this ion.",
         indent=6)
    wrap("From there a grade falls one further notch for a second reason. For most records"
         " that reason is the outlier check, which cannot be evaluated for an ion this corpus"
         " has not measured. For this one it is that the submitted value sits just outside the"
         " range the correction was fitted over.", indent=6)
    wrap("So `weak` here means a small corpus and one unanswerable check - not a bad"
         " measurement, and not a broken pipeline.", indent=6)


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
    lo, hi = correction.robust.slope_interval
    print(f"   robust slope        {correction.robust.slope:.5f}"
          f"   interval ({lo:.4f}, {hi:.4f})")
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
    places = band.decimals_supported
    low, high = outward(band.low, band.high, places)
    shift = 100 * (result.harmonized_ccs - result.original_ccs) / result.original_ccs
    width_percent = 100 * band.width / band.at_value

    print(f"   ORIGINAL            {result.original_ccs:8.3f} A^2   <- returned unchanged, always")
    print(f"   harmonized          {to_places(result.harmonized_ccs, places):8.{places}f} A^2"
          f"   ({result.basis.value})")
    print(f"      slope-derived    {to_places(result.slope_derived_ccs, places):8.{places}f}")
    print(f"      median-derived   {to_places(result.median_derived_ccs, places):8.{places}f}")
    print(f"      robust-derived   {to_places(result.robust_derived_ccs, places):8.{places}f}")
    print(f"   shift               {shift:+.3f}%")
    print(f"   interval            {low:.{places}f} to {high:.{places}f} A^2"
          f"   (width {band.width:.3f}, {width_percent:.2f}% of CCS)")
    print(f"   coverage            nominal {band.nominal_coverage:.0%},"
          f" GUARANTEED {band.guaranteed_coverage:.0%}   [jackknife+ proves 1-2*alpha]")
    print(f"   informative?        {band.interval_is_informative}"
          f"   tail ratio {band.tail_ratio:.2f}")
    wrap("tail ratio: the largest leave-one-out residual over the one that sets the interval"
         " width, so at 1.00 a single ion would be setting it, and"
         f" {'one is' if band.driven_by_one_ion else 'none is'}.", indent=6)
    print(f"   The model holds more digits than these. {places} decimal(s) is what an interval"
          f" {band.width:.2f} wide")
    print("   supports, and it is what the API serves; the rest would be arithmetic rather than")
    print("   evidence. The ORIGINAL above is never rounded.")
    print()

    # WHY A CORRECTION SMALLER THAN ITS OWN INTERVAL IS WORTH APPLYING. The obvious
    # objection to this whole platform, answered on the screen that provokes it rather
    # than in a document nobody opens.
    wrap(f"The shift is {abs(shift):.3f}% and the interval is {width_percent:.2f}% wide, so the"
         " correction is smaller than the uncertainty around it. It is still worth applying,"
         " because the two measure different things. The shift is a SYSTEMATIC offset between"
         " two platforms - the same direction for every ion in this stratum - and leaving it in"
         " biases every comparison anybody makes with these numbers. The interval is the SPREAD"
         " of individual ions around that offset, and correcting the offset does not reduce it."
         " A bias you know about does not average out over many measurements. Scatter does.")
    print()
    print(f"   grade               {result.confidence.grade.value}")
    for demotion in result.confidence.demotions:
        wrap(f"{demotion.rule}: {demotion.detail}", indent=6)
    for note in result.confidence.not_checked:
        wrap(f"NOT CHECKED - {note}", indent=6)
    print()
    explain_the_grade(model, result)

    # --- 7. WHAT IT MAY BE CALLED -----------------------------------------------------------
    head(7, "SCOPE: what this number is, and what it is not")
    wrap(result.scope.caveat())
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
