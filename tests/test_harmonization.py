"""M4 on the real corpus: what it corrects, what it refuses, and what it will not claim.

This runs on the steroid data rather than on fixtures, because M4's whole justification
is that it was not built on synthetic bias. Where a fixture is used it is because the
case being tested does not occur in the real corpus, and each one says so.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest

from wmxccs.contracts import CorrectionBasis
from wmxccs.harmonization import (
    JackknifePlusInterval,
    outward,
    to_places,
    PRIMARY_REFERENCE,
    corpus_digest,
    fit_harmonization,
    harmonize,
    jackknife_plus,
    measure_coverage,
)
from wmxccs.loader import load_measurements_file
from wmxccs.matching import build_matched_ions
from wmxccs.readiness import DataMaturity, smallest_informative_calibration_set
from wmxccs.scope import Claim, ComparisonScope, ScopeExceededError, assert_may_be_quoted_as
from wmxccs.statistics import compare_platforms

SEED = Path(__file__).resolve().parents[1] / "data" / "seed" / "steroid_jasms2022.csv"


@pytest.fixture(scope="module")
def report():
    return load_measurements_file(SEED)


@pytest.fixture(scope="module")
def model(report):
    return fit_harmonization(compare_platforms(build_matched_ions(report.cleared)))


# --- what was fitted ----------------------------------------------------------------------


def test_eighteen_strata_are_fitted_and_only_nine_are_ever_applied(model):
    """Six platform pairs by three adducts, and only the ones anchored on the primary apply.

    A correction between two CALIBRATED platforms anchors nothing: both sides already rest
    on somebody else's reference values. Those nine are computed and published and never
    applied to a measurement, and there is no transitive composition to route around it.
    """
    assert len(model.corrections) == 18
    assert len(model.applied) == 9
    assert len(model.diagnostics) == 9
    assert {c.reference_platform for c in model.applied} == {PRIMARY_REFERENCE}
    assert PRIMARY_REFERENCE not in {c.reference_platform for c in model.diagnostics}


def test_no_stratum_pools_two_calibration_groups(model):
    """One correction per platform pair per adduct, never one across adducts."""
    keys = [(c.reference_platform, c.other_platform, c.stratum.other_group) for c in model.corrections]
    assert len(keys) == len(set(keys)), "a stratum appears twice, so something was pooled"
    adducts = Counter(c.stratum.points[0].ion.key.adduct for c in model.applied)
    assert dict(adducts) == {"[M+H]+": 3, "[M+Na]+": 3, "[M-H]-": 3}


def test_the_headline_basis_is_derived_from_the_slope_interval_not_chosen(model):
    """ROBUST_SLOPE exactly where the rank interval excludes 1, MEDIAN everywhere else.

    Not a preference: where the slope cannot be told from 1, a slope-derived correction and
    a constant offset describe the same data and the offset is the one that does not
    extrapolate.
    """
    for correction in model.corrections:
        expected = (
            CorrectionBasis.ROBUST_SLOPE
            if correction.robust is not None and correction.robust.slope_distinguishable_from_unity
            else CorrectionBasis.MEDIAN
        )
        assert correction.basis is expected
    assert model.basis_counts == {"robust_slope": 10, "median": 8}


def test_all_three_corrections_are_always_reported(model):
    """Both the slope-derived and the median-derived value travel with every estimate.

    They diverge exactly when a few ions are levering the fit, which is the case a reader
    most needs to see, so neither may be dropped for being the one not chosen.
    """
    for correction in model.applied:
        probe = correction.typical_other_ccs
        assert correction.slope_derived(probe) is not None
        assert correction.median_derived(probe) is not None
        assert correction.robust_derived(probe) is not None


# --- the interval -----------------------------------------------------------------------------


def test_every_applied_interval_is_informative_rather_than_the_full_observed_range(model):
    """The reason jackknife+ was chosen over a split.

    A split would have had to divide strata of 23 to 41 into a fit half and a calibration
    half, and a 90 per cent conformal interval is the FULL OBSERVED RANGE below 19
    calibration points. Jackknife+ uses every ion for both, so n_cal = n = 23, 29 or 41,
    all above the informative floor.
    """
    floor = smallest_informative_calibration_set(0.10)
    for correction in model.applied:
        band = correction.interval_at_a_typical_ion
        assert band is not None
        assert correction.n >= floor, f"n={correction.n} is below the informative floor {floor}"
        assert band.interval_is_informative is True
        assert band.readiness_refusal is None


def test_the_interval_states_what_it_guarantees_not_only_what_it_is_nominally(model):
    """Jackknife+ proves 1-2alpha. A 90 per cent interval is guaranteed at 80."""
    for correction in model.applied:
        band = correction.interval_at_a_typical_ion
        assert band.nominal_coverage == pytest.approx(0.90)
        assert band.guaranteed_coverage == pytest.approx(0.80)
        assert band.guaranteed_coverage < band.nominal_coverage


def test_the_interval_depends_on_the_value_it_is_computed_for(model):
    """It is a prediction interval, not a band. A stratum-wide band would be the spread.

    The first implementation made exactly that mistake and produced intervals of 170 to
    217 square angstrom - the whole range of steroid cross sections - for every stratum.
    """
    correction = model.applied[0]
    small = correction.interval_for(160.0)
    large = correction.interval_for(260.0)
    assert small.low != large.low
    assert small.high != large.high
    assert small.at_value == 160.0 and large.at_value == 260.0
    assert small.high < large.low, "an interval must move with the value it describes"


def test_an_interval_width_is_quoted_at_a_named_value(model):
    """Because a jackknife+ interval has no single width, any width must name its query."""
    for correction in model.applied:
        band = correction.interval_at_a_typical_ion
        assert band.at_value == pytest.approx(correction.typical_other_ccs)
        assert 0 < band.width


def test_a_stratum_too_small_to_leave_one_out_is_refused(model):
    """Fewer compounds than a refit needs is a refusal, not a smaller interval."""
    correction = model.applied[0]
    result = jackknife_plus(correction.stratum.points[:3])
    assert isinstance(result, str)
    assert "leave one out" in result


# --- coverage, and what it is not ----------------------------------------------------------------


def test_leave_one_out_coverage_meets_the_guarantee_on_every_applied_stratum(model):
    for correction in model.applied:
        check = measure_coverage(correction.stratum)
        assert not isinstance(check, str), check
        assert check.meets_guarantee


def test_coverage_is_arithmetically_pinned_and_says_so(model):
    """The finding that kills coverage as validation, asserted rather than described.

    The interval is built from the quantile of the residual set, so the fraction inside it
    is forced to the quantile index over n. Measured on all nine applied strata the
    observed rate equals the pinned rate to the last digit, which means it could not have
    come out otherwise and is therefore not evidence about the fit.
    """
    for correction in model.applied:
        check = measure_coverage(correction.stratum)
        assert check.is_pinned, "coverage differed from its forced value, which would be worth reading"
        assert check.is_evidence_of_calibration is False


def test_coverage_still_falsifies_a_grossly_wrong_interval(model):
    """It has no teeth against miscalibration and real teeth against nonsense.

    Shrinking every residual must take coverage below the guarantee. If it did not, the
    check would be incapable of failing and should not be reported at all.
    """
    from wmxccs.harmonization import LeaveOneOutFits, _compound_of

    correction = model.applied[0]
    points = list(correction.stratum.points)
    covered = 0
    for held in points:
        rest = [p for p in points if _compound_of(p) != _compound_of(held)]
        fits = jackknife_plus(rest)
        shrunk = LeaveOneOutFits(
            slopes=fits.slopes,
            intercepts=fits.intercepts,
            residuals=tuple(r * 0.25 for r in fits.residuals),
            groups=fits.groups,
            alpha=fits.alpha,
        )
        band = shrunk.interval_for(held.other_ccs)
        covered += band.low <= held.reference_ccs <= band.high
    assert covered / len(points) < 0.80, "a quarter-width interval must fail the guarantee"


# --- the scope, structurally ----------------------------------------------------------------------


def test_every_correction_carries_a_within_study_scope(model):
    for correction in model.corrections:
        assert correction.scope.scope is ComparisonScope.WITHIN_STUDY
        assert correction.scope.studies == ("doi:10.1021/jasms.2c00196",)
        assert "NOT interlaboratory reproducibility" in correction.scope.caveat()


def test_the_model_refuses_to_be_published_as_an_interlaboratory_figure(model):
    from wmxccs.harmonization import assert_may_be_published

    assert_may_be_published(model, Claim.PLATFORM_DIFFERENCE_WITHIN_A_STUDY)
    with pytest.raises(ScopeExceededError):
        assert_may_be_published(model, Claim.INTERLABORATORY_REPRODUCIBILITY)
    with pytest.raises(ScopeExceededError):
        assert_may_be_published(model, Claim.PLATFORM_DIFFERENCE_ACROSS_STUDIES)


def test_the_maturity_is_provisional_and_validated_is_not_reachable(model):
    """Validation needs data the model was not fitted on, and one study has none."""
    assert model.maturity.data_maturity is DataMaturity.PROVISIONAL
    assert model.maturity.matched_ion_count == 142


# --- harmonizing a measurement --------------------------------------------------------------------


def test_the_original_is_returned_untouched_and_by_identity(report, model):
    """CLAUDE.md constraint 3. Not a copy that happens to be equal: the same object."""
    corrected = 0
    for record in report.cleared:
        result = harmonize(model, record)
        assert result.original is record
        assert result.original_ccs == record.ccs
        corrected += result.was_corrected
    assert corrected == 417


def test_a_measurement_already_on_the_primary_platform_is_not_corrected(report, model):
    """There is nothing to refer it to, and that is the right answer rather than a failure."""
    primary = [r for r in report.cleared if str(r.dtims_method or "") == "stepped_field"]
    assert primary, "the corpus must hold primary records or this test proves nothing"
    for record in primary:
        result = harmonize(model, record)
        assert result.was_corrected is False
        assert any("already on the primary" in reason for reason in result.refusals)


def test_a_corrected_value_carries_all_three_bases_an_interval_and_a_grade(report, model):
    corrected = next(h for h in (harmonize(model, r) for r in report.cleared) if h.was_corrected)
    assert corrected.harmonized_ccs is not None
    assert corrected.slope_derived_ccs is not None
    assert corrected.median_derived_ccs is not None
    assert corrected.robust_derived_ccs is not None
    assert corrected.interval is not None
    assert corrected.confidence is not None
    assert corrected.scope is not None
    assert corrected.basis in (CorrectionBasis.ROBUST_SLOPE, CorrectionBasis.MEDIAN)


def test_the_harmonized_value_is_the_chosen_basis_and_not_a_fourth_number(report, model):
    for record in report.cleared:
        result = harmonize(model, record)
        if not result.was_corrected:
            continue
        expected = (
            result.robust_derived_ccs
            if result.basis is CorrectionBasis.ROBUST_SLOPE
            else result.median_derived_ccs
        )
        assert result.harmonized_ccs == pytest.approx(expected, abs=1e-12)


def test_a_correction_moves_the_value_by_about_the_stratum_offset(report, model):
    """A sanity check with teeth: the correction must be small and in the right direction.

    A harmonized value differing from its original by more than a few per cent would mean
    the inversion is wrong, which no unit test on synthetic points would catch.
    """
    for record in report.cleared:
        result = harmonize(model, record)
        if not result.was_corrected:
            continue
        shift = 100.0 * (result.harmonized_ccs - result.original_ccs) / result.original_ccs
        assert abs(shift) < 5.0, f"a {shift:.2f}% correction is not a harmonization"


# --- the two guards the real corpus cannot exercise -------------------------------------------
#
# Both of these were caught by the mutation sweep rather than by a failing test, and both
# are the shape LIMITATIONS 4.5 names: a guard whose case does not occur in the data it
# was written against, so every test passes whether it works or not.


class _Key:
    def __init__(self, analyte):
        self.analyte = analyte


class _Ion:
    def __init__(self, analyte):
        self.key = _Key(analyte)


class _Point:
    """The attributes jackknife_plus reads, and nothing else."""

    def __init__(self, analyte, reference_ccs, other_ccs):
        self.ion = _Ion(analyte)
        self.reference_ccs = reference_ccs
        self.other_ccs = other_ccs


def _stratum_where_one_compound_has_two_conformers():
    """Eight compounds, one of which contributes TWO points, as a pair of conformers would.

    This case does not occur in the steroid corpus - there every compound contributes
    exactly one ion per stratum, which is why leave-one-compound-out and leave-one-ion-out
    are indistinguishable on it and why the grouping needs a fixture of its own.
    """
    points = []
    for i in range(8):
        x = 150.0 + 10.0 * i
        points.append(_Point(f"compound-{i}", x, 2.0 + 1.03 * x))
    # the second conformer of compound-0: a different value, the SAME compound, and
    # deliberately far off the line so that leaking it into a fit is visible
    points.append(_Point("compound-0", 152.0, 2.0 + 1.03 * 152.0 + 25.0))
    return points


def test_the_leave_one_out_groups_are_compounds_so_a_conformer_cannot_calibrate_its_own_fit():
    """Leaving out an ION would leave its sibling conformer in the fit it is calibrating.

    With compound grouping, holding out compound-0 removes BOTH of its conformers, so the
    residual of the off-line one is measured against a fit that never saw it. Grouping by
    ion instead leaves the sibling in, the fit is pulled towards it, and the residual - and
    therefore the interval - comes out smaller than it should.
    """
    points = _stratum_where_one_compound_has_two_conformers()
    fits = jackknife_plus(points)
    assert not isinstance(fits, str), fits

    # NINE POINTS, EIGHT REFITS. That is the whole assertion: one refit per COMPOUND, one
    # residual per POINT. Grouping by ion would give nine refits, and nine is what the
    # mutation that swaps the grouping produces.
    assert fits.groups == 8, "one refit per compound, not one per ion"
    assert len(fits.residuals) == 9, "one residual per point, including both conformers"
    assert fits.refits == 9, "both conformers carry a (refit, residual) pair"

    # And the off-line conformer's residual is measured against a fit that never saw its
    # compound at all, so it is the full offset rather than a shrunken one.
    assert max(fits.residuals) > 10.0


def test_a_measurement_is_corrected_on_its_own_adduct_and_not_another(report, model):
    """Matching on the platform alone would pool two calibration groups.

    Every adduct's correction is small, so a wrongly-matched one still produces a plausible
    number - which is exactly why this has to be asserted on the group rather than on the
    size of the shift.
    """
    checked = 0
    for record in report.cleared:
        result = harmonize(model, record)
        if not result.was_corrected:
            continue
        checked += 1
        assert result.correction.stratum.other_group == str(record.calibration_group), (
            "a correction was applied from a different calibration group"
        )
        assert result.correction.stratum.points[0].ion.key.adduct == record.adduct
    assert checked == 417


def test_two_adducts_of_one_platform_get_different_corrections(model):
    """The consequence, stated as a fact about the model rather than about one record."""
    twims = [c for c in model.applied if c.other_platform == "TWIMS"]
    assert len(twims) == 3
    offsets = {c.stratum.points[0].ion.key.adduct: c.median_offset_percent for c in twims}
    assert len(set(offsets.values())) == 3, "three adducts must give three different offsets"
    # and they differ by more than rounding: -0.44, -0.08, -0.76
    assert max(offsets.values()) - min(offsets.values()) > 0.5


# --- the digest internals ------------------------------------------------------------------
#
# These are function-level on purpose. The end-to-end versioning tests change a measurement,
# which moves the records, the fit and the residuals together - so they cannot establish
# which of those the digest actually depends on. All three properties below survived the
# mutation sweep until they were tested directly.


def test_the_digest_does_not_depend_on_the_order_the_lines_arrive_in():
    """Sorted before hashing, or one model hashes two ways on two runs.

    Dict and set iteration order is stable within a run and not guaranteed across changes to
    the code that builds them. A digest that moved with it would report a model change on
    every refactor.
    """
    from wmxccs.harmonization import _digest

    lines = ["gamma|3", "alpha|1", "beta|2"]
    assert _digest(lines) == _digest(list(reversed(lines)))
    assert _digest(lines) == _digest(sorted(lines))
    # and it still depends on the CONTENT
    assert _digest(lines) != _digest(lines + ["delta|4"])


def test_the_corpus_digest_does_not_depend_on_record_order(report):
    records = list(report.cleared)
    assert corpus_digest(records) == corpus_digest(list(reversed(records)))
    assert corpus_digest(records) != corpus_digest(records[:-1])


def test_floats_enter_the_digest_exactly_and_are_not_rounded():
    """A rounded digest calls two different fits identical.

    repr round-trips exactly in Python 3. Formatting to a few decimal places would make any
    change below that precision invisible - and a slope differing in the ninth decimal is a
    different fit, however little it matters to an answer.
    """
    from wmxccs.harmonization import _canonical

    assert _canonical(1.0) != _canonical(1.000000001)
    assert _canonical(0.1 + 0.2) != _canonical(0.3)
    assert float(_canonical(1.2345678901234567)) == 1.2345678901234567


def test_a_change_far_below_three_decimal_places_still_moves_the_corpus_digest(report):
    """The end-to-end test changes a value by 0.001, which survives rounding. This does not."""
    records = list(report.cleared)
    nudged = [records[0].model_copy(update={"ccs": records[0].ccs + 1e-9})] + records[1:]
    assert corpus_digest(nudged) != corpus_digest(records)


def test_the_parameters_digest_covers_the_residuals_that_set_every_interval(model):
    """Two models with identical fits and different residuals are two models.

    The residuals decide every interval, so a digest omitting them would call a model with
    twice the interval width identical to this one. Nothing in the end-to-end tests can see
    that, because changing the data changes the slopes too.
    """
    import dataclasses

    from wmxccs.harmonization import parameters_digest

    corrections = list(model.applied)
    original = parameters_digest(corrections, model.alpha)

    first = corrections[0]
    widened = dataclasses.replace(
        first, loo=dataclasses.replace(first.loo, residuals=tuple(r * 2 for r in first.loo.residuals))
    )
    altered = [widened] + corrections[1:]

    assert parameters_digest(altered, model.alpha) != original


def test_the_parameters_digest_covers_alpha(model):
    """The coverage level is part of what an answer is, so it is part of the digest."""
    from wmxccs.harmonization import parameters_digest

    corrections = list(model.applied)
    assert parameters_digest(corrections, 0.10) != parameters_digest(corrections, 0.05)


def test_the_parameters_digest_does_not_depend_on_stratum_order(model):
    from wmxccs.harmonization import parameters_digest

    corrections = list(model.applied)
    assert parameters_digest(corrections, model.alpha) == parameters_digest(
        list(reversed(corrections)), model.alpha
    )


# --- the precision a served number is allowed to claim ------------------------------------------
#
# THESE ARE FUNCTION-LEVEL ON PURPOSE. An end-to-end test that submits a measurement and
# reads the wire exercises this rule but cannot isolate WHICH input it read: on most of the
# seed corpus, reading the full width instead of the half width lands on the same answer, so
# an end-to-end assertion passes either way. Three mutations survived a sweep on exactly that
# before these were written - the same lesson the versioning digest taught in M5.


def an_interval(low: float, high: float) -> JackknifePlusInterval:
    """A minimal interval. Only `low` and `high` matter to the rule under test."""
    return JackknifePlusInterval(
        low=low,
        high=high,
        at_value=(low + high) / 2.0,
        nominal_coverage=0.9,
        guaranteed_coverage=0.8,
        groups_left_out=1,
        refits=20,
        quantile_index_low=1,
        quantile_index_high=20,
        largest_residual=1.0,
        quantile_residual=1.0,
        interval_is_informative=True,
    )


def test_the_precision_is_read_off_the_half_width_not_the_whole_width():
    """A width of 1.09 is a half-width of 0.54, and those support different precisions.

    The interval brackets the value on BOTH sides, so what the number can claim is set by
    how far it may be from the truth in one direction - the half-width. Reading the full
    width understates the precision by a factor of two and lands on the same answer for
    most of the seed corpus, which is why only a function-level test can tell them apart.
    """
    assert an_interval(100.0, 101.09).decimals_supported == 2
    # and the same figure read as a half-width gives the other answer, which is the point
    assert an_interval(100.0, 102.18).decimals_supported == 1


def test_a_wider_interval_supports_fewer_decimals_than_a_tighter_one():
    """The direction of the rule, asserted as a direction rather than at one point."""
    tight = an_interval(100.0, 100.02).decimals_supported      # half 0.01
    middling = an_interval(100.0, 101.2).decimals_supported     # half 0.6
    wide = an_interval(100.0, 111.0).decimals_supported         # half 5.5
    assert tight > middling > wide
    assert (tight, middling, wide) == (3, 2, 1)


def test_the_precision_never_exceeds_the_three_decimals_the_corpus_itself_carries():
    """The upper clamp. No seed interval reaches it, so nothing else would notice it going.

    The seed measurements carry at most three decimal places. An interval tight enough to
    justify six would be claiming to know the answer better than any measurement behind it
    was ever recorded, which is arithmetic outrunning its evidence.
    """
    assert an_interval(100.0, 100.000002).decimals_supported == 3
    assert an_interval(100.0, 100.0).decimals_supported == 3


def test_the_precision_never_falls_below_one_decimal():
    """The lower clamp, and the other direction: a huge interval still says something."""
    assert an_interval(100.0, 400.0).decimals_supported == 1
    assert an_interval(100.0, 100000.0).decimals_supported == 1


def test_an_interval_rounds_outward_so_it_is_never_narrower_than_it_was():
    # Chosen so that NEAREST would move BOTH bounds inward: it gives (1.24, 6.78), an
    # interval 0.02 narrower than the one that was computed. Outward gives (1.23, 6.79).
    low, high = outward(1.2371, 6.7849, 2)
    assert (low, high) == (1.23, 6.79)
    assert low <= 1.2371 and high >= 6.7849
    assert low < round(1.2371, 2) and high > round(6.7849, 2)
    assert high - low > round(6.7849, 2) - round(1.2371, 2)


def test_a_value_rounds_to_nearest_because_it_sits_inside_its_interval():
    assert to_places(166.56725251726084, 1) == 166.6
    assert to_places(166.56725251726084, 2) == 166.57
