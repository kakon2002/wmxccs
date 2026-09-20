"""Confidence grading: a scheme of demotions, asserted in both directions.

THE PROPERTY THE WHOLE DESIGN RESTS ON IS THAT A GRADE ONLY EVER FALLS, and the
grade is the WORST demotion found rather than a score, an average or a count. So
the first two tests here are a matched pair: one severe rule beside several mild
ones must still grade at the severe one, and two mild rules must NOT add up to a
severe one. A suite holding only the first would pass a module that returned
"unsupported" for everything; a suite holding only the second would pass one that
never escalated at all.

Every rule below is asserted twice, once where it must fire and once where it must
not. A test that only reads the happy path kills nothing: a rule that always
demoted and a rule that never demoted would each satisfy half of this file.

THE SINGLE MOST IMPORTANT TEST HERE is the one asserting that the leverage
diagnostic does not remove the outliers it measures. "Refit without the outliers
to see how much they matter" is one edit away from "refit without the outliers",
and nothing but a test stands between the two. That test therefore checks the
stratum's points and its reported Deming slope before and after the call, not
only the number the call returns.

The numbers the benchmark stratum carries are read off fixtures.benchmark_corpus,
whose biases are declared in code rather than measured: twenty-four matched ions,
a centre near the injected +2.0 per cent, and exactly three ions carrying the
injected seven per cent that do not transfer. No figure in this file describes any
instrument.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import pytest

from wmxccs import fixtures
from wmxccs.grading import (
    EXTRAPOLATION_MARGIN,
    LEVERAGE_LIMIT_PERCENT,
    Confidence,
    ConfidenceGrade,
    Demotion,
    correction_driven_by_outliers,
    flagged_as_an_outlier,
    grade_correction,
    grade_rules,
    outside_calibration_range,
    slope_leverage_percent,
    the_outlier_rule_can_run_for,
    thinly_populated,
    unverified_source_behind_it,
)
from wmxccs.matching import build_matched_ions
from wmxccs.readiness import TARGET_MATCHED_IONS
from wmxccs.reuse import ReuseStatus
from wmxccs.statistics import (
    MIN_POINTS_FOR_DEMING,
    MIN_POINTS_FOR_LIMITS_OF_AGREEMENT,
    compare_platforms,
    deming_slope_intercept,
)

# Platform labels as statistics.py spells them.
DTIMS = "DTIMS/stepped_field"
TWIMS = "TWIMS"
TIMS = "TIMS"

# The benchmark TWIMS stratum, read off fixtures.benchmark_corpus. Asserted in
# test_the_benchmark_stratum_is_the_stratum_these_tests_believe_it_is below, so
# that every later assertion rests on a shape that was checked rather than assumed.
BENCHMARK_N = 24
BENCHMARK_OUTLIERS = 3
BENCHMARK_CENTRE_PERCENT = 2.0
BENCHMARK_LOW_CCS = 150.0
BENCHMARK_HIGH_CCS = 380.0
BENCHMARK_LEVERAGE_PERCENT = 0.61
BENCHMARK_DEMING_SLOPE = 1.03098

# Just past the top of the fitted range but inside the extrapolation margin of
# 10 per cent of the span, which is 23 square angstrom here.
JUST_OUTSIDE_CCS = BENCHMARK_HIGH_CCS + 1.0
# An antibody-sized ion. Nothing in a 150 to 380 fit resembles it.
FAR_OUTSIDE_CCS = 7000.0


# --- building strata ------------------------------------------------------------------------


def stratum_of(records):
    """The one stratum a corpus of one platform pair yields, built as the pipeline builds it."""
    report = compare_platforms(build_matched_ions(records))
    assert len(report.pairs) == 1, f"expected one platform pair, got {len(report.pairs)}"
    assert len(report.pairs[0].strata) == 1, f"expected one stratum, got {len(report.pairs[0].strata)}"
    return report.pairs[0].strata[0]


def paired_corpus(count: int, *, tag: str, bias_percent: float = 2.0):
    """`count` ions on stepped-field DTIMS and TWIMS, biased by `bias_percent`.

    Every ion gets its own molecule. Two ions sharing one would share a
    matched-ion key, merge into a single matched set, and quietly shrink the
    stratum below the size the test asked for.
    """
    records = []
    for index in range(count):
        analyte = fixtures.small_molecule(fixtures.synthetic_inchikey(f"{tag}{index:03d}"))
        reference_ccs = 150.0 + 10.0 * index
        records.append(
            fixtures.measurement(
                analyte=analyte,
                platform="dtims",
                ccs=reference_ccs,
                source="fixture DTIMS",
                doi=fixtures.DOI_A,
            )
        )
        records.append(
            fixtures.measurement(
                analyte=analyte,
                platform="twims",
                # A little scatter so the points are not collinear, and far too
                # little for any ion to be flagged against its own centre.
                ccs=round(reference_ccs * (1 + bias_percent / 100.0) + (index % 3) * 0.1, 4),
                source="fixture TWIMS",
                doi=fixtures.DOI_B,
            )
        )
    return tuple(records)


def stratum_of_size(count: int, *, tag: str):
    stratum = stratum_of(paired_corpus(count, tag=tag))
    assert stratum.n == count, f"asked for {count} matched ions, built {stratum.n}"
    return stratum


def with_first_point_recorded_as(stratum, record):
    """A copy of `stratum` whose first point carries `record` on the comparison side.

    Built here rather than through the pipeline on purpose. The licence gate keeps
    a record nobody may use out of a matched set before any stratum is formed, so a
    stratum holding one cannot be produced by compare_platforms at all. The rule
    exists for exactly that case: a grade that trusted the gate to have run would
    be trusting something it cannot see.
    """
    first, *rest = stratum.points
    return replace(stratum, points=(replace(first, other=record), *rest))


@dataclass(frozen=True)
class ARecordWhoseReuseStatusCannotBeRead:
    """A stand-in for a record carrying something that is not a reuse status at all.

    No CCSMeasurement can be built in this state, which is the point: the rule
    reads reuse_status off whatever it is handed, and what it does with something
    unreadable has to be asserted against something unreadable.
    """

    ccs: float
    source: str
    reuse_status: object = None


# --- the corpora the assertions below rest on -------------------------------------------------


@pytest.fixture(scope="module")
def benchmark():
    return compare_platforms(build_matched_ions(fixtures.benchmark_corpus()))


@pytest.fixture(scope="module")
def twims_stratum(benchmark):
    """TWIMS against stepped-field DTIMS: 24 ions, centre near +2%, three flagged."""
    found = [
        stratum
        for pair in benchmark.pairs
        if (pair.reference_platform, pair.other_platform) == (DTIMS, TWIMS)
        for stratum in pair.strata
    ]
    assert len(found) == 1, f"expected one TWIMS against DTIMS stratum, got {len(found)}"
    return found[0]


@pytest.fixture(scope="module")
def tims_stratum(benchmark):
    """TIMS against stepped-field DTIMS: 24 ions and nothing flagged, so leverage is unmeasurable."""
    found = [
        stratum
        for pair in benchmark.pairs
        if (pair.reference_platform, pair.other_platform) == (DTIMS, TIMS)
        for stratum in pair.strata
    ]
    assert len(found) == 1, f"expected one TIMS against DTIMS stratum, got {len(found)}"
    return found[0]


@pytest.fixture(scope="module")
def outlying_point(twims_stratum):
    """One of the three ions the stratum itself reports as not transferring."""
    assert twims_stratum.outliers, "nothing was flagged, so the outlier tests would prove nothing"
    return twims_stratum.outliers[0]


@pytest.fixture(scope="module")
def ordinary_point(twims_stratum):
    """An ion sitting at its stratum's centre. The control for every outlier assertion."""
    flagged = {id(point) for point in twims_stratum.outliers}
    ordinary = [point for point in twims_stratum.points if id(point) not in flagged]
    assert len(ordinary) == BENCHMARK_N - BENCHMARK_OUTLIERS
    return ordinary[0]


def test_the_benchmark_stratum_is_the_stratum_these_tests_believe_it_is(twims_stratum):
    # Every figure below is read off this stratum, so its shape is pinned here
    # once rather than assumed in twenty places.
    assert twims_stratum.n == BENCHMARK_N
    assert len(twims_stratum.outliers) == BENCHMARK_OUTLIERS
    assert twims_stratum.centre_percent == pytest.approx(BENCHMARK_CENTRE_PERCENT, abs=0.05)
    values = [point.reference_ccs for point in twims_stratum.points]
    assert min(values) == pytest.approx(BENCHMARK_LOW_CCS, abs=1e-9)
    assert max(values) == pytest.approx(BENCHMARK_HIGH_CCS, abs=1e-9)


# --- 1. worst wins, and a count of mild concerns is not a severe one ---------------------------


def test_one_severe_rule_beside_several_mild_ones_still_grades_at_the_severe_one(
    twims_stratum, outlying_point
):
    """A flagged ion in a thinly populated, outlier-levered stratum grades unsupported.

    Three rules fire at once here: the ion is flagged (unsupported), the
    correction is levered by ions that do not transfer (weak), and 24 matched ions
    is below the target (qualified). If the grade were an average or a count, the
    two mild concerns would dilute the severe one into something usable.
    """
    confidence = grade_correction(
        outlying_point.reference_ccs, outlying_point.ion.key, twims_stratum
    )
    grades = {demotion.grade for demotion in confidence.demotions}
    assert grades == {
        ConfidenceGrade.UNSUPPORTED,
        ConfidenceGrade.WEAK,
        ConfidenceGrade.QUALIFIED,
    }, f"this test needs one severe rule and several mild ones, got {confidence.summary()}"
    assert confidence.grade is ConfidenceGrade.UNSUPPORTED
    assert confidence.usable is False


def test_two_mild_rules_do_not_add_up_to_a_severe_grade(twims_stratum, ordinary_point):
    """An extrapolated, levered, thinly populated correction is weak and no worse.

    Three demotions, none of them severe. A scheme that scored or counted would
    have reached unsupported somewhere along the way; this one cannot, because the
    grade is the worst single demotion and the worst of these is weak.
    """
    confidence = grade_correction(JUST_OUTSIDE_CCS, ordinary_point.ion.key, twims_stratum)
    assert len(confidence.demotions) >= 2, f"this test needs several mild rules, got {confidence.summary()}"
    assert ConfidenceGrade.UNSUPPORTED not in {demotion.grade for demotion in confidence.demotions}
    assert confidence.grade is ConfidenceGrade.WEAK
    assert confidence.usable is True


def test_the_grade_is_the_worst_demotion_and_never_softened_by_the_rest(
    twims_stratum, ordinary_point, outlying_point
):
    # The same assertion made twice over, through two different severe rules, so
    # that a change to one scenario cannot quietly retire the property.
    worst = ConfidenceGrade.UNSUPPORTED
    for ccs, key in (
        (FAR_OUTSIDE_CCS, ordinary_point.ion.key),
        (outlying_point.reference_ccs, outlying_point.ion.key),
    ):
        confidence = grade_correction(ccs, key, twims_stratum)
        assert confidence.grade is worst
        assert len(confidence.demotions) > 1, "a single demotion would not test the reduction"


# --- 2. outside the calibration range ----------------------------------------------------------


def test_an_ion_inside_the_fitted_range_is_not_demoted_for_extrapolation(twims_stratum):
    inside = (BENCHMARK_LOW_CCS + BENCHMARK_HIGH_CCS) / 2
    assert outside_calibration_range(inside, twims_stratum) is None
    assert outside_calibration_range(BENCHMARK_LOW_CCS, twims_stratum) is None
    assert outside_calibration_range(BENCHMARK_HIGH_CCS, twims_stratum) is None


def test_an_ion_just_outside_the_fitted_range_falls_to_weak_and_the_message_names_the_range(
    twims_stratum, ordinary_point
):
    demotion = outside_calibration_range(JUST_OUTSIDE_CCS, twims_stratum)
    assert demotion is not None
    assert demotion.grade is ConfidenceGrade.WEAK
    assert demotion.rule == "outside the calibration range"
    assert "150" in demotion.detail and "380" in demotion.detail
    assert f"{EXTRAPOLATION_MARGIN:.0%}" in demotion.detail
    assert grade_correction(JUST_OUTSIDE_CCS, ordinary_point.ion.key, twims_stratum).usable is True


def test_an_antibody_sized_ion_far_outside_the_fitted_range_falls_to_unsupported(
    twims_stratum, ordinary_point
):
    """7000 square angstrom against a fit over 150 to 380. Not a caveat but a refusal."""
    demotion = outside_calibration_range(FAR_OUTSIDE_CCS, twims_stratum)
    assert demotion is not None
    assert demotion.grade is ConfidenceGrade.UNSUPPORTED
    assert demotion.rule == "far outside the calibration range"
    assert "150" in demotion.detail and "380" in demotion.detail
    assert "7000" in demotion.detail
    confidence = grade_correction(FAR_OUTSIDE_CCS, ordinary_point.ion.key, twims_stratum)
    assert confidence.grade is ConfidenceGrade.UNSUPPORTED
    assert confidence.usable is False


def test_the_margin_is_the_line_between_extrapolating_and_asserting(twims_stratum):
    # The span is 230, so the margin is 23. Either side of 403 the rule must give
    # a different answer, or the margin is not doing anything.
    span = BENCHMARK_HIGH_CCS - BENCHMARK_LOW_CCS
    edge = BENCHMARK_HIGH_CCS + span * EXTRAPOLATION_MARGIN
    assert outside_calibration_range(edge - 0.5, twims_stratum).grade is ConfidenceGrade.WEAK
    assert outside_calibration_range(edge + 0.5, twims_stratum).grade is ConfidenceGrade.UNSUPPORTED


# --- 3. thinly populated calibration groups ------------------------------------------------------


def test_a_stratum_below_the_points_a_slope_needs_is_refused_rather_than_graded():
    stratum = stratum_of_size(MIN_POINTS_FOR_DEMING - 1, tag="slopefew")
    demotion = thinly_populated(stratum)
    assert demotion is not None
    assert demotion.grade is ConfidenceGrade.UNSUPPORTED
    assert demotion.rule == "no correction can be fitted"
    assert str(MIN_POINTS_FOR_DEMING) in demotion.detail
    confidence = grade_correction(stratum.points[0].reference_ccs, stratum.points[0].ion.key, stratum)
    assert confidence.grade is ConfidenceGrade.UNSUPPORTED
    assert confidence.usable is False


def test_a_stratum_below_the_points_limits_of_agreement_need_falls_to_weak():
    stratum = stratum_of_size(MIN_POINTS_FOR_LIMITS_OF_AGREEMENT - 5, tag="loafew")
    demotion = thinly_populated(stratum)
    assert demotion is not None
    assert demotion.grade is ConfidenceGrade.WEAK
    assert demotion.rule == "thinly populated calibration group"
    assert str(MIN_POINTS_FOR_LIMITS_OF_AGREEMENT) in demotion.detail


def test_a_stratum_below_the_target_matched_ions_falls_to_qualified():
    stratum = stratum_of_size(MIN_POINTS_FOR_LIMITS_OF_AGREEMENT + 2, tag="belowtgt")
    demotion = thinly_populated(stratum)
    assert demotion is not None
    assert demotion.grade is ConfidenceGrade.QUALIFIED
    assert demotion.rule == "thinly populated calibration group"
    assert str(TARGET_MATCHED_IONS) in demotion.detail


def test_a_stratum_at_the_target_matched_ions_is_not_demoted_for_population():
    """The other direction. A rule that demoted every stratum would pass the three above."""
    stratum = stratum_of_size(TARGET_MATCHED_IONS, tag="attarget")
    assert stratum.n >= TARGET_MATCHED_IONS
    assert thinly_populated(stratum) is None
    confidence = grade_correction(stratum.points[0].reference_ccs, stratum.points[0].ion.key, stratum)
    # NOT DEMOTED FOR POPULATION, which is the whole of what this test is about. Asserting
    # the grade itself would make it fail for a reason it does not test: since 20 September
    # 2026 an unevaluable rule demotes one notch, and leverage is unmeasurable on a fixture
    # this clean. So pin the RULES, which is the claim, rather than the grade, which is a
    # consequence of every rule at once.
    assert [demotion.rule for demotion in confidence.demotions] == [
        "a rule of the scheme could not be evaluated for this ion"
    ], "a stratum at the target population is demoted for something other than the unevaluable rule"
    assert confidence.grade is ConfidenceGrade.QUALIFIED


def test_the_three_population_rungs_are_ordered_worst_first():
    # A rung that demoted no further than the one below it would leave the
    # thresholds doing nothing, and each test above would still pass alone.
    rungs = [
        thinly_populated(stratum_of_size(count, tag=tag))
        for count, tag in (
            (MIN_POINTS_FOR_DEMING - 1, "rungOne"),
            (MIN_POINTS_FOR_DEMING + 1, "rungTwo"),
            (MIN_POINTS_FOR_LIMITS_OF_AGREEMENT + 1, "rungThree"),
        )
    ]
    assert [rung.grade for rung in rungs] == [
        ConfidenceGrade.UNSUPPORTED,
        ConfidenceGrade.WEAK,
        ConfidenceGrade.QUALIFIED,
    ]


# --- 4. flagged as an outlier in its own stratum ---------------------------------------------


def test_an_ion_its_own_stratum_flagged_grades_unsupported_and_the_message_quotes_both_figures(
    twims_stratum, outlying_point
):
    """The message has to carry the ion's own difference AND the stratum's centre.

    Either number alone is unreadable: seven per cent is only alarming beside a
    centre of two, and a centre of two says nothing about this ion.
    """
    demotion = flagged_as_an_outlier(outlying_point.ion.key, twims_stratum)
    assert demotion is not None
    assert demotion.grade is ConfidenceGrade.UNSUPPORTED
    assert demotion.rule == "flagged as an outlier in its own stratum"
    assert f"{outlying_point.difference_percent:+.2f}%" in demotion.detail
    assert f"{twims_stratum.centre_percent:+.2f}%" in demotion.detail
    confidence = grade_correction(outlying_point.reference_ccs, outlying_point.ion.key, twims_stratum)
    assert confidence.grade is ConfidenceGrade.UNSUPPORTED
    assert confidence.usable is False


def test_an_ordinary_ion_in_the_same_stratum_is_not_flagged(twims_stratum, ordinary_point):
    assert flagged_as_an_outlier(ordinary_point.ion.key, twims_stratum) is None
    confidence = grade_correction(
        ordinary_point.reference_ccs, ordinary_point.ion.key, twims_stratum
    )
    assert "flagged as an outlier in its own stratum" not in {
        demotion.rule for demotion in confidence.demotions
    }
    assert confidence.usable is True


def test_an_ion_that_is_not_in_the_stratum_at_all_is_not_flagged_by_it(twims_stratum):
    # A key nothing in the stratum holds must not match the first outlier by
    # accident, which a comparison that had stopped looking at the key would do.
    assert flagged_as_an_outlier("not a key this stratum holds", twims_stratum) is None


# --- 5. a source that may not be used ---------------------------------------------------------


def test_a_stratum_holding_an_unverified_record_is_unsupported_and_the_message_names_the_source(
    twims_stratum,
):
    offending_source = "a source nobody has read the terms of"
    spoiled = with_first_point_recorded_as(
        twims_stratum,
        fixtures.measurement(
            ccs=twims_stratum.points[0].other_ccs,
            source=offending_source,
            reuse_status=ReuseStatus.UNVERIFIED,
        ),
    )
    demotion = unverified_source_behind_it(spoiled)
    assert demotion is not None
    assert demotion.grade is ConfidenceGrade.UNSUPPORTED
    assert demotion.rule == "a source behind the correction may not be used"
    assert offending_source in demotion.detail
    assert ReuseStatus.UNVERIFIED.value in demotion.detail
    confidence = grade_correction(
        spoiled.points[0].reference_ccs, spoiled.points[0].ion.key, spoiled
    )
    assert confidence.grade is ConfidenceGrade.UNSUPPORTED
    assert confidence.usable is False


def test_a_record_whose_reuse_status_cannot_be_read_at_all_counts_against_the_correction(
    twims_stratum,
):
    """Unreadable is not the same as cleared, and must not be waved through.

    A status that cannot be parsed is the case where the gate has already failed
    in some way nobody predicted, so it is counted against rather than assumed
    harmless.
    """
    offending_source = "a source whose reuse status is not a status"
    spoiled = with_first_point_recorded_as(
        twims_stratum,
        ARecordWhoseReuseStatusCannotBeRead(
            ccs=twims_stratum.points[0].other_ccs, source=offending_source
        ),
    )
    demotion = unverified_source_behind_it(spoiled)
    assert demotion is not None
    assert demotion.grade is ConfidenceGrade.UNSUPPORTED
    assert offending_source in demotion.detail
    assert "unreadable reuse status" in demotion.detail


def test_a_stratum_where_every_record_is_cleared_does_not_trip_the_licence_rule(twims_stratum):
    """The other direction. Without it a rule that always fired would pass the two above."""
    assert unverified_source_behind_it(twims_stratum) is None
    confidence = grade_correction(
        twims_stratum.points[0].reference_ccs, twims_stratum.points[0].ion.key, twims_stratum
    )
    assert "a source behind the correction may not be used" not in {
        demotion.rule for demotion in confidence.demotions
    }


# --- 6. leverage: a diagnostic, and only a diagnostic ------------------------------------------


def test_the_benchmark_strata_leverage_is_the_figure_the_rule_was_written_against(twims_stratum):
    leverage = slope_leverage_percent(twims_stratum)
    assert leverage is not None
    assert leverage == pytest.approx(BENCHMARK_LEVERAGE_PERCENT, abs=0.01)
    assert leverage > LEVERAGE_LIMIT_PERCENT


def test_a_correction_levered_past_the_limit_by_ions_that_do_not_transfer_falls_to_weak(
    twims_stratum,
):
    demotion = correction_driven_by_outliers(twims_stratum)
    assert demotion is not None
    assert demotion.grade is ConfidenceGrade.WEAK
    assert demotion.rule == "the correction is driven by ions that do not transfer"
    assert f"{BENCHMARK_LEVERAGE_PERCENT:.2f}%" in demotion.detail
    assert f"{LEVERAGE_LIMIT_PERCENT:g}%" in demotion.detail
    assert "NOT removed" in demotion.detail


def test_a_stratum_with_no_outliers_has_no_leverage_to_measure_and_trips_nothing(tims_stratum):
    assert tims_stratum.outliers == ()
    assert slope_leverage_percent(tims_stratum) is None
    assert correction_driven_by_outliers(tims_stratum) is None


def test_the_leverage_diagnostic_does_not_remove_the_outliers_it_measures(twims_stratum):
    """THE MOST IMPORTANT TEST IN THIS FILE.

    Measuring how much the outliers move the correction means refitting without
    them, and a refit that was kept rather than thrown away would be the one thing
    this repository refuses hardest. So the stratum is checked before and after:
    the same points, in the same order, the same three flagged ions, and the same
    reported Deming slope. If the diagnostic ever starts writing its refit back,
    the slope here moves off 1.031 towards the injected 1.020 and this fails.
    """
    points_before = [id(point) for point in twims_stratum.points]
    outliers_before = [id(point) for point in twims_stratum.outliers]
    slope_before = twims_stratum.agreement.deming_slope
    intercept_before = twims_stratum.agreement.deming_intercept
    centre_before = twims_stratum.centre_percent

    assert slope_leverage_percent(twims_stratum) == pytest.approx(
        BENCHMARK_LEVERAGE_PERCENT, abs=0.01
    )
    assert correction_driven_by_outliers(twims_stratum) is not None
    grade_correction(twims_stratum.points[0].reference_ccs, twims_stratum.points[0].ion.key, twims_stratum)

    assert [id(point) for point in twims_stratum.points] == points_before
    assert twims_stratum.n == BENCHMARK_N == len(points_before)
    assert [id(point) for point in twims_stratum.outliers] == outliers_before
    assert len(twims_stratum.outliers) == BENCHMARK_OUTLIERS
    assert twims_stratum.agreement.n == BENCHMARK_N
    assert twims_stratum.agreement.deming_slope == slope_before
    assert twims_stratum.agreement.deming_intercept == intercept_before
    assert twims_stratum.centre_percent == centre_before
    # The reported slope is still the levered one, fitted on all 24 ions.
    assert twims_stratum.agreement.deming_slope == pytest.approx(BENCHMARK_DEMING_SLOPE, abs=5e-4)


def test_leverage_is_measured_at_a_typical_ion_rather_than_at_a_cross_section_of_zero(
    twims_stratum,
):
    """The figure is a percentage OF CCS at the middle of the range, not at the intercept.

    Two lines differing by a hair at a real ion can differ by a great deal
    extrapolated back to zero, so a comparison made at the intercept would report a
    correction as levered when nothing an actual ion receives has moved. The refit
    is redone here, independently, and the two evaluation points are compared: the
    figure the module returns must be the one taken at the median ion.
    """
    from statistics import median

    flagged = {id(point) for point in twims_stratum.outliers}
    kept = [point for point in twims_stratum.points if id(point) not in flagged]
    without = deming_slope_intercept(
        [point.reference_ccs for point in kept],
        [point.other_ccs for point in kept],
        twims_stratum.agreement.lambda_used,
    )
    assert without is not None
    slope_without, intercept_without = without
    slope = twims_stratum.agreement.deming_slope
    intercept = twims_stratum.agreement.deming_intercept

    def moved_at(size: float) -> float:
        with_outliers = slope * size + intercept
        without_outliers = slope_without * size + intercept_without
        return 100.0 * abs(with_outliers - without_outliers) / size

    middle = median([point.reference_ccs for point in twims_stratum.points])
    assert middle == pytest.approx(265.0, abs=1e-9)
    assert slope_leverage_percent(twims_stratum) == pytest.approx(moved_at(middle), abs=1e-9)
    # The same two lines compared at a cross section no ion in this stratum has
    # give a materially different answer, which is why the evaluation point is
    # stated rather than left at the intercept.
    assert moved_at(1.0) != pytest.approx(moved_at(middle), abs=0.05)
    assert slope_leverage_percent(twims_stratum) != pytest.approx(moved_at(1.0), abs=0.05)


# --- 7. not checked is not the same as checked and fine ----------------------------------------


def test_a_stratum_with_no_outliers_reports_the_leverage_check_as_not_checked(tims_stratum):
    confidence = grade_correction(
        tims_stratum.points[0].reference_ccs, tims_stratum.points[0].ion.key, tims_stratum
    )
    leverage_notes = [note for note in confidence.not_checked if note.startswith("leverage:")]
    assert len(leverage_notes) == 1, f"expected one leverage note, got {confidence.not_checked}"
    assert "no outliers to exclude" in leverage_notes[0]


def test_a_check_that_could_not_be_run_never_appears_among_the_demotions(tims_stratum):
    confidence = grade_correction(
        tims_stratum.points[0].reference_ccs, tims_stratum.points[0].ion.key, tims_stratum
    )
    assert confidence.not_checked, "this test needs an unevaluable check to be meaningful"
    for note in confidence.not_checked:
        assert note not in {demotion.detail for demotion in confidence.demotions}
    assert "the correction is driven by ions that do not transfer" not in {
        demotion.rule for demotion in confidence.demotions
    }
    rendered = confidence.summary()
    assert "[not checked]" in rendered
    assert f"[{ConfidenceGrade.SUPPORTED}]" not in rendered


def test_a_stratum_whose_leverage_could_be_measured_does_not_report_it_as_unchecked(
    twims_stratum,
):
    """The other direction. A module that noted everything as unchecked would pass the two above."""
    confidence = grade_correction(
        twims_stratum.points[0].reference_ccs, twims_stratum.points[0].ion.key, twims_stratum
    )
    assert not [note for note in confidence.not_checked if note.startswith("leverage:")]
    assert "the correction is driven by ions that do not transfer" in {
        demotion.rule for demotion in confidence.demotions
    }


# --- 7b. a rule that could not be evaluated is named, and it demotes ---------------------------
#
# THE FAILURE THIS SECTION EXISTS FOR was live in every release up to 20 September 2026 and was
# invisible because it had the shape of a passing check. `flagged_as_an_outlier` is a CORPUS
# LOOKUP: it asks whether the submitted ion is among the ions this stratum already recorded as
# not transferring. For an ion that is not in the corpus at all the answer is not "no" - there
# is nothing to look in. It returned None either way, grade_correction read None as "checked
# and fine", and `not_checked` came back empty, which this API uses everywhere as a positive
# claim that every rule ran.
#
# The population that reaches this service is precisely the ions NOT already in its corpus. So
# the only rule in the scheme that speaks about the submitted ion rather than about the stratum
# around it was inert for every real caller, while the response asserted it had been applied.
#
# Every test below is written in both directions on purpose. A module that reported everything
# as unevaluable, or that demoted everything, would satisfy half of them.


def a_key_this_stratum_has_never_seen():
    """A real matched-ion key belonging to a different corpus.

    Taken from another stratum rather than hand-built, so it is a key of exactly the kind
    the code meets in production: a well-formed ion that simply is not here.
    """
    return stratum_of_size(5, tag="foreignion").points[0].ion.key


def test_the_outlier_rule_cannot_run_for_an_ion_the_stratum_has_never_seen(twims_stratum):
    assert the_outlier_rule_can_run_for(a_key_this_stratum_has_never_seen(), twims_stratum) is False


def test_the_outlier_rule_can_run_for_an_ion_the_stratum_holds(twims_stratum):
    """The other direction. A predicate that returned False always would pass the test above."""
    assert the_outlier_rule_can_run_for(twims_stratum.points[0].ion.key, twims_stratum) is True


def test_a_novel_ion_has_the_outlier_check_reported_as_not_checked(twims_stratum):
    confidence = grade_correction(
        twims_stratum.points[0].reference_ccs, a_key_this_stratum_has_never_seen(), twims_stratum
    )
    notes = [note for note in confidence.not_checked if note.startswith("outlier check:")]
    assert len(notes) == 1, f"expected the outlier lookup to be named, got {confidence.not_checked}"
    assert "not in the stratum" in notes[0]
    # and it is NOT reported as a rule that fired, which would be the opposite lie
    assert "flagged as an outlier in its own stratum" not in {d.rule for d in confidence.demotions}


def test_an_ion_the_stratum_holds_has_every_rule_evaluated(twims_stratum):
    """The other direction, and the claim an empty not_checked is making.

    This stratum measures its own leverage and holds this ion, so nothing is unevaluable. If
    this ever reports a note, an empty `not_checked` has stopped meaning what the rest of this
    file and the API both read it as meaning.
    """
    flagged = {id(point) for point in twims_stratum.outliers}
    clean = [point for point in twims_stratum.points if id(point) not in flagged][0]
    confidence = grade_correction(clean.reference_ccs, clean.ion.key, twims_stratum)
    assert confidence.not_checked == ()
    assert "a rule of the scheme could not be evaluated for this ion" not in {
        d.rule for d in confidence.demotions
    }


def test_an_unevaluable_rule_demotes_by_exactly_one_notch(twims_stratum):
    """Not to unsupported, and not nowhere. The ion may be perfectly fine."""
    flagged = {id(point) for point in twims_stratum.outliers}
    clean = [point for point in twims_stratum.points if id(point) not in flagged][0]
    checked = grade_correction(clean.reference_ccs, clean.ion.key, twims_stratum)
    unchecked = grade_correction(
        clean.reference_ccs, a_key_this_stratum_has_never_seen(), twims_stratum
    )
    assert checked.not_checked == () and unchecked.not_checked != ()
    assert checked.grade is ConfidenceGrade.WEAK
    assert unchecked.grade is ConfidenceGrade.UNSUPPORTED
    assert "a rule of the scheme could not be evaluated for this ion" in {
        d.rule for d in unchecked.demotions
    }


def test_two_unevaluable_rules_still_demote_only_one_notch():
    """The demotion is for grading on an incomplete scheme, not a tally of missing rules.

    A per-rule penalty would reach unsupported on any stratum with two gaps, reporting a
    confident refusal where the honest answer is only that less was checked.
    """
    stratum = stratum_of_size(TARGET_MATCHED_IONS, tag="onenotch")
    one = grade_correction(stratum.points[0].reference_ccs, stratum.points[0].ion.key, stratum)
    two = grade_correction(
        stratum.points[0].reference_ccs, a_key_this_stratum_has_never_seen(), stratum
    )
    assert len(one.not_checked) == 1 and len(two.not_checked) == 2
    assert one.grade is ConfidenceGrade.QUALIFIED
    assert two.grade is ConfidenceGrade.QUALIFIED


def test_the_demotion_for_an_unevaluable_rule_says_which_rule_could_not_run(twims_stratum):
    # A demotion whose reason a caller cannot read is a demotion they cannot argue with.
    confidence = grade_correction(
        twims_stratum.points[0].reference_ccs, a_key_this_stratum_has_never_seen(), twims_stratum
    )
    detail = [
        d.detail
        for d in confidence.demotions
        if d.rule == "a rule of the scheme could not be evaluated for this ion"
    ]
    assert len(detail) == 1
    assert "outlier check" in detail[0]
    assert "fewer checks than the scheme advertises" in detail[0]


def test_an_unevaluable_rule_cannot_demote_below_unsupported():
    """Unsupported has no notch beneath it, and the clamp is arithmetic, so assert it."""
    stratum = stratum_of_size(MIN_POINTS_FOR_DEMING - 1, tag="clamped")
    already = grade_correction(stratum.points[0].reference_ccs, stratum.points[0].ion.key, stratum)
    assert already.grade is ConfidenceGrade.UNSUPPORTED
    confidence = grade_correction(
        stratum.points[0].reference_ccs, a_key_this_stratum_has_never_seen(), stratum
    )
    assert confidence.grade is ConfidenceGrade.UNSUPPORTED


def test_the_scope_limit_on_the_outlier_rule_is_published_not_only_per_response():
    """A caller deciding whether to trust the service reads the scheme before sending anything.

    So the limit has to be legible there, not only in the not_checked of a response they have
    already received. This is the only rule in the scheme that is not "any ion".
    """
    published = {rule["rule"]: rule for rule in grade_rules()}
    applies = published["flagged as an outlier in its own stratum"]["applies_to"]
    assert "ONLY IONS ALREADY IN THIS CORPUS" in applies
    assert "new ion" in applies
    for name, rule in published.items():
        if name != "flagged as an outlier in its own stratum":
            assert "ONLY IONS ALREADY IN THIS CORPUS" not in rule["applies_to"], (
                f"{name}: if every rule carries the caveat, the caveat marks nothing"
            )


# --- 7c. the ceiling: what the top grade needs, and who can never reach it ----------------------
#
# No record in the seed corpus grades `supported` - 186 qualified, 216 weak, 15 unsupported,
# of 417. That is a fact about corpus size and not about the measurements: every applied
# stratum holds 23 to 46 matched ions against a TARGET_MATCHED_IONS of 100, so the population
# rule fires on all 417.
#
# Which makes it worth asserting, in both directions, that the top grade is LIVE CODE and that
# what stands between a caller and it is data. A grade nothing can ever return is a grade that
# should not be published on /confidence/rules, and a grade that turns out to be reachable for
# the wrong reason would be worse.


def a_stratum_of(count: int, *, tag: str, with_an_outlier: bool):
    """`count` matched ions, optionally with one that does not transfer.

    Leverage is only measurable where at least one ion is flagged, so a perfectly clean
    stratum reports the leverage check as unevaluable and can never reach `supported`
    however large it is. That is a real property of the scheme and this makes it visible.
    """
    records = []
    for index in range(count):
        analyte = fixtures.small_molecule(fixtures.synthetic_inchikey(f"{tag}{index:03d}"))
        reference = 150.0 + 10.0 * index
        bias = 1.07 if (with_an_outlier and index == 0) else 1.02
        records.append(
            fixtures.measurement(
                analyte=analyte, platform="dtims", ccs=reference, source="fixture DTIMS",
                doi=fixtures.DOI_A,
            )
        )
        records.append(
            fixtures.measurement(
                analyte=analyte, platform="twims", ccs=reference * bias,
                source="fixture TWIMS", doi=fixtures.DOI_A,
            )
        )
    return stratum_of(records)


def test_the_top_grade_is_reachable_so_publishing_it_is_not_a_fiction():
    """At the target population, with leverage measurable, an ion in the stratum grades supported.

    /confidence/rules publishes a four-grade scale. If the top one were unreachable by
    construction rather than by corpus size, publishing it would be advertising something
    the service cannot do - and the honest fix would be to remove it, not to explain it.
    """
    stratum = a_stratum_of(TARGET_MATCHED_IONS, tag="ceilA", with_an_outlier=True)
    assert stratum.n == TARGET_MATCHED_IONS
    assert thinly_populated(stratum) is None
    assert slope_leverage_percent(stratum) is not None, "this test needs leverage to be measurable"
    clean = [
        point for point in stratum.points
        if id(point) not in {id(outlier) for outlier in stratum.outliers}
    ][0]
    confidence = grade_correction(clean.reference_ccs, clean.ion.key, stratum)
    assert confidence.not_checked == ()
    assert confidence.demotions == ()
    assert confidence.grade is ConfidenceGrade.SUPPORTED


def test_below_the_target_population_the_ceiling_is_qualified_however_good_the_ion():
    """The other direction, and the reason no seed record grades supported.

    Same stratum, same clean ion, one fewer matched ion than the target. Nothing about the
    measurement changed; the grade falls because the corpus is small.
    """
    stratum = a_stratum_of(TARGET_MATCHED_IONS - 1, tag="ceilB", with_an_outlier=True)
    clean = [
        point for point in stratum.points
        if id(point) not in {id(outlier) for outlier in stratum.outliers}
    ][0]
    confidence = grade_correction(clean.reference_ccs, clean.ion.key, stratum)
    assert confidence.grade is ConfidenceGrade.QUALIFIED
    assert [d.rule for d in confidence.demotions] == ["thinly populated calibration group"]


def test_an_ion_the_corpus_has_never_seen_tops_out_at_qualified_at_any_corpus_size():
    """THE CEILING FOR THE POPULATION THIS PLATFORM EXISTS TO SERVE, and it does not move.

    Growing the corpus lifts the population demotion, so an ion already in it can reach
    `supported`. A genuinely new ion cannot: the outlier rule is a corpus lookup, it can
    never answer for them, and an unevaluable rule always costs one notch. So `qualified`
    is their ceiling however large this corpus becomes - a consequence of the scheme rather
    than of the data, and asserted here so it cannot quietly stop being true or quietly
    start being worse.
    """
    novel = a_stratum_of(5, tag="ceilnew", with_an_outlier=False).points[0].ion.key
    for count in (TARGET_MATCHED_IONS, TARGET_MATCHED_IONS * 3):
        stratum = a_stratum_of(count, tag=f"ceilC{count}", with_an_outlier=True)
        assert thinly_populated(stratum) is None
        confidence = grade_correction(stratum.points[0].reference_ccs, novel, stratum)
        assert confidence.grade is ConfidenceGrade.QUALIFIED, (
            f"a new ion graded {confidence.grade.value} on a stratum of {count}"
        )
        assert [d.rule for d in confidence.demotions] == [
            "a rule of the scheme could not be evaluated for this ion"
        ]


def test_the_population_threshold_is_not_quietly_lowered_to_make_the_top_grade_appear():
    """Adjusting a check until it passes is fitting the check to the data.

    The seed corpus's largest applied stratum holds 46. A threshold edged down to 46, or to
    anything a stratum already reaches, would make `supported` appear without a single new
    measurement - and would be the platform grading its own corpus as sufficient because it
    is the corpus it has.
    """
    assert TARGET_MATCHED_IONS == 100
    stratum = a_stratum_of(46, tag="ceilD", with_an_outlier=True)
    assert thinly_populated(stratum) is not None, (
        "the largest applied stratum in the seed corpus no longer triggers the population rule,"
        " which means the threshold moved"
    )


# --- 8. usable ----------------------------------------------------------------------------------


def test_only_an_unsupported_grade_is_unusable():
    for grade in ConfidenceGrade:
        confidence = Confidence(grade=grade)
        assert confidence.usable is (grade is not ConfidenceGrade.UNSUPPORTED), grade
    assert Confidence(grade=ConfidenceGrade.WEAK).usable is True
    assert Confidence(grade=ConfidenceGrade.UNSUPPORTED).usable is False


def test_usable_follows_the_graded_result_and_not_the_number_of_demotions(
    twims_stratum, ordinary_point
):
    weak = grade_correction(ordinary_point.reference_ccs, ordinary_point.ion.key, twims_stratum)
    unsupported = grade_correction(FAR_OUTSIDE_CCS, ordinary_point.ion.key, twims_stratum)
    assert len(weak.demotions) >= 2 and weak.usable is True
    assert unsupported.usable is False


# --- 9. the published scheme cannot drift from the code ------------------------------------------

# The demotion each published rule covers. Declared here rather than derived, so
# that a new rule in grading.py with no published entry, or a published entry no
# rule produces, fails this file instead of shipping. The two severe variants are
# separate demotions with their own wording and are published under the rule they
# are the severe end of.
PUBLISHED_BY_DEMOTION = {
    "outside the calibration range": "outside the calibration range",
    "far outside the calibration range": "outside the calibration range",
    "thinly populated calibration group": "thinly populated calibration group",
    "no correction can be fitted": "thinly populated calibration group",
    "flagged as an outlier in its own stratum": "flagged as an outlier in its own stratum",
    "a source behind the correction may not be used": "a source behind the correction may not be used",
    "the correction is driven by ions that do not transfer": (
        "the correction is driven by ions that do not transfer"
    ),
    "a rule of the scheme could not be evaluated for this ion": (
        "a rule of the scheme could not be evaluated for this ion"
    ),
}


def every_demotion_grade_correction_can_produce(twims_stratum):
    """Every (rule, grade) pair the rules can actually emit, gathered by driving them."""
    produced: set[tuple[str, ConfidenceGrade]] = set()
    flagged = {id(point) for point in twims_stratum.outliers}
    key = [point for point in twims_stratum.points if id(point) not in flagged][0].ion.key
    scenarios = [
        (JUST_OUTSIDE_CCS, key, twims_stratum),
        (FAR_OUTSIDE_CCS, key, twims_stratum),
        (twims_stratum.outliers[0].reference_ccs, twims_stratum.outliers[0].ion.key, twims_stratum),
    ]
    for count, tag in (
        (MIN_POINTS_FOR_DEMING - 1, "scenone"),
        (MIN_POINTS_FOR_DEMING + 1, "scentwo"),
        (MIN_POINTS_FOR_LIMITS_OF_AGREEMENT + 1, "scenthr"),
    ):
        stratum = stratum_of_size(count, tag=tag)
        scenarios.append((stratum.points[0].reference_ccs, stratum.points[0].ion.key, stratum))
    # A WELL-POPULATED STRATUM WHOSE ONLY FAULT IS AN UNEVALUABLE RULE. Without this the
    # scenarios never leave the grade at supported with a non-empty not_checked, so the
    # meta-rule's first rung - supported down to qualified - is never driven and the
    # falls_to test reads a genuine publication as an over-publication.
    at_target = stratum_of_size(TARGET_MATCHED_IONS, tag="scenfour")
    scenarios.append((at_target.points[0].reference_ccs, at_target.points[0].ion.key, at_target))
    spoiled = with_first_point_recorded_as(
        twims_stratum,
        fixtures.measurement(
            ccs=twims_stratum.points[0].other_ccs,
            source="a source nobody has read the terms of",
            reuse_status=ReuseStatus.UNVERIFIED,
        ),
    )
    scenarios.append((spoiled.points[0].reference_ccs, spoiled.points[0].ion.key, spoiled))
    for ccs, matched_ion_key, stratum in scenarios:
        for demotion in grade_correction(ccs, matched_ion_key, stratum).demotions:
            produced.add((demotion.rule, demotion.grade))
    return produced


def test_every_published_rule_carries_a_why_a_threshold_and_a_basis():
    rules = grade_rules()
    assert rules, "the published scheme is empty"
    seen = set()
    for rule in rules:
        assert set(rule) == {"rule", "falls_to", "why", "threshold", "basis", "applies_to"}
        for field in ("rule", "why", "threshold", "basis", "applies_to"):
            assert isinstance(rule[field], str) and rule[field].strip(), f"{rule['rule']}: empty {field}"
        assert rule["rule"] not in seen, f"{rule['rule']} is published twice"
        seen.add(rule["rule"])
        assert rule["falls_to"], f"{rule['rule']}: publishes no grade to fall to"
        for grade in rule["falls_to"]:
            assert ConfidenceGrade(grade) in ConfidenceGrade
            assert ConfidenceGrade(grade) is not ConfidenceGrade.SUPPORTED, (
                f"{rule['rule']}: a rule that falls to supported is not a demotion"
            )


def test_the_published_scheme_names_exactly_the_rules_the_code_can_produce(twims_stratum):
    """If the two ever drift, a reader is arguing with a scheme the code does not run."""
    produced = every_demotion_grade_correction_can_produce(twims_stratum)
    produced_rules = {rule for rule, _grade in produced}
    assert produced_rules == set(PUBLISHED_BY_DEMOTION), (
        "grading.py produces a demotion this file does not know about, or stopped producing one it did"
    )
    published = {rule["rule"] for rule in grade_rules()}
    assert published == set(PUBLISHED_BY_DEMOTION.values())


def test_every_published_rule_falls_exactly_as_far_as_the_code_makes_it_fall(twims_stratum):
    produced = every_demotion_grade_correction_can_produce(twims_stratum)
    observed: dict[str, set[str]] = {}
    for rule, grade in produced:
        observed.setdefault(PUBLISHED_BY_DEMOTION[rule], set()).add(grade.value)
    for rule in grade_rules():
        assert set(rule["falls_to"]) == observed[rule["rule"]], (
            f"{rule['rule']}: publishes {sorted(rule['falls_to'])} and produces {sorted(observed[rule['rule']])}"
        )


def test_a_demotion_renders_with_its_grade_and_its_reason(twims_stratum):
    # The published scheme is only as useful as the demotion a reader actually
    # receives, so the rendering carries both halves.
    demotion = correction_driven_by_outliers(twims_stratum)
    assert isinstance(demotion, Demotion)
    rendered = str(demotion)
    assert rendered.startswith(f"[{ConfidenceGrade.WEAK}]")
    assert demotion.rule in rendered
    assert demotion.detail in rendered


# --- 5b. the licence question this demotion actually asks --------------------------------------
#
# The same defect as in matching.blockers, in the same words and found at the same
# time: this rule asked `can_train_commercial` where it meant "may this platform
# use the record". The test above did not catch it because `unverified` is refused
# under both questions. academic_only is the branch where they differ.


def test_an_academic_only_source_does_not_demote_the_correction(twims_stratum):
    """A source this platform may use is not a source behind the correction that may not be used.

    Had this rule kept asking the commercial question, every correction derived
    from the steroid interplatform study would have graded UNSUPPORTED - not
    because anything about it is weak, but because a predicate two modules away
    answered a question nobody asked it.
    """
    academic = with_first_point_recorded_as(
        twims_stratum,
        fixtures.measurement(
            ccs=twims_stratum.points[0].other_ccs,
            source="Feuerstein et al., J. Am. Soc. Mass Spectrom. 2022",
            reuse_status=ReuseStatus.ACADEMIC_ONLY,
        ),
    )
    assert unverified_source_behind_it(academic) is None


@pytest.mark.parametrize(
    "status",
    [
        ReuseStatus.UNVERIFIED,
        ReuseStatus.OPEN_SHARE_ALIKE,
        ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES,
        ReuseStatus.EXCLUDED,
    ],
    ids=lambda status: status.value,
)
def test_a_source_this_platform_may_not_use_still_demotes_the_correction(twims_stratum, status):
    """The other side: widening the predicate once did not widen it to everything."""
    spoiled = with_first_point_recorded_as(
        twims_stratum,
        fixtures.measurement(
            ccs=twims_stratum.points[0].other_ccs,
            source="a source this platform may not use",
            reuse_status=status,
        ),
    )
    demotion = unverified_source_behind_it(spoiled)
    assert demotion is not None
    assert demotion.grade is ConfidenceGrade.UNSUPPORTED
    assert status.value in demotion.detail
