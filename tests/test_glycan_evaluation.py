"""Scoring CCS predictions: the maths, and what it refuses to score.

Every array here is hand-computed. None of these numbers came from running the
code under test, and none of them is a real measurement.
"""

import dataclasses
import math
import statistics

import pytest

from wmxglycan.contracts import DataMaturity, MaturityStamp
from wmxglycan.evaluation import (
    MIN_GROUP_RECORDS,
    MIN_GROUP_STRUCTURES,
    NO_PAIRS,
    EvaluationReport,
    GroupResult,
    Metrics,
    ModeConfusionError,
    PooledError,
    group_refusal,
    pool,
)
from wmxglycan.splits import SplitMode
from wmxglycan.models import CalibrationGroup, Derivatisation, DriftGas, DTIMSMethod, IMSType, ReducingEndLabel

# Worked by hand: relative errors 1%, 2%, 1%, 10%.
OBSERVED = [100.0, 200.0, 300.0, 400.0]
PREDICTED = [101.0, 204.0, 297.0, 440.0]

PROVISIONAL = MaturityStamp(data_maturity=DataMaturity.PROVISIONAL, training_record_count=0)


def a_group(calibrant="fixture calibrant"):
    return CalibrationGroup(
        IMSType.TWIMS, None, DriftGas.N2, calibrant, "[M+2H]2+", ReducingEndLabel.NATIVE, Derivatisation.UNDERIVATISED
    )


# --- the maths, against numbers computed by hand ------------------------------


def test_the_metrics_on_a_hand_computed_array():
    metrics = Metrics.of(OBSERVED, PREDICTED)
    assert metrics.n == 4
    # errors 1, 2, 1, 10 -> sorted 1, 1, 2, 10 -> median is the mean of the middle two
    assert metrics.median_relative_error == pytest.approx(1.5)
    assert metrics.mean_absolute_error == pytest.approx((1 + 4 + 3 + 40) / 4)
    assert metrics.rmse == pytest.approx(math.sqrt((1 + 16 + 9 + 1600) / 4))
    assert metrics.within_1_percent == pytest.approx(50.0)
    assert metrics.within_2_percent == pytest.approx(75.0)


def test_an_even_number_of_errors_takes_the_mean_of_the_middle_two():
    # Pinned because the median is the headline number and the convention is a choice.
    assert statistics.median([1, 2, 3, 4]) == 2.5
    metrics = Metrics.of([100.0, 100.0], [101.0, 103.0])
    assert metrics.median_relative_error == pytest.approx(2.0)


def test_relative_error_is_relative_to_the_measurement_not_the_prediction():
    # Predicting 200 for a measured 100 is 100 per cent out, not 50.
    assert Metrics.of([100.0], [200.0]).median_relative_error == pytest.approx(100.0)


def test_a_perfect_prediction_scores_zero():
    metrics = Metrics.of(OBSERVED, list(OBSERVED))
    assert metrics.median_relative_error == 0.0
    assert metrics.within_1_percent == 100.0


# --- what it refuses to score -------------------------------------------------


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf, 0.0, -1.0])
def test_a_prediction_that_is_not_a_positive_finite_number_is_refused_not_dropped(bad):
    with pytest.raises(ValueError) as caught:
        Metrics.of([100.0, 200.0], [101.0, bad])
    message = str(caught.value)
    assert "index 1" in message  # names where
    assert repr(bad) in message  # and what
    # Nothing was clamped, skipped or absolute-valued into plausibility.
    assert "not a positive finite number" in message


@pytest.mark.parametrize("bad", [math.nan, 0.0, -5.0])
def test_a_measurement_that_cannot_be_divided_by_is_refused(bad):
    with pytest.raises(ValueError, match="observed value at index 0"):
        Metrics.of([bad], [100.0])


def test_no_pairs_gives_our_refusal_not_a_statistics_error():
    with pytest.raises(ValueError) as caught:
        Metrics.of([], [])
    assert NO_PAIRS in str(caught.value)
    # StatisticsError IS a ValueError, so reaching statistics.median first would
    # let a caller's row-level except clause read "no data" as "one bad row".
    assert not isinstance(caught.value, statistics.StatisticsError)


def test_a_length_mismatch_names_both_lengths():
    with pytest.raises(ValueError) as caught:
        Metrics.of([100.0, 200.0], [101.0])
    assert "2 observed" in str(caught.value) and "1 predicted" in str(caught.value)


def test_a_boolean_is_not_a_measurement():
    # True == 1 numerically, which would otherwise sail through a > 0 check.
    with pytest.raises(ValueError):
        Metrics.of([True], [100.0])


# --- square angstroms never leave their calibration group ---------------------


def test_square_angstrom_numbers_are_never_pooled_across_calibration_groups():
    # The absence of these fields is the enforcement, so the test is on the shape
    # of the type rather than on any value it happens to hold.
    fields = {f.name for f in dataclasses.fields(PooledError)}
    assert {"mean_absolute_error", "rmse"}.isdisjoint(fields)
    assert "record_weighted_mean_of_group_medians" in fields  # dimensionless, so it does pool


def test_pooled_coverage_is_never_rounded_into_a_false_100_per_cent():
    # 996 of 1000 is 99.6 per cent. One decimal place in the POOLED summary as
    # well as in a group's, because rounding it up would claim perfect coverage.
    # The group-level version of this rule has always been tested and this one
    # never was: a single mutation anchor matched both sites, mutated whichever
    # came first, and reported a kill it had earned somewhere else entirely.
    pooled = PooledError(
        n=1000,
        groups=3,
        record_weighted_mean_of_group_medians=0.5,
        within_1_percent=99.6,
        within_2_percent=99.9,
    )
    text = pooled.summary()
    assert "within 1% 99.6%" in text and "within 2% 99.9%" in text
    assert "100%" not in text


def test_the_pooled_figure_is_not_called_a_median_because_it_is_not_one():
    # A mean of medians is bounded by the extreme group medians, so it cannot see
    # the pooled tail. The field is named for what it computes.
    fields = {f.name for f in dataclasses.fields(PooledError)}
    assert "median_relative_error" not in fields
    thirty = GroupResult(a_group("C"), n=30, structures=10, metrics=Metrics.of([100.0] * 30, [100.1] * 30))
    thirty_one = GroupResult(a_group("D"), n=31, structures=10, metrics=Metrics.of([100.0] * 31, [110.0] * 31))
    pooled = pool([thirty, thirty_one])
    # The true median of the 61 pooled errors is 0.1; this figure is over fifty
    # times that, which is exactly why it must not be labelled a median.
    assert pooled.record_weighted_mean_of_group_medians > 5.0
    assert "median relative error" not in pooled.summary()
    assert "record-weighted mean of group medians" in pooled.summary()


def test_pooling_weights_each_measurement_equally():
    one = GroupResult(a_group("A"), n=10, structures=10, metrics=Metrics.of([100.0] * 10, [101.0] * 10))
    two = GroupResult(a_group("B"), n=30, structures=10, metrics=Metrics.of([100.0] * 30, [103.0] * 30))
    pooled = pool([one, two])
    assert pooled.groups == 2 and pooled.n == 40
    assert pooled.record_weighted_mean_of_group_medians == pytest.approx((1.0 * 10 + 3.0 * 30) / 40)


def test_pooling_nothing_gives_nothing_rather_than_zero():
    refused = GroupResult(a_group(), n=3, structures=1, refusal="too small")
    assert pool([refused]) is None
    assert pool([]) is None


# --- a group too small to score is refused, and says why ----------------------


def test_a_calibration_group_below_the_floor_cannot_be_scored():
    reason = group_refusal(records=MIN_GROUP_RECORDS - 1, structures=MIN_GROUP_STRUCTURES)
    assert reason is not None
    assert str(MIN_GROUP_RECORDS) in reason and str(MIN_GROUP_STRUCTURES) in reason
    assert group_refusal(records=MIN_GROUP_RECORDS - 1, structures=MIN_GROUP_STRUCTURES - 1) is not None
    assert group_refusal(records=MIN_GROUP_RECORDS, structures=MIN_GROUP_STRUCTURES) is None


def test_remeasuring_one_glycan_is_not_analyte_diversity():
    # Plenty of records, but they are the same few molecules over and over.
    assert group_refusal(records=100, structures=2) is not None


def test_every_reported_number_carries_its_record_count():
    # n must be the count the metrics actually scored, or the number a reader
    # checks is not the number the summary prints.
    scored = GroupResult(a_group("A"), n=25, structures=12, metrics=Metrics.of([100.0] * 25, [101.0] * 25))
    refused = GroupResult(a_group("B"), n=3, structures=1, refusal=group_refusal(3, 1))
    for result in (scored, refused):
        if result.metrics is not None:
            assert result.n >= MIN_GROUP_RECORDS
            assert result.metrics.n == result.n
            assert result.refusal is None
        else:
            assert result.refusal  # a group without a number says why, rather than showing nothing
            assert "not scored" in result.summary()


def test_a_group_result_cannot_disagree_with_its_own_metrics():
    with pytest.raises(ValueError, match="scored"):
        GroupResult(a_group(), n=25, structures=12, metrics=Metrics.of(OBSERVED, PREDICTED))


def test_a_group_result_carries_a_score_or_a_refusal_and_never_both():
    with pytest.raises(ValueError, match="never both"):
        GroupResult(a_group(), n=4, structures=2, metrics=Metrics.of(OBSERVED, PREDICTED), refusal="too small")
    with pytest.raises(ValueError, match="never both"):
        GroupResult(a_group(), n=0, structures=0)  # neither would render "not scored: None"


# --- the report cannot be quoted without its context --------------------------


def a_report(**overrides):
    fields = dict(
        analyte_level="composition",
        by_study=False,
        baseline="group_median",
        records_in=0,
        records_not_scored=0,
        groups=(),
        pooled=None,
        maturity=PROVISIONAL,
    )
    fields.update(overrides)
    return EvaluationReport(**fields)


def test_the_report_is_stamped_provisional_and_m3_never_validates():
    report = a_report()
    assert report.maturity.data_maturity is DataMaturity.PROVISIONAL
    assert report.maturity.training_record_count == 0


def test_a_report_cannot_claim_validated_in_this_milestone():
    validated = MaturityStamp(data_maturity=DataMaturity.VALIDATED, training_record_count=5000)
    with pytest.raises(ValueError, match="provisional"):
        a_report(maturity=validated)


def test_the_summary_names_the_level_and_the_maturity_on_the_first_line():
    # A number quoted without its grouping level is a number that will be over-read.
    first = a_report(records_in=7, records_not_scored=7).summary().splitlines()[0]
    assert "composition" in first and "provisional" in first
    assert "group_median" in first


def test_a_report_with_nothing_scored_says_so_rather_than_printing_zero():
    assert "nothing was scored" in a_report().summary()


# --- the guardrail: an isomer result cannot speak as performance ----------------


def a_scored_report(**overrides):
    scored = GroupResult(a_group("A"), n=25, structures=12, metrics=Metrics.of([100.0] * 25, [101.0] * 25))
    pooled = pool([scored])
    return a_report(groups=(scored,), pooled=pooled, records_in=25, **overrides)


def test_a_deployment_result_yields_its_headline():
    report = a_scored_report()
    assert report.is_performance
    assert report.headline.n == 25


def test_an_isomer_result_refuses_to_yield_a_headline():
    # The number would look right, which is exactly the danger: nothing
    # downstream questions a plausible error rate.
    report = a_scored_report(mode=SplitMode.ISOMER_DISCRIMINATION)
    assert not report.is_performance
    with pytest.raises(ModeConfusionError, match="not deployment performance"):
        report.headline


def test_the_refusal_cannot_be_swallowed_by_a_broad_except():
    # A reporting path wrapped in `except (ValueError, TypeError)` must not be
    # able to catch this and print the diagnostic as though it were performance.
    assert not issubclass(ModeConfusionError, (ValueError, TypeError))


def test_an_isomer_summary_never_renders_the_pooled_figure():
    report = a_scored_report(mode=SplitMode.ISOMER_DISCRIMINATION)
    text = report.summary()
    assert "DIAGNOSTIC, NOT PERFORMANCE" in text
    assert "optimistic bound" in text
    # The figure itself is absent, so there is nothing in the text to quote.
    assert "record-weighted mean of group medians" not in text
    assert "pooled over" not in text
    # ...whereas the deployment rendering of the same numbers does show it.
    assert "pooled over" in a_scored_report().summary()


def test_a_scored_deployment_report_with_nothing_pooled_still_refuses_a_headline():
    with pytest.raises(ModeConfusionError, match="nothing was scored"):
        a_report().headline


def test_the_mode_is_a_field_so_it_travels_with_the_result():
    assert "mode" in {f.name for f in dataclasses.fields(EvaluationReport)}
    assert a_report().mode is SplitMode.DEPLOYMENT  # the default is the reportable one


def test_the_maturity_is_a_field_not_a_footnote():
    # It cannot be dropped by a caller who renders only the metrics.
    assert "maturity" in {f.name for f in dataclasses.fields(EvaluationReport)}


def test_a_maturity_stamp_cannot_be_flipped_to_validated_after_the_check():
    # Stamps are routinely shared between reports, so a mutable one could be
    # changed after the check that refused it and every holder would re-render.
    report = a_report()
    with pytest.raises(Exception):
        report.maturity.data_maturity = DataMaturity.VALIDATED
    with pytest.raises(Exception):
        report.maturity.training_record_count = 5000
    assert report.maturity.data_maturity is DataMaturity.PROVISIONAL
    assert "provisional" in report.summary()


def test_a_report_cannot_claim_a_pooled_figure_it_did_not_earn():
    scored = GroupResult(a_group("A"), n=25, structures=12, metrics=Metrics.of([100.0] * 25, [101.0] * 25))
    # A scored group with no pooled figure prints "nothing was scored" above real metrics.
    with pytest.raises(ValueError, match="cannot be absent"):
        a_report(groups=(scored,), pooled=None)
    # A pooled figure must agree with the groups it claims to pool.
    with pytest.raises(ValueError, match="claims"):
        a_report(
            groups=(scored,),
            pooled=PooledError(n=25, groups=42, record_weighted_mean_of_group_medians=1.0,
                               within_1_percent=100.0, within_2_percent=100.0),
        )


def test_a_report_cannot_lose_more_records_than_it_was_given():
    with pytest.raises(ValueError, match="went unscored"):
        a_report(records_in=10, records_not_scored=999)


def test_coverage_is_never_rounded_up_into_a_false_hundred_per_cent():
    # 996 of 1000 is 99.6 per cent. Printing "100%" would claim perfect coverage.
    nearly = Metrics.of([100.0] * 1000, [100.5] * 996 + [110.0] * 4)
    assert nearly.within_1_percent == pytest.approx(99.6)
    assert "100%" not in nearly.summary()
    assert "99.6%" in nearly.summary()
