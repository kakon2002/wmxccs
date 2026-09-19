"""Passing-Bablok: that it is exact where it should be, and immovable where it must be.

The estimator exists for one reason - a few ions must not decide the correction applied
to every other ion - so the tests that matter are the ones that break it deliberately and
check that it does not move. A robustness claim tested only on well-behaved data is a
claim nobody checked.
"""

from __future__ import annotations

import math

import pytest

from wmxccs.robust import (
    MIN_POINTS_FOR_ROBUST_SLOPE,
    ORDER_REVERSAL_LIMIT,
    THEIL_SEN_BREAKDOWN,
    passing_bablok,
)
from wmxccs.statistics import deming_slope_intercept

# A clean line, and a realistic one: cross sections of this size, spaced as steroids are.
CLEAN_X = [150.0, 158.0, 166.0, 174.0, 182.0, 190.0, 198.0, 206.0, 212.0, 220.0, 228.0]
SLOPE = 1.04
INTERCEPT = -3.0
CLEAN_Y = [INTERCEPT + SLOPE * x for x in CLEAN_X]


def fit(xs, ys):
    result = passing_bablok(xs, ys)
    assert not isinstance(result, tuple), f"refused: {result}"
    return result


# --- exactness -------------------------------------------------------------------------------


def test_an_exact_line_is_recovered_exactly():
    """A median of pairwise slopes on collinear points is the slope, to floating point."""
    result = fit(CLEAN_X, CLEAN_Y)
    assert result.slope == pytest.approx(SLOPE, abs=1e-12)
    assert result.intercept == pytest.approx(INTERCEPT, abs=1e-9)
    assert result.negative_slopes == 0
    assert result.shift_k == 0
    assert result.tied_x_pairs == 0


def test_the_identity_relation_gives_slope_one_and_intercept_zero():
    result = fit(CLEAN_X, list(CLEAN_X))
    assert result.slope == pytest.approx(1.0, abs=1e-12)
    assert result.intercept == pytest.approx(0.0, abs=1e-9)
    assert result.slope_distinguishable_from_unity is False, "an identity relation is not a size dependence"


# --- robustness, which is the entire point ------------------------------------------------------


def test_one_catastrophic_ion_moves_the_robust_slope_not_at_all_and_deming_by_a_factor_of_six():
    """The measurement that justifies this module existing.

    Deming's breakdown point is zero: one point placed far enough determines the answer.
    Here it takes the slope from 1.04 to over 6 while Passing-Bablok does not move at all.
    """
    wrecked = list(CLEAN_Y)
    wrecked[-1] = 900.0

    robust_clean = fit(CLEAN_X, CLEAN_Y).slope
    robust_wrecked = fit(CLEAN_X, wrecked).slope
    deming_clean = deming_slope_intercept(CLEAN_X, CLEAN_Y, 1.0)[0]
    deming_wrecked = deming_slope_intercept(CLEAN_X, wrecked, 1.0)[0]

    assert robust_wrecked == pytest.approx(robust_clean, abs=1e-12), "the robust slope must not move"
    assert deming_wrecked > 5 * deming_clean, "and Deming must, or this test proves nothing"


def test_the_robust_slope_survives_its_stated_breakdown_budget():
    """Corrupt as many points as the breakdown point allows and the slope must hold.

    The budget is not a guess: floor(0.293 * n). This corrupts exactly that many and
    requires the slope to be unmoved, then corrupts well past it and requires that it IS
    moved - because an estimator that never moved would not be measuring anything.
    """
    result = fit(CLEAN_X, CLEAN_Y)
    budget = result.ions_that_may_be_arbitrarily_bad
    assert budget == int(THEIL_SEN_BREAKDOWN * len(CLEAN_X))
    assert budget >= 1

    within = list(CLEAN_Y)
    for i in range(budget):
        within[i] = 1000.0 + 100.0 * i
    assert fit(CLEAN_X, within).slope == pytest.approx(SLOPE, abs=1e-9)

    # Past the budget the estimator must either move or refuse, and it must not
    # quietly return the old answer. Corrupting more than half the points inverts most
    # pairs, at which point the Passing-Bablok shift runs past the end of the slope list
    # and the fit is REFUSED rather than clamped to the largest pairwise slope.
    beyond = list(CLEAN_Y)
    for i in range(len(CLEAN_Y) // 2 + 1):
        beyond[i] = 1000.0 + 100.0 * i
    result = passing_bablok(CLEAN_X, beyond)
    if isinstance(result, tuple):
        assert "no single line describes this" in result[1][0]
    else:
        assert result.slope != pytest.approx(SLOPE, abs=1e-6)


def test_a_stratum_whose_pairs_are_all_inverted_is_refused_and_never_clamped():
    """The shift can run past the end of the slope list, and the largest slope is not a median.

    UNCONDITIONAL. An earlier version of this test branched on the result and asserted
    something true in either branch, so a clamped answer satisfied it and the mutation that
    clamps survived the sweep. Every pairwise slope here is far below -1, so K equals the
    number of slopes and the shifted index is past the end by construction.
    """
    xs = [100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0]
    ys = [900.0, 700.0, 500.0, 300.0, 100.0, 50.0, 10.0]
    result = passing_bablok(xs, ys)

    assert isinstance(result, tuple), "a stratum with every pair inverted must be REFUSED, not fitted"
    assert result[0] is None
    assert "pairwise slopes are below -1" in result[1][0]
    assert "no single line describes this" in result[1][0]


def test_the_refusal_names_the_rank_it_would_have_needed():
    """So a reader can check the arithmetic rather than trust the refusal."""
    xs = [100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0]
    ys = [900.0, 700.0, 500.0, 300.0, 100.0, 50.0, 10.0]
    message = passing_bablok(xs, ys)[1][0]
    pairs = len(xs) * (len(xs) - 1) // 2
    assert f"{pairs} pairwise slopes" in message
    assert f"of {pairs}" in message


# --- the rank interval, which is the basis selector ----------------------------------------------


def test_a_genuine_size_dependence_gives_an_interval_that_excludes_one():
    result = fit(CLEAN_X, CLEAN_Y)
    low, high = result.slope_interval
    assert low <= SLOPE <= high
    assert not (low <= 1.0 <= high)
    assert result.slope_distinguishable_from_unity is True


def test_a_pure_offset_gives_an_interval_that_contains_one():
    """Where the platforms differ by a constant, the slope must not be claimed to differ.

    This is the case that selects the MEDIAN basis, and getting it wrong would apply a
    size-dependent correction to data that shows no size dependence.
    """
    offset = [x + 4.0 for x in CLEAN_X]
    result = fit(CLEAN_X, offset)
    low, high = result.slope_interval
    assert low <= 1.0 <= high
    assert result.slope_distinguishable_from_unity is False


def test_an_absent_interval_reads_as_cannot_tell_rather_than_as_no():
    """The selector must not treat a missing interval as evidence of no size dependence."""
    scattered_x = [100.0, 100.5, 101.0, 100.2, 100.8, 100.1]
    scattered_y = [200.0, 190.0, 210.0, 195.0, 205.0, 188.0]
    result = fit(scattered_x, scattered_y)
    if result.slope_interval is None:
        assert result.slope_distinguishable_from_unity is False
        assert any("UNKNOWN rather than no" in w for w in result.warnings)


# --- refusals and counted exclusions -------------------------------------------------------------


@pytest.mark.parametrize("n", [0, 1, 2, 3, 4])
def test_too_few_points_is_refused_rather_than_fitted(n):
    xs = [100.0 + i for i in range(n)]
    result = passing_bablok(xs, list(xs))
    assert isinstance(result, tuple)
    assert result[0] is None
    assert str(MIN_POINTS_FOR_ROBUST_SLOPE) in result[1][0]


def test_every_point_sharing_one_reference_value_is_refused():
    result = passing_bablok([180.0] * 8, [1.0, 2, 3, 4, 5, 6, 7, 8])
    assert isinstance(result, tuple)
    assert "one abscissa" in result[1][0]


def test_tied_reference_values_are_dropped_and_counted_rather_than_silently_lost():
    """A drop nobody counts is a stratum working on less data than its n suggests."""
    xs = CLEAN_X + [CLEAN_X[0]]
    ys = CLEAN_Y + [CLEAN_Y[0] + 1.0]
    result = fit(xs, ys)
    assert result.tied_x_pairs == 1
    assert result.n_points == len(xs)
    assert result.pairwise_slopes == len(xs) * (len(xs) - 1) // 2 - result.tied_x_pairs


def test_the_order_reversal_warning_is_proportional_and_does_not_fire_on_clean_data():
    """It fires on a stratum with no relationship, and not on one with a clean one.

    A warning that fires on every stratum is a warning nobody reads: over the eighteen
    real strata the negative-slope share runs 0.5 to 4.3 per cent, all below the limit.
    """
    assert fit(CLEAN_X, CLEAN_Y).warnings == () or all(
        "OPPOSITE" not in w for w in fit(CLEAN_X, CLEAN_Y).warnings
    )
    noise_x = [100.0, 101, 102, 103, 104, 105, 106, 107]
    noise_y = [200.0, 199, 203, 198, 204, 197, 205, 196]
    noisy = fit(noise_x, noise_y)
    assert noisy.order_reversed_share > ORDER_REVERSAL_LIMIT
    assert any("OPPOSITE" in w for w in noisy.warnings)


# --- inverting the fit, which is what a correction does --------------------------------------------


def test_the_correction_inverts_the_fit_exactly():
    result = fit(CLEAN_X, CLEAN_Y)
    for x in CLEAN_X:
        predicted = result.predict(x)
        assert result.correct(predicted) == pytest.approx(x, abs=1e-9)


def test_the_implied_lambda_is_the_slope_squared():
    """What Passing-Bablok assumes about the two platforms' error variances.

    Reported so it can be compared against the lambda Deming measured from the records:
    where the two disagree the estimators are answering under incompatible assumptions.
    """
    result = fit(CLEAN_X, CLEAN_Y)
    assert result.implied_lambda == pytest.approx(SLOPE**2, abs=1e-12)


def test_the_breakdown_point_is_derived_and_not_a_written_down_number():
    assert THEIL_SEN_BREAKDOWN == pytest.approx(1 - 1 / math.sqrt(2), abs=1e-15)
    assert 0.29 < THEIL_SEN_BREAKDOWN < 0.30
