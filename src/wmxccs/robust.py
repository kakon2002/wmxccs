"""Passing-Bablok regression: a slope a few bad ions cannot move.

WHY THIS EXISTS, in one measurement
-----------------------------------
M3 reported that three outliers in twenty-four moved a Deming slope from the injected
1.020 to 1.031, and the correction applied to every well-behaved ion was therefore in
part the work of the three that were not. That was synthetic. On the real steroid
corpus the same effect is far larger, and it is concentrated in negative mode:

    stratum                          Deming slope   Passing-Bablok
    DTIMS/stepped vs TWIMS  [M-H]-       1.19024          1.05052
    DTIMS/stepped vs TIMS   [M-H]-       1.14698          1.00160
    DTIMS/single  vs TIMS   [M-H]-       1.10959          1.01609

Those Deming slopes carry intercepts of -31 to -40 square angstrom to compensate. A
slope of 1.147 between two platforms whose mean bias is -0.3 per cent is not a
description of those platforms; it is two leveraged ions and an intercept. The same
strata give Lin's concordance above 0.97, so nothing in the association statistics
flags it.

WHAT PASSING-BABLOK IS
----------------------
The slope is the median of the n(n-1)/2 pairwise slopes between all point pairs,
shifted by the number of slopes below -1. Because it is a median of pairwise slopes it
is the Theil-Sen estimator with a sign correction, and it inherits Theil-Sen's
breakdown point of 1 - 1/sqrt(2) = 0.293: up to 29 per cent of the points can be
arbitrarily bad before the slope can be dragged anywhere. Deming's breakdown point is
zero - one point placed far enough determines the answer.

It also suits this problem for a reason beyond robustness. Deming needs lambda, the
ratio of the two platforms' error variances, and on this corpus lambda is MEASURED in
some strata and ASSUMED EQUAL in others depending on whether every paired record
happens to report a convertible uncertainty. An estimator whose answer depends on that
accident is hard to defend. Passing-Bablok assumes only that the errors are
proportional to the measurement, which for a cross section reported to three
significant figures is the more defensible assumption, and it needs no lambda at all.

WHAT IT CANNOT DO, and these are not small
------------------------------------------
- It is NOT unbiased for lambda far from the slope squared. Where lambda is genuinely
  large and known, Deming is the better estimator and this is the worse one. Both are
  reported, always, for that reason.
- The rank confidence interval rests on the Kendall tau variance, which is a normal
  approximation. At n = 23 it is indicative, not exact, and this module says so.
- It says nothing about whether a LINEAR correction is the right model. A slope and an
  intercept fitted robustly are still a slope and an intercept.

Pure Python by deliberate choice, as the rest of this package is. `statistics` from the
standard library supplies the normal quantile, so the 95 per cent constant is derived
here rather than written down.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from statistics import NormalDist, median
from typing import Sequence

# The breakdown point of a median of pairwise slopes. Derived, not chosen: a single
# arbitrarily bad point among n contaminates n-1 of the n(n-1)/2 pairwise slopes, which
# is a fraction 2/n of them, and the median moves only once more than half the slopes
# have moved. Theil and Sen give 1 - 1/sqrt(2).
THEIL_SEN_BREAKDOWN = 1.0 - 1.0 / math.sqrt(2.0)

# Below this a median of pairwise slopes is not a robust statistic in any useful sense:
# with n = 4 a single bad point already contaminates half the slopes.
MIN_POINTS_FOR_ROBUST_SLOPE = 5

# The rank interval uses the Kendall tau variance, which is a normal approximation.
# Below this it is reported with a warning rather than withheld, because the interval is
# the BASIS SELECTOR in harmonization.py and a missing one is treated as "cannot tell".
MIN_POINTS_FOR_A_TRUSTWORTHY_SLOPE_INTERVAL = 10

SLOPE_INTERVAL_COVERAGE = 0.95

# POLICY, and marked as such because it is not derived from anything.
#
# Some negative pairwise slopes are inevitable rather than diagnostic: within one adduct
# the cross sections cluster tightly, so a pair of neighbouring ions has a tiny
# denominator and the measurement noise decides the sign. Measured over the eighteen real
# strata of the steroid corpus the fraction runs from 0.49 to 4.30 per cent, so a warning
# on "any negative slope" fires on all eighteen and means nothing.
#
# Five per cent is above every value observed and below anything that would suggest the
# two platforms genuinely disagree about the ORDER of these ions. It does not fire on the
# present corpus, which is the honest state: this is a guard for a stratum noisier than
# any yet seen, and its test uses a fixture built for it.
ORDER_REVERSAL_LIMIT = 0.05

TOO_FEW_FOR_ROBUST = (
    "a median of pairwise slopes needs at least {needed} points to be robust and this has {n}:"
    " with fewer, one bad point already contaminates half the pairwise slopes and the median moves with it"
)
ALL_X_TIED = (
    "every pair of points shares its reference value, so there is no pairwise slope to take a median of."
    " Two platforms cannot be related by a line through one abscissa"
)
INTERVAL_APPROXIMATE = (
    "the rank interval on the slope rests on the Kendall tau variance, a normal approximation, and n is"
    " {n}: it is indicative rather than exact below {needed} and is reported for direction, not for width"
)
SLOPE_INTERVAL_UNAVAILABLE = (
    "no rank interval could be formed on the slope, so whether it differs from 1 is UNKNOWN rather than no"
)
SHIFT_RUNS_OFF_THE_END = (
    "{k} of {total} pairwise slopes are below -1, so the Passing-Bablok shift puts the median at rank"
    " {index} of {total} and there is no such slope. More than half the point pairs are ordered inversely"
    " by the two platforms: no single line describes this, and returning the largest pairwise slope as"
    " though it were a median would be an answer rather than the truth"
)
ORDER_REVERSED = (
    "{negative} of {total} pairwise slopes are negative ({share:.1%}, above the {limit:.0%} this package"
    " treats as ordinary): for those pairs the two platforms rank the two ions in OPPOSITE orders. A single"
    " line through this stratum describes it poorly whatever estimator fits it, and a correction derived"
    " from that line will be confidently wrong for some ions in it"
)


@dataclass(frozen=True)
class PassingBablokFit:
    """A robust slope and intercept, with everything needed to distrust them.

    `slope_interval` is load-bearing rather than decorative: harmonization.py selects
    between a slope-derived and a median-derived correction on whether this interval
    excludes 1, so a None here means "cannot tell" and is treated as such.
    """

    slope: float
    intercept: float
    n_points: int
    pairwise_slopes: int
    shift_k: int  # slopes below -1; the Passing-Bablok correction
    negative_slopes: int
    tied_x_pairs: int  # dropped, and counted, because a drop nobody counts is a loss
    slopes_equal_minus_one: int  # excluded by the definition of the estimator
    even_count_interpolated: bool
    slope_interval: tuple[float, float] | None = None
    intercept_interval: tuple[float, float] | None = None
    slope_interval_coverage: float = SLOPE_INTERVAL_COVERAGE
    breakdown_point: float = THEIL_SEN_BREAKDOWN
    warnings: tuple[str, ...] = ()

    @property
    def order_reversed_share(self) -> float:
        """The fraction of pairwise slopes that are negative. 0.005 to 0.043 on the real corpus."""
        return self.negative_slopes / self.pairwise_slopes if self.pairwise_slopes else 0.0

    @property
    def ions_that_may_be_arbitrarily_bad(self) -> int:
        """How many points could be anywhere without moving the slope. The robustness budget."""
        return int(self.breakdown_point * self.n_points)

    @property
    def slope_distinguishable_from_unity(self) -> bool:
        """Whether the data supports a SIZE-DEPENDENT correction at all.

        THE BASIS SELECTOR, and the reason it lives on the fit rather than in a caller:
        it is a property of the estimator's own interval, not a policy. Where the
        interval contains 1 the slope is indistinguishable from unity, so a slope-derived
        correction and a constant offset are indistinguishable too - and of the two, the
        offset is the one that does not extrapolate.

        False when there is no interval. "Cannot tell" is not "yes".
        """
        if self.slope_interval is None:
            return False
        low, high = self.slope_interval
        return not (low <= 1.0 <= high)

    @property
    def implied_lambda(self) -> float:
        """The error-variance ratio this fit implicitly assumes: the slope squared.

        Passing-Bablok assumes errors proportional to the measurement, which is
        lambda = slope^2. Reported so it can be compared against the lambda Deming
        measured from the records: where the two disagree wildly, the two estimators are
        answering under incompatible assumptions and neither is obviously right.
        """
        return self.slope**2

    def predict(self, reference_ccs: float) -> float:
        """Forward: what this fit says the OTHER platform would report."""
        return self.intercept + self.slope * reference_ccs

    def correct(self, other_ccs: float) -> float | None:
        """Inverse: refer a value measured on the other platform BACK to the reference.

        None where the slope is zero, which no fit on real cross sections produces and
        which would otherwise divide by zero.
        """
        if self.slope == 0:
            return None
        return (other_ccs - self.intercept) / self.slope


def _pairwise_slopes(xs: Sequence[float], ys: Sequence[float]) -> tuple[list[float], int, int]:
    """Every pairwise slope, with the two exclusions counted rather than silent."""
    slopes: list[float] = []
    tied = 0
    minus_one = 0
    n = len(xs)
    for i in range(n):
        for j in range(i + 1, n):
            dx = xs[j] - xs[i]
            if dx == 0:
                # Vertical: no slope exists. Passing-Bablok drops these; this counts them,
                # because a stratum where many points share an abscissa is one where the
                # estimator is working on less data than its n suggests.
                tied += 1
                continue
            slope = (ys[j] - ys[i]) / dx
            if slope == -1.0:
                # Excluded by the definition of the estimator: the shift K counts slopes
                # strictly below -1, and a slope of exactly -1 belongs to neither side.
                minus_one += 1
                continue
            slopes.append(slope)
    return slopes, tied, minus_one


def _shifted_median(slopes: list[float], k: int) -> tuple[float, bool] | None:
    """The Passing-Bablok slope: the median of the sorted slopes, offset by k.

    The offset is what distinguishes this from plain Theil-Sen. Slopes below -1
    correspond to pairs the second platform orders inversely, and Passing and Bablok
    showed the median must be taken k places up the sorted list to stay consistent.

    None where the shift runs past the end of the list, which happens when most pairs are
    inverted. Not clamped: the largest pairwise slope is not a median.
    """
    count = len(slopes)
    if count % 2 == 1:
        index = (count + 1) // 2 - 1 + k
        if index >= count:
            return None
        return slopes[index], False
    if count // 2 + k >= count:
        return None
    low = slopes[count // 2 - 1 + k]
    high = slopes[count // 2 + k]
    if low > 0 and high > 0:
        # Geometric mean: a slope is a ratio, so the ratio scale is the natural one to
        # interpolate on. At slopes near unity this differs from the arithmetic mean by
        # less than 1e-6, which a test pins - the choice is principled, not numerical.
        return math.sqrt(low * high), True
    return (low + high) / 2.0, True


def passing_bablok(
    xs: Sequence[float], ys: Sequence[float], coverage: float = SLOPE_INTERVAL_COVERAGE
) -> PassingBablokFit | tuple[None, tuple[str, ...]]:
    """Fit a robust slope and intercept, or refuse and say why.

    Returns a `PassingBablokFit`, or `(None, refusals)` where no median of pairwise
    slopes exists. Refusing is a tuple rather than an exception because a stratum that
    cannot be fitted is an ordinary reportable state in this pipeline, not an error.
    """
    if len(xs) != len(ys):
        raise ValueError(f"{len(xs)} reference values against {len(ys)} other values")
    n = len(xs)
    if n < MIN_POINTS_FOR_ROBUST_SLOPE:
        return None, (TOO_FEW_FOR_ROBUST.format(needed=MIN_POINTS_FOR_ROBUST_SLOPE, n=n),)

    slopes, tied, minus_one = _pairwise_slopes(xs, ys)
    if not slopes:
        return None, (ALL_X_TIED,)
    slopes.sort()

    k = sum(1 for slope in slopes if slope < -1.0)
    negative = sum(1 for slope in slopes if slope < 0.0)
    shifted = _shifted_median(slopes, k)
    if shifted is None:
        return None, (
            SHIFT_RUNS_OFF_THE_END.format(
                k=k, total=len(slopes), index=(len(slopes) + 1) // 2 + k
            ),
        )
    slope, interpolated = shifted
    intercept = median([ys[i] - slope * xs[i] for i in range(n)])

    warnings: list[str] = []
    share = negative / len(slopes)
    if share > ORDER_REVERSAL_LIMIT:
        warnings.append(
            ORDER_REVERSED.format(
                negative=negative, total=len(slopes), share=share, limit=ORDER_REVERSAL_LIMIT
            )
        )

    # The rank interval. C is the square root of the Kendall tau variance, so the bounds
    # are ranks in the sorted slope list rather than a standard error around the slope.
    slope_interval: tuple[float, float] | None = None
    intercept_interval: tuple[float, float] | None = None
    total = len(slopes)
    spread = math.sqrt(n * (n - 1) * (2 * n + 5) / 18.0)
    c = NormalDist().inv_cdf(1.0 - (1.0 - coverage) / 2.0) * spread
    lower_rank = round((total - c) / 2.0)
    upper_rank = total - lower_rank + 1
    if 1 <= lower_rank and upper_rank + k <= total and lower_rank + k >= 1:
        low = slopes[lower_rank + k - 1]
        high = slopes[upper_rank + k - 1]
        if low <= high:
            slope_interval = (low, high)
            # The intercept bounds invert: a steeper slope gives a lower intercept.
            intercept_interval = (
                median([ys[i] - high * xs[i] for i in range(n)]),
                median([ys[i] - low * xs[i] for i in range(n)]),
            )
    if slope_interval is None:
        warnings.append(SLOPE_INTERVAL_UNAVAILABLE)
    elif n < MIN_POINTS_FOR_A_TRUSTWORTHY_SLOPE_INTERVAL:
        warnings.append(
            INTERVAL_APPROXIMATE.format(n=n, needed=MIN_POINTS_FOR_A_TRUSTWORTHY_SLOPE_INTERVAL)
        )

    return PassingBablokFit(
        slope=slope,
        intercept=intercept,
        n_points=n,
        pairwise_slopes=total,
        shift_k=k,
        negative_slopes=negative,
        tied_x_pairs=tied,
        slopes_equal_minus_one=minus_one,
        even_count_interpolated=interpolated,
        slope_interval=slope_interval,
        intercept_interval=intercept_interval,
        slope_interval_coverage=coverage,
        warnings=tuple(warnings),
    )
