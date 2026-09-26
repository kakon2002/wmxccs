"""Cross-platform statistics: how two platforms agree, and how they merely correlate.

THE ONE THING THIS MODULE EXISTS TO KEEP APART
-----------------------------------------------
ASSOCIATION is whether two platforms move together. AGREEMENT is whether they
give the same number. They are different questions, they have different answers,
and the second one is the one this platform was built to answer.

A platform can correlate at r = 0.995 with DTIMS and still run two per cent high
on every single ion. The correlation is not wrong; it is answering a question
nobody asked. The published interplatform figures make the point better than any
argument: TWIMS against DTIMS correlates at 0.9949, and 95 per cent of ions sit
within 2 per cent bias rather than within nothing.

So association and agreement are computed separately, reported in separate
blocks, and there is no combined score anywhere in this module. There is nothing
to average, because averaging them would produce a number that answers neither
question.

WHY THERE IS NO ORDINARY LEAST SQUARES HERE
-------------------------------------------
OLS assumes the x axis is known exactly and all the error is in y. Between two
ion mobility platforms that assumption is simply false: both axes are
measurements, both carry error, and an OLS slope is therefore biased towards zero
by an amount that depends on how noisy the x platform happens to be. Regressing
DTIMS on TWIMS and TWIMS on DTIMS would give two different answers, and neither
would be the relationship between the instruments.

Deming regression is the correct tool: it minimises the perpendicular distance
weighted by the ratio of the two error variances, and it gives one answer
whichever way round the platforms are put. OLS is therefore NOT OFFERED by this
module - not computed, not reported, not available behind a flag. If it is ever
added, it must be labelled as OLS at the point of use and must stay out of every
headline figure, because a reader who sees "slope" will assume the slope means
what the Deming slope means.

WHAT IS NEVER POOLED
--------------------
Statistics are computed per platform pair AND per calibration-group pair within
it. Two ions compared between the same two platforms but calibrated against
different reference sets, measured against different gases, or run at different
cyclic pass counts are not the same comparison, and a figure pooling them
describes no instrument pair that exists. A platform-pair figure is produced only
when that pair holds exactly one calibration-group stratum; where it holds more,
the pooled figure is REFUSED and the strata are listed instead.

OUTLIERS ARE REPORTED AND NEVER REMOVED
---------------------------------------
Published interplatform work finds the bulk of ions within 1 to 2 per cent and a
small fraction - under 1.5 per cent of them - disagreeing by as much as 7. Those
ions are the result. They are where a harmonization model will be wrong, they are
what a confidence grade exists to warn about, and a model that reduced average
disagreement while hiding them would be worse than no model. Nothing in this
module drops a point, and there is no parameter that would let it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from statistics import fmean, median
from typing import Mapping, Sequence

from .matching import MatchedIon, MatchingReport
from .models import DTIMSMethod, IMSType, UncertaintyType

# --- thresholds ---------------------------------------------------------------------------

# DERIVED, not chosen. With two points Pearson r is exactly +1 or -1 whatever the
# points are, so it carries no information at all; a third point is the first that
# can disagree with a line.
MIN_POINTS_FOR_CORRELATION = 3

# DERIVED. Deming fits a slope and an intercept, so three points leave one
# residual degree of freedom. At two the line passes through both and the residual
# is zero by construction, which reads as a perfect fit.
MIN_POINTS_FOR_DEMING = 3

# POLICY, and the loosest figure in this module. Limits of agreement are the mean
# difference plus and minus 1.96 standard deviations of the differences, so they
# are only as stable as that standard deviation. Ten is the conventional floor
# below which an SD of differences is too unstable to quote as a limit; it is not
# derived from anything in this repository and is marked as policy so that
# somebody with a real corpus can argue with it.
MIN_POINTS_FOR_LIMITS_OF_AGREEMENT = 10

# How far from ITS OWN STRATUM'S median an ion must sit to be FLAGGED, never
# dropped.
#
# Relative to the stratum median rather than to zero, and the difference matters
# more than it looks. A platform pair has a systematic offset - that is the whole
# finding, and it is what a harmonization model exists to correct. An ion sitting
# at the offset is not interesting; it is typical. The interesting ions are the
# ones far from it, because those are still wrong AFTER the systematic part is
# removed, and they are where a model will be confidently wrong.
#
# An absolute threshold gets this backwards. Tried first, at 2 per cent from zero,
# it flagged 14 of 24 ions in a corpus whose injected bias was exactly 2 per cent:
# every ordinary ion tripped it and the three that genuinely did not transfer were
# buried among them. A flag that fires on the typical case is not a flag.
#
# The margin is 2 per cent because that is the width of the published envelope for
# the worst-behaved pair - 95 per cent of TWIMS ions within 2 per cent of DTIMS -
# so an ion more than that from its own pair's centre is outside the spread its
# pair is expected to hold. Those figures are quoted in CONTEXT.md and have not
# been read from the paper here.
OUTLIER_MARGIN_PERCENT = 2.0

# Coverage IS reported against zero, at these two bands, because that is the
# published comparison: "95 per cent of ions within 1 per cent for TIMS and 2 per
# cent for TWIMS relative to DTIMS". Coverage answers "how close is this platform
# to the reference", which is an absolute question; outlier flagging answers "which
# ions do not behave like their neighbours", which is a relative one. Both are
# wanted, and they are not the same question.
COVERAGE_BANDS_PERCENT = (1.0, 2.0)

# 95 per cent limits, the Bland-Altman convention.
LOA_MULTIPLIER = 1.96

TOO_FEW_FOR_CORRELATION = (
    "{n} paired ion(s): a correlation needs at least {needed}, because with two points r is exactly"
    " +1 or -1 whatever the points are"
)
TOO_FEW_FOR_DEMING = (
    "{n} paired ion(s): a Deming slope needs at least {needed}. A slope from fewer has no residual"
    " left to disagree with it, and a two-point line fits perfectly by construction"
)
TOO_FEW_FOR_LOA = (
    "{n} paired ion(s): limits of agreement rest on the standard deviation of the differences, which"
    " is too unstable to quote below {needed}"
)
NO_VARIATION = (
    "every {side} value is the same number, so there is no variation to correlate or to fit a slope"
    " through"
)
POOLED_REFUSED = (
    "this platform pair holds {strata} calibration-group strata, so there is no single figure for it:"
    " values calibrated against different reference sets, gases or pass counts are not the same"
    " comparison. The strata are reported separately"
)
NO_PLATFORM_PAIR = (
    "every measurement in this corpus is on a single platform ({platforms}), so no ion is measured on"
    " two and there is nothing to compare. More values on the same platform will not change that"
)
NO_SHARED_ION = (
    "{platforms} platforms are represented, but no ion is held on more than one of them: the analyte"
    " sets do not overlap, or they differ in adduct, charge, gas or structural state"
)
ALL_SETS_UNUSABLE = (
    "{sets} matched set(s) exist and none may be used. A matched set is only as usable as its least"
    " usable member; the matched-ion report names which member and why"
)
NOTHING_MATCHED = "there are no matched ions at all, so there is nothing to compare"

SYNTHETIC_STATISTICS = (
    "every figure here was computed over synthetic fixtures. They are records declared in code, and"
    " no number from them describes any instrument. This report may be read and may not be quoted"
)


class NotQuotableError(Exception):
    """Statistics computed over synthetic fixtures were asked to be quoted as a result."""


# --- the arithmetic, kept small and readable ----------------------------------------------


def _moments(xs: Sequence[float], ys: Sequence[float]) -> tuple[float, float, float, float, float]:
    """(mean x, mean y, Sxx, Syy, Sxy), all as sums of squares about the mean."""
    mx, my = fmean(xs), fmean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return mx, my, sxx, syy, sxy


def pearson_r(xs: Sequence[float], ys: Sequence[float]) -> float | None:
    """Pearson correlation, or None where one axis does not vary at all."""
    _mx, _my, sxx, syy, sxy = _moments(xs, ys)
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def deming_slope_intercept(
    xs: Sequence[float], ys: Sequence[float], lam: float
) -> tuple[float, float] | None:
    """Deming slope and intercept of y on x at error-variance ratio `lam`.

    `lam` is var(error in y) / var(error in x). At lam = 1 this is orthogonal
    regression. Unlike OLS it is symmetric: fitting x on y gives the reciprocal
    slope and the same line.

    Returns None where the fit is not defined - no covariance at all, which means
    the points form a horizontal or vertical cloud with no direction to fit.
    """
    if lam <= 0:
        raise ValueError(f"the error-variance ratio must be positive, not {lam!r}")
    _mx, _my, sxx, syy, sxy = _moments(xs, ys)
    if sxy == 0:
        return None
    mx, my = fmean(xs), fmean(ys)
    discriminant = (syy - lam * sxx) ** 2 + 4 * lam * sxy**2
    slope = (syy - lam * sxx + math.sqrt(discriminant)) / (2 * sxy)
    return slope, my - slope * mx


def lins_concordance(xs: Sequence[float], ys: Sequence[float]) -> float | None:
    """Lin's concordance correlation coefficient.

    Pearson r asks whether the points lie on A line. This asks whether they lie on
    THE line - y = x - so it falls when a platform is precise and biased, which is
    exactly the case correlation alone hides. Population moments, as Lin defines it.
    """
    n = len(xs)
    if n < 2:
        return None
    mx, my, sxx, syy, sxy = _moments(xs, ys)
    var_x, var_y, cov = sxx / n, syy / n, sxy / n
    denominator = var_x + var_y + (mx - my) ** 2
    if denominator == 0:
        return None
    return 2 * cov / denominator


def _standard_deviation_of(record: object) -> float | None:
    """The record's uncertainty expressed as one standard deviation, or None.

    This is what `uncertainty_type` was carried for. A spread is only usable here
    once it is known what kind of spread it is, and two of the five kinds cannot
    be converted without something the record may not have:

    - SD converts to itself, TWO_SD is halved;
    - SEM is SD / sqrt(replicates), so it converts only where the replicate count
      is recorded, and is refused where it is not;
    - CI95 is NOT converted, deliberately. A 95 per cent interval may be reported
      as a half-width or as a full width, and this repository has never fixed
      which. Dividing by 1.96 on the assumption that it is a half-width would be
      a guess that silently halves or doubles every weight built on it;
    - UNKNOWN converts to nothing, which is the whole point of the value.
    """
    kind = getattr(record, "uncertainty_type", None)
    spread = getattr(record, "ccs_uncertainty", None)
    if kind is None or spread is None:
        return None
    if kind is UncertaintyType.SD:
        return float(spread)
    if kind is UncertaintyType.TWO_SD:
        return float(spread) / 2
    if kind is UncertaintyType.SEM:
        replicates = getattr(record, "replicates", None)
        if not replicates:
            return None
        return float(spread) * math.sqrt(replicates)
    return None


def error_variance_ratio(reference: Sequence[object], other: Sequence[object]) -> tuple[float, str]:
    """(lambda, how it was arrived at) for the Deming fit.

    Where both sides report a convertible uncertainty, the ratio is measured from
    the data. Where they do not, it falls back to 1 - orthogonal regression - and
    SAYS SO, because assuming equal error variances is an assumption and a reader
    has to know it was made rather than measured.
    """
    x_sd = [sd for sd in (_standard_deviation_of(record) for record in reference) if sd is not None]
    y_sd = [sd for sd in (_standard_deviation_of(record) for record in other) if sd is not None]
    if len(x_sd) == len(reference) and len(y_sd) == len(other) and x_sd and y_sd:
        var_x = fmean([sd**2 for sd in x_sd])
        var_y = fmean([sd**2 for sd in y_sd])
        if var_x > 0:
            return var_y / var_x, "measured from the reported uncertainties of both platforms"
    return 1.0, "ASSUMED EQUAL: not every paired record reports a convertible uncertainty, so the fit is orthogonal"


# --- one paired ion -----------------------------------------------------------------------


@dataclass(frozen=True)
class PairedPoint:
    """One ion measured on both platforms of a comparison."""

    ion: MatchedIon
    reference: object  # the measurement on the reference platform
    other: object

    @property
    def reference_ccs(self) -> float:
        return float(self.reference.ccs)

    @property
    def other_ccs(self) -> float:
        return float(self.other.ccs)

    @property
    def difference(self) -> float:
        """Other minus reference, in square angstrom. Signed, always."""
        return self.other_ccs - self.reference_ccs

    @property
    def difference_percent(self) -> float:
        """Signed difference as a percentage OF THE REFERENCE.

        Of the reference rather than of the pair mean, because that is what the
        published figures this will be compared against are relative to: "within
        1 per cent for TIMS and 2 per cent for TWIMS RELATIVE TO DTIMS".
        """
        return 100.0 * self.difference / self.reference_ccs

    def is_outlier_against(self, centre: float) -> bool:
        """Whether this ion sits more than the margin from its stratum's centre.

        Takes the centre as an argument rather than reading a constant, because a
        point cannot know its own stratum and a point that guessed would be
        comparing against the wrong centre.
        """
        return abs(self.difference_percent - centre) > OUTLIER_MARGIN_PERCENT

    def __str__(self) -> str:
        return (
            f"{self.ion.key}  {self.reference_ccs:g} -> {self.other_ccs:g}"
            f"  ({self.difference_percent:+.2f}%)"
        )


# --- the two blocks, kept apart ------------------------------------------------------------


@dataclass(frozen=True)
class Association:
    """Whether two platforms move together. NOT whether they agree."""

    n: int
    pearson_r: float | None = None
    refusals: tuple[str, ...] = ()

    @property
    def r_squared(self) -> float | None:
        """The square of Pearson r.

        The square of a CORRELATION, not a regression R squared. It says how much
        of the variation in one platform is shared with the other, and it says
        nothing whatever about whether the two give the same number.
        """
        return None if self.pearson_r is None else self.pearson_r**2

    def summary(self) -> list[str]:
        lines = [f"  association  (do they move together?)   n = {self.n}"]
        if self.pearson_r is None:
            for refusal in self.refusals:
                lines.append(f"    refused: {refusal}")
            return lines
        lines.append(f"    Pearson r    {self.pearson_r:+.4f}")
        lines.append(f"    r squared    {self.r_squared:.4f}   (correlation squared, NOT a regression fit)")
        return lines


@dataclass(frozen=True)
class Agreement:
    """Whether two platforms give the same number. The question that matters here."""

    n: int
    deming_slope: float | None = None
    deming_intercept: float | None = None
    lambda_used: float = 1.0
    lambda_basis: str = ""
    bias: float | None = None  # mean signed difference, square angstrom
    bias_sd: float | None = None
    limits_of_agreement: tuple[float, float] | None = None
    mean_difference_percent: float | None = None
    median_difference_percent: float | None = None
    mae: float | None = None
    mape: float | None = None
    rmse: float | None = None
    lins_ccc: float | None = None
    coverage: Mapping[float, float] = field(default_factory=dict)
    refusals: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    def summary(self) -> list[str]:
        lines = [f"  agreement    (do they give the same number?)   n = {self.n}"]
        if self.bias is not None:
            lines.append(
                f"    bias         {self.bias:+.3f} A^2"
                f"   mean {self.mean_difference_percent:+.3f}%   median {self.median_difference_percent:+.3f}%"
            )
        if self.limits_of_agreement is not None:
            low, high = self.limits_of_agreement
            lines.append(f"    95% limits   {low:+.3f} to {high:+.3f} A^2   (Bland-Altman)")
        if self.deming_slope is not None:
            lines.append(
                f"    Deming       slope {self.deming_slope:.5f}   intercept {self.deming_intercept:+.3f}"
                f"   lambda {self.lambda_used:.4g}"
            )
            lines.append(f"                 lambda {self.lambda_basis}")
        if self.lins_ccc is not None:
            lines.append(f"    Lin's CCC    {self.lins_ccc:+.4f}   (concordance, not correlation)")
        if self.mae is not None:
            lines.append(f"    MAE {self.mae:.3f} A^2   MAPE {self.mape:.3f}%   RMSE {self.rmse:.3f} A^2")
        if self.coverage:
            covered = "   ".join(f"within {band:g}%: {share:.1f}%" for band, share in sorted(self.coverage.items()))
            lines.append(f"    coverage     {covered}")
        for refusal in self.refusals:
            lines.append(f"    refused: {refusal}")
        for warning in self.warnings:
            lines.append(f"    weak: {warning}")
        return lines


@dataclass(frozen=True)
class ResidualTrend:
    """Whether the disagreement grows with the size of the ion.

    MASS IS NOT AVAILABLE. The brief asks for residuals against mass as well as
    against CCS, and no record in this repository carries a mass: CCSMeasurement
    has no mass or m/z field, and the one place a mass could have been derived -
    the glycan residue-mass table - was deliberately left behind with the rest of
    the glycan chemistry. So this is the CCS half only, and CCS is a reasonable
    size proxy for exactly the purpose intended. Adding mass means adding a field
    and a source for it, which is a schema decision for a later milestone.
    """

    n: int
    correlation_with_reference_ccs: float | None = None

    def summary(self) -> list[str]:
        if self.correlation_with_reference_ccs is None:
            return ["    residuals    not computable: too few points, or no variation in CCS"]
        return [
            f"    residuals    r = {self.correlation_with_reference_ccs:+.4f} against reference CCS"
            "   (size-dependent bias; mass is not held by any record, see ResidualTrend)"
        ]


@dataclass(frozen=True)
class StratumStatistics:
    """One comparison: two platforms, one calibration-group pair, and never pooled with another."""

    reference_platform: str
    other_platform: str
    reference_group: str
    other_group: str
    points: tuple[PairedPoint, ...]
    association: Association
    agreement: Agreement
    residuals: ResidualTrend
    synthetic: bool = False

    @property
    def n(self) -> int:
        return len(self.points)

    @property
    def centre_percent(self) -> float:
        """The stratum's own systematic offset, as the median signed difference.

        The median and not the mean, deliberately: the mean is dragged by the very
        ions this is used to find, so a few badly transferring ions would move the
        centre towards themselves and hide themselves behind it.
        """
        return median([point.difference_percent for point in self.points]) if self.points else 0.0

    @property
    def outliers(self) -> tuple[PairedPoint, ...]:
        """Ions far from their own stratum's offset. Reported, never removed.

        There is no parameter anywhere that removes them, and adding one would be
        a bug rather than a feature: these ions are the result. They are where a
        harmonization model will be wrong, and a model that reduced average
        disagreement while hiding them would be worse than no model at all.
        """
        centre = self.centre_percent
        return tuple(point for point in self.points if point.is_outlier_against(centre))

    def summary(self) -> list[str]:
        lines = [
            f"{self.other_platform} against {self.reference_platform}   (reference first)",
            f"  calibration groups: {self.reference_group}  |  {self.other_group}",
        ]
        lines += self.association.summary()
        lines += self.agreement.summary()
        lines += self.residuals.summary()
        if self.outliers:
            lines.append(
                f"    OUTLIERS     {len(self.outliers)} of {self.n} ion(s) more than"
                f" {OUTLIER_MARGIN_PERCENT:g}% from this stratum's centre of"
                f" {self.centre_percent:+.3f}%, reported and kept:"
            )
            for point in sorted(self.outliers, key=lambda p: -abs(p.difference_percent)):
                lines.append(f"      {point}")
        return lines


# --- building the comparisons ---------------------------------------------------------------


def _is_primary(record: object) -> bool:
    """Stepped-field DTIMS: the only CCS that stands on first principles."""
    return (
        getattr(record, "ims_type", None) is IMSType.DTIMS
        and getattr(record, "dtims_method", None) is DTIMSMethod.STEPPED_FIELD
    )


def _platform_of(record: object) -> str:
    ims = getattr(record, "ims_type", None)
    method = getattr(record, "dtims_method", None)
    return str(ims) if method is None else f"{ims}/{method}"


def choose_reference(one: str, other: str, records_by_platform: Mapping[str, object]) -> tuple[str, str]:
    """Which platform is the reference, and why it is not arbitrary.

    A primary value - stepped-field DTIMS - is the reference whenever one is
    present, because it is the only CCS in the comparison that does not itself
    rest on somebody else's reference values. Comparing a calibrated platform
    against a primary one measures the calibrated platform. Comparing two
    calibrated platforms measures their difference and anchors nothing, which is
    worth knowing when reading the result.

    With no primary side, the reference is the alphabetically first platform.
    That IS arbitrary, and it is arbitrary in a way that does not change any
    conclusion: the Deming slope is symmetric, and the sign of the bias flips
    with the choice and is reported with the platforms named either way round.
    """
    primary = [name for name in (one, other) if _is_primary(records_by_platform[name])]
    if len(primary) == 1:
        reference = primary[0]
        return reference, other if reference == one else one
    first, second = sorted((one, other))
    return first, second


def _paired_points(ion: MatchedIon, reference_group: str, other_group: str) -> PairedPoint | None:
    """The one point this ion contributes, or None if it contributes none.

    An ion contributes AT MOST ONE point to a stratum. Where a side holds more
    than one measurement in the same calibration group - replicates - the ion is
    skipped rather than averaged or than having one silently chosen: averaging
    would merge two originals, which this platform does not do anywhere, and
    picking the first would make the result depend on file order.
    """
    on_reference = [m.record for m in ion.measurements if str(m.record.calibration_group) == reference_group]
    on_other = [m.record for m in ion.measurements if str(m.record.calibration_group) == other_group]
    if len(on_reference) != 1 or len(on_other) != 1:
        return None
    return PairedPoint(ion=ion, reference=on_reference[0], other=on_other[0])


def _association(xs: Sequence[float], ys: Sequence[float]) -> Association:
    n = len(xs)
    if n < MIN_POINTS_FOR_CORRELATION:
        return Association(
            n=n, refusals=(TOO_FEW_FOR_CORRELATION.format(n=n, needed=MIN_POINTS_FOR_CORRELATION),)
        )
    r = pearson_r(xs, ys)
    if r is None:
        side = "reference" if len(set(xs)) == 1 else "comparison"
        return Association(n=n, refusals=(NO_VARIATION.format(side=side),))
    return Association(n=n, pearson_r=r)


def _agreement(points: Sequence[PairedPoint]) -> Agreement:
    n = len(points)
    xs = [point.reference_ccs for point in points]
    ys = [point.other_ccs for point in points]
    differences = [point.difference for point in points]
    percents = [point.difference_percent for point in points]

    refusals: list[str] = []
    warnings: list[str] = []

    # Bias and the difference figures need only one point to be arithmetically
    # defined, and they are reported from one, because a single disagreement is a
    # real observation. What a single point cannot support is a LIMIT, a SLOPE or
    # a CORRELATION, and each of those refuses separately below.
    bias = fmean(differences)
    bias_sd = None
    limits = None
    if n >= 2:
        bias_sd = math.sqrt(sum((d - bias) ** 2 for d in differences) / (n - 1))
    if n >= MIN_POINTS_FOR_LIMITS_OF_AGREEMENT and bias_sd is not None:
        limits = (bias - LOA_MULTIPLIER * bias_sd, bias + LOA_MULTIPLIER * bias_sd)
    else:
        refusals.append(TOO_FEW_FOR_LOA.format(n=n, needed=MIN_POINTS_FOR_LIMITS_OF_AGREEMENT))

    slope = intercept = None
    lam, basis = error_variance_ratio([p.reference for p in points], [p.other for p in points])
    if n < MIN_POINTS_FOR_DEMING:
        refusals.append(TOO_FEW_FOR_DEMING.format(n=n, needed=MIN_POINTS_FOR_DEMING))
    else:
        fitted = deming_slope_intercept(xs, ys, lam)
        if fitted is None:
            refusals.append(
                "the paired values have no covariance at all, so there is no direction to fit a slope through"
            )
        else:
            slope, intercept = fitted

    ccc = lins_concordance(xs, ys) if n >= MIN_POINTS_FOR_CORRELATION else None

    coverage = {
        band: 100.0 * sum(1 for p in percents if abs(p) <= band) / n for band in COVERAGE_BANDS_PERCENT
    }
    if n < MIN_POINTS_FOR_LIMITS_OF_AGREEMENT:
        warnings.append(
            f"{n} ion(s) is a small comparison: every figure here is reported because it is arithmetically"
            " defined, not because it is settled"
        )

    return Agreement(
        n=n,
        deming_slope=slope,
        deming_intercept=intercept,
        lambda_used=lam,
        lambda_basis=basis,
        bias=bias,
        bias_sd=bias_sd,
        limits_of_agreement=limits,
        mean_difference_percent=fmean(percents),
        median_difference_percent=median(percents),
        mae=fmean([abs(d) for d in differences]),
        mape=fmean([abs(p) for p in percents]),
        rmse=math.sqrt(fmean([d**2 for d in differences])),
        lins_ccc=ccc,
        coverage=coverage,
        refusals=tuple(refusals),
        warnings=tuple(warnings),
    )


def _residuals(points: Sequence[PairedPoint]) -> ResidualTrend:
    n = len(points)
    if n < MIN_POINTS_FOR_CORRELATION:
        return ResidualTrend(n=n)
    sizes = [point.reference_ccs for point in points]
    residuals = [point.difference_percent for point in points]
    return ResidualTrend(n=n, correlation_with_reference_ccs=pearson_r(sizes, residuals))


@dataclass(frozen=True)
class PlatformPair:
    """One pair of platforms, and every calibration-group stratum inside it.

    THE POOLED FIGURE IS PRODUCED ONLY WHEN THERE IS NOTHING TO POOL. With one
    stratum, the pair figure and the stratum figure are the same numbers and the
    pair is reported. With more than one, there is no honest single figure: the
    strata were calibrated against different reference sets, or measured against
    different gases, or run at different cyclic pass counts, and averaging them
    would describe an instrument pair that does not exist. The refusal names how
    many strata there are, and the strata are reported one by one.
    """

    reference_platform: str
    other_platform: str
    strata: tuple[StratumStatistics, ...]

    @property
    def n(self) -> int:
        return sum(stratum.n for stratum in self.strata)

    @property
    def pooled(self) -> StratumStatistics | None:
        return self.strata[0] if len(self.strata) == 1 else None

    @property
    def pooling_refusal(self) -> str | None:
        return None if len(self.strata) <= 1 else POOLED_REFUSED.format(strata=len(self.strata))

    @property
    def outliers(self) -> tuple[PairedPoint, ...]:
        return tuple(point for stratum in self.strata for point in stratum.outliers)

    def summary(self) -> list[str]:
        lines = [
            "",
            f"=== {self.other_platform} against {self.reference_platform} ===",
            f"  {self.n} paired ion(s) in {len(self.strata)} calibration-group stratum/strata",
        ]
        if self.pooling_refusal is not None:
            lines.append(f"  NO POOLED FIGURE: {self.pooling_refusal}")
        for stratum in self.strata:
            lines.append("")
            lines += [f"  {line}" for line in stratum.summary()]
        return lines


@dataclass(frozen=True)
class ComparisonReport:
    """Every cross-platform comparison the corpus supports, and every one it refuses."""

    pairs: tuple[PlatformPair, ...] = ()
    ions_considered: int = 0
    ions_skipped_for_replicates: tuple[str, ...] = ()
    synthetic: bool = False
    unusable_sets_skipped: tuple[str, ...] = ()
    # Pairs refused because one member's calibration traces to the other member. Kept on the
    # report because a refusal nobody can count is indistinguishable from an absence of data.
    pairs_refused_as_circular: tuple[str, ...] = ()
    # Why there is no comparison, where there is none. A report that says only
    # "nothing here" leaves a reader unable to tell a corpus that is too small
    # from one that can never work however much of it arrives, and those call for
    # completely different decisions.
    no_comparison_reason: str | None = None

    @property
    def n_points(self) -> int:
        return sum(pair.n for pair in self.pairs)

    @property
    def outliers(self) -> tuple[PairedPoint, ...]:
        return tuple(point for pair in self.pairs for point in pair.outliers)

    @property
    def quotable(self) -> bool:
        """Whether any figure here may be presented as describing real instruments."""
        return not self.synthetic

    def refusal(self) -> str | None:
        return SYNTHETIC_STATISTICS if self.synthetic else None

    def summary(self) -> str:
        lines = [
            "Cross-platform statistics",
            f"  matched ions considered   {self.ions_considered}",
            f"  paired points used        {self.n_points}",
            f"  platform pairs            {len(self.pairs)}",
            f"  ions far from their pair  {len(self.outliers)}"
            f"  (more than {OUTLIER_MARGIN_PERCENT:g}% from their own stratum's centre; reported, never removed)",
        ]
        if self.unusable_sets_skipped:
            lines.append(
                f"  matched sets not used     {len(self.unusable_sets_skipped)}"
                "  (a member may not be used; see the matched-ion report for which)"
            )
        if self.ions_skipped_for_replicates:
            lines.append(
                f"  ions held for replicates  {len(self.ions_skipped_for_replicates)}"
                "  (more than one value on one side of a comparison; averaging would merge two originals)"
            )
        if self.synthetic:
            lines.append("")
            lines.append(f"  SYNTHETIC: {SYNTHETIC_STATISTICS}")
        if not self.pairs:
            lines.append("")
            lines.append("  NO CROSS-PLATFORM COMPARISON IS POSSIBLE FROM THIS CORPUS.")
            if self.no_comparison_reason:
                lines.append(f"  {self.no_comparison_reason}")
        for pair in self.pairs:
            lines += pair.summary()
        return "\n".join(lines)


CIRCULAR_PAIR = (
    "refused as circular: the {calibrated} value was calibrated against reference values that"
    " {origin}, which is what the {reference} value in this pair is. Comparing them measures how"
    " well the calibration reproduces its own reference set, not how the two platforms differ."
    " {detail}"
)


def circularity_between(one, other) -> str | None:
    """Why this pair must not be compared, or None if it may be.

    THE RISK IS NOT HYPOTHETICAL AND THE SCHEMA WAS BUILT FOR IT. A TWIMS cross section is
    calibrated against reference values, those reference values were usually measured by
    DTIMS, and if the comparison then puts that TWIMS value against the same DTIMS values it
    calibrated from, the agreement reported is an artefact of the calibration. It is not a
    weak comparison, it is a different measurement entirely: it measures the calibration.

    `CalibrationReference` has recorded the lineage since M0 - the reference set by name, the
    publication its values came from, and the platform they were measured on - and until
    26 September 2026 NOTHING READ IT at comparison time. The circularity was detectable and
    undetected, which is the failure class LIMITATIONS 4.5 collects. This is the detection.

    A PAIR IS REFUSED, NOT FLAGGED. A flag on a figure is read by whoever is looking for it;
    the Bush Lab MicroSource records are the case in hand, and their reference values are
    their own laboratory's nitrogen DTIMS measurements, so a Bush Lab DTIMS record arriving
    later must not quietly form a stratum with them.

    WHAT THIS CATCHES is provable circularity: the reference values were published in the
    very paper being compared against, matched by DOI, on the platform the lineage names.
    WHAT IT DOES NOT CATCH is a same-laboratory chain across two different papers, where
    whether the values are independent is a judgement about the two papers rather than a fact
    in either record. That is stated here rather than implied by silence, because a guard that
    reads broader than it is, is the thing this repository keeps finding.
    """
    for calibrated, reference, who in ((one, other, "first"), (other, one, "second")):
        lineage = getattr(calibrated, "calibration_reference", None)
        if lineage is None:
            continue
        if lineage.platform is None or lineage.platform != reference.ims_type:
            continue
        if lineage.method is not None and lineage.method != reference.dtims_method:
            continue
        # THE TEST IS THE DOI, and only the DOI. An earlier draft of this also tried to
        # match on the reference-set NAME, which was both broken - it compared a string to
        # a CalibrationReference object, so it could never fire - and wrong in principle:
        # two records naming one calibrant share a calibrant, which is not circularity.
        # Circularity is that THIS record's reference values were published in THAT record's
        # paper, and the DOI is the only thing in either record that can establish it.
        if not lineage.doi:
            continue
        # THE PUBLICATION IS MATCHED IN TWO PLACES, and the second is not sloppiness. Where a
        # value was obtained from a compiler rather than from its publisher - a database whose
        # terms were read, citing a paper whose licence was not - the licence gate correctly
        # refuses a reuse claim attached to the paper's DOI, so `doi` is left null and the
        # publication is named in `source_locator` instead. That is the decision recorded for
        # CCSbase in LIMITATIONS 7E and it is how the Bush Lab records are seeded. A
        # circularity test that looked only at `doi` would therefore never fire on exactly the
        # records it was written for.
        locator = getattr(reference, "source_locator", None) or ""
        if not (lineage.doi == getattr(reference, "doi", None) or lineage.doi in locator):
            continue
        return CIRCULAR_PAIR.format(
            calibrated=calibrated.ims_type,
            reference=reference.ims_type,
            origin=f"were published in doi:{lineage.doi}",
            detail=(
                f"Reference set {lineage.reference_set!r}, measured on"
                f" {lineage.platform}{'/' + str(lineage.method) if lineage.method else ''}."
                f" The {who} record of the pair is the calibrated one."
            ),
        )
    return None


def compare_platforms(matching: MatchingReport) -> ComparisonReport:
    """Every cross-platform comparison the matched ions support.

    Only sets that are USABLE are compared: a matched set holding a member nobody
    may use is skipped and counted, not quietly included. The exception is the
    synthetic blocker, which is carried through to the report instead - a synthetic
    corpus must still be computable, or the arithmetic could never be exercised at
    all, and `quotable` is what stops the result being presented as real.
    """
    considered = 0
    skipped_replicates: list[str] = []
    unusable: list[str] = []
    refused_circular: list[str] = []
    synthetic = bool(matching.synthetic_groups)

    # stratum key -> points
    strata: dict[tuple[str, str, str, str], list[PairedPoint]] = {}

    for ion in matching.matched:
        considered += 1
        blocking = [problem for problem in ion.blockers if not problem.startswith("this set is synthetic")]
        if blocking:
            unusable.append(str(ion.key))
            continue
        by_platform = {m.platform: m.record for m in ion.measurements}
        names = sorted(by_platform)
        for index, one in enumerate(names):
            for other in names[index + 1 :]:
                reference_platform, other_platform = choose_reference(one, other, by_platform)
                # REFUSED BEFORE A STRATUM EXISTS, so a circular pair cannot reach a fit even
                # as one point among many. Counted on the report rather than dropped: the
                # gate here raises, it does not thin the corpus quietly.
                circular = circularity_between(
                    by_platform[reference_platform], by_platform[other_platform]
                )
                if circular:
                    refused_circular.append(f"{ion.key}: {circular}")
                    continue
                reference_group = str(by_platform[reference_platform].calibration_group)
                other_group = str(by_platform[other_platform].calibration_group)
                point = _paired_points(ion, reference_group, other_group)
                if point is None:
                    skipped_replicates.append(str(ion.key))
                    continue
                key = (reference_platform, other_platform, reference_group, other_group)
                strata.setdefault(key, []).append(point)

    built: dict[tuple[str, str], list[StratumStatistics]] = {}
    for (reference_platform, other_platform, reference_group, other_group), points in strata.items():
        xs = [point.reference_ccs for point in points]
        ys = [point.other_ccs for point in points]
        stratum = StratumStatistics(
            reference_platform=reference_platform,
            other_platform=other_platform,
            reference_group=reference_group,
            other_group=other_group,
            points=tuple(points),
            association=_association(xs, ys),
            agreement=_agreement(points),
            residuals=_residuals(points),
            synthetic=synthetic,
        )
        built.setdefault((reference_platform, other_platform), []).append(stratum)

    pairs = tuple(
        PlatformPair(
            reference_platform=reference_platform,
            other_platform=other_platform,
            strata=tuple(sorted(found, key=lambda s: (s.reference_group, s.other_group))),
        )
        for (reference_platform, other_platform), found in sorted(built.items())
    )

    reason = None
    if not pairs:
        platforms = sorted({m.platform for group in matching.all_groups for m in group.measurements})
        if unusable and not considered - len(unusable):
            reason = ALL_SETS_UNUSABLE.format(sets=len(unusable))
        elif not matching.all_groups:
            reason = NOTHING_MATCHED
        elif len(platforms) < 2:
            reason = NO_PLATFORM_PAIR.format(platforms=", ".join(platforms) or "none")
        else:
            reason = NO_SHARED_ION.format(platforms=len(platforms))

    return ComparisonReport(
        pairs=pairs,
        ions_considered=considered,
        ions_skipped_for_replicates=tuple(dict.fromkeys(skipped_replicates)),
        pairs_refused_as_circular=tuple(dict.fromkeys(refused_circular)),
        synthetic=synthetic,
        unusable_sets_skipped=tuple(dict.fromkeys(unusable)),
        no_comparison_reason=reason,
    )


def assert_quotable(report: ComparisonReport) -> None:
    """Raise unless every figure in `report` describes real instruments.

    The call that stands between a synthetic fixture and a number in a document.
    Whatever reports a statistic calls this first.
    """
    refusal = report.refusal()
    if refusal is not None:
        raise NotQuotableError(refusal)
