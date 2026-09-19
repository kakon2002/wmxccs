"""M4: the harmonization model. A corrected cross section BESIDE the original, never instead.

WHAT THIS DOES
--------------
For each platform pair and each calibration group separately - never pooled - it fits
the relationship between the two platforms three ways, decides which of them the data
actually supports, and produces a corrected value with a prediction interval and a
confidence grade. The original measurement is returned untouched alongside it.

THE THREE CORRECTIONS, AND WHY THE HEADLINE IS CHOSEN RATHER THAN FIXED
----------------------------------------------------------------------
- SLOPE: from the Deming fit already in statistics.py. Symmetric, errors in both
  variables, and with a breakdown point of zero.
- ROBUST_SLOPE: from Passing-Bablok in robust.py. A median of pairwise slopes, breakdown
  point 0.293.
- MEDIAN: the stratum's median signed difference, applied as a constant per-cent offset.
  Cannot express a size-dependent correction, and cannot be moved by a few ions.

All three are always computed and always reported. The HEADLINE is chosen by a rule
derived from the estimator rather than by preference:

    basis = ROBUST_SLOPE   if the Passing-Bablok rank interval on the slope EXCLUDES 1
            MEDIAN         otherwise

The reasoning: where the slope cannot be distinguished from 1, a slope-derived
correction and a constant offset are two descriptions of the same data, and of the two
the offset is the one that does not extrapolate. Measured on the real corpus this picks
ROBUST_SLOPE in 10 of 18 strata and MEDIAN in 8.

That rule was not the first choice. Preferring the robust slope unconditionally fails on
this corpus for a reason worth recording: in negative mode the Deming slopes reach
1.19 with intercepts of -40 square angstrom, and Passing-Bablok correctly flattens them
towards 1.00 - at which point the slope is doing nothing and an affine correction fitted
through it can change SIGN inside its own fitted range. The median offset does not.

THE INTERVAL: grouped leave-one-compound-out jackknife+
------------------------------------------------------
The strata hold 23 to 41 ions. A split-conformal design would have to divide that into a
fit half and a calibration half, and readiness.py establishes that a 90 per cent
conformal interval is the FULL OBSERVED RANGE for any calibration set below 19 - so a
split would either starve the fit or produce an interval that excludes nothing. Neither
is acceptable.

Jackknife+ (Barber, Candes, Ramdas and Tibshirani, 2021) avoids the split entirely: it
refits leaving out one group at a time and combines the leave-one-out predictions with
their own residuals, so every ion serves as both fit and calibration. n_cal = n = 23, 29
or 41, all of which are above 19, so every interval here is informative. The price is the
guarantee: jackknife+ proves coverage of 1 - 2*alpha rather than 1 - alpha. That is
stated in the output rather than quietly nominal - `guaranteed_coverage` is 0.80 where
`nominal_coverage` is 0.90, and both travel with the interval.

The leave-one-out groups are COMPOUNDS, not ions. In this corpus each compound appears
once per stratum so the two coincide, but a compound with two conformers would otherwise
have one conformer calibrating a fit the other trained.

WHAT IS APPLIED AND WHAT IS ONLY PUBLISHED
------------------------------------------
Of the 18 strata, only the 9 whose reference platform is stepped-field DTIMS are ever
APPLIED to a measurement. Stepped-field DTIMS is the only primary method here - the only
CCS that does not itself rest on somebody else's reference values - so it is the only
sensible thing to refer a value TO. The other 9 strata (calibrated platform against
calibrated platform) are computed, reported and never applied, because a correction
between two calibrated platforms anchors nothing. There is NO transitive composition:
this package will not refer TWIMS to TIMS and then TIMS to the drift tube.

THE SCOPE CAVEAT IS STRUCTURAL
------------------------------
Every output of this module carries a `ScopeStamp`, which is REQUIRED and not defaulted,
and whose scope is DERIVED from the provenance of the records behind the fit. See
scope.py. A figure from this module cannot be serialised without its scope, and the claim
"interlaboratory reproducibility" is not representable anywhere in this package.

WHAT THIS IS NOT
----------------
It is not validated. Every corpus it has seen comes from one study, so the maturity stamp
is provisional and `DataMaturity.VALIDATED` is unreachable while the scope is
within-study. Coverage measured by leave-one-out on the same data that set the quantile
is not evidence of calibration - it is arithmetically pinned near k/n - so this module
reports interval WIDTH and the TAIL RATIO instead, and says so.
"""

from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from dataclasses import dataclass, field
from statistics import median
from typing import Iterable, Mapping, Sequence

from .contracts import CorrectionBasis
from .grading import Confidence, grade_correction
from .readiness import (
    CONFORMAL_ALPHA,
    MaturityStamp,
    conformal_quantile_index,
    interval_readiness,
    smallest_informative_calibration_set,
)
from .robust import MIN_POINTS_FOR_ROBUST_SLOPE, PassingBablokFit, passing_bablok
from .scope import Claim, ScopeStamp, assert_may_be_quoted_as, stamp_over, stamp_over_stamps
from .statistics import ComparisonReport, PairedPoint, StratumStatistics

# Jackknife+ proves coverage of 1 - 2*alpha, not 1 - alpha. Both are reported: the
# nominal level the quantile is taken at, and the level actually guaranteed.
NOMINAL_COVERAGE = 1.0 - CONFORMAL_ALPHA
GUARANTEED_COVERAGE = 1.0 - 2.0 * CONFORMAL_ALPHA

# A leave-one-out refit needs enough points left to still be a robust fit.
MIN_POINTS_FOR_JACKKNIFE = MIN_POINTS_FOR_ROBUST_SLOPE + 1

NO_ROBUST_FIT = "no robust slope could be fitted, so no correction is offered from this stratum: {why}"
NOT_APPLIED_NOT_PRIMARY = (
    "this stratum compares two CALIBRATED platforms ({reference} against {other}), so a correction derived"
    " from it anchors nothing: both sides already rest on somebody else's reference values. It is computed"
    " and reported, and never applied to a measurement. Only a stratum whose reference is stepped-field"
    " DTIMS - the one primary method here - is applied"
)
NO_INTERVAL = (
    "no prediction interval could be formed: {why}. A corrected value without an interval is not offered,"
    " because the interval is the part that says how much to trust it"
)
COVERAGE_IS_NOT_EVIDENCE = (
    "leave-one-out coverage on the same points that set the quantile is pinned near k/n by construction and"
    " is NOT evidence of calibration. Reported here for completeness only; the informative figures are the"
    " interval WIDTH and the TAIL RATIO below"
)
PRIMARY_REFERENCE = "DTIMS/stepped_field"


@dataclass(frozen=True)
class LeaveOneOutFits:
    """The refits and residuals a jackknife+ interval is built from.

    STORED RATHER THAN REDUCED TO A BAND, because a jackknife+ interval is a property of
    the QUERY as well as of the stratum: each left-out refit predicts the query value, and
    each prediction is offset by that refit's own residual. Collapsing this to one band
    per stratum produces the spread of the stratum instead of an interval around anything.
    """

    slopes: tuple[float, ...]
    intercepts: tuple[float, ...]
    residuals: tuple[float, ...]
    groups: int
    alpha: float

    @property
    def refits(self) -> int:
        return len(self.slopes)

    def interval_for(self, other_ccs: float) -> JackknifePlusInterval | str:
        """The interval around the correction of one value, or why there is none."""
        total = len(self.residuals)
        if not total:
            return "no leave-one-out residuals, so no interval"
        lows: list[float] = []
        highs: list[float] = []
        for slope, intercept, residual in zip(self.slopes, self.intercepts, self.residuals):
            if slope == 0:
                return "a leave-one-out refit has a zero slope, so its correction does not invert"
            predicted = (other_ccs - intercept) / slope
            lows.append(predicted - residual)
            highs.append(predicted + residual)
        lows.sort()
        highs.sort()
        low_index = max(1, math.floor(self.alpha * (total + 1)))
        high_index = min(total, math.ceil((1.0 - self.alpha) * (total + 1)))
        refusal, warning = interval_readiness(total, self.alpha)
        ordered = sorted(self.residuals)
        quantile_residual = ordered[min(total, conformal_quantile_index(total, self.alpha)) - 1]
        return JackknifePlusInterval(
            low=lows[low_index - 1],
            high=highs[high_index - 1],
            at_value=other_ccs,
            nominal_coverage=1.0 - self.alpha,
            guaranteed_coverage=1.0 - 2.0 * self.alpha,
            groups_left_out=self.groups,
            refits=self.refits,
            quantile_index_low=low_index,
            quantile_index_high=high_index,
            largest_residual=max(self.residuals),
            quantile_residual=quantile_residual,
            interval_is_informative=total >= smallest_informative_calibration_set(self.alpha),
            readiness_warning=warning,
            readiness_refusal=refusal,
        )


@dataclass(frozen=True)
class JackknifePlusInterval:
    """A prediction interval from leave-one-group-out refits, and what it is worth.

    `nominal_coverage` is the level the quantiles are taken at; `guaranteed_coverage` is
    what jackknife+ actually proves. They differ by construction and both are carried,
    because an interval quoted at its nominal level overstates itself by a factor of two
    in alpha.
    """

    low: float
    high: float
    at_value: float  # the OTHER platform's value this interval was computed around
    nominal_coverage: float
    guaranteed_coverage: float
    groups_left_out: int
    refits: int
    quantile_index_low: int
    quantile_index_high: int
    largest_residual: float
    quantile_residual: float
    interval_is_informative: bool
    readiness_warning: str | None = None
    readiness_refusal: str | None = None

    @property
    def width(self) -> float:
        return self.high - self.low

    @property
    def tail_ratio(self) -> float:
        """Largest residual over the quantile residual. 1.0 means one ion sets the width."""
        if self.quantile_residual == 0:
            return math.inf
        return self.largest_residual / self.quantile_residual

    @property
    def driven_by_one_ion(self) -> bool:
        """Whether the width is set by the single worst residual.

        True exactly when the quantile IS the maximum, which readiness.py shows happens
        for every calibration set below 19. It is reported because an interval set by one
        ion is a statement about that ion.
        """
        return self.tail_ratio <= 1.0


@dataclass(frozen=True)
class StratumCorrection:
    """One platform pair, one calibration group, three corrections and a scope.

    Never pooled with another. There is no method here that combines two of these.
    """

    stratum: StratumStatistics
    scope: ScopeStamp
    median_offset_percent: float
    robust: PassingBablokFit | None = None
    robust_refusals: tuple[str, ...] = ()
    loo: LeaveOneOutFits | None = None
    warnings: tuple[str, ...] = ()
    refusals: tuple[str, ...] = ()

    @property
    def n(self) -> int:
        return self.stratum.n

    @property
    def reference_platform(self) -> str:
        return self.stratum.reference_platform

    @property
    def other_platform(self) -> str:
        return self.stratum.other_platform

    @property
    def basis(self) -> CorrectionBasis:
        """The headline basis, DERIVED from the robust fit's own interval.

        Not a preference and not a setting. Where the rank interval on the slope excludes
        1 the data supports a size-dependent correction and the robust slope is the
        headline; otherwise the slope is indistinguishable from unity, a slope-derived
        correction and a constant offset describe the same data, and the offset is the one
        that does not extrapolate.
        """
        if self.robust is not None and self.robust.slope_distinguishable_from_unity:
            return CorrectionBasis.ROBUST_SLOPE
        return CorrectionBasis.MEDIAN

    @property
    def is_applied(self) -> bool:
        """Whether a measurement may be corrected using this stratum.

        Only strata anchored on the one primary method. See the module docstring.
        """
        return self.reference_platform == PRIMARY_REFERENCE and self.loo is not None

    def interval_for(self, other_ccs: float) -> JackknifePlusInterval | None:
        """The prediction interval around one corrected value, or None where there is none."""
        if self.loo is None:
            return None
        band = self.loo.interval_for(other_ccs)
        return None if isinstance(band, str) else band

    @property
    def typical_other_ccs(self) -> float | None:
        """The median value on the other platform: what a per-stratum width is quoted AT.

        A jackknife+ interval has no single width - it is a property of the query - so any
        width quoted for a whole stratum has to name the value it was computed at, and
        this is that value.
        """
        if not self.stratum.points:
            return None
        return median([p.other_ccs for p in self.stratum.points])

    @property
    def interval_at_a_typical_ion(self) -> JackknifePlusInterval | None:
        """The interval at the stratum's median ion. For reporting a representative width."""
        typical = self.typical_other_ccs
        return None if typical is None else self.interval_for(typical)

    def slope_derived(self, other_ccs: float) -> float | None:
        """The Deming correction: invert slope and intercept."""
        slope = self.stratum.agreement.deming_slope
        intercept = self.stratum.agreement.deming_intercept
        if slope in (None, 0) or intercept is None:
            return None
        return (other_ccs - intercept) / slope

    def median_derived(self, other_ccs: float) -> float:
        """The constant-offset correction, applied multiplicatively.

        The stratum's median is a percentage OF THE REFERENCE, so inverting it is a
        division rather than a subtraction: other = reference * (1 + m/100).
        """
        return other_ccs / (1.0 + self.median_offset_percent / 100.0)

    def robust_derived(self, other_ccs: float) -> float | None:
        return None if self.robust is None else self.robust.correct(other_ccs)

    def headline(self, other_ccs: float) -> float | None:
        """The corrected value on the chosen basis, or None where that basis cannot compute."""
        if self.basis is CorrectionBasis.ROBUST_SLOPE:
            return self.robust_derived(other_ccs)
        return self.median_derived(other_ccs)

    def summary(self) -> list[str]:
        lines = [
            f"  {self.other_platform} referred to {self.reference_platform}   n = {self.n}",
            f"    calibration groups {self.stratum.reference_group}  |  {self.stratum.other_group}",
            f"    basis        {self.basis.value}"
            f"   {'APPLIED' if self.is_applied else 'reported only, never applied'}",
            f"    median offset {self.median_offset_percent:+.3f}%",
        ]
        if self.robust is not None:
            interval = self.robust.slope_interval
            shown = f" ({interval[0]:.4f}, {interval[1]:.4f})" if interval else " (no rank interval)"
            lines.append(
                f"    robust slope {self.robust.slope:.5f}{shown}"
                f"   intercept {self.robust.intercept:+.3f}"
                f"   differs from 1: {'yes' if self.robust.slope_distinguishable_from_unity else 'no'}"
            )
            lines.append(
                f"    Deming slope {self.stratum.agreement.deming_slope:.5f}"
                f"   (breakdown 0; robust breakdown {self.robust.breakdown_point:.3f},"
                f" {self.robust.ions_that_may_be_arbitrarily_bad} ions may be anywhere)"
            )
        band = self.interval_at_a_typical_ion
        if band is not None:
            lines.append(
                f"    interval at a typical ion ({band.at_value:.1f} A^2 on {self.other_platform}):"
                f" {band.low:.2f} to {band.high:.2f} A^2"
            )
            lines.append(
                f"                 width {band.width:.3f} A^2 ({100 * band.width / band.at_value:.2f}% of CCS)"
                f"   nominal {band.nominal_coverage:.0%}, GUARANTEED {band.guaranteed_coverage:.0%}"
            )
            lines.append(
                f"                 tail ratio {band.tail_ratio:.2f}"
                f"   {'SET BY ONE ION' if band.driven_by_one_ion else 'not set by one ion'}"
                f"   {'informative' if band.interval_is_informative else 'DEGENERATE'}"
            )
        for note in self.refusals:
            lines.append(f"    refused: {note}")
        for note in self.warnings:
            lines.append(f"    weak: {note}")
        lines.append(f"    scope        {self.scope.caveat()}")
        return lines


def _compound_of(point: PairedPoint) -> str:
    """The grouping unit for the leave-one-out split: the compound, never the ion."""
    return str(point.ion.key.analyte)


def jackknife_plus(points: Sequence[PairedPoint], alpha: float = CONFORMAL_ALPHA) -> LeaveOneOutFits | str:
    """The leave-one-compound-out refits an interval is built from, or why there are none.

    Returns `LeaveOneOutFits`, from which `interval_for(value)` gives the interval around
    one corrected value. It does NOT return a band: see LeaveOneOutFits.

    The construction, following Barber et al. 2021 rather than a simplification of it:
    for each group g, refit on every point not in g, and for each left-out point i in g
    record the leave-one-out residual R_i = |corrected_i - reference_i|. The interval
    around a new prediction is then taken from the quantiles of the leave-one-out
    predictions offset by their own residuals - which is what gives the 1-2alpha
    guarantee, as opposed to a single symmetric band around a full-data fit, which
    guarantees nothing.

    Here the correction is applied to the OTHER platform's value to recover the
    reference, so the residual is in the reference platform's units and the interval is an
    interval on the corrected value.
    """
    groups: dict[str, list[PairedPoint]] = defaultdict(list)
    for point in points:
        groups[_compound_of(point)].append(point)
    if len(groups) < MIN_POINTS_FOR_JACKKNIFE:
        return (
            f"leave-one-compound-out needs at least {MIN_POINTS_FOR_JACKKNIFE} compounds to leave one out and"
            f" still fit robustly, and this stratum has {len(groups)}"
        )

    slopes: list[float] = []
    intercepts: list[float] = []
    residuals: list[float] = []
    for held_out, held_points in groups.items():
        kept = [p for p in points if _compound_of(p) != held_out]
        fit = passing_bablok([p.reference_ccs for p in kept], [p.other_ccs for p in kept])
        if isinstance(fit, tuple):
            return f"a leave-one-out refit refused: {fit[1][0]}"
        for point in held_points:
            corrected = fit.correct(point.other_ccs)
            if corrected is None:
                return "a leave-one-out refit produced a zero slope, so its correction does not invert"
            # One (refit, residual) pair per LEFT-OUT POINT, so a compound contributing two
            # conformers contributes two - the refit is shared and the residuals are not.
            slopes.append(fit.slope)
            intercepts.append(fit.intercept)
            residuals.append(abs(corrected - point.reference_ccs))

    if not residuals:
        return "no leave-one-out residual could be computed"

    return LeaveOneOutFits(
        slopes=tuple(slopes),
        intercepts=tuple(intercepts),
        residuals=tuple(residuals),
        groups=len(groups),
        alpha=alpha,
    )


@dataclass(frozen=True)
class CoverageCheck:
    """Nested leave-one-compound-out coverage, and an honest account of its worth.

    READ `is_evidence_of_calibration` BEFORE READING `rate`. It is False, always, and the
    reason is arithmetic rather than a shortcoming of this corpus: the interval is built
    from the quantile of the residual set, so the fraction of points inside it is pinned at
    the quantile index over n. Measured over the nine applied strata, `rate` equals
    `pinned_rate` to the last digit in all nine.

    What the check DOES establish is that the interval is not grossly too narrow. Scaling
    every residual by one half takes coverage below the guarantee, so a broken interval
    fails this. A merely miscalibrated one passes.
    """

    n: int
    covered: int
    guaranteed_coverage: float
    nominal_coverage: float
    mean_width_percent: float

    @property
    def rate(self) -> float:
        return self.covered / self.n if self.n else 0.0

    @property
    def pinned_rate(self) -> float:
        """What the rate must be by construction: the quantile index over n."""
        return math.ceil(self.n * self.nominal_coverage) / self.n if self.n else 0.0

    @property
    def is_pinned(self) -> bool:
        """Whether the observed rate is the arithmetically forced one. True on every real stratum."""
        return abs(self.rate - self.pinned_rate) < 1e-12

    @property
    def is_evidence_of_calibration(self) -> bool:
        """Always False. Kept as a property so the answer is in the code, not only in prose."""
        return False

    @property
    def meets_guarantee(self) -> bool:
        return self.rate >= self.guaranteed_coverage

    def summary(self) -> list[str]:
        return [
            f"    coverage     {self.covered}/{self.n} = {self.rate:.3f}"
            f"   guarantee {self.guaranteed_coverage:.2f}"
            f"   {'meets it' if self.meets_guarantee else 'FAILS THE GUARANTEE'}",
            f"                 pinned at {self.pinned_rate:.3f} by the quantile index"
            f"   {'(identical: not evidence)' if self.is_pinned else '(differs: worth reading)'}",
            f"                 mean width {self.mean_width_percent:.2f}% of CCS",
        ]


def measure_coverage(
    stratum: StratumStatistics, alpha: float = CONFORMAL_ALPHA
) -> CoverageCheck | str:
    """Hold each compound out entirely, rebuild the interval, and count what it covers.

    Honest out-of-sample in the only sense available on one study: the interval for a
    compound is built without that compound. See CoverageCheck for why the number this
    produces is not evidence of calibration.
    """
    points = list(stratum.points)
    if not points:
        return "no points to measure coverage over"
    covered = 0
    widths: list[float] = []
    measured = 0
    for held in points:
        rest = [p for p in points if _compound_of(p) != _compound_of(held)]
        fits = jackknife_plus(rest, alpha)
        if isinstance(fits, str):
            return f"a held-out refit refused: {fits}"
        band = fits.interval_for(held.other_ccs)
        if isinstance(band, str):
            return f"a held-out interval refused: {band}"
        measured += 1
        if band.low <= held.reference_ccs <= band.high:
            covered += 1
        widths.append(100.0 * band.width / band.at_value)
    return CoverageCheck(
        n=measured,
        covered=covered,
        guaranteed_coverage=1.0 - 2.0 * alpha,
        nominal_coverage=1.0 - alpha,
        mean_width_percent=sum(widths) / len(widths),
    )


def fit_stratum(stratum: StratumStatistics, alpha: float = CONFORMAL_ALPHA) -> StratumCorrection:
    """Fit one stratum three ways. Never pooled with another stratum, by construction:
    this function cannot see another one."""
    records = [p.reference for p in stratum.points] + [p.other for p in stratum.points]
    scope = stamp_over(records)
    offset = stratum.agreement.median_difference_percent
    if offset is None:
        offset = median([p.difference_percent for p in stratum.points]) if stratum.points else 0.0

    warnings: list[str] = []
    refusals: list[str] = []

    fit = passing_bablok([p.reference_ccs for p in stratum.points], [p.other_ccs for p in stratum.points])
    robust: PassingBablokFit | None
    if isinstance(fit, tuple):
        robust = None
        robust_refusals = fit[1]
        refusals.append(NO_ROBUST_FIT.format(why="; ".join(fit[1])))
    else:
        robust = fit
        robust_refusals = ()
        warnings.extend(fit.warnings)

    fits = jackknife_plus(stratum.points, alpha)
    loo: LeaveOneOutFits | None
    if isinstance(fits, str):
        loo = None
        refusals.append(NO_INTERVAL.format(why=fits))
    else:
        loo = fits
        # Checked at the stratum's median ion: the readiness verdict and the tail ratio
        # depend only on the residual set, so they are the same at every query value.
        probe = fits.interval_for(median([p.other_ccs for p in stratum.points]))
        if isinstance(probe, str):
            refusals.append(NO_INTERVAL.format(why=probe))
            loo = None
        else:
            if probe.readiness_warning:
                warnings.append(probe.readiness_warning)
            if probe.readiness_refusal:
                refusals.append(probe.readiness_refusal)
            if probe.driven_by_one_ion:
                warnings.append(
                    f"the interval width is set by the single largest leave-one-out residual"
                    f" ({probe.largest_residual:.3f} A^2), so it describes that ion as much as this stratum"
                )

    if stratum.reference_platform != PRIMARY_REFERENCE:
        warnings.append(
            NOT_APPLIED_NOT_PRIMARY.format(
                reference=stratum.reference_platform, other=stratum.other_platform
            )
        )

    return StratumCorrection(
        stratum=stratum,
        scope=scope,
        median_offset_percent=offset,
        robust=robust,
        robust_refusals=robust_refusals,
        loo=loo,
        warnings=tuple(warnings),
        refusals=tuple(refusals),
    )


def _canonical(value: object) -> str:
    """One spelling per value, so the same model always hashes the same.

    Floats go through repr rather than format: repr round-trips exactly in Python 3, and a
    rounded digest would call two different fits identical.
    """
    if value is None:
        return "~"
    if isinstance(value, float):
        return repr(value)
    return str(value)


def _digest(lines: Sequence[str]) -> str:
    """sha256 over SORTED lines, so the digest does not depend on iteration order."""
    joined = "\n".join(sorted(lines))
    return hashlib.sha256(joined.encode("utf-8")).hexdigest()


def corpus_digest(records: Iterable[object]) -> str:
    """A digest of the measurements BEHIND THE FIT. Not of every file on disk.

    Over the RECORDS rather than the files, so it is meaningful however the model was built
    and so that two deployments reading the same data by different paths agree.

    WHAT IT DOES NOT COVER, said plainly because the name suggests otherwise: a record that
    feeds no correction is not in it. Loading a seed file whose measurements pair nothing -
    the Struwe files are one platform, so they pair nothing - leaves this digest unchanged,
    and so does deleting it. That is deliberate. The digest moves exactly when an ANSWER
    could move, which is the promise being made; hashing every file instead would report a
    change whenever any data moved, making two identical answers look like two models.
    """
    lines = []
    for record in records:
        analyte = getattr(record, "analyte", None)
        identity = analyte.identity_key() if analyte is not None else None
        lines.append(
            "|".join(
                _canonical(part)
                for part in (
                    identity,
                    getattr(record, "adduct", None),
                    getattr(record, "charge", None),
                    getattr(record, "drift_gas", None),
                    getattr(record, "ims_type", None),
                    getattr(record, "dtims_method", None),
                    getattr(record, "calibrant", None),
                    getattr(record, "ccs", None),
                    getattr(record, "ccs_uncertainty", None),
                    getattr(record, "uncertainty_type", None),
                    getattr(record, "conformer", None),
                    getattr(record, "source", None),
                    getattr(record, "doi", None),
                )
            )
        )
    return _digest(lines)


def parameters_digest(corrections: Sequence["StratumCorrection"], alpha: float) -> str:
    """A digest of everything that decides what a correction returns.

    Includes the leave-one-out slopes, intercepts and residuals, because those set every
    interval - a digest covering only the headline slope would call two models identical
    while their intervals differed.
    """
    lines = [f"alpha|{_canonical(alpha)}"]
    for correction in corrections:
        robust = correction.robust
        agreement = correction.stratum.agreement
        loo = correction.loo
        lines.append(
            "|".join(
                _canonical(part)
                for part in (
                    correction.reference_platform,
                    correction.other_platform,
                    correction.stratum.reference_group,
                    correction.stratum.other_group,
                    correction.n,
                    correction.basis.value,
                    correction.is_applied,
                    correction.median_offset_percent,
                    None if robust is None else robust.slope,
                    None if robust is None else robust.intercept,
                    None if robust is None else robust.slope_interval,
                    agreement.deming_slope,
                    agreement.deming_intercept,
                    None if loo is None else loo.slopes,
                    None if loo is None else loo.intercepts,
                    None if loo is None else loo.residuals,
                )
            )
        )
    return _digest(lines)


@dataclass(frozen=True)
class ModelFingerprint:
    """What this model is, so that two answers can be told apart or told the same.

    Carried on `/health` and on every harmonized estimate. An answer without one cannot be
    reproduced: the model is refitted at every startup, and nothing else records that the
    data or the code behind it moved.

    TWO DIGESTS, AND THE PAIR IS THE USEFUL PART:

      corpus      over the records BEHIND THE FIT - not the files they came from, and not
                  every record loaded. A measurement that feeds no correction is not in it,
                  so data that cannot change an answer does not change this. It moves
                  exactly when an answer could.
      parameters  over everything that decides an answer: which stratum applies, on what
                  basis, with which slopes, intercepts and offsets, and the leave-one-out
                  residuals that set every interval.

    A different corpus with the same parameters means the data moved without moving the fit.
    The same corpus with different parameters means the code did. One combined hash would
    say only that something had.
    """

    corpus: str
    parameters: str

    @property
    def short(self) -> str:
        """Twelve hex characters of each, for logs and for a human comparing two responses."""
        return f"{self.corpus[:12]}/{self.parameters[:12]}"

    def same_corpus_as(self, other: "ModelFingerprint") -> bool:
        return self.corpus == other.corpus

    def same_parameters_as(self, other: "ModelFingerprint") -> bool:
        return self.parameters == other.parameters


@dataclass(frozen=True)
class HarmonizationModel:
    """Every stratum correction, and the scope of the corpus they were fitted on."""

    corrections: tuple[StratumCorrection, ...]
    scope: ScopeStamp
    fingerprint: ModelFingerprint
    alpha: float = CONFORMAL_ALPHA

    @property
    def applied(self) -> tuple[StratumCorrection, ...]:
        """The corrections that may actually be applied to a measurement."""
        return tuple(c for c in self.corrections if c.is_applied)

    @property
    def diagnostics(self) -> tuple[StratumCorrection, ...]:
        """Computed, reported, never applied."""
        return tuple(c for c in self.corrections if not c.is_applied)

    @property
    def maturity(self) -> MaturityStamp:
        """Always provisional while the corpus is one study.

        `DataMaturity.VALIDATED` is not reachable from here: it would require held-out
        matched ions from a source the model was not fitted on, and a within-study scope
        means by definition that there is no such source.
        """
        from .readiness import DataMaturity

        ions = len({id(p.ion) for c in self.corrections for p in c.stratum.points})
        return MaturityStamp(data_maturity=DataMaturity.PROVISIONAL, matched_ion_count=ions)

    @property
    def basis_counts(self) -> Mapping[str, int]:
        counts: dict[str, int] = defaultdict(int)
        for correction in self.corrections:
            counts[correction.basis.value] += 1
        return dict(counts)

    def correction_for(self, platform: str, reference_group: str, other_group: str) -> StratumCorrection | None:
        """The one correction that applies to a measurement, or None.

        Matched on the calibration groups rather than the platform alone, because two
        adducts of one platform are two different comparisons and pooling them is what
        this whole package refuses.
        """
        for correction in self.applied:
            if (
                correction.other_platform == platform
                and correction.stratum.other_group == other_group
                and correction.stratum.reference_group == reference_group
            ):
                return correction
        return None

    def summary(self) -> str:
        lines = [
            "Harmonization model",
            f"  strata fitted              {len(self.corrections)}",
            f"  applied (anchored on {PRIMARY_REFERENCE})  {len(self.applied)}",
            f"  reported only              {len(self.diagnostics)}",
            f"  headline basis             {self.basis_counts}",
            f"  maturity                   {self.maturity.data_maturity.value},"
            f" {self.maturity.matched_ion_count} matched ions",
            f"  fingerprint                {self.fingerprint.short}"
            f"   (corpus/parameters, sha256)",
            f"  SCOPE                      {self.scope.caveat()}",
            f"  {COVERAGE_IS_NOT_EVIDENCE}",
            "",
        ]
        for correction in self.corrections:
            lines += correction.summary()
            lines.append("")
        return "\n".join(lines)


def fit_harmonization(comparison: ComparisonReport, alpha: float = CONFORMAL_ALPHA) -> HarmonizationModel:
    """Fit every stratum in a comparison report.

    Takes the M3 report rather than raw records, so the licence gate, the matched-ion
    construction and the no-pooling rule have all already been applied and this cannot
    bypass them.

    Note on synthetic corpora: this deliberately does NOT refuse them. `compare_platforms`
    computes over synthetic data on purpose - the arithmetic has to be exercisable - and
    refuses at the point a figure is quoted. Putting the synthetic gate here would make
    every fixture unfittable and force an `allow_synthetic` parameter, which is the
    "floor that can be passed in" this package refuses elsewhere. The gate is on the
    publication path: see `assert_may_be_published`.
    """
    corrections = tuple(
        fit_stratum(stratum, alpha) for pair in comparison.pairs for stratum in pair.strata
    )
    if corrections:
        scope = stamp_over_stamps([c.scope for c in corrections])
    else:
        scope = ScopeStamp(studies=(), instruments=(), platforms=(), records_behind_it=0)
    # Over the DISTINCT records behind the fit. A record paired into several strata is one
    # measurement and is counted once, or the digest would depend on how many comparisons
    # happened to use it.
    records = {
        id(record): record
        for correction in corrections
        for point in correction.stratum.points
        for record in (point.reference, point.other)
    }
    fingerprint = ModelFingerprint(
        corpus=corpus_digest(records.values()),
        parameters=parameters_digest(corrections, alpha),
    )
    return HarmonizationModel(
        corrections=corrections, scope=scope, fingerprint=fingerprint, alpha=alpha
    )


def assert_may_be_published(model: HarmonizationModel, claim: Claim) -> None:
    """The publication gate. Raises unless the corpus supports the claim being made."""
    assert_may_be_quoted_as(model.scope, claim)


NO_CORRECTION_FOR_THIS_PLATFORM = (
    "no correction is applied for {platform} in calibration group {group}. Either this platform was not"
    " compared against the primary method ({primary}) in this corpus, or the comparison exists and is"
    " reported without being applied. A value is not corrected on a stratum that anchors nothing"
)
ORIGINAL_IS_THE_PRIMARY = (
    "this measurement is already on the primary method ({primary}), so there is nothing to refer it to."
    " It is returned unchanged, which is the correct answer rather than a failure"
)


@dataclass(frozen=True)
class Harmonized:
    """One measurement, corrected, with the ORIGINAL beside it and never in place of it.

    `original` is the object that was passed in, by identity. Nothing in this module
    constructs a modified copy of a measurement, so there is no code path by which a
    caller could receive a changed one.
    """

    original: object
    correction: StratumCorrection | None
    harmonized_ccs: float | None
    slope_derived_ccs: float | None
    median_derived_ccs: float | None
    robust_derived_ccs: float | None
    basis: CorrectionBasis | None
    interval: JackknifePlusInterval | None
    confidence: Confidence | None
    scope: ScopeStamp | None
    refusals: tuple[str, ...] = ()

    @property
    def was_corrected(self) -> bool:
        return self.harmonized_ccs is not None

    @property
    def original_ccs(self) -> float:
        """The value as submitted. Read from the original object every time, never cached."""
        return float(self.original.ccs)


def harmonize(model: HarmonizationModel, measurement: object) -> Harmonized:
    """Refer one measurement to the primary platform, or say why it was not.

    THE ORIGINAL IS NEVER MODIFIED. It is carried through by reference and returned; this
    function has no code that writes to it, and `Harmonized` is frozen.
    """
    platform = _platform_of_record(measurement)
    group = str(getattr(measurement, "calibration_group", ""))

    if platform == PRIMARY_REFERENCE:
        return Harmonized(
            original=measurement,
            correction=None,
            harmonized_ccs=None,
            slope_derived_ccs=None,
            median_derived_ccs=None,
            robust_derived_ccs=None,
            basis=None,
            interval=None,
            confidence=None,
            scope=model.scope,
            refusals=(ORIGINAL_IS_THE_PRIMARY.format(primary=PRIMARY_REFERENCE),),
        )

    correction = None
    for candidate in model.applied:
        if candidate.other_platform == platform and candidate.stratum.other_group == group:
            correction = candidate
            break
    if correction is None:
        return Harmonized(
            original=measurement,
            correction=None,
            harmonized_ccs=None,
            slope_derived_ccs=None,
            median_derived_ccs=None,
            robust_derived_ccs=None,
            basis=None,
            interval=None,
            confidence=None,
            scope=model.scope,
            refusals=(
                NO_CORRECTION_FOR_THIS_PLATFORM.format(
                    platform=platform, group=group or "(none)", primary=PRIMARY_REFERENCE
                ),
            ),
        )

    observed = float(measurement.ccs)
    headline = correction.headline(observed)
    band = correction.interval_for(observed)
    grade = grade_correction(observed, getattr(measurement, "matched_ion_key", None), correction.stratum)
    refusals: list[str] = []
    if headline is None:
        refusals.append("the chosen basis could not produce a value")
    if band is None:
        refusals.append("no prediction interval could be formed for this value")
    return Harmonized(
        original=measurement,
        correction=correction,
        harmonized_ccs=headline if band is not None else None,
        slope_derived_ccs=correction.slope_derived(observed),
        median_derived_ccs=correction.median_derived(observed),
        robust_derived_ccs=correction.robust_derived(observed),
        basis=correction.basis,
        interval=band,
        confidence=grade,
        scope=correction.scope,
        refusals=tuple(refusals),
    )


def _platform_of_record(record: object) -> str:
    ims = getattr(record, "ims_type", None)
    method = getattr(record, "dtims_method", None)
    if ims is None:
        return ""
    return f"{ims}/{method}" if method else str(ims)
