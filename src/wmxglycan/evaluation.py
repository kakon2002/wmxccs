"""How a CCS prediction is scored, and the report that carries the score.

The maths here is stdlib only and takes plain numbers, so it can be exercised
in full against hand-computed arrays while no model and no measurement exist.

Two rules shape the whole module.

A CCS value means nothing outside its calibration group. Values measured on
different platforms, in different gases, against different calibrants or on
differently labelled analytes are not the same quantity, so an absolute error
in square angstrom cannot be averaged across them. Per-group results carry
absolute errors; the pooled result structurally cannot, because PooledError has
no square-angstrom field at all. The absence is the enforcement.

A number without its record count, and without the maturity of the data behind
it, is a number that will be quoted out of context. Both travel as fields of
the report rather than as a footnote.

Nothing here fits anything, and no report this milestone produces is ever
stamped validated. That flag belongs to the milestone that checks a trained
model against held-out measurements, and there is no such model.
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

from .contracts import DataMaturity, MaturityStamp
from .models import CalibrationGroup

from .splits import SplitMode

if TYPE_CHECKING:  # imported for typing only, so the dependency arrow stays one-way
    from .splits import AnalyteLevel, SplitReport
    from .training import TrainingSet

# Before a calibration group may be scored at all. They live here, beside the
# code that applies them, because training.py imports evaluation and not the
# other way round. Both are POLICY, chosen rather than measured: with fewer than
# twenty measurements a median relative error is a statement about four or five
# glycans, and with fewer than ten distinct analytes it largely measures how
# often one glycan was remeasured.
MIN_GROUP_RECORDS = 20
MIN_GROUP_STRUCTURES = 10

# Refusal wording is an interface: written once here, asserted verbatim in tests.
NO_PAIRS = "there are no observed and predicted pairs to score"
LENGTH_MISMATCH = "there are {observed} observed values and {predicted} predicted ones"
BAD_PREDICTION = "prediction at index {index} is {value!r}, which is not a positive finite number"
BAD_OBSERVATION = "observed value at index {index} is {value!r}, which is not a positive finite number"
GROUP_TOO_SMALL = (
    "this calibration group holds {records} measurements over {structures} analyte groups;"
    " scoring needs at least {min_records} over {min_structures}"
)
ISOMER_IS_NOT_PERFORMANCE = (
    "this result was measured in isomer-discrimination mode, which deliberately places isomers of one"
    " composition on both sides of the split. It answers whether any learnable CCS signal separates"
    " isomers at all, and it is an optimistic bound and a diagnostic. It is not deployment performance"
    " and must not be quoted as one: re-measure in deployment mode for a number that may be reported"
)
NOTHING_TO_HEADLINE = "nothing was scored, so there is no figure to report"


class ModeConfusionError(Exception):
    """An isomer-discrimination result was asked to speak as deployment performance.

    Deliberately not a ValueError. A number that looks right is the dangerous
    kind, and a caller wrapping its reporting in a broad `except (ValueError,
    TypeError)` would swallow this refusal and print the diagnostic as though it
    were performance.
    """


def _positive_finite(value: float) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) and value > 0


@dataclass(frozen=True)
class Metrics:
    """How far predictions sat from measurements, within one calibration group."""

    n: int
    median_relative_error: float  # per cent; how the ion mobility literature reports CCS
    mean_absolute_error: float  # square angstrom
    rmse: float  # square angstrom
    within_1_percent: float  # per cent of predictions inside 1 per cent of the measurement
    within_2_percent: float

    @classmethod
    def of(cls, observed: Sequence[float], predicted: Sequence[float]) -> Metrics:
        """Score predictions against measurements.

        Nothing is clamped, skipped or absolute-valued into plausibility: a
        prediction that is not a positive finite number is a defect in whatever
        produced it, and it is named by index and value rather than quietly
        dropped, which would silently improve the score.

        The empty case is refused here, before statistics.median is reached.
        StatisticsError is a ValueError, so a caller's row-level except clause
        would swallow it and read the absence of data as a bad row.
        """
        if len(observed) != len(predicted):
            raise ValueError(LENGTH_MISMATCH.format(observed=len(observed), predicted=len(predicted)))
        if not observed:
            raise ValueError(NO_PAIRS)
        for index, value in enumerate(observed):
            if not _positive_finite(value):
                raise ValueError(BAD_OBSERVATION.format(index=index, value=value))
        for index, value in enumerate(predicted):
            if not _positive_finite(value):
                raise ValueError(BAD_PREDICTION.format(index=index, value=value))

        # Relative to the MEASUREMENT, not the prediction: predicting 200 for a
        # measured 100 is 100 per cent out, not 50.
        relative = [100.0 * abs(p - o) / o for o, p in zip(observed, predicted)]
        absolute = [abs(p - o) for o, p in zip(observed, predicted)]
        n = len(observed)
        return cls(
            n=n,
            # statistics.median, so an even count takes the mean of the middle two.
            # Pinned by a test, because this is the headline number.
            median_relative_error=statistics.median(relative),
            mean_absolute_error=sum(absolute) / n,
            rmse=math.sqrt(sum(value * value for value in absolute) / n),
            within_1_percent=100.0 * sum(1 for value in relative if value <= 1.0) / n,
            within_2_percent=100.0 * sum(1 for value in relative if value <= 2.0) / n,
        )

    def summary(self) -> str:
        return (
            f"n={self.n}  median relative error {self.median_relative_error:.2f}%"
            f"  MAE {self.mean_absolute_error:.2f} A^2  RMSE {self.rmse:.2f} A^2"
            # One decimal place, deliberately: 996 of 1000 is 99.6 per cent, and
            # rounding it to "100%" would be a claim of perfect coverage.
            f"  within 1% {self.within_1_percent:.1f}%  within 2% {self.within_2_percent:.1f}%"
        )


@dataclass(frozen=True)
class PooledError:
    """One figure across calibration groups, in the only units that survive pooling.

    There is deliberately no mean_absolute_error and no rmse here. A relative
    error is dimensionless and comparable between a He-referenced TWIMS value
    and an N2 drift-tube one; a square-angstrom error is not, because the values
    themselves are not comparable. Adding such a field would make an incorrect
    number expressible, so the field does not exist.
    """

    n: int
    groups: int
    # NOT a median. Each group's median relative error, weighted by that group's
    # record count. A mean of medians is bounded by the extreme group medians, so
    # it cannot see the pooled tail: two groups of 30 at 0.1% and 31 at 10% give
    # 5.13 here against a true pooled median of 0.10. Named for what it is,
    # because the honest alternative would be to pool the per-record errors, and
    # that is a different statistic rather than a correction to this one.
    record_weighted_mean_of_group_medians: float
    # These two DO pool exactly: a record-weighted mean of per-group percentages
    # is the true pooled percentage.
    within_1_percent: float
    within_2_percent: float

    def summary(self) -> str:
        return (
            f"pooled over {self.groups} calibration groups, n={self.n}:"
            f" record-weighted mean of group medians"
            f" {self.record_weighted_mean_of_group_medians:.2f}%"
            f"  within 1% {self.within_1_percent:.1f}%  within 2% {self.within_2_percent:.1f}%"
        )


@dataclass(frozen=True)
class GroupResult:
    """One calibration group's result, or the reason it could not be scored."""

    group: CalibrationGroup
    n: int
    structures: int
    metrics: Metrics | None = None
    refusal: str | None = None

    def __post_init__(self) -> None:
        # A result that carries both a score and a refusal renders only the
        # score, so the refusal would be invisible; one carrying neither renders
        # the literal "not scored: None".
        if (self.metrics is None) == (self.refusal is None):
            raise ValueError("a group result carries either metrics or a refusal, never both and never neither")
        if self.metrics is not None and self.metrics.n != self.n:
            # Otherwise the count a reader checks and the count the summary
            # prints are different numbers.
            raise ValueError(f"this group says it holds {self.n} records, but its metrics scored {self.metrics.n}")
        if self.n < 0 or self.structures < 0:
            raise ValueError("a group cannot hold a negative number of records or analytes")

    @property
    def scored(self) -> bool:
        return self.metrics is not None

    def summary(self) -> str:
        head = f"{self.group}"
        if self.metrics is None:
            return f"{head}\n    not scored: {self.refusal}"
        return f"{head}\n    {self.metrics.summary()}  ({self.structures} analyte groups)"


def group_refusal(records: int, structures: int) -> str | None:
    """Why a calibration group cannot be scored, or None if it can."""
    if records < MIN_GROUP_RECORDS or structures < MIN_GROUP_STRUCTURES:
        return GROUP_TOO_SMALL.format(
            records=records,
            structures=structures,
            min_records=MIN_GROUP_RECORDS,
            min_structures=MIN_GROUP_STRUCTURES,
        )
    return None


def pool(results: Sequence[GroupResult]) -> PooledError | None:
    """Pool the scored groups, weighting each measurement equally. None if none scored."""
    scored = [result for result in results if result.metrics is not None]
    if not scored:
        return None
    n = sum(result.metrics.n for result in scored)
    weighted = sum(result.metrics.median_relative_error * result.metrics.n for result in scored)
    within1 = sum(result.metrics.within_1_percent * result.metrics.n for result in scored)
    within2 = sum(result.metrics.within_2_percent * result.metrics.n for result in scored)
    return PooledError(
        n=n,
        groups=len(scored),
        record_weighted_mean_of_group_medians=weighted / n,
        within_1_percent=within1 / n,
        within_2_percent=within2 / n,
    )


@dataclass(frozen=True)
class EvaluationReport:
    """What a baseline scored, under which split, on data of which maturity.

    analyte_level and maturity are fields rather than commentary: a number
    cannot be read off this report without the grouping level it was measured
    under and the maturity of the data behind it.
    """

    analyte_level: AnalyteLevel
    by_study: bool
    baseline: str
    records_in: int
    records_not_scored: int
    groups: tuple[GroupResult, ...]
    pooled: PooledError | None
    maturity: MaturityStamp
    split: SplitReport | None = None
    # Which question this result answers. Only a deployment result may be
    # reported as performance, and `headline` is the only way to get the figure.
    mode: SplitMode = SplitMode.DEPLOYMENT

    def __post_init__(self) -> None:
        if self.maturity.data_maturity is not DataMaturity.PROVISIONAL:
            # M6 assigns the validated flag, once a model has been checked
            # against held-out measurements. Nothing in M3 may set it.
            raise ValueError("an evaluation report from this milestone is provisional; nothing here validates a model")
        if self.records_in < 0 or self.records_not_scored < 0:
            raise ValueError("a report cannot hold a negative record count")
        if self.records_not_scored > self.records_in:
            raise ValueError(
                f"{self.records_not_scored} records went unscored out of {self.records_in} that went in"
            )
        scored = [result for result in self.groups if result.metrics is not None]
        if self.pooled is None and scored:
            # Otherwise the report prints "nothing was scored" directly above a
            # group carrying full metrics.
            raise ValueError(f"{len(scored)} group(s) were scored, so the pooled figure cannot be absent")
        if self.pooled is not None:
            if self.pooled.groups != len(scored):
                raise ValueError(f"the pooled figure claims {self.pooled.groups} groups; {len(scored)} were scored")
            if self.pooled.n != sum(result.metrics.n for result in scored):
                raise ValueError("the pooled record count disagrees with the groups it claims to pool")

    @property
    def headline(self) -> PooledError:
        """The figure that may be reported as performance, or a refusal.

        The single door to the number, so an isomer-mode result cannot reach a
        report by any path that does not pass this check.
        """
        if self.mode is not SplitMode.DEPLOYMENT:
            raise ModeConfusionError(ISOMER_IS_NOT_PERFORMANCE)
        if self.pooled is None:
            raise ModeConfusionError(NOTHING_TO_HEADLINE)
        return self.pooled

    @property
    def is_performance(self) -> bool:
        return self.mode is SplitMode.DEPLOYMENT

    def summary(self) -> str:
        lines = [
            f"baseline {self.baseline!r} in {self.mode} mode, grouped at {self.analyte_level} level"
            f"{', split across studies' if self.by_study else ''}"
            f"  [{self.maturity.data_maturity}, trained on {self.maturity.training_record_count} records]",
            f"  {self.records_in} records in, {self.records_not_scored} not scored",
        ]
        if self.mode is not SplitMode.DEPLOYMENT:
            # Never renders the pooled figure. A reader skimming for a number
            # must find a sentence about what this is instead of one to quote.
            lines.append("  DIAGNOSTIC, NOT PERFORMANCE. " + ISOMER_IS_NOT_PERFORMANCE)
            for result in self.groups:
                lines.append("  " + result.summary().replace("\n", "\n  "))
            return "\n".join(lines)
        lines.append(f"  {self.pooled.summary()}" if self.pooled else "  nothing was scored")
        for result in self.groups:
            lines.append("  " + result.summary().replace("\n", "\n  "))
        return "\n".join(lines)


def evaluate_baseline(
    training_set: TrainingSet, baseline: str, level: AnalyteLevel | None = None, *, synthetic_ok: bool = False
):
    """Score one baseline. Refuses today, because fitting it refuses today.

    The refusal comes from the training gate rather than from here, so there is
    one place that decides whether a fit may happen at all. `synthetic_ok` is
    passed through to it: a set of synthetic fixtures is scored only when the
    call says it is a test.
    """
    from .training import fit_ccs_baseline  # deferred: training imports evaluation, so this closes the cycle

    fit_ccs_baseline(training_set, baseline, synthetic_ok=synthetic_ok)  # raises; the next line would mean a model exists
    raise NotImplementedError("scoring continues here the day a baseline can be fitted")


def m3_report(training_set: TrainingSet) -> str:
    """The report this milestone can actually produce: readiness, not accuracy.

    No median relative error is produced, because none can be: there are no
    measurements to score against. What can be shown honestly is what is held,
    what cleared the licence gate, and what each shortfall is.
    """
    return training_set.readiness.summary()
