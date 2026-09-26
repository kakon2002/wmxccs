"""What must be true before a model may be fitted, and what to report while it is not.

Nothing here fits anything. There are no CCS training records, so every path to
a fit refuses, and the refusal is the deliverable.

ONE rule refuses here: the result could not be evaluated at all.

Everything else warns. A fit on forty records goes ahead and is reported as what
it is, with the record count travelling in the maturity stamp beside every
number. A weak result that reports its own weakness is useful; a confident
interval drawn from a calibration set too small to support it is not.

So MIN_TRAINING_RECORDS is a WARNING threshold, not a refusal. The brief's
earlier wording said to stop under about two hundred records; the owner has
superseded that with the inevaluability rule, and this module follows the later
framing. Reinstating it as a blocker would refuse the corpus we actually expect
to get, which is when the guard would matter most.

A blocker means one of these, and nothing else: there is nothing to fit; the
corpus cannot be divided into grouped folds at all; or no calibration group
holds enough to be scored within. See also the derived conformal floor below,
which is the sharpest case of an interval that cannot be formed rather than one
that is merely wide.

The sharpest case is the conformal calibration set, and its floor is derivable
rather than chosen. Split conformal takes the ceil((n+1)(1-alpha))-th smallest
conformity score as its quantile. When that index exceeds n there is no such
score, so no finite interval exists and any number printed in its place is
invented. The smallest workable calibration set is therefore ceil(1/alpha) - 1:
nine for 90 per cent coverage, nineteen for 95 per cent. At exactly that size
the quantile is the largest observed score, so the interval spans the full
observed range: valid, and uninformative. That is weakness, not impossibility,
so it is reported rather than refused.

Floors are module constants and never function parameters. A floor that can be
passed in is a floor that can be passed a zero.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

from .contracts import DataMaturity, MaturityStamp
from .evaluation import MIN_GROUP_RECORDS, MIN_GROUP_STRUCTURES
from .features import has_structure
from .licensing import (
    LicenceGateError,
    TrainingGateError,
    UnbackedClaimError,
    assert_trainable,
    declares_synthetic,
)
from .models import CalibrationGroup
from .splits import MIN_ANALYTE_GROUPS, AnalyteLevel, analyte_key, group_records, provenance_atoms

# A WARNING threshold, not a blocker. The brief's original rule was to stop under
# about two hundred usable records; the owner has since superseded it, because a
# forty-record corpus is perfectly evaluable and refusing it would suppress a
# weak but honest result rather than a misleading one.
MIN_TRAINING_RECORDS = 200
# POLICY, chosen rather than measured. With five folds this puts roughly twenty
# glycans in a held-out fold; whether twenty suffices to resolve the difference
# this platform exists to resolve is itself unknown until real records exist.
TARGET_ANALYTE_GROUPS = 100
# Two before leave-one-source-out produces anything at all, three before any
# cross-study claim is worth making.
MIN_SOURCES_FOR_CROSS_STUDY = 3

# Nominal coverage the conformal interval is quoted at, and the calibration size
# it requires. Both derived from the quantile index, not chosen: see the module
# docstring. Recomputed here rather than hard-coded so the derivation is the code.
CONFORMAL_ALPHA = 0.10


def smallest_calibration_set(alpha: float = CONFORMAL_ALPHA) -> int:
    """The fewest calibration points at which the nominal coverage is attainable.

    Below this, ceil((n+1)(1-alpha)) exceeds n, no finite conformal quantile
    exists, and an interval quoted at that coverage is fabricated.
    """
    if not 0 < alpha < 1:
        raise ValueError(f"alpha must lie strictly between 0 and 1, not {alpha!r}")
    return math.ceil(1 / alpha) - 1


MIN_CALIBRATION_RECORDS = smallest_calibration_set()

NOTHING_TO_FIT = "there is nothing to fit: the training set holds no records"
NO_RESOLVED_STRUCTURES = (
    "none of the {held} records carries a resolved structure, so a structure-feature fit has nothing to"
    " learn from: a composition-only measurement says nothing about isomers. This is a curation gap, not"
    " a licence one - a published CCS row has to be tied to a specific structure by hand"
)
TOO_FEW_RECORDS = "the training set holds {held} records; fitting needs at least {needed}"
TOO_FEW_GROUPS_TO_FIT = "the training set covers {groups} analyte groups; fitting needs at least {needed}"
UNKNOWN_BASELINE = "{baseline!r} is not a known baseline; the registry holds {known}"
NOT_A_TRAINING_SET = (
    "fitting takes a TrainingSet, not a {given}. A bare sequence has not been through the licence gate,"
    " so it may hold records that must never train"
)
SYNTHETIC_FIXTURES = (
    "{synthetic} of {held} records are synthetic fixtures, declared as such by the record or by a structure"
    " inside it: this set is a test, and no number from it describes real data"
)
SYNTHETIC_IN_FIT = (
    "the training set holds {synthetic} synthetic fixture(s) of {held} records, counting a declaration on a"
    " structure inside a record: a fit on them is a test, not a model. A test that means to reach the fit"
    " passes synthetic_ok=True; nothing else should"
)
CALIBRATION_TOO_SMALL = (
    "a calibration set of {held} cannot support {coverage:.0%} coverage: the conformal quantile would need"
    " the {index} smallest score of {held}, which does not exist, so no finite interval can be formed."
    " At least {needed} are required"
)
CALIBRATION_DEGENERATE = (
    "the calibration set holds {held}, the fewest at which {coverage:.0%} coverage is attainable, so the"
    " interval is the full observed range: honest, and uninformative"
)

BASELINES = ("group_median", "mass_power_law", "gradient_boosted")


class InsufficientTrainingDataError(TrainingGateError):
    """There is not enough data to fit, or not enough to evaluate what was fitted.

    A subclass of the training gate rather than of the licence gate: a data gap
    is not a licence problem, and conflating them would let one be reported as
    the other.
    """


@dataclass(frozen=True)
class GroupCount:
    """How much one calibration group holds."""

    records: int
    structures: int

    @property
    def scorable(self) -> bool:
        return self.records >= MIN_GROUP_RECORDS and self.structures >= MIN_GROUP_STRUCTURES


def calibration_warning(held: int, alpha: float = CONFORMAL_ALPHA) -> str | None:
    """A warning where a conformal interval would be valid but uninformative, or None."""
    needed = smallest_calibration_set(alpha)
    if held == needed:
        return CALIBRATION_DEGENERATE.format(held=held, coverage=1 - alpha)
    return None


def _ordinal(number: int) -> str:
    """1st, 2nd, 3rd, 4th, 11th. Written out because "1th" in a refusal reads as a bug."""
    if 10 <= number % 100 <= 20:
        return f"{number}th"
    return f"{number}{ {1: 'st', 2: 'nd', 3: 'rd'}.get(number % 10, 'th') }"


def calibration_refusal(held: int, alpha: float = CONFORMAL_ALPHA) -> str | None:
    """Why a conformal interval cannot be formed at this calibration size, or None."""
    if held < 0:
        raise ValueError(f"a calibration set cannot hold {held} records")
    needed = smallest_calibration_set(alpha)
    if held < needed:
        return CALIBRATION_TOO_SMALL.format(
            held=held, coverage=1 - alpha, index=_ordinal(math.ceil((held + 1) * (1 - alpha))), needed=needed
        )
    return None


@dataclass(frozen=True)
class Readiness:
    """What is held, what cleared the gate, and every reason a fit is not yet possible.

    Producible at zero records, and prints an honest zero rather than an empty
    table. This is the report this milestone can actually deliver.
    """

    records_held: int
    records_cleared: int
    analyte_groups: int
    distinct_sources: int
    level: AnalyteLevel
    # Counted apart from everything else, because being blocked on curation is a
    # different problem from being blocked on a licence and is fixed differently.
    records_with_structure: int = 0
    by_group: Mapping[str, GroupCount] = field(default_factory=dict)
    # Each record's OWN status, so this tally sums to records_cleared. A status
    # on a structure INSIDE a record is not here; it is in synthetic_records,
    # which asks the whole record rather than its outermost layer.
    by_reuse_status: Mapping[str, int] = field(default_factory=dict)
    # Records that declare themselves invented, at any depth. Counted separately
    # because a declaration on a component would otherwise be counted nowhere,
    # and a set of fixtures would print as data.
    synthetic_records: int = 0

    @property
    def blockers(self) -> tuple[str, ...]:
        """Reasons the result could not be evaluated at all. These, and only these, refuse."""
        problems: list[str] = []
        if self.records_cleared == 0:
            problems.append(NOTHING_TO_FIT)
            return tuple(problems)
        if self.analyte_groups < MIN_ANALYTE_GROUPS:
            # Too few groups to divide into folds, so there is nothing to hold out.
            problems.append(TOO_FEW_GROUPS_TO_FIT.format(groups=self.analyte_groups, needed=MIN_ANALYTE_GROUPS))
        if not any(count.scorable for count in self.by_group.values()):
            problems.append(
                f"no calibration group holds {MIN_GROUP_RECORDS} measurements over"
                f" {MIN_GROUP_STRUCTURES} analyte groups, so nothing could be scored within one"
            )
        if self.records_with_structure == 0:
            problems.append(NO_RESOLVED_STRUCTURES.format(held=self.records_cleared))
        return tuple(problems)

    @property
    def records_without_structure(self) -> int:
        """Cleared the licence gate, but carry only a composition."""
        return self.records_cleared - self.records_with_structure

    @property
    def warnings(self) -> tuple[str, ...]:
        """Reasons the result will be weak. A weak result that says so is still useful."""
        notes: list[str] = []
        if self.synthetic_records:
            notes.append(SYNTHETIC_FIXTURES.format(synthetic=self.synthetic_records, held=self.records_cleared))
        if 0 < self.records_cleared < MIN_TRAINING_RECORDS:
            notes.append(TOO_FEW_RECORDS.format(held=self.records_cleared, needed=MIN_TRAINING_RECORDS))
        if self.records_cleared and self.analyte_groups < TARGET_ANALYTE_GROUPS:
            notes.append(
                f"{self.analyte_groups} analyte groups is below the {TARGET_ANALYTE_GROUPS} at which the"
                " number is worth quoting, so expect wide intervals"
            )
        if self.records_cleared and self.distinct_sources < MIN_SOURCES_FOR_CROSS_STUDY:
            notes.append(
                f"{self.distinct_sources} source(s): under {MIN_SOURCES_FOR_CROSS_STUDY} no cross-study claim"
                " can be made, so every number is a within-study number"
            )
        if self.records_without_structure and self.records_with_structure:
            notes.append(
                f"{self.records_without_structure} of {self.records_cleared} records carry only a composition"
                " and are refused from a structure-feature fit; they are blocked on curation, not on licence"
            )
        return tuple(notes)

    @property
    def shortfalls(self) -> tuple[str, ...]:
        """Everything worth saying, blockers first. Kept for callers that want one list."""
        return self.blockers + self.warnings

    @property
    def ready(self) -> bool:
        """Whether a fit may proceed. Weakness does not stop it; inevaluability does."""
        return not self.blockers

    @property
    def refusal(self) -> str | None:
        return "; ".join(self.blockers) or None

    @property
    def maturity(self) -> MaturityStamp:
        # Always provisional here. The validated flag belongs to the milestone
        # that checks a trained model against held-out measurements.
        return MaturityStamp(data_maturity=DataMaturity.PROVISIONAL, training_record_count=self.records_cleared)

    def summary(self) -> str:
        lines = [
            f"CCS training readiness, grouped at {self.level} level"
            f"  [{self.maturity.data_maturity}, {self.maturity.training_record_count} records behind any number]",
            f"  records held               {self.records_held}",
            *(
                [
                    f"  SYNTHETIC FIXTURES         {self.synthetic_records}  (the record, or a structure in it,"
                    " declared in code by a test; nothing here describes real data)"
                ]
                if self.synthetic_records
                else []
            ),
            f"  cleared the licence gate   {self.records_cleared}  (of {MIN_TRAINING_RECORDS} needed)",
            f"  carry a resolved structure {self.records_with_structure}"
            f"  ({self.records_without_structure} composition-only, blocked on curation not licence)",
            f"  analyte groups             {self.analyte_groups}  (of {MIN_ANALYTE_GROUPS} needed,"
            f" {TARGET_ANALYTE_GROUPS} before the number is worth quoting)",
            f"  distinct sources           {self.distinct_sources}"
            f"  (of {MIN_SOURCES_FOR_CROSS_STUDY} needed for any cross-study claim)",
            f"  conformal calibration      needs at least {MIN_CALIBRATION_RECORDS} held-out records for"
            f" {1 - CONFORMAL_ALPHA:.0%} coverage to be attainable at all",
        ]
        if self.by_reuse_status:
            lines.append("  by reuse status:")
            for status, count in sorted(self.by_reuse_status.items()):
                lines.append(f"    {status}: {count}")
        if self.by_group:
            lines.append("  by calibration group:")
            for group, count in sorted(self.by_group.items()):
                mark = "scorable" if count.scorable else "too small to score"
                lines.append(f"    {group}: {count.records} records, {count.structures} analytes ({mark})")
        lines.append("  ready to fit: " + ("yes" if self.ready else "no"))
        for blocker in self.blockers:
            lines.append(f"    cannot be evaluated: {blocker}")
        for warning in self.warnings:
            lines.append(f"    will be weak: {warning}")
        return "\n".join(lines)


@dataclass(frozen=True)
class TrainingSet:
    """Records that have been through the licence gate, and may be trained on.

    The gate runs in __post_init__, so no instance can hold an ungated record
    however it was constructed. There is no classmethod to route around.
    """

    records: tuple[object, ...] = ()
    level: AnalyteLevel = AnalyteLevel.COMPOSITION

    def __post_init__(self) -> None:
        object.__setattr__(self, "records", tuple(self.records))
        faults: list[str] = []
        licence_fault = False
        unbacked_only = True
        for index, record in enumerate(self.records):
            try:
                assert_trainable(record)
            except TrainingGateError as exc:
                faults.append(f"record {index}: {exc}")
                licence_fault = licence_fault or isinstance(exc, LicenceGateError)
                unbacked_only = unbacked_only and isinstance(exc, UnbackedClaimError)
        if faults:
            problem = "; ".join(faults)
            # The most specific class the faults allow: a caller that catches
            # UnbackedClaimError around a TrainingSet must not be told "licence
            # not cleared" for a record whose licence is fine and whose CLAIM is
            # what nothing backs, because that points at the wrong remedy.
            error = UnbackedClaimError if unbacked_only else LicenceGateError if licence_fault else TrainingGateError
            raise error(f"a training set cannot hold {len(faults)} record(s) that may not train: {problem}")

    def __len__(self) -> int:
        return len(self.records)

    @property
    def readiness(self) -> Readiness:
        by_group: dict[str, GroupCount] = {}
        analytes: dict[str, set[str]] = {}
        statuses: dict[str, int] = {}
        for record in self.records:
            group = str(getattr(record, "calibration_group", "unknown"))
            analytes.setdefault(group, set()).add(analyte_key(record, self.level))
            existing = by_group.get(group)
            by_group[group] = GroupCount(
                records=(existing.records if existing else 0) + 1,
                structures=len(analytes[group]),
            )
            status = str(getattr(record, "reuse_status", "unknown"))
            statuses[status] = statuses.get(status, 0) + 1
        # Studies, not source strings. Twenty-four rows transcribed from one table
        # are one source however their source text varies row to row, and the
        # DOI-preferring provenance atom is what the splitter groups by, so the
        # readiness figure and the split agree on what a source is.
        sources = {atom for record in self.records for atom in provenance_atoms(record)}
        # Counted the way the splitter counts, as connected components rather
        # than as distinct keys. The two diverge exactly where an identifier
        # merges two analytes, and reporting keys here would certify a corpus
        # ready that grouped_folds then refuses as too few groups.
        analyte_groups = len(group_records(self.records, level=self.level).groups) if self.records else 0
        return Readiness(
            records_held=len(self.records),
            # Every record here has already cleared the gate, by construction.
            records_cleared=len(self.records),
            records_with_structure=sum(1 for record in self.records if has_structure(record)),
            analyte_groups=analyte_groups,
            distinct_sources=len(sources),
            level=self.level,
            by_group=by_group,
            by_reuse_status=statuses,
            # Asked of the whole record, not of its outermost status: a structure
            # invented by a test inside a measurement claiming a real status is
            # still a fixture, and a set of them is still a test.
            synthetic_records=sum(1 for record in self.records if declares_synthetic(record)),
        )


def fit_ccs_baseline(training_set: TrainingSet, baseline: str = "gradient_boosted", *, synthetic_ok: bool = False):
    """Fit a CCS baseline. Refuses today, because there are no records to fit on.

    Every check below is an explicit raise rather than an assert, so it survives
    python -O, and every one sits lexically above the estimator import, so no
    empty dataset can reach it.

    A set holding synthetic fixtures is refused unless `synthetic_ok` is passed:
    the gate lets a fixture through so that a test can build a training set,
    but a fit on fixtures is a test and not a model, and the call must say so.
    """
    if not isinstance(training_set, TrainingSet):
        # A bare list is the natural accidental bypass, and a TypeError would be
        # swallowed by a caller skipping malformed rows.
        raise TrainingGateError(NOT_A_TRAINING_SET.format(given=type(training_set).__name__))
    readiness = training_set.readiness
    if readiness.refusal is not None:
        raise InsufficientTrainingDataError(readiness.refusal)
    if baseline not in BASELINES:
        raise InsufficientTrainingDataError(UNKNOWN_BASELINE.format(baseline=baseline, known=", ".join(BASELINES)))
    if readiness.synthetic_records and not synthetic_ok:
        raise TrainingGateError(SYNTHETIC_IN_FIT.format(synthetic=readiness.synthetic_records, held=len(training_set)))

    from sklearn.ensemble import HistGradientBoostingRegressor  # noqa: F401  (below every guard, deliberately)

    raise NotImplementedError(
        "fitting continues here the day readiness clears; no model has ever been trained in this repository"
    )
