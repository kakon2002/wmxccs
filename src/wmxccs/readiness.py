"""What must be true before anything may be fitted, and what to report while it is not.

Nothing here fits anything. There are no cross-platform matched ions, so every
path to a fit refuses, and the refusal is the deliverable.

ONE CATEGORY REFUSES HERE: the result could not be evaluated at all. Everything
else warns. A fit on forty matched ions goes ahead and is reported as what it is,
with the count travelling in the maturity stamp beside every number. A weak
result that reports its own weakness is useful; a confident interval drawn from a
calibration set too small to support it is not.

WHAT WAS RETHRESHOLDED, AND WHY
-------------------------------
The glycan platform counted TRAINING RECORDS and ANALYTE GROUPS, because it was
fitting a model from structure to CCS. This platform compares platforms, so its
binding count is different and smaller in a way that matters:

- the unit is the MATCHED ION measured on MORE THAN ONE PLATFORM. A corpus of a
  hundred thousand values on one platform yields nothing at all. This is the
  count everything downstream rests on, and it has its own blocker.
- PLATFORM COUNT is a blocker of its own, and it is arithmetic rather than
  policy: with fewer than two platforms represented there are zero possible
  pairs, now and after any amount of further loading. The glycan platform had no
  such threshold because it never compared platforms, and inheriting its source
  count alone would have scored the current corpus as "2 sources, a warning"
  while the real problem - one platform, therefore no pair can ever exist - had
  nothing watching it at all.
- the per-platform-pair floor is derived from the fold design, not chosen: a
  grouped k-fold needs at least one group per fold, so a pair carrying fewer
  matched ions than there are folds cannot be cross-validated at all.

MIN_RECORDS_HELD stays a WARNING threshold, not a refusal, exactly as it was.
Refusing a small corpus would suppress a weak but honest result rather than a
misleading one, which is when the guard would matter least and cost most.

Floors are module constants and never function parameters. A floor that can be
passed in is a floor that can be passed a zero.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Mapping

from pydantic import BaseModel, ConfigDict, Field

from .licensing import (
    LicenceGateError,
    TrainingGateError,
    UnbackedClaimError,
    assert_trainable,
    declares_synthetic,
)

# --- the maturity stamp, ported from the glycan platform's contracts.py -------------------


class DataMaturity(StrEnum):
    """How far the basis of a number has been taken.

    M0 ships "provisional": no harmonization model has been fitted, let alone
    checked against held-out measurements. The stamp travels with every output so
    the caveat cannot be separated from the number.
    """

    PROVISIONAL = "provisional"
    VALIDATED = "validated"


class MaturityStamp(BaseModel):
    """The maturity of the data behind a number, carried in the number's own report.

    Frozen. A stamp is routinely shared between reports, so a mutable one could
    be flipped to validated after the check that refused it, and every report
    holding it would re-render the new claim.
    """

    model_config = ConfigDict(frozen=True)

    data_maturity: DataMaturity = Field(
        description="'provisional' until a harmonization model has been checked against held-out matched"
        " ions; 'validated' only after that."
    )
    matched_ion_count: int = Field(
        ge=0,
        description="Cross-platform matched ions behind any number quoted. 0 while none exist, which is why"
        " M0 is provisional.",
    )


# --- thresholds ---------------------------------------------------------------------------

# ARITHMETIC, not policy. A comparison needs two things to compare. Below this
# the matched-ion table is empty by construction and no amount of further data
# on the same platform changes it.
MIN_PLATFORMS = 2

# The fold design, ported from the glycan platform's splitter. Grouped k-fold:
# one analyte never spans train and test.
N_FOLDS = 5

# DERIVED from N_FOLDS, not chosen: a grouped k-fold needs at least one group per
# fold, so a platform pair carrying fewer matched ions than there are folds
# cannot be cross-validated at all. Below it, any error figure quoted for that
# pair is an in-sample figure wearing a cross-validated label.
MIN_MATCHED_IONS_PER_PAIR = N_FOLDS

# DERIVED: Deming regression fits a slope and an intercept, so three points leave
# one residual degree of freedom. At two, the line passes through both points and
# the residual is zero by construction, which reads as a perfect fit.
MIN_DEMING_POINTS = 3

# POLICY, chosen rather than measured, and carried over unchanged in spirit from
# the glycan platform's TARGET_ANALYTE_GROUPS. Whether a hundred matched ions
# suffice to characterise a platform pair is itself unknown until real matched
# ions exist.
TARGET_MATCHED_IONS = 100

# A WARNING threshold. Kept at the glycan platform's figure because the argument
# for it did not change: a forty-record corpus is perfectly evaluable, and
# refusing it would suppress an honest weak result.
MIN_RECORDS_HELD = 200

# Two before leave-one-source-out produces anything at all, three before any
# cross-study claim is worth making. Ported unchanged - but note it counts
# STUDIES, and for this platform the binding count is PLATFORMS, which has its
# own blocker above. Both are reported; neither stands in for the other.
MIN_SOURCES_FOR_CROSS_STUDY = 3

# Nominal coverage a conformal interval is quoted at, and the calibration size it
# requires. Both derived from the quantile index, not chosen. Recomputed here
# rather than hard-coded so the derivation is the code.
CONFORMAL_ALPHA = 0.10


def smallest_calibration_set(alpha: float = CONFORMAL_ALPHA) -> int:
    """The fewest calibration points at which the nominal coverage is attainable.

    Split conformal takes the ceil((n+1)(1-alpha))-th smallest conformity score
    as its quantile. When that index exceeds n there is no such score, so no
    finite interval exists and any number printed in its place is invented. The
    smallest workable calibration set is therefore ceil(1/alpha) - 1: nine for 90
    per cent coverage, nineteen for 95. At exactly that size the quantile is the
    largest observed score, so the interval spans the full observed range: valid,
    and uninformative.
    """
    if not 0 < alpha < 1:
        raise ValueError(f"alpha must lie strictly between 0 and 1, not {alpha!r}")
    return math.ceil(1 / alpha) - 1


MIN_CALIBRATION_RECORDS = smallest_calibration_set()

NOTHING_HELD = "there is nothing to work with: no record cleared the licence gate"
ONE_PLATFORM_ONLY = (
    "every cleared record comes from a single platform ({platforms}), so no ion is measured on two"
    " platforms and the matched-ion table is empty by construction. This is not a shortage of data."
    " More values on the same platform will not produce one pair: a second platform is required"
)
NO_MATCHED_IONS = (
    "{keys} distinct matched-ion keys are held across {platforms} platforms, but none of them is measured"
    " on more than one, so there is no pair to compare and no bias to estimate"
)
NO_EVALUABLE_PAIR = (
    "no platform pair holds {needed} matched ions, the fewest that {folds}-fold grouped cross-validation"
    " can divide, so no pair could be scored even if one were fitted"
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
NOT_A_MATCHED_ION_SET = (
    "this takes a MatchedIonSet, not a {given}. A bare sequence has not been through the licence gate,"
    " so it may hold records that must never train"
)


class InsufficientDataError(TrainingGateError):
    """There is not enough data to fit, or not enough to evaluate what was fitted.

    A subclass of the training gate rather than of the licence gate: a data gap
    is not a licence problem, and conflating them would let one be reported as
    the other.
    """


def _ordinal(number: int) -> str:
    """1st, 2nd, 3rd, 4th, 11th. Written out because "1th" in a refusal reads as a bug."""
    if 10 <= number % 100 <= 20:
        return f"{number}th"
    return f"{number}{ {1: 'st', 2: 'nd', 3: 'rd'}.get(number % 10, 'th') }"


def calibration_warning(held: int, alpha: float = CONFORMAL_ALPHA) -> str | None:
    """A warning where a conformal interval would be valid but uninformative, or None."""
    needed = smallest_calibration_set(alpha)
    if held == needed:
        return CALIBRATION_DEGENERATE.format(held=held, coverage=1 - alpha)
    return None


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


def interval_readiness(calibration_size: int, alpha: float = CONFORMAL_ALPHA) -> tuple[str | None, str | None]:
    """(refusal, warning) for quoting a conformal interval at this calibration size.

    THE ONE CALLER THAT MUST EXIST AND DOES NOT YET. In the glycan platform the
    two functions above were derived correctly, tested, and then never called by
    anything that refuses: the floor reached a human only as a descriptive line
    in a summary. That is the failure shape this repository keeps meeting - a
    guard that is present, correct and not wired to anything.

    It is wired here, and it is reachable and tested. The interval-quoting path
    itself does not exist until the harmonization model does, so the obligation
    is recorded in LIMITATIONS.md rather than left implicit: whatever computes a
    prediction interval calls this first, and refuses on a refusal.
    """
    return calibration_refusal(calibration_size, alpha), calibration_warning(calibration_size, alpha)


@dataclass(frozen=True)
class PairCount:
    """How many matched ions one platform pair holds, and whether it could be scored."""

    matched_ions: int

    @property
    def fittable(self) -> bool:
        """Enough points for a two-parameter fit to have a residual at all."""
        return self.matched_ions >= MIN_DEMING_POINTS

    @property
    def scorable(self) -> bool:
        """Enough matched ions to divide into grouped folds."""
        return self.matched_ions >= MIN_MATCHED_IONS_PER_PAIR


@dataclass(frozen=True)
class Readiness:
    """What is held, what cleared the gate, and every reason a comparison is not yet possible.

    Producible at zero records, and prints an honest zero rather than an empty
    table. For M0 this IS the deliverable: the expected result is that the corpus
    yields no matched pairs at all, and a report that says so clearly, with the
    reason, is worth more than a table of numbers that do not exist.
    """

    records_held: int
    records_cleared: int
    matched_ion_keys: int
    matched_ions_multi_platform: int
    platforms: tuple[str, ...]
    distinct_sources: int
    # Records whose charge carrier the source never named. Each holds a key unique
    # to itself and can never pair with anything, so they are counted apart from
    # matched_ion_keys rather than inflating it.
    #
    # ZERO HERE TODAY, AND DELIBERATELY KEPT ANYWAY. An unstated carrier is a
    # training blocker, so the gate refuses such a record before it can reach a
    # MatchedIonSet, and this figure is a second line of defence rather than the
    # only one. The place the count actually bites is the loader's held figures,
    # where such records ARE present: see MeasurementLoadReport.unmatchable_held.
    # If the blocker is ever relaxed, this reports the consequence instead of the
    # corpus quietly gaining keys that can never pair.
    unmatchable_records: int = 0
    by_platform: Mapping[str, int] = field(default_factory=dict)
    by_pair: Mapping[str, PairCount] = field(default_factory=dict)
    # Each record's OWN status, so this tally sums to records_cleared. A status
    # on an analyte INSIDE a record is not here; it is in synthetic_records,
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
            return (NOTHING_HELD,)
        if len(self.platforms) < MIN_PLATFORMS:
            problems.append(ONE_PLATFORM_ONLY.format(platforms=", ".join(self.platforms) or "none"))
            # Everything below is downstream of this and would only restate it.
            return tuple(problems)
        if self.matched_ions_multi_platform == 0:
            problems.append(NO_MATCHED_IONS.format(keys=self.matched_ion_keys, platforms=len(self.platforms)))
            return tuple(problems)
        if not any(count.scorable for count in self.by_pair.values()):
            problems.append(NO_EVALUABLE_PAIR.format(needed=MIN_MATCHED_IONS_PER_PAIR, folds=N_FOLDS))
        return tuple(problems)

    @property
    def warnings(self) -> tuple[str, ...]:
        """Reasons the result will be weak. A weak result that says so is still useful."""
        notes: list[str] = []
        if self.synthetic_records:
            notes.append(
                f"{self.synthetic_records} of {self.records_cleared} records are synthetic fixtures, declared"
                " as such by the record or by an analyte inside it: this set is a test, and no number from it"
                " describes real data"
            )
        if 0 < self.records_cleared < MIN_RECORDS_HELD:
            notes.append(
                f"the corpus holds {self.records_cleared} cleared records; {MIN_RECORDS_HELD} is the point at"
                " which a figure stops being provisional"
            )
        if self.matched_ions_multi_platform and self.matched_ions_multi_platform < TARGET_MATCHED_IONS:
            notes.append(
                f"{self.matched_ions_multi_platform} cross-platform matched ions is below the"
                f" {TARGET_MATCHED_IONS} at which the number is worth quoting, so expect wide intervals"
            )
        if self.records_cleared and self.distinct_sources < MIN_SOURCES_FOR_CROSS_STUDY:
            notes.append(
                f"{self.distinct_sources} source(s): under {MIN_SOURCES_FOR_CROSS_STUDY} no cross-study claim"
                " can be made, so every number is a within-study number"
            )
        if self.unmatchable_records:
            notes.append(
                f"{self.unmatchable_records} of {self.records_cleared} records do not name the charge carrier"
                " of their ion, so none of them can be matched with anything, including each other: two"
                " sources both reporting a 24+ ion may mean protons and ammonium, which are not the same ion."
                " Resolving them means reading each source's methods, never assuming protons"
            )
        for name, count in sorted(self.by_pair.items()):
            if count.matched_ions and not count.scorable:
                notes.append(
                    f"platform pair {name} holds {count.matched_ions} matched ion(s), under the"
                    f" {MIN_MATCHED_IONS_PER_PAIR} needed for grouped cross-validation; any figure for this"
                    " pair is an in-sample figure"
                )
        return tuple(notes)

    @property
    def shortfalls(self) -> tuple[str, ...]:
        """Everything worth saying, blockers first. Kept for callers that want one list."""
        return self.blockers + self.warnings

    @property
    def ready(self) -> bool:
        """Whether a comparison may proceed. Weakness does not stop it; inevaluability does."""
        return not self.blockers

    @property
    def refusal(self) -> str | None:
        return "; ".join(self.blockers) or None

    @property
    def maturity(self) -> MaturityStamp:
        # Always provisional here. The validated flag belongs to the milestone
        # that checks a fitted model against held-out matched ions.
        return MaturityStamp(
            data_maturity=DataMaturity.PROVISIONAL, matched_ion_count=self.matched_ions_multi_platform
        )

    def summary(self) -> str:
        lines = [
            "Cross-platform CCS readiness"
            f"  [{self.maturity.data_maturity}, {self.maturity.matched_ion_count} matched ions behind any number]",
            f"  records held               {self.records_held}",
        ]
        if self.synthetic_records:
            lines.append(
                f"  SYNTHETIC FIXTURES         {self.synthetic_records}  (the record, or an analyte in it,"
                " declared in code by a test; nothing here describes real data)"
            )
        lines += [
            f"  cleared the licence gate  {self.records_cleared}",
            f"  distinct matched ions     {self.matched_ion_keys}",
            *(
                [
                    f"  UNMATCHABLE               {self.unmatchable_records}  (charge carrier not stated by"
                    " the source; each can never pair with anything, including another such record)"
                ]
                if self.unmatchable_records
                else []
            ),
            f"  measured on >1 platform   {self.matched_ions_multi_platform}"
            f"  (of {TARGET_MATCHED_IONS} before the number is worth quoting)",
            f"  platforms represented     {len(self.platforms)}  (of {MIN_PLATFORMS} needed for any pair"
            " to be possible at all)",
            f"  distinct sources          {self.distinct_sources}"
            f"  (of {MIN_SOURCES_FOR_CROSS_STUDY} needed for any cross-study claim)",
            f"  conformal calibration     needs at least {MIN_CALIBRATION_RECORDS} held-out matched ions for"
            f" {1 - CONFORMAL_ALPHA:.0%} coverage to be attainable at all",
        ]
        if self.by_platform:
            lines.append("  by platform:")
            for name, count in sorted(self.by_platform.items()):
                lines.append(f"    {name}: {count} records")
        if self.by_reuse_status:
            lines.append("  by reuse status:")
            for status, count in sorted(self.by_reuse_status.items()):
                lines.append(f"    {status}: {count}")
        if self.by_pair:
            lines.append("  by platform pair:")
            for name, count in sorted(self.by_pair.items()):
                mark = "scorable" if count.scorable else "too small to score"
                lines.append(f"    {name}: {count.matched_ions} matched ions ({mark})")
        else:
            lines.append("  by platform pair:  none - no ion is held on two platforms")
        lines.append("  ready to compare: " + ("yes" if self.ready else "no"))
        for blocker in self.blockers:
            lines.append(f"    cannot be evaluated: {blocker}")
        for warning in self.warnings:
            lines.append(f"    will be weak: {warning}")
        return "\n".join(lines)


def _provenance_atom(record: object) -> str:
    """What counts as one study: the DOI where there is one, else the source text.

    Studies, not source strings. Twenty-four rows transcribed from one table are
    one source however their source text varies row to row.
    """
    doi = getattr(record, "doi", None)
    if doi:
        return f"doi:{str(doi).strip().casefold()}"
    return f"source:{str(getattr(record, 'source', '') or '').strip()}"


def _platform_of(record: object) -> str:
    ims = getattr(record, "ims_type", None)
    method = getattr(record, "dtims_method", None)
    return str(ims) if method is None else f"{ims}/{method}"


@dataclass(frozen=True)
class MatchedIonSet:
    """Records that have been through the licence gate, and the matched ions they form.

    The gate runs in __post_init__, so no instance can hold an ungated record
    however it was constructed. There is no classmethod to route around.

    Grouping is by the matched-ion key exactly as identity.py defines it, so the
    readiness figure and any later pairing agree on what one ion is. Where two
    records name one molecule through different identifiers they do NOT group
    here; merging those is matched-ion construction's job in M2, and doing it
    here would mean choosing silently between identifiers.
    """

    records: tuple[object, ...] = ()

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
            # UnbackedClaimError must not be told "licence not cleared" for a
            # record whose licence is fine and whose CLAIM is what nothing
            # backs, because that points at the wrong remedy.
            error = UnbackedClaimError if unbacked_only else LicenceGateError if licence_fault else TrainingGateError
            raise error(f"a matched-ion set cannot hold {len(faults)} record(s) that may not train: {problem}")

    def __len__(self) -> int:
        return len(self.records)

    def by_key(self) -> dict[object, set[str]]:
        """Matched-ion key -> the set of platforms it was measured on.

        Unmatchable keys are included here, because they are held and a caller
        asking what the corpus contains should see them. They are counted apart
        in the readiness report, and no pair can ever form from one: each carries
        its own record's provenance, so it is equal to nothing but itself.
        """
        found: dict[object, set[str]] = {}
        for record in self.records:
            key = getattr(record, "matched_ion_key", None)
            if key is None:
                continue
            found.setdefault(key, set()).add(_platform_of(record))
        return found

    @property
    def readiness(self) -> Readiness:
        keys = self.by_key()
        platforms = Counter(_platform_of(record) for record in self.records)
        pairs: Counter[str] = Counter()
        for measured_on in keys.values():
            for first in sorted(measured_on):
                for second in sorted(measured_on):
                    if first < second:
                        pairs[f"{first} vs {second}"] += 1
        statuses: dict[str, int] = {}
        for record in self.records:
            status = str(getattr(record, "reuse_status", "unknown"))
            statuses[status] = statuses.get(status, 0) + 1
        sources = {_provenance_atom(record) for record in self.records}
        matchable = {key: seen for key, seen in keys.items() if getattr(key, "matchable", True)}
        return Readiness(
            records_held=len(self.records),
            # Every record here has already cleared the gate, by construction.
            records_cleared=len(self.records),
            matched_ion_keys=len(matchable),
            unmatchable_records=len(keys) - len(matchable),
            matched_ions_multi_platform=sum(1 for measured_on in matchable.values() if len(measured_on) > 1),
            platforms=tuple(sorted(platforms)),
            distinct_sources=len(sources),
            by_platform=dict(platforms),
            by_pair={name: PairCount(matched_ions=count) for name, count in pairs.items()},
            by_reuse_status=statuses,
            # Asked of the whole record, not of its outermost status: an analyte
            # invented by a test inside a measurement claiming a real status is
            # still a fixture, and a set of them is still a test.
            synthetic_records=sum(1 for record in self.records if declares_synthetic(record)),
        )


def require_ready(matched_ion_set: MatchedIonSet) -> Readiness:
    """The readiness of a gated set, or InsufficientDataError naming why it cannot be evaluated.

    Every check is an explicit raise rather than an assert, so it survives
    python -O. A bare list is the natural accidental bypass, and a TypeError
    would be swallowed by a caller skipping malformed rows, so it is refused as
    a gate error instead.
    """
    if not isinstance(matched_ion_set, MatchedIonSet):
        raise TrainingGateError(NOT_A_MATCHED_ION_SET.format(given=type(matched_ion_set).__name__))
    readiness = matched_ion_set.readiness
    if readiness.refusal is not None:
        raise InsufficientDataError(readiness.refusal)
    return readiness
