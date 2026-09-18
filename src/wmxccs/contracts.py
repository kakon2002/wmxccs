"""Request and response bodies for the HTTP API.

THE RESPONSE CONTRACT, AND THE ONE PROMISE IT MAKES
---------------------------------------------------
Every harmonization response carries FOUR things, and the first of them is the
reason the other three can be trusted:

1. THE ORIGINALS, UNTOUCHED. Every measurement submitted comes back exactly as it
   was sent. Not rounded, not normalised, not re-expressed. A harmonized value is
   an ADDITIONAL number, never a replacement, and a caller who disagrees with the
   correction must be able to reach past it to the measurement it was computed
   from without asking anybody for it.
2. THE HARMONIZED ESTIMATE, with its interval and the kind of interval it is.
3. THE CONFIDENCE GRADE, with every reason it is not higher.
4. THE PROVENANCE OF EVERY INPUT: where it came from, on what terms, and who
   recorded those terms.

NO PLACEHOLDER NUMBERS, ANYWHERE
--------------------------------
While no model exists, `/harmonize` answers 501 and the body carries no numbers at
all - no zero, no null standing in for a value, no interval of infinite width. The
harmonized estimate is absent rather than empty. That is the difference between
"we cannot do this yet" and "here is a number we made up", and every field below
is shaped so that the second is not expressible.

WHY THE ESTIMATE CARRIES TWO CORRECTIONS
----------------------------------------
M3 found that three outlying ions in twenty-four moved a fitted Deming slope from
the injected 1.020 to 1.031 - so a slope-derived correction applied to every
well-behaved ion was partly the work of the three that were not. The contract
therefore has room for BOTH the slope-derived and the median-derived correction,
and says which one the headline value came from. A caller comparing the two learns
something a single number would have hidden.
"""

from __future__ import annotations

from datetime import date
from enum import StrEnum
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from .grading import ConfidenceGrade
from .models import CCSMeasurement
from .readiness import DataMaturity, MaturityStamp
from .reuse import ReuseStatus


class HealthResponse(BaseModel):
    status: Literal["ok"] = "ok"
    version: str
    model_loaded: bool = Field(
        description="Whether a harmonization model is loaded. False until one has been fitted on real"
        " cross-platform matched ions, which has not happened."
    )
    matched_ions_available: int = Field(
        ge=0,
        description="Cross-platform matched ions the loaded model rests on. Zero while none exist, which"
        " is why no model does.",
    )


class SourceProvenance(BaseModel):
    """Where one input came from and on what terms, carried back with the answer.

    Not a convenience. A harmonized CCS inherits the licence of every measurement
    behind it, and a caller who cannot see those terms cannot know what they are
    allowed to do with the number they have been given.
    """

    model_config = ConfigDict(frozen=True)

    source: str = Field(description="The source as the record states it. Free text, and not a citation.")
    doi: str | None = Field(default=None, description="Bare DOI, where the record has one.")
    reuse_status: ReuseStatus = Field(
        description="What the source's terms allow. 'unverified' is the default and means nobody has read them."
    )
    licence: str | None = Field(
        default=None, description="The licence as the registry records it, where this source is registered."
    )
    licence_reported_by: str | None = Field(
        default=None, description="Who reported reading those terms. Not a claim that this code verified anything."
    )
    licence_reported_on: date | None = Field(default=None, description="When that report was received.")
    registered: bool = Field(
        description="Whether the source is in the licence registry at all. False means the claim on the record"
        " is backed by nothing."
    )


class CorrectionBasis(StrEnum):
    """Which correction the headline harmonized value came from.

    Both are reported when both exist. They differ exactly when the fit is levered
    by ions that do not transfer, which is the case a caller most needs to see.
    """

    # From the Deming slope and intercept.
    SLOPE = "slope"
    # From the stratum's median signed difference. Unmoved by a few bad ions,
    # and unable to express a correction that varies with size.
    MEDIAN = "median"
    # A slope fitted so that a few ions cannot dominate it.
    ROBUST_SLOPE = "robust_slope"


class IntervalKind(StrEnum):
    """What the interval IS. It travels with the numbers, like an uncertainty type.

    An interval without its kind is read as whatever the reader assumes, which is
    the same failure uncertainty_type exists to prevent one level down.
    """

    CONFORMAL_PREDICTION = "conformal_prediction"
    LIMITS_OF_AGREEMENT = "limits_of_agreement"


class HarmonizedEstimate(BaseModel):
    """The corrected value. ABSENT, never empty, while no model exists."""

    model_config = ConfigDict(frozen=True)

    ccs: float = Field(gt=0, description="The harmonized cross section, square angstrom.")
    basis: CorrectionBasis = Field(description="Which correction the value above came from.")
    slope_derived_ccs: float | None = Field(
        default=None,
        gt=0,
        description="What the Deming slope and intercept give. Reported alongside the median-derived value"
        " because the two diverge exactly when a few ions are levering the fit.",
    )
    median_derived_ccs: float | None = Field(
        default=None,
        gt=0,
        description="What the stratum's median signed difference gives. Unmoved by outliers, and unable to"
        " express a correction that varies with size.",
    )
    interval_low: float = Field(gt=0)
    interval_high: float = Field(gt=0)
    interval_coverage: float = Field(
        gt=0, lt=1, description="Nominal coverage, e.g. 0.90. Never 1, which no finite interval attains."
    )
    interval_kind: IntervalKind
    reference_platform: str = Field(description="The platform the correction refers the value to.")
    matched_ions_behind_it: int = Field(
        ge=0, description="How many matched ions the correction was fitted on. Travels with the number."
    )


class ConfidenceReason(BaseModel):
    model_config = ConfigDict(frozen=True)

    rule: str
    falls_to: ConfidenceGrade
    detail: str


class ConfidenceReport(BaseModel):
    """The grade, and every reason it is not higher."""

    model_config = ConfigDict(frozen=True)

    grade: ConfidenceGrade
    reasons: tuple[ConfidenceReason, ...] = ()
    not_checked: tuple[str, ...] = Field(
        default=(),
        description="Checks that could not be run, with why. Kept apart from the reasons: not checked is a"
        " different statement from checked and fine.",
    )


class HarmonizeRequest(BaseModel):
    """One measurement or a set of them.

    A set rather than only a single record because the unit of comparison is the
    matched ion, and a caller holding two platforms' values for one ion should be
    able to submit both and be told what they say about each other.
    """

    model_config = ConfigDict(extra="forbid")

    measurements: tuple[CCSMeasurement, ...] = Field(
        min_length=1, description="The measurements to harmonize. Returned untouched in the response."
    )


class HarmonizedMeasurement(BaseModel):
    """One input, its provenance, and what could be said about it."""

    model_config = ConfigDict(frozen=True)

    original: CCSMeasurement = Field(
        description="EXACTLY as submitted. Never rounded, normalised or re-expressed, and never replaced by"
        " the harmonized value."
    )
    provenance: SourceProvenance
    harmonized: HarmonizedEstimate | None = Field(
        default=None,
        description="Absent, not empty, where no correction can be made. There is no placeholder value.",
    )
    confidence: ConfidenceReport | None = Field(
        default=None, description="Absent where there is no harmonized value to grade."
    )


class HarmonizeResponse(BaseModel):
    """The 200 body. Not reachable until a model exists; specified so callers can build against it."""

    model_config = ConfigDict(frozen=True)

    measurements: tuple[HarmonizedMeasurement, ...]
    maturity: MaturityStamp


class HarmonizationUnavailable(BaseModel):
    """The 501 body.

    It carries the provenance of every input and the maturity stamp, and NO
    NUMBERS. A caller gets back exactly what they sent, plus the terms it came on,
    plus a statement that no model exists. They do not get a value, an interval or
    a grade, because there is nothing behind any of the three.
    """

    model_config = ConfigDict(frozen=True)

    error: Literal["not_implemented"] = "not_implemented"
    detail: str
    measurements: tuple[HarmonizedMeasurement, ...] = Field(
        description="The originals and their provenance, returned unchanged. Every harmonized estimate and"
        " confidence report inside is absent."
    )
    maturity: MaturityStamp


class ConfidenceRule(BaseModel):
    """One grading rule, published so the scheme can be argued with."""

    model_config = ConfigDict(frozen=True)

    rule: str
    falls_to: tuple[str, ...]
    why: str
    threshold: str
    basis: str = Field(description="Whether the threshold is derived from something or is policy.")


class ConfidenceRulesResponse(BaseModel):
    model_config = ConfigDict(frozen=True)

    grades: tuple[str, ...] = Field(description="Every grade, worst last.")
    rules: tuple[ConfidenceRule, ...]
    note: str


__all__ = [
    "ConfidenceReason",
    "ConfidenceReport",
    "ConfidenceRule",
    "ConfidenceRulesResponse",
    "CorrectionBasis",
    "DataMaturity",
    "HarmonizationUnavailable",
    "HarmonizeRequest",
    "HarmonizeResponse",
    "HarmonizedEstimate",
    "HarmonizedMeasurement",
    "HealthResponse",
    "IntervalKind",
    "MaturityStamp",
    "SourceProvenance",
]
