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

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .grading import ConfidenceGrade
from .models import CCSMeasurement
from .scope import ComparisonScope
from .readiness import DataMaturity, MaturityStamp
from .reuse import ReuseStatus


class ModelVersion(BaseModel):
    """Which model produced an answer, so the answer can be reproduced or told apart.

    TWO DIGESTS rather than one, and the pair is the useful part. A different corpus with
    the same parameters means the data moved without moving the fit; the same corpus with
    different parameters means the code did. A single combined hash would say only that
    something had.

    The model is refitted from the seed files at every startup, which is fine. What is not
    fine is a changed seed file changing the answers with nothing recording it - two
    deployments could give different numbers for one input and neither response would say
    so. This is what makes an answer repeatable.
    """

    model_config = ConfigDict(frozen=True)

    corpus_sha256: str = Field(
        min_length=64,
        max_length=64,
        description="Over the records BEHIND THE FIT - not the files they came from, and not every record"
        " loaded. A measurement that feeds no correction is not in it, so data that cannot change an"
        " answer does not change this digest. It moves exactly when an answer could.",
    )
    parameters_sha256: str = Field(
        min_length=64,
        max_length=64,
        description="Over everything that decides an answer: which stratum applies, on what basis, with"
        " which slopes and offsets, and the leave-one-out residuals that set every interval.",
    )

    @property
    def short(self) -> str:
        return f"{self.corpus_sha256[:12]}/{self.parameters_sha256[:12]}"

    @classmethod
    def of(cls, fingerprint) -> "ModelVersion":
        return cls(
            corpus_sha256=fingerprint.corpus, parameters_sha256=fingerprint.parameters
        )


class HealthResponse(BaseModel):
    status: Literal["ok"] = "ok"
    version: str
    model_loaded: bool = Field(
        description="Whether a harmonization model is loaded. False until one has been fitted on real"
        " cross-platform matched ions, which has not happened."
    )
    matched_ions_available: int = Field(
        ge=0,
        description="Cross-platform matched ions the loaded model rests on. Zero where none is loaded.",
    )
    model_version: ModelVersion | None = Field(
        default=None,
        description="Absent only where no model is loaded. Compare it against the version on an estimate"
        " to know whether two answers came from the same model.",
    )


class SourceProvenance(BaseModel):
    """Where one input came from and on what terms, carried back with the answer.

    Not a convenience. A harmonized CCS inherits the licence of every measurement
    behind it, and a caller who cannot see those terms cannot know what they are
    allowed to do with the number they have been given.

    WHAT A SUBMITTED RECORD'S LICENCE STATUS DOES NOT DO, decided deliberately on
    20 September 2026 and stated here rather than left to be inferred:

    It does NOT gate the response. A measurement claiming `unverified`, `excluded` or
    `non_commercial_no_derivatives` is still corrected, and the reasons are these:

    - The licence gate governs what may enter a FIT. The model is already fitted, on
      records that passed the gate. A submitted measurement enters no fit and changes
      no parameter.
    - The status is SELF-ASSERTED AND UNVERIFIABLE. A caller refused for `excluded`
      edits the field to `open_attribution` and resubmits. A gate on a field the caller
      supplies protects nothing and obstructs only honest callers.
    - The caller's relationship to their own data is not knowable here. They may be its
      author.

    RESPONSIBILITY FOR THE INPUT'S TERMS THEREFORE REMAINS THE CALLER'S, and that is
    the disclosure rather than an assumption. It matters most for
    `non_commercial_no_derivatives`: a harmonized value IS a derivative of the input, so
    a caller holding data under a no-derivatives licence is the party making a
    derivative of it. This service does not and cannot check that.

    `synthetic_fixture` is the exception and is REFUSED, because it is not a licence
    claim at all - it is a declaration that the record was built in code and is not a
    measurement. Refusing it protects the integrity of the output, which is this
    service's responsibility, rather than a licence, which is not.
    """

    model_config = ConfigDict(frozen=True)

    source: str = Field(description="The source as the record states it. Free text, and not a citation.")
    doi: str | None = Field(default=None, description="Bare DOI, where the record has one.")
    reuse_status_claimed: ReuseStatus = Field(
        description="What the SUBMITTED RECORD says its source's terms allow. A claim by whoever sent it,"
        " self-asserted and unverifiable by this service. Named 'claimed' because it used to be named"
        " 'reuse_status', which read as established fact - so a response could show a caller's claim of"
        " 'excluded' beside the registry's open-access licence text and look like it had verified both."
    )
    reuse_status_in_registry: ReuseStatus | None = Field(
        default=None,
        description="What the LICENCE REGISTRY holds for this DOI, where the DOI is registered. None means"
        " nobody has recorded reading the terms. This is the half that is actually backed by somebody.",
    )
    claim_backed_by_registry: bool | None = Field(
        default=None,
        description="Whether the claim and the registry agree. None where the source is not registered, so"
        " there is nothing to agree with. FALSE IS THE INTERESTING CASE: the record asserts terms nobody"
        " recorded, and this field is what makes that visible instead of leaving two fields to be compared"
        " by a reader who may not notice.",
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
    # Distinct from the above, and the distinction is the guarantee rather than the
    # family. Split conformal proves 1-alpha; jackknife+ proves 1-2*alpha, so a
    # nominally 90 per cent jackknife+ interval is guaranteed at 80. Labelling one as
    # the other would be read as the stronger claim, which is the failure this enum
    # exists to prevent one level up.
    JACKKNIFE_PLUS = "jackknife_plus"
    LIMITS_OF_AGREEMENT = "limits_of_agreement"


class ScopeReport(BaseModel):
    """What the correction behind a number is entitled to be quoted as.

    REQUIRED on every harmonized estimate, and there is no default. A number whose scope
    could be omitted is a number that will be quoted without it.

    `scope` is not settable independently: it is filled from `scope.ScopeStamp`, where it is
    a property of the studies rather than a field, and a validator here re-derives it from
    `studies` so that a hand-built report cannot widen the claim without naming a second
    study. There is deliberately no member for interlaboratory reproducibility anywhere in
    this package - see scope.ComparisonScope.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    scope: ComparisonScope
    studies: tuple[str, ...] = Field(
        min_length=1, description="The studies behind the fit, as DOIs where known. At least one."
    )
    platforms: tuple[str, ...] = Field(min_length=1)
    records_behind_it: int = Field(ge=1)
    caveat: str = Field(
        min_length=1,
        description="The scope in words, to be printed WITH the number and never instead of it.",
    )

    @model_validator(mode="after")
    def scope_matches_the_studies(self) -> "ScopeReport":
        expected = (
            ComparisonScope.CROSS_STUDY if len(set(self.studies)) > 1 else ComparisonScope.WITHIN_STUDY
        )
        if self.scope is not expected:
            raise ValueError(
                f"scope {self.scope.value!r} does not follow from {len(set(self.studies))} study(ies):"
                f" it must be {expected.value!r}. The scope is DERIVED from the provenance, so widening it"
                f" means naming another study, which is data somebody has to produce"
            )
        return self

    @classmethod
    def of(cls, stamp) -> "ScopeReport":
        """Build from a scope.ScopeStamp. The only intended construction path."""
        return cls(
            scope=ComparisonScope(stamp.scope.value),
            studies=tuple(stamp.studies),
            platforms=tuple(stamp.platforms),
            records_behind_it=stamp.records_behind_it,
            caveat=stamp.caveat(),
        )


class HarmonizedEstimate(BaseModel):
    """The corrected value. ABSENT, never empty, while no model exists.

    THE FIELD ORDER IS PART OF THE CONTRACT, and it is ordered by what a reader must not
    be allowed to miss rather than by what the code computes first.

    A response is read from the top, and truncated from the bottom - by a log line, a
    console, a table of the "main" columns, a person skimming. So the first four fields are
    the number, the two bounds it sits between, and the scope caveat saying what it may be
    quoted as; then the two coverage figures, adjacent, with the GUARANTEED one first
    because it is the weaker and therefore the honest one. A reader who takes only the
    opening of this object still has a number, its uncertainty, and the sentence that stops
    them calling it interlaboratory reproducibility.

    It used to open with the number, the basis and two alternative corrections, and put the
    scope caveat ten fields down and the guaranteed coverage thirteen. Anyone who stopped
    reading early got three cross sections and no idea what any of them meant.
    """

    model_config = ConfigDict(frozen=True)

    # --- what the number is, and what it is worth ----------------------------------------
    ccs: float = Field(gt=0, description="The harmonized cross section, square angstrom.")
    interval_low: float = Field(gt=0)
    interval_high: float = Field(gt=0)
    scope: ScopeReport = Field(
        description="REQUIRED, no default. What this correction may be quoted as. A harmonized cross"
        " section cannot be serialised without it, and it is placed HERE, immediately under the number"
        " and its interval, because a caveat that arrives ten fields later has already been skipped."
    )
    guaranteed_coverage: float = Field(
        gt=0,
        lt=1,
        description="WHAT THE INTERVAL PROVES, and the figure to rely on. Jackknife+ guarantees"
        " 1-2*alpha, so a nominally 90% interval is guaranteed at 80%: about one value in five may fall"
        " outside its stated range rather than one in ten. First of the two coverage figures, and"
        " directly beside the interval it describes, because it is the weaker of the two.",
    )
    interval_coverage: float = Field(
        gt=0,
        lt=1,
        description="Nominal coverage, e.g. 0.90 - the level the quantile is taken at, not a proved"
        " guarantee. Never 1, which no finite interval attains. Adjacent to `guaranteed_coverage` so"
        " that the two can never be read apart; quoting this one alone overstates the interval.",
    )
    interval_kind: IntervalKind
    interval_is_informative: bool = Field(
        description="False where the interval is the full observed range - honest, and excluding nothing."
        " See readiness.smallest_informative_calibration_set."
    )
    interval_accounts_for_submitted_uncertainty: bool = Field(
        default=False,
        description="FALSE, always, and stated because the response invites the opposite reading: the"
        " submitted measurement's own ccs_uncertainty is echoed on `original`, and the two are unrelated."
        " The interval is MODEL-DERIVED - it comes from the spread of leave-one-out refits of the stratum"
        " - and is identical whether the submitted uncertainty is 0.0001, 50, or absent entirely."
        " Propagating the caller's uncertainty into it is not implemented; this field exists so that is a"
        " stated fact rather than something a reader has to discover.",
    )

    # --- how it was arrived at -------------------------------------------------------------
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
    reference_platform: str = Field(description="The platform the correction refers the value to.")
    matched_ions_behind_it: int = Field(
        ge=0, description="How many matched ions the correction was fitted on. Travels with the number."
    )
    model_version: ModelVersion = Field(
        description="REQUIRED, no default. Which model produced this number. Without it the same input"
        " could give two different answers on two days with nothing saying which was which."
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


MAX_MEASUREMENTS_PER_REQUEST = 1_000
"""A ceiling on one request, to go with the floor that was always there.

Not a rate limit and not a security control - there is no authentication on this service,
and LIMITATIONS 7G records both as deliberate omissions with what each would take. This is
the narrower thing: a single request cannot ask for unbounded work by accident. 1000 is
roughly twice the whole seed corpus, so no honest use of this endpoint meets it.
"""


class HarmonizeRequest(BaseModel):
    """One measurement or a set of them.

    A set rather than only a single record because the unit of comparison is the
    matched ion, and a caller holding two platforms' values for one ion should be
    able to submit both and be told what they say about each other.
    """

    model_config = ConfigDict(extra="forbid")

    measurements: tuple[CCSMeasurement, ...] = Field(
        min_length=1,
        max_length=MAX_MEASUREMENTS_PER_REQUEST,
        description="The measurements to harmonize. Returned untouched in the response. Bounded at"
        f" {MAX_MEASUREMENTS_PER_REQUEST}: the floor was here from the start and the ceiling was not, so"
        " a request could carry an unbounded number of bounded records - the same problem one level up"
        " from the text fields. Split a larger set across requests; nothing here is stateful and the"
        " model is identical between them, as the model_version on each estimate shows.",
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
        default=None,
        description="The grade. Present even where `harmonized` is absent IF a grade was computed and the"
        " number withheld because of it - a caller is owed the reason the value is missing.",
    )
    not_harmonized_because: str | None = Field(
        default=None,
        description="REQUIRED whenever `harmonized` is absent, and absent whenever it is present. Which of"
        " those two it is decides what a caller should do next: a platform the model does not cover is a"
        " different problem from an ion outside the range its correction was fitted over.",
    )

    @model_validator(mode="after")
    def a_missing_value_carries_its_reason(self) -> "HarmonizedMeasurement":
        if self.harmonized is None and not (self.not_harmonized_because or "").strip():
            raise ValueError(
                "a measurement with no harmonized value must say why. An absent value with no reason is"
                " indistinguishable from an oversight, and a caller cannot tell whether to fix their input,"
                " wait for more data, or stop asking"
            )
        if self.harmonized is not None and self.not_harmonized_because:
            raise ValueError(
                "a measurement cannot both carry a harmonized value and a reason there is none:"
                f" {self.not_harmonized_because!r}"
            )
        return self


class HarmonizeResponse(BaseModel):
    """The 200 body.

    Carries a validator rather than only a docstring: a validated maturity beside a
    within-study scope is incoherent, and an incoherent response is the artefact somebody
    would quote.
    """

    model_config = ConfigDict(frozen=True)

    measurements: tuple[HarmonizedMeasurement, ...]
    maturity: MaturityStamp

    @model_validator(mode="after")
    def validated_maturity_needs_more_than_one_study(self) -> "HarmonizeResponse":
        if self.maturity.data_maturity is not DataMaturity.VALIDATED:
            return self
        within = [
            estimate.scope.studies
            for measurement in self.measurements
            if (estimate := measurement.harmonized) is not None
            and estimate.scope.scope is ComparisonScope.WITHIN_STUDY
        ]
        if within:
            raise ValueError(
                "a response may not report maturity 'validated' while a correction in it is within-study."
                " Validation means checked against data the model was not fitted on, and a single-study"
                f" corpus has none by definition. {len(within)} estimate(s) here are within-study"
            )
        return self


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
    """One grading rule, published so the scheme can be argued with.

    `extra="forbid"` IS LOAD-BEARING. The endpoint builds this straight from the dict
    grade_rules() returns, so a key added to the scheme and not added here fails loudly at
    the first request instead of being dropped on the floor. It was dropped on the floor:
    `applies_to` existed in the scheme for the length of one afternoon on 20 September 2026
    while the endpoint hand-copied five fields by name, so the scope limit it carries was
    published to nobody. That is the same failure as the one it was written to describe - a
    thing that looks complete because nothing checks the join.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    rule: str
    falls_to: tuple[str, ...]
    why: str
    threshold: str
    basis: str = Field(description="Whether the threshold is derived from something or is policy.")
    applies_to: str = Field(
        description="WHICH IONS THIS RULE CAN ANSWER FOR. Most rules are properties of the"
        " calibration group and answer for any ion. One is not: the outlier rule is a lookup"
        " against ions already in this corpus, so it cannot answer for a new one. Where a rule"
        " cannot run it is named in the response's not_checked and the grade is demoted a notch;"
        " it is never counted as passed."
    )


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
