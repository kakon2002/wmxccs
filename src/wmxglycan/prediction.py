"""The served shapes: what a prediction is, what may be attached to it, and what a comparison says.

WHAT A PREDICTION IS HERE, AND WHAT IT IS NOT

It is a frozen, ranked set of candidate structures for one composition and one ion, with the
evidence behind each and the reasons the platform gives for asking for experimental validation.
IT CONTAINS NO PREDICTED CROSS SECTION, because there is no glycan CCS model and V1 will not
have one. Anything that looks like a CCS in a response is a MEASURED value - either one the
platform holds as a reference, or one the caller attached - and the field names say which.

THE NUMBER IS NEVER CALLED A PROBABILITY. Ruled 27 September 2026: the spec asks for a
probability-LIKE confidence, so the served field is `evidence_share` and no field in this module
or in `api.py` is named `probability`, `p_correct` or `likelihood`. A test walks every response
model and fails on any such name, because the name is what a consumer reads and a caveat two
files away does not travel with it. The declared POLICY prior and all four standard alternatives
ride in the same object as the number.

EVERY RESPONSE CARRIES THE STAMP. `ModelStamp` holds the package version, the two pipeline
digests and the data snapshot the answer was made against. It is on every response shape here,
including the errors, so a prediction read back in a year can be checked against the pipeline
that produced it rather than against whatever is installed then.

THE DECISION FIELD IS CONSTANT, AND THAT IS THE ANSWER RATHER THAN A DEFECT. Every prediction
lands on IM_VALIDATION_REQUIRED, for two reasons that hold for every possible input. The rules
are published per response with what each applies to, so the constancy is readable from the
response and not only from the documentation. See docs/GLYCAN_LIMITATIONS.md.
"""

from __future__ import annotations

from enum import StrEnum
from typing import Mapping, Self

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .ranking import Calibration, CoverageState, Decision, RankedSet, SetCompleteness

# POLICY, CHOSEN RATHER THAN MEASURED. How far an attached measurement may sit from a reference
# before the comparison calls it a disagreement, as a percentage of the reference. It is the same
# kind of choice as the CCS core's limits-of-agreement threshold and is declared the same way:
# named here, published in the response, and never passed in as an argument, because a limit that
# can be passed in is a limit that can be passed a number chosen to make a case pass.
AGREEMENT_LIMIT_PERCENT = 2.0

NO_PREDICTED_VALUE = (
    "there is no predicted cross section to compare against, and there will not be one in V1."
    " This platform ranks candidate STRUCTURES on evidence; it does not predict a CCS, because"
    " no glycan CCS model exists here. So a delta against a prediction is not merely unmeasured,"
    " it is undefined - and an interval coverage figure would be a statement about an interval"
    " that was never produced. What CAN be compared is reported beside this: an attached"
    " measurement against a reference the platform holds, and attached measurements against each"
    " other"
)
INTERVAL_UNEVALUABLE = (
    "no interval was predicted, so coverage cannot be evaluated. An interval requires a predicted"
    " value and an uncertainty on it; this platform produces neither"
)


class ModelStamp(BaseModel):
    """What a response was produced by. On every shape in this module, without exception."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    package: str = "wmxglycan"
    version: str
    # `<parameters>/<corpus>`, the same form the CCS core prints, so the two read alike.
    fingerprint: str
    fingerprint_parameters: str
    fingerprint_corpus: str
    data_snapshot: str
    dataset: str
    dataset_licence: str
    dataset_attribution: str
    structures_indexed: int
    distinct_structures: int
    compositions_in_corpus: int

    @classmethod
    def of(cls, fingerprint) -> ModelStamp:
        snapshot = fingerprint.snapshot
        return cls(
            version=fingerprint.version,
            fingerprint=fingerprint.short,
            fingerprint_parameters=fingerprint.parameters,
            fingerprint_corpus=fingerprint.corpus,
            data_snapshot=snapshot.digest,
            dataset=snapshot.dataset,
            dataset_licence=snapshot.licence,
            dataset_attribution=snapshot.attribution,
            structures_indexed=snapshot.structures_indexed,
            distinct_structures=snapshot.distinct_structures,
            compositions_in_corpus=snapshot.compositions,
        )


class ReasonOut(BaseModel):
    """One piece of the reasoning behind a candidate, with whatever cites it."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    kind: str = Field(description='"enzyme", "rule" or "ordering".')
    text: str
    enzymes: tuple[str, ...] = ()
    reference: str | None = None


class CandidateOut(BaseModel):
    """One enumerated structure, and the evidence for it."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    canonical_key: str
    iupac_condensed: str
    class_id: str = Field(description="The indistinguishable class this candidate belongs to.")
    attested: bool
    attesting_accessions: tuple[str, ...] = Field(
        default=(),
        description="GlyTouCan accessions attesting this structure, where the records carry one."
        " An EMPTY tuple does not mean unattested: 8,809 of 12,664 reference records carry an"
        " accession, so `attested` is the authority.",
    )
    reference_rows: int = Field(
        default=0,
        description="Reference ROWS behind this structure. 2 or 3 means the source spelled one"
        " structure twice; reported so the discrepancy is visible, and never scored.",
    )
    rationale: tuple[ReasonOut, ...] = ()
    rules_applicable: int = 0
    ordering_caveats: int = 0


class ClassOut(BaseModel):
    """Candidates no feature the platform computes can separate, and their pooled evidence.

    THE UNIT OF RANKING. `members` is a sorted list FOR RENDERING ONLY; the order carries no
    meaning and `is_a_tie` says so explicitly, because a JSON array is ordered whatever the
    producer intends and a consumer will read the first element as the best one unless told.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    class_id: str
    members: tuple[str, ...]
    members_order_is_not_meaningful: bool = True
    is_a_tie: bool
    attested_structures: int
    attesting_accessions: tuple[str, ...] = ()
    reference_rows: int = 0
    # NOT a probability. See the module docstring and the `confidence` block on the response.
    evidence_share: float | None = None
    evidence_share_under_priors: Mapping[str, float | None] = Field(default_factory=dict)


class BandOut(BaseModel):
    """Classes holding equal evidence. Unordered within, ordered between."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    rank: int
    classes: tuple[str, ...]
    classes_order_is_not_meaningful: bool = True
    candidates: int
    attested_structures: int
    share_per_class: float | None = None
    band_total_share: float | None = None
    share_under_priors: Mapping[str, float | None] = Field(default_factory=dict)


class ConfidenceOut(BaseModel):
    """The number, what it means, and the lever it moves on.

    Named `confidence` in the response and NEVER `probability`. Everything that qualifies the
    number travels in the same object, so a consumer cannot read the figure without the prior it
    was computed under sitting beside it.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    prior: str = Field(description="Which of `priors` the served evidence_share uses.")
    prior_pseudocount: float
    priors: tuple[str, ...] = Field(
        description="Every prior the response reports a share under. Four standard"
        " non-informative Dirichlet concentrations; the choice is POLICY and it is the single"
        " largest lever in the output."
    )
    calibration: Calibration
    what_would_calibrate_it: str
    what_the_number_means: str

    # --- THE SHARES DO NOT SUM TO 1, AND THE ANSWER IS HERE RATHER THAN IN A DOCUMENT ------------
    # Somebody who adds up the per-class shares gets about 0.69 for Hex5HexNAc4Fuc1 and has to be
    # told why WHERE THEY ARE LOOKING. The three fields that DO account to 1 are named as such, the
    # sum is served so nobody has to compute it, and the one field that is a SUBSET rather than a
    # fourth term says so in its own description - adding that one too is the obvious next mistake.
    candidate_shares_sum: float | None = Field(
        default=None,
        description="Every share served over the candidates shown, added up. It is LESS THAN 1 by"
        " design and `mass_not_on_any_candidate` is the difference.",
    )
    mass_not_on_any_candidate: float | None = Field(
        default=None,
        description="1 - candidate_shares_sum. Evidence mass on hypotheses that are not among the"
        " candidates at all, which is why a candidate set can never be treated as exhaustive.",
    )
    mass_on_reference_structures_not_enumerated: float | None = Field(
        default=None,
        description="Part of `mass_not_on_any_candidate`: structures the reference corpus holds"
        " for this composition that the enumerator did not propose.",
    )
    mass_on_a_structure_nobody_proposed: float | None = Field(
        default=None,
        description="Part of `mass_not_on_any_candidate`: the catch-all hypothesis, a structure"
        " neither proposed here nor deposited anywhere. It is never zero, which is what stops the"
        " shares summing to 1.",
    )
    mass_on_unattested_classes: float | None = Field(
        default=None,
        description="A SUBSET OF `candidate_shares_sum`, NOT a fourth term: mass on candidates"
        " that ARE shown but which no reference structure attests. Do not add this to the three"
        " that account to 1 - it is already inside `candidate_shares_sum`.",
    )
    shares_account_to_one: bool | None = Field(
        default=None,
        description="Whether candidate_shares_sum + mass_on_reference_structures_not_enumerated +"
        " mass_on_a_structure_nobody_proposed == 1, to within floating point. Served as a checked"
        " identity rather than a promise; a consumer can verify the accounting from the response.",
    )
    why_the_shares_do_not_sum_to_one: str = Field(
        default="",
        description="In one paragraph, for a reader who has just added up the candidate shares.",
    )
    mass_under_priors: Mapping[str, Mapping[str, float | None]] = Field(default_factory=dict)


class CoverageOut(BaseModel):
    """Whether the candidate set is even capable of containing the answer."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    composition: str
    state: CoverageState
    completeness: SetCompleteness = Field(
        description="INCOMPLETE where the corpus proves something is missing, otherwise"
        " NOT_ESTABLISHABLE. There is no COMPLETE value: a corpus records what has been"
        " deposited, not what exists."
    )
    truth_may_not_be_in_the_candidate_set: bool
    reference_rows: int
    reference_structures: int
    attested_by_a_candidate: int
    not_enumerated: int
    not_enumerated_examples: tuple[str, ...] = ()
    summary: str


class DecisionRuleOut(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    fires: bool
    applies_to: str
    because: str
    decision_if_it_fires: Decision


class DecisionOut(BaseModel):
    """What the platform says should happen next, and every rule behind it.

    CURRENTLY CONSTANT at IM_VALIDATION_REQUIRED for every input, and that is the true state of
    a platform with no CCS model rather than a broken field. `constant_today` and
    `why_it_is_constant` are served so a consumer meets that fact in the response.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    decision: Decision
    rules: tuple[DecisionRuleOut, ...]
    constant_today: bool = True
    why_it_is_constant: str = (
        "Two of the rules below fire on every possible input, so IM_VALIDATION_REQUIRED is the"
        " only reachable value: no validated model exists, and the completeness of a candidate"
        " set cannot be established from a corpus that records depositions rather than existence."
        " This is not a defect to be fixed. What would move it is a model validated against"
        " independently known structures, and a cross section that resolves to one structure"
        " rather than to a composition and an ion."
    )


class RuleAccountingOut(BaseModel):
    """The three different numbers people mean by "rules satisfied", named apart."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    rules_in_scheme: int
    rules_violated: int
    verified_by_running_the_check: bool
    rules_applicable_range: tuple[int, ...] = Field(
        default=(),
        description="The distinct counts of rules that actually BEAR on a candidate. It varies,"
        " and it is deliberately not a score: it is largely a count of how branched a candidate"
        " is, and no curated rule says a more branched structure is more likely.",
    )
    summary: str


class PredictionRequest(BaseModel):
    """Create a prediction from a composition and its ion metadata."""

    model_config = ConfigDict(extra="forbid")

    composition: str = Field(
        min_length=1,
        max_length=256,
        description='A composition such as "Hex5HexNAc4Fuc1". Residues in any order; dHex reads'
        " as Fuc. A string that does not parse completely is refused rather than partly read.",
    )
    adduct: str = Field(min_length=1, max_length=64, description='For example "[M+H]+".')
    charge: int = Field(description="Signed charge state, e.g. 1 or -1.")
    client_reference: str | None = Field(
        default=None,
        max_length=128,
        description="The caller's own idempotency handle. A second create with the same"
        " reference is REFUSED with 409 and the id of the prediction that stands, because"
        " returning a new prediction would silently duplicate and merging would mutate.",
    )

    @model_validator(mode="after")
    def the_charge_is_not_zero(self) -> Self:
        if self.charge == 0:
            raise ValueError("charge 0 is not an ion; give a signed charge such as 1 or -1")
        return self


class PredictionResponse(BaseModel):
    """A frozen prediction. Immutable once served: the same id returns the same bytes."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    prediction_id: str
    created_at: str
    composition: str
    adduct: str
    charge: int
    client_reference: str | None = None
    frozen: bool = True
    candidates_total: int
    classes_total: int
    tied_candidates: int
    largest_indistinguishable_class: int
    enumerator_built: int = 0
    enumerator_rejected: int = 0
    enumerator_truncated: bool = False
    bands: tuple[BandOut, ...] = ()
    classes: tuple[ClassOut, ...] = ()
    candidates: tuple[CandidateOut, ...] = ()
    confidence: ConfidenceOut
    coverage: CoverageOut
    rules: RuleAccountingOut
    decision: DecisionOut
    ccs_evidence: Mapping[str, object] = Field(
        default_factory=dict,
        description="What the platform knows about a cross section for this composition and ion."
        " COMPOSITION- or ION-level, so it cannot separate candidates and says so. Never a"
        " prediction.",
    )
    is_a_ranking: bool = Field(
        description="False when the platform refuses to order the set. The candidates are still"
        " returned; what is withheld is an order, not an answer."
    )
    refusal: str | None = None
    weaknesses: tuple[str, ...] = ()
    stamp: ModelStamp


class MeasurementIn(BaseModel):
    """One experimental measurement a caller attaches. Never merged into the prediction."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    ccs: float = Field(gt=0, description="The measured cross section, in square angstroms.")
    uncertainty: float = Field(
        ge=0,
        description="As measured. REQUIRED: hard constraint 7 refuses a number without a type,"
        " and a type without a number is the same gap wearing a label.",
    )
    uncertainty_type: str = Field(
        min_length=1, description='SD, 2SD, SEM or CI95. Never "unknown".'
    )
    adduct: str = Field(min_length=1, max_length=64)
    charge: int
    ims_type: str = Field(min_length=1, description="DTIMS, TWIMS, TIMS or cyclic.")
    drift_gas: str = Field(
        min_length=1, description='The gas the value refers to. Never "UNSTATED".'
    )
    calibrant: str | None = None
    instrument: str | None = None
    source: str = Field(min_length=1, description="Who measured it, or which run it came from.")
    measured_on: str | None = None
    replicates: int | None = Field(default=None, ge=1)

    @model_validator(mode="after")
    def the_number_carries_its_type_and_its_gas(self) -> Self:
        if self.uncertainty_type.strip().casefold() in {"", "unknown", "none", "n/a"}:
            raise ValueError(
                f"uncertainty_type {self.uncertainty_type!r} is not a type. Hard constraint 7"
                " refuses a number without one, and a comparison of intervals cannot be made"
                " from a value whose spread is unstated"
            )
        if self.drift_gas.strip().casefold() in {"", "unstated", "unknown", "n/a"}:
            raise ValueError(
                f"drift_gas {self.drift_gas!r} is not a gas. A cross section measured against an"
                " unstated gas is not comparable with one measured against a stated gas, which"
                " is hard constraint 4"
            )
        if self.charge == 0:
            raise ValueError("charge 0 is not an ion")
        return self


class ValidationRequest(BaseModel):
    """Experimental measurements to attach to a frozen prediction."""

    model_config = ConfigDict(extra="forbid")

    measurements: tuple[MeasurementIn, ...] = Field(min_length=1, max_length=200)
    note: str | None = Field(default=None, max_length=2000)


class ValidationResponse(BaseModel):
    """What was attached, and proof that the prediction was not touched."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    validation_id: str
    prediction_id: str
    created_at: str
    measurements_attached: int
    validations_on_this_prediction: int = Field(
        description="Including this one. Attaching is APPEND ONLY: a second attach adds a"
        " record, it does not replace the first."
    )
    prediction_digest_before: str
    prediction_digest_after: str
    prediction_unchanged: bool = Field(
        description="Read from the stored payload before and after the attach, not asserted."
    )
    stamp: ModelStamp


class ComparisonState(StrEnum):
    """Why a comparison could or could not be made. Four reasons, never collapsed into one."""

    COMPARED = "compared"
    NO_PREDICTED_VALUE = "no_predicted_value"
    REFERENCE_HELD_NOT_RELEASABLE = "reference_held_not_releasable"
    NO_REFERENCE_IN_CORPUS = "no_reference_in_corpus"
    REFERENCE_NOT_CONSULTED = "reference_not_consulted"


class MeasurementComparison(BaseModel):
    """One attached measurement, against whatever could be compared with it."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    ccs: float
    uncertainty: float
    uncertainty_type: str
    adduct: str
    charge: int
    ims_type: str
    drift_gas: str
    source: str
    against_reference: ComparisonState
    reference_ccs: float | None = None
    delta_ccs: float | None = None
    delta_percent: float | None = None
    within_agreement_limit: bool | None = None
    unevaluable_because: str | None = None


class AttachedSpread(BaseModel):
    """How far the caller's own measurements sit from each other, per ion.

    The one comparison this platform can always make, and it is worth making: it is the
    caller's own reproducibility, and it needs no model and no reference.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    adduct: str
    charge: int
    ims_type: str
    drift_gas: str
    measurements: int
    lowest: float
    highest: float
    spread: float
    spread_percent: float
    within_agreement_limit: bool


class ComparisonResponse(BaseModel):
    """Delta CCS, interval coverage and status - each with the reason it is what it is."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    prediction_id: str
    composition: str
    adduct: str
    charge: int
    validations: int
    measurements: int
    status: ComparisonState
    # Always NO_PREDICTED_VALUE, and it says why rather than leaving a null.
    against_prediction: ComparisonState = ComparisonState.NO_PREDICTED_VALUE
    against_prediction_because: str = NO_PREDICTED_VALUE
    interval_coverage: str | None = None
    interval_coverage_because: str = INTERVAL_UNEVALUABLE
    agreement_limit_percent: float = AGREEMENT_LIMIT_PERCENT
    agreement_limit_is_policy: bool = True
    per_measurement: tuple[MeasurementComparison, ...] = ()
    among_attached: tuple[AttachedSpread, ...] = ()
    stamp: ModelStamp


class RunOut(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    run_id: str
    created_at: str
    kind: str
    prediction_id: str | None = None
    composition: str | None = None
    fingerprint: str
    data_snapshot: str
    detail: str


class RunsResponse(BaseModel):
    """A page of history. `total` is what the page is a page OF, so a caller can page honestly."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    total: int
    returned: int
    limit: int
    offset: int
    filters: Mapping[str, object] = Field(default_factory=dict)
    runs: tuple[RunOut, ...] = ()
    stamp: ModelStamp


class DomainOut(BaseModel):
    """The applicability domain: what this platform can answer for, and what it refuses."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    glycan_class: str = "N-linked"
    builds_on: str = (
        "the complete branched Man3GlcNAc2 core. A truncated core cannot be represented, so"
        " paucimannose species, endoglycosidase products keeping one core GlcNAc, and"
        " degradation products are out of reach BY CONSTRUCTION rather than by rule"
    )
    residues_supported: tuple[str, ...] = ()
    biosynthetic_rules: int = Field(
        default=0,
        description="The SIZE OF THE CURATED TABLE the enumerator loads, which is 15. Eleven of"
        " those reach N-glycan enumeration and four never do, so this is the scheme size and NOT"
        " the number of rules that can bear on a candidate. See `rules_govern`.",
    )
    rules_govern: str = (
        "ELEVEN of the 15 curated rules reach N-glycan enumeration: nine MGAT rules on branching"
        " order and bisecting interference, FUT8 on core fucosylation, and one class-agnostic blood"
        " group rule. The other four are O-glycan core rules and never apply here. NOT ONE of the"
        " eleven constrains galactosylation type, fucose position, chain extension or LacdiNAc -"
        " which is what the candidates for one composition differ in, and is why the curated rules"
        " cannot order a candidate set. (This field said all fifteen govern MGAT branching and"
        " bisecting until 27 September 2026. It was the last of five copies of that sentence and"
        " the only one served over HTTP.)"
    )
    species_assumption: str = (
        "mammalian. Insect and plant N-glycans carrying a core alpha1,3-fucose are outside it,"
        " and the enumerator places core fucose only at alpha1,6"
    )
    predicts_ccs: bool = False
    predicts_ccs_note: str = (
        "No. This platform ranks candidate structures on evidence. There is no glycan CCS model"
        " and V1 will not have one, so every cross section in a response is a MEASURED value"
    )
    holds_a_measured_cross_section_for_any_candidate: bool = False
    measured_reference_is_unreachable_note: str = (
        "AND IT HOLDS NONE. Stated plainly because the previous field only rules out PREDICTED"
        " values: the measured-reference path is UNREACHABLE for every candidate set in this"
        " release, not merely thin. All 24 glycan cross sections that clear every gate the loader"
        " applies are milk oligosaccharides recording no composition, so none of them can be keyed"
        " to a candidate set, and the 89 Struwe 2015 N-glycan values are held on a drift gas and an"
        " uncertainty type absent from the source. The LICENCE gate refuses none of them - all 117"
        " rows are open_attribution - so permission is not what is missing. Every response therefore carries a STATED"
        " ABSENCE of CCS evidence rather than a value, and that absence is the honest output"
        " rather than a gap in the wiring. What would change it: a cleared record carrying a"
        " composition. What would make such a value SEPARATE candidates is a further step - a"
        " measurement resolving to one structure rather than to a composition and an ion"
    )
    known_coverage_limit: str = (
        "the enumerator misses between a third and four fifths of the fully-resolved reference"
        " structures of a given composition - 17 of 34 for Hex5HexNAc4Fuc1, 22 of 28 for"
        " Hex5HexNAc2. Diagnosed: the mammalian scope, structures the curated rules deliberately"
        " reject, and five monolinks the placement vocabulary cannot make"
    )


class CurrentModelResponse(BaseModel):
    """What is serving, what it was built from, and what it will answer for."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    stamp: ModelStamp
    ccs_model_fitted: bool = False
    ccs_model_note: str = (
        "There is no fitted glycan CCS model. The fingerprint above is of a DETERMINISTIC"
        " PIPELINE - two curated tables, one reference corpus and one declared policy constant -"
        " which is what makes a past prediction reconstructible"
    )
    calibration: Calibration = Calibration.NEVER_CALIBRATED
    what_would_calibrate_it: str
    decision_values: tuple[str, ...] = ()
    decision_reachable_today: tuple[str, ...] = Field(
        default=(),
        description="Which of `decision_values` any input can actually produce. Derived from the"
        " gates in `ranking.decision_reachability()`, not listed by hand.",
    )
    decision_unreachable_today: Mapping[str, str] = Field(
        default_factory=dict,
        description="The values nothing can reach today, each against what WOULD reach it. Served"
        " so a consumer who finds two enum values that never appear is told why here rather than"
        " concluding the field is a stub.",
    )
    domain: DomainOut
    predictions_frozen: int = 0
