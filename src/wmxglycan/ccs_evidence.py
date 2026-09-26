"""CCS as evidence, never as a prediction, and never attached to a structure it was not measured on.

THE OWNER'S DECISION OF 26 SEPTEMBER 2026: there is no glycan CCS model and V1 will
not have one. So a cross section may enter a ranked result in exactly two ways -
as a measured reference we actually hold, or as a stated absence - and in no
other way at all. Nothing here predicts anything.

WHY EVIDENCE ATTACHES TO THE SET AND NOT TO A CANDIDATE

Every candidate for a composition IS that composition and that ion: they are
isomers. A reference cross section keyed on composition and ion therefore
matches all of them equally, and attaching it to each would render ONE
measurement as up to 167 confirmations of 167 different structures. That is a
value presented as belonging to something it was not measured on, which is the
thing hard constraint 3 exists to forbid, and it would make a candidate set the
platform cannot order look evidence-backed.

So `evidence_level` is part of the finding, and it is checked rather than
documented: a COMPOSITION-level or ION-level finding cannot be attached to an
individual candidate, and `for_candidate` raises if asked. STRUCTURE-level is
the only level that may, and it requires the cleared record itself to be fully
resolved with a canonical key equal to that candidate's.

WHAT THE FIVE STATES ARE FOR

The point of five rather than two is that "we looked and found nothing" is a
scientific claim about the literature, and three quite different things get
mistaken for it:

  MEASURED_REFERENCE      we hold a cleared value; here it is, with its provenance
  HELD_NOT_RELEASABLE     values exist for this ion and are blocked; the count and
                          the blockers travel, THE VALUE NEVER DOES
  NONE_IN_CORPUS_SEARCHED we looked, in a named corpus, and there is nothing
  NOT_CONSULTED           no evidence source was wired in. THE DEFAULT, so a ranker
                          without an adapter cannot report an absence it never looked for
  LOOKUP_FAILED           the lookup raised, or the key did not normalise. NEVER an
                          absence: a wiring bug must not read as a finding about the
                          literature

`NONE_IN_CORPUS_SEARCHED` cannot be constructed without naming the corpus and how
many records were consulted, because an absence with no denominator is not a
finding. That is the same asymmetry as instance Eleven in LIMITATIONS 4.5: a
field meaning "nothing to report" needs the same evidence as any other claim.

WHY A PROVENANCE GUARD AND NOT AN `is_measured_reference` FLAG

The first design carried a `Literal[True]` field so that "a prediction cannot be
placed here". A flag the producer sets about itself forbids only the one spelling
nobody would use, and it is satisfied by its own default - the same shape as
`extra="forbid"` rejecting unknown keys while missing ones sail through. So the
provenance IS the guard: a reference requires a source, a DOI, an uncertainty
value with a type that is not "unknown", a licence status, and the full
measurement conditions. A prediction cannot supply a DOI and an uncertainty
type, so the schema does the work and no flag is trusted.
"""

from __future__ import annotations

from enum import StrEnum
from typing import Protocol, Self, runtime_checkable

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .composition import Composition
from .reuse import (
    ReuseStatus,
    can_redistribute,
    can_train_commercial,
    is_inference_only,
)

# Refusal wording is an interface: written once, imported by tests rather than matched as prose.
WRONG_LEVEL_FOR_A_CANDIDATE = (
    "this CCS finding is {level} evidence and cannot be attached to an individual candidate."
    " Every candidate of a composition is that composition and that ion, so attaching one"
    " measurement to each would present it as a separate confirmation of each of them - a value"
    " claimed for something it was not measured on. Only STRUCTURE-level evidence, whose own"
    " canonical key equals the candidate's, may be attached to one candidate"
)
NO_CORPUS_NAMED = (
    "a searched-absence finding must name the corpus it searched and how many records it"
    " consulted, and the count must be above zero. An absence with no denominator is not a"
    " finding about the literature, it is an unexamined field, and the two read identically"
    " downstream; a denominator of zero is the same thing wearing a number"
)
NO_PERMITTED_USE = (
    "a measured reference cannot be shown under reuse status {status!r}, which has NO PERMITTED"
    " USE AT ALL - not training, not redistribution, and not the inference-side reference use"
    " this field is for. The three such statuses are `unverified` (nobody has read the terms),"
    " `excluded` (read, and deliberately kept out) and `non_commercial_no_derivatives` (the"
    " no-derivatives clause forbids sharing adapted material). The gate is default-deny and"
    " attaching a value to a ranked candidate is showing it"
)
FLAG_WITHOUT_THE_LEVEL = (
    "this finding is {level} evidence and cannot also claim to discriminate between candidates."
    " Every candidate of a composition is that composition and that ion, so a value keyed on"
    " either is shared by all of them equally. Only STRUCTURE-level evidence can separate"
    " candidates, and the flag is refused here rather than trusted downstream because a flag a"
    " producer sets about itself is exactly the guard that was removed from CCSReference"
)
UNCERTAINTY_TYPE_UNKNOWN = (
    "a measured reference needs an uncertainty with a stated type. This one reports"
    " uncertainty_type={given!r}, and a number without a type is refused - it is hard"
    " constraint 7, and it is the reason the 89 Struwe 2015 N-glycan values are held rather"
    " than shown"
)


class CCSEvidenceState(StrEnum):
    """What is known about a cross section for this composition and ion."""

    MEASURED_REFERENCE = "measured_reference"
    HELD_NOT_RELEASABLE = "held_not_releasable"
    NONE_IN_CORPUS_SEARCHED = "none_in_corpus_searched"
    NOT_CONSULTED = "not_consulted"
    LOOKUP_FAILED = "lookup_failed"


class EvidenceLevel(StrEnum):
    """What the evidence is evidence ABOUT. Decides whether it may touch one candidate."""

    COMPOSITION = "composition"
    ION = "ion"
    STRUCTURE = "structure"


class CCSReference(BaseModel):
    """A measured cross section we hold and may show, with everything that makes it mean something.

    Every field below is required because a prediction cannot supply them. That is the guard:
    there is no flag saying "this is measured", there is a set of obligations only a real
    measurement can meet.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    ccs: float = Field(description="The value as published. Never a harmonized or corrected value.")
    uncertainty: float = Field(description="As published. Required: hard constraint 7.")
    uncertainty_type: str = Field(description='SD, 2SD, SEM or CI95. Never "unknown".')
    adduct: str
    charge: int
    polarity: str
    ims_type: str
    drift_gas: str = Field(description='The gas the value refers to. Never "UNSTATED".')
    calibrant: str | None = None
    source: str
    doi: str = Field(description="Required. A value with no citable source is not a reference.")
    source_locator: str
    reuse_status: ReuseStatus
    # The structure the measurement was actually made on, where the source resolved one. None
    # means the source gave a composition and an ion and no structure, which is the usual case
    # and is exactly why evidence is normally COMPOSITION-level.
    measured_on_canonical_key: str | None = None

    @model_validator(mode="after")
    def the_provenance_is_the_guard(self) -> Self:
        if self.uncertainty_type.strip().casefold() in {"", "unknown", "none"}:
            raise ValueError(UNCERTAINTY_TYPE_UNKNOWN.format(given=self.uncertainty_type))
        if self.drift_gas.strip().casefold() in {"", "unstated", "unknown"}:
            raise ValueError(
                "a measured reference needs the gas its value refers to; an unstated gas makes the"
                " number incomparable, which is hard constraint 4"
            )
        if not self.doi.strip():
            raise ValueError("a measured reference needs a DOI")
        # THE LICENCE GATE APPLIES HERE, and it did not in the first version: every status was
        # accepted, including `excluded` and the `unverified` default. Requiring the FIELD to be
        # present is not the same as requiring it to permit anything, which is the asymmetry
        # this project keeps meeting - a schema that rejects the one spelling nobody uses.
        permitted = (
            can_train_commercial(self.reuse_status)
            or can_redistribute(self.reuse_status)
            or is_inference_only(self.reuse_status)
        )
        if not permitted:
            raise ValueError(NO_PERMITTED_USE.format(status=self.reuse_status.value))
        return self


class HeldValues(BaseModel):
    """Values that exist for this ion and may not be shown, with what blocks them.

    The count and the blockers travel. THE VALUES DO NOT. This is not squeamishness: the 89
    Struwe 2015 N-glycan values are blocked because the drift gas is UNSTATED and the
    uncertainty type is `unknown` with no value, so there is no number here that could be
    shown without inventing the thing that makes it mean something.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    records: int = Field(gt=0, description="How many records exist for this ion and are blocked.")
    blockers: tuple[str, ...] = Field(min_length=1)
    # Whether any single blocker, cleared alone, would release the values. Recorded because a
    # single-item blocker list invites the reading "resolve this one thing and you have data",
    # and for Struwe 2015 that reading is false: reading the paper for the gas leaves the
    # uncertainty exactly where it is, because it is absent from the source.
    any_blocker_resolves_alone: bool = False
    doi: str | None = None
    what_would_release_them: str


class CCSEvidence(BaseModel):
    """What the platform knows about a cross section for one composition and ion.

    Attached ONCE to a ranked result, never to a candidate. See the module docstring.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    state: CCSEvidenceState = CCSEvidenceState.NOT_CONSULTED
    level: EvidenceLevel = EvidenceLevel.COMPOSITION
    composition: str | None = None
    adduct: str | None = None
    charge: int | None = None
    reference: CCSReference | None = None
    held: HeldValues | None = None
    # Required by NONE_IN_CORPUS_SEARCHED: an absence needs a denominator.
    corpus_searched: str | None = None
    records_consulted: int | None = None
    # Required by LOOKUP_FAILED: what broke and what was asked for.
    failure: str | None = None
    key_attempted: str | None = None
    # Always true of COMPOSITION- and ION-level evidence, and said rather than left to be
    # inferred: a number shared by every candidate cannot separate them.
    discriminates_between_candidates: bool = False

    @model_validator(mode="after")
    def each_state_carries_what_it_claims(self) -> Self:
        state = self.state
        if state is CCSEvidenceState.MEASURED_REFERENCE and self.reference is None:
            raise ValueError("a measured-reference finding must carry the reference")
        if state is not CCSEvidenceState.MEASURED_REFERENCE and self.reference is not None:
            raise ValueError(f"a {state.value} finding must not carry a reference value")
        if state is CCSEvidenceState.HELD_NOT_RELEASABLE and self.held is None:
            raise ValueError("a held finding must say how many records are held and what blocks them")
        if state is not CCSEvidenceState.HELD_NOT_RELEASABLE and self.held is not None:
            raise ValueError(f"a {state.value} finding must not carry held values")
        if state is CCSEvidenceState.NONE_IN_CORPUS_SEARCHED and (
            not self.corpus_searched
            or self.records_consulted is None
            or self.records_consulted <= 0
        ):
            # A DENOMINATOR OF ZERO IS NOT A DENOMINATOR. The first version checked only for
            # None, so "we searched a corpus of 0 records and found nothing" was accepted as a
            # finding about the literature - which is exactly what an unsearched corpus looks
            # like too. `HeldValues.records` next door was already pinned gt=0; this was the
            # asymmetry between the two.
            raise ValueError(NO_CORPUS_NAMED)
        if state is CCSEvidenceState.LOOKUP_FAILED and not self.failure:
            raise ValueError(
                "a failed lookup must say what failed. An exception that degrades to an absence"
                " turns a wiring bug into a finding about the literature"
            )
        if (
            state is CCSEvidenceState.MEASURED_REFERENCE
            and self.level is EvidenceLevel.STRUCTURE
            and (self.reference is None or self.reference.measured_on_canonical_key is None)
        ):
            raise ValueError(
                "structure-level evidence must name the structure the measurement was made on"
            )
        if self.discriminates_between_candidates and self.level is not EvidenceLevel.STRUCTURE:
            raise ValueError(FLAG_WITHOUT_THE_LEVEL.format(level=self.level.value.upper()))
        return self

    def for_candidate(self, canonical_key: str) -> CCSReference:
        """The reference, if it really is evidence about THIS structure. Otherwise it raises.

        Deliberately a raise and not a None. A None would be absorbed by a caller that then
        showed nothing, and the mistake this prevents - one measurement rendered as a per
        structure claim - is silent by nature.
        """
        if self.level is not EvidenceLevel.STRUCTURE:
            raise ValueError(WRONG_LEVEL_FOR_A_CANDIDATE.format(level=self.level.value.upper()))
        if self.reference is None or self.reference.measured_on_canonical_key != canonical_key:
            raise ValueError(
                f"this reference was measured on {self.reference.measured_on_canonical_key!r}"
                f" and not on {canonical_key!r}"
            )
        return self.reference

    def summary(self) -> str:
        if self.state is CCSEvidenceState.MEASURED_REFERENCE and self.reference is not None:
            return (
                f"{self.level.value}-level measured reference: {self.reference.ccs} A^2"
                f" +/- {self.reference.uncertainty} ({self.reference.uncertainty_type}),"
                f" {self.reference.ims_type} in {self.reference.drift_gas},"
                f" {self.reference.source} doi:{self.reference.doi}"
            )
        if self.state is CCSEvidenceState.HELD_NOT_RELEASABLE and self.held is not None:
            return (
                f"{self.held.records} value(s) exist for this ion and are HELD, not shown:"
                f" {'; '.join(self.held.blockers)}. {self.held.what_would_release_them}"
            )
        if self.state is CCSEvidenceState.NONE_IN_CORPUS_SEARCHED:
            return (
                f"no cross section for this ion in {self.corpus_searched}"
                f" ({self.records_consulted} record(s) consulted)"
            )
        if self.state is CCSEvidenceState.LOOKUP_FAILED:
            return f"the CCS lookup FAILED and this is not an absence: {self.failure}"
        return (
            "no CCS evidence source was consulted. This is not a statement that none exists;"
            " nothing looked"
        )


@runtime_checkable
class CCSEvidenceLookup(Protocol):
    """What the ranker needs from whatever holds cross sections.

    A Protocol, so that `wmxglycan` states its requirement without importing anything that
    satisfies it. The CCS core is on the other side of a wall this package does not cross;
    an adapter outside both packages wires the two together. See tools/glycan_ccs_evidence.py.
    """

    def evidence_for(
        self, composition: Composition, adduct: str, charge: int
    ) -> CCSEvidence:  # pragma: no cover - a protocol has no body
        ...
