"""Ranking candidate N-glycans on evidence, and refusing to order what cannot be ordered.

WHAT THIS RANKS ON, AND WHAT IT REFUSES TO RANK ON

There is no glycan CCS model and V1 will not have one, so nothing here predicts a
cross section. Candidates are separated by ONE thing: whether the reference
corpus attests them. Everything else the platform knows is reported and is
deliberately kept out of the number, and each exclusion has a reason:

- THE CURATED RULES CANNOT ORDER THE SET. All fifteen govern MGAT branching order
  and bisecting interference, and the enumerator already rejected every candidate
  that breaks one. What is constant across survivors is the count of rules
  VIOLATED, which is zero - and this module VERIFIES that by running the check
  rather than trusting that the enumerator did. "Rules satisfied" is a different
  quantity and it measurably varies: 1 to 5 rules bear on a candidate, across
  41, 31, 50, 39 and 6 of the 167 Hex5HexNAc4Fuc1 candidates, and zero for all six
  Hex5HexNAc2 candidates. (An earlier draft of this paragraph said 1 to 3 across
  54, 83 and 30. That was measured BEFORE the ordering-note context gate was fixed,
  when 108 candidates carried a caveat instead of a satisfied rule, and it was
  carried in here without being re-measured. No test pinned it, which is why it
  survived; one does now.) It is not plausibility: the count of applicable rules
  is largely a count of how branched a candidate is, and no curated rule says a
  more branched structure is more likely. Reported in three separately named
  quantities, and none of them is a score.

- THE ORDERING CAVEATS ARE REPORTED AND NOT SCORED. Four of the nine FORBIDS
  rules describe assembly order, so they travel with a candidate as a caveat
  instead of excluding it. Turning "this needs a particular assembly order" into
  a number needs a weight, and the curated table gives none - FUT8's own
  rationale says the reverse order is "permitted and common", which is evidence
  AGAINST treating it as a penalty.

THE UNIT OF RANKING IS THE INDISTINGUISHABLE CLASS, NOT THE CANDIDATE

Two candidates whose 24 structure feature columns are identical cannot be
separated by anything this platform computes. Measured on Hex5HexNAc4Fuc1: 167
candidates form 61 such classes, 46 of them holding more than one candidate and
the largest holding 10, so 152 of 167 candidates have a twin. Ordering inside a
class would be arbitrary, and an arbitrary order rendered as a ranking is the
failure this project is built against.

So the class is the unit and attestation is POOLED to it. That is not a
concession, it is the correct arithmetic: members are indistinguishable, so an
observation attesting one member is an observation of the class. Pooling also
restores the discrimination that per-candidate scoring destroys - Hex5HexNAc4Fuc1
goes from a flat 17-against-150 split to bands of 1, 13 and 47 classes, the top
band holding four distinct attested structures.

A band is a `frozenset` of class ids and a class holds a `frozenset` of candidate
keys, so an order WITHIN either is unrepresentable rather than merely
undocumented. The bands themselves are a tuple, because between bands there is a
real ordering.

THE KEY IS THE None-BEARING ROW, NEVER `as_row()`

`FeatureVector.as_row()` maps an absent value to the module-level `math.nan`
singleton, and `(1.0, math.nan) == (1.0, math.nan)` is True only because tuple
comparison short-circuits on identity: `float("nan")` compares False against
itself, as does a numpy NaN, and a JSON round trip produces a fresh one. Keyed on
`as_row()`, every class would silently become a singleton the moment anyone wrote
`float("nan")` or routed a row through numpy, and the platform would start
ordering 152 of 167 candidates it cannot distinguish with no test failing. So the
key is built from the None-bearing mapping, where `None == None` is identity-true
and robust. This is a latent defect in the ported featuriser, not in this module;
see docs/GLYCAN_LIMITATIONS.md.

THE NUMBER, AND WHY IT IS NOT A PROBABILITY

The owner asked for a normalised confidence, and this module returns one. It is
NOT a probability that a candidate is correct, and three separate things stop it
being one:

1. THE HYPOTHESIS SPACE IS NOT THE CANDIDATE SET. Normalising over the candidates
   alone asserts the answer is among them, and that is measured false between a
   third and four fifths of the time: 17 of 34 fully-resolved reference
   structures of Hex5HexNAc4Fuc1 are not among the 167 candidates, and 22 of 28
   for Hex5HexNAc2. So every reference structure of the composition that no
   candidate matches is carried as its own hypothesis, and the share of the mass
   sitting on them is reported as `share_not_enumerated`. For Hex5HexNAc2 that is
   about four fifths, which is the honest headline for that composition.

   AND THAT WAS NOT ENOUGH, WHICH IS WORTH RECORDING BECAUSE THE FIRST VERSION LOOKED
   FINISHED. Reserving mass only for reference structures the corpus KNOWS about still
   closed the world whenever the corpus happened to know of none that were missed:
   Hex6HexNAc3Fuc3 has exactly one fully-resolved reference structure, that structure IS
   among its candidates, and the shares therefore summed to exactly 1.0 with
   share_not_enumerated published as a measured 0.0%. On the strength of one deposited
   structure the response asserted the answer was in the set.

   So the space carries one further hypothesis, ALWAYS: A STRUCTURE NEITHER PROPOSED NOR
   IN THE CORPUS. It has no observations, so its weight is the prior pseudocount and
   nothing else, and its share is published as `share_not_proposed`. The corpus can show
   a candidate set is INCOMPLETE; it can never show one is COMPLETE, and a hypothesis
   space that omits the unproposed is a space that claims it can.

2. THE PRIOR IS A CHOSEN CONSTANT AND IS PUBLISHED AS ONE. `PRIOR_PSEUDOCOUNT` is
   POLICY. An earlier draft of this module claimed add-one smoothing was "forced"
   by the uniform Dirichlet prior; that is wrong, because uniformity on the
   simplex is a choice of coordinates, and Jeffreys (1/2), Perks (1/K) and
   Haldane (0) are each standard and each published. The choice is also the single
   largest lever in the output. So the share is reported under ALL FOUR, in
   `share_under_priors`, in every response - a field that ships four numbers for
   one quantity cannot be mistaken for the probability.

3. THE EVENT IS A LITERATURE DEPOSITION, NOT A MOLECULE IN A SAMPLE. SugarBase
   records that a structure has been reported, with no abundance, no tissue and
   for most rows no species. `SHARE_MEANS` states this and travels with the
   result.

A share is UNDEFINED rather than zero under Haldane's prior when nothing in the space
is attested: the concentration and the observations are both zero, so the denominator
is. `None` is published for that prior in that case. A 0.0 there would make one column
appear to sum to nothing while the other three sum to one.

CALIBRATION HAS NEVER HAPPENED and cannot happen in this repository: it needs
compositions whose true structure is independently known, scored blind, and this
repository holds none.

REFUSE ON INEVALUABILITY, WARN ON WEAKNESS - hard constraint 8, applied exactly.
One band means nothing distinguishes any class from any other, so the result is
REFUSED as a ranking (Hex5HexNAc2 is the standing case: six candidates, six
classes, every one attested by exactly one structure, one band). Bands that exist
but rest on a handful of structures are a WEAKNESS, reported and not refused.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Iterable, Mapping, Sequence

from .attestation import AttestationIndex, default_attestation_index
from .ccs_evidence import CCSEvidence, CCSEvidenceLookup, CCSEvidenceState, EvidenceLevel
from .composition import Composition
from .enumeration import Enumerator, EnumerationResult, GlycanCandidate, is_order_constraint
from .features import extract_structure
from .glycan_graph import GlycanGraph

# --- policy, declared as policy ------------------------------------------------------------------

# The prior pseudocount. POLICY, CHOSEN RATHER THAN MEASURED, and the reason it is named here
# rather than written inline is that it is the one free parameter in the whole module and it
# moves the answer more than the data does. Every response publishes the share under all four
# standard alternatives so the lever is visible; see `PRIORS`.
PRIOR_PSEUDOCOUNT = 1.0

# The four standard non-informative Dirichlet concentrations. Perks is 1/K and so depends on the
# size of the hypothesis space, which is why it is a callable rather than a number.
PRIORS: Mapping[str, object] = {
    "haldane_0": 0.0,
    "perks_1_over_k": "1/K",
    "jeffreys_0.5": 0.5,
    "uniform_1": 1.0,
}

# Which of the four PRIORS the headline `evidence_share` uses. Named, and carried on every
# response, because identifying it by comparing floats against the four published columns is
# ambiguous the moment two of them coincide - which they do whenever the space is small.
POLICY_PRIOR = "uniform_1"

SHARE_MEANS = (
    "evidence_share is the share of a weighted hypothesis space held by this indistinguishable"
    " class, where the space is every candidate class PLUS every fully-resolved reference"
    " structure of this composition that no candidate matches, and the weight is the number of"
    " distinct reference structures attesting it plus a declared prior pseudocount. THE EVENT IT"
    " DESCRIBES IS A LITERATURE DEPOSITION, NOT A MOLECULE IN YOUR SAMPLE: SugarBase records that"
    " a structure has been reported, with no abundance, no tissue and for most rows no species."
    " It is NOT the probability that this class is the structure you measured, and it has never"
    " been calibrated. It also MOVES WITH THE PRIOR, and by how much depends on the row: the"
    " best-attested class of Hex5HexNAc4Fuc1 ranges over about a factor of three across the four"
    " standard priors, while an UNATTESTED class ranges from zero upwards, which is an unbounded"
    " ratio rather than a factor of three. All four are reported in share_under_priors so the"
    " size of the lever can be read off rather than taken on trust."
)
WHAT_WOULD_CALIBRATE_IT = (
    "a set of compositions whose true structure is independently known, enumerated and scored"
    " blind, with the realised frequency of the true structure landing in each band compared"
    " against that band's share. This repository holds no such set, and nothing in it can"
    " construct one: the 24 cleared glycan CCS records are milk oligosaccharides on a lactose"
    " core and the 89 N-glycan records are held. Until then the number is uncalibrated, and that"
    " is a statement about the evidence and not a promise about a later release"
)

# Refusal wording is an interface: written once, imported by tests rather than matched as prose.
ONE_BAND = (
    "RANKING REFUSED: all {classes} distinguishable class(es) of {composition} carry the same"
    " evidence ({structures} distinct attesting reference structure(s) each), so nothing"
    " separates any of them from any other. An order over them would be arbitrary, and an"
    " arbitrary order presented as a ranking is the failure this platform is built against. The"
    " candidates are returned as an unordered set"
)
NOTHING_TO_RANK = "RANKING REFUSED: the enumerator produced no candidate for {composition}: {reason}"
RULES_BROKEN = (
    "RANKING REFUSED: {violations} curated rule violation(s) were found among the candidates for"
    " {composition} by running the rule check over them. The enumerator is supposed to have"
    " rejected every violator, so this is not a statement about the glycans - it is the"
    " enumerator and the rule check disagreeing, and nothing should be ranked while they do"
)
ONLY_ONE_CANDIDATE = (
    "NOT A RANKING: the enumerator produced exactly one candidate for {composition}, so there is"
    " no order to report and none is being withheld. The candidate and its evidence are"
    " returned. One candidate is not a claim that the structure is settled: see the coverage"
    " block, which cannot establish that the enumerator proposed everything it should have"
)


class Decision(StrEnum):
    """What the platform says should happen next. The spec's three values.

    **TWO OF THESE THREE ARE UNREACHABLE TODAY, AND THAT IS NOT A BUG IN THIS ENUM.** Read this
    before concluding the decision field is a stub: `IM_VALIDATION_REQUIRED` is the only value any
    input can produce, because two of the published rules in `decide()` fire on EVERY input.

        AI_ONLY                    UNREACHABLE. Two independent gates, and both would have to
                                   change: no model here has been validated against known truth
                                   (`a_validated_model_exists()`), and the completeness of a
                                   candidate set cannot be established from a corpus of
                                   depositions (`CORPUS_CAN_ESTABLISH_COMPLETENESS`).
        IM_VALIDATION_RECOMMENDED  UNREACHABLE. Blocked by the completeness gate alone.
        IM_VALIDATION_REQUIRED     The value served for every request in this release.

    **What would reach them** is in `decision_reachability()` below, which DERIVES the answer by
    reading the same two gates `decide()` reads instead of restating it in prose. So the day either
    gate changes, the reachability this module publishes changes with it - and the sentence above
    cannot quietly become false, because `tests/test_glycan_ranking.py` patches both gates and
    asserts the other two values then appear. That is also how the branches are shown to be LIVE
    CODE rather than an enum member nothing references.

    Do not "fix" this by relaxing a threshold. A platform with no validated CCS model asking for
    instrument validation on every answer is the correct output, and it is what the specification
    means by identifying predictions that require experimental validation.
    """

    AI_ONLY = "AI_ONLY"
    IM_VALIDATION_RECOMMENDED = "IM_VALIDATION_RECOMMENDED"
    IM_VALIDATION_REQUIRED = "IM_VALIDATION_REQUIRED"


class Calibration(StrEnum):
    NEVER_CALIBRATED = "never_calibrated"
    CALIBRATED = "calibrated"  # reachable in principle; see WHAT_WOULD_CALIBRATE_IT


class CoverageState(StrEnum):
    """Whether the completeness of the candidate set could be checked at all."""

    MEASURED = "measured"
    UNEVALUABLE = "unevaluable"


class SetCompleteness(StrEnum):
    """What the corpus can say about whether the answer is among the candidates.

    THERE IS NO `COMPLETE` MEMBER, and its absence is the point. The corpus can show that a
    candidate set MISSES a structure; nothing in it can show that a set misses none, because
    the corpus is a record of what has been deposited and not a census of what exists. A
    version of this enum with a third member would let a composition with one deposited
    structure - Hex6HexNAc3Fuc3 has exactly one - be reported as settled.
    """

    # Measured: the corpus holds fully-resolved structures of this composition that no
    # candidate matches. The strongest statement available, and it is a negative one.
    INCOMPLETE = "incomplete"
    # Either nothing is missing from what the corpus knows, or the corpus knows nothing of
    # this composition. Neither establishes completeness.
    NOT_ESTABLISHABLE = "not_establishable"


# THE FIRST OF THE TWO GATES THAT MAKE AI_ONLY UNREACHABLE. A named module-level fact rather than
# a condition buried in a branch, so that a test can patch it and prove the branch is LIVE CODE -
# "AI_ONLY never appears" is equally satisfied by dead code, a misspelled comparison or an enum
# member nothing references, and this project has met that shape repeatedly (LIMITATIONS 4.5
# counts them; the count is kept in one place on purpose).
VALIDATED_MODEL: object | None = None

# THE SECOND GATE, and the one that also blocks IM_VALIDATION_RECOMMENDED. Named here for the same
# reason as the first: `Coverage.truth_may_not_be_in_the_candidate_set` READS THIS rather than
# returning a bare True, so the fact lives in one place and patching it moves both that property
# and `decision_reachability()` together. A corpus records what has been deposited; nothing in it
# can show that a candidate set misses nothing. See SetCompleteness, which has no COMPLETE member
# for the same reason.
CORPUS_CAN_ESTABLISH_COMPLETENESS = False


def a_validated_model_exists() -> bool:
    """False, and it has never been anything else. See VALIDATED_MODEL."""
    return VALIDATED_MODEL is not None


@dataclass(frozen=True)
class DecisionReachability:
    """Whether one `Decision` value can be produced today, and what would change that."""

    decision: Decision
    reachable_today: bool
    blocked_by: tuple[str, ...]
    what_would_reach_it: str | None


def decision_reachability() -> tuple[DecisionReachability, ...]:
    """Which of the three decisions any input can produce, DERIVED from the two gates.

    WHY DERIVED. The served list of reachable values was a hand-written tuple in `api.py` until
    27 September 2026, which is the shape LIMITATIONS 4.5 collects: a statement about what the
    code can do, maintained by hand next to the code that does it. This reads the same two
    predicates `decide()` reads, so a gate that opens cannot leave a published claim behind.

    THE LOGIC, which is just `decide()` read backwards. Every published rule maps to
    IM_VALIDATION_REQUIRED and `decide()` returns it the moment ANY rule fires. Two rules fire on
    every input while the gates hold: "the candidate set may not contain the answer" (the
    completeness gate) and "no validated model exists". So while the completeness gate holds,
    nothing else is reachable at all; if it opened, a single-class set with structure-level
    evidence would reach IM_VALIDATION_RECOMMENDED, and AI_ONLY additionally needs a validated
    model.
    """
    validated = a_validated_model_exists()
    completeness = CORPUS_CAN_ESTABLISH_COMPLETENESS

    no_model = (
        "no model here has been validated against independently known structures"
        " (ranking.VALIDATED_MODEL is None)"
    )
    no_completeness = (
        "the completeness of a candidate set is not establishable from a corpus of depositions"
        " (ranking.CORPUS_CAN_ESTABLISH_COMPLETENESS is False)"
    )

    ai_only_blockers = tuple(
        reason
        for reason, holds in ((no_model, not validated), (no_completeness, not completeness))
        if holds
    )
    recommended_blockers = (no_completeness,) if not completeness else ()

    return (
        DecisionReachability(
            decision=Decision.AI_ONLY,
            reachable_today=not ai_only_blockers,
            blocked_by=ai_only_blockers,
            what_would_reach_it=(
                None
                if not ai_only_blockers
                else (
                    "a model validated against structures known by another method, AND a way to"
                    " establish that a candidate set is complete. Both, not either: each blocks"
                    " this value on its own"
                )
            ),
        ),
        DecisionReachability(
            decision=Decision.IM_VALIDATION_RECOMMENDED,
            reachable_today=not recommended_blockers,
            blocked_by=recommended_blockers,
            what_would_reach_it=(
                None
                if not recommended_blockers
                else (
                    "a source that establishes completeness rather than recording depositions."
                    " A validated model is NOT required for this value, only for AI_ONLY"
                )
            ),
        ),
        DecisionReachability(
            decision=Decision.IM_VALIDATION_REQUIRED,
            reachable_today=True,
            blocked_by=(),
            what_would_reach_it=None,
        ),
    )


# --- the pieces of a result -----------------------------------------------------------------------


@dataclass(frozen=True)
class RuleAccounting:
    """The three different numbers people mean by "rules satisfied", named apart.

    Written as three fields because the first version of this design reported one scalar and
    called it constant, and the underlying quantity measurably partitions the candidate set. A
    field documented as constant is a field a reader is told not to check.
    """

    rules_in_scheme: int
    # VERIFIED, not assumed: computed by running the rule check over every candidate rather
    # than trusting that the enumerator rejected the violators.
    rules_violated: int
    rules_violated_verified_by_running_the_check: bool = True
    # canonical key -> how many curated rules actually bear on that candidate. Varies.
    rules_applicable: Mapping[str, int] = field(default_factory=dict)
    ordering_caveats: Mapping[str, int] = field(default_factory=dict)

    @property
    def applicable_varies(self) -> bool:
        return len(set(self.rules_applicable.values())) > 1

    def summary(self) -> str:
        if not self.rules_violated_verified_by_running_the_check:
            return (
                "the rule check was NOT run: there were no candidates to run it over, so these"
                " counts are absent rather than zero"
            )
        counts = sorted(set(self.rules_applicable.values()))
        return (
            f"{self.rules_in_scheme} curated rules in the scheme; {self.rules_violated} violated"
            f" (verified by running the check, not inferred); rules bearing on a candidate:"
            f" {counts or [0]}"
            + ("  <- VARIES, and is deliberately not a score" if self.applicable_varies else "")
        )


@dataclass(frozen=True)
class IndistinguishableClass:
    """Candidates no feature this platform computes can separate, and their pooled evidence."""

    class_id: str
    members: frozenset[str]
    # Pooled: how many DISTINCT reference structures attest any member. Members are
    # indistinguishable, so an observation of one is an observation of the class.
    attested_structures: int
    attesting_accessions: tuple[str, ...] = ()
    # Rows behind those structures. 2 or 3 means SugarBase spelled one structure twice; reported
    # so the discrepancy stays visible, never scored.
    reference_rows: int = 0
    evidence_share: float | None = None
    share_under_priors: Mapping[str, float] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.members)

    @property
    def is_a_tie(self) -> bool:
        """More than one candidate, none of which can be put before another."""
        return len(self.members) > 1

    @property
    def rendered_members(self) -> tuple[str, ...]:
        """Members in a stable order FOR DISPLAY ONLY. Sorted, so it cannot be read as a ranking."""
        return tuple(sorted(self.members))


@dataclass(frozen=True)
class Band:
    """Classes holding equal evidence. Unordered within; ordered between."""

    rank: int
    classes: frozenset[str]
    attested_structures: int
    # PER CLASS, and named so. The first version called this `evidence_share`, the same name
    # the class carries, and printed it straight after the band's class and candidate counts -
    # so a consumer reading a band line took a per-class number as the band's mass and
    # understated band 3 of Hex5HexNAc4Fuc1 by a factor of 47.
    share_per_class: float | None = None
    band_total_share: float | None = None
    share_under_priors: Mapping[str, float | None] = field(default_factory=dict)
    candidates: int = 0

    @property
    def rendered_classes(self) -> tuple[str, ...]:
        return tuple(sorted(self.classes))


@dataclass(frozen=True)
class Coverage:
    """Whether the candidate set is even capable of containing the answer.

    The measured refutation of the one claim a normalised score would otherwise make.
    """

    composition: str
    state: CoverageState
    reference_rows: int = 0
    reference_structures: int = 0
    attested_by_a_candidate: int = 0
    not_enumerated: int = 0
    not_enumerated_examples: tuple[str, ...] = ()

    @property
    def completeness(self) -> SetCompleteness:
        """INCOMPLETE where the corpus proves something is missing; otherwise NOT_ESTABLISHABLE.

        Never COMPLETE. See SetCompleteness.
        """
        if self.not_enumerated > 0:
            return SetCompleteness.INCOMPLETE
        return SetCompleteness.NOT_ESTABLISHABLE

    @property
    def truth_may_not_be_in_the_candidate_set(self) -> bool:
        """ALWAYS TRUE, and it is a finding rather than a placeholder.

        The first version of this returned False whenever every reference structure of the
        composition happened to be among the candidates - which for Hex6HexNAc3Fuc3 meant
        asserting the answer was in the set on the strength of ONE deposited structure. The
        published decision rule then reported "the candidate set may not contain the answer:
        fires=False", the coverage weakness vanished from the report, and the shares summed to
        exactly 1.0 over the candidates.

        Completeness is not a thing this corpus can establish, so the honest value is a
        constant, and `completeness` carries the distinction that actually varies.

        IT READS `CORPUS_CAN_ESTABLISH_COMPLETENESS` rather than returning a bare True, so this
        fact and the reachability `decision_reachability()` publishes come from one place. A bare
        True here and a sentence about unreachability elsewhere is two copies of one claim.
        """
        return not CORPUS_CAN_ESTABLISH_COMPLETENESS

    def summary(self) -> str:
        if self.state is CoverageState.UNEVALUABLE:
            return (
                f"COVERAGE UNEVALUABLE for {self.composition}: the reference corpus holds no"
                " fully-resolved structure of this composition, so whether the candidates could"
                " contain the answer cannot be checked. This is NOT evidence that they do"
            )
        if self.not_enumerated == 0:
            return (
                f"all {self.reference_structures} fully-resolved reference structure(s) of"
                f" {self.composition} are among the candidates, and NOTHING FOLLOWS FROM THAT"
                " about completeness: the corpus records what has been deposited, not what"
                " exists, so the answer may still be a structure nobody has deposited"
            )
        return (
            f"{self.attested_by_a_candidate} of {self.reference_structures} fully-resolved"
            f" reference structures of {self.composition} are among the candidates;"
            f" {self.not_enumerated} are NOT, so the answer may not be in the set"
        )


@dataclass(frozen=True)
class DecisionRule:
    """One published reason the platform asks for ion-mobility validation.

    Published, with `applies_to`, for the same reason the CCS confidence rules are: a rule the
    caller cannot read is a rule the caller cannot check, and one of the CCS rules turned out to
    be inert for every ion not already in the corpus without anybody noticing.
    """

    name: str
    fires: bool
    applies_to: str
    because: str
    decision_if_it_fires: Decision


@dataclass(frozen=True)
class RankedSet:
    """The ranker's answer for one composition: bands, evidence, coverage and what to do next."""

    composition: str
    candidates_total: int
    classes_total: int
    bands: tuple[Band, ...]
    classes: Mapping[str, IndistinguishableClass]
    candidates: Mapping[str, GlycanCandidate]
    coverage: Coverage
    rule_accounting: RuleAccounting
    ccs_evidence: CCSEvidence
    decision: Decision
    decision_rules: tuple[DecisionRule, ...]
    calibration: Calibration
    share_means: str
    what_would_calibrate_it: str
    attestation_provenance: str
    # The licence and the attribution travel with the evidence, because the whole score derives
    # from SugarBase and a consumer rendering a ranked set otherwise has no licence statement
    # for the thing the ranking rests on. `AttestationIndex.summary()` already printed both.
    attestation_licence: str = ""
    attestation_attribution: str = ""
    enumerator_built: int = 0
    enumerator_rejected: int = 0
    # The enumerator hit its tree budget and stopped, so the candidate set is TRUNCATED and
    # every figure computed from it describes a fragment. Carried because the first version
    # dropped it: a capped result was published with a coverage line reading "N reference
    # structures are NOT among the candidates" as though the enumerator had considered and
    # not produced them, when it had simply stopped.
    enumerator_capped: bool = False
    # The named policy prior, so a consumer does not have to identify it by float comparison
    # against the four published columns - which is ambiguous whenever two priors coincide.
    prior: str = ""
    prior_pseudocount: float = 0.0
    share_not_enumerated: float | None = None
    # The catch-all: a structure neither proposed by the enumerator nor deposited in the
    # corpus. Never zero, because nothing can rule it out.
    share_not_proposed: float | None = None
    # Both masses under every published prior, because three of the four columns did not sum
    # to one while only the policy prior's complement was published.
    mass_under_priors: Mapping[str, Mapping[str, float | None]] = field(default_factory=dict)
    # POOLED OVER CLASSES, not candidates, and named so. 47 unattested CLASSES hold 127 of the
    # 167 Hex5HexNAc4Fuc1 candidates, so "share on unattested candidates: 42.0%" invited the
    # reading that three quarters of the set held 42% of the mass.
    share_on_unattested_classes: float | None = None
    weaknesses: tuple[str, ...] = ()
    refusal: str | None = None

    @property
    def is_a_ranking(self) -> bool:
        return self.refusal is None

    @property
    def tied_candidates(self) -> int:
        """Candidates sharing a class with at least one other, so no order over them exists."""
        return sum(len(one) for one in self.classes.values() if one.is_a_tie)

    @property
    def largest_tie(self) -> int:
        return max((len(one) for one in self.classes.values()), default=0)

    def summary(self) -> str:
        lines = [
            f"{self.composition}: {self.candidates_total} candidates in"
            f" {self.classes_total} indistinguishable class(es)",
            f"  {self.coverage.summary()}",
            f"  {self.rule_accounting.summary()}",
            f"  CCS: {self.ccs_evidence.summary()}",
            f"  decision: {self.decision.value}",
        ]
        if self.refusal:
            lines.append(f"  {self.refusal}")
        else:
            for band in self.bands:
                each = "-" if band.share_per_class is None else f"{band.share_per_class:7.3%}"
                total = "-" if band.band_total_share is None else f"{band.band_total_share:7.3%}"
                lines.append(
                    f"  band {band.rank}: {len(band.classes):3} class(es),"
                    f" {band.candidates:4} candidate(s), attested by"
                    f" {band.attested_structures} distinct structure(s),"
                    f" {each} per class, {total} for the band"
                )
        if self.share_not_enumerated is not None:
            unattested = (
                "undefined"
                if self.share_on_unattested_classes is None
                else f"{self.share_on_unattested_classes:.1%}"
            )
            lines.append(
                f"  mass on reference structures NOT enumerated: {self.share_not_enumerated:.1%}"
                f"   |   on a structure nobody proposed: {self.share_not_proposed:.2%}"
                f"   |   on unattested CLASSES: {unattested}"
            )
        lines.append(f"  prior: {self.prior} (pseudocount {self.prior_pseudocount})")
        lines.append(
            f"  evidence: {self.attestation_provenance}"
            f"  [{self.attestation_licence}, {self.attestation_attribution}]"
        )
        for weakness in self.weaknesses:
            lines.append(f"  WEAKNESS: {weakness}")
        lines.append(f"  calibration: {self.calibration.value}")
        return "\n".join(lines)


# --- building one ---------------------------------------------------------------------------------


def class_key_for(candidate: GlycanCandidate) -> tuple[tuple[str, object], ...]:
    """The 24 structure feature columns, in a form where two equal rows really compare equal.

    Built from the None-bearing mapping and never from `FeatureVector.as_row()`, which maps an
    absent value to `math.nan`. See the module docstring: `as_row()` equality holds by CPython
    identity and `float("nan")`, numpy and JSON each break it silently.
    """
    values = extract_structure(candidate.iupac_condensed)
    row = tuple(sorted((name, value) for name, value in values.items()))
    if any(isinstance(value, float) and value != value for _name, value in row):
        raise ValueError(
            "a class key contains a NaN, so two identical rows would compare unequal and every"
            " indistinguishable class would silently become a singleton. Absent values must"
            " arrive as None"
        )
    return row


def rank(
    result: EnumerationResult,
    *,
    index: AttestationIndex | None = None,
    enumerator: Enumerator | None = None,
    ccs: CCSEvidenceLookup | None = None,
    adduct: str | None = None,
    charge: int | None = None,
) -> RankedSet:
    """Band the candidates of one composition by the evidence that attests them.

    `ccs` is optional and its absence is reported as NOT_CONSULTED rather than as an absence of
    evidence; a lookup that raises is reported as LOOKUP_FAILED and never as an absence either.
    """
    index = default_attestation_index() if index is None else index
    composition_text = result.composition.canonical

    evidence = _ccs_evidence(result.composition, ccs, adduct, charge)

    if not result.candidates:
        return _refused(
            result,
            composition_text,
            NOTHING_TO_RANK.format(
                composition=composition_text, reason=result.refusal or "no reason given"
            ),
            evidence,
            index,
        )

    # --- the classes, and the evidence pooled to them ---------------------------------------------
    grouped: dict[tuple, list[GlycanCandidate]] = {}
    for candidate in sorted(result.candidates, key=lambda c: c.canonical_key):
        grouped.setdefault(class_key_for(candidate), []).append(candidate)

    classes: dict[str, IndistinguishableClass] = {}
    for members in grouped.values():
        keys = frozenset(c.canonical_key for c in members)
        class_id = min(keys)  # derived from the members, so it cannot depend on enumerator order
        attested = {key for key in keys if index.attests(key)}
        classes[class_id] = IndistinguishableClass(
            class_id=class_id,
            members=keys,
            attested_structures=len(attested),
            attesting_accessions=tuple(
                sorted({a for key in sorted(attested) for a in index.accessions_for(key)})
            ),
            reference_rows=sum(index.rows_for(key) for key in sorted(attested)),
        )

    # --- coverage: could the answer even be in here ------------------------------------------------
    candidate_keys = {c.canonical_key for c in result.candidates}
    reference_keys = index.structures_of_composition(result.composition)
    coverage = Coverage(
        composition=composition_text,
        state=CoverageState.MEASURED if reference_keys else CoverageState.UNEVALUABLE,
        reference_rows=index.rows_of_composition(result.composition),
        reference_structures=len(reference_keys),
        attested_by_a_candidate=len(reference_keys & candidate_keys),
        not_enumerated=len(reference_keys - candidate_keys),
        not_enumerated_examples=tuple(sorted(reference_keys - candidate_keys))[:5],
    )

    # --- the shares, over a hypothesis space that admits the answer may be elsewhere ---------------
    counts = [one.attested_structures for one in classes.values()]
    # Every reference structure no candidate matches is its own hypothesis, attested once;
    # plus ONE catch-all for a structure neither proposed nor deposited, which carries no
    # observation and exists so that nothing here can sum to 1 over the candidates alone.
    missing = coverage.not_enumerated
    space = len(counts) + missing + 1
    observations = sum(counts) + missing
    shares = {
        name: _shares(counts, missing, space, observations, alpha) for name, alpha in PRIORS.items()
    }
    policy = _shares(counts, missing, space, observations, PRIOR_PSEUDOCOUNT)

    classes = {
        class_id: IndistinguishableClass(
            class_id=one.class_id,
            members=one.members,
            attested_structures=one.attested_structures,
            attesting_accessions=one.attesting_accessions,
            reference_rows=one.reference_rows,
            evidence_share=policy["by_count"][one.attested_structures],
            share_under_priors={
                name: value["by_count"][one.attested_structures] for name, value in shares.items()
            },
        )
        for class_id, one in classes.items()
    }

    # --- bands ------------------------------------------------------------------------------------
    by_count: dict[int, list[IndistinguishableClass]] = {}
    for one in classes.values():
        by_count.setdefault(one.attested_structures, []).append(one)
    bands = tuple(
        Band(
            rank=position,
            classes=frozenset(one.class_id for one in group),
            attested_structures=count,
            share_per_class=policy["by_count"][count],
            band_total_share=(
                None
                if policy["by_count"][count] is None
                else policy["by_count"][count] * len(group)
            ),
            share_under_priors={name: value["by_count"][count] for name, value in shares.items()},
            candidates=sum(len(one) for one in group),
        )
        for position, (count, group) in enumerate(
            sorted(by_count.items(), key=lambda item: -item[0]), start=1
        )
    )

    accounting = _rule_accounting(result.candidates, enumerator)
    weaknesses = _weaknesses(bands, coverage, policy, classes, capped=result.capped)

    refusal = None
    if accounting.rules_violated:
        # MEASURED AND THEN ACTED ON. The first version computed this, published it, and did
        # nothing with it: a candidate set breaking curated biosynthetic rules was banded and
        # returned with no refusal and no weakness. The enumerator is supposed to have rejected
        # every violator, so a non-zero count here does not mean the candidates are unusual -
        # it means the enumerator and the rule check disagree, and nothing downstream should
        # act on a set in that state.
        refusal = RULES_BROKEN.format(
            violations=accounting.rules_violated, composition=composition_text
        )
    elif len(classes) == 1 and not any(one.is_a_tie for one in classes.values()):
        # One candidate, so there is nothing to order and nothing is being withheld. Refusing
        # it with the one-band wording would have said "nothing separates any of them from any
        # other" about a set of one, which reads as a failure where there is none.
        refusal = ONLY_ONE_CANDIDATE.format(composition=composition_text)
    elif len(bands) == 1:
        refusal = ONE_BAND.format(
            classes=len(classes),
            composition=composition_text,
            structures=bands[0].attested_structures,
        )

    decision, rules = _decide(
        classes=classes, bands=bands, coverage=coverage, evidence=evidence, refused=refusal is not None
    )

    return RankedSet(
        composition=composition_text,
        candidates_total=len(result.candidates),
        classes_total=len(classes),
        bands=bands,
        classes=classes,
        # SORTED. The only ordered container in the result whose order the design does not
        # assert, and it arrived in the enumerator's depth-first order - so a consumer taking
        # the first entry got whichever branch the generator happened to walk first, beside a
        # refusal saying the candidates are returned as an unordered set.
        candidates={
            c.canonical_key: c
            for c in sorted(result.candidates, key=lambda one: one.canonical_key)
        },
        coverage=coverage,
        rule_accounting=accounting,
        ccs_evidence=evidence,
        decision=decision,
        decision_rules=rules,
        calibration=Calibration.NEVER_CALIBRATED,
        share_means=SHARE_MEANS,
        what_would_calibrate_it=WHAT_WOULD_CALIBRATE_IT,
        attestation_provenance=index.provenance,
        attestation_licence=index.dataset_version.licence,
        attestation_attribution=index.dataset_version.attribution,
        enumerator_built=result.built,
        enumerator_rejected=result.rejected,
        enumerator_capped=result.capped,
        prior=POLICY_PRIOR,
        prior_pseudocount=PRIOR_PSEUDOCOUNT,
        share_not_enumerated=policy["not_enumerated"],
        share_not_proposed=policy["not_proposed"],
        mass_under_priors={
            name: {
                "not_enumerated": value["not_enumerated"],
                "not_proposed": value["not_proposed"],
            }
            for name, value in shares.items()
        },
        # THREE CASES, AND THEY ARE DIFFERENT ANSWERS. No unattested class is a measured
        # ZERO - Hex5HexNAc2 has none, because all six of its classes are attested. An
        # unattested class whose share is undefined (Haldane with nothing observed) is NONE. A
        # first version collapsed the first case into the second and summary() then crashed
        # formatting it, which no test caught because no test had rendered a refused set.
        share_on_unattested_classes=_unattested_mass(policy["by_count"], counts),
        weaknesses=weaknesses,
        refusal=refusal,
    )


def _shares(
    counts: Sequence[int], missing: int, space: int, observations: int, alpha: object
) -> dict:
    """Shares under one prior concentration, over the WHOLE hypothesis space.

    `space` counts the candidate classes, the reference structures no candidate matches, AND
    the one catch-all hypothesis that the answer is a structure neither proposed nor deposited.
    That last one is why nothing here can sum to 1 over the candidates alone.
    """
    a = 1.0 / space if alpha == "1/K" else float(alpha)  # type: ignore[arg-type]
    total = observations + a * space
    if total <= 0:
        # Haldane's prior with nothing attested: every weight is zero and the share is
        # UNDEFINED, not zero. Publishing 0.0 would make this column read as a measured
        # near-impossibility for every hypothesis at once.
        return {
            "by_count": {count: None for count in set(counts) | {1}},
            "not_enumerated": None,
            "not_proposed": None,
        }
    return {
        "by_count": {count: (count + a) / total for count in set(counts) | {1}},
        # Every not-enumerated reference structure is attested exactly once.
        "not_enumerated": missing * (1.0 + a) / total,
        # The catch-all carries the prior and nothing else: nothing has ever been observed
        # that is neither a candidate nor in the corpus, which is not the same as nothing
        # being there.
        "not_proposed": a / total,
    }


def _unattested_mass(by_count: Mapping[int, float | None], counts: Sequence[int]) -> float | None:
    """Mass held by the classes nothing attests: 0.0 when there are none, None when undefined."""
    unattested = sum(1 for count in counts if count == 0)
    if unattested == 0:
        return 0.0
    each = by_count.get(0)
    return None if each is None else each * unattested


def _rule_accounting(
    candidates: Iterable[GlycanCandidate], enumerator: Enumerator | None
) -> RuleAccounting:
    candidates = tuple(candidates)
    working = enumerator if enumerator is not None else Enumerator()
    violated = 0
    applicable: dict[str, int] = {}
    caveats: dict[str, int] = {}
    for candidate in candidates:
        graph = GlycanGraph.from_iupac_condensed(candidate.iupac_condensed)
        # VERIFIED by running the check. The enumerator is supposed to have rejected every
        # violator; "supposed to" is not a measurement.
        violated += len(working.broken_rules(graph))
        applicable[candidate.canonical_key] = sum(
            1 for reason in candidate.rationale if reason.kind == "rule"
        )
        caveats[candidate.canonical_key] = sum(
            1 for reason in candidate.rationale if reason.kind == "ordering"
        )
    return RuleAccounting(
        rules_in_scheme=len(working.constraints),
        rules_violated=violated,
        rules_applicable=applicable,
        ordering_caveats=caveats,
    )


def _weaknesses(
    bands: Sequence[Band],
    coverage: Coverage,
    policy: Mapping[str, object],
    classes: Mapping[str, IndistinguishableClass],
    capped: bool = False,
) -> tuple[str, ...]:
    found: list[str] = []
    if capped:
        found.append(
            "THE CANDIDATE SET IS TRUNCATED: the enumerator hit its tree budget and stopped, so"
            " every figure here describes a fragment of the candidate space and the coverage"
            " line below is not a statement about what the rules permit"
        )
    top = bands[0] if bands else None
    if top is not None and coverage.state is CoverageState.MEASURED:
        found.append(
            f"the top band rests on {top.attested_structures} distinct reference structure(s) out"
            f" of {coverage.reference_structures} known for this composition. That is the whole"
            " of the discrimination"
        )
    if coverage.truth_may_not_be_in_the_candidate_set:
        found.append(coverage.summary())
    tied = sum(len(one) for one in classes.values() if one.is_a_tie)
    if tied:
        largest = max(len(one) for one in classes.values())
        found.append(
            f"{tied} candidate(s) share an indistinguishable class with at least one other and"
            f" cannot be ordered at all; the largest such class holds {largest}"
        )
    share = policy.get("not_enumerated")
    if isinstance(share, float) and share > 0:
        found.append(
            f"{share:.1%} of the evidence mass sits on reference structures the enumerator did"
            " not propose, so it is not on any candidate returned here"
        )
    return tuple(found)


def _ccs_evidence(
    composition: Composition, ccs: CCSEvidenceLookup | None, adduct: str | None, charge: int | None
) -> CCSEvidence:
    if ccs is not None and (adduct is None or charge is None):
        # A WIRING BUG, AND IT RAISES. The caller supplied an evidence source and then did not
        # say which ion to look up, so nothing was consulted - and reporting that as
        # NOT_CONSULTED would deliver a bug as the honest default, indistinguishable from a
        # deployment with no adapter at all. That is the one-value-for-two-reasons shape this
        # repository has met before.
        raise ValueError(
            "a CCS evidence source was supplied without both an adduct and a charge, so there"
            " was no ion to look up. Refused rather than reported as NOT_CONSULTED, which would"
            " present a wiring mistake as a deliberate absence"
        )
    if ccs is None:
        # NOT_CONSULTED, never an absence. A ranker with no evidence source has not looked, and
        # saying "nothing found" would be a claim about the literature that nothing supports.
        return CCSEvidence(
            state=CCSEvidenceState.NOT_CONSULTED, composition=composition.canonical
        )
    try:
        return ccs.evidence_for(composition, adduct, charge)
    except Exception as failure:  # noqa: BLE001 - any failure, and it must not become an absence
        return CCSEvidence(
            state=CCSEvidenceState.LOOKUP_FAILED,
            composition=composition.canonical,
            adduct=adduct,
            charge=charge,
            failure=f"{type(failure).__name__}: {failure}",
            key_attempted=f"{composition.canonical} {adduct} {charge:+d}",
        )


def decide(
    *,
    classes: Mapping[str, IndistinguishableClass],
    bands: Sequence[Band],
    coverage: Coverage,
    evidence: CCSEvidence,
    refused: bool,
) -> tuple[Decision, tuple[DecisionRule, ...]]:
    """Public: the decision and the published rules behind it.

    Exists because the service layer looks the CCS evidence up ONCE, freezes it into the
    prediction body, and must re-derive the decision from the evidence it froze - one of the
    published rules reads that evidence, so keeping the earlier decision would serve rules that
    disagree with the evidence printed beside them. Reaching into a private for that would make
    the service depend on something this module never promised to keep.
    """
    return _decide(
        classes=classes, bands=bands, coverage=coverage, evidence=evidence, refused=refused
    )


def _decide(
    *,
    classes: Mapping[str, IndistinguishableClass],
    bands: Sequence[Band],
    coverage: Coverage,
    evidence: CCSEvidence,
    refused: bool,
) -> tuple[Decision, tuple[DecisionRule, ...]]:
    """The decision and the published rules behind it. The worst outcome any rule reaches wins."""
    # THE LEVEL IS CHECKED HERE, not only the flag. The first version trusted
    # `discriminates_between_candidates`, which is a field a producer sets about itself - the
    # same shape as the `is_measured_reference` flag that was removed from CCSReference for
    # being self-asserted. A COMPOSITION-level finding with that flag set to True was accepted
    # by the validator and made the "no cross section is held for this structure" rule report
    # fires=False, which is the one rule standing between a shared measurement and a
    # per-structure claim. Both ends are fixed: the validator now forbids the combination and
    # this reads the level.
    structure_level = (
        evidence.state is CCSEvidenceState.MEASURED_REFERENCE
        and evidence.level is EvidenceLevel.STRUCTURE
        and evidence.discriminates_between_candidates
    )
    rules = (
        DecisionRule(
            name="more than one structure is possible",
            fires=len(classes) > 1 or any(one.is_a_tie for one in classes.values()),
            applies_to="EVERY COMPOSITION THE ENUMERATOR RETURNS MORE THAN ONE STRUCTURE FOR",
            because=(
                "the platform cannot name one structure, so any downstream use of a single named"
                " structure would be choosing one of several the evidence does not separate"
            ),
            decision_if_it_fires=Decision.IM_VALIDATION_REQUIRED,
        ),
        DecisionRule(
            name="no cross section is held for this structure",
            fires=not structure_level,
            applies_to="EVERY ION, because no measurement in this corpus resolves to a structure",
            because=(
                "a cross section keyed on composition and ion is shared by every isomer, so it"
                " cannot confirm one of them. Structure-level evidence is what would, and this"
                " repository holds none: the 24 cleared glycan records are milk oligosaccharides"
                " and the 89 N-glycan records are held"
            ),
            decision_if_it_fires=Decision.IM_VALIDATION_REQUIRED,
        ),
        DecisionRule(
            name="the candidate set may not contain the answer",
            fires=coverage.truth_may_not_be_in_the_candidate_set,
            applies_to=(
                "EVERY COMPOSITION WITHOUT EXCEPTION, and it fires on every one. The corpus can"
                " show a candidate set misses a structure and can never show it misses none, so"
                " there is no input for which this rule does not fire. It is kept as a rule"
                f" rather than folded into the prose because `completeness` reports which half"
                " applies: measured-incomplete ({completeness}) or not establishable"
            ).format(completeness=coverage.completeness.value),
            because=(
                "validating one of these candidates cannot settle the question if the real"
                " structure was never proposed"
            ),
            decision_if_it_fires=Decision.IM_VALIDATION_REQUIRED,
        ),
        DecisionRule(
            name="the bands do not separate",
            fires=refused,
            applies_to="A SET WHOSE CLASSES ALL CARRY THE SAME EVIDENCE",
            because="there is no ranking to act on, only an unordered set",
            decision_if_it_fires=Decision.IM_VALIDATION_REQUIRED,
        ),
        DecisionRule(
            name="no validated model exists",
            fires=not a_validated_model_exists(),
            applies_to="EVERY REQUEST, FOR AS LONG AS THERE IS NO VALIDATED MODEL",
            because=(
                "AI_ONLY asserts the platform's own answer needs no instrument. Nothing here has"
                " ever been validated against known truth, so that assertion has no basis. This"
                " is the rule that makes AI_ONLY unreachable, and it is one named predicate"
                " rather than a condition spread across branches so that it can be patched and"
                " the branch shown to be live code"
            ),
            decision_if_it_fires=Decision.IM_VALIDATION_REQUIRED,
        ),
    )
    if any(rule.fires for rule in rules):
        return Decision.IM_VALIDATION_REQUIRED, rules
    # One structure, structure-level evidence for it, a validated model, complete coverage.
    if structure_level and a_validated_model_exists():
        return Decision.AI_ONLY, rules
    return Decision.IM_VALIDATION_RECOMMENDED, rules


def _refused(
    result: EnumerationResult,
    composition_text: str,
    reason: str,
    evidence: CCSEvidence,
    index: AttestationIndex,
) -> RankedSet:
    # THE INDEX IS CONSULTED. The first version hardcoded an UNEVALUABLE coverage with
    # reference_structures=0, so a composition the enumerator declines was reported as one the
    # corpus knows nothing about. 49 of the 687 indexed compositions fail the enumerator's
    # complete-core gate, and Hex2HexNAc2Fuc2 is one of them - the corpus holds 15 distinct
    # fully-resolved structures of it, and the response said none. That is a fabricated
    # absence, and the fact it is attached to a refusal makes it more misleading rather than
    # less: a reader takes it as the reason nothing could be ranked.
    reference_keys = index.structures_of_composition(composition_text)
    coverage = Coverage(
        composition=composition_text,
        state=CoverageState.MEASURED if reference_keys else CoverageState.UNEVALUABLE,
        reference_rows=index.rows_of_composition(composition_text),
        reference_structures=len(reference_keys),
        attested_by_a_candidate=0,
        not_enumerated=len(reference_keys),
        not_enumerated_examples=tuple(sorted(reference_keys))[:5],
    )
    decision, rules = _decide(
        classes={}, bands=(), coverage=coverage, evidence=evidence, refused=True
    )
    return RankedSet(
        composition=composition_text,
        candidates_total=0,
        classes_total=0,
        bands=(),
        classes={},
        candidates={},
        coverage=coverage,
        # NOT ZERO, AND NOT "VERIFIED". The refused path used to publish
        # RuleAccounting(rules_in_scheme=0, rules_violated=0) with the verified flag left at
        # its default True, so a response read "0 curated rules in the scheme; 0 violated
        # (verified by running the check)" when nothing had been checked and the scheme holds
        # fifteen. There were no candidates to run the check over, so the counts are ABSENT.
        rule_accounting=RuleAccounting(
            rules_in_scheme=0,
            rules_violated=0,
            rules_violated_verified_by_running_the_check=False,
        ),
        ccs_evidence=evidence,
        decision=decision,
        decision_rules=rules,
        calibration=Calibration.NEVER_CALIBRATED,
        share_means=SHARE_MEANS,
        what_would_calibrate_it=WHAT_WOULD_CALIBRATE_IT,
        attestation_provenance=index.provenance,
        attestation_licence=index.dataset_version.licence,
        attestation_attribution=index.dataset_version.attribution,
        refusal=reason,
    )
