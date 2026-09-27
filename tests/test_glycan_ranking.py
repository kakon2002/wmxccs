"""The ranker: what it orders, what it refuses to order, and what its number does not mean.

Written against the real corpus wherever the figure is the point. Several of these tests exist
because a reviewer measured the first version of this design and found it wrong; where that is
so, the comment says what the wrong version did, because the shape is more useful than the fix.
"""

from __future__ import annotations

import ast
import inspect
import json
import math
import pathlib

import pytest

from wmxglycan import ranking
from wmxglycan.attestation import default_attestation_index
from wmxglycan.ccs_evidence import (
    CCSEvidence,
    CCSEvidenceState,
    CCSReference,
    EvidenceLevel,
    HeldValues,
)
from wmxglycan.composition import parse_composition
from wmxglycan.enumeration import Enumerator, is_order_constraint, _nodes_matching
from wmxglycan.features import extract_structure
from wmxglycan.glycan_graph import GlycanGraph
from wmxglycan.ranking import (
    ONE_BAND,
    PRIOR_PSEUDOCOUNT,
    PRIORS,
    SHARE_MEANS,
    Band,
    Calibration,
    Coverage,
    CoverageState,
    Decision,
    IndistinguishableClass,
    class_key_for,
    rank,
)
from wmxglycan.reuse import ReuseStatus

G2F = "Hex5HexNAc4Fuc1"
MAN5 = "Hex5HexNAc2"
G0F = "Hex3HexNAc4Fuc1"
G1F = "Hex4HexNAc4Fuc1"
ALL_FOUR = (G0F, G1F, G2F, MAN5)

# Measured in this repository against glycowork 1.10.0. Pinned so a drift on either side of the
# match - the enumerator's or the corpus's - fails rather than quietly re-ranking.
CLASSES = {G0F: 9, G1F: 32, G2F: 61, MAN5: 6}
BANDS = {G0F: 2, G1F: 3, G2F: 3, MAN5: 1}
# composition -> {attested structures: how many classes}
POOLED = {
    G0F: {1: 6, 0: 3},
    G1F: {2: 2, 1: 12, 0: 18},
    G2F: {4: 1, 1: 13, 0: 47},
    MAN5: {1: 6},
}
TIED = {G0F: 2, G1F: 51, G2F: 152, MAN5: 0}
LARGEST_TIE = {G0F: 2, G1F: 4, G2F: 10, MAN5: 1}


@pytest.fixture(scope="module")
def enumerator():
    return Enumerator()


@pytest.fixture(scope="module")
def index():
    return default_attestation_index()


@pytest.fixture(scope="module")
def ranked(enumerator, index):
    return {
        composition: rank(
            enumerator.enumerate(composition), index=index, enumerator=enumerator
        )
        for composition in ALL_FOUR
    }


# --- the unit of ranking is the indistinguishable class ----------------------------------------


@pytest.mark.parametrize("composition", ALL_FOUR)
def test_the_classes_are_derived_and_the_count_is_pinned(ranked, composition):
    result = ranked[composition]
    assert result.classes_total == CLASSES[composition]
    assert sum(len(one) for one in result.classes.values()) == result.candidates_total


def test_the_grouping_is_neither_empty_nor_saturated(ranked):
    # BOTH DIRECTIONS OF VACUITY. An empty grouping - every candidate its own class - would let
    # the platform order 152 of 167 G2F candidates it cannot distinguish, which is the failure
    # this project is built against. A saturated grouping - all 167 in one class, which is what
    # happens if a candidate is featurised without its structure and every column comes back
    # None - would make the isomer guard pass for the wrong reason.
    result = ranked[G2F]
    assert result.classes_total == 61
    assert 1 < result.classes_total < result.candidates_total
    assert result.tied_candidates == TIED[G2F] == 152
    assert result.largest_tie == LARGEST_TIE[G2F] == 10
    multi = [one for one in result.classes.values() if one.is_a_tie]
    assert len(multi) == 46


@pytest.mark.parametrize("composition", ALL_FOUR)
def test_tie_figures_are_pinned(ranked, composition):
    assert ranked[composition].tied_candidates == TIED[composition]
    assert ranked[composition].largest_tie == LARGEST_TIE[composition]


def test_a_class_holds_a_frozenset_so_an_order_within_it_is_unrepresentable(ranked):
    for one in ranked[G2F].classes.values():
        assert isinstance(one.members, frozenset)


def test_a_band_holds_a_frozenset_of_classes(ranked):
    for band in ranked[G2F].bands:
        assert isinstance(band.classes, frozenset)


def test_every_member_of_a_class_really_has_the_same_structure_row(ranked, enumerator):
    result = ranked[G2F]
    for one in result.classes.values():
        rows = {
            tuple(sorted(extract_structure(result.candidates[key].iupac_condensed).items()))
            for key in one.members
        }
        assert len(rows) == 1, f"class {one.class_id} holds candidates with different rows"


def test_candidates_in_one_class_are_genuinely_different_structures(ranked):
    """Indistinguishable is not identical: the members are different molecules.

    THIS ASSERTED NOTHING AND IS KEPT AS A LESSON. It read
    `len(one.members) == len(set(one.members))` where `members` is a frozenset, so
    `set(frozenset)` has the same length by construction and the assertion held for any code at
    all. What the property actually needs is that the members differ as STRUCTURES while
    agreeing as feature rows - which is the whole point of a class - and that they are drawn
    from the candidate set rather than invented.
    """
    result = ranked[G2F]
    checked = 0
    for one in result.classes.values():
        if not one.is_a_tie:
            continue
        checked += 1
        structures = {result.candidates[key].iupac_condensed for key in one.members}
        assert len(structures) == len(one.members), (
            f"class {one.class_id} holds two candidates with the same structure string, so the"
            " enumerator failed to deduplicate and the class is measuring the wrong thing"
        )
        rows = {
            tuple(sorted(extract_structure(result.candidates[key].iupac_condensed).items()))
            for key in one.members
        }
        assert len(rows) == 1, "the members must agree on every feature column"
    assert checked == 46, "the 46 multi-member classes must all be examined"


# --- the class key must not depend on a NaN identity accident ----------------------------------


def test_the_class_key_contains_no_nan(ranked, enumerator):
    for candidate in enumerator.enumerate(G2F).candidates[:40]:
        for _name, value in class_key_for(candidate):
            assert not (isinstance(value, float) and math.isnan(value))


def test_the_class_key_survives_a_json_round_trip(enumerator):
    # FeatureVector.as_row() maps an absent value to the math.nan singleton, and
    # (1.0, math.nan) == (1.0, math.nan) holds only because tuple comparison short-circuits on
    # identity. float("nan"), a numpy float64 and a JSON round trip each break it, at which
    # point every class becomes a singleton and nothing fails. The key is built from the
    # None-bearing mapping instead, so this round trip preserves grouping.
    candidate = enumerator.enumerate(G2F).candidates[0]
    key = class_key_for(candidate)
    again = tuple(tuple(pair) for pair in json.loads(json.dumps([list(p) for p in key])))
    assert again == key


def test_the_identity_accident_is_real_so_the_choice_of_key_matters():
    # Pins the reason rather than the fix, so that someone "simplifying" the key learns why.
    assert (1.0, math.nan) == (1.0, math.nan)
    assert (1.0, float("nan")) != (1.0, float("nan"))
    assert (1.0, None) == (1.0, None)


def test_a_key_carrying_a_nan_is_refused(monkeypatch):
    monkeypatch.setattr(ranking, "extract_structure", lambda _text: {"bonds": float("nan")})
    candidate = Enumerator().enumerate(MAN5).candidates[0]
    with pytest.raises(ValueError, match="contains a NaN"):
        class_key_for(candidate)


# --- attestation is pooled to the class, and deduplicated ---------------------------------------


@pytest.mark.parametrize("composition", ALL_FOUR)
def test_the_pooled_attestation_histogram_is_pinned(ranked, composition):
    # The top G2F class holds 4 distinct attested structures, which is the one genuine
    # discriminating signal in a set of 167. Before deduplication and pooling, the top of this
    # ranking was three candidates SugarBase happens to spell twice.
    result = ranked[composition]
    histogram: dict[int, int] = {}
    for one in result.classes.values():
        histogram[one.attested_structures] = histogram.get(one.attested_structures, 0) + 1
    assert histogram == POOLED[composition]


def test_pooled_attestation_equals_the_distinct_attested_members(ranked, index):
    for one in ranked[G2F].classes.values():
        assert one.attested_structures == sum(1 for key in one.members if index.attests(key))


def test_the_row_count_is_reported_beside_the_structure_count(ranked):
    # So a future release that dedupes or duplicates differently is visible rather than
    # silently re-ranking. Rows may exceed structures; they may never decide a band.
    top = max(ranked[G2F].classes.values(), key=lambda one: one.attested_structures)
    assert top.attested_structures == 4
    assert top.reference_rows >= top.attested_structures


# --- bands, and that no band splits a class -----------------------------------------------------


@pytest.mark.parametrize("composition", ALL_FOUR)
def test_the_band_count_is_pinned(ranked, composition):
    assert len(ranked[composition].bands) == BANDS[composition]


def test_bands_are_ordered_by_evidence_and_cover_every_class(ranked):
    for composition in ALL_FOUR:
        result = ranked[composition]
        counts = [band.attested_structures for band in result.bands]
        assert counts == sorted(counts, reverse=True)
        assert sum(len(band.classes) for band in result.bands) == result.classes_total
        assert sum(band.candidates for band in result.bands) == result.candidates_total


def test_no_class_is_split_across_bands(ranked):
    for composition in ALL_FOUR:
        seen: set[str] = set()
        for band in ranked[composition].bands:
            assert not (seen & band.classes)
            seen |= band.classes


def test_every_candidate_in_one_class_shares_one_band(ranked):
    # The product rule, asserted as an equal band index rather than as a populated twin list.
    # A separate field a consumer must join against the bands is a guard that relies on someone
    # reading it, and this repository's recurring failure class is exactly that.
    #
    # THIS TEST WAS A TAUTOLOGY AND IS KEPT AS A LESSON. It read:
    #     ranks = {band_of[one.class_id]}
    #     assert len(ranks) == 1
    # A one-element set literal has length one whatever the code does; only a KeyError could
    # have reddened it. The property it names is that every CANDIDATE in a class sits in the
    # same band, so the band must be looked up per MEMBER and the set built across members.
    result = ranked[G2F]
    band_of = {
        class_id: band.rank for band in result.bands for class_id in band.classes
    }
    member_to_class = {
        key: one.class_id for one in result.classes.values() for key in one.members
    }
    assert len(member_to_class) == result.candidates_total
    for one in result.classes.values():
        ranks = {band_of[member_to_class[key]] for key in one.members}
        assert len(ranks) == 1, f"class {one.class_id} straddles bands {sorted(ranks)}"
    # And the same property stated over the whole set, so a class of one cannot satisfy it
    # trivially: every candidate's band is the band of its own class, for all 167.
    assert {band_of[member_to_class[key]] for key in member_to_class} == {
        band.rank for band in result.bands
    }


# --- the ranker never consumes enumerator order -------------------------------------------------


def test_the_result_does_not_depend_on_the_order_candidates_arrive_in(enumerator, index):
    import dataclasses

    result = enumerator.enumerate(G1F)
    forward = rank(result, index=index, enumerator=enumerator)
    backward = rank(
        dataclasses.replace(result, candidates=tuple(reversed(result.candidates))),
        index=index,
        enumerator=enumerator,
    )
    assert forward.classes_total == backward.classes_total
    assert {b.rank: (b.classes, b.attested_structures) for b in forward.bands} == {
        b.rank: (b.classes, b.attested_structures) for b in backward.bands
    }
    assert set(forward.classes) == set(backward.classes)


def test_a_class_id_is_derived_from_its_members(ranked):
    for one in ranked[G2F].classes.values():
        assert one.class_id == min(one.members)


# --- refuse on inevaluability, warn on weakness -------------------------------------------------


def test_man5_is_refused_because_nothing_separates_its_candidates(ranked):
    # THE STANDING REGRESSION CASE. Six candidates, six classes, every one attested by exactly
    # one distinct structure: one band, no discrimination. Under the first design's row counting
    # one of them would have won on the strength of a structure entered twice, once with
    # accession G83351GR and once with none - and Man5 is the most studied N-glycan there is.
    result = ranked[MAN5]
    assert not result.is_a_ranking
    assert result.refusal is not None
    assert "RANKING REFUSED" in result.refusal
    assert "arbitrary order presented as a ranking" in result.refusal
    assert result.candidates_total == 6
    assert len(result.bands) == 1


def test_the_refusal_trigger_is_the_band_structure_and_not_zero_attestation(ranked):
    # Man5 has attestation on every candidate, so a trigger keyed on "no candidate is attested"
    # would let it through as a ranking. The trigger is one band.
    result = ranked[MAN5]
    assert all(one.attested_structures == 1 for one in result.classes.values())
    assert result.refusal is not None


def test_a_set_whose_bands_separate_is_not_refused(ranked):
    # The positive control, so the refusal is not vacuously always on.
    for composition in (G0F, G1F, G2F):
        assert ranked[composition].is_a_ranking, composition
        assert ranked[composition].refusal is None


def test_an_empty_candidate_set_is_refused_with_the_enumerators_reason(enumerator, index):
    result = enumerator.enumerate("Hex2HexNAc2")
    assert not result.candidates
    ranked_empty = rank(result, index=index, enumerator=enumerator)
    assert not ranked_empty.is_a_ranking
    assert "no candidate" in ranked_empty.refusal
    assert ranked_empty.decision is Decision.IM_VALIDATION_REQUIRED


def test_the_one_band_wording_is_importable(ranked):
    assert "{composition}" in ONE_BAND


@pytest.mark.parametrize("composition", ALL_FOUR)
def test_every_result_reports_its_weaknesses(ranked, composition):
    assert ranked[composition].weaknesses


def test_the_weaknesses_name_the_thinness_of_the_discrimination(ranked):
    joined = " ".join(ranked[G2F].weaknesses)
    assert "4 distinct reference structure(s) out of 34" in joined
    assert "152 candidate(s) share an indistinguishable class" in joined


# --- coverage: the candidate set may not contain the answer --------------------------------------


@pytest.mark.parametrize("composition", ALL_FOUR)
def test_coverage_is_measured_and_says_the_answer_may_be_elsewhere(ranked, composition):
    coverage = ranked[composition].coverage
    assert coverage.state is CoverageState.MEASURED
    assert coverage.not_enumerated > 0
    assert coverage.truth_may_not_be_in_the_candidate_set
    assert coverage.attested_by_a_candidate + coverage.not_enumerated == coverage.reference_structures


def test_coverage_is_unevaluable_when_the_corpus_holds_nothing_for_the_composition(
    enumerator, index
):
    # And UNEVALUABLE must still say the answer may not be in the set: no reference structure of
    # a composition is not evidence that the enumerator is complete for it.
    # Hex4HexNAc4Fuc1NeuAc2 is enumerable - 160 candidates - and the corpus holds no
    # fully-resolved structure of it. An earlier version of this test used a composition the
    # enumerator refuses outright, so it passed on the empty-set path and never reached the
    # coverage path it names. That is the fixture-where-the-wrong-answer-is-the-right-one shape.
    novel = "Hex4HexNAc4Fuc1NeuAc2"
    result = rank(enumerator.enumerate(novel), index=index, enumerator=enumerator)
    assert result.candidates_total == 160, "the premise is that this composition IS enumerable"
    assert result.coverage.reference_structures == 0
    assert result.coverage.state is CoverageState.UNEVALUABLE
    assert result.coverage.truth_may_not_be_in_the_candidate_set
    assert "NOT evidence that they do" in result.coverage.summary()


# --- the number, and what stops it being a probability -------------------------------------------


def test_nothing_sums_to_one_over_the_candidate_set_alone(ranked):
    # Normalising over the candidates asserts the answer is among them, which is measured false
    # for a large and highly variable fraction of compositions. Two reserved masses make the shares honest:
    # the reference structures no candidate matches, and the catch-all for a structure nobody
    # proposed at all.
    for composition in (G0F, G1F, G2F):
        result = ranked[composition]
        over_classes = sum(
            one.evidence_share for one in result.classes.values() if one.evidence_share
        )
        assert over_classes < 1.0
        assert result.share_not_enumerated > 0
        assert result.share_not_proposed > 0
        assert over_classes + result.share_not_enumerated + result.share_not_proposed == (
            pytest.approx(1.0, abs=1e-9)
        )


def test_the_catch_all_hypothesis_is_never_zero(ranked):
    # THE REGRESSION GUARD FOR THE CLOSED-WORLD BUG. Reserving mass only for reference
    # structures the corpus knows about still closed the world whenever it knew of none that
    # were missed. The catch-all cannot be ruled out by any corpus, so it always holds mass.
    for composition in ALL_FOUR:
        result = ranked[composition]
        assert result.share_not_proposed is not None and result.share_not_proposed > 0


def test_a_composition_whose_every_reference_structure_is_enumerated_still_does_not_sum_to_one(
    enumerator, index
):
    # Hex6HexNAc3Fuc3 is the case that broke the first version: it has exactly ONE
    # fully-resolved reference structure, that structure IS among its candidates, so
    # not_enumerated was 0 and the class shares summed to exactly 1.0 - the response asserting
    # the answer was in the set on the strength of one deposited structure.
    result = rank(enumerator.enumerate("Hex6HexNAc3Fuc3"), index=index, enumerator=enumerator)
    assert result.coverage.reference_structures == 1
    assert result.coverage.not_enumerated == 0
    assert result.share_not_enumerated == 0.0
    over_classes = sum(
        one.evidence_share for one in result.classes.values() if one.evidence_share
    )
    assert over_classes < 1.0, "the shares must not close the world"
    assert result.share_not_proposed > 0
    assert over_classes + result.share_not_proposed == pytest.approx(1.0, abs=1e-9)


def test_completeness_is_never_reported_as_established(ranked, enumerator, index):
    # There is no COMPLETE member of SetCompleteness, and that is the design. A corpus records
    # what has been deposited, not what exists.
    from wmxglycan.ranking import SetCompleteness

    assert set(SetCompleteness) == {
        SetCompleteness.INCOMPLETE,
        SetCompleteness.NOT_ESTABLISHABLE,
    }
    for composition in ALL_FOUR:
        assert ranked[composition].coverage.truth_may_not_be_in_the_candidate_set is True
        assert ranked[composition].coverage.completeness is SetCompleteness.INCOMPLETE
    # And the nothing-missing case is NOT_ESTABLISHABLE rather than complete.
    nothing_missing = rank(
        enumerator.enumerate("Hex6HexNAc3Fuc3"), index=index, enumerator=enumerator
    )
    assert nothing_missing.coverage.completeness is SetCompleteness.NOT_ESTABLISHABLE
    assert nothing_missing.coverage.truth_may_not_be_in_the_candidate_set is True
    assert "NOTHING FOLLOWS FROM THAT" in nothing_missing.coverage.summary()


def test_every_published_prior_has_its_own_reserved_masses(ranked):
    # Only the policy prior's complement was published at first, so three of the four columns
    # did not sum to one and a reader taking the Haldane column got 0.804.
    for composition in (G0F, G1F, G2F):
        result = ranked[composition]
        assert set(result.mass_under_priors) == set(PRIORS)
        for name in PRIORS:
            masses = result.mass_under_priors[name]
            over_classes = sum(
                one.share_under_priors[name]
                for one in result.classes.values()
                if one.share_under_priors[name] is not None
            )
            if masses["not_enumerated"] is None:
                continue
            total = over_classes + masses["not_enumerated"] + masses["not_proposed"]
            assert total == pytest.approx(1.0, abs=1e-9), f"{composition} under {name}: {total}"


def test_the_policy_prior_is_named_in_the_response(ranked):
    # Identifying it by comparing floats against four columns is ambiguous whenever two of them
    # coincide, which happens as soon as the space is small.
    from wmxglycan.ranking import POLICY_PRIOR

    for composition in ALL_FOUR:
        result = ranked[composition]
        assert result.prior == POLICY_PRIOR
        assert result.prior_pseudocount == PRIOR_PSEUDOCOUNT
        assert result.prior in PRIORS


def test_an_undefined_share_is_published_as_none_and_never_as_zero(enumerator, index):
    # Haldane's prior with nothing attested: the denominator is zero and the share does not
    # exist. A 0.0 there would make one column read as a measured near-impossibility for every
    # hypothesis while the other three summed to one.
    result = rank(
        enumerator.enumerate("Hex4HexNAc4Fuc1NeuAc2"), index=index, enumerator=enumerator
    )
    assert result.coverage.reference_structures == 0
    for one in result.classes.values():
        assert one.share_under_priors["haldane_0"] is None
    assert result.mass_under_priors["haldane_0"]["not_enumerated"] is None


def test_the_share_on_structures_nobody_enumerated_is_a_headline(ranked):
    # Slightly below the figures before the catch-all hypothesis was added, because the space
    # grew by one: 22/28 of Man5's reference structures are unreachable by the enumerator.
    assert ranked[MAN5].share_not_enumerated == pytest.approx(0.7719, abs=1e-3)
    assert ranked[G2F].share_not_enumerated == pytest.approx(0.3009, abs=1e-3)


def test_the_share_is_published_under_all_four_standard_priors(ranked):
    for composition in (G0F, G1F, G2F):
        for one in ranked[composition].classes.values():
            assert set(one.share_under_priors) == set(PRIORS)
        for band in ranked[composition].bands:
            assert set(band.share_under_priors) == set(PRIORS)


def test_the_prior_is_a_lever_and_the_response_shows_it(ranked):
    # The refutation of the claim the first design made, that add-one smoothing was "forced".
    # If the four priors gave the same answer there would be nothing to declare; they do not.
    top = max(ranked[G2F].classes.values(), key=lambda one: one.attested_structures)
    values = set(round(value, 6) for value in top.share_under_priors.values())
    assert len(values) == 4, top.share_under_priors
    assert max(values) > 2 * min(values)


def test_the_policy_share_is_the_declared_prior(ranked):
    top = max(ranked[G2F].classes.values(), key=lambda one: one.attested_structures)
    assert PRIOR_PSEUDOCOUNT == 1.0
    assert top.evidence_share == pytest.approx(top.share_under_priors["uniform_1"])


def test_the_meaning_of_the_number_travels_with_the_result(ranked):
    # In the response, not only in LIMITATIONS, because the consumer reads the response.
    for composition in ALL_FOUR:
        result = ranked[composition]
        assert result.share_means is SHARE_MEANS
        assert "NOT A MOLECULE IN YOUR SAMPLE" in result.share_means
        assert "never been calibrated" in result.share_means
        assert result.calibration is Calibration.NEVER_CALIBRATED
        assert "independently known" in result.what_would_calibrate_it


def test_no_field_is_called_probability_or_confidence():
    source = pathlib.Path(ranking.__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    names = {
        target.target.id if isinstance(target, ast.AnnAssign) else ""
        for node in ast.walk(tree)
        if isinstance(node, ast.ClassDef)
        for target in node.body
        if isinstance(target, ast.AnnAssign) and isinstance(target.target, ast.Name)
    }
    for forbidden in ("probability", "confidence", "likelihood", "p_correct"):
        assert not any(forbidden in name for name in names), forbidden


def test_a_band_reports_both_its_per_class_share_and_its_total(ranked):
    # The first version gave the band a field called `evidence_share` - the same name the class
    # carries - and printed it beside the band's class and candidate counts, so a consumer read
    # a per-class number as the band's mass and understated band 3 by a factor of 47.
    result = ranked[G2F]
    for band in result.bands:
        assert band.share_per_class is not None
        assert band.band_total_share == pytest.approx(band.share_per_class * len(band.classes))
    assert not hasattr(result.bands[0], "evidence_share")
    totals = sum(band.band_total_share for band in result.bands)
    assert totals + result.share_not_enumerated + result.share_not_proposed == pytest.approx(
        1.0, abs=1e-9
    )
    assert "per class" in result.summary() and "for the band" in result.summary()


def test_the_unattested_mass_is_labelled_as_classes_and_not_candidates(ranked):
    # 47 unattested CLASSES hold 127 of the 167 candidates, so calling the pooled class mass a
    # candidate mass invited the reading that three quarters of the set held 42%.
    result = ranked[G2F]
    unattested_classes = [one for one in result.classes.values() if one.attested_structures == 0]
    unattested_candidates = sum(len(one) for one in unattested_classes)
    assert len(unattested_classes) == 47
    assert unattested_candidates == 127
    assert result.share_on_unattested_classes == pytest.approx(
        len(unattested_classes) * unattested_classes[0].evidence_share
    )
    assert "unattested CLASSES" in result.summary()
    assert not hasattr(result, "share_on_unattested")


def test_the_candidates_mapping_is_sorted_and_not_in_enumerator_order(ranked):
    # The only ordered container whose order the design does not assert, and it arrived in the
    # enumerator's depth-first order - beside a refusal saying the candidates are unordered.
    for composition in ALL_FOUR:
        keys = list(ranked[composition].candidates)
        assert keys == sorted(keys), composition


def test_a_truncated_candidate_set_says_so(enumerator, index):
    # A capped enumeration was published with a coverage line reading "N reference structures
    # are NOT among the candidates" as though the enumerator had considered and declined them,
    # when it had simply stopped.
    result = enumerator.enumerate("Hex6HexNAc5Fuc1")
    assert result.capped, "this composition must exceed the tree budget or the test proves nothing"
    ranked_capped = rank(result, index=index, enumerator=enumerator)
    assert ranked_capped.enumerator_capped is True
    assert any("TRUNCATED" in weakness for weakness in ranked_capped.weaknesses)


def test_an_uncapped_candidate_set_does_not_claim_truncation(ranked):
    for composition in ALL_FOUR:
        assert ranked[composition].enumerator_capped is False
        assert not any("TRUNCATED" in w for w in ranked[composition].weaknesses)


def test_a_single_candidate_is_not_refused_with_the_one_band_wording(enumerator, index):
    # Man3GlcNAc2 gives one candidate. The one-band refusal said "nothing separates any of them
    # from any other" about a set of one, which reads as a failure where there is none.
    result = rank(enumerator.enumerate("Hex3HexNAc2"), index=index, enumerator=enumerator)
    assert result.candidates_total == 1
    assert result.refusal is not None
    assert "NOT A RANKING" in result.refusal
    assert "exactly one candidate" in result.refusal
    assert "nothing separates any of them" not in result.refusal
    assert "none is being withheld" in result.refusal


def test_a_rule_violation_among_the_candidates_refuses_rather_than_being_reported_and_ignored(
    enumerator, index, monkeypatch
):
    # The first version measured rules_violated, published it, and did nothing with it: a
    # candidate set breaking curated rules was banded and returned with no refusal and no
    # weakness. A non-zero count means the enumerator and the rule check disagree.
    monkeypatch.setattr(
        Enumerator, "broken_rules", lambda _self, _graph: ("MGAT9: a fabricated violation",)
    )
    result = rank(enumerator.enumerate(G2F), index=index, enumerator=Enumerator())
    assert result.rule_accounting.rules_violated == 167
    assert result.refusal is not None
    assert "curated rule violation" in result.refusal
    assert not result.is_a_ranking


def test_a_composition_the_enumerator_declines_reports_the_corpus_it_actually_has(
    enumerator, index
):
    # The refused path hardcoded reference_structures=0, so a composition the enumerator
    # declines was reported as one the corpus knows nothing about. It holds 15 structures of
    # this one, and the claim was attached to a refusal, where a reader takes it as the reason.
    result = rank(enumerator.enumerate("Hex2HexNAc2Fuc2"), index=index, enumerator=enumerator)
    assert result.candidates_total == 0
    assert result.coverage.reference_structures == 15
    assert result.coverage.not_enumerated == 15
    assert result.coverage.state is CoverageState.MEASURED
    assert "15 are NOT" in result.coverage.summary()


def test_a_class_share_is_not_divided_among_its_members(ranked):
    # Splitting a class's share across its members would dilute evidence across cells the
    # platform cannot separate, which is the dilution pooling exists to stop.
    result = ranked[G2F]
    by_count = {one.attested_structures: one.evidence_share for one in result.classes.values()}
    for one in result.classes.values():
        assert one.evidence_share == by_count[one.attested_structures]


# --- the rules: three quantities, named apart ---------------------------------------------------


@pytest.mark.parametrize("composition", ALL_FOUR)
def test_no_candidate_violates_a_rule_and_it_is_verified_by_running_the_check(ranked, composition):
    accounting = ranked[composition].rule_accounting
    assert accounting.rules_violated == 0
    assert accounting.rules_violated_verified_by_running_the_check
    assert accounting.rules_in_scheme == 15


def test_the_applicable_rule_count_varies_and_is_not_the_advertised_constant(ranked):
    # The first design reported one scalar and called it constant. Measured, the number of rules
    # bearing on a candidate varies, and a field documented as constant is a field a reader is
    # told not to check.
    assert ranked[G2F].rule_accounting.applicable_varies
    assert len(set(ranked[G2F].rule_accounting.rules_applicable.values())) > 1
    assert "VARIES" in ranked[G2F].rule_accounting.summary()


def test_no_curated_rule_bears_on_any_man5_candidate(ranked):
    # And so "15 rules satisfied" would have been false for all six of them.
    assert set(ranked[MAN5].rule_accounting.rules_applicable.values()) == {0}
    assert not ranked[MAN5].rule_accounting.applicable_varies


def test_the_ordering_caveats_are_counted_and_kept_out_of_the_number(ranked):
    result = ranked[G2F]
    caveats = result.rule_accounting.ordering_caveats
    assert sum(caveats.values()) == 6
    with_caveat = {key for key, count in caveats.items() if count}
    assert len(with_caveat) == 6
    assert with_caveat <= set(result.candidates)

    # THE STATED PROPERTY, NOW ACTUALLY ASSERTED. The old version built a `bands` dict, checked
    # it was truthy, and threw it away - so "the caveat cannot have moved a candidate's band"
    # was a comment rather than a test. A candidate's band must be a function of its class's
    # attestation and of nothing else, so two candidates with the same attestation must share a
    # band whether or not one of them carries a caveat.
    member_to_class = {
        key: one.class_id for one in result.classes.values() for key in one.members
    }
    band_of = {cid: band.rank for band in result.bands for cid in band.classes}
    attested_of = {
        one.class_id: one.attested_structures for one in result.classes.values()
    }
    by_attestation: dict[int, set[int]] = {}
    for key in result.candidates:
        cls = member_to_class[key]
        by_attestation.setdefault(attested_of[cls], set()).add(band_of[cls])
    for count, ranks in by_attestation.items():
        assert len(ranks) == 1, f"attestation {count} spans bands {sorted(ranks)}"
    # And a caveat-carrying candidate is not concentrated in the top band, which is what a
    # silent penalty or bonus would look like.
    caveat_bands = {band_of[member_to_class[key]] for key in with_caveat}
    assert caveat_bands, "the six caveat-carrying candidates must land in some band"


def test_every_ordering_caveat_is_true_of_the_candidate_that_carries_it(ranked, enumerator):
    # THE PORTED DEFECT THIS FIXES. Before the fix, 114 of 167 G2F candidates carried an
    # ordering note and only 18 contained the bisecting GlcNAc the note names as its context,
    # so 108 published a sentence carrying an enzyme and a Schachter 1986 citation that was not
    # true of them. Carrying rationale through unchanged would have turned one enumerator bug
    # into 108 false cited claims per query.
    order_contexts = [
        contexts for rule, _product, contexts in enumerator._rules if is_order_constraint(rule)
    ]
    assert order_contexts, "there must be order-constrained rules or this proves nothing"
    for key, candidate in ranked[G2F].candidates.items():
        carries = [reason for reason in candidate.rationale if reason.kind == "ordering"]
        if not carries:
            continue
        graph = GlycanGraph.from_iupac_condensed(candidate.iupac_condensed)
        assert any(
            _nodes_matching(graph, context) for contexts in order_contexts for context in contexts
        ), f"{key} publishes an ordering note without the context it names"


def test_the_rationale_and_its_citation_reach_the_ranked_candidate(ranked):
    result = ranked[G2F]
    cited = [
        reason
        for candidate in result.candidates.values()
        for reason in candidate.rationale
        if reason.reference
    ]
    assert cited
    for reason in cited[:10]:
        assert reason.enzymes
        assert reason.reference.strip()


# --- CCS evidence is a property of the set, never of a candidate ---------------------------------


def _reference(**overrides) -> CCSReference:
    fields = dict(
        ccs=700.0,
        uncertainty=3.0,
        uncertainty_type="two_sd",
        adduct="[M+H]+",
        charge=1,
        polarity="positive",
        ims_type="TWIMS",
        drift_gas="N2",
        calibrant="dextran",
        source="a synthetic fixture",
        doi="10.0000/fixture",
        source_locator="Table S1",
        reuse_status=ReuseStatus.SYNTHETIC_FIXTURE,
    )
    fields.update(overrides)
    return CCSReference(**fields)


def test_with_no_evidence_source_the_state_is_not_consulted_and_never_an_absence(ranked):
    for composition in ALL_FOUR:
        evidence = ranked[composition].ccs_evidence
        assert evidence.state is CCSEvidenceState.NOT_CONSULTED
        assert "nothing looked" in evidence.summary()
        assert evidence.state is not CCSEvidenceState.NONE_IN_CORPUS_SEARCHED


def test_a_lookup_that_raises_is_a_failure_and_never_an_absence(enumerator, index):
    class Broken:
        def evidence_for(self, composition, adduct, charge):
            raise RuntimeError("the adapter is misconfigured")

    result = rank(
        enumerator.enumerate(MAN5),
        index=index,
        enumerator=enumerator,
        ccs=Broken(),
        adduct="[M-H]-",
        charge=-1,
    )
    evidence = result.ccs_evidence
    assert evidence.state is CCSEvidenceState.LOOKUP_FAILED
    assert "RuntimeError" in evidence.failure
    assert evidence.key_attempted == "Hex5HexNAc2 [M-H]- -1"
    assert "not an absence" in evidence.summary()


def test_composition_level_evidence_cannot_be_attached_to_one_candidate(enumerator, index):
    class Composition_level:
        def evidence_for(self, composition, adduct, charge):
            return CCSEvidence(
                state=CCSEvidenceState.MEASURED_REFERENCE,
                level=EvidenceLevel.COMPOSITION,
                composition=composition.canonical,
                adduct=adduct,
                charge=charge,
                reference=_reference(),
            )

    result = rank(
        enumerator.enumerate(MAN5),
        index=index,
        enumerator=enumerator,
        ccs=Composition_level(),
        adduct="[M+H]+",
        charge=1,
    )
    assert result.ccs_evidence.state is CCSEvidenceState.MEASURED_REFERENCE
    assert not result.ccs_evidence.discriminates_between_candidates
    a_candidate = next(iter(result.candidates))
    with pytest.raises(ValueError, match="cannot be attached to an individual candidate"):
        result.ccs_evidence.for_candidate(a_candidate)


def test_a_measured_reference_shared_by_every_candidate_does_not_make_the_decision_softer(
    enumerator, index
):
    class Composition_level:
        def evidence_for(self, composition, adduct, charge):
            return CCSEvidence(
                state=CCSEvidenceState.MEASURED_REFERENCE,
                level=EvidenceLevel.COMPOSITION,
                composition=composition.canonical,
                adduct=adduct,
                charge=charge,
                reference=_reference(),
            )

    result = rank(
        enumerator.enumerate(G2F),
        index=index,
        enumerator=enumerator,
        ccs=Composition_level(),
        adduct="[M+H]+",
        charge=1,
    )
    assert result.decision is Decision.IM_VALIDATION_REQUIRED


# --- the decision field, and AI_ONLY proved unreachable rather than merely never seen -------------


@pytest.mark.parametrize("composition", ALL_FOUR)
def test_every_real_composition_lands_on_validation_required(ranked, composition):
    assert ranked[composition].decision is Decision.IM_VALIDATION_REQUIRED


def test_the_decision_rules_are_published_with_what_they_apply_to(ranked):
    rules = ranked[G2F].decision_rules
    assert len(rules) == 5
    for rule in rules:
        assert rule.applies_to.strip()
        assert rule.because.strip()
        assert rule.decision_if_it_fires is Decision.IM_VALIDATION_REQUIRED
    assert any(rule.fires for rule in rules)


def test_the_validated_model_gate_is_one_named_predicate():
    assert ranking.VALIDATED_MODEL is None
    assert ranking.a_validated_model_exists() is False


def test_ai_only_is_unreachable_across_the_whole_decision_input_space():
    # "AI_ONLY never appears" is equally satisfied by dead code, a misspelled comparison or an
    # enum member nothing references, so the sweep is only half the proof; the positive control
    # below is the other half.
    coverages = (
        Coverage(composition="X", state=CoverageState.MEASURED, reference_structures=3,
                 attested_by_a_candidate=3, not_enumerated=0),
        Coverage(composition="X", state=CoverageState.MEASURED, reference_structures=3,
                 attested_by_a_candidate=1, not_enumerated=2),
        Coverage(composition="X", state=CoverageState.UNEVALUABLE),
    )
    one = IndistinguishableClass(class_id="a", members=frozenset({"a"}), attested_structures=1)
    two = IndistinguishableClass(class_id="b", members=frozenset({"b", "c"}), attested_structures=1)
    class_sets = ({}, {"a": one}, {"b": two}, {"a": one, "b": two})
    # EVERY STATE AT EVERY LEVEL IT CAN LEGALLY CARRY, and the flag only where the validator
    # permits it - which is STRUCTURE level alone, since the validator now refuses a
    # composition-level finding that claims to discriminate.
    evidences = []
    for state in CCSEvidenceState:
        for level in EvidenceLevel:
            for discriminates in (False, True):
                if discriminates and level is not EvidenceLevel.STRUCTURE:
                    continue
                extra: dict = {}
                if state is CCSEvidenceState.MEASURED_REFERENCE:
                    extra["reference"] = _reference(
                        measured_on_canonical_key="a" if level is EvidenceLevel.STRUCTURE else None
                    )
                elif state is CCSEvidenceState.HELD_NOT_RELEASABLE:
                    extra["held"] = HeldValues(
                        records=1, blockers=("x",), what_would_release_them="y"
                    )
                elif state is CCSEvidenceState.NONE_IN_CORPUS_SEARCHED:
                    extra["corpus_searched"] = "c"
                    # Above zero: a denominator of nothing is not a denominator, and the
                    # validator now refuses it.
                    extra["records_consulted"] = 7
                elif state is CCSEvidenceState.LOOKUP_FAILED:
                    extra["failure"] = "boom"
                evidences.append(
                    CCSEvidence(
                        state=state,
                        level=level,
                        discriminates_between_candidates=discriminates,
                        **extra,
                    )
                )
    assert len(evidences) >= 15, "the sweep must cover a real space, not an empty one"

    # WHAT THIS SWEEP ESTABLISHES, AND WHAT IT CANNOT.
    #
    # It establishes that no COMBINATION OF INPUTS a caller can construct yields AI_ONLY. That
    # is the threat model worth guarding: a caller cannot talk the platform into AI_ONLY by
    # arranging its evidence, its classes or its coverage.
    #
    # It cannot establish that the AI_ONLY branch is live code, and an earlier version of this
    # test pretended otherwise. Two of the five decision rules fire on EVERY input there is -
    # "no validated model exists" and "the candidate set may not contain the answer" - so no
    # input reaches the code past the rule check, and `assert AI_ONLY not in seen` would have
    # passed just as well with the branch deleted. An attempt to fix that by holding the
    # coverage gate open failed for the same reason one level down: the model gate is itself one
    # of the rules, so the sweep still could not get there.
    #
    # So the division of labour is explicit. THE SWEEP asserts that every input lands on
    # REQUIRED and that the two always-firing rules are the reason. THE POSITIVE CONTROL below
    # patches both gates and reaches AI_ONLY, which is what shows the branch exists at all.
    seen = set()
    always_fired: set[str] = set()
    first = True
    combinations = 0
    for classes in class_sets:
        for coverage in coverages:
            for evidence in evidences:
                for refused in (False, True):
                    combinations += 1
                    decision, rules = ranking._decide(
                        classes=classes, bands=(), coverage=coverage,
                        evidence=evidence, refused=refused,
                    )
                    seen.add(decision)
                    fired = {rule.name for rule in rules if rule.fires}
                    assert fired, "some rule must fire, or the sweep says nothing"
                    always_fired = fired if first else (always_fired & fired)
                    first = False

    assert combinations >= 200, f"the sweep must be a real space, not a token one: {combinations}"
    assert seen == {Decision.IM_VALIDATION_REQUIRED}, seen
    assert Decision.AI_ONLY not in seen
    # The two rules that fire on every single input, named. If either stopped being universal
    # this fails, which is the signal that the reachability of AI_ONLY has changed.
    assert always_fired == {
        "no validated model exists",
        "the candidate set may not contain the answer",
    }, sorted(always_fired)


def test_the_decision_field_has_exactly_one_reachable_value_today(ranked):
    """And it is a fact about the evidence rather than a tuning choice.

    The owner asked for three values and expected "almost everything" to land on REQUIRED. The
    honest figure is stronger and is stated rather than softened: today it is EVERYTHING, for
    two independent reasons that hold for every possible input - no validated model exists, and
    the completeness of a candidate set is not establishable from a corpus of depositions.

    Nothing was tuned to make some cases look better. What would move it is a validated model
    and a cross section that resolves to one structure, and this repository has neither.
    """
    assert {result.decision for result in ranked.values()} == {Decision.IM_VALIDATION_REQUIRED}
    universal = [
        rule.name
        for rule in ranked[G2F].decision_rules
        if rule.fires and "validated model" in rule.name or "may not contain" in rule.name
    ]
    assert len(universal) == 2


def test_ai_only_needs_two_independent_gates_opened_and_is_live_code(monkeypatch):
    """THE POSITIVE CONTROL, and it now has to open TWO gates.

    AI_ONLY was unreachable for one reason when this was written - no validated model exists.
    Fixing the coverage overclaim added a second, independent one: completeness can never be
    established from a corpus that records depositions rather than existence, so the rule "the
    candidate set may not contain the answer" fires on every input there is.

    Both are patched here, because a positive control that opens one gate and finds the branch
    still shut proves nothing about whether the branch is live. Opening both and reaching
    AI_ONLY is what shows this is real code rather than an enum member nothing references.
    """
    monkeypatch.setattr(ranking, "VALIDATED_MODEL", object())
    monkeypatch.setattr(
        ranking.Coverage, "truth_may_not_be_in_the_candidate_set", property(lambda _self: False)
    )
    assert ranking.a_validated_model_exists() is True
    one = IndistinguishableClass(class_id="a", members=frozenset({"a"}), attested_structures=1)
    decision, rules = ranking._decide(
        classes={"a": one},
        bands=(),
        coverage=Coverage(
            composition="X", state=CoverageState.MEASURED, reference_structures=1,
            attested_by_a_candidate=1, not_enumerated=0,
        ),
        evidence=CCSEvidence(
            state=CCSEvidenceState.MEASURED_REFERENCE,
            level=EvidenceLevel.STRUCTURE,
            reference=_reference(measured_on_canonical_key="a"),
            discriminates_between_candidates=True,
        ),
        refused=False,
    )
    assert decision is Decision.AI_ONLY
    assert not any(rule.fires for rule in rules)


def test_opening_only_the_model_gate_leaves_ai_only_shut(monkeypatch):
    # Which is the point of the test above needing two. If this passed by opening one gate, the
    # second would not be doing anything.
    monkeypatch.setattr(ranking, "VALIDATED_MODEL", object())
    one = IndistinguishableClass(class_id="a", members=frozenset({"a"}), attested_structures=1)
    decision, rules = ranking._decide(
        classes={"a": one},
        bands=(),
        coverage=Coverage(
            composition="X", state=CoverageState.MEASURED, reference_structures=1,
            attested_by_a_candidate=1, not_enumerated=0,
        ),
        evidence=CCSEvidence(
            state=CCSEvidenceState.MEASURED_REFERENCE,
            level=EvidenceLevel.STRUCTURE,
            reference=_reference(measured_on_canonical_key="a"),
            discriminates_between_candidates=True,
        ),
        refused=False,
    )
    assert decision is Decision.IM_VALIDATION_REQUIRED
    fired = [rule.name for rule in rules if rule.fires]
    assert fired == ["the candidate set may not contain the answer"]


def test_the_ai_only_literal_appears_in_exactly_one_returning_branch():
    """Exactly one place DECIDES AI_ONLY, so there is one branch a test can reach and patch.

    SHARPENED 27 September 2026. This walked every node under every `return` and counted any
    mention of AI_ONLY, which caught `decision_reachability()` - a function that returns a
    DESCRIPTION of AI_ONLY's reachability rather than the decision itself, where the name appears
    as a keyword argument to a constructor. Broadening the allowance to two would have weakened the
    guard permanently, so what it counts is narrowed to what it was always about instead: AI_ONLY
    appearing as the returned value, or as an element of a returned tuple, rather than anywhere
    inside a nested call. `return Decision.AI_ONLY, rules` counts; `return (Thing(decision=...),)`
    does not, and a second real decision site still fails this.
    """
    source = pathlib.Path(ranking.__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)

    def decided_here(node: ast.Return) -> bool:
        value = node.value
        if value is None:
            return False
        returned = value.elts if isinstance(value, ast.Tuple) else [value]
        return any(
            isinstance(one, ast.Attribute) and one.attr == "AI_ONLY" for one in returned
        )

    returns = [node for node in ast.walk(tree) if isinstance(node, ast.Return) and decided_here(node)]
    assert len(returns) == 1, "AI_ONLY must be DECIDED in exactly one place"

    # The floor: the narrowing above must not have made this count nothing at all, which would
    # leave the assertion passing on a module that never mentions AI_ONLY.
    assert "AI_ONLY" in source
    decide = inspect.getsource(ranking._decide)
    assert "a_validated_model_exists()" in decide


def test_removing_the_model_gate_alone_does_not_open_ai_only(monkeypatch):
    # Belt and braces: with a validated model but more than one possible structure, the answer
    # must still be REQUIRED. Otherwise the gate would be the only thing holding the field shut.
    monkeypatch.setattr(ranking, "VALIDATED_MODEL", object())
    one = IndistinguishableClass(class_id="a", members=frozenset({"a"}), attested_structures=1)
    two = IndistinguishableClass(class_id="b", members=frozenset({"b"}), attested_structures=0)
    decision, _rules = ranking._decide(
        classes={"a": one, "b": two},
        bands=(),
        coverage=Coverage(composition="X", state=CoverageState.MEASURED, reference_structures=1,
                          attested_by_a_candidate=1, not_enumerated=0),
        evidence=CCSEvidence(
            state=CCSEvidenceState.MEASURED_REFERENCE,
            level=EvidenceLevel.STRUCTURE,
            reference=_reference(measured_on_canonical_key="a"),
            discriminates_between_candidates=True,
        ),
        refused=False,
    )
    assert decision is Decision.IM_VALIDATION_REQUIRED


# --- provenance ------------------------------------------------------------------------------------


def test_the_result_carries_the_live_dataset_version(ranked):
    from importlib.metadata import version

    for composition in ALL_FOUR:
        provenance = ranked[composition].attestation_provenance
        assert version("glycowork") in provenance
        assert "SugarBase v12" in provenance


def test_the_candidate_total_is_what_is_returned(ranked):
    for composition in ALL_FOUR:
        result = ranked[composition]
        assert result.candidates_total == len(result.candidates)
        assert result.candidates_total == sum(len(one) for one in result.classes.values())


def test_the_enumerators_own_accounting_travels(ranked):
    assert ranked[G2F].enumerator_built == 1012
    assert ranked[G2F].enumerator_rejected == 845


def test_the_rule_check_is_actually_run_and_not_assumed(enumerator, index, monkeypatch):
    # Without this, `rules_violated == 0` passes whether the check ran or not: zero is both the
    # right answer and what an unrun check returns. So the check is forced to report a violation
    # and the accounting must carry it. This is the vacuity that made the FIRST version of this
    # design claim a constant it had never measured.
    monkeypatch.setattr(
        Enumerator, "broken_rules", lambda _self, _graph: ("MGAT9: a fabricated violation",)
    )
    result = rank(enumerator.enumerate(MAN5), index=index, enumerator=Enumerator())
    assert result.rule_accounting.rules_violated == 6, (
        "the accounting must report what the check returns, not what the enumerator promised"
    )


def test_the_share_is_refused_rather_than_invented_when_there_is_no_hypothesis_space(
    enumerator, index
):
    # A composition the corpus knows nothing about: every class unattested, no not-enumerated
    # mass. The shares are then pure prior and carry no information, and the one-band refusal
    # is what stops them being read as a ranking.
    result = rank(enumerator.enumerate("Hex4HexNAc4Fuc1NeuAc2"), index=index, enumerator=enumerator)
    assert result.coverage.state is CoverageState.UNEVALUABLE
    assert all(one.attested_structures == 0 for one in result.classes.values())
    assert len(result.bands) == 1
    assert not result.is_a_ranking
    assert result.share_not_enumerated == 0.0


def test_a_composition_the_enumerator_refuses_is_a_different_refusal(enumerator, index):
    # Kept beside the one above so the two paths cannot be conflated again: no candidates at all
    # is not the same answer as candidates that nothing distinguishes.
    result = rank(enumerator.enumerate("Hex5HexNAc4NeuGc2"), index=index, enumerator=enumerator)
    assert result.candidates_total == 0
    assert result.bands == ()
    assert "no candidate" in result.refusal


# --- rendering. Barely tested until it crashed on a refused set. -------------------------------------


@pytest.mark.parametrize("composition", ALL_FOUR)
def test_the_summary_renders_for_every_composition_including_a_refused_one(ranked, composition):
    # THIS FOUND A CRASH. `summary()` raised TypeError on Hex5HexNAc2 because
    # share_on_unattested_classes was None - Man5 has no unattested class, so the honest value
    # is a measured 0.0 rather than "undefined" - and no test had rendered a refused set. A
    # method nothing calls is a method nothing protects.
    text = ranked[composition].summary()
    assert composition in text
    assert "calibration: never_calibrated" in text
    assert "prior: uniform_1" in text
    assert "MIT" in text and "Daniel Bojar" in text


def test_the_unattested_mass_is_zero_when_there_are_no_unattested_classes(ranked):
    # Distinguished from "undefined", which is what an unsatisfiable denominator gives.
    assert all(one.attested_structures == 1 for one in ranked[MAN5].classes.values())
    assert ranked[MAN5].share_on_unattested_classes == 0.0
    assert ranked[G2F].share_on_unattested_classes > 0


def test_a_refused_set_still_returns_its_candidates_and_its_evidence(ranked):
    # Refusing to ORDER is not refusing to ANSWER: the candidates, their attestation and the
    # coverage all travel, which is what makes the refusal usable rather than a dead end.
    result = ranked[MAN5]
    assert not result.is_a_ranking
    assert len(result.candidates) == 6
    assert len(result.classes) == 6
    assert result.coverage.reference_structures == 28
    assert result.attestation_licence == "MIT"


def test_the_decision_reads_the_evidence_level_even_when_the_validator_is_bypassed():
    """Defence in depth, and THE SWEEP FOUND IT UNTESTED.

    `_decide` checks `evidence.level is EvidenceLevel.STRUCTURE` as well as the
    `discriminates_between_candidates` flag, and the `CCSEvidence` validator forbids that
    combination outright - so no ordinary construction can tell the two versions of `_decide`
    apart, and the mutation that removes the level check SURVIVED. A guard present and
    unprotected is the failure class this project keeps meeting.

    So the input is built with `model_construct`, which skips every validator, exactly as the
    CCS core's suite does for a record that reached it without validation. That is the real
    threat model too: the validator guards the constructor, and a producer using
    `model_construct` or `model_copy(update=...)` reaches `_decide` without passing it.
    """
    smuggled = CCSEvidence.model_construct(
        state=CCSEvidenceState.MEASURED_REFERENCE,
        level=EvidenceLevel.COMPOSITION,
        reference=_reference(),
        discriminates_between_candidates=True,
    )
    # The bypass really did produce the illegal object the validator refuses.
    assert smuggled.level is EvidenceLevel.COMPOSITION
    assert smuggled.discriminates_between_candidates is True
    with pytest.raises(ValueError, match="cannot also claim to discriminate"):
        CCSEvidence(
            state=CCSEvidenceState.MEASURED_REFERENCE,
            level=EvidenceLevel.COMPOSITION,
            reference=_reference(),
            discriminates_between_candidates=True,
        )

    one = IndistinguishableClass(class_id="a", members=frozenset({"a"}), attested_structures=1)
    _decision, rules = ranking._decide(
        classes={"a": one},
        bands=(),
        coverage=Coverage(
            composition="X", state=CoverageState.MEASURED, reference_structures=1,
            attested_by_a_candidate=1, not_enumerated=0,
        ),
        evidence=smuggled,
        refused=False,
    )
    held = [r for r in rules if r.name == "no cross section is held for this structure"]
    assert len(held) == 1
    assert held[0].fires is True, (
        "a composition-level measurement is shared by every candidate, so the rule must still"
        " fire however the finding labels itself"
    )
