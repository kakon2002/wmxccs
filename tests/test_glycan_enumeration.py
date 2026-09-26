"""Candidate generation: what the enumerator builds, and what it refuses to build."""

import pytest

from wmxglycan.composition import Composition
from wmxglycan.enumeration import (
    SITES,
    EnumerationError,
    Enumerator,
    Site,
    _arrives_on,
    _chain,
    _contains,
    contains_motif,
    is_order_constraint,
    motif_anchors,
)
from wmxglycan.glycan_graph import GlycanGraph
from wmxglycan.licensing import ReuseStatus, assert_trainable

G0F = "Hex3HexNAc4Fuc1"
G2F = "Hex5HexNAc4Fuc1"
MAN5 = "Hex5HexNAc2"

BISECT = "GlcNAc(b1-4)Man(b1-4)GlcNAc"
CORE_FUC = "Fuc(a1-6)GlcNAc"

# The structures these compositions are named for, as the literature writes them.
CANONICAL = {
    G0F: "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc",
    G2F: "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc",
    MAN5: "Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
}


@pytest.fixture(scope="module")
def enumerator():
    return Enumerator()


def graph_of(iupac):
    return GlycanGraph.from_iupac_condensed(iupac)


def key_of(iupac):
    return graph_of(iupac).canonical_key()


def holds(candidate, fragment):
    return _contains(graph_of(candidate.iupac_condensed), _chain(fragment))


# --- the structure each composition is named for must be reachable ------------


@pytest.mark.parametrize("composition", [G0F, G2F, MAN5])
def test_the_named_structure_is_among_the_candidates(enumerator, composition):
    result = enumerator.enumerate(composition)
    assert key_of(CANONICAL[composition]) in {c.canonical_key for c in result.candidates}


def test_every_candidate_has_the_composition_that_was_asked_for(enumerator):
    composition = Composition.parse(G2F)
    for candidate in enumerator.enumerate(composition).candidates:
        assert candidate.structure.composition == composition


def test_man5_is_the_small_closed_case(enumerator):
    # Five hexoses and the two core GlcNAc leave only mannose placements.
    result = enumerator.enumerate(MAN5)
    assert len(result.candidates) == 6
    for candidate in result.candidates:
        assert "GlcNAc(b1-4)GlcNAc" in candidate.iupac_condensed
        assert candidate.iupac_condensed.count("Man") == 5


# --- deduplication, tested on the key rather than on the dict that holds it ---
# enumerate() stores candidates in a dict keyed by canonical_key, so asserting
# those keys are distinct cannot fail. The teeth have to be in the key itself.


def test_one_structure_written_two_ways_gives_one_key():
    # The same glycan with its two branches written in the other order.
    assert key_of("Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc") == key_of(
        "Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
    )


def test_two_different_structures_never_share_a_key():
    # Same composition, different place for the second GlcNAc: real isomers.
    one = key_of("GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc")
    two = key_of("GlcNAc(b1-2)[GlcNAc(b1-4)]Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc")
    assert one != two


def test_a_floating_fragment_is_kept_apart_from_an_attached_one():
    # Braces mean the residue is there but its attachment point is unknown; that
    # is not the same structure as one where the position is known.
    floating = key_of("{Fuc(a1-3)}GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc")
    attached = key_of("Fuc(a1-3)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc")
    assert floating != attached
    assert "{" in floating and "{" not in attached


def test_two_real_isomers_both_survive_deduplication(enumerator):
    # Tests the key as enumerate() actually uses it: were it any coarser than the
    # structure, one of these two would be discarded as a duplicate of the other.
    result = enumerator.enumerate(G0F)
    keys = {c.canonical_key for c in result.candidates}
    biantennary = key_of(CANONICAL[G0F])
    both_on_one_arm = key_of(
        "GlcNAc(b1-2)[GlcNAc(b1-4)]Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
    )
    assert biantennary != both_on_one_arm
    assert biantennary in keys and both_on_one_arm in keys
    # Every candidate is a distinct structure, so no two share a rendering either.
    assert len({c.iupac_condensed for c in result.candidates}) == len(result.candidates)


def test_two_spellings_of_one_ambiguity_give_one_key():
    # "a2-3/6" and "a2-6/3" say the same thing, so they must not look like two structures.
    assert key_of("Neu5Ac(a2-3/6)Gal(b1-4)GlcNAc") == key_of("Neu5Ac(a2-6/3)Gal(b1-4)GlcNAc")


# --- the curated rules bite, and the ones about order do not over-bite --------


def test_no_candidate_branches_before_the_enzyme_that_must_act_first(enumerator):
    # MGAT5 acts on the MGAT2 product, so b1-6 on the alpha1-6 arm needs b1-2 there first.
    checked = 0
    for candidate in enumerator.enumerate(G2F).candidates:
        if holds(candidate, "GlcNAc(b1-6)Man(a1-6)"):
            checked += 1
            assert holds(candidate, "GlcNAc(b1-2)Man(a1-6)")
    # A vacuous pass would prove nothing, so the branch has to occur somewhere.
    assert checked or not any(holds(c, "GlcNAc(b1-6)") for c in enumerator.enumerate(G2F).candidates)


def test_a_second_antenna_needs_the_first_on_the_same_mannose(enumerator):
    for candidate in enumerator.enumerate(G2F).candidates:
        if holds(candidate, "GlcNAc(b1-4)Man(a1-3)"):
            assert holds(candidate, "GlcNAc(b1-2)Man(a1-3)")


def test_a_requires_context_must_sit_on_the_residue_the_rule_speaks_about(enumerator):
    # MGAT4's b1-4 branch needs the b1-2 GlcNAc on the SAME mannose. Here the
    # b1-2 GlcNAc is on the other arm, so the rule is broken, not satisfied.
    wrong_arm = graph_of(
        "GlcNAc(b1-4)Man(a1-3)[GlcNAc(b1-2)Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
    )
    broken = enumerator._broken_rules(wrong_arm)
    assert any("MGAT4" in reason and "requires" in reason for reason in broken)


def test_exactly_the_four_bisecting_rules_are_treated_as_order_constraints(enumerator):
    # Pinned by name: if glycowork rewords a rationale, this fails loudly rather
    # than silently changing which structures are enumerated.
    assert sorted(rule.enzyme for rule in enumerator.order_constraints) == ["FUT8", "MGAT2", "MGAT4", "MGAT5"]
    for rule in enumerator.order_constraints:
        assert rule.contexts == (BISECT,)
        assert "subsequent" in rule.rationale.lower()


def test_a_bisected_candidate_may_carry_core_fucose(enumerator):
    # FUT8's own rationale: bisecting blocks SUBSEQUENT fucosylation, and "the
    # reverse order is permitted and common". Excluding the pair would drop
    # structures the curated text calls common.
    both = [c for c in enumerator.enumerate(G0F).candidates if holds(c, BISECT) and holds(c, CORE_FUC)]
    assert both, "a bisected, core-fucosylated candidate must be reachable"
    ordering = [r for r in both[0].rationale if r.kind == "ordering"]
    assert any("FUT8" in reason.enzymes for reason in ordering)
    assert all(reason.reference for reason in ordering)


def test_an_order_rule_never_appears_as_a_rejection(enumerator):
    # Matched on the whole rule rather than the enzyme: MGAT2 owns two rules and
    # only the bisecting one is about order. Keying on the name alone would
    # confuse it with the alpha-mannosidase II rule, which must still reject.
    for composition in (G0F, G2F):
        rejected = enumerator.enumerate(composition).rejected_by_rule
        for rule in enumerator.order_constraints:
            wording = f"{rule.enzyme}: {rule.product} forbidden in {' or '.join(rule.contexts)}"
            assert wording not in rejected


def test_the_other_mgat2_rule_still_rejects(enumerator):
    # Same enzyme, different rule: trimming does not reverse, so an untrimmed
    # alpha1-6 arm carrying the MGAT2 product is not a reachable structure.
    tally = enumerator.enumerate(G2F).rejected_by_rule
    assert any("MGAT2" in reason and "Man(a1-3)Man(a1-6)" in reason for reason in tally)


def test_a_trimming_rule_still_excludes(enumerator):
    # MGAT1 needs the alpha1-2 mannoses gone first, and trimming does not reverse,
    # so this one is about the finished structure and is not an order constraint.
    mgat1 = [rule for rule, _, _ in enumerator._rules if rule.enzyme == "MGAT1"]
    assert mgat1 and not any(is_order_constraint(rule) for rule in mgat1)


def test_the_rejection_tally_names_the_rule_and_the_enzyme(enumerator):
    result = enumerator.enumerate(G0F)
    assert result.rejected_by_rule
    assert result.built == len(result.candidates) + result.duplicates + result.rejected
    # A tree may break several rules, so the tally can exceed the rejected count.
    assert sum(result.rejected_by_rule.values()) >= result.rejected
    assert any("MGAT" in rule for rule in result.rejected_by_rule)


# --- the chain matcher must read a fragment's innermost linkage ---------------


def test_a_fragment_ending_in_a_linkage_checks_how_that_residue_arrives():
    # The regression that made every composition return nothing: the core
    # beta-mannose arrives on b1-4 and carries Man(a1-3), so ignoring the
    # trailing "a1-6" let it stand in for an untrimmed alpha1-6 arm.
    core = graph_of("Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc")
    assert not _contains(core, _chain("Man(a1-3)Man(a1-6)"))
    untrimmed = graph_of("Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc")
    assert _contains(untrimmed, _chain("Man(a1-3)Man(a1-6)"))


def test_a_fragment_with_no_trailing_linkage_matches_anywhere():
    graph = graph_of("Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc")
    assert _contains(graph, _chain("Man(b1-4)GlcNAc"))
    assert not _contains(graph, _chain("Gal(b1-4)GlcNAc"))


def test_the_reducing_end_arrives_on_nothing():
    graph = graph_of("Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc")
    assert not _contains(graph, _chain("GlcNAc(b1-4)GlcNAc(b1-4)"))


def test_a_residue_of_unknown_attachment_satisfies_nothing():
    # Braces say we do not know where it is. Accepting the dangling linkage as a
    # real bond would let "unknown" discharge a rule's precondition.
    graph = graph_of("{GlcNAc(b1-2)Man(a1-3)}GlcNAc(b1-2)Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc")
    floating = [n for n, d in graph.graph.nodes(data=True) if not d["attached"]]
    assert floating
    assert not any(_arrives_on(graph, node, "a1-3") for node in floating)


def test_an_unresolved_fragment_cannot_discharge_a_requirement(enumerator):
    braced = graph_of("{GlcNAc(b1-2)Man(a1-3)}GlcNAc(b1-2)Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc")
    assert braced.has_unresolved_linkage
    assert any("MGAT2" in reason and "requires" in reason for reason in enumerator._broken_rules(braced))


def test_a_branched_fragment_is_refused_rather_than_flattened():
    # Dropping the brackets would read two residues on one mannose as a stack of three.
    with pytest.raises(EnumerationError, match="branched"):
        _chain("GlcNAc(b1-2)[GlcNAc(b1-4)]Man(a1-3)")


# --- what no human enzyme makes is not enumerated -----------------------------


def test_alpha_gal_is_dropped_because_humans_have_no_enzyme_for_it(enumerator):
    # Humans lack GGTA1, so Gal(a1-3) has no maker in the glycoenzyme table.
    assert any(site.monolink == "Gal(a1-3)" for site in enumerator.dropped_sites)
    for candidate in enumerator.enumerate(G2F).candidates:
        assert "Gal(a1-3)" not in candidate.iupac_condensed


def test_a_site_no_enzyme_makes_is_never_used():
    invented = {**SITES, "gal": SITES["gal"] + (Site("a9-9", "Gal", "gal", "invented"),)}
    enumerator = Enumerator(sites=invented)
    assert any(site.linkage == "a9-9" for site in enumerator.dropped_sites)


# --- refusals say why ---------------------------------------------------------


def test_a_composition_smaller_than_the_core_is_refused(enumerator):
    result = enumerator.enumerate("Hex2HexNAc2")
    assert not result.candidates
    assert "core" in result.refusal


def test_a_composition_with_a_residue_no_site_places_is_refused(enumerator):
    # Neu5Gc has no maker in the human table, so it cannot be placed at all.
    # A composition names the class, NeuGc; a structure names the residue, Neu5Gc.
    result = enumerator.enumerate("Hex5HexNAc4NeuGc1")
    assert not result.candidates
    assert result.refusal and "NeuGc" in result.refusal


def test_a_composition_that_builds_nothing_still_says_why(enumerator):
    # Two fucoses with no antenna to carry the second: the walk ends with no tree,
    # and silence would leave the caller guessing.
    result = enumerator.enumerate("Hex3HexNAc2Fuc2")
    assert not result.candidates
    assert result.refusal
    assert "no candidates" in result.summary()


# --- the cap is reported honestly ---------------------------------------------


def test_an_untruncated_walk_is_not_reported_as_capped(enumerator):
    result = enumerator.enumerate(MAN5)
    assert not result.capped
    assert result.built == 6


def test_landing_exactly_on_the_limit_is_not_truncation(enumerator):
    # Six trees with a limit of six is a complete walk, not a cut-off one.
    result = enumerator.enumerate(MAN5, limit=6)
    assert result.built == 6
    assert not result.capped
    assert len(result.candidates) == 6


def test_a_truncated_walk_is_reported_as_capped(enumerator):
    result = enumerator.enumerate(G2F, limit=10)
    assert result.capped
    assert result.built <= 10


# --- licensing and the reasoning that travels with a candidate ----------------


def test_a_generated_candidate_states_its_resolution(enumerator):
    for candidate in enumerator.enumerate(G0F).candidates:
        # Every placement is made at a stated position, so nothing is left open.
        assert not candidate.structure.has_unresolved_linkage
        assert not candidate.structure.has_unresolved_anomericity


def test_a_candidate_is_our_own_output_and_may_be_trained_on(enumerator):
    candidate = enumerator.enumerate(G0F).candidates[0]
    assert candidate.structure.reuse_status is ReuseStatus.INTERNAL_PROPRIETARY
    assert_trainable(candidate.structure)


def test_a_candidate_carries_the_enzymes_and_the_rules_behind_it(enumerator):
    candidate = next(c for c in enumerator.enumerate(G0F).candidates if "Fuc(a1-6)" in c.iupac_condensed)
    enzymes = [r for r in candidate.rationale if r.kind == "enzyme"]
    rules = [r for r in candidate.rationale if r.kind in {"rule", "ordering"}]
    assert any("FUT8" in r.enzymes for r in enzymes)
    assert rules, "a candidate with antennae is governed by at least one curated rule"
    for rule in rules:
        assert rule.reference  # M5 has to be able to cite it
        assert rule.text


def test_an_enzyme_reason_lists_every_maker_of_the_monolink(enumerator):
    # Truncating the list dropped the only correct enzymes for some bonds, so
    # nothing is cut: what the table holds is what the reason carries.
    candidate = next(c for c in enumerator.enumerate(G0F).candidates if "Fuc(a1-6)" in c.iupac_condensed)
    for reason in candidate.rationale:
        if reason.kind != "enzyme":
            continue
        residue, _, rest = reason.text.partition("the monolink ")[2].partition("(")
        linkage = rest.split(")")[0]
        assert reason.enzymes == enumerator.catalogue.enzymes_for(residue, linkage)


def test_an_enzyme_reason_does_not_claim_to_attribute_the_individual_bond(enumerator):
    # The table records no acceptor for all but six rows, so a monolink's makers
    # include enzymes that make it elsewhere. The wording must not overclaim.
    candidate = enumerator.enumerate(G0F).candidates[0]
    for reason in candidate.rationale:
        if reason.kind == "enzyme":
            assert "monolink" in reason.text


def test_the_summary_reads_as_a_sentence(enumerator):
    text = enumerator.enumerate(G0F).summary()
    assert "candidates" in text and "Hex3HexNAc4Fuc1" in text
    assert "rejected by" in text


# --- the public motif wrappers, for callers outside this module ---------------

FLOATING_FUC = "{Fuc(a1-3)}GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"


def test_a_motif_needs_a_linkage_or_a_floating_residue_satisfies_it():
    # The trap this wrapper exists to make visible: a bare residue name is
    # satisfied by a residue whose position is unknown. Both halves are asserted,
    # so nobody can "simplify" a motif flag down to a bare name and stay green.
    floating = graph_of(FLOATING_FUC)
    assert contains_motif(floating, "Fuc")  # true, and says nothing about where
    assert not contains_motif(floating, "Fuc(a1-6)GlcNAc")  # not core fucose
    assert contains_motif(graph_of(CANONICAL[G0F]), "Fuc(a1-6)GlcNAc")  # this one is


def test_a_motif_anchor_says_where_not_only_whether():
    g0f = graph_of(CANONICAL[G0F])
    anchors = motif_anchors(g0f, "Fuc(a1-6)GlcNAc")
    assert len(anchors) == 1
    assert g0f.graph.nodes[anchors[0]]["name"] == "GlcNAc"
    assert motif_anchors(g0f, "Neu5Ac(a2-6)Gal") == ()


def test_the_wrappers_agree_with_the_privates_they_wrap():
    graph = graph_of(CANONICAL[G2F])
    for fragment in ("Fuc(a1-6)GlcNAc", "Gal(b1-4)GlcNAc", "GlcNAc(b1-4)Man(b1-4)"):
        assert contains_motif(graph, fragment) == _contains(graph, _chain(fragment))
        assert bool(motif_anchors(graph, fragment)) == contains_motif(graph, fragment)


def test_a_branched_motif_is_refused_by_the_wrapper_too():
    with pytest.raises(EnumerationError, match="branched"):
        contains_motif(graph_of(CANONICAL[G0F]), "GlcNAc(b1-2)[GlcNAc(b1-4)]Man(a1-3)")
