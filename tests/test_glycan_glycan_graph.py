"""IUPAC-condensed strings into graphs: residues as nodes, linkages as edges, unknowns kept unknown."""

import pytest

from wmxglycan import glycan_graph
from wmxglycan.composition import Composition, UnsupportedResidueError
from wmxglycan.glycan_graph import GlycanGraph, GlycanGraphError, Linkage, is_linkage_label

M3 = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
G0F_CORE = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
FLOATING_RESIDUE = "{Man(a1-2)}Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
FLOATING_FRAGMENT = "{Gal(b1-3/4)GlcNAc(b1-3/4/6)}Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"


@pytest.mark.parametrize(
    "raw, anomer, donor, acceptors, anomer_known, position_known",
    [
        ("a1-3", "a", 1, (3,), True, True),
        ("b1-4", "b", 1, (4,), True, True),
        ("a2-3/6", "a", 2, (3, 6), True, False),  # ambiguous: both positions kept, neither chosen
        ("?1-?", None, 1, (), False, False),
        ("a1-?", "a", 1, (), True, False),
        ("1-5", None, 1, (5,), False, True),  # written without an anomer
        ("a2-11", "a", 2, (11,), True, True),  # two-digit position
        ("b6-2", "b", 6, (2,), True, True),  # donor other than C1 or C2
        ("a1-3/?", "a", 1, (3,), True, False),  # one position offered, one left open
        ("a1-?/3", "a", 1, (3,), True, False),
    ],
)
def test_linkage_parsing(raw, anomer, donor, acceptors, anomer_known, position_known):
    linkage = Linkage.parse(raw)
    assert (linkage.anomer, linkage.donor, linkage.acceptors) == (anomer, donor, acceptors)
    assert linkage.anomer_known is anomer_known
    assert linkage.position_known is position_known
    assert str(linkage) == raw


def test_an_open_alternative_keeps_the_position_open():
    # "a1-3/?" says the position may be 3 or may be unknown, which is not a resolved position.
    assert Linkage.parse("a1-3/?").acceptor_unknown is True
    assert Linkage.parse("a1-3/?").position_known is False
    assert Linkage.parse("a1-3/6").acceptor_unknown is False
    assert Linkage.parse("a1-3").acceptor_unknown is False


@pytest.mark.parametrize("label", ["a1-3", "b1-4", "?1-?", "1-5", "a2-3/6", "b1-2/4/6"])
def test_linkage_labels_are_recognised(label):
    assert is_linkage_label(label)


@pytest.mark.parametrize(
    "label", ["Man", "GlcNAc", "Neu5Ac", "Gal3Me", "GlcNAc6S", "Rib5P-ol", "2,5-Anhydro-Man-ol", "D-Rha4NFo"]
)
def test_residue_names_are_not_mistaken_for_linkages(label):
    assert not is_linkage_label(label)


def test_a_branched_structure_becomes_residues_and_linkages():
    graph = GlycanGraph.from_iupac_condensed(G0F_CORE)
    assert len(graph) == 6
    assert sorted(graph.residue_names) == ["Fuc", "GlcNAc", "GlcNAc", "Man", "Man", "Man"]
    assert graph.composition() == Composition(hex=3, hexnac=2, fuc=1)
    assert graph.unattached_residues == ()
    assert graph.has_unresolved_linkage is False
    assert graph.has_unresolved_anomericity is False


def test_edges_run_from_the_reducing_end_outwards():
    graph = GlycanGraph.from_iupac_condensed(M3)
    names = graph.graph.nodes
    root = graph.root
    assert names[root]["name"] == "GlcNAc"  # the reducing end carries the rest
    assert graph.graph.in_degree(root) == 0
    joined = {
        (names[parent]["name"], names[child]["name"], str(data["linkage"]))
        for parent, child, data in graph.graph.edges(data=True)
    }
    assert ("GlcNAc", "GlcNAc", "b1-4") in joined
    assert ("Man", "Man", "a1-3") in joined and ("Man", "Man", "a1-6") in joined


def test_a_floating_residue_is_unattached_but_still_counted():
    graph = GlycanGraph.from_iupac_condensed(FLOATING_RESIDUE)
    assert graph.unattached_residues == ("Man",)
    assert graph.composition() == Composition(hex=4, hexnac=2)  # the floating residue is part of the molecule
    assert graph.has_unresolved_linkage is True  # we do not know where it sits
    assert graph.has_unresolved_anomericity is False
    assert graph.graph.nodes[graph.root]["attached"] is True


def test_a_floating_fragment_is_unattached_as_a_whole():
    graph = GlycanGraph.from_iupac_condensed(FLOATING_FRAGMENT)
    # Both residues of "{Gal(b1-3/4)GlcNAc(b1-3/4/6)}" hang off the structure at an unknown point.
    assert sorted(graph.unattached_residues) == ["Gal", "GlcNAc"]
    assert graph.root is not None and graph.graph.nodes[graph.root]["name"] == "GlcNAc"
    assert graph.has_unresolved_linkage is True
    internal = {
        (graph.graph.nodes[p]["name"], graph.graph.nodes[c]["name"])
        for p, c in graph.graph.edges
        if not graph.graph.nodes[p]["attached"]
    }
    assert ("GlcNAc", "Gal") in internal  # the fragment keeps the linkage it does know


@pytest.mark.parametrize(
    "text, unresolved_linkage, unresolved_anomer",
    [
        (M3, False, False),
        ("Neu5Ac(a2-3/6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc", True, False),
        ("Man(?1-?)Man(b1-4)GlcNAc", True, True),
        ("Man(a1-?)Man(b1-4)GlcNAc", True, False),
    ],
)
def test_what_the_string_leaves_open_stays_open(text, unresolved_linkage, unresolved_anomer):
    graph = GlycanGraph.from_iupac_condensed(text)
    assert graph.has_unresolved_linkage is unresolved_linkage
    assert graph.has_unresolved_anomericity is unresolved_anomer


def test_a_composition_is_refused_rather_than_rounded_off():
    graph = GlycanGraph.from_iupac_condensed("Xyl(b1-2)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc")
    with pytest.raises(UnsupportedResidueError) as caught:
        graph.composition()
    assert caught.value.names == ("Xyl",)
    assert "Xyl" in str(caught.value)


@pytest.mark.parametrize("text", ["", "   ", None])
def test_an_empty_string_is_refused(text):
    with pytest.raises(GlycanGraphError):
        GlycanGraph.from_iupac_condensed(text)


def test_a_parser_failure_becomes_our_own_error(monkeypatch):
    def explode(text):
        raise RuntimeError("glycowork fell over")

    monkeypatch.setattr(glycan_graph, "_glycowork_parser", lambda: explode)
    with pytest.raises(GlycanGraphError, match="glycowork fell over"):
        GlycanGraph.from_iupac_condensed(M3)


def test_a_bad_linkage_label_is_refused():
    with pytest.raises(GlycanGraphError, match="cannot read"):
        Linkage.parse("sideways")


def test_repr_names_the_structure():
    assert M3 in repr(GlycanGraph.from_iupac_condensed(M3))
