"""A biological claim must come from the graph, never from a string.

Three defects of one shape have been caught in this project: a chitinase
credited with making a bond, a salivary amylase credited the same way, and a
motif pattern that called antennary fucose core fucose. Each was a string
standing in for a claim about chemistry.

The rule: any flag, column or label asserting something biological must be
derived by inspecting the graph at the node it is about. Not by matching a
residue-and-linkage pattern, and not by reading how a name is spelt.

This file enforces the rule structurally, so a future edit that reaches for a
substring has to fail here rather than ship a wrong claim. The behavioural
counterparts live beside the features they test.
"""

import ast
import inspect
from pathlib import Path

import pytest

import wmxglycan.enzymes as enzymes_module
import wmxglycan.features as features_module
from wmxglycan.enumeration import contains_motif, motif_anchors
from wmxglycan.features import extract_structure
from wmxglycan.glycan_graph import GlycanGraph

SOURCE_DIR = Path(__file__).resolve().parent.parent / "src" / "wmxglycan"

# Positional claims: each names a residue AND the place it must sit.
POSITIONAL_FLAGS = ("core_fucose", "bisecting_glcnac")

# The matcher walks the graph, but it matches a CHAIN, so it can only say what the
# chain spells. It cannot express "at the reducing end" or "on the beta-mannose",
# which is exactly how antennary fucose came to be scored as core fucose.
CHAIN_MATCHERS = {"contains_motif", "motif_anchors", "_contains", "_nodes_matching"}

CORE_FUCOSE = "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
ANTENNARY_FUCOSE = "Fuc(a1-6)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
FLOATING_PAIR = "{Fuc(a1-6)GlcNAc(b1-4)}Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"


def called_names(module) -> set[str]:
    tree = ast.parse(inspect.getsource(module))
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            called = node.func
            if isinstance(called, ast.Name):
                names.add(called.id)
            elif isinstance(called, ast.Attribute):
                names.add(called.attr)
    return names


# --- the rule, enforced structurally ------------------------------------------


def test_the_featuriser_never_asks_a_chain_matcher_for_a_biological_flag():
    assert not (called_names(features_module) & CHAIN_MATCHERS)


def test_the_featuriser_does_not_even_import_the_chain_matcher():
    # Removing the import is what makes the rule hard to breach by accident.
    tree = ast.parse(inspect.getsource(features_module))
    imported = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    }
    assert not (imported & CHAIN_MATCHERS)


@pytest.mark.parametrize("flag", POSITIONAL_FLAGS)
def test_every_positional_flag_is_anchored_to_a_node(flag):
    # Anchored means the computation names the node the claim is about, rather
    # than searching the whole structure for a pattern.
    source = inspect.getsource(features_module)
    line = next(text for text in source.splitlines() if f'"{flag}"' in text and "float(" in text)
    assert "_child_on(" in line, f"{flag} is not anchored to a node: {line.strip()}"


def test_no_biological_flag_is_read_off_the_structure_string():
    # A substring test on the IUPAC text would be the crudest form of the bug.
    source = inspect.getsource(features_module)
    for banned in ('in iupac', 'in text', '.startswith("Fuc', '.startswith("GlcNAc'):
        assert banned not in source, banned


def test_the_enzyme_table_exposes_no_activity_it_cannot_know():
    # glycoclass in the shipped table is a gene-symbol family root, not an
    # activity class, so nothing may present enzyme activity as though it were
    # known. There is no API here that claims to.
    public = {name for name in dir(enzymes_module.EnzymeCatalogue) if not name.startswith("_")}
    for claim in ("is_transferase", "is_hydrolase", "activity", "enzyme_class", "glycoclass", "ec_number"):
        assert claim not in public, claim


# --- and the behaviour the rule exists to protect -------------------------------


def test_the_chain_matcher_really_cannot_express_a_positional_claim():
    # Not a hypothetical. This is why the rule exists, asserted directly, so
    # nobody re-adopts the matcher for a flag believing it would be safe.
    antennary = GlycanGraph.from_iupac_condensed(ANTENNARY_FUCOSE)
    assert contains_motif(antennary, "Fuc(a1-6)GlcNAc")  # the pattern says yes...
    assert extract_structure(ANTENNARY_FUCOSE)["core_fucose"] == 0.0  # ...the anchored flag says no
    assert motif_anchors(antennary, "Fuc(a1-6)GlcNAc")  # it even reports an anchor


def test_a_fragment_of_unknown_placement_makes_no_positional_claim():
    floating = GlycanGraph.from_iupac_condensed(FLOATING_PAIR)
    assert contains_motif(floating, "Fuc(a1-6)GlcNAc")
    assert extract_structure(FLOATING_PAIR)["core_fucose"] == 0.0


def test_the_anchored_flag_still_finds_the_real_thing():
    # A rule that only ever says no would be useless.
    assert extract_structure(CORE_FUCOSE)["core_fucose"] == 1.0
