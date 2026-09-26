"""Biosynthetic constraints: loaded whole, with the reasoning that has to travel with them."""

import pytest

from wmxglycan.constraints import BiosyntheticConstraint, ConstraintKind, ConstraintSet, GlycanClass

HEADER = "﻿product,context,kind,glycan_class,enzyme,rationale,reference\n"  # the shipped file carries a BOM
ROWS = (
    "GlcNAc(b1-2)Man(a1-3),Man(a1-2)Man,forbids,N,MGAT1,alpha1-2 mannoses must be trimmed first,Schachter 1986\n"
    "GlcNAc(b1-6)GalNAc,Gal(b1-3)GalNAc|GlcNAc(b1-3)GalNAc,requires,O,GCNT1/GCNT3,core 2 and core 4 branch,Brockhausen 2022\n"
    "GalNAc(a1-3)Gal,Fuc(a1-2)Gal,requires,,ABO,the A transferase needs the H antigen,Yamamoto 2004\n"
)
COLUMNS = ("product", "context", "kind", "glycan_class", "enzyme", "rationale", "reference")


def one_row(**overrides):
    fields = {
        "product": "GlcNAc(b1-2)Man(a1-3)",
        "context": "Man(a1-2)Man",
        "kind": "forbids",
        "glycan_class": "N",
        "enzyme": "MGAT1",
        "rationale": "why",
        "reference": "somewhere",
    }
    fields.update(overrides)
    return HEADER + ",".join(fields[name] for name in COLUMNS) + "\n"


def test_a_rule_keeps_its_enzyme_rationale_and_reference():
    rule = ConstraintSet.from_csv(HEADER + ROWS).for_enzyme("MGAT1")[0]
    assert rule.product == "GlcNAc(b1-2)Man(a1-3)"
    assert rule.contexts == ("Man(a1-2)Man",)
    assert rule.kind is ConstraintKind.FORBIDS
    assert rule.glycan_class is GlycanClass.N_LINKED
    assert rule.rationale.startswith("alpha1-2 mannoses")
    assert rule.reference == "Schachter 1986"


def test_alternative_contexts_are_split_on_the_bar():
    rule = ConstraintSet.from_csv(HEADER + ROWS).for_product("GlcNAc(b1-6)GalNAc")[0]
    assert rule.contexts == ("Gal(b1-3)GalNAc", "GlcNAc(b1-3)GalNAc")


def test_a_rule_without_a_class_applies_to_every_class():
    rule = ConstraintSet.from_csv(HEADER + ROWS).for_enzyme("ABO")[0]
    assert rule.glycan_class is None
    assert rule.applies_to(GlycanClass.N_LINKED) and rule.applies_to(GlycanClass.O_LINKED)


def test_queries():
    rules = ConstraintSet.from_csv(HEADER + ROWS)
    assert len(rules) == 3
    assert len(rules.of_kind(ConstraintKind.REQUIRES)) == 2
    assert {rule.enzyme for rule in rules.for_glycan_class(GlycanClass.N_LINKED)} == {"MGAT1", "ABO"}
    assert {rule.enzyme for rule in rules.for_glycan_class(GlycanClass.O_LINKED)} == {"GCNT1/GCNT3", "ABO"}
    assert rules.enzymes == ("ABO", "GCNT1/GCNT3", "MGAT1")
    assert len(list(rules)) == 3


# --- the sentence shown beside a ranked candidate ----------------------------


def test_a_forbids_rule_reads_as_a_bar_on_the_product():
    rule = ConstraintSet.from_csv(HEADER + ROWS).for_enzyme("MGAT1")[0]
    assert rule.explain().startswith("GlcNAc(b1-2)Man(a1-3) cannot form where Man(a1-2)Man is present")


def test_a_requires_rule_reads_as_a_precondition_not_an_obligation():
    # The context is what must already be there; it does not force the product wherever it appears.
    rule = ConstraintSet.from_csv(HEADER + ROWS).for_enzyme("GCNT1/GCNT3")[0]
    explanation = rule.explain()
    assert explanation.startswith(
        "GlcNAc(b1-6)GalNAc can only form where Gal(b1-3)GalNAc or GlcNAc(b1-3)GalNAc is already present"
    )
    assert "is required in the context" not in explanation


def test_explain_shows_the_reasoning_a_ranked_candidate_will_need():
    explanation = ConstraintSet.from_csv(HEADER + ROWS).for_enzyme("MGAT1")[0].explain()
    assert "alpha1-2 mannoses must be trimmed first" in explanation
    assert "MGAT1" in explanation and "Schachter 1986" in explanation


# --- a rule that cannot explain itself does not load -------------------------


@pytest.mark.parametrize("column", ["enzyme", "rationale", "reference", "product", "context", "kind"])
def test_a_rule_that_lost_part_of_itself_stops_the_load(column):
    with pytest.raises(ValueError, match="line 2"):
        ConstraintSet.from_csv(one_row(**{column: ""}))


@pytest.mark.parametrize("placeholder", ["-", "N/A", "unknown", "?", "none", "n.d."])
@pytest.mark.parametrize("column", ["enzyme", "rationale", "reference"])
def test_a_placeholder_is_not_reasoning(column, placeholder):
    # "N/A" beside a ranked candidate would be an invented citation, not an explanation.
    with pytest.raises(ValueError, match="line 2"):
        ConstraintSet.from_csv(one_row(**{column: placeholder}))


def test_an_unknown_kind_stops_the_load():
    with pytest.raises(ValueError, match="line 2"):
        ConstraintSet.from_csv(one_row(kind="encourages"))


def test_missing_columns_are_an_error():
    with pytest.raises(ValueError, match="missing the columns"):
        ConstraintSet.from_csv("product,context\nx,y\n")


def test_an_empty_file_is_an_error():
    with pytest.raises(ValueError):
        ConstraintSet.from_csv(HEADER)


def test_a_constraint_needs_at_least_one_context():
    with pytest.raises(ValueError):
        BiosyntheticConstraint(
            product="p", contexts=(), kind="forbids", glycan_class=None, enzyme="e", rationale="r", reference="ref"
        )


# --- the rules that actually ship in the pinned glycowork release ------------


def test_the_shipped_rules_load_whole():
    rules = ConstraintSet.from_glycowork()
    assert len(rules) == 15
    assert len(rules.of_kind(ConstraintKind.FORBIDS)) == 9
    assert len(rules.of_kind(ConstraintKind.REQUIRES)) == 6
    assert len(rules.for_glycan_class(GlycanClass.N_LINKED)) == 11  # ten N rules plus the class-agnostic one
    assert len(rules.for_glycan_class(GlycanClass.O_LINKED)) == 5
    assert {"MGAT1", "MGAT2", "MGAT3", "MGAT4", "MGAT5", "FUT8"} <= set(rules.enzymes)


def test_every_shipped_rule_carries_its_reasoning():
    for rule in ConstraintSet.from_glycowork():
        assert rule.enzyme and rule.rationale and rule.reference
        assert rule.contexts and all(rule.contexts)
        assert rule.rationale in rule.explain() and rule.reference in rule.explain()


def test_no_shipped_rule_is_explained_backwards():
    for rule in ConstraintSet.from_glycowork():
        explanation = rule.explain()
        assert explanation.startswith(rule.product)  # the product stays the subject
        expected = "cannot form where" if rule.kind is ConstraintKind.FORBIDS else "can only form where"
        assert expected in explanation
