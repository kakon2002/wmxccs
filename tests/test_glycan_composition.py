"""Composition parsing, canonical form, mass, the Man3GlcNAc2 core check and N-glycan rules."""

import itertools

import pytest

from wmxglycan.composition import (
    RESIDUE_MASS,
    WATER_MASS,
    Composition,
    CompositionError,
    Residue,
    canonical_composition,
    has_full_man3glcnac2_core,
    is_plausible_n_glycan,
    man3glcnac2_core_warnings,
    n_glycan_implausibility_reasons,
    parse_composition,
)


def test_parse():
    comp = parse_composition("Hex5HexNAc4Fuc1")
    assert comp == Composition(hex=5, hexnac=4, fuc=1)
    assert comp.counts() == {"Hex": 5, "HexNAc": 4, "Fuc": 1}


def test_every_residue_order_gives_one_composition():
    spellings = {"".join(p) for p in itertools.permutations(["Hex5", "HexNAc4", "Fuc1", "NeuAc2"])}
    assert len(spellings) == 24
    assert {parse_composition(s) for s in spellings} == {Composition(hex=5, hexnac=4, fuc=1, neuac=2)}
    assert {canonical_composition(s) for s in spellings} == {"Hex5HexNAc4Fuc1NeuAc2"}


def test_dhex_is_read_as_fuc():
    assert parse_composition("Hex5HexNAc4dHex1") == parse_composition("Hex5HexNAc4Fuc1")
    assert canonical_composition("dHex1HexNAc4Hex5") == "Hex5HexNAc4Fuc1"


def test_canonical_order_and_zero_counts():
    assert canonical_composition("NeuGc1NeuAc1Fuc1HexNAc4Hex5") == "Hex5HexNAc4Fuc1NeuAc1NeuGc1"
    assert canonical_composition("Hex5HexNAc2Fuc0") == "Hex5HexNAc2"


def test_canonical_round_trips():
    comp = parse_composition("NeuAc2dHex1HexNAc4Hex5")
    assert parse_composition(str(comp)) == comp
    assert hash(parse_composition(str(comp))) == hash(comp)


def test_outer_whitespace_is_ignored():
    assert parse_composition("  Hex5HexNAc2\n") == Composition(hex=5, hexnac=2)


@pytest.mark.parametrize(
    "text",
    [
        "Hex5HexNac4",  # wrong case in HexNAc
        "hex5hexnac4",  # lower case
        "Hex5HexNAc4Fuk1",  # misspelt residue
        "Hex5HexNAc4Fuc",  # missing count
        "Hex5 HexNAc4",  # internal space
        "Hex5,HexNAc4",  # separator
        "Hex5HexNAc4Fuc1x",  # trailing junk
        "5Hex",  # count before name
        "Hex-1",  # negative count
        "Hex3.5HexNAc2",  # fractional count
        "Hex05HexNAc4",  # leading zero
        "Hex3Hex2HexNAc4",  # repeated residue
        "Fuc1dHex1Hex3HexNAc2",  # repeated through the dHex alias
        "H5N4F1",  # short notation is not supported
        "Neu5Ac1Hex5HexNAc4",  # Neu5Ac spelling is not supported
        "Hex５HexNAc4",  # full-width digit five
        "Hex0",  # no residues
        "",
        "   ",
    ],
)
def test_malformed_compositions_raise(text):
    with pytest.raises(CompositionError):
        parse_composition(text)


def test_a_typo_is_not_partially_parsed():
    # "Hex5HexNAc4" is a valid prefix; it must not come back as the answer.
    with pytest.raises(CompositionError, match="Fuk"):
        parse_composition("Hex5HexNAc4Fuk1")


def test_error_suggests_the_right_case():
    with pytest.raises(CompositionError, match="did you mean 'HexNAc'"):
        parse_composition("Hex5HexNac4")


@pytest.mark.parametrize("value", [None, 5, b"Hex5"])
def test_non_string_input_is_a_type_error(value):
    with pytest.raises(TypeError):
        parse_composition(value)


@pytest.mark.parametrize("counts", [{"hex": -1, "hexnac": 2}, {"hex": True}, {"hex": 2.0}, {}])
def test_direct_construction_is_validated(counts):
    with pytest.raises(CompositionError):
        Composition(**counts)


# Published free-glycan monoisotopic masses, quoted to 2 decimals, so the
# tolerance is half the last digit. A missing water would be off by 18.01.
@pytest.mark.parametrize(
    "text, published",
    [
        ("Hex3HexNAc4Fuc1", 1462.54),  # G0F
        ("Hex5HexNAc4Fuc1", 1786.65),  # G2F
        ("Hex5HexNAc2", 1234.43),  # Man5
    ],
)
def test_mass_matches_published_values(text, published):
    assert parse_composition(text).monoisotopic_mass == pytest.approx(published, abs=0.005)


def test_mass_includes_exactly_one_water():
    comp = parse_composition("Hex5HexNAc4Fuc1")
    assert comp.monoisotopic_mass - comp.residue_mass == pytest.approx(WATER_MASS, abs=1e-9)


# Independent check of the residue mass constants from elemental formulas.
# Monoisotopic atomic masses from the AME2020 evaluation; C-12 is exact by definition.
ATOM = {"C": 12.0, "H": 1.00782503223, "N": 14.00307400443, "O": 15.99491461957}
RESIDUE_FORMULA = {  # monosaccharide minus H2O
    Residue.HEX: {"C": 6, "H": 10, "O": 5},
    Residue.HEXNAC: {"C": 8, "H": 13, "N": 1, "O": 5},
    Residue.FUC: {"C": 6, "H": 10, "O": 4},
    Residue.NEUAC: {"C": 11, "H": 17, "N": 1, "O": 8},
    Residue.NEUGC: {"C": 11, "H": 17, "N": 1, "O": 9},
}


@pytest.mark.parametrize("residue", list(Residue))
def test_residue_masses_agree_with_elemental_formulas(residue):
    from_atoms = sum(ATOM[element] * n for element, n in RESIDUE_FORMULA[residue].items())
    assert RESIDUE_MASS[residue] == pytest.approx(from_atoms, abs=1e-5)


def test_water_mass_agrees_with_elements():
    assert WATER_MASS == pytest.approx(2 * ATOM["H"] + ATOM["O"], abs=1e-5)


# --- Man3GlcNAc2 core and N-glycan rules ------------------------------------


@pytest.mark.parametrize(
    "text",
    [
        "Hex3HexNAc4Fuc1",  # G0F
        "Hex4HexNAc4Fuc1",  # G1F
        "Hex5HexNAc4Fuc1",  # G2F
        "Hex5HexNAc2",  # Man5
        "Hex3HexNAc2Fuc1",  # core plus core fucose
        "Hex5HexNAc4Fuc1NeuAc2",  # sialylated, with antennae
        "Hex5HexNAc4Fuc2",  # second fucose, with antennae
        "Hex4HexNAc3NeuAc1",  # one antenna carrying a sialic acid
    ],
)
def test_plausible_n_glycans_with_the_full_core(text):
    assert has_full_man3glcnac2_core(text)
    assert man3glcnac2_core_warnings(text) == []
    assert is_plausible_n_glycan(text)
    assert n_glycan_implausibility_reasons(text) == []


@pytest.mark.parametrize(
    "text, fragment",
    [
        ("Hex2HexNAc2Fuc1", "paucimannosidic"),  # two mannoses and a core fucose
        ("Hex1HexNAc2", "paucimannosidic"),
        ("Hex5HexNAc1", "endoglycosidase"),  # one core GlcNAc left
    ],
)
def test_a_missing_core_is_a_warning_not_a_rejection(text, fragment):
    warnings = man3glcnac2_core_warnings(text)
    assert len(warnings) == 1
    assert fragment in warnings[0] and "Man3GlcNAc2" in warnings[0]
    assert not has_full_man3glcnac2_core(text)
    assert is_plausible_n_glycan(text)
    assert n_glycan_implausibility_reasons(text) == []


@pytest.mark.parametrize("text", ["HexNAc2Fuc1", "Hex2HexNAc4", "Hex2HexNAc1"])
def test_paucimannose_is_named_only_for_one_or_two_mannoses_on_the_chitobiose_core(text):
    warnings = man3glcnac2_core_warnings(text)
    assert warnings
    assert not any("paucimannosidic" in warning for warning in warnings)


@pytest.mark.parametrize("text", ["Hex5", "Hex9", "Fuc1", "Hex3Fuc1", "Hex1Fuc2NeuAc1"])
def test_a_composition_with_no_hexnac_is_not_an_n_glycan(text):
    reasons = n_glycan_implausibility_reasons(text)
    assert len(reasons) == 1 and reasons[0].startswith("HexNAc0")
    assert not is_plausible_n_glycan(text)
    # No endoglycosidase leaves a glycan without any GlcNAc, so none is offered as an explanation.
    assert not any("endoglycosidase" in warning for warning in man3glcnac2_core_warnings(text))


@pytest.mark.parametrize(
    "text, fragment",
    [
        ("Hex5HexNAc2NeuAc1", "sialic acid"),  # rule a
        ("Hex5HexNAc2NeuGc1", "sialic acid"),  # rule a, NeuGc
        ("Hex5HexNAc2Fuc2", "at most Fuc1"),  # rule b
    ],
)
def test_each_mammalian_rule_gives_one_readable_reason(text, fragment):
    reasons = n_glycan_implausibility_reasons(text)
    assert len(reasons) == 1
    assert fragment in reasons[0]
    assert not is_plausible_n_glycan(text)
    assert has_full_man3glcnac2_core(text)


def test_rules_and_core_warnings_are_reported_separately_and_completely():
    text = "Hex1HexNAc1Fuc2NeuAc1"
    assert len(n_glycan_implausibility_reasons(text)) == 2
    assert len(man3glcnac2_core_warnings(text)) == 2


def test_core_and_rule_checks_reject_typos_too():
    for check in (is_plausible_n_glycan, has_full_man3glcnac2_core):
        with pytest.raises(CompositionError):
            check("Hex5HexNac4")
