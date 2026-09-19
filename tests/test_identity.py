"""The analyte union and the matched-ion key: when two measurements are the same ion.

Every test here holds one promise. A match is made on a DECLARED IDENTIFIER, an
adduct, a signed charge, the gas the value refers to and the structural state,
and never on what a paper happens to print in its compound column. Two failure
shapes are equally bad and both are tested for:

- keying APART two records that are the same ion, which reports a clean run over
  half the pairs it should have found;
- keying TOGETHER two records that are not the same ion, which reports the
  difference between two molecules, two charge states or a folded and an
  unfolded protein as inter-platform bias.

So every guard is tested in both directions: the valid case passes, and the
invalid case raises with a message a human can act on.
"""

from __future__ import annotations

import warnings

import pytest
from pydantic import TypeAdapter, ValidationError

from conftest import (
    INCHIKEY,
    LNH_IUPAC,
    LNNH_IUPAC,
    OTHER_INCHIKEY,
    adc,
    antibody,
    antibody_identity,
    glycan,
    measurement,
    peptide,
    protein,
    small_molecule,
)
from wmxccs.identity import (
    Analyte,
    AntibodyIdentity,
    Composition,
    CompositionError,
    Derivatisation,
    FoldingState,
    GlycanAnalyte,
    PeptideAnalyte,
    ProteinAnalyte,
    ReducingEndLabel,
    SmallMoleculeAnalyte,
    UnverifiedFormatWarning,
    adduct_carrier_is_unstated,
    adduct_components,
    canonical_composition,
    parse_adduct,
    parse_composition,
)

# Two distinct WURCS strings. Neither is asserted to describe any particular
# sugar; they are here as two different structure identifiers.
WURCS_ONE = "WURCS=2.0/1,1,0/[a2122h-1b_1-5]/1/"
WURCS_TWO = "WURCS=2.0/2,2,1/[a2122h-1b_1-5][a1122h-1b_1-5]/1-2/a4-b1"

# The composition LNH and LNnH share. Given for both, so that a platform keying
# on composition alone would call the two one molecule.
LNH_COMPOSITION = "Hex4HexNAc2"


# --- the four assertions the brief requires -------------------------------------------


def test_a_native_24_plus_and_a_denatured_40_plus_ion_of_one_antibody_are_different_matched_ions():
    native_24 = measurement(
        analyte=antibody(folding_state=FoldingState.NATIVE), adduct="[M+24H]24+", charge=24, ccs=7000.0
    )
    denatured_40 = measurement(
        analyte=antibody(folding_state=FoldingState.DENATURED), adduct="[M+40H]40+", charge=40, ccs=11000.0
    )
    assert native_24.matched_ion_key != denatured_40.matched_ion_key


def test_a_native_24_plus_and_a_native_40_plus_antibody_ion_differ_by_charge_alone():
    """Proves the charge half of the assertion above. Without it, that test passes
    even when only the folding state reaches the key."""
    native_24 = measurement(
        analyte=antibody(folding_state=FoldingState.NATIVE), adduct="[M+24H]24+", charge=24, ccs=7000.0
    )
    native_40 = measurement(
        analyte=antibody(folding_state=FoldingState.NATIVE), adduct="[M+40H]40+", charge=40, ccs=7000.0
    )
    assert native_24.matched_ion_key != native_40.matched_ion_key


def test_a_native_40_plus_and_a_denatured_40_plus_antibody_ion_differ_by_folding_state_alone():
    """The other half: one charge state, one adduct, one gas, and still two ions,
    because a folded antibody and an unfolded one have different cross sections."""
    native_40 = measurement(
        analyte=antibody(folding_state=FoldingState.NATIVE), adduct="[M+40H]40+", charge=40, ccs=7000.0
    )
    denatured_40 = measurement(
        analyte=antibody(folding_state=FoldingState.DENATURED), adduct="[M+40H]40+", charge=40, ccs=11000.0
    )
    assert native_40.matched_ion_key != denatured_40.matched_ion_key


def test_two_analytes_sharing_a_display_name_are_not_the_same_matched_ion():
    one = measurement(analyte=small_molecule(display_name="the standard", inchikey=INCHIKEY))
    other = measurement(analyte=small_molecule(display_name="the standard", inchikey=OTHER_INCHIKEY))
    assert one.matched_ion_key != other.matched_ion_key


def test_two_antibodies_sharing_a_display_name_are_not_the_same_matched_ion():
    one = measurement(
        analyte=antibody(display_name="the mAb", antibody=antibody_identity(inn="trastuzumab")),
        adduct="[M+24H]24+",
        charge=24,
        ccs=7000.0,
    )
    other = measurement(
        analyte=antibody(display_name="the mAb", antibody=antibody_identity(inn="rituximab")),
        adduct="[M+24H]24+",
        charge=24,
        ccs=7000.0,
    )
    assert one.matched_ion_key != other.matched_ion_key


@pytest.mark.parametrize(
    "build",
    [
        pytest.param(lambda name: small_molecule(display_name=name), id="small_molecule"),
        pytest.param(lambda name: peptide(display_name=name), id="peptide"),
        pytest.param(lambda name: glycan(display_name=name), id="glycan"),
        pytest.param(lambda name: protein(display_name=name), id="protein"),
    ],
)
def test_changing_the_display_name_alone_does_not_change_the_matched_ion_key(build):
    printed_one_way = measurement(analyte=build("Man5"))
    printed_another_way = measurement(analyte=build("high-mannose 5"))
    assert printed_one_way.matched_ion_key == printed_another_way.matched_ion_key


@pytest.mark.parametrize(
    "one_spelling, other_spelling, charge",
    [
        ("[M+Na+H]2+", "[M+H+Na]2+", 2),
        ("[M+1H]1+", "[M+H]+", 1),
    ],
)
def test_two_spellings_of_one_adduct_give_one_matched_ion_key(one_spelling, other_spelling, charge):
    one = measurement(adduct=one_spelling, charge=charge)
    other = measurement(adduct=other_spelling, charge=charge)
    assert one.matched_ion_key == other.matched_ion_key


# --- the adduct ------------------------------------------------------------------------


@pytest.mark.parametrize(
    "written, canonical, charge",
    [
        ("[M+H]+", "[M+H]+", 1),
        ("[M+1H]1+", "[M+H]+", 1),
        ("[M+Na+H]2+", "[M+H+Na]2+", 2),
        ("[M+H+Na]2+", "[M+H+Na]2+", 2),
        ("[M-H]-", "[M-H]-", -1),
        ("[M+2H]2+", "[M+2H]2+", 2),
        ("[2M+Na]+", "[2M+Na]+", 1),
        ("[M+H-H2O]+", "[M+H-H2O]+", 1),
    ],
)
def test_an_adduct_canonicalises_to_one_spelling_and_a_signed_charge(written, canonical, charge):
    assert parse_adduct(written) == (canonical, charge)


def test_an_adduct_naming_one_species_twice_is_refused_rather_than_merged():
    """Merging +H+H into +2H would guess that a transcription fault was a count."""
    with pytest.raises(ValueError, match="more than once"):
        parse_adduct("[M+H+H]2+")


@pytest.mark.parametrize("written", ["M+H", "[M+H]", "[M+H]2", "[X+H]+", "", "[M+h]+", "[M+H]+2"])
def test_an_adduct_not_in_bracket_notation_is_refused(written):
    with pytest.raises(ValueError):
        parse_adduct(written)


def test_adduct_components_come_back_in_canonical_order_with_their_counts():
    assert adduct_components("[M+Na+H]2+") == (("+", 1, "H"), ("+", 1, "Na"))
    assert adduct_components("[M+2H]2+") == (("+", 2, "H"),)
    assert adduct_components("[M+H-H2O]+") == (("+", 1, "H"), ("-", 1, "H2O"))


def test_an_unstated_charge_carrier_parses_and_keeps_its_charge():
    assert parse_adduct("[M+24?]24+") == ("[M+24?]24+", 24)
    assert adduct_components("[M+24?]24+") == (("+", 24, "?"),)


@pytest.mark.parametrize(
    "written, unstated",
    [
        ("[M+24?]24+", True),
        ("[M+24H]24+", False),
        ("[M+H]+", False),
        ("[M-H]-", False),
        ("[M+H+?]2+", True),
    ],
)
def test_adduct_carrier_is_unstated_sees_only_the_unnamed_carrier(written, unstated):
    assert adduct_carrier_is_unstated(written) is unstated


def test_an_unstated_charge_carrier_keys_apart_from_a_protonated_ion_of_the_same_charge():
    """Twenty-four unnamed charges are not twenty-four protons: the mass and the
    cross section both depend on which, so the two may never be pooled."""
    unstated = measurement(
        analyte=antibody(folding_state=FoldingState.NATIVE), adduct="[M+24?]24+", charge=24, ccs=7000.0
    )
    protonated = measurement(
        analyte=antibody(folding_state=FoldingState.NATIVE), adduct="[M+24H]24+", charge=24, ccs=7000.0
    )
    assert unstated.matched_ion_key != protonated.matched_ion_key


# --- text that is actually a missing value ---------------------------------------------

# Pinned here rather than imported from the module, so that deleting a spelling
# from the module's list fails a test instead of quietly narrowing the guard.
PLACEHOLDER_SPELLINGS = [
    "-",
    "?",
    "#n/a",
    "<na>",
    "missing",
    "n.a.",
    "n.d.",
    "n/a",
    "na",
    "nan",
    "nat",
    "nd",
    "none",
    "not available",
    "not determined",
    "not reported",
    "null",
    "tbd",
    "unknown",
]


@pytest.mark.parametrize("spelling", PLACEHOLDER_SPELLINGS)
def test_every_placeholder_spelling_is_refused_where_text_is_required(spelling):
    with pytest.raises(ValidationError, match="placeholder"):
        small_molecule(display_name=spelling)


@pytest.mark.parametrize("spelling", ["N/A", "  n/a  ", "Unknown", "NOT REPORTED"])
def test_a_placeholder_is_refused_however_it_is_cased_or_padded(spelling):
    with pytest.raises(ValidationError, match="placeholder"):
        small_molecule(display_name=spelling)


@pytest.mark.parametrize("spelling", ["--", "...", "???", "—"])
def test_a_string_with_no_letter_or_digit_is_refused_as_a_placeholder(spelling):
    with pytest.raises(ValidationError, match="placeholder"):
        small_molecule(display_name=spelling)


def test_a_placeholder_is_refused_in_a_required_text_field_too():
    with pytest.raises(ValidationError, match="placeholder"):
        small_molecule(source="not reported")


@pytest.mark.parametrize("name", ["Man5", "LNH", "trastuzumab", "2-AB"])
def test_ordinary_text_is_kept_and_stripped(name):
    assert small_molecule(display_name=f"  {name} ").display_name == name


# --- a record that skipped validation ---------------------------------------------------


def test_a_record_built_with_model_construct_is_caught_as_unvalidated():
    """model_construct skips every validator, so these modifications stay unsorted
    and the record is not what validation would have stored. The gate must notice."""
    smuggled = PeptideAnalyte.model_construct(
        sequence="PEPTIDE", modifications=("Phospho@S5", "Acetyl@K2"), source="a test"
    )
    blockers = smuggled.training_blockers()
    assert blockers, "a record built without validation was trusted"
    assert any("without validation" in blocker for blocker in blockers)


def test_a_record_whose_text_was_never_stripped_is_caught_as_unvalidated():
    smuggled = ProteinAnalyte.model_construct(accession="P01857", source="  a test  ")
    assert any("without validation" in blocker for blocker in smuggled.training_blockers())


def test_a_properly_built_record_has_no_revalidation_blocker():
    assert peptide(modifications=("Acetyl@K2", "Phospho@S5")).training_blockers() == []
    assert protein(folding_state=FoldingState.NATIVE).training_blockers() == []


# --- small molecules --------------------------------------------------------------------


def test_a_small_molecule_without_an_inchikey_cannot_be_built():
    with pytest.raises(ValidationError):
        SmallMoleculeAnalyte(source="a test", display_name="the standard")


@pytest.mark.parametrize(
    "written",
    [
        "RYYVLZVUVIJVGH-UHFFFAOYSA",
        "RYYVLZVUVIJV-UHFFFAOYSA-N",
        "ryyvlzvuvijvgh-uhfffaoysa-n",
        "RYYVLZVUVIJVGH UHFFFAOYSA N",
        "RYYVLZVUVIJVGH-UHFFFAOYSA-NN",
        "InChIKey=RYYVLZVUVIJVGH-UHFFFAOYSA-N",
        "C1=CC=CC=C1",
    ],
)
def test_an_inchikey_of_the_wrong_shape_is_refused(written):
    with pytest.raises(ValidationError, match="is not an InChIKey"):
        small_molecule(inchikey=written)


def test_a_well_formed_inchikey_is_the_whole_identity_of_a_small_molecule():
    assert INCHIKEY in small_molecule().identity_key()
    assert small_molecule().identity_key() != small_molecule(inchikey=OTHER_INCHIKEY).identity_key()


def test_smiles_is_carried_but_is_never_the_identity():
    """One molecule has many valid SMILES and one InChIKey, so keying on SMILES
    would split one ion into as many keys as there are ways of writing it."""
    one_way = small_molecule(smiles="C1=CC=CC=C1")
    another_way = small_molecule(smiles="c1ccccc1")
    assert one_way.smiles == "C1=CC=CC=C1"
    assert one_way.identity_key() == another_way.identity_key()


# --- peptides ----------------------------------------------------------------------------


@pytest.mark.parametrize("code", ["B", "Z", "J", "X"])
def test_a_peptide_sequence_holding_an_ambiguity_code_is_refused(code):
    with pytest.raises(ValidationError, match="ambiguity codes"):
        peptide(sequence=f"PEPT{code}DE")


@pytest.mark.parametrize("residue", ["U", "O"])
def test_selenocysteine_and_pyrrolysine_are_accepted_as_real_residues(residue):
    assert peptide(sequence=f"PEPT{residue}DE").sequence == f"PEPT{residue}DE"


@pytest.mark.parametrize("written", ["PEPTIDE1", "PEP TIDE", "PEPTIDE-2", ""])
def test_a_peptide_sequence_that_is_not_amino_acids_is_refused(written):
    with pytest.raises(ValidationError):
        peptide(sequence=written)


def test_a_lower_case_sequence_gives_the_same_identity_as_an_upper_case_one():
    assert peptide(sequence="peptide").identity_key() == peptide(sequence="PEPTIDE").identity_key()


def test_modifications_are_held_in_a_canonical_order_so_one_peptide_has_one_key():
    one_way = peptide(modifications=("Phospho@S5", "Acetyl@K2"))
    another_way = peptide(modifications=("Acetyl@K2", "Phospho@S5"))
    assert one_way.identity_key() == another_way.identity_key()


def test_a_modification_named_twice_is_refused_rather_than_merged():
    """Two identical entries are either a transcription fault or a claim about two
    sites, and choosing between them would be inventing data."""
    with pytest.raises(ValidationError, match="more than once"):
        peptide(modifications=("Phospho@S5", "Phospho@S5"))


def test_a_modified_peptide_is_a_different_ion_from_the_unmodified_one():
    assert peptide().identity_key() != peptide(modifications=("Phospho@S5",)).identity_key()


# --- glycan composition -------------------------------------------------------------------


def test_a_composition_is_canonicalised_to_one_spelling():
    assert canonical_composition("HexNAc4dHex1Hex5") == "Hex5HexNAc4Fuc1"


@pytest.mark.parametrize(
    "written",
    ["Hex5HexNAc4Fuc1", "HexNAc4Hex5Fuc1", "Fuc1Hex5HexNAc4", "HexNAc4dHex1Hex5", " Hex5HexNAc4Fuc1 "],
)
def test_every_spelling_of_one_composition_canonicalises_alike(written):
    assert canonical_composition(written) == "Hex5HexNAc4Fuc1"
    assert parse_composition(written) == Composition(hex=5, hexnac=4, fuc=1)


def test_two_spellings_of_one_composition_give_one_glycan_identity_key():
    assert glycan(composition="HexNAc4dHex1Hex5").identity_key() == glycan(composition="Hex5HexNAc4Fuc1").identity_key()


def test_dhex_and_fuc_are_one_residue_so_a_composition_naming_both_is_refused():
    with pytest.raises(CompositionError, match="given more than once"):
        parse_composition("Hex5dHex1Fuc1")


@pytest.mark.parametrize(
    "written, problem",
    [
        ("Hex5Foo2", "unknown residue"),
        ("hex5", "unknown residue"),
        ("HEXNAC2", "unknown residue"),
        ("Hex", "is not followed by a count"),
        ("Hex5HexNAc", "is not followed by a count"),
        ("Hex5-HexNAc2", "not a residue name"),
        ("Hex05", "leading zero"),
        ("Hex2Hex3", "given more than once"),
        ("Hex0", "at least one residue"),
        ("", "empty composition"),
        ("   ", "empty composition"),
    ],
)
def test_every_malformed_composition_is_refused_with_a_reason(written, problem):
    with pytest.raises(CompositionError, match=problem):
        parse_composition(written)


def test_a_malformed_composition_is_refused_by_the_analyte_too():
    with pytest.raises(ValidationError, match="unknown residue"):
        glycan(composition="Hex5Foo2")


# --- glycan identity ------------------------------------------------------------------------


def test_a_glycan_identified_by_nothing_at_all_is_refused():
    with pytest.raises(ValidationError, match="a glycan needs a composition"):
        GlycanAnalyte(source="a test", display_name="Man5")


def test_lnh_and_lnnh_key_apart_although_one_composition_is_given_for_both():
    """The case the whole rule exists for: two milk oligosaccharides of one
    composition whose deprotonated cross sections differ by 11.4 per cent."""
    lnh = glycan(composition=LNH_COMPOSITION, iupac_condensed=LNH_IUPAC)
    lnnh = glycan(composition=LNH_COMPOSITION, iupac_condensed=LNNH_IUPAC)
    assert lnh.identity_key() != lnnh.identity_key()


def test_lnh_and_lnnh_are_different_matched_ions_under_one_adduct():
    lnh = measurement(analyte=glycan(composition=LNH_COMPOSITION, iupac_condensed=LNH_IUPAC))
    lnnh = measurement(analyte=glycan(composition=LNH_COMPOSITION, iupac_condensed=LNNH_IUPAC))
    assert lnh.matched_ion_key != lnnh.matched_ion_key


def test_a_glycan_keys_on_its_wurcs_rather_than_its_composition():
    one = glycan(composition="Hex5HexNAc2", wurcs=WURCS_ONE)
    other = glycan(composition="Hex5HexNAc2", wurcs=WURCS_TWO)
    assert one.identity_key() != other.identity_key()
    assert WURCS_ONE in one.identity_key()


def test_a_glycan_keys_on_its_glytoucan_accession_rather_than_its_composition():
    one = glycan(composition="Hex5HexNAc2", glytoucan_ac="G12345AB")
    other = glycan(composition="Hex5HexNAc2", glytoucan_ac="G54321BA")
    assert one.identity_key() != other.identity_key()


def test_a_glycan_with_only_a_composition_keys_on_that_composition():
    assert glycan(composition="Hex5HexNAc2").identity_key() == glycan(composition="HexNAc2Hex5").identity_key()


@pytest.mark.parametrize("claim", ["has_unresolved_linkage", "has_unresolved_anomericity"])
@pytest.mark.parametrize(
    "identifiers",
    [
        pytest.param({"composition": "Hex5HexNAc2"}, id="composition_only"),
        pytest.param({"composition": "Hex5HexNAc2", "glytoucan_ac": "not-an-accession"}, id="malformed_accession"),
        pytest.param({"composition": "Hex5HexNAc2", "iupac_condensed": "ManManGlcNAc"}, id="no_linkage_stated"),
    ],
)
def test_a_resolved_linkage_claim_needs_a_real_structure_identifier(claim, identifiers):
    with pytest.raises(ValidationError, match="cannot claim resolved linkage or anomericity"):
        glycan(**identifiers, **{claim: False})


@pytest.mark.parametrize(
    "identifiers",
    [
        pytest.param({"wurcs": WURCS_ONE}, id="wurcs"),
        pytest.param({"glytoucan_ac": "G12345AB"}, id="well_formed_accession"),
        pytest.param({"iupac_condensed": "Man(a1-3)Man(b1-4)GlcNAc"}, id="linkage_stated"),
        pytest.param({"iupac_condensed": LNH_IUPAC}, id="lnh"),
    ],
)
def test_a_record_with_a_structure_identifier_may_claim_resolved_linkage(identifiers):
    resolved = glycan(composition=None, has_unresolved_linkage=False, **identifiers)
    assert resolved.has_unresolved_linkage is False


# Was an xfail against a real defect: a glycan stating only an accession of an
# unverified form, or a structure naming no linkage, validated and then raised
# AttributeError the moment anything asked for its key. identity_key now keys on
# what such a record actually states, tagged with why it is weak.
@pytest.mark.parametrize(
    "identifiers, expected",
    [
        pytest.param(
            {"glytoucan_ac": "GXX"},
            ("glycan", "glytoucan_unverified_form", "GXX"),
            id="accession_of_unverified_form",
        ),
        pytest.param(
            {"iupac_condensed": "ManManGlcNAc"},
            ("glycan", "iupac_no_linkage_stated", "ManManGlcNAc"),
            id="structure_with_no_linkage_stated",
        ),
    ],
)
def test_a_glycan_that_can_be_built_always_has_a_matched_ion_key(identifiers, expected):
    # The key must name WHAT THE RECORD ACTUALLY STATES, not merely be non-empty.
    # Asserting only that a key exists lets a key built from the wrong field pass:
    # a record identified by an accession would key on a null structure string,
    # which is a key every such record would share.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UnverifiedFormatWarning)
        keyless = glycan(composition=None, **identifiers)
    assert keyless.identity_key() == expected
    # And it is tagged, so it can never be read as a resolved structure.
    assert "unverified" in expected[1] or "no_linkage" in expected[1]


def test_two_glycans_identified_only_by_different_weak_identifiers_still_key_apart():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UnverifiedFormatWarning)
        one = glycan(composition=None, glytoucan_ac="GXX")
        other = glycan(composition=None, glytoucan_ac="GYY")
    assert one.identity_key() != other.identity_key()


def test_a_glytoucan_accession_of_an_unexpected_form_warns_rather_than_being_refused():
    """The form was inferred from accessions seen in use, not from documentation,
    so an unverified rule may warn but may not reject a real record."""
    with pytest.warns(UnverifiedFormatWarning, match="does not match the expected form"):
        odd = glycan(composition="Hex5HexNAc2", glytoucan_ac="GXX")
    assert odd.glytoucan_ac == "GXX"
    assert any("does not match the expected form" in note for note in odd.validation_warnings)


def test_a_well_formed_glytoucan_accession_warns_about_nothing():
    assert glycan(composition="Hex5HexNAc2", glytoucan_ac="G12345AB").validation_warnings == []


# --- what was done to a glycan before it was measured --------------------------------------


def test_a_glycans_structural_state_carries_both_the_label_and_the_derivatisation():
    labelled = glycan(reducing_end_label=ReducingEndLabel.TWO_AB, derivatisation=Derivatisation.UNDERIVATISED)
    permethylated = glycan(reducing_end_label=ReducingEndLabel.TWO_AB, derivatisation=Derivatisation.PERMETHYLATION)
    assert labelled.structural_state() != permethylated.structural_state()


def test_two_labels_of_one_glycan_are_two_matched_ions():
    native_label = measurement(analyte=glycan(reducing_end_label=ReducingEndLabel.NATIVE))
    two_ab = measurement(analyte=glycan(reducing_end_label=ReducingEndLabel.TWO_AB))
    assert native_label.matched_ion_key != two_ab.matched_ion_key


def test_two_derivatisations_of_one_glycan_are_two_matched_ions():
    underivatised = measurement(analyte=glycan(derivatisation=Derivatisation.UNDERIVATISED))
    permethylated = measurement(analyte=glycan(derivatisation=Derivatisation.PERMETHYLATION))
    assert underivatised.matched_ion_key != permethylated.matched_ion_key


def test_an_open_reducing_end_label_blocks_training_as_an_open_bucket():
    """'other' blocks for its own reason: it would pool unrelated label chemistries
    under one key, and the fix is to add the label to the enum, not to record 'other'."""
    blockers = glycan(
        reducing_end_label=ReducingEndLabel.OTHER, derivatisation=Derivatisation.UNDERIVATISED
    ).training_blockers()
    assert any("open bucket" in blocker for blocker in blockers)


def test_an_unknown_reducing_end_label_blocks_training():
    blockers = glycan(
        reducing_end_label=ReducingEndLabel.UNKNOWN, derivatisation=Derivatisation.UNDERIVATISED
    ).training_blockers()
    assert any("not defined" in blocker for blocker in blockers)


def test_an_unknown_derivatisation_blocks_training():
    blockers = glycan(
        reducing_end_label=ReducingEndLabel.NATIVE, derivatisation=Derivatisation.UNKNOWN
    ).training_blockers()
    assert any("not defined" in blocker for blocker in blockers)


@pytest.mark.parametrize("label", sorted(set(ReducingEndLabel) - {ReducingEndLabel.OTHER, ReducingEndLabel.UNKNOWN}))
def test_a_stated_label_with_a_stated_derivatisation_leaves_no_blocker(label):
    assert glycan(reducing_end_label=label, derivatisation=Derivatisation.UNDERIVATISED).training_blockers() == []


@pytest.mark.parametrize(
    "derivatisation",
    [Derivatisation.SIALIC_ACID_AMIDATION, Derivatisation.SIALIC_ACID_ESTERIFICATION],
)
def test_a_sialic_acid_derivatisation_needs_a_sialic_acid_in_the_composition(derivatisation):
    with pytest.raises(ValidationError, match="needs a sialic acid to act on"):
        glycan(composition="Hex5HexNAc2", derivatisation=derivatisation)


@pytest.mark.parametrize(
    "derivatisation",
    [Derivatisation.SIALIC_ACID_AMIDATION, Derivatisation.SIALIC_ACID_ESTERIFICATION],
)
def test_a_sialic_acid_derivatisation_is_refused_where_no_composition_can_check_it(derivatisation):
    """Unchecked is not the same as checked and passed: a claim naming a residue
    must be made against a record that says whether the residue is there."""
    with pytest.raises(ValidationError, match="gives no composition"):
        glycan(composition=None, iupac_condensed=LNH_IUPAC, derivatisation=derivatisation)


def test_a_sialic_acid_derivatisation_on_a_sialylated_glycan_is_stored_but_blocks_training():
    sialylated = glycan(
        composition="Hex5HexNAc4NeuAc2",
        reducing_end_label=ReducingEndLabel.NATIVE,
        derivatisation=Derivatisation.SIALIC_ACID_AMIDATION,
    )
    assert sialylated.derivatisation is Derivatisation.SIALIC_ACID_AMIDATION
    assert any("open bucket" in blocker for blocker in sialylated.training_blockers())


def test_a_neugc_composition_also_satisfies_the_sialic_acid_rule():
    sialylated = glycan(composition="Hex5HexNAc4NeuGc1", derivatisation=Derivatisation.SIALIC_ACID_ESTERIFICATION)
    assert sialylated.composition.neugc == 1


# --- proteins and their folding state -------------------------------------------------------


def test_a_protein_identified_only_by_a_display_name_is_refused():
    with pytest.raises(ValidationError, match="needs an accession or a sequence"):
        ProteinAnalyte(source="a test", display_name="the light chain")


def test_a_protein_may_be_identified_by_a_sequence_alone():
    assert protein(accession=None, sequence="PEPTIDE").identity_key() != protein().identity_key()


def test_a_subunit_is_part_of_a_proteins_identity():
    """A light chain and a heavy chain of one accession are not one ion."""
    light = protein(subunit="light chain")
    heavy = protein(subunit="heavy chain")
    assert light.identity_key() != heavy.identity_key()


@pytest.mark.parametrize(
    "build",
    [
        pytest.param(protein, id="protein"),
        pytest.param(antibody, id="intact_antibody"),
        pytest.param(adc, id="adc"),
    ],
)
def test_an_unstated_folding_state_blocks_training(build):
    blockers = build(folding_state=FoldingState.UNSTATED).training_blockers()
    assert any("folding_state" in blocker for blocker in blockers)


@pytest.mark.parametrize(
    "build",
    [
        pytest.param(protein, id="protein"),
        pytest.param(antibody, id="intact_antibody"),
        pytest.param(adc, id="adc"),
    ],
)
@pytest.mark.parametrize("state", [FoldingState.NATIVE, FoldingState.DENATURED])
def test_a_stated_folding_state_leaves_no_folding_blocker(build, state):
    assert build(folding_state=state).training_blockers() == []


def test_an_unstated_folding_state_keys_apart_from_both_native_and_denatured():
    """UNSTATED is a positive record of what the source does not say, so it may
    not quietly join either real state."""
    unstated = protein(folding_state=FoldingState.UNSTATED).structural_state()
    native = protein(folding_state=FoldingState.NATIVE).structural_state()
    denatured = protein(folding_state=FoldingState.DENATURED).structural_state()
    assert len({unstated, native, denatured}) == 3


def test_a_native_and_a_denatured_protein_ion_of_one_accession_are_two_matched_ions():
    native = measurement(analyte=protein(folding_state=FoldingState.NATIVE))
    denatured = measurement(analyte=protein(folding_state=FoldingState.DENATURED))
    assert native.matched_ion_key != denatured.matched_ion_key


# --- antibodies ------------------------------------------------------------------------------


def test_an_antibody_identified_by_nothing_at_all_is_refused():
    with pytest.raises(ValidationError, match="needs an INN, an accession or a sequence"):
        AntibodyIdentity()


@pytest.mark.parametrize("written", ["Trastuzumab", "TRASTUZUMAB", "trastuzumab"])
def test_two_spellings_of_one_inn_give_one_antibody_key(written):
    assert AntibodyIdentity(inn=written).key() == AntibodyIdentity(inn="trastuzumab").key()


def test_two_spellings_of_one_inn_give_one_matched_ion_key():
    one = measurement(
        analyte=antibody(antibody=antibody_identity(inn="Trastuzumab"), folding_state=FoldingState.NATIVE),
        adduct="[M+24H]24+",
        charge=24,
        ccs=7000.0,
    )
    other = measurement(
        analyte=antibody(antibody=antibody_identity(inn="trastuzumab"), folding_state=FoldingState.NATIVE),
        adduct="[M+24H]24+",
        charge=24,
        ccs=7000.0,
    )
    assert one.matched_ion_key == other.matched_ion_key


def test_an_antibodys_glycoform_is_part_of_its_identity():
    """G0F/G0F and G2F/G2F differ by four hexoses, so they are different masses
    and different cross sections."""
    g0f = antibody(glycoform="G0F/G0F")
    g2f = antibody(glycoform="G2F/G2F")
    assert g0f.identity_key() != g2f.identity_key()


def test_an_unstated_glycoform_keys_apart_from_a_stated_one():
    assert antibody(glycoform=None).identity_key() != antibody(glycoform="G0F/G0F").identity_key()


def test_an_antibodys_ciu_state_is_part_of_its_structural_state():
    """A collision-induced-unfolding step is not a replicate of the compact form."""
    compact = antibody(folding_state=FoldingState.NATIVE, ciu_state=None)
    unfolded = antibody(folding_state=FoldingState.NATIVE, ciu_state="CIU state 2")
    assert compact.structural_state() != unfolded.structural_state()


def test_two_ciu_states_of_one_antibody_ion_are_two_matched_ions():
    compact = measurement(
        analyte=antibody(folding_state=FoldingState.NATIVE, ciu_state="CIU state 1"),
        adduct="[M+24H]24+",
        charge=24,
        ccs=7000.0,
    )
    unfolded = measurement(
        analyte=antibody(folding_state=FoldingState.NATIVE, ciu_state="CIU state 2"),
        adduct="[M+24H]24+",
        charge=24,
        ccs=8000.0,
    )
    assert compact.matched_ion_key != unfolded.matched_ion_key


# --- ADCs -------------------------------------------------------------------------------------


def test_an_adc_whose_loading_is_unrecorded_is_refused():
    with pytest.raises(ValidationError, match="needs a resolved DAR or a named conjugation state"):
        adc(dar=None, conjugation_state=None)


def test_an_adc_may_state_a_conjugation_state_where_the_dar_is_unresolved():
    assert adc(dar=None, conjugation_state="unconjugated").identity_key() != adc(dar=0).identity_key()


def test_an_adcs_payload_class_is_part_of_its_identity():
    """One antibody at one DAR carrying two different payloads is two molecules."""
    vc_mmae = adc(linker_payload_class="vc-MMAE")
    mc_mmaf = adc(linker_payload_class="mc-MMAF")
    assert vc_mmae.identity_key() != mc_mmaf.identity_key()


def test_an_adcs_dar_is_part_of_its_identity():
    assert adc(dar=2).identity_key() != adc(dar=4).identity_key()


def test_a_negative_dar_is_refused():
    with pytest.raises(ValidationError):
        adc(dar=-1)


def test_two_payload_classes_of_one_antibody_are_two_matched_ions():
    vc_mmae = measurement(
        analyte=adc(linker_payload_class="vc-MMAE", folding_state=FoldingState.NATIVE),
        adduct="[M+24H]24+",
        charge=24,
        ccs=7200.0,
    )
    mc_mmaf = measurement(
        analyte=adc(linker_payload_class="mc-MMAF", folding_state=FoldingState.NATIVE),
        adduct="[M+24H]24+",
        charge=24,
        ccs=7200.0,
    )
    assert vc_mmae.matched_ion_key != mc_mmaf.matched_ion_key


# --- the union itself ---------------------------------------------------------------------------


def test_a_dumped_protein_comes_back_as_a_protein():
    adapter = TypeAdapter(Analyte)
    restored = adapter.validate_python(protein(sequence="PEPTIDE").model_dump())
    assert isinstance(restored, ProteinAnalyte)


def test_an_analyte_payload_with_no_kind_tag_is_refused_rather_than_read_as_another_kind():
    """A protein described only by a sequence is also a valid peptide payload. Without
    an explicit discriminator the union picks an arm silently, and a protein can come
    back a peptide while the two dumps still compare equal."""
    adapter = TypeAdapter(Analyte)
    ambiguous = {"sequence": "PEPTIDE", "source": "a test"}
    with pytest.raises(ValidationError):
        adapter.validate_python(ambiguous)


@pytest.mark.parametrize(
    "build, expected",
    [
        pytest.param(small_molecule, SmallMoleculeAnalyte, id="small_molecule"),
        pytest.param(peptide, PeptideAnalyte, id="peptide"),
        pytest.param(glycan, GlycanAnalyte, id="glycan"),
        pytest.param(protein, ProteinAnalyte, id="protein"),
    ],
)
def test_every_analyte_kind_survives_a_dump_and_revalidation_as_its_own_kind(build, expected):
    adapter = TypeAdapter(Analyte)
    assert isinstance(adapter.validate_python(build().model_dump()), expected)


def test_two_kinds_of_analyte_stating_the_same_text_are_not_the_same_ion():
    """A peptide sequence and a protein sequence of the same letters are different
    claims about what was measured, so the kind is in the key."""
    assert peptide(sequence="PEPTIDE").identity_key() != protein(accession=None, sequence="PEPTIDE").identity_key()


# --- a compound its dataset names and nothing else identifies ------------------------------
#
# Real published CCS tables very often carry no InChIKey. The steroid interplatform
# study gives a compound name, a commercial name, a formula and an m/z, and nothing
# that identifies a structure. The InChIKey was REQUIRED until that data arrived,
# and keeping it required would have meant either refusing the only cross-platform
# dataset this project has, or resolving 87 names against a structure database -
# a lookup that fails silently and wrongly on exactly the compounds that matter
# here, which are isomers with similar names.
#
# So a dataset-scoped id is accepted as a fallback, and the tests below are about
# keeping it NARROW. It must pair inside its own source and nowhere else.


def test_a_small_molecule_identified_by_nothing_at_all_is_refused():
    """Loosening the InChIKey requirement must not loosen it to a display name.

    A display name is not an identity: two laboratories spell one compound two ways,
    and matching on the spelling would pair ions that are not the same ion.
    """
    with pytest.raises(ValidationError, match="dataset-scoped compound id"):
        small_molecule(inchikey=None, display_name="trenbolone")


def test_a_dataset_scoped_id_identifies_a_compound_when_no_inchikey_exists():
    analyte = small_molecule(inchikey=None, dataset_compound_id="steroid_jasms2022:trenbolone")
    assert analyte.identity_key() == ("small_molecule", "dataset_compound", "steroid_jasms2022:trenbolone")
    assert analyte.identity_atoms() == frozenset({"dataset_compound:steroid_jasms2022:trenbolone"})


def test_two_rows_of_one_dataset_naming_one_compound_are_the_same_compound():
    """Within a dataset the id IS a real identity, on that table's own authority.

    That authority is a better warrant than any inference this code could make, and
    it is what lets 142 ions pair across four platforms.
    """
    one = small_molecule(inchikey=None, dataset_compound_id="steroid_jasms2022:trenbolone")
    two = small_molecule(inchikey=None, dataset_compound_id="steroid_jasms2022:trenbolone")
    assert one.identity_key() == two.identity_key()


def test_the_same_compound_name_in_two_datasets_does_not_pair():
    """The namespace is the point, and this is correct rather than unfortunate.

    Nobody has established that these two rows are the same compound. A name is
    what this platform refuses to match on, and two tables both listing
    "trenbolone" are two tables using a name.
    """
    here = small_molecule(inchikey=None, dataset_compound_id="steroid_jasms2022:trenbolone")
    there = small_molecule(inchikey=None, dataset_compound_id="ccsbase:trenbolone")
    assert here.identity_key() != there.identity_key()
    assert not (here.identity_atoms() & there.identity_atoms())


def test_an_inchikey_wins_over_a_dataset_id_so_a_resolved_compound_pairs_across_sources():
    """The fallback is a fallback. Resolving a dataset id is the act that makes a
    compound pair across sources, and once resolved the structure is the key."""
    resolved = small_molecule(inchikey=INCHIKEY, dataset_compound_id="steroid_jasms2022:trenbolone")
    plain = small_molecule(inchikey=INCHIKEY)
    assert resolved.identity_key() == plain.identity_key()
    assert resolved.identity_key()[1] == "inchikey"
    # Both atoms are carried, so the shared-peak check sees the dataset id too and a
    # resolved record still licenses a merge with its unresolved former self.
    assert resolved.identity_atoms() == frozenset(
        {f"inchikey:{INCHIKEY}", "dataset_compound:steroid_jasms2022:trenbolone"}
    )


def test_a_dataset_scoped_id_is_an_identity_atom_the_shared_peak_check_can_see():
    """Asserted here as well as in the loader, because the two surfaces must agree.

    `identity_atoms` and `identity_key` drifting apart has already happened once in
    this repository - ProteinAnalyte dropped the subunit from its atoms while keeping
    it in its key, so two records could share an atom, licensing a merge, while their
    keys differed. A new identity atom is exactly when that can happen again.
    """
    analyte = small_molecule(inchikey=None, dataset_compound_id="steroid_jasms2022:trenbolone")
    assert analyte.identity_atoms()
    for atom in analyte.identity_atoms():
        assert str(analyte.identity_key()[2]) in atom


# --- the namespace on a dataset-scoped id, which IS the safety property ---------------------
#
# Added 19 September 2026 after finding it documented and not enforced. The class
# docstring explained at length that the namespace is what stops two tables listing one
# compound name from pairing - and a bare name was accepted, so two bare names from two
# adapters matched each other. Name bridging, through the field built to prevent it,
# with nothing failing anywhere.
#
# The owner instruction is explicit: do not bridge on compound names. These tests are
# that instruction made structural rather than written down.


def test_a_bare_compound_name_is_refused_as_a_dataset_scoped_id():
    """The hole that existed. A bare name matches every other bare name."""
    with pytest.raises(ValidationError, match="has no namespace"):
        small_molecule(inchikey=None, dataset_compound_id="trenbolone")


@pytest.mark.parametrize(
    "bad",
    ["Trenbolone:17-beta", "PC(18:1/16:0)", "X:y", "17:beta", "CCSBASE:trenbolone"],
)
def test_a_namespace_that_is_not_a_dataset_tag_is_refused(bad):
    """So a compound name containing a colon cannot pose as a namespaced id.

    Real compound names do carry colons - the lipid nomenclature in CCSbase is full of
    them, PC(18:1/16:0) and worse - so "contains a colon" is not enough to make an id
    namespaced. An uppercase tag is refused too, because two spellings of one namespace
    would split one dataset into two.
    """
    with pytest.raises(ValidationError):
        small_molecule(inchikey=None, dataset_compound_id=bad)


def test_a_namespace_with_no_compound_after_it_is_refused():
    with pytest.raises(ValidationError, match="names a dataset and no compound"):
        small_molecule(inchikey=None, dataset_compound_id="ccsbase:")


def test_the_real_dataset_tags_are_accepted():
    """The three sources this repository has actually met. Not a hypothetical pattern."""
    for good in (
        "steroid_jasms2022:4-androstene-17-methyl-17-ol-3-one",
        "ccsbase:CCSBASE_A4F2E9AA6E",
        "bushlab:melittin",
    ):
        analyte = small_molecule(inchikey=None, dataset_compound_id=good)
        assert analyte.identity_key()[2] == good


def test_a_compound_name_containing_a_colon_survives_inside_the_identifier():
    """The split is on the FIRST colon, so a lipid name keeps its own."""
    analyte = small_molecule(inchikey=None, dataset_compound_id="ccsbase:PC(18:1/16:0)")
    assert analyte.identity_key()[2] == "ccsbase:PC(18:1/16:0)"


def test_one_compound_name_in_two_datasets_still_does_not_pair():
    """The whole point, asserted now that the namespace is enforced.

    Why it matters concretely: "androstenedione" names both 4-androstene-3,17-dione and
    5-androstene-3,17-dione, and betamethasone and dexamethasone are C16 epimers. 30 of
    the 87 steroid compounds share a commercial name with a CCSbase compound, so this is
    a live 30-wide temptation rather than a hypothetical one.
    """
    here = small_molecule(inchikey=None, dataset_compound_id="steroid_jasms2022:androstenedione")
    there = small_molecule(inchikey=None, dataset_compound_id="ccsbase:androstenedione")
    assert here.identity_key() != there.identity_key()
    assert not (here.identity_atoms() & there.identity_atoms())


def test_a_display_name_is_never_part_of_the_identity_surface():
    """Neither the key nor the atoms may carry it, so nothing can match on it.

    Both surfaces are checked because they license different things: the key decides
    what PAIRS, and a shared ATOM licenses a MERGE in the shared-peak check. A name
    reaching either one would bridge two compounds.
    """
    analyte = small_molecule(
        inchikey=None, dataset_compound_id="ccsbase:CCSBASE_1", display_name="Betamethasone"
    )
    assert "Betamethasone" not in str(analyte.identity_key())
    assert not any("Betamethasone" in atom for atom in analyte.identity_atoms())
    # and two records sharing ONLY a display name are two compounds
    other = small_molecule(
        inchikey=None,
        dataset_compound_id="steroid_jasms2022:9-fluoro-16-methylprednisolone",
        display_name="Betamethasone",
    )
    assert analyte.identity_key() != other.identity_key()
    assert not (analyte.identity_atoms() & other.identity_atoms())
