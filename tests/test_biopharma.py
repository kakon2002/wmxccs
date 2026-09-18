"""The biopharmaceutical analyte types: antibody, subunit, ADC and glycopeptide.

M1 items 2 and 3. These four kinds are the ones this platform exists to serve,
and they are the ones whose identity is easiest to get wrong, because every
distinction that matters here is invisible in a compound column:

- a native ion and a denatured ion of ONE antibody are two molecules' worth of
  cross section apart, and both are printed as "trastuzumab";
- two ADCs differing only in how much payload is attached are two masses and two
  cross sections, and both are printed as the same conjugate;
- a glycopeptide is neither its backbone nor its glycan, and a paper naming the
  glycan and the peptide separately has still reported one molecule that is
  neither of them.

So every rule below is tested in BOTH directions: the records that must key
apart do, and the records that must key together do. A file that only proved
things differ would pass just as happily over a key that made everything
unique, which finds no pairs at all and reports a clean run.

The gate regression at the top is not decoration. Listing AntibodyIdentity as a
component record made the licence gate refuse every intact antibody and every
ADC measurement whatever its licence, which closed this entire layer once
already. It is asserted for all four kinds rather than for one, because the
next component added will be added to one class and inherited by the rest.
"""

from __future__ import annotations

import pytest
from pydantic import TypeAdapter, ValidationError

from conftest import (
    FIXTURE,
    SOURCE,
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
    ADCAnalyte,
    Analyte,
    AnalyteKind,
    AntibodyIdentity,
    Derivatisation,
    FoldingState,
    GlycanAnalyte,
    GlycopeptideAnalyte,
    IntactAntibodyAnalyte,
    PeptideAnalyte,
    ProteinAnalyte,
    SmallMoleculeAnalyte,
)
from wmxccs.licensing import LicenceGateError, TrainingGateError, assert_trainable
from wmxccs.models import DTIMSMethod, IMSType
from wmxccs.reuse import ReuseStatus

# One backbone, two asparagines on it, and glycans that differ. Nothing here is
# asserted to be a real glycopeptide; they are four identifiers that must not be
# confused with one another.
BACKBONE = "NLTK"
OTHER_BACKBONE = "NKTL"
SITE = "N297"
OTHER_SITE = "N100"
NEUTRAL_GLYCAN = "Hex5HexNAc2"
OTHER_GLYCAN = "Hex5HexNAc4"
SIALYLATED_GLYCAN = "Hex5HexNAc4NeuAc2"
NEUGC_GLYCAN = "Hex5HexNAc4NeuGc1"
GLYCAN_WURCS = "WURCS=2.0/1,1,0/[a2122h-1b_1-5]/1/"
GLYCAN_IUPAC = "Gal(b1-4)GlcNAc"  # states a linkage, so it counts as a structure
GLYCAN_GLYTOUCAN = "G12345AB"  # the expected form: G, five digits, two capitals


def glycopeptide(**overrides) -> GlycopeptideAnalyte:
    """A localised, underivatised glycopeptide: the minimal record that does NOT block training.

    Underivatised rather than the field's own default of 'unknown', for the
    reason conftest's cyclic_settings gives: a test that wants a blocker has to
    ask for one, rather than getting one by accident from a default and passing
    for the wrong reason.

    Not in conftest, because no other file needs it yet.
    """
    fields = dict(
        sequence=BACKBONE,
        attachment_site=SITE,
        glycan_composition=NEUTRAL_GLYCAN,
        derivatisation=Derivatisation.UNDERIVATISED,
        source=SOURCE,
        reuse_status=FIXTURE,
    )
    fields.update(overrides)
    return GlycopeptideAnalyte(**fields)


# --- one measurement per biopharmaceutical kind -----------------------------------------
#
# Each states its folding state where it has one, and each carries an adduct
# whose charge is the charge a paper would report for that kind of ion. They
# take overrides so that a test can change ONE thing and nothing else.


def antibody_record(**overrides):
    fields = dict(
        analyte=antibody(folding_state=FoldingState.NATIVE), adduct="[M+24H]24+", charge=24, ccs=7000.0
    )
    fields.update(overrides)
    return measurement(**fields)


def adc_record(**overrides):
    fields = dict(analyte=adc(folding_state=FoldingState.NATIVE), adduct="[M+24H]24+", charge=24, ccs=7200.0)
    fields.update(overrides)
    return measurement(**fields)


def subunit_record(**overrides):
    fields = dict(
        analyte=protein(subunit="light chain", folding_state=FoldingState.DENATURED),
        adduct="[M+11H]11+",
        charge=11,
        ccs=2400.0,
    )
    fields.update(overrides)
    return measurement(**fields)


def glycopeptide_record(**overrides):
    fields = dict(analyte=glycopeptide(), adduct="[M+3H]3+", charge=3, ccs=600.0)
    fields.update(overrides)
    return measurement(**fields)


BIOPHARMA_RECORDS = [
    pytest.param(antibody_record, id="intact_antibody"),
    pytest.param(adc_record, id="adc"),
    pytest.param(subunit_record, id="protein_subunit"),
    pytest.param(glycopeptide_record, id="glycopeptide"),
]

FOLDED_RECORDS = [
    pytest.param(antibody_record, id="intact_antibody"),
    pytest.param(adc_record, id="adc"),
    pytest.param(subunit_record, id="protein_subunit"),
]


def as_stepped_field_dtims(build, **overrides):
    """The same ion on the one primary platform. Nothing about the ION changes."""
    return build(ims_type=IMSType.DTIMS, dtims_method=DTIMSMethod.STEPPED_FIELD, calibrant=None, **overrides)


# --- (a) one antibody, two ions ---------------------------------------------------------


def test_a_native_24_plus_and_a_denatured_40_plus_ion_of_one_antibody_are_different_matched_ions():
    """The assertion the brief names. Do not delete: pooling these two would report
    the difference between a folded and an unfolded molecule as inter-platform bias."""
    native_24 = antibody_record()
    denatured_40 = antibody_record(
        analyte=antibody(folding_state=FoldingState.DENATURED), adduct="[M+40H]40+", charge=40, ccs=11000.0
    )
    assert native_24.matched_ion_key != denatured_40.matched_ion_key


def test_a_native_24_plus_and_a_native_40_plus_antibody_ion_differ_by_charge_alone():
    """The charge half of the assertion above, which would otherwise pass on the
    folding state alone. The adduct changes with the charge because it has to: an
    adduct whose charge contradicts the charge field is refused."""
    native_24 = antibody_record()
    native_40 = antibody_record(adduct="[M+40H]40+", charge=40)
    assert native_24.matched_ion_key != native_40.matched_ion_key


def test_a_native_40_plus_and_a_denatured_40_plus_antibody_ion_differ_by_folding_state_alone():
    """The folding half: one charge, one adduct, one gas, one glycoform, two ions."""
    native_40 = antibody_record(adduct="[M+40H]40+", charge=40)
    denatured_40 = antibody_record(
        analyte=antibody(folding_state=FoldingState.DENATURED), adduct="[M+40H]40+", charge=40
    )
    assert native_40.matched_ion_key != denatured_40.matched_ion_key


# --- (b) two ADCs, one difference -------------------------------------------------------


def test_two_adcs_differing_only_in_dar_are_different_matched_ions():
    """The assertion the brief names. One antibody, one payload class, one glycoform,
    one charge state, one gas: only the drug load differs, and that is two molecules."""
    dar_two = adc_record(analyte=adc(folding_state=FoldingState.NATIVE, dar=2))
    dar_four = adc_record(analyte=adc(folding_state=FoldingState.NATIVE, dar=4))
    assert dar_two.matched_ion_key != dar_four.matched_ion_key


def test_two_adcs_differing_only_in_dar_have_different_identity_keys():
    assert adc(dar=2).identity_key() != adc(dar=4).identity_key()


def test_two_adcs_differing_only_in_dar_share_no_identity_atom():
    """M2 merges two analytes that share an atom into one molecule. A DAR 2 and a
    DAR 4 conjugate must not be mergeable by any route."""
    assert not (adc(dar=2).identity_atoms() & adc(dar=4).identity_atoms())


def test_two_adcs_of_one_dar_are_the_same_matched_ion():
    """The other direction. Without this, a key that made every record unique would
    pass every test above and find no pairs at all."""
    assert adc_record().matched_ion_key == adc_record().matched_ion_key


# --- (c) a glycopeptide is neither of its halves -----------------------------------------


def test_a_glycopeptide_and_its_bare_peptide_backbone_are_different_matched_ions():
    """The assertion the brief names. Same sequence, same adduct, same charge, same
    gas; one carries a glycan and the other does not, and their masses differ by it."""
    glycosylated = glycopeptide_record()
    backbone = glycopeptide_record(analyte=peptide(sequence=BACKBONE), ccs=480.0)
    assert glycosylated.matched_ion_key != backbone.matched_ion_key


def test_a_glycopeptide_and_its_bare_peptide_backbone_share_no_identity_atom():
    """Sharing one atom would let matched-ion construction merge the two as one
    molecule in M2, whatever their keys say."""
    assert not (glycopeptide().identity_atoms() & peptide(sequence=BACKBONE).identity_atoms())


def test_a_glycopeptide_and_the_free_glycan_of_its_composition_are_different_matched_ions():
    glycosylated = glycopeptide_record()
    free = glycopeptide_record(
        analyte=glycan(composition=NEUTRAL_GLYCAN, derivatisation=Derivatisation.UNDERIVATISED), ccs=280.0
    )
    assert glycosylated.matched_ion_key != free.matched_ion_key


def test_a_glycopeptide_and_the_free_glycan_of_its_composition_share_no_identity_atom():
    assert not (glycopeptide().identity_atoms() & glycan(composition=NEUTRAL_GLYCAN).identity_atoms())


def test_every_identity_atom_of_a_glycopeptide_is_namespaced_to_the_glycopeptide():
    """The seam M2 merges on. An atom of a glycopeptide must not be spelt like a
    peptide's or a glycan's, or the merge step would join a glycopeptide to a
    molecule it is not, on a string comparison nobody would look at twice."""
    atoms = glycopeptide(glycan_wurcs=GLYCAN_WURCS, glycan_glytoucan_ac=GLYCAN_GLYTOUCAN).identity_atoms()
    assert atoms
    assert all(atom.startswith("glycopeptide:") for atom in atoms)


def test_two_glycopeptides_differing_only_in_attachment_site_share_no_identity_atom():
    """The site is part of the molecule, so it must be part of every atom as well as
    of the key. One glycan on Asn297 and the same glycan on another asparagine of the
    same backbone are different shapes, and that difference is most of this field."""
    here = glycopeptide(attachment_site=SITE)
    there = glycopeptide(attachment_site=OTHER_SITE)
    assert not (here.identity_atoms() & there.identity_atoms())


# --- M1 item 3: the gate regression, for every biopharmaceutical kind --------------------


@pytest.mark.parametrize("build", BIOPHARMA_RECORDS)
def test_a_valid_biopharma_record_clears_the_training_gate(build):
    """The M0 fix, asserted for all four kinds. AntibodyIdentity was listed as a
    component record, the gate refused every component that cannot state a reuse
    status, and so every antibody and every ADC was refused whatever its licence."""
    assert assert_trainable(build()) is None


@pytest.mark.parametrize("build", BIOPHARMA_RECORDS)
def test_the_gate_still_refuses_a_biopharma_record_whose_licence_is_unverified(build):
    """The other direction, so the four tests above cannot pass on a gate that
    stopped refusing anything at all."""
    with pytest.raises(LicenceGateError, match="nobody has checked its terms yet"):
        assert_trainable(build(reuse_status=ReuseStatus.UNVERIFIED))


def test_the_antibody_identity_is_not_a_component_record_of_the_measurement():
    """Why the four records above clear: the antibody identity is the structured form
    of the analyte's own identity, not a separately sourced record, and it carries no
    reuse status for the gate to read. Listing it is what closed this layer."""
    record = adc_record()
    identity = record.analyte.antibody
    assert not hasattr(identity, "reuse_status")
    assert all(component is not identity for component in record.component_records())


@pytest.mark.parametrize("build", FOLDED_RECORDS)
def test_a_biopharma_record_that_does_not_state_its_folding_state_is_still_refused(build):
    """A separate guard, and a correct one: an ion of unrecorded conformation cannot
    be compared with one whose conformation is recorded. The gate fix must not have
    relaxed it."""
    analyte = build().analyte.model_copy(update={"folding_state": FoldingState.UNSTATED})
    unstated = build(analyte=type(analyte).model_validate(analyte.model_dump()))
    with pytest.raises(TrainingGateError, match="folding_state is 'UNSTATED'"):
        assert_trainable(unstated)


def test_a_glycopeptide_whose_derivatisation_is_unrecorded_is_refused():
    """The glycopeptide's equivalent guard, for the same reason."""
    with pytest.raises(TrainingGateError, match="derivatisation is 'unknown'"):
        assert_trainable(glycopeptide_record(analyte=glycopeptide(derivatisation=Derivatisation.UNKNOWN)))


# --- the key carries no platform, for these kinds too ------------------------------------


@pytest.mark.parametrize("build", BIOPHARMA_RECORDS)
def test_one_biopharma_ion_measured_on_two_platforms_lands_on_one_matched_ion_key(build):
    """The whole point of the key. If a biopharmaceutical record keyed by platform
    there would be nothing to compare and the pipeline would report zero pairs over
    a corpus full of them."""
    travelling_wave = build()
    stepped_field = as_stepped_field_dtims(build, ccs=build().ccs * 1.01)
    assert travelling_wave.matched_ion_key == stepped_field.matched_ion_key


@pytest.mark.parametrize("build", BIOPHARMA_RECORDS)
def test_the_calibration_group_of_those_two_platforms_is_not_one_group(build):
    """The complement of the test above, and the reason the two keys exist: what may
    be compared is not what may be pooled."""
    assert build().calibration_group != as_stepped_field_dtims(build).calibration_group


# --- ADC identity ------------------------------------------------------------------------


def test_an_adc_needs_a_resolved_dar_or_a_named_conjugation_state():
    with pytest.raises(ValidationError, match="needs a resolved DAR or a named conjugation state"):
        adc(dar=None, conjugation_state=None)


@pytest.mark.parametrize(
    "loading",
    [
        pytest.param(dict(dar=4, conjugation_state=None), id="resolved_dar"),
        pytest.param(dict(dar=0, conjugation_state=None), id="dar_zero_is_a_species"),
        pytest.param(dict(dar=None, conjugation_state="unconjugated"), id="named_conjugation_state"),
    ],
)
def test_an_adc_that_states_its_loading_either_way_can_be_built(loading):
    assert adc(**loading).identity_key()


def test_a_resolved_dar_of_zero_and_a_named_unconjugated_species_are_not_one_ion():
    """A record saying "DAR 0" claims the loading was resolved and came out zero. A
    record saying "unconjugated" says what the paper called the peak. They are not
    the same claim, so they are not merged into one."""
    assert adc(dar=0).identity_key() != adc(dar=None, conjugation_state="unconjugated").identity_key()


def test_an_average_dar_is_refused_because_a_mixture_is_not_an_ion():
    with pytest.raises(ValidationError):
        adc(dar=3.5)


def test_an_adc_without_a_linker_payload_class_cannot_be_built():
    """Required, not optional: two conjugates of one antibody at one DAR are
    different molecules if the payload differs."""
    with pytest.raises(ValidationError, match="linker_payload_class"):
        ADCAnalyte(antibody=antibody_identity(), dar=2, source=SOURCE, reuse_status=FIXTURE)


def test_two_adcs_differing_only_in_payload_class_have_different_identity_keys():
    assert adc(linker_payload_class="vc-MMAE").identity_key() != adc(linker_payload_class="mc-MMAF").identity_key()


def test_two_adcs_differing_only_in_payload_class_are_different_matched_ions():
    vc_mmae = adc_record(analyte=adc(folding_state=FoldingState.NATIVE, linker_payload_class="vc-MMAE"))
    mc_mmaf = adc_record(analyte=adc(folding_state=FoldingState.NATIVE, linker_payload_class="mc-MMAF"))
    assert vc_mmae.matched_ion_key != mc_mmaf.matched_ion_key


def test_two_adcs_of_one_payload_class_written_alike_are_the_same_matched_ion():
    assert adc_record().matched_ion_key == adc_record().matched_ion_key


def test_two_adcs_differing_only_in_glycoform_are_different_matched_ions():
    stated = adc_record(analyte=adc(folding_state=FoldingState.NATIVE, glycoform="G0F/G0F"))
    unstated = adc_record(analyte=adc(folding_state=FoldingState.NATIVE, glycoform=None))
    assert stated.matched_ion_key != unstated.matched_ion_key


def test_an_adcs_ciu_state_is_part_of_its_structural_state():
    compact = adc(folding_state=FoldingState.NATIVE, ciu_state=None)
    unfolded = adc(folding_state=FoldingState.NATIVE, ciu_state="CIU state 2")
    assert compact.structural_state() != unfolded.structural_state()


def test_two_adcs_of_one_antibody_differing_in_nothing_but_the_antibody_are_two_ions():
    trastuzumab = adc(antibody=antibody_identity(inn="trastuzumab"))
    rituximab = adc(antibody=antibody_identity(inn="rituximab"))
    assert trastuzumab.identity_key() != rituximab.identity_key()


# --- intact antibody identity -------------------------------------------------------------


def test_an_antibodys_glycoform_is_part_of_its_identity():
    """G0F/G0F and G2F/G2F differ by four hexoses: two masses, two cross sections."""
    assert antibody(glycoform="G0F/G0F").identity_key() != antibody(glycoform="G2F/G2F").identity_key()


def test_an_unstated_glycoform_keys_apart_from_a_stated_one():
    """Optional, because public intact-antibody data almost never resolves it, and
    keyed apart, because "not resolved" is not the same claim as "G0F/G0F"."""
    assert antibody(glycoform=None).identity_key() != antibody(glycoform="G0F/G0F").identity_key()


def test_two_antibodies_of_one_stated_glycoform_are_the_same_matched_ion():
    one = antibody_record(analyte=antibody(folding_state=FoldingState.NATIVE, glycoform="G0F/G0F"))
    other = antibody_record(analyte=antibody(folding_state=FoldingState.NATIVE, glycoform="G0F/G0F"))
    assert one.matched_ion_key == other.matched_ion_key


def test_an_antibodys_ciu_state_is_part_of_its_structural_state():
    compact = antibody(folding_state=FoldingState.NATIVE, ciu_state=None)
    unfolded = antibody(folding_state=FoldingState.NATIVE, ciu_state="CIU state 2")
    assert compact.structural_state() != unfolded.structural_state()


def test_two_ciu_states_of_one_antibody_ion_are_two_matched_ions():
    """An activation ramp gives several cross sections for one charge state. They are
    not replicates of one number."""
    first = antibody_record(analyte=antibody(folding_state=FoldingState.NATIVE, ciu_state="CIU state 1"))
    second = antibody_record(
        analyte=antibody(folding_state=FoldingState.NATIVE, ciu_state="CIU state 2"), ccs=8200.0
    )
    assert first.matched_ion_key != second.matched_ion_key


def test_an_antibody_needs_an_inn_an_accession_or_a_sequence():
    with pytest.raises(ValidationError, match="needs an INN, an accession or a sequence"):
        AntibodyIdentity()


def test_an_antibody_identity_will_not_even_hold_a_display_name():
    """Stronger than keeping the name out of the key: there is no field to put it in,
    so an identity cannot be built out of what a paper calls the molecule. The
    printed name belongs to the enclosing analyte."""
    with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
        AntibodyIdentity(display_name="the mAb")
    assert antibody(display_name="the mAb").identity_key() == antibody().identity_key()


def test_two_spellings_of_one_inn_give_one_matched_ion_key():
    one = antibody_record(analyte=antibody(antibody=antibody_identity(inn="Trastuzumab"), folding_state=FoldingState.NATIVE))
    other = antibody_record(analyte=antibody(antibody=antibody_identity(inn="trastuzumab"), folding_state=FoldingState.NATIVE))
    assert one.matched_ion_key == other.matched_ion_key


def test_an_intact_antibody_and_an_adc_of_one_antibody_are_not_the_same_ion():
    """A naked antibody and a conjugate of it share an INN and nothing else."""
    assert antibody(folding_state=FoldingState.NATIVE).identity_key() != adc(
        folding_state=FoldingState.NATIVE
    ).identity_key()


def test_an_intact_antibody_and_an_adc_of_one_antibody_share_no_identity_atom():
    naked = antibody(folding_state=FoldingState.NATIVE)
    conjugate = adc(folding_state=FoldingState.NATIVE)
    assert not (naked.identity_atoms() & conjugate.identity_atoms())


# --- the unstated charge carrier, where native MS actually leaves it ----------------------


def test_two_antibody_records_whose_charge_carrier_is_unstated_never_match_each_other():
    """"The 24+ ion" from two papers is not one ion: twenty-four protons and
    twenty-four ammonium adducts differ by 408 Da. The key is made unique to the
    record, so the two can never be paired, and this is asserted for an antibody
    because native MS is where the carrier is most often left unsaid."""
    one = antibody_record(adduct="[M+24?]24+", source="one paper")
    other = antibody_record(adduct="[M+24?]24+", source="another paper")
    assert one.matched_ion_key != other.matched_ion_key
    assert not one.matched_ion_key.matchable


def test_an_unstated_carrier_keys_apart_from_the_protonated_ion_of_one_antibody():
    assert antibody_record(adduct="[M+24?]24+").matched_ion_key != antibody_record().matched_ion_key


def test_an_antibody_record_whose_charge_carrier_is_unstated_is_refused_by_the_gate():
    with pytest.raises(TrainingGateError, match="does not name the charge carrier"):
        assert_trainable(antibody_record(adduct="[M+24?]24+"))


# --- protein and subunit identity ----------------------------------------------------------


def test_a_protein_identified_only_by_a_display_name_is_refused():
    """"The light chain" identifies nothing across two laboratories, and a subunit in
    particular is whatever the digest produced."""
    with pytest.raises(ValidationError, match="needs an accession or a sequence"):
        protein(accession=None, display_name="the light chain")


@pytest.mark.parametrize(
    "identifiers",
    [
        pytest.param(dict(accession="P01857", sequence=None), id="accession"),
        pytest.param(dict(accession=None, sequence="PEPTIDE"), id="sequence"),
        pytest.param(dict(accession="P01857", sequence="PEPTIDE"), id="both"),
    ],
)
def test_a_protein_stating_an_accession_or_a_sequence_can_be_built(identifiers):
    assert protein(**identifiers).identity_key()


def test_a_light_chain_and_a_heavy_chain_of_one_accession_are_different_matched_ions():
    """The subunit is part of the protein key, so one accession covering two chains
    does not collapse them into one ion."""
    light = subunit_record(
        analyte=protein(subunit="light chain", folding_state=FoldingState.DENATURED), ccs=2400.0
    )
    heavy = subunit_record(
        analyte=protein(subunit="heavy chain", folding_state=FoldingState.DENATURED), ccs=3600.0
    )
    assert light.matched_ion_key != heavy.matched_ion_key


def test_an_unstated_subunit_keys_apart_from_a_stated_one():
    """An intact protein record and a record of one of its chains are not the same
    molecule, and the intact record is the one that names no subunit."""
    assert protein(subunit=None).identity_key() != protein(subunit="light chain").identity_key()


def test_two_light_chains_of_one_accession_are_the_same_matched_ion():
    assert subunit_record().matched_ion_key == subunit_record().matched_ion_key


def test_a_native_and_a_denatured_subunit_of_one_accession_are_two_matched_ions():
    native = subunit_record(analyte=protein(subunit="light chain", folding_state=FoldingState.NATIVE))
    denatured = subunit_record(analyte=protein(subunit="light chain", folding_state=FoldingState.DENATURED))
    assert native.matched_ion_key != denatured.matched_ion_key


# --- glycopeptide identity -------------------------------------------------------------------


def test_a_glycopeptide_with_no_glycan_at_all_is_refused():
    """Without a glycan this is a peptide record claiming to be a glycopeptide, and it
    would key as a molecule nobody has identified."""
    with pytest.raises(ValidationError, match="a glycopeptide needs a glycan"):
        glycopeptide(
            glycan_composition=None,
            glycan_wurcs=None,
            glycan_glytoucan_ac=None,
            glycan_iupac_condensed=None,
        )


@pytest.mark.parametrize(
    "identifier",
    [
        pytest.param(dict(glycan_composition=NEUTRAL_GLYCAN), id="composition"),
        pytest.param(dict(glycan_wurcs=GLYCAN_WURCS), id="wurcs"),
        pytest.param(dict(glycan_glytoucan_ac=GLYCAN_GLYTOUCAN), id="glytoucan"),
        pytest.param(dict(glycan_iupac_condensed=GLYCAN_IUPAC), id="iupac"),
    ],
)
def test_a_glycopeptide_stating_any_one_glycan_identifier_can_be_built_and_has_a_key(identifier):
    """The other direction of the rule above: any one of the four is enough, and each
    of them produces a key rather than an AttributeError the first time M2 asks."""
    blank = dict(
        glycan_composition=None, glycan_wurcs=None, glycan_glytoucan_ac=None, glycan_iupac_condensed=None
    )
    assert glycopeptide(**(blank | identifier)).identity_key()


def test_one_glycan_on_two_attachment_sites_gives_two_matched_ions():
    """The site is identity, not metadata. A great deal of biopharmaceutical analysis
    is about exactly this difference."""
    here = glycopeptide_record(analyte=glycopeptide(attachment_site=SITE))
    there = glycopeptide_record(analyte=glycopeptide(attachment_site=OTHER_SITE), ccs=590.0)
    assert here.matched_ion_key != there.matched_ion_key


def test_one_glycan_on_one_site_gives_one_matched_ion():
    assert glycopeptide_record().matched_ion_key == glycopeptide_record().matched_ion_key


def test_an_unlocalised_glycopeptide_keys_apart_from_a_localised_one():
    """A glycopeptide of unknown site is not known to be the same molecule as one of
    a stated site, so it is not pooled with it."""
    assert glycopeptide(attachment_site=None).identity_key() != glycopeptide(attachment_site=SITE).identity_key()


def test_two_glycans_on_one_backbone_and_one_site_give_two_matched_ions():
    """The glycan is in the key, so one backbone does not pool every glycoform of
    itself into a single ion."""
    smaller = glycopeptide_record(analyte=glycopeptide(glycan_composition=NEUTRAL_GLYCAN))
    larger = glycopeptide_record(analyte=glycopeptide(glycan_composition=OTHER_GLYCAN), ccs=650.0)
    assert smaller.matched_ion_key != larger.matched_ion_key


def test_two_backbones_carrying_one_glycan_give_two_matched_ions():
    one = glycopeptide_record(analyte=glycopeptide(sequence=BACKBONE))
    other = glycopeptide_record(analyte=glycopeptide(sequence=OTHER_BACKBONE))
    assert one.matched_ion_key != other.matched_ion_key


def test_a_glycopeptides_derivatisation_is_part_of_its_structural_state():
    underivatised = glycopeptide(derivatisation=Derivatisation.UNDERIVATISED)
    permethylated = glycopeptide(derivatisation=Derivatisation.PERMETHYLATION)
    assert underivatised.structural_state() != permethylated.structural_state()


def test_two_derivatisations_of_one_glycopeptide_are_two_matched_ions():
    underivatised = glycopeptide_record()
    permethylated = glycopeptide_record(
        analyte=glycopeptide(derivatisation=Derivatisation.PERMETHYLATION), ccs=640.0
    )
    assert underivatised.matched_ion_key != permethylated.matched_ion_key


@pytest.mark.parametrize(
    "derivatisation",
    [Derivatisation.SIALIC_ACID_AMIDATION, Derivatisation.SIALIC_ACID_ESTERIFICATION],
)
def test_a_sialic_acid_derivatisation_on_a_glycopeptide_needs_a_sialic_acid_to_act_on(derivatisation):
    with pytest.raises(ValidationError, match="needs a sialic acid to act on"):
        glycopeptide(glycan_composition=NEUTRAL_GLYCAN, derivatisation=derivatisation)


@pytest.mark.parametrize(
    "derivatisation",
    [Derivatisation.SIALIC_ACID_AMIDATION, Derivatisation.SIALIC_ACID_ESTERIFICATION],
)
def test_a_sialic_acid_derivatisation_on_a_glycopeptide_is_refused_where_nothing_can_check_it(derivatisation):
    """A structure with no composition cannot answer whether a sialic acid is there,
    so the claim is refused rather than waved through."""
    with pytest.raises(ValidationError, match="names a sialic acid and this record gives no glycan composition"):
        glycopeptide(
            glycan_composition=None, glycan_wurcs=GLYCAN_WURCS, derivatisation=derivatisation
        )


@pytest.mark.parametrize("composition", [SIALYLATED_GLYCAN, NEUGC_GLYCAN])
def test_a_sialic_acid_derivatisation_on_a_sialylated_glycopeptide_is_stored_but_blocks_training(composition):
    """The other direction, and it stops halfway: the record is accepted, because it
    is what the source says, and it is refused for training, because the bucket
    covers chemistries of different mass."""
    analyte = glycopeptide(
        glycan_composition=composition, derivatisation=Derivatisation.SIALIC_ACID_AMIDATION
    )
    assert analyte.identity_key()
    with pytest.raises(TrainingGateError, match="an open bucket covering chemistries of different mass"):
        assert_trainable(glycopeptide_record(analyte=analyte))


def test_a_glycopeptide_stating_its_derivatisation_leaves_no_blocker():
    assert glycopeptide(derivatisation=Derivatisation.UNDERIVATISED).training_blockers() == []


# --- the union discriminates ----------------------------------------------------------------


def test_an_analyte_payload_valid_for_two_arms_is_refused_rather_than_read_as_one_of_them():
    """A subunit identified by a sequence alone is also a valid peptide payload. Left
    to a smart union the arm is picked silently, and a protein subunit can come back a
    peptide while the two dumps still compare equal, so the revalidation guard that
    compares dumps would never notice."""
    adapter = TypeAdapter(Analyte)
    ambiguous = {"sequence": "PEPTIDE", "source": SOURCE, "reuse_status": FIXTURE}
    with pytest.raises(ValidationError, match="kind_tag"):
        adapter.validate_python(ambiguous)


@pytest.mark.parametrize(
    "kind, expected",
    [
        pytest.param(AnalyteKind.PROTEIN, ProteinAnalyte, id="protein"),
        pytest.param(AnalyteKind.PEPTIDE, PeptideAnalyte, id="peptide"),
    ],
)
def test_that_same_payload_reaches_the_arm_its_kind_tag_names(kind, expected):
    """The other direction: the payload is fine, and it is the missing discriminator
    that is refused above."""
    adapter = TypeAdapter(Analyte)
    payload = {"kind_tag": kind.value, "sequence": "PEPTIDE", "source": SOURCE, "reuse_status": FIXTURE}
    assert type(adapter.validate_python(payload)) is expected


# --- every kind has a builder -----------------------------------------------------------------

ANALYTE_BUILDERS = {
    AnalyteKind.SMALL_MOLECULE: small_molecule,
    AnalyteKind.PEPTIDE: peptide,
    AnalyteKind.GLYCAN: glycan,
    AnalyteKind.GLYCOPEPTIDE: glycopeptide,
    AnalyteKind.PROTEIN: protein,
    AnalyteKind.INTACT_ANTIBODY: antibody,
    AnalyteKind.ADC: adc,
}

ANALYTE_CLASSES = {
    AnalyteKind.SMALL_MOLECULE: SmallMoleculeAnalyte,
    AnalyteKind.PEPTIDE: PeptideAnalyte,
    AnalyteKind.GLYCAN: GlycanAnalyte,
    AnalyteKind.GLYCOPEPTIDE: GlycopeptideAnalyte,
    AnalyteKind.PROTEIN: ProteinAnalyte,
    AnalyteKind.INTACT_ANTIBODY: IntactAntibodyAnalyte,
    AnalyteKind.ADC: ADCAnalyte,
}


def test_every_analyte_kind_has_a_builder_so_a_new_kind_cannot_arrive_untested():
    """Set coverage, not a count. A kind added to the enum with no builder fails here
    rather than shipping with no test behind it."""
    assert set(ANALYTE_BUILDERS) == set(AnalyteKind)
    assert set(ANALYTE_CLASSES) == set(AnalyteKind)


@pytest.mark.parametrize("kind", list(AnalyteKind), ids=[kind.value for kind in AnalyteKind])
def test_every_analyte_kind_builds_keys_and_survives_the_union_as_its_own_kind(kind):
    built = ANALYTE_BUILDERS[kind]()
    assert type(built) is ANALYTE_CLASSES[kind]
    assert built.kind is kind
    assert built.identity_key()
    assert built.identity_atoms()
    restored = TypeAdapter(Analyte).validate_python(built.model_dump())
    assert type(restored) is ANALYTE_CLASSES[kind]
    assert restored.identity_key() == built.identity_key()


@pytest.mark.parametrize("kind", list(AnalyteKind), ids=[kind.value for kind in AnalyteKind])
def test_no_two_analyte_kinds_share_an_identity_key_prefix(kind):
    """The kind leads every identity key, so two kinds stating the same text - a
    peptide sequence and a subunit sequence, a free glycan composition and a
    glycopeptide's - are never one ion."""
    assert ANALYTE_BUILDERS[kind]().identity_key()[0] == kind.value


# --- the two identity surfaces must not contradict each other ---------------------------


def test_a_proteins_sequence_atom_carries_the_subunit_like_its_key_does():
    """identity_key and identity_atoms must agree about what is one molecule.

    The key holds two subunit labels apart. An atom that dropped the label would
    say the same two records are one molecule and may be merged, and matched-ion
    construction is told that sharing an atom licenses a merge. The merge would
    then join records the key deliberately separates, and the disagreement would
    surface only as a pair nobody could explain.
    """
    from wmxccs.identity import ProteinAnalyte
    from wmxccs.reuse import ReuseStatus

    fixture = ReuseStatus.SYNTHETIC_FIXTURE
    light = ProteinAnalyte(sequence="PEPTIDE", subunit="light chain", source="a test", reuse_status=fixture)
    heavy = ProteinAnalyte(sequence="PEPTIDE", subunit="heavy chain", source="a test", reuse_status=fixture)

    assert light.identity_key() != heavy.identity_key()
    assert not (light.identity_atoms() & heavy.identity_atoms())


def test_two_proteins_of_one_sequence_and_one_subunit_do_still_share_an_atom():
    """The negative direction. Without it the test above passes on a protein whose
    atoms are unique per record, which would stop M2 ever merging anything."""
    from wmxccs.identity import ProteinAnalyte
    from wmxccs.reuse import ReuseStatus

    fixture = ReuseStatus.SYNTHETIC_FIXTURE
    one = ProteinAnalyte(sequence="PEPTIDE", subunit="light chain", source="one", reuse_status=fixture)
    other = ProteinAnalyte(sequence="PEPTIDE", subunit="light chain", source="another", reuse_status=fixture)

    assert one.identity_key() == other.identity_key()
    assert one.identity_atoms() & other.identity_atoms()
