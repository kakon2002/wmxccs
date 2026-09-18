"""What a measurement record refuses, and the two keys it has to keep apart.

THE TWO ASSERTIONS THIS FILE EXISTS FOR, because getting either backwards breaks
the platform without breaking a single test that only reads the happy path:

  the MATCHED-ION KEY carries NO platform, so a DTIMS value and a TWIMS value of
  one ion land on one key. A key carrying the platform can never pair anything,
  and the pipeline then reports a clean run of zero pairs over a corpus full of
  them.

  the CALIBRATION GROUP carries NO analyte, so two analytes measured the same way
  land in one group. A group as fine as the ion holds exactly one record, and the
  shared-peak check that compares two laboratories can then never fire. A check
  that can never fire passes every test it is given.

They are close to complements, they sit a few lines apart in the source, and
swapping them produces a platform that looks like it works. The tests below are
named so that nobody can swap them later without the names becoming false.
"""

from __future__ import annotations

from datetime import date, datetime

import pytest
from pydantic import ValidationError

from wmxccs.licensing import TrainingGateError, assert_trainable

from conftest import (
    OTHER_INCHIKEY,
    adc,
    antibody,
    glycan,
    measurement,
    peptide,
    primary,
    protein,
    small_molecule,
)
from wmxccs.identity import (
    DriftGas,
    FoldingState,
    Polarity,
    ReducingEndLabel,
    UnverifiedFormatWarning,
)
from wmxccs.models import (
    CalibrationLineage,
    CalibrationReference,
    DTIMSMethod,
    IMSType,
    MeasurementConditions,
    UncertaintyType,
)

# Every way a platform can be stated, each with the calibrant its method requires.
# These are five different ways of PRODUCING a number, and not one of them may
# change WHICH ION the number is of.
PLATFORMS = {
    "twims": dict(ims_type=IMSType.TWIMS, calibrant="dextran"),
    "tims": dict(ims_type=IMSType.TIMS, calibrant="dextran"),
    "cyclic": dict(ims_type=IMSType.CYCLIC, calibrant="dextran"),
    "dtims_single_field": dict(
        ims_type=IMSType.DTIMS, dtims_method=DTIMSMethod.SINGLE_FIELD, calibrant="dextran"
    ),
    "dtims_stepped_field": dict(
        ims_type=IMSType.DTIMS, dtims_method=DTIMSMethod.STEPPED_FIELD, calibrant=None
    ),
}
CALIBRATED_PLATFORMS = [name for name in PLATFORMS if name != "dtims_stepped_field"]
STATED_UNCERTAINTY_TYPES = [kind for kind in UncertaintyType if kind is not UncertaintyType.UNKNOWN]


def reference(**overrides) -> CalibrationReference:
    fields = dict(reference_set="dextran DTCCS ladder")
    fields.update(overrides)
    return CalibrationReference(**fields)


def native_antibody_ion(**overrides):
    """A native intact antibody ion, 24+ and protonated unless the test says otherwise."""
    fields = dict(
        analyte=antibody(folding_state=FoldingState.NATIVE),
        adduct="[M+24H]24+",
        charge=24,
        ccs=7000.0,
    )
    fields.update(overrides)
    return measurement(**fields)


# --- (a) the matched-ion key carries NO platform ----------------------------------------


@pytest.mark.parametrize("platform", sorted(PLATFORMS), ids=sorted(PLATFORMS))
def test_the_matched_ion_key_carries_no_platform_so_a_dtims_and_a_twims_value_of_one_ion_pair(platform):
    on_twims = measurement(**PLATFORMS["twims"])
    elsewhere = measurement(**PLATFORMS[platform])
    assert elsewhere.matched_ion_key == on_twims.matched_ion_key
    assert str(elsewhere.matched_ion_key) == str(on_twims.matched_ion_key)


def test_the_platforms_that_share_one_matched_ion_key_really_are_different_platforms():
    # Without this the parametrize above could pass by comparing a record with itself.
    records = [measurement(**fields) for fields in PLATFORMS.values()]
    assert len({(record.ims_type, record.dtims_method) for record in records}) == len(PLATFORMS)


def test_the_matched_ion_key_carries_neither_the_calibrant_nor_the_instrument_nor_the_source():
    # All provenance. If any of it entered the key, no two laboratories could ever
    # produce a matched pair at all.
    plain = measurement()
    elsewhere = measurement(
        calibrant="polyalanine",
        instrument="a different instrument",
        source="another laboratory",
        doi="10.1038/s41467-021-00001-2",
        source_locator="Table S2",
        replicates=5,
        measured_on="2021-05-04",
    )
    assert elsewhere.matched_ion_key == plain.matched_ion_key


def test_the_matched_ion_key_uses_the_gas_the_value_refers_to_and_not_the_gas_in_the_cell():
    # A TWIMS value measured in nitrogen but calibrated against helium references is
    # a helium value: it pools with helium values and not with nitrogen ones.
    helium_value_from_a_nitrogen_cell = measurement(drift_gas=DriftGas.HE, cell_gas=DriftGas.N2)
    helium_throughout = measurement(drift_gas=DriftGas.HE, cell_gas=DriftGas.HE)
    nitrogen_throughout = measurement(drift_gas=DriftGas.N2, cell_gas=DriftGas.N2)
    assert helium_value_from_a_nitrogen_cell.matched_ion_key == helium_throughout.matched_ion_key
    assert helium_value_from_a_nitrogen_cell.matched_ion_key != nitrogen_throughout.matched_ion_key


def test_every_drift_gas_including_unstated_keys_apart_from_every_other():
    keys = {measurement(drift_gas=gas).matched_ion_key for gas in DriftGas}
    assert len(keys) == len(DriftGas)


def test_the_matched_ion_key_states_the_signed_charge():
    # No two valid records can differ in charge alone, because the adduct fixes the
    # charge. So a key that dropped the charge shows only on the key itself.
    assert measurement().matched_ion_key.charge == 1
    assert measurement(adduct="[M-H]-", charge=-1, polarity=Polarity.NEGATIVE).matched_ion_key.charge == -1
    assert native_antibody_ion().matched_ion_key.charge == 24
    assert "+24" in str(native_antibody_ion().matched_ion_key)


def test_two_charge_states_of_one_antibody_are_different_matched_ions():
    forty_plus = measurement(
        analyte=antibody(folding_state=FoldingState.NATIVE), adduct="[M+40H]40+", charge=40, ccs=9000.0
    )
    assert forty_plus.matched_ion_key != native_antibody_ion().matched_ion_key


def test_a_native_and_a_denatured_ion_of_one_antibody_are_different_matched_ions():
    # Same antibody, same adduct, same charge, same gas: the folding state is the
    # ONLY difference, so this is the pair that catches a key which drops the state.
    native = native_antibody_ion()
    denatured = native_antibody_ion(analyte=antibody(folding_state=FoldingState.DENATURED))
    assert native.matched_ion_key != denatured.matched_ion_key


def test_a_native_24_plus_and_a_denatured_40_plus_ion_of_one_antibody_are_different_matched_ions():
    denatured_40 = measurement(
        analyte=antibody(folding_state=FoldingState.DENATURED), adduct="[M+40H]40+", charge=40, ccs=11000.0
    )
    assert denatured_40.matched_ion_key != native_antibody_ion().matched_ion_key


def test_every_folding_state_including_unstated_keys_apart_for_a_protein():
    keys = {measurement(analyte=protein(folding_state=state)).matched_ion_key for state in FoldingState}
    assert len(keys) == len(FoldingState)


def test_a_labelled_and_an_unlabelled_glycan_are_different_matched_ions():
    two_ab = measurement(analyte=glycan(reducing_end_label=ReducingEndLabel.TWO_AB))
    unlabelled = measurement(analyte=glycan(reducing_end_label=ReducingEndLabel.NATIVE))
    assert two_ab.matched_ion_key != unlabelled.matched_ion_key


def test_two_spellings_of_one_adduct_give_one_matched_ion_key():
    assert (
        measurement(adduct="[M+Na+H]2+", charge=2).matched_ion_key
        == measurement(adduct="[M+H+Na]2+", charge=2).matched_ion_key
    )


def test_an_adduct_is_stored_canonically_so_the_key_and_the_group_agree():
    assert measurement(adduct="[M+Na+H]2+", charge=2).adduct == "[M+H+Na]2+"


def test_a_sodiated_ion_is_not_a_protonated_one_of_the_same_charge():
    assert measurement(adduct="[M+Na]+").matched_ion_key != measurement(adduct="[M+H]+").matched_ion_key


def test_a_charge_carrier_the_source_never_named_keys_apart_from_protons():
    unstated = native_antibody_ion(adduct="[M+24?]24+")
    assert unstated.matched_ion_key != native_antibody_ion().matched_ion_key


def test_a_display_name_is_in_no_key_so_a_compound_name_alone_never_makes_a_match():
    one_name_two_molecules = (
        measurement(analyte=small_molecule(display_name="the compound")),
        measurement(analyte=small_molecule(display_name="the compound", inchikey=OTHER_INCHIKEY)),
    )
    assert one_name_two_molecules[0].matched_ion_key != one_name_two_molecules[1].matched_ion_key
    two_names_one_molecule = (
        measurement(analyte=small_molecule(display_name="one spelling")),
        measurement(analyte=small_molecule(display_name="another spelling")),
    )
    assert two_names_one_molecule[0].matched_ion_key == two_names_one_molecule[1].matched_ion_key


# Was an xfail against a real defect in identity.py GlycanAnalyte.identity_key,
# which fell through to self.composition.canonical on a record whose composition
# is None. Every record the package accepts can now produce its key.
def test_every_record_the_package_accepts_can_produce_its_matched_ion_key():
    # A structure string stating no linkages is a real shape in published tables:
    # the source names a residue and nothing about how anything is joined.
    unlinked_structure = glycan(composition=None, iupac_condensed="GlcNAc")
    assert measurement(analyte=unlinked_structure).matched_ion_key is not None


def test_the_matched_ion_key_ignores_the_conformer_index_so_two_peaks_of_one_ion_stay_one_ion():
    # Two arrival-time peaks under one set of conditions are a legitimate pair, not
    # two ions. Putting the index in the key would hide that they are one.
    first = measurement(conformer=1, conformers_total=2)
    second = measurement(conformer=2, conformers_total=2, ccs=330.0)
    assert first.matched_ion_key == second.matched_ion_key


# --- (b) the calibration group carries NO analyte ----------------------------------------


def test_the_calibration_group_carries_no_analyte_so_two_analytes_measured_alike_share_a_group():
    one_glycan = measurement(analyte=glycan(composition="Hex5HexNAc2"))
    another_glycan = measurement(analyte=glycan(composition="Hex3HexNAc2"))
    assert one_glycan.matched_ion_key != another_glycan.matched_ion_key  # really two analytes
    assert one_glycan.calibration_group == another_glycan.calibration_group
    assert str(one_glycan.calibration_group) == str(another_glycan.calibration_group)


def test_two_analyte_kinds_measured_alike_share_a_calibration_group():
    # A small molecule and a peptide have no structural state, so nothing but the
    # analyte itself could separate them, and the analyte is not in the group.
    assert measurement(analyte=small_molecule()).calibration_group == measurement(analyte=peptide()).calibration_group


def test_the_calibration_group_is_built_without_ever_asking_the_analyte_who_it_is():
    """The strongest form of (b): a record whose analyte cannot produce an identity
    key at all must still produce a calibration group, because the group never asks.

    The analyte here is a glycan identified only by an accession that does not match
    the expected GlyTouCan form, which the package accepts with a warning. Its
    identity key cannot be built; see the suspected defects reported with this file.
    """
    with pytest.warns(UnverifiedFormatWarning):
        unidentifiable = glycan(composition=None, glytoucan_ac="not-a-form")
    assert measurement(analyte=unidentifiable).calibration_group == measurement().calibration_group


def test_the_calibration_group_carries_the_calibrant_because_two_ladders_are_two_ways_to_a_number():
    assert measurement(calibrant="dextran").calibration_group != measurement(calibrant="polyalanine").calibration_group


def test_the_calibration_group_carries_the_platform_which_the_matched_ion_key_does_not():
    twims = measurement(**PLATFORMS["twims"])
    single_field = measurement(**PLATFORMS["dtims_single_field"])
    assert twims.matched_ion_key == single_field.matched_ion_key
    assert twims.calibration_group != single_field.calibration_group


def test_the_calibration_group_carries_the_gas_the_value_refers_to():
    assert measurement(drift_gas=DriftGas.HE).calibration_group != measurement(drift_gas=DriftGas.N2).calibration_group


def test_the_calibration_group_carries_the_structural_state():
    two_ab = measurement(analyte=glycan(reducing_end_label=ReducingEndLabel.TWO_AB))
    unlabelled = measurement(analyte=glycan(reducing_end_label=ReducingEndLabel.NATIVE))
    assert two_ab.calibration_group != unlabelled.calibration_group


def test_a_calibrant_may_not_be_a_placeholder_so_the_dash_in_a_group_string_cannot_collide():
    # The group prints "-" where there is no calibrant. That is safe only because a
    # calibrant spelt "-" is refused as a placeholder.
    assert "|-|" in str(primary().calibration_group)
    with pytest.raises(ValidationError, match="placeholder"):
        measurement(calibrant="-")


# --- the platform validators -------------------------------------------------------------


def test_a_dtims_record_must_say_whether_it_is_stepped_field_or_single_field():
    with pytest.raises(ValidationError, match="must state dtims_method"):
        measurement(ims_type=IMSType.DTIMS)
    assert measurement(**PLATFORMS["dtims_single_field"]).dtims_method is DTIMSMethod.SINGLE_FIELD


@pytest.mark.parametrize("ims_type", [IMSType.TWIMS, IMSType.TIMS, IMSType.CYCLIC])
def test_only_a_dtims_record_may_state_a_dtims_method(ims_type):
    with pytest.raises(ValidationError, match="applies only to DTIMS"):
        measurement(ims_type=ims_type, dtims_method=DTIMSMethod.STEPPED_FIELD)
    assert measurement(ims_type=ims_type).dtims_method is None


@pytest.mark.parametrize("platform", CALIBRATED_PLATFORMS, ids=CALIBRATED_PLATFORMS)
def test_a_calibrated_platform_with_no_calibrant_recorded_is_refused(platform):
    with pytest.raises(ValidationError, match="the calibrant is part of the measurement"):
        measurement(**{**PLATFORMS[platform], "calibrant": None})
    assert measurement(**PLATFORMS[platform]).ccs_is_calibrated is True


def test_a_primary_stepped_field_value_may_not_name_a_calibrant():
    with pytest.raises(ValidationError, match="no calibrant enters the value"):
        primary(calibrant="dextran")
    assert primary().calibrant is None
    assert primary().ccs_is_calibrated is False
    assert primary().calibration_lineage is CalibrationLineage.PRIMARY


def test_single_field_dtims_is_calibrated_and_not_primary():
    # Single-field DTIMS reads its CCS off a calibration curve exactly as TWIMS does.
    # Treating it as primary would let a calibrated value stand as first principles.
    single_field = measurement(**PLATFORMS["dtims_single_field"])
    assert single_field.ccs_is_calibrated is True
    assert single_field.calibration_lineage is CalibrationLineage.DERIVED
    assert primary().calibration_lineage is CalibrationLineage.PRIMARY


@pytest.mark.parametrize("platform", CALIBRATED_PLATFORMS, ids=CALIBRATED_PLATFORMS)
def test_every_platform_but_stepped_field_dtims_is_derived(platform):
    assert measurement(**PLATFORMS[platform]).calibration_lineage is CalibrationLineage.DERIVED


# --- the ion validators ------------------------------------------------------------------


def test_a_charge_of_zero_is_refused_because_ion_mobility_measures_ions():
    with pytest.raises(ValidationError, match="charge cannot be 0"):
        measurement(adduct="[M+H]+", charge=0)


def test_a_charge_must_agree_with_the_stated_polarity():
    with pytest.raises(ValidationError, match="does not match positive polarity"):
        measurement(adduct="[M-H]-", charge=-1, polarity=Polarity.POSITIVE)
    with pytest.raises(ValidationError, match="does not match negative polarity"):
        measurement(adduct="[M+H]+", charge=1, polarity=Polarity.NEGATIVE)
    assert measurement(adduct="[M-H]-", charge=-1, polarity=Polarity.NEGATIVE).charge == -1


def test_the_adduct_must_carry_the_charge_the_record_states():
    with pytest.raises(ValidationError, match="carries charge"):
        measurement(adduct="[M+2H]2+", charge=1)
    with pytest.raises(ValidationError, match="carries charge"):
        measurement(adduct="[M+H]+", charge=2)
    assert measurement(adduct="[M+2H]2+", charge=2).charge == 2


def test_a_primary_value_must_refer_to_the_gas_it_was_actually_measured_in():
    with pytest.raises(ValidationError, match="must equal drift_gas"):
        primary(cell_gas=DriftGas.HE)
    assert primary(cell_gas=DriftGas.N2).cell_gas is DriftGas.N2
    # A calibrated value may refer to a gas other than the one in the cell.
    assert measurement(drift_gas=DriftGas.HE, cell_gas=DriftGas.N2).cell_gas is DriftGas.N2


# --- every uncertainty carries its type --------------------------------------------------


def test_a_spread_with_no_statement_of_what_it_is_is_refused():
    with pytest.raises(ValidationError, match="ccs_uncertainty needs uncertainty_type"):
        measurement(ccs_uncertainty=5.0)
    assert measurement(ccs_uncertainty=5.0, uncertainty_type=UncertaintyType.SD).ccs_uncertainty == 5.0


def test_uncertainty_type_unknown_is_the_one_type_that_may_stand_alone():
    # It does not describe a number; it records that the source states no spread.
    record = measurement(uncertainty_type=UncertaintyType.UNKNOWN)
    assert record.ccs_uncertainty is None
    assert record.uncertainty_type is UncertaintyType.UNKNOWN


@pytest.mark.parametrize("kind", STATED_UNCERTAINTY_TYPES, ids=[k.value for k in STATED_UNCERTAINTY_TYPES])
def test_every_uncertainty_type_other_than_unknown_needs_a_number_to_describe(kind):
    with pytest.raises(ValidationError, match="may stand alone"):
        measurement(uncertainty_type=kind)
    assert measurement(ccs_uncertainty=5.0, uncertainty_type=kind).uncertainty_type is kind


# --- the conformer rules ------------------------------------------------------------------


def test_a_conformer_index_with_no_total_beside_it_is_refused():
    with pytest.raises(ValidationError, match="says nothing without conformers_total"):
        measurement(conformer=1)
    assert measurement(conformer=1, conformers_total=1).conformer == 1


def test_a_row_claiming_several_conformers_must_say_which_one_it_is():
    with pytest.raises(ValidationError, match="must say which of them it is"):
        measurement(conformers_total=3)
    assert measurement(conformer=3, conformers_total=3).conformer == 3
    assert measurement(conformers_total=1).conformers_total == 1


def test_a_conformer_index_may_not_exceed_the_total_the_source_resolved():
    with pytest.raises(ValidationError, match="the index exceeds the total"):
        measurement(conformer=3, conformers_total=2)
    assert measurement(conformer=2, conformers_total=2).conformer == 2


# --- calibration reference lineage --------------------------------------------------------


def test_a_primary_value_may_not_carry_a_calibration_reference():
    with pytest.raises(ValidationError, match="cannot have a calibration reference"):
        primary(calibration_reference=reference())
    assert measurement(calibration_reference=reference()).calibration_reference == reference()


@pytest.mark.parametrize(
    "platform", [None, IMSType.TWIMS, IMSType.TIMS, IMSType.CYCLIC], ids=["unstated", "twims", "tims", "cyclic"]
)
def test_a_reference_set_may_state_a_dtims_method_only_where_it_is_itself_dtims(platform):
    with pytest.raises(ValidationError, match="applies only to DTIMS"):
        reference(platform=platform, method=DTIMSMethod.STEPPED_FIELD)
    assert reference(platform=IMSType.DTIMS, method=DTIMSMethod.STEPPED_FIELD).method is DTIMSMethod.STEPPED_FIELD


@pytest.mark.parametrize(
    "reference_fields, expected",
    [
        pytest.param(None, None, id="a_derived_value_whose_reference_set_is_unrecorded"),
        pytest.param({}, None, id="a_reference_set_that_does_not_say_what_it_was_measured_on"),
        pytest.param(dict(platform=IMSType.DTIMS), None, id="a_dtims_reference_set_that_does_not_say_which_method"),
        pytest.param(dict(platform=IMSType.TWIMS), False, id="a_twims_reference_set_is_itself_calibrated"),
        pytest.param(dict(platform=IMSType.TIMS), False, id="a_tims_reference_set_is_itself_calibrated"),
        pytest.param(
            dict(platform=IMSType.DTIMS, method=DTIMSMethod.SINGLE_FIELD),
            False,
            id="a_single_field_dtims_reference_set_is_itself_calibrated",
        ),
        pytest.param(
            dict(platform=IMSType.DTIMS, method=DTIMSMethod.STEPPED_FIELD),
            True,
            id="a_stepped_field_dtims_reference_set_is_primary",
        ),
    ],
)
def test_traces_to_primary_answers_the_circularity_question_and_says_unknown_when_it_is_unknown(
    reference_fields, expected
):
    # None is a real answer here and not a missing one: it is the state most
    # published TWIMS values are in, and reporting it as True would hide exactly the
    # circularity this field was added to expose.
    reference_set = None if reference_fields is None else reference(**reference_fields)
    assert measurement(calibration_reference=reference_set).traces_to_primary is expected
    if reference_set is not None:
        assert reference_set.traces_to_primary is expected


def test_a_stepped_field_value_traces_to_primary_with_no_reference_set_at_all():
    assert primary().traces_to_primary is True


def test_a_calibration_reference_records_a_bare_doi_or_nothing():
    assert reference(doi="10.1038/s41467-021-00001-2").doi == "10.1038/s41467-021-00001-2"
    assert reference().doi is None
    with pytest.raises(ValidationError, match="is not a DOI"):
        reference(doi="https://doi.org/10.1038/s41467-021-00001-2")


# --- what blocks training ------------------------------------------------------------------


def test_a_clean_record_has_no_training_blockers(glycan_measurement):
    assert glycan_measurement.training_blockers() == []


def test_an_unstated_drift_gas_blocks_training():
    blockers = measurement(drift_gas=DriftGas.UNSTATED).training_blockers()
    assert any("drift_gas" in blocker for blocker in blockers)
    assert measurement(drift_gas=DriftGas.N2).training_blockers() == []


def test_an_uncertainty_of_unknown_kind_blocks_training():
    blockers = measurement(uncertainty_type=UncertaintyType.UNKNOWN).training_blockers()
    assert any("uncertainty_type" in blocker for blocker in blockers)


@pytest.mark.parametrize("kind", STATED_UNCERTAINTY_TYPES, ids=[k.value for k in STATED_UNCERTAINTY_TYPES])
def test_a_spread_whose_kind_is_stated_does_not_block_training(kind):
    assert measurement(ccs_uncertainty=5.0, uncertainty_type=kind).training_blockers() == []


def test_an_adduct_that_does_not_name_its_charge_carrier_blocks_training():
    unstated = native_antibody_ion(adduct="[M+24?]24+")
    assert any("charge carrier" in blocker for blocker in unstated.training_blockers())
    assert native_antibody_ion().training_blockers() == []


def test_a_record_reports_every_blocker_it_has_and_not_only_the_first():
    record = measurement(drift_gas=DriftGas.UNSTATED, uncertainty_type=UncertaintyType.UNKNOWN)
    assert len(record.training_blockers()) >= 2


# --- what the licence gate is offered --------------------------------------------------------


def test_component_records_offers_the_analyte_to_the_licence_gate(glycan_measurement):
    # The analyte carries its own reuse status: a CCS value from an open paper does
    # not clear an identity taken from a restricted one.
    assert glycan_measurement.component_records() == (glycan_measurement.analyte,)


@pytest.mark.parametrize("build", [antibody, adc], ids=["intact_antibody", "adc"])
def test_the_nested_antibody_identity_is_not_offered_as_a_separately_licensed_record(build):
    # It was, and the consequence was severe enough to be worth a test of its own.
    # AntibodyIdentity carries no reuse status, the gate refuses anything in
    # component_records() that cannot state one, and so EVERY intact-antibody and
    # ADC measurement was refused whatever its licence - closing the whole
    # biopharmaceutical layer, which is the platform's stated differentiator.
    #
    # It is the structured form of the enclosing analyte's identity, not a record
    # sourced from somewhere else, and the analyte's own source and reuse_status
    # already say where that identity came from and on what terms.
    analyte = build()
    parts = measurement(analyte=analyte, ccs=7000.0).component_records()
    assert parts == (analyte,)
    assert analyte.antibody not in parts


@pytest.mark.parametrize("build", [antibody, adc], ids=["intact_antibody", "adc"])
def test_an_antibody_measurement_can_clear_the_licence_gate(build):
    # The negative direction of the test above, and the one that would have caught
    # the defect: it is not enough that the nested identity is absent from the
    # tuple, the record has to actually pass the gate.
    #
    # The folding state has to be stated for it to pass, which is the point of
    # the other guard and not an inconvenience: an ion of unrecorded conformation
    # is not defined well enough to train on.
    analyte = build(folding_state=FoldingState.NATIVE)
    assert_trainable(measurement(analyte=analyte, ccs=7000.0))


@pytest.mark.parametrize("build", [antibody, adc], ids=["intact_antibody", "adc"])
def test_an_antibody_whose_folding_state_is_unstated_is_still_refused(build):
    with pytest.raises(TrainingGateError, match="folding_state"):
        assert_trainable(measurement(analyte=build(), ccs=7000.0))


def test_a_measurement_reports_its_analytes_warnings_without_blocking_on_them():
    with pytest.warns(UnverifiedFormatWarning):
        odd_accession = glycan(glytoucan_ac="not-a-form")
    warnings_reported = measurement(analyte=odd_accession).validation_warnings
    assert warnings_reported and all(warning.startswith("analyte: ") for warning in warnings_reported)
    assert measurement().validation_warnings == []


# --- a record says what it was given, and nothing else ----------------------------------------


def test_a_measurement_cannot_be_changed_after_validation(glycan_measurement):
    # An original value is never replaced. A record that could be edited in place
    # would let a harmonized number overwrite the measurement it came from.
    with pytest.raises(ValidationError):
        glycan_measurement.ccs = 400.0


def test_a_misspelt_field_name_is_refused_rather_than_silently_ignored():
    with pytest.raises(ValidationError, match="ccs_uncertianty"):
        measurement(ccs_uncertianty=5.0)


@pytest.mark.parametrize("value", [0.0, -1.0, float("inf"), float("nan")])
def test_a_ccs_must_be_a_real_positive_number(value):
    with pytest.raises(ValidationError):
        measurement(ccs=value)


def test_a_ccs_written_as_text_is_not_quietly_read_as_a_number():
    with pytest.raises(ValidationError):
        measurement(ccs="300.0")


@pytest.mark.parametrize("filler", ["n/a", "N.D.", "-", "unknown", "not reported", "   "])
def test_a_placeholder_is_refused_where_a_value_is_required(filler):
    with pytest.raises(ValidationError, match="placeholder|at least 1 character"):
        measurement(source=filler)


@pytest.mark.parametrize("written", ["1462060800", 1462060800, "04/05/2021", "2021-5-4", "May 2021"])
def test_a_measurement_date_is_read_only_as_written_so_a_count_never_becomes_a_date(written):
    with pytest.raises(ValidationError, match="written YYYY-MM-DD"):
        measurement(measured_on=written)


def test_a_measurement_date_is_a_date_and_not_a_datetime():
    # Cutting a time of day off would be a silent change to what the source claimed.
    with pytest.raises(ValidationError, match="not a datetime"):
        measurement(measured_on=datetime(2021, 5, 4, 12, 0))
    assert measurement(measured_on="2021-05-04").measured_on == date(2021, 5, 4)
    assert measurement(measured_on=date(2021, 5, 4)).measured_on == date(2021, 5, 4)


def test_a_measurement_doi_is_the_bare_doi_and_not_a_url():
    assert measurement(doi="10.1038/s41467-021-00001-2").doi == "10.1038/s41467-021-00001-2"
    with pytest.raises(ValidationError, match="is not a DOI"):
        measurement(doi="doi:10.1038/s41467-021-00001-2")


# --- the conditions on their own --------------------------------------------------------------


def conditions(**overrides) -> MeasurementConditions:
    fields = dict(
        adduct="[M+H]+",
        charge=1,
        polarity=Polarity.POSITIVE,
        ims_type=IMSType.TWIMS,
        drift_gas=DriftGas.N2,
        calibrant="dextran",
    )
    fields.update(overrides)
    return MeasurementConditions(**fields)


def test_a_bare_set_of_conditions_is_held_to_the_same_rules_as_a_measurement():
    # Conditions are shared with harmonization queries later. A query typed more
    # loosely than a record would be answered with records it does not match.
    assert conditions().ccs_is_calibrated is True
    assert conditions().calibration_lineage is CalibrationLineage.DERIVED
    with pytest.raises(ValidationError, match="must state dtims_method"):
        conditions(ims_type=IMSType.DTIMS)
    with pytest.raises(ValidationError, match="charge cannot be 0"):
        conditions(charge=0)
    with pytest.raises(ValidationError, match="carries charge"):
        conditions(adduct="[M+2H]2+")
    with pytest.raises(ValidationError, match="no calibrant enters the value"):
        conditions(ims_type=IMSType.DTIMS, dtims_method=DTIMSMethod.STEPPED_FIELD)


# --- an ion whose charge carrier the source never named ---------------------------------
#
# Added after the first mutation sweep, which is when this behaviour was added.
# Native-MS papers routinely report a charge state and nothing else - "the 24+
# ion" - without saying whether those charges are protons, sodium or ammonium.
# The Bush Lab protein data is expected to be full of them, so this is the common
# case in the biopharmaceutical layer rather than an exotic one.


def unstated_carrier_ion(**overrides):
    fields = dict(
        analyte=protein(folding_state=FoldingState.NATIVE),
        adduct="[M+24?]24+",
        charge=24,
        ccs=7000.0,
    )
    fields.update(overrides)
    return measurement(**fields)


def test_an_ion_whose_charge_carrier_is_unstated_can_be_recorded_at_all():
    # The alternative to this form is writing [M+24H]24+, which asserts protons,
    # which the paper does not say. Refusing to record the measurement and
    # inventing the chemistry are both worse than recording the ambiguity.
    record = unstated_carrier_ion()
    assert record.adduct == "[M+24?]24+"
    assert record.charge == 24


def test_an_unstated_charge_carrier_blocks_training_like_the_other_unknowns():
    blockers = unstated_carrier_ion().training_blockers()
    assert any("charge carrier" in blocker for blocker in blockers)
    with pytest.raises(TrainingGateError, match="charge carrier"):
        assert_trainable(unstated_carrier_ion())


def test_two_ions_whose_charge_carrier_is_unstated_do_not_match_each_other():
    # THE POINT OF THE WHOLE FORM. Two papers both reporting "the 24+ ion" of one
    # protein have not reported the same ion: one may be twenty-four protons and
    # the other twenty-four ammonium adducts, which differ by 408 Da and do not
    # have the same cross section. Pairing them would report the difference
    # between two different ions as inter-platform bias.
    one = unstated_carrier_ion(source="one laboratory", ccs=7000.0)
    other = unstated_carrier_ion(source="another laboratory", ccs=7010.0)
    assert one.matched_ion_key != other.matched_ion_key
    assert not one.matched_ion_key.matchable
    assert not other.matched_ion_key.matchable
    assert "unmatchable" in str(one.matched_ion_key)


def test_an_unstated_carrier_ion_does_not_match_the_same_ion_with_a_named_carrier():
    named = measurement(
        analyte=protein(folding_state=FoldingState.NATIVE), adduct="[M+24H]24+", charge=24, ccs=7000.0
    )
    assert unstated_carrier_ion().matched_ion_key != named.matched_ion_key


def test_an_ion_whose_carrier_IS_named_keeps_an_ordinary_matchable_key():
    # The negative direction: the unmatchable component must be null for every
    # ordinary record, or nothing would ever pair with anything.
    named = measurement(
        analyte=protein(folding_state=FoldingState.NATIVE), adduct="[M+24H]24+", charge=24, ccs=7000.0
    )
    assert named.matched_ion_key.matchable
    assert named.matched_ion_key.unmatchable is None
    elsewhere = measurement(
        analyte=protein(folding_state=FoldingState.NATIVE),
        adduct="[M+24H]24+",
        charge=24,
        ccs=7100.0,
        source="another laboratory",
        ims_type=IMSType.TIMS,
    )
    assert named.matched_ion_key == elsewhere.matched_ion_key
