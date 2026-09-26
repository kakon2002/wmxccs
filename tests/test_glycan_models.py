"""Record-level constraints.

Every record built here is a synthetic fixture for exercising validation.
None of the values is a real measurement or a real structure assignment.
"""

import os
import subprocess
import sys
import warnings
from enum import StrEnum

import pytest
from pydantic import ValidationError

import wmxglycan.composition
import wmxglycan.licensing
import wmxglycan.models
from wmxglycan.composition import Composition
from wmxglycan.licensing import LicenceGateError, ReuseStatus, TrainingGateError, assert_trainable
from wmxglycan.models import (
    _PLACEHOLDERS,  # the spellings a loader might use for a missing value
    _is_placeholder,
    adduct_components,
    CalibrationGroup,
    CCSMeasurement,
    Derivatisation,
    DriftGas,
    DTIMSMethod,
    GlycanStructure,
    IMSType,
    ReducingEndLabel,
    UncertaintyType,
    UnverifiedFormatWarning,
)

FIXTURE_SOURCE = "synthetic test fixture, not a real record"
FIXTURE_CCS = 100.0  # arbitrary; only validation is under test
FIXTURE_CALIBRANT = "fixture calibrant"
FIXTURE_WURCS = "WURCS=2.0/synthetic-fixture"
MISSING = object()
RESOLVED = {"has_unresolved_linkage": False, "has_unresolved_anomericity": False}
SIALYLATED = "Hex5HexNAc4Fuc1NeuAc2"

STEPPED = {"ims_type": "DTIMS", "dtims_method": "stepped_field", "calibrant": MISSING}
SINGLE = {"ims_type": "DTIMS", "dtims_method": "single_field"}

LABEL_VALUES = {"native", "reduced", "PA", "2-AB", "2-AA", "procainamide", "RapiFluor-MS", "APTS", "other", "unknown"}
DERIVATISATION_VALUES = {"underivatised", "permethylation", "sialic_acid_amidation", "sialic_acid_esterification", "unknown"}
SIALIC_ACID_DERIVATISATIONS = ["sialic_acid_amidation", "sialic_acid_esterification"]
TRAINABLE_LABELS = sorted(LABEL_VALUES - {"other", "unknown"})
TRAINABLE_DERIVATISATIONS = sorted(DERIVATISATION_VALUES - {"unknown", *SIALIC_ACID_DERIVATISATIONS})

# Fillers a literature or dataframe loader might write for a missing cell. None of them is a value.
LOADER_FILLERS = ["n.d.", "not reported", "#N/A", "<NA>", "NaT", "--", "—", "–", "..."]


class _Ambiguous:
    def __bool__(self):
        raise TypeError("boolean value of NA is ambiguous")


class _AmbiguousArray:
    def __bool__(self):
        raise ValueError("the truth value of an array with more than one element is ambiguous")


class NALike:
    """Stands in for pandas.NA: hashable, but comparing it gives something with no truth value."""

    __hash__ = object.__hash__

    def __eq__(self, other):
        return _Ambiguous()


class ArrayLike:
    """Stands in for a numpy array: unhashable, and its comparison has no single truth value."""

    __hash__ = None

    def __eq__(self, other):
        return _AmbiguousArray()


class CollidesWithOther(NALike):
    """NA-like, with the hash of 'other', so a set lookup has to compare it with ReducingEndLabel.OTHER."""

    def __hash__(self):
        return hash("other")


def make_glycan(**overrides):
    return GlycanStructure(**({"composition": "Hex5HexNAc4Fuc1", "source": FIXTURE_SOURCE} | overrides))


def measurement_fields(**overrides):
    fields = {
        "glycan": make_glycan(),
        "ccs": FIXTURE_CCS,
        "adduct": "[M+2H]2+",
        "charge": 2,
        "polarity": "positive",
        "reducing_end_label": "native",
        "derivatisation": "underivatised",
        "ims_type": "TWIMS",
        "drift_gas": "N2",
        "calibrant": FIXTURE_CALIBRANT,
        "source": FIXTURE_SOURCE,
    }
    fields |= overrides
    return {k: v for k, v in fields.items() if v is not MISSING}


def make_measurement(**overrides):
    return CCSMeasurement(**measurement_fields(**overrides))


def cleared(**overrides):
    """A measurement whose licences both permit training.

    Declared synthetic, on both the value and the structure: the gate accepts
    that without a licence record, and the file loader refuses it, so nothing
    built here can be mistaken for data.
    """
    return make_measurement(
        **({"reuse_status": "synthetic_fixture", "glycan": make_glycan(reuse_status="synthetic_fixture")} | overrides)
    )


# --- GlycanStructure --------------------------------------------------------


def test_glycan_composition_is_parsed_and_serialised_canonically():
    glycan = make_glycan(composition="HexNAc4dHex1Hex5")
    assert glycan.composition == Composition(hex=5, hexnac=4, fuc=1)
    assert glycan.model_dump()["composition"] == "Hex5HexNAc4Fuc1"
    assert GlycanStructure.model_validate_json(glycan.model_dump_json()) == glycan


def test_glycan_defaults_to_unverified_and_unresolved():
    glycan = make_glycan()
    assert glycan.reuse_status is ReuseStatus.UNVERIFIED
    assert glycan.has_unresolved_linkage is True
    assert glycan.has_unresolved_anomericity is True


def test_glycan_rejects_a_typo_composition():
    with pytest.raises(ValidationError, match="HexNac"):
        make_glycan(composition="Hex5HexNac4")


@pytest.mark.parametrize("flag", ["has_unresolved_linkage", "has_unresolved_anomericity"])
def test_composition_only_record_cannot_claim_to_be_resolved(flag):
    with pytest.raises(ValidationError, match="only a composition"):
        make_glycan(**{flag: False})


@pytest.mark.parametrize("identifier", [{"wurcs": FIXTURE_WURCS}, {"glytoucan_ac": "G00000AA"}])
def test_a_structure_identifier_allows_resolved_flags(identifier):
    glycan = make_glycan(**RESOLVED, **identifier)
    assert glycan.has_unresolved_linkage is False


@pytest.mark.parametrize("value", ["G1234AB", "g00000aa", "G00000A", "00000AA", "G00000AAA"])
def test_an_unexpected_glytoucan_format_warns_but_is_kept(value):
    with pytest.warns(UnverifiedFormatWarning, match="GlyTouCan"):
        glycan = make_glycan(glytoucan_ac=value)
    assert glycan.glytoucan_ac == value
    assert len(glycan.validation_warnings) == 1


@pytest.mark.parametrize("value", ["G1234AB", "pending", "see Table S2", "Hex5HexNAc4Fuc1"])
def test_an_unverified_accession_cannot_make_a_composition_only_record_resolved(value):
    # Unresolved stays unresolved: text that is not a well-formed accession is not structural evidence.
    with pytest.warns(UnverifiedFormatWarning), pytest.raises(ValidationError, match="cannot claim resolved linkage"):
        make_glycan(glytoucan_ac=value, **RESOLVED)


def test_an_unverified_accession_is_fine_alongside_a_wurcs():
    with pytest.warns(UnverifiedFormatWarning):
        glycan = make_glycan(glytoucan_ac="G1234AB", wurcs=FIXTURE_WURCS, **RESOLVED)
    assert glycan.has_unresolved_linkage is False


def test_the_expected_glytoucan_format_is_silent():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        glycan = make_glycan(glytoucan_ac="G00000AA")
    assert glycan.validation_warnings == []


def test_the_glytoucan_warning_is_shown_by_default():
    assert issubclass(UnverifiedFormatWarning, UserWarning)
    code = (
        "from wmxglycan.models import GlycanStructure\n"
        "GlycanStructure(composition='Hex5HexNAc4Fuc1', source='synthetic test fixture', glytoucan_ac='G1234AB')\n"
    )
    env = {k: v for k, v in os.environ.items() if k != "PYTHONWARNINGS"}
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
    assert result.returncode == 0, result.stderr
    assert "UnverifiedFormatWarning" in result.stderr


@pytest.mark.parametrize("value", ["", "   ", "unknown", "N/A", *LOADER_FILLERS])
def test_a_blank_or_placeholder_glytoucan_accession_is_still_rejected(value):
    with pytest.raises(ValidationError):
        make_glycan(glytoucan_ac=value)


def test_a_measurement_surfaces_its_glycans_warnings():
    with pytest.warns(UnverifiedFormatWarning):
        glycan = make_glycan(glytoucan_ac="G1234AB")
    assert make_measurement(glycan=glycan).validation_warnings == [f"glycan: {glycan.validation_warnings[0]}"]


@pytest.mark.parametrize("value", ["", "wurcs=2.0/x", "WURCS=2.0/a b", "G00000AA"])
def test_wurcs_format(value):
    with pytest.raises(ValidationError, match="WURCS"):
        make_glycan(wurcs=value)


@pytest.mark.parametrize("value", ["no", "false", 0])
def test_resolution_flags_must_be_real_booleans(value):
    with pytest.raises(ValidationError):
        make_glycan(has_unresolved_linkage=value)


@pytest.mark.parametrize("source", [MISSING, "", "   ", "unknown", "N/A", "none", *LOADER_FILLERS])
def test_every_record_needs_a_real_source(source):
    fields = {"composition": "Hex5HexNAc4Fuc1"}
    if source is not MISSING:
        fields["source"] = source
    with pytest.raises(ValidationError):
        GlycanStructure(**fields)


def test_records_are_frozen():
    with pytest.raises(ValidationError):
        make_glycan().reuse_status = ReuseStatus.OPEN_ATTRIBUTION
    with pytest.raises(ValidationError):
        make_measurement().calibrant = None
    with pytest.raises(ValidationError):
        make_measurement().reducing_end_label = ReducingEndLabel.PA
    with pytest.raises(ValidationError):
        make_measurement().derivatisation = Derivatisation.PERMETHYLATION


def test_misspelt_fields_are_rejected_not_ignored():
    with pytest.raises(ValidationError):
        make_glycan(reuse_stauts="open_attribution")
    with pytest.raises(ValidationError):
        make_measurement(calibrent="x")


@pytest.mark.parametrize(
    "glycan_fields",
    [{"source": FIXTURE_SOURCE}, {"source": FIXTURE_SOURCE, "composition": "Hex5HexNAc4Fuc1"}],
    ids=["no-composition", "unparsed-composition"],
)
def test_a_measurement_of_a_glycan_built_without_validation_is_rejected_cleanly(glycan_fields):
    # A raw AttributeError would escape a loader that catches ValidationError to report bad rows.
    with pytest.raises(ValidationError, match="built without validation"):
        make_measurement(glycan=GlycanStructure.model_construct(**glycan_fields))


# --- CCSMeasurement ---------------------------------------------------------


def test_valid_twims_measurement():
    measurement = make_measurement()
    assert measurement.ims_type is IMSType.TWIMS
    assert measurement.dtims_method is None
    assert measurement.calibrant == FIXTURE_CALIBRANT
    assert measurement.reuse_status is ReuseStatus.UNVERIFIED


def test_twims_record_without_calibrant_is_rejected():
    with pytest.raises(ValidationError, match="TWIMS gives calibrated CCS"):
        make_measurement(ims_type="TWIMS", calibrant=MISSING)


CALIBRATED = [{"ims_type": "TWIMS"}, {"ims_type": "TIMS"}, {"ims_type": "CYCLIC"}, SINGLE]
NO_CALIBRANT = [MISSING, None, "", "   ", "none", "N/A", "unknown", *LOADER_FILLERS]


@pytest.mark.parametrize("platform", CALIBRATED, ids=["TWIMS", "TIMS", "CYCLIC", "DTIMS-single-field"])
@pytest.mark.parametrize("calibrant", NO_CALIBRANT, ids=[f"calibrant{i}" for i in range(len(NO_CALIBRANT))])
def test_calibrated_platforms_need_a_real_calibrant(platform, calibrant):
    with pytest.raises(ValidationError):
        make_measurement(**(platform | {"calibrant": calibrant}))


@pytest.mark.parametrize("filler", LOADER_FILLERS)
def test_a_loader_filler_is_not_an_instrument_either(filler):
    with pytest.raises(ValidationError, match="placeholder"):
        make_measurement(instrument=filler)


@pytest.mark.parametrize(
    "platform, calibrated",
    [({"ims_type": "TWIMS"}, True), ({"ims_type": "TIMS"}, True), ({"ims_type": "CYCLIC"}, True), (SINGLE, True), (STEPPED, False)],
    ids=["TWIMS", "TIMS", "CYCLIC", "DTIMS-single-field", "DTIMS-stepped-field"],
)
def test_which_platforms_give_calibrated_ccs(platform, calibrated):
    assert make_measurement(**platform).ccs_is_calibrated is calibrated


def test_dtims_methods_are_exactly_the_agreed_pair():
    assert {method.value for method in DTIMSMethod} == {"stepped_field", "single_field"}


def test_a_dtims_record_must_state_its_method():
    with pytest.raises(ValidationError, match="dtims_method"):
        make_measurement(ims_type="DTIMS", calibrant=MISSING)
    with pytest.raises(ValidationError, match="dtims_method"):
        make_measurement(ims_type="DTIMS")


@pytest.mark.parametrize("method", ["stepped_field", "single_field"])
@pytest.mark.parametrize("platform", ["TWIMS", "TIMS", "CYCLIC"])
def test_dtims_method_is_rejected_for_other_platforms(platform, method):
    with pytest.raises(ValidationError, match="only to DTIMS"):
        make_measurement(ims_type=platform, dtims_method=method)


@pytest.mark.parametrize("method", ["stepped", "single-field", "STEPPED_FIELD", "unknown", ""])
def test_unknown_dtims_methods_are_rejected(method):
    with pytest.raises(ValidationError):
        make_measurement(ims_type="DTIMS", dtims_method=method)


def test_stepped_field_dtims_needs_no_calibrant():
    measurement = make_measurement(**STEPPED)
    assert measurement.dtims_method is DTIMSMethod.STEPPED_FIELD
    assert measurement.calibrant is None


def test_stepped_field_dtims_may_not_record_a_calibrant():
    with pytest.raises(ValidationError, match="primary CCS"):
        make_measurement(ims_type="DTIMS", dtims_method="stepped_field")


def test_single_field_dtims_without_calibrant_is_rejected():
    with pytest.raises(ValidationError, match="single_field gives calibrated CCS"):
        make_measurement(**SINGLE, calibrant=MISSING)


def test_single_field_dtims_with_calibrant():
    measurement = make_measurement(**SINGLE)
    assert measurement.dtims_method is DTIMSMethod.SINGLE_FIELD
    assert measurement.calibrant == FIXTURE_CALIBRANT


@pytest.mark.parametrize("charge, adduct, polarity", [(2, "[M+2H]2+", "negative"), (-1, "[M-H]-", "positive")])
def test_charge_sign_must_match_polarity(charge, adduct, polarity):
    with pytest.raises(ValidationError, match="polarity"):
        make_measurement(charge=charge, adduct=adduct, polarity=polarity)


def test_negative_mode_measurement():
    assert make_measurement(charge=-1, adduct="[M-H]-", polarity="negative").charge == -1


def test_zero_charge_is_rejected():
    with pytest.raises(ValidationError, match="charge cannot be 0"):
        make_measurement(charge=0)


@pytest.mark.parametrize("charge", [True, 2.0, "2"])
def test_charge_must_be_an_int(charge):
    with pytest.raises(ValidationError):
        make_measurement(charge=charge)


def test_adduct_charge_must_match_charge():
    with pytest.raises(ValidationError, match="carries charge"):
        make_measurement(adduct="[M+H]+", charge=2)


@pytest.mark.parametrize("adduct", ["M+2H", "[M+2H]", "[M+2H]+2", "[M+2H]++", "[2H]2+", "[M+2H]02+", "[M+2H] 2+", "[0M+H]+"])
def test_malformed_adducts_are_rejected(adduct):
    with pytest.raises(ValidationError, match="bracket notation"):
        make_measurement(adduct=adduct)


# --- adduct normalisation ---------------------------------------------------


@pytest.mark.parametrize(
    "spelling, canonical, charge",
    [
        ("[M+Na+H]2+", "[M+H+Na]2+", 2),
        ("[M+H+Na]2+", "[M+H+Na]2+", 2),
        ("[M-2H+Na]-", "[M+Na-2H]-", -1),
        ("[M-H2O+H]+", "[M+H-H2O]+", 1),
        ("[M-H2O-H]-", "[M-H-H2O]-", -1),  # losses are sorted too
        ("[M+1H]1+", "[M+H]+", 1),
        ("[1M+Na]+", "[M+Na]+", 1),
        ("[2M+Na]+", "[2M+Na]+", 1),
        ("[M+10H]10+", "[M+10H]10+", 10),  # counts of ten or more
        ("[M+NH4]+", "[M+NH4]+", 1),
        ("[M+HCOO]-", "[M+HCOO]-", -1),
    ],
)
def test_adduct_components_are_sorted_into_one_canonical_spelling(spelling, canonical, charge):
    polarity = "positive" if charge > 0 else "negative"
    assert make_measurement(adduct=spelling, charge=charge, polarity=polarity).adduct == canonical


def test_reordered_adducts_share_a_calibration_group():
    a = make_measurement(adduct="[M+Na+H]2+")
    b = make_measurement(adduct="[M+H+Na]2+")
    assert a.adduct == b.adduct == "[M+H+Na]2+"
    assert a.calibration_group == b.calibration_group


@pytest.mark.parametrize("adduct", ["[M+H+H]2+", "[M+Na-Na+H]+", "[M+2H+H]3+"])
def test_an_adduct_naming_a_species_twice_is_rejected(adduct):
    with pytest.raises(ValidationError, match="more than once"):
        make_measurement(adduct=adduct)


@pytest.mark.parametrize("adduct", ["[M|H]+", "[M+h]+", "[M+0H]+", "[M+02H]2+", "[M++H]+", "[M+]+", "[M+H-]+"])
def test_malformed_adduct_components_are_rejected(adduct):
    with pytest.raises(ValidationError, match="expected components"):
        make_measurement(adduct=adduct)


def test_an_overlong_adduct_is_rejected_before_it_is_parsed():
    with pytest.raises(ValidationError, match="at most 100 characters"):
        make_measurement(adduct="[M" + "+H" * 60 + "]+")


@pytest.mark.parametrize("ccs", [0.0, -1.0, float("nan"), float("inf"), True, "100.0"])
def test_ccs_must_be_a_positive_finite_number(ccs):
    with pytest.raises(ValidationError):
        make_measurement(ccs=ccs)


@pytest.mark.parametrize(
    "field, value", [("ims_type", "FAIMS"), ("ims_type", "twims"), ("drift_gas", "nitrogen"), ("polarity", "+")]
)
def test_unknown_enum_values_are_rejected(field, value):
    with pytest.raises(ValidationError):
        make_measurement(**{field: value})


# --- reducing-end label and derivatisation ----------------------------------


def test_reducing_end_labels_are_exactly_the_agreed_list():
    assert {label.value for label in ReducingEndLabel} == LABEL_VALUES


def test_derivatisations_are_exactly_the_agreed_list():
    assert {derivatisation.value for derivatisation in Derivatisation} == DERIVATISATION_VALUES


@pytest.mark.parametrize("label", sorted(LABEL_VALUES))
@pytest.mark.parametrize("derivatisation", sorted(DERIVATISATION_VALUES))
def test_label_and_derivatisation_are_independent(label, derivatisation):
    # A sialylated glycan, so every derivatisation, sialic-acid ones included, applies.
    measurement = make_measurement(
        reducing_end_label=label, derivatisation=derivatisation, glycan=make_glycan(composition=SIALYLATED)
    )
    assert (measurement.reducing_end_label, measurement.derivatisation) == (label, derivatisation)


@pytest.mark.parametrize("label", ["2AB", "2-ab", "pa", "Rapifluor-MS", "permethylated", "labelled", ""])
def test_misspelt_or_retired_labels_are_rejected(label):
    with pytest.raises(ValidationError):
        make_measurement(reducing_end_label=label)


@pytest.mark.parametrize("derivatisation", ["none", "permethylated", "Permethylation", "amidation", "esterified", "other", ""])
def test_unagreed_or_retired_derivatisations_are_rejected(derivatisation):
    with pytest.raises(ValidationError):
        make_measurement(derivatisation=derivatisation)


def test_label_and_derivatisation_default_to_unknown():
    measurement = make_measurement(reducing_end_label=MISSING, derivatisation=MISSING)
    assert measurement.reducing_end_label is ReducingEndLabel.UNKNOWN
    assert measurement.derivatisation is Derivatisation.UNKNOWN


def test_the_old_reducing_end_field_name_is_rejected():
    with pytest.raises(ValidationError):
        make_measurement(reducing_end="native")
    # On its own too, so it would fail if the old name ever became an alias for the new one.
    with pytest.raises(ValidationError):
        make_measurement(reducing_end_label=MISSING, reducing_end="native")


# --- sialic-acid derivatisation ---------------------------------------------


@pytest.mark.parametrize("derivatisation", SIALIC_ACID_DERIVATISATIONS)
def test_sialic_acid_derivatisation_needs_a_sialic_acid(derivatisation):
    with pytest.raises(ValidationError, match="needs a sialic acid"):
        make_measurement(derivatisation=derivatisation, glycan=make_glycan(composition="Hex5HexNAc4Fuc1"))
    for composition in ("Hex5HexNAc4Fuc1NeuAc1", "Hex5HexNAc4Fuc1NeuGc1"):
        make_measurement(derivatisation=derivatisation, glycan=make_glycan(composition=composition))


def test_the_rejection_points_to_unknown_when_something_else_may_have_reacted():
    # A 2-AA label carries a carboxyl that the same reagents can modify, so 'underivatised' may be false.
    with pytest.raises(ValidationError, match="'unknown' if something else may have reacted"):
        make_measurement(reducing_end_label="2-AA", derivatisation="sialic_acid_esterification")


@pytest.mark.parametrize("derivatisation", ["underivatised", "permethylation", "unknown"])
def test_other_derivatisations_need_no_sialic_acid(derivatisation):
    make_measurement(derivatisation=derivatisation, glycan=make_glycan(composition="Hex5HexNAc4Fuc1"))


@pytest.mark.parametrize("derivatisation", SIALIC_ACID_DERIVATISATIONS)
def test_sialic_acid_derivatisations_block_training(derivatisation):
    glycan = make_glycan(composition=SIALYLATED, reuse_status="synthetic_fixture")
    with pytest.raises(TrainingGateError, match=f"derivatisation is '{derivatisation}', an open bucket") as caught:
        assert_trainable(cleared(derivatisation=derivatisation, glycan=glycan))
    assert not isinstance(caught.value, LicenceGateError)


# --- placeholders -----------------------------------------------------------

# Every StrEnum in the package, found rather than listed, so a new enum is covered automatically.
PACKAGE_ENUMS = sorted(
    {
        obj
        for module in (wmxglycan.composition, wmxglycan.licensing, wmxglycan.models)
        for obj in vars(module).values()
        if isinstance(obj, type) and issubclass(obj, StrEnum) and obj is not StrEnum
    },
    key=lambda enum: enum.__name__,
)
# The only values spelt like a placeholder, allowed because a filler that lands on them blocks training.
PLACEHOLDER_SPELT_BUT_BLOCKING = {ReducingEndLabel.UNKNOWN, Derivatisation.UNKNOWN, UncertaintyType.UNKNOWN}


def test_the_enum_audit_sees_every_enum():
    names = {enum.__name__ for enum in PACKAGE_ENUMS}
    expected = {"Derivatisation", "DriftGas", "DTIMSMethod", "IMSType", "Polarity", "ReducingEndLabel", "Residue", "ReuseStatus"}
    assert expected <= names


@pytest.mark.parametrize("enum", PACKAGE_ENUMS, ids=lambda enum: enum.__name__)
def test_no_enum_value_is_spelt_like_a_placeholder_unless_it_blocks_training(enum):
    for member in enum:
        if _is_placeholder(member.value):
            assert member in PLACEHOLDER_SPELT_BUT_BLOCKING, f"{enum.__name__}.{member.name} = {member.value!r}"


@pytest.mark.parametrize("placeholder", sorted(_PLACEHOLDERS))
@pytest.mark.parametrize("field", ["reducing_end_label", "derivatisation"])
def test_a_placeholder_in_an_analyte_field_never_trains(field, placeholder):
    # A loader that fills empty cells with a placeholder must be rejected, or land on a value that blocks training.
    try:
        record = cleared(**{field: placeholder})
    except ValidationError:
        return
    with pytest.raises(TrainingGateError):
        assert_trainable(record)


@pytest.mark.parametrize("placeholder", sorted(_PLACEHOLDERS))
@pytest.mark.parametrize("field", ["ims_type", "drift_gas", "cell_gas", "reuse_status"])
def test_a_placeholder_in_any_other_enum_field_is_rejected(field, placeholder):
    with pytest.raises(ValidationError):
        make_measurement(**{field: placeholder})


# --- drift gas and cell gas -------------------------------------------------


def test_drift_gas_is_defined_as_the_gas_the_ccs_value_refers_to():
    description = CCSMeasurement.model_fields["drift_gas"].description
    assert "refers to" in description
    assert "not necessarily the gas in the cell" in description


def test_a_he_referenced_twims_value_does_not_pool_with_n2_values():
    # Both measured in an N2 cell; one calibrated against He reference values, so it is a He value.
    he_referenced = make_measurement(ims_type="TWIMS", cell_gas="N2", drift_gas="He").calibration_group
    n2_referenced = make_measurement(ims_type="TWIMS", cell_gas="N2", drift_gas="N2").calibration_group
    assert he_referenced != n2_referenced


def test_cell_gas_defaults_to_not_reported():
    assert make_measurement().cell_gas is None


def test_cell_gas_is_provenance_and_kept_out_of_the_calibration_group():
    assert "provenance" in CCSMeasurement.model_fields["cell_gas"].description
    groups = {make_measurement(cell_gas=gas).calibration_group for gas in (MISSING, "N2", "He")}
    assert groups == {make_measurement().calibration_group}


def test_a_primary_ccs_refers_to_the_gas_it_was_measured_in():
    with pytest.raises(ValidationError, match="must equal drift_gas"):
        make_measurement(**STEPPED, cell_gas="He", drift_gas="N2")
    assert make_measurement(**STEPPED, cell_gas="N2", drift_gas="N2").cell_gas is DriftGas.N2


# --- calibration group ------------------------------------------------------


def test_calibration_group_key():
    group = make_measurement().calibration_group
    assert group == CalibrationGroup(
        IMSType.TWIMS,
        None,
        DriftGas.N2,
        FIXTURE_CALIBRANT,
        "[M+2H]2+",
        ReducingEndLabel.NATIVE,
        Derivatisation.UNDERIVATISED,
    )
    assert str(group) == "TWIMS|N2|fixture calibrant|[M+2H]2+|native|underivatised"


def test_calibration_group_ignores_analyte_value_and_provenance():
    a = make_measurement()
    b = make_measurement(
        glycan=make_glycan(composition="Hex3HexNAc4Fuc1"),
        ccs=FIXTURE_CCS * 1.1,
        instrument="another fixture instrument",
        source="another fixture source",
    )
    assert len({a.calibration_group, b.calibration_group}) == 1


@pytest.mark.parametrize(
    "change",
    [
        {"ims_type": "TIMS"},
        {"drift_gas": "He"},
        {"calibrant": "other fixture calibrant"},
        {"adduct": "[M+H+Na]2+"},
        {"reducing_end_label": "PA"},
        {"derivatisation": "permethylation"},
        SINGLE,
    ],
    ids=["platform", "gas", "calibrant", "adduct", "reducing-end-label", "derivatisation", "dtims-method"],
)
def test_each_component_splits_the_calibration_group(change):
    assert make_measurement(**change).calibration_group != make_measurement().calibration_group


def test_pa_labelled_and_native_glycans_are_different_groups():
    pa = make_measurement(reducing_end_label="PA").calibration_group
    native = make_measurement(reducing_end_label="native").calibration_group
    assert pa != native
    assert str(pa) == "TWIMS|N2|fixture calibrant|[M+2H]2+|PA|underivatised"


def test_reduced_permethylated_and_permethylated_glycans_are_different_groups():
    reduced = make_measurement(reducing_end_label="reduced", derivatisation="permethylation").calibration_group
    free_end = make_measurement(reducing_end_label="native", derivatisation="permethylation").calibration_group
    assert reduced != free_end
    assert str(reduced) == "TWIMS|N2|fixture calibrant|[M+2H]2+|reduced|permethylation"


def test_stepped_and_single_field_dtims_are_different_groups():
    stepped = make_measurement(**STEPPED).calibration_group
    single = make_measurement(**SINGLE).calibration_group
    assert stepped != single
    assert str(stepped) == "DTIMS/stepped_field|N2|-|[M+2H]2+|native|underivatised"
    assert str(single) == "DTIMS/single_field|N2|fixture calibrant|[M+2H]2+|native|underivatised"


# --- training gate on records -----------------------------------------------


def test_academic_only_record_cannot_enter_training():
    record = make_measurement(reuse_status="academic_only", glycan=make_glycan(reuse_status="academic_only"))
    with pytest.raises(LicenceGateError, match="academic_only"):
        assert_trainable(record)

    # Trying harder: relabelling the record in place is refused, because records are frozen...
    with pytest.raises(ValidationError):
        record.reuse_status = ReuseStatus.OPEN_ATTRIBUTION
    assert record.reuse_status is ReuseStatus.ACADEMIC_ONLY
    with pytest.raises(LicenceGateError):
        assert_trainable(record)

    # ...and building it with validation skipped does not skip the gate.
    unvalidated = CCSMeasurement.model_construct(
        reuse_status="academic_only",
        glycan=make_glycan(reuse_status="synthetic_fixture"),
        reducing_end_label="native",
        derivatisation="underivatised",
    )
    with pytest.raises(LicenceGateError, match="academic_only"):
        assert_trainable(unvalidated)


def test_a_cleared_measurement_of_a_restricted_structure_is_refused():
    record = make_measurement(reuse_status="synthetic_fixture", glycan=make_glycan(reuse_status="academic_only"))
    with pytest.raises(LicenceGateError, match="GlycanStructure .* 'academic_only'"):
        assert_trainable(record)


def test_share_alike_record_cannot_enter_training():
    with pytest.raises(LicenceGateError, match="open_share_alike"):
        assert_trainable(cleared(reuse_status="open_share_alike"))
    with pytest.raises(LicenceGateError, match="open_share_alike"):
        assert_trainable(cleared(glycan=make_glycan(reuse_status="open_share_alike")))


def test_a_record_left_at_its_default_is_refused():
    with pytest.raises(LicenceGateError, match="unverified"):
        assert_trainable(make_measurement())


def test_a_fully_cleared_record_passes():
    assert_trainable(cleared())


def test_a_record_that_skipped_validation_is_judged_by_value():
    # model_construct leaves raw strings; a complete record with defined values must still train.
    fields = measurement_fields(reuse_status="synthetic_fixture", glycan=make_glycan(reuse_status="synthetic_fixture"))
    assert_trainable(CCSMeasurement.model_construct(**fields))


# --- a trainable status is a claim the gate checks; a synthetic fixture is not a claim ---


def test_an_open_claim_on_a_hand_built_record_is_refused_without_a_record():
    # No loader involved. The record says open_attribution and nothing backs it.
    with pytest.raises(LicenceGateError, match="no DOI"):
        assert_trainable(make_measurement(reuse_status="open_attribution", glycan=make_glycan(reuse_status="synthetic_fixture")))
    with pytest.raises(LicenceGateError, match="neither the row's DOI nor a registered dataset"):
        assert_trainable(cleared(glycan=make_glycan(reuse_status="open_attribution")))


def test_a_synthetic_fixture_is_trainable_and_says_what_it_is():
    record = cleared()
    assert record.reuse_status is ReuseStatus.SYNTHETIC_FIXTURE
    assert record.glycan.reuse_status is ReuseStatus.SYNTHETIC_FIXTURE
    assert_trainable(record)
    assert_trainable(cleared(doi="10.1000/synthetic-fixture"))  # an invented DOI on an invented record


def test_a_synthetic_fixture_cannot_cite_the_recorded_paper():
    from wmxglycan.sources import STRUWE_2016

    with pytest.raises(LicenceGateError, match="the record governs"):
        assert_trainable(cleared(doi=STRUWE_2016.doi))


@pytest.mark.parametrize(
    "update, message",
    [
        ({"calibrant": None}, "does not pass validation"),  # TWIMS without a calibrant
        ({"calibrant": "not reported"}, "does not pass validation"),  # a loader filler in place of a calibrant
        ({"derivatisation": "permethylated"}, "does not pass validation"),  # not an agreed value
        ({"derivatisation": "sialic_acid_esterification"}, "does not pass validation"),  # no sialic acid to act on
        ({"adduct": "[M+Na+H]2+"}, "differs from its validated form"),  # unsorted adduct, a split key
        ({"calibrant": "  fixture calibrant  "}, "differs from its validated form"),  # unstripped text
    ],
    ids=[
        "calibrant-removed",
        "calibrant-filler",
        "unagreed-derivatisation",
        "sialic-without-sialic-acid",
        "unsorted-adduct",
        "unstripped-text",
    ],
)
def test_a_record_changed_without_validation_cannot_train(update, message):
    # Records are frozen, so model_copy(update=...) is how one gets "corrected", and it skips validation.
    with pytest.raises(TrainingGateError, match=message):
        assert_trainable(cleared().model_copy(update=update))


def test_a_glycan_changed_without_validation_cannot_train():
    glycan = make_glycan(reuse_status="synthetic_fixture").model_copy(update=RESOLVED)  # resolved, with no identifier
    # Nesting it in a new measurement re-validates it...
    with pytest.raises(ValidationError, match="only a composition"):
        cleared(glycan=glycan)
    # ...but copying it into an existing one does not, so the gate has to catch it.
    with pytest.raises(TrainingGateError, match="GlycanStructure .*does not pass validation"):
        assert_trainable(cleared().model_copy(update={"glycan": glycan}))


def test_the_calibration_group_prints_for_a_record_changed_without_validation():
    record = cleared().model_copy(update={"reducing_end_label": "PA", "derivatisation": "underivatised"})
    assert str(record.calibration_group) == "TWIMS|N2|fixture calibrant|[M+2H]2+|PA|underivatised"


def test_an_unknown_label_blocks_training():
    with pytest.raises(TrainingGateError, match="reducing-end label is 'unknown'") as caught:
        assert_trainable(cleared(reducing_end_label="unknown"))
    assert not isinstance(caught.value, LicenceGateError)  # a data gap, not a licence problem


def test_an_other_label_blocks_training():
    with pytest.raises(TrainingGateError, match="reducing-end label is 'other'.*add the label") as caught:
        assert_trainable(cleared(reducing_end_label="other"))
    assert not isinstance(caught.value, LicenceGateError)


def test_an_unknown_derivatisation_blocks_training():
    with pytest.raises(TrainingGateError, match="derivatisation is 'unknown'") as caught:
        assert_trainable(cleared(derivatisation="unknown"))
    assert not isinstance(caught.value, LicenceGateError)


def test_both_analyte_gaps_are_reported_together():
    with pytest.raises(TrainingGateError) as caught:
        assert_trainable(cleared(reducing_end_label="other", derivatisation="unknown"))
    assert "reducing-end label is 'other'" in str(caught.value)
    assert "derivatisation is 'unknown'" in str(caught.value)


def test_licence_and_analyte_problems_are_reported_together():
    with pytest.raises(LicenceGateError) as caught:
        assert_trainable(cleared(reducing_end_label="unknown", reuse_status="academic_only"))
    assert "academic_only" in str(caught.value)
    assert "reducing-end label" in str(caught.value)


@pytest.mark.parametrize("label", TRAINABLE_LABELS)
@pytest.mark.parametrize("derivatisation", TRAINABLE_DERIVATISATIONS)
def test_every_defined_label_and_derivatisation_may_train(label, derivatisation):
    assert_trainable(cleared(reducing_end_label=label, derivatisation=derivatisation))


@pytest.mark.parametrize(
    "field, message", [("reducing_end_label", "reducing-end label"), ("derivatisation", "derivatisation is")]
)
@pytest.mark.parametrize(
    "value",
    [MISSING, "unknown", "not a value", ["native"], {"value": "native"}, NALike(), ArrayLike(), CollidesWithOther()],
    ids=["omitted", "unknown", "garbage", "list", "dict", "na-like", "array-like", "hash-collides-with-other"],
)
def test_an_undefined_analyte_blocks_training_even_when_validation_is_skipped(field, message, value):
    fields = {
        "reuse_status": "synthetic_fixture",
        "glycan": make_glycan(reuse_status="synthetic_fixture"),
        "reducing_end_label": "native",
        "derivatisation": "underivatised",
    }
    if value is MISSING:
        del fields[field]
    else:
        fields[field] = value
    with pytest.raises(TrainingGateError, match=message):
        assert_trainable(CCSMeasurement.model_construct(**fields))
    # A raw ValueError or TypeError would let a loader that skips malformed rows swallow a licence refusal.
    with pytest.raises(LicenceGateError, match="academic_only"):
        assert_trainable(CCSMeasurement.model_construct(**(fields | {"reuse_status": "academic_only"})))


# --- an IUPAC-condensed structure as a structure identifier ------------------


def test_a_structure_that_states_its_linkages_is_a_structure_identifier():
    glycan = make_glycan(iupac_condensed="Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc", **RESOLVED)
    assert glycan.has_unresolved_linkage is False
    assert glycan.wurcs is None and glycan.glytoucan_ac is None


@pytest.mark.parametrize("text", ["not-a-glycan", "x", "Man", "GlcNAc"])
def test_a_string_with_no_linkage_in_it_is_not_a_structure_identifier(text):
    # However structure-shaped it looks, it says nothing about how residues are joined.
    with pytest.raises(ValidationError, match="cannot claim resolved linkage"):
        make_glycan(iupac_condensed=text, **RESOLVED)


def test_a_structure_string_is_still_stored_when_it_states_no_linkage():
    assert make_glycan(iupac_condensed="Man").iupac_condensed == "Man"


@pytest.mark.parametrize("text", ["not a glycan", "   ", "N/A", ""])
def test_a_blank_or_placeholder_structure_string_is_rejected(text):
    with pytest.raises(ValidationError):
        make_glycan(iupac_condensed=text)


# --- provenance per value ----------------------------------------------------


def test_provenance_fields_default_to_not_reported():
    # Null, never a filler: "not reported" is a fact about the source, and a
    # replicate count of 1 would be an invented one.
    measurement = make_measurement()
    assert measurement.doi is None
    assert measurement.source_locator is None
    assert measurement.replicates is None
    assert measurement.measured_on is None


def test_a_recorded_doi_locator_replicates_and_date():
    measurement = make_measurement(
        doi="10.1038/s41467-021-00001-2",
        source_locator="Supplementary Data 3",
        replicates=5,
        measured_on="2021-06-30",
    )
    assert measurement.doi == "10.1038/s41467-021-00001-2"
    assert measurement.source_locator == "Supplementary Data 3"
    assert measurement.replicates == 5
    assert measurement.measured_on.year == 2021 and measurement.measured_on.month == 6


@pytest.mark.parametrize(
    "doi",
    [
        "https://doi.org/10.1038/s41467-021-00001-2",  # a URL is not a DOI
        "doi:10.1038/s41467-021-00001-2",  # nor is a prefixed citation
        "10.1038",  # no suffix
        "10.10/x",  # registrant code too short
        "Smith et al. 10.1038/s41467-021-00001-2",  # a citation containing one
        "n/a",
    ],
)
def test_a_doi_that_is_not_a_bare_doi_is_rejected(doi):
    with pytest.raises(ValidationError):
        make_measurement(doi=doi)


@pytest.mark.parametrize("replicates", [0, -1, 2.0, True, "3"])
def test_a_replicate_count_must_be_a_whole_number_of_at_least_one(replicates):
    with pytest.raises(ValidationError):
        make_measurement(replicates=replicates)


@pytest.mark.parametrize("when", ["not a date", "2021-13-01", "", "n/a"])
def test_a_measurement_date_must_be_a_date(when):
    with pytest.raises(ValidationError):
        make_measurement(measured_on=when)


@pytest.mark.parametrize(
    "when",
    [
        "1462060800",  # pydantic would read this as a Unix timestamp: 2016-05-01
        "1462060800.0",
        1462060800,
        "2016-05-01T00:00:00",  # a midnight datetime is not how a date is written
        "20160501",
        "2016-5-1",
    ],
)
def test_a_measurement_date_is_only_ever_written_as_a_date(when):
    # A cell that held a count must not become a plausible date with nothing to show it.
    with pytest.raises(ValidationError, match="YYYY-MM-DD"):
        make_measurement(measured_on=when)


def test_a_datetime_is_not_a_measurement_date():
    from datetime import date, datetime

    with pytest.raises(ValidationError, match="not a datetime"):
        make_measurement(measured_on=datetime(2016, 5, 1, 0, 0))
    assert make_measurement(measured_on=date(2016, 5, 1)).measured_on == date(2016, 5, 1)
    assert make_measurement(measured_on="2016-05-01").measured_on == date(2016, 5, 1)


def test_provenance_is_not_part_of_the_calibration_group():
    # Two measurements of one analyte from different papers still pool: what a
    # value means is fixed by the conditions, not by who reported it.
    a = make_measurement(doi="10.1038/s41467-021-00001-2", replicates=3, measured_on="2021-06-30")
    b = make_measurement(doi="10.1000/other-paper", source_locator="Table S2")
    assert a.calibration_group == b.calibration_group == make_measurement().calibration_group


@pytest.mark.parametrize("filler", LOADER_FILLERS)
def test_a_loader_filler_is_not_a_locator(filler):
    with pytest.raises(ValidationError, match="placeholder"):
        make_measurement(source_locator=filler)


# --- the uncertainty travels with its type ------------------------------------


def test_uncertainty_types_are_exactly_the_agreed_list():
    assert {kind.value for kind in UncertaintyType} == {"sd", "two_sd", "sem", "ci95", "unknown"}


def test_an_uncertainty_is_stored_with_its_type():
    measurement = make_measurement(ccs_uncertainty=2.4, uncertainty_type="two_sd")
    assert measurement.ccs_uncertainty == 2.4
    assert measurement.uncertainty_type is UncertaintyType.TWO_SD


def test_an_uncertainty_without_a_type_is_rejected():
    # A spread of unknown kind is read as whatever the reader assumes. Loading a
    # two-standard-deviation spread into a field read as one halves every interval
    # built on it, and nothing downstream could tell. So the type is mandatory.
    with pytest.raises(ValidationError, match="needs uncertainty_type"):
        make_measurement(ccs_uncertainty=2.4)


def test_a_type_without_an_uncertainty_is_rejected():
    with pytest.raises(ValidationError, match="not given"):
        make_measurement(uncertainty_type="sd")


@pytest.mark.parametrize("value", [0.0, -1.0, float("nan"), float("inf"), True, "2.4"])
def test_an_uncertainty_must_be_a_positive_finite_number(value):
    with pytest.raises(ValidationError):
        make_measurement(ccs_uncertainty=value, uncertainty_type="sd")


@pytest.mark.parametrize("kind", ["2sd", "SD", "stdev", "95ci", "other", ""])
def test_unagreed_uncertainty_types_are_rejected(kind):
    with pytest.raises(ValidationError):
        make_measurement(ccs_uncertainty=1.0, uncertainty_type=kind)


def test_an_unknown_uncertainty_type_is_stored_but_blocks_training():
    record = cleared(ccs_uncertainty=1.0, uncertainty_type="unknown")
    assert record.uncertainty_type is UncertaintyType.UNKNOWN
    with pytest.raises(TrainingGateError, match="uncertainty_type is 'unknown'") as caught:
        assert_trainable(record)
    assert not isinstance(caught.value, LicenceGateError)  # a data gap, not a licence problem


@pytest.mark.parametrize("kind", ["sd", "two_sd", "sem", "ci95"])
def test_every_stated_uncertainty_type_may_train(kind):
    assert_trainable(cleared(ccs_uncertainty=1.0, uncertainty_type=kind))


def test_no_uncertainty_at_all_is_not_a_blocker():
    # Not reported is not the same as reported with an unknown kind.
    assert_trainable(cleared())


def test_the_uncertainty_is_not_part_of_the_calibration_group():
    a = make_measurement(ccs_uncertainty=0.4, uncertainty_type="two_sd")
    b = make_measurement(ccs_uncertainty=3.0, uncertainty_type="sd")
    assert a.calibration_group == b.calibration_group == make_measurement().calibration_group


def test_unknown_may_stand_alone_because_it_describes_no_number():
    # It records that the source states no spread, which is a fact worth keeping
    # and which blocks training. Every other type describes a number.
    record = cleared(uncertainty_type="unknown")
    assert record.uncertainty_type is UncertaintyType.UNKNOWN and record.ccs_uncertainty is None
    with pytest.raises(TrainingGateError, match="uncertainty_type is 'unknown'"):
        assert_trainable(record)


@pytest.mark.parametrize("kind", ["sd", "two_sd", "sem", "ci95"])
def test_every_type_that_describes_a_number_still_needs_one(kind):
    with pytest.raises(ValidationError, match="only 'unknown' may stand alone"):
        make_measurement(uncertainty_type=kind)


# --- an unstated drift gas is a fact, and it blocks -----------------------------


def test_an_unstated_drift_gas_is_stored_and_blocks_training():
    # The source did not say which gas its values refer to. Recorded, not guessed.
    record = cleared(drift_gas="UNSTATED", cell_gas="UNSTATED")
    assert record.drift_gas is DriftGas.UNSTATED
    with pytest.raises(TrainingGateError, match="drift_gas is 'UNSTATED'") as caught:
        assert_trainable(record)
    assert "calibration reference" in str(caught.value)
    assert not isinstance(caught.value, LicenceGateError)  # a data gap, not a licence problem


def test_an_unstated_gas_is_not_spelt_like_a_placeholder():
    # Deliberately not "unknown", which this package refuses in every enum field.
    assert not _is_placeholder(DriftGas.UNSTATED.value)
    with pytest.raises(ValidationError):
        make_measurement(drift_gas="unknown")


def test_an_unstated_drift_gas_keeps_its_own_calibration_group():
    # It must not pool with helium or with nitrogen: the quantity is undefined.
    unstated = make_measurement(drift_gas="UNSTATED").calibration_group
    assert unstated != make_measurement(drift_gas="He").calibration_group
    assert unstated != make_measurement(drift_gas="N2").calibration_group
    assert "UNSTATED" in str(unstated)


def test_an_unreported_cell_gas_is_not_stricter_than_an_unstated_one():
    # cell_gas is nullable and null does not block, so an explicit "the paper did
    # not say" must not block either; drift_gas is where the quantity is defined.
    assert_trainable(cleared(cell_gas="UNSTATED"))
    assert_trainable(cleared())


# --- conformers: one structure, more than one peak ------------------------------


def test_a_conformer_pair_is_stored_with_its_index_and_total():
    record = cleared(conformer=2, conformers_total=2)
    assert record.conformer == 2 and record.conformers_total == 2


def test_a_single_conformer_needs_no_index():
    assert cleared(conformers_total=1).conformer is None


def test_an_index_without_a_total_is_rejected():
    with pytest.raises(ValidationError, match="says nothing without conformers_total"):
        make_measurement(conformer=1)


def test_a_total_above_one_demands_an_index():
    with pytest.raises(ValidationError, match="must say which of them it is"):
        make_measurement(conformers_total=2)


def test_an_index_beyond_the_total_is_rejected():
    with pytest.raises(ValidationError, match="exceeds the total"):
        make_measurement(conformer=3, conformers_total=2)


@pytest.mark.parametrize("value", [0, -1, 1.0, True, "1"])
def test_a_conformer_index_must_be_a_whole_number_of_at_least_one(value):
    with pytest.raises(ValidationError):
        make_measurement(conformer=value, conformers_total=2)


def test_conformers_do_not_split_the_calibration_group_or_block_training():
    # Two peaks of one structure under one set of conditions are one group.
    a = cleared(conformer=1, conformers_total=2)
    b = cleared(conformer=2, conformers_total=2, ccs=FIXTURE_CCS + 10.0)
    assert a.calibration_group == b.calibration_group == cleared().calibration_group
    assert_trainable(a)
    assert_trainable(b)


# --- the calibrant reference is provenance for the calibration curve -----------


def test_a_calibrant_reference_is_stored_for_a_calibrated_value():
    measurement = make_measurement(calibrant_reference="helium CCS reference values")
    assert measurement.calibrant_reference == "helium CCS reference values"


def test_a_primary_value_may_not_carry_a_calibrant_reference():
    # Stepped-field DTIMS has no calibrant, so nothing was referenced against anything.
    with pytest.raises(ValidationError, match="cannot have a calibrant reference"):
        make_measurement(**STEPPED, calibrant_reference="helium CCS reference values")


def test_the_calibrant_reference_is_kept_out_of_the_calibration_group():
    # It records where the curve came from; drift_gas already says which gas the value refers to.
    a = make_measurement(calibrant_reference="helium CCS reference values")
    b = make_measurement(calibrant_reference="nitrogen CCS reference values")
    assert a.calibration_group == b.calibration_group == make_measurement().calibration_group


@pytest.mark.parametrize("filler", LOADER_FILLERS)
def test_a_loader_filler_is_not_a_calibrant_reference(filler):
    with pytest.raises(ValidationError, match="placeholder"):
        make_measurement(calibrant_reference=filler)


# --- adduct components ------------------------------------------------------


@pytest.mark.parametrize(
    "adduct, parts",
    [
        ("[M+H]+", (("+", 1, "H"),)),
        ("[M+2H]2+", (("+", 2, "H"),)),
        ("[M+Na+H]2+", (("+", 1, "H"), ("+", 1, "Na"))),  # canonical order, additions by species
        ("[M-H]-", (("-", 1, "H"),)),
        ("[M-H2O+H]+", (("+", 1, "H"), ("-", 1, "H2O"))),  # additions before losses
        ("[M+NH4]+", (("+", 1, "NH4"),)),
        ("[M+10H]10+", (("+", 10, "H"),)),
        ("[2M+Na]+", (("+", 1, "Na"),)),  # the multiplier is not a component
    ],
)
def test_adduct_components_are_the_parts_the_canonical_spelling_discards(adduct, parts):
    assert adduct_components(adduct) == parts


def test_two_spellings_of_one_ion_give_one_component_sequence():
    assert adduct_components("[M+Na+H]2+") == adduct_components("[M+H+Na]2+")


def test_a_sodiated_ion_is_not_a_protonated_one_of_the_same_charge():
    # The whole reason the parts are needed: these two share a charge and differ
    # in mass, so a featuriser that saw only the charge would call them equal.
    assert adduct_components("[M+2H]2+") != adduct_components("[M+H+Na]2+")


@pytest.mark.parametrize("adduct", ["M+2H", "[M+2H]", "[M+h]+", "[M+H+H]2+", ""])
def test_a_malformed_adduct_has_no_components(adduct):
    with pytest.raises(ValueError):
        adduct_components(adduct)


@pytest.mark.parametrize("adduct", ["[M+H]+", "[M+2H]2+", "[M+Na+H]2+", "[M-H]-", "[M-H2O+H]+", "[M+NH4]+", "[2M+Na]+"])
def test_the_components_reconstruct_the_canonical_spelling(adduct):
    # The real invariant. Counting additions against losses is not one: a
    # neutral loss such as water changes the mass without changing the charge,
    # and losing a proton makes the charge negative rather than smaller.
    canonical, _ = wmxglycan.models._parse_adduct(adduct)
    body = "".join(
        f"{sign}{'' if count == 1 else count}{species}" for sign, count, species in adduct_components(adduct)
    )
    assert f"M{body}]" in canonical


@pytest.mark.parametrize("adduct, charge", [("[M+H]+", 1), ("[M+2H]2+", 2), ("[M+10H]10+", 10)])
def test_a_purely_protonated_ion_carries_its_charge_in_protons(adduct, charge):
    protons = sum(count for sign, count, species in adduct_components(adduct) if sign == "+" and species == "H")
    assert protons == charge
