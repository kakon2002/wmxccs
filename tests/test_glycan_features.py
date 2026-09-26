"""Turning a record into a row of numbers, and what the row refuses to say.

Every record here is a synthetic fixture. No CCS value in this file is a real
measurement, and none is ever read by the extractor under test.
"""

import inspect
import math
import re

import pytest

import wmxglycan.features as features_module
from wmxglycan.contracts import PredictRequest
from wmxglycan.features import (
    CATEGORICAL_FEATURES,
    DEPLOYMENT_SUPPLIED_FEATURES,
    EXCLUDED,
    FEATURE_NAMES,
    FITTED_ANALYTE_FEATURES,
    FITTED_CONDITION_FEATURES,
    FITTED_FEATURES,
    REQUIRES_STRUCTURE,
    STRUCTURE_FEATURES,
    UNSUPPORTED_RESIDUES,
    ExtractionReport,
    FeatureVector,
    extract,
    extract_all,
    extract_structure,
    has_structure,
    is_fully_resolved,
    resolution_gaps,
)
from wmxglycan.glycan_graph import GlycanGraph
from wmxglycan.models import CCSMeasurement, Derivatisation, GlycanStructure, IMSType, ReducingEndLabel

FIXTURE_SOURCE = "synthetic test fixture, not a real record"
FIXTURE_CCS = 100.0  # arbitrary; the extractor must never look at it

G0F = "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
# Same composition as G0F's galactosylated form, differing only in which arm carries the galactose.
G1F_A3 = "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
G1F_A6 = "GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
FLOATING_FUC = "{Fuc(a1-3)}GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
SIALYL_AMBIGUOUS = "Neu5Ac(a2-3/6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
SIALYL_A2_3 = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
NO_ANOMER = "Man(1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"

RESOLVED = {"has_unresolved_linkage": False, "has_unresolved_anomericity": False}

PREDICT_BODY = {
    "composition": "Hex5HexNAc4Fuc1",
    "adduct": "[M+2H]2+",
    "charge": 2,
    "polarity": "positive",
    "reducing_end_label": "native",
    "derivatisation": "underivatised",
    "ims_type": "TWIMS",
    "drift_gas": "N2",
    "calibrant": "fixture calibrant",
}


def a_glycan(composition="Hex5HexNAc4Fuc1", **overrides):
    return GlycanStructure(**({"composition": composition, "source": FIXTURE_SOURCE} | overrides))


def a_measurement(**overrides):
    fields = {
        "glycan": a_glycan(),
        "ccs": FIXTURE_CCS,
        "adduct": "[M+2H]2+",
        "charge": 2,
        "polarity": "positive",
        "reducing_end_label": "native",
        "derivatisation": "underivatised",
        "ims_type": "TWIMS",
        "drift_gas": "N2",
        "calibrant": "fixture calibrant",
        "source": FIXTURE_SOURCE,
    }
    return CCSMeasurement(**(fields | overrides))


def with_structure(iupac, composition="Hex5HexNAc4Fuc1"):
    return a_measurement(glycan=a_glycan(composition, iupac_condensed=iupac, **RESOLVED))


def value_of(vector, name):
    return vector.as_mapping()[name]


# --- the shape of the row never moves ----------------------------------------


def test_every_feature_name_is_classified():
    assert FEATURE_NAMES == FITTED_ANALYTE_FEATURES + FITTED_CONDITION_FEATURES + STRUCTURE_FEATURES
    assert len(FEATURE_NAMES) == len(set(FEATURE_NAMES))  # no name appears twice
    assert FITTED_FEATURES == FEATURE_NAMES  # all three blocks are fitted
    assert DEPLOYMENT_SUPPLIED_FEATURES == FITTED_ANALYTE_FEATURES + FITTED_CONDITION_FEATURES
    assert REQUIRES_STRUCTURE <= set(FEATURE_NAMES)
    assert CATEGORICAL_FEATURES <= set(FITTED_CONDITION_FEATURES)


def test_the_row_is_pinned_name_by_name():
    # Pinned literally, not derived. A width check that reads len(FEATURE_NAMES)
    # follows the constant it is meant to police, so dropping a feature would
    # narrow the row from 44 to 43 in silence.
    assert len(FEATURE_NAMES) == 44
    assert FITTED_ANALYTE_FEATURES == (
        "hex",
        "hexnac",
        "fuc",
        "neuac",
        "neugc",
        "residue_total",
        "monoisotopic_mass",
        "antenna_hexnac",
        "has_full_core",
    )
    assert FITTED_CONDITION_FEATURES == (
        "charge_magnitude",
        "polarity_positive",
        "adduct_protons",
        "adduct_sodium",
        "adduct_other_additions",
        "adduct_losses",
        "ccs_is_calibrated",
        "drift_gas_he",
        "ims_type",
        "reducing_end_label",
        "derivatisation",
    )
    assert STRUCTURE_FEATURES == (
        "residues_attached",
        "unplaced_residues",
        "bonds",
        "max_depth",
        "mean_depth",
        "branch_points",
        "max_out_degree",
        "terminal_residues",
        "wiener_index",
        "mean_path_length",
        "arm_a3_size",
        "arm_a6_size",
        "arm_asymmetry",
        "core_fucose",
        "bisecting_glcnac",
        "linkages_total",
        "alpha_linkages",
        "beta_linkages",
        "anomer_unknown",
        "position_unknown",
        "flexible_1_6_linkages",
        "sia_a2_3",
        "sia_a2_6",
        "sia_position_ambiguous",
    )
    assert (len(FITTED_ANALYTE_FEATURES), len(FITTED_CONDITION_FEATURES), len(STRUCTURE_FEATURES)) == (9, 11, 24)


@pytest.mark.parametrize(
    "record",
    [
        a_measurement(),  # composition only
        with_structure(G0F, "Hex3HexNAc4Fuc1"),
        with_structure(FLOATING_FUC, "Hex3HexNAc3Fuc1"),
        with_structure(SIALYL_AMBIGUOUS, "Hex4HexNAc3NeuAc1"),
    ],
    ids=["composition-only", "resolved", "floating", "ambiguous"],
)
def test_the_vector_width_and_order_never_change(record):
    vector = extract(record)
    assert len(vector.values) == len(FEATURE_NAMES)
    assert list(vector.as_mapping()) == list(FEATURE_NAMES)


def test_a_vector_of_the_wrong_width_is_refused():
    with pytest.raises(ValueError, match="feature vector holds"):
        FeatureVector(values=(1.0, 2.0))


def test_the_categorical_codes_cover_every_enum_member():
    # An upstream enum addition must fail here rather than silently take a NaN.
    for enum, name in [
        (IMSType, "ims_type"),
        (ReducingEndLabel, "reducing_end_label"),
        (Derivatisation, "derivatisation"),
    ]:
        order = getattr(features_module, f"_{ {'ims_type': 'IMS', 'reducing_end_label': 'LABEL', 'derivatisation': 'DERIVATISATION'}[name] }_ORDER")
        assert set(order) == set(enum), name
        assert len(order) == len(set(order)), name


# --- missing is None, and never zero -----------------------------------------


def test_a_missing_value_is_none_and_never_zero():
    # The single most important property here: zero means the extractor looked
    # and found none; None means it could not look.
    vector = extract(a_measurement())  # no structure at all
    mapping = vector.as_mapping()
    for name in STRUCTURE_FEATURES:
        assert mapping[name] is None, name
        assert mapping[name] != 0 and mapping[name] != 0.0, name

    row = vector.as_row(FEATURE_NAMES)
    for index, name in enumerate(FEATURE_NAMES):
        if name in REQUIRES_STRUCTURE:
            assert math.isnan(row[index]), name  # NaN only at the wire
        else:
            assert not math.isnan(row[index]), name

    mask = vector.known_mask()
    known = {name for name, ok in zip(FEATURE_NAMES, mask) if ok}
    assert known == set(FEATURE_NAMES) - REQUIRES_STRUCTURE


def test_a_zero_count_is_a_real_zero():
    # G0F carries no sialic acid: that is a measured zero, not an absence.
    assert value_of(extract(with_structure(G0F, "Hex3HexNAc4Fuc1")), "sia_a2_3") == 0.0
    assert value_of(extract(a_measurement()), "sia_a2_3") is None


# --- the fitted row is exactly what a prediction request can supply -----------


def test_the_whole_row_is_fitted_including_the_structure_block():
    # The model scores an enumerated candidate, not the request, so a structure
    # is available at inference. Without these columns nothing separates two
    # isomers of one composition and a ranking layer has nothing to rank on.
    assert FITTED_FEATURES == FEATURE_NAMES
    assert REQUIRES_STRUCTURE <= set(FITTED_FEATURES)


def test_a_request_alone_is_no_longer_a_complete_model_input():
    # It becomes one once a candidate has been enumerated for it. The
    # deployment-supplied block is exactly the part a request can fill by itself.
    vector = extract(PredictRequest(**PREDICT_BODY))
    for name, value in zip(DEPLOYMENT_SUPPLIED_FEATURES, vector.as_row(DEPLOYMENT_SUPPLIED_FEATURES)):
        assert not math.isnan(value), name
    # ...and the structure half is still absent, which is why a candidate is needed.
    whole = dict(zip(FEATURE_NAMES, vector.as_row(FEATURE_NAMES)))
    assert all(math.isnan(whole[name]) for name in REQUIRES_STRUCTURE)


def test_a_candidate_structure_completes_the_fitted_row():
    vector = extract(with_structure(G0F, "Hex3HexNAc4Fuc1"))
    for name, value in zip(FITTED_FEATURES, vector.as_row(FITTED_FEATURES)):
        assert not math.isnan(value), name


def test_a_composition_only_record_cannot_teach_a_structure_feature():
    assert not has_structure(a_measurement())
    assert has_structure(with_structure(G0F, "Hex3HexNAc4Fuc1"))


def test_a_request_and_a_measurement_of_one_analyte_agree_on_what_a_request_supplies():
    request = PredictRequest(**PREDICT_BODY)
    measurement = a_measurement()
    supplied = DEPLOYMENT_SUPPLIED_FEATURES
    assert extract(request).as_row(supplied) == extract(measurement).as_row(supplied)


def test_the_analyte_row_is_the_molecule_and_not_the_conditions():
    # Two measurements of one glycan under different adducts: same analyte row,
    # different fitted row. This is what the leakage law is checked over.
    one = extract(a_measurement())
    two = extract(a_measurement(adduct="[M+H+Na]2+"))
    assert one.analyte_row() == two.analyte_row()
    assert one.as_row(FITTED_FEATURES) != two.as_row(FITTED_FEATURES)


def test_two_isomers_of_one_composition_share_an_analyte_row():
    # The reason the split cannot be keyed on structure: the model cannot tell
    # these apart from the fitted row, so a split that separates them leaks.
    a3 = extract(with_structure(G1F_A3, "Hex4HexNAc4Fuc1"))
    a6 = extract(with_structure(G1F_A6, "Hex4HexNAc4Fuc1"))
    assert a3.analyte_row() == a6.analyte_row()
    # ...and yet they really are two molecules, which the structure block sees.
    assert value_of(a3, "arm_a3_size") != value_of(a6, "arm_a3_size")


# --- no record influences another --------------------------------------------


def test_a_feature_vector_does_not_depend_on_the_other_records():
    # What forbids target encoding, frequency encoding and any fitted transform,
    # and what licenses featurising before the split.
    alone = extract(a_measurement())
    batch = extract_all([a_measurement(), a_measurement(ccs=222.0), a_measurement(adduct="[M+H]+", charge=1)])
    assert batch.vectors[0].values == alone.values


# --- the target and the provenance are never read -----------------------------


class ExplodingCCS:
    """Carries every condition a featuriser needs, and detonates if the value is read."""

    glycan = None
    adduct = "[M+2H]2+"
    charge = 2
    polarity = "positive"
    reducing_end_label = ReducingEndLabel.NATIVE
    derivatisation = Derivatisation.UNDERIVATISED
    ims_type = IMSType.TWIMS
    drift_gas = "N2"
    ccs_is_calibrated = True
    composition = None

    def __init__(self, composition):
        self.composition = composition

    @property
    def ccs(self):
        raise AssertionError("the extractor read the measured value")


def test_the_extractor_never_reads_the_measured_value():
    # Trying harder than a name check: this fails loudly rather than silently.
    vector = extract(ExplodingCCS(a_glycan().composition))
    assert value_of(vector, "hex") == 5.0


BANNED_READS = ["ccs", "source", "instrument", "cell_gas", "reuse_status", "glytoucan_ac", "wurcs"]


def test_no_excluded_field_reaches_the_features():
    source = inspect.getsource(features_module)
    for name in BANNED_READS:
        # A word boundary, so the legitimate derived flag ccs_is_calibrated is
        # not caught while a read of the target itself would be.
        assert not re.search(rf"\.{name}\b", source), f"features.py reads .{name}"


def test_every_banned_field_carries_the_reason_it_is_banned():
    for name in BANNED_READS + ["iupac_condensed", "calibrant"]:
        assert name in EXCLUDED, name
        assert EXCLUDED[name].strip(), name
        assert name not in FEATURE_NAMES, name


# --- unresolved stays unresolved ----------------------------------------------


def test_an_ambiguous_sialic_linkage_is_not_resolved_to_either_position():
    vector = extract(with_structure(SIALYL_AMBIGUOUS, "Hex4HexNAc3NeuAc1"))
    assert value_of(vector, "sia_a2_3") == 0.0
    assert value_of(vector, "sia_a2_6") == 0.0
    assert value_of(vector, "sia_position_ambiguous") == 1.0
    assert value_of(vector, "position_unknown") >= 1.0
    # A resolved one is counted where it belongs, so the fixture is discriminating.
    resolved = extract(with_structure(SIALYL_A2_3, "Hex4HexNAc3NeuAc1"))
    assert value_of(resolved, "sia_a2_3") == 1.0
    assert value_of(resolved, "sia_position_ambiguous") == 0.0


SIA_BETA = "Neu5Ac(b2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
SIA_NO_ANOMER = "Neu5Ac(2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"


@pytest.mark.parametrize("iupac", [SIA_BETA, SIA_NO_ANOMER], ids=["stated-beta", "unstated-anomer"])
def test_a_sialic_acid_that_is_not_alpha_is_never_counted_in_an_alpha_column(iupac):
    # sia_a2_3 names an alpha linkage. A stated beta, or an anomer nobody wrote
    # down, is not one, and folding either into it would be inventing chemistry.
    vector = extract(with_structure(iupac, "Hex4HexNAc3NeuAc1"))
    assert value_of(vector, "sia_a2_3") == 0.0
    assert value_of(vector, "sia_a2_6") == 0.0
    assert value_of(vector, "sia_position_ambiguous") == 1.0


def test_a_flexible_linkage_needs_a_known_position():
    # "a1-6/?" says the bond might be 1-6, not that it is.
    maybe = "Man(a1-3)[Man(a1-6/3)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
    vector = extract(with_structure(maybe, "Hex3HexNAc2"))
    assert value_of(vector, "position_unknown") >= 1.0
    assert value_of(vector, "flexible_1_6_linkages") == 0.0


def test_an_unknown_anomer_is_counted_as_unknown():
    vector = extract(with_structure(NO_ANOMER, "Hex3HexNAc2"))
    total = value_of(vector, "linkages_total")
    alpha = value_of(vector, "alpha_linkages")
    beta = value_of(vector, "beta_linkages")
    unknown = value_of(vector, "anomer_unknown")
    assert unknown >= 1.0
    # The invariant is what stops a future edit folding the unknown into either bucket.
    assert alpha + beta + unknown == total


def test_a_resolution_gap_is_reported_on_the_vector():
    assert extract(with_structure(SIALYL_AMBIGUOUS, "Hex4HexNAc3NeuAc1")).resolution_gaps
    assert not extract(with_structure(G0F, "Hex3HexNAc4Fuc1")).resolution_gaps
    assert is_fully_resolved(GlycanGraph.from_iupac_condensed(G0F))
    assert not is_fully_resolved(GlycanGraph.from_iupac_condensed(SIALYL_AMBIGUOUS))
    assert resolution_gaps(GlycanGraph.from_iupac_condensed(FLOATING_FUC))


# --- a motif is a placement claim, not a substring -----------------------------


def test_a_floating_residue_counts_toward_the_analyte_but_never_toward_shape():
    vector = extract(with_structure(FLOATING_FUC, "Hex3HexNAc3Fuc1"))
    # Seven residues in the string, one of them floating, so six are placed.
    assert value_of(vector, "residue_total") == 7.0  # the analyte is all seven
    assert value_of(vector, "unplaced_residues") == 1.0
    assert value_of(vector, "residues_attached") == 6.0
    # An a1-3 fucose of unknown attachment is not core fucose.
    assert value_of(vector, "core_fucose") == 0.0
    placed = extract(with_structure(G0F, "Hex3HexNAc4Fuc1"))
    assert value_of(placed, "core_fucose") == 1.0
    assert value_of(placed, "unplaced_residues") == 0.0


ANTENNARY_FUC = "Fuc(a1-6)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
FLOATING_PAIR = "{Fuc(a1-6)GlcNAc(b1-4)}Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"


def test_core_fucose_means_the_reducing_end_and_not_any_glcnac():
    # A fragment match cannot express this: "Fuc(a1-6)GlcNAc" leaves the
    # innermost GlcNAc with no arrival linkage, so it is satisfied by an
    # antennary one. The flag is anchored to the reducing end instead.
    antennary = extract(with_structure(ANTENNARY_FUC, "Hex3HexNAc3Fuc1"))
    assert value_of(antennary, "core_fucose") == 0.0
    assert value_of(extract(with_structure(G0F, "Hex3HexNAc4Fuc1")), "core_fucose") == 1.0


def test_a_multi_residue_brace_fragment_cannot_satisfy_a_placement_flag():
    # The single-node floating fixture has no internal edge and so cannot catch
    # this; a two-residue fragment carries the bond a substring match would find.
    vector = extract(with_structure(FLOATING_PAIR, "Hex3HexNAc3Fuc1"))
    assert value_of(vector, "unplaced_residues") == 2.0
    assert value_of(vector, "core_fucose") == 0.0


def test_the_bisecting_flag_needs_the_bond_not_the_residue():
    assert value_of(extract(with_structure(G0F, "Hex3HexNAc4Fuc1")), "bisecting_glcnac") == 0.0
    bisected = "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)][GlcNAc(b1-4)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
    assert value_of(extract(with_structure(bisected, "Hex3HexNAc5")), "bisecting_glcnac") == 1.0


def test_the_arm_features_are_none_when_there_is_no_branched_core():
    # None, not zero: a linear fragment has no arms to measure, which is not the
    # same as having arms of size zero.
    vector = extract(with_structure("Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc", "Hex2HexNAc2"))
    assert value_of(vector, "arm_a3_size") is None
    assert value_of(vector, "arm_asymmetry") is None


# --- residues the composition vocabulary cannot express ------------------------


def test_a_residue_outside_the_five_classes_is_reported_not_rounded_off():
    sulfated = "GlcNAc6S(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
    record = with_structure(sulfated, "Hex3HexNAc3")
    vector = extract(record)
    # Block A still comes from the record's own validated composition...
    assert value_of(vector, "hexnac") == 3.0
    # ...and the disagreement is reported rather than absorbed.
    assert UNSUPPORTED_RESIDUES in vector.notes
    report = extract_all([record])
    assert report.failure_counts.get(UNSUPPORTED_RESIDUES) == 1


# --- the batch report accounts for every row -----------------------------------


def test_every_row_is_accounted_for():
    records = [a_measurement(), with_structure(G0F, "Hex3HexNAc4Fuc1"), with_structure(SIALYL_A2_3, "Hex4HexNAc3NeuAc1")]
    report = extract_all(records)
    assert report.rows_in == 3
    assert report.vectors_built == 3
    assert report.records_with_structure == 2
    assert report.structures_parsed == 2
    assert "rows in" in report.summary()


def test_an_empty_batch_is_a_caller_mistake_not_a_clean_run():
    with pytest.raises(ValueError, match="no records"):
        extract_all([])


def test_the_report_is_frozen_and_counts_are_readable():
    report = extract_all([a_measurement()])
    assert isinstance(report, ExtractionReport)
    assert report.rows_failed == 0
