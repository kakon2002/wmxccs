"""Loading SugarBase from the installed glycowork package into licence-tagged records."""

from importlib.metadata import version

import pytest

from wmxglycan.licensing import ReuseStatus, assert_trainable
from wmxglycan.sugarbase import (
    COMPOSITION_DISAGREES,
    PARSE_FAILED,
    UNSUPPORTED_RESIDUES,
    dataset_version,
    load_sugarbase,
)

M3 = "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
ROW = {"glycan": M3, "glycan_type": "N", "glytoucan_id": "G00000AA", "Composition": "{'Hex': 3, 'HexNAc': 2}"}


def row(**overrides):
    return {**ROW, **overrides}


def test_the_dataset_version_names_the_installed_package():
    dataset = dataset_version()
    assert dataset.package_version == version("glycowork")
    assert str(dataset) == f"SugarBase v12 via glycowork {version('glycowork')}"
    assert dataset.file in dataset.provenance
    assert dataset.licence == "MIT"
    assert "Daniel Bojar" in dataset.attribution


def test_a_row_becomes_a_record():
    report = load_sugarbase(rows=[row()])
    (structure,) = report.structures
    assert structure.composition.canonical == "Hex3HexNAc2"
    assert structure.iupac_condensed == M3
    assert structure.glytoucan_ac == "G00000AA"
    assert structure.has_unresolved_linkage is False
    assert structure.has_unresolved_anomericity is False
    assert structure.reuse_status is ReuseStatus.OPEN_ATTRIBUTION
    assert str(report.dataset_version) in structure.source
    assert report.dataset_version.file in structure.source


def test_a_record_may_enter_training_because_the_licence_permits_it():
    (structure,) = load_sugarbase(rows=[row()]).structures
    assert_trainable(structure)  # MIT: open attribution


@pytest.mark.parametrize("missing", [None, "", "   ", 42, float("nan")])
def test_a_row_without_an_accession_still_becomes_a_record(missing):
    report = load_sugarbase(rows=[row(glytoucan_id=missing)])
    (structure,) = report.structures
    assert structure.glytoucan_ac is None
    assert report.records_with_accession == 0


def test_unresolved_linkages_are_carried_into_the_record():
    (structure,) = load_sugarbase(
        rows=[row(glycan="Man(a1-?)Man(b1-4)GlcNAc", Composition="{'Hex': 2, 'HexNAc': 1}")]
    ).structures
    assert structure.has_unresolved_linkage is True
    assert structure.has_unresolved_anomericity is False


def test_residues_we_cannot_represent_are_reported_not_rounded_off():
    report = load_sugarbase(rows=[row(glycan="Xyl(b1-2)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc")])
    assert report.structures == ()
    assert report.failure_counts == {UNSUPPORTED_RESIDUES: 1}
    assert "Xyl" in report.failure_examples[UNSUPPORTED_RESIDUES][0]


def test_a_composition_that_disagrees_with_sugarbase_is_reported():
    report = load_sugarbase(rows=[row(Composition="{'Hex': 9, 'HexNAc': 2}")])
    assert report.structures == ()
    assert report.failure_counts == {COMPOSITION_DISAGREES: 1}


def test_a_structure_that_will_not_parse_is_reported():
    report = load_sugarbase(rows=[row(glycan="")])
    assert report.structures == ()
    assert report.failure_counts == {PARSE_FAILED: 1}


def test_rows_are_filtered_by_glycan_type_and_the_rest_are_still_counted():
    rows = [row(), row(glycan_type="O"), row(glycan_type=None)]
    n_linked = load_sugarbase(rows=rows)
    assert n_linked.rows_in_dataset == 3  # every row was read
    assert n_linked.rows_selected == 1  # only the N-linked one was accounted for
    everything = load_sugarbase(glycan_types=None, rows=rows)
    assert everything.rows_in_dataset == 3
    assert everything.rows_selected == 3
    assert everything.rows_n_linked == 1
    assert everything.records_n_linked == 1
    assert everything.records_built == 3


def test_limit_stops_early():
    assert load_sugarbase(rows=[row()] * 5, limit=2).records_built == 2


def test_every_selected_row_is_accounted_for():
    rows = [row(), row(glycan="Xyl(b1-2)Man(b1-4)GlcNAc"), row(glycan=""), row(Composition="{'Hex': 9}")]
    report = load_sugarbase(rows=rows)
    assert report.rows_selected == 4
    assert report.records_built + report.rows_failed == report.rows_selected
    assert report.graphs_parsed == 3  # the empty string never became a graph


def test_the_summary_states_the_licence_the_attribution_and_both_row_counts():
    summary = load_sugarbase(rows=[row(), row(glycan_type="O")]).summary()
    assert "MIT" in summary and "Daniel Bojar" in summary
    assert "rows in the dataset      2" in summary
    assert "selected by type         1" in summary
    assert "records built" in summary


# --- against the SugarBase release that actually ships -----------------------


@pytest.fixture(scope="module")
def n_linked():
    return load_sugarbase(glycan_types=("N",))


def test_the_shipped_n_linked_load(n_linked):
    # Counts are tied to the pinned glycowork release; they move only when the pin moves.
    assert n_linked.rows_in_dataset == 50461
    assert n_linked.rows_selected == 15712
    assert n_linked.rows_n_linked == 15712
    assert n_linked.graphs_parsed == 15712  # every N-linked structure parses
    assert n_linked.failure_counts.get(PARSE_FAILED, 0) == 0
    assert n_linked.records_built == 12664
    assert n_linked.records_n_linked == 12664
    assert n_linked.records_with_accession == 8809
    assert n_linked.failure_counts == {UNSUPPORTED_RESIDUES: 3048}
    assert n_linked.records_built + n_linked.rows_failed == n_linked.rows_selected


def test_the_shipped_records_keep_what_the_strings_left_open(n_linked):
    fully_resolved = [
        s for s in n_linked.structures if not s.has_unresolved_linkage and not s.has_unresolved_anomericity
    ]
    assert len(fully_resolved) == 4001
    assert sum(1 for s in n_linked.structures if s.has_unresolved_linkage) == 8659
    assert sum(1 for s in n_linked.structures if s.has_unresolved_anomericity) == 1402


def test_every_shipped_record_is_tagged_and_traceable(n_linked):
    for structure in n_linked.structures[:200]:
        assert structure.reuse_status is ReuseStatus.OPEN_ATTRIBUTION
        assert structure.iupac_condensed
        assert "glycowork" in structure.source
