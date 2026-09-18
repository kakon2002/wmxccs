"""What the seeded corpus actually holds, asserted by loading it.

This file is the evidence behind the M0 report. Every number in it is computed
from the files on disk - the delivered transcriptions, the converted seed files
and what the loader makes of them - and compared against another number computed
the same way. Where a literal appears (28, 89, 117) it is the SECOND half of an
assertion whose first half came from a file, so that a conversion which lost the
same row from both sides is still caught.

Three things here are load-bearing and easy to lose quietly:

- the 2016 transcription is delivered under a name with a SPACE in it. Tidying
  the name breaks the provenance trail, and a loader pointed at the tidy name
  finds nothing and reports a clean run over zero rows.
- the four held 2016 rows are held because the loader DERIVED a suspected shared
  peak from the values themselves. Nothing in the file marks them. If the
  detector ever stops firing, those four rows clear in silence and a model learns
  that two different oligosaccharides share a cross section.
- the 89 held 2015 rows have TWO blockers, not one. CONTEXT.md records only the
  unstated drift gas. Resolving the gas alone would unblock nothing, and the day
  somebody resolves it this suite has to say so rather than report a surprise.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from tools.seed_struwe import FILE_2015, FILE_2016, _registry_status, _resolution_flags
from wmxccs.identity import DriftGas
from wmxccs.licensing import TrainingGateError, assert_trainable
from wmxccs.loader import (
    GATE_ANALYTE,
    GATE_LONE_CONFORMER,
    GATE_SHARED_PEAK,
    load_measurements_file,
)
from wmxccs.models import IMSType
from wmxccs.reuse import ReuseStatus
from wmxccs.sources import (
    SEED_STRUWE_2015,
    SEED_STRUWE_2016,
    STRUWE_2015,
    STRUWE_2016,
    dataset_for,
)

REPO = Path(__file__).resolve().parents[1]
SEED_DIR = REPO / "data" / "seed"
AS_DELIVERED = SEED_DIR / "as_delivered"

SEED_FILE_2016 = SEED_DIR / "struwe2016_chemcommun.csv"
SEED_FILE_2015 = SEED_DIR / "struwe2015_analyst.csv"

# The two delivered transcriptions, named exactly as they were delivered. The
# 2016 one really does carry a space before the extension.
DELIVERED_2016 = AS_DELIVERED / "struwe2016_ccs .csv"
DELIVERED_2015 = AS_DELIVERED / "struwe2015_analyst_ccs.csv"

# Counted from the delivered files in the row-conservation tests before they are
# used as anything. They appear as the second half of a comparison, never as the
# only statement of how many rows there are.
EXPECTED_ROWS_2016 = 28
EXPECTED_ROWS_2015 = 89
EXPECTED_ROWS_TOGETHER = 117

# The two CCS values the 2016 ESI reports against both LNH and LNnH.
SHARED_PEAK_VALUES = {(228.9, "[M+H]+"), (245.0, "[M+Cl]-")}


def data_rows(path: Path) -> list[dict[str, str]]:
    """Every data row of a CSV, counted with the csv module rather than by eye."""
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def held_records(report):
    """The records a report holds back: built, kept, and not cleared to train."""
    cleared = {id(record) for record in report.cleared}
    return [record for record in report.records if id(record) not in cleared]


@pytest.fixture(scope="module")
def report_2016():
    return load_measurements_file(SEED_FILE_2016)


@pytest.fixture(scope="module")
def report_2015():
    return load_measurements_file(SEED_FILE_2015)


@pytest.fixture(scope="module")
def corpus(report_2016, report_2015):
    """Both seed files loaded into one list, cleared or not.

    The corpus is what the platform HOLDS. Pairing is asked of everything held,
    not only of what may train, or a file that clears nothing would look like a
    file with nothing in it.
    """
    return list(report_2016.records) + list(report_2015.records)


# --- the provenance trail ---------------------------------------------------------------


def test_the_2016_transcription_is_kept_under_the_name_it_was_delivered_with_space_and_all():
    # A tidy-up that renames this file breaks the trail back to what the owner
    # delivered, and the converter, which looks the file up by that exact name,
    # would then find nothing and have nothing to convert.
    assert DELIVERED_2016.is_file(), f"the delivered 2016 transcription is missing from {AS_DELIVERED}"
    assert DELIVERED_2016.name == "struwe2016_ccs .csv"
    assert " .csv" in DELIVERED_2016.name, "the space before the extension is part of the delivered name"
    tidied = AS_DELIVERED / "struwe2016_ccs.csv"
    assert not tidied.exists(), "the transcription has been renamed; the delivered name is the provenance"


def test_the_2015_transcription_is_kept_under_the_name_it_was_delivered_with():
    assert DELIVERED_2015.is_file(), f"the delivered 2015 transcription is missing from {AS_DELIVERED}"
    assert DELIVERED_2015.name == "struwe2015_analyst_ccs.csv"


def test_the_converter_looks_for_the_transcriptions_under_the_names_they_are_kept_under():
    # The converter names its inputs as constants. If those constants and the
    # files in as_delivered ever disagree, the conversion cannot be re-run and
    # the seed files stop being reproducible from the transcriptions.
    assert FILE_2016 == DELIVERED_2016.name
    assert FILE_2015 == DELIVERED_2015.name


# --- row conservation, file against file -------------------------------------------------


@pytest.mark.parametrize(
    "delivered, seeded, expected",
    [
        pytest.param(DELIVERED_2016, SEED_FILE_2016, EXPECTED_ROWS_2016, id="2016"),
        pytest.param(DELIVERED_2015, SEED_FILE_2015, EXPECTED_ROWS_2015, id="2015"),
    ],
)
def test_every_row_of_the_delivered_transcription_survives_into_the_seed_file(delivered, seeded, expected):
    # Both counts are read off the files. The literal is asserted separately, so
    # that a conversion which lost a row from BOTH sides is still caught.
    delivered_rows = data_rows(delivered)
    seeded_rows = data_rows(seeded)
    assert len(seeded_rows) == len(delivered_rows), (
        f"{seeded.name} holds {len(seeded_rows)} rows against {len(delivered_rows)} in {delivered.name}:"
        " the conversion lost or invented a row"
    )
    assert len(delivered_rows) == expected


@pytest.mark.parametrize(
    "seeded, delivered",
    [
        pytest.param(SEED_FILE_2016, DELIVERED_2016, id="2016"),
        pytest.param(SEED_FILE_2015, DELIVERED_2015, id="2015"),
    ],
)
def test_the_loader_reads_exactly_the_rows_the_seed_file_holds(seeded, delivered):
    report = load_measurements_file(seeded)
    assert report.rows_read == len(data_rows(seeded))
    assert report.rows_read == len(data_rows(delivered))
    # The loader's own conservation law: nothing is dropped between reading and building.
    assert report.records_built + report.rows_failed == report.rows_read
    assert report.records_cleared + report.records_refused == report.records_built


def test_the_two_seed_files_together_hold_one_hundred_and_seventeen_transcribed_values(report_2016, report_2015):
    delivered = len(data_rows(DELIVERED_2016)) + len(data_rows(DELIVERED_2015))
    loaded = report_2016.rows_read + report_2015.rows_read
    built = report_2016.records_built + report_2015.records_built
    assert loaded == delivered
    assert built == delivered
    assert delivered == EXPECTED_ROWS_TOGETHER


# --- the 2016 file: 28 rows, 24 cleared, 4 held for a derived shared peak -------------------


def test_the_2016_seed_builds_a_record_from_every_row_and_fails_none(report_2016):
    assert report_2016.rows_read == EXPECTED_ROWS_2016
    assert report_2016.records_built == EXPECTED_ROWS_2016
    assert report_2016.rows_failed == 0
    assert report_2016.failure_counts == {}


def test_the_2016_seed_clears_twenty_four_records_and_holds_four(report_2016):
    assert report_2016.records_cleared == 24
    assert report_2016.records_refused == 4
    assert len(held_records(report_2016)) == 4


def test_the_four_held_2016_rows_are_held_for_the_shared_peak_reason_and_no_other(report_2016):
    assert set(report_2016.gate_counts) == {GATE_SHARED_PEAK}
    assert report_2016.gate_counts[GATE_SHARED_PEAK] == 4
    for example in report_2016.gate_examples[GATE_SHARED_PEAK]:
        assert "is reported for 2 analytes" in example


def test_the_held_2016_rows_are_the_LNH_and_LNnH_pairs_at_228_9_and_245_0(report_2016):
    """The finding CONTEXT.md records, re-derived from the file.

    This hold is DERIVED by the loader's shared-peak check from the values
    themselves; nothing in the file marks these rows. If the detector stops
    firing, these four rows clear in silence and a fit is taught that two
    different oligosaccharides share a cross section. That is what this test is
    here to notice, so it asserts which rows are held, not merely how many.
    """
    held = held_records(report_2016)
    assert {(record.ccs, record.adduct) for record in held} == SHARED_PEAK_VALUES

    for value, adduct in sorted(SHARED_PEAK_VALUES):
        at_value = [record for record in held if (record.ccs, record.adduct) == (value, adduct)]
        assert len(at_value) == 2, f"CCS {value} {adduct} should be held against both compounds"
        # Two rows, two DIFFERENT analytes. One analyte reported twice is a
        # replicate and must not be held; the point of the check is a collision.
        assert len({record.analyte.identity_key() for record in at_value}) == 2
        assert {record.analyte.display_name for record in at_value} == {"LNH", "LNnH"}
        # And they really are in one calibration group, which is what makes the
        # identical value a collision rather than two values from two methods.
        assert len({record.calibration_group for record in at_value}) == 1

    # The negative direction: the file's other 24 rows clear. A detector that
    # held everything would satisfy the assertions above and be useless.
    cleared_values = {(record.ccs, record.adduct) for record in report_2016.cleared}
    assert cleared_values & SHARED_PEAK_VALUES == set()
    assert (225.5, "[M+H]+") in cleared_values


def test_no_row_of_either_seed_file_carries_a_curation_flag_so_every_hold_is_derived():
    # If a seed row carried a curation flag the loader would hold it for that
    # reason instead, and the shared-peak test above would pass without the
    # detector doing anything at all.
    for path in (SEED_FILE_2016, SEED_FILE_2015):
        flags = {row.get("curation_flag", "") for row in data_rows(path)}
        assert flags <= {""}, f"{path.name} marks rows for curation; the holds would no longer be derived"


def test_the_2016_values_carry_a_stated_uncertainty_type_so_nothing_is_held_for_the_spread(report_2016):
    # The 2016 transcription's ccs_2sd column is loaded as two standard
    # deviations on the owner's instruction. A spread stored without its type,
    # or with 'unknown', would block all 28 rows rather than 4.
    assert {str(record.uncertainty_type) for record in report_2016.records} == {"two_sd"}
    assert {record.ccs_uncertainty is not None for record in report_2016.records} == {True}
    assert set(report_2016.uncertainty_types) == {"two_sd"}
    assert sum(report_2016.uncertainty_types.values()) == report_2016.records_cleared


# --- the 2015 file: 89 rows, nothing clears, and TWO blockers ------------------------------


def test_the_2015_seed_builds_a_record_from_every_row_and_fails_none(report_2015):
    assert report_2015.rows_read == EXPECTED_ROWS_2015
    assert report_2015.records_built == EXPECTED_ROWS_2015
    assert report_2015.rows_failed == 0
    assert report_2015.failure_counts == {}


def test_the_2015_seed_clears_nothing_and_holds_every_one_of_its_eighty_nine_records(report_2015):
    assert report_2015.records_cleared == 0
    assert report_2015.records_refused == EXPECTED_ROWS_2015
    assert len(held_records(report_2015)) == EXPECTED_ROWS_2015
    # Held, not dropped. Storing a record is not the same as training on it.
    assert len(report_2015.records) == EXPECTED_ROWS_2015


def test_every_held_2015_row_is_held_for_the_analyte_gate_reason_and_no_other(report_2015):
    assert set(report_2015.gate_counts) == {GATE_ANALYTE}
    assert report_2015.gate_counts[GATE_ANALYTE] == EXPECTED_ROWS_2015


def test_the_2015_hold_names_both_the_unstated_drift_gas_and_the_unknown_uncertainty_type(report_2015):
    """CONTEXT.md says the gas blocks these rows. It is not the only thing that does.

    Every 2015 record carries two blockers: the drift gas the ESI never states,
    and an uncertainty column with no stated type. A reader who fixes only the
    one CONTEXT.md names would expect 89 rows to clear and get none.
    """
    for record in report_2015.records:
        blockers = record.training_blockers()
        assert any("drift_gas is 'UNSTATED'" in blocker for blocker in blockers), blockers
        assert any("uncertainty_type is 'unknown'" in blocker for blocker in blockers), blockers
        assert len(blockers) >= 2

    # And the gate itself says both, in the message a person would read.
    with pytest.raises(TrainingGateError) as refusal:
        assert_trainable(report_2015.records[0])
    assert "drift_gas is 'UNSTATED'" in str(refusal.value)
    assert "uncertainty_type is 'unknown'" in str(refusal.value)


def test_the_loaders_own_report_of_the_2015_hold_names_both_blockers_and_not_just_the_gas(report_2015):
    # The example text is what a person reads in the M0 report, and it is the
    # thing CONTEXT.md got wrong. The loader truncates a gate example at 400
    # characters and the second blocker currently ends ON that boundary, so a few
    # more words in the first blocker's wording would leave the report naming the
    # gas alone. If that happens this fires, which is the right outcome: the
    # report has to say both.
    for example in report_2015.gate_examples[GATE_ANALYTE]:
        assert "drift_gas is 'UNSTATED'" in example
        assert "uncertainty_type is 'unknown'" in example, (
            "the gate example no longer names the second blocker; it is being cut off by the 400-character"
            f" truncation in loader._first_line. Example was: {example}"
        )


@pytest.mark.xfail(
    reason="SUSPECTED DEFECT in loader.MeasurementLoadReport: uncertainty_types is computed over `cleared`"
    " only, and there is no uncertainty_types_held beside adducts_held, gas_split_held, platforms_held and"
    " analyte_kinds_held. The 2015 seed clears nothing, so its summary says nothing at all about the spreads"
    " it holds - and an uncertainty type of 'unknown' is one of the two things blocking all 89 rows. The"
    " report's own docstring says the _held four exist 'so that a report of a fully held file still says what"
    " is in it'; this is the fifth, and it is missing.",
    strict=True,
)
def test_the_summary_of_a_fully_held_file_says_what_its_values_uncertainty_types_are(report_2015):
    # Every one of the 89 records reports an uncertainty type, and it is one of
    # the two reasons none of them may train.
    assert {str(record.uncertainty_type) for record in report_2015.records} == {"unknown"}
    written = [line.strip() for line in report_2015.summary().splitlines()]
    assert any(line.startswith("uncertainty") for line in written), (
        "the summary of a file that clears nothing never says what its values' spreads are, so the reader"
        " cannot see the second blocker at all"
    )


def test_resolving_the_drift_gas_alone_would_not_clear_a_single_2015_row(report_2015):
    # The second blocker, pinned from the direction somebody will actually come
    # at it: the day the calibration reference is read and the gas is settled,
    # this says there is still an uncertainty type to resolve.
    for record in report_2015.records:
        fields = record.model_dump()
        fields["drift_gas"] = DriftGas.HE
        with_gas_resolved = type(record)(**fields)
        assert with_gas_resolved.training_blockers(), "the gas was the only blocker after all"
        with pytest.raises(TrainingGateError) as refusal:
            assert_trainable(with_gas_resolved)
        assert "uncertainty_type is 'unknown'" in str(refusal.value)
        assert "drift_gas" not in str(refusal.value)


def test_no_seed_row_is_held_for_an_incomplete_conformer_set(report_2015, report_2016):
    # The 2015 file declares eleven two-conformer sets. All twenty-two rows are
    # present, so none is held for a missing sibling; if a transcription ever
    # loses one, the reason a row is held changes and this fires.
    declaring = [record for record in report_2015.records if (record.conformers_total or 1) > 1]
    assert len(declaring) == 22, "the 2015 file should declare eleven complete conformer pairs"
    assert GATE_LONE_CONFORMER not in report_2015.gate_counts
    assert GATE_LONE_CONFORMER not in report_2016.gate_counts


def test_the_2015_locator_names_the_sample_origin_so_two_origins_do_not_pool(report_2015):
    # Sample origin is folded into the locator, and the locator is what keeps one
    # origin's conformer pair from satisfying another origin's missing sibling.
    # Give origin a column of its own and the separation disappears.
    locators = {record.source_locator for record in report_2015.records}
    assert len(locators) > 1, "one locator for the whole file would pool every sample origin"
    assert all(locator.startswith("ESI Table S1 (") and locator.endswith(")") for locator in locators), locators


# --- zero matched pairs, and why it is over-determined -------------------------------------


def test_no_matched_ion_key_in_the_seeded_corpus_is_held_on_more_than_one_platform(corpus, report_2016, report_2015):
    """The expected M0 result: zero cross-platform matched pairs.

    Asserted over the corpus as a whole rather than per file, because a pair
    between the two files is exactly what somebody would hope for.
    """
    platforms_by_key: dict[object, set[str]] = {}
    for record in corpus:
        platforms_by_key.setdefault(record.matched_ion_key, set()).add(str(record.ims_type))
    paired = {key for key, platforms in platforms_by_key.items() if len(platforms) > 1}
    assert paired == set(), f"{len(paired)} matched ions are claimed on two platforms"

    # The keys are disjoint between the files too, so the corpus holds exactly as
    # many distinct ions as the two files hold separately: nothing was matched.
    assert len(platforms_by_key) == report_2016.matched_ion_keys_held + report_2015.matched_ion_keys_held
    # Not a vacuous zero: there are ions there to pair, they simply never pair.
    assert len(platforms_by_key) > 0


def test_both_seed_files_are_travelling_wave_so_there_is_no_second_platform_to_pair_with(corpus):
    # Reason one, on its own. Even if every other difference were resolved, one
    # platform cannot produce a cross-platform pair.
    assert {record.ims_type for record in corpus} == {IMSType.TWIMS}
    assert {str(record.ims_type) for record in corpus} == {"TWIMS"}


def test_the_two_seed_files_share_no_analyte_identity_at_all(report_2016, report_2015):
    # Reason two, on its own. The 2016 file identifies four milk
    # oligosaccharides by structure; the 2015 file identifies seven high-mannose
    # N-glycans by composition. No identifier is shared, so no ion could match
    # even on one platform.
    atoms_2016 = {atom for record in report_2016.records for atom in record.analyte.identity_atoms()}
    atoms_2015 = {atom for record in report_2015.records for atom in record.analyte.identity_atoms()}
    assert atoms_2016, "the 2016 analytes state no identifier at all"
    assert atoms_2015, "the 2015 analytes state no identifier at all"
    assert atoms_2016 & atoms_2015 == set()

    keys_2016 = {record.analyte.identity_key() for record in report_2016.records}
    keys_2015 = {record.analyte.identity_key() for record in report_2015.records}
    assert len(keys_2016) == 4
    assert len(keys_2015) == 7
    assert keys_2016 & keys_2015 == set()


def test_the_gas_alone_would_keep_the_two_seed_files_apart(report_2016, report_2015):
    # Reason three, on its own. The gas the value refers to is part of the
    # matched-ion key, so helium values and values of an unstated gas cannot land
    # on one key however alike everything else is.
    gases_2016 = {record.drift_gas for record in report_2016.records}
    gases_2015 = {record.drift_gas for record in report_2015.records}
    assert gases_2016 == {DriftGas.HE}
    assert gases_2015 == {DriftGas.UNSTATED}
    assert gases_2016 & gases_2015 == set()
    # And the gas really is in the key, so the disjointness above bites.
    for record in list(report_2016.records) + list(report_2015.records):
        assert record.matched_ion_key.drift_gas == record.drift_gas


# --- what backs the seeded analytes --------------------------------------------------------


@pytest.mark.parametrize(
    "seed_file, dataset",
    [
        pytest.param(SEED_FILE_2016, SEED_STRUWE_2016, id="2016"),
        pytest.param(SEED_FILE_2015, SEED_STRUWE_2015, id="2015"),
    ],
)
def test_every_seeded_analyte_names_the_dataset_provenance_the_registry_records(seed_file, dataset):
    report = load_measurements_file(seed_file)
    sources = {record.analyte.source for record in report.records}
    assert sources == {dataset.provenance}
    # Read from the file as text as well: the gate matches this string character
    # for character, so a stray edit to the column is a broken licence claim.
    assert {row["analyte_source"] for row in data_rows(seed_file)} == {dataset.provenance}
    assert dataset_for(dataset.provenance) is dataset


def test_a_provenance_string_that_is_merely_similar_backs_nothing():
    # The negative direction of the check above. No fuzzy matching: a dataset
    # claim is matched exactly or it is not matched.
    assert dataset_for("Struwe 2016 Chem Commun ESI Table S1") is None
    assert dataset_for(SEED_STRUWE_2016.provenance.replace("wmxccs", "wmxglycan")) is None
    assert dataset_for(None) is None
    # Surrounding whitespace is not a different dataset, though.
    assert dataset_for(f"  {SEED_STRUWE_2016.provenance}  ") is SEED_STRUWE_2016


@pytest.mark.parametrize(
    "seed_file, delivered, entry",
    [
        pytest.param(SEED_FILE_2016, DELIVERED_2016, STRUWE_2016, id="2016"),
        pytest.param(SEED_FILE_2015, DELIVERED_2015, STRUWE_2015, id="2015"),
    ],
)
def test_the_seeded_reuse_status_came_from_the_registry_and_not_from_the_licence_cell(seed_file, delivered, entry):
    # The transcription carries its own licence cell. A licence claim in a data
    # cell is not a licence, so what is seeded is the registry's status, read by
    # a named person on a stated date.
    cells = {row["licence"] for row in data_rows(delivered)}
    seeded = {row["reuse_status"] for row in data_rows(seed_file)}
    assert seeded == {entry.reuse_status.value}
    assert seeded == {ReuseStatus.OPEN_ATTRIBUTION.value}
    assert seeded & cells == set(), "the licence cell was copied through instead of being looked up"
    assert entry.reported_by and entry.reported_on, "a licence status needs who read the terms and when"


def test_the_licence_cell_in_a_row_decides_nothing_and_an_unknown_doi_stops_the_conversion():
    # Both directions of the rule the converter states: the registry's status is
    # what is written whatever the cell says, and a DOI with no record refuses.
    assert _registry_status(STRUWE_2016.doi, "public domain, honestly", "a test") is STRUWE_2016.reuse_status
    assert _registry_status(STRUWE_2016.doi, "", "a test") is STRUWE_2016.reuse_status
    with pytest.raises(SystemExit) as refusal:
        _registry_status("10.0000/not-a-real-doi", "CC-BY-4.0", "a test")
    assert "no licence record" in str(refusal.value)


# --- how each file's analytes are identified ------------------------------------------------


def test_the_2016_analytes_are_identified_by_structure_and_carry_no_composition(report_2016):
    # No composition is derived from an IUPAC string: that needs a table mapping
    # monosaccharide names onto residue classes, which is glycan chemistry and
    # does not come across. The structure is the finer identifier anyway.
    for record in report_2016.records:
        assert record.analyte.composition is None
        assert record.analyte.iupac_condensed
        assert record.analyte.identity_key()[1] == "iupac"


def test_the_2015_analytes_are_identified_by_composition_and_carry_no_structure(report_2015):
    # The structures sit in the paper's Figure 1, which nobody here has read.
    # Nothing is filled in from memory, so these key on composition.
    for record in report_2015.records:
        assert record.analyte.composition is not None
        assert record.analyte.iupac_condensed is None
        assert record.analyte.wurcs is None
        assert record.analyte.glytoucan_ac is None
        assert record.analyte.identity_key()[1] == "composition"


def test_the_2015_analytes_do_not_claim_a_resolution_their_records_cannot_support(report_2015):
    # A record holding only a composition cannot claim resolved linkage or
    # anomericity. Claiming it would make two isomers of one composition look
    # like one analyte with two cross sections.
    for record in report_2015.records:
        assert record.analyte.has_unresolved_linkage is True
        assert record.analyte.has_unresolved_anomericity is True


def test_the_2016_resolution_flags_are_read_off_the_structure_strings_in_the_file():
    # The flags in the seed file are derived from the structure the row carries,
    # not asserted by the transcriber. Re-derive them and compare.
    for row in data_rows(SEED_FILE_2016):
        linkage, anomericity = _resolution_flags(row["analyte_iupac_condensed"])
        assert row["analyte_has_unresolved_linkage"] == ("true" if linkage else "false")
        assert row["analyte_has_unresolved_anomericity"] == ("true" if anomericity else "false")
    # All four 2016 structures state their linkages and anomeric configurations.
    assert {row["analyte_has_unresolved_linkage"] for row in data_rows(SEED_FILE_2016)} == {"false"}


@pytest.mark.parametrize(
    "iupac, expected",
    [
        pytest.param("Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)Glc", (False, False), id="fully-resolved"),
        pytest.param("Gal(?1-4)Glc", (False, True), id="anomericity-unknown"),
        pytest.param("Gal(b1-?)Glc", (True, False), id="linkage-position-unknown"),
        pytest.param("Gal(?1-?)Glc", (True, True), id="neither-stated"),
        pytest.param("Hex5HexNAc2", (True, True), id="no-linkage-at-all"),
        pytest.param("", (True, True), id="nothing-stated"),
    ],
)
def test_a_structure_string_stating_nothing_is_not_evidence_that_its_linkages_are_known(iupac, expected):
    # The direction that matters: absence of a linkage in the string must read as
    # unresolved, never as resolved. A blank or composition-shaped string that
    # came back "resolved" would let a composition-only record claim a structure.
    assert _resolution_flags(iupac) == expected


# --- the calibration reference lineage the objective document has no field for ---------------


@pytest.mark.parametrize(
    "seed_file, reference_set",
    [
        pytest.param(SEED_FILE_2016, "helium CCS reference values", id="2016"),
        pytest.param(SEED_FILE_2015, "DTCCS dextran ladder, gas not stated in ESI", id="2015"),
    ],
)
def test_every_seeded_value_is_derived_and_records_the_reference_set_it_was_calibrated_against(
    seed_file, reference_set
):
    # Both files are travelling-wave, so every value is calibration-derived
    # rather than primary, and each carries the reference set its source names.
    report = load_measurements_file(seed_file)
    for record in report.records:
        assert record.ccs_is_calibrated is True
        assert str(record.calibration_lineage) == "derived"
        assert record.calibration_reference is not None
        assert record.calibration_reference.reference_set == reference_set


def test_no_seeded_value_can_yet_be_traced_back_to_a_primary_measurement(report_2016, report_2015):
    # The honest answer is that nobody here knows: resolving it means reading the
    # paper the calibration cites. None, not False, and never True by assumption
    # - this is the reference-circularity risk the objective document names.
    for record in list(report_2016.records) + list(report_2015.records):
        assert record.traces_to_primary is None
