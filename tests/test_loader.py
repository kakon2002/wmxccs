"""The loader loses nothing and invents nothing.

Every row of a file ends up in exactly one bucket, and the two conservation laws
in MeasurementLoadReport say so arithmetically:

    records_built + rows_failed == rows_read
    records_cleared + records_refused == records_built

A row that vanishes between those counts is the failure this module exists to
prevent, and it is the failure that looks most like success: a clean report over
fewer rows than the file holds.

CSV text is built from wmxccs.loader.COLUMNS rather than written out by hand, so
a new column added to the format does not rewrite this file. The record builders
in conftest cannot serve here: they build records, and the whole point of the
loader is the path from text to record, which has to start as text.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path

import pytest

from conftest import INCHIKEY, LNH_IUPAC, LNNH_IUPAC
from wmxccs.identity import AnalyteKind
from wmxccs.loader import (
    ANALYTE_REJECTED,
    COLUMNS,
    GATE_ANALYTE,
    GATE_CURATION,
    GATE_LICENCE,
    GATE_LONE_CONFORMER,
    GATE_SHARED_PEAK,
    GATE_UNBACKED_CLAIM,
    MEASUREMENT_REJECTED,
    NOT_A_BOOLEAN,
    NOT_A_DATE,
    NOT_A_NUMBER,
    SYNTHETIC_STATUS_IN_FILE,
    UNKNOWN_ANALYTE_KIND,
    WRONG_CELL_COUNT,
    load_measurements,
    load_measurements_file,
    read_table,
)
from wmxccs.models import UncertaintyType
from wmxccs.reuse import ReuseStatus
from wmxccs.sources import SEED_STRUWE_2016, STRUWE_2016

LABEL = "a test file"
SEED = Path(__file__).resolve().parent.parent / "data" / "seed"

# The licence half of a row that may train: a value whose DOI is on record, and
# an analyte identity whose source names a registered dataset. Taken from the
# registry rather than typed, so a change to either entry is felt here.
BACKED_DOI = STRUWE_2016.doi
BACKED_PROVENANCE = SEED_STRUWE_2016.provenance

# One row that builds a record AND clears the gate. Everything else in this file
# is this row with something changed, so a test says what it is about by what it
# overrides.
CLEAN: dict[str, str] = {
    "analyte_kind": "glycan",
    "analyte_display_name": "LNH",
    "analyte_source": BACKED_PROVENANCE,
    "analyte_reuse_status": "open_attribution",
    "analyte_iupac_condensed": LNH_IUPAC,
    "analyte_has_unresolved_linkage": "false",
    "analyte_has_unresolved_anomericity": "false",
    "analyte_reducing_end_label": "native",
    "analyte_derivatisation": "underivatised",
    "ccs": "225.5",
    "ccs_uncertainty": "0.4",
    "uncertainty_type": "two_sd",
    "adduct": "[M+H]+",
    "charge": "1",
    "polarity": "positive",
    "ims_type": "TWIMS",
    "drift_gas": "He",
    "calibrant": "dextran",
    "cell_gas": "N2",
    "instrument": "Synapt G2-S HDMS",
    "source": "struwe2016",
    "doi": BACKED_DOI,
    "source_locator": "ESI Table S1",
    "replicates": "2",
    "reuse_status": "open_attribution",
}

# The analyte columns that say WHICH analyte. Blanked wholesale when a row builds
# a kind other than the glycan CLEAN describes, because the models forbid unknown
# fields and a leftover glycan column would refuse a peptide for the wrong reason.
IDENTITY_COLUMNS = tuple(
    name
    for name in COLUMNS
    if name.startswith("analyte_")
    and name not in ("analyte_kind", "analyte_source", "analyte_reuse_status", "analyte_display_name")
)


def _line(values) -> str:
    """One physical CSV line, quoted by the csv module so a cell may hold a comma."""
    buffer = io.StringIO()
    csv.writer(buffer, lineterminator="\n").writerow(list(values))
    return buffer.getvalue()


def header(columns: tuple[str, ...] = COLUMNS) -> str:
    return _line(columns)


def row(columns: tuple[str, ...] = COLUMNS, **overrides: object) -> str:
    """One data row as CSV text: CLEAN with `overrides` applied, in `columns` order."""
    unknown = sorted(set(overrides) - set(COLUMNS) - set(columns))
    if unknown:  # a typo in a test must not quietly become a row with a default in it
        raise AssertionError(f"{unknown} are not columns of the format")
    cells = dict.fromkeys(COLUMNS, "") | CLEAN | {name: str(value) for name, value in overrides.items()}
    return _line(cells.get(name, "") for name in columns)


def analyte_row(kind: str, **overrides: object) -> str:
    """A row building `kind`, with every identity column blanked but the ones given."""
    blanked: dict[str, object] = {name: "" for name in IDENTITY_COLUMNS}
    blanked.update(overrides)
    return row(analyte_kind=kind, **blanked)


def table(*rows: str, columns: tuple[str, ...] = COLUMNS) -> str:
    return header(columns) + "".join(rows)


def load(text: str):
    return load_measurements(text, label=LABEL)


# --- the file as a whole: faults that raise rather than count --------------------------


def test_a_file_with_no_header_is_refused_whole_rather_than_read_as_an_empty_run():
    with pytest.raises(ValueError, match="no header row"):
        read_table("", label=LABEL)


def test_a_header_with_an_unnamed_column_is_refused_whole():
    # An unnamed column cannot be matched to a field, so every cell under it is
    # silently discarded. Refused, because the discard would not show anywhere.
    with pytest.raises(ValueError, match="unnamed column"):
        read_table("ccs,,adduct\n1,2,3\n", label=LABEL)


def test_a_header_naming_one_column_twice_is_refused_rather_than_losing_that_columns_cells():
    # csv.DictReader keeps only the last cell under a repeated name, so one of
    # the two columns would vanish without a count.
    with pytest.raises(ValueError, match="more than once"):
        read_table("ccs,adduct,ccs\n1,2,3\n", label=LABEL)


def test_an_unclosed_quotation_mark_is_refused_whole_rather_than_folding_rows_together():
    with pytest.raises(ValueError, match="malformed CSV"):
        read_table('ccs,adduct\n"225.5,[M+H]+\n', label=LABEL)


def test_a_data_row_spanning_two_physical_lines_is_refused_whole():
    # Two stray quotation marks fencing a whole row is the one quoting fault the
    # csv module reads happily, one row short, with nothing to say it happened.
    with pytest.raises(ValueError, match="runs on to line"):
        read_table('ccs,adduct\n"225.5\n[M+H]+",x\n', label=LABEL)


def test_a_header_spanning_two_physical_lines_is_refused_whole():
    with pytest.raises(ValueError, match="runs on to line"):
        read_table('"ccs\nvalue",adduct\n1,2\n', label=LABEL)


def test_a_missing_required_column_raises_rather_than_reporting_a_clean_run_over_nothing():
    # The dangerous shape: without the raise, every row fails for a reason that
    # reads like bad data rather than like the wrong file.
    columns = tuple(name for name in COLUMNS if name != "ccs")
    with pytest.raises(ValueError, match="missing the required columns"):
        load(table(row(columns=columns), columns=columns))


def test_a_header_with_no_data_rows_is_a_valid_empty_load_and_not_an_error():
    report = load(table())
    assert report.rows_read == 0
    assert report.records_built == 0
    assert report.records_cleared == 0
    assert report.readiness is not None  # the thresholds are visible from the first row, and before it


def test_a_blank_line_is_skipped_and_drops_nothing():
    report = load(header() + row(ccs="225.5") + "\n" + row(ccs="233.6"))
    assert report.rows_read == 2
    assert report.records_built == 2
    assert report.rows_failed == 0


def test_a_utf8_byte_order_mark_is_stripped_rather_than_hiding_the_first_column():
    # With the mark left on, "analyte_kind" reads as "﻿analyte_kind": a
    # required column missing, and a whole file refused for an invisible reason.
    report = load("﻿" + table(row()))
    assert report.rows_read == 1
    assert report.records_built == 1
    assert report.unknown_columns == ()


def test_a_column_the_format_does_not_know_is_reported_and_its_cells_ignored():
    columns = COLUMNS + ("collision_energy",)
    report = load(table(row(columns=columns, collision_energy="30 V"), columns=columns))
    assert report.unknown_columns == ("collision_energy",)
    assert report.records_built == 1


# --- the conservation laws -------------------------------------------------------------


def test_the_conservation_laws_hold_over_a_file_holding_every_outcome_at_once():
    """The centrepiece: one file with a row of every kind, and nothing lost between them.

    A clean row, a row with the wrong cell count, a row with a bad number, a row
    the analyte validator refuses, a row the measurement validator refuses, a row
    a person held for curation, and a row the licence gate refuses. Seven rows,
    three records, one of them trainable, and the two laws hold across all of it.
    """
    text = table(
        row(ccs="225.5"),  # clean: built and cleared
        "glycan,LNH\n",  # two cells against a full header: read at all
        row(ccs="1_000"),  # a number in a spelling this format does not take
        row(ccs="231.1", analyte_iupac_condensed="", analyte_composition=""),  # a glycan identified by nothing
        row(ccs="232.2", polarity="negative"),  # charge +1 against negative polarity
        row(ccs="233.3", curation_flag="checking the ESI against the table"),  # built, held by a person
        row(ccs="234.4", reuse_status="unverified"),  # built, refused by the gate
    )
    report = load(text)

    assert report.rows_read == 7
    assert report.records_built + report.rows_failed == report.rows_read
    assert report.records_cleared + report.records_refused == report.records_built

    assert report.records_built == 3
    assert report.rows_failed == 4
    assert report.records_cleared == 1
    assert report.records_refused == 2

    assert report.failure_counts == {
        WRONG_CELL_COUNT: 1,
        NOT_A_NUMBER: 1,
        ANALYTE_REJECTED: 1,
        MEASUREMENT_REJECTED: 1,
    }
    assert report.gate_counts == {GATE_CURATION: 1, GATE_LICENCE: 1}


@pytest.mark.parametrize("name", ["struwe2016_chemcommun.csv", "struwe2015_analyst.csv"])
def test_the_conservation_laws_hold_over_each_seeded_transcription(name):
    """The laws again, over the real files, through the real file entry point.

    No count is written down here: the counts belong in the milestone report, and
    a test that pinned them would fail the day a row is corrected. What must hold
    whatever the file says is that every row is accounted for and that a shipped
    transcription still parses.
    """
    path = SEED / name
    if not path.exists():
        pytest.skip(f"{name} is not in data/seed")
    report = load_measurements_file(path)
    assert report.label == name
    assert report.rows_read > 0
    assert report.records_built + report.rows_failed == report.rows_read
    assert report.records_cleared + report.records_refused == report.records_built
    assert report.failure_counts == {}  # a seeded file that stopped parsing is a regression, not data


def test_a_row_with_the_wrong_cell_count_is_counted_against_its_own_line_and_not_skipped():
    report = load(table(row(), "glycan,LNH\n", row(ccs="233.6")))
    assert report.rows_read == 3
    assert report.failure_counts[WRONG_CELL_COUNT] == 1
    assert report.records_built == 2
    # The fault stays inside its own line: the example names it, and the rows
    # either side of it still load.
    assert "line 3" in report.failure_examples[WRONG_CELL_COUNT][0]


def test_a_record_the_gate_refuses_is_kept_in_the_file_and_only_kept_out_of_training():
    # Storing a record is not the same as training on it. A gate that dropped the
    # row would make the file's own contents unreportable.
    report = load(table(row(reuse_status="unverified")))
    assert report.records_built == 1
    assert report.records_cleared == 0
    assert report.records_refused == 1
    assert len(report.records) == 1
    assert report.records[0].reuse_status is ReuseStatus.UNVERIFIED
    # And the held record is still described, which is what the _held counts are for.
    assert report.matched_ion_keys_held == 1
    assert report.matched_ion_keys == 0


# --- statuses, placeholders and holds --------------------------------------------------


@pytest.mark.parametrize("column", ["reuse_status", "analyte_reuse_status"])
@pytest.mark.parametrize("spelling", ["synthetic_fixture", "SYNTHETIC_FIXTURE", "Synthetic_Fixture"])
def test_synthetic_fixture_in_either_status_column_is_refused_before_any_record_exists(column, spelling):
    # The status is a test's declaration that a record was built in code. A row
    # in a file is a real record, so it never carries one, however it is spelt.
    report = load(table(row(**{column: spelling})))
    assert report.records_built == 0
    assert report.failure_counts == {SYNTHETIC_STATUS_IN_FILE: 1}
    assert report.gate_counts == {}


def test_a_synthetic_fixture_row_is_counted_under_its_own_reason_not_under_the_licence_gate():
    report = load(table(row(reuse_status="synthetic_fixture")))
    assert SYNTHETIC_STATUS_IN_FILE in report.failure_counts
    assert "synthetic_fixture" in report.failure_examples[SYNTHETIC_STATUS_IN_FILE][0]


@pytest.mark.parametrize(
    ("column", "filler", "reason"),
    [
        ("analyte_display_name", "n/a", ANALYTE_REJECTED),
        ("analyte_display_name", "n.d.", ANALYTE_REJECTED),
        ("analyte_display_name", "unknown", ANALYTE_REJECTED),
        ("instrument", "n/a", MEASUREMENT_REJECTED),
        ("instrument", "n.d.", MEASUREMENT_REJECTED),
        ("calibrant", "-", MEASUREMENT_REJECTED),
        ("source_locator", "unknown", MEASUREMENT_REJECTED),
    ],
)
def test_a_placeholder_cell_reaches_validation_and_is_refused_there_not_mapped_to_null(column, filler, reason):
    # A loader that mapped fillers to null would absorb a data-quality fault
    # instead of showing it: the row would load, and nothing would say that the
    # instrument or the calibrant had been thrown away.
    report = load(table(row(**{column: filler})))
    assert report.records_built == 0
    assert report.failure_counts == {reason: 1}


def test_a_blank_cell_is_absent_so_the_models_own_default_applies():
    report = load(table(row(replicates="", instrument="", measured_on="")))
    record = report.records[0]
    assert record.replicates is None
    assert record.instrument is None
    assert record.measured_on is None


def test_uncertainty_type_unknown_is_a_real_value_that_is_kept_and_then_held():
    # 'unknown' is spelt like a placeholder but is a stated enum member here: it
    # records that the source reports no spread. Kept, and it blocks training on
    # its own, which is a different outcome from being dropped.
    report = load(table(row(ccs_uncertainty="", uncertainty_type="unknown")))
    assert report.records_built == 1
    assert report.records[0].uncertainty_type is UncertaintyType.UNKNOWN
    assert report.gate_counts == {GATE_ANALYTE: 1}


def test_a_curation_flag_holds_an_otherwise_clean_row_without_dropping_it():
    report = load(table(row(curation_flag="transcription not yet checked against the ESI")))
    assert report.records_built == 1
    assert report.records_cleared == 0
    assert report.gate_counts == {GATE_CURATION: 1}
    assert "transcription not yet checked" in report.gate_examples[GATE_CURATION][0]


def test_a_blank_curation_flag_holds_nothing():
    report = load(table(row(curation_flag="")))
    assert report.records_cleared == 1
    assert report.gate_counts == {}


def test_the_gate_counts_its_three_refusals_separately():
    # An unbacked claim, a refused licence and an undefined ion are three
    # different things to fix, and one bucket for all three would say only that
    # something was wrong.
    text = table(
        row(ccs="225.5", doi="10.9999/nobody-has-read-this"),
        row(ccs="226.6", reuse_status="unverified"),
        row(ccs="227.7", drift_gas="UNSTATED"),
    )
    report = load(text)
    assert report.records_built == 3
    assert report.records_cleared == 0
    assert report.gate_counts == {GATE_UNBACKED_CLAIM: 1, GATE_LICENCE: 1, GATE_ANALYTE: 1}


def test_a_row_whose_licence_record_backs_its_claim_clears_the_gate():
    report = load(table(row()))
    assert report.records_cleared == 1
    assert report.gate_counts == {}


# --- the shared-peak check -------------------------------------------------------------


def test_two_different_analytes_with_one_ccs_in_one_calibration_group_are_both_held():
    # LNH and LNnH share a composition and not a cross section. One identical
    # value against both is almost certainly one unresolved peak transcribed
    # twice, and training on it would teach that two molecules share a value.
    text = table(
        row(analyte_display_name="LNH", analyte_iupac_condensed=LNH_IUPAC, ccs="228.9"),
        row(analyte_display_name="LNnH", analyte_iupac_condensed=LNNH_IUPAC, ccs="228.9"),
    )
    report = load(text)
    assert report.records_built == 2
    assert report.records_cleared == 0
    assert report.gate_counts == {GATE_SHARED_PEAK: 2}


def test_the_same_analyte_measured_twice_at_one_ccs_is_a_replicate_and_is_not_held():
    # The negative direction of the check above: two rows of one analyte carrying
    # one value are a replicate, and holding them would cost data for nothing.
    text = table(
        row(ccs="228.9", source_locator="ESI Table S1"),
        row(ccs="228.9", source_locator="ESI Table S2"),
    )
    report = load(text)
    assert report.records_built == 2
    assert report.records_cleared == 2
    assert report.gate_counts == {}


def test_two_analytes_sharing_a_ccs_in_different_calibration_groups_are_not_held():
    # The check is per calibration group: two values produced different ways are
    # not one peak however alike the numbers are.
    text = table(
        row(analyte_iupac_condensed=LNH_IUPAC, ccs="228.9", adduct="[M+H]+", charge="1"),
        row(analyte_iupac_condensed=LNNH_IUPAC, ccs="228.9", adduct="[M+Na]+", charge="1"),
    )
    report = load(text)
    assert report.records_cleared == 2
    assert report.gate_counts == {}


def test_a_row_sharing_a_value_with_a_row_that_did_not_validate_is_held_fail_closed():
    # The failed row's calibration group could not be computed, so it could not
    # be compared. Held until that row is fixed, rather than cleared on the
    # strength of a comparison nobody could make.
    text = table(
        row(ccs="300.0"),
        row(
            ccs="300.0",
            analyte_iupac_condensed="",
            analyte_composition="Hex9HexNAc2",
            analyte_has_unresolved_linkage="true",
            analyte_has_unresolved_anomericity="true",
            polarity="negative",  # against charge +1: the measurement validator refuses it
        ),
    )
    report = load(text)
    assert report.failure_counts == {MEASUREMENT_REJECTED: 1}
    assert report.records_built == 1
    assert report.records_cleared == 0
    assert report.gate_counts == {GATE_SHARED_PEAK: 1}
    assert "did not validate" in report.gate_examples[GATE_SHARED_PEAK][0]


def test_two_rows_identified_by_nothing_are_read_as_two_analytes_and_not_merged_into_one():
    """Fail closed where identity is unclear. Pinned on the helper, because nothing else can reach it.

    No file can produce this case today: every analyte that validates states at
    least one identifier, and the fail-closed comparison only ever pairs ONE
    built row with ONE failed row, so two rows identified by nothing never meet
    inside a single call. The rule the helper encodes is still the one that
    decides the safe direction - two rows identified by nothing are DIFFERENT
    analytes, so an identical value holds them rather than clearing them as a
    replicate - and a shared "unidentified" bucket would merge every such row
    into one analyte and clear the lot. Pinned now rather than left to whichever
    caller reaches it first: M2 pairs matched ions from these same atoms.
    """
    from wmxccs.loader import _components

    assert len(_components([("line 2", frozenset()), ("line 3", frozenset())])) == 2
    # The other direction: rows sharing any atom are one analyte, and a shared
    # value between them is a replicate rather than a suspected shared peak.
    shared = frozenset({f"inchikey:{INCHIKEY}"})
    assert len(_components([("line 2", shared), ("line 3", shared)])) == 1


def test_a_row_sharing_a_value_with_a_failed_row_of_the_same_analyte_is_not_held():
    # Fail-closed does not mean hold everything: the failed row names the same
    # analyte, so the two are a replicate and not a suspected collision.
    text = table(
        row(ccs="300.0"),
        row(ccs="300.0", replicates="two"),  # a bad number, and the same analyte
    )
    report = load(text)
    assert report.failure_counts == {NOT_A_NUMBER: 1}
    assert report.records_cleared == 1
    assert report.gate_counts == {}


# --- the lone-conformer check ----------------------------------------------------------


def test_a_conformer_pair_whose_members_are_both_in_the_file_clears():
    text = table(
        row(ccs="250.0", conformer="1", conformers_total="2"),
        row(ccs="260.0", conformer="2", conformers_total="2"),
    )
    report = load(text)
    assert report.records_built == 2
    assert report.records_cleared == 2
    assert report.gate_counts == {}


def test_a_conformer_whose_sibling_is_not_in_the_file_is_held():
    # "Conformer 1 of 2" alone is either a transcription that lost a row or a
    # value that is not what it claims. Training on it would teach that this
    # molecule has a single peak at this value.
    report = load(table(row(ccs="250.0", conformer="1", conformers_total="2")))
    assert report.records_built == 1
    assert report.records_cleared == 0
    assert report.gate_counts == {GATE_LONE_CONFORMER: 1}
    assert "declares 2 conformers" in report.gate_examples[GATE_LONE_CONFORMER][0]


def test_two_samples_each_missing_a_sibling_do_not_pool_across_source_locators():
    # Sample A kept conformer 1, sample B kept conformer 2. Pooled on analyte
    # alone the indices 1 and 2 are both present and the set looks complete, so
    # two lost rows would clear each other. The locator is what keeps the two
    # samples apart.
    text = table(
        row(ccs="250.0", conformer="1", conformers_total="2", source_locator="ESI Table S1 (sample A)"),
        row(ccs="260.0", conformer="2", conformers_total="2", source_locator="ESI Table S1 (sample B)"),
    )
    report = load(text)
    assert report.records_built == 2
    assert report.records_cleared == 0
    assert report.gate_counts == {GATE_LONE_CONFORMER: 2}


def test_a_row_reporting_a_single_peak_claims_no_siblings_and_is_not_held():
    report = load(table(row(ccs="250.0", conformer="1", conformers_total="1")))
    assert report.records_cleared == 1
    assert report.gate_counts == {}


# --- numbers, booleans and dates, in one spelling each ---------------------------------


@pytest.mark.parametrize(
    "spelling",
    [
        "1_000",  # Python's float() takes digit-group underscores; a transcriber does not write them
        "٣٠٠",  # Arabic-Indic digits: float() takes every Unicode digit
        "nan",
        "inf",
        "-inf",
        "1,000",
        "225.5 A^2",
    ],
)
def test_a_ccs_written_in_a_second_spelling_is_refused_as_not_a_number(spelling):
    report = load(table(row(ccs=spelling)))
    assert report.records_built == 0
    assert report.failure_counts == {NOT_A_NUMBER: 1}


@pytest.mark.parametrize("spelling", ["2_0", "٣", "2.0", "two"])
def test_an_integer_column_written_in_a_second_spelling_is_refused_as_not_a_number(spelling):
    report = load(table(row(replicates=spelling)))
    assert report.records_built == 0
    assert report.failure_counts == {NOT_A_NUMBER: 1}


def test_a_plain_ascii_number_is_read_as_the_number_it_shows():
    report = load(table(row(ccs="225.5", ccs_uncertainty="0.4", replicates="2")))
    record = report.records[0]
    assert record.ccs == 225.5
    assert record.ccs_uncertainty == 0.4
    assert record.replicates == 2


@pytest.mark.parametrize(
    "spelling",
    [
        "1462060800",  # a Unix timestamp; date.fromisoformat reads it as 1462-06-08
        "01/02/2026",
        "2026-W38-5",  # an ISO week date: a date, but not this format's one spelling
        "18-09-2026",
        "2026-9-18",
        "September 2026",
    ],
)
def test_a_measurement_date_not_written_yyyy_mm_dd_is_refused_as_not_a_date(spelling):
    report = load(table(row(measured_on=spelling)))
    assert report.records_built == 0
    assert report.failure_counts == {NOT_A_DATE: 1}


def test_a_date_written_yyyy_mm_dd_is_read_as_the_date_it_shows():
    from datetime import date

    report = load(table(row(measured_on="2026-09-18")))
    assert report.records[0].measured_on == date(2026, 9, 18)


@pytest.mark.parametrize("spelling", ["true", "TRUE", "False", "yes", "no", "1", "0"])
def test_a_boolean_column_takes_the_spellings_the_format_lists(spelling):
    report = load(table(row(analyte_has_unresolved_linkage=spelling, analyte_has_unresolved_anomericity="true")))
    assert report.records_built == 1


@pytest.mark.parametrize("spelling", ["t", "y", "unresolved", "2"])
def test_a_boolean_column_written_any_other_way_is_refused(spelling):
    report = load(table(row(analyte_has_unresolved_linkage=spelling)))
    assert report.records_built == 0
    assert report.failure_counts == {NOT_A_BOOLEAN: 1}


def test_a_semicolon_separated_list_column_is_read_as_a_list_of_its_parts():
    report = load(
        table(analyte_row("peptide", analyte_sequence="PEPTIDE", analyte_modifications="Phospho@S5; Acetyl@K2"))
    )
    assert report.records_built == 1, report.failure_counts
    assert report.records[0].analyte.modifications == ("Acetyl@K2", "Phospho@S5")


# --- the analyte union: a row builds the kind it NAMES ---------------------------------

KIND_ROWS = {
    AnalyteKind.SMALL_MOLECULE: {"analyte_inchikey": INCHIKEY},
    AnalyteKind.PEPTIDE: {"analyte_sequence": "PEPTIDE"},
    AnalyteKind.GLYCAN: {"analyte_composition": "Hex5HexNAc2"},
    AnalyteKind.GLYCOPEPTIDE: {
        "analyte_sequence": "NLTK",
        "analyte_glycan_composition": "Hex5HexNAc2",
        "analyte_attachment_site": "N297",
    },
    AnalyteKind.PROTEIN: {"analyte_accession": "P01857", "analyte_folding_state": "native"},
    AnalyteKind.INTACT_ANTIBODY: {"analyte_inn": "trastuzumab", "analyte_folding_state": "native"},
    AnalyteKind.ADC: {
        "analyte_inn": "trastuzumab",
        "analyte_linker_payload_class": "vc-MMAE",
        "analyte_dar": "2",
        "analyte_folding_state": "native",
    },
}


@pytest.mark.parametrize("kind", list(AnalyteKind), ids=lambda kind: kind.value)
def test_a_row_builds_the_analyte_kind_it_names(kind):
    report = load(table(analyte_row(kind.value, **KIND_ROWS[kind])))
    assert report.records_built == 1, report.failure_counts
    assert report.records[0].analyte.kind is kind
    assert report.records[0].analyte.identity_key()[0] == kind.value


def test_every_analyte_kind_of_the_union_has_a_row_in_this_file():
    # Guards the parametrize above: a seventh kind added to the union must arrive
    # here with a row of its own rather than raising a KeyError nobody reads.
    assert set(KIND_ROWS) == set(AnalyteKind)


def test_a_row_naming_no_kind_at_all_is_refused_rather_than_guessed_from_its_columns():
    # A filled composition column is not a declaration that the row is a glycan.
    # Guessing would build a record the file never claimed.
    report = load(table(analyte_row("", analyte_composition="Hex5HexNAc2")))
    assert report.records_built == 0
    assert report.rows_failed == 1
    assert report.records_built + report.rows_failed == report.rows_read


def test_a_row_naming_a_kind_whose_required_fields_are_absent_is_counted_not_guessed_elsewhere():
    # It says peptide and carries an InChIKey. That is a fault to report, not a
    # small molecule to infer.
    report = load(table(analyte_row("peptide", analyte_inchikey=INCHIKEY)))
    assert report.records_built == 0
    assert report.failure_counts == {ANALYTE_REJECTED: 1}


def test_a_kind_the_union_does_not_hold_is_refused_and_the_message_names_the_kinds_it_could_be():
    report = load(table(analyte_row("protein_complex", analyte_accession="P01857")))
    assert report.records_built == 0
    assert report.rows_failed == 1
    example = next(iter(report.failure_examples.values()))[0]
    assert "protein_complex" in example
    for kind in AnalyteKind:
        assert kind.value in example


# Was an xfail against a real defect: the call site caught bare Exception, so a
# row naming a kind the union does not hold was counted as failed validation
# rather than under its own reason. _Coercion is now caught separately.
@pytest.mark.parametrize("named", ["", "protein_complex"])
def test_an_analyte_kind_the_union_does_not_hold_is_counted_under_its_own_reason(named):
    report = load(table(analyte_row(named, analyte_composition="Hex5HexNAc2")))
    assert report.failure_counts == {UNKNOWN_ANALYTE_KIND: 1}


def test_a_native_and_a_denatured_ion_of_one_antibody_are_two_matched_ions_in_one_file():
    # The loader must carry the distinction through from the file: folding state
    # is identity, so these two rows are not replicates of one number.
    text = table(
        analyte_row(
            "intact_antibody",
            analyte_inn="trastuzumab",
            analyte_folding_state="native",
            adduct="[M+24H]24+",
            charge="24",
            ccs="7000.0",
        ),
        analyte_row(
            "intact_antibody",
            analyte_inn="trastuzumab",
            analyte_folding_state="denatured",
            adduct="[M+40H]40+",
            charge="40",
            ccs="12000.0",
        ),
    )
    report = load(text)
    assert report.records_built == 2, report.failure_counts
    assert report.matched_ion_keys_held == 2


# --- ions that can never pair -----------------------------------------------------------
#
# Added after the first mutation sweep, with the behaviour itself.


def test_records_whose_charge_carrier_is_unstated_are_counted_apart_from_matched_ions():
    # Folded into the matched-ion count they would read as ions waiting for a
    # partner. A native-MS protein file would then report hundreds of distinct
    # matched ions and no possible pair, which reads as a puzzle rather than as
    # the one-line problem it is.
    unstated = analyte_row(
        "protein",
        analyte_accession="P01857",
        analyte_folding_state="native",
        adduct="[M+24?]24+",
        charge="24",
        polarity="positive",
        ccs="7000.0",
    )
    named = analyte_row(
        "protein",
        analyte_accession="P01857",
        analyte_folding_state="native",
        adduct="[M+24H]24+",
        charge="24",
        polarity="positive",
        ccs="7000.0",
    )
    report = load(table(unstated, named))
    assert report.records_built == 2
    assert report.unmatchable_held == 1
    # The named-carrier record is the only one that counts as a matched ion.
    assert report.matched_ion_keys_held == 1
    assert "UNMATCHABLE" in report.summary()
