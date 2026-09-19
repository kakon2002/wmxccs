"""The first real data: the steroid interplatform database, end to end.

Everything else in this suite runs on a corpus this repository invented. This file
runs on published measurements of real instruments, and it is the only place where
a number failing an assertion might mean the DATA is not what we think rather than
the code.

So the assertions are split deliberately:

- WHAT THE COMMITTED SEED FILE HOLDS. These run in any clone. They pin the counts,
  the coverage, the platforms, the matched ions, and which records are held from
  training and why. If the conversion changes, one of these goes red with a number
  a reviewer can compare against the paper.
- THAT NO VALUE WAS TRANSCRIBED BY HAND. Every CCS in the seed file is checked
  against the sheet dumped verbatim beside it. This is the assertion that matters
  most: a hand-typed digit is the one error that no amount of downstream statistics
  would reveal, because the wrong number is still a plausible cross section.
- THAT THE ADAPTER REFUSES A SHEET IT DOES NOT RECOGNISE. It reads columns by
  position, so this is the guard standing between a revised supporting file and a
  silent swap of one platform's values for another's.

The paper's own headline is 142 CCS values across the three technologies, which is
what `EXPECTED_DATA_ROWS` below is, and it is the one number here that was read
from the publication rather than computed from the file.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import pytest

from wmxccs.loader import load_measurements_file
from wmxccs.matching import build_matched_ions
from wmxccs.models import UNSTATED_CALIBRANT
from wmxccs.reuse import ReuseStatus
from wmxccs.sources import REGISTRY, STEROID_INTERPLATFORM_2022
from wmxccs.statistics import compare_platforms

REPO = Path(__file__).resolve().parents[1]
SEED = REPO / "data" / "seed" / "steroid_jasms2022.csv"
AS_DELIVERED = REPO / "data" / "seed" / "as_delivered" / "js2c00196_si_003_S2.csv"
DOI = "10.1021/jasms.2c00196"

# From the publication, not from the file.
EXPECTED_DATA_ROWS = 142

# Per platform, counted from the file and checked against it. The two drift-tube
# columns are short, and that unevenness is the point: it is the "present on one
# platform only" case arriving in the first real dataset rather than in theory.
EXPECTED_COVERAGE = {
    ("TWIMS", ""): 142,
    ("TIMS", ""): 142,
    ("DTIMS", "single_field"): 135,
    ("DTIMS", "stepped_field"): 102,
}
EXPECTED_RECORDS = sum(EXPECTED_COVERAGE.values())  # 521


@pytest.fixture(scope="module")
def rows() -> list[dict[str, str]]:
    with SEED.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


@pytest.fixture(scope="module")
def report():
    return load_measurements_file(SEED)


# --- what the committed seed file holds --------------------------------------------------


def test_the_seed_file_is_committed_alongside_the_sheet_it_came_from() -> None:
    """Both, or the hand-transcription check below cannot run in a clone.

    data/raw is gitignored, so the verbatim dump in data/seed/as_delivered IS the
    input as far as this repository is concerned.
    """
    assert SEED.exists(), "run tools/ingest_steroid.py"
    assert AS_DELIVERED.exists(), "run tools/ingest_steroid.py"


def test_the_file_holds_one_record_per_published_value_and_no_others(rows) -> None:
    assert len(rows) == EXPECTED_RECORDS
    coverage = Counter((row["ims_type"], row["dtims_method"]) for row in rows)
    assert dict(coverage) == EXPECTED_COVERAGE
    # One row of the sheet per compound-ion, each contributing between two and four
    # records. Conservation in the direction that catches a duplicated row.
    assert len({row["source_locator"] for row in rows}) == EXPECTED_DATA_ROWS


def test_every_row_carries_the_same_source_string_as_its_analyte(rows) -> None:
    """Not cosmetic. A component's licence claim is backed by the enclosing record's
    DOI only when the component names the enclosing record's own source, so these
    two drifting apart would refuse all 521 rows with a message about licences."""
    for row in rows:
        assert row["analyte_source"] == row["source"]
        assert row["doi"] == DOI


def test_the_registry_backs_what_every_row_claims(rows) -> None:
    assert DOI in REGISTRY
    assert STEROID_INTERPLATFORM_2022.reuse_status is ReuseStatus.ACADEMIC_ONLY
    assert {row["reuse_status"] for row in rows} == {ReuseStatus.ACADEMIC_ONLY.value}
    # And the entry records who decided this platform may use an academic-only
    # source, which is the whole basis on which these rows load at all.
    assert STEROID_INTERPLATFORM_2022.context_basis.strip()


def test_the_calibrants_are_recorded_as_the_supporting_information_states_them(rows) -> None:
    """Three different statements, recorded as three different things.

    Two platforms name a calibrant. One states only its MASS calibration, so its
    CCS calibrant is UNSTATED. One is primary and uses none, which is an empty
    calibrant and NOT the same claim as UNSTATED.
    """
    by_platform = {
        (row["ims_type"], row["dtims_method"]): row["calibrant"] for row in rows
    }
    assert by_platform[("TWIMS", "")] == "Waters Major Mix"
    assert by_platform[("DTIMS", "single_field")] == "Agilent ESI-L tune mix (G1969-85000)"
    assert by_platform[("TIMS", "")] == UNSTATED_CALIBRANT
    assert by_platform[("DTIMS", "stepped_field")] == ""


def test_the_gas_the_values_refer_to_is_stated_and_the_gas_in_the_cell_is_not(rows) -> None:
    """TWCCSN2, TIMCCSN2, DTCCSN2: the column names state the reference gas.

    What gas was in each cell they do not state, and the two are different facts.
    Filling cell_gas from the reference gas would be inventing a method detail.
    """
    assert {row["drift_gas"] for row in rows} == {"N2"}
    assert {row["cell_gas"] for row in rows} == {""}


def test_no_compound_carries_an_inchikey_and_every_one_carries_a_dataset_scoped_id(rows) -> None:
    """The sheet identifies nothing structurally, and nothing here pretends otherwise.

    Resolving 87 names against a structure database is a curation act with a
    provenance trail. Until somebody does it, these compounds pair inside this
    source and nowhere else - which is correct, not a shortfall.
    """
    assert {row["analyte_inchikey"] for row in rows} == {""}
    for row in rows:
        assert row["analyte_dataset_compound_id"].startswith("steroid_jasms2022:")


# --- that no value was transcribed by hand ------------------------------------------------


def test_every_ccs_in_the_seed_file_appears_in_the_sheet_it_came_from(rows) -> None:
    """The assertion that matters most in this file.

    A mistyped digit produces a number that is still a plausible cross section, so
    no statistic downstream would reveal it - the bias would simply be wrong, and
    wrong in a way that looks like a finding. This compares every value against the
    sheet dumped cell for cell.

    Matched on the text openpyxl read, at the precision the adapter writes, so that
    a value which merely ROUNDS to a real cell does not pass: the comparison set is
    built by formatting the delivered cells the same way.
    """
    delivered: set[str] = set()
    with AS_DELIVERED.open(newline="", encoding="utf-8") as handle:
        for line in csv.reader(handle):
            for cell in line:
                try:
                    delivered.add(f"{float(cell):.6g}")
                except ValueError:
                    continue
    assert delivered, "the verbatim dump holds no numbers"
    for row in rows:
        assert row["ccs"] in delivered, f"{row['ccs']} is in no cell of the sheet"
        if row["ccs_uncertainty"]:
            assert row["ccs_uncertainty"] in delivered


def test_every_compound_name_in_the_seed_file_appears_in_the_sheet(rows) -> None:
    """The same check for the identities, including the conformation asterisk.

    The asterisk is kept rather than stripped: two compounds carry it, and removing
    it would merge a starred compound with an unstarred one of the same name if the
    sheet ever listed both.
    """
    delivered = set()
    with AS_DELIVERED.open(newline="", encoding="utf-8") as handle:
        for line in csv.reader(handle):
            delivered.update(cell.strip() for cell in line)
    starred = 0
    for row in rows:
        name = row["analyte_dataset_compound_id"].removeprefix("steroid_jasms2022:")
        assert name in delivered, f"{name!r} is in no cell of the sheet"
        starred += "*" in name
    assert starred > 0, "the conformation asterisk is no longer being carried through"


# --- that the adapter refuses a sheet it does not recognise --------------------------------


def test_the_adapter_refuses_a_sheet_whose_columns_have_moved() -> None:
    """The guard standing between a revised supporting file and a silent platform swap.

    Columns are read by POSITION. If a future version of the supporting information
    inserts one, every travelling-wave value would be recorded as a trapped-ion one
    and the resulting bias figure would look like a result.
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location("ingest_steroid", REPO / "tools" / "ingest_steroid.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    header = [None, None, None, None, None] + [None] * 13
    for column, text in module.IDENTITY_HEADERS:
        header[column - 1] = text
    for entry in module.PLATFORMS:
        header[entry[0] - 1] = entry[2]
    header[module.NOT_INGESTED[0] - 1] = module.NOT_INGESTED[1]
    module.check_headers(tuple(header))  # the real headers pass

    # One column shifted: refused, and the message says which.
    moved = list(header)
    moved.insert(6, "an inserted column")
    with pytest.raises(module.SheetChanged) as refusal:
        module.check_headers(tuple(moved))
    assert "not where this adapter reads them" in str(refusal.value)
    assert "will not guess" in str(refusal.value)


# --- end to end: load, match, compare -----------------------------------------------------


def test_every_row_loads_and_nothing_is_lost(report) -> None:
    assert report.rows_read == EXPECTED_RECORDS
    assert report.rows_failed == 0
    assert report.failure_counts == {}
    assert report.records_built == EXPECTED_RECORDS
    assert report.unknown_columns == ()


def test_the_trapped_ion_records_are_held_from_training_by_the_unstated_calibrant(report) -> None:
    """142 real values, held rather than refused, and held for one stated reason.

    This is the single largest thing standing between this repository and a
    three-technology harmonization, and it is one paragraph of somebody's methods
    away. Recorded as a held record with a reason, not dropped.
    """
    held = [record for record in report.records if record not in report.cleared]
    tims = [record for record in report.records if record.ims_type.value == "TIMS"]
    assert len(tims) == 142
    assert all(record.calibrant == UNSTATED_CALIBRANT for record in tims)
    assert not set(tims) & set(report.cleared), "no unstated-calibrant record may train"
    assert len(held) == 146  # the 142, plus 4 the shared-peak check flagged


def test_the_shared_peak_check_fires_on_real_data(report) -> None:
    """It does, and on two genuinely different steroids reporting one cross section.

    Not a defect in the check and not necessarily one in the paper: two isomers can
    agree to two decimal places. It is held for a person to look at, which is what
    the check is for.
    """
    reasons = " | ".join(report.gate_counts)
    assert "shared peak" in reasons
    assert report.gate_counts["suspected shared peak: identical CCS for different analytes in one calibration group"] == 4


def test_the_real_matched_ions(report) -> None:
    """THE FIRST REAL RESULT THIS PROJECT HAS. 142 cross-platform matched ions.

    Every one is a genuine pairing of one ion measured on two to four technologies,
    asserted by the paper itself rather than inferred by us from a compound name.
    Nothing here is synthetic and nothing is held: `usable` is all 142.
    """
    matching = build_matched_ions(report.records)
    assert len(matching.matched) == 142
    assert matching.single_platform == ()
    assert matching.unmatchable == ()
    assert matching.duplicates_collapsed == 0
    assert len(matching.usable) == 142, "an academic_only member must not block its set"
    assert matching.widest_match == 4
    assert dict(sorted(Counter(len(ion.platforms) for ion in matching.matched).items())) == {2: 2, 3: 43, 4: 97}
    assert matching.quotable and matching.refusal() is None


def test_the_matched_ions_that_may_actually_train(report) -> None:
    """Fewer, because the trapped-ion values are held: 140 across three platforms.

    Two ions drop to a single platform once the trapped-ion value is held, which is
    why this is asserted separately from the count above. Both numbers are true and
    they answer different questions.
    """
    matching = build_matched_ions(report.cleared)
    assert len(matching.matched) == 140
    assert len(matching.single_platform) == 2
    assert len(matching.usable) == 140
    assert matching.widest_match == 3


def test_the_comparison_is_stratified_by_calibration_group_and_never_pooled(report) -> None:
    """Three platform pairs, each split by adduct, and no pooled figure for any of them.

    The adducts land in different calibration groups, so there is no single number
    for "travelling wave against drift tube" here - and the bias really does differ
    between them, which is why refusing to pool was worth building.
    """
    comparison = compare_platforms(build_matched_ions(report.cleared))
    assert comparison.quotable, "these are real instruments and may be quoted"
    assert comparison.refusal() is None
    assert len(comparison.pairs) == 3
    assert comparison.ions_considered == 140
    for pair in comparison.pairs:
        assert len(pair.strata) == 3, "one stratum per adduct"
        assert pair.pooled is None, "values calibrated differently are not one comparison"


def test_the_outliers_are_reported_and_kept(report) -> None:
    """Five, every one kept in the data, and the two worst are the same two compounds twice.

    An ester and a diglucuronide - the largest and most flexible ions in the set -
    disagree between platforms by 2% and 7%. That is a real finding about where
    harmonization is hard, and it is exactly what removing outliers would erase.
    """
    comparison = compare_platforms(build_matched_ions(report.cleared))
    assert len(comparison.outliers) == 5
    named = {point.ion.key.analyte for point in comparison.outliers}
    assert len(named) == 3, "the same few compounds recur across strata"
    # and nothing was dropped to achieve that
    assert comparison.n_points == 326


def test_the_uncertainty_columns_hold_four_different_things_and_all_are_accounted_for(rows) -> None:
    """Messy reality, pinned so it cannot be flattened away later.

    Every standard-deviation cell of an ingested value is one of four things, and
    three of them become "no uncertainty reported":

      a number     a real spread over the stated n
      exactly 0    6 cells. Either the replicates agreed to the reported precision
                   or nothing was computed; the sheet does not say which, and the
                   model refuses 0 as a spread in any case
      "n.d."       23 cells, every one a negative-mode stepped-field row whose n is
                   1. Not determined, and not determinable from one measurement -
                   a coherent absence rather than a gap in the data
      blank        none in this file

    The counts are asserted because a future change that started recording a zero as
    a spread, or that let an "n.d." through as text, would otherwise be invisible.
    """
    missing = [row for row in rows if not row["ccs_uncertainty"]]
    assert len(missing) == 29, "6 zeros plus 23 'n.d.'"
    # The 23 are exactly the negative-mode stepped-field rows, where n is 1.
    nd = [row for row in missing if row["dtims_method"] == "stepped_field" and row["polarity"] == "negative"]
    assert len(nd) == 23
    assert {row["replicates"] for row in nd} == {"1"}, "n=1 cannot have a standard deviation"
    # Everything that DOES report a spread reports it as a standard deviation, and
    # none of them is zero: the model refuses that, so a zero reaching it would be a
    # loader failure rather than a held row.
    reported = [row for row in rows if row["ccs_uncertainty"]]
    assert len(reported) == 492
    assert {row["uncertainty_type"] for row in reported} == {"sd"}
    assert all(float(row["ccs_uncertainty"]) > 0 for row in reported)
