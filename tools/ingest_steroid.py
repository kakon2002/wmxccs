"""Convert the steroid interplatform database into the canonical wmxccs seed format.

THE SOURCE
----------
Feuerstein et al., "Comparability of Steroid Collision Cross Sections Using Three
Different IM-HRMS Technologies: An Interplatform Study", J. Am. Soc. Mass
Spectrom. 2022, DOI 10.1021/jasms.2c00196. Supporting file js2c00196_si_003.xlsx,
sheet "S2_Interplatform CCS Database".

Downloaded 19 September 2026 from Europe PMC's supplementaryFiles REST endpoint
for PMC9545150, which is a documented open API. The publisher's own copy answers
403 to an automated request, and PMC serves a proof-of-work challenge for binary
downloads; neither was worked around. Nothing here was transcribed by hand.

WHY THIS DATASET AND NOT ANOTHER
--------------------------------
It is the only source identified so far that is cross-platform BY DESIGN: one row
per ion, with that ion's cross section from several technologies side by side. Every
other candidate is one platform deep, so pairing across it would depend on
resolving compound identities between papers, which nobody has done. Here the
paper itself asserts that the values in a row are of the same ion, and that
assertion is the warrant - not an inference of ours.

WHAT THE SHEET ACTUALLY LOOKS LIKE, and it is not tidy
------------------------------------------------------
- TWO header rows. Row 1 groups the columns by technology AND doubles as a
  category banner ("ANDROGENS AND RELATED COMPOUNDS"); row 2 holds the column
  names. Row 2 has text in the m/z column, so "has an m/z" alone does not
  identify a data row.
- A FOOTNOTE row at 145, then four empty rows. All five are accounted for below
  rather than skipped by a truthiness test that would also have swallowed a real
  row with a blank cell.
- FIVE platform columns, not three. Travelling wave appears twice - once as one
  laboratory's value and once as an average over four instruments - and the drift
  tube appears twice, single-field and stepped-field. That the two drift-tube
  columns are different platforms is not a detail: only the stepped-field value is
  primary, and this repository already treats them as two platforms for exactly
  that reason.
- COVERAGE IS UNEVEN. Travelling wave and trapped ion cover all 142 ions;
  single-field drift tube covers 135; stepped-field covers 102. Forty ions have no
  primary value at all, which is the "present on one platform only" problem
  arriving immediately rather than theoretically.

WHAT THIS PRODUCED
------------------
142 cross-platform matched ions - 2 across two platforms, 43 across three, 97 across
all four - from 521 records. The first result in this project that describes real
instruments. Of the 521, 375 may train: the 142 trapped-ion values are held by their
unstated calibrant and 4 more by the shared-peak check, which fired on real data for
the first time on two isomers reporting an identical cross section.

THE THREE JUDGEMENT CALLS, all reversible, none silent
------------------------------------------------------
1. THE SINGLE-LABORATORY TRAVELLING-WAVE COLUMN IS NOT INGESTED. Ingesting both
   travelling-wave columns would put two values of one ion on one platform in one
   calibration group, which is a replicate pair - and matched-ion construction
   holds an ion with replicates on one side rather than averaging or picking one.
   Every ion would be held and the dataset would yield nothing. The
   cross-laboratory column is ingested because it is the one the paper presents as
   its interplatform reference, and it carries its own n and standard deviation.
   The other column is counted and reported, never silently dropped.
2. THE TRAPPED-ION CALIBRANT IS RECORDED AS UNSTATED, and the paper does not merely
   omit it - it refers to a calibrant it never names. The supporting information
   states the travelling-wave calibrant (Waters Major Mix) and the single-field
   drift-tube calibrant (Agilent ESI-L G1969-85000) explicitly. For trapped ion
   mobility, page 7 gives only the MASS calibration - "The TOFMS was calibrated
   using 10 mM sodium formate and a 7th order high-performance calibration" - and
   then says:

       "In addition to the external calibration, each sample was automatically
        post-run calibrated by injecting a 1:1 mixture of BOTH calibrants."

   Both calibrants. One is named. The sentence presupposes a second one, which is
   the mobility calibrant, and the methods never identify it.

   The likely answer is not a mystery: the background section calls an Agilent tune
   mix typical for trapped ion mobility, and this same laboratory used exactly that
   mix for its single-field drift-tube work. But "typical" with an "e.g." is not a
   statement about this experiment, and an obvious inference written into a record as
   though the paper stated it is indistinguishable, downstream, from a fact. So the
   calibrant is UNSTATED: those 142 rows are held out of any fit, stay visible and
   countable, and key apart from every named calibrant so they cannot pool with one.

   This costs more than it appears to. Computed descriptively, the trapped-ion values
   agree with the drift tube BETTER than the travelling-wave ones do. The platform
   currently excluded is the one that agrees best, and one paragraph of somebody's
   methods would admit it. Recorded in the registry entry's what_to_check.
3. COMPOUNDS ARE IDENTIFIED BY THIS DATASET'S OWN NAMES. The sheet carries a
   compound name, a commercial name, a formula and an m/z, and nothing that
   identifies a structure. No InChIKey is invented and none is looked up: a
   name-to-structure lookup fails silently and wrongly on exactly the compounds
   that matter here, which are isomers with similar names. So each compound is
   recorded as "steroid_jasms2022:<name>", which is a real identity inside this
   dataset and matches nothing outside it until a person resolves it.

Run: python tools/ingest_steroid.py
"""

from __future__ import annotations

import csv
import pathlib
import sys
from collections import Counter

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import openpyxl  # noqa: E402

from wmxccs.loader import COLUMNS  # noqa: E402
from wmxccs.models import UNSTATED_CALIBRANT  # noqa: E402
from wmxccs.reuse import ReuseStatus  # noqa: E402
from wmxccs.sources import REGISTRY, STEROID_INTERPLATFORM_2022  # noqa: E402

SOURCE_FILE = REPO / "data" / "raw" / "steroid_jasms2022" / "js2c00196_si_003.xlsx"
SHEET = "S2_Interplatform CCS Database"
TARGET = REPO / "data" / "seed" / "steroid_jasms2022.csv"

# The sheet dumped cell for cell, no interpretation, committed beside the converted
# file. data/raw is gitignored precisely so that the delivered data reaches the
# repository as a reviewable artefact in data/seed instead, and this is that
# artefact for a source that arrives as a binary workbook: it makes the conversion
# diffable against its input, and it lets the check that no value here was
# transcribed by hand run from a clone with no raw directory at all.
AS_DELIVERED = REPO / "data" / "seed" / "as_delivered" / "js2c00196_si_003_S2.csv"

# Taken FROM the registry rather than retyped beside it. The gate backs a record's
# licence claim by looking its DOI up in the registry, so a DOI written here by
# hand that drifted from the entry by one character would refuse all 421 rows for
# a reason that looks like a licence problem and is a typo.
DOI = STEROID_INTERPLATFORM_2022.doi
STATUS = STEROID_INTERPLATFORM_2022.reuse_status.value
assert DOI in REGISTRY, DOI
assert STEROID_INTERPLATFORM_2022.reuse_status is ReuseStatus.ACADEMIC_ONLY

# The analyte carries this too, and the two must be the SAME STRING. A component's
# licence claim is backed by the enclosing record's DOI only when the component
# names the enclosing record's own source; otherwise any matching-looking string
# would let one record's registry entry vouch for another's analyte.
CITATION = "Feuerstein et al., J. Am. Soc. Mass Spectrom. 2022"

# The calibrants, each as the supporting information states it. See the module
# docstring for why the trapped-ion one is UNSTATED.
MAJOR_MIX = "Waters Major Mix"
ESI_L = "Agilent ESI-L tune mix (G1969-85000)"

# (ccs column, sd column, header the ccs column must carry, ims_type, dtims_method,
#  calibrant, replicates, what it is)
#
# COLUMNS ARE READ BY POSITION, so each entry also carries the header text that
# column must hold, and convert() refuses the sheet if it does not. Without that a
# revised supporting file with one column inserted would be ingested silently, with
# travelling-wave values recorded as trapped-ion ones and a bias figure to match.
# Nothing about the failure would look like a failure. The header substring is also
# where the replicate counts below come from - "(n=12; four TWIMS platforms)",
# "(n=3-6)", "(n=3)", "(n=3 pos; n=1 neg)" - so a changed n cannot slip past either.
PLATFORMS = (
    (
        7, 8, "Average TWCCSN2 (Å2) (n=12; four TWIMS platforms)",
        "TWIMS", "", MAJOR_MIX, 12, "cross-laboratory average over four Waters TWIM-MS systems",
    ),
    (
        10, 11, "Average TIMCCSN2 (Å2) (n=3-6)",
        # n is a RANGE in the header, so no single replicate count is true of every
        # row. Left absent rather than reduced to one end of it.
        "TIMS", "", UNSTATED_CALIBRANT, None, "timsTOF pro",
    ),
    (
        13, 14, "Average DTCCSN2 (Å2) (n=3)",
        "DTIMS", "single_field", ESI_L, 3, "Agilent DTIM-MS, single-field calibration, 4-bit",
    ),
    (
        16, 17, "Average DTCCSN2 (Å2) (n=3 pos; n=1 neg)",
        # No calibrant: the stepped-field method is primary and is not calibrated
        # against a reference set. An empty calibrant is a different statement from
        # UNSTATED - this source says what this method does, and it does not use one.
        "DTIMS", "stepped_field", "", None, "Agilent DTIM-MS, stepped-field method",
    ),
)

# Not ingested, and reported. See judgement call 1. Its header is checked too, so
# that the column being skipped is verifiably the one meant to be skipped.
NOT_INGESTED = (6, "TWCCSN2 (Å2)a", "TWIMS single-laboratory database (one Synapt G2-S HDMS)")

# The identifying columns, checked for the same reason.
IDENTITY_HEADERS = ((1, "Compound name"), (2, "Commercial name"), (3, "Formula"), (4, "Ion"), (5, "m/z"))


class SheetChanged(SystemExit):
    """The sheet is not the one this adapter was written against."""


def check_headers(header_row: tuple) -> None:
    """Refuse a sheet whose columns have moved. Positional reads demand this."""
    expected = [(col, text) for col, text in IDENTITY_HEADERS]
    expected += [(entry[0], entry[2]) for entry in PLATFORMS]
    expected.append((NOT_INGESTED[0], NOT_INGESTED[1]))
    wrong = []
    for column, text in expected:
        found = header_row[column - 1] if column - 1 < len(header_row) else None
        # Compared on stripped text, because the sheet pads some headers with a
        # trailing space, and a trailing space is not a changed column.
        if str(found or "").strip() != text.strip():
            wrong.append(f"column {column}: expected {text!r}, found {found!r}")
    if wrong:
        raise SheetChanged(
            "the sheet's columns are not where this adapter reads them, so it will not guess:\n  "
            + "\n  ".join(wrong)
            + "\nRe-read the sheet and correct PLATFORMS, rather than loosening this check."
        )

POLARITY = {"+": "positive", "-": "negative"}

# The line ending, named once and built without a backslash escape, so that both
# writers below use the same one and neither can drift from the other.
NEWLINE = chr(10)


def _polarity_of(adduct: str) -> str:
    return POLARITY[adduct.strip()[-1]]


def _charge_of(adduct: str) -> int:
    """Signed charge from the bracket suffix. These are all singly charged."""
    body = adduct.strip()
    sign = 1 if body.endswith("+") else -1
    magnitude = body[body.rindex("]") + 1 : -1]
    return sign * (int(magnitude) if magnitude else 1)


def _blank_row() -> dict[str, str]:
    return {column: "" for column in COLUMNS}


def _number(value: object) -> str:
    """A cell as text, or empty. Only a real number becomes a number."""
    if value is None or isinstance(value, bool):
        return ""
    if isinstance(value, (int, float)):
        return f"{value:.6g}"
    return ""


def convert() -> tuple[list[dict[str, str]], dict[str, int], list[str]]:
    workbook = openpyxl.load_workbook(SOURCE_FILE, data_only=True)
    sheet = workbook[SHEET]
    rows = list(sheet.iter_rows(values_only=True))

    counts: Counter[str] = Counter()
    notes: list[str] = []
    out: list[dict[str, str]] = []

    for line, row in enumerate(rows, start=1):
        counts["rows in the sheet"] += 1
        # Account for every row explicitly. A truthiness test would also swallow a
        # real row whose first cell happened to be blank.
        if line == 1:
            counts["header: platform grouping and category banner"] += 1
            continue
        if line == 2:
            check_headers(row)
            counts["header: column names, all checked to be where this adapter reads them"] += 1
            continue
        if all(cell is None for cell in row):
            counts["blank"] += 1
            continue
        if row[4] is None:
            counts["footnote or note row"] += 1
            notes.append(f"line {line}: {str(row[0])[:100]}")
            continue

        name, commercial, formula, adduct, mz = (row[0], row[1], row[2], row[3], row[4])
        if not (name and adduct):
            counts["skipped: no compound name or no ion"] += 1
            continue
        counts["data rows"] += 1

        compound = str(name).strip()
        # The footnote at line 145 says an asterisk marks a compound with multiple
        # conformations observed. The marker is kept in the recorded identity, because
        # stripping it would silently merge a starred compound with an unstarred one
        # if both appear, and reported here so nobody has to guess what it means.
        if "*" in compound:
            counts["data rows whose compound name carries the conformation asterisk"] += 1

        polarity = _polarity_of(str(adduct))
        charge = _charge_of(str(adduct))

        for ccs_col, sd_col, _header, ims_type, method, calibrant, replicates, what in PLATFORMS:
            ccs = row[ccs_col - 1]
            if ccs is None:
                counts[f"no value on {ims_type}{'/' + method if method else ''}"] += 1
                continue
            if not isinstance(ccs, (int, float)) or isinstance(ccs, bool):
                # Prospective. Every CCS cell in this file is a number; if one were
                # not, _number would quietly return "" and the row would leave here
                # with an empty required field, to be refused downstream with a
                # message about a missing number rather than about a changed sheet.
                raise SheetChanged(
                    f"line {line}, column {ccs_col}: a cross section that is not a number: {ccs!r}."
                    " Read the sheet and decide what it means; do not let it become an absent value."
                )
            sd = row[sd_col - 1]
            record = _blank_row()
            record.update(
                {
                    "analyte_kind": "small_molecule",
                    "analyte_dataset_compound_id": f"steroid_jasms2022:{compound}",
                    "analyte_display_name": str(commercial).strip() if commercial else compound,
                    "analyte_source": CITATION,
                    "analyte_reuse_status": STATUS,
                    "ccs": _number(ccs),
                    "adduct": str(adduct).strip(),
                    "charge": str(charge),
                    "polarity": polarity,
                    "ims_type": ims_type,
                    "dtims_method": method,
                    # Stated in the column names themselves: TWCCSN2, TIMCCSN2,
                    # DTCCSN2. The gas in the cell is not stated, so cell_gas is
                    # left blank rather than assumed to be the same.
                    "drift_gas": "N2",
                    "calibrant": calibrant,
                    "instrument": what,
                    "source": CITATION,
                    "doi": DOI,
                    "source_locator": f"SI_3 sheet {SHEET}, line {line}",
                    "reuse_status": STATUS,
                }
            )
            # THE STANDARD DEVIATION COLUMNS HOLD FOUR DIFFERENT THINGS, and all but
            # the first become "no uncertainty reported". Which one it was is counted
            # rather than flattened, because they do not mean the same:
            #   a number   - a real spread over the stated n
            #   exactly 0  - 6 cells. Either the replicates agreed to the reported
            #                precision or nothing was computed; the sheet does not
            #                say which, and the model refuses 0 as a spread anyway
            #   "n.d."     - 23 cells, all of them the negative-mode stepped-field
            #                rows, whose n is 1. Not determined, and not determinABLE
            #                from one measurement: a coherent absence, not a gap
            #   blank      - none in this file, counted in case a revision has some
            platform = f"{ims_type}{'/' + method if method else ''}"
            if isinstance(sd, bool) or sd is None:
                counts[f"{platform}: no SD reported (blank cell)"] += 1
                spread = ""
            elif isinstance(sd, (int, float)):
                if sd > 0:
                    spread = _number(sd)
                else:
                    counts[f"{platform}: SD reported as exactly {sd!r}, recorded as no spread"] += 1
                    spread = ""
            else:
                counts[f"{platform}: SD reported as the text {str(sd).strip()!r}, recorded as no spread"] += 1
                spread = ""
            if spread:
                record["ccs_uncertainty"] = spread
                record["uncertainty_type"] = "sd"
            if replicates is not None:
                record["replicates"] = str(replicates)
            elif ims_type == "DTIMS" and method == "stepped_field":
                # The column header states n=3 for positive mode and n=1 for negative.
                record["replicates"] = "3" if polarity == "positive" else "1"
            out.append(record)
            counts[f"records emitted for {ims_type}{'/' + method if method else ''}"] += 1

    counts[f"column {NOT_INGESTED[0]} ({NOT_INGESTED[1]}) NOT ingested: {NOT_INGESTED[2]}"] = counts["data rows"]
    return out, dict(counts), notes


def dump_as_delivered() -> int:
    """Write the sheet out verbatim. Every cell as read, nothing dropped, nothing parsed."""
    sheet = openpyxl.load_workbook(SOURCE_FILE, data_only=True)[SHEET]
    rows = list(sheet.iter_rows(values_only=True))
    width = max(len(row) for row in rows)
    AS_DELIVERED.parent.mkdir(parents=True, exist_ok=True)
    with AS_DELIVERED.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator=NEWLINE)
        for row in rows:
            padded = list(row) + [None] * (width - len(row))
            writer.writerow(["" if cell is None else str(cell) for cell in padded])
    return len(rows)


def main() -> int:
    delivered = dump_as_delivered()
    records, counts, notes = convert()
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    with TARGET.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(COLUMNS), lineterminator=NEWLINE)
        writer.writeheader()
        writer.writerows(records)

    # Row conservation, asserted rather than hoped for.
    with TARGET.open(newline="", encoding="utf-8") as handle:
        written = list(csv.DictReader(handle))
    if len(written) != len(records):
        raise SystemExit(f"wrote {len(written)} rows, built {len(records)}")

    print(f"{SOURCE_FILE.name} sheet {SHEET!r}")
    print(f"  -> {AS_DELIVERED.relative_to(REPO)}  ({delivered} rows, verbatim)")
    print(f"  -> {TARGET.relative_to(REPO)}")
    for label, count in sorted(counts.items()):
        print(f"  {label}: {count}")
    if notes:
        print("  note rows carried through to this report rather than dropped:")
        for note in notes:
            print(f"    {note}")
    print(f"  records written: {len(written)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
