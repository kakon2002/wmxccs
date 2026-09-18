"""Convert the two transcribed Struwe CSVs into the canonical wmxccs seed format.

WHY A CONVERSION STEP AND NOT A SECOND LOADER
---------------------------------------------
loader.py says it in its own docstring: one format, one loader, one set of
checks, so a value cannot arrive by a route that skips the licence gate or the
validators. A per-source adapter reading its own columns straight into records is
exactly that second route. So each source is converted ONCE, by this script, into
a file a person can read and review, and that file then goes through the same
loader as everything else.

The conversion is the thing that has to be trusted, so it is made checkable:

- row conservation is asserted, not hoped for. Rows in must equal rows out, and
  the script fails loudly rather than writing a short file. A swallowed row that
  reports success is one of the two failure classes this project keeps meeting.
- the raw transcriptions are copied into data/seed verbatim, filename and all,
  beside the converted files. The transcription as delivered stays reviewable,
  and the conversion stays diffable against it.
- THE LICENCE COLUMN IN THE FILE IS NOT BELIEVED. A licence claim in a data cell
  is not a licence. The DOI is looked up in the registry, the registry's status
  is what is written, and a row whose DOI is not on record stops the conversion.
  The `licence` cell is compared against the registry only so that a disagreement
  is reported; it never decides anything.

WHAT IS DELIBERATELY NOT DERIVED
--------------------------------
- composition for the 2016 file. Those rows carry an IUPAC-condensed structure
  and no composition column. Deriving one means mapping monosaccharide names onto
  residue classes, which is glycan chemistry and does not come across. The
  structure is the finer identifier anyway, so composition stays null.
- the gas the 2015 values refer to. The ESI says the travelling-wave values were
  calibrated against a dextran ladder of known DTCCS and never says which gas
  those reference values were in. It is recorded as UNSTATED, not inferred from
  the same group's 2016 paper.
- the platform the calibration references were measured on. "DTCCS" in the 2015
  cell plainly suggests a drift tube, and that is still a reading of an
  abbreviation rather than of a paper, so calref_platform stays null on both
  files. Resolving it means reading Hofmann 2014, Anal. Chem. 86, 10789.

Run: python tools/seed_struwe.py
"""

from __future__ import annotations

import csv
import re
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from wmxccs.loader import COLUMNS  # noqa: E402
from wmxccs.reuse import ReuseStatus  # noqa: E402
from wmxccs.sources import SEED_STRUWE_2015, SEED_STRUWE_2016, licence_for  # noqa: E402

SOURCE_DIR = Path(r"C:/Users/User/Project2/data/raw")
SEED_DIR = REPO / "data" / "seed"
AS_DELIVERED = SEED_DIR / "as_delivered"

# As delivered, space included. The 2016 file really is named with a space
# before the extension, and normalising it silently would mean a later loader
# looks for a file that is not there and reports a clean run over zero rows.
FILE_2016 = "struwe2016_ccs .csv"
FILE_2015 = "struwe2015_analyst_ccs.csv"

DOI_2016 = "10.1039/c6cc06247d"
DOI_2015 = "10.1039/c5an01092f"

# A linkage is fully resolved only when its anomeric configuration is stated and
# nothing in it is a question mark. "Gal(b1-4)Glc" is resolved; "Gal(?1-4)Glc"
# and "Gal(1-4)Glc" are not. Used to decide the two resolution flags from the
# string itself rather than asserting them.
_LINKAGE = re.compile(r"\(([ab?]?)([0-9?]+)-([0-9?]+)(?:/[0-9?]+)*\)")


def _resolution_flags(iupac: str) -> tuple[bool, bool]:
    """(has_unresolved_linkage, has_unresolved_anomericity) read off the string.

    Returns "unresolved" for a string with no linkage in it at all, because a
    string that states no linkage is not evidence that the linkages are known.
    """
    found = _LINKAGE.findall(iupac or "")
    if not found:
        return True, True
    unresolved_linkage = any("?" in positions or "?" in target for _anomer, positions, target in found)
    unresolved_anomericity = any(anomer not in ("a", "b") for anomer, _positions, _target in found)
    return unresolved_linkage, unresolved_anomericity


def _registry_status(doi: str, licence_cell: str, where: str) -> ReuseStatus:
    """The registry's status for this DOI. The row's own licence cell decides nothing."""
    entry = licence_for(doi)
    if entry is None:
        raise SystemExit(
            f"{where}: DOI {doi!r} has no licence record. A transcription cannot be seeded on the strength of"
            f" a {licence_cell!r} cell: record the licence in sources.py, read from the publisher's page,"
            " with who read it and when"
        )
    return entry.reuse_status


def _blank_row() -> dict[str, str]:
    return {column: "" for column in COLUMNS}


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [{k: (v or "").strip() for k, v in row.items()} for row in csv.DictReader(handle)]


def convert_2016(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for index, cell in enumerate(rows, start=2):
        where = f"{FILE_2016} line {index}"
        if cell["doi"].strip().casefold() != DOI_2016:
            raise SystemExit(f"{where}: DOI {cell['doi']!r} is not this file's DOI {DOI_2016!r}")
        status = _registry_status(DOI_2016, cell.get("licence", ""), where)
        unresolved_linkage, unresolved_anomericity = _resolution_flags(cell["structure_iupac"])
        row = _blank_row()
        row.update(
            {
                "analyte_kind": "glycan",
                "analyte_display_name": cell["compound"],
                # The dataset provenance, exactly as sources.py records it. The
                # licence gate matches this string character for character.
                "analyte_source": SEED_STRUWE_2016.provenance,
                "analyte_reuse_status": status.value,
                "analyte_iupac_condensed": cell["structure_iupac"],
                # No composition column in this file, and none is derived. See the module docstring.
                "analyte_composition": "",
                "analyte_has_unresolved_linkage": "true" if unresolved_linkage else "false",
                "analyte_has_unresolved_anomericity": "true" if unresolved_anomericity else "false",
                "analyte_reducing_end_label": cell["reducing_end_label"],
                "analyte_derivatisation": cell["derivatisation"],
                "ccs": cell["ccs_A2"],
                # The owner's instruction: this column reports two standard deviations.
                "ccs_uncertainty": cell["ccs_2sd"],
                "uncertainty_type": "two_sd" if cell["ccs_2sd"] else "",
                "adduct": cell["adduct"],
                "charge": cell["charge"],
                "polarity": cell["polarity"],
                "ims_type": cell["ims_type"],
                "drift_gas": cell["drift_gas"],
                "cell_gas": cell["cell_gas"],
                "calibrant": cell["calibrant"],
                "calref_reference_set": cell["calibrant_reference"],
                "instrument": cell["instrument"],
                "replicates": cell["replicates"],
                "source": cell["source_id"],
                "doi": cell["doi"],
                "source_locator": cell["source_locator"],
                "reuse_status": status.value,
            }
        )
        out.append(row)
    return out


def convert_2015(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for index, cell in enumerate(rows, start=2):
        where = f"{FILE_2015} line {index}"
        if cell["doi"].strip().casefold() != DOI_2015:
            raise SystemExit(f"{where}: DOI {cell['doi']!r} is not this file's DOI {DOI_2015!r}")
        status = _registry_status(DOI_2015, cell.get("licence", ""), where)
        # Sample origin is folded into the locator, exactly as the glycan
        # platform's adapter did, and for a reason worth restating: the eleven
        # conformer pairs in this file are the same four compounds measured from
        # four sample origins, and the ONLY thing keeping one origin's pair from
        # pooling with another's is the locator they are compared on. Give sample
        # origin a column of its own and stop folding it in, and the lone-conformer
        # check silently stops separating them. There is a test for this.
        locator = cell["source_locator"]
        if cell.get("sample_origin"):
            locator = f"{locator} ({cell['sample_origin']})"
        row = _blank_row()
        row.update(
            {
                "analyte_kind": "glycan",
                "analyte_display_name": cell["compound"],
                "analyte_source": SEED_STRUWE_2015.provenance,
                "analyte_reuse_status": status.value,
                "analyte_composition": cell["composition"],
                # structure_iupac is empty in every row of this file: the structures
                # sit in the paper's Figure 1, which nobody has read here, and the
                # structure_status column says so. Nothing is filled in from memory.
                "analyte_iupac_condensed": cell["structure_iupac"],
                "analyte_has_unresolved_linkage": "true",
                "analyte_has_unresolved_anomericity": "true",
                "analyte_reducing_end_label": cell["reducing_end_label"],
                "analyte_derivatisation": cell["derivatisation"],
                "ccs": cell["ccs_A2"],
                "ccs_uncertainty": cell["ccs_uncertainty"],
                "uncertainty_type": cell["uncertainty_type"],
                "conformer": cell["conformer"],
                "conformers_total": cell["conformers_total"],
                "adduct": cell["adduct"],
                "charge": cell["charge"],
                "polarity": cell["polarity"],
                "ims_type": cell["ims_type"],
                "drift_gas": cell["drift_gas"],
                "cell_gas": cell["cell_gas"],
                "calibrant": cell["calibrant"],
                "calref_reference_set": cell["calibrant_reference"],
                "instrument": cell["instrument"],
                "source": cell["source_id"],
                "doi": cell["doi"],
                "source_locator": locator,
                "reuse_status": status.value,
            }
        )
        out.append(row)
    return out


def write(rows: list[dict[str, str]], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(COLUMNS), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    SEED_DIR.mkdir(parents=True, exist_ok=True)
    AS_DELIVERED.mkdir(parents=True, exist_ok=True)

    for name, converter, out_name in (
        (FILE_2016, convert_2016, "struwe2016_chemcommun.csv"),
        (FILE_2015, convert_2015, "struwe2015_analyst.csv"),
    ):
        source = SOURCE_DIR / name
        if not source.exists():
            raise SystemExit(f"missing transcription: {source}")
        # Verbatim, filename and all, so the delivered transcription stays reviewable.
        shutil.copy2(source, AS_DELIVERED / name)

        rows = _read(source)
        converted = converter(rows)
        # Row conservation, asserted rather than hoped for.
        if len(converted) != len(rows):
            raise SystemExit(
                f"{name}: read {len(rows)} rows and produced {len(converted)}. A conversion that loses a row"
                " and exits zero is the failure this check exists for"
            )
        target = SEED_DIR / out_name
        write(converted, target)
        # And again from the file that was actually written, in case the writer lost one.
        written = _read(target)
        if len(written) != len(rows):
            raise SystemExit(f"{name}: wrote {len(written)} rows to {target.name}, expected {len(rows)}")
        print(f"{name}: {len(rows)} rows in -> {target.name}, {len(written)} rows out")
    print(f"raw transcriptions copied verbatim to {AS_DELIVERED.relative_to(REPO)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
