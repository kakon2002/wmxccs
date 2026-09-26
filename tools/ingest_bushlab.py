"""Bush Lab CCS database -> flat conversion, and a refusal where the platform is unknown.

Run: .venv/Scripts/python tools/ingest_bushlab.py    (.venv/bin/python on macOS and Linux)

WHAT THIS DOES AND WHAT IT DELIBERATELY WILL NOT DO
---------------------------------------------------
It reads all eight sheets of the downloaded workbook, resolves every structural trap in the
file, and writes the result to data/seed/as_delivered/ as one flat table: one row per CCS
VALUE, carrying the sheet, the cited paper, the analyte, the charge, the gas, the NATIVE
unit, the NATIVE value and the converted value in square angstroms. That file is the audit
trail for the conversion: every number the platform would ever use can be traced back to a
cell, a unit and a factor.

IT DOES NOT WRITE A SEED FILE, and that is the point of this adapter rather than a gap in
it. A seed row needs an `ims_type`, and THE WORKBOOK STATES NO PLATFORM ANYWHERE - not on
any sheet, not in any cell, not in any header. It states the gas, the charge, the cross
section and a per-row `Ref` naming one of nine papers. The platform is in those papers.

So this script refuses to emit a seed file until it is handed a per-paper platform map, and
the refusal is the feature: filling `ims_type` from the database page's sentence about the
collection being "primarily traveling-wave" would be inventing a per-row fact out of a
general one, and `primarily` is not a value any row can carry. sources.BUSH_LAB_CCS records
the same conclusion and the reason for it.

Pass a map and it writes:

    python tools/ingest_bushlab.py --platform-map bushlab_platforms.json

where the map names, per cited paper, the platform, the DTIMS method where the platform is
DTIMS, and the page and citation it was read from - so that the resolution is evidence and
not a preference.

THE FIVE TRAPS, ALL VERIFIED AGAINST THE FILE BEFORE THIS WAS WRITTEN
--------------------------------------------------------------------
1. UNITS DIFFER BY A FACTOR OF 100. Four sheets report nm^2 and four report A^2. The unit is
   read from the column header, per column, and never inferred from the magnitude - a
   protein at 19.1 nm^2 and a small molecule at 124.5 A^2 are both four-digit-ish numbers
   and the magnitudes overlap once charge states spread out. The native unit and the native
   value are both carried through so the conversion can be checked rather than trusted.
2. TWO IDENTICALLY LABELLED HELIUM COLUMNS. On 'Native-Like Protein Cations and', columns 5
   and 6 are both headed 'O(He) / nm^2' and are distinguished ONLY by a merged banner one
   row above reading 'Cations' and 'Anions'. Reading them as helium-against-nitrogen would
   silently turn an anion measurement into a nitrogen measurement of a cation. This adapter
   reads the banner, and if the banner is not exactly the pair it expects it REFUSES the
   sheet rather than guessing which column is which.
3. CHARGE WITHOUT AN ADDUCT. Every protein, peptide and polymer sheet gives a charge state
   and never says what the charges are. That is the unstated-carrier case the schema already
   carries as `[M+n?]n+`, which keys uniquely per record, so two such records never match
   each other. Recorded, not guessed.
4. TWO-ROW HEADERS on three sheets, where a merged banner sits above the column names. The
   header depth is declared per sheet in LAYOUT and asserted against the file on every run.
5. MICROSOURCE IS THE ONLY SHEET WITH AN ADDUCT COLUMN, is nitrogen only, and is the only
   sheet carrying a per-row standard deviation. It needs none of the carrier handling and
   all of the unit handling.

A SIXTH, FOUND WHILE WRITING THIS
'Small Molecular Ions' HAS NO CHARGE COLUMN AT ALL. It gives a name, a class, an m/z and two
cross sections. A CCSMeasurement requires a charge, so those 23 rows cannot be built even
once the platform is known, and they are reported separately rather than folded in.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

WORKBOOK = REPO / "data" / "raw" / "bushlab" / "bushlab_ccs_database.xlsx"
CONVERSION = REPO / "data" / "seed" / "as_delivered" / "bushlab_ccs_database_converted.csv"
SEED = REPO / "data" / "seed" / "bushlab_microsource.csv"

SOURCE = "Bush Lab CCS database"
SOURCE_URL = "https://biophysicalms.org/ccsdatabase"

# 1 nm^2 = 100 A^2. Stated as a named constant because it is the whole of trap 1 and a bare
# 100 in an expression is the kind of thing that gets "simplified" away.
A2_PER_NM2 = 100.0
UNITS = {"nm^2": A2_PER_NM2, "A^2": 1.0}


class SheetRefused(Exception):
    """A sheet whose structure is not what this adapter verified. Never a warning."""


# Per sheet: header depth, analyte column, charge column (None where absent), Ref column, and
# the CCS columns as (column, gas, unit, charge_rule). `charge_rule` is an int where the
# column group fixes the charge, '+'/'-' where a banner fixes the SIGN of the sheet's charge
# column, or None where the charge column is taken as written.
LAYOUT: dict[str, tuple] = {
    "Native-Like Protein Cations": (1, 0, 3, 8, (
        (4, "He", "nm^2", None), (5, "N2", "nm^2", None))),
    "Native-Like Protein Cations and": (2, 0, 3, 7, (
        (4, "He", "nm^2", "+"), (5, "He", "nm^2", "-"))),
    "Denatured Protein Cations": (1, 0, 2, 6, (
        (3, "He", "nm^2", None), (4, "N2", "nm^2", None))),
    "Polyalaine Cations": (2, 0, None, 7, (
        (1, "He", "A^2", 1), (2, "N2", "A^2", 1),
        (3, "He", "A^2", 2), (4, "N2", "A^2", 2),
        (5, "He", "A^2", 3), (6, "N2", "A^2", 3))),
    "Anionic Homopolymers": (2, 0, None, 5, (
        (1, "He", "A^2", -1), (2, "N2", "A^2", -1),
        (3, "He", "A^2", -2), (4, "N2", "A^2", -2))),
    "Other Peptides": (1, 0, 2, 6, (
        (3, "He", "nm^2", None), (4, "N2", "nm^2", None))),
    "Small Molecular Ions": (1, 0, None, 5, (
        (3, "He", "A^2", None), (4, "N2", "A^2", None))),
    "MicroSource Collection": (1, 0, 7, 11, (
        (8, "N2", "A^2", None),)),
}

# The banner that makes trap 2 readable. If this is not what the file says, the sheet is
# refused: two columns with one label and no banner cannot be told apart at all.
BANNER_SHEET = "Native-Like Protein Cations and"
BANNER_EXPECTED = {4: "Cations", 5: "Anions"}

# Which sheets carry no adduct, so every row takes an unstated carrier.
NO_ADDUCT_SHEETS = frozenset(LAYOUT) - {"MicroSource Collection"}

FIELDS = (
    "sheet", "cited_paper", "analyte", "analyte_kind", "charge", "polarity", "adduct",
    "drift_gas", "native_unit", "native_value", "ccs_a2", "conversion_factor",
    "uncertainty_native", "uncertainty_a2", "uncertainty_type", "formula", "cas",
    "bioactivity", "mass_kda", "n_subunits", "origin", "supplier", "source", "source_url",
    "source_locator", "ims_type", "dtims_method", "platform_evidence", "blocked_because",
)


def rows_of(sheet, header_rows: int) -> tuple[list, list]:
    """Return (header rows, data rows), asserting the declared header depth against the file."""
    rows = [r for r in sheet.iter_rows(values_only=True)
            if any(c not in (None, "") for c in r)]
    first_numeric = next(
        (i for i, r in enumerate(rows) if any(isinstance(c, (int, float)) for c in r)),
        len(rows),
    )
    if first_numeric != header_rows:
        raise SheetRefused(
            f"{sheet.title!r}: declared {header_rows} header row(s), the file has"
            f" {first_numeric}. The layout in this adapter no longer matches the workbook;"
            f" re-read the sheet rather than adjusting the number to make it pass."
        )
    return rows[:header_rows], rows[header_rows:]


def unit_from_header(header: list, column: int, declared: str) -> str:
    """Read the unit out of the column header, and refuse if it is not what was declared.

    This is trap 1's guard. The declared unit in LAYOUT is a claim about the file; this
    checks it against the file on every run, so a workbook that changes its units cannot be
    converted by a factor that is no longer right.
    """
    text = str(header[column]) if column < len(header) and header[column] else ""
    found = [name for name in UNITS if name.replace("^", "^") in text.replace("Å", "A")]
    if declared not in found:
        raise SheetRefused(
            f"column {column} header is {text!r}, which does not state the declared unit"
            f" {declared!r}. The unit must be read from the header, never inferred from the"
            f" magnitude of the values."
        )
    return declared


def check_banner(header_rows: list) -> None:
    """Trap 2. Two columns, one label, told apart only by the row above."""
    banner = header_rows[0] if header_rows else []
    seen = {
        i: str(banner[i]).strip()
        for i in BANNER_EXPECTED
        if i < len(banner) and banner[i] not in (None, "")
    }
    if seen != BANNER_EXPECTED:
        raise SheetRefused(
            f"{BANNER_SHEET!r}: columns 4 and 5 are both headed as helium and are"
            f" distinguished only by the banner above them, which should read"
            f" {BANNER_EXPECTED} and reads {seen}. Without that banner the two columns"
            f" cannot be told apart, and reading them as helium against nitrogen would"
            f" turn an anion into a nitrogen measurement of a cation. Refusing the sheet."
        )


def unstated_carrier_adduct(charge: int) -> str:
    """`[M+3?]3+` - a well-formed ion whose charge carrier the source never named.

    Keys uniquely per record, so two such records never match each other. That is the
    schema's answer to trap 3 and it was built for exactly this source.
    """
    sign = "+" if charge > 0 else "-"
    return f"[M{sign}{abs(charge)}?]{abs(charge)}{sign}"


def adduct_problem(adduct: str) -> str:
    """Why the schema would refuse this adduct, or "" if it accepts it.

    Checked at conversion time rather than left to fail at load time, so the count is
    reported beside the other blockers instead of surfacing as a validation error on a row
    somebody has to go and find.
    """
    from wmxccs.identity import parse_adduct

    try:
        parse_adduct(adduct)
    except ValueError as refusal:
        return str(refusal)
    return ""


def analyte_kind_for(sheet: str) -> str:
    if sheet in ("Native-Like Protein Cations", "Native-Like Protein Cations and",
                 "Denatured Protein Cations"):
        return "protein"
    if sheet == "Other Peptides":
        return "peptide"
    return "small_molecule"


def convert(workbook_path: pathlib.Path, platform_map: dict | None) -> tuple[list, dict]:
    import openpyxl

    book = openpyxl.load_workbook(workbook_path, read_only=True, data_only=True)
    out: list[dict] = []
    tally: collections.Counter = collections.Counter()
    refused_sheets: dict[str, str] = {}

    for name, (header_rows, a_col, z_col, ref_col, ccs_cols) in LAYOUT.items():
        sheet = book[name]
        try:
            headers, data = rows_of(sheet, header_rows)
            column_header = headers[-1] if headers else []
            for column, _gas, declared, _rule in ccs_cols:
                unit_from_header(column_header, column, declared)
            if name == BANNER_SHEET:
                check_banner(headers)
        except SheetRefused as refusal:
            refused_sheets[name] = str(refusal)
            tally["sheets refused"] += 1
            continue

        tally["data rows read"] += len(data)
        for line, row in enumerate(data, start=header_rows + 1):
            analyte = str(row[a_col]).strip() if a_col < len(row) and row[a_col] else ""
            if not analyte:
                tally["rows with no analyte"] += 1
                continue
            paper = str(row[ref_col]).strip() if ref_col < len(row) and row[ref_col] else ""

            for column, gas, unit, rule in ccs_cols:
                if column >= len(row) or not isinstance(row[column], (int, float)):
                    continue
                native = float(row[column])
                tally["ccs values read"] += 1

                # --- charge, per trap 3 and the sheet's own shape ------------------------
                blocked: list[str] = []
                if isinstance(rule, int):
                    charge = rule
                elif rule in ("+", "-"):
                    raw = row[z_col] if z_col is not None and z_col < len(row) else None
                    charge = int(raw) if isinstance(raw, (int, float)) else None
                    if charge is not None and rule == "-":
                        charge = -charge
                elif z_col is not None and z_col < len(row) and isinstance(row[z_col], (int, float)):
                    charge = int(row[z_col])
                else:
                    charge = None
                if charge is None:
                    blocked.append(
                        "no charge stated: this sheet has no charge column, and a"
                        " measurement without a charge is not a measurement of an ion"
                    )
                    tally["values with no charge"] += 1

                # --- adduct: stated on MicroSource, an unstated carrier everywhere else ---
                if name in NO_ADDUCT_SHEETS:
                    adduct = unstated_carrier_adduct(charge) if charge is not None else ""
                    if charge is not None:
                        tally["values with an unstated charge carrier"] += 1
                else:
                    adduct = str(row[4]).strip() if row[4] else ""
                    if not adduct:
                        blocked.append("no adduct stated on a sheet that has an adduct column")
                    else:
                        # THE ADDUCT IS VALIDATED HERE AND NEVER REPAIRED. Three MicroSource
                        # values carry strings this schema refuses: 'M+' twice, with no
                        # brackets, and '[M+3H]+3' once, with the charge suffix reversed.
                        # '[M+3H]+3' is almost certainly '[M+3H]3+' and that is exactly why
                        # it is not corrected here - "almost certainly" is a guess about
                        # somebody else's data, and a silently repaired adduct is a fact
                        # invented in a column that decides which ions match.
                        problem = adduct_problem(adduct)
                        if problem:
                            blocked.append(f"adduct not in a form the schema accepts: {problem}")
                            tally["values with an unusable adduct"] += 1

                factor = UNITS[unit]
                spread = row[9] if name == "MicroSource Collection" and len(row) > 9 else None
                spread = float(spread) if isinstance(spread, (int, float)) else None

                # --- the platform, which is the whole blocker --------------------------
                ims_type = dtims_method = evidence = ""
                if platform_map and paper in platform_map:
                    entry = platform_map[paper]
                    ims_type = entry.get("ims_type", "")
                    dtims_method = entry.get("dtims_method", "")
                    evidence = entry.get("read_from", "")
                else:
                    blocked.append(
                        f"no platform stated: the workbook names no instrument or method"
                        f" anywhere, and the paper it cites ({paper or 'unnamed'}) has not"
                        f" been read. A seed row needs an ims_type and it cannot be"
                        f" inferred from this file"
                    )
                    tally["values with no platform"] += 1

                out.append({
                    "sheet": name,
                    "cited_paper": paper,
                    "analyte": analyte,
                    "analyte_kind": analyte_kind_for(name),
                    "charge": "" if charge is None else charge,
                    "polarity": "" if charge is None else ("positive" if charge > 0 else "negative"),
                    "adduct": adduct,
                    "drift_gas": gas,
                    "native_unit": unit,
                    "native_value": repr(native),
                    "ccs_a2": repr(native * factor),
                    "conversion_factor": repr(factor),
                    "uncertainty_native": "" if spread is None else repr(spread),
                    "uncertainty_a2": "" if spread is None else repr(spread * factor),
                    "uncertainty_type": "sd" if spread is not None else "unknown",
                    "formula": str(row[1]).strip() if name == "MicroSource Collection" and row[1] else "",
                    "cas": str(row[2]).strip() if name == "MicroSource Collection" and row[2] else "",
                    "bioactivity": str(row[3]).strip() if name == "MicroSource Collection" and row[3] else "",
                    "mass_kda": repr(row[2]) if name != "MicroSource Collection" and isinstance(row[2] if 2 < len(row) else None, (int, float)) else "",
                    "n_subunits": repr(row[1]) if name.startswith("Native-Like") and isinstance(row[1], (int, float)) else "",
                    "origin": str(row[6]).strip() if name.startswith("Native-Like") and 6 < len(row) and row[6] else "",
                    "supplier": "",
                    "source": SOURCE,
                    "source_url": SOURCE_URL,
                    "source_locator": f"sheet {name!r}, line {line}, column {column}",
                    "ims_type": ims_type,
                    "dtims_method": dtims_method,
                    "platform_evidence": evidence,
                    "blocked_because": " | ".join(blocked),
                })
    book.close()
    return out, {"tally": tally, "refused_sheets": refused_sheets}


SEED_HEADER_FROM = REPO / "data" / "seed" / "steroid_jasms2022.csv"


def seed_fieldnames() -> list[str]:
    """The seed schema, taken from an existing seed file rather than retyped.

    Sixty-three columns, and a hand-copied list of them is a list that will drift. Reading it
    from a committed seed file means a schema change surfaces here as a KeyError on the next
    run rather than as a silently empty column.
    """
    with SEED_HEADER_FROM.open(encoding="utf-8-sig", newline="") as handle:
        return next(csv.reader(handle))


def conformer_assignments(rows: list[dict], named: set[str]) -> dict[str, tuple[int, int, str]]:
    """Which converted rows are conformers of one ion, by locator.

    THE PAPER LUMPS THREE SITUATIONS AND THE SHEET DOES NOT SEPARATE THEM. Hines 2017 says 16
    drugs "display two peaks, had two major adducts, or were mixtures for which we have
    reported a CCS value of each component". Four of the sixteen differ by ADDUCT and are
    ordinary separate ions - nothing to do here. The other twelve share an adduct, so their two
    values are either two peaks of one ion or two components of a mixture.

    A PROTOMER IS A CONFORMER. A MIXTURE COMPONENT IS A DIFFERENT ANALYTE. Calling the second
    a conformer would merge two compounds into one ion, so the twelve are not treated alike:
    the six the paper NAMES - the fluoroquinolone protomers and its new cephalosporin finding -
    are conformers on the paper's own authority, and the other six are recorded as conformers
    AND carry a curation flag saying the sheet does not distinguish the two cases. Neither is
    averaged and neither is dropped, which is what the instruction required; what is not done
    is asserting protomerism for six compounds the paper never called protomers.
    """
    by_compound: dict[tuple[str, str], list[dict]] = {}
    for row in rows:
        by_compound.setdefault((row["analyte"].casefold(), row["adduct"]), []).append(row)

    out: dict[str, tuple[int, int, str]] = {}
    for (compound, _adduct), group in by_compound.items():
        if len(group) < 2:
            continue
        # stable order: by the cross section itself, so the numbering is reproducible and
        # does not depend on how the sheet happened to be read
        ordered = sorted(group, key=lambda r: float(r["ccs_a2"]))
        flag = "" if compound in named else (
            "MAY BE A MIXTURE RATHER THAN TWO CONFORMERS, and the sheet does not say which."
            " Hines 2017 reports that sixteen compounds 'display two peaks, had two major"
            " adducts, or were mixtures for which we have reported a CCS value of each"
            " component', and it NAMES only the fluoroquinolone protomers and its new"
            " cephalosporin finding. This compound is not among those, so its two values at"
            " one adduct are recorded as two conformers on the evidence available and may"
            " instead be two components of a mixture - which would make them two ANALYTES,"
            " not two conformers of one ion. The clearest case is antimycin A, which this"
            " sheet itself labels 'antimycin a (a1 shown)' and which is a mixture of"
            " antimycins A1 to A4. Resolving this means reading the paper's figures; it is"
            " deliberately not resolved, because these records can form no matched ion and"
            " the stakes are a reference-library entry rather than a correction."
        )
        for index, row in enumerate(ordered, start=1):
            out[row["source_locator"]] = (index, len(ordered), flag)
    return out


def seed_rows(rows: list[dict], platform_map: dict) -> tuple[list[dict], dict]:
    """Converted rows -> seed rows, for the papers the map resolves. Blocked rows are skipped."""
    fields = seed_fieldnames()
    emitted: list[dict] = []
    tally: collections.Counter = collections.Counter()

    eligible = [r for r in rows if r["cited_paper"] in platform_map and not r["blocked_because"]]
    named_per_paper = {
        paper: {n.casefold() for n in entry.get("protomer_pairs_named_by_the_paper", ())}
        for paper, entry in platform_map.items()
        if isinstance(entry, dict)
    }
    conformers: dict[str, tuple[int, int, str]] = {}
    for paper in {r["cited_paper"] for r in eligible}:
        conformers.update(
            conformer_assignments(
                [r for r in eligible if r["cited_paper"] == paper],
                named_per_paper.get(paper, set()),
            )
        )

    for r in eligible:
        entry = platform_map[r["cited_paper"]]
        row = {name: "" for name in fields}
        index, total, flag = conformers.get(r["source_locator"], (None, None, ""))
        row.update({
            "analyte_kind": r["analyte_kind"],
            "analyte_display_name": r["analyte"],
            "analyte_source": SOURCE,
            "analyte_reuse_status": "academic_only",
            # DATASET-SCOPED, so it pairs inside this source and matches nothing outside it.
            # The sheet also gives a CAS number, which would be a real cross-source identifier -
            # the seed schema has no column for one, so it stays in the as_delivered conversion
            # and this id deliberately bridges nothing.
            "analyte_dataset_compound_id": f"bushlab_microsource:{r['analyte']}",
            "ccs": r["ccs_a2"],
            "ccs_uncertainty": r["uncertainty_a2"],
            "uncertainty_type": r["uncertainty_type"],
            "conformer": "" if index is None else index,
            "conformers_total": "" if total is None else total,
            "adduct": r["adduct"],
            "charge": r["charge"],
            "polarity": r["polarity"],
            "ims_type": entry["ims_type"],
            "dtims_method": entry.get("dtims_method", ""),
            "drift_gas": entry.get("drift_gas", r["drift_gas"]),
            "cell_gas": entry.get("cell_gas", ""),
            "calibrant": entry.get("calibrant", ""),
            "instrument": entry.get("instrument", ""),
            "source": SOURCE,
            # `doi` IS LEFT BLANK ON PURPOSE, and the licence gate is why. Permission to use
            # these values comes from the Bush Lab database page, whose terms were read and
            # recorded; it does not come from the ACS paper, whose licence nobody has read -
            # Europe PMC flags it open access and states no licence, and an open-access flag
            # is not a licence. Naming the paper's DOI here would assert `academic_only` for a
            # publication with no licence record, and the gate refuses that, correctly. This
            # is the decision already taken for CCSbase and recorded in LIMITATIONS 7E: the
            # source is where the value was obtained and whose terms permit the use, with the
            # primary paper named in the locator so the real origin stays visible.
            "doi": "",
            "source_locator": f"{r['source_locator']}; published in doi:{entry.get('doi', '')}",
            "replicates": entry.get("replicates", ""),
            "reuse_status": "academic_only",
            # CONDITION 1: the lineage that makes these records dangerous to compare.
            "calref_reference_set": entry.get("calref_reference_set", ""),
            "calref_doi": entry.get("calref_doi", ""),
            "calref_platform": entry.get("calref_platform", ""),
            "calref_method": entry.get("calref_method", ""),
            "curation_flag": flag,
        })
        emitted.append(row)
        tally[entry["ims_type"]] += 1
        if index is not None:
            tally["conformer rows"] += 1
        if flag:
            tally["conformer rows carrying a curation flag"] += 1
    return emitted, tally


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python tools/ingest_bushlab.py",
        description="Convert the Bush Lab CCS workbook; refuse to seed it without a platform.",
    )
    parser.add_argument(
        "--platform-map",
        type=pathlib.Path,
        help="JSON: {\"<cited paper>\": {\"ims_type\": ..., \"dtims_method\": ..., "
        "\"read_from\": \"page and citation\"}}. Without it, no seed file is written.",
    )
    args = parser.parse_args(argv)

    if not WORKBOOK.is_file():
        print(f"workbook not found: {WORKBOOK}", file=sys.stderr)
        return 2

    platform_map = None
    if args.platform_map:
        platform_map = json.loads(args.platform_map.read_text(encoding="utf-8"))

    try:
        rows, report = convert(WORKBOOK, platform_map)
    except SheetRefused as refusal:
        print(f"REFUSED: {refusal}", file=sys.stderr)
        return 1

    CONVERSION.parent.mkdir(parents=True, exist_ok=True)
    with CONVERSION.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    tally = report["tally"]
    print(f"workbook   {WORKBOOK.relative_to(REPO)}")
    print(f"conversion {CONVERSION.relative_to(REPO)}   ({len(rows)} rows written)")
    print()
    for key in ("data rows read", "ccs values read", "values with an unstated charge carrier",
                "values with no charge", "values with an unusable adduct",
                "values with no platform", "rows with no analyte", "sheets refused"):
        print(f"   {key:42} {tally[key]}")
    for name, why in report["refused_sheets"].items():
        print(f"   REFUSED SHEET {name!r}: {why}")

    by_sheet = collections.Counter(r["sheet"] for r in rows)
    by_gas = collections.Counter(r["drift_gas"] for r in rows)
    by_unit = collections.Counter(r["native_unit"] for r in rows)
    print()
    print("   values by sheet:")
    for name, count in by_sheet.items():
        print(f"      {name:38} {count:5}")
    print(f"   by gas  {dict(by_gas)}")
    print(f"   by unit {dict(by_unit)}   (nm^2 converted at x{A2_PER_NM2:g})")

    blocked = [r for r in rows if r["blocked_because"]]

    # --- the seed, for whatever the map resolves -------------------------------------------
    if platform_map:
        seeded, seed_tally = seed_rows(rows, platform_map)
        if seeded:
            with SEED.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=seed_fieldnames())
                writer.writeheader()
                writer.writerows(seeded)
            print()
            print(f"   SEED WRITTEN: {SEED.relative_to(REPO)}   ({len(seeded)} records)")
            for key, count in sorted(seed_tally.items()):
                print(f"      {key:44} {count}")
            resolved = sorted(k for k in platform_map if not k.startswith("_"))
            print(f"      resolved papers: {resolved}")

    print()
    if blocked:
        print(f"   {len(blocked)} of {len(rows)} values remain blocked.")
        reasons = collections.Counter(
            part.split(":")[0] for r in blocked for part in r["blocked_because"].split(" | ")
        )
        for reason, count in reasons.most_common():
            print(f"      {count:5}  {reason}")
        papers = sorted({r["cited_paper"] for r in rows if r["cited_paper"]})
        print()
        print(f"   The platform is answerable per paper. {len(papers)} papers are cited:")
        for paper in papers:
            n = sum(1 for r in rows if r["cited_paper"] == paper)
            print(f"      {paper:20} {n:5} value(s)")
        print()
        print("   Supply --platform-map with an ims_type per paper and the page it was read")
        print("   from, and this writes their rows to the seed file. It will not guess one.")
        return 0

    print(f"   {len(rows) - len(blocked)} value(s) are not blocked and can be seeded.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
