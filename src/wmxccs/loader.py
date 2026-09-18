"""CCS measurements, read from a file of one row per measured value.

This is the path every measurement enters by, whether a table transcribed from a
paper or an export from an instrument. One format, one loader, one set of checks,
so a value cannot arrive by a route that skips the licence gate or the
validators. A source whose own columns differ is converted to this format first,
and the conversion is a separate, auditable step that writes a file somebody can
read - not a second entrance.

The format is a CSV with one row per measurement and one column per field, with
the analyte's fields flattened in under an `analyte_` prefix and an
`analyte_kind` column deciding which variant of the union a row builds. Every row
carries its own provenance - DOI, locator, reuse status - rather than inheriting
it from a file header, because a file is a convenience and a row is the record. A
row is one line: no cell of this format may contain a line break, and a file
whose rows cannot be told apart is refused whole rather than read as best it can
be. A transcription with a stray quotation mark in one cell would otherwise fold
the next row into that cell, and the report would count one row fewer with
nothing to say that it had.

Nothing is guessed and nothing is dropped quietly. A row that does not validate
is counted under the reason it failed, with an example naming the line in the
file. A row that validates but may not train is KEPT, and the reason it may not
train is counted separately: storing a record is not the same as training on it,
and the two counts are the two things a reader needs.

A blank cell means "not reported" and becomes null; a filler such as "n/a" is a
data-quality fault and is refused by validation rather than mapped to null here,
so it is seen rather than absorbed. Numbers and dates are read in one spelling
each - plain ASCII digits, and YYYY-MM-DD - so that a cell can only become the
value it visibly is.

Both of a row's reuse-status claims - the measurement's and the analyte's - are
checked against the registry in sources.py by the licence gate itself, so the
check is the same for a record from a file and for one built in code.
"open_attribution" typed into a cell is not a licence; a licence record is.

One status may never come from a file at all. synthetic_fixture is a test's
declaration that a record was built in code; a row in a file is a real record, so
a row carrying that status in either reuse-status column is refused before any
record exists, and counted under its own reason. That is what keeps the status
inside the tests.

Two more things keep a record out of training without dropping it. A row may
carry a `curation_flag` naming a reason a person wants it reviewed. And the
loader itself flags a SUSPECTED SHARED PEAK: two different analytes in one
calibration group reporting an identical CCS may be one unresolved peak
transcribed against both compounds, and training on it would teach the model that
two molecules share a value. Such rows are held back until a person confirms they
are separate measurements. The check fails closed: a row that shares a value with
a row that did not validate is held too, because the partner's conditions could
not be compared. Holding back is the safe direction: it costs data, not honesty.
"""

from __future__ import annotations

import csv
import io
import math
import re
from collections import Counter
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Mapping

from .identity import ANALYTE_TYPES, AnalyteKind
from .licensing import LicenceGateError, ReuseStatus, TrainingGateError, UnbackedClaimError, assert_trainable
from .models import CCSMeasurement, CalibrationReference
from .readiness import MatchedIonSet, Readiness

# The analyte's fields, flattened under a prefix. `analyte_source` and
# `analyte_reuse_status` are prefixed like the rest because the measurement
# carries a `source` and `reuse_status` of its own: where the ANALYTE IDENTITY
# was assigned is not where the VALUE came from.
ANALYTE_PREFIX = "analyte_"
ANALYTE_COLUMNS = (
    "analyte_kind",
    "analyte_display_name",
    "analyte_source",
    "analyte_reuse_status",
    # small molecule
    "analyte_inchikey",
    "analyte_smiles",
    # peptide, protein, antibody
    "analyte_sequence",
    "analyte_modifications",
    "analyte_accession",
    "analyte_subunit",
    # glycan
    "analyte_composition",
    "analyte_wurcs",
    "analyte_glytoucan_ac",
    "analyte_iupac_condensed",
    "analyte_has_unresolved_linkage",
    "analyte_has_unresolved_anomericity",
    "analyte_reducing_end_label",
    "analyte_derivatisation",
    # protein, antibody, ADC
    "analyte_folding_state",
    # antibody and ADC
    "analyte_inn",
    "analyte_glycoform",
    "analyte_ciu_state",
    # ADC
    "analyte_linker_payload_class",
    "analyte_dar",
    "analyte_conjugation_state",
)

# Which flattened columns belong to the nested AntibodyIdentity rather than to
# the analyte itself. Nested because one antibody identity is shared by the
# intact-antibody and ADC variants, and duplicating its validation would let the
# two drift apart.
ANTIBODY_COLUMNS = ("analyte_inn", "analyte_accession", "analyte_sequence")

CALIBRATION_REFERENCE_COLUMNS = (
    "calref_reference_set",
    "calref_doi",
    "calref_platform",
    "calref_method",
)

MEASUREMENT_COLUMNS = (
    "ccs",
    "ccs_uncertainty",
    "uncertainty_type",
    "conformer",
    "conformers_total",
    "adduct",
    "charge",
    "polarity",
    "ims_type",
    "dtims_method",
    "drift_gas",
    "calibrant",
    "cell_gas",
    "instrument",
    "source",
    "doi",
    "source_locator",
    "replicates",
    "measured_on",
    "reuse_status",
)

# Read by the loader and never passed to a record: a reason a person wants the
# row reviewed before it may train. Free text; any non-blank value holds the row.
LOADER_COLUMNS = ("curation_flag",)
COLUMNS = ANALYTE_COLUMNS + CALIBRATION_REFERENCE_COLUMNS + MEASUREMENT_COLUMNS + LOADER_COLUMNS

# Without these a row cannot even be attempted. Everything else is optional and
# falls to the model's own defaults, which for the licence means "unverified",
# which blocks training.
REQUIRED_COLUMNS = (
    "analyte_kind",
    "analyte_source",
    "ccs",
    "adduct",
    "charge",
    "polarity",
    "ims_type",
    "drift_gas",
    "source",
)

_BOOLEAN_COLUMNS = frozenset({"analyte_has_unresolved_linkage", "analyte_has_unresolved_anomericity"})
_INTEGER_COLUMNS = frozenset({"charge", "replicates", "conformer", "conformers_total", "analyte_dar"})
_FLOAT_COLUMNS = frozenset({"ccs", "ccs_uncertainty"})
_DATE_COLUMNS = frozenset({"measured_on"})
_STATUS_COLUMNS = frozenset({"reuse_status", "analyte_reuse_status"})
# Semicolon-separated, because a comma is the field separator and a modification
# name may legitimately contain neither.
_LIST_COLUMNS = frozenset({"analyte_modifications"})
_BOOLEANS = {"true": True, "false": False, "1": True, "0": False, "yes": True, "no": False}
# One spelling each, ASCII digits only. Python's int() and float() also take
# digit-group underscores and every Unicode digit, and float() takes "nan" and
# "inf"; none of those is a number a transcriber wrote on purpose.
_INTEGER = re.compile(r"[+-]?[0-9]+")
_FLOAT = re.compile(r"[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")
_ISO_DATE = re.compile(r"[0-9]{4}-[0-9]{2}-[0-9]{2}")

# Why a row did not become a record.
NOT_A_NUMBER = "a numeric column does not hold a number"
NOT_A_BOOLEAN = "a boolean column does not hold true or false"
NOT_A_DATE = "a date column does not hold a date written YYYY-MM-DD"
WRONG_CELL_COUNT = "the row does not have one cell per header column"
UNKNOWN_ANALYTE_KIND = "analyte_kind is not one of the analyte kinds"
SYNTHETIC_STATUS_IN_FILE = (
    "a file may not carry synthetic_fixture: that status is a test's declaration of a record built in code,"
    " and a row in a file is a real record"
)
ANALYTE_REJECTED = "the analyte failed validation"
MEASUREMENT_REJECTED = "the measurement failed validation"

# Why a built record may not train.
GATE_UNBACKED_CLAIM = "reuse-status claim not backed by a licence record"
GATE_LICENCE = "licence not cleared"
GATE_ANALYTE = "not defined well enough to train: the analyte, the ion, the reported uncertainty, or the record itself"
GATE_CURATION = "held for curation review"
GATE_SHARED_PEAK = "suspected shared peak: identical CCS for different analytes in one calibration group"
GATE_LONE_CONFORMER = "a conformer set is incomplete: a row claims siblings that are not in the file"


class _Coercion(ValueError):
    def __init__(self, reason: str, detail: str) -> None:
        super().__init__(detail)
        self.reason = reason


class _Union:
    """Union-find over identity atoms. Ported; the 'smallest root wins' rule matters.

    It makes a component id stable and order-independent, which is what lets the
    shared-peak check give the same answer however the rows are ordered in the
    file.
    """

    def __init__(self) -> None:
        self._parent: dict[str, str] = {}

    def find(self, node: str) -> str:
        parent = self._parent.setdefault(node, node)
        while parent != node:
            node, parent = parent, self._parent.setdefault(parent, parent)
        return node

    def union(self, a: str, b: str) -> None:
        root_a, root_b = self.find(a), self.find(b)
        if root_a != root_b:
            low, high = sorted((root_a, root_b))
            self._parent[high] = low


# --- reading a table -----------------------------------------------------------------


@dataclass(frozen=True)
class Row:
    """One data row of a CSV table, by the physical line it sits on."""

    line: int  # 1-based line in the file; the header is line 1
    cells: Mapping[str, str]  # column -> stripped text, blank where the cell is blank
    fault: str | None = None  # why the row cannot be read at all, if it cannot

    @property
    def where(self) -> str:
        return f"line {self.line}"


def read_table(text: str, *, label: str) -> tuple[tuple[str, ...], tuple[Row, ...]]:
    """The header and every data row of a CSV, or a ValueError naming what is wrong with the file.

    Raised, not counted, because each of these makes the row boundaries or the
    column meanings uncertain for the whole file, and a partial report would look
    like a clean run over fewer rows:

    - no header, an unnamed header cell, or a column name that appears twice
      (csv.DictReader would keep only the last cell under a repeated name);
    - a quoting fault, such as an unclosed quotation mark, which the csv module
      would otherwise repair by reading the following rows into the open cell;
    - a row that spans more than one physical line, which is the one form of that
      fault the csv module cannot see (two stray quotation marks fencing a whole
      row) and which no cell of this format may legitimately produce.

    A row with the wrong number of cells is returned with its `fault` set, for the
    caller to count: that fault stays inside its own line.
    """
    reader = csv.reader(io.StringIO(text.lstrip("﻿")), strict=True)
    try:
        header = next(reader, None)
        if not header or not any(name.strip() for name in header):
            raise ValueError(f"{label}: the file has no header row")
        if reader.line_num != 1:
            raise ValueError(f"{label}: the header row runs on to line {reader.line_num}; a cell contains a line break")
        fields = tuple(name.strip() for name in header)
        if any(not name for name in fields):
            raise ValueError(f"{label}: the header has an unnamed column (cell {fields.index('') + 1})")
        repeated = sorted(name for name, count in Counter(fields).items() if count > 1)
        if repeated:
            raise ValueError(f"{label}: the header names a column more than once: {repeated}")

        rows: list[Row] = []
        while True:
            start = reader.line_num + 1
            cells = next(reader, None)
            if cells is None:
                break
            if not cells:
                continue  # a blank line holds nothing, so skipping it drops nothing
            end = reader.line_num
            if end != start:
                raise ValueError(
                    f"{label}: the row at line {start} runs on to line {end}; a cell contains a line break,"
                    " which no cell of this format may (an unclosed quotation mark?)"
                )
            if len(cells) != len(fields):
                fault = f"{len(cells)} cells against {len(fields)} header columns"
                rows.append(Row(line=start, cells={}, fault=fault))
                continue
            rows.append(Row(line=start, cells={name: cell.strip() for name, cell in zip(fields, cells)}))
    except csv.Error as exc:
        raise ValueError(f"{label}: malformed CSV at line {reader.line_num}: {exc}") from exc
    return fields, tuple(rows)


# --- the load ------------------------------------------------------------------------


@dataclass(frozen=True)
class MeasurementLoadReport:
    """What a file produced, what it could not, and how far it gets us.

    Two conservation laws hold, and tests pin both:
    records_built + rows_failed == rows_read, and
    records_cleared + records_refused == records_built.
    """

    label: str
    records: tuple[CCSMeasurement, ...]  # every record that validated, cleared or not
    cleared: tuple[CCSMeasurement, ...]  # the subset that may train
    rows_read: int
    failure_counts: Mapping[str, int] = field(default_factory=dict)
    failure_examples: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    gate_counts: Mapping[str, int] = field(default_factory=dict)
    gate_examples: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    unknown_columns: tuple[str, ...] = ()
    readiness: Readiness | None = None

    @property
    def records_built(self) -> int:
        return len(self.records)

    @property
    def records_cleared(self) -> int:
        return len(self.cleared)

    @property
    def rows_failed(self) -> int:
        return sum(self.failure_counts.values())

    @property
    def records_refused(self) -> int:
        return sum(self.gate_counts.values())

    # These four describe what may TRAIN. When rows are held they describe less
    # and less of the file, and a load where nothing clears reports zero of
    # everything while holding a table full of it. The `_held` four below
    # describe what the file HOLDS, cleared or not, so that a report of a fully
    # held file still says what is in it.

    @property
    def matched_ion_keys(self) -> int:
        return len({record.matched_ion_key for record in self.cleared})

    @property
    def calibration_groups(self) -> int:
        return len({record.calibration_group for record in self.cleared})

    @property
    def adducts(self) -> Mapping[str, int]:
        return dict(Counter(record.adduct for record in self.cleared))

    @property
    def gas_split(self) -> Mapping[str, int]:
        """Cleared records by (drift gas the value refers to, gas in the cell).

        The two are different questions. A travelling-wave value calibrated
        against helium reference values is a helium CCS whatever gas was in the
        cell, and it must pool with helium values, not with nitrogen ones.
        """
        return dict(Counter(_gas_pair(record) for record in self.cleared))

    @property
    def matched_ion_keys_held(self) -> int:
        """Distinct matched ions the file holds, EXCLUDING ones that can never pair."""
        return len({record.matched_ion_key for record in self.records if record.matched_ion_key.matchable})

    @property
    def unmatchable_held(self) -> int:
        """Records whose charge carrier the source never named.

        Counted apart, because each carries a key unique to itself and can never
        pair with anything - including another record of the same shape. Folded
        into matched_ion_keys_held they would read as ions waiting for a partner,
        and a native-MS protein file would report hundreds of distinct matched
        ions and no possible pair, which reads as a puzzle rather than as the
        one-line problem it is.
        """
        return sum(1 for record in self.records if not record.matched_ion_key.matchable)

    @property
    def calibration_groups_held(self) -> int:
        return len({record.calibration_group for record in self.records})

    @property
    def adducts_held(self) -> Mapping[str, int]:
        return dict(Counter(record.adduct for record in self.records))

    @property
    def gas_split_held(self) -> Mapping[str, int]:
        return dict(Counter(_gas_pair(record) for record in self.records))

    @property
    def platforms_held(self) -> Mapping[str, int]:
        return dict(Counter(_platform_of(record) for record in self.records))

    @property
    def analyte_kinds_held(self) -> Mapping[str, int]:
        return dict(Counter(str(record.analyte.kind) for record in self.records))

    @property
    def uncertainty_types(self) -> Mapping[str, int]:
        """Cleared records by what their reported spread is, with 'not reported' counted."""
        return _uncertainty_tally(self.cleared)

    @property
    def uncertainty_types_held(self) -> Mapping[str, int]:
        """The same tally over everything the file holds, cleared or not.

        The fifth of the `_held` figures, and it was missing. A file where nothing
        clears said nothing at all about its spreads - and for the 2015 seed an
        uncertainty type of 'unknown' is one of the two things blocking all 89
        rows, so the report was silent about half of its own headline.
        """
        return _uncertainty_tally(self.records)

    def summary(self) -> str:
        lines = [
            f"{self.label}",
            f"  rows read                 {self.rows_read}",
            f"  records built             {self.records_built}",
            f"  rows not loaded           {self.rows_failed}",
        ]
        for reason, count in sorted(self.failure_counts.items(), key=lambda item: -item[1]):
            lines.append(f"    {reason}: {count}")
            for example in self.failure_examples.get(reason, ())[:2]:
                lines.append(f"      e.g. {example}")
        lines += [
            f"  cleared the licence gate  {self.records_cleared}",
            f"  kept but may not train    {self.records_refused}",
        ]
        for reason, count in sorted(self.gate_counts.items(), key=lambda item: -item[1]):
            lines.append(f"    {reason}: {count}")
            for example in self.gate_examples.get(reason, ())[:2]:
                lines.append(f"      e.g. {example}")
        lines += [
            f"  distinct matched ions     {self.matched_ion_keys}  (of the cleared)",
            f"  calibration groups        {self.calibration_groups}",
        ]
        if self.adducts:
            lines.append("  adducts                   " + ", ".join(f"{a} x{n}" for a, n in sorted(self.adducts.items())))
        if self.gas_split:
            lines.append("  gas                       " + ", ".join(f"{g} x{n}" for g, n in sorted(self.gas_split.items())))
        if self.uncertainty_types:
            lines.append(
                "  uncertainty reported as   "
                + ", ".join(f"{k} x{n}" for k, n in sorted(self.uncertainty_types.items()))
            )
        if self.records_built != self.records_cleared:
            lines += [
                "  in the file, cleared or not:",
                f"    distinct matched ions   {self.matched_ion_keys_held}",
                *(
                    [
                        f"    UNMATCHABLE             {self.unmatchable_held}  (charge carrier not stated;"
                        " can never pair, not even with each other)"
                    ]
                    if self.unmatchable_held
                    else []
                ),
                f"    calibration groups      {self.calibration_groups_held}",
            ]
            if self.platforms_held:
                lines.append(
                    "    platforms               " + ", ".join(f"{p} x{n}" for p, n in sorted(self.platforms_held.items()))
                )
            if self.analyte_kinds_held:
                lines.append(
                    "    analyte kinds           "
                    + ", ".join(f"{k} x{n}" for k, n in sorted(self.analyte_kinds_held.items()))
                )
            if self.adducts_held:
                lines.append(
                    "    adducts                 " + ", ".join(f"{a} x{n}" for a, n in sorted(self.adducts_held.items()))
                )
            if self.gas_split_held:
                lines.append(
                    "    gas                     " + ", ".join(f"{g} x{n}" for g, n in sorted(self.gas_split_held.items()))
                )
            if self.uncertainty_types_held:
                lines.append(
                    "    uncertainty reported as "
                    + ", ".join(f"{k} x{n}" for k, n in sorted(self.uncertainty_types_held.items()))
                )
        if self.unknown_columns:
            lines.append(f"  columns not recognised    {', '.join(self.unknown_columns)}  (ignored, and reported)")
        if self.readiness is not None:
            lines.append("")
            lines.append(self.readiness.summary())
        return "\n".join(lines)


@dataclass(frozen=True)
class _Failed:
    """What a row that produced no record still said, for the shared-peak check."""

    where: str
    ccs: float | None
    atoms: frozenset[str]


def _uncertainty_tally(records) -> Mapping[str, int]:
    return dict(
        Counter(
            str(record.uncertainty_type) if record.uncertainty_type is not None else "not reported"
            for record in records
        )
    )


def _gas_pair(record: CCSMeasurement) -> str:
    cell = record.cell_gas if record.cell_gas is not None else "not reported"
    return f"drift {record.drift_gas} / cell {cell}"


def _platform_of(record: CCSMeasurement) -> str:
    return str(record.ims_type) if record.dtims_method is None else f"{record.ims_type}/{record.dtims_method}"


def load_measurements(text: str, *, label: str = "measurements") -> MeasurementLoadReport:
    """Read a CSV of measurements into gated records, accounting for every row.

    Raises ValueError for a file that cannot be read as a table (see read_table)
    or whose header lacks a required column: that is a caller mistake, not a data
    gap, and must not look like a clean run with nothing in it. A header with no
    rows is an empty load and is reported as one.
    """
    fields, rows = read_table(text, label=label)
    missing = [column for column in REQUIRED_COLUMNS if column not in fields]
    if missing:
        raise ValueError(f"{label}: the header is missing the required columns {missing}")
    unknown = tuple(name for name in fields if name not in COLUMNS)

    records: list[CCSMeasurement] = []
    built: list[tuple[str, CCSMeasurement, str | None]] = []  # (where, record, curation flag)
    failed: list[_Failed] = []
    failures: Counter[str] = Counter()
    failure_examples: dict[str, list[str]] = {}
    gate: Counter[str] = Counter()
    gate_examples: dict[str, list[str]] = {}

    for row in rows:
        where = row.where
        if row.fault is not None:
            _note(failures, failure_examples, WRONG_CELL_COUNT, f"{where}: {row.fault}")
            continue  # no cell can be trusted to be in its column, so nothing is read from it
        try:
            analyte_fields, calref_fields, measurement_fields, flag = _coerce(row.cells)
        except _Coercion as exc:
            _note(failures, failure_examples, exc.reason, f"{where}: {exc}")
            failed.append(_as_failed(where, row.cells))
            continue
        try:
            analyte = _build_analyte(analyte_fields)
        except _Coercion as exc:
            # Its own bucket. A row naming a kind the union does not hold is a
            # wrong COLUMN, and counting it as failed validation makes it read as
            # bad data, which is a different fix by a different person.
            _note(failures, failure_examples, exc.reason, f"{where}: {exc}")
            failed.append(_as_failed(where, row.cells))
            continue
        except Exception as exc:
            _note(failures, failure_examples, ANALYTE_REJECTED, f"{where}: {_first_line(exc)}")
            failed.append(_as_failed(where, row.cells))
            continue
        try:
            if calref_fields:
                measurement_fields["calibration_reference"] = CalibrationReference(**calref_fields)
            record = CCSMeasurement(analyte=analyte, **measurement_fields)
        except Exception as exc:
            _note(failures, failure_examples, MEASUREMENT_REJECTED, f"{where}: {_first_line(exc)}")
            failed.append(_as_failed(where, row.cells))
            continue
        records.append(record)
        built.append((where, record, flag))

    # A shared peak can only be seen once every row is in: the same CCS under the
    # same conditions, reported against two different analytes.
    shared = _suspected_shared_peaks(built, failed)
    lone = _lone_conformers(built)

    cleared: list[CCSMeasurement] = []
    for where, record, flag in built:
        # Kept either way. The gate decides what may TRAIN, not what may be
        # stored. Order: a person's hold first, then the loader's own suspicion,
        # then the licence claims against the registry, and only then the gate.
        if flag:
            _note(gate, gate_examples, GATE_CURATION, f"{where}: {flag}")
            continue
        if where in shared:
            _note(gate, gate_examples, GATE_SHARED_PEAK, f"{where}: {shared[where]}")
            continue
        if where in lone:
            _note(gate, gate_examples, GATE_LONE_CONFORMER, f"{where}: {lone[where]}")
            continue
        try:
            assert_trainable(record)
        except UnbackedClaimError as exc:
            _note(gate, gate_examples, GATE_UNBACKED_CLAIM, f"{where}: {_first_line(exc)}")
        except LicenceGateError as exc:
            _note(gate, gate_examples, GATE_LICENCE, f"{where}: {_first_line(exc)}")
        except TrainingGateError as exc:
            _note(gate, gate_examples, GATE_ANALYTE, f"{where}: {_first_line(exc)}")
        else:
            cleared.append(record)

    return MeasurementLoadReport(
        label=label,
        records=tuple(records),
        cleared=tuple(cleared),
        rows_read=len(rows),
        failure_counts=dict(failures),
        failure_examples={reason: tuple(found) for reason, found in failure_examples.items()},
        gate_counts=dict(gate),
        gate_examples={reason: tuple(found) for reason, found in gate_examples.items()},
        unknown_columns=unknown,
        # Always produced, at zero records too, so the thresholds are visible from the first row.
        readiness=MatchedIonSet(records=tuple(cleared)).readiness,
    )


def load_measurements_file(path: str | Path) -> MeasurementLoadReport:
    path = Path(path)
    return load_measurements(path.read_text(encoding="utf-8"), label=path.name)


# --- building one analyte ---------------------------------------------------------------


def _build_analyte(fields: dict[str, object]):
    """Build the right variant of the analyte union from the flattened columns.

    The kind is read from the row and never guessed from which columns happen to
    be filled: a row that fills `analyte_inchikey` and says it is a peptide is a
    fault to be reported, not a small molecule to be inferred.
    """
    raw_kind = fields.pop("kind", None)
    try:
        kind = AnalyteKind(str(raw_kind))
    except ValueError:
        known = ", ".join(k.value for k in AnalyteKind)
        raise _Coercion(UNKNOWN_ANALYTE_KIND, f"analyte_kind={raw_kind!r}; expected one of: {known}") from None
    antibody_fields = {name: fields.pop(name) for name in ("inn", "accession", "sequence") if name in fields}
    model = ANALYTE_TYPES[kind]
    if kind in (AnalyteKind.INTACT_ANTIBODY, AnalyteKind.ADC):
        from .identity import AntibodyIdentity

        fields["antibody"] = AntibodyIdentity(**antibody_fields)
    else:
        fields.update(antibody_fields)
    return model(**fields)


# --- the shared-peak check -------------------------------------------------------------


def _suspected_shared_peaks(
    built: list[tuple[str, CCSMeasurement, str | None]], failed: list[_Failed]
) -> dict[str, str]:
    """Rows whose CCS is identical to a DIFFERENT analyte's under the same conditions.

    "Different" is judged by identity atoms: two rows are one analyte if they
    share any identifier, and are otherwise taken as different. So two spellings
    of one molecule are one molecule and not a shared peak (that is a replicate,
    which is fine), while two rows identified only by different accessions are
    two molecules. Where identity is unclear, "different" is the safe reading: it
    holds a row rather than clearing one.

    The check fails closed over rows that produced no record. A built row is held
    if a row that failed validation carries the same value against an analyte it
    cannot be identified with, because the failed row's calibration group could
    not be computed and so could not be compared. That hold lifts when the failed
    row is fixed.
    """
    flagged: dict[str, str] = {}

    by_key: dict[tuple, list[tuple[str, CCSMeasurement]]] = {}
    for where, record, _ in built:
        by_key.setdefault((record.calibration_group, record.ccs), []).append((where, record))
    for (group, ccs), members in by_key.items():
        components = _components([(where, record.analyte.identity_atoms()) for where, record in members])
        if len(components) < 2:
            continue
        names = ", ".join(sorted(_analyte_name(record.analyte) for _, record in members))
        for where, _ in members:
            flagged[where] = f"CCS {ccs} in {group} is reported for {len(components)} analytes: {names}"

    for where, record, _ in built:
        if where in flagged:
            continue
        for other in failed:
            if other.ccs is None or other.ccs != record.ccs:
                continue
            if len(_components([(where, record.analyte.identity_atoms()), (other.where, other.atoms)])) < 2:
                continue
            flagged[where] = (
                f"CCS {record.ccs} is also reported at {other.where} for a row that did not validate, so its"
                " calibration group could not be compared; held until that row is fixed"
            )
            break
    return flagged


def _lone_conformers(built: list[tuple[str, CCSMeasurement, str | None]]) -> dict[str, str]:
    """Rows claiming a conformer set whose other members are not in the file.

    One molecule can give more than one arrival-time peak under a single set of
    conditions. Those rows are a legitimate pair rather than a collision, and the
    shared-peak check never sees them: it compares rows carrying the SAME value,
    and conformers carry different ones. So this is not that check relaxed - it is
    the opposite check, and it is about completeness.

    A row saying "conformer 1 of 2" whose partner is absent is either a
    transcription that lost a row or a value that is not what it claims, and
    training on it would teach the model that this molecule has a single peak at
    this value. Such rows are held, never averaged and never dropped.

    Siblings are rows agreeing on everything that ought to be identical between
    them: the calibration group, the analyte, the declared total, and the source
    locator, which is what says which sample and which table row a value came
    from. Two samples of one molecule therefore do not pool, provided their
    locators differ.

    KNOWN LIMITATION, carried over and recorded rather than quietly inherited:
    where a file gives no locator the rows DO pool, and the check weakens to
    requiring that the indices present cover the set and that each appears
    equally often. Two samples of one molecule that have each lost a different
    sibling then satisfy that count and are not caught. See LIMITATIONS.md.
    """
    by_key: dict[tuple, list[tuple[str, int | None]]] = {}
    for where, record, _ in built:
        total = record.conformers_total
        if total is None or total < 2:
            continue  # a single conformer claims no siblings
        key = (record.calibration_group, record.analyte.identity_key(), total, record.source_locator)
        by_key.setdefault(key, []).append((where, record.conformer))

    flagged: dict[str, str] = {}
    for (group, analyte, total, _locator), members in by_key.items():
        counts = Counter(index for _, index in members)
        if set(counts) == set(range(1, total + 1)) and len(set(counts.values())) == 1:
            continue
        present = ", ".join(f"{index}x{count}" for index, count in sorted(counts.items()))
        for where, _ in members:
            flagged[where] = (
                f"{analyte} in {group} declares {total} conformers, but the rows present are {present};"
                " every conformer of a set must be in the file before any of them may train"
            )
    return flagged


def _components(members: list[tuple[str, frozenset[str]]]) -> dict[str, list[str]]:
    """Rows joined into analytes by any shared identity atom."""
    union = _Union()
    for where, atoms in members:
        union.find(where)
        for atom in atoms or (f"unidentified:{where}",):
            union.union(where, atom)
    components: dict[str, list[str]] = {}
    for where, _ in members:
        components.setdefault(union.find(where), []).append(where)
    return components


def _analyte_name(analyte) -> str:
    """A readable name for a report, never for a key."""
    key = analyte.identity_key()
    return ":".join(str(part) for part in key[1:] if part)


def _as_failed(where: str, cells: Mapping[str, str]) -> _Failed:
    """As much of a failed row as the shared-peak check can use: its value and its identifiers."""
    ccs = None
    text = cells.get("ccs", "")
    if _FLOAT.fullmatch(text) and math.isfinite(float(text)):
        ccs = float(text)
    atoms = set()
    for column, namespace in (
        ("analyte_inchikey", "inchikey"),
        ("analyte_wurcs", "wurcs"),
        ("analyte_glytoucan_ac", "glytoucan"),
        ("analyte_iupac_condensed", "iupac"),
        ("analyte_composition", "composition"),
        ("analyte_accession", "accession"),
        ("analyte_sequence", "sequence"),
        ("analyte_inn", "inn"),
    ):
        if value := cells.get(column):
            atoms.add(f"{namespace}:{value}")
    return _Failed(where=where, ccs=ccs, atoms=frozenset(atoms))


# --- one row ---------------------------------------------------------------------------


def _coerce(cells: Mapping[str, str]) -> tuple[dict[str, object], dict[str, object], dict[str, object], str | None]:
    """Split one row into the analyte's fields, the calibration reference's, the measurement's, and a flag.

    A blank cell is absent, so the model's default applies. Numbers, booleans and
    dates are read in one spelling each and refused otherwise. Any other text is
    passed through untouched, so a filler like "n/a" reaches validation and is
    refused there with its own message, which is where a reader will look for it.
    """
    analyte: dict[str, object] = {}
    calref: dict[str, object] = {}
    measurement: dict[str, object] = {}
    flag: str | None = None
    for column, text in cells.items():
        if not text:
            continue
        value: object = text
        if column in _STATUS_COLUMNS and text.casefold() == ReuseStatus.SYNTHETIC_FIXTURE.value:
            # Refused before a record exists: a synthetic fixture is built in code, never read from a file.
            raise _Coercion(SYNTHETIC_STATUS_IN_FILE, f"{column}={text!r}")
        if column in _BOOLEAN_COLUMNS:
            key = text.casefold()
            if key not in _BOOLEANS:
                raise _Coercion(NOT_A_BOOLEAN, f"{column}={text!r}")
            value = _BOOLEANS[key]
        elif column in _INTEGER_COLUMNS:
            if not _INTEGER.fullmatch(text):
                raise _Coercion(NOT_A_NUMBER, f"{column}={text!r}")
            value = int(text)
        elif column in _FLOAT_COLUMNS:
            if not _FLOAT.fullmatch(text) or not math.isfinite(float(text)):
                raise _Coercion(NOT_A_NUMBER, f"{column}={text!r}")
            value = float(text)
        elif column in _DATE_COLUMNS:
            # Left as text, pydantic would read a run of digits as a Unix timestamp.
            if not _ISO_DATE.fullmatch(text):
                raise _Coercion(NOT_A_DATE, f"{column}={text!r}")
            try:
                value = date.fromisoformat(text)
            except ValueError:
                raise _Coercion(NOT_A_DATE, f"{column}={text!r}") from None
        elif column in _LIST_COLUMNS:
            value = tuple(part.strip() for part in text.split(";") if part.strip())

        if column in ANALYTE_COLUMNS:
            analyte[column[len(ANALYTE_PREFIX) :]] = value
        elif column in CALIBRATION_REFERENCE_COLUMNS:
            calref[column[len("calref_") :]] = value
        elif column in MEASUREMENT_COLUMNS:
            measurement[column] = value
        elif column == "curation_flag":
            flag = text
        # an unrecognised column was already reported from the header; its cells are ignored
    return analyte, calref, measurement, flag


def _first_line(exc: BaseException) -> str:
    """The reason, not the count.

    A pydantic error says "1 validation error" on its first line and puts the
    message underneath, so the first line alone would tell a reader nothing about
    what to fix.
    """
    errors = getattr(exc, "errors", None)
    if callable(errors):
        try:
            first = errors()[0]
            where = ".".join(str(part) for part in first.get("loc", ())) or "record"
            return f"{where}: {first.get('msg', '')}"[:400]
        except Exception:
            pass
    text = str(exc).strip().splitlines()
    return (text[0] if text else type(exc).__name__)[:400]


def _note(counts: Counter[str], examples: dict[str, list[str]], reason: str, example: str) -> None:
    counts[reason] += 1
    kept = examples.setdefault(reason, [])
    if len(kept) < 3:
        kept.append(example)
