"""SugarBase, read from the installed glycowork package.

No data is copied into this repository. glycowork ships SugarBase, pyproject
pins the glycowork version, and that version goes into the dataset version
string that every record's provenance carries, so a record can always be traced
back to the release it came from.

glycowork and its bundled data are MIT licensed (Copyright (c) 2021 Daniel
Bojar), which permits commercial training and redistribution with attribution,
so records load with reuse status open_attribution.
"""

from __future__ import annotations

import ast
from collections import Counter
from dataclasses import dataclass, field
from importlib.metadata import version
from typing import Iterable, Iterator, Mapping

from .composition import Composition, UnsupportedResidueError
from .glycan_graph import GlycanGraph, GlycanGraphError
from .models import GlycanStructure
from .reuse import ReuseStatus  # from below the gate, which imports the registry, which imports this loader

SUGARBASE_RELEASE = "v12"
SUGARBASE_FILE = "glycan_data/v12_sugarbase.json"
SUGARBASE_LICENCE = "MIT"
SUGARBASE_ATTRIBUTION = "Copyright (c) 2021 Daniel Bojar"
N_LINKED = "N"

# Why a row did not become a record.
PARSE_FAILED = "the structure string did not parse"
UNSUPPORTED_RESIDUES = "residues outside the five composition classes"
COMPOSITION_DISAGREES = "our composition disagrees with the one SugarBase records"
RECORD_REJECTED = "the record failed validation"

# SugarBase's own composition keys, in our terms.
_COMPOSITION_KEYS = {"Hex": "hex", "HexNAc": "hexnac", "dHex": "fuc", "Neu5Ac": "neuac", "Neu5Gc": "neugc"}


@dataclass(frozen=True)
class DatasetVersion:
    """Which release of which dataset a record came from, and on what terms."""

    name: str
    release: str
    package: str
    package_version: str
    file: str
    licence: str
    attribution: str

    def __str__(self) -> str:
        return f"{self.name} {self.release} via {self.package} {self.package_version}"

    @property
    def provenance(self) -> str:
        """What goes in a record's source field."""
        return f"{self} ({self.file})"


def dataset_version() -> DatasetVersion:
    """The dataset version for the glycowork release that is actually installed."""
    return DatasetVersion(
        name="SugarBase",
        release=SUGARBASE_RELEASE,
        package="glycowork",
        package_version=version("glycowork"),
        file=SUGARBASE_FILE,
        licence=SUGARBASE_LICENCE,
        attribution=SUGARBASE_ATTRIBUTION,
    )


@dataclass(frozen=True)
class LoadReport:
    """What a load produced, and what it could not."""

    dataset_version: DatasetVersion
    structures: tuple[GlycanStructure, ...]
    rows_in_dataset: int = 0  # every row the dataset holds
    rows_selected: int = 0  # those the glycan_types filter kept, and so the only ones this report accounts for
    rows_n_linked: int = 0
    graphs_parsed: int = 0
    records_built: int = 0
    records_n_linked: int = 0
    records_with_accession: int = 0
    failure_counts: Mapping[str, int] = field(default_factory=dict)
    failure_examples: Mapping[str, tuple[str, ...]] = field(default_factory=dict)

    @property
    def rows_failed(self) -> int:
        return sum(self.failure_counts.values())

    def summary(self) -> str:
        lines = [
            f"{self.dataset_version}  [{self.dataset_version.licence}, {self.dataset_version.attribution}]",
            f"  rows in the dataset      {self.rows_in_dataset}",
            f"  selected by type         {self.rows_selected}",
            f"  N-linked rows            {self.rows_n_linked}",
            f"  parsed into graphs       {self.graphs_parsed}",
            f"  records built            {self.records_built}",
            f"    of which N-linked      {self.records_n_linked}",
            f"    with a GlyTouCan id    {self.records_with_accession}",
            f"  rows not loaded          {self.rows_failed}",
        ]
        for reason, count in sorted(self.failure_counts.items(), key=lambda item: -item[1]):
            lines.append(f"    {reason}: {count}")
            for example in self.failure_examples.get(reason, ())[:2]:
                lines.append(f"      e.g. {example}")
        return "\n".join(lines)


def load_sugarbase(
    glycan_types: Iterable[str] | None = (N_LINKED,),
    limit: int | None = None,
    rows: Iterable[Mapping[str, object]] | None = None,
) -> LoadReport:
    """Load SugarBase structures as GlycanStructure records.

    `glycan_types` filters on SugarBase's own typing ("N", "O", ...); None
    loads every row. `rows` accepts an iterable of row mappings instead of the
    installed package, which is how the tests drive the loader without pandas.

    Nothing is guessed: a structure whose residues have no place in a
    Hex/HexNAc/Fuc/NeuAc/NeuGc composition is reported as a failure rather than
    rounded off, and the composition read from the structure is checked against
    the one SugarBase records.
    """
    wanted = set(glycan_types) if glycan_types is not None else None
    dataset = dataset_version()
    structures: list[GlycanStructure] = []
    failures: Counter[str] = Counter()
    examples: dict[str, list[str]] = {}
    rows_in_dataset = rows_selected = rows_n_linked = graphs_parsed = records_n_linked = with_accession = 0

    for row in _rows(rows):
        rows_in_dataset += 1
        glycan_type = _text(row.get("glycan_type"))
        if wanted is not None and glycan_type not in wanted:
            continue
        rows_selected += 1
        if glycan_type == N_LINKED:
            rows_n_linked += 1
        glycan = _text(row.get("glycan")) or ""

        try:
            graph = GlycanGraph.from_iupac_condensed(glycan)
        except GlycanGraphError as exc:
            _note(failures, examples, PARSE_FAILED, f"{glycan[:60]}: {exc}")
            continue
        graphs_parsed += 1

        try:
            composition = graph.composition()
        except UnsupportedResidueError as exc:
            _note(failures, examples, UNSUPPORTED_RESIDUES, f"{glycan[:60]}: {', '.join(exc.names)}")
            continue

        recorded = _recorded_composition(row.get("Composition"))
        if recorded is not None and recorded != composition:
            _note(failures, examples, COMPOSITION_DISAGREES, f"{glycan[:60]}: {composition} vs {recorded}")
            continue

        accession = _text(row.get("glytoucan_id"))
        try:
            structure = GlycanStructure(
                composition=composition,
                glytoucan_ac=accession,
                iupac_condensed=glycan,
                has_unresolved_linkage=graph.has_unresolved_linkage,
                has_unresolved_anomericity=graph.has_unresolved_anomericity,
                source=dataset.provenance,
                reuse_status=ReuseStatus.OPEN_ATTRIBUTION,
            )
        except Exception as exc:
            _note(failures, examples, RECORD_REJECTED, f"{glycan[:60]}: {str(exc)[:120]}")
            continue

        structures.append(structure)
        if glycan_type == N_LINKED:
            records_n_linked += 1
        if structure.glytoucan_ac is not None:
            with_accession += 1
        if limit is not None and len(structures) >= limit:
            break

    return LoadReport(
        dataset_version=dataset,
        structures=tuple(structures),
        rows_in_dataset=rows_in_dataset,
        rows_selected=rows_selected,
        rows_n_linked=rows_n_linked,
        graphs_parsed=graphs_parsed,
        records_built=len(structures),
        records_n_linked=records_n_linked,
        records_with_accession=with_accession,
        failure_counts=dict(failures),
        failure_examples={reason: tuple(found) for reason, found in examples.items()},
    )


def _rows(rows: Iterable[Mapping[str, object]] | None) -> Iterator[Mapping[str, object]]:
    if rows is not None:
        yield from rows
        return
    from glycowork.glycan_data.loader import df_glycan  # imported here: it pulls in pandas

    for record in df_glycan.to_dict(orient="records"):
        yield record


def _text(value: object) -> str | None:
    """A trimmed string, or None for anything a table uses to mean "not there"."""
    if value is None or not isinstance(value, str):
        return None
    trimmed = value.strip()
    return trimmed or None


def _recorded_composition(value: object) -> Composition | None:
    """SugarBase's own composition for the row, when it maps onto our five classes."""
    if isinstance(value, str):
        try:
            value = ast.literal_eval(value)
        except (ValueError, SyntaxError):
            return None
    if not isinstance(value, Mapping):
        return None
    if set(value) - set(_COMPOSITION_KEYS):
        return None  # holds residues we cannot represent; the residue check already reports those
    try:
        return Composition(**{_COMPOSITION_KEYS[key]: int(count) for key, count in value.items()})
    except Exception:
        return None


def _note(counts: Counter[str], examples: dict[str, list[str]], reason: str, example: str) -> None:
    counts[reason] += 1
    kept = examples.setdefault(reason, [])
    if len(kept) < 3:
        kept.append(example)
