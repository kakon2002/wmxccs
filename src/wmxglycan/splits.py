"""Dividing records into folds so that a reported error means something.

The rule the whole module serves: if a model cannot tell two records apart, they
must not sit on opposite sides of a split. A test record answered by an
identical training record measures memorisation, and the error it produces is
not an error at all.

That fixes the grouping key, and it fixes it lower than intuition suggests. A
prediction request carries measurement conditions and a composition, and no
structure at all, so the analyte half of any row a model can ever see at request
time is a function of composition alone. Two isomers of one composition are
byte-identical to it. Composition is therefore the coarsest key that is still
correct, and every structure-derived key is FINER than it: keying on a canonical
structure key or an accession does not merely fail to help, it manufactures the
leak, by scattering the isomers of one composition across folds.

This is also why there is no fallback ladder. "Use the accession, else the
structure, else the composition" gives two records of one glycan different keys
whenever they carry different kinds of evidence, and lets that glycan straddle.
Every record emits every atom it can support, and records sharing any atom are
one group: the group is a connected component, and adding an atom can only ever
merge, never split.

Two axes compose as the JOIN of their partitions, which is what connected
components compute. A compound (analyte, source) tuple is the MEET, strictly
finer than either axis, and so prevents neither leak while looking rigorous.

Floors and tolerances are module constants and never function parameters. A
floor that can be passed in is a floor that can be passed a zero.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from enum import StrEnum
from fractions import Fraction
from typing import Iterable, Mapping, Sequence

from .composition import Composition
from .features import FeatureVector, extract
from .glycan_graph import GlycanGraph, GlycanGraphError
from .licensing import LicenceGateError, TrainingGateError, UnbackedClaimError, assert_trainable

N_FOLDS = 5
MIN_GROUPS_PER_FOLD = 5  # POLICY, chosen rather than measured: fewer than five groups in a fold
# makes a fold's error a statement about a handful of glycans.
MIN_ANALYTE_GROUPS = N_FOLDS * MIN_GROUPS_PER_FOLD  # derived, not chosen
# How far a realised fold may sit from an equal share. An exact fraction, not a
# float: 1/5 - 0.05 evaluates to 0.15000000000000002 in binary floating point, so
# a fold holding exactly 15 per cent would be refused by a message reading
# "15.0% of the records, outside 15.0% to 25.0%", which denies itself.
FOLD_TOLERANCE = Fraction(1, 20)
_EQUAL_SHARE = Fraction(1, N_FOLDS)
_FOLD_LOW = _EQUAL_SHARE - FOLD_TOLERANCE
_FOLD_HIGH = _EQUAL_SHARE + FOLD_TOLERANCE
REPORTED_LEVELS: tuple[AnalyteLevel, ...]  # filled in below, once the enum exists

_DOI = re.compile(r"10\.\d{4,9}/[^\s\"'<>,;]+")

# Refusal wording is an interface: written once, imported by tests rather than matched as prose.
NOTHING_TO_SPLIT = "there are no records to split"
TOO_FEW_GROUPS = "{groups} {unit} is too few to split into {folds} folds; at least {needed} are needed{detail}"
WELDED_BY_SOURCES = (
    ". The corpus holds {keys} analyte keys, welded into {groups} components by shared sources;"
    " that is what a cross-study split costs here"
)
GROUP_OVER_BUDGET = (
    "the largest group {group!r} holds {share:.1%} of the records, over the {budget:.1%} a fold can take;"
    " a group is never divided to hit a fraction"
)
FOLDS_UNEVEN = "fold {fold} holds {share:.1%} of the records, outside {low:.1%} to {high:.1%}"
FEATURES_STRADDLE = (
    "records {first} and {second} have the same analyte features but are in groups {group_a!r} and {group_b!r};"
    " a model cannot tell them apart, so a split that separates them measures memorisation"
)
GROUP_STRADDLES = "group {group!r} appears in folds {folds}"
UNFEATURISABLE_RECORD = "record {index} cannot be featurised, so it cannot be checked for leakage: {reason}"
NO_STRUCTURE_FOR_ISOMER_MODE = (
    "record {index} carries no structure, so it cannot take part in an isomer-discrimination split:"
    " the split groups on the structure, and a composition-only record has none. Either resolve it or"
    " run the deployment split, which groups on composition"
)
NO_CROSS_STUDY_SPLIT = (
    "this corpus cannot support a cross-study split: {groups} group(s) over {sources} source(s)."
    " The number beside it is a within-study number. {bridge}Use leave_one_source_out for the practical"
    " cross-study estimate"
)


class AnalyteLevel(StrEnum):
    """How coarsely one analyte is defined for the purpose of grouping."""

    # The featuriser's own kernel, forced by what a prediction request can carry.
    COMPOSITION = "composition"
    # The elaboration residues dropped, so the fucosylation and sialylation
    # series of one backbone stay together. A judgement about analogue leakage,
    # reported beside the headline rather than instead of it.
    BACKBONE = "backbone"


REPORTED_LEVELS = (AnalyteLevel.COMPOSITION, AnalyteLevel.BACKBONE)


class SplitMode(StrEnum):
    """Which question a split answers. The two are not interchangeable.

    DEPLOYMENT is the real task: an unseen composition arrives, candidates are
    enumerated, each is scored. Records are grouped on composition, so no test
    composition was ever trained on. A number measured this way is the only one
    that may be reported as performance.

    ISOMER_DISCRIMINATION asks a different question: is there any learnable CCS
    signal separating isomers at all. It deliberately puts isomers of one
    composition on both sides, grouping on the structure instead, so the model is
    asked to tell apart molecules whose near relatives it has already seen. That
    is an optimistic bound and a diagnostic, never a performance estimate, and
    the reporting path refuses to let it be quoted as one.
    """

    DEPLOYMENT = "deployment"
    ISOMER_DISCRIMINATION = "isomer_discrimination"


class AtomKind(StrEnum):
    ANALYTE = "analyte"
    GLYTOUCAN = "glytoucan"
    WURCS = "wurcs"
    STRUCTURE = "structure"
    DOI = "doi"
    STUDY = "study"


class SplitError(TrainingGateError):
    """A split cannot be made, or the one that was made is not honest."""


class SplitRefused(SplitError):
    """The data cannot support a split that would mean anything."""


class LeakageError(SplitError):
    """A split that was made, or would be made, puts indistinguishable records on both sides."""


# --- atoms --------------------------------------------------------------------


def _composition_of(subject: object) -> Composition:
    glycan = getattr(subject, "glycan", None)
    if glycan is not None:
        return glycan.composition
    if isinstance(subject, Composition):
        return subject
    composition = getattr(subject, "composition", None)
    if isinstance(composition, Composition):
        return composition
    raise TypeError(f"{type(subject).__name__} carries no composition to group on")


def analyte_key(subject: object, level: AnalyteLevel = AnalyteLevel.COMPOSITION) -> str:
    """The analyte this record is about, at the requested coarseness.

    Reads the composition and nothing else: not the conditions, not the measured
    value, not the provenance, not the licence. Total on every composition that
    can exist, so no record is ever dropped for want of a key.
    """
    composition = _composition_of(subject)
    if level is AnalyteLevel.BACKBONE:
        return f"Hex{composition.hex}HexNAc{composition.hexnac}"
    return composition.canonical


def analyte_atoms(record: object, level: AnalyteLevel = AnalyteLevel.COMPOSITION) -> tuple[str, ...]:
    """Exactly one, always. This is the floor every record stands on."""
    return (f"{AtomKind.ANALYTE}:{analyte_key(record, level)}",)


def identity_atoms(record: object, notes: list[str] | None = None) -> tuple[str, ...]:
    """Structure identifiers, as a merge-only safety net.

    An identifier can only ever join records; it can never separate them. Where
    two records of genuinely different analytes share one identifier, merging is
    the safe direction, because it costs statistical power rather than honesty.
    Such a merge is reported as a conflict and never silently corrected:
    correcting an accession is a data decision with a source behind it.
    """
    glycan = getattr(record, "glycan", None) or record
    atoms: list[str] = []
    accession = getattr(glycan, "glytoucan_ac", None)
    if accession:
        # Kept even when malformed, and upper-cased, because merging is safe.
        atoms.append(f"{AtomKind.GLYTOUCAN}:{accession.strip().upper()}")
    wurcs = getattr(glycan, "wurcs", None)
    if wurcs:
        atoms.append(f"{AtomKind.WURCS}:{wurcs}")
    text = getattr(glycan, "iupac_condensed", None)
    if text:
        try:
            atoms.append(f"{AtomKind.STRUCTURE}:{GlycanGraph.from_iupac_condensed(text).canonical_key()}")
        except GlycanGraphError as exc:
            if notes is not None:
                notes.append(f"structure did not parse, so it contributes no atom: {exc}")
    return tuple(atoms)


def _trim_doi(found: str) -> str:
    """Drop punctuation a citation left on the end of a DOI.

    A DOI written in prose arrives as "(10.1000/paper-a)" or "10.1000/paper-a.",
    and keeping the bracket or the full stop makes one paper look like several
    sources. That is the unsafe direction: it FAILS to merge, so a corpus can
    clear the group floor on study diversity that does not exist.

    A closing bracket is dropped only when it is unmatched, because a DOI suffix
    may legitimately contain one.
    """
    changed = True
    while changed and found:
        changed = False
        if found[-1] in ".,;:]}>":
            found, changed = found[:-1], True
        elif found.endswith(")") and found.count(")") > found.count("("):
            found, changed = found[:-1], True
    return found


def grouping_atoms(
    record: object,
    index: int,
    level: AnalyteLevel = AnalyteLevel.COMPOSITION,
    mode: SplitMode = SplitMode.DEPLOYMENT,
    notes: list[str] | None = None,
) -> tuple[str, ...]:
    """The atoms that decide which group a record lands in, for this mode.

    In deployment mode the floor is the composition, and identity atoms may only
    merge on top of it. In isomer mode the floor is the STRUCTURE, because the
    whole point is to place isomers of one composition on opposite sides; adding
    the composition atom there would weld them straight back together.
    """
    if mode is SplitMode.ISOMER_DISCRIMINATION:
        structures = [
            atom for atom in identity_atoms(record, notes) if atom.startswith(f"{AtomKind.STRUCTURE}:")
        ]
        if not structures:
            raise SplitRefused(NO_STRUCTURE_FOR_ISOMER_MODE.format(index=index))
        return (structures[0],)
    return analyte_atoms(record, level) + identity_atoms(record, notes)


def provenance_atoms(record: object, aliases: Mapping[str, str] | None = None) -> tuple[str, ...]:
    """Where the MEASUREMENT came from, as a study atom.

    The glycan's own source is deliberately not used: it records where the
    structure assignment came from, so counting it would put every record whose
    structure came from one database into a single group.

    A cohort alias may only merge, and only from a human-curated file with a
    citation per entry. Nothing is inferred from how a source string is spelt.
    """
    # A recorded DOI is a stated fact and is preferred outright. Scraping one out
    # of free text is the fallback, and it is the weaker path: punctuation left by
    # a citation used to split one paper into several sources.
    stated = getattr(record, "doi", None)
    if stated:
        return (f"{AtomKind.DOI}:{stated.strip().casefold()}",)
    source = getattr(record, "source", None)
    if not source:
        return ()
    text = source.strip()
    doi = _DOI.search(text)
    if doi:
        return (f"{AtomKind.DOI}:{_trim_doi(doi.group(0)).casefold()}",)
    key = text.casefold()
    if aliases and key in aliases:
        return (f"{AtomKind.STUDY}:{aliases[key]}",)
    return (f"{AtomKind.STUDY}:{key}",)


# --- union-find ----------------------------------------------------------------


class _Union:
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
            # Smallest root wins, so the component id is stable and order-independent.
            low, high = sorted((root_a, root_b))
            self._parent[high] = low


# --- grouping -------------------------------------------------------------------


@dataclass(frozen=True)
class Grouping:
    """Which records belong together, and why they were joined."""

    level: AnalyteLevel
    by_study: bool
    records: tuple[object, ...]
    group_of: tuple[str, ...]  # parallel to records
    groups: Mapping[str, tuple[int, ...]]
    sources: Mapping[str, int]
    bridges: Mapping[str, int]  # analyte atom -> how many distinct sources it joins
    identity_conflicts: Mapping[str, tuple[str, ...]]
    notes: tuple[str, ...] = ()
    mode: SplitMode = SplitMode.DEPLOYMENT

    @property
    def largest_share(self) -> float:
        if not self.records:
            return 0.0
        return max(len(members) for members in self.groups.values()) / len(self.records)

    def summary(self) -> str:
        lines = [
            f"{len(self.records)} records in {len(self.groups)} groups at {self.level} level"
            f"{', study axis on' if self.by_study else ''}",
            f"  largest group holds {self.largest_share:.1%}",
            f"  {len(self.sources)} distinct sources",
        ]
        for atom, count in sorted(self.identity_conflicts.items()):
            lines.append(f"  identity conflict on {atom}: {', '.join(count)}")
        return "\n".join(lines)


def _gate(records: Sequence[object]) -> None:
    """Run the licence gate over every record, reporting every fault at once.

    Before grouping, not after: an ungated record could otherwise be the bridge
    that silently changes which sources merge.
    """
    faults: list[str] = []
    licence_fault = False
    unbacked_only = True
    for index, record in enumerate(records):
        try:
            assert_trainable(record)
        except TrainingGateError as exc:
            faults.append(f"record {index}: {exc}")
            licence_fault = licence_fault or isinstance(exc, LicenceGateError)
            unbacked_only = unbacked_only and isinstance(exc, UnbackedClaimError)
    if faults:
        problem = "; ".join(faults)
        # The most specific class the faults allow, so an unbacked claim is not
        # reported as a licence refusal, which points at the wrong remedy.
        error = UnbackedClaimError if unbacked_only else LicenceGateError if licence_fault else TrainingGateError
        raise error(f"{len(faults)} record(s) may not be used for training: {problem}")


def group_records(
    records: Iterable[object],
    level: AnalyteLevel = AnalyteLevel.COMPOSITION,
    by_study: bool = False,
    aliases: Mapping[str, str] | None = None,
    mode: SplitMode = SplitMode.DEPLOYMENT,
) -> Grouping:
    """Partition records into groups that must not be split across folds."""
    records = tuple(records)
    # Emptiness is checked here and never delegated to the licence gate: a
    # wrapper with a trainable status and no component records passes the gate
    # silently, so "the gate passed" must never be read as "there is data".
    if not records:
        raise SplitRefused(NOTHING_TO_SPLIT)
    _gate(records)

    union = _Union()
    notes: list[str] = []
    atoms_of: list[tuple[str, ...]] = []
    analyte_of: list[str] = []
    for index, record in enumerate(records):
        atoms = grouping_atoms(record, index, level, mode, notes)
        if by_study:
            atoms += provenance_atoms(record, aliases)
        atoms_of.append(atoms)
        analyte_of.append(atoms[0])
        for atom in atoms:
            union.union(atoms[0], atom)

    members: dict[str, list[int]] = defaultdict(list)
    group_of: list[str] = []
    for index, atoms in enumerate(atoms_of):
        group = union.find(atoms[0])
        group_of.append(group)
        members[group].append(index)

    # An identifier that joins two different analyte keys means two records
    # disagree about the composition of one glycan: a transcription error.
    conflicts: dict[str, tuple[str, ...]] = {}
    seen: dict[str, set[str]] = defaultdict(set)
    for index, atoms in enumerate(atoms_of):
        for atom in atoms[1:]:
            if atom.startswith((f"{AtomKind.GLYTOUCAN}:", f"{AtomKind.WURCS}:", f"{AtomKind.STRUCTURE}:")):
                seen[atom].add(analyte_of[index])
    for atom, keys in seen.items():
        if len(keys) > 1:
            conflicts[atom] = tuple(sorted(keys))

    sources = Counter(
        atom for atoms in atoms_of for atom in atoms if atom.startswith((f"{AtomKind.DOI}:", f"{AtomKind.STUDY}:"))
    )
    bridges: dict[str, int] = {}
    if by_study:
        joined: dict[str, set[str]] = defaultdict(set)
        for index, atoms in enumerate(atoms_of):
            for atom in atoms:
                if atom.startswith((f"{AtomKind.DOI}:", f"{AtomKind.STUDY}:")):
                    joined[analyte_of[index]].add(atom)
        bridges = {atom: len(found) for atom, found in joined.items() if len(found) > 1}

    return Grouping(
        level=level,
        by_study=by_study,
        mode=mode,
        records=records,
        group_of=tuple(group_of),
        groups={group: tuple(indices) for group, indices in sorted(members.items())},
        sources=dict(sources),
        bridges=bridges,
        identity_conflicts=conflicts,
        notes=tuple(notes),
    )


# --- the guards ------------------------------------------------------------------


def assert_features_do_not_straddle(grouping: Grouping, rows: Sequence[tuple[float, ...]]) -> None:
    """The law: records a model cannot tell apart must be in one group.

    The law does not change between modes; what the model can SEE does, so the
    caller passes the row the model actually gets. In deployment mode that is the
    analyte row, because the fitted structure columns are the same for every
    candidate of one composition at grouping time. In isomer mode it is the whole
    fitted row, which distinguishes isomers, so identical structures still cannot
    straddle while isomers of one composition deliberately do.
    """
    if len(rows) != len(grouping.records):
        raise ValueError(f"{len(rows)} rows for {len(grouping.records)} records")
    first_seen: dict[tuple[float, ...], int] = {}
    for index, row in enumerate(rows):
        earlier = first_seen.setdefault(row, index)
        if grouping.group_of[earlier] != grouping.group_of[index]:
            raise LeakageError(
                FEATURES_STRADDLE.format(
                    first=earlier,
                    second=index,
                    group_a=grouping.group_of[earlier],
                    group_b=grouping.group_of[index],
                )
            )


def assert_no_group_straddles(folds: Sequence[Sequence[int]], group_of: Sequence[str]) -> None:
    """No group may appear in more than one fold."""
    where: dict[str, set[int]] = defaultdict(set)
    for number, members in enumerate(folds):
        for index in members:
            where[group_of[index]].add(number)
    for group, found in sorted(where.items()):
        if len(found) > 1:
            raise LeakageError(GROUP_STRADDLES.format(group=group, folds=sorted(found)))


# --- folds -------------------------------------------------------------------------


@dataclass(frozen=True)
class SplitReport:
    """The folds, and everything needed to judge whether they mean anything."""

    level: AnalyteLevel
    by_study: bool
    mode: SplitMode
    folds: tuple[tuple[int, ...], ...]
    assignment: tuple[int, ...]  # parallel to records: which fold each is in
    group_of: tuple[str, ...]
    group_sizes: Mapping[str, int]
    fold_sizes: tuple[int, ...]
    records: tuple[object, ...]
    distinct_compositions: int
    distinct_canonical_keys: int
    sources: Mapping[str, int]
    bridges: Mapping[str, int]
    identity_conflicts: Mapping[str, tuple[str, ...]]
    refusals: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    @property
    def fold_shares(self) -> tuple[float, ...]:
        total = len(self.records) or 1
        return tuple(size / total for size in self.fold_sizes)

    @property
    def largest_group(self) -> tuple[str, int]:
        return max(self.group_sizes.items(), key=lambda item: (item[1], item[0]))

    def summary(self) -> str:
        group, size = self.largest_group
        lines = [
            f"{len(self.records)} records, {len(self.group_sizes)} groups, {N_FOLDS} folds"
            f" in {self.mode} mode at {self.level} level{', study axis on' if self.by_study else ''}",
            f"  fold sizes {list(self.fold_sizes)}  shares {[f'{s:.1%}' for s in self.fold_shares]}",
            f"  largest group {group!r} holds {size} records",
            f"  {self.distinct_compositions} compositions, {self.distinct_canonical_keys} distinct structures",
        ]
        for refusal in self.refusals:
            lines.append(f"  finding: {refusal}")
        return "\n".join(lines)


def _pack(groups: Mapping[str, Sequence[int]]) -> tuple[list[list[int]], list[str]]:
    """Largest group first into the emptiest fold. Deterministic, and seedless.

    No RNG and no hash anywhere: an auditor re-running this gets byte-identical
    folds, and a flattering split cannot be produced by reseeding.
    """
    folds: list[list[int]] = [[] for _ in range(N_FOLDS)]
    order = sorted(groups.items(), key=lambda item: (-len(item[1]), item[0]))
    placed: list[str] = []
    for group, members in order:
        target = min(range(N_FOLDS), key=lambda number: (len(folds[number]), number))
        folds[target].extend(members)
        placed.append(group)
    return folds, placed


def grouped_folds(
    records: Iterable[object],
    level: AnalyteLevel = AnalyteLevel.COMPOSITION,
    by_study: bool = False,
    aliases: Mapping[str, str] | None = None,
    mode: SplitMode = SplitMode.DEPLOYMENT,
) -> SplitReport:
    """Folds that no group straddles, or a refusal saying why none can be made."""
    grouping = group_records(records, level=level, by_study=by_study, aliases=aliases, mode=mode)
    total = len(grouping.records)

    # Assessed BEFORE the group-count guard. Turning the study axis on welds
    # analyte keys together through shared sources, so the component count can
    # fall below the floor while the corpus holds plenty of distinct analytes.
    # Reporting the component count alone would state something untrue of the data.
    analyte_keys = len({analyte_key(record, level) for record in grouping.records})
    welded = by_study and len(grouping.groups) < analyte_keys
    findings: list[str] = []
    if welded:
        bridge = ""
        if grouping.bridges:
            worst = max(grouping.bridges.items(), key=lambda item: (item[1], item[0]))
            bridge = f"The analyte {worst[0]!r} alone joins {worst[1]} sources. "
        findings.append(
            NO_CROSS_STUDY_SPLIT.format(
                groups=len(grouping.groups), sources=len(grouping.sources), bridge=bridge
            )
        )

    if len(grouping.groups) < MIN_ANALYTE_GROUPS:
        raise SplitRefused(
            TOO_FEW_GROUPS.format(
                groups=len(grouping.groups),
                unit="grouping components" if by_study else "analyte groups",
                folds=N_FOLDS,
                needed=MIN_ANALYTE_GROUPS,
                detail=WELDED_BY_SOURCES.format(keys=analyte_keys, groups=len(grouping.groups)) if welded else "",
            )
        )
    biggest = max(grouping.groups.items(), key=lambda item: (len(item[1]), item[0]))
    if Fraction(len(biggest[1]), total) > _FOLD_HIGH:
        raise SplitRefused(
            GROUP_OVER_BUDGET.format(
                group=biggest[0], share=len(biggest[1]) / total, budget=float(_FOLD_HIGH)
            )
        )

    rows: list[tuple[float, ...]] = []
    for index, record in enumerate(grouping.records):
        try:
            vector = extract(record)
            rows.append(vector.as_row() if mode is SplitMode.ISOMER_DISCRIMINATION else vector.analyte_row())
        except Exception as exc:  # a record the featuriser cannot render is named, never made a singleton
            raise SplitRefused(UNFEATURISABLE_RECORD.format(index=index, reason=exc)) from exc
    assert_features_do_not_straddle(grouping, rows)

    folds, _ = _pack(grouping.groups)
    for number, members in enumerate(folds):
        # Currently unreachable, and kept deliberately. The per-group budget above
        # refuses any group over 1/N + tolerance, and with every group inside that
        # bound and at least MIN_ANALYTE_GROUPS of them, largest-first packing has
        # never produced a fold outside the window: 400,000 random corpora and a
        # set of adversarial ones (one group sitting exactly at the budget, the
        # rest as lumpy as allowed) all balanced. A mutation deleting this check
        # therefore survives the suite, which is a statement about reachability
        # rather than a missing test. It stays because it is the guard that would
        # catch a future change to the budget, the fold count or the packing
        # order; if the budget guard is ever loosened, this is what starts firing.
        # Exact rational comparison, so a fold sitting on the boundary is inside it.
        if not _FOLD_LOW <= Fraction(len(members), total) <= _FOLD_HIGH:
            raise SplitRefused(
                FOLDS_UNEVEN.format(
                    fold=number, share=len(members) / total, low=float(_FOLD_LOW), high=float(_FOLD_HIGH)
                )
            )
    assert_no_group_straddles(folds, grouping.group_of)

    assignment = [0] * total
    for number, members in enumerate(folds):
        for index in members:
            assignment[index] = number

    keys = {atom for record in grouping.records for atom in identity_atoms(record) if atom.startswith("structure:")}
    return SplitReport(
        level=level,
        by_study=by_study,
        mode=mode,
        folds=tuple(tuple(sorted(members)) for members in folds),
        assignment=tuple(assignment),
        group_of=grouping.group_of,
        group_sizes={group: len(members) for group, members in grouping.groups.items()},
        fold_sizes=tuple(len(members) for members in folds),
        records=grouping.records,
        distinct_compositions=len({analyte_key(record) for record in grouping.records}),
        distinct_canonical_keys=len(keys),
        sources=grouping.sources,
        bridges=grouping.bridges,
        identity_conflicts=grouping.identity_conflicts,
        refusals=tuple(findings),
        notes=grouping.notes,
    )


@dataclass(frozen=True)
class SourceHoldout:
    """One source held out, and how much of it the training side already contains.

    The overlap is reported rather than refused. Analyte overlap between studies
    is intrinsic to cross-study cross-validation, because every CCS paper
    measures the common N-glycans; refusing on it would make this function
    useless on every real corpus. Returning bare folds and saying nothing is the
    one option that is not honest, so the contamination travels with the split.
    """

    source: str
    held_out: tuple[int, ...]
    trained_on: tuple[int, ...]
    overlapping_records: tuple[int, ...]  # held-out records whose analyte is also in training
    overlapping_analytes: tuple[str, ...]

    @property
    def clean(self) -> bool:
        return not self.overlapping_records

    def summary(self) -> str:
        share = len(self.overlapping_records) / len(self.held_out) if self.held_out else 0.0
        return (
            f"{self.source}: {len(self.held_out)} held out, {len(self.trained_on)} trained on;"
            f" {len(self.overlapping_records)} held-out record(s) ({share:.0%}) share an analyte with"
            f" training, over {len(self.overlapping_analytes)} analyte(s)"
        )


def leave_one_source_out(
    records: Iterable[object],
    level: AnalyteLevel = AnalyteLevel.COMPOSITION,
    aliases: Mapping[str, str] | None = None,
) -> tuple[SourceHoldout, ...]:
    """The practical cross-study estimate when the combined scheme collapses.

    One entry per source, each carrying the analyte overlap it could not avoid.
    """
    records = tuple(records)
    if not records:
        raise SplitRefused(NOTHING_TO_SPLIT)
    _gate(records)
    by_source: dict[str, list[int]] = defaultdict(list)
    for index, record in enumerate(records):
        for atom in provenance_atoms(record, aliases) or ("study:unknown",):
            by_source[atom].append(index)
    if len(by_source) < 2:
        raise SplitRefused(f"leave-one-source-out needs at least 2 sources; this corpus has {len(by_source)}")

    keys = [analyte_key(record, level) for record in records]
    everything = set(range(len(records)))
    holdouts: list[SourceHoldout] = []
    for source, held in sorted(by_source.items()):
        rest = sorted(everything - set(held))
        trained_keys = {keys[index] for index in rest}
        overlapping = tuple(sorted(index for index in held if keys[index] in trained_keys))
        holdouts.append(
            SourceHoldout(
                source=source,
                held_out=tuple(sorted(held)),
                trained_on=tuple(rest),
                overlapping_records=overlapping,
                overlapping_analytes=tuple(sorted({keys[index] for index in overlapping})),
            )
        )
    return tuple(holdouts)
