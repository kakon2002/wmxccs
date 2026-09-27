"""How well a structure is attested in the reference corpus, and what that is worth.

This is the ONLY component of the ranker's score that can discriminate between
candidates, and it is worth being exact about why.

The curated biosynthetic rules cannot score. ELEVEN of the fifteen reach N-glycan
enumeration - nine MGAT rules on branching order and bisecting interference, FUT8 on
core fucosylation, and one class-agnostic blood group rule - and the other four are
O-glycan core rules. (Corrected 27 September 2026 from "all fifteen govern MGAT
branching order and bisecting interference".) The enumerator has already rejected every
candidate that breaks one of the eleven, so what is constant across a candidate set is the
number of rules VIOLATED, which is zero. "Rules satisfied" is a different
quantity and it measurably varies: for Hex5HexNAc4Fuc1 the number of rules that
actually bear on a candidate is 1, 2 or 3 across 54, 83 and 30 candidates, and
for Hex5HexNAc2 it is zero for all six. None of that is plausibility - the count
of applicable rules is largely a count of how branched a candidate is, and no
curated rule says a more branched structure is more likely. So the rules are
reported, in three separately named quantities, and none of them enters a score.

So attestation carries the whole discriminating load, and it is thin:

- ONLY FULLY RESOLVED REFERENCE STRUCTURES ATTEST, and resolvedness is DERIVED
  FROM THE PARSE, never read off the record's own flags. `GlycanStructure`
  validates only that SOME structure identifier exists; it never compares the
  flags to the string, so a record can declare itself resolved while its string
  says `Man(a2-3/6)...GlcNAc(?1-?)GlcNAc`. The SugarBase loader happens to set
  the flags from the graph, so a test driven by the loader would pass while a
  fixture or a second ingest path walked straight through. A record whose
  declared flags disagree with its parsed string is REFUSED BY NAME rather than
  silently preferred either way.

- THE UNIT IS THE DISTINCT STRUCTURE, NOT THE ROW. Measured here: 4,001 fully
  resolved records index to 3,640 distinct canonical keys; 352 keys carry more
  than one row (343 twice, 9 three times). Not one of those collisions is two
  identical IUPAC strings - every one is the same structure written with its
  branches in a different order, and several carry the SAME GlyTouCan accession
  on both rows, which is the source's own statement that they are one structure.
  Counting rows would have made every multi-row candidate outrank every other,
  so the entire top of a ranking would have been decided by how SugarBase
  happens to spell things. `attests()` and `structures_for()` are the authority;
  `rows_for()` exists so the discrepancy stays visible rather than being
  quietly corrected away.

- A COUNT HERE IS A COUNT OF DEPOSITED STRUCTURES, NOT AN ABUNDANCE. SugarBase
  carries no abundance, no tissue, and for most rows no species. Nothing here
  supports a statement about how common a structure is in a sample.

- ABSENCE IS NOT EVIDENCE OF ABSENCE. For Hex5HexNAc4Fuc1, 150 of 167
  candidates have no attestation at all.

- THE CANDIDATE SET IS NOT KNOWN TO CONTAIN THE TRUTH, and this is why
  `structures_of_composition` exists. For Hex5HexNAc4Fuc1, 17 of the 34
  fully-resolved reference structures of that composition are not among the 167
  candidates; for Hex5HexNAc2 it is 22 of 28. Any number normalised over the
  candidate set alone would assert the truth is in the set, which is measured
  false for a large and highly variable fraction of compositions: 0 to 100 per cent across those
  holding at least three reference structures, so no range is quoted.

WHY AN EMPTY INDEX IS REFUSED RATHER THAN RETURNED

An index that loaded nothing would give every candidate a count of zero, and a
count of zero is indistinguishable from "this candidate is unattested". Every
candidate would come out equally weighted, the ranker would return one band, and
the output would look exactly like the honest answer for a composition nothing
attests. So the index refuses to exist empty, and it carries the counts it was
built from so that a caller can tell a thin corpus from a broken load.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from typing import Mapping

from .composition import Composition
from .glycan_graph import GlycanGraph, GlycanGraphError
from .sugarbase import DatasetVersion, LoadReport, load_sugarbase

# Refusal wording is an interface: written once, imported by tests rather than matched as prose.
NOTHING_INDEXED = (
    "the attestation index would hold no structure at all. Every candidate would score zero,"
    " which is indistinguishable from every candidate being unattested, so the ranker would"
    " return a plausible-looking flat result built on nothing. Refused rather than returned:"
    " check that glycowork is installed and that load_sugarbase() is returning records"
)
FLAGS_DISAGREE = (
    "reference record {name!r} declares has_unresolved_linkage={declared_linkage} and"
    " has_unresolved_anomericity={declared_anomer}, but its own structure string parses as"
    " linkage={parsed_linkage} and anomericity={parsed_anomer}: {structure!r}. Resolvedness"
    " decides whether this record may attest a specific linkage assignment, so a record that"
    " disagrees with itself about it is refused rather than resolved in either direction."
    " GlycanStructure does not check this - it only requires that some structure identifier"
    " exists - so the disagreement has to be caught here"
)


@dataclass(frozen=True)
class AttestationIndex:
    """Canonical structure keys the reference corpus attests, and how often.

    Keys are `GlycanGraph.canonical_key()`, so a candidate and a reference structure are
    compared through the same canonicalisation rather than by string equality of two IUPAC
    spellings. The candidate side is already deduplicated on that key by the enumerator;
    this index deduplicates the reference side on the same key, which is the symmetry the
    first version of this module did not have.
    """

    dataset_version: DatasetVersion
    # canonical key -> how many ROWS carry it. Kept as rows because the row/structure
    # discrepancy is a measured property of the release worth reporting; `attests` and
    # `structures_for` are what a score may use.
    rows: Mapping[str, int]
    # GlyTouCan accessions per key, so a ranked candidate can cite what attests it. Not every
    # record has one (8,809 of 12,664 do), so an empty tuple here does NOT mean unattested.
    accessions: Mapping[str, tuple[str, ...]]
    # canonical composition string -> the distinct canonical keys of that composition. This is
    # what makes "the truth may not be in the candidate set" a measured number rather than a
    # caveat in prose.
    keys_by_composition: Mapping[str, frozenset[str]]
    rows_by_composition: Mapping[str, int]
    structures_considered: int
    structures_indexed: int
    structures_unresolved: int
    structures_unparsed: int = 0

    def __post_init__(self) -> None:
        if self.structures_indexed <= 0 or not self.rows:
            raise ValueError(NOTHING_INDEXED)

    # --- what a score may use: the distinct structure ------------------------------------------

    def attests(self, canonical_key: str) -> bool:
        return canonical_key in self.rows

    def structures_for(self, canonical_key: str) -> int:
        """0 or 1. The unit of evidence: one distinct reference structure, however spelled."""
        return 1 if canonical_key in self.rows else 0

    # --- measured facts about the release, reported and never scored ---------------------------

    def rows_for(self, canonical_key: str) -> int:
        """How many reference ROWS carry this structure. 2 or 3 means a re-spelling, not evidence."""
        return self.rows.get(canonical_key, 0)

    def accessions_for(self, canonical_key: str) -> tuple[str, ...]:
        return self.accessions.get(canonical_key, ())

    def structures_of_composition(self, composition: Composition | str) -> frozenset[str]:
        """Every distinct fully-resolved reference structure of that composition.

        The denominator for coverage. A candidate set that misses members of this set is a
        candidate set that may not contain the answer.
        """
        text = composition if isinstance(composition, str) else composition.canonical
        return self.keys_by_composition.get(text, frozenset())

    def rows_of_composition(self, composition: Composition | str) -> int:
        text = composition if isinstance(composition, str) else composition.canonical
        return self.rows_by_composition.get(text, 0)

    # --- the shape of the index itself ----------------------------------------------------------

    @property
    def distinct_structures(self) -> int:
        return len(self.rows)

    @property
    def duplicate_rows(self) -> int:
        """Rows beyond the first for a structure. Re-spellings, measured: 361 in SugarBase v12."""
        return self.structures_indexed - self.distinct_structures

    @property
    def keys_with_more_than_one_row(self) -> int:
        return sum(1 for count in self.rows.values() if count > 1)

    @property
    def provenance(self) -> str:
        """What a ranked candidate cites as the source of its attestation."""
        return self.dataset_version.provenance

    def summary(self) -> str:
        return (
            f"{self.dataset_version}  [{self.dataset_version.licence},"
            f" {self.dataset_version.attribution}]\n"
            f"  records built              {self.structures_considered}\n"
            f"  fully resolved, indexed   {self.structures_indexed}"
            "   (resolvedness derived from the parse, not from the record's flags)\n"
            f"  carrying an unknown       {self.structures_unresolved}"
            "   (cannot attest a linkage assignment)\n"
            f"  resolved but unparsed     {self.structures_unparsed}\n"
            f"  DISTINCT structures       {self.distinct_structures}\n"
            f"  duplicate rows            {self.duplicate_rows}"
            f"   over {self.keys_with_more_than_one_row} structures, all of them re-spellings\n"
            f"  distinct compositions     {len(self.keys_by_composition)}"
        )


def build_attestation_index(report: LoadReport | None = None) -> AttestationIndex:
    """Index the fully-resolved reference structures by canonical structure key.

    `report` is taken rather than loaded when given, so a test can index a handful of
    structures without reading the whole dataset.
    """
    report = load_sugarbase() if report is None else report
    rows: Counter[str] = Counter()
    accessions: dict[str, list[str]] = {}
    keys_by_composition: dict[str, set[str]] = {}
    rows_by_composition: Counter[str] = Counter()
    unresolved = 0
    unparsed = 0

    for record in report.structures:
        text = record.iupac_condensed
        if not text:
            # A record can be fully resolved on the strength of a WURCS or an accession rather
            # than an IUPAC string. It states a structure this module cannot key, which is not
            # the same as stating none, so it is counted separately and never as unresolved.
            if not (record.has_unresolved_linkage or record.has_unresolved_anomericity):
                unparsed += 1
            else:
                unresolved += 1
            continue
        try:
            graph = GlycanGraph.from_iupac_condensed(text)
        except GlycanGraphError:
            unparsed += 1
            continue

        # RESOLVEDNESS COMES FROM THE PARSE. The record's own flags are compared against it and
        # a disagreement stops the load; see the module docstring.
        parsed_linkage = graph.has_unresolved_linkage
        parsed_anomer = graph.has_unresolved_anomericity
        if (record.has_unresolved_linkage, record.has_unresolved_anomericity) != (
            parsed_linkage,
            parsed_anomer,
        ):
            raise ValueError(
                FLAGS_DISAGREE.format(
                    name=record.glytoucan_ac or record.source,
                    declared_linkage=record.has_unresolved_linkage,
                    declared_anomer=record.has_unresolved_anomericity,
                    parsed_linkage=parsed_linkage,
                    parsed_anomer=parsed_anomer,
                    structure=text,
                )
            )
        if parsed_linkage or parsed_anomer:
            unresolved += 1
            continue

        key = graph.canonical_key()
        rows[key] += 1
        if record.glytoucan_ac:
            accessions.setdefault(key, []).append(record.glytoucan_ac)
        composition = graph.composition().canonical
        keys_by_composition.setdefault(composition, set()).add(key)
        rows_by_composition[composition] += 1

    return AttestationIndex(
        dataset_version=report.dataset_version,
        rows=dict(rows),
        accessions={key: tuple(found) for key, found in accessions.items()},
        keys_by_composition={
            composition: frozenset(keys) for composition, keys in keys_by_composition.items()
        },
        rows_by_composition=dict(rows_by_composition),
        structures_considered=len(report.structures),
        structures_indexed=sum(rows.values()),
        structures_unresolved=unresolved,
        structures_unparsed=unparsed,
    )


@lru_cache(maxsize=1)
def default_attestation_index() -> AttestationIndex:
    """The index over the installed glycowork release, built once per process.

    Cached because building it reads 50,461 rows and re-parses 4,001 structures, which is
    about thirteen seconds. The result is frozen and its mappings are never handed out
    mutable, so the cache is shared safely.
    """
    return build_attestation_index()
