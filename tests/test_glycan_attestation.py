"""The attestation index: what may count as evidence, and what a zero must not hide.

The index is the only discriminating input the ranker has, so almost every way it could be
silently wrong is a way the ranker silently stops meaning anything. These tests are written
against the real SugarBase release rather than fixtures wherever the figure is the point,
because a fixture-only suite here would be the classic case of the wrong answer being the
right one: a canonicalisation that matched nothing would score every candidate zero, which is
exactly what an honest thin corpus looks like.
"""

from __future__ import annotations

import pytest

from wmxglycan.attestation import (
    FLAGS_DISAGREE,
    NOTHING_INDEXED,
    AttestationIndex,
    build_attestation_index,
    default_attestation_index,
)
from wmxglycan.composition import parse_composition
from wmxglycan.enumeration import Enumerator
from wmxglycan.glycan_graph import GlycanGraph
from wmxglycan.models import GlycanStructure
from wmxglycan.reuse import ReuseStatus
from wmxglycan.sugarbase import LoadReport, dataset_version

# Measured against glycowork 1.10.0 in this repository. Pinned, because every one of these is a
# figure a later change could move without anything else noticing.
INDEXED = 4001
DISTINCT = 3640
UNRESOLVED = 8663
DUPLICATE_ROWS = INDEXED - DISTINCT
KEYS_WITH_MORE_THAN_ONE_ROW = 352

# composition -> (candidates, reference rows, distinct reference structures, attested, not enumerated)
MEASURED = {
    "Hex3HexNAc4Fuc1": (10, 26, 22, 6, 16),
    "Hex4HexNAc4Fuc1": (63, 27, 26, 16, 10),
    "Hex5HexNAc4Fuc1": (167, 38, 34, 17, 17),
    "Hex5HexNAc2": (6, 30, 28, 6, 22),
}


@pytest.fixture(scope="module")
def index():
    return default_attestation_index()


@pytest.fixture(scope="module")
def enumerator():
    return Enumerator()


# --- the index refuses to exist empty --------------------------------------------------------


def test_an_empty_index_is_refused_rather_than_returned():
    # A zero count and "the index loaded nothing" are the same number downstream, and the second
    # produces output indistinguishable from the honest answer for an unattested composition.
    with pytest.raises(ValueError, match="no structure at all"):
        AttestationIndex(
            dataset_version=dataset_version(),
            rows={},
            accessions={},
            keys_by_composition={},
            rows_by_composition={},
            structures_considered=0,
            structures_indexed=0,
            structures_unresolved=0,
        )


def test_the_refusal_wording_is_importable_rather_than_matched_as_prose():
    assert "Refused rather than returned" in NOTHING_INDEXED


def test_a_report_with_no_usable_structure_is_refused():
    report = LoadReport(dataset_version=dataset_version(), structures=())
    with pytest.raises(ValueError, match="no structure at all"):
        build_attestation_index(report)


# --- resolvedness is derived from the parse, never self-certified ------------------------------


def _record(**overrides) -> GlycanStructure:
    fields = dict(
        composition="Hex5HexNAc2",
        iupac_condensed="Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
        has_unresolved_linkage=False,
        has_unresolved_anomericity=False,
        source="a test",
        reuse_status=ReuseStatus.SYNTHETIC_FIXTURE,
    )
    fields.update(overrides)
    return GlycanStructure(**fields)


def test_a_record_that_declares_itself_resolved_while_its_string_is_not_is_refused_by_name():
    # GlycanStructure validates only that SOME structure identifier exists; it never compares
    # the flags to the string. The SugarBase loader happens to set them from the graph, so a
    # loader-driven test passes and a fixture or a second ingest path walks straight through.
    ambiguous = _record(
        composition="Hex3HexNAc2",
        iupac_condensed="Man(a2-3/6)[Man(a1-6)]Man(b1-4)GlcNAc(?1-?)GlcNAc",
        has_unresolved_linkage=False,
        has_unresolved_anomericity=False,
    )
    graph = GlycanGraph.from_iupac_condensed(ambiguous.iupac_condensed)
    assert graph.has_unresolved_linkage or graph.has_unresolved_anomericity, (
        "the fixture must actually be ambiguous, or this test proves nothing"
    )
    report = LoadReport(dataset_version=dataset_version(), structures=(ambiguous,))
    with pytest.raises(ValueError) as caught:
        build_attestation_index(report)
    assert "declares has_unresolved_linkage" in str(caught.value)
    assert "refused rather than resolved in either direction" in str(caught.value)


def test_the_disagreement_wording_names_both_the_declaration_and_the_parse():
    assert "{declared_linkage}" in FLAGS_DISAGREE and "{parsed_linkage}" in FLAGS_DISAGREE


def test_a_record_that_declares_itself_unresolved_and_is_unresolved_is_simply_skipped():
    honest = _record(
        composition="Hex3HexNAc2",
        iupac_condensed="Man(a2-3/6)[Man(a1-6)]Man(b1-4)GlcNAc(?1-?)GlcNAc",
        has_unresolved_linkage=True,
        has_unresolved_anomericity=True,
    )
    resolved = _record()
    report = LoadReport(dataset_version=dataset_version(), structures=(honest, resolved))
    built = build_attestation_index(report)
    assert built.structures_indexed == 1
    assert built.structures_unresolved == 1


# --- the unit is the distinct structure, not the row -------------------------------------------


def test_the_real_release_holds_duplicate_rows_for_one_structure(index):
    # THE FINDING THIS MODULE WAS REWRITTEN FOR. Counting rows would have decided the top of
    # every ranking by how SugarBase happens to spell a structure.
    assert index.structures_indexed == INDEXED
    assert index.distinct_structures == DISTINCT
    assert index.duplicate_rows == DUPLICATE_ROWS == 361
    assert index.keys_with_more_than_one_row == KEYS_WITH_MORE_THAN_ONE_ROW


def test_a_structure_spelled_twice_attests_once(index):
    # structures_for is 0 or 1 for every key in the corpus, however many rows carry it.
    multi = [key for key, rows in index.rows.items() if rows > 1]
    assert multi, "the corpus must hold at least one re-spelled structure or this proves nothing"
    for key in multi[:20]:
        assert index.rows_for(key) > 1
        assert index.structures_for(key) == 1
        assert index.attests(key)


def test_structures_for_is_never_more_than_one(index):
    assert set(index.structures_for(key) for key in list(index.rows)[:500]) == {1}
    assert index.structures_for("NOT A REAL KEY") == 0


def test_two_rows_of_one_structure_really_are_the_same_structure_written_differently(index):
    # The claim behind the dedup: collisions are re-spellings, not two different molecules.
    # Both sides of the collision canonicalise to one key, which is the definition, so what is
    # asserted here is that the row count exceeds the structure count for those keys.
    doubled = [key for key, rows in index.rows.items() if rows == 2]
    tripled = [key for key, rows in index.rows.items() if rows == 3]
    assert len(doubled) == 343
    assert len(tripled) == 9
    assert len(doubled) + len(tripled) == KEYS_WITH_MORE_THAN_ONE_ROW


# --- the positive control: the canonicalisation really matches something -----------------------


@pytest.mark.parametrize("composition", sorted(MEASURED))
def test_the_index_matches_the_measured_number_of_candidates(index, enumerator, composition):
    # THE MOST IMPORTANT TEST IN THIS FILE. A canonicalisation that silently stopped matching -
    # a changed residue name, a changed linkage spelling, a reversed branch sort - would score
    # every candidate zero, and zero is what an honestly unattested composition looks like. So
    # the attested count is pinned per composition, not merely asserted to be non-zero.
    candidates, rows, distinct, attested, not_enumerated = MEASURED[composition]
    keys = {c.canonical_key for c in enumerator.enumerate(composition).candidates}
    assert len(keys) == candidates
    references = index.structures_of_composition(composition)
    assert index.rows_of_composition(composition) == rows
    assert len(references) == distinct
    assert len(keys & references) == attested
    assert len(references - keys) == not_enumerated


def test_the_off_support_accounting_balances(index, enumerator):
    # attested + not-enumerated must equal the distinct structures of the composition. A
    # coverage figure computed from the wrong denominator still looks plausible.
    for composition, (_c, _rows, distinct, attested, not_enumerated) in MEASURED.items():
        keys = {c.canonical_key for c in enumerator.enumerate(composition).candidates}
        references = index.structures_of_composition(composition)
        assert len(keys & references) + len(references - keys) == distinct
        assert attested + not_enumerated == distinct


def test_every_candidates_stored_key_round_trips_against_its_own_structure(enumerator):
    # The field the whole match depends on. A key drifted from its structure scores zero
    # attestation and is indistinguishable from a genuinely unattested candidate.
    for composition in sorted(MEASURED):
        for candidate in enumerator.enumerate(composition).candidates:
            graph = GlycanGraph.from_iupac_condensed(candidate.iupac_condensed)
            assert graph.canonical_key() == candidate.canonical_key


def test_an_unknown_composition_reports_an_empty_reference_set_rather_than_raising(index):
    # And the caller must be able to tell that apart from "no structures matched", which the
    # ranker does via CoverageState.UNEVALUABLE.
    assert index.structures_of_composition(parse_composition("Hex2HexNAc7NeuGc4")) == frozenset()
    assert index.rows_of_composition("Hex2HexNAc7NeuGc4") == 0


# --- provenance travels -----------------------------------------------------------------------


def test_the_index_carries_a_live_dataset_version_not_a_literal(index):
    from importlib.metadata import version

    assert index.dataset_version.package_version == version("glycowork")
    assert index.dataset_version.licence == "MIT"
    assert "Daniel Bojar" in index.dataset_version.attribution
    assert "SugarBase v12" in index.provenance


def test_the_unresolved_records_are_counted_rather_than_dropped(index):
    assert index.structures_unresolved == UNRESOLVED
    assert index.structures_considered == index.structures_indexed + UNRESOLVED + index.structures_unparsed


def test_the_summary_says_resolvedness_is_derived(index):
    assert "derived from the parse" in index.summary()
    assert "re-spellings" in index.summary()
