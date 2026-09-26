"""Grouping records into folds, and the leaks the grouping exists to prevent.

Every measurement here is a synthetic fixture. No CCS value is real.

The tests that matter most are the ones asserting that two records a model
cannot tell apart stay together. Several of them assert first that the fixture
really does contain two different molecules, because a leakage test that passes
on two identical inputs proves nothing.
"""

import inspect
import os
import subprocess
import sys
from fractions import Fraction

import pytest

from wmxglycan.features import STRUCTURE_PARSE_FAILED, extract, extract_all
from wmxglycan.glycan_graph import GlycanGraph
from wmxglycan.licensing import LicenceGateError, TrainingGateError
from wmxglycan.models import CCSMeasurement, GlycanStructure
from wmxglycan.splits import (
    FOLD_TOLERANCE,
    MIN_ANALYTE_GROUPS,
    MIN_GROUPS_PER_FOLD,
    N_FOLDS,
    NOTHING_TO_SPLIT,
    AnalyteLevel,
    Grouping,
    LeakageError,
    SourceHoldout,
    SplitError,
    SplitMode,
    SplitRefused,
    analyte_key,
    assert_features_do_not_straddle,
    assert_no_group_straddles,
    group_records,
    grouped_folds,
    identity_atoms,
    leave_one_source_out,
    provenance_atoms,
)

FIXTURE_SOURCE = "synthetic test fixture, not a real record"
RESOLVED = {"has_unresolved_linkage": False, "has_unresolved_anomericity": False}

G1F_A3 = "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
G1F_A6 = "GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc"
SIA_A2_3 = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
SIA_A2_6 = "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
SIA_AMBIGUOUS = "Neu5Ac(a2-3/6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
SIA_COMPOSITION = "Hex4HexNAc3NeuAc1"


def a_glycan(composition="Hex5HexNAc4Fuc1", **overrides):
    fields = {"composition": composition, "source": FIXTURE_SOURCE, "reuse_status": "synthetic_fixture"}
    return GlycanStructure(**(fields | overrides))


def a_measurement(**overrides):
    fields = {
        "glycan": a_glycan(),
        "ccs": 100.0,
        "adduct": "[M+2H]2+",
        "charge": 2,
        "polarity": "positive",
        "reducing_end_label": "native",
        "derivatisation": "underivatised",
        "ims_type": "TWIMS",
        "drift_gas": "N2",
        "calibrant": "fixture calibrant",
        "source": FIXTURE_SOURCE,
        "reuse_status": "synthetic_fixture",
    }
    return CCSMeasurement(**(fields | overrides))


def structured(iupac, composition, **overrides):
    return a_measurement(glycan=a_glycan(composition, iupac_condensed=iupac, **RESOLVED), **overrides)


def many(count=MIN_ANALYTE_GROUPS, per_group=1):
    """Enough distinct compositions to make a legal split, one record each."""
    records = []
    for index in range(count):
        composition = f"Hex{3 + index // 4}HexNAc{2 + index % 4}"
        for repeat in range(per_group):
            records.append(a_measurement(glycan=a_glycan(composition), ccs=100.0 + index + repeat))
    return records


def key_of(iupac):
    return GlycanGraph.from_iupac_condensed(iupac).canonical_key()


def group_ids(records, **kwargs):
    grouping = group_records(records, **kwargs)
    return set(grouping.group_of), grouping


# --- the law: what a model cannot tell apart stays together --------------------


def test_records_the_model_cannot_tell_apart_share_a_group_and_a_fold():
    records = many(per_group=2)
    report = grouped_folds(records)
    rows = [extract(record).analyte_row() for record in records]
    buckets = {}
    for index, row in enumerate(rows):
        buckets.setdefault(row, []).append(index)
    collisions = [members for members in buckets.values() if len(members) > 1]
    assert collisions, "the fixture must contain at least one feature-identical pair"
    for members in collisions:
        assert len({report.group_of[i] for i in members}) == 1
        assert len({report.assignment[i] for i in members}) == 1


def test_two_isomers_of_one_composition_land_in_one_group():
    # Assert first that these really are two different structures, or the test
    # would pass on two identical inputs.
    assert key_of(SIA_A2_3) != key_of(SIA_A2_6)
    ids, _ = group_ids([structured(SIA_A2_3, SIA_COMPOSITION), structured(SIA_A2_6, SIA_COMPOSITION)])
    assert len(ids) == 1


def test_the_two_arm_isomers_are_never_separated():
    # The isomer axis this platform exists for: same composition, different arm.
    a3, a6 = structured(G1F_A3, "Hex4HexNAc4Fuc1"), structured(G1F_A6, "Hex4HexNAc4Fuc1")
    assert key_of(G1F_A3) != key_of(G1F_A6)
    assert extract(a3).as_mapping()["arm_a3_size"] != extract(a6).as_mapping()["arm_a3_size"]
    ids, _ = group_ids([a3, a6])
    assert len(ids) == 1


def test_a_structure_keyed_split_would_have_separated_them():
    # Encodes the refutation, so a later "simplification" to the obvious
    # structure key fails with a message that explains itself.
    records = [structured(G1F_A3, "Hex4HexNAc4Fuc1"), structured(G1F_A6, "Hex4HexNAc4Fuc1")]
    structure_keyed = {key_of(r.glycan.iupac_condensed) for r in records}
    assert len(structure_keyed) == 2, "a structure key splits these two"
    ids, _ = group_ids(records)
    assert len(ids) == 1, "the analyte key keeps them together"


def test_the_same_glycan_from_two_papers_is_one_group():
    a = a_measurement(source="10.1000/paper-a")
    b = a_measurement(source="10.1000/paper-b")
    assert a != b and len({a, b}) == 2  # equality includes provenance, so a set dedupes nothing
    ids, _ = group_ids([a, b])
    assert len(ids) == 1


def test_a_compound_key_of_analyte_and_source_would_split_them():
    # The join, not the meet. A tuple key is strictly finer than either axis and
    # so prevents neither leak while looking rigorous.
    a = a_measurement(source="10.1000/paper-p")
    b = a_measurement(source="10.1000/paper-q")
    c = a_measurement(glycan=a_glycan("Hex3HexNAc4Fuc1"), source="10.1000/paper-q")
    tuple_keyed = {(analyte_key(r), r.source) for r in (a, b, c)}
    assert len(tuple_keyed) == 3
    ids, _ = group_ids([a, b, c], by_study=True)
    assert len(ids) == 1


def test_the_same_glycan_under_two_adducts_is_not_split():
    a = a_measurement()
    b = a_measurement(adduct="[M+H+Na]2+")
    assert a.calibration_group != b.calibration_group  # genuinely different conditions
    ids, _ = group_ids([a, b])
    assert len(ids) == 1


def test_the_leakage_guard_reads_the_analyte_row_not_the_whole_row():
    records = [a_measurement(), a_measurement(adduct="[M+H+Na]2+")]
    grouping = group_records(records)
    rows = [extract(record).analyte_row() for record in records]
    assert rows[0] == rows[1]  # one analyte
    assert extract(records[0]).as_row() != extract(records[1]).as_row()  # two conditions
    assert_features_do_not_straddle(grouping, rows)  # must not raise


def test_the_leakage_guard_bites_when_a_group_is_made_too_fine():
    records = [a_measurement(), a_measurement(ccs=222.0)]
    grouping = group_records(records)
    broken = Grouping(
        level=grouping.level,
        by_study=grouping.by_study,
        records=grouping.records,
        group_of=("one", "two"),  # hand-corrupted: the same analyte in two groups
        groups={"one": (0,), "two": (1,)},
        sources=grouping.sources,
        bridges=grouping.bridges,
        identity_conflicts=grouping.identity_conflicts,
    )
    rows = [extract(record).analyte_row() for record in records]
    with pytest.raises(LeakageError, match="cannot tell them apart"):
        assert_features_do_not_straddle(broken, rows)


def test_an_unresolved_linkage_is_not_resolved_to_make_a_key():
    ambiguous = structured(SIA_AMBIGUOUS, SIA_COMPOSITION)
    resolved = structured(SIA_A2_6, SIA_COMPOSITION)
    # Nothing was completed: the key still says 3/6, and differs from the resolved one.
    assert "3/6" in key_of(SIA_AMBIGUOUS)
    assert key_of(SIA_AMBIGUOUS) != key_of(SIA_A2_6)
    # And they are nonetheless one group, through the shared analyte atom.
    ids, _ = group_ids([ambiguous, resolved])
    assert len(ids) == 1


# --- the key reads the analyte and nothing else --------------------------------


def test_the_key_ignores_the_value_the_provenance_and_the_licence():
    variants = [
        a_measurement(),
        a_measurement(ccs=999.0),
        a_measurement(source="another fixture source"),
        a_measurement(instrument="a fixture instrument"),
        a_measurement(cell_gas="He"),
        a_measurement(adduct="[M+H]+", charge=1),
        a_measurement(ims_type="TIMS"),
        a_measurement(reducing_end_label="PA"),
    ]
    assert len({analyte_key(record) for record in variants}) == 1


@pytest.mark.parametrize("composition", ["Hex4HexNAc4Fuc1", "Hex5HexNAc5Fuc1", "Hex5HexNAc4", "Hex5HexNAc4Fuc1NeuAc1"])
def test_each_residue_count_splits_the_key(composition):
    assert analyte_key(a_measurement(glycan=a_glycan(composition))) != analyte_key(a_measurement())


def test_the_calibration_group_and_the_instrument_are_not_grouping_atoms():
    a = a_measurement(glycan=a_glycan("Hex3HexNAc2"), instrument="one machine")
    b = a_measurement(glycan=a_glycan("Hex9HexNAc2"), instrument="one machine")
    assert a.calibration_group == b.calibration_group  # they share everything but the analyte
    ids, _ = group_ids([a, b])
    assert len(ids) == 2


def test_the_glycans_own_source_is_not_a_provenance_atom():
    # It records where the structure assignment came from, so counting it would
    # put every record from one database in a single group.
    a = a_measurement(glycan=a_glycan("Hex3HexNAc2"), source="10.1000/paper-a")
    b = a_measurement(glycan=a_glycan("Hex9HexNAc2"), source="10.1000/paper-b")
    assert a.glycan.source == b.glycan.source
    ids, _ = group_ids([a, b], by_study=True)
    assert len(ids) == 2


# --- backbone level -------------------------------------------------------------


def test_backbones_one_residue_apart_are_still_different_groups():
    # The key is a projection, not a transitive closure. Closing over "one
    # residue apart" would put nearly the whole census in one component.
    a = a_measurement(glycan=a_glycan("Hex5HexNAc4"))
    b = a_measurement(glycan=a_glycan("Hex6HexNAc4"))
    ids, _ = group_ids([a, b], level=AnalyteLevel.BACKBONE)
    assert len(ids) == 2


def test_the_fucosylation_series_merges_only_at_backbone_level():
    records = [
        a_measurement(glycan=a_glycan("Hex5HexNAc4")),
        a_measurement(glycan=a_glycan("Hex5HexNAc4Fuc1")),
        a_measurement(glycan=a_glycan("Hex5HexNAc4NeuAc1")),
    ]
    at_composition, _ = group_ids(records)
    at_backbone, _ = group_ids(records, level=AnalyteLevel.BACKBONE)
    assert len(at_composition) == 3  # too coarse a level would fail here
    assert len(at_backbone) == 1  # too fine a level would fail here


# --- identifiers merge, and say when they disagree -------------------------------


def test_an_identifier_shared_across_two_analyte_keys_merges_and_is_reported():
    a = a_measurement(glycan=a_glycan("Hex5HexNAc4Fuc1", glytoucan_ac="G00000AA"))
    b = a_measurement(glycan=a_glycan("Hex3HexNAc2", glytoucan_ac="G00000AA"))
    grouping = group_records([a, b])
    assert len(set(grouping.group_of)) == 1  # merging is the safe direction
    assert grouping.identity_conflicts  # but it is never silent
    # The conflict names both analyte atoms, so a reader can see what disagreed.
    conflict = next(iter(grouping.identity_conflicts.values()))
    assert any("Hex5HexNAc4Fuc1" in atom for atom in conflict)
    assert any("Hex3HexNAc2" in atom for atom in conflict)


def test_an_accession_is_normalised_before_it_is_compared():
    atoms = identity_atoms(a_measurement(glycan=a_glycan(glytoucan_ac="G00000AA")))
    assert any(atom == "glytoucan:G00000AA" for atom in atoms)


def test_a_structure_that_does_not_parse_still_groups_on_its_composition():
    # 'Gal(b1-4' really is accepted by GlycanStructure: validation refuses an
    # unparseable string only when the resolved-linkage flags are also claimed.
    # An unknown residue name is no route here, because 'NotASugar(b1-4)GlcNAc'
    # parses perfectly well.
    broken = a_measurement(glycan=a_glycan("Hex5HexNAc4Fuc1", iupac_condensed="Gal(b1-4"))
    assert identity_atoms(broken) == ()  # it contributes no structure atom
    grouping = group_records([a_measurement(), broken])
    assert grouping.notes and any("did not parse" in note for note in grouping.notes)
    # ...and it is still grouped on its composition rather than dropped.
    assert len(set(grouping.group_of)) == 1
    assert grouping.group_of[1] == f"analyte:{analyte_key(broken)}"


def test_a_structure_that_does_not_parse_is_counted_by_the_featuriser():
    broken = a_measurement(glycan=a_glycan("Hex5HexNAc4Fuc1", iupac_condensed="Gal(b1-4"))
    report = extract_all([broken])
    assert report.rows_in == 1
    assert report.vectors_built == 0
    assert report.failure_counts.get(STRUCTURE_PARSE_FAILED) == 1


# --- the degenerate case, so the machinery earns its place -----------------------


def test_with_the_study_axis_off_the_components_are_exactly_the_analyte_classes():
    records = many()
    grouping = group_records(records)
    for record, group in zip(records, grouping.group_of):
        assert group == f"analyte:{analyte_key(record)}"


# --- folds ------------------------------------------------------------------------


def test_every_record_is_in_exactly_one_group_and_one_fold():
    records = many(per_group=2)
    report = grouped_folds(records)
    assert sum(report.fold_sizes) == len(records) == len(report.assignment)
    assert sum(report.group_sizes.values()) == len(records)
    seen = [index for fold in report.folds for index in fold]
    assert sorted(seen) == list(range(len(records)))


def test_no_group_appears_in_two_folds_and_the_checker_bites_when_one_does():
    records = many(per_group=2)
    report = grouped_folds(records)
    assert_no_group_straddles(report.folds, report.group_of)  # the clean case
    corrupted = [list(fold) for fold in report.folds]
    moved = corrupted[0].pop()
    corrupted[1].append(moved)
    with pytest.raises(LeakageError, match="appears in folds"):
        assert_no_group_straddles(corrupted, report.group_of)


def test_the_folds_do_not_depend_on_the_input_order():
    records = many(per_group=2)
    forward = grouped_folds(records)
    backward = grouped_folds(list(reversed(records)))
    assert sorted(forward.fold_sizes) == sorted(backward.fold_sizes)
    assert {group: size for group, size in forward.group_sizes.items()} == dict(backward.group_sizes)


def test_the_folds_do_not_depend_on_the_hash_seed():
    code = (
        "from tests.test_glycan_splits import many\n"
        "from wmxglycan.splits import grouped_folds\n"
        "print(grouped_folds(many(per_group=2)).fold_sizes)\n"
    )
    outputs = set()
    for seed in ("0", "1", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=seed, PYTHONPATH=os.getcwd())
        result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
        assert result.returncode == 0, result.stderr
        outputs.add(result.stdout.strip())
    assert len(outputs) == 1, outputs


def test_the_signature_takes_no_seed():
    # A lucky split must not be reachable by reseeding.
    for function in (grouped_folds, group_records):
        names = set(inspect.signature(function).parameters)
        assert not {name for name in names if any(word in name for word in ("seed", "random", "shuffle"))}


def test_no_floor_is_a_parameter():
    # A floor that can be passed in is a floor that can be passed a zero.
    for function in (grouped_folds, group_records, leave_one_source_out):
        names = set(inspect.signature(function).parameters)
        assert not {n for n in names if any(w in n for w in ("min", "max", "threshold", "floor", "tolerance"))}


def test_the_floors_are_pinned_to_their_literal_values():
    # Asserted literally, not derived. Every fixture here sizes itself from
    # MIN_ANALYTE_GROUPS, so lowering a floor would make the fixtures follow it
    # and the suite would stay green while accepting five glycans for 5-fold CV.
    assert (N_FOLDS, MIN_GROUPS_PER_FOLD, MIN_ANALYTE_GROUPS) == (5, 5, 25)
    assert MIN_ANALYTE_GROUPS == N_FOLDS * MIN_GROUPS_PER_FOLD
    # An exact fraction: the float 0.05 is not one twentieth, which is the whole
    # reason the bound is rational.
    assert FOLD_TOLERANCE == Fraction(1, 20)
    assert float(FOLD_TOLERANCE) == 0.05


def lopsided():
    """A corpus whose groups differ in size, so the packing order actually matters."""
    records = []
    for index in range(MIN_ANALYTE_GROUPS):
        composition = f"Hex{3 + index // 4}HexNAc{2 + index % 4}"
        for repeat in range(1 + index % 4):  # sizes 1, 2, 3, 4, 1, 2, ...
            records.append(a_measurement(glycan=a_glycan(composition), ccs=100.0 + index + repeat / 10))
    return records


def test_groups_of_differing_sizes_are_still_packed_within_tolerance():
    # Every group in `many()` is the same size, so round-robin packing would pass
    # there. Heterogeneous sizes are what require largest-first.
    report = grouped_folds(lopsided())
    total = len(report.records)
    for size in report.fold_sizes:
        assert Fraction(1, N_FOLDS) - FOLD_TOLERANCE <= Fraction(size, total) <= Fraction(1, N_FOLDS) + FOLD_TOLERANCE
    assert len(set(report.group_sizes.values())) > 1, "the fixture must hold groups of different sizes"


# --- refusals ----------------------------------------------------------------------


def test_the_split_refuses_an_empty_corpus_and_says_so():
    with pytest.raises(SplitRefused, match=NOTHING_TO_SPLIT) as caught:
        grouped_folds([])
    # An absence of data is not a licence problem.
    assert not isinstance(caught.value, LicenceGateError)


def test_too_few_groups_is_refused_rather_than_packed():
    with pytest.raises(SplitRefused, match="too few") as caught:
        grouped_folds([a_measurement(), a_measurement(glycan=a_glycan("Hex3HexNAc2"))])
    assert str(MIN_ANALYTE_GROUPS) in str(caught.value)
    assert not isinstance(caught.value, LicenceGateError)


def test_a_group_holding_most_of_the_data_is_refused_and_named():
    records = many() + [a_measurement()] * 60  # one composition swamping the rest
    with pytest.raises(SplitRefused, match="never divided") as caught:
        grouped_folds(records)
    assert "Hex5HexNAc4Fuc1" in str(caught.value)


def test_the_refusal_cannot_be_swallowed_by_a_row_level_except_clause():
    assert issubclass(SplitError, TrainingGateError)
    assert not issubclass(SplitError, (ValueError, TypeError))
    assert issubclass(SplitRefused, SplitError) and issubclass(LeakageError, SplitError)


def test_the_licence_gate_runs_before_any_grouping_and_reports_every_fault():
    unverified = a_measurement(reuse_status="unverified")
    with pytest.raises(LicenceGateError) as caught:
        group_records([a_measurement(), unverified])
    assert "unverified" in str(caught.value)
    two_bad = a_measurement(reuse_status="academic_only")
    with pytest.raises(LicenceGateError) as both:
        group_records([unverified, two_bad])
    assert "record 0" in str(both.value) and "record 1" in str(both.value)


def test_the_licence_gate_runs_before_a_leave_one_source_out_split_too():
    # group_records has had its gate tested since the beginning; this entry point
    # never did, and a mutation removing the call survived every sweep unnoticed
    # because one anchor matched both call sites and only ever mutated the first.
    # An ungated record reaching a held-out split is the same fault wherever it
    # gets in, so the gate has to run at both doors.
    ungated = a_measurement(reuse_status="unverified", source="10.1000/paper-b")
    with pytest.raises(LicenceGateError, match="unverified"):
        leave_one_source_out([a_measurement(), ungated])


def test_an_unbacked_claim_keeps_its_class_through_the_splitter():
    # The splitter's gate is the other admission point, and it must not report an
    # unbacked claim as a plain licence refusal: that tells the reader to relabel
    # a record whose licence is fine, instead of to record the licence.
    from wmxglycan.licensing import UnbackedClaimError

    unbacked = a_measurement(reuse_status="open_attribution")  # no DOI, so nothing backs it
    with pytest.raises(UnbackedClaimError):
        group_records([unbacked])
    # A genuine licence refusal alongside it is the more serious fault and wins.
    with pytest.raises(LicenceGateError) as caught:
        group_records([unbacked, a_measurement(reuse_status="academic_only")])
    assert not isinstance(caught.value, UnbackedClaimError)


def test_a_cleared_measurement_wrapping_a_restricted_structure_is_refused():
    record = a_measurement(glycan=a_glycan(reuse_status="academic_only"))
    with pytest.raises(LicenceGateError, match="academic_only"):
        group_records([record])


# --- the study axis is always attempted, and its refusal is reported ---------------


def test_a_single_source_corpus_reports_that_it_cannot_support_a_cross_study_split():
    records = [a_measurement(source="10.1000/only-paper") for _ in range(3)]
    grouping = group_records(records, by_study=True)
    assert len(set(grouping.group_of)) == 1
    with pytest.raises(SplitRefused):
        leave_one_source_out(records)


def test_leave_one_source_out_holds_out_whole_sources():
    records = [
        a_measurement(source="10.1000/paper-a"),
        a_measurement(glycan=a_glycan("Hex3HexNAc2"), source="10.1000/paper-b"),
    ]
    splits = leave_one_source_out(records)
    assert len(splits) == 2
    for holdout in splits:
        assert isinstance(holdout, SourceHoldout)
        assert set(holdout.held_out) & set(holdout.trained_on) == set()
        assert len(holdout.held_out) + len(holdout.trained_on) == len(records)
    # Different analytes in each paper, so nothing is contaminated here.
    assert all(holdout.clean for holdout in splits)


def test_leave_one_source_out_reports_the_overlap_it_cannot_avoid():
    # The same glycan measured by two papers. Refusing on this would make the
    # function useless on every real CCS corpus, because every paper measures the
    # common N-glycans; returning bare folds and saying nothing would be the
    # dishonest option. So the contamination travels with the split.
    records = [
        a_measurement(source="10.1000/paper-a"),
        a_measurement(source="10.1000/paper-b", ccs=222.0),
    ]
    splits = leave_one_source_out(records)
    assert len(splits) == 2
    for holdout in splits:
        assert not holdout.clean
        assert holdout.overlapping_records == holdout.held_out  # every held-out record has a twin
        assert holdout.overlapping_analytes == ("Hex5HexNAc4Fuc1",)
        assert "share an analyte with training" in holdout.summary()


def test_a_recorded_doi_is_preferred_over_one_scraped_from_prose():
    # A stated DOI is a fact; scraping one out of free text is the fallback that
    # used to split a single paper into several sources on trailing punctuation.
    stated = a_measurement(source="Smith et al., see the supplement", doi="10.1000/paper-a")
    scraped = a_measurement(source="Smith et al. (10.1000/paper-a)")
    assert provenance_atoms(stated) == provenance_atoms(scraped) == ("doi:10.1000/paper-a",)


def test_a_recorded_doi_wins_when_the_prose_names_a_different_one():
    # A methods paper cited in the source text must not become the study identity.
    record = a_measurement(source="method of 10.1000/methods-paper", doi="10.1000/the-actual-study")
    assert provenance_atoms(record) == ("doi:10.1000/the-actual-study",)


# --- the two split modes -------------------------------------------------------


def arm_isomers(length):
    """One composition, two structures: the same chain on the a1-3 arm or the a1-6 arm.

    Chosen because the feature set can actually see this difference, through
    arm_a3_size and arm_a6_size. Not every pair of isomers is separable by these
    features - see the blind-spot test below - so a fixture built from an
    invisible pair would be asking the split to do something the model cannot.
    """
    chain = "Man(a1-2)" * length
    on_a3 = f"{chain}Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
    on_a6 = f"Man(a1-3)[{chain}Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc"
    return f"Hex{3 + length}HexNAc2", on_a3, on_a6


def linear(hexes, hexnacs, arm="a1-2"):
    """A linear chain of exactly that composition. `arm` changes the linkage only."""
    return "Man(%s)" % arm * (hexes - 1) + "Man(b1-4)" + "GlcNAc(b1-4)" * (hexnacs - 1) + "GlcNAc"


def isomer_pairs(compositions=13):
    """Two feature-distinguishable structures for each of several compositions."""
    records = []
    for index in range(1, compositions + 1):
        composition, on_a3, on_a6 = arm_isomers(index)
        for iupac, ccs in ((on_a3, 100.0), (on_a6, 200.0)):
            records.append(
                a_measurement(
                    glycan=a_glycan(composition, iupac_condensed=iupac, **RESOLVED), ccs=ccs + index
                )
            )
    return records


def test_isomers_of_one_composition_are_separated_only_in_isomer_mode():
    a, b = isomer_pairs(compositions=1)
    assert key_of(a.glycan.iupac_condensed) != key_of(b.glycan.iupac_condensed)  # genuinely two molecules
    assert extract(a).as_row() != extract(b).as_row()  # and the features can see the difference
    # Deployment mode keeps them together: the deployment question is about an
    # unseen composition, so a test composition must never have been trained on.
    deployment, _ = group_ids([a, b])
    assert len(deployment) == 1
    # Isomer mode deliberately places them apart. That is the whole point of it.
    isomer, _ = group_ids([a, b], mode=SplitMode.ISOMER_DISCRIMINATION)
    assert len(isomer) == 2


def test_isomer_mode_refuses_isomers_the_features_cannot_tell_apart():
    # Fitting the structure block does not buy isomer discrimination outright: it
    # buys it only for the isomer classes these 24 features can see. Two linear
    # chains differing only in linkage position are invisible to them - same
    # topology, same residue counts, same alpha and beta totals, and no branched
    # core so the arm features are absent. Splitting such a pair would measure
    # memorisation, and the guard refuses rather than letting it look like signal.
    composition = "Hex3HexNAc2"
    a = a_measurement(glycan=a_glycan(composition, iupac_condensed=linear(3, 2, "a1-2"), **RESOLVED))
    b = a_measurement(glycan=a_glycan(composition, iupac_condensed=linear(3, 2, "a1-3"), **RESOLVED), ccs=222.0)
    assert key_of(a.glycan.iupac_condensed) != key_of(b.glycan.iupac_condensed)  # two molecules...
    rows = [extract(a).as_row(), extract(b).as_row()]
    assert rows[0] == rows[1]  # ...that the model cannot tell apart
    grouping = group_records([a, b], mode=SplitMode.ISOMER_DISCRIMINATION)
    assert len(set(grouping.group_of)) == 2  # the split would separate them
    with pytest.raises(LeakageError, match="cannot tell them apart"):
        assert_features_do_not_straddle(grouping, rows)


def test_isomer_mode_groups_on_the_structure_and_still_keeps_duplicates_together():
    a, b = isomer_pairs(compositions=1)
    same_again = a_measurement(glycan=a.glycan, ccs=333.0)  # one more measurement of structure a
    ids, grouping = group_ids([a, b, same_again], mode=SplitMode.ISOMER_DISCRIMINATION)
    assert len(ids) == 2
    assert grouping.group_of[0] == grouping.group_of[2]  # the repeat is not a new group
    assert grouping.mode is SplitMode.ISOMER_DISCRIMINATION


def test_isomer_mode_refuses_a_record_with_no_structure():
    # It groups on the structure, and a composition-only record has none. Making
    # it a singleton would silently place it wherever the packing put it.
    a, _ = isomer_pairs(compositions=1)
    with pytest.raises(SplitRefused, match="no structure"):
        group_records([a, a_measurement()], mode=SplitMode.ISOMER_DISCRIMINATION)


def test_an_isomer_mode_split_carries_its_mode_on_the_report():
    report = grouped_folds(isomer_pairs(), mode=SplitMode.ISOMER_DISCRIMINATION)
    assert report.mode is SplitMode.ISOMER_DISCRIMINATION
    assert "isomer_discrimination" in report.summary()
    assert len(report.group_sizes) >= MIN_ANALYTE_GROUPS


def test_the_default_mode_is_deployment():
    report = grouped_folds(many())
    assert report.mode is SplitMode.DEPLOYMENT
    assert "deployment" in report.summary()


def test_a_doi_written_two_ways_is_one_study_atom():
    a = provenance_atoms(a_measurement(source="https://doi.org/10.1000/Paper-A"))
    b = provenance_atoms(a_measurement(source="doi:10.1000/paper-a"))
    assert a == b


@pytest.mark.parametrize(
    "spelling",
    [
        "10.1000/paper-a",
        "doi:10.1000/paper-a",
        "https://doi.org/10.1000/paper-a",
        "Smith et al. (10.1000/paper-a)",
        "Smith et al., 10.1000/paper-a.",
        "[10.1000/paper-a]",
        "see 10.1000/paper-a;",
        "10.1000/paper-a,",
    ],
)
def test_one_paper_cited_eight_ways_is_one_source(spelling):
    # Trailing punctuation left by a citation is the UNSAFE direction: it fails
    # to merge, so one paper counts as several sources and a corpus can clear the
    # group floor on study diversity that does not exist.
    assert provenance_atoms(a_measurement(source=spelling)) == ("doi:10.1000/paper-a",)


def test_a_bracket_that_belongs_to_the_doi_is_kept():
    # Only an unmatched closing bracket is dropped, because a DOI suffix may
    # legitimately contain one.
    assert provenance_atoms(a_measurement(source="10.1000/paper(a)")) == ("doi:10.1000/paper(a)",)


def test_a_single_source_corpus_is_refused_with_a_count_that_is_true_of_it():
    # The corpus holds 30 analyte keys. Reporting "1 analyte groups" would be a
    # false claim about the data: 1 is the component count after the study axis
    # welded every key together through the shared source.
    records = [
        a_measurement(glycan=a_glycan(f"Hex{3 + i // 6}HexNAc{2 + i % 6}"), source="10.1000/only-paper")
        for i in range(30)
    ]
    with pytest.raises(SplitRefused) as caught:
        grouped_folds(records, by_study=True)
    message = str(caught.value)
    assert "grouping components" in message  # not "analyte groups"
    assert "30 analyte keys" in message  # and it says what the corpus really holds
    assert "cross-study" in message
    # The same corpus splits perfectly well within-study.
    assert grouped_folds(records).refusals == ()
