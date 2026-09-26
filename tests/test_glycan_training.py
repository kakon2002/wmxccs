"""What must be true before a model may be fitted, and what is reported while it is not.

No model is fitted anywhere in this file, deliberately. Exercising the ready
branch on synthetic records would produce exactly the thing this milestone
refuses to produce, so the tests that reach readiness stop there.

Every record here is a synthetic fixture, and says so: each carries the reuse
status synthetic_fixture, which the gate accepts without a licence record and
the file loader refuses. No CCS value is real.
"""

import inspect
import math
import subprocess
import sys
from pathlib import Path

import pytest

from wmxglycan.contracts import DataMaturity
from wmxglycan.evaluation import MIN_GROUP_RECORDS, MIN_GROUP_STRUCTURES
from wmxglycan.licensing import LicenceGateError, TrainingGateError, UnbackedClaimError
from wmxglycan.models import CCSMeasurement, GlycanStructure
from wmxglycan.splits import MIN_ANALYTE_GROUPS, AnalyteLevel, SplitRefused, analyte_key, grouped_folds
from wmxglycan.training import (
    BASELINES,
    CONFORMAL_ALPHA,
    MIN_CALIBRATION_RECORDS,
    MIN_SOURCES_FOR_CROSS_STUDY,
    MIN_TRAINING_RECORDS,
    TARGET_ANALYTE_GROUPS,
    TOO_FEW_GROUPS_TO_FIT,
    TOO_FEW_RECORDS,
    GroupCount,
    InsufficientTrainingDataError,
    Readiness,
    TrainingSet,
    calibration_refusal,
    calibration_warning,
    fit_ccs_baseline,
    smallest_calibration_set,
)

FIXTURE_SOURCE = "synthetic test fixture, not a real record"


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


def linear(hexes, hexnacs, arm="a1-2"):
    """A parseable glycan holding exactly that composition.

    Fixtures have to carry real structures now: a composition-only record is
    refused from a structure-feature fit, because it has nothing to teach about
    isomers. Varying `arm` gives a second structure of the same composition.
    """
    return "Man(%s)" % arm * (hexes - 1) + "Man(b1-4)" + "GlcNAc(b1-4)" * (hexnacs - 1) + "GlcNAc"


def a_corpus(records=MIN_TRAINING_RECORDS, compositions=TARGET_ANALYTE_GROUPS, with_structure=True):
    """Enough cleared records, over enough analytes, to clear every readiness floor."""
    made = []
    for index in range(records):
        which = index % compositions
        hexes, hexnacs = 3 + which // 6, 2 + which % 6
        extra = {"iupac_condensed": linear(hexes, hexnacs)} if with_structure else {}
        made.append(
            a_measurement(glycan=a_glycan(f"Hex{hexes}HexNAc{hexnacs}", **extra), ccs=100.0 + index)
        )
    return made


# --- the conformal floor is derived, not chosen --------------------------------


@pytest.mark.parametrize(
    "alpha, smallest",
    [(0.32, 3), (0.20, 4), (0.10, 9), (0.05, 19), (0.02, 49), (0.01, 99)],
)
def test_the_smallest_calibration_set_is_where_the_quantile_first_exists(alpha, smallest):
    assert smallest_calibration_set(alpha) == smallest
    # Below it, the index the interval needs is larger than the set itself.
    assert math.ceil((smallest - 1 + 1) * (1 - alpha)) > smallest - 1
    # At it, the index fits.
    assert math.ceil((smallest + 1) * (1 - alpha)) <= smallest


def test_the_floor_matches_the_closed_form():
    for alpha in (0.32, 0.20, 0.10, 0.05, 0.02, 0.01):
        assert smallest_calibration_set(alpha) == math.ceil(1 / alpha) - 1


def test_the_default_floor_is_nine_for_ninety_per_cent():
    assert CONFORMAL_ALPHA == 0.10
    assert MIN_CALIBRATION_RECORDS == 9
    assert MIN_CALIBRATION_RECORDS == smallest_calibration_set(CONFORMAL_ALPHA)


def test_the_conformal_floor_is_derived_rather_than_written_down():
    # Comparing values cannot tell a derivation from a literal here, because 9 IS
    # the right answer for alpha=0.10: a hard-coded 9 passes the equality above.
    # The two only diverge if CONFORMAL_ALPHA moves, at which point a literal
    # would silently stop matching the coverage it claims to support. So the
    # binding is what gets pinned.
    import wmxglycan.training as training_module

    source = inspect.getsource(training_module)
    assert "MIN_CALIBRATION_RECORDS = smallest_calibration_set()" in source
    # And the derivation genuinely moves with alpha, so the binding is not idle.
    assert smallest_calibration_set(0.05) != smallest_calibration_set(0.10)


@pytest.mark.parametrize("held, ordinal", [(0, "1st"), (1, "2nd"), (2, "3rd"), (3, "4th")])
def test_the_refusal_counts_in_readable_ordinals(held, ordinal):
    # "the 1th smallest score" reads as a bug and undermines the number beside it.
    reason = calibration_refusal(held)
    assert f"the {ordinal} smallest score" in reason


def test_a_negative_calibration_count_is_refused_rather_than_described():
    # Otherwise it reports "would need the -3th smallest score of -5".
    with pytest.raises(ValueError, match="cannot hold"):
        calibration_refusal(-5)


@pytest.mark.parametrize("held", [0, 1, 8])
def test_a_calibration_set_too_small_for_the_coverage_is_refused(held):
    reason = calibration_refusal(held)
    assert reason is not None
    assert "does not exist" in reason  # no finite quantile, so no honest interval
    assert str(MIN_CALIBRATION_RECORDS) in reason


def test_at_exactly_the_floor_the_interval_is_valid_and_uninformative():
    # Weakness, not impossibility, so it is reported rather than refused.
    assert calibration_refusal(MIN_CALIBRATION_RECORDS) is None
    warning = calibration_warning(MIN_CALIBRATION_RECORDS)
    assert warning is not None and "full observed range" in warning


def test_above_the_floor_there_is_neither_refusal_nor_warning():
    assert calibration_refusal(MIN_CALIBRATION_RECORDS + 1) is None
    assert calibration_warning(MIN_CALIBRATION_RECORDS + 1) is None


@pytest.mark.parametrize("alpha", [0.0, 1.0, -0.1, 1.5])
def test_a_coverage_outside_zero_to_one_is_refused(alpha):
    with pytest.raises(ValueError, match="strictly between"):
        smallest_calibration_set(alpha)


# --- readiness at zero records --------------------------------------------------


def test_readiness_is_producible_at_zero_records_and_prints_an_honest_zero():
    readiness = TrainingSet().readiness
    assert readiness.records_held == 0
    assert readiness.records_cleared == 0
    assert not readiness.ready
    text = readiness.summary()
    assert "records held               0" in text
    assert "ready to fit: no" in text
    assert str(MIN_TRAINING_RECORDS) in text  # says what it would need


def test_readiness_is_always_provisional_and_carries_the_count():
    readiness = TrainingSet().readiness
    assert readiness.maturity.data_maturity is DataMaturity.PROVISIONAL
    assert readiness.maturity.training_record_count == 0


def test_the_summary_names_every_floor_so_the_cost_is_on_the_same_page():
    text = TrainingSet().readiness.summary()
    for floor in (MIN_TRAINING_RECORDS, MIN_ANALYTE_GROUPS, TARGET_ANALYTE_GROUPS, MIN_SOURCES_FOR_CROSS_STUDY):
        assert str(floor) in text
    assert str(MIN_CALIBRATION_RECORDS) in text


def test_a_short_but_evaluable_corpus_goes_ahead_and_says_it_will_be_weak():
    # The rule is inevaluability, not weakness. Thirty records can be split and
    # scored, so a fit proceeds; refusing it would suppress an honest weak result.
    readiness = TrainingSet(records=tuple(a_corpus(records=30, compositions=30))).readiness
    assert readiness.ready
    assert readiness.refusal is None
    assert readiness.blockers == ()
    assert TOO_FEW_RECORDS.format(held=30, needed=MIN_TRAINING_RECORDS) in readiness.warnings
    text = readiness.summary()
    assert "ready to fit: yes" in text
    assert "will be weak" in text and "cannot be evaluated" not in text


def test_the_record_count_travels_with_the_weak_result():
    # The caveat rides on the number rather than sitting in a report elsewhere.
    readiness = TrainingSet(records=tuple(a_corpus(records=30, compositions=30))).readiness
    assert readiness.maturity.training_record_count == 30
    assert readiness.maturity.data_maturity is DataMaturity.PROVISIONAL


def test_a_corpus_of_one_glycan_is_short_of_analyte_diversity():
    readiness = TrainingSet(records=tuple(a_corpus(records=MIN_TRAINING_RECORDS, compositions=1))).readiness
    assert not readiness.ready
    # Asserted as the formatted rule, not as the substring "analyte groups",
    # which the calibration-group shortfall also contains: matching the substring
    # let the whole diversity check be deleted without any test failing.
    assert TOO_FEW_GROUPS_TO_FIT.format(groups=1, needed=MIN_ANALYTE_GROUPS) in readiness.refusal


def test_enough_records_over_too_few_analytes_is_not_ready():
    # 200 records is plenty; twelve glycans is not.
    readiness = TrainingSet(records=tuple(a_corpus(records=MIN_TRAINING_RECORDS, compositions=12))).readiness
    assert readiness.records_cleared == MIN_TRAINING_RECORDS
    assert readiness.analyte_groups == 12
    assert not readiness.ready
    assert TOO_FEW_GROUPS_TO_FIT.format(groups=12, needed=MIN_ANALYTE_GROUPS) in readiness.refusal


def test_readiness_counts_groups_the_way_the_splitter_counts_them():
    # An accession shared by two different compositions merges them, so the
    # corpus holds 25 analyte keys but only 24 grouping components. Counting keys
    # here would certify it ready and let grouped_folds refuse it afterwards.
    records = a_corpus(records=MIN_TRAINING_RECORDS, compositions=MIN_ANALYTE_GROUPS)
    shared = "G00000AA"
    # Both compositions are already in the corpus, and several other records
    # carry each, so the key count is unchanged and only the components merge.
    # Substituting a new composition instead would add a key as it merged one,
    # leaving the totals level and hiding the very divergence under test.
    records[0] = a_measurement(glycan=a_glycan("Hex3HexNAc2", glytoucan_ac=shared))
    records[1] = a_measurement(glycan=a_glycan("Hex3HexNAc3", glytoucan_ac=shared))
    readiness = TrainingSet(records=tuple(records)).readiness
    keys = len({analyte_key(record) for record in records})
    assert keys == readiness.analyte_groups + 1  # the merge is real
    assert readiness.analyte_groups < MIN_ANALYTE_GROUPS
    assert not readiness.ready
    # And the splitter agrees, rather than contradicting the readiness report.
    with pytest.raises(SplitRefused):
        grouped_folds(records)


def test_a_calibration_group_too_small_to_score_is_named():
    readiness = TrainingSet(records=tuple(a_corpus(records=MIN_TRAINING_RECORDS))).readiness
    assert readiness.by_group
    for group, count in readiness.by_group.items():
        assert isinstance(count, GroupCount)
    assert "by calibration group" in readiness.summary()


def test_group_counts_need_both_records_and_analyte_diversity():
    assert GroupCount(records=MIN_GROUP_RECORDS, structures=MIN_GROUP_STRUCTURES).scorable
    assert not GroupCount(records=MIN_GROUP_RECORDS - 1, structures=MIN_GROUP_STRUCTURES).scorable
    # Remeasuring one glycan a hundred times is not diversity.
    assert not GroupCount(records=100, structures=2).scorable


# --- the gate cannot be routed around -------------------------------------------


def test_a_training_set_cannot_hold_an_ungated_record():
    # Constructed directly, so there is no classmethod to skip.
    with pytest.raises(LicenceGateError, match="unverified"):
        TrainingSet(records=(a_measurement(reuse_status="unverified"),))


def test_a_training_set_reports_every_fault_at_once():
    with pytest.raises(LicenceGateError) as caught:
        TrainingSet(records=(a_measurement(reuse_status="unverified"), a_measurement(reuse_status="academic_only")))
    assert "record 0" in str(caught.value) and "record 1" in str(caught.value)


def test_a_cleared_measurement_wrapping_a_restricted_structure_is_refused():
    with pytest.raises(LicenceGateError, match="academic_only"):
        TrainingSet(records=(a_measurement(glycan=a_glycan(reuse_status="academic_only")),))


def test_a_training_set_refuses_a_record_whose_open_claim_nothing_backs():
    # The bypass the registry check in the loader alone left open: a record
    # constructed directly, claiming a licence nobody recorded. The gate stops
    # it wherever the record was built.
    with pytest.raises(LicenceGateError, match="no DOI"):
        TrainingSet(records=(a_measurement(reuse_status="open_attribution"),))
    with pytest.raises(LicenceGateError, match="no licence record"):
        TrainingSet(records=(a_measurement(reuse_status="open_attribution", doi="10.1000/no-record"),))
    with pytest.raises(LicenceGateError, match="neither the row's DOI nor a registered dataset"):
        TrainingSet(records=(a_measurement(glycan=a_glycan(reuse_status="open_attribution")),))


def test_synthetic_fixtures_are_counted_warned_and_printed():
    readiness = TrainingSet(records=tuple(a_corpus(records=30, compositions=30))).readiness
    assert readiness.synthetic_records == 30
    assert readiness.by_reuse_status == {"synthetic_fixture": 30}
    assert readiness.ready  # trainable: a test can build a set and reach every threshold
    assert any("synthetic fixtures" in w and "no number from it describes real data" in w for w in readiness.warnings)
    assert "SYNTHETIC FIXTURES         30" in readiness.summary()


def a_real_measurement(**overrides):
    """Wellmatix's own record: trainable, backed without any registry entry, and not a fixture."""
    fields = {
        "reuse_status": "internal_proprietary",
        "glycan": a_glycan(reuse_status="internal_proprietary"),
    }
    return a_measurement(**(fields | overrides))


def a_real_corpus(records=MIN_TRAINING_RECORDS, compositions=TARGET_ANALYTE_GROUPS):
    """A ready corpus of records that are NOT fixtures: our own data, with structures."""
    made = []
    for index in range(records):
        which = index % compositions
        hexes, hexnacs = 3 + which // 6, 2 + which % 6
        made.append(
            a_measurement(
                glycan=a_glycan(
                    f"Hex{hexes}HexNAc{hexnacs}",
                    iupac_condensed=linear(hexes, hexnacs),
                    reuse_status="internal_proprietary",
                ),
                reuse_status="internal_proprietary",
                ccs=100.0 + index,
            )
        )
    return made


def a_mixed_corpus(records=MIN_TRAINING_RECORDS, compositions=TARGET_ANALYTE_GROUPS):
    """Half the structures invented by a test, half our own. Every record's own status is real."""
    made = []
    for index in range(records):
        which = index % compositions
        hexes, hexnacs = 3 + which // 6, 2 + which % 6
        made.append(
            a_measurement(
                glycan=a_glycan(
                    f"Hex{hexes}HexNAc{hexnacs}",
                    iupac_condensed=linear(hexes, hexnacs),
                    reuse_status="synthetic_fixture" if index % 2 == 0 else "internal_proprietary",
                ),
                reuse_status="internal_proprietary",
                ccs=100.0 + index,
            )
        )
    return made


def test_a_ready_set_of_real_records_reaches_the_fit_without_being_called_a_test():
    # Every corpus in this file was synthetic, so a guard keyed on the wrong
    # count - refusing any non-empty set rather than a set holding fixtures -
    # would have gone unnoticed, and no real corpus could ever have been fitted
    # without the caller declaring it a test. A real one must get through alone.
    ready = TrainingSet(records=tuple(a_real_corpus()))
    assert ready.readiness.ready and ready.readiness.synthetic_records == 0
    with pytest.raises(NotImplementedError, match="no model has ever been trained"):
        fit_ccs_baseline(ready)


def test_the_synthetic_warning_names_the_synthetic_count_not_the_total():
    readiness = TrainingSet(records=tuple(a_mixed_corpus(records=10, compositions=10))).readiness
    assert readiness.synthetic_records == 5
    warning = next(note for note in readiness.warnings if "synthetic fixtures" in note)
    assert "5 of 10" in warning  # not "10 of 10": a report about honesty must count honestly


def test_the_fits_refusal_names_the_synthetic_count_not_the_whole_set():
    ready = TrainingSet(records=tuple(a_mixed_corpus()))
    half = MIN_TRAINING_RECORDS // 2
    assert ready.readiness.synthetic_records == half
    with pytest.raises(TrainingGateError) as caught:
        fit_ccs_baseline(ready)
    assert f"{half} synthetic fixture(s) of {MIN_TRAINING_RECORDS} records" in str(caught.value)


def test_a_set_of_real_records_reports_no_synthetic_records():
    # A real set, not an empty one: emptiness would satisfy these assertions on its own.
    readiness = TrainingSet(records=tuple(a_real_measurement(ccs=100.0 + i) for i in range(30))).readiness
    assert len(readiness.by_reuse_status) and readiness.by_reuse_status == {"internal_proprietary": 30}
    assert readiness.synthetic_records == 0
    assert "SYNTHETIC" not in readiness.summary()
    assert not any("synthetic" in warning for warning in readiness.warnings)


def test_a_fixture_declared_on_the_structure_is_still_a_fixture():
    # The defect the review found: counting only each record's OWN status let a
    # measurement claiming a real status carry an invented structure, and the
    # whole set printed and fitted as data with nothing marking it.
    records = tuple(
        a_real_measurement(
            glycan=a_glycan(f"Hex{3 + i // 6}HexNAc{2 + i % 6}", iupac_condensed=linear(3 + i // 6, 2 + i % 6)),
            ccs=100.0 + i,
        )
        for i in range(30)
    )
    assert {str(r.glycan.reuse_status) for r in records} == {"synthetic_fixture"}  # a_glycan's default
    assert {str(r.reuse_status) for r in records} == {"internal_proprietary"}
    readiness = TrainingSet(records=records).readiness
    assert readiness.synthetic_records == 30
    assert any("structure inside it" in warning for warning in readiness.warnings)
    assert "SYNTHETIC FIXTURES" in readiness.summary()


def test_a_fixture_on_the_structure_cannot_reach_a_fit_unannounced():
    records = tuple(
        a_real_measurement(
            glycan=a_glycan(f"Hex{3 + i // 6}HexNAc{2 + i % 6}", iupac_condensed=linear(3 + i // 6, 2 + i % 6)),
            ccs=100.0 + i,
        )
        for i in range(MIN_TRAINING_RECORDS)
    )
    ready = TrainingSet(records=records)
    assert ready.readiness.ready
    with pytest.raises(TrainingGateError, match="synthetic fixture"):
        fit_ccs_baseline(ready)
    with pytest.raises(NotImplementedError):
        fit_ccs_baseline(ready, synthetic_ok=True)


def test_the_status_tally_counts_each_record_once_and_the_synthetic_count_looks_deeper():
    # The two counts answer different questions, and this is the case that
    # separates them: every record's OWN status here is a real one, so
    # by_reuse_status says "ten real records" and sums to the records held -
    # while half of them are built from a structure invented in a test, which
    # only the synthetic count can see.
    real = tuple(a_real_measurement(ccs=100.0 + i) for i in range(5))
    hiding_a_fixture = tuple(a_real_measurement(glycan=a_glycan(), ccs=200.0 + i) for i in range(5))
    readiness = TrainingSet(records=real + hiding_a_fixture).readiness
    assert sum(readiness.by_reuse_status.values()) == readiness.records_cleared == 10
    assert readiness.by_reuse_status == {"internal_proprietary": 10}
    assert readiness.synthetic_records == 5
    assert "SYNTHETIC FIXTURES         5" in readiness.summary()


def test_an_unbacked_claim_keeps_its_class_through_a_training_set():
    # The loader distinguishes an unbacked claim from a licence refusal; a caller
    # that catches UnbackedClaimError around a TrainingSet must get the same
    # answer, or it is told to relabel a record whose licence is fine.
    unbacked = a_measurement(reuse_status="open_attribution", glycan=a_glycan(reuse_status="internal_proprietary"))
    with pytest.raises(UnbackedClaimError):
        TrainingSet(records=(unbacked,))
    # A real licence refusal in the same set is the more serious fault and wins.
    with pytest.raises(LicenceGateError) as caught:
        TrainingSet(records=(unbacked, a_measurement(reuse_status="academic_only")))
    assert not isinstance(caught.value, UnbackedClaimError)


def test_an_analyte_gap_is_not_reported_as_a_licence_problem():
    with pytest.raises(TrainingGateError) as caught:
        TrainingSet(records=(a_measurement(reducing_end_label="unknown"),))
    assert not isinstance(caught.value, LicenceGateError)


# --- fitting refuses --------------------------------------------------------------


def test_an_empty_training_set_refuses_to_fit():
    with pytest.raises(InsufficientTrainingDataError, match="nothing to fit") as caught:
        fit_ccs_baseline(TrainingSet())
    assert not isinstance(caught.value, LicenceGateError)  # a data gap, not a licence problem
    assert not isinstance(caught.value, (ValueError, TypeError))  # not swallowable by a row-level except


@pytest.mark.parametrize("held", [30, MIN_TRAINING_RECORDS - 1])
def test_a_corpus_below_the_volume_floor_is_warned_about_and_not_refused(held):
    # Both of these used to refuse. The volume floor is a warning now, because a
    # corpus this size is perfectly evaluable and a weak result that reports its
    # own weakness is useful.
    readiness = TrainingSet(records=tuple(a_corpus(records=held, compositions=min(held, 40)))).readiness
    assert readiness.ready
    assert any(str(held) in warning and str(MIN_TRAINING_RECORDS) in warning for warning in readiness.warnings)


@pytest.mark.parametrize(
    "records, compositions, blocked",
    [
        (0, 0, True),  # nothing to fit
        (MIN_TRAINING_RECORDS, 1, True),  # cannot be divided into folds
        (30, 30, False),  # small, but evaluable
        (MIN_TRAINING_RECORDS, 100, False),
    ],
)
def test_only_inevaluability_blocks_a_fit(records, compositions, blocked):
    corpus = a_corpus(records=records, compositions=compositions) if records else []
    readiness = TrainingSet(records=tuple(corpus)).readiness
    assert bool(readiness.blockers) is blocked
    assert readiness.ready is not blocked


def test_a_corpus_that_can_be_split_but_scored_nowhere_is_blocked():
    # Thirty analytes spread across thirty calibration groups: the split is fine,
    # and no group holds enough to be scored within. CCS values are not
    # comparable across calibration groups, so there is nowhere to compute an
    # error at all. That is inevaluability, not weakness, so it blocks.
    records = [
        a_measurement(glycan=a_glycan(f"Hex{3 + i // 6}HexNAc{2 + i % 6}"), calibrant=f"fixture calibrant {i}")
        for i in range(30)
    ]
    readiness = TrainingSet(records=tuple(records)).readiness
    assert readiness.analyte_groups == 30  # the split itself is not the problem
    assert len(readiness.by_group) == 30
    assert not any(count.scorable for count in readiness.by_group.values())
    assert not readiness.ready
    blocker = "; ".join(readiness.blockers)
    assert "nothing could be scored" in blocker
    assert str(MIN_GROUP_RECORDS) in blocker and str(MIN_GROUP_STRUCTURES) in blocker
    # And the same records in ONE calibration group are scorable, so the fixture
    # is discriminating rather than merely small.
    together = TrainingSet(records=tuple(a_corpus(records=30, compositions=30))).readiness
    assert together.ready


def test_a_corpus_with_no_resolved_structures_cannot_train_structure_features():
    # Composition-only records have nothing to teach about isomers, and the
    # structure block is the only signal that separates them.
    readiness = TrainingSet(records=tuple(a_corpus(with_structure=False))).readiness
    assert readiness.records_with_structure == 0
    assert readiness.records_without_structure == MIN_TRAINING_RECORDS
    assert not readiness.ready
    blocker = "; ".join(readiness.blockers)
    assert "no resolved structure" in blocker or "resolved structure" in blocker
    # Named as a curation gap, so it is not mistaken for a licence problem.
    assert "curation" in blocker and "licence" in blocker


def test_records_blocked_on_curation_are_counted_apart_from_the_rest():
    mixed = a_corpus(records=MIN_TRAINING_RECORDS, compositions=50)
    for index in range(20):
        mixed[index] = a_measurement(glycan=a_glycan("Hex3HexNAc2"))  # composition only
    readiness = TrainingSet(records=tuple(mixed)).readiness
    assert readiness.records_without_structure == 20
    assert readiness.records_with_structure == MIN_TRAINING_RECORDS - 20
    # Some structures remain, so this is weakness rather than inevaluability.
    assert readiness.ready
    assert any("blocked on curation, not on licence" in warning for warning in readiness.warnings)
    assert "composition-only" in readiness.summary()


def test_sources_are_counted_as_studies_not_as_source_strings():
    # Twenty-four rows transcribed from one table carry one DOI and, plausibly,
    # twenty-four different source texts. That is one source. Counting strings
    # would silence the warning that no cross-study claim can be made. The DOI
    # is invented with the rest of the record: a synthetic fixture citing a
    # paper whose licence is on record would be refused as that paper's record.
    records = [
        a_measurement(
            glycan=a_glycan("Hex3HexNAc2", iupac_condensed=linear(3, 2)),
            source=f"{FIXTURE_SOURCE}, row {i}",
            doi="10.1000/one-invented-paper",
            ccs=100.0 + i,
        )
        for i in range(24)
    ]
    readiness = TrainingSet(records=tuple(records)).readiness
    assert readiness.distinct_sources == 1
    assert any("no cross-study claim" in warning for warning in readiness.warnings)


def test_without_a_doi_a_source_is_its_text():
    # No DOI to prefer, so the text is the study identity, as the splitter treats it.
    records = [
        a_measurement(glycan=a_glycan("Hex3HexNAc2", iupac_condensed=linear(3, 2)), source=f"lab {i}", ccs=100.0 + i)
        for i in range(3)
    ]
    assert TrainingSet(records=tuple(records)).readiness.distinct_sources == 3


def test_the_volume_floor_is_never_a_blocker():
    # Pinned as a property rather than an example: no corpus is ever blocked for
    # holding too few records alone.
    for held in (1, 30, 199):
        readiness = TrainingSet(records=tuple(a_corpus(records=held, compositions=min(held, 40)))).readiness
        assert not any(str(MIN_TRAINING_RECORDS) in blocker for blocker in readiness.blockers)


def test_handing_the_fit_a_bare_list_is_refused_by_the_gate_not_by_a_type_error():
    with pytest.raises(TrainingGateError, match="not a list") as caught:
        fit_ccs_baseline([a_measurement()])
    assert not isinstance(caught.value, TypeError)


def test_an_unknown_baseline_is_refused_and_names_the_registry():
    with pytest.raises(InsufficientTrainingDataError) as caught:
        fit_ccs_baseline(TrainingSet(), "a_baseline_that_does_not_exist")
    # The emptiness refusal comes first, which is itself the right order.
    assert "nothing to fit" in str(caught.value)


@pytest.mark.parametrize("baseline", BASELINES)
def test_every_baseline_goes_through_the_same_guard(baseline):
    with pytest.raises(InsufficientTrainingDataError):
        fit_ccs_baseline(TrainingSet(), baseline)


def test_no_floor_is_a_parameter():
    # A floor that can be passed in is a floor that can be passed a zero.
    names = set(inspect.signature(fit_ccs_baseline).parameters)
    assert not {n for n in names if any(w in n for w in ("min", "max", "threshold", "floor", "tolerance", "seed"))}


def test_the_refusal_survives_python_optimised_mode():
    # Catches any bare assert standing in for a check.
    code = (
        "from wmxglycan.training import TrainingSet, fit_ccs_baseline, InsufficientTrainingDataError\n"
        "try:\n"
        "    fit_ccs_baseline(TrainingSet())\n"
        "except InsufficientTrainingDataError:\n"
        "    print('refused')\n"
    )
    result = subprocess.run([sys.executable, "-O", "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "refused" in result.stdout


# --- no estimator library is loaded, and nothing fits -----------------------------


# DERIVED FROM THE PACKAGE, not written out here. The original was a hand-written list of
# eleven module names, and the port is what exposed it: the list named `wmxglycan.api`,
# which this repository deliberately does not carry, so the guard failed on a straight copy
# for a reason that had nothing to do with estimators. It would have failed the same way the
# day anyone added a module and forgot the list - which is the failure class this project
# keeps finding, a guard that looks tested and is not. Derived, it covers whatever modules
# exist, including ones added after this test was written. The floor below stops an empty
# glob from making it vacuous, which would be the same defect moved one level up.
GLYCAN_PACKAGE = Path(__file__).resolve().parent.parent / "src" / "wmxglycan"
PACKAGE_MODULES = [
    f"wmxglycan.{path.stem}" for path in sorted(GLYCAN_PACKAGE.glob("*.py")) if path.stem != "__init__"
]


def test_the_derived_module_list_actually_covers_the_package():
    # Without this, a glob returning nothing would leave the estimator guard green while
    # importing nothing at all. The named six are the modules the port was asked for.
    assert len(PACKAGE_MODULES) >= 15, PACKAGE_MODULES
    for name in ("composition", "enumeration", "features", "splits", "evaluation", "training"):
        assert f"wmxglycan.{name}" in PACKAGE_MODULES


def test_no_estimator_library_is_imported_by_the_package():
    # Importing the package must not drag in an estimator. Any accidental fit
    # anywhere in the suite would make this fail.
    code = (
        "import sys\n"
        + "".join(f"import {name}\n" for name in PACKAGE_MODULES)
        + "print(sorted(m for m in ('sklearn', 'torch', 'xgboost', 'lightgbm') if m in sys.modules))\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "[]"


def test_nothing_in_the_package_fits_anything_today():
    # An AST walk, not a substring scan: getattr(estimator, "fit")(...) contains
    # no ".fit(" and would sail past a text search.
    import ast

    source_dir = Path(__file__).resolve().parent.parent / "src" / "wmxglycan"
    fitting = set()
    for path in sorted(source_dir.glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            called = node.func
            if isinstance(called, ast.Attribute) and called.attr in {"fit", "partial_fit", "fit_transform"}:
                fitting.add(f"{path.name}:{node.lineno}")
            if (
                isinstance(called, ast.Call)
                and isinstance(called.func, ast.Name)
                and called.func.id == "getattr"
                and len(called.args) >= 2
                and isinstance(called.args[1], ast.Constant)
                and called.args[1].value in {"fit", "partial_fit"}
            ):
                fitting.add(f"{path.name}:{node.lineno}")
    assert not fitting, f"no model is fitted in this milestone: {sorted(fitting)}"


# --- the entry points past the readiness gate, which a green suite otherwise never reaches ---


def test_a_ready_training_set_still_refuses_because_no_model_exists():
    # Without this, replacing the terminal raise with a real fit, or with a
    # fabricated accuracy figure, leaves the whole suite green. The corpus is
    # synthetic and the call says so; that is the only way past the fixture guard.
    ready = TrainingSet(records=tuple(a_corpus()))
    assert ready.readiness.ready
    with pytest.raises(NotImplementedError, match="no model has ever been trained"):
        fit_ccs_baseline(ready, synthetic_ok=True)


def test_a_fit_refuses_synthetic_fixtures_unless_the_call_says_it_is_a_test():
    # The gate lets a fixture through so a test can build a set; a fit on
    # fixtures is a test and not a model, and must be asked for as one.
    ready = TrainingSet(records=tuple(a_corpus()))
    assert ready.readiness.ready  # not a readiness refusal: the set is evaluable, it is just not data
    with pytest.raises(TrainingGateError, match="synthetic fixture") as caught:
        fit_ccs_baseline(ready)
    assert not isinstance(caught.value, (LicenceGateError, InsufficientTrainingDataError))
    assert "synthetic_ok=True" in str(caught.value)


def test_a_ready_set_with_an_unknown_baseline_is_refused_by_name():
    ready = TrainingSet(records=tuple(a_corpus()))
    with pytest.raises(InsufficientTrainingDataError, match="not a known baseline") as caught:
        fit_ccs_baseline(ready, "a_baseline_that_does_not_exist")
    for known in BASELINES:
        assert known in str(caught.value)


def test_evaluating_a_baseline_never_returns_a_number_in_this_milestone():
    from wmxglycan.evaluation import evaluate_baseline

    with pytest.raises(InsufficientTrainingDataError, match="nothing to fit"):
        evaluate_baseline(TrainingSet(), "group_median")
    ready = TrainingSet(records=tuple(a_corpus()))
    with pytest.raises(TrainingGateError, match="synthetic fixture"):
        evaluate_baseline(ready, "group_median")  # the fixture guard is the fit's, reached through scoring too
    with pytest.raises(NotImplementedError):
        evaluate_baseline(ready, "group_median", synthetic_ok=True)


def test_the_m3_report_is_exactly_the_readiness_summary_and_adds_no_number():
    # Equality, not containment: an appended "median relative error 1.8%" would
    # pass a containment check while inventing a benchmark.
    from wmxglycan.evaluation import m3_report

    for training_set in (TrainingSet(), TrainingSet(records=tuple(a_corpus(records=30, compositions=30)))):
        assert m3_report(training_set) == training_set.readiness.summary()


def test_the_estimator_import_sits_below_every_guard():
    source = inspect.getsource(fit_ccs_baseline)
    import_line = source.index("from sklearn")
    # Every guard, by name. The list is the weak point: the synthetic-fixture
    # guard was added and not listed here, so the invariant had a hole exactly
    # its shape and a mutation moving it below the import survived. Add new
    # guards here when you add them there.
    for guard in ("NOT_A_TRAINING_SET", "readiness.refusal", "UNKNOWN_BASELINE", "synthetic_records"):
        assert source.index(guard) < import_line, guard


def test_the_validated_flag_belongs_to_a_later_milestone():
    source_dir = Path(__file__).resolve().parent.parent / "src" / "wmxglycan"
    for name in ("training.py", "evaluation.py", "splits.py", "features.py"):
        assert "VALIDATED" not in (source_dir / name).read_text(encoding="utf-8"), name
