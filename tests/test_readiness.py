"""What readiness refuses, what it merely warns about, and what it prints when it holds nothing.

M0's deliverable IS the refusal. The expected result for the corpus actually
loaded is that no ion is held on two platforms, so the report that says so has to
be exactly right: right about the count, right about WHICH rung of the ladder is
blocking, and right that more data on the one platform already held would not
help. A report that quietly said "ready" over an empty table, or that blamed the
wrong rung, would be worse than no report.

Every test here builds its corpus through the conftest builders, so every record
carries the synthetic-fixture status and every set is visibly a test.
"""

from __future__ import annotations

import pytest

from conftest import measurement, peptide, primary, protein
from wmxccs.identity import FoldingState
from wmxccs.licensing import LicenceGateError, TrainingGateError, UnbackedClaimError
from wmxccs.models import UncertaintyType
from wmxccs.readiness import (
    CONFORMAL_ALPHA,
    MIN_CALIBRATION_RECORDS,
    MIN_DEMING_POINTS,
    MIN_MATCHED_IONS_PER_PAIR,
    MIN_PLATFORMS,
    MIN_RECORDS_HELD,
    MIN_SOURCES_FOR_CROSS_STUDY,
    N_FOLDS,
    TARGET_MATCHED_IONS,
    DataMaturity,
    InsufficientDataError,
    MatchedIonSet,
    PairCount,
    Readiness,
    calibration_refusal,
    calibration_warning,
    conformal_quantile_index,
    smallest_informative_calibration_set,
    interval_readiness,
    require_ready,
    smallest_calibration_set,
)
from wmxccs.reuse import ReuseStatus

# Distinct peptide sequences, one per matched ion a corpus needs. Peptides rather
# than glycans because a glycan needs its reducing-end label and its derivatisation
# stated before it clears the gate, and none of that is what these tests are about.
_SEQUENCES = ("PEPTIDEA", "PEPTIDEC", "PEPTIDED", "PEPTIDEE", "PEPTIDEF", "PEPTIDEG")


def _twims_ion(sequence: str):
    """One calibrated TWIMS measurement of the named peptide."""
    return measurement(analyte=peptide(sequence=sequence))


def _dtims_ion(sequence: str):
    """The same ion on a second platform: stepped-field DTIMS, therefore primary."""
    return primary(analyte=peptide(sequence=sequence))


def _cross_platform_corpus(ion_count: int) -> MatchedIonSet:
    """`ion_count` distinct ions, each measured on BOTH TWIMS and stepped-field DTIMS."""
    records = []
    for sequence in _SEQUENCES[:ion_count]:
        records.append(_twims_ion(sequence))
        records.append(_dtims_ion(sequence))
    return MatchedIonSet(tuple(records))


def _tidy(report: str) -> list[str]:
    """The summary's lines with their column padding collapsed, so a test reads a value not a layout."""
    return [" ".join(line.split()) for line in report.splitlines()]


# --- the thresholds, pinned ---------------------------------------------------------------
#
# Asserted against literals and never against the constants themselves, because
# `assert MIN_PLATFORMS == MIN_PLATFORMS` passes whatever the floor has quietly
# been moved to. These tests exist to notice the move.


def test_the_platform_floor_is_two_because_a_comparison_needs_two_things_to_compare():
    assert MIN_PLATFORMS == 2


def test_the_fold_count_and_the_two_counts_derived_from_it_are_at_their_stated_values():
    assert N_FOLDS == 5
    # One group per fold, so a pair under this cannot be cross-validated at all.
    assert MIN_MATCHED_IONS_PER_PAIR == 5
    # A slope and an intercept leave a residual only from the third point on.
    assert MIN_DEMING_POINTS == 3


def test_the_warning_thresholds_are_at_their_stated_values():
    assert TARGET_MATCHED_IONS == 100
    assert MIN_RECORDS_HELD == 200
    assert MIN_SOURCES_FOR_CROSS_STUDY == 3


def test_the_conformal_calibration_floor_is_nine_held_out_ions():
    assert MIN_CALIBRATION_RECORDS == 9
    assert CONFORMAL_ALPHA == 0.10


# --- zero records -------------------------------------------------------------------------


def test_a_readiness_report_is_producible_when_nothing_at_all_is_held():
    held = MatchedIonSet(()).readiness

    assert isinstance(held, Readiness)
    assert held.records_held == 0
    assert held.records_cleared == 0
    assert held.matched_ion_keys == 0
    assert held.matched_ions_multi_platform == 0
    assert held.platforms == ()
    assert held.by_pair == {}
    assert held.ready is False


def test_the_blocker_for_an_empty_corpus_is_that_nothing_is_held_and_nothing_else():
    held = MatchedIonSet(()).readiness

    assert len(held.blockers) == 1
    assert "no record cleared the licence gate" in held.blockers[0]
    assert held.refusal is not None
    assert "nothing to work with" in held.refusal


def test_an_empty_corpus_prints_an_honest_zero_rather_than_an_empty_table():
    # The failure this catches is a report that renders a blank where a count
    # belongs, which a reader takes for "not applicable" rather than "none".
    lines = _tidy(MatchedIonSet(()).readiness.summary())

    assert "records held 0" in lines
    assert "cleared the licence gate 0" in lines
    assert "distinct matched ions 0" in lines
    assert "by platform pair: none - no ion is held on two platforms" in lines
    assert "ready to compare: no" in lines
    # No per-platform table at all, rather than a heading over nothing.
    assert "by platform:" not in lines
    assert any("cannot be evaluated" in line and "nothing to work with" in line for line in lines)


# --- the blocker ladder, one rung at a time -----------------------------------------------


def test_a_corpus_from_one_platform_is_blocked_because_a_second_platform_is_required():
    # The single most important distinction this report draws: it is not a
    # shortage of data, and loading more TWIMS will not move it.
    one_platform = MatchedIonSet((_twims_ion("PEPTIDEA"), _twims_ion("PEPTIDEA")))
    held = one_platform.readiness

    assert held.platforms == ("TWIMS",)
    assert held.ready is False
    blocker = held.blockers[0]
    assert "a second platform is required" in blocker
    assert "More values on the same platform will not produce one pair" in blocker
    assert "not a shortage of data" in blocker
    assert "TWIMS" in blocker


def test_the_blocker_ladder_short_circuits_at_the_one_platform_rung():
    # Every rung below this one is downstream of it and would only restate it.
    # A reader handed three blockers cannot tell which one to act on.
    held = MatchedIonSet((_twims_ion("PEPTIDEA"), _twims_ion("PEPTIDEC"))).readiness

    assert len(held.blockers) == 1
    only = held.blockers[0]
    assert "a second platform is required" in only
    assert not any("no pair to compare" in blocker for blocker in held.blockers)
    assert not any("could be scored even if one were fitted" in blocker for blocker in held.blockers)


def test_two_platforms_holding_no_common_ion_are_blocked_for_having_no_pair():
    # A native protein ion on one platform and a denatured one on the other are
    # two ions, not a pair, so a second platform on its own buys nothing.
    native = measurement(
        analyte=protein(folding_state=FoldingState.NATIVE), adduct="[M+11H]11+", charge=11, ccs=3000.0
    )
    denatured = primary(
        analyte=protein(folding_state=FoldingState.DENATURED), adduct="[M+11H]11+", charge=11, ccs=4000.0
    )
    held = MatchedIonSet((native, denatured)).readiness

    assert len(held.platforms) == 2
    assert held.matched_ion_keys == 2
    assert held.matched_ions_multi_platform == 0
    assert held.ready is False
    assert len(held.blockers) == 1
    assert "none of them is measured on more than one" in held.blockers[0]
    assert "no pair to compare" in held.blockers[0]


def test_a_pair_that_exists_but_is_too_small_to_cross_validate_is_blocked_not_scored():
    held = _cross_platform_corpus(MIN_MATCHED_IONS_PER_PAIR - 1).readiness

    assert held.matched_ions_multi_platform == 4
    assert held.ready is False
    assert len(held.blockers) == 1
    assert "no platform pair holds 5 matched ions" in held.blockers[0]
    assert "5-fold grouped cross-validation" in held.blockers[0]


# --- the happy direction, so that a Readiness which always refused would fail --------------


def test_one_ion_measured_on_two_platforms_is_counted_and_named_as_a_pair():
    held = _cross_platform_corpus(1).readiness

    assert held.matched_ion_keys == 1
    assert held.matched_ions_multi_platform == 1
    assert len(held.by_pair) == 1
    name, count = next(iter(held.by_pair.items()))
    assert "TWIMS" in name and "DTIMS" in name
    assert count.matched_ions == 1


def test_five_matched_ions_on_two_platforms_clear_every_blocker_but_still_warn_about_size():
    held = _cross_platform_corpus(MIN_MATCHED_IONS_PER_PAIR).readiness

    assert held.matched_ions_multi_platform == 5
    assert held.blockers == ()
    assert held.ready is True
    assert held.refusal is None
    # Ready is not the same as strong, and the report has to keep saying so.
    assert any("200" in warning for warning in held.warnings)
    assert any("100" in warning for warning in held.warnings)
    assert any("cross-study" in warning for warning in held.warnings)


def test_require_ready_returns_the_readiness_of_a_corpus_that_can_be_evaluated():
    held = require_ready(_cross_platform_corpus(MIN_MATCHED_IONS_PER_PAIR))

    assert held.ready is True
    assert held.matched_ions_multi_platform == 5


def test_require_ready_refuses_an_unevaluable_corpus_as_insufficient_data():
    with pytest.raises(InsufficientDataError, match="nothing to work with"):
        require_ready(MatchedIonSet(()))


# --- counting ------------------------------------------------------------------------------


def test_the_multi_platform_count_counts_ions_on_two_platforms_and_not_keys_held():
    # Three keys, only one of which is on both platforms. A count taken from the
    # key total would report three cross-platform ions where one exists.
    records = (
        _twims_ion("PEPTIDEA"),
        _dtims_ion("PEPTIDEA"),
        _twims_ion("PEPTIDEC"),
        _dtims_ion("PEPTIDED"),
    )
    held = MatchedIonSet(records).readiness

    assert held.matched_ion_keys == 3
    assert held.matched_ions_multi_platform == 1


def test_the_maturity_stamp_is_provisional_and_carries_the_multi_platform_count():
    records = (
        _twims_ion("PEPTIDEA"),
        _dtims_ion("PEPTIDEA"),
        _twims_ion("PEPTIDEC"),
        _dtims_ion("PEPTIDED"),
    )
    stamp = MatchedIonSet(records).readiness.maturity

    assert stamp.data_maturity is DataMaturity.PROVISIONAL
    assert stamp.matched_ion_count == 1


def test_the_maturity_stamp_on_an_empty_corpus_is_provisional_at_zero_and_says_so_in_the_report():
    held = MatchedIonSet(()).readiness

    assert held.maturity.data_maturity is DataMaturity.PROVISIONAL
    assert held.maturity.matched_ion_count == 0
    assert "provisional" in held.summary()
    assert "validated" not in held.summary()


def test_a_maturity_stamp_cannot_be_edited_after_the_check_that_accepted_it():
    stamp = MatchedIonSet(()).readiness.maturity

    with pytest.raises(Exception):
        stamp.data_maturity = DataMaturity.VALIDATED


# --- the gate runs on construction ----------------------------------------------------------


def test_a_matched_ion_set_cannot_be_constructed_around_an_unverified_record():
    unverified = measurement(analyte=peptide(), reuse_status=ReuseStatus.UNVERIFIED)

    with pytest.raises(LicenceGateError, match="may not train"):
        MatchedIonSet((unverified,))


def test_an_unverified_record_is_a_licence_fault_and_not_reported_as_an_unbacked_claim():
    # The two point at different remedies: read the terms, versus record who read
    # them. Naming the wrong one sends somebody to the wrong place.
    unverified = measurement(analyte=peptide(), reuse_status=ReuseStatus.UNVERIFIED)

    with pytest.raises(LicenceGateError) as raised:
        MatchedIonSet((unverified,))
    assert not isinstance(raised.value, UnbackedClaimError)
    assert "nobody has checked its terms" in str(raised.value)


def test_a_trainable_claim_that_no_registry_entry_backs_raises_the_unbacked_claim_error():
    # The status is fine and the licence may well be fine; what is missing is the
    # registry entry naming who read the terms. Relabelling the record is not the fix.
    claimed_open = measurement(analyte=peptide(), reuse_status=ReuseStatus.OPEN_ATTRIBUTION)

    with pytest.raises(UnbackedClaimError):
        MatchedIonSet((claimed_open,))


def test_a_non_licence_training_blocker_raises_the_plain_gate_error_and_not_a_licence_error():
    # An uninterpretable uncertainty type is a data fault, not a licence fault.
    # Reporting it as one would send somebody to read a licence that is not the problem.
    unreadable_spread = measurement(analyte=peptide(), uncertainty_type=UncertaintyType.UNKNOWN)

    with pytest.raises(TrainingGateError) as raised:
        MatchedIonSet((unreadable_spread,))
    assert not isinstance(raised.value, LicenceGateError)
    assert "uncertainty_type" in str(raised.value)


def test_a_matched_ion_set_names_how_many_records_may_not_train_and_which_ones():
    good = _twims_ion("PEPTIDEA")
    bad = measurement(analyte=peptide(), reuse_status=ReuseStatus.UNVERIFIED)

    with pytest.raises(LicenceGateError) as raised:
        MatchedIonSet((good, bad))
    assert "1 record(s)" in str(raised.value)
    assert "record 1" in str(raised.value)


def test_require_ready_refuses_a_bare_list_with_a_gate_error_rather_than_a_type_error():
    # A TypeError would be swallowed by any caller that skips malformed rows,
    # and a bare list is the natural way to route around the gate by accident.
    with pytest.raises(TrainingGateError, match="MatchedIonSet") as raised:
        require_ready([_twims_ion("PEPTIDEA")])
    assert not isinstance(raised.value, TypeError)
    assert "has not been through the licence gate" in str(raised.value)


@pytest.mark.parametrize("not_a_set", [[], (), None, 0, "records"])
def test_require_ready_refuses_anything_that_is_not_a_gated_set(not_a_set):
    with pytest.raises(TrainingGateError, match="MatchedIonSet"):
        require_ready(not_a_set)


# --- synthetic fixtures ---------------------------------------------------------------------


def test_a_corpus_of_fixtures_is_counted_and_warned_about_rather_than_presented_as_data():
    held = _cross_platform_corpus(1).readiness

    assert held.synthetic_records == 2
    assert any("synthetic fixtures" in warning for warning in held.warnings)
    assert any("no number from it describes real data" in warning for warning in held.warnings)
    assert any("SYNTHETIC FIXTURES" in line for line in _tidy(held.summary()))


def test_a_fixture_declared_only_on_the_analyte_inside_a_record_is_still_counted():
    # The record's own status says internal_proprietary, which is a real status.
    # Counting the outermost layer alone would print a test as in-house data.
    disguised = measurement(analyte=peptide(), reuse_status=ReuseStatus.INTERNAL_PROPRIETARY)
    held = MatchedIonSet((disguised,)).readiness

    assert held.by_reuse_status == {"internal_proprietary": 1}
    assert held.synthetic_records == 1
    assert any("synthetic fixtures" in warning for warning in held.warnings)


# --- the conformal calibration floor -----------------------------------------------------------


def test_the_smallest_calibration_set_is_nine_at_ninety_per_cent_and_nineteen_at_ninety_five():
    assert smallest_calibration_set(0.10) == 9
    assert smallest_calibration_set(0.05) == 19
    assert smallest_calibration_set() == 9


@pytest.mark.parametrize("alpha", [0.0, 1.0, -0.1, 1.5, 2, -1])
def test_an_alpha_outside_the_open_unit_interval_is_refused(alpha):
    with pytest.raises(ValueError, match="strictly between 0 and 1"):
        smallest_calibration_set(alpha)


def test_a_calibration_set_below_the_floor_is_refused_and_names_the_quantile_that_is_missing():
    refusal = calibration_refusal(8)

    assert refusal is not None
    assert "cannot support 90% coverage" in refusal
    assert "9th smallest score of 8" in refusal
    assert "At least 9 are required" in refusal


@pytest.mark.parametrize("held", [9, 10, 40])
def test_a_calibration_set_at_or_above_the_floor_is_not_refused(held):
    assert calibration_refusal(held) is None


@pytest.mark.parametrize("held", [0, 1, 8])
def test_every_calibration_set_under_the_floor_is_refused(held):
    assert calibration_refusal(held) is not None


def test_a_calibration_set_exactly_at_the_floor_warns_that_the_interval_is_uninformative():
    warning = calibration_warning(9)

    assert warning is not None
    assert "uninformative" in warning
    assert "full observed range" in warning


# THE FLOOR AND THE INFORMATIVE SIZE ARE DIFFERENT NUMBERS, and this test file used to
# assume they were the same. An interval EXISTS from nine calibration points; it is
# narrower than the observed data only from nineteen. Between those the quantile is the
# largest observed score, so the interval spans every residual seen and excludes nothing.
#
# This was wrong until 19 September 2026 and a test pinned it: `calibration_warning(10)`
# was asserted to be None. Found while sizing M4's strata, whose grouped splits land
# calibration sets at 11, 14 and 20 - two of the three inside the silent range.


@pytest.mark.parametrize("held", list(range(9, 19)), ids=lambda n: f"n{n}")
def test_every_calibration_size_whose_quantile_is_the_largest_score_warns(held):
    """Nine through eighteen, not nine alone."""
    warning = calibration_warning(held)

    assert warning is not None, f"n={held} gives a full-range interval and must say so"
    assert "full observed range" in warning
    assert "LARGEST" in warning
    assert conformal_quantile_index(held) == held


@pytest.mark.parametrize("held", [19, 20, 41, 200], ids=lambda n: f"n{n}")
def test_a_calibration_set_large_enough_to_be_informative_does_not_warn(held):
    """From nineteen the quantile is no longer the maximum, so the interval says something."""
    assert calibration_warning(held) is None
    assert conformal_quantile_index(held) < held


def test_the_two_calibration_sizes_are_derived_from_alpha_and_not_tabulated():
    """Both move with the coverage asked for, and the derivation is the code.

    90 per cent: an interval exists from 9 and is informative from 19.
    95 per cent: an interval exists from 19 and is informative from 39.
    """
    assert smallest_calibration_set(0.10) == 9
    assert smallest_informative_calibration_set(0.10) == 19
    assert smallest_calibration_set(0.05) == 19
    assert smallest_informative_calibration_set(0.05) == 39
    # and the informative size is always the larger of the two
    for alpha in (0.01, 0.05, 0.10, 0.20, 0.5):
        assert smallest_informative_calibration_set(alpha) > smallest_calibration_set(alpha)


def test_the_boundary_between_refused_and_degenerate_is_the_quantile_index():
    """Refused when the score does not exist, degenerate when it is the last one."""
    assert conformal_quantile_index(8) == 9 > 8  # no 9th score of 8: refused
    assert calibration_refusal(8) is not None
    assert calibration_warning(8) is None, "a refused set is not a degenerate one"
    assert conformal_quantile_index(9) == 9  # the 9th of 9: the maximum
    assert calibration_refusal(9) is None
    assert calibration_warning(9) is not None


def test_a_negative_calibration_set_is_refused_as_impossible():
    with pytest.raises(ValueError, match="cannot hold"):
        calibration_refusal(-1)


@pytest.mark.parametrize(
    "held,expect_refusal,expect_warning",
    [
        (8, True, False),  # no such score: refused
        (9, False, True),  # the score exists and is the largest: degenerate
        (10, False, True),  # STILL the largest. This case read (False, False) until 19 Sep 2026
        (18, False, True),  # the last degenerate size at 90 per cent
        (19, False, False),  # the first size that says something
    ],
)
def test_interval_readiness_returns_the_refusal_and_the_warning_together(held, expect_refusal, expect_warning):
    refusal, warning = interval_readiness(held)

    assert (refusal is not None) is expect_refusal
    assert (warning is not None) is expect_warning


def test_interval_readiness_at_ninety_five_per_cent_moves_both_sizes_not_only_the_floor():
    """Tighter coverage costs calibration points twice over, and both costs are real.

    At 95 per cent an interval EXISTS from 19 - that is the floor moving - and it is
    narrower than the observed data only from 39. Asking for tighter coverage on a
    calibration set that cannot support it does not fail loudly; it returns the full
    observed range, which is why the warning has to reach that far.
    """
    assert interval_readiness(18, alpha=0.05)[0] is not None, "below the floor: refused"
    assert interval_readiness(19, alpha=0.05)[0] is None, "at the floor: an interval exists"
    assert interval_readiness(19, alpha=0.05)[1] is not None, "and it is the full observed range"
    # 20 is above the floor and still degenerate: this asserted (None, None) until the
    # derivation was corrected.
    assert interval_readiness(20, alpha=0.05)[0] is None
    assert interval_readiness(20, alpha=0.05)[1] is not None
    assert interval_readiness(38, alpha=0.05)[1] is not None, "the last degenerate size at 95 per cent"
    assert interval_readiness(39, alpha=0.05) == (None, None), "the first informative size"


# --- pair counts ------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "matched_ions,fittable,scorable",
    [(0, False, False), (2, False, False), (3, True, False), (4, True, False), (5, True, True), (9, True, True)],
)
def test_a_pair_is_fittable_from_three_points_and_scorable_only_from_five(matched_ions, fittable, scorable):
    count = PairCount(matched_ions=matched_ions)

    assert count.fittable is fittable
    assert count.scorable is scorable


def test_a_pair_under_the_fold_floor_is_warned_about_by_name_in_the_report():
    # The warning is per pair, so a reader can see WHICH pair would give an
    # in-sample figure wearing a cross-validated label.
    held = _cross_platform_corpus(1).readiness

    assert any("in-sample figure" in warning for warning in held.warnings)
    assert any("TWIMS" in warning and "DTIMS" in warning for warning in held.warnings)


def test_a_pair_at_the_fold_floor_draws_no_in_sample_warning():
    held = _cross_platform_corpus(MIN_MATCHED_IONS_PER_PAIR).readiness

    assert held.ready is True
    assert not any("in-sample figure" in warning for warning in held.warnings)


# --- what does not fit yet ----------------------------------------------------------------------


# Was an xfail against a real defect: no antibody or ADC measurement could enter
# a MatchedIonSet, because the nested AntibodyIdentity was offered to the gate as
# a component record and carries no reuse status.
def test_an_antibody_ion_can_be_held_in_a_matched_ion_set_and_reported_on():
    from conftest import antibody

    native = measurement(
        analyte=antibody(folding_state=FoldingState.NATIVE), adduct="[M+24H]24+", charge=24, ccs=7000.0
    )
    held = MatchedIonSet((native,)).readiness

    assert held.records_held == 1
