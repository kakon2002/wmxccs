"""Cyclic ion mobility: what a cyclic record must state, and why it never pools with TWIMS.

THE ASSERTION THIS FILE EXISTS FOR, and it has two halves that must both hold:

  a six-pass value and a single-pass value of ONE ION share a MATCHED-ION KEY,
  because they are the same ion and comparing them is the whole point;

  and they have DIFFERENT CALIBRATION GROUPS, because the calibration that turns
  an arrival time into a cross section is specific to the path the ion took, so
  averaging six passes with one would average two different quantities.

Asserting only the second half passes even when the pass count has wrongly
entered the matched-ion key, and a cyclic value whose key carries its pass count
can never pair with a DTIMS or a TWIMS value of the same ion. The platform would
then report a clean run of zero matches over a corpus full of them. So every
pooling test below asserts both halves, and is named so that nobody can later
delete one half without the name becoming false.

The second thing this file holds is the shape of the record itself: a cyclic
record must state its cyclic settings, no other platform may carry them, and the
settings' training blockers have to reach the measurement that owns them.
"""

from __future__ import annotations

import re

import pytest
from pydantic import ValidationError

from wmxccs.licensing import LicenceGateError, TrainingGateError, assert_trainable

from conftest import cyclic_measurement, cyclic_settings, measurement, small_molecule
from wmxccs.models import CyclicSettings, DTIMSMethod, IMSType, PassMode
from wmxccs.reuse import ReuseStatus

# Every platform that is NOT cyclic, each with the calibrant its method requires.
# None of them has a pass count, so none of them may carry cyclic settings.
NON_CYCLIC_PLATFORMS = {
    "twims": dict(ims_type=IMSType.TWIMS, calibrant="dextran"),
    "tims": dict(ims_type=IMSType.TIMS, calibrant="dextran"),
    "dtims_single_field": dict(
        ims_type=IMSType.DTIMS, dtims_method=DTIMSMethod.SINGLE_FIELD, calibrant="dextran"
    ),
    "dtims_stepped_field": dict(
        ims_type=IMSType.DTIMS, dtims_method=DTIMSMethod.STEPPED_FIELD, calibrant=None
    ),
}

# A cyclic run with every field a source could state actually stated. Used for
# the round trip, and deliberately not used anywhere a default matters.
FULLY_STATED = dict(
    passes=6,
    pass_mode=PassMode.MULTIPASS,
    effective_path_length_m=6.5,
    tw_velocity_m_per_s=375.0,
    tw_height_v=22.0,
    arrival_time_correction="dead time measured on the calibrant and subtracted",
    wrap_around=False,
)


def multipass(**overrides) -> CyclicSettings:
    """Six passes with wrap-around ruled out: the multipass twin of conftest's cyclic_settings().

    Six rather than two so that a count read as a flag, or a flag read as a
    count, shows up in the printed group as well as in the comparison.
    """
    fields = dict(passes=6, pass_mode=PassMode.MULTIPASS, wrap_around=False)
    fields.update(overrides)
    return cyclic_settings(**fields)


def licensed_cyclic(**overrides):
    """A cyclic record whose ANALYTE has no blockers of its own, for gate tests.

    The default glycan fixture leaves its reducing-end label and derivatisation
    unknown, and those block training for reasons that have nothing to do with
    cyclic settings. A gate test that used it would pass whether the cyclic
    layer worked or not.
    """
    fields = dict(analyte=small_molecule())
    fields.update(overrides)
    return cyclic_measurement(**fields)


# --- the assertion this file exists for -------------------------------------------------


def test_a_six_pass_and_a_single_pass_value_of_one_ion_share_a_key_and_do_not_share_a_group():
    # BOTH halves, in one test, because each without the other is misleading.
    # Same key: they are the same ion, and a cyclic value that could not pair
    # with anything would be dead weight in the corpus.
    # Different group: they were calibrated against different effective path
    # lengths, so nothing may average them.
    one_pass = cyclic_measurement(cyclic=cyclic_settings())
    six_passes = cyclic_measurement(cyclic=multipass())

    assert one_pass.matched_ion_key == six_passes.matched_ion_key
    assert one_pass.calibration_group != six_passes.calibration_group


def test_the_pass_count_is_the_only_thing_that_differs_between_those_two_calibration_groups():
    # Which pins WHERE the difference comes from. If the groups differed for some
    # other reason the test above would pass while the pass count sat outside the
    # group entirely, and six passes would still pool with one.
    one_pass = cyclic_measurement(cyclic=cyclic_settings()).calibration_group
    six_passes = cyclic_measurement(cyclic=multipass()).calibration_group

    assert one_pass.cyclic_passes != six_passes.cyclic_passes
    assert one_pass._replace(cyclic_passes=None) == six_passes._replace(cyclic_passes=None)


def test_two_cyclic_values_with_no_count_but_different_modes_share_a_key_and_do_not_pool():
    # A source may report "multipass" without saying how many. The mode alone
    # still decides pooling, so the group carries the mode as well as the count.
    stated_multipass = cyclic_measurement(
        cyclic=CyclicSettings(passes=None, pass_mode=PassMode.MULTIPASS, wrap_around=False)
    )
    stated_single = cyclic_measurement(
        cyclic=CyclicSettings(passes=None, pass_mode=PassMode.SINGLE_PASS, wrap_around=False)
    )

    assert stated_multipass.matched_ion_key == stated_single.matched_ion_key
    assert stated_multipass.calibration_group != stated_single.calibration_group


def test_a_cyclic_value_and_a_single_pass_twims_value_of_one_ion_share_a_key_and_do_not_pool():
    # The M1 headline: cyclic is held apart from single-pass TWIMS. Apart in the
    # calibration group, together in the matched-ion key, which is what "held
    # apart without being made incomparable" means.
    travelling_wave = measurement()
    cyclic = cyclic_measurement()

    assert cyclic.matched_ion_key == travelling_wave.matched_ion_key
    assert cyclic.calibration_group != travelling_wave.calibration_group
    assert travelling_wave.calibration_group.cyclic_passes is None
    assert cyclic.calibration_group.cyclic_passes is not None


def test_two_cyclic_values_produced_the_same_way_do_share_a_calibration_group():
    # The positive direction, and the one that catches a group made so fine that
    # it holds one record. A group that can never hold two records makes every
    # within-platform agreement check pass without ever firing.
    one = cyclic_measurement(ccs=300.0)
    other = cyclic_measurement(ccs=305.0, source="another laboratory")

    assert one.calibration_group == other.calibration_group


def test_the_calibration_group_prints_its_pass_count_so_two_cyclic_groups_can_be_told_apart():
    one_pass = str(cyclic_measurement(cyclic=cyclic_settings()).calibration_group)
    six_passes = str(cyclic_measurement(cyclic=multipass()).calibration_group)
    no_count = str(
        cyclic_measurement(
            cyclic=CyclicSettings(passes=None, pass_mode=PassMode.MULTIPASS, wrap_around=False)
        ).calibration_group
    )

    assert "1 pass/single_pass" in one_pass
    assert "6 pass/multipass" in six_passes
    # A count the source never gave prints as a question mark, not as a one.
    assert "? pass/multipass" in no_count
    assert len({one_pass, six_passes, no_count}) == 3
    assert "pass" not in str(measurement().calibration_group)


def test_a_cyclic_settings_difference_outside_the_pass_count_does_not_split_the_group():
    # arrival_time_correction is free text, deliberately, and is NOT in the
    # group: two laboratories describing the same correction in different words
    # would otherwise land in two groups and never be compared.
    described_one_way = cyclic_measurement(
        cyclic=multipass(arrival_time_correction="dead time subtracted")
    )
    described_another = cyclic_measurement(
        cyclic=multipass(arrival_time_correction="time outside the separation region removed")
    )

    assert described_one_way.calibration_group == described_another.calibration_group


# --- what a cyclic record stores --------------------------------------------------------


@pytest.mark.parametrize("field", sorted(FULLY_STATED), ids=sorted(FULLY_STATED))
def test_every_cyclic_field_a_source_states_round_trips_through_the_measurement(field):
    settings = CyclicSettings(**FULLY_STATED)
    record = cyclic_measurement(cyclic=settings)

    assert getattr(settings, field) == FULLY_STATED[field]
    assert getattr(record.cyclic, field) == FULLY_STATED[field]
    assert record.cyclic == settings


def test_a_cyclic_setting_the_source_did_not_state_is_left_null_rather_than_filled_in():
    # Guessing a travelling-wave velocity from what is usual would be inventing
    # an instrument setting, which reads downstream exactly like a measured one.
    bare = CyclicSettings()

    assert bare.passes is None
    assert bare.pass_mode is PassMode.UNSTATED
    assert bare.effective_path_length_m is None
    assert bare.tw_velocity_m_per_s is None
    assert bare.tw_height_v is None
    assert bare.arrival_time_correction is None
    assert bare.wrap_around is None


@pytest.mark.parametrize(
    "field", ["effective_path_length_m", "tw_velocity_m_per_s", "tw_height_v"]
)
@pytest.mark.parametrize("value", [0.0, -1.0], ids=["zero", "negative"])
def test_a_cyclic_instrument_setting_that_cannot_be_a_measurement_is_refused(field, value):
    with pytest.raises(ValidationError, match="greater than 0"):
        cyclic_settings(**{field: value})
    assert getattr(cyclic_settings(**{field: 12.5}), field) == 12.5


@pytest.mark.parametrize(
    "field", ["effective_path_length_m", "tw_velocity_m_per_s", "tw_height_v"]
)
def test_a_cyclic_instrument_setting_may_not_be_nan_or_infinite(field):
    for value in (float("nan"), float("inf")):
        with pytest.raises(ValidationError):
            cyclic_settings(**{field: value})


@pytest.mark.parametrize(
    "field, value",
    [
        pytest.param("passes", "6", id="a_pass_count_written_as_text"),
        pytest.param("wrap_around", "true", id="a_wrap_around_flag_written_as_text"),
        pytest.param("effective_path_length_m", "6.5", id="a_path_length_written_as_text"),
    ],
)
def test_a_cyclic_field_is_not_silently_coerced_from_text(field, value):
    # A path length that arrives as text has come from somewhere unvalidated, and
    # a value read loosely here is a value nothing downstream can question.
    with pytest.raises(ValidationError):
        CyclicSettings(**{**FULLY_STATED, field: value})


def test_a_misspelt_cyclic_field_is_refused_rather_than_ignored():
    with pytest.raises(ValidationError, match="extra_forbidden|not permitted|Extra inputs"):
        CyclicSettings(passes=1, pass_mode=PassMode.SINGLE_PASS, wrap_arround=False)


# --- cyclic settings belong to cyclic records, and to nothing else -----------------------


def test_a_cyclic_record_that_states_no_cyclic_settings_is_refused():
    # Without them a cyclic value cannot be told apart from a single-pass
    # travelling-wave value, and the two would silently pool.
    with pytest.raises(ValidationError, match="must state its cyclic settings"):
        measurement(ims_type=IMSType.CYCLIC)
    assert cyclic_measurement().cyclic is not None


@pytest.mark.parametrize("platform", sorted(NON_CYCLIC_PLATFORMS), ids=sorted(NON_CYCLIC_PLATFORMS))
def test_cyclic_settings_are_refused_on_a_platform_that_has_no_passes(platform):
    fields = NON_CYCLIC_PLATFORMS[platform]
    with pytest.raises(ValidationError, match="cyclic settings apply only to CYCLIC records"):
        measurement(**fields, cyclic=cyclic_settings())
    assert measurement(**fields).cyclic is None


@pytest.mark.parametrize("count", [0, -1], ids=["zero", "negative"])
def test_a_pass_count_below_one_is_refused_because_passes_are_counted_from_one(count):
    # A zero-pass cyclic measurement is not a measurement, and a count that can
    # be zero makes "passes > 1" the only thing standing between one pass and six.
    with pytest.raises(ValidationError, match="greater than or equal to 1"):
        cyclic_settings(passes=count)
    assert cyclic_settings(passes=1).passes == 1


# --- the count and the mode have to agree -----------------------------------------------


@pytest.mark.parametrize(
    "passes, mode",
    [
        pytest.param(1, PassMode.MULTIPASS, id="one_pass_called_multipass"),
        pytest.param(6, PassMode.SINGLE_PASS, id="six_passes_called_single_pass"),
        pytest.param(2, PassMode.SINGLE_PASS, id="two_passes_called_single_pass"),
    ],
)
def test_a_pass_count_and_a_pass_mode_that_contradict_each_other_are_refused(passes, mode):
    # One of the two is wrong and there is no way to tell which. Picking one
    # would change what the value means, so the record is refused instead.
    with pytest.raises(ValidationError, match="one of the two is wrong"):
        CyclicSettings(passes=passes, pass_mode=mode, wrap_around=False)


@pytest.mark.parametrize(
    "passes, mode",
    [
        pytest.param(1, PassMode.SINGLE_PASS, id="one_pass_called_single_pass"),
        pytest.param(6, PassMode.MULTIPASS, id="six_passes_called_multipass"),
        pytest.param(2, PassMode.MULTIPASS, id="two_passes_called_multipass"),
        pytest.param(None, PassMode.SINGLE_PASS, id="a_mode_with_no_count"),
        pytest.param(None, PassMode.MULTIPASS, id="a_multipass_mode_with_no_count"),
    ],
)
def test_a_pass_count_and_a_pass_mode_that_agree_are_accepted(passes, mode):
    settings = CyclicSettings(passes=passes, pass_mode=mode, wrap_around=False)
    assert settings.passes == passes
    assert settings.pass_mode is mode


@pytest.mark.parametrize("count", [1, 2, 6], ids=["one", "two", "six"])
def test_an_unstated_pass_mode_is_allowed_with_any_count_because_the_count_then_decides(count):
    # UNSTATED is a positive record of what the paper does not say, not a filler,
    # and it contradicts nothing: a stated count already answers the question.
    settings = CyclicSettings(passes=count, pass_mode=PassMode.UNSTATED, wrap_around=False)

    assert settings.pass_mode is PassMode.UNSTATED
    assert settings.is_multipass is (count > 1)


# --- is_multipass -----------------------------------------------------------------------


@pytest.mark.parametrize(
    "passes, mode, expected",
    [
        pytest.param(1, PassMode.SINGLE_PASS, False, id="one_pass_stated_both_ways_is_not_multipass"),
        pytest.param(1, PassMode.UNSTATED, False, id="one_pass_with_no_mode_is_not_multipass"),
        pytest.param(None, PassMode.SINGLE_PASS, False, id="single_pass_by_mode_alone"),
        pytest.param(None, PassMode.UNSTATED, False, id="neither_stated_is_not_a_claim_of_multipass"),
        pytest.param(2, PassMode.MULTIPASS, True, id="two_passes_is_multipass"),
        pytest.param(6, PassMode.UNSTATED, True, id="the_count_decides_when_the_mode_is_unstated"),
        pytest.param(None, PassMode.MULTIPASS, True, id="the_mode_decides_when_there_is_no_count"),
    ],
)
def test_is_multipass_is_true_by_count_and_by_mode_and_false_for_one_pass(passes, mode, expected):
    assert CyclicSettings(passes=passes, pass_mode=mode, wrap_around=False).is_multipass is expected


# --- training blockers, in both directions ----------------------------------------------

# Each case: settings that block, and a fragment of the reason they give. Used
# for the blocker itself, for the prefix, for the measurement's own list, and for
# the gate, so that a blocker cannot be reported in one place and lost in another.
BLOCKING_SETTINGS = {
    "neither_a_count_nor_a_mode": (
        dict(passes=None, pass_mode=PassMode.UNSTATED, wrap_around=False),
        "neither a pass count nor a pass mode",
    ),
    "multipass_silent_on_wrap_around": (
        dict(passes=6, pass_mode=PassMode.MULTIPASS, wrap_around=None),
        "whether wrap-around occurred",
    ),
    "wrap_around_with_no_correction": (
        dict(passes=6, pass_mode=PassMode.MULTIPASS, wrap_around=True),
        "no arrival-time correction is recorded",
    ),
}


def test_a_cyclic_value_stating_neither_a_pass_count_nor_a_pass_mode_blocks_training():
    silent = CyclicSettings(passes=None, pass_mode=PassMode.UNSTATED, wrap_around=False)
    assert any("neither a pass count nor a pass mode" in b for b in silent.training_blockers())
    # Either one alone is enough to say how far the ion travelled.
    assert CyclicSettings(passes=1, pass_mode=PassMode.UNSTATED, wrap_around=False).training_blockers() == []
    assert CyclicSettings(passes=None, pass_mode=PassMode.SINGLE_PASS, wrap_around=False).training_blockers() == []


def test_a_multipass_value_that_never_says_whether_ions_lapped_blocks_training():
    # On a multipass value that silence is a gap, not a no: the longer the
    # flight, the likelier a fast ion laps a slow one.
    silent = multipass(wrap_around=None)
    assert any("whether wrap-around occurred" in b for b in silent.training_blockers())
    assert multipass(wrap_around=False).training_blockers() == []


def test_a_single_pass_value_that_does_not_mention_wrap_around_does_not_block():
    # One pass cannot lap anything, so silence there is not a gap. This is also
    # the assertion that fails if one pass is ever read as several.
    settings = cyclic_settings(wrap_around=None)
    assert settings.is_multipass is False
    assert settings.training_blockers() == []


def test_wrap_around_with_no_arrival_time_correction_blocks_training():
    uncorrected = multipass(wrap_around=True)
    assert any("no arrival-time correction is recorded" in b for b in uncorrected.training_blockers())


def test_wrap_around_with_an_arrival_time_correction_recorded_does_not_block():
    # The negative direction: wrap-around is not itself disqualifying. What
    # disqualifies is wrap-around that nobody corrected for.
    corrected = multipass(
        wrap_around=True, arrival_time_correction="unwrapped against a single-pass reference"
    )
    assert corrected.training_blockers() == []


def test_a_clean_single_pass_cyclic_record_has_no_cyclic_blocker_at_all():
    record = cyclic_measurement()
    assert record.cyclic.training_blockers() == []
    assert record.training_blockers() == []


@pytest.mark.parametrize("case", sorted(BLOCKING_SETTINGS), ids=sorted(BLOCKING_SETTINGS))
def test_a_cyclic_settings_blocker_reaches_the_measurements_own_training_blockers(case):
    # A blocker that only exists on CyclicSettings is a blocker nothing consults:
    # the gate reads the measurement, not the settings hanging off it.
    fields, fragment = BLOCKING_SETTINGS[case]
    record = cyclic_measurement(cyclic=CyclicSettings(**fields))

    assert any(fragment in blocker for blocker in record.training_blockers())


@pytest.mark.parametrize("case", sorted(BLOCKING_SETTINGS), ids=sorted(BLOCKING_SETTINGS))
def test_a_cyclic_settings_blocker_is_prefixed_so_a_reader_can_see_where_it_came_from(case):
    fields, fragment = BLOCKING_SETTINGS[case]
    record = cyclic_measurement(cyclic=CyclicSettings(**fields))

    reported = [blocker for blocker in record.training_blockers() if fragment in blocker]
    assert reported
    assert all(blocker.startswith("cyclic settings: ") for blocker in reported)


@pytest.mark.parametrize("case", sorted(BLOCKING_SETTINGS), ids=sorted(BLOCKING_SETTINGS))
def test_the_gate_refuses_a_cyclic_record_whose_settings_block_and_says_which_settings(case):
    fields, fragment = BLOCKING_SETTINGS[case]
    record = licensed_cyclic(cyclic=CyclicSettings(**fields))

    # Both halves of the message: that it came from the cyclic settings, and
    # which of them. "cyclic settings" alone would not say what to go and read.
    with pytest.raises(TrainingGateError, match="cyclic settings"):
        assert_trainable(record)
    with pytest.raises(TrainingGateError, match=re.escape(fragment)):
        assert_trainable(record)


# --- cyclic settings are not a component record ------------------------------------------


def test_cyclic_settings_are_not_a_component_record_because_they_can_state_no_reuse_status():
    # THE ANTIBODY GATE BUG, in the place it would come back. AntibodyIdentity
    # was listed in component_records(), the gate asked it for a reuse status it
    # has no field for, and every antibody and ADC record was refused: the whole
    # biopharmaceutical layer, closed by one tuple entry. CyclicSettings has the
    # same shape, conditions of a measurement rather than a separately sourced
    # record, so it must stay out of that tuple for the same reason.
    record = cyclic_measurement()
    parts = record.component_records()

    assert record.cyclic is not None
    assert record.cyclic not in parts
    assert not any(isinstance(part, CyclicSettings) for part in parts)
    assert parts == (record.analyte,)


def test_the_gate_would_refuse_cyclic_settings_outright_which_is_why_they_are_not_a_component():
    # The other half, and the one that says what the damage would be rather than
    # only that the tuple is short. If CyclicSettings ever enters
    # component_records(), THIS is the refusal every cyclic record would get.
    settings = cyclic_settings()

    assert not hasattr(settings, "reuse_status")
    with pytest.raises(LicenceGateError, match="cannot tell the reuse status of a CyclicSettings"):
        assert_trainable(settings)


def test_a_cyclic_record_with_a_valid_licence_clears_the_gate():
    # The assertion the antibody bug needed and did not have: it is not enough
    # that the settings are absent from the tuple, a clean cyclic record has to
    # actually pass.
    assert_trainable(licensed_cyclic())
    assert_trainable(licensed_cyclic(cyclic=multipass()))
    assert_trainable(
        licensed_cyclic(cyclic=multipass(wrap_around=True, arrival_time_correction="unwrapped"))
    )


def test_a_cyclic_record_whose_licence_is_unverified_is_still_refused():
    # Default deny does not stop applying because the record is cyclic.
    with pytest.raises(LicenceGateError, match="unverified"):
        assert_trainable(licensed_cyclic(reuse_status=ReuseStatus.UNVERIFIED))
