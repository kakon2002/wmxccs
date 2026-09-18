"""The licence gate: what it refuses, how it refuses, and what it refuses to assume.

Three properties are worth more than any amount of coverage here, and they are
what most of this file is about.

IT RAISES, IT DOES NOT RETURN FALSE. A loader that filters rows on a boolean
drops refused rows silently and the loss only surfaces if somebody later audits
row counts. An exception stops the load on the offending row, so every test
below asserts a raise and never a return value.

THE ERROR CLASS IS PART OF THE ANSWER. A licence fault, an unbacked claim and
an ordinary blocker have three different remedies - read the terms, record the
evidence, fix the record - so the gate raises three different classes and the
worst one wins. Reporting an unbacked claim as a licence refusal would send
somebody to relabel a record when what was missing was a registry entry.

IT FAILS CLOSED. A record whose status, blockers or components cannot be READ is
refused rather than assumed. None of those four branches can be reached with a
real CCSMeasurement, because pydantic will not build one with an accessor that
raises, so they are tested with deliberately broken stubs. They are the reason
the gate can be trusted with a file nobody has looked at, and without these
tests they are the least exercised code in the package.
"""

from __future__ import annotations

import pytest
from conftest import SOURCE, antibody, glycan, measurement

from wmxccs.identity import Derivatisation, DriftGas, FoldingState, ReducingEndLabel
from wmxccs.licensing import (
    LicenceGateError,
    ReuseStatus,
    TrainingGateError,
    UnbackedClaimError,
    as_reuse_status,
    assert_trainable,
    can_train_commercial,
    declares_synthetic,
)

# A DOI the registry holds, with the status it is recorded at, and the exact
# provenance string the seed transcription writes into an analyte's source.
# Taken from sources.py rather than invented: a DOI nobody recorded would test
# the unbacked path while pretending to test the backed one.
REGISTERED_DOI = "10.1039/c6cc06247d"
REGISTERED_SEED = "wmxccs seed transcription: Struwe 2016 Chem Commun ESI Table S1"
# Recorded as non_commercial_no_derivatives, so a row claiming an open licence
# on it claims more than the registry says.
RESTRICTED_DOI = "10.1038/s41467-025-67069-w"
UNREGISTERED_DOI = "10.1000/nobody-has-read-this"

# Every status the gate must refuse, with a phrase the refusal has to carry.
# Hard-coded rather than read back out of reuse.py: a list computed from the
# tier the gate consults would shrink in step with a mistake there and take
# these tests with it.
REFUSED_STATUSES = (
    (ReuseStatus.OPEN_SHARE_ALIKE, "by policy"),
    (ReuseStatus.NON_COMMERCIAL, "do not permit commercial use"),
    (ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES, "each clause blocks training on its own"),
    (ReuseStatus.ACADEMIC_ONLY, "academic or research purposes"),
    (ReuseStatus.UNVERIFIED, "nobody has checked its terms yet"),
    (ReuseStatus.EXCLUDED, "deliberately kept out of the platform"),
)


def defined_glycan(**overrides):
    """conftest's glycan() with its label and derivatisation stated.

    The builder leaves both 'unknown', which blocks training on its own and is
    right for the loader's tests. A gate test needs an analyte whose only
    interesting property is its licence, or every refusal below would pass for
    the wrong reason.
    """
    fields = dict(reducing_end_label=ReducingEndLabel.NATIVE, derivatisation=Derivatisation.UNDERIVATISED)
    fields.update(overrides)
    return glycan(**fields)


def clean(**overrides):
    """A measurement with nothing wrong with it: synthetic, fully defined, no blockers."""
    fields = dict(analyte=defined_glycan())
    fields.update(overrides)
    return measurement(**fields)


def published(**overrides):
    """A record standing on the registry: an open DOI, and an analyte from the seed it was transcribed into."""
    fields = dict(
        analyte=defined_glycan(source=REGISTERED_SEED, reuse_status=ReuseStatus.OPEN_ATTRIBUTION),
        source="Struwe et al., Chem. Commun. 2016",
        doi=REGISTERED_DOI,
        reuse_status=ReuseStatus.OPEN_ATTRIBUTION,
    )
    fields.update(overrides)
    return measurement(**fields)


# --- the deliberately broken stubs ------------------------------------------------------
#
# Not CCSMeasurements, and they cannot be: a validated record cannot be built
# with an accessor that raises, which is exactly why the gate's defensive
# branches would otherwise never run.


class Stub:
    """The smallest thing the gate accepts: a status it can read, and a source to name it by."""

    reuse_status = ReuseStatus.SYNTHETIC_FIXTURE
    source = "a stub in the licence gate tests"


class StublessOfStatus:
    """A record with no reuse_status at all. Deliberately not a subclass of Stub."""

    source = "a stub in the licence gate tests"


class StatusBomb(Stub):
    @property
    def reuse_status(self):  # type: ignore[override]
        raise RuntimeError("the status column was never loaded")


class BlockerBomb(Stub):
    def training_blockers(self):
        raise RuntimeError("the blocker table could not be read")


class ComponentBomb(Stub):
    def component_records(self):
        raise RuntimeError("the component list could not be read")


class DoiBomb(Stub):
    @property
    def doi(self):
        raise RuntimeError("the DOI column could not be read")


class Refusing(Stub):
    reuse_status = ReuseStatus.UNVERIFIED
    source = "a component nobody has checked the terms of"


class Holder(Stub):
    """A record built from whatever it is handed."""

    def __init__(self, *parts: object) -> None:
        self.parts = parts

    def component_records(self):
        return self.parts


class Loop:
    """A record that is built from another record that is built from it."""

    source = "a stub in the licence gate tests"

    def __init__(self, reuse_status: ReuseStatus = ReuseStatus.OPEN_ATTRIBUTION) -> None:
        self.reuse_status = reuse_status
        self.other: object = None

    def component_records(self):
        return (self.other,)


def loop_of_two(*statuses: ReuseStatus) -> Loop:
    first, second = (Loop(status) for status in statuses)
    first.other, second.other = second, first
    return first


# --- the baseline, so that a refusal below cannot pass for the wrong reason ---------------


def test_a_record_with_nothing_wrong_with_it_passes_the_gate() -> None:
    """The anchor for the whole file.

    Every other test asserts a raise. If the builder produced a record the gate
    refused for some unrelated reason, all of them would pass while testing
    nothing, so this one goes first.
    """
    assert assert_trainable(clean()) is None


def test_a_record_whose_open_licence_the_registry_backs_may_enter_training() -> None:
    assert assert_trainable(published()) is None


@pytest.mark.parametrize("status", (ReuseStatus.OPEN_ATTRIBUTION, ReuseStatus.INTERNAL_PROPRIETARY, ReuseStatus.SYNTHETIC_FIXTURE))
def test_a_bare_trainable_status_passes_because_a_status_alone_makes_no_claim(status: ReuseStatus) -> None:
    assert assert_trainable(status) is None
    assert assert_trainable(status.value) is None


# --- it raises, and what it raises ------------------------------------------------------


@pytest.mark.parametrize(("status", "phrase"), REFUSED_STATUSES)
def test_the_gate_raises_on_every_untrainable_status_and_says_why(status: ReuseStatus, phrase: str) -> None:
    with pytest.raises(LicenceGateError) as refusal:
        assert_trainable(clean(reuse_status=status))
    assert status.value in str(refusal.value)
    assert phrase in str(refusal.value)


def test_a_refusal_tells_the_reader_not_to_filter_the_row_out_quietly() -> None:
    """The remedy belongs in the message: the failure mode this gate exists to stop
    is somebody catching the refusal and dropping the row."""
    with pytest.raises(TrainingGateError, match="do not filter it out quietly"):
        assert_trainable(clean(reuse_status=ReuseStatus.UNVERIFIED))


def test_the_gate_refuses_by_raising_rather_than_by_returning_false() -> None:
    """A boolean answer is the whole failure mode: it can be ignored, and is.

    Asserted as a property of the call, not of the class hierarchy: nothing the
    gate returns for a refused record may be usable as a filter.
    """
    refused = clean(reuse_status=ReuseStatus.EXCLUDED)
    with pytest.raises(TrainingGateError):
        assert_trainable(refused)
    assert assert_trainable(clean()) is None


def test_a_gate_refusal_is_not_a_value_error() -> None:
    """A loader skipping malformed rows catches ValueError. It must not catch these."""
    assert not issubclass(TrainingGateError, ValueError)
    assert issubclass(LicenceGateError, TrainingGateError)
    assert issubclass(UnbackedClaimError, LicenceGateError)
    with pytest.raises(TrainingGateError):
        try:
            assert_trainable(clean(reuse_status=ReuseStatus.UNVERIFIED))
        except ValueError as exc:  # pragma: no cover - only runs if the hierarchy breaks
            raise AssertionError("a licence refusal was catchable as a malformed row") from exc


# --- which class, and in which order ----------------------------------------------------


def test_a_licence_fault_gives_a_licence_gate_error() -> None:
    with pytest.raises(LicenceGateError, match="academic or research purposes"):
        assert_trainable(clean(reuse_status=ReuseStatus.ACADEMIC_ONLY))


@pytest.mark.parametrize(
    ("doi", "phrase"),
    (
        (None, "with no DOI"),
        (UNREGISTERED_DOI, "which has no licence record"),
        (RESTRICTED_DOI, "the record governs"),
    ),
)
def test_an_unbacked_claim_alone_gives_an_unbacked_claim_error(doi: str | None, phrase: str) -> None:
    """Three ways a trainable claim can be unbacked, and all three point at the registry.

    Told apart from a licence fault because the remedy differs: here the status
    may well be right and the evidence for it is missing, so the fix is to read
    the licence and record it, never to relabel the record.
    """
    record = published(doi=doi, reuse_status=ReuseStatus.OPEN_ATTRIBUTION)
    with pytest.raises(UnbackedClaimError) as refusal:
        assert_trainable(record)
    assert phrase in str(refusal.value)


def test_an_unbacked_claim_on_the_analyte_alone_gives_an_unbacked_claim_error() -> None:
    record = published(
        analyte=defined_glycan(source="a spreadsheet somebody sent", reuse_status=ReuseStatus.OPEN_ATTRIBUTION)
    )
    with pytest.raises(UnbackedClaimError, match="neither the row's DOI nor a registered dataset backs it"):
        assert_trainable(record)


@pytest.mark.parametrize(
    ("overrides", "phrase"),
    (
        (dict(drift_gas=DriftGas.UNSTATED), "the gas this CCS refers to is not defined"),
        (dict(adduct="[M+24?]24+", charge=24, ccs=7000.0), "does not name the charge carrier"),
    ),
)
def test_a_blocker_that_is_not_a_licence_gives_a_plain_training_gate_error(overrides: dict, phrase: str) -> None:
    """An undefined gas or an unnamed charge carrier is not a licence problem.

    Raising LicenceGateError here would send somebody to read a licence for a
    record whose licence is fine and whose ION is not defined.
    """
    with pytest.raises(TrainingGateError) as refusal:
        assert_trainable(clean(**overrides))
    assert not isinstance(refusal.value, LicenceGateError)
    assert phrase in str(refusal.value)


def test_a_licence_fault_outranks_an_unbacked_claim_and_an_ordinary_blocker() -> None:
    """All three at once. The licence is the one that must be reported."""
    record = published(
        analyte=defined_glycan(source="a spreadsheet somebody sent", reuse_status=ReuseStatus.NON_COMMERCIAL),
        doi=UNREGISTERED_DOI,
        drift_gas=DriftGas.UNSTATED,
    )
    with pytest.raises(LicenceGateError) as refusal:
        assert_trainable(record)
    assert not isinstance(refusal.value, UnbackedClaimError)
    assert "non_commercial" in str(refusal.value)


def test_an_unbacked_claim_outranks_an_ordinary_blocker() -> None:
    record = published(doi=UNREGISTERED_DOI, drift_gas=DriftGas.UNSTATED)
    with pytest.raises(UnbackedClaimError) as refusal:
        assert_trainable(record)
    assert "which has no licence record" in str(refusal.value)
    # The blocker is still reported; it is the CLASS that is decided by the worst fault.
    assert "the gas this CCS refers to is not defined" in str(refusal.value)


# --- fail closed: what cannot be read is refused, never assumed --------------------------


def test_a_record_whose_reuse_status_cannot_be_read_is_refused_rather_than_assumed() -> None:
    """A raising accessor must not leave the gate as a raw RuntimeError either.

    It would travel straight past every `except TrainingGateError` a caller has
    and abort a whole load instead of refusing one row.
    """
    with pytest.raises(LicenceGateError) as refusal:
        assert_trainable(StatusBomb())
    assert "refused rather than assumed" in str(refusal.value)
    assert "RuntimeError" in str(refusal.value)


def test_a_record_with_no_reuse_status_at_all_is_refused() -> None:
    with pytest.raises(LicenceGateError, match="cannot tell the reuse status"):
        assert_trainable(StublessOfStatus())


@pytest.mark.parametrize("subject", (None, "open_sesame", ""))
def test_a_subject_that_is_no_status_at_all_is_refused(subject: object) -> None:
    with pytest.raises(LicenceGateError, match="unrecognised reuse status"):
        assert_trainable(subject)


def test_a_record_whose_training_blockers_raise_is_refused_as_an_ordinary_blocker() -> None:
    """Unreadable blockers are a refusal, and not a licence one: nothing is known
    to be wrong with the licence."""
    with pytest.raises(TrainingGateError) as refusal:
        assert_trainable(BlockerBomb())
    assert not isinstance(refusal.value, LicenceGateError)
    assert "training blockers could not be checked" in str(refusal.value)


def test_a_record_whose_components_cannot_be_listed_is_refused_as_a_licence_fault() -> None:
    """The subtree is abandoned, so its licences are unchecked, so it is a licence refusal.

    Clearing it instead would let a record launder anything it is built from by
    failing to list it.
    """
    with pytest.raises(LicenceGateError) as refusal:
        assert_trainable(ComponentBomb())
    assert "could not be listed" in str(refusal.value)
    assert "licences are unchecked" in str(refusal.value)


def test_abandoning_one_unreadable_subtree_does_not_abandon_its_siblings() -> None:
    """The walk gives up on the branch it cannot read, not on the record.

    A sibling component with a refusing licence must still be reported, or one
    unreadable part of a record would hide every other fault in it.
    """
    with pytest.raises(LicenceGateError) as refusal:
        assert_trainable(Holder(ComponentBomb(), Refusing()))
    assert "could not be listed" in str(refusal.value)
    assert "unverified" in str(refusal.value)


def test_a_claim_that_cannot_be_checked_is_refused_rather_than_waved_through() -> None:
    """The claim check is wrapped like the rest: an exception inside it is a refusal.

    A record whose DOI cannot be read has a claim nobody can trace, which is the
    same position as a claim nothing backs.
    """
    with pytest.raises(UnbackedClaimError) as refusal:
        assert_trainable(DoiBomb())
    assert "claim could not be checked" in str(refusal.value)
    assert "RuntimeError" in str(refusal.value)


# --- components: a record is only as clear as what it is built from ----------------------


def test_a_measurement_is_refused_when_its_analyte_carries_a_refusing_status() -> None:
    """The record's own licence is fine and it is refused anyway.

    A CCS value from an open paper does not clear an identity taken from a
    restricted one, and a gate that only read the outermost status would admit
    exactly that.
    """
    record = clean(analyte=defined_glycan(reuse_status=ReuseStatus.NON_COMMERCIAL))
    assert can_train_commercial(record.reuse_status) is True
    with pytest.raises(LicenceGateError) as refusal:
        assert_trainable(record)
    assert "GlycanAnalyte" in str(refusal.value)
    assert "non_commercial" in str(refusal.value)


@pytest.mark.parametrize("status", (ReuseStatus.UNVERIFIED, ReuseStatus.EXCLUDED, ReuseStatus.ACADEMIC_ONLY))
def test_every_refusing_status_on_a_component_refuses_the_whole_record(status: ReuseStatus) -> None:
    with pytest.raises(LicenceGateError, match=status.value):
        assert_trainable(clean(analyte=defined_glycan(reuse_status=status)))


def test_a_blocker_on_a_component_refuses_the_whole_record() -> None:
    """conftest's glycan() leaves its label unknown, which is what a real
    unreported record looks like, and the measurement around it is spotless."""
    with pytest.raises(TrainingGateError) as refusal:
        assert_trainable(measurement())
    assert not isinstance(refusal.value, LicenceGateError)
    assert "reducing-end label is 'unknown'" in str(refusal.value)


# --- declares_synthetic: a fixture must never be countable as data -----------------------


def test_a_synthetic_declaration_is_found_on_a_component_not_only_at_the_top_level() -> None:
    """The count that matters: a record claiming a real licence whose analyte was
    built by a test. Reading the outermost status alone would present it as data."""
    record = published(analyte=defined_glycan())  # the builder's analyte is a synthetic fixture
    assert as_reuse_status(record.reuse_status) is ReuseStatus.OPEN_ATTRIBUTION
    assert declares_synthetic(record) is True


def test_a_record_that_declares_itself_synthetic_is_counted_at_the_top_level_too() -> None:
    assert declares_synthetic(clean()) is True
    assert declares_synthetic(ReuseStatus.SYNTHETIC_FIXTURE) is True
    assert declares_synthetic(ReuseStatus.SYNTHETIC_FIXTURE.value) is True


def test_a_record_with_no_synthetic_declaration_anywhere_in_it_is_not_counted() -> None:
    assert declares_synthetic(published()) is False
    assert declares_synthetic(ReuseStatus.UNVERIFIED) is False
    assert declares_synthetic("open_sesame") is False
    assert declares_synthetic(None) is False


def test_a_record_whose_components_cannot_be_listed_counts_as_synthetic() -> None:
    """Fail closed in the direction that matters here: label it rather than
    present something unreadable as real data."""
    unreadable = ComponentBomb()
    unreadable.reuse_status = ReuseStatus.OPEN_ATTRIBUTION  # type: ignore[misc]
    assert declares_synthetic(unreadable) is True


def test_a_cycle_between_two_records_does_not_hang_the_synthetic_check() -> None:
    assert declares_synthetic(loop_of_two(ReuseStatus.OPEN_ATTRIBUTION, ReuseStatus.OPEN_ATTRIBUTION)) is False
    assert declares_synthetic(loop_of_two(ReuseStatus.OPEN_ATTRIBUTION, ReuseStatus.SYNTHETIC_FIXTURE)) is True


# --- the re-export: the gate and the statuses are one interface --------------------------


def test_the_statuses_are_importable_from_the_gate_and_are_the_same_objects() -> None:
    """licensing re-exports reuse deliberately. A second ReuseStatus type would
    make `is` comparisons across the two modules silently false."""
    from wmxccs import reuse

    assert ReuseStatus is reuse.ReuseStatus
    assert as_reuse_status is reuse.as_reuse_status
    assert can_train_commercial is reuse.can_train_commercial


# --- a refusal this file believes is wrong -----------------------------------------------


# Was an xfail against a real defect: AntibodyIdentity was listed as a component
# record and carries no reuse status, so the gate refused every antibody and ADC
# measurement whatever its licence. It is the structured form of the enclosing
# analyte's identity, not a separately sourced record, and is no longer offered.
def test_a_native_antibody_measurement_that_declares_itself_synthetic_may_enter_training() -> None:
    """The biopharmaceutical layer is the platform's differentiator and no record
    in it can currently pass the gate. Reported as a suspected defect in identity.py
    rather than worked around here."""
    record = measurement(
        analyte=antibody(folding_state=FoldingState.NATIVE, source=SOURCE),
        adduct="[M+24H]24+",
        charge=24,
        ccs=7000.0,
    )
    assert assert_trainable(record) is None
