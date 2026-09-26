"""Licence gate: default deny, and a refused record stops the load instead of vanishing."""

import subprocess
import sys

import pytest

from wmxglycan.licensing import (
    DEFAULT_REUSE_STATUS,
    LicenceGateError,
    ReuseStatus,
    TrainingGateError,
    UnbackedClaimError,
    assert_trainable,
    can_redistribute,
    can_train_commercial,
    is_inference_only,
)
from wmxglycan.sources import STRUWE_2016
from wmxglycan.sugarbase import dataset_version

S = ReuseStatus

# (trainable commercially, redistributable, inference-only) for every status.
POLICY = {
    S.OPEN_ATTRIBUTION: (True, True, False),
    S.OPEN_SHARE_ALIKE: (False, True, True),  # kept out of training in V1 by policy
    S.NON_COMMERCIAL: (False, False, True),
    # No permitted use here at all, so no tier: the no-derivatives clause blocks
    # even the inference-side reference use that plain non-commercial allows.
    S.NON_COMMERCIAL_NO_DERIVATIVES: (False, False, False),
    S.ACADEMIC_ONLY: (False, False, True),
    S.INTERNAL_PROPRIETARY: (True, False, False),
    S.SYNTHETIC_FIXTURE: (True, False, False),  # a test's declaration; never a real record, see below
    S.UNVERIFIED: (False, False, False),
    S.EXCLUDED: (False, False, False),
}
TRAINABLE = [s for s in ReuseStatus if POLICY[s][0]]
NOT_TRAINABLE = [s for s in ReuseStatus if not POLICY[s][0]]


def test_default_is_unverified():
    assert DEFAULT_REUSE_STATUS is S.UNVERIFIED


def test_every_status_has_a_policy_row():
    # A new status added without a decision here fails this test.
    assert set(POLICY) == set(ReuseStatus)


@pytest.mark.parametrize("status", list(ReuseStatus))
def test_policy(status):
    assert (can_train_commercial(status), can_redistribute(status), is_inference_only(status)) == POLICY[status]


@pytest.mark.parametrize("status", list(ReuseStatus))
def test_tiers_are_consistent(status):
    # At most one usable tier, and nothing without a permitted use is passed on.
    assert not (can_train_commercial(status) and is_inference_only(status))
    if can_redistribute(status):
        assert can_train_commercial(status) or is_inference_only(status)


def test_no_derivatives_is_stricter_than_non_commercial_and_says_which_clauses():
    # Folding CC BY-NC-ND into non_commercial would lose the ND clause, which is
    # the stricter of the two. The refusal names both, so relabelling is visibly
    # not a remedy.
    assert not can_train_commercial(S.NON_COMMERCIAL_NO_DERIVATIVES)
    assert not is_inference_only(S.NON_COMMERCIAL_NO_DERIVATIVES)  # unlike plain non_commercial
    assert is_inference_only(S.NON_COMMERCIAL)
    with pytest.raises(LicenceGateError) as caught:
        assert_trainable(S.NON_COMMERCIAL_NO_DERIVATIVES)
    message = str(caught.value)
    assert "derivative works" in message and "commercial use" in message
    assert "not cleared for training" not in message  # not the unlisted-status fallback


def test_share_alike_is_not_trainable_in_v1():
    assert not can_train_commercial(S.OPEN_SHARE_ALIKE)
    with pytest.raises(LicenceGateError, match="open_share_alike"):
        assert_trainable(S.OPEN_SHARE_ALIKE)


def test_share_alike_refusal_blames_policy_not_the_licence():
    # A message blaming the licence would be false, and would invite relabelling the record.
    with pytest.raises(LicenceGateError) as caught:
        assert_trainable(S.OPEN_SHARE_ALIKE)
    message = str(caught.value)
    assert "by policy" in message and "allows commercial use" in message
    assert "do not permit" not in message


@pytest.mark.parametrize("status", NOT_TRAINABLE)
def test_every_refused_status_gives_its_own_reason(status):
    with pytest.raises(LicenceGateError) as caught:
        assert_trainable(status)
    assert "not cleared for training" not in str(caught.value)  # the fallback for an unlisted status


@pytest.mark.parametrize("status", TRAINABLE)
def test_gate_passes_trainable_statuses(status):
    assert_trainable(status)
    assert_trainable(status.value)


@pytest.mark.parametrize("status", NOT_TRAINABLE)
def test_gate_raises_for_every_other_status(status):
    with pytest.raises(LicenceGateError, match=status.value):
        assert_trainable(status)
    with pytest.raises(LicenceGateError, match=status.value):
        assert_trainable(status.value)


@pytest.mark.parametrize("junk", ["", "ACADEMIC_ONLY", "cc-by", None, True, 1, object(), {"reuse_status": "open_attribution"}])
def test_gate_raises_when_status_is_unrecognisable(junk):
    with pytest.raises(LicenceGateError):
        assert_trainable(junk)


@pytest.mark.parametrize("junk", ["", "ACADEMIC_ONLY", "cc-by"])
def test_predicates_reject_unknown_strings(junk):
    for predicate in (can_train_commercial, can_redistribute, is_inference_only):
        with pytest.raises(ValueError):
            predicate(junk)


@pytest.mark.parametrize("junk", [None, True, 1])
def test_predicates_reject_non_strings(junk):
    for predicate in (can_train_commercial, can_redistribute, is_inference_only):
        with pytest.raises(TypeError):
            predicate(junk)


def test_gate_errors_are_not_swallowed_by_row_level_except_clauses():
    # Loaders commonly `except ValueError` to skip malformed rows.
    assert issubclass(LicenceGateError, TrainingGateError)
    for error in (TrainingGateError, LicenceGateError):
        assert not issubclass(error, (ValueError, TypeError))


def test_gate_still_raises_under_python_O():
    code = (
        "from wmxglycan.licensing import LicenceGateError, assert_trainable\n"
        "try:\n"
        "    assert_trainable('academic_only')\n"
        "except LicenceGateError:\n"
        "    raise SystemExit(0)\n"
        "raise SystemExit(1)\n"
    )
    result = subprocess.run([sys.executable, "-O", "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


class StubRecord:
    """A record whose gate hooks can be made to fail."""

    source = "synthetic test fixture"

    def __init__(self, reuse_status, blockers=(), components=()):
        self.reuse_status = reuse_status
        self._blockers = blockers
        self._components = components

    def training_blockers(self):
        if isinstance(self._blockers, Exception):
            raise self._blockers
        return list(self._blockers)

    def component_records(self):
        if isinstance(self._components, Exception):
            raise self._components
        return tuple(self._components)


@pytest.mark.parametrize("failure", [ValueError("ambiguous truth value"), TypeError("unhashable"), KeyError("x")])
def test_a_record_whose_blockers_cannot_be_checked_is_refused_not_crashed(failure):
    # A raw ValueError or TypeError could be swallowed by a loader that skips malformed rows.
    with pytest.raises(TrainingGateError, match="could not be checked") as caught:
        assert_trainable(StubRecord(S.SYNTHETIC_FIXTURE, blockers=failure))
    assert not isinstance(caught.value, LicenceGateError)
    with pytest.raises(LicenceGateError, match="academic_only"):
        assert_trainable(StubRecord(S.ACADEMIC_ONLY, blockers=failure))


def test_a_record_whose_components_cannot_be_listed_is_refused_on_licence():
    with pytest.raises(LicenceGateError, match="could not be listed"):
        assert_trainable(StubRecord(S.SYNTHETIC_FIXTURE, components=AttributeError("glycan")))


def test_components_are_checked_through_the_hook():
    with pytest.raises(LicenceGateError, match="academic_only"):
        assert_trainable(StubRecord(S.SYNTHETIC_FIXTURE, components=[StubRecord(S.ACADEMIC_ONLY)]))
    assert_trainable(StubRecord(S.SYNTHETIC_FIXTURE, components=[StubRecord(S.INTERNAL_PROPRIETARY)]))


# --- a trainable status is a claim, and the gate checks it against the registry ------


class Measured(StubRecord):
    """A record with provenance of its own, the way a CCSMeasurement has a DOI and a source."""

    def __init__(self, reuse_status, doi=None, source="synthetic test fixture", components=()):
        super().__init__(reuse_status, components=components)
        self.doi = doi
        self.source = source


def a_structure(reuse_status, source="synthetic test fixture"):
    structure = StubRecord(reuse_status)
    structure.source = source
    return structure


def test_an_open_claim_with_no_doi_is_refused_as_unbacked():
    # The bypass: a record built in code, never through the loader, claiming a
    # licence nobody recorded. The gate refuses it; no convention is involved.
    with pytest.raises(UnbackedClaimError, match="no DOI") as caught:
        assert_trainable(Measured(S.OPEN_ATTRIBUTION))
    assert isinstance(caught.value, LicenceGateError)


def test_an_open_claim_for_an_unrecorded_doi_is_refused():
    with pytest.raises(UnbackedClaimError, match="no licence record"):
        assert_trainable(Measured(S.OPEN_ATTRIBUTION, doi="10.1000/no-record"))


def test_an_open_claim_for_a_recorded_doi_passes():
    assert_trainable(Measured(S.OPEN_ATTRIBUTION, doi=STRUWE_2016.doi))


def test_an_internal_claim_needs_no_doi_and_no_record():
    assert_trainable(Measured(S.INTERNAL_PROPRIETARY))


def test_a_synthetic_fixture_needs_no_record():
    # A test's declaration that the record is invented. No DOI enters the registry for it.
    assert_trainable(Measured(S.SYNTHETIC_FIXTURE))
    assert_trainable(Measured(S.SYNTHETIC_FIXTURE, doi="10.1000/invented-with-the-rest"))
    assert_trainable(StubRecord(S.SYNTHETIC_FIXTURE))  # no provenance at all
    assert_trainable(Measured(S.SYNTHETIC_FIXTURE, components=[a_structure(S.SYNTHETIC_FIXTURE)]))


def test_a_synthetic_fixture_cannot_cite_a_recorded_paper():
    # A record citing a paper whose licence is on record is that paper's record,
    # and the record governs. Synthetic does not launder a real DOI.
    with pytest.raises(UnbackedClaimError, match="the record governs"):
        assert_trainable(Measured(S.SYNTHETIC_FIXTURE, doi=STRUWE_2016.doi))


def test_a_bare_status_has_no_claim_to_check():
    # A status alone is not a record; nothing can be trained from it. Its tier is all there is.
    assert_trainable(S.OPEN_ATTRIBUTION)
    assert_trainable("open_attribution")


def test_a_component_is_backed_by_the_enclosing_doi_only_when_it_names_the_same_source():
    paper = "the paper's table"
    assert_trainable(
        Measured(S.OPEN_ATTRIBUTION, doi=STRUWE_2016.doi, source=paper, components=[a_structure(S.OPEN_ATTRIBUTION, paper)])
    )
    # A registered DOI on the record does not launder a structure that names something else.
    with pytest.raises(UnbackedClaimError, match="named as the measurement's own source"):
        assert_trainable(
            Measured(
                S.OPEN_ATTRIBUTION,
                doi=STRUWE_2016.doi,
                source=paper,
                components=[a_structure(S.OPEN_ATTRIBUTION, "some other database")],
            )
        )


def test_a_component_whose_own_doi_is_absent_falls_to_the_structure_rule():
    # A component carrying a doi attribute set to None names no citation, so it
    # must be judged as a structure - by its source - rather than refused for
    # having no DOI of its own. A real GlycanStructure has no doi field at all,
    # but the gate accepts any record, so the branch has to be right for those too.
    inner = a_structure(S.OPEN_ATTRIBUTION, dataset_version().provenance)
    inner.doi = None
    assert_trainable(Measured(S.OPEN_ATTRIBUTION, doi=STRUWE_2016.doi, source="the paper", components=[inner]))


def test_a_licence_refusal_outranks_an_unbacked_claim():
    # Both can hold at once: the measurement's status is refused outright while
    # the structure's claim is merely unbacked. A caller catching
    # UnbackedClaimError must not be handed a record whose licence forbids
    # training, because the remedy that class implies - go and record the
    # licence - is not the remedy.
    record = Measured(
        S.ACADEMIC_ONLY, source="the paper", components=[a_structure(S.OPEN_ATTRIBUTION, "nowhere")]
    )
    with pytest.raises(LicenceGateError) as caught:
        assert_trainable(record)
    assert not isinstance(caught.value, UnbackedClaimError)
    assert "academic_only" in str(caught.value)


def test_a_structure_standing_alone_is_backed_by_its_dataset_or_by_nothing():
    assert_trainable(a_structure(S.OPEN_ATTRIBUTION, dataset_version().provenance))
    with pytest.raises(UnbackedClaimError, match="neither the row's DOI nor a registered dataset"):
        assert_trainable(a_structure(S.OPEN_ATTRIBUTION))


def test_a_claim_that_cannot_be_checked_is_refused_not_crashed():
    # The registry lookup would raise on a DOI that is not text; that becomes a refusal, not an escape.
    with pytest.raises(UnbackedClaimError, match="could not be checked"):
        assert_trainable(Measured(S.OPEN_ATTRIBUTION, doi=object()))


def test_an_unbacked_claim_is_a_licence_gate_error_and_not_swallowable():
    assert issubclass(UnbackedClaimError, LicenceGateError)
    assert not issubclass(UnbackedClaimError, (ValueError, TypeError))


# --- the whole record is asked what it is, not just its outermost layer ---------


def test_declares_synthetic_reaches_a_declaration_on_a_component():
    from wmxglycan.licensing import declares_synthetic

    assert declares_synthetic(Measured(S.SYNTHETIC_FIXTURE))
    assert declares_synthetic(S.SYNTHETIC_FIXTURE) and declares_synthetic("synthetic_fixture")
    # The case the review found: a real-looking record built from an invented one.
    assert declares_synthetic(Measured(S.INTERNAL_PROPRIETARY, components=[a_structure(S.SYNTHETIC_FIXTURE)]))
    # And two levels down, because component_records is a tree, not a pair.
    deep = Measured(S.INTERNAL_PROPRIETARY, components=[StubRecord(S.INTERNAL_PROPRIETARY, components=[a_structure(S.SYNTHETIC_FIXTURE)])])
    assert declares_synthetic(deep)


def test_a_record_with_nothing_synthetic_in_it_is_not_declared_synthetic():
    from wmxglycan.licensing import declares_synthetic

    assert not declares_synthetic(Measured(S.INTERNAL_PROPRIETARY, components=[a_structure(S.INTERNAL_PROPRIETARY)]))
    assert not declares_synthetic(S.OPEN_ATTRIBUTION)


def test_a_record_whose_components_cannot_be_listed_counts_as_synthetic():
    from wmxglycan.licensing import declares_synthetic

    # Fail closed toward labelling: presenting something unreadable as real data
    # is the failure that matters. (The gate refuses such a record anyway.)
    assert declares_synthetic(StubRecord(S.INTERNAL_PROPRIETARY, components=AttributeError("glycan")))


def test_declares_synthetic_survives_a_cycle():
    from wmxglycan.licensing import declares_synthetic

    looped = StubRecord(S.INTERNAL_PROPRIETARY)
    looped._components = [looped]
    assert not declares_synthetic(looped)


# --- nothing leaves the gate except a gate error ---------------------------------


def test_a_record_whose_source_cannot_be_read_is_refused_not_crashed():
    # getattr's default absorbs only AttributeError. A raw exception here would
    # pass every `except TrainingGateError` a caller has and abort a whole load.
    class Exploding:
        reuse_status = S.SYNTHETIC_FIXTURE

        @property
        def source(self):
            raise RuntimeError("source blew up")

    with pytest.raises(LicenceGateError, match="source unreadable"):
        assert_trainable(Exploding())


def test_a_record_whose_status_cannot_be_read_is_refused_not_crashed():
    class Exploding:
        source = "synthetic test fixture"

        @property
        def reuse_status(self):
            raise RuntimeError("status blew up")

    with pytest.raises(LicenceGateError, match="could not be read"):
        assert_trainable(Exploding())


def test_a_component_is_checked_on_its_own_doi_when_it_names_one():
    # A matching source string must not let the enclosing record's citation back
    # a component that cites something else: that is the laundering the
    # structure rule exists to prevent, reached from the other side.
    paper = "the paper's table"
    inner = a_structure(S.OPEN_ATTRIBUTION, paper)
    inner.doi = "10.1000/never-registered"
    with pytest.raises(UnbackedClaimError, match="no licence record"):
        assert_trainable(Measured(S.OPEN_ATTRIBUTION, doi=STRUWE_2016.doi, source=paper, components=[inner]))
    # With its own DOI on record, the same component is backed.
    inner.doi = STRUWE_2016.doi
    assert_trainable(Measured(S.OPEN_ATTRIBUTION, doi=STRUWE_2016.doi, source=paper, components=[inner]))


def test_the_gate_imports_the_registry_without_a_cycle():
    # The check must be in the gate, so the gate imports the registry, which
    # imports the dataset loader, the models, and the statuses. Importing the
    # gate first, in a fresh interpreter, must not trip over that chain.
    code = "import wmxglycan.licensing as g; g.assert_trainable('synthetic_fixture'); print('ok')"
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0 and "ok" in result.stdout, result.stderr
