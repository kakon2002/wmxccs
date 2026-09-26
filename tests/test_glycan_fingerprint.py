"""The pipeline fingerprint: what moves it, what does not, and the honesty of the difference.

A fingerprint is only worth carrying if its readers trust it, and the way to lose that is to
have it move for reasons that cannot change an answer. So this file pins BOTH directions: which
changes must move it, and which must not.
"""

from __future__ import annotations

import dataclasses

import pytest

from wmxglycan import __version__
from wmxglycan.attestation import build_attestation_index, default_attestation_index
from wmxglycan.enumeration import Enumerator
from wmxglycan.fingerprint import (
    DataSnapshot,
    default_fingerprint,
    fingerprint_of,
    parameters_of,
)
from wmxglycan.sugarbase import LoadReport, dataset_version


@pytest.fixture(scope="module")
def enumerator():
    return Enumerator()


@pytest.fixture(scope="module")
def index():
    return default_attestation_index()


@pytest.fixture(scope="module")
def fingerprint(enumerator, index):
    return fingerprint_of(enumerator, index)


# --- it is stable ---------------------------------------------------------------------------------


def test_the_fingerprint_is_the_same_on_two_builds(enumerator, index):
    assert fingerprint_of(enumerator, index).short == fingerprint_of(enumerator, index).short


def test_the_fingerprint_does_not_depend_on_dict_insertion_order(enumerator, index):
    # The digest is taken over a JSON rendering with sorted keys, so this holds. Without
    # sort_keys it would move between runs under a different hash seed.
    one = parameters_of(enumerator)
    other = {key: one[key] for key in reversed(list(one))}
    from wmxglycan.fingerprint import _digest

    assert _digest(one) == _digest(other)


def test_the_cached_default_matches_a_fresh_build(fingerprint):
    assert default_fingerprint().short == fingerprint.short


def test_the_two_halves_are_reported_in_full_as_well_as_short(fingerprint):
    assert len(fingerprint.parameters) == 64
    assert len(fingerprint.corpus) == 64
    assert fingerprint.short == f"{fingerprint.parameters[:12]}/{fingerprint.corpus[:12]}"
    assert fingerprint.version == __version__


# --- what MUST move it ----------------------------------------------------------------------------


def test_a_changed_prior_moves_the_parameters_digest(enumerator, index, monkeypatch):
    before = fingerprint_of(enumerator, index).parameters
    monkeypatch.setattr("wmxglycan.fingerprint.PRIOR_PSEUDOCOUNT", 0.5)
    assert fingerprint_of(enumerator, index).parameters != before


def test_a_changed_policy_prior_moves_the_parameters_digest(enumerator, index, monkeypatch):
    before = fingerprint_of(enumerator, index).parameters
    monkeypatch.setattr("wmxglycan.fingerprint.POLICY_PRIOR", "jeffreys_0.5")
    assert fingerprint_of(enumerator, index).parameters != before


def test_a_changed_feature_set_moves_the_parameters_digest(enumerator, index, monkeypatch):
    before = fingerprint_of(enumerator, index).parameters
    monkeypatch.setattr("wmxglycan.fingerprint.FEATURE_NAMES", ("only_one_column",))
    assert fingerprint_of(enumerator, index).parameters != before


def test_a_dropped_curated_rule_moves_the_parameters_digest(enumerator, index):
    from wmxglycan.constraints import ConstraintSet

    before = fingerprint_of(enumerator, index).parameters
    fewer = Enumerator(
        catalogue=enumerator.catalogue,
        constraints=ConstraintSet(tuple(enumerator.constraints)[:-1]),
    )
    assert fingerprint_of(fewer, index).parameters != before


def test_a_different_corpus_moves_the_corpus_digest(enumerator, index):
    from wmxglycan.models import GlycanStructure
    from wmxglycan.reuse import ReuseStatus

    before = fingerprint_of(enumerator, index).corpus
    one = GlycanStructure(
        composition="Hex5HexNAc2",
        iupac_condensed="Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
        has_unresolved_linkage=False,
        has_unresolved_anomericity=False,
        source="a test",
        reuse_status=ReuseStatus.SYNTHETIC_FIXTURE,
    )
    smaller = build_attestation_index(
        LoadReport(dataset_version=dataset_version(), structures=(one,))
    )
    assert fingerprint_of(enumerator, smaller).corpus != before


def test_two_corpora_with_the_same_counts_but_different_structures_differ(enumerator):
    """Counts alone would not distinguish them, which is why the keys are digested."""
    from wmxglycan.models import GlycanStructure
    from wmxglycan.reuse import ReuseStatus

    def one(iupac: str, composition: str):
        return GlycanStructure(
            composition=composition,
            iupac_condensed=iupac,
            has_unresolved_linkage=False,
            has_unresolved_anomericity=False,
            source="a test",
            reuse_status=ReuseStatus.SYNTHETIC_FIXTURE,
        )

    first = build_attestation_index(
        LoadReport(
            dataset_version=dataset_version(),
            structures=(one("Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc", "Hex3HexNAc2"),),
        )
    )
    second = build_attestation_index(
        LoadReport(
            dataset_version=dataset_version(),
            structures=(
                one("Man(a1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc", "Hex4HexNAc2"),
            ),
        )
    )
    assert first.structures_indexed == second.structures_indexed == 1
    assert first.distinct_structures == second.distinct_structures == 1
    assert DataSnapshot.of(first).digest != DataSnapshot.of(second).digest


# --- what must NOT move it, and the cost of that choice, stated ------------------------------------


def test_the_fingerprint_does_not_digest_module_source(enumerator, index):
    """DELIBERATE, AND THE COST IS STATED RATHER THAN HIDDEN.

    A fingerprint that moves on a comment trains its readers to ignore it, and the CCS core's has
    held across five releases precisely because it digests inputs and parameters rather than
    code. The consequence is real: A LOGIC CHANGE THAT TOUCHES NO TABLE, NO CORPUS AND NO POLICY
    CONSTANT WILL NOT MOVE THIS DIGEST. The mutation catalogue is what guards that, not this.
    """
    from wmxglycan import fingerprint as module

    text = parameters_of(enumerator)
    # Nothing in the digested payload is source text.
    rendered = repr(text)
    assert "def " not in rendered
    assert "import " not in rendered
    assert module.__doc__ is not None
    assert "A LOGIC CHANGE THAT ALTERS NO TABLE" in module.__doc__


def test_the_parameters_payload_names_what_it_covers(enumerator):
    covered = parameters_of(enumerator)
    assert set(covered) == {
        "constraints",
        "monolinks",
        "enzyme_count",
        "sites",
        "sites_dropped",
        "features",
        "prior_policy",
        "prior_pseudocount",
        "priors",
    }
    # Returned rather than digested directly so a test can say WHICH parameter moved.
    assert len(covered["constraints"]) == len(tuple(enumerator.constraints))
    assert len(covered["features"]) == 44


def test_the_snapshot_carries_the_licence_and_the_attribution(index):
    snapshot = DataSnapshot.of(index)
    assert snapshot.licence == "MIT"
    assert "Daniel Bojar" in snapshot.attribution
    assert "SugarBase v12" in snapshot.dataset
    assert snapshot.structures_indexed == 4001
    assert snapshot.distinct_structures == 3640
    assert snapshot.duplicate_rows == 361
    assert snapshot.compositions == 687


def test_the_summary_prints_both_digests_and_the_licence(fingerprint):
    text = fingerprint.summary()
    assert fingerprint.parameters in text
    assert fingerprint.corpus in text
    assert "MIT" in text and "Daniel Bojar" in text
    assert "4001 resolved" in text


def test_a_snapshot_is_frozen_so_it_cannot_drift_after_being_stamped(fingerprint):
    with pytest.raises(dataclasses.FrozenInstanceError):
        fingerprint.snapshot.structures_indexed = 0  # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        fingerprint.parameters = "x"  # type: ignore[misc]
