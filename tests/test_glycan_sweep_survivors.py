"""The six behaviours the 27 September mutation sweep found with no test behind them.

WHY THIS FILE EXISTS. The sweep of `c145eaf` ran 120 mutations and six survived - four in
`api.py` and two in `store.py`. A surviving mutation is not a note here: it is a behaviour the
suite would not have noticed breaking, and this repository treats that as a failure. All six were
in code written the same week, which is the honest version of the pattern: the tests written
alongside new code assert what the code does on the happy path, and the sweep is what finds the
guard nobody broke on purpose.

WHAT MADE THEM HARD TO REACH, AND WHY THAT IS NOT AN EXCUSE. Four of the six only misbehave in a
state the rest of the system prevents - a prediction that moved under an attach, a held value
carrying a reference, a page smaller than its result. That is precisely the shape of a
defence-in-depth check, and a defence-in-depth check with no test is decoration. Each test below
therefore manufactures the state the outer layer prevents, and says which outer layer it is
bypassing and why that is legitimate.

EACH TEST FIRST PROVES IT COULD FAIL. Four tests earlier this week could not fail - a one-element
set literal, a `len(frozenset) == len(set(frozenset))`, an unasserted stated property, and a sweep
that never reached the branch it asserted about. So where a test asserts that two readings of the
same thing differ, it asserts that they really do differ for this input, rather than only
asserting the one it prefers.
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from wmxglycan.api import create_app
from wmxglycan.attestation import default_attestation_index
from wmxglycan.ccs_evidence import (
    CCSEvidence,
    CCSEvidenceState,
    CCSReference,
    EvidenceLevel,
    HeldValues,
)
from wmxglycan.enumeration import Enumerator
from wmxglycan.fingerprint import default_fingerprint
from wmxglycan.prediction import ComparisonState
from wmxglycan.reuse import ReuseStatus
from wmxglycan.store import Store, StoreError

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

G2F = "Hex5HexNAc4Fuc1"
MAN5 = "Hex5HexNAc2"

MEASUREMENT = {
    "ccs": 712.4,
    "uncertainty": 2.1,
    "uncertainty_type": "SD",
    "adduct": "[M+H]+",
    "charge": 1,
    "ims_type": "TWIMS",
    "drift_gas": "N2",
    "source": "lab run 1",
}

FROZEN = {
    "composition": "Hex5HexNAc2",
    "adduct": "[M-H]-",
    "charge": -1,
    "pipeline_version": "0.1.0",
    "fingerprint": "aaaaaaaaaaaa/bbbbbbbbbbbb",
    "snapshot_digest": "c" * 64,
}


@pytest.fixture(scope="module")
def shared():
    """The enumerator, index and fingerprint, built once: each costs seconds."""
    return {
        "enumerator": Enumerator(),
        "index": default_attestation_index(),
        "fingerprint": default_fingerprint(),
    }


def created(client, composition=MAN5, adduct="[M-H]-", charge=-1, **extra):
    body = {"composition": composition, "adduct": adduct, "charge": charge, **extra}
    response = client.post("/v1/predictions", json=body)
    assert response.status_code == 201, response.text
    return response.json()


def a_reference(**over) -> CCSReference:
    return CCSReference(
        **{
            "ccs": 700.0,
            "uncertainty": 3.0,
            "uncertainty_type": "SD",
            "adduct": "[M+H]+",
            "charge": 1,
            "polarity": "positive",
            "ims_type": "TWIMS",
            "drift_gas": "N2",
            "source": "a synthetic fixture",
            "doi": "10.0000/fixture",
            "source_locator": "Table S1",
            "reuse_status": ReuseStatus.SYNTHETIC_FIXTURE,
            **over,
        }
    )


# --- [S] store.py ---------------------------------------------------------------------------------


def test_the_attach_guard_notices_a_prediction_that_moved_under_it(tmp_path):
    """`store.attach` re-reads the prediction after the insert and refuses if its digest moved.

    THE STATE BEING MANUFACTURED. Nothing in this module can move a frozen payload - there is no
    UPDATE and no DELETE statement in it, asserted elsewhere by reading the source. So the only
    way to exercise the guard is to make the re-read disagree, which is what this subclass does on
    the SECOND read inside `attach` - the one after the INSERT. The guard exists for a future edit
    that starts touching the prediction here, and that future edit is what this stands in for.
    """

    class MovesDuringAttach(Store):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._inside_attach = False
            self._reads = 0

        def attach(self, **kwargs):
            self._inside_attach, self._reads = True, 0
            try:
                return super().attach(**kwargs)
            finally:
                self._inside_attach = False

        def get(self, prediction_id):
            row = super().get(prediction_id)
            if self._inside_attach and row is not None:
                self._reads += 1
                if self._reads >= 2:
                    return dataclasses.replace(row, payload_digest="0" * 64)
            return row

    store = MovesDuringAttach(tmp_path / "moved.sqlite3")
    row = store.freeze(
        prediction_id="pred_moved", created_at="2026-09-27T00:00:00Z",
        payload={"body": "the frozen answer"}, **FROZEN,
    )

    with pytest.raises(StoreError, match="must never touch the prediction"):
        store.attach(
            prediction_id=row.prediction_id,
            validation_id="val_moved",
            created_at="2026-09-27T00:01:00Z",
            fingerprint=FROZEN["fingerprint"],
            snapshot_digest=FROZEN["snapshot_digest"],
            payload={"measurements": [MEASUREMENT]},
        )

    # AND IT ROLLED BACK. Refusing after keeping the row would be the worse half of the failure:
    # a validation recorded against a prediction the caller was told not to trust.
    assert store.validations(row.prediction_id) == ()
    store.close()


def test_the_run_total_counts_the_whole_result_and_not_one_page(tmp_path):
    """`count_runs` is what a paged response calls itself a page OF, so it must not be paged."""
    store = Store(tmp_path / "runs.sqlite3")
    for index in range(55):
        store.record_run(
            run_id=f"run_{index:03d}",
            created_at=f"2026-09-27T00:{index:02d}:00Z",
            kind="prediction.created",
            fingerprint=FROZEN["fingerprint"],
            snapshot_digest=FROZEN["snapshot_digest"],
            detail="one run",
            composition="Hex5HexNAc2",
        )

    # The page really is smaller than the result, which is what makes the two readings differ and
    # this test able to fail. Without this line the assertion below would hold at any page size.
    assert len(store.runs()) == 50
    assert store.count_runs() == 55
    assert store.count_runs(kind="prediction.created") == 55
    assert store.count_runs(composition="Hex5HexNAc2") == 55
    assert store.count_runs(kind="validation.attached") == 0
    store.close()


# --- [A] api.py -----------------------------------------------------------------------------------


def test_attestation_is_reported_per_candidate_and_never_read_off_its_class(shared):
    """A class of 10 holding 4 attested structures has 4 attested candidates, not 10.

    This is the defect the code comment at the site describes, and it was fixed without a test.
    For G2F the two readings differ by more than a factor of two, which is asserted here so the
    test cannot pass by the readings happening to agree.
    """
    app = create_app(store=Store(), **shared)
    with TestClient(app) as client:
        body = created(client, G2F, "[M+H]+", 1)

    index = shared["index"]
    per_candidate = [one for one in body["candidates"] if one["attested"]]
    from_the_index = [one for one in body["candidates"] if index.attests(one["canonical_key"])]
    assert len(per_candidate) == len(from_the_index) > 0
    assert {one["canonical_key"] for one in per_candidate} == {
        one["canonical_key"] for one in from_the_index
    }

    attested_classes = {
        one["class_id"] for one in body["classes"] if one["attested_structures"] > 0
    }
    read_off_the_class = [
        one for one in body["candidates"] if one["class_id"] in attested_classes
    ]
    assert len(read_off_the_class) > len(per_candidate), (
        "the two readings agree for this input, so this test would pass on the wrong one;"
        " pick a composition with a partly attested class"
    )

    # And the class keeps its own count, which is a statement about the CLASS and still true.
    for one in body["classes"]:
        members = [c for c in body["candidates"] if c["class_id"] == one["class_id"]]
        assert one["attested_structures"] <= len(members)


def test_the_attach_response_computes_whether_the_prediction_moved_rather_than_asserting_it(shared):
    """`prediction_unchanged` is the evidence offered to a caller who should not have to trust us.

    THE STATE BEING MANUFACTURED. The store refuses an attach that moves the prediction, so in a
    working system this field is always True - which is exactly why a hardcoded True would never
    be noticed. The store here reports a moved digest only on the read the ENDPOINT makes after
    `attach` has returned, so the store's own guard is not what is under test.
    """

    class MovedAfterAttach(Store):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._attached = False

        def attach(self, **kwargs):
            record = super().attach(**kwargs)
            self._attached = True
            return record

        def get(self, prediction_id):
            row = super().get(prediction_id)
            if self._attached and row is not None:
                return dataclasses.replace(row, payload_digest="f" * 64)
            return row

    app = create_app(store=MovedAfterAttach(), **shared)
    with TestClient(app) as client:
        body = created(client, MAN5, "[M-H]-", -1)
        response = client.post(
            f"/v1/predictions/{body['prediction_id']}/validation",
            json={"measurements": [MEASUREMENT]},
        )
        assert response.status_code == 201, response.text
        attached = response.json()

    assert attached["prediction_digest_before"] != attached["prediction_digest_after"], (
        "the two digests are equal, so this test cannot tell a computed field from a constant"
    )
    assert attached["prediction_unchanged"] is False


def test_the_decision_is_re_derived_from_the_evidence_frozen_beside_it(shared):
    """One published decision rule reads the CCS evidence, so attaching evidence must recompute it.

    `rank()` runs without an evidence source and the lookup happens once in the endpoint. If the
    endpoint sets the evidence field without re-deriving, the response serves rules that disagree
    with the evidence printed next to them - a prediction whose own reasons contradict its data.
    """
    # A REAL candidate key, not a placeholder string. `look_up` in api.py swallows every exception
    # from an evidence source into LOOKUP_FAILED - correctly, so a wiring bug cannot read as an
    # absence - which also means a fixture that fails validation would silently become the state
    # this test is trying to move away from. A first version of this test passed a candidate OBJECT
    # where a key was required and spent its evidence on LOOKUP_FAILED without saying so.
    structure_key = shared["enumerator"].enumerate(G2F).candidates[0].canonical_key

    class Structural:
        def evidence_for(self, composition, adduct, charge):
            return CCSEvidence(
                state=CCSEvidenceState.MEASURED_REFERENCE,
                level=EvidenceLevel.STRUCTURE,
                composition=composition.canonical,
                adduct=adduct,
                charge=charge,
                discriminates_between_candidates=True,
                reference=a_reference(
                    adduct=adduct, charge=charge, measured_on_canonical_key=structure_key
                ),
            )

    app = create_app(store=Store(), evidence=Structural(), **shared)
    with TestClient(app) as client:
        body = created(client, G2F, "[M+H]+", 1)

    rules = {one["name"]: one for one in body["decision"]["rules"]}
    held = rules["no cross section is held for this structure"]
    assert held["fires"] is False, (
        "structure-level evidence was frozen into this prediction and the rule that reads it"
        " still reports fires=True, so the rules were not re-derived from the evidence"
    )

    # The decision itself does not move, and that is the point of the ruling of 27 September: two
    # other rules fire on every input, so IM_VALIDATION_REQUIRED stays even with a real reference.
    assert body["decision"]["decision"] == "IM_VALIDATION_REQUIRED"
    assert body["decision"]["constant_today"] is True
    # Without this the assertion above could hold because NO rule reads the evidence at all.
    assert any(one["fires"] for one in rules.values())


def test_a_held_value_cannot_reach_the_comparison_even_if_the_evidence_carries_one(shared):
    """The licence gate, at the last layer that could leak it.

    THE OUTER LAYER BEING BYPASSED. `CCSEvidence`'s validator forbids a non-MEASURED_REFERENCE
    finding from carrying a reference at all, so this state is unconstructible through the normal
    door and `model_construct` is used to walk past it. That is deliberate: the comparison's own
    check on the state is defence in depth behind that validator, and defence in depth that
    nothing tests is decoration. If the validator is ever relaxed - and it is one line - this is
    the test that stops a held number being served as a delta.
    """
    held = CCSEvidence.model_construct(
        state=CCSEvidenceState.HELD_NOT_RELEASABLE,
        level=EvidenceLevel.ION,
        composition="Hex5HexNAc2",
        adduct="[M-H]-",
        charge=-1,
        reference=a_reference(ccs=1234.5, adduct="[M-H]-", charge=-1),
        held=HeldValues(
            records=8,
            blockers=("the drift gas is UNSTATED",),
            any_blocker_resolves_alone=False,
            doi="10.1039/c5an01092f",
            what_would_release_them="nothing available today",
        ),
        corpus_searched=None,
        records_consulted=None,
        failure=None,
        key_attempted=None,
        discriminates_between_candidates=False,
    )
    # The fixture really is the forbidden combination, or this test proves nothing.
    assert held.state is not CCSEvidenceState.MEASURED_REFERENCE
    assert held.reference is not None

    class Leaky:
        def evidence_for(self, composition, adduct, charge):
            return held

    app = create_app(store=Store(), evidence=Leaky(), **shared)
    with TestClient(app) as client:
        body = created(client, MAN5, "[M-H]-", -1)
        prediction_id = body["prediction_id"]
        client.post(
            f"/v1/predictions/{prediction_id}/validation",
            json={"measurements": [{**MEASUREMENT, "adduct": "[M-H]-", "charge": -1}]},
        )
        response = client.get(f"/v1/predictions/{prediction_id}/comparison")
        assert response.status_code == 200, response.text
        comparison = response.json()

    # The reference axis is reported PER MEASUREMENT; the top-level status is the prediction axis,
    # which is undefined here for the separate reason that no cross section was predicted at all.
    assert comparison["against_prediction"] == ComparisonState.NO_PREDICTED_VALUE.value
    assert comparison["per_measurement"], "nothing was compared, so nothing below is asserted"
    for measurement in comparison["per_measurement"]:
        assert measurement["against_reference"] == (
            ComparisonState.REFERENCE_HELD_NOT_RELEASABLE.value
        )
        assert measurement["reference_ccs"] is None
        assert measurement["delta_ccs"] is None
        assert measurement["delta_percent"] is None

    # And the number itself appears NOWHERE in the response, not merely not in the delta fields.
    assert "1234.5" not in response.text
