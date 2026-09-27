"""The six endpoints, over the real ranker and the real corpus.

No endpoint is mocked. The app a test builds is the app a deployment builds, with the store and
the evidence source injected - which is the reason `create_app` takes them rather than reaching
for module-level globals.

THE THREE RULINGS OF 27 SEPTEMBER 2026 ARE ENFORCED HERE, not documented:

  - nothing is overwritten, and a second write to a frozen prediction is refused rather than
    merged. Tested at the API as well as at the store, because a store that refuses and an
    endpoint that quietly retries would together produce the thing the rule forbids.
  - no served field is named `probability`. Walked over every response model, so a field added
    later cannot slip past by not being in a hand-written list.
  - the decision field is constant, and the response says so. Asserted as a property of the
    responses rather than of the documentation.
"""

from __future__ import annotations

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
from wmxglycan.prediction import AGREEMENT_LIMIT_PERCENT, ComparisonState
from wmxglycan.reuse import ReuseStatus
from wmxglycan.store import Store

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


@pytest.fixture(scope="module")
def shared():
    """The enumerator, index and fingerprint, built once: each costs seconds."""
    return {
        "enumerator": Enumerator(),
        "index": default_attestation_index(),
        "fingerprint": default_fingerprint(),
    }


@pytest.fixture
def client(shared):
    app = create_app(store=Store(), **shared)
    with TestClient(app) as made:
        yield made


@pytest.fixture
def held_client(shared):
    """A client whose evidence source reports HELD values, as the seed corpus does for Man5."""

    class Held:
        def evidence_for(self, composition, adduct, charge):
            return CCSEvidence(
                state=CCSEvidenceState.HELD_NOT_RELEASABLE,
                level=EvidenceLevel.ION,
                composition=composition.canonical,
                adduct=adduct,
                charge=charge,
                held=HeldValues(
                    records=8,
                    blockers=("the drift gas is UNSTATED", "uncertainty_type is unknown"),
                    any_blocker_resolves_alone=False,
                    doi="10.1039/c5an01092f",
                    what_would_release_them="nothing available today",
                ),
            )

    app = create_app(store=Store(), evidence=Held(), **shared)
    with TestClient(app) as made:
        yield made


def created(client, composition=MAN5, adduct="[M-H]-", charge=-1, **extra):
    body = {"composition": composition, "adduct": adduct, "charge": charge, **extra}
    response = client.post("/v1/predictions", json=body)
    assert response.status_code == 201, response.text
    return response.json()


# --- 1. POST /v1/predictions ---------------------------------------------------------------------


def test_creating_a_prediction_returns_201_and_a_frozen_answer(client):
    body = created(client, G2F, "[M+H]+", 1)
    assert body["frozen"] is True
    assert body["composition"] == G2F
    assert body["candidates_total"] == 167
    assert body["classes_total"] == 61
    assert len(body["bands"]) == 3
    assert body["tied_candidates"] == 152
    assert body["largest_indistinguishable_class"] == 10
    assert body["is_a_ranking"] is True


def test_a_composition_that_does_not_parse_is_refused_with_422(client):
    response = client.post(
        "/v1/predictions", json={"composition": "Hex5HexNAcFuc1", "adduct": "[M+H]+", "charge": 1}
    )
    assert response.status_code == 422
    assert "could not be read" in response.json()["detail"]
    assert "Nothing is partly parsed here" in response.json()["detail"]


def test_charge_zero_is_refused(client):
    response = client.post(
        "/v1/predictions", json={"composition": MAN5, "adduct": "[M]", "charge": 0}
    )
    assert response.status_code == 422


def test_a_composition_the_enumerator_declines_is_still_a_prediction_with_a_reason(client):
    # Not an error: the platform answering "no candidate, and here is why" is an answer.
    body = created(client, "Hex2HexNAc2Fuc2", "[M+H]+", 1)
    assert body["candidates_total"] == 0
    assert body["is_a_ranking"] is False
    assert "no candidate" in body["refusal"]
    # And the coverage still reports the corpus it actually has, rather than asserting an absence.
    assert body["coverage"]["reference_structures"] == 15


def test_man5_is_served_as_a_refusal_to_order_rather_than_an_order(client):
    body = created(client)
    assert body["candidates_total"] == 6
    assert body["is_a_ranking"] is False
    assert "RANKING REFUSED" in body["refusal"]
    # The candidates are STILL RETURNED. What is withheld is an order, not an answer.
    assert len(body["candidates"]) == 6
    assert len(body["classes"]) == 6


# --- the immutability ruling, at the API ---------------------------------------------------------


def test_a_second_create_under_the_same_client_reference_is_refused_with_409(client):
    first = created(client, G2F, "[M+H]+", 1, client_reference="run-A")
    again = client.post(
        "/v1/predictions",
        json={
            "composition": G2F,
            "adduct": "[M+H]+",
            "charge": 1,
            "client_reference": "run-A",
        },
    )
    assert again.status_code == 409
    detail = again.json()["detail"]
    assert detail["error"] == "already_frozen"
    # THE ID THAT STANDS is named, so a caller fetches it rather than guessing.
    assert detail["prediction_id"] == first["prediction_id"]
    assert "NOT modified" in detail["detail"]
    # And exactly one prediction exists, not two.
    runs = client.get("/v1/runs", params={"kind": "prediction.created"}).json()
    assert runs["total"] == 1


def test_a_refused_second_create_does_not_change_the_first(client):
    first = created(client, MAN5, "[M-H]-", -1, client_reference="run-B")
    client.post(
        "/v1/predictions",
        json={
            "composition": G2F,
            "adduct": "[M+H]+",
            "charge": 1,
            "client_reference": "run-B",
        },
    )
    again = client.get(f"/v1/predictions/{first['prediction_id']}").json()
    assert again == first, "the standing prediction must be byte-identical after a refused write"


def test_two_creates_without_a_client_reference_are_two_predictions(client):
    one = created(client, MAN5)
    two = created(client, MAN5)
    assert one["prediction_id"] != two["prediction_id"]
    # Creating twice is not a mutation. What is refused is a second write to ONE prediction.
    assert one["candidates_total"] == two["candidates_total"]


# --- 2. GET /v1/predictions/{id} -----------------------------------------------------------------


def test_reading_a_prediction_back_returns_exactly_what_was_served(client):
    body = created(client, G2F, "[M+H]+", 1)
    again = client.get(f"/v1/predictions/{body['prediction_id']}")
    assert again.status_code == 200
    assert again.json() == body


def test_reading_a_prediction_does_not_recompute_it(client, shared):
    """THE POINT OF FREEZING. A read must not re-run the enumerator, because then a change to
    the tables or the corpus would silently alter what a past prediction says."""
    body = created(client, MAN5)
    # Swap in an enumerator that would raise if it were consulted. A read must not touch it.
    class Exploding:
        def enumerate(self, *args, **kwargs):
            raise AssertionError("a read recomputed the prediction instead of serving the bytes")

        constraints = shared["enumerator"].constraints

    client.app.state.enumerator = Exploding()
    again = client.get(f"/v1/predictions/{body['prediction_id']}")
    assert again.status_code == 200
    assert again.json() == body


def test_an_unknown_prediction_is_404_and_says_ids_are_never_deleted(client):
    response = client.get("/v1/predictions/pred_nothing")
    assert response.status_code == 404
    assert "was never issued" in response.json()["detail"]


# --- 3. POST /v1/predictions/{id}/validation -----------------------------------------------------


def test_attaching_measurements_returns_201_and_proves_the_prediction_was_not_touched(client):
    body = created(client, G2F, "[M+H]+", 1)
    pid = body["prediction_id"]
    response = client.post(
        f"/v1/predictions/{pid}/validation", json={"measurements": [MEASUREMENT]}
    )
    assert response.status_code == 201
    attached = response.json()
    assert attached["measurements_attached"] == 1
    assert attached["validations_on_this_prediction"] == 1
    # READ FROM THE STORE BEFORE AND AFTER, not asserted: the evidence is in the reply.
    assert attached["prediction_digest_before"] == attached["prediction_digest_after"]
    assert attached["prediction_unchanged"] is True
    assert client.get(f"/v1/predictions/{pid}").json() == body


def test_a_second_attach_appends_and_never_merges(client):
    body = created(client, G2F, "[M+H]+", 1)
    pid = body["prediction_id"]
    for n in range(3):
        response = client.post(
            f"/v1/predictions/{pid}/validation",
            json={"measurements": [{**MEASUREMENT, "ccs": 710.0 + n, "source": f"run {n}"}]},
        )
        assert response.status_code == 201
        assert response.json()["validations_on_this_prediction"] == n + 1
    # Three validations, three measurements, and the prediction still exactly as frozen.
    comparison = client.get(f"/v1/predictions/{pid}/comparison").json()
    assert comparison["validations"] == 3
    assert comparison["measurements"] == 3
    assert client.get(f"/v1/predictions/{pid}").json() == body


def test_a_measurement_without_an_uncertainty_type_is_refused(client):
    body = created(client, G2F, "[M+H]+", 1)
    for bad in ("unknown", "", "n/a"):
        response = client.post(
            f"/v1/predictions/{body['prediction_id']}/validation",
            json={"measurements": [{**MEASUREMENT, "uncertainty_type": bad}]},
        )
        assert response.status_code == 422, bad


def test_a_measurement_against_an_unstated_gas_is_refused(client):
    body = created(client, G2F, "[M+H]+", 1)
    response = client.post(
        f"/v1/predictions/{body['prediction_id']}/validation",
        json={"measurements": [{**MEASUREMENT, "drift_gas": "UNSTATED"}]},
    )
    assert response.status_code == 422


def test_attaching_to_an_unknown_prediction_is_404(client):
    response = client.post(
        "/v1/predictions/pred_nothing/validation", json={"measurements": [MEASUREMENT]}
    )
    assert response.status_code == 404


def test_an_empty_measurement_list_is_refused(client):
    body = created(client, G2F, "[M+H]+", 1)
    response = client.post(
        f"/v1/predictions/{body['prediction_id']}/validation", json={"measurements": []}
    )
    assert response.status_code == 422


# --- 4. GET /v1/predictions/{id}/comparison -------------------------------------------------------


def test_the_comparison_says_there_is_no_predicted_value_and_why(client):
    body = created(client, G2F, "[M+H]+", 1)
    pid = body["prediction_id"]
    client.post(f"/v1/predictions/{pid}/validation", json={"measurements": [MEASUREMENT]})
    comparison = client.get(f"/v1/predictions/{pid}/comparison").json()
    assert comparison["against_prediction"] == ComparisonState.NO_PREDICTED_VALUE.value
    assert "it does not exist" in comparison["against_prediction_because"] or (
        "undefined" in comparison["against_prediction_because"]
    )
    assert comparison["interval_coverage"] is None
    assert "no interval was predicted" in comparison["interval_coverage_because"]
    assert comparison["agreement_limit_percent"] == AGREEMENT_LIMIT_PERCENT
    assert comparison["agreement_limit_is_policy"] is True


def test_with_no_evidence_source_a_reference_delta_is_not_consulted_rather_than_absent(client):
    body = created(client, G2F, "[M+H]+", 1)
    pid = body["prediction_id"]
    client.post(f"/v1/predictions/{pid}/validation", json={"measurements": [MEASUREMENT]})
    one = client.get(f"/v1/predictions/{pid}/comparison").json()["per_measurement"][0]
    assert one["against_reference"] == ComparisonState.REFERENCE_NOT_CONSULTED.value
    assert one["delta_ccs"] is None
    assert "nothing looked" in one["unevaluable_because"]


def test_a_held_reference_is_reported_as_held_and_never_as_a_delta(held_client):
    body = created(held_client, MAN5)
    pid = body["prediction_id"]
    held_client.post(f"/v1/predictions/{pid}/validation", json={"measurements": [MEASUREMENT]})
    one = held_client.get(f"/v1/predictions/{pid}/comparison").json()["per_measurement"][0]
    assert one["against_reference"] == ComparisonState.REFERENCE_HELD_NOT_RELEASABLE.value
    assert one["delta_ccs"] is None
    assert one["reference_ccs"] is None, "a held value must never reach the response"
    assert "8 value(s) exist" in one["unevaluable_because"]


def test_a_releasable_reference_produces_a_real_delta(shared):
    """THE POSITIVE CONTROL. Every state above is an unevaluable one against the real corpus, so
    without this the delta arithmetic would be unreachable and the endpoint would be a refusal
    generator whose only working path nothing exercises."""

    class Measured:
        def evidence_for(self, composition, adduct, charge):
            return CCSEvidence(
                state=CCSEvidenceState.MEASURED_REFERENCE,
                level=EvidenceLevel.COMPOSITION,
                composition=composition.canonical,
                adduct=adduct,
                charge=charge,
                reference=CCSReference(
                    ccs=700.0,
                    uncertainty=3.0,
                    uncertainty_type="SD",
                    adduct=adduct,
                    charge=charge,
                    polarity="positive",
                    ims_type="TWIMS",
                    drift_gas="N2",
                    source="a synthetic fixture",
                    doi="10.0000/fixture",
                    source_locator="Table S1",
                    reuse_status=ReuseStatus.SYNTHETIC_FIXTURE,
                ),
            )

    app = create_app(store=Store(), evidence=Measured(), **shared)
    with TestClient(app) as client:
        body = created(client, G2F, "[M+H]+", 1)
        pid = body["prediction_id"]
        client.post(
            f"/v1/predictions/{pid}/validation",
            json={"measurements": [{**MEASUREMENT, "ccs": 707.0}, {**MEASUREMENT, "ccs": 750.0}]},
        )
        comparison = client.get(f"/v1/predictions/{pid}/comparison").json()
        assert comparison["status"] == ComparisonState.COMPARED.value
        near, far = comparison["per_measurement"]
        assert near["delta_ccs"] == pytest.approx(7.0)
        assert near["delta_percent"] == pytest.approx(1.0)
        assert near["within_agreement_limit"] is True
        assert far["delta_percent"] == pytest.approx(50.0 / 7.0, abs=1e-6)
        assert far["within_agreement_limit"] is False
        # Even with a real delta, the prediction axis is still undefined.
        assert comparison["against_prediction"] == ComparisonState.NO_PREDICTED_VALUE.value


def test_the_spread_among_the_callers_own_measurements_is_always_computable(client):
    # The one comparison that needs no model and no reference, so it is the one that works.
    body = created(client, G2F, "[M+H]+", 1)
    pid = body["prediction_id"]
    client.post(
        f"/v1/predictions/{pid}/validation",
        json={
            "measurements": [
                {**MEASUREMENT, "ccs": 710.0},
                {**MEASUREMENT, "ccs": 713.0},
                {**MEASUREMENT, "ccs": 716.0},
            ]
        },
    )
    spreads = client.get(f"/v1/predictions/{pid}/comparison").json()["among_attached"]
    assert len(spreads) == 1
    only = spreads[0]
    assert only["measurements"] == 3
    assert only["lowest"] == 710.0 and only["highest"] == 716.0
    assert only["spread"] == pytest.approx(6.0)
    assert only["within_agreement_limit"] is True


def test_measurements_on_different_gases_are_not_spread_against_each_other(client):
    # Hard constraint 4: values measured against different reference gases are not
    # interchangeable, so they are grouped apart rather than averaged into one spread.
    body = created(client, G2F, "[M+H]+", 1)
    pid = body["prediction_id"]
    client.post(
        f"/v1/predictions/{pid}/validation",
        json={
            "measurements": [
                {**MEASUREMENT, "ccs": 710.0, "drift_gas": "N2"},
                {**MEASUREMENT, "ccs": 712.0, "drift_gas": "N2"},
                {**MEASUREMENT, "ccs": 300.0, "drift_gas": "He"},
                {**MEASUREMENT, "ccs": 302.0, "drift_gas": "He"},
            ]
        },
    )
    spreads = client.get(f"/v1/predictions/{pid}/comparison").json()["among_attached"]
    assert {one["drift_gas"] for one in spreads} == {"N2", "He"}
    assert all(one["measurements"] == 2 for one in spreads)
    assert all(one["spread"] == pytest.approx(2.0) for one in spreads)


def test_a_single_measurement_yields_no_spread_rather_than_a_zero(client):
    # A spread of one value is not zero agreement, it is no comparison. Reporting 0.0 would read
    # as perfect reproducibility from a single run.
    body = created(client, G2F, "[M+H]+", 1)
    pid = body["prediction_id"]
    client.post(f"/v1/predictions/{pid}/validation", json={"measurements": [MEASUREMENT]})
    assert client.get(f"/v1/predictions/{pid}/comparison").json()["among_attached"] == []


def test_comparing_with_nothing_attached_is_404_and_says_what_to_do(client):
    body = created(client, G2F, "[M+H]+", 1)
    response = client.get(f"/v1/predictions/{body['prediction_id']}/comparison")
    assert response.status_code == 404
    assert "POST to /v1/predictions" in response.json()["detail"]


# --- 5. GET /v1/runs -----------------------------------------------------------------------------


def test_the_run_log_records_creation_and_attachment(client):
    body = created(client, G2F, "[M+H]+", 1)
    client.post(
        f"/v1/predictions/{body['prediction_id']}/validation",
        json={"measurements": [MEASUREMENT]},
    )
    runs = client.get("/v1/runs").json()
    assert runs["total"] == 2
    assert [one["kind"] for one in runs["runs"]] == ["validation.attached", "prediction.created"]
    assert all(one["fingerprint"] for one in runs["runs"])
    assert all(one["data_snapshot"] for one in runs["runs"])


def test_the_run_filters_compose(client):
    a = created(client, G2F, "[M+H]+", 1)
    created(client, MAN5)
    client.post(f"/v1/predictions/{a['prediction_id']}/validation", json={"measurements": [MEASUREMENT]})
    assert client.get("/v1/runs", params={"kind": "prediction.created"}).json()["total"] == 2
    assert client.get("/v1/runs", params={"composition": MAN5}).json()["total"] == 1
    assert (
        client.get("/v1/runs", params={"prediction_id": a["prediction_id"]}).json()["total"] == 2
    )
    assert (
        client.get(
            "/v1/runs", params={"kind": "validation.attached", "composition": G2F}
        ).json()["total"]
        == 1
    )


def test_a_composition_filter_is_canonicalised(client):
    created(client, MAN5)
    # A caller searching "HexNAc2Hex5" must find runs stored as "Hex5HexNAc2".
    assert client.get("/v1/runs", params={"composition": "HexNAc2Hex5"}).json()["total"] == 1


def test_an_unparseable_composition_filter_returns_nothing_rather_than_an_error(client):
    created(client, MAN5)
    response = client.get("/v1/runs", params={"composition": "not a composition"})
    assert response.status_code == 200
    assert response.json()["total"] == 0


def test_paging_reports_the_total_before_the_page(client):
    for _ in range(5):
        created(client, MAN5)
    page = client.get("/v1/runs", params={"limit": 2, "offset": 0}).json()
    assert page["total"] == 5
    assert page["returned"] == 2
    assert page["limit"] == 2 and page["offset"] == 0
    last = client.get("/v1/runs", params={"limit": 2, "offset": 4}).json()
    assert last["returned"] == 1 and last["total"] == 5


def test_the_filters_that_were_applied_are_echoed(client):
    created(client, MAN5)
    body = client.get("/v1/runs", params={"kind": "prediction.created"}).json()
    assert body["filters"] == {"kind": "prediction.created"}


# --- 6. GET /v1/models/current -------------------------------------------------------------------


def test_the_current_model_reports_the_pipeline_and_says_no_ccs_model_is_fitted(client, shared):
    body = client.get("/v1/models/current").json()
    assert body["ccs_model_fitted"] is False
    assert "no fitted glycan CCS model" in body["ccs_model_note"]
    assert body["stamp"]["fingerprint"] == shared["fingerprint"].short
    assert body["stamp"]["data_snapshot"] == shared["fingerprint"].snapshot.digest
    assert body["calibration"] == "never_calibrated"
    assert "independently known" in body["what_would_calibrate_it"]


def test_the_domain_says_what_the_platform_will_not_answer_for(client, shared):
    """THE SHAPE of the domain block. Its CLAIMS are verified in test_glycan_served_claims.py.

    TWO ASSERTIONS WERE REMOVED FROM HERE ON 27 SEPTEMBER 2026, and they are the reason that file
    exists. They read:

        assert set(domain["residues_supported"]) == {"Hex", "HexNAc", "Fuc", "NeuAc", "NeuGc"}
        assert "misses between a third and four fifths" in domain["known_coverage_limit"]

    The first pinned a field that advertised NeuGc, which the enumerator refuses every composition
    for. The second pinned a range that 76 of 145 measured compositions fall outside. Both passed
    because they compared the served string to a copy of itself: a test that matches text can only
    notice that someone changed it, never that it is false. They are replaced by measurements -
    `test_every_residue_the_domain_advertises_can_actually_be_placed` submits a composition needing
    each residue, and the coverage figures are recomputed from the index.
    """
    domain = client.get("/v1/models/current").json()["domain"]
    assert domain["glycan_class"] == "N-linked"
    assert domain["predicts_ccs"] is False
    assert "complete branched Man3GlcNAc2 core" in domain["builds_on"]
    assert "mammalian" in domain["species_assumption"]
    # DERIVED, not listed: what the enumerator can place is what the domain must advertise.
    assert set(domain["residues_supported"]) == {
        residue.value for residue in shared["enumerator"].placeable_residues
    }
    assert domain["biosynthetic_rules"] == len(shared["enumerator"].constraints)
    # The block must still carry a coverage statement; what it says is checked by measurement.
    assert domain["known_coverage_limit"]


def test_the_reachable_decisions_are_reported_as_one(client):
    body = client.get("/v1/models/current").json()
    assert set(body["decision_values"]) == {
        "AI_ONLY",
        "IM_VALIDATION_RECOMMENDED",
        "IM_VALIDATION_REQUIRED",
    }
    assert body["decision_reachable_today"] == ["IM_VALIDATION_REQUIRED"]


def test_the_count_of_frozen_predictions_is_live(client):
    assert client.get("/v1/models/current").json()["predictions_frozen"] == 0
    created(client, MAN5)
    created(client, G2F, "[M+H]+", 1)
    assert client.get("/v1/models/current").json()["predictions_frozen"] == 2


# --- provenance on every response ----------------------------------------------------------------


def test_every_endpoint_carries_the_fingerprint_and_the_data_snapshot(client, shared):
    body = created(client, G2F, "[M+H]+", 1)
    pid = body["prediction_id"]
    client.post(f"/v1/predictions/{pid}/validation", json={"measurements": [MEASUREMENT]})
    paths = [
        f"/v1/predictions/{pid}",
        f"/v1/predictions/{pid}/comparison",
        "/v1/runs",
        "/v1/models/current",
    ]
    stamps = [body["stamp"]] + [client.get(path).json()["stamp"] for path in paths]
    # And the attach response too, which is a POST.
    stamps.append(
        client.post(
            f"/v1/predictions/{pid}/validation", json={"measurements": [MEASUREMENT]}
        ).json()["stamp"]
    )
    for stamp in stamps:
        assert stamp["fingerprint"] == shared["fingerprint"].short
        assert stamp["fingerprint_parameters"] == shared["fingerprint"].parameters
        assert stamp["fingerprint_corpus"] == shared["fingerprint"].corpus
        assert stamp["data_snapshot"] == shared["fingerprint"].snapshot.digest
        assert stamp["dataset_licence"] == "MIT"
        assert "Daniel Bojar" in stamp["dataset_attribution"]
        assert stamp["version"]


def test_a_frozen_prediction_keeps_the_stamp_it_was_made_under(client, shared):
    """So a past prediction can be reconstructed. If a read re-stamped it with whatever is
    installed now, the one thing the stamp is for would be lost."""
    body = created(client, MAN5)
    moved = "deadbeefdead/beefdeadbeef"
    client.app.state.fingerprint = type(shared["fingerprint"])(
        version="9.9.9",
        parameters="deadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeef",
        corpus="beefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdead",
        snapshot=shared["fingerprint"].snapshot,
    )
    again = client.get(f"/v1/predictions/{body['prediction_id']}").json()
    assert again["stamp"]["fingerprint"] == shared["fingerprint"].short
    assert again["stamp"]["fingerprint"] != moved
    assert again == body


# --- the naming ruling, enforced -----------------------------------------------------------------


def test_no_served_field_is_named_probability():
    """RULED 27 September 2026. Walked over every response model rather than a hand-written list,
    so a field added later cannot slip past by not being on the list - which is instance Twelve's
    shape and would be especially poor here, in the test that exists to prevent it."""
    import wmxglycan.prediction as contracts
    from pydantic import BaseModel

    models = [
        value
        for value in vars(contracts).values()
        if isinstance(value, type) and issubclass(value, BaseModel) and value is not BaseModel
    ]
    assert len(models) >= 15, f"the walk found only {len(models)} models"
    forbidden = ("probability", "p_correct", "likelihood", "prob_")
    for model in models:
        for name in model.model_fields:
            for word in forbidden:
                assert word not in name.lower(), f"{model.__name__}.{name}"


def test_the_number_is_served_as_confidence_with_its_prior_beside_it(client):
    body = created(client, G2F, "[M+H]+", 1)
    confidence = body["confidence"]
    assert confidence["prior"] == "uniform_1"
    assert confidence["prior_pseudocount"] == 1.0
    assert set(confidence["priors"]) == {
        "haldane_0",
        "jeffreys_0.5",
        "perks_1_over_k",
        "uniform_1",
    }
    assert confidence["calibration"] == "never_calibrated"
    assert "NOT A MOLECULE IN YOUR SAMPLE" in confidence["what_the_number_means"]
    # And every class carries its share under all four, so the lever is visible per row.
    for one in body["classes"]:
        assert set(one["evidence_share_under_priors"]) == set(confidence["priors"])


def test_the_reserved_masses_travel_so_the_number_cannot_read_as_closed(client):
    confidence = created(client, G2F, "[M+H]+", 1)["confidence"]
    assert confidence["mass_on_reference_structures_not_enumerated"] > 0
    assert confidence["mass_on_a_structure_nobody_proposed"] > 0
    assert set(confidence["mass_under_priors"]) == set(confidence["priors"])


# --- the constant decision, enforced -------------------------------------------------------------


@pytest.mark.parametrize(
    "composition,adduct,charge",
    [(G2F, "[M+H]+", 1), (MAN5, "[M-H]-", -1), ("Hex3HexNAc4Fuc1", "[M+Na]+", 1)],
)
def test_every_prediction_lands_on_validation_required(client, composition, adduct, charge):
    body = created(client, composition, adduct, charge)
    assert body["decision"]["decision"] == "IM_VALIDATION_REQUIRED"


def test_the_response_says_the_decision_is_constant_and_why(client):
    decision = created(client, G2F, "[M+H]+", 1)["decision"]
    assert decision["constant_today"] is True
    assert "only reachable value" in decision["why_it_is_constant"]
    assert "not a defect to be fixed" in decision["why_it_is_constant"]
    # And the rules are published with what each applies to, so the constancy is readable from
    # the response rather than only from the documentation.
    assert len(decision["rules"]) == 5
    universal = [
        rule["name"]
        for rule in decision["rules"]
        if rule["fires"] and ("validated model" in rule["name"] or "may not contain" in rule["name"])
    ]
    assert len(universal) == 2
    for rule in decision["rules"]:
        assert rule["applies_to"].strip()
        assert rule["because"].strip()


def test_a_held_reference_does_not_soften_the_decision(held_client):
    body = created(held_client, MAN5)
    assert body["ccs_evidence"]["state"] == "held_not_releasable"
    assert body["decision"]["decision"] == "IM_VALIDATION_REQUIRED"


# --- the wall ------------------------------------------------------------------------------------


def test_the_api_module_imports_nothing_from_the_ccs_core():
    import ast

    source = Path(__import__("wmxglycan.api", fromlist=["x"]).__file__).read_text(encoding="utf-8")
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            assert not node.module.startswith("wmxccs"), node.module
        if isinstance(node, ast.Import):
            for alias in node.names:
                assert not alias.name.startswith("wmxccs"), alias.name


def test_the_default_app_reports_not_consulted_rather_than_an_absence(client):
    # A deployment nobody wired must not make a claim about the literature.
    body = created(client, MAN5)
    assert body["ccs_evidence"]["state"] == "not_consulted"
    assert "nothing looked" in body["ccs_evidence"]["summary"]


def test_an_evidence_source_that_raises_is_reported_as_a_failure_not_an_absence(shared):
    class Broken:
        def evidence_for(self, composition, adduct, charge):
            raise RuntimeError("the adapter is misconfigured")

    app = create_app(store=Store(), evidence=Broken(), **shared)
    with TestClient(app) as client:
        body = created(client, MAN5)
        assert body["ccs_evidence"]["state"] == "lookup_failed"
        assert "RuntimeError" in body["ccs_evidence"]["failure"]


def test_the_composition_root_wires_both_packages_and_finds_the_held_values():
    """The only place the two packages meet, exercised against the real seed corpus."""
    from tools.glycan_service import build

    app = build(database=":memory:")
    with TestClient(app) as client:
        body = created(client, MAN5, "[M-H]-", -1)
        evidence = body["ccs_evidence"]
        assert evidence["state"] == "held_not_releasable"
        assert evidence["held"]["records"] == 8
        assert len(evidence["held"]["blockers"]) == 2
        assert evidence["held"]["any_blocker_resolves_alone"] is False
        assert evidence["reference"] is None
