"""The wired /harmonize: what it answers, what it refuses, and what it never extrapolates.

Separate from test_api.py, which tests the contract shapes and the provenance echo. This
file tests the endpoint against a REAL model fitted from the seed corpus, because the
behaviour being checked - which measurements are covered and which are not - is a property
of the data rather than of the routing.
"""

from __future__ import annotations

import json

import pytest
from fastapi.testclient import TestClient

from wmxccs.api import build_default_model, create_app
from wmxccs.contracts import CorrectionBasis, IntervalKind
from wmxccs.grading import ConfidenceGrade
from wmxccs.loader import load_measurements_file
from wmxccs.readiness import DataMaturity
from wmxccs.scope import ComparisonScope

from pathlib import Path

SEED = Path(__file__).resolve().parents[1] / "data" / "seed" / "steroid_jasms2022.csv"
# Named from THIS FILE rather than left to the package default, which resolves relative to
# the installed package. That default is right for a deployment and wrong for a test: the
# mutation harness runs this suite against a copy of the package in a temporary directory,
# where the default finds no data at all.
SEED_DIRECTORY = SEED.parent


@pytest.fixture(scope="module")
def records():
    return load_measurements_file(SEED).cleared


@pytest.fixture(scope="module")
def served():
    model = build_default_model(SEED_DIRECTORY)
    assert model is not None, "the seed corpus must yield a model or this whole file is vacuous"
    return TestClient(create_app(model))


def body_for(*records) -> dict:
    return {"measurements": [json.loads(record.model_dump_json()) for record in records]}


def a_covered_record(records):
    """A travelling-wave record: a platform the model compares against the primary."""
    return next(r for r in records if r.ims_type.value == "TWIMS")


def a_primary_record(records):
    return next(r for r in records if str(r.dtims_method or "") == "stepped_field")


# --- the model is actually loaded ----------------------------------------------------------


def test_health_reports_the_model_the_deployment_serves(served):
    payload = served.get("/health").json()
    assert payload["model_loaded"] is True
    assert payload["matched_ions_available"] == 142


def test_the_default_model_is_fitted_from_the_seed_corpus_on_disk():
    """The deployment wiring. `uvicorn wmxccs.api:app` must serve a model, not a skeleton."""
    model = build_default_model(SEED_DIRECTORY)
    assert len(model.applied) == 9
    assert model.maturity.data_maturity is DataMaturity.PROVISIONAL


def test_a_directory_with_no_seed_data_yields_no_model_rather_than_an_empty_one(tmp_path):
    """So an installation with no data answers 501 instead of 200-with-nothing."""
    assert build_default_model(tmp_path) is None


def test_data_that_loads_but_pairs_nothing_also_yields_no_model(tmp_path):
    """The branch the empty-directory test does NOT reach, and the sweep found it uncovered.

    An empty directory returns None one branch earlier, at "no cleared records". The branch
    that matters is: records load, clear the gate, and still produce no APPLICABLE
    correction. Serving an empty model in that case would make /harmonize answer 200 with
    every value absent, which is a refusal wearing a success code.

    The Struwe seed files are exactly that data - all travelling-wave from one laboratory,
    so nothing pairs across platforms and there is no stratum to anchor on the primary.
    """
    import shutil

    source = SEED_DIRECTORY / "struwe2016_chemcommun.csv"
    assert source.exists(), "this test needs the one-platform seed file"
    shutil.copy(source, tmp_path / source.name)

    # It loads and it clears: this is not an empty-directory case in disguise.
    report = load_measurements_file(tmp_path / source.name)
    assert report.records_cleared > 0

    assert build_default_model(tmp_path) is None


# --- 200: what a harmonized answer carries --------------------------------------------------


def test_a_covered_measurement_comes_back_harmonized_with_everything_that_must_travel(served, records):
    response = served.post("/harmonize", json=body_for(a_covered_record(records)))
    assert response.status_code == 200
    measurement = response.json()["measurements"][0]
    estimate = measurement["harmonized"]

    assert estimate is not None
    assert measurement["not_harmonized_because"] is None
    assert estimate["ccs"] > 0
    assert estimate["basis"] in {CorrectionBasis.MEDIAN.value, CorrectionBasis.ROBUST_SLOPE.value}
    # Both corrections, always, because they diverge where the fit is levered.
    assert estimate["slope_derived_ccs"] > 0
    assert estimate["median_derived_ccs"] > 0
    assert estimate["interval_low"] < estimate["ccs"] < estimate["interval_high"]
    assert measurement["confidence"]["grade"] in {grade.value for grade in ConfidenceGrade}
    assert measurement["provenance"]["doi"] == "10.1021/jasms.2c00196"


def test_both_coverage_figures_travel_and_the_interval_names_its_method(served, records):
    """The one that would be easiest to get quietly wrong.

    Jackknife+ proves 1-2*alpha. A response carrying only the nominal 0.90 would be read
    as a 90 per cent guarantee, and labelling the interval `conformal_prediction` would
    imply split conformal's stronger claim.
    """
    estimate = served.post("/harmonize", json=body_for(a_covered_record(records))).json()[
        "measurements"
    ][0]["harmonized"]
    assert estimate["interval_coverage"] == pytest.approx(0.90)
    assert estimate["guaranteed_coverage"] == pytest.approx(0.80)
    assert estimate["interval_kind"] == IntervalKind.JACKKNIFE_PLUS.value
    assert estimate["interval_is_informative"] is True


def test_every_harmonized_estimate_carries_its_scope_and_the_scope_says_within_study(served, records):
    estimate = served.post("/harmonize", json=body_for(a_covered_record(records))).json()[
        "measurements"
    ][0]["harmonized"]
    scope = estimate["scope"]
    assert scope["scope"] == ComparisonScope.WITHIN_STUDY.value
    assert scope["studies"] == ["doi:10.1021/jasms.2c00196"]
    assert "NOT interlaboratory reproducibility" in scope["caveat"]


def test_the_original_is_echoed_exactly_beside_the_harmonized_value(served, records):
    record = a_covered_record(records)
    measurement = served.post("/harmonize", json=body_for(record)).json()["measurements"][0]
    assert measurement["original"]["ccs"] == record.ccs
    assert measurement["original"]["adduct"] == record.adduct
    assert measurement["harmonized"]["ccs"] != record.ccs, "a correction that changes nothing is not one"


def test_the_maturity_on_a_200_is_still_provisional(served, records):
    payload = served.post("/harmonize", json=body_for(a_covered_record(records))).json()
    assert payload["maturity"]["data_maturity"] == DataMaturity.PROVISIONAL.value
    assert payload["maturity"]["matched_ion_count"] == 142


# --- 501: what is refused, and why ------------------------------------------------------------


def test_a_measurement_already_on_the_primary_platform_is_refused_with_its_reason(served, records):
    response = served.post("/harmonize", json=body_for(a_primary_record(records)))
    assert response.status_code == 501
    measurement = response.json()["measurements"][0]
    assert measurement["harmonized"] is None
    assert "already on the primary" in measurement["not_harmonized_because"]
    assert measurement["original"]["ccs"] == a_primary_record(records).ccs


def test_a_request_the_model_covers_nothing_in_is_a_501_and_says_so(served, records):
    response = served.post("/harmonize", json=body_for(a_primary_record(records)))
    payload = response.json()
    assert response.status_code == 501
    assert payload["error"] == "not_implemented"
    assert "covers none" in payload["detail"]


def test_a_mixed_request_is_a_200_and_each_absent_value_says_why(served, records):
    """The case a per-request status code alone could not express.

    One covered measurement and one not: 200, because something was answered, with the
    uncovered one carrying its reason rather than a number.
    """
    response = served.post(
        "/harmonize", json=body_for(a_covered_record(records), a_primary_record(records))
    )
    assert response.status_code == 200
    covered, uncovered = response.json()["measurements"]
    assert covered["harmonized"] is not None
    assert covered["not_harmonized_because"] is None
    assert uncovered["harmonized"] is None
    assert uncovered["not_harmonized_because"].strip()


def test_an_app_with_no_model_refuses_everything_and_returns_the_originals(records):
    client = TestClient(create_app(None))
    response = client.post("/harmonize", json=body_for(a_covered_record(records)))
    assert response.status_code == 501
    payload = response.json()
    assert "No harmonization model is loaded" in payload["detail"]
    measurement = payload["measurements"][0]
    assert measurement["harmonized"] is None
    assert measurement["original"]["ccs"] == a_covered_record(records).ccs
    assert payload["maturity"]["matched_ion_count"] == 0


def test_nothing_is_extrapolated_a_value_graded_unsupported_is_withheld(served, records):
    """A grade of `unsupported` means do not use this number, so the number is not returned.

    Handing over a value while saying not to use it is a contradiction a caller resolves in
    favour of the value. Built with an ion far outside the fitted range, because the real
    corpus is fitted over its own range and produces no such case on its own.
    """
    far = a_covered_record(records).model_copy(update={"ccs": 900.0})
    response = served.post("/harmonize", json=body_for(far))
    assert response.status_code == 501
    measurement = response.json()["measurements"][0]
    assert measurement["harmonized"] is None
    assert measurement["confidence"]["grade"] == ConfidenceGrade.UNSUPPORTED.value
    assert "do not use the number" in measurement["not_harmonized_because"]
    # the grade's reasons are still returned, so a caller learns WHY it was withheld
    assert any("calibration range" in reason["rule"] for reason in measurement["confidence"]["reasons"])


def test_a_malformed_request_is_still_a_422_and_not_a_refusal(served):
    """Validation happens before any of this. A 501 swallowing bad input would mislead."""
    assert served.post("/harmonize", json={"measurements": []}).status_code == 422
    assert served.post("/harmonize", json={"nonsense": 1}).status_code == 422
