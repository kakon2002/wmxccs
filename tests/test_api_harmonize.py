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


def test_the_estimate_states_that_the_interval_ignores_the_submitted_uncertainty(served, records):
    """Half of the ruling of 20 September 2026, and the half that is only a disclosure.

    The response carries the caller's own `ccs_uncertainty` in `original` a few fields above
    an interval that has nothing to do with it, which invites exactly the wrong reading. The
    field says so rather than leaving it to be discovered.
    """
    estimate = served.post("/harmonize", json=body_for(a_covered_record(records))).json()[
        "measurements"
    ][0]["harmonized"]
    assert estimate["interval_accounts_for_submitted_uncertainty"] is False


def test_the_interval_really_is_independent_of_the_submitted_uncertainty(served, records):
    """The field above ASSERTS something; this proves the assertion is true.

    A disclosure nothing checks is the same class of defect as the unevaluable rule it was
    ruled on beside: a statement in the response that no code is holding to. If propagation
    is ever implemented, this test fails and the field must stop saying False - which is the
    point of writing it this way round.
    """
    record = a_covered_record(records)
    intervals = []
    for uncertainty in (0.0001, 50.0):
        body = json.loads(record.model_dump_json())
        body["ccs_uncertainty"] = uncertainty
        body["uncertainty_type"] = "sd"
        estimate = served.post("/harmonize", json={"measurements": [body]}).json()[
            "measurements"
        ][0]["harmonized"]
        assert estimate is not None, "this test needs a harmonized value to compare"
        assert estimate["interval_accounts_for_submitted_uncertainty"] is False
        intervals.append((estimate["interval_low"], estimate["ccs"], estimate["interval_high"]))
    assert intervals[0] == intervals[1], (
        "the interval moved with the submitted uncertainty, so the response's claim that it"
        f" does not is false: {intervals[0]} against {intervals[1]}"
    )


def test_a_relabelled_outlier_is_served_but_the_response_stops_claiming_it_was_checked(
    served, records
):
    """A KNOWN AND ACCEPTED RESIDUAL, pinned so it cannot change without somebody deciding to.

    The outlier rule is the only rule keyed on the analyte's identity, and identity is
    supplied by the caller. So a record the service refuses under its own name can still be
    resubmitted under a new one and receive a value - the ruling of 20 September 2026 is
    explicit that an unevaluable rule demotes one notch rather than refusing, because a
    genuinely new ion may be fine and refusing every new ion refuses the platform's purpose.

    What the fix changed is the CLAIM, not the number. This test asserts both halves: the
    value still comes back, AND the response no longer reports the outlier check as having
    passed. If the first half ever starts failing, the platform has begun refusing novel
    ions. If the second half ever starts failing, the original defect is back.
    """
    refused = None
    for record in records:
        measurement = served.post("/harmonize", json=body_for(record)).json()["measurements"][0]
        confidence = measurement.get("confidence")
        if confidence and confidence["grade"] == ConfidenceGrade.UNSUPPORTED.value:
            flagged = "flagged as an outlier in its own stratum" in {
                reason["rule"] for reason in confidence["reasons"]
            }
            if flagged:
                refused = record
                break
    assert refused is not None, "this test needs a record the corpus refuses as an outlier"
    assert refused.analyte.dataset_compound_id is not None

    body = body_for(refused)
    body["measurements"][0]["analyte"]["dataset_compound_id"] = "customerlab:brand-new-compound"
    body["measurements"][0]["analyte"]["display_name"] = "never seen before"
    renamed = served.post("/harmonize", json=body).json()["measurements"][0]

    # the residual: a value the service just refused, served under another name
    assert renamed["harmonized"] is not None
    # and the disclosure that is the whole of the fix
    assert renamed["confidence"]["grade"] != ConfidenceGrade.UNSUPPORTED.value
    assert [note for note in renamed["confidence"]["not_checked"] if note.startswith("outlier check:")]
    assert "a rule of the scheme could not be evaluated for this ion" in {
        reason["rule"] for reason in renamed["confidence"]["reasons"]
    }
    assert "flagged as an outlier in its own stratum" not in {
        reason["rule"] for reason in renamed["confidence"]["reasons"]
    }


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


# --- model versioning: what makes an answer repeatable ----------------------------------------
#
# The model is refitted from the seed files at every startup. That is fine. What is not fine
# is a changed seed file changing the answers with nothing recording it: two deployments
# could give different numbers for one input and neither response would say so. These tests
# are the difference between a demonstration and something somebody can rely on twice.


def test_the_fingerprint_is_identical_across_two_fits_of_the_same_data():
    """Deterministic, or it identifies nothing.

    sha256 over a sorted canonical form with floats written by repr. No paths, no timestamps,
    no dict ordering - all three would make two fits of one corpus look like two models.
    """
    first = build_default_model(SEED_DIRECTORY)
    second = build_default_model(SEED_DIRECTORY)
    assert first.fingerprint == second.fingerprint
    assert first.fingerprint.corpus == second.fingerprint.corpus
    assert first.fingerprint.parameters == second.fingerprint.parameters


def test_the_same_data_reached_by_a_different_path_gives_the_same_fingerprint(tmp_path):
    """Over the RECORDS, not the files, so a copy of the corpus is the same corpus."""
    import shutil

    for source in SEED_DIRECTORY.glob("*.csv"):
        shutil.copy(source, tmp_path / source.name)
    assert build_default_model(tmp_path).fingerprint == build_default_model(SEED_DIRECTORY).fingerprint


def test_changing_one_cross_section_by_a_thousandth_changes_both_digests(tmp_path):
    """The property the whole mechanism exists for.

    A digest that did not move on a changed measurement would be decoration. Both move: the
    corpus because a record changed, the parameters because the fit did.
    """
    import csv
    import shutil

    for source in SEED_DIRECTORY.glob("*.csv"):
        shutil.copy(source, tmp_path / source.name)
    before = build_default_model(tmp_path).fingerprint

    target = tmp_path / "steroid_jasms2022.csv"
    rows = list(csv.DictReader(target.open(newline="", encoding="utf-8")))
    rows[0]["ccs"] = str(float(rows[0]["ccs"]) + 0.001)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    after = build_default_model(tmp_path).fingerprint
    assert not after.same_corpus_as(before), "a changed measurement must change the corpus digest"
    assert not after.same_parameters_as(before), "and the fit it feeds"


def test_dropping_measurements_that_feed_a_correction_changes_the_fingerprint(tmp_path):
    """Removing data the fit uses is a change, and a digest that missed it would be worse than none."""
    import csv
    import shutil

    for source in SEED_DIRECTORY.glob("*.csv"):
        shutil.copy(source, tmp_path / source.name)
    full = build_default_model(tmp_path).fingerprint

    target = tmp_path / "steroid_jasms2022.csv"
    rows = list(csv.DictReader(target.open(newline="", encoding="utf-8")))
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator=chr(10))
        writer.writeheader()
        writer.writerows(rows[:-40])

    assert build_default_model(tmp_path).fingerprint != full


def test_removing_data_that_feeds_no_correction_leaves_the_fingerprint_alone(tmp_path):
    """Deliberate, and the opposite of what the name "corpus" suggests.

    The Struwe seed files are all travelling-wave from one laboratory, so they pair nothing,
    enter no stratum and can change no answer. Removing them therefore changes nothing a
    caller could observe, and the digest says so.

    This was written the other way round first, on the assumption that a digest called
    "corpus" covers everything on disk. It does not, and the weaker promise would be worse:
    hashing every file would make two identical answers look like they came from different
    models whenever unrelated data moved.
    """
    import shutil

    for source in SEED_DIRECTORY.glob("*.csv"):
        shutil.copy(source, tmp_path / source.name)
    full = build_default_model(tmp_path).fingerprint

    struwe = tmp_path / "struwe2016_chemcommun.csv"
    assert struwe.exists()
    struwe.unlink()

    assert build_default_model(tmp_path).fingerprint == full


def test_the_two_digests_answer_different_questions():
    """Which is why there are two rather than one combined hash.

    A different corpus with the same parameters means the data moved without moving the fit;
    the same corpus with different parameters means the code did. One hash says only that
    something changed.
    """
    model = build_default_model(SEED_DIRECTORY)
    assert model.fingerprint.corpus != model.fingerprint.parameters
    assert len(model.fingerprint.corpus) == 64
    assert len(model.fingerprint.parameters) == 64
    assert model.fingerprint.short.count("/") == 1


def test_health_serves_the_version_of_the_model_it_is_running(served):
    payload = served.get("/health").json()
    expected = build_default_model(SEED_DIRECTORY).fingerprint
    assert payload["model_version"]["corpus_sha256"] == expected.corpus
    assert payload["model_version"]["parameters_sha256"] == expected.parameters


def test_health_reports_no_version_where_no_model_is_loaded():
    """Absent, not a placeholder hash. There is no model to identify."""
    payload = TestClient(create_app(None)).get("/health").json()
    assert payload["model_loaded"] is False
    assert payload["model_version"] is None


def test_every_harmonized_estimate_carries_the_version_that_produced_it(served, records):
    """Required with no default: an answer that cannot say which model made it cannot be repeated."""
    estimate = served.post("/harmonize", json=body_for(a_covered_record(records))).json()[
        "measurements"
    ][0]["harmonized"]
    expected = build_default_model(SEED_DIRECTORY).fingerprint
    assert estimate["model_version"]["corpus_sha256"] == expected.corpus
    assert estimate["model_version"]["parameters_sha256"] == expected.parameters


def test_the_version_on_an_estimate_matches_the_one_on_health(served, records):
    """So a caller can tell whether the answer they hold came from the model now running."""
    health = served.get("/health").json()["model_version"]
    estimate = served.post("/harmonize", json=body_for(a_covered_record(records))).json()[
        "measurements"
    ][0]["harmonized"]["model_version"]
    assert health == estimate


# --- the synthetic refusal must not depend on whether a model is loaded -----------------------
#
# Caught by the mutation sweep, not by a test. The synthetic-fixture test in test_api.py runs
# against an app with NO model, where `model is None` is already true - so a mutation making
# the refusal conditional on that changed nothing there and survived. The objection is to the
# RECORD, so it has to be asserted where a model IS loaded.


def test_a_synthetic_fixture_is_refused_even_when_a_model_is_loaded(served, records):
    """With a model present there is a correction available, and it is still refused.

    This is the half test_api.py cannot reach: against a model-less app the record would be
    refused anyway for want of a model, so that test cannot tell the two reasons apart.
    """
    record = a_covered_record(records)
    body = body_for(record)
    body["measurements"][0]["reuse_status"] = "synthetic_fixture"

    response = served.post("/harmonize", json=body)
    assert response.status_code == 501
    measurement = response.json()["measurements"][0]
    assert measurement["harmonized"] is None
    assert "may not carry synthetic_fixture" in measurement["not_harmonized_because"]
    assert "integrity of the answer" in measurement["not_harmonized_because"]

    # And the identical record WITHOUT that status is harmonized by this same app, so the
    # refusal is attributable to the status rather than to anything else about the record.
    clean = served.post("/harmonize", json=body_for(record))
    assert clean.status_code == 200
    assert clean.json()["measurements"][0]["harmonized"] is not None
