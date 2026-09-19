"""The HTTP contract: what comes back untouched, and what does not come back at all.

TWO PROMISES ARE ASSERTED HERE AND NOTHING ELSE MATTERS AS MUCH.

The first is that THE ORIGINALS COME BACK EXACTLY AS THEY WERE SENT. Not rounded,
not normalised, not re-expressed, and never replaced by a corrected value. So the
round trip below is asserted on the WHOLE record rather than on its cross section:
a response that kept `ccs` and quietly dropped a conformer index, a measurement
date or an uncertainty's type would satisfy a test that only read the number, and
would have lost the parts that say what the number means.

The second is that WHILE NO MODEL EXISTS, /harmonize RETURNS NO NUMBER IT WAS NOT
GIVEN. Not a zero cross section, not an interval of infinite width, not a default
grade computed against nothing. That is asserted by walking the whole response and
requiring that every number in it is either an echo of the input or the count of
matched ions that exist, which is zero because none do. An API that returned a
plausible value with a caveat in a field would be used and the caveat would not be
read.

Provenance is the third thing, and the pair of tests on it is the point rather
than either half: the API reports what the REGISTRY says about a source AND what
the record CLAIMS, so a reader can see when the two disagree. A record may claim
any reuse status it likes; only the registry entry names somebody who read the
terms and when.
"""

from __future__ import annotations

import json
from copy import deepcopy

import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError

from wmxccs import sources
from wmxccs.api import CONFIDENCE_NOTE, create_app
from wmxccs.scope import ComparisonScope
from wmxccs.contracts import (
    ScopeReport,
    ConfidenceReport,
    CorrectionBasis,
    HarmonizationUnavailable,
    HarmonizedEstimate,
    HarmonizedMeasurement,
    HarmonizeRequest,
    HarmonizeResponse,
    IntervalKind,
    SourceProvenance,
)
from wmxccs.grading import ConfidenceGrade, grade_rules
from wmxccs.models import CCSMeasurement, UncertaintyType
from wmxccs.readiness import DataMaturity
from wmxccs.reuse import ReuseStatus

from conftest import measurement

# A DOI the licence registry holds, with somebody's name and date against it.
REGISTERED_DOI = "10.1039/c6cc06247d"
# Not allocated to any registrant, so it can never be in the registry.
UNREGISTERED_DOI = "10.9999/nobody.has.registered.this"


@pytest.fixture
def app():
    return create_app()


@pytest.fixture
def client(app):
    return TestClient(app)


def awkward_record(**overrides) -> CCSMeasurement:
    """A record whose every field would show rounding or re-expression if either happened.

    The cross section carries nine decimal places, the uncertainty carries its
    type, the row says which of three conformers it is, and it has a measurement
    date. A response that re-expressed any of them comes back unequal.
    """
    fields = dict(
        ccs=187.123456789,
        ccs_uncertainty=0.987654321,
        uncertainty_type=UncertaintyType.TWO_SD,
        conformer=2,
        conformers_total=3,
        measured_on="2026-03-04",
        source_locator="Table S2",
        replicates=5,
        instrument="a fixture instrument",
        doi=UNREGISTERED_DOI,
    )
    fields.update(overrides)
    return measurement(**fields)


def body_for(*records) -> dict:
    """A request body built from real records, the way a caller would build one."""
    return {"measurements": [json.loads(record.model_dump_json()) for record in records]}


def numbers_in(value, path=()):
    """Every number in a decoded JSON tree, with the path it sits at."""
    if isinstance(value, bool):
        return
    if isinstance(value, (int, float)):
        yield path, value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from numbers_in(item, path + (key,))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from numbers_in(item, path + (index,))


# --- 1. the originals come back untouched --------------------------------------------------------


def test_a_posted_measurement_comes_back_as_an_equal_record_and_not_merely_an_equal_ccs(client):
    """The third hard constraint, and the reason this endpoint exists in this state."""
    record = awkward_record()
    response = client.post("/harmonize", json=body_for(record))
    assert response.status_code == 501
    returned = response.json()["measurements"]
    assert len(returned) == 1
    round_tripped = CCSMeasurement.model_validate(returned[0]["original"])
    assert round_tripped == record
    # Spelled out as well as compared, because an equality that started returning
    # True for everything would pass the line above in silence.
    assert round_tripped.ccs == record.ccs == 187.123456789
    assert round_tripped.ccs_uncertainty == record.ccs_uncertainty
    assert round_tripped.uncertainty_type is UncertaintyType.TWO_SD
    assert round_tripped.conformer == 2 and round_tripped.conformers_total == 3
    assert round_tripped.measured_on == record.measured_on
    assert round_tripped.analyte == record.analyte


def test_the_returned_original_is_not_rounded_or_re_expressed_in_the_json_itself(client):
    # The comparison above goes through the model, which would forgive a value
    # re-expressed as a string. This one reads the wire.
    record = awkward_record()
    original = client.post("/harmonize", json=body_for(record)).json()["measurements"][0]["original"]
    assert original["ccs"] == 187.123456789
    assert original["ccs_uncertainty"] == 0.987654321
    assert original["uncertainty_type"] == UncertaintyType.TWO_SD.value
    assert original["measured_on"] == "2026-03-04"
    assert original["conformer"] == 2
    assert original == json.loads(record.model_dump_json())


def test_every_measurement_posted_comes_back_in_the_order_it_was_sent(client):
    first = awkward_record(ccs=187.123456789)
    second = awkward_record(ccs=211.987654321, conformer=1, conformers_total=3)
    response = client.post("/harmonize", json=body_for(first, second))
    returned = response.json()["measurements"]
    assert len(returned) == 2
    assert [CCSMeasurement.model_validate(item["original"]) for item in returned] == [first, second]


def test_the_harmonized_value_is_an_additional_field_and_never_replaces_the_original(client):
    record = awkward_record()
    returned = client.post("/harmonize", json=body_for(record)).json()["measurements"][0]
    assert set(returned) == {"original", "provenance", "harmonized", "confidence"}
    assert returned["harmonized"] is None
    assert returned["original"]["ccs"] == record.ccs


# --- 2. the refusal carries no number it was not given --------------------------------------------


def test_harmonize_refuses_with_501_and_names_the_refusal(client):
    response = client.post("/harmonize", json=body_for(awkward_record()))
    assert response.status_code == 501
    payload = response.json()
    assert payload["error"] == "not_implemented"
    assert payload["detail"].strip()


def test_the_refusal_carries_no_harmonized_value_and_no_grade(client):
    payload = client.post("/harmonize", json=body_for(awkward_record(), awkward_record(ccs=300.5))).json()
    for item in payload["measurements"]:
        # Absent, not empty. An empty estimate would be a shape a caller could
        # read a zero out of.
        assert item["harmonized"] is None
        assert item["confidence"] is None


def test_no_number_in_the_refusal_is_anything_but_an_echo_of_the_input(client):
    """No placeholder anywhere: not a zero cross section, not an infinite interval.

    The originals are the only thing in the body that may carry a number, so they
    are stripped and what remains is required to hold none, apart from the count
    of matched ions that exist, which is zero because there are none.
    """
    payload = client.post("/harmonize", json=body_for(awkward_record())).json()
    stripped = deepcopy(payload)
    for item in stripped["measurements"]:
        del item["original"]
    found = sorted(numbers_in(stripped))
    assert found == [(("maturity", "matched_ion_count"), 0)], f"an unexplained number appeared: {found}"
    assert not any(isinstance(value, float) for _path, value in found)


def test_no_default_grade_and_no_interval_creeps_into_the_refusal(client):
    payload = client.post("/harmonize", json=body_for(awkward_record())).json()
    stripped = deepcopy(payload)
    for item in stripped["measurements"]:
        del item["original"]
    rendered = json.dumps(stripped)
    for grade in ConfidenceGrade:
        assert grade.value not in rendered, f"a {grade.value} grade appeared with nothing behind it"
    for word in ("Infinity", "-Infinity", "NaN", "interval_low", "interval_high", "interval_coverage"):
        assert word not in rendered
    assert "harmonized" in json.dumps(payload)  # the key is there; the value is null


# --- 3 and 4. the refusal does not swallow a bad request -------------------------------------------


def test_a_malformed_measurement_is_refused_with_422_rather_than_501(client):
    """A 501 that swallowed a bad body would tell a caller their record was fine."""
    response = client.post("/harmonize", json={"measurements": [{"ccs": "not a number"}]})
    assert response.status_code == 422
    assert response.status_code != 501


def test_a_body_that_is_not_a_request_at_all_is_refused_with_422(client):
    assert client.post("/harmonize", json={"nothing": "useful"}).status_code == 422
    assert client.post("/harmonize", json=[]).status_code == 422


def test_a_measurement_carrying_a_spread_with_no_stated_type_is_refused_with_422(client):
    # The seventh hard constraint reaching the wire: a number without its type is
    # refused at the edge rather than stored and puzzled over later.
    body = body_for(awkward_record())
    body["measurements"][0]["uncertainty_type"] = None
    assert client.post("/harmonize", json=body).status_code == 422


def test_a_request_with_no_measurements_at_all_is_refused(client):
    """min_length=1. A request to harmonize nothing is a mistake, not an empty answer."""
    response = client.post("/harmonize", json={"measurements": []})
    assert response.status_code == 422
    with pytest.raises(ValidationError):
        HarmonizeRequest(measurements=())


def test_a_request_holding_one_measurement_is_accepted(client):
    """The other direction, so the refusal above is the length rule and not a broken endpoint."""
    assert client.post("/harmonize", json=body_for(awkward_record())).status_code == 501
    assert len(HarmonizeRequest(measurements=(awkward_record(),)).measurements) == 1


# --- 5. provenance comes from the registry, and the record's claim is reported beside it ------------


def test_provenance_for_a_registered_source_is_read_from_the_registry(client):
    record = awkward_record(doi=REGISTERED_DOI)
    provenance = client.post("/harmonize", json=body_for(record)).json()["measurements"][0]["provenance"]
    assert provenance["registered"] is True
    assert provenance["licence"] == sources.STRUWE_2016.licence
    assert provenance["licence_reported_by"] == sources.STRUWE_2016.reported_by
    assert provenance["licence_reported_on"] == sources.STRUWE_2016.reported_on.isoformat()
    assert provenance["doi"] == REGISTERED_DOI
    assert provenance["source"] == record.source


def test_provenance_for_an_unregistered_source_says_so_and_invents_no_licence(client):
    record = awkward_record(doi=UNREGISTERED_DOI)
    assert sources.licence_for(UNREGISTERED_DOI) is None
    provenance = client.post("/harmonize", json=body_for(record)).json()["measurements"][0]["provenance"]
    assert provenance["registered"] is False
    assert provenance["licence"] is None
    assert provenance["licence_reported_by"] is None
    assert provenance["licence_reported_on"] is None


def test_the_records_own_reuse_claim_is_echoed_whether_or_not_the_registry_agrees(client):
    """The pair is the point: a reader can see the claim and the registry disagree.

    The record below claims to be a synthetic fixture while the registry says the
    DOI it cites is CC BY 3.0, read by a named person on a stated date. Both are
    reported, and neither is allowed to overwrite the other.
    """
    claimed = ReuseStatus.SYNTHETIC_FIXTURE
    registered = awkward_record(doi=REGISTERED_DOI, reuse_status=claimed)
    unregistered = awkward_record(doi=UNREGISTERED_DOI, reuse_status=claimed)
    returned = client.post("/harmonize", json=body_for(registered, unregistered)).json()["measurements"]
    assert returned[0]["provenance"]["reuse_status"] == claimed.value
    assert returned[1]["provenance"]["reuse_status"] == claimed.value
    assert returned[0]["provenance"]["registered"] is True
    assert returned[1]["provenance"]["registered"] is False
    assert returned[0]["provenance"]["licence"] != returned[1]["provenance"]["licence"]
    assert sources.STRUWE_2016.reuse_status is not claimed


def test_a_record_with_no_doi_at_all_is_reported_as_unregistered(client):
    provenance = client.post("/harmonize", json=body_for(awkward_record(doi=None))).json()[
        "measurements"
    ][0]["provenance"]
    assert provenance["doi"] is None
    assert provenance["registered"] is False
    assert provenance["licence"] is None


# --- 6. nothing claims to be more mature than it is ------------------------------------------------


def test_health_reports_no_model_and_no_matched_ions(client):
    payload = client.get("/health").json()
    assert payload["status"] == "ok"
    assert payload["model_loaded"] is False
    assert payload["matched_ions_available"] == 0
    assert payload["version"]


def test_health_would_report_a_model_if_one_were_loaded(app):
    """The other direction, so the flag is reading the app rather than returning a constant."""
    app.state.model = object()
    app.state.matched_ion_count = 7
    payload = TestClient(app).get("/health").json()
    assert payload["model_loaded"] is True
    assert payload["matched_ions_available"] == 7


def test_the_maturity_stamp_on_the_refusal_is_provisional_over_zero_matched_ions(client):
    maturity = client.post("/harmonize", json=body_for(awkward_record())).json()["maturity"]
    assert maturity["data_maturity"] == DataMaturity.PROVISIONAL.value
    assert maturity["matched_ion_count"] == 0


def test_the_maturity_stamp_only_says_validated_once_something_has_been_checked(app):
    """The other direction. A stamp that always said provisional would pass the test above."""
    app.state.model_validated = True
    maturity = TestClient(app).post("/harmonize", json=body_for(awkward_record())).json()["maturity"]
    assert maturity["data_maturity"] == DataMaturity.VALIDATED.value


# --- 7. the grading scheme is published and needs no model -------------------------------------------


def test_the_confidence_rules_answer_on_a_fresh_app_with_no_model(app):
    """The endpoint exists precisely because these rules need no training data."""
    assert app.state.model is None
    assert app.state.matched_ion_count == 0
    response = TestClient(app).get("/confidence/rules")
    assert response.status_code == 200


def test_the_confidence_rules_list_every_grade_and_every_rule(client):
    payload = client.get("/confidence/rules").json()
    assert payload["grades"] == [grade.value for grade in ConfidenceGrade]
    assert len(payload["rules"]) == len(grade_rules())
    assert [rule["rule"] for rule in payload["rules"]] == [rule["rule"] for rule in grade_rules()]
    for served, declared in zip(payload["rules"], grade_rules()):
        assert served["falls_to"] == list(declared["falls_to"])
        assert served["why"] == declared["why"]
        assert served["threshold"] == declared["threshold"]
        assert served["basis"] == declared["basis"]


def test_the_published_note_says_the_grade_only_ever_falls(client):
    note = client.get("/confidence/rules").json()["note"]
    assert note == CONFIDENCE_NOTE
    assert "only ever falls" in note
    assert "worst demotion" in note
    assert "never an average" in note


# --- 8. the 200 shape is specified and is not what any endpoint returns -------------------------------


def test_the_eventual_200_shape_is_specified_for_callers_building_ahead_of_the_model():
    assert set(HarmonizeResponse.model_fields) == {"measurements", "maturity"}
    assert set(HarmonizedMeasurement.model_fields) == {
        "original",
        "provenance",
        "harmonized",
        "confidence",
    }
    assert set(ConfidenceReport.model_fields) == {"grade", "reasons", "not_checked"}
    assert {"ccs", "basis", "interval_low", "interval_high", "interval_coverage", "interval_kind"} <= set(
        HarmonizedEstimate.model_fields
    )
    assert set(SourceProvenance.model_fields) == {
        "source",
        "doi",
        "reuse_status",
        "licence",
        "licence_reported_by",
        "licence_reported_on",
        "registered",
    }


def test_no_endpoint_returns_the_200_shape_today(app):
    """/harmonize answers with the refusal shape, which is a different shape on purpose.

    Declaring HarmonizeResponse here would let a caller write code against a
    harmonized value that is never present, and the schema would say the value was
    merely optional rather than absent by policy.
    """
    harmonize = [route for route in app.routes if getattr(route, "path", None) == "/harmonize"]
    assert len(harmonize) == 1
    assert harmonize[0].response_model is HarmonizationUnavailable
    assert harmonize[0].response_model is not HarmonizeResponse
    assert harmonize[0].status_code == 501
    declared = {getattr(route, "response_model", None) for route in app.routes}
    assert HarmonizeResponse not in declared


# --- 9. the contract refuses a placeholder estimate at the model level ---------------------------------


def good_scope(**overrides) -> dict:
    """The scope every estimate must carry. Within-study, because that is what the corpus is."""
    fields = dict(
        scope=ComparisonScope.WITHIN_STUDY,
        studies=("doi:10.1021/jasms.2c00196",),
        platforms=("DTIMS/stepped_field", "TWIMS"),
        records_behind_it=62,
        caveat="WITHIN ONE STUDY. Not interlaboratory reproducibility.",
    )
    fields.update(overrides)
    return fields


def good_estimate(**overrides) -> dict:
    """A well-formed estimate, so each refusal below is about the one field it changes."""
    fields = dict(
        ccs=187.5,
        basis=CorrectionBasis.MEDIAN,
        interval_low=185.0,
        interval_high=190.0,
        interval_coverage=0.90,
        interval_kind=IntervalKind.LIMITS_OF_AGREEMENT,
        reference_platform="DTIMS/stepped_field",
        matched_ions_behind_it=24,
        scope=ScopeReport(**good_scope()),
        interval_is_informative=True,
        guaranteed_coverage=0.80,
    )
    fields.update(overrides)
    return fields


# --- the scope is REQUIRED and cannot be widened by declaring it -------------------------------


def test_a_harmonized_estimate_without_a_scope_is_refused():
    """No default. A number whose scope could be omitted is a number quoted without it."""
    fields = good_estimate()
    del fields["scope"]
    with pytest.raises(ValidationError):
        HarmonizedEstimate(**fields)


def test_a_scope_report_may_not_claim_more_than_its_studies_support():
    """The scope is DERIVED from the provenance, so declaring a wider one is refused.

    This is the structural half of the requirement: widening the claim means naming a
    second study, which is data somebody has to produce, not a field somebody can set.
    """
    with pytest.raises(ValidationError, match="does not follow from"):
        ScopeReport(**good_scope(scope=ComparisonScope.CROSS_STUDY))
    # and naming two studies makes it legitimate
    across = ScopeReport(
        **good_scope(
            scope=ComparisonScope.CROSS_STUDY,
            studies=("doi:10.1021/jasms.2c00196", "doi:10.1021/acs.analchem.9b05247"),
        )
    )
    assert across.scope is ComparisonScope.CROSS_STUDY
    # the other direction too: two studies may not be called within-study
    with pytest.raises(ValidationError, match="does not follow from"):
        ScopeReport(
            **good_scope(
                scope=ComparisonScope.WITHIN_STUDY,
                studies=("doi:a/1", "doi:b/2"),
            )
        )


def test_a_scope_report_needs_at_least_one_study_and_one_platform():
    with pytest.raises(ValidationError):
        ScopeReport(**good_scope(studies=()))
    with pytest.raises(ValidationError):
        ScopeReport(**good_scope(platforms=()))


def test_no_scope_member_means_interlaboratory_reproducibility():
    """The claim is not representable, so it cannot be recorded or returned."""
    assert {member.value for member in ComparisonScope} == {"within_study", "cross_study"}


def test_a_well_formed_harmonized_estimate_is_accepted():
    """The control. Without it every refusal below would be satisfied by a model that refused all."""
    estimate = HarmonizedEstimate(**good_estimate())
    assert estimate.ccs == pytest.approx(187.5, abs=1e-9)
    assert estimate.interval_coverage == pytest.approx(0.90, abs=1e-9)
    assert estimate.matched_ions_behind_it == 24


@pytest.mark.parametrize("ccs", [0.0, -1.0, -187.5])
def test_a_harmonized_cross_section_of_zero_or_less_is_refused(ccs):
    with pytest.raises(ValidationError):
        HarmonizedEstimate(**good_estimate(ccs=ccs))


def test_an_interval_coverage_of_exactly_one_is_refused_because_no_finite_interval_attains_it():
    with pytest.raises(ValidationError):
        HarmonizedEstimate(**good_estimate(interval_coverage=1.0))
    with pytest.raises(ValidationError):
        HarmonizedEstimate(**good_estimate(interval_coverage=0.0))
    # And the other direction: a coverage strictly inside the interval is fine.
    assert HarmonizedEstimate(**good_estimate(interval_coverage=0.999)).interval_coverage == pytest.approx(
        0.999, abs=1e-9
    )


def test_a_harmonized_estimate_must_say_how_many_matched_ions_are_behind_it():
    """The count travels with the number. A default of zero would be a number nobody chose."""
    fields = good_estimate()
    del fields["matched_ions_behind_it"]
    with pytest.raises(ValidationError):
        HarmonizedEstimate(**fields)
    with pytest.raises(ValidationError):
        HarmonizedEstimate(**good_estimate(matched_ions_behind_it=-1))
    assert HarmonizedEstimate(**good_estimate(matched_ions_behind_it=0)).matched_ions_behind_it == 0
