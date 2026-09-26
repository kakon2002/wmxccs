"""The dashboard: a real client of the six endpoints, and it must not manufacture certainty.

A browser cannot be run here, so what is tested is everything short of rendering:

  - the page is served by the service itself, same origin as the endpoints it calls;
  - EVERY field the JavaScript reads off a response is a real field of a real response model,
    derived from the models on one side and from the script on the other. A mistyped field name
    in a client renders `undefined` in a browser and fails nothing, which is the worst kind of
    silence: the page looks finished and one number is missing;
  - the candidate view never numbers a candidate, because 152 of 167 Hex5HexNAc4Fuc1 candidates
    sit in a class the platform's own features cannot separate and a list running 1 to 167 would
    present an arbitrary order as a ranking;
  - an unevaluable comparison axis shows its REASON and never a blank. A dash in a delta column
    reads as "about zero" to anyone skimming.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from pydantic import BaseModel

import wmxglycan.ccs_evidence as evidence_models
import wmxglycan.prediction as contracts
from wmxglycan.api import DASHBOARD, create_app
from wmxglycan.attestation import default_attestation_index
from wmxglycan.enumeration import Enumerator
from wmxglycan.fingerprint import default_fingerprint
from wmxglycan.store import Store

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

PAGES = (
    "New Prediction",
    "Prediction Result",
    "Experimental Validation",
    "Comparison",
    "Run History",
    "Model Monitor",
)
# The prototype's design tokens. This is wiring, not a redesign, so they must all survive.
PROTOTYPE_TOKENS = (
    "--ink:#153c55",
    "--blue:#2a6f97",
    "--sky:#eaf3f7",
    "--line:#d9e4ea",
    "--bg:#f6f9fb",
    "--good:#247a52",
    "--warn:#a86b16",
    "--bad:#a33b3b",
    "--muted:#667985",
)
PROTOTYPE_CLASSES = (
    ".top", ".brand", ".tag", ".shell", "aside", ".nav", ".main", ".hero", ".badge", ".grid",
    ".card", ".section-label", ".row", ".field", ".help", ".btns", ".btn", ".primary",
    ".secondary", ".outline", ".metric-grid", ".metric", ".candidate", ".cand-head", ".rank",
    ".score", ".bar", ".chips", ".chip", ".status", ".compare", ".big", ".arrow", ".hidden",
    ".footer-note", ".tabs", ".tab",
)


@pytest.fixture(scope="module")
def shared():
    return {
        "enumerator": Enumerator(),
        "index": default_attestation_index(),
        "fingerprint": default_fingerprint(),
    }


@pytest.fixture(scope="module")
def page() -> str:
    return DASHBOARD.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def script(page: str) -> str:
    return page.split("<script>", 1)[1].split("</script>", 1)[0]


@pytest.fixture
def client(shared):
    with TestClient(create_app(store=Store(), **shared)) as made:
        yield made


# --- it is served by the service, same origin -------------------------------------------------


def test_the_service_serves_the_dashboard_at_the_root(client):
    response = client.get("/")
    assert response.status_code == 200
    assert response.headers["content-type"].startswith("text/html")
    assert "WELLMATIX" in response.text


def test_the_dashboard_calls_relative_paths_so_it_is_same_origin(script):
    # No base URL, no CORS to configure, and no second process that could drift out of step with
    # the API it is a client of.
    assert 'const API = ""' in script
    for path in ("/v1/predictions", "/v1/runs", "/v1/models/current"):
        assert path in script


def test_the_dashboard_is_one_self_contained_file(page):
    # No build step, no bundler, no CDN: the same property the prototype had, kept.
    assert "<script src=" not in page
    assert "<link" not in page
    assert "http://" not in page.replace("http://127.0.0.1:8010/", "")
    assert "https://" not in page


def test_it_is_a_client_rather_than_a_demonstration(script):
    assert script.count("fetch(") >= 1
    # The prototype's fake behaviour is gone: no alert, and no hard-coded predicted CCS.
    assert "alert(" not in script
    assert "244.7" not in script


# --- the prototype's design survives ------------------------------------------------------------


@pytest.mark.parametrize("token", PROTOTYPE_TOKENS)
def test_every_design_token_from_the_prototype_survives(page, token):
    assert token in page


@pytest.mark.parametrize("selector", PROTOTYPE_CLASSES)
def test_every_layout_class_from_the_prototype_survives(page, selector):
    assert selector + "{" in page or selector + "," in page or selector + " " in page


@pytest.mark.parametrize("name", PAGES)
def test_all_six_pages_are_present_and_navigable(page, name):
    assert name in page
    assert f'>{name}<' in page


def test_the_sidebar_offers_all_six_pages(page):
    views = re.findall(r'data-view="([a-z]+)"', page)
    assert views == ["predict", "results", "validate", "comparison", "history", "monitor"]


# --- the contract between the client and the server ---------------------------------------------


def _model_fields() -> set[str]:
    """Every field name of every response model, derived from the models themselves."""
    found: set[str] = set()
    for module in (contracts, evidence_models):
        for value in vars(module).values():
            if isinstance(value, type) and issubclass(value, BaseModel) and value is not BaseModel:
                found |= set(value.model_fields)
    return found


def test_every_field_the_client_reads_is_a_real_field_of_a_real_response(script):
    """DERIVED ON BOTH SIDES, so neither can drift without this failing.

    snake_case is the discriminator: the API's fields are snake_case and the browser's own
    APIs are camelCase, so a `.something_like_this` in the script is a field read off a
    response. A typo produces `undefined` in a browser and fails nothing at all, which is why
    this is worth deriving rather than reviewing.
    """
    fields = _model_fields()
    assert len(fields) > 100, "the model walk found too little to be a real check"
    reads = {token for token in re.findall(r"\.([a-z][a-z0-9_]*)\b", script) if "_" in token}
    assert len(reads) > 50, "the script walk found too little to be a real check"
    unknown = sorted(reads - fields)
    assert unknown == [], f"the client reads fields no response has: {unknown}"


def test_the_contract_check_would_notice_a_typo(script):
    # The floor under the test above: if the walk could not see a bad field, it would pass on
    # anything. A deliberately wrong name must be caught.
    fields = _model_fields()
    broken = script.replace(".delta_percent", ".delta_pcent")
    reads = {token for token in re.findall(r"\.([a-z][a-z0-9_]*)\b", broken) if "_" in token}
    assert "delta_pcent" in sorted(reads - fields)


def test_every_endpoint_the_client_calls_exists_on_the_service(client, script):
    paths = {route.path for route in client.app.routes}
    assert "/v1/predictions" in paths
    assert "/v1/predictions/{prediction_id}" in paths
    assert "/v1/predictions/{prediction_id}/validation" in paths
    assert "/v1/predictions/{prediction_id}/comparison" in paths
    assert "/v1/runs" in paths
    assert "/v1/models/current" in paths
    assert "/" in paths
    # And the client calls all six, not a subset.
    for fragment in (
        '"/v1/predictions"',
        '"/v1/predictions/" + encodeURIComponent(id)',
        '"/validation"',
        '"/comparison"',
        '"/v1/runs?"',
        '"/v1/models/current"',
    ):
        assert fragment in script, fragment


# --- the candidate view must not manufacture an order -------------------------------------------


def test_the_client_never_numbers_a_candidate(script):
    """THE RULE THAT MATTERS MOST ON THIS PAGE.

    The prototype rendered `#1 Candidate GLY-ISO-001` with a score of 0.54. For
    Hex5HexNAc4Fuc1 that shape would number 167 candidates of which 152 sit in a class the
    platform's own 24 structure columns cannot separate - an arbitrary order presented as a
    ranking, which is the failure this project is built against.
    """
    assert "#1 Candidate" not in script
    assert "GLY-ISO" not in script
    # No index-based numbering of candidates or classes anywhere.
    assert "forEach((c, i)" not in script
    assert "index + 1" not in script
    assert ", i) =>" not in script
    # Bands ARE ordered and are the only thing given a position.
    assert "Band ' + band.rank" in script
    assert '" of " + total' in script


def test_the_client_says_out_loud_that_nothing_inside_a_band_is_ordered(script):
    assert "Nothing inside a band is ordered" in script
    assert "would present an arbitrary order as a ranking" in script
    assert "They are shown unordered because nothing separates them." in script


def test_a_tied_class_is_labelled_as_tied_and_says_the_platform_cannot_choose(script):
    assert "tied &middot; not separable" in script or "tied · not separable" in script
    assert "candidates sharing this position" in script
    assert "The platform cannot tell you which of them it is" in script


def test_the_client_reports_how_many_share_each_position(client, script):
    # The count, per class, from the response - not a rank.
    assert "c.members.length" in script
    assert "p.tied_candidates" in script
    assert "p.largest_indistinguishable_class" in script


def test_the_response_the_client_renders_really_does_carry_the_class_structure(client):
    """So the view above is not describing a shape the API does not produce."""
    body = client.post(
        "/v1/predictions",
        json={"composition": "Hex5HexNAc4Fuc1", "adduct": "[M+H]+", "charge": 1},
    ).json()
    assert body["candidates_total"] == 167
    assert body["classes_total"] == 61
    assert body["tied_candidates"] == 152
    assert body["largest_indistinguishable_class"] == 10
    assert len(body["bands"]) == 3
    # Every class says whether it is a tie, and every member list says its order is meaningless.
    assert all("is_a_tie" in one for one in body["classes"])
    assert all(one["members_order_is_not_meaningful"] for one in body["classes"])
    assert sum(len(one["members"]) for one in body["classes"]) == 167
    # And the multi-member classes really are the majority of the set.
    tied = [one for one in body["classes"] if one["is_a_tie"]]
    assert len(tied) == 46
    assert sum(len(one["members"]) for one in tied) == 152


def test_a_refused_set_is_rendered_as_not_ordered(client, script):
    body = client.post(
        "/v1/predictions", json={"composition": "Hex5HexNAc2", "adduct": "[M-H]-", "charge": -1}
    ).json()
    assert body["is_a_ranking"] is False
    assert "RANKING REFUSED" in body["refusal"]
    # The client has a branch for it that says so in capitals rather than showing a quiet list.
    assert "THIS SET IS NOT ORDERED." in script
    assert 'p.is_a_ranking ? "FROZEN · BANDED" : "FROZEN · NOT ORDERED"' in script


# --- no empty deltas: an unevaluable axis shows its reason ---------------------------------------


def test_an_unevaluable_delta_shows_the_reason_and_never_a_blank(script):
    """The owner's second rule. A dash in a delta column reads as "about zero"."""
    assert "No delta, and not because it is zero:" in script
    assert "m.unevaluable_because" in script
    # The em-dash placeholder is used for a missing COMPOSITION in the run table, never for a
    # delta. Assert the delta branch has no such fallback.
    delta_branch = script.split("THE REASON, NOT A BLANK", 1)[1].split("}).join", 1)[0]
    assert "unevaluable_because" in delta_branch
    assert '"—"' not in delta_branch and "'—'" not in delta_branch


def test_all_three_axes_are_rendered_with_their_own_state(script):
    assert "Axis 1 &mdash; against the AI prediction" in script
    assert "Axis 2 &mdash; against a reference the platform holds" in script
    assert "Axis 3 &mdash; your own measurements against each other" in script
    # Axis 1 is undefined, always, and carries the service's own reason.
    assert "c.against_prediction_because" in script
    assert '<span class="pill bad">undefined</span>' in script


def test_interval_coverage_is_never_rendered_as_yes_or_no(script):
    """The prototype showed `Interval coverage: YES`. There is no interval, so that was false."""
    assert "c.interval_coverage_because" in script
    assert '"YES"' not in script.split("Interval coverage", 1)[1].split("</div>", 1)[0]
    assert "not evaluable" in script.split("Interval coverage", 1)[1][:400]


def test_the_agreement_limit_travels_with_any_verdict_it_supports(script):
    assert "c.agreement_limit_percent" in script
    assert "c.agreement_limit_is_policy" in script
    assert "POLICY, chosen rather than measured" in script


def test_the_comparison_the_client_renders_really_carries_the_three_axes(client):
    created = client.post(
        "/v1/predictions", json={"composition": "Hex5HexNAc4Fuc1", "adduct": "[M+H]+", "charge": 1}
    ).json()
    pid = created["prediction_id"]
    client.post(
        f"/v1/predictions/{pid}/validation",
        json={
            "measurements": [
                {"ccs": 712.4, "uncertainty": 2.1, "uncertainty_type": "SD", "adduct": "[M+H]+",
                 "charge": 1, "ims_type": "TWIMS", "drift_gas": "N2", "source": "run A"},
                {"ccs": 715.9, "uncertainty": 2.6, "uncertainty_type": "SD", "adduct": "[M+H]+",
                 "charge": 1, "ims_type": "TWIMS", "drift_gas": "N2", "source": "run B"},
            ]
        },
    )
    body = client.get(f"/v1/predictions/{pid}/comparison").json()
    assert body["against_prediction"] == "no_predicted_value"
    assert body["against_prediction_because"]
    assert body["interval_coverage"] is None
    assert body["interval_coverage_because"]
    # Axis 2: a reason per measurement, never a null with nothing beside it.
    for one in body["per_measurement"]:
        assert one["delta_ccs"] is None
        assert one["unevaluable_because"], "a null delta must arrive with its reason"
    # Axis 3: the one that works.
    assert len(body["among_attached"]) == 1
    assert body["among_attached"][0]["spread"] == pytest.approx(3.5)


# --- nothing false is left over from the prototype ------------------------------------------------


def test_no_demo_value_survives(page):
    for stale in (
        "244.7",            # the prototype's predicted CCS
        "241.9",            # its interval
        "WMX-PRED-DEMO-001",
        "GLY-ISO-001",
        "DEMO VALUES",
        "glycan-ccs-v0.1-demo",
        "PROVISIONAL PASS",
        "Illustrative front-end behavior",
        "illustrative UI values",
    ):
        assert stale not in page, stale


def test_the_stale_ccs_corpus_figures_are_gone(page):
    """The prototype hard-coded 2,075 / 541 / 142 on the monitor page and flagged them itself as
    describing the CCS corpus rather than glycan coverage. They now come from the service.

    Searched in the BODY and as whole tokens. A first version of this test searched the whole
    file for "142" and failed on `max-width:1420px` in the prototype's own stylesheet - a test
    that reads a stylesheet as a scientific claim.
    """
    body = page.split("</style>", 1)[1]
    for stale in ("2,075", "2075", "541", "142"):
        assert not re.search(r"(?<![\d,.])" + re.escape(stale) + r"(?![\d,.])", body), stale
    assert "m.stamp.structures_indexed" in page
    assert "m.predictions_frozen" in page


def test_the_page_does_not_claim_a_predicted_ccs_anywhere(page):
    assert "Predicted CCS" not in page
    assert "No predicted CCS" in page or "no predicted CCS" in page
    assert "There is no glycan CCS model" in page or "no glycan CCS model" in page


def test_fields_the_service_does_not_accept_are_shown_as_not_sent(page):
    # The prototype offered drift gas, instrument, antibody class, site, retention time and
    # MS/MS on the prediction form. A field that looks submitted and is not is the same class of
    # false certainty as a fake rank, so they are named as not sent rather than quietly removed.
    assert "Not sent, and not silently dropped" in page
    for field in ("drift gas", "instrument / method", "antibody class", "retention time", "MS/MS evidence"):
        assert field in page, field


def test_the_client_states_what_the_service_actually_takes(page):
    assert "takes <b>composition, adduct and charge</b> and nothing else" in page


# --- provenance is on the page, not only in the response -----------------------------------------


def test_the_stamp_is_rendered_on_every_page_that_shows_a_response(script):
    assert script.count("stampCard(") >= 4
    for field in ("fingerprint_parameters", "fingerprint_corpus", "data_snapshot", "dataset_licence"):
        assert field in script, field


def test_the_dashboard_says_the_fingerprint_does_not_digest_source(script):
    # The one honest caveat about the fingerprint, carried to the reader rather than left in a
    # docstring they will never open.
    assert "not the source text" in script


def test_the_decision_rules_are_rendered_with_what_each_applies_to(script):
    assert "r.applies_to" in script
    assert "d.why_it_is_constant" in script
    assert "currently constant, and that is the answer rather than a defect" in script


def test_a_missing_dashboard_file_is_an_error_rather_than_an_empty_page(shared, tmp_path, monkeypatch):
    monkeypatch.setattr("wmxglycan.api.DASHBOARD", tmp_path / "gone.html")
    with TestClient(create_app(store=Store(), **shared)) as client:
        response = client.get("/")
    assert response.status_code == 500
    assert "ships as package data" in response.json()["detail"]


def test_the_dashboard_ships_as_package_data():
    # Or an installed wheel serves six working endpoints and no page.
    pyproject = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    assert "[tool.setuptools.package-data]" in pyproject
    assert 'wmxglycan = ["static/*.html"]' in pyproject
    assert DASHBOARD.is_file()
    assert DASHBOARD.parent.name == "static"


def test_the_page_is_still_about_the_size_of_the_prototype(page):
    # The prototype was 15 KB with zero fetch calls. Wiring it costs something; a redesign
    # would cost far more, and this is the cheapest available check that it stayed wiring.
    size = len(page.encode("utf-8"))
    assert 30_000 < size < 90_000, size


# --- the page's own JavaScript, executed against real responses ------------------------------------


def _render_with_node(prediction: dict, comparison: dict) -> str:
    """Run the SHIPPED page's rendering functions over real responses, via node.

    Not a re-implementation: the functions are cut out of dashboard.html, so what runs is what
    ships. Everything above this point tests the page as text; this tests it as code, which is
    the only way to find a rendering function that throws on a shape the API really returns.
    """
    import shutil
    import subprocess
    import tempfile

    node = shutil.which("node")
    if not node:
        pytest.skip("node is not on PATH, so the page's JavaScript cannot be executed here")
    script = DASHBOARD.read_text(encoding="utf-8").split("<script>", 1)[1].split("</script>", 1)[0]
    wanted = ("esc", "pct", "num", "kv", "metric", "stampCard", "candidatesCard", "shapeCard",
              "confidenceCard", "coverageCard", "decisionCard", "evidenceCard", "rulesCard")
    blocks = []
    for name in wanted:
        for pattern in (rf"\nfunction {name}\(.*?\n\}}\n", rf"\nconst {name} = .*?;\n"):
            found = re.search(pattern, script, re.S)
            if found:
                blocks.append(found.group(0))
                break
        else:
            raise AssertionError(f"could not extract {name} from the shipped page")
    assert len(blocks) == len(wanted)
    harness = "\n".join(blocks) + """
const fs = require("fs");
const p = JSON.parse(fs.readFileSync(process.argv[2], "utf-8"));
const c = JSON.parse(fs.readFileSync(process.argv[3], "utf-8"));
process.stdout.write(candidatesCard(p) + shapeCard(p) + confidenceCard(p) + coverageCard(p)
  + decisionCard(p) + evidenceCard(p) + rulesCard(p) + stampCard(p.stamp));
"""
    work = Path(tempfile.mkdtemp(prefix="wmx-render-"))
    try:
        (work / "r.js").write_text(harness, encoding="utf-8")
        (work / "p.json").write_text(json.dumps(prediction), encoding="utf-8")
        (work / "c.json").write_text(json.dumps(comparison), encoding="utf-8")
        result = subprocess.run(
            [node, str(work / "r.js"), str(work / "p.json"), str(work / "c.json")],
            capture_output=True, text=True, encoding="utf-8",
        )
        assert result.returncode == 0, result.stderr[-3000:]
        return result.stdout
    finally:
        shutil.rmtree(work, ignore_errors=True)


@pytest.fixture(scope="module")
def rendered(shared) -> str:
    with TestClient(create_app(store=Store(), **shared)) as client:
        prediction = client.post(
            "/v1/predictions",
            json={"composition": "Hex5HexNAc4Fuc1", "adduct": "[M+H]+", "charge": 1},
        ).json()
        pid = prediction["prediction_id"]
        client.post(
            f"/v1/predictions/{pid}/validation",
            json={"measurements": [
                {"ccs": 712.4, "uncertainty": 2.1, "uncertainty_type": "SD", "adduct": "[M+H]+",
                 "charge": 1, "ims_type": "TWIMS", "drift_gas": "N2", "source": "run A"},
                {"ccs": 715.9, "uncertainty": 2.6, "uncertainty_type": "SD", "adduct": "[M+H]+",
                 "charge": 1, "ims_type": "TWIMS", "drift_gas": "N2", "source": "run B"}]},
        )
        comparison = client.get(f"/v1/predictions/{pid}/comparison").json()
    return _render_with_node(prediction, comparison)


def test_the_rendered_page_gives_a_position_only_to_a_band(rendered):
    """THE RULE, CHECKED ON THE OUTPUT AND NOT ONLY ON THE SOURCE.

    Every `.rank` label the page emits for a 167-candidate set must be a BAND position. If a
    candidate were ever numbered it would appear here.
    """
    labels = sorted(set(re.findall(r'class="rank">([^<]*)', rendered)))
    assert labels == ["Band 1 of 3", "Band 2 of 3", "Band 3 of 3"], labels
    # And nothing anywhere in the output numbers a candidate.
    assert not re.search(r"#\s*\d+\s*(Candidate|candidate)", rendered)
    assert "Candidate GLY" not in rendered


def test_the_rendered_page_shows_every_tied_group_and_its_size(rendered):
    sizes = [int(n) for n in re.findall(r"<b>(\d+) candidates sharing this position</b>", rendered)]
    assert len(sizes) == 46, "all 46 multi-member classes must be rendered"
    assert sum(sizes) == 152, "the rendered tied candidates must total the 152 the API reports"
    assert max(sizes) == 10
    assert rendered.count("alone in its class") == 15, "and the 15 singletons say they are alone"
    assert rendered.count("tied") >= 46


def test_the_rendered_page_says_the_platform_cannot_separate_them(rendered):
    for phrase in (
        "cannot be ordered at all",
        "The platform cannot tell you which of them it is",
        "would present an arbitrary order as a ranking",
        "shown unordered because nothing separates them",
    ):
        assert phrase in rendered, phrase


def test_the_rendered_page_carries_no_predicted_cross_section(rendered):
    assert "Predicted CCS" not in rendered
    assert "244.7" not in rendered
    assert "no glycan CCS model" in rendered or "NOT A MOLECULE IN YOUR SAMPLE" in rendered


def test_the_rendered_page_carries_the_provenance(rendered):
    assert "Pipeline fingerprint" in rendered
    assert "Data snapshot" in rendered
    assert "MIT" in rendered and "Daniel Bojar" in rendered
    assert "not the source text" in rendered


def test_the_rendered_page_reports_every_reserved_mass(rendered):
    # So the number cannot be read as closed over the candidates shown.
    for label in (
        "On reference structures not enumerated",
        "On a structure nobody proposed",
        "On unattested classes",
    ):
        assert label in rendered, label
    assert "This is not a probability." in rendered
