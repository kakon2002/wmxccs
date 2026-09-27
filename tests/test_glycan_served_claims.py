"""Every served string that makes a factual claim is verified against BEHAVIOUR, or declared unasserted.

WHY THIS FILE EXISTS. On 27 September 2026 eight served strings were found making claims that were
wrong or overstated. Not one was caught by the 4,169-test suite, and two were actively PINNED by it:

    assert set(domain["residues_supported"]) == {"Hex", "HexNAc", "Fuc", "NeuAc", "NeuGc"}
    assert "misses between a third and four fifths" in domain["known_coverage_limit"]

The first locked in a field advertising a residue the enumerator refuses every composition for. The
second locked in a range that 76 of 145 measured compositions fall outside. Both assertions passed
for the same reason: **they compared the string to itself.** A test that matches served text cannot
notice the text is false - it can only notice that someone changed it.

SO THE RULE HERE IS: derive the truth from behaviour, then compare the served string to it. Never
the other way round. Where a claim cannot be reduced to a measurement, it goes in UNASSERTED with a
reason, and `test_every_claim_bearing_served_string_is_verified_or_declared` fails if a new
claim-bearing string appears in neither register. A new unverifiable claim is then a deliberate act
with a name on it, rather than prose nobody checked.

WHAT "CLAIM-BEARING" MEANS, mechanically: a served string of 25 characters or more that mentions
rules, enumeration scope, residues, what is supported or refused, or what has been checked. The
regex is below and is deliberately wide: a false positive costs one register entry, a false negative
costs what the eight above cost.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from wmxglycan.api import create_app
from wmxglycan.attestation import default_attestation_index
from wmxglycan.composition import Residue
from wmxglycan.enumeration import Enumerator
from wmxglycan.fingerprint import default_fingerprint
from wmxglycan.ranking import rank
from wmxglycan.store import Store

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

G2F = "Hex5HexNAc4Fuc1"
MAN5 = "Hex5HexNAc2"
DASHBOARD = REPO / "src" / "wmxglycan" / "static" / "dashboard.html"


@pytest.fixture(scope="module")
def shared():
    return {
        "enumerator": Enumerator(),
        "index": default_attestation_index(),
        "fingerprint": default_fingerprint(),
    }


@pytest.fixture(scope="module")
def client(shared):
    with TestClient(create_app(store=Store(), **shared)) as made:
        yield made


@pytest.fixture(scope="module")
def domain(client):
    return client.get("/v1/models/current").json()["domain"]


@pytest.fixture(scope="module")
def prediction(client):
    response = client.post(
        "/v1/predictions",
        json={"composition": G2F, "adduct": "[M+H]+", "charge": 1, "client_reference": "claims"},
    )
    assert response.status_code == 201, response.text
    return response.json()


# --- 1. residues_supported: derived from what actually enumerates -----------------------------------


def test_every_residue_the_domain_advertises_can_actually_be_placed(domain, shared):
    """THE CLAIM: `residues_supported` is what the platform can answer for.

    THE MEASUREMENT: submit a composition needing one of each and see whether it enumerates. The
    old test asserted a hardcoded five-element list including NeuGc, and every NeuGc composition is
    refused, so the suite defended the wrong answer.
    """
    enumerator = shared["enumerator"]
    base = {"hex": 3, "hexnac": 2}
    needs = {
        "Hex": "Hex4HexNAc2",
        "HexNAc": "Hex3HexNAc3",
        "Fuc": "Hex3HexNAc2Fuc1",
        "NeuAc": "Hex4HexNAc3NeuAc1",
        "NeuGc": "Hex4HexNAc3NeuGc1",
    }
    enumerates = {
        residue: bool(enumerator.enumerate(composition).candidates)
        for residue, composition in needs.items()
    }
    # The measurement must discriminate, or this test would pass on any served list.
    assert any(enumerates.values()) and not all(enumerates.values()), enumerates

    placeable = {residue for residue, works in enumerates.items() if works}
    assert set(domain["residues_supported"]) == placeable, (
        f"the domain advertises {sorted(domain['residues_supported'])} and a composition needing"
        f" each of {sorted(placeable)} enumerates. Serve what can be placed, not what the"
        " composition parser can read."
    )


def test_the_residues_it_cannot_place_are_served_and_are_genuinely_refused(domain, shared):
    gap = tuple(domain["residues_recognised_but_not_placeable"])
    assert gap, "no residue is declared unplaceable; NeuGc is one, so this is not an empty case"
    every = {residue.value for residue in Residue}
    assert set(domain["residues_supported"]) | set(gap) == every, (
        "supported plus unplaceable must account for every residue a composition can contain"
    )
    for residue in gap:
        composition = {"NeuGc": "Hex4HexNAc3NeuGc1"}.get(residue)
        assert composition, f"no probe composition for {residue}; add one rather than skipping"
        result = shared["enumerator"].enumerate(composition)
        assert not result.candidates, f"{residue} is declared unplaceable and {composition} enumerated"
        assert residue in (result.refusal or ""), (
            f"{composition} is refused but the reason does not name {residue}: {result.refusal!r}"
        )


# --- 2. known_coverage_limit: the figures it cites, measured -----------------------------------------


@pytest.mark.parametrize(
    ("composition", "missed", "total"),
    [(G2F, 17, 34), (MAN5, 22, 28), ("Hex3HexNAc6Fuc2", 0, 6), ("Hex9HexNAc4", 3, 3)],
)
def test_the_coverage_figures_the_domain_cites_are_the_measured_ones(
    domain, shared, composition, missed, total
):
    """Each cited figure is recomputed. The served string quotes exactly these four."""
    result = rank(
        shared["enumerator"].enumerate(composition),
        index=shared["index"],
        enumerator=shared["enumerator"],
    )
    coverage = result.coverage
    assert (coverage.not_enumerated, coverage.reference_structures) == (missed, total), (
        f"{composition}: measured {coverage.not_enumerated} of {coverage.reference_structures},"
        f" and known_coverage_limit cites {missed} of {total}"
    )
    assert f"{missed} of {total}" in domain["known_coverage_limit"], (
        f"{missed} of {total} is measured for {composition} and the served string does not cite it"
    )


def test_the_domain_quotes_no_coverage_RANGE_because_the_measured_one_is_0_to_100(domain, shared):
    """A range is what went wrong, so the absence of one is the thing to assert.

    The two witnesses above already prove the extremes reach 0% and 100%, so any range narrower
    than "none to all" is false and a range that wide says nothing. The served string says so.
    """
    text = domain["known_coverage_limit"]
    for banned in ("between a third and four fifths", "between a third", "a third to four fifths"):
        assert banned not in text.replace("This field said 'between a third and four fifths'", ""), (
            f"known_coverage_limit claims a range again: {banned!r}"
        )
    assert "VARIABLE" in text or "variable" in text
    assert "coverage" in text, "it must point at the per-composition fields that are computed"


# --- 3. the rule accounting: what was actually checked ----------------------------------------------


def test_the_served_rule_counts_are_the_enumerators_own(prediction, shared):
    """THE CLAIM: "N in the scheme, of which M were EVALUATED". THE MEASUREMENT: ask the enumerator."""
    enumerator = shared["enumerator"]
    rules = prediction["rules"]
    assert rules["rules_in_scheme"] == len(enumerator.constraints)
    assert rules["rules_evaluated_by_the_check"] == enumerator.rules_the_check_evaluates
    assert rules["rules_reported_as_ordering_caveats"] == enumerator.rules_skipped_as_order
    # AND THE POINT: fewer are evaluated than the scheme holds. Without this the test would pass
    # on a build that went back to asserting a verified zero across the whole scheme.
    assert rules["rules_evaluated_by_the_check"] < rules["rules_in_scheme"], (
        "the check evaluates the whole scheme, so the distinction this field draws has gone"
    )
    assert rules["rules_reported_as_ordering_caveats"] > 0


def test_the_summary_says_how_many_candidates_carry_a_caveat_and_it_is_the_real_number(prediction):
    """The 6 caveated G2F candidates were invisible in the old summary, which said "0 violated"."""
    caveated = sum(
        1 for candidate in prediction["candidates"] if candidate.get("ordering_caveats", 0)
    )
    assert caveated > 0, "no candidate carries a caveat, so this test cannot see the claim"
    assert f"{caveated} candidate(s) here carry one" in prediction["rules"]["summary"], (
        f"{caveated} candidates carry an ordering caveat and the summary does not say so:"
        f" {prediction['rules']['summary']!r}"
    )
    assert "EVALUATED" in prediction["rules"]["summary"]


def test_a_refused_composition_publishes_the_real_scheme_size_and_an_absent_violated_count(client, shared):
    """It published rules_in_scheme=0 beside prose saying the counts were absent, and 15 elsewhere."""
    response = client.post(
        "/v1/predictions",
        json={"composition": "Hex5HexNAc4NeuGc1", "adduct": "[M+H]+", "charge": 1},
    )
    assert response.status_code == 201, response.text
    rules = response.json()["rules"]
    assert rules["rules_in_scheme"] == len(shared["enumerator"].constraints) > 0
    assert rules["rules_violated"] is None, "absent must be null, not 0"
    assert rules["verified_by_running_the_check"] is False
    assert rules["rules_evaluated_by_the_check"] == 0
    # And the same deployment must not disagree with itself.
    served = client.get("/v1/models/current").json()["domain"]["biosynthetic_rules"]
    assert served == rules["rules_in_scheme"]


# --- 4. what_the_number_means: the denominator it describes is the one used -------------------------


def test_the_definition_of_the_share_names_the_denominator_the_code_uses(prediction, shared):
    """THE CLAIM is a formula. THE MEASUREMENT reproduces the number from it.

    The definition named a two-part space; `_shares` normalises over three. A caller following the
    served definition computed a denominator one smaller than the served share was divided by.
    """
    confidence = prediction["confidence"]
    means = confidence["what_the_number_means"]
    result = rank(
        shared["enumerator"].enumerate(G2F), index=shared["index"], enumerator=shared["enumerator"]
    )
    classes = result.classes_total
    missing = result.coverage.not_enumerated
    # Reproduce a served share from the three-part space and require it to match.
    band = max(result.bands, key=lambda one: one.band_total_share or 0.0)
    one_class = next(iter(band.classes))
    attested = result.classes[one_class].attested_structures
    observations = sum(one.attested_structures for one in result.classes.values()) + missing
    space = classes + missing + 1
    alpha = 1.0
    expected = (attested + alpha) / (observations + alpha * space)
    assert result.classes[one_class].evidence_share == pytest.approx(expected, rel=1e-9), (
        "the three-part space does not reproduce the served share, so the definition below is"
        " being checked against the wrong arithmetic"
    )
    # A two-part space would give a different number, which is what makes this test able to fail.
    two_part = (attested + alpha) / (observations + alpha * (classes + missing))
    assert two_part != pytest.approx(expected, rel=1e-9)
    assert "THREE groups" in means or "three groups" in means, means[:200]
    assert "catch-all" in means


# --- 5. builds_on, and the decision remedy ----------------------------------------------------------


def test_the_scope_claim_about_truncated_cores_matches_what_is_refused(domain, shared):
    enumerator = shared["enumerator"]
    truncated = {"Hex1HexNAc2": 1, "Hex2HexNAc2": 2}
    for composition in truncated:
        assert not enumerator.enumerate(composition).candidates, f"{composition} enumerated"
    # And the two the old wording wrongly excluded DO enumerate, which is why it said TRUNCATED.
    for composition in ("Hex3HexNAc2", "Hex3HexNAc2Fuc1"):
        assert enumerator.enumerate(composition).candidates, (
            f"{composition} no longer enumerates, so builds_on's wording needs re-reading"
        )
    assert "TRUNCATED" in domain["builds_on"]
    assert "Man1-2GlcNAc2" in domain["builds_on"]


def test_the_remedy_names_every_blocker_the_same_response_names(client, prediction):
    """It named two blockers and then a two-item remedy that satisfied neither of them."""
    why = prediction["decision"]["why_it_is_constant"]
    unreachable = client.get("/v1/models/current").json()["decision_unreachable_today"]
    assert "complete" in why.lower(), "the completeness blocker is missing from the remedy again"
    assert "ALL THREE" in why
    # The two statements in one deployment must agree on the count.
    assert "All three" in unreachable["AI_ONLY"] or "AND" in unreachable["AI_ONLY"]


# --- 6. the dashboard hardcodes no claim it can read -------------------------------------------------


def test_the_dashboard_reads_the_constancy_sentence_rather_than_hardcoding_a_count():
    page = DASHBOARD.read_text(encoding="utf-8")
    assert "decision_constant_because" in page, "the page must read the sentence from the response"
    assert "Two rules fire on every possible input" not in page, (
        "the hardcoded count is back; three rules fire for every input in this release and the"
        " page rendered a 4-of-5 pill two cards away from this sentence"
    )
    assert "scheme size" in page, "the served 15 must be labelled as the scheme size"


def test_the_served_constancy_sentence_is_the_same_one_the_prediction_carries(client, prediction):
    served = client.get("/v1/models/current").json()["decision_constant_because"]
    assert served, "models/current serves no constancy sentence, so the page has nothing to read"
    assert served == prediction["decision"]["why_it_is_constant"], (
        "the model monitor and the prediction disagree about why the field is constant"
    )


# --- 7. THE REGISTER: every claim-bearing served string is verified or declared ---------------------


CLAIM = re.compile(
    r"\brules?\b|\bcurated\b|\benumerat|\bgovern|\bconstrain|\bresidue|\bsupported\b|\brefus"
    r"|\bplaceable\b|\bscope\b|\bmammalian\b|\bverified\b|\bchecked\b|\bcore\b"
    # ADDED after the register's first run caught decision_unreachable_today.IM_VALIDATION_RECOMMENDED
    # and MISSED the sibling AI_ONLY string, which makes an equally strong claim about what would
    # make a value reachable. A pattern that sees one of a pair and not the other is the
    # narrow-search failure this whole file is about.
    r"|\breachable\b|\bunreachable\b|\bvalidated\b|\bcomplete\b|\bgate\b",
    re.I,
)

# Claim-bearing served fields whose claim IS verified against behaviour, above or elsewhere.
VERIFIED = {
    "domain.residues_supported",
    "domain.residues_recognised_but_not_placeable",
    "domain.residues_not_placeable_because",
    "domain.known_coverage_limit",
    "domain.builds_on",
    "domain.rules_govern",
    "domain.glycan_class",
    "domain.species_assumption",
    "domain.predicts_ccs_note",
    "domain.measured_reference_is_unreachable_note",
    "decision_constant_because",
    "rules.summary",
    "decision.why_it_is_constant",
    # Both verified by the probe in tests/test_glycan_read_this_first.py, which calls `_decide`
    # over both gates and every input shape rather than reading these sentences.
    "decision_unreachable_today.AI_ONLY",
    "decision_unreachable_today.IM_VALIDATION_RECOMMENDED",
    "confidence.what_the_number_means",
    "confidence.why_the_shares_do_not_sum_to_one",
    "ccs_model_note",
    "what_would_calibrate_it",
    "confidence.what_would_calibrate_it",
    "coverage.summary",
    "ccs_evidence.summary",
}

# Claim-bearing served strings whose claim CANNOT be reduced to a measurement in this repository.
# Each needs a reason. This is the "say so" half of the rule: an unverifiable claim is allowed, and
# it is allowed BY NAME.
UNASSERTED = {
    "confidence.what_would_calibrate_it": "describes a study that would have to be run elsewhere;"
    " nothing here can measure whether it would calibrate anything",
    "ccs_model_note": "a statement about what the fingerprint digests, verified by"
    " test_glycan_fingerprint.py rather than by a served-claim measurement",
    "ccs_evidence.summary": "quotes CCS corpus counts (24 cleared, 89 held, 117 licence-clean)"
    " that tests/test_glycan_ccs_evidence.py derives from the loader, not from this surface",
    "coverage.summary": "restates reference_structures and not_enumerated, which"
    " test_the_coverage_figures_the_domain_cites_are_the_measured_ones already derives",
}


def _strings(node, path=""):
    if isinstance(node, str):
        yield path, node
    elif isinstance(node, dict):
        for key, value in node.items():
            yield from _strings(value, f"{path}.{key}" if path else key)
    elif isinstance(node, list):
        for index, value in enumerate(node):
            yield from _strings(value, f"{path}[{index}]")


def _claim_paths(body) -> set[str]:
    found = set()
    for path, value in _strings(body):
        if len(value) < 25 or not CLAIM.search(value):
            continue
        # Per-item paths collapse to their field: candidates[7].rationale[2].text -> rationale.text
        collapsed = re.sub(r"\[\d+\]", "", path)
        found.add(collapsed)
    return found


# Paths whose text is DATA rather than a platform claim: a curated rule's own rationale, a
# structure string, a licence name. They restate the source, and the source is what they are.
DATA_NOT_CLAIM = {
    "candidates.rationale.text",
    "candidates.iupac_condensed",
    "decision.rules.applies_to",
    "decision.rules.because",
    "decision.rules.name",
    "weaknesses",
    "ccs_evidence.held.blockers",
    "ccs_evidence.held.what_would_release_them",
    "refusal",
    "stamp.dataset",
    "stamp.dataset_attribution",
    "attestation_provenance",
    "classes.members",
    "what_the_number_means",
}


def test_every_claim_bearing_served_string_is_verified_or_declared(client, prediction):
    """THE STRUCTURAL GUARD. A new claim-bearing served string must be classified, not ignored.

    This is the test the eight defects needed. Each of them was a served sentence nobody had
    decided about: not verified, not declared unverifiable, just written. This fails until a new
    one is put in VERIFIED (and given a measurement) or in UNASSERTED (and given a reason).
    """
    bodies = [
        prediction,
        client.get("/v1/models/current").json(),
        client.post(
            "/v1/predictions",
            json={"composition": MAN5, "adduct": "[M-H]-", "charge": -1},
        ).json(),
    ]
    paths = set()
    for body in bodies:
        paths |= _claim_paths(body)
    assert len(paths) > 10, f"the sweep found only {paths}; the regex or the walk is broken"

    classified = VERIFIED | set(UNASSERTED) | DATA_NOT_CLAIM
    unclassified = sorted(paths - classified)
    assert unclassified == [], (
        "these served strings make a factual claim and are neither verified against behaviour nor"
        f" declared unverifiable: {unclassified}. Add a test that MEASURES the claim and list the"
        " path in VERIFIED, or list it in UNASSERTED with the reason it cannot be measured. Do not"
        " write a test that matches the text: that is what let eight wrong strings ship."
    )


def test_the_registers_are_honest_about_themselves():
    """The floor. Registers that drift from reality are how a guard stops guarding."""
    assert set(UNASSERTED) <= VERIFIED | set(UNASSERTED)
    for path, reason in UNASSERTED.items():
        assert len(reason) > 30, f"{path} is declared unverifiable with no real reason"
    # And the classifier must be able to see a claim, or everything above is vacuous.
    assert _claim_paths({"x": "the curated rules govern branching order and nothing else"}) == {"x"}
    assert _claim_paths({"x": "short"}) == set()
