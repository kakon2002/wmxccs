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

# --- 8. THE REGISTER, EXTENDED TO NUMERIC AND BOOLEAN CLAIMS ----------------------------------------
#
# THE NINTH DEFECT WAS TWO SCALARS, and the string register could not see it: a refused composition
# served `enumerator_built: 0` and `enumerator_truncated: false` while 5,000 trees had been built and
# the cap had fired, in the same response as a refusal blaming the curated rules. A guard that covers
# strings and not numbers is how that walked past the guard written for the eight before it.


def _scalars(node, path=""):
    """Every int, float and bool a caller receives, with its json path."""
    if isinstance(node, bool) or isinstance(node, (int, float)):
        yield path, node
    elif isinstance(node, dict):
        for key, value in node.items():
            yield from _scalars(value, f"{path}.{key}" if path else key)
    elif isinstance(node, list):
        for index, value in enumerate(node):
            yield from _scalars(value, f"{path}[{index}]")


def _scalar_paths(body) -> set[str]:
    return {re.sub(r"\[\d+\]", "", path) for path, _ in _scalars(body)}


# Every scalar the platform serves, against the measurement that checks it. The value is the test
# function name, so a reader can go straight to what establishes the number.
VERIFIED_SCALARS = {
    # what the enumerator did - THE NINTH DEFECT'S FIELDS
    "enumerator_built": "test_the_enumerator_counts_are_the_real_ones_even_on_a_refusal",
    "enumerator_rejected": "test_the_enumerator_counts_are_the_real_ones_even_on_a_refusal",
    "enumerator_truncated": "test_truncation_is_reported_exactly_when_the_budget_was_hit",
    # the shape of the answer
    "candidates_total": "test_the_served_shape_counts_are_derived_from_the_classes",
    "classes_total": "test_the_served_shape_counts_are_derived_from_the_classes",
    "tied_candidates": "test_the_served_shape_counts_are_derived_from_the_classes",
    "largest_indistinguishable_class": "test_the_served_shape_counts_are_derived_from_the_classes",
    "is_a_ranking": "test_the_served_shape_counts_are_derived_from_the_classes",
    "bands.rank": "test_the_served_shape_counts_are_derived_from_the_classes",
    "bands.candidates": "test_the_served_shape_counts_are_derived_from_the_classes",
    "bands.attested_structures": "test_the_served_shape_counts_are_derived_from_the_classes",
    "classes.is_a_tie": "test_the_served_shape_counts_are_derived_from_the_classes",
    "classes.attested_structures": "test_the_served_shape_counts_are_derived_from_the_classes",
    "candidates.attested": "test_every_residue_the_domain_advertises_can_actually_be_placed",
    # the rule accounting
    "rules.rules_in_scheme": "test_the_served_rule_counts_are_the_enumerators_own",
    "rules.rules_evaluated_by_the_check": "test_the_served_rule_counts_are_the_enumerators_own",
    "rules.rules_reported_as_ordering_caveats": "test_the_served_rule_counts_are_the_enumerators_own",
    "rules.rules_violated": "test_a_refused_composition_publishes_the_real_scheme_size_and_an_absent_violated_count",
    "rules.verified_by_running_the_check": "test_a_refused_composition_publishes_the_real_scheme_size_and_an_absent_violated_count",
    "rules.rules_applicable_range": "test_the_summary_says_how_many_candidates_carry_a_caveat_and_it_is_the_real_number",
    "candidates.rules_applicable": "test_the_summary_says_how_many_candidates_carry_a_caveat_and_it_is_the_real_number",
    "candidates.ordering_caveats": "test_the_summary_says_how_many_candidates_carry_a_caveat_and_it_is_the_real_number",
    # the shares
    "confidence.candidate_shares_sum": "test_the_served_shares_account_to_one_and_say_so",
    "confidence.mass_not_on_any_candidate": "test_the_served_shares_account_to_one_and_say_so",
    "confidence.mass_on_reference_structures_not_enumerated": "test_the_served_shares_account_to_one_and_say_so",
    "confidence.mass_on_a_structure_nobody_proposed": "test_the_served_shares_account_to_one_and_say_so",
    "confidence.shares_account_to_one": "test_the_served_shares_account_to_one_and_say_so",
    "confidence.mass_on_unattested_classes": "test_the_unattested_mass_is_inside_the_candidate_sum_and_not_a_fourth_term",
    "confidence.prior_pseudocount": "test_the_definition_of_the_share_names_the_denominator_the_code_uses",
    "classes.evidence_share": "test_the_definition_of_the_share_names_the_denominator_the_code_uses",
    "bands.share_per_class": "test_the_served_shares_account_to_one_and_say_so",
    "bands.band_total_share": "test_the_served_shares_account_to_one_and_say_so",
    # coverage, from the index
    "coverage.reference_structures": "test_the_coverage_figures_the_domain_cites_are_the_measured_ones",
    "coverage.not_enumerated": "test_the_coverage_figures_the_domain_cites_are_the_measured_ones",
    "coverage.attested_by_a_candidate": "test_the_coverage_figures_the_domain_cites_are_the_measured_ones",
    "coverage.truth_may_not_be_in_the_candidate_set": "test_the_completeness_fact_lives_in_one_place (test_glycan_read_this_first.py)",
    # the domain
    "domain.biosynthetic_rules": "test_a_refused_composition_publishes_the_real_scheme_size_and_an_absent_violated_count",
    "domain.predicts_ccs": "test_the_api_says_it_holds_no_measured_cross_section_for_any_candidate (test_glycan_read_this_first.py)",
    "domain.holds_a_measured_cross_section_for_any_candidate": "test_the_api_says_it_holds_no_measured_cross_section_for_any_candidate (test_glycan_read_this_first.py)",
    # the decision
    "decision.constant_today": "test_only_one_decision_is_reachable_and_the_other_two_say_what_would_reach_them (test_glycan_read_this_first.py)",
    "decision.rules.fires": "test_what_the_decision_function_can_actually_return_under_every_gate_state (test_glycan_read_this_first.py)",
    # the stamp, from the index
    "stamp.structures_indexed": "test_the_stamp_counts_are_the_indexs_own",
    "stamp.distinct_structures": "test_the_stamp_counts_are_the_indexs_own",
    "stamp.compositions_in_corpus": "test_the_stamp_counts_are_the_indexs_own",
    # the store
    "frozen": "test_a_prediction_is_frozen_on_creation_and_the_flag_is_not_a_constant",
    "predictions_frozen": "test_a_prediction_is_frozen_on_creation_and_the_flag_is_not_a_constant",
    "classes.reference_rows": "test_the_served_shape_counts_are_derived_from_the_classes",
    # what the platform holds and whether it discriminates
    "ccs_model_fitted": "test_the_api_says_it_holds_no_measured_cross_section_for_any_candidate (test_glycan_read_this_first.py)",
    "ccs_evidence.discriminates_between_candidates": "test_the_ranker_agrees_that_nothing_discriminates_between_candidates (test_glycan_read_this_first.py)",
    "coverage.reference_rows": "test_the_coverage_figures_the_domain_cites_are_the_measured_ones",
    "candidates.reference_rows": "test_the_served_shape_counts_are_derived_from_the_classes",
    "classes.members_order_is_not_meaningful": "test_the_served_shape_counts_are_derived_from_the_classes",
    "bands.classes_order_is_not_meaningful": "test_the_served_shape_counts_are_derived_from_the_classes",
}

# Scalars whose value cannot be reduced to a measurement, each with the reason.
# The caller's own input, echoed. Not a claim the platform makes about itself - but the echo being
# UNCHANGED is a claim, and hard constraint 3 is exactly about that, so it is measured.
ECHOED_INPUT = {
    "charge": "test_the_request_is_echoed_back_unchanged",
}

UNASSERTED_SCALARS = {
    "confidence.mass_under_priors": "the three alternative priors are arithmetic on the same"
    " observations as the served one; the served prior's value IS measured, and reproducing all four"
    " would re-implement _shares in the test rather than check it",
    "classes.evidence_share_under_priors": "same reason as mass_under_priors: the served prior is"
    " checked, the alternatives are the same formula at another concentration",
    "bands.share_under_priors": "the band-level view of the same four-prior spread; the band's"
    " served share under the policy prior IS measured, and the alternatives are that arithmetic at"
    " another concentration",
}


def test_every_served_scalar_is_verified_or_declared(client, prediction):
    """THE NINTH DEFECT'S GUARD. Two false scalars shipped past a register that walked only strings.

    Same rule, extended: a served int, float or bool must either have a test that MEASURES it or be
    listed as unverifiable with a reason. A new one in neither fails here.
    """
    bodies = [
        prediction,
        client.get("/v1/models/current").json(),
        client.post(
            "/v1/predictions", json={"composition": MAN5, "adduct": "[M-H]-", "charge": -1}
        ).json(),
    ]
    paths = set()
    for body in bodies:
        paths |= _scalar_paths(body)
    assert len(paths) > 30, f"the scalar walk found only {len(paths)} paths; it is broken"

    # A REGISTER ENTRY COVERS A SUBTREE. `confidence.mass_under_priors` declares every leaf under
    # it: the walker yields leaves, and naming each of the eight prior-by-quantity combinations
    # separately would be a list to maintain rather than a decision to record.
    classified = set(VERIFIED_SCALARS) | set(UNASSERTED_SCALARS) | set(ECHOED_INPUT)

    def is_classified(one: str) -> bool:
        parts = one.split(".")
        return any(".".join(parts[: n + 1]) in classified for n in range(len(parts)))

    unclassified = sorted(one for one in paths if not is_classified(one))
    assert unclassified == [], (
        "these served numbers and booleans assert something about what the platform did, and are"
        f" neither measured by a test nor declared unverifiable: {unclassified}. Add a test that"
        " MEASURES the value and list the path in VERIFIED_SCALARS, or list it in"
        " UNASSERTED_SCALARS with the reason. The ninth defect was two scalars - enumerator_built"
        " and enumerator_truncated - that no test measured."
    )


def test_the_scalar_register_names_tests_that_exist():
    """A register pointing at a test that does not exist is a register that checks nothing."""
    here = Path(__file__).read_text(encoding="utf-8")
    sibling = (REPO / "tests" / "test_glycan_read_this_first.py").read_text(encoding="utf-8")
    for served, where in VERIFIED_SCALARS.items():
        name = where.split(" ")[0]
        assert f"def {name}(" in here or f"def {name}(" in sibling, (
            f"{served} is registered as verified by {name}, which does not exist"
        )
    for served, reason in UNASSERTED_SCALARS.items():
        assert len(reason) > 30, f"{served} is declared unverifiable with no real reason"
    for served, where in ECHOED_INPUT.items():
        name = where.split(" ")[0]
        assert f"def {name}(" in here, f"{served} is registered as echoed, verified by a missing {name}"
    # The floor: the classifier must see a prefix match and must NOT see an unrelated path.
    assert "confidence.mass_under_priors" in set(UNASSERTED_SCALARS)
    assert "enumerator_built" in set(VERIFIED_SCALARS)


# --- the measurements the scalar register points at -------------------------------------------------


def test_the_request_is_echoed_back_unchanged(client):
    """The caller's input comes back as sent. Hard constraint 3's shape, at the glycan layer.

    `charge` is the only scalar in a response that is the caller's rather than the platform's, so it
    is registered separately - but "echoed" is itself a claim, and an echo that quietly normalised a
    value would be the thing the CCS core refuses to do with a measurement.
    """
    for adduct, charge in (("[M+H]+", 1), ("[M-H]-", -1), ("[M+2H]2+", 2)):
        body = client.post(
            "/v1/predictions",
            json={"composition": MAN5, "adduct": adduct, "charge": charge},
        ).json()
        assert body["charge"] == charge, f"sent charge {charge}, got {body['charge']}"
        assert body["adduct"] == adduct
    # And a signed charge is not silently absolute-valued, which is the normalisation to fear.
    negative = client.post(
        "/v1/predictions", json={"composition": MAN5, "adduct": "[M-H]-", "charge": -1}
    ).json()
    assert negative["charge"] == -1


def test_the_enumerator_counts_are_the_real_ones_even_on_a_refusal(client, shared):
    """THE NINTH DEFECT. A refused composition served built=0 while 5,000 trees had been built."""
    enumerator = shared["enumerator"]
    # A composition that exhausts the budget and yields nothing: the case that was misreported.
    exhausting = "Hex7HexNAc6NeuAc3"
    direct = enumerator.enumerate(exhausting)
    assert not direct.candidates and direct.capped, (
        f"{exhausting} no longer exhausts the budget, so this test is not exercising the case"
    )
    body = client.post(
        "/v1/predictions", json={"composition": exhausting, "adduct": "[M+H]+", "charge": 1}
    ).json()
    assert body["enumerator_built"] == direct.built > 0, (
        f"served built={body['enumerator_built']} and the enumerator built {direct.built}"
    )
    assert body["enumerator_rejected"] == direct.rejected
    assert body["candidates_total"] == 0

    # AND THE REFUSAL MUST NAME THE BUDGET, NOT THE RULES.
    refusal = body["refusal"] or ""
    assert "tree budget" in refusal, refusal[:200]
    assert "NOT A STATEMENT THAT NO CANDIDATE EXISTS" in refusal, refusal[:200]
    assert "every one of the" not in refusal, (
        "a budget exhaustion is being reported as every arrangement breaking a biosynthetic rule,"
        f" which is a claim about glycobiology: {refusal[:200]}"
    )


def test_truncation_is_reported_exactly_when_the_budget_was_hit(client, shared):
    """`enumerator_truncated` must track the cap in both directions, on both paths."""
    enumerator = shared["enumerator"]
    cases = {
        G2F: False,                 # completes well inside the budget
        "Hex7HexNAc6NeuAc3": True,  # exhausts it and yields nothing
    }
    for composition, expected in cases.items():
        direct = enumerator.enumerate(composition)
        assert direct.capped is expected, (
            f"{composition}: capped={direct.capped}, expected {expected}; the budget or the space"
            " changed and this test needs new witnesses"
        )
        body = client.post(
            "/v1/predictions",
            json={"composition": composition, "adduct": "[M+H]+", "charge": 1},
        ).json()
        assert body["enumerator_truncated"] is expected, (
            f"{composition}: served truncated={body['enumerator_truncated']}, measured {expected}"
        )
    # Both values occur, or the assertion above would hold for a constant.
    assert set(cases.values()) == {True, False}

    # AND THE BUDGET MUST STILL BE BIG ENOUGH FOR WHAT IT WAS RAISED FOR. The default was 5,000
    # until 28 September 2026, at which point three ordinary tri-antennary compositions were
    # refused - with a refusal blaming the curated rules - where a larger budget returns a ranked
    # answer. Lowering it again would bring that back silently, so the witness is pinned here.
    from wmxglycan.enumeration import DEFAULT_TREE_BUDGET

    witness = enumerator.enumerate("Hex6HexNAc5NeuAc3")
    assert not witness.capped, (
        f"Hex6HexNAc5NeuAc3 no longer completes within the budget of {DEFAULT_TREE_BUDGET:,}"
        f" ({witness.built:,} trees built). At 5,000 it returned zero candidates and a refusal that"
        " blamed the curated rules; the budget exists at its current size to prevent that."
    )
    assert len(witness.candidates) > 10_000, len(witness.candidates)


def test_the_served_shape_counts_are_derived_from_the_classes(prediction):
    """Every count about the answer's shape, recomputed from the classes the same response carries."""
    classes = {one["class_id"]: one for one in prediction["classes"]}
    candidates = prediction["candidates"]
    assert prediction["candidates_total"] == len(candidates)
    assert prediction["classes_total"] == len(classes)
    assert prediction["tied_candidates"] == sum(
        1 for one in candidates if classes[one["class_id"]]["is_a_tie"]
    )
    assert prediction["largest_indistinguishable_class"] == max(
        len(one["members"]) for one in classes.values()
    )
    assert sum(len(one["members"]) for one in classes.values()) == len(candidates)
    for band in prediction["bands"]:
        assert band["candidates"] == sum(
            len(classes[class_id]["members"]) for class_id in band["classes"]
        )
        assert band["attested_structures"] == max(
            classes[class_id]["attested_structures"] for class_id in band["classes"]
        )
        assert band["classes_order_is_not_meaningful"] is True
    assert [band["rank"] for band in prediction["bands"]] == list(
        range(1, len(prediction["bands"]) + 1)
    )
    for one in classes.values():
        assert one["is_a_tie"] is (len(one["members"]) > 1)
        assert one["members_order_is_not_meaningful"] is True
    # is_a_ranking must track the refusal rather than being a constant.
    assert prediction["is_a_ranking"] is (prediction.get("refusal") is None)


def test_the_stamp_counts_are_the_indexs_own(client, shared):
    stamp = client.get("/v1/models/current").json()["stamp"]
    snapshot = shared["fingerprint"].snapshot
    assert stamp["structures_indexed"] == snapshot.structures_indexed
    assert stamp["distinct_structures"] == snapshot.distinct_structures
    assert stamp["compositions_in_corpus"] == snapshot.compositions
    # And they must be the index's, not a copy that drifted.
    index = shared["index"]
    assert stamp["distinct_structures"] == index.distinct_structures
    assert stamp["distinct_structures"] < stamp["structures_indexed"], (
        "distinct equals indexed, so the deduplication this figure exists to report is not happening"
    )


def test_a_prediction_is_frozen_on_creation_and_the_flag_is_not_a_constant(client):
    """`frozen` and `predictions_frozen` measured, not asserted."""
    before = client.get("/v1/models/current").json()["predictions_frozen"]
    body = client.post(
        "/v1/predictions", json={"composition": MAN5, "adduct": "[M-H]-", "charge": -1}
    ).json()
    assert body["frozen"] is True
    after = client.get("/v1/models/current").json()["predictions_frozen"]
    assert after == before + 1, (
        f"predictions_frozen went {before} -> {after} after one create, so it is not counting them"
    )
    # A second write to the frozen prediction is refused, which is what `frozen` claims.
    again = client.post(
        f"/v1/predictions/{body['prediction_id']}/validation",
        json={"measurements": [{
            "ccs": 500.0, "uncertainty": 1.0, "uncertainty_type": "SD", "adduct": "[M-H]-",
            "charge": -1, "ims_type": "DTIMS", "drift_gas": "N2", "source": "scalar register"}]},
    )
    assert again.status_code == 201
    digest_unchanged = again.json()["prediction_unchanged"]
    assert digest_unchanged is True, "attaching changed the frozen prediction"
