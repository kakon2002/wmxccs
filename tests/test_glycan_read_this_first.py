"""The four claims of 27 September 2026 that were moved out of documents and into the product.

The owner's instruction was that three of the seven Tuesday misreadings should be met BEFORE a
reader can go wrong rather than after, and that the fifth should be stated without hedging. A
document alone cannot do that, so each of them landed in code:

  item 1  the two maturities, in the first lines of whatever a reader opens first
  item 3  the shares not summing to 1, next to the number in the response and on the dashboard
  item 4  the two unreachable decisions, at the enum, with what would reach them
  item 5  the measured-reference path stated as UNREACHABLE, and counted rather than claimed

WHAT THESE TESTS DO AND DO NOT ASSERT. Asserting that a docstring contains a sentence is a weak
test and it is the kind this repository has been burned by, so the checks here are mostly on
BEHAVIOUR: the numbers the API serves, the identity they satisfy, and the fact that the published
reachability MOVES when the gate it is derived from is patched. Where a test does read prose - the
banner in README, the seven sections of READ_THIS_FIRST - it is checking that a specific
navigational promise made elsewhere is kept, and it says so.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

import wmxglycan.ranking as ranking
from wmxglycan.api import create_app
from wmxglycan.attestation import default_attestation_index
from wmxglycan.ccs_evidence import CCSEvidenceState
from wmxglycan.enumeration import Enumerator
from wmxglycan.fingerprint import default_fingerprint
from wmxglycan.ranking import Decision, decision_reachability, rank
from wmxglycan.store import Store

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

READ_THIS_FIRST = REPO / "docs" / "READ_THIS_FIRST.md"
G2F = "Hex5HexNAc4Fuc1"
MAN5 = "Hex5HexNAc2"
# Four compositions with deliberately different shapes: a small set, two large ones, and one the
# enumerator refuses to order. The accounting identity has to hold for all of them, including the
# refusal - a summary that crashed on the refused set is how this was found once already.
EVERY_SHAPE = ("Hex3HexNAc4Fuc1", "Hex4HexNAc4Fuc1", G2F, MAN5)


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


def created(client, composition, adduct, charge):
    response = client.post(
        "/v1/predictions",
        json={"composition": composition, "adduct": adduct, "charge": charge},
    )
    assert response.status_code == 201, response.text
    return response.json()


# --- item 1: the maturity wall, in the first lines --------------------------------------------------


@pytest.mark.parametrize(
    ("relative", "lines"),
    [("README.md", 12), ("docs/HANDOVER.md", 14), ("docs/READ_THIS_FIRST.md", 10)],
)
def test_the_opening_lines_say_two_maturities_and_no_cross_section_prediction(relative, lines):
    """A reader meets this before they can carry the CCS core's maturity across the wall.

    It reads only the OPENING of each file, because a statement in section 9 does not prevent a
    misreading formed in section 1. The line budget is what makes this test able to fail: moving
    the banner further down breaks it.
    """
    opening = "\n".join((REPO / relative).read_text(encoding="utf-8").splitlines()[:lines])
    assert "wmxglycan" in opening and "wmxccs" in opening, opening
    assert re.search(r"two\s+(?:things|packages)", opening, re.I), opening
    assert re.search(r"maturit", opening, re.I), opening
    assert re.search(r"no\s+(?:fitted\s+)?(?:glycan\s+)?(?:CCS\s+)?model", opening, re.I), opening
    assert re.search(r"not\s+predict|does not predict|no glycan CCS model", opening, re.I), opening


def test_readme_and_handover_both_point_a_new_reader_at_read_this_first():
    for relative in ("README.md", "docs/HANDOVER.md"):
        text = (REPO / relative).read_text(encoding="utf-8")
        assert "READ_THIS_FIRST.md" in text, f"{relative} does not point at the orientation document"


def test_read_this_first_covers_all_seven_and_names_the_files_that_enforce_them():
    text = READ_THIS_FIRST.read_text(encoding="utf-8")
    headings = re.findall(r"^## (\d)\. ", text, re.M)
    assert headings == ["1", "2", "3", "4", "5", "6", "7"], headings
    # Each item points at something checkable rather than only asserting. These are the files the
    # document tells a reader to prefer over itself, so a rename must not leave the advice dangling.
    for named in (
        "tests/test_glycan_boundary.py",
        "src/wmxglycan/ranking.py",
        "src/wmxglycan/ccs_evidence.py",
        "docs/GLYCAN_PORT.md",
        "docs/GLYCAN_LIMITATIONS.md",
        "LIMITATIONS.md",
    ):
        assert named in text, f"READ_THIS_FIRST names no {named}"
        assert (REPO / named).exists(), f"READ_THIS_FIRST points at {named}, which does not exist"


# --- item 3: the shares, next to the number ---------------------------------------------------------


@pytest.mark.parametrize("composition", EVERY_SHAPE)
def test_the_served_shares_account_to_one_and_say_so(client, composition):
    """The identity, checked against the real numbers for every shape of set the platform makes."""
    adduct, charge = ("[M+H]+", 1) if composition != MAN5 else ("[M-H]-", -1)
    confidence = created(client, composition, adduct, charge)["confidence"]

    total = confidence["candidate_shares_sum"]
    not_enumerated = confidence["mass_on_reference_structures_not_enumerated"]
    not_proposed = confidence["mass_on_a_structure_nobody_proposed"]
    assert None not in (total, not_enumerated, not_proposed), confidence

    assert total + not_enumerated + not_proposed == pytest.approx(1.0, abs=1e-9)
    assert confidence["shares_account_to_one"] is True
    assert confidence["mass_not_on_any_candidate"] == pytest.approx(1.0 - total, abs=1e-12)

    # AND IT IS STRICTLY BELOW 1, which is the thing a reader trips over. Without this the test
    # would pass on a build whose shares summed to exactly 1 over the candidates - the defect this
    # whole arrangement exists to prevent.
    assert total < 1.0, "the candidate shares sum to 1, so the catch-all hypothesis has vanished"
    assert not_proposed > 0.0, "the catch-all is zero, so the set is being treated as exhaustive"


@pytest.mark.parametrize("composition", EVERY_SHAPE)
def test_the_unattested_mass_is_inside_the_candidate_sum_and_not_a_fourth_term(client, composition):
    adduct, charge = ("[M+H]+", 1) if composition != MAN5 else ("[M-H]-", -1)
    confidence = created(client, composition, adduct, charge)["confidence"]
    unattested = confidence["mass_on_unattested_classes"]
    assert unattested is not None
    assert unattested <= confidence["candidate_shares_sum"] + 1e-12, confidence


def test_the_explanation_travels_with_the_number_rather_than_living_in_a_document(client):
    """Item 3's instruction: the answer must be where somebody adding the shares up is looking."""
    confidence = created(client, G2F, "[M+H]+", 1)["confidence"]
    why = confidence["why_the_shares_do_not_sum_to_one"]
    assert why, "the response carries no explanation of why the shares fall short of 1"
    # It has to name the two sinks and the subset trap, or it does not answer the question asked.
    # Case-insensitive: the served text puts HYPOTHESIS SPACE in capitals for emphasis, and a test
    # that pins the casing of prose fails on an edit that improves it.
    for owed in ("hypothesis", "catch-all", "mass_on_unattested_classes"):
        assert owed in why.lower(), f"the explanation does not mention {owed!r}"
    # Ruling 2 of 27 September: the number is never named a probability. Checked as a SUBSTRING of
    # every key, because `probability_like` would pass an equality test.
    assert not [key for key in confidence if "probability" in key.lower()], list(confidence)


def test_the_dashboard_answers_it_on_the_same_card_as_the_number():
    """The dashboard is the other place a reader adds them up, so it carries the same answer."""
    page = (REPO / "src" / "wmxglycan" / "static" / "dashboard.html").read_text(encoding="utf-8")
    assert "why_the_shares_do_not_sum_to_one" in page
    assert "candidate_shares_sum" in page
    assert "shares_account_to_one" in page, "the page cannot report a broken accounting"
    assert "not a fourth share" in page, "the subset trap is not called out where the number is"
    # The old heading filed the unattested mass under "not on any candidate shown", which was
    # wrong: those candidates ARE shown. If it comes back, this fails.
    assert "Evidence mass that is not on any candidate shown" not in page


# --- item 4: the two unreachable decisions, derived ------------------------------------------------


def test_only_one_decision_is_reachable_and_the_other_two_say_what_would_reach_them():
    reachable = [one for one in decision_reachability() if one.reachable_today]
    blocked = [one for one in decision_reachability() if not one.reachable_today]
    assert [one.decision for one in reachable] == [Decision.IM_VALIDATION_REQUIRED]
    assert {one.decision for one in blocked} == {
        Decision.AI_ONLY,
        Decision.IM_VALIDATION_RECOMMENDED,
    }
    for one in blocked:
        assert one.blocked_by, f"{one.decision} is unreachable and names nothing that blocks it"
        assert one.what_would_reach_it, f"{one.decision} does not say what would reach it"
    # AI_ONLY is blocked TWICE and RECOMMENDED once. If that ever collapses to one reason for
    # both, the claim that AI_ONLY has two independent gates has stopped being true.
    by_decision = {one.decision: one for one in blocked}
    assert len(by_decision[Decision.AI_ONLY].blocked_by) == 2
    assert len(by_decision[Decision.IM_VALIDATION_RECOMMENDED].blocked_by) == 1


def test_the_published_reachability_moves_when_a_gate_is_patched(monkeypatch):
    """DERIVED, not listed - and this is what says so.

    A hand-written tuple would pass every assertion above and would not move here. That tuple is
    what api.py held until 27 September 2026, under a comment claiming it was derived.
    """
    # raising=True is the default and is load-bearing here: a typo in the attribute name would
    # otherwise patch nothing, the assertions below would test the unpatched module, and this test
    # would quietly become one that cannot fail.
    monkeypatch.setattr(ranking, "CORPUS_CAN_ESTABLISH_COMPLETENESS", True)
    opened = {one.decision for one in decision_reachability() if one.reachable_today}
    assert Decision.IM_VALIDATION_RECOMMENDED in opened, opened
    assert Decision.AI_ONLY not in opened, "AI_ONLY needs the model gate too, and it is still shut"

    monkeypatch.setattr(ranking, "VALIDATED_MODEL", object())
    both = {one.decision for one in decision_reachability() if one.reachable_today}
    assert both == set(Decision), both


def test_the_completeness_fact_lives_in_one_place(monkeypatch, shared):
    """`Coverage.truth_may_not_be_in_the_candidate_set` and `decision_reachability()` read ONE fact.

    WHY A TEST AND NOT A MUTATION. Replacing the property's body with a bare `return True` would be
    behaviour-identical today, so no test could kill that mutation and it would survive the sweep
    meaning nothing. The single-source link is therefore asserted directly: patch the constant and
    require BOTH readers to move. If the property ever goes back to a hardcoded True, this fails.
    """
    result = rank(
        shared["enumerator"].enumerate(MAN5),
        index=shared["index"],
        enumerator=shared["enumerator"],
    )
    assert result.coverage.truth_may_not_be_in_the_candidate_set is True

    monkeypatch.setattr(ranking, "CORPUS_CAN_ESTABLISH_COMPLETENESS", True)
    assert result.coverage.truth_may_not_be_in_the_candidate_set is False, (
        "the coverage property did not move with the constant, so the fact is written twice and"
        " one of the two copies will eventually be wrong"
    )
    assert Decision.IM_VALIDATION_RECOMMENDED in {
        one.decision for one in decision_reachability() if one.reachable_today
    }


def test_the_api_serves_the_derivation_rather_than_a_list(client):
    body = client.get("/v1/models/current").json()
    assert body["decision_reachable_today"] == [Decision.IM_VALIDATION_REQUIRED.value]
    unreachable = body["decision_unreachable_today"]
    assert set(unreachable) == {Decision.AI_ONLY.value, Decision.IM_VALIDATION_RECOMMENDED.value}
    for value, what_would in unreachable.items():
        assert what_would, f"{value} is served as unreachable with no route out of it"
    assert set(body["decision_values"]) == {one.value for one in Decision}


def test_the_enum_itself_carries_the_warning_a_reader_would_otherwise_have_to_look_for():
    """Item 4's instruction was "at the enum definition itself", so the docstring is the artefact.

    Weak on its own, which is why every test above it is about behaviour. It is here because a
    developer who greps for AI_ONLY lands on the enum and nowhere else.
    """
    text = Decision.__doc__ or ""
    assert "UNREACHABLE" in text
    assert "decision_reachability()" in text
    assert "stub" in text.lower()


def test_the_dashboard_names_the_unreachable_values_and_what_would_reach_them():
    page = (REPO / "src" / "wmxglycan" / "static" / "dashboard.html").read_text(encoding="utf-8")
    assert "decision_unreachable_today" in page
    assert "Unreachable today, and what would reach it" in page


# --- item 5: unreachable, and counted ---------------------------------------------------------------


def test_the_measured_reference_state_is_stated_as_unreachable_at_the_enum():
    text = CCSEvidenceState.__doc__ or ""
    assert "UNREACHABLE" in text, "the enum does not say the measured-reference path is unreachable"
    assert "milk oligosaccharide" in text
    assert "composition" in text


def test_the_api_says_it_holds_no_measured_cross_section_for_any_candidate(client):
    domain = client.get("/v1/models/current").json()["domain"]
    assert domain["predicts_ccs"] is False
    assert domain["holds_a_measured_cross_section_for_any_candidate"] is False
    note = domain["measured_reference_is_unreachable_note"]
    assert "UNREACHABLE" in note
    # The distinction the field exists to draw: not predicting is one claim, holding nothing is
    # another, and the first was being read as covering the second.
    assert "not merely thin" in note or "not merely" in note


def test_the_glycan_limitations_states_it_without_hedging():
    text = (REPO / "docs" / "GLYCAN_LIMITATIONS.md").read_text(encoding="utf-8")
    assert "UNREACHABLE" in text
    assert re.search(r"not\s+thin|not merely|Not thin", text), "the plain statement is missing"


def test_every_composition_returns_a_stated_absence_rather_than_omitting_the_field(client):
    """The behavioural half of item 5, and the part a document cannot promise."""
    for composition in EVERY_SHAPE:
        adduct, charge = ("[M+H]+", 1) if composition != MAN5 else ("[M-H]-", -1)
        body = created(client, composition, adduct, charge)
        assert "ccs_evidence" in body, body.keys()
        evidence = body["ccs_evidence"]
        assert evidence is not None, "the field was omitted, which reads as 'not applicable'"
        assert evidence["state"] != CCSEvidenceState.MEASURED_REFERENCE.value, (
            "a measured reference reached a candidate set, so the unreachability claim in"
            " READ_THIS_FIRST item 5 and at CCSEvidenceState is now false and must be re-read"
        )


def test_the_ranker_agrees_that_nothing_discriminates_between_candidates(shared):
    """Belt and braces on the same claim, one layer below the API."""
    for composition in EVERY_SHAPE:
        result = rank(
            shared["enumerator"].enumerate(composition),
            index=shared["index"],
            enumerator=shared["enumerator"],
        )
        assert result.ccs_evidence.discriminates_between_candidates is False
        assert result.ccs_evidence.state is not CCSEvidenceState.MEASURED_REFERENCE
