"""CCS as evidence: the provenance guard, the five states, and the adapter across the wall.

The states exist so that "we looked and found nothing" cannot be confused with "we could not
look", "we hold it and may not show it", or "the lookup broke". Every one of those has been a
real bug somewhere in this repository's history in some other clothes, which is why each is a
separate state with its own obligations rather than a nullable value.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from wmxglycan.ccs_evidence import (
    NO_CORPUS_NAMED,
    WRONG_LEVEL_FOR_A_CANDIDATE,
    CCSEvidence,
    CCSEvidenceLookup,
    CCSEvidenceState,
    CCSReference,
    EvidenceLevel,
    HeldValues,
)
from wmxglycan.composition import parse_composition
from wmxglycan.reuse import ReuseStatus

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tools.glycan_ccs_evidence import (  # noqa: E402  (after the path fix, deliberately)
    CORPUS_NAME,
    STRUWE_2015_BLOCKERS,
    SeedCorpusEvidence,
)


def reference(**overrides) -> CCSReference:
    fields = dict(
        ccs=700.0,
        uncertainty=3.0,
        uncertainty_type="two_sd",
        adduct="[M+H]+",
        charge=1,
        polarity="positive",
        ims_type="TWIMS",
        drift_gas="N2",
        calibrant="dextran",
        source="a synthetic fixture",
        doi="10.0000/fixture",
        source_locator="Table S1",
        reuse_status=ReuseStatus.SYNTHETIC_FIXTURE,
    )
    fields.update(overrides)
    return CCSReference(**fields)


# --- the provenance IS the guard, and no flag is trusted ------------------------------------------


def test_a_reference_needs_an_uncertainty_type_that_is_not_unknown():
    # Hard constraint 7, and the reason the 89 Struwe 2015 values are held rather than shown.
    for given in ("unknown", "UNKNOWN", "", "   ", "none"):
        with pytest.raises(ValueError, match="uncertainty with a stated type"):
            reference(uncertainty_type=given)


def test_a_reference_needs_the_gas_its_value_refers_to():
    for given in ("UNSTATED", "unstated", "", "unknown"):
        with pytest.raises(ValueError, match="the gas its value refers to"):
            reference(drift_gas=given)


def test_a_reference_needs_a_doi():
    with pytest.raises(ValueError, match="needs a DOI"):
        reference(doi="   ")


@pytest.mark.parametrize(
    "missing", ["ccs", "uncertainty", "uncertainty_type", "adduct", "charge", "polarity",
                "ims_type", "drift_gas", "source", "doi", "source_locator", "reuse_status"]
)
def test_every_provenance_field_is_required(missing):
    # A prediction cannot supply these. That is the guard: not a Literal[True] flag the producer
    # sets about itself, which forbids only the one spelling nobody would use and is satisfied
    # by its own default - the same asymmetry as extra="forbid" letting missing keys through.
    fields = {
        "ccs": 700.0, "uncertainty": 3.0, "uncertainty_type": "two_sd", "adduct": "[M+H]+",
        "charge": 1, "polarity": "positive", "ims_type": "TWIMS", "drift_gas": "N2",
        "source": "x", "doi": "10.0/x", "source_locator": "T1",
        "reuse_status": ReuseStatus.SYNTHETIC_FIXTURE,
    }
    fields.pop(missing)
    with pytest.raises(Exception):
        CCSReference(**fields)


def test_a_reference_cannot_be_built_from_a_bare_number():
    with pytest.raises(Exception):
        CCSReference(ccs=700.0)


def test_no_field_claims_the_value_is_measured():
    # The flag that used to be here was removed on review. If it comes back, so does the
    # weakness: a self-asserted boolean that nothing checks.
    assert not any(
        "is_measured" in name or "is_prediction" in name for name in CCSReference.model_fields
    )


# --- each state carries what it claims -------------------------------------------------------------


def test_the_default_state_is_not_consulted():
    # So a ranker with no adapter cannot report an absence it never looked for.
    assert CCSEvidence().state is CCSEvidenceState.NOT_CONSULTED
    assert "nothing looked" in CCSEvidence().summary()


def test_a_measured_reference_must_carry_the_reference():
    with pytest.raises(ValueError, match="must carry the reference"):
        CCSEvidence(state=CCSEvidenceState.MEASURED_REFERENCE)


def test_no_other_state_may_carry_a_reference_value():
    for state in CCSEvidenceState:
        if state is CCSEvidenceState.MEASURED_REFERENCE:
            continue
        with pytest.raises(ValueError, match="must not carry a reference value"):
            CCSEvidence(state=state, reference=reference(), failure="x", corpus_searched="c",
                        records_consulted=0)


def test_a_held_finding_carries_its_count_and_blockers_and_never_a_value():
    held = CCSEvidence(
        state=CCSEvidenceState.HELD_NOT_RELEASABLE,
        held=HeldValues(records=8, blockers=("gas unstated",), what_would_release_them="nothing"),
    )
    assert held.reference is None
    assert "8 value(s) exist" in held.summary()
    with pytest.raises(ValueError, match="must not carry held values"):
        CCSEvidence(
            state=CCSEvidenceState.NONE_IN_CORPUS_SEARCHED,
            corpus_searched="c",
            records_consulted=0,
            held=HeldValues(records=1, blockers=("x",), what_would_release_them="y"),
        )


def test_held_values_need_at_least_one_blocker_and_a_positive_count():
    with pytest.raises(Exception):
        HeldValues(records=0, blockers=("x",), what_would_release_them="y")
    with pytest.raises(Exception):
        HeldValues(records=1, blockers=(), what_would_release_them="y")


def test_a_searched_absence_cannot_be_constructed_without_a_denominator():
    # An absence with no corpus and no count is an unexamined field, and the two read
    # identically downstream.
    with pytest.raises(ValueError, match="must name the corpus it searched"):
        CCSEvidence(state=CCSEvidenceState.NONE_IN_CORPUS_SEARCHED)
    with pytest.raises(ValueError, match="must name the corpus it searched"):
        CCSEvidence(state=CCSEvidenceState.NONE_IN_CORPUS_SEARCHED, corpus_searched="c")
    assert NO_CORPUS_NAMED.startswith("a searched-absence finding")


def test_a_failed_lookup_must_say_what_failed():
    with pytest.raises(ValueError, match="must say what failed"):
        CCSEvidence(state=CCSEvidenceState.LOOKUP_FAILED)
    failed = CCSEvidence(state=CCSEvidenceState.LOOKUP_FAILED, failure="KeyError: nope")
    assert "not an absence" in failed.summary()


def test_structure_level_evidence_must_name_the_structure_it_was_measured_on():
    with pytest.raises(ValueError, match="must name the structure"):
        CCSEvidence(
            state=CCSEvidenceState.MEASURED_REFERENCE,
            level=EvidenceLevel.STRUCTURE,
            reference=reference(),
        )


# --- scope: one measurement is never many confirmations ---------------------------------------------


def test_composition_and_ion_level_evidence_cannot_reach_a_candidate():
    for level in (EvidenceLevel.COMPOSITION, EvidenceLevel.ION):
        evidence = CCSEvidence(
            state=CCSEvidenceState.MEASURED_REFERENCE, level=level, reference=reference()
        )
        with pytest.raises(ValueError, match="cannot be attached to an individual candidate"):
            evidence.for_candidate("any key")
    assert "was not measured on" in WRONG_LEVEL_FOR_A_CANDIDATE


def test_structure_level_evidence_reaches_only_the_structure_it_names():
    evidence = CCSEvidence(
        state=CCSEvidenceState.MEASURED_REFERENCE,
        level=EvidenceLevel.STRUCTURE,
        reference=reference(measured_on_canonical_key="the-right-one"),
    )
    assert evidence.for_candidate("the-right-one").ccs == 700.0
    with pytest.raises(ValueError, match="and not on"):
        evidence.for_candidate("a-different-structure")


def test_a_raise_rather_than_a_none_because_the_mistake_is_silent():
    # A None would be absorbed by a caller that then showed nothing, and one measurement
    # rendered as 167 per-structure claims is silent by nature.
    evidence = CCSEvidence(
        state=CCSEvidenceState.MEASURED_REFERENCE,
        level=EvidenceLevel.COMPOSITION,
        reference=reference(),
    )
    with pytest.raises(ValueError):
        evidence.for_candidate("x")


def test_composition_level_evidence_says_it_cannot_discriminate():
    evidence = CCSEvidence(
        state=CCSEvidenceState.MEASURED_REFERENCE,
        level=EvidenceLevel.COMPOSITION,
        reference=reference(),
    )
    assert evidence.discriminates_between_candidates is False


# --- the adapter, against the real corpus ------------------------------------------------------------


@pytest.fixture(scope="module")
def adapter():
    return SeedCorpusEvidence()


def test_the_adapter_satisfies_the_protocol_the_glycan_layer_declares(adapter):
    assert isinstance(adapter, CCSEvidenceLookup)


def test_a_stub_defined_in_the_glycan_layers_own_tests_satisfies_the_protocol():
    # The Protocol must be satisfiable without importing anything from the CCS core, or the
    # wall is decorative.
    class Stub:
        def evidence_for(self, composition, adduct, charge):
            return CCSEvidence()

    assert isinstance(Stub(), CCSEvidenceLookup)


def test_man5_returns_held_values_from_the_real_corpus(adapter):
    # 8 rows, not 19: 19 is Hex5HexNAc2 across every adduct, and this is the singly-charged
    # deprotonated ion. Pinned because a key normalisation that widened the match would look
    # like more evidence.
    found = adapter.evidence_for(parse_composition("Hex5HexNAc2"), "[M-H]-", -1)
    assert found.state is CCSEvidenceState.HELD_NOT_RELEASABLE
    assert found.held.records == 8
    assert found.reference is None
    assert found.level is EvidenceLevel.ION


def test_both_blockers_travel_and_neither_resolves_alone(adapter):
    found = adapter.evidence_for(parse_composition("Hex5HexNAc2"), "[M-H]-", -1)
    assert len(found.held.blockers) == 2
    assert found.held.any_blocker_resolves_alone is False
    joined = " ".join(found.held.blockers)
    assert "drift gas is UNSTATED" in joined
    assert "absent from the source rather than unread" in joined
    # A single-blocker list would invite "read Hofmann 2014 and you have data", which is false.
    assert "cannot be resolved by reading anything" in found.held.what_would_release_them


def test_the_conformer_rows_are_not_presented_as_independent_measurements(adapter):
    found = adapter.evidence_for(parse_composition("Hex5HexNAc2"), "[M-H]-", -1)
    assert "numbered conformers of one ion" in found.held.what_would_release_them


def test_the_fucosylated_compositions_return_a_searched_absence_with_a_denominator(adapter):
    for text in ("Hex3HexNAc4Fuc1", "Hex4HexNAc4Fuc1", "Hex5HexNAc4Fuc1"):
        found = adapter.evidence_for(parse_composition(text), "[M-H]-", -1)
        assert found.state is CCSEvidenceState.NONE_IN_CORPUS_SEARCHED
        assert found.corpus_searched == CORPUS_NAME
        assert found.records_consulted == 89


def test_a_renormalised_composition_spelling_still_finds_the_same_ion(adapter):
    # The classic fixture-only failure: a key mismatch reads as a scientific absence. Both
    # spellings must reach the held values.
    one = adapter.evidence_for(parse_composition("Hex5HexNAc2"), "[M-H]-", -1)
    other = adapter.evidence_for(parse_composition("HexNAc2Hex5"), "[M-H]-", -1)
    assert one.state is other.state is CCSEvidenceState.HELD_NOT_RELEASABLE
    assert one.held.records == other.held.records == 8


def test_an_ion_the_corpus_does_not_hold_is_an_absence_and_not_a_failure(adapter):
    found = adapter.evidence_for(parse_composition("Hex5HexNAc2"), "[M+2H]2+", 2)
    assert found.state is CCSEvidenceState.NONE_IN_CORPUS_SEARCHED


def test_the_reachable_states_are_counted_rather_than_asserted(adapter):
    reachable = adapter.reachable_states()
    assert "REACHABLE" in reachable[CCSEvidenceState.HELD_NOT_RELEASABLE.value]
    assert "REACHABLE" in reachable[CCSEvidenceState.LOOKUP_FAILED.value]
    measured = reachable[CCSEvidenceState.MEASURED_REFERENCE.value]
    assert measured.startswith("NOT REACHABLE")
    assert "0 of 541 cleared record(s) carry a composition" in measured


def test_the_unreachability_claim_fails_the_moment_a_cleared_record_carries_a_composition(adapter):
    # Held to the same standard as AI_ONLY: derived by counting, so it changes by itself rather
    # than being a sentence that quietly becomes false.
    with_composition = sum(
        1
        for record in adapter._cleared
        if getattr(getattr(record, "analyte", None), "composition", None) is not None
    )
    assert with_composition == 0, (
        "a cleared glycan record now carries a composition, so MEASURED_REFERENCE is reachable"
        " and the published reachability statement must be re-read"
    )


def test_every_cleared_glycan_record_is_a_milk_oligosaccharide_without_a_composition(adapter):
    # The measured reason MEASURED_REFERENCE cannot fire. Not a code gap: the two seed files are
    # complementary, and this is the half with the values.
    glycans = [
        record
        for record in adapter._cleared
        if type(getattr(record, "analyte", None)).__name__ == "GlycanAnalyte"
    ]
    assert len(glycans) == 24
    assert all(record.analyte.composition is None for record in glycans)
    assert all(record.analyte.iupac_condensed for record in glycans)


def test_the_blocker_wording_is_importable_rather_than_matched_as_prose():
    assert len(STRUWE_2015_BLOCKERS) == 2
    assert "hard constraint 7" in STRUWE_2015_BLOCKERS[1]


# --- the licence gate applies to a reference we are about to SHOW ------------------------------------


@pytest.mark.parametrize(
    "status",
    [ReuseStatus.UNVERIFIED, ReuseStatus.EXCLUDED, ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES],
)
def test_a_status_with_no_permitted_use_cannot_carry_a_reference(status):
    # The first version required the FIELD to be present and accepted every value, including
    # the `unverified` DEFAULT and `excluded`, which reuse.py says have no permitted use at all.
    # Requiring a field to exist is not requiring it to permit anything - the same asymmetry as
    # extra="forbid" rejecting unknown keys while missing ones sail through.
    with pytest.raises(ValueError, match="NO PERMITTED USE"):
        reference(reuse_status=status)


@pytest.mark.parametrize(
    "status",
    [
        ReuseStatus.OPEN_ATTRIBUTION,
        ReuseStatus.OPEN_SHARE_ALIKE,
        ReuseStatus.NON_COMMERCIAL,
        ReuseStatus.ACADEMIC_ONLY,
        ReuseStatus.INTERNAL_PROPRIETARY,
        ReuseStatus.SYNTHETIC_FIXTURE,
    ],
)
def test_a_status_with_some_permitted_use_is_accepted(status):
    # The positive control: if the guard rejected everything it would be indistinguishable from
    # a guard that works, and the six usable tiers would be unreachable.
    assert reference(reuse_status=status).reuse_status is status


def test_the_three_refused_statuses_are_exactly_those_with_no_permitted_use():
    # Derived from the predicates rather than listed here, so the guard cannot drift from the
    # tiers it is enforcing.
    from wmxglycan.reuse import can_redistribute, can_train_commercial, is_inference_only

    refused = {
        status
        for status in ReuseStatus
        if not (
            can_train_commercial(status) or can_redistribute(status) or is_inference_only(status)
        )
    }
    assert refused == {
        ReuseStatus.UNVERIFIED,
        ReuseStatus.EXCLUDED,
        ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES,
    }
    for status in refused:
        with pytest.raises(ValueError):
            reference(reuse_status=status)


# --- the discriminates flag needs the level that could justify it ------------------------------------


@pytest.mark.parametrize("level", [EvidenceLevel.COMPOSITION, EvidenceLevel.ION])
def test_evidence_below_structure_level_cannot_claim_to_discriminate(level):
    # `_decide` used to trust this flag without reading the level, so a composition-level
    # finding with the flag set made the rule "no cross section is held for this structure"
    # report fires=False - the one rule between a shared measurement and a per-structure claim.
    # Both ends are guarded now: the level is read there, and the combination is refused here.
    with pytest.raises(ValueError, match="cannot also claim to discriminate"):
        CCSEvidence(
            state=CCSEvidenceState.MEASURED_REFERENCE,
            level=level,
            reference=reference(),
            discriminates_between_candidates=True,
        )


def test_structure_level_evidence_may_claim_to_discriminate():
    found = CCSEvidence(
        state=CCSEvidenceState.MEASURED_REFERENCE,
        level=EvidenceLevel.STRUCTURE,
        reference=reference(measured_on_canonical_key="k"),
        discriminates_between_candidates=True,
    )
    assert found.discriminates_between_candidates is True
    assert found.for_candidate("k").ccs == 700.0
