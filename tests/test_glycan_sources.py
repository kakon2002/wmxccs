"""The licence registry: a reuse-status claim traces back to something a person read."""

from datetime import date

import pytest

from wmxglycan.licensing import ReuseStatus
from wmxglycan.sources import (
    CLAIM_EXCEEDS_RECORD,
    DATASETS,
    NO_DOI_FOR_PUBLISHED_CLAIM,
    NO_RECORD_FOR_DOI,
    REGISTRY,
    STRUCTURE_CLAIM_EXCEEDS_DATASET,
    STRUWE_2016,
    SUGARBASE,
    DatasetLicence,
    SourceLicence,
    claim_problem,
    dataset_for,
    licence_for,
    structure_claim_problem,
)
from wmxglycan.sugarbase import dataset_version


def an_entry(**overrides):
    fields = dict(
        doi="10.1000/synthetic-fixture",
        citation="a synthetic fixture, not a real paper",
        title="a synthetic fixture",
        licence="CC BY 4.0",
        reuse_status=ReuseStatus.OPEN_ATTRIBUTION,
        evidence=("read on a fixture page",),
        reported_by="a fixture reader",
        reported_on=date(2026, 1, 1),
        attribution="cite the fixture",
    )
    fields.update(overrides)
    return SourceLicence(**fields)


# --- the entry that is actually there -------------------------------------------


def test_struwe_2016_is_recorded_with_its_evidence():
    assert STRUWE_2016.doi == "10.1039/c6cc06247d"
    assert STRUWE_2016.reuse_status is ReuseStatus.OPEN_ATTRIBUTION
    assert STRUWE_2016.permits_training
    assert "CC BY 3.0" in STRUWE_2016.licence
    # Three places on the publisher's site, each named, as the owner reported them.
    assert len(STRUWE_2016.evidence) == 3
    assert any("article page" in item for item in STRUWE_2016.evidence)
    assert any("issue listing" in item for item in STRUWE_2016.evidence)
    assert any("PDF header" in item for item in STRUWE_2016.evidence)
    assert STRUWE_2016.reported_by and STRUWE_2016.reported_on == date(2026, 9, 15)


def test_the_entry_names_who_read_the_page():
    # Stated by the owner directly, not inferred from the git identity:
    # "reported_by is Shawon Chakrabarty Kakon. I read the RSC page."
    assert STRUWE_2016.reported_by == "Shawon Chakrabarty Kakon"
    # Nothing in the evidence claims to quote the page verbatim.
    assert not any('"' in item for item in STRUWE_2016.evidence)


def test_the_stale_third_party_flag_is_recorded_settled_and_attributed():
    # A discrepancy that was met is written down with how it was resolved and
    # who found what, rather than left for the next reader to rediscover.
    note = next(note for note in STRUWE_2016.notes if "Europe PMC" in note)
    assert "stale" in note
    assert "Claude found" in note and "the owner read" in note  # who observed which


def test_the_record_says_none_of_the_values_were_recalled():
    assert any("reconstructed from memory" in note for note in STRUWE_2016.notes)
    assert any("awaited from a person" in note for note in STRUWE_2016.notes)


def test_the_supported_finding_is_the_narrow_one_and_marked_pending():
    # Comparable as [M+Na]+, over 7 per cent apart as [M-H]-. Nothing wider.
    finding = next(note for note in STRUWE_2016.notes if "LNH" in note)
    assert "[M+Na]+" in finding and "[M-H]-" in finding and "7 per cent" in finding
    assert "pending the paper" in finding
    assert "corrected" in finding  # the earlier wider summary is on record as corrected by the owner
    assert "as the owner summarised" in finding  # attributed, not asserted as fact


def test_the_registry_is_keyed_by_normalised_doi():
    assert licence_for("10.1039/C6CC06247D") is STRUWE_2016
    assert licence_for("  10.1039/c6cc06247d  ") is STRUWE_2016
    assert licence_for("10.1000/not-there") is None
    assert licence_for(None) is None and licence_for("") is None
    assert set(REGISTRY) == {"10.1039/c6cc06247d", "10.1039/c5an01092f", "10.1038/s41467-025-67069-w"}


def test_struwe_2015_is_recorded_with_the_footer_its_licence_was_read_from():
    from wmxglycan.sources import STRUWE_2015

    assert STRUWE_2015.doi == "10.1039/c5an01092f"
    assert STRUWE_2015.reuse_status is ReuseStatus.OPEN_ATTRIBUTION and STRUWE_2015.permits_training
    assert "CC BY 3.0" in STRUWE_2015.licence
    assert any("PDF footer" in item for item in STRUWE_2015.evidence)
    assert STRUWE_2015.reported_by == "Shawon Chakrabarty Kakon"
    # The title was not reported, and is not invented to fill the field.
    assert "not recorded" in STRUWE_2015.title
    # The two gaps a reader must not rediscover: the unstated gas, and the
    # composition that no structure can check.
    notes = " ".join(STRUWE_2015.notes)
    assert "UNSTATED" in notes and "Hofmann 2014" in notes and "NOT inferred" in notes
    assert "taken as transcribed" in notes


def test_the_non_commercial_no_derivatives_paper_is_recorded_so_nobody_rechecks_it():
    from wmxglycan.sources import SASTRE_TORANO_2025

    assert SASTRE_TORANO_2025.doi == "10.1038/s41467-025-67069-w"
    assert SASTRE_TORANO_2025.reuse_status is ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES
    assert not SASTRE_TORANO_2025.permits_training
    # The clause itself is stored verbatim, which is the point of the entry.
    quoted = " ".join(SASTRE_TORANO_2025.evidence)
    assert "You do not have permission under this licence to share adapted material" in quoted
    assert "Rights and permissions" in quoted and "15 September 2026" in quoted
    # And why this status rather than plain non_commercial.
    assert "no permitted use here at all" in " ".join(SASTRE_TORANO_2025.notes)


def test_a_claim_of_open_licence_on_the_nc_nd_paper_is_refused_by_its_record():
    problem = claim_problem("10.1038/s41467-025-67069-w", ReuseStatus.OPEN_ATTRIBUTION)
    assert problem is not None and "the record governs" in problem


# --- an entry cannot exist without its evidence --------------------------------------


def test_an_entry_with_no_evidence_is_refused():
    with pytest.raises(ValueError, match="records no evidence"):
        an_entry(evidence=())
    with pytest.raises(ValueError, match="records no evidence"):
        an_entry(evidence=("", "  "))


def test_an_entry_naming_nobody_is_refused():
    with pytest.raises(ValueError, match="names nobody"):
        an_entry(reported_by="   ")


def test_an_entry_with_no_doi_is_refused():
    with pytest.raises(ValueError, match="needs a DOI"):
        an_entry(doi="  ")


def test_an_entry_naming_no_licence_is_refused():
    # "open_attribution" with no licence named is a claim with nothing behind it.
    with pytest.raises(ValueError, match="names no licence"):
        an_entry(licence="  ")


def test_an_entry_without_a_real_date_is_refused():
    with pytest.raises(ValueError, match="no date"):
        an_entry(reported_on="not a date")


def test_an_entry_whose_status_is_a_bare_string_is_refused():
    with pytest.raises(ValueError, match="not a ReuseStatus"):
        an_entry(reuse_status="open_attribution")


def test_every_registered_entry_satisfies_its_own_invariants():
    for entry in REGISTRY.values():
        assert entry.evidence and entry.reported_by and entry.doi == entry.doi.casefold()
        assert isinstance(entry.reported_on, date) and isinstance(entry.reuse_status, ReuseStatus)


# --- datasets: a structure drawn from a shipped database ---------------------------


def test_sugarbase_is_registered_under_the_exact_provenance_its_loader_writes():
    assert SUGARBASE.provenance == dataset_version().provenance
    assert SUGARBASE.reuse_status is ReuseStatus.OPEN_ATTRIBUTION
    assert SUGARBASE.licence == "MIT"
    assert set(DATASETS) == {SUGARBASE.provenance}


def test_a_dataset_is_found_only_by_its_exact_provenance():
    # No fuzzy matching: inferring a licence from how a source is spelt is the
    # failure this registry exists to prevent.
    assert dataset_for(dataset_version().provenance) is SUGARBASE
    assert dataset_for("SugarBase") is None
    assert dataset_for("sugarbase v12 via glycowork") is None
    assert dataset_for(None) is None and dataset_for("") is None


def test_a_dataset_licence_needs_a_provenance_and_a_licence():
    with pytest.raises(ValueError, match="needs a provenance"):
        DatasetLicence(provenance=" ", licence="MIT", reuse_status=ReuseStatus.OPEN_ATTRIBUTION, attribution="x")
    with pytest.raises(ValueError, match="not a ReuseStatus"):
        DatasetLicence(provenance="p", licence="MIT", reuse_status="open_attribution", attribution="x")


# --- a structure's claim is backed by its paper or its database -------------------


def test_a_structure_from_the_registered_paper_is_backed_by_the_paper():
    # The structure names the paper by giving the measurement's own source.
    assert (
        structure_claim_problem(
            "10.1039/c6cc06247d", "the paper", ReuseStatus.OPEN_ATTRIBUTION, measurement_source="the paper"
        )
        is None
    )


def test_a_row_doi_does_not_back_a_structure_that_names_another_source():
    # Before this, any registered DOI on the row laundered any structure at all:
    # a structure from an unregistered database cleared on the measurement's licence.
    problem = structure_claim_problem(
        "10.1039/c6cc06247d", "some other database", ReuseStatus.OPEN_ATTRIBUTION, measurement_source="the paper"
    )
    assert problem is not None and "named as the measurement's own source" in problem
    # And with no measurement source given at all, the DOI backs nothing.
    assert structure_claim_problem("10.1039/c6cc06247d", "the paper", ReuseStatus.OPEN_ATTRIBUTION) is not None


def test_a_structure_from_sugarbase_is_backed_by_the_dataset_without_a_doi():
    provenance = dataset_version().provenance
    assert structure_claim_problem(None, provenance, ReuseStatus.OPEN_ATTRIBUTION) is None


def test_a_structure_of_unknown_origin_claiming_open_licence_is_unbacked():
    problem = structure_claim_problem(None, "somewhere", ReuseStatus.OPEN_ATTRIBUTION)
    assert problem is not None and "neither the row's DOI nor a registered dataset" in problem


def test_a_structure_with_a_doi_that_has_no_record_is_unbacked():
    problem = structure_claim_problem(
        "10.1000/no-record", "the paper", ReuseStatus.OPEN_ATTRIBUTION, measurement_source="the paper"
    )
    assert problem == NO_RECORD_FOR_DOI.format(claimed=ReuseStatus.OPEN_ATTRIBUTION, doi="10.1000/no-record")


def test_an_internal_structure_needs_no_backing():
    assert structure_claim_problem(None, "our lab", ReuseStatus.INTERNAL_PROPRIETARY) is None


def test_an_internal_claim_on_a_dataset_structure_is_refused_by_the_dataset_record():
    # The structure came from SugarBase; SugarBase's record says what it is.
    provenance = dataset_version().provenance
    problem = structure_claim_problem(None, provenance, ReuseStatus.INTERNAL_PROPRIETARY)
    assert problem == STRUCTURE_CLAIM_EXCEEDS_DATASET.format(
        claimed=ReuseStatus.INTERNAL_PROPRIETARY,
        source=provenance,
        recorded=ReuseStatus.OPEN_ATTRIBUTION,
        licence=SUGARBASE.licence,
    )


def test_an_empty_doi_is_no_doi_at_all():
    # Not reachable through a record, because the model refuses an empty DOI.
    # These are module functions with callers of their own, though, and "" must
    # not be mistaken for a stated DOI that merely has no licence record: the
    # two send a reader to different remedies.
    assert claim_problem("", ReuseStatus.OPEN_ATTRIBUTION) == NO_DOI_FOR_PUBLISHED_CLAIM.format(
        claimed=ReuseStatus.OPEN_ATTRIBUTION
    )
    assert claim_problem("", ReuseStatus.INTERNAL_PROPRIETARY) is None  # internal needs no DOI
    problem = structure_claim_problem("", "the paper", ReuseStatus.OPEN_ATTRIBUTION, measurement_source="the paper")
    assert problem is not None and "neither the row's DOI nor a registered dataset" in problem


def test_a_structure_source_is_compared_with_its_whitespace_stripped():
    # A transcription that pads a cell still names the same paper. Without the
    # strip the structure would be refused as unbacked for a trailing space.
    assert (
        structure_claim_problem(
            "10.1039/c6cc06247d", "  the paper  ", ReuseStatus.OPEN_ATTRIBUTION, measurement_source="the paper"
        )
        is None
    )


def test_a_synthetic_fixture_needs_no_record_and_may_carry_an_invented_doi():
    # Not a claim about any source: nothing to look up, nothing to record.
    assert claim_problem(None, ReuseStatus.SYNTHETIC_FIXTURE) is None
    assert claim_problem("10.1000/invented", ReuseStatus.SYNTHETIC_FIXTURE) is None
    assert structure_claim_problem(None, "anywhere", ReuseStatus.SYNTHETIC_FIXTURE) is None
    assert structure_claim_problem("10.1000/invented", "the paper", ReuseStatus.SYNTHETIC_FIXTURE, measurement_source="the paper") is None


def test_a_synthetic_fixture_cannot_cite_a_recorded_paper_or_dataset():
    # A DOI that is on record makes the record that paper's; the record governs.
    problem = claim_problem("10.1039/c6cc06247d", ReuseStatus.SYNTHETIC_FIXTURE)
    assert problem is not None and "the record governs" in problem
    provenance = dataset_version().provenance
    assert "the record governs" in structure_claim_problem(None, provenance, ReuseStatus.SYNTHETIC_FIXTURE)
    assert "the record governs" in structure_claim_problem(
        "10.1039/c6cc06247d", "the paper", ReuseStatus.SYNTHETIC_FIXTURE, measurement_source="the paper"
    )


def test_an_internal_claim_on_a_structure_from_the_paper_is_refused_by_the_paper_record():
    # Same rule as for the measurement's own claim: the record governs.
    problem = structure_claim_problem(
        "10.1039/c6cc06247d", "the paper", ReuseStatus.INTERNAL_PROPRIETARY, measurement_source="the paper"
    )
    assert problem is not None and "the record governs" in problem


def test_the_registries_are_read_only():
    # An entry is added by editing sources.py with its evidence, and by no other route.
    with pytest.raises(TypeError):
        REGISTRY["10.1000/injected"] = an_entry(doi="10.1000/injected")  # type: ignore[index]
    with pytest.raises(TypeError):
        DATASETS["anything"] = SUGARBASE  # type: ignore[index]
    assert licence_for("10.1000/injected") is None


def test_a_non_trainable_structure_claim_is_not_checked():
    assert structure_claim_problem(None, "anywhere", ReuseStatus.UNVERIFIED) is None
    assert structure_claim_problem(None, "anywhere", ReuseStatus.ACADEMIC_ONLY) is None


def test_a_bare_string_claim_is_refused_with_the_right_reason():
    # Not silently coerced and not compared by identity to an enum it can never equal.
    assert "not a ReuseStatus" in claim_problem(None, "internal_proprietary")
    assert "not a ReuseStatus" in structure_claim_problem(None, "x", "open_attribution")


# --- a row's claim is checked against the registry -------------------------------------


def test_a_matching_claim_is_backed():
    assert claim_problem("10.1039/C6CC06247D", ReuseStatus.OPEN_ATTRIBUTION) is None


def test_a_non_trainable_claim_needs_no_record():
    # Saying "unverified" is the safe direction; the gate refuses it anyway.
    # Every status the gate refuses is here, share-alike and excluded included.
    from wmxglycan.licensing import can_train_commercial

    refused = [status for status in ReuseStatus if not can_train_commercial(status)]
    assert ReuseStatus.OPEN_SHARE_ALIKE in refused and len(refused) >= 4
    for status in refused:
        assert claim_problem("10.1000/no-record", status) is None
        assert claim_problem(None, status) is None
        assert structure_claim_problem(None, "anywhere", status) is None


def test_an_internal_measurement_needs_no_doi():
    # Wellmatix's own data has no publication to point at.
    assert claim_problem(None, ReuseStatus.INTERNAL_PROPRIETARY) is None


def test_an_open_claim_with_no_doi_cannot_be_checked():
    problem = claim_problem(None, ReuseStatus.OPEN_ATTRIBUTION)
    assert problem == NO_DOI_FOR_PUBLISHED_CLAIM.format(claimed=ReuseStatus.OPEN_ATTRIBUTION)


def test_an_open_claim_for_an_unrecorded_doi_is_not_backed():
    problem = claim_problem("10.1000/no-record", ReuseStatus.OPEN_ATTRIBUTION)
    assert problem == NO_RECORD_FOR_DOI.format(claimed=ReuseStatus.OPEN_ATTRIBUTION, doi="10.1000/no-record")


def test_a_claim_wider_than_the_record_is_not_backed():
    # Claiming internal ownership of a CC BY paper is a different status from the recorded one.
    problem = claim_problem("10.1039/c6cc06247d", ReuseStatus.INTERNAL_PROPRIETARY)
    assert problem is not None
    assert "the record governs" in problem
    assert problem == CLAIM_EXCEEDS_RECORD.format(
        claimed=ReuseStatus.INTERNAL_PROPRIETARY,
        doi="10.1039/c6cc06247d",
        recorded=ReuseStatus.OPEN_ATTRIBUTION,
        licence=STRUWE_2016.licence,
    )
