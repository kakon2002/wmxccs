"""The registry is evidence, and these tests are what make it evidence.

There are two jobs here and they are different jobs.

The first is CONTENT. M0 asks for the registry seeded with every source named in
CONTEXT.md, at the status stated there, with the reader and the date given. A
registry that is well shaped but missing a source, or holding one at a status
nobody actually read, is worse than no registry at all: it reads as a licence
check that has happened. So every source CONTEXT.md names is pinned here by its
key and by its status, and the set is pinned closed, so that a twelfth entry
cannot arrive without a test saying so.

The second job is the GUARDS. Each one is asserted in both directions: the valid
record is built and kept, and the invalid one raises with a message a person can
act on. A guard only ever exercised in the passing direction catches nothing,
and a claim check that never refuses is the licence gate failing open.
"""

from __future__ import annotations

from datetime import date

import pytest

from wmxccs import sources
from wmxccs.reuse import ReuseStatus, can_train_commercial
from wmxccs.sources import (
    DatasetLicence,
    SourceLicence,
    claim_problem,
    component_claim_problem,
    dataset_for,
    licence_for,
    source_for,
    unverified_sources,
)

STRUWE_2016_DOI = "10.1039/c6cc06247d"
STRUWE_2015_DOI = "10.1039/c5an01092f"
STEROID_DOI = "10.1021/jasms.2c00196"
BAYESIAN_DOI = "10.1021/acs.analchem.5c06667"

# A DOI that is deliberately not in the registry, for every "no record" case.
UNRECORDED_DOI = "10.9999/nobody.has.read.this"

SEED_2016 = sources.SEED_STRUWE_2016.provenance

# An analyte identity whose origin nobody recorded. The case the component check
# exists for: "open_attribution" typed against something of unknown provenance.
UNKNOWN_ORIGIN = "an identity of unknown origin"


def licence_record(**overrides) -> SourceLicence:
    """A minimal VALID registry entry, so each guard test changes exactly one field.

    Not a real source and not registered: it exists only so that the failing
    direction of every schema guard can be reached from a record that is
    otherwise beyond reproach.
    """
    fields = dict(
        citation="A Reader et al., Journal of Nothing 2026, 1, 1",
        title="A paper whose licence was read",
        licence="Creative Commons Attribution 4.0 International (CC BY 4.0)",
        reuse_status=ReuseStatus.OPEN_ATTRIBUTION,
        evidence=("Licence badge on the article page, read 18 September 2026",),
        reported_by="A Reader",
        reported_on=date(2026, 9, 18),
        attribution="Correct acknowledgement.",
        doi=UNRECORDED_DOI,
    )
    fields.update(overrides)
    return SourceLicence(**fields)


def dataset_record(**overrides) -> DatasetLicence:
    fields = dict(
        provenance="a test dataset",
        licence="Creative Commons Attribution 4.0 International (CC BY 4.0)",
        reuse_status=ReuseStatus.OPEN_ATTRIBUTION,
        attribution="Correct acknowledgement.",
    )
    fields.update(overrides)
    return DatasetLicence(**fields)


# =========================================================================================
# CONTENT: every source CONTEXT.md names, at the status it states
# =========================================================================================

# (readable name, the entry, the key it is found by, the status CONTEXT.md states)
CONTEXT_SOURCES = (
    ("Struwe 2016 Chem Commun", sources.STRUWE_2016, STRUWE_2016_DOI, ReuseStatus.OPEN_ATTRIBUTION),
    ("Struwe 2015 Analyst", sources.STRUWE_2015, STRUWE_2015_DOI, ReuseStatus.OPEN_ATTRIBUTION),
    ("steroid interplatform 2022", sources.STEROID_INTERPLATFORM_2022, STEROID_DOI, ReuseStatus.UNVERIFIED),
    ("METLIN-CCS", sources.METLIN_CCS, "https://metlin.scripps.edu/", ReuseStatus.UNVERIFIED),
    ("Bush Lab CCS database", sources.BUSH_LAB_CCS, "https://biophysicalms.org/ccsdatabase", ReuseStatus.UNVERIFIED),
    ("Bayesian harmonization 2026", sources.BAYESIAN_HARMONIZATION_2026, BAYESIAN_DOI, ReuseStatus.UNVERIFIED),
    ("CCSbase", sources.CCSBASE, "https://ccsbase.net/", ReuseStatus.ACADEMIC_ONLY),
    (
        "Sastre Torano 2025",
        sources.SASTRE_TORANO_2025,
        "10.1038/s41467-025-67069-w",
        ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES,
    ),
    (
        "Manabe 2022",
        sources.MANABE_2022,
        "https://pubs.acs.org/journal/jamsef",
        ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES,
    ),
    ("AllCCS2", sources.ALLCCS2, "https://allccs.zhulab.cn/", ReuseStatus.EXCLUDED),
    ("GlycoMob", sources.GLYCOMOB, "http://www.glycomob.org/", ReuseStatus.EXCLUDED),
)

CONTEXT_IDS = [name for name, _entry, _key, _status in CONTEXT_SOURCES]


@pytest.mark.parametrize(("name", "entry", "key", "status"), CONTEXT_SOURCES, ids=CONTEXT_IDS)
def test_every_source_named_in_context_is_in_the_registry_at_the_status_context_states(name, entry, key, status):
    assert entry.reuse_status is status, f"{name} is recorded at the wrong status"
    assert entry.key == key, f"{name} is not findable by the key it should be findable by"
    assert source_for(key) is entry, f"{name} is not reachable through the registry lookup"


def test_the_registry_holds_exactly_the_sources_context_names_and_no_others():
    # Pinned closed on purpose. An entry added without its own row above is an
    # entry nobody asserted the status of, which is the failure this file exists
    # to prevent.
    assert set(sources.REGISTRY) == {key for _name, _entry, key, _status in CONTEXT_SOURCES}


@pytest.mark.parametrize(("name", "entry", "key", "status"), CONTEXT_SOURCES, ids=CONTEXT_IDS)
def test_every_registry_entry_names_a_reader_and_a_date_for_its_report(name, entry, key, status):
    # Constraint 2: a licence claim is only valid when backed by an entry naming
    # who read the terms and when. An entry with neither is not evidence.
    assert entry.reported_by.strip(), f"{name} names nobody who read it"
    assert isinstance(entry.reported_on, date), f"{name} has no date for its report"


def test_struwe_2016_records_the_reader_the_date_and_the_rsc_page_it_was_read_on():
    entry = sources.STRUWE_2016
    assert entry.reported_by == "Shawon Chakrabarty Kakon"
    assert entry.reported_on == date(2026, 9, 15)
    assert "Attribution 3.0" in entry.licence
    joined = " ".join(entry.evidence)
    assert "RSC article page" in joined, "the evidence does not say which page was read"
    assert "Creative Commons Attribution 3.0" in joined


def test_struwe_2015_records_its_own_reader_and_a_later_date_than_the_2016_paper():
    entry = sources.STRUWE_2015
    assert entry.reported_by == "Shawon Chakrabarty Kakon"
    assert entry.reported_on == date(2026, 9, 17)
    assert entry.evidence, "the 2015 paper's entry records no evidence"


def test_only_the_two_struwe_papers_are_recorded_as_trainable():
    # Constraint 2 in the direction that matters commercially: nothing else in
    # the registry may enter a fit, and CCSbase in particular is blocked while
    # the platform is treated as commercial.
    trainable = {entry.key for entry in sources.REGISTRY.values() if entry.permits_training}
    assert trainable == {STRUWE_2016_DOI, STRUWE_2015_DOI}
    assert sources.CCSBASE.permits_training is False


def test_unverified_sources_returns_exactly_the_four_sources_nobody_has_read_the_terms_for():
    assert {entry.key for entry in unverified_sources()} == {
        STEROID_DOI,
        "https://metlin.scripps.edu/",
        "https://biophysicalms.org/ccsdatabase",
        BAYESIAN_DOI,
    }


def test_unverified_sources_lists_only_entries_whose_status_is_unverified():
    assert [entry for entry in unverified_sources() if entry.reuse_status is not ReuseStatus.UNVERIFIED] == []


@pytest.mark.parametrize("entry", unverified_sources(), ids=lambda e: e.key)
def test_every_unverified_entry_says_what_would_settle_it_so_nobody_re_reads_the_same_page(entry):
    assert entry.what_to_check.strip(), (
        f"{entry.key} is unverified and does not say what would settle it, so the next person"
        " spends the same day on the same page"
    )


def test_the_registry_is_a_read_only_mapping_a_caller_cannot_add_to():
    with pytest.raises(TypeError):
        sources.REGISTRY[UNRECORDED_DOI] = licence_record()


def test_the_registry_is_a_read_only_mapping_a_caller_cannot_overwrite_an_entry_in():
    with pytest.raises(TypeError):
        sources.REGISTRY[STRUWE_2016_DOI] = licence_record()


def test_the_dataset_registry_is_a_read_only_mapping_a_caller_cannot_add_to():
    with pytest.raises(TypeError):
        sources.DATASETS["a dataset nobody registered"] = dataset_record()


def test_the_two_seed_transcriptions_are_registered_as_datasets_under_their_papers_licence():
    assert set(sources.DATASETS) == {
        sources.SEED_STRUWE_2016.provenance,
        sources.SEED_STRUWE_2015.provenance,
    }
    assert sources.SEED_STRUWE_2016.reuse_status is ReuseStatus.OPEN_ATTRIBUTION
    assert sources.SEED_STRUWE_2016.derived_from_doi == STRUWE_2016_DOI
    assert sources.SEED_STRUWE_2015.derived_from_doi == STRUWE_2015_DOI


# =========================================================================================
# SCHEMA GUARDS, each in the failing direction
# =========================================================================================


def test_a_record_with_a_doi_and_a_reader_and_evidence_is_accepted():
    entry = licence_record()
    assert entry.key == UNRECORDED_DOI
    assert entry.permits_training is True


def test_an_entry_with_neither_a_doi_nor_a_url_is_refused_because_nobody_can_re_check_it():
    with pytest.raises(ValueError, match="neither a DOI nor a URL"):
        licence_record(doi=None, url=None)


def test_an_entry_with_only_a_url_is_accepted_and_is_found_by_that_url():
    entry = licence_record(doi=None, url="https://example.invalid/terms")
    assert entry.key == "https://example.invalid/terms"


@pytest.mark.parametrize("blank", ["", "   "])
def test_an_entry_naming_no_licence_at_all_is_refused(blank):
    with pytest.raises(ValueError, match="names no licence"):
        licence_record(licence=blank)


def test_an_entry_with_no_evidence_behind_it_is_refused():
    with pytest.raises(ValueError, match="records no evidence"):
        licence_record(evidence=())


def test_an_entry_whose_evidence_is_blank_text_is_refused_as_firmly_as_one_with_none():
    # A whitespace string is the shape a placeholder takes when somebody wants
    # the record to construct. It is not evidence.
    with pytest.raises(ValueError, match="records no evidence"):
        licence_record(evidence=("Licence badge read", "   "))


def test_an_entry_with_nobody_named_as_having_read_it_is_refused():
    with pytest.raises(ValueError, match="names nobody who reported it"):
        licence_record(reported_by="  ")


@pytest.mark.parametrize("not_a_date", ["2026-09-18", 20260918, None])
def test_an_entry_whose_report_date_is_not_a_date_is_refused(not_a_date):
    with pytest.raises(ValueError, match="no date for its report"):
        licence_record(reported_on=not_a_date)


@pytest.mark.parametrize("not_a_status", ["open_attribution", "made up", None])
def test_an_entry_whose_reuse_status_is_not_a_reuse_status_is_refused(not_a_status):
    with pytest.raises(ValueError, match="not a ReuseStatus"):
        licence_record(reuse_status=not_a_status)


def test_an_unverified_entry_that_does_not_say_what_would_settle_it_is_refused():
    with pytest.raises(ValueError, match="does not say what would settle it"):
        licence_record(reuse_status=ReuseStatus.UNVERIFIED, what_to_check="")


def test_an_unverified_entry_is_accepted_once_it_says_which_page_would_settle_it():
    entry = licence_record(
        reuse_status=ReuseStatus.UNVERIFIED,
        what_to_check="Read the licence badge on the publisher's page and record it verbatim.",
    )
    assert entry.reuse_status is ReuseStatus.UNVERIFIED
    assert entry.permits_training is False


@pytest.mark.parametrize("status", [s for s in ReuseStatus if can_train_commercial(s)], ids=lambda s: s.value)
def test_an_entry_that_permits_training_and_names_no_attribution_is_refused(status):
    # The obligation a permissive licence asks for in return. An entry that
    # clears training and names none is an obligation nobody can discharge.
    with pytest.raises(ValueError, match="names no attribution"):
        licence_record(reuse_status=status, attribution="   ")


def test_an_entry_that_permits_nothing_may_leave_attribution_blank():
    # The guard is scoped to trainable entries: a blocked source has no
    # attribution to record, and requiring one would push somebody to invent one.
    entry = licence_record(reuse_status=ReuseStatus.EXCLUDED, attribution="")
    assert entry.permits_training is False


def test_evidence_passed_as_a_list_is_held_as_a_tuple_so_a_frozen_record_holds_nothing_mutable():
    entry = licence_record(evidence=["the licence badge", "the PDF footer"])
    assert isinstance(entry.evidence, tuple)
    assert entry.evidence == ("the licence badge", "the PDF footer")
    hash(entry)  # a record holding a list is unhashable, and every later test on it trips


def test_notes_passed_as_a_list_are_held_as_a_tuple():
    entry = licence_record(notes=["one discrepancy, settled"])
    assert isinstance(entry.notes, tuple)
    assert entry.notes == ("one discrepancy, settled",)
    hash(entry)


# Was an xfail against a real defect: evidence was coerced with tuple(), so a
# bare string became one piece of evidence per character and the entry read as
# thoroughly evidenced while holding nothing. A bare string is now refused.
def test_evidence_given_as_one_bare_string_is_refused_rather_than_split_into_characters():
    # Both spellings must be refused, and for the SAME reason. A string with no
    # spaces used to pass every check as one piece of evidence per character,
    # leaving the entry reading as thoroughly evidenced while holding nothing;
    # a string with spaces failed for the wrong reason, "records no evidence",
    # which points the reader at a missing field rather than a wrong type.
    for spelling in ("RSCarticlepage", "RSC article page"):
        with pytest.raises(ValueError, match="not one string"):
            licence_record(evidence=spelling)


def test_notes_given_as_one_bare_string_are_refused_the_same_way():
    with pytest.raises(ValueError, match="not one string"):
        licence_record(notes="settled on the publisher's page")


def test_evidence_given_as_a_list_is_still_coerced_to_a_tuple():
    # The coercion itself is worth keeping: a frozen record holding a mutable,
    # unhashable field trips every later membership test.
    entry = licence_record(evidence=["RSC article page", "the issue listing"])
    assert entry.evidence == ("RSC article page", "the issue listing")


@pytest.mark.parametrize("dataset", list(sources.DATASETS.values()), ids=lambda d: d.provenance)
def test_no_registered_dataset_claims_more_than_the_paper_it_was_transcribed_from(dataset):
    # derived_from_doi is the only thing tying a seeded corpus back to a licence
    # somebody read. Nothing in the module consults it, so this holds the
    # invariant from outside: a dataset entry may not out-claim its own paper.
    if dataset.derived_from_doi is None:
        pytest.skip("this dataset was not transcribed from a publication")
    paper = licence_for(dataset.derived_from_doi)
    assert paper is not None, f"{dataset.provenance} derives from a DOI with no licence record"
    assert dataset.reuse_status is paper.reuse_status, (
        f"{dataset.provenance} claims {dataset.reuse_status} but the record for the paper it was"
        f" transcribed from says {paper.reuse_status}"
    )


@pytest.mark.parametrize("dataset", list(sources.DATASETS.values()), ids=lambda d: d.provenance)
def test_every_registered_dataset_that_permits_training_names_an_attribution(dataset):
    # The same obligation SourceLicence enforces on a publication. DatasetLicence
    # does not check it, so it is checked here on the entries that exist.
    if not can_train_commercial(dataset.reuse_status):
        return
    assert dataset.attribution.strip(), (
        f"{dataset.provenance} permits training and names no attribution, which is a licence"
        " obligation nobody can discharge"
    )


def test_a_dataset_entry_with_no_provenance_string_is_refused():
    with pytest.raises(ValueError, match="needs a provenance string and a licence"):
        dataset_record(provenance="   ")


def test_a_dataset_entry_naming_no_licence_is_refused():
    with pytest.raises(ValueError, match="needs a provenance string and a licence"):
        dataset_record(licence="")


def test_a_dataset_entry_whose_reuse_status_is_not_a_reuse_status_is_refused():
    with pytest.raises(ValueError, match="not a ReuseStatus"):
        dataset_record(reuse_status="open_attribution")


# =========================================================================================
# LOOKUP
# =========================================================================================


@pytest.mark.parametrize(
    "spelling",
    ["10.1039/c6cc06247d", "10.1039/C6CC06247D", "  10.1039/C6cc06247D  ", "10.1039/C6CC06247d"],
)
def test_a_doi_finds_its_record_however_it_is_spelled(spelling):
    # A DOI is case-insensitive. Matching one case-sensitively silently loses the
    # record, and a lost record reads exactly like a source nobody ever read.
    assert licence_for(spelling) is sources.STRUWE_2016


def test_a_doi_with_no_record_finds_nothing_rather_than_the_nearest_thing():
    assert licence_for(UNRECORDED_DOI) is None
    assert licence_for(None) is None
    assert licence_for("") is None


def test_a_url_keyed_source_is_found_by_its_url():
    assert source_for("https://biophysicalms.org/ccsdatabase") is sources.BUSH_LAB_CCS
    assert source_for("https://BIOPHYSICALMS.org/ccsdatabase") is sources.BUSH_LAB_CCS
    assert source_for(None) is None


def test_a_dataset_is_found_only_by_its_exact_provenance_string():
    assert dataset_for(SEED_2016) is sources.SEED_STRUWE_2016
    assert dataset_for("  " + SEED_2016 + "  ") is sources.SEED_STRUWE_2016


@pytest.mark.parametrize(
    "near_miss",
    [
        "wmxccs seed transcription",
        "wmxccs seed transcription: Struwe 2016",
        "Struwe 2016 Chem Commun ESI Table S1",
        "wmxccs seed transcription: Struwe 2016 Chem Commun ESI Table S1, row 3",
        "WMXCCS SEED TRANSCRIPTION: STRUWE 2016 CHEM COMMUN ESI TABLE S1",
    ],
)
def test_a_provenance_that_merely_resembles_a_registered_dataset_matches_nothing(near_miss):
    # No fuzzy and no substring matching. Loose matching here would let any
    # string that happens to sit inside a registered provenance inherit that
    # dataset's licence.
    assert dataset_for(near_miss) is None


def test_a_provenance_of_none_matches_nothing():
    assert dataset_for(None) is None
    assert dataset_for("") is None


# =========================================================================================
# claim_problem: is a row's own claim backed?
# =========================================================================================


@pytest.mark.parametrize("status", [s for s in ReuseStatus if not can_train_commercial(s)], ids=lambda s: s.value)
def test_a_claim_that_cannot_train_is_never_checked_because_it_errs_in_the_safe_direction(status):
    assert claim_problem(None, status) is None
    assert claim_problem(UNRECORDED_DOI, status) is None
    assert claim_problem(STRUWE_2016_DOI, status) is None


def test_a_claim_that_is_not_a_reuse_status_at_all_is_refused():
    problem = claim_problem(None, "open_attribution")
    assert problem is not None
    assert "not a ReuseStatus" in problem


def test_an_open_claim_on_a_doi_whose_record_says_the_same_thing_is_backed():
    assert claim_problem(STRUWE_2016_DOI, ReuseStatus.OPEN_ATTRIBUTION) is None
    assert claim_problem(STRUWE_2015_DOI, ReuseStatus.OPEN_ATTRIBUTION) is None


def test_an_open_claim_is_backed_whatever_case_the_doi_is_written_in():
    assert claim_problem("10.1039/C6CC06247D", ReuseStatus.OPEN_ATTRIBUTION) is None


def test_an_open_claim_with_no_doi_at_all_is_refused_because_nothing_can_be_checked():
    problem = claim_problem(None, ReuseStatus.OPEN_ATTRIBUTION)
    assert problem is not None
    assert "with no DOI" in problem
    assert "open_attribution" in problem


def test_an_open_claim_on_a_doi_with_no_licence_record_is_refused():
    problem = claim_problem(UNRECORDED_DOI, ReuseStatus.OPEN_ATTRIBUTION)
    assert problem is not None
    assert "has no licence record" in problem
    assert UNRECORDED_DOI in problem


@pytest.mark.parametrize(
    ("doi", "recorded"),
    [
        (STEROID_DOI, "unverified"),
        (BAYESIAN_DOI, "unverified"),
        ("https://ccsbase.net/", "academic_only"),
    ],
)
def test_a_claim_that_exceeds_what_the_record_says_is_refused_and_the_record_governs(doi, recorded):
    # The whole point of the registry: a row may type anything it likes, and the
    # thing a person actually read is what decides.
    problem = claim_problem(doi, ReuseStatus.OPEN_ATTRIBUTION)
    assert problem is not None
    assert "the record governs" in problem
    assert recorded in problem


def test_an_internal_claim_on_a_paper_that_is_on_record_is_refused_as_exceeding_that_record():
    problem = claim_problem(STRUWE_2016_DOI, ReuseStatus.INTERNAL_PROPRIETARY)
    assert problem is not None
    assert "the record governs" in problem


def test_an_internal_measurement_with_no_doi_needs_no_licence_record():
    assert claim_problem(None, ReuseStatus.INTERNAL_PROPRIETARY) is None
    assert claim_problem("", ReuseStatus.INTERNAL_PROPRIETARY) is None


def test_an_internal_claim_that_cites_a_doi_nobody_recorded_is_refused():
    # Citing a publication makes it a claim about that publication, and a claim
    # about a publication needs the page read.
    problem = claim_problem(UNRECORDED_DOI, ReuseStatus.INTERNAL_PROPRIETARY)
    assert problem is not None
    assert "has no licence record" in problem


def test_a_synthetic_fixture_needs_no_licence_record_with_or_without_an_invented_doi():
    assert claim_problem(None, ReuseStatus.SYNTHETIC_FIXTURE) is None
    assert claim_problem(UNRECORDED_DOI, ReuseStatus.SYNTHETIC_FIXTURE) is None


def test_a_record_citing_a_doi_that_is_on_record_may_not_call_itself_synthetic():
    # It is that paper's record, whatever the row says about itself.
    problem = claim_problem(STRUWE_2016_DOI, ReuseStatus.SYNTHETIC_FIXTURE)
    assert problem is not None
    assert "the record governs" in problem


# =========================================================================================
# component_claim_problem: is an analyte identity's claim backed?
# =========================================================================================


def test_an_analyte_drawn_from_a_registered_dataset_is_backed_by_that_dataset():
    assert component_claim_problem(None, SEED_2016, ReuseStatus.OPEN_ATTRIBUTION) is None


def test_an_analyte_claim_that_exceeds_its_datasets_licence_is_refused():
    problem = component_claim_problem(None, SEED_2016, ReuseStatus.INTERNAL_PROPRIETARY)
    assert problem is not None
    assert "the record governs" in problem
    assert SEED_2016 in problem


def test_an_analyte_from_a_registered_dataset_may_not_call_itself_synthetic_either():
    # The dataset's record governs, an internal or a synthetic claim included:
    # the identity came from somewhere, and that somewhere is on record.
    problem = component_claim_problem(None, SEED_2016, ReuseStatus.SYNTHETIC_FIXTURE)
    assert problem is not None
    assert "the record governs" in problem


def test_an_analyte_that_names_the_measurements_own_source_is_backed_by_the_rows_doi():
    assert (
        component_claim_problem(
            STRUWE_2016_DOI,
            "Struwe 2016 Chem Commun",
            ReuseStatus.OPEN_ATTRIBUTION,
            measurement_source="Struwe 2016 Chem Commun",
        )
        is None
    )


def test_an_analyte_naming_the_same_source_is_still_held_to_what_that_doi_record_says():
    problem = component_claim_problem(
        STEROID_DOI,
        "the steroid paper",
        ReuseStatus.OPEN_ATTRIBUTION,
        measurement_source="the steroid paper",
    )
    assert problem is not None
    assert "the record governs" in problem


def test_a_rows_doi_does_not_back_an_analyte_that_names_a_different_source():
    # The laundering case. If a valid DOI backed any identity attached to the
    # row, the DOI would launder an identity of unknown origin into an open
    # licence, which is precisely what the component check exists to stop.
    problem = component_claim_problem(
        STRUWE_2016_DOI,
        UNKNOWN_ORIGIN,
        ReuseStatus.OPEN_ATTRIBUTION,
        measurement_source="Struwe 2016 Chem Commun",
    )
    assert problem is not None
    assert "neither the row's DOI nor a registered dataset" in problem
    assert UNKNOWN_ORIGIN in problem


def test_a_rows_doi_does_not_back_an_analyte_when_the_measurement_names_no_source_at_all():
    problem = component_claim_problem(STRUWE_2016_DOI, UNKNOWN_ORIGIN, ReuseStatus.OPEN_ATTRIBUTION)
    assert problem is not None
    assert "neither the row's DOI nor a registered dataset" in problem


def test_an_analyte_whose_source_is_merely_a_fragment_of_a_registered_dataset_is_not_backed():
    # Exact dataset matching, seen from the gate's side: a provenance that sits
    # inside a registered one must not inherit that dataset's licence.
    problem = component_claim_problem(None, "wmxccs seed transcription", ReuseStatus.OPEN_ATTRIBUTION)
    assert problem is not None
    assert "neither the row's DOI nor a registered dataset" in problem


def test_an_analyte_claim_backed_by_nothing_at_all_is_refused():
    problem = component_claim_problem(None, UNKNOWN_ORIGIN, ReuseStatus.OPEN_ATTRIBUTION)
    assert problem is not None
    assert "neither the row's DOI nor a registered dataset" in problem


def test_an_analyte_with_no_source_recorded_at_all_is_refused_for_an_open_claim():
    problem = component_claim_problem(None, None, ReuseStatus.OPEN_ATTRIBUTION)
    assert problem is not None
    assert "neither the row's DOI nor a registered dataset" in problem


@pytest.mark.parametrize(
    "status",
    [ReuseStatus.INTERNAL_PROPRIETARY, ReuseStatus.SYNTHETIC_FIXTURE],
    ids=lambda s: s.value,
)
def test_an_in_house_or_invented_analyte_identity_needs_no_backing(status):
    assert component_claim_problem(None, UNKNOWN_ORIGIN, status) is None
    assert component_claim_problem(None, None, status) is None


@pytest.mark.parametrize("status", [s for s in ReuseStatus if not can_train_commercial(s)], ids=lambda s: s.value)
def test_an_analyte_claim_that_cannot_train_is_never_checked(status):
    assert component_claim_problem(None, UNKNOWN_ORIGIN, status) is None
    assert component_claim_problem(STRUWE_2016_DOI, SEED_2016, status) is None


def test_an_analyte_claim_that_is_not_a_reuse_status_at_all_is_refused():
    problem = component_claim_problem(None, SEED_2016, "open_attribution")
    assert problem is not None
    assert "not a ReuseStatus" in problem
