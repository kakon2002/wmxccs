"""Matched-ion construction: what pairs, what only looks like it pairs, and what may be quoted.

Every statistic in a later milestone is computed over what matching.py calls a
matched set, so both directions of every rule are tested here. A matcher that
pairs everything passes a suite built only from things that ought to pair, which
is why for each rule there is a case that MUST group and a case that MUST NOT.

Two things are load-bearing and easy to lose quietly:

- the corpus these tests run on is invented. It says so - every record carries
  ReuseStatus.SYNTHETIC_FIXTURE - and the point of the synthetic section below is
  that the guards stopping an invented number from being quoted as a result
  actually fire. Those tests prove nothing on their own, so there are also
  NON-synthetic records here (internal, no DOI, so they clear the gate) whose
  report IS quotable. If the synthetic guard were broken to always refuse, or
  always allow, one half or the other goes red.
- the real corpus pairs nothing, and it does so for THREE independent reasons.
  They are asserted separately, so that resolving one does not let the suite
  claim pairs exist on the strength of the other two still holding.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from wmxccs import fixtures
from wmxccs.fixtures import CASES
from wmxccs.identity import DriftGas, SmallMoleculeAnalyte
from wmxccs.licensing import assert_trainable, declares_synthetic
from wmxccs.loader import load_measurements_file
from wmxccs.matching import (
    MatchedIon,
    MatchingReport,
    NotQuotableError,
    assert_quotable,
    build_matched_ions,
    deduplicate,
)
from wmxccs.reuse import ReuseStatus

REPO = Path(__file__).resolve().parents[1]
SEED_FILE_2016 = REPO / "data" / "seed" / "struwe2016_chemcommun.csv"
SEED_FILE_2015 = REPO / "data" / "seed" / "struwe2015_analyst.csv"


def report_for(case: str) -> MatchingReport:
    """The report for one named case in fixtures.CASES, built in isolation."""
    return build_matched_ions(CASES[case]())


def only_group(report: MatchingReport) -> MatchedIon:
    assert len(report.all_groups) == 1, f"expected one ion, got {len(report.all_groups)}"
    return report.all_groups[0]


# --- records that are NOT synthetic -------------------------------------------------------
#
# Internal data with no DOI: trainable, needing no licence record, and declaring
# nothing synthetic anywhere in the record or the analyte inside it. Without
# these the synthetic tests below would be satisfied by a report that simply
# refused everything.

HOUSE = ReuseStatus.INTERNAL_PROPRIETARY
HOUSE_SOURCE = "an in-house run"


def house_record(tag: str, platform: str, ccs: float, **overrides):
    """One real (non-synthetic) measurement: our own data, unpublished, no DOI."""
    analyte = SmallMoleculeAnalyte(
        inchikey=fixtures.synthetic_inchikey(tag), source=HOUSE_SOURCE, reuse_status=HOUSE
    )
    fields = dict(
        analyte=analyte, platform=platform, ccs=ccs, source=HOUSE_SOURCE, doi=None, reuse_status=HOUSE
    )
    fields.update(overrides)
    return fixtures.measurement(**fields)


def house_pair(tag: str = "house") -> tuple:
    """A genuine cross-platform pair of our own, with nothing invented in it."""
    return (house_record(tag, "dtims", 180.0), house_record(tag, "twims", 183.6))


def test_the_non_synthetic_records_this_file_relies_on_really_are_not_synthetic():
    # If these ever started declaring themselves synthetic, every "is quotable"
    # test below would pass for the wrong reason and prove nothing.
    for record in house_pair():
        assert not declares_synthetic(record)
        assert_trainable(record)  # must not raise: internal data needs no licence record


# --- rule 1: a matched set is 2+ measurements from 2+ DIFFERENT platforms ------------------


def test_one_ion_on_three_platforms_is_one_set_of_three_and_not_three_pairs():
    report = report_for("a three-way match across DTIMS, TWIMS and TIMS")
    assert len(report.matched) == 1
    group = report.matched[0]
    assert len(group.platforms) == 3
    assert len(group.measurements) == 3
    assert group.is_cross_platform
    # Flattening a three-way set into pairs would count one compound three times.
    assert len(group.platform_pairs) == 3
    assert set(group.platform_pairs) == {
        "DTIMS/stepped_field vs TIMS",
        "DTIMS/stepped_field vs TWIMS",
        "TIMS vs TWIMS",
    }
    assert report.widest_match == 3
    assert report.platform_pair_counts == {pair: 1 for pair in group.platform_pairs}


def test_one_ion_on_two_platforms_is_a_matched_set():
    report = report_for("a two-way match")
    assert len(report.matched) == 1
    assert len(report.matched[0].platforms) == 2
    assert report.single_platform == ()
    assert report.widest_match == 2


def test_two_values_of_one_ion_on_one_platform_are_replicates_and_do_not_match():
    # The negative half of rule 1: two measurements, one platform. A set that
    # counted measurements rather than platforms would call this a match.
    report = report_for("replicates on one platform")
    assert report.matched == ()
    assert len(report.single_platform) == 1
    group = report.single_platform[0]
    assert len(group.measurements) == 2
    assert len(group.platforms) == 1
    assert not group.is_cross_platform
    assert group.platform_pairs == ()
    assert report.widest_match == 0
    assert report.platform_pair_counts == {}


def test_two_replicates_on_one_platform_are_still_two_measurements_not_a_duplicate():
    report = report_for("replicates on one platform")
    assert report.duplicates_collapsed == 0
    assert report.measurements == 2


# --- rule 2: deduplication ----------------------------------------------------------------


def test_one_measurement_republished_under_two_dois_is_one_measurement_not_a_pair():
    report = report_for("one measurement republished in a review")
    assert report.records_in == 2
    assert report.duplicates_collapsed == 1
    assert report.measurements == 1
    assert report.matched == ()
    assert len(report.single_platform) == 1


def test_the_second_doi_of_a_republished_measurement_is_kept_in_its_provenance():
    # Dropping the citation would lose the fact that two papers carry the number.
    records = CASES["one measurement republished in a review"]()
    measurements, collapsed = deduplicate(records)
    assert collapsed == 1
    (kept,) = measurements
    assert kept.republished
    assert kept.also_published_as == (("a review reprinting it", fixtures.DOI_REVIEW),)
    assert kept.sources == ("fixture laboratory A", "a review reprinting it")
    assert set(kept.dois) == {fixtures.DOI_A, fixtures.DOI_REVIEW}


def test_deduplication_keeps_the_first_record_itself_rather_than_merging_a_new_one():
    # Constraint 3: nothing about an original measurement is rewritten here.
    records = CASES["one measurement republished in a review"]()
    (kept,), _collapsed = deduplicate(records)
    assert kept.record is records[0]
    assert kept.record.ccs == records[0].ccs


def test_a_republished_measurement_carries_both_dois_up_into_the_group():
    report = report_for("one measurement republished in a review")
    group = only_group(report)
    assert set(group.dois) == {fixtures.DOI_A, fixtures.DOI_REVIEW}
    assert set(group.sources) == {"fixture laboratory A", "a review reprinting it"}


def test_two_platforms_reporting_the_identical_value_are_a_match_and_not_a_duplicate():
    # The counterweight: collapsing on the value alone would delete exactly the
    # evidence this platform exists to find.
    report = report_for("two laboratories agreeing exactly")
    assert report.duplicates_collapsed == 0
    assert report.measurements == 2
    assert len(report.matched) == 1
    assert len(report.matched[0].platforms) == 2


def test_deduplication_does_not_key_on_the_compound_name():
    # Two different molecules, one printed name, one identical number, one
    # platform. Nothing here may collapse: the name is not identity.
    def named(inchikey: str):
        analyte = SmallMoleculeAnalyte(
            inchikey=inchikey,
            display_name="compound 7",
            source="a fixture",
            reuse_status=ReuseStatus.SYNTHETIC_FIXTURE,
        )
        return fixtures.measurement(analyte=analyte, platform="dtims", ccs=180.0)

    records = (named(fixtures.synthetic_inchikey("nameone")), named(fixtures.synthetic_inchikey("nametwo")))
    assert records[0].analyte.display_name == records[1].analyte.display_name
    assert records[0].ccs == records[1].ccs

    report = build_matched_ions(records)
    assert report.duplicates_collapsed == 0
    assert report.measurements == 2
    assert len(report.all_groups) == 2
    assert report.matched == ()


def test_two_different_molecules_on_two_platforms_do_not_match():
    report = report_for("two different molecules")
    assert report.matched == ()
    assert len(report.single_platform) == 2


def test_two_glycan_isomers_of_one_composition_do_not_match():
    # LNH and LNnH: one composition, two structures, two cross sections.
    report = report_for("two glycan isomers of one composition")
    assert report.matched == ()
    assert len(report.single_platform) == 2


# --- rule 3: conformers -------------------------------------------------------------------


def test_two_conformers_on_one_platform_stay_two_groups_and_are_neither_duplicates_nor_a_match():
    report = report_for("two conformers on one platform")
    assert report.duplicates_collapsed == 0
    assert report.measurements == 2
    assert report.matched == ()
    assert len(report.single_platform) == 2
    assert {group.conformer for group in report.single_platform} == {1, 2}
    assert len(report.conformer_groups) == 2


def test_the_conformer_index_is_part_of_the_grouping_key_in_its_own_right():
    # Two conformers reported at the SAME rounded value. If the index left the
    # grouping key these would become one set; the differing value cannot be what
    # holds them apart, because here there is no differing value.
    def peak(index: int):
        return fixtures.measurement(
            analyte=fixtures.small_molecule(fixtures.synthetic_inchikey("samevalue")),
            platform="dtims",
            ccs=180.0,
            conformer=index,
            conformers_total=2,
        )

    report = build_matched_ions((peak(1), peak(2)))
    assert len(report.all_groups) == 2
    assert {group.conformer for group in report.all_groups} == {1, 2}


def test_the_conformer_index_is_part_of_the_duplicate_key_in_its_own_right():
    # The same two peaks through deduplication alone: same key, same platform,
    # same value, different conformer is two measurements, not one republished.
    def peak(index: int):
        return fixtures.measurement(
            analyte=fixtures.small_molecule(fixtures.synthetic_inchikey("samevalue")),
            platform="dtims",
            ccs=180.0,
            conformer=index,
            conformers_total=2,
        )

    measurements, collapsed = deduplicate((peak(1), peak(2)))
    assert collapsed == 0
    assert len(measurements) == 2


def test_the_same_conformer_index_on_two_platforms_groups_but_is_refused_as_source_local():
    report = report_for("the same conformer index on two platforms")
    assert len(report.matched) == 1
    group = report.matched[0]
    assert group.conformer == 1
    assert len(group.platforms) == 2
    blockers = " | ".join(group.blockers)
    assert "conformer numbering is source-local" in blockers
    assert "conformer 1" in blockers
    assert not group.usable


def test_a_cross_platform_set_with_no_conformer_index_carries_no_conformer_blocker():
    # The negative direction of the conformer blocker: it must fire on the
    # conformer-resolved set and only on it.
    group = report_for("a two-way match").matched[0]
    assert group.conformer is None
    assert not any("source-local" in blocker for blocker in group.blockers)


# --- rule 4: an unstated charge carrier can never match ------------------------------------


def test_two_ions_with_an_unstated_charge_carrier_never_match_even_each_other():
    report = report_for("two unstated charge carriers")
    assert report.matched == ()
    assert report.single_platform == ()
    assert len(report.unmatchable) == 2
    for group in report.unmatchable:
        assert not group.is_cross_platform
        assert len(group.measurements) == 1
        assert not group.key.matchable


def test_an_unstated_carrier_does_not_match_the_same_ion_with_a_named_carrier():
    # Same protein, same charge, same platforms: one row names its carrier and
    # one does not. Knowing it is 24 protons is not the same as not knowing.
    unstated = CASES["two unstated charge carriers"]()[0]
    named = unstated.model_copy(update={"adduct": "[M+24H]24+"})
    other_platform = named.model_copy(update={"source": "fixture laboratory B"})

    report = build_matched_ions((unstated, named, other_platform))
    assert len(report.unmatchable) == 1
    assert report.matched == ()
    assert unstated.matched_ion_key != named.matched_ion_key


def test_two_named_carrier_records_of_that_same_ion_do_match_so_the_carrier_is_what_blocks_it():
    # The control for the test above: with the carrier stated, the pair forms.
    unstated = CASES["two unstated charge carriers"]()[0]
    named = unstated.model_copy(update={"adduct": "[M+24H]24+"})
    partner = CASES["two unstated charge carriers"]()[1].model_copy(update={"adduct": "[M+24H]24+"})

    report = build_matched_ions((named, partner))
    assert report.unmatchable == ()
    assert len(report.matched) == 1
    assert len(report.matched[0].platforms) == 2


# --- the biopharmaceutical and gas cases: same-looking records that are different ions -----


def test_a_native_and_a_denatured_antibody_do_not_match():
    report = report_for("a native and a denatured antibody")
    assert report.matched == ()
    assert len(report.single_platform) == 2


def test_two_adcs_differing_only_in_dar_do_not_match():
    report = report_for("two ADCs differing only in DAR")
    assert report.matched == ()
    assert len(report.single_platform) == 2


def test_a_glycopeptide_and_its_bare_backbone_do_not_match():
    report = report_for("a glycopeptide and its bare backbone")
    assert report.matched == ()
    assert len(report.single_platform) == 2


def test_one_ion_referred_to_two_gases_does_not_match():
    # Constraint 4: values referred to different gases never pool.
    report = report_for("one ion referred to two gases")
    assert report.matched == ()
    assert len(report.single_platform) == 2


def test_a_cyclic_multipass_value_still_matches_a_dtims_value_of_the_same_ion():
    # Pass count belongs to the calibration group, not to the ion. It must not
    # reach the matched-ion key, or a cyclic value could never pair with anything.
    report = report_for("a cyclic pair at different pass counts")
    assert len(report.matched) == 1
    assert set(report.matched[0].platforms) == {"CYCLIC", "DTIMS/stepped_field"}


# --- rule 5: provenance and licence --------------------------------------------------------


def test_a_matched_set_exposes_its_sources_dois_and_reuse_statuses():
    group = report_for("a three-way match across DTIMS, TWIMS and TIMS").matched[0]
    assert set(group.sources) == {
        "fixture laboratory A",
        "fixture laboratory B",
        "fixture laboratory C",
    }
    assert set(group.dois) == {fixtures.DOI_A, fixtures.DOI_B, "10.9999/fixture.c"}
    assert group.reuse_statuses == {ReuseStatus.SYNTHETIC_FIXTURE.value: 3}


def test_a_matched_set_with_an_unverified_member_is_refused_and_the_refusal_names_that_member():
    report = report_for("a match with an unusable member")
    assert len(report.matched) == 1
    group = report.matched[0]
    assert not group.usable
    assert report.usable == ()
    blockers = " | ".join(group.blockers)
    assert "a source nobody has checked" in blockers
    assert ReuseStatus.UNVERIFIED.value in blockers
    assert "least usable member" in blockers
    # and the member that is fine is not the one blamed
    assert "fixture laboratory A" not in blockers
    assert group.reuse_statuses[ReuseStatus.UNVERIFIED.value] == 1


def test_a_matched_set_with_an_unreadable_reuse_status_is_refused_rather_than_assumed():
    # Default deny: a status nobody can parse is not a status anybody cleared.
    good, other = CASES["a two-way match"]()
    unreadable = other.model_copy(update={"reuse_status": "not a status at all"})
    group = build_matched_ions((good, unreadable)).matched[0]
    blockers = " | ".join(group.blockers)
    assert "unreadable reuse status" in blockers
    assert "fixture laboratory B" in blockers
    assert not group.usable
    assert group.reuse_statuses.get("unreadable") == 1


def test_a_matched_set_whose_members_are_all_fine_carries_only_the_synthetic_blocker():
    group = report_for("a two-way match").matched[0]
    assert len(group.blockers) == 1
    assert "this set is synthetic" in group.blockers[0]
    assert not group.usable


def test_a_matched_set_of_real_records_that_are_all_fine_is_usable():
    # The other direction of the same rule, without the synthetic blocker in the
    # way: a real pair with clean licences enters a fit.
    report = build_matched_ions(house_pair())
    assert len(report.matched) == 1
    assert report.matched[0].blockers == ()
    assert report.matched[0].usable
    assert report.usable == report.matched


def test_a_blocked_matched_set_is_still_reported_rather_than_dropped():
    # It IS the same ion. Dropping it would report a smaller corpus with nothing
    # to say a source is waiting on a licence check.
    report = report_for("a match with an unusable member")
    assert len(report.matched) == 1
    assert report.usable == ()
    assert "may not enter a fit" in report.summary()
    assert "a source nobody has checked" in report.summary()


def test_a_real_matched_set_with_one_unverified_member_is_refused_and_named():
    # The licence rule tested clear of the synthetic blocker, so that removing
    # the synthetic guard cannot make this one pass by accident.
    good, bad = house_pair("housemixed")
    bad = bad.model_copy(update={"reuse_status": ReuseStatus.UNVERIFIED, "source": "a source nobody read"})
    report = build_matched_ions((good, bad))
    group = report.matched[0]
    assert not group.usable
    assert report.usable == ()
    assert "a source nobody read" in " | ".join(group.blockers)
    assert report.quotable  # not synthetic; unusable and unquotable are different questions


# --- rule 6: synthetic safety ---------------------------------------------------------------


@pytest.mark.parametrize("case", sorted(CASES))
def test_every_record_of_every_fixture_case_declares_itself_synthetic(case: str):
    # Asked of the WHOLE record, because one case deliberately carries an
    # unverified status on the measurement to exercise the licence rule. Its
    # analyte still declares the fixture, which is what keeps that case out of a
    # quotable report; a case declaring nothing anywhere would be a real record
    # in a file of invented ones.
    records = CASES[case]()
    assert records, f"case {case!r} built no records"
    for record in records:
        assert declares_synthetic(record)
        statuses = {record.reuse_status, *(part.reuse_status for part in record.component_records())}
        assert ReuseStatus.SYNTHETIC_FIXTURE in statuses


def test_every_group_built_from_the_fixture_corpus_declares_synthetic():
    report = build_matched_ions(fixtures.corpus())
    assert report.all_groups
    for group in report.all_groups:
        assert group.declares_synthetic
        assert group.synthetic_measurements == len(group.measurements)


def test_the_fixture_corpus_report_is_not_quotable_and_says_why():
    report = build_matched_ions(fixtures.corpus())
    assert report.quotable is False
    refusal = report.refusal()
    assert refusal is not None
    assert "synthetic" in refusal
    assert "may not be quoted as a result" in refusal
    assert str(report.measurements) in refusal


def test_assert_quotable_refuses_a_synthetic_report():
    report = build_matched_ions(fixtures.corpus())
    with pytest.raises(NotQuotableError, match="may not be quoted as a result"):
        assert_quotable(report)


def test_not_quotable_error_is_not_a_value_error():
    # A caller catching ValueError around a numeric routine must not swallow this.
    assert issubclass(NotQuotableError, Exception)
    assert not issubclass(NotQuotableError, ValueError)


def test_a_report_over_real_records_is_quotable_and_assert_quotable_returns_cleanly():
    # Without this the synthetic tests prove nothing: a gate that refused every
    # report would satisfy them all.
    report = build_matched_ions(house_pair())
    assert report.synthetic_groups == ()
    assert report.quotable is True
    assert report.refusal() is None
    assert assert_quotable(report) is None


def test_a_synthetic_declaration_on_the_analyte_alone_still_makes_the_set_synthetic():
    # The status can sit on a component. Counting only the outermost layer would
    # present an invented identity as real data.
    analyte = fixtures.small_molecule(fixtures.synthetic_inchikey("insidesynth"))
    assert analyte.reuse_status is ReuseStatus.SYNTHETIC_FIXTURE
    record = fixtures.measurement(
        analyte=analyte, platform="dtims", ccs=180.0, source=HOUSE_SOURCE, doi=None, reuse_status=HOUSE
    )
    assert record.reuse_status is HOUSE

    report = build_matched_ions((record,))
    group = only_group(report)
    assert group.synthetic_measurements == 1
    assert group.declares_synthetic
    assert report.quotable is False
    with pytest.raises(NotQuotableError):
        assert_quotable(report)


def test_a_report_is_not_quotable_when_only_a_single_platform_group_is_synthetic():
    # quotable asks about every group, not only the matched ones. A report mixing
    # real and invented records cannot be quoted without saying which is which.
    lonely = fixtures.measurement(
        analyte=fixtures.small_molecule(fixtures.synthetic_inchikey("lonely")), platform="tims", ccs=200.0
    )
    report = build_matched_ions((*house_pair(), lonely))
    assert len(report.matched) == 1
    assert not report.matched[0].declares_synthetic
    assert len(report.single_platform) == 1
    assert report.single_platform[0].declares_synthetic
    assert report.quotable is False
    with pytest.raises(NotQuotableError):
        assert_quotable(report)


def test_a_report_is_not_quotable_when_only_an_unmatchable_group_is_synthetic():
    # The same point from the other side of the report: an unmatchable group is
    # still a group, and a synthetic one still poisons the whole report.
    report = build_matched_ions((*house_pair(), *CASES["two unstated charge carriers"]()))
    assert len(report.matched) == 1
    assert not report.matched[0].declares_synthetic
    assert len(report.unmatchable) == 2
    assert all(group.declares_synthetic for group in report.unmatchable)
    assert report.quotable is False
    assert report.refusal() is not None


def test_a_synthetic_matched_set_is_never_counted_as_usable():
    report = build_matched_ions(fixtures.corpus())
    assert report.matched  # there ARE matched sets in the fixture corpus
    assert report.usable == ()
    for group in report.matched:
        assert not group.usable
        assert any("may not be quoted as a result" in blocker for blocker in group.blockers)


def test_the_summary_of_a_synthetic_report_says_so_in_words():
    summary = build_matched_ions(fixtures.corpus()).summary()
    assert "SYNTHETIC" in summary
    assert "may not be quoted as a result" in summary


# --- conservation ---------------------------------------------------------------------------


def test_records_in_equals_measurements_plus_duplicates_collapsed_on_a_corpus_with_duplicates():
    corpus = fixtures.corpus()
    report = build_matched_ions(corpus)
    assert report.duplicates_collapsed > 0, "this law is only worth testing where something collapsed"
    assert report.records_in == len(corpus)
    assert report.records_in == report.measurements + report.duplicates_collapsed


def test_every_measurement_lands_in_exactly_one_group():
    report = build_matched_ions(fixtures.corpus())
    in_groups = sum(len(group.measurements) for group in report.all_groups)
    assert in_groups == report.measurements
    assert len(report.all_groups) == len(report.matched) + len(report.single_platform) + len(
        report.unmatchable
    )


def test_an_empty_corpus_reports_nothing_rather_than_failing():
    report = build_matched_ions(())
    assert report.records_in == 0
    assert report.measurements == 0
    assert report.all_groups == ()
    assert report.widest_match == 0
    assert report.quotable is True  # nothing invented, because nothing at all
    assert report.refusal() is None


# --- the real corpus: zero matched sets, for three independent reasons -----------------------


def real_corpus() -> tuple:
    """Every record the two seed files hold, cleared or not. 117 of them."""
    held = []
    for path in (SEED_FILE_2016, SEED_FILE_2015):
        held.extend(load_measurements_file(path).records)
    return tuple(held)


def test_the_real_corpus_holds_no_matched_sets_at_all():
    corpus = real_corpus()
    assert len(corpus) == 117
    report = build_matched_ions(corpus)
    assert report.records_in == 117
    assert report.matched == ()
    assert report.usable == ()
    assert report.widest_match == 0
    assert report.platform_pair_counts == {}
    # and nothing was lost on the way to that answer
    assert report.records_in == report.measurements + report.duplicates_collapsed


def test_the_real_corpus_pairs_nothing_because_only_one_platform_is_represented():
    platforms = {(record.ims_type, record.dtims_method) for record in real_corpus()}
    assert len(platforms) == 1, f"more than one platform is now held: {platforms}"


def test_the_real_corpus_pairs_nothing_because_the_two_files_share_no_identity_atom():
    of_2016 = {record.matched_ion_key.analyte for record in load_measurements_file(SEED_FILE_2016).records}
    of_2015 = {record.matched_ion_key.analyte for record in load_measurements_file(SEED_FILE_2015).records}
    assert of_2016 and of_2015
    assert of_2016.isdisjoint(of_2015), f"the files now share an analyte: {of_2016 & of_2015}"


def test_the_real_corpus_pairs_nothing_because_the_two_files_refer_to_different_gases():
    of_2016 = {record.drift_gas for record in load_measurements_file(SEED_FILE_2016).records}
    of_2015 = {record.drift_gas for record in load_measurements_file(SEED_FILE_2015).records}
    assert of_2016 == {DriftGas.HE}
    assert of_2015 == {DriftGas.UNSTATED}
    assert of_2016.isdisjoint(of_2015)


def test_no_real_record_is_ever_reported_as_synthetic():
    # The gate has to run in both directions: if declares_synthetic said yes to a
    # real record, the whole seeded corpus would become unquotable in silence.
    report = build_matched_ions(real_corpus())
    assert report.synthetic_groups == ()
    assert report.quotable is True
    assert assert_quotable(report) is None


# --- two decisions pinned after review, both flagged during M2 ---------------------------


def test_stepped_field_and_single_field_dtims_count_as_two_platforms():
    """A primary value and a calibrated one are a genuine cross-platform comparison.

    Stepped-field DTIMS follows from first principles; single-field DTIMS reads
    off a calibration curve and inherits its reference set's bias. Comparing them
    measures exactly what this pipeline exists to measure, and the published
    reproducibility figures treat them as separate numbers. Folding them together
    would also hide the primary-versus-derived distinction from any later fit,
    whose DTIMS arm would then be regressing against a moving anchor.
    """
    report = build_matched_ions(fixtures.a_primary_and_a_calibrated_dtims_pair())
    assert len(report.matched) == 1
    assert report.matched[0].platforms == ("DTIMS/single_field", "DTIMS/stepped_field")
    assert report.matched[0].is_cross_platform


def test_two_dtims_values_of_the_same_method_are_still_only_one_platform():
    """The negative direction. Without it the test above passes on a matcher that
    calls every record its own platform, which would make replicates matches."""
    report = build_matched_ions(fixtures.replicates_on_one_platform())
    assert report.matched == ()
    assert len(report.single_platform) == 1


def test_two_conformers_reported_at_the_same_value_are_still_two_measurements():
    """The case where the conformer index in the deduplication key is load-bearing.

    The ordinary conformer pair carries two different CCS values, so the value
    alone keeps those rows apart and the index never has to do any work. Here the
    two peaks share a value, and the index is the only thing standing between a
    legitimate pair and being collapsed into one measurement.
    """
    records = fixtures.two_conformers_sharing_one_value()
    assert records[0].ccs == records[1].ccs
    report = build_matched_ions(records)
    assert report.duplicates_collapsed == 0
    assert report.measurements == 2
    assert len(report.single_platform) == 2
    assert report.records_in == report.measurements + report.duplicates_collapsed
