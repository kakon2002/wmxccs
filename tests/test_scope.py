"""The scope caveat is structural, and these tests are what "structural" means.

The requirement is that a correction fitted on one study cannot be quoted as an
interlaboratory reproducibility figure. A sentence in a docstring does not achieve that.
What is tested here is that the claim is UNREPRESENTABLE and that the scope cannot be
set - only derived - so widening it requires producing data rather than passing a
keyword.
"""

from __future__ import annotations

import dataclasses

import pytest

from wmxccs.scope import (
    Claim,
    ComparisonScope,
    ScopeExceededError,
    ScopeStamp,
    assert_may_be_quoted_as,
    provenance_atom,
    stamp_over,
    stamp_over_stamps,
)


class FakeRecord:
    """A record-shaped object. Only the attributes the scope reads."""

    def __init__(self, doi=None, source="a laboratory", instrument="an instrument", platform="TWIMS"):
        self.doi = doi
        self.source = source
        self.instrument = instrument
        self.platform = platform
        self.ccs = 187.5


ONE_STUDY = [FakeRecord(doi="10.1021/jasms.2c00196") for _ in range(4)]
TWO_STUDIES = ONE_STUDY + [FakeRecord(doi="10.1021/acs.analchem.9b05247")]


# --- the scope is derived, not declared -----------------------------------------------------


def test_scope_is_not_a_field_so_there_is_no_argument_to_set_it_with():
    """The heart of it. A frozen dataclass with a scope field still has a scope= keyword."""
    fields = {field.name for field in dataclasses.fields(ScopeStamp)}
    assert "scope" not in fields, "scope must be derived from studies, never stored"
    with pytest.raises(TypeError):
        ScopeStamp(  # type: ignore[call-arg]
            studies=("doi:a/1",),
            instruments=(),
            platforms=("TWIMS",),
            records_behind_it=1,
            scope=ComparisonScope.CROSS_STUDY,
        )


def test_widening_the_scope_requires_naming_a_second_study():
    """Which is data somebody has to produce, not a flag somebody can flip."""
    one = stamp_over(ONE_STUDY)
    assert one.scope is ComparisonScope.WITHIN_STUDY
    two = stamp_over(TWO_STUDIES)
    assert two.scope is ComparisonScope.CROSS_STUDY
    assert len(two.studies) == 2


def test_the_real_corpus_is_within_study():
    stamp = stamp_over(ONE_STUDY)
    assert stamp.scope is ComparisonScope.WITHIN_STUDY
    assert stamp.studies == ("doi:10.1021/jasms.2c00196",)
    assert stamp.records_behind_it == 4
    assert "NOT interlaboratory reproducibility" in stamp.caveat()


# --- what cannot be claimed ------------------------------------------------------------------


def test_there_is_no_scope_member_for_interlaboratory_reproducibility():
    """Unrepresentable, so it cannot be recorded, serialised or returned."""
    assert {member.value for member in ComparisonScope} == {"within_study", "cross_study"}


def test_a_single_study_may_not_be_quoted_across_studies():
    stamp = stamp_over(ONE_STUDY)
    assert_may_be_quoted_as(stamp, Claim.PLATFORM_DIFFERENCE_WITHIN_A_STUDY)
    with pytest.raises(ScopeExceededError, match="may not be quoted"):
        assert_may_be_quoted_as(stamp, Claim.PLATFORM_DIFFERENCE_ACROSS_STUDIES)


@pytest.mark.parametrize("records", [ONE_STUDY, TWO_STUDIES], ids=["one_study", "two_studies"])
def test_reproducibility_is_refused_however_many_studies_there_are(records):
    """The branch that no amount of data unlocks, and the reason it exists.

    A comparison BETWEEN PLATFORMS is not a reproducibility figure for either of them,
    whatever its provenance: reproducibility is a distribution of one platform's values
    over many laboratories. So this refusal is unreachable-by-data on purpose - it exists
    so the request has an answer rather than a silent success.
    """
    stamp = stamp_over(records)
    with pytest.raises(ScopeExceededError) as refusal:
        assert_may_be_quoted_as(stamp, Claim.INTERLABORATORY_REPRODUCIBILITY)
    assert "DISTRIBUTION" in str(refusal.value)
    assert Claim.INTERLABORATORY_REPRODUCIBILITY not in stamp.supports


def test_supports_never_includes_reproducibility():
    for records in (ONE_STUDY, TWO_STUDIES):
        assert Claim.INTERLABORATORY_REPRODUCIBILITY not in stamp_over(records).supports


# --- failing closed --------------------------------------------------------------------------


def test_a_stamp_over_no_records_is_refused_rather_than_defaulted():
    """An empty stamp would read as the narrowest scope and pass the narrowest check."""
    with pytest.raises(ScopeExceededError, match="no records at all"):
        stamp_over([])
    with pytest.raises(ScopeExceededError):
        stamp_over_stamps([])


def test_a_source_string_rather_than_a_doi_is_flagged_as_possibly_overstating_the_scope():
    """Two spellings of one citation would count as two studies and grant CROSS_STUDY wrongly."""
    spelled_two_ways = [
        FakeRecord(doi=None, source="Feuerstein et al 2022"),
        FakeRecord(doi=None, source="Feuerstein et al., JASMS 2022"),
    ]
    stamp = stamp_over(spelled_two_ways)
    assert stamp.scope is ComparisonScope.CROSS_STUDY
    assert stamp.scope_is_understated_risk is True, "a source-identified study must say the scope may be wrong"
    # and a DOI-identified corpus carries no such risk
    assert stamp_over(ONE_STUDY).scope_is_understated_risk is False


def test_a_doi_identifies_a_study_and_a_citation_spelling_does_not():
    assert provenance_atom(FakeRecord(doi="10.1/A")) == "doi:10.1/a"
    assert provenance_atom(FakeRecord(doi="10.1/a")) == provenance_atom(FakeRecord(doi="10.1/A"))
    assert provenance_atom(FakeRecord(doi=None, source="A Lab")).startswith("source:")


# --- combining ---------------------------------------------------------------------------------


def test_combining_stamps_unions_the_provenance_and_never_narrows_it():
    """Two within-study stamps from different studies combine to cross-study, not to one."""
    a = stamp_over([FakeRecord(doi="10.1/a")])
    b = stamp_over([FakeRecord(doi="10.1/b")])
    assert a.scope is ComparisonScope.WITHIN_STUDY
    assert b.scope is ComparisonScope.WITHIN_STUDY
    combined = stamp_over_stamps([a, b])
    assert combined.scope is ComparisonScope.CROSS_STUDY
    assert combined.records_behind_it == 2
    assert set(combined.studies) == {"doi:10.1/a", "doi:10.1/b"}


def test_combining_stamps_from_one_study_stays_within_study():
    a = stamp_over([FakeRecord(doi="10.1/a")])
    b = stamp_over([FakeRecord(doi="10.1/a")])
    assert stamp_over_stamps([a, b]).scope is ComparisonScope.WITHIN_STUDY
