"""What a correction is entitled to be quoted as, derived from the data behind it.

THE PROBLEM THIS SOLVES
-----------------------
Every matched ion in this repository comes from one paper. A bias figure fitted on it
describes the difference between the instruments that paper used. It is NOT an
interlaboratory reproducibility figure, and the difference is not pedantic: published
interlaboratory reproducibility for stepped-field DTIMS is 0.29 per cent RSD over many
laboratories, and the numbers this repository produces are single-comparison offsets of
0.1 to 1.1 per cent from one study. Quoting the second as the first would overstate
what is known by a category, not by a margin.

A sentence in a docstring does not prevent that. The precedent this module follows is
the synthetic guard: `declares_synthetic` is a FUNCTION OVER THE RECORDS, so a fixture
cannot pretend to be real data by setting a field, and `assert_quotable` raises at the
point a figure would be published rather than trusting a reader.

HOW THIS IS STRUCTURAL RATHER THAN WRITTEN DOWN
-----------------------------------------------
1. `ScopeStamp` stores PROVENANCE - which studies, which instruments, which platforms -
   and `scope` is a COMPUTED PROPERTY of it. There is no scope field, so there is no
   scope argument, so there is nothing to set. Claiming a wider scope requires naming a
   second study, which is data somebody would have to produce.
2. `ComparisonScope` has no member meaning "interlaboratory reproducibility". The claim
   is not representable, so it cannot be recorded, serialised or returned.
3. `assert_may_be_quoted_as` raises `ScopeExceededError` for any claim the provenance
   does not support. It is called on the publication path.
4. The stamp is REQUIRED on every harmonization output - not optional with a default -
   so a correction cannot be constructed without one.

WHAT THIS DOES NOT DO
---------------------
It cannot stop somebody reading a number out of a report and writing it into a slide
under a different heading. Nothing in a library can. What it stops is this repository
producing, serialising or returning a figure that CLAIMS a scope its data does not
support, and it makes the true scope travel beside every number so that the slide is at
least contradicted by its own source.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Iterable, Sequence


class ComparisonScope(StrEnum):
    """What the values behind a figure span.

    TWO MEMBERS, and the absence of a third is the design. There is no
    INTERLABORATORY_REPRODUCIBILITY, because no corpus this repository has ever held
    could support it: that would need the same ions measured independently by many
    laboratories, with a distribution across them rather than one aggregate value. A
    member that can never be reached honestly is a member somebody will eventually set.
    """

    WITHIN_STUDY = "within_study"
    CROSS_STUDY = "cross_study"


class Claim(StrEnum):
    """What a caller wants to say with a number. Checked against what the data supports."""

    # "these two platforms differed by X on the instruments of one study"
    PLATFORM_DIFFERENCE_WITHIN_A_STUDY = "platform_difference_within_a_study"
    # "these two platforms differ by X across studies"
    PLATFORM_DIFFERENCE_ACROSS_STUDIES = "platform_difference_across_studies"
    # "X is the reproducibility of this platform between laboratories"
    INTERLABORATORY_REPRODUCIBILITY = "interlaboratory_reproducibility"


class ScopeError(Exception):
    """Base for every refusal in this module."""


class ScopeExceededError(ScopeError):
    """A claim wider than the provenance behind it supports."""


SCOPE_TOO_NARROW = (
    "this figure rests on {studies} stud{plural} ({named}) and may not be quoted as {claim!r}."
    " {why}"
)
NEVER_REPRODUCIBILITY = (
    "Interlaboratory reproducibility is a DISTRIBUTION of one platform's values over many laboratories"
    " measuring the same ions independently. This is a single comparison between platforms, so there is no"
    " such distribution here however many studies it spans, and no corpus this repository has held could"
    " supply one. The figure that answers that question is in the methods literature - Stow et al. 2017"
    " report 0.29% RSD for stepped-field DTIMS - and it is not this figure"
)
WITHIN_ONE_STUDY = (
    "Values from a single study describe the instruments that study used. A difference between two"
    " laboratories cannot be separated from a difference between two instruments when there is only one"
    " laboratory in the corpus"
)
NO_PROVENANCE = (
    "a scope stamp was built over no records at all, so there is no provenance to reason about and no claim"
    " is supported. This is refused rather than defaulted, because an empty stamp would otherwise read as"
    " the narrowest scope and pass the narrowest check"
)


def provenance_atom(record: object) -> str:
    """The study a record came from, as one comparable string.

    The DOI where there is one, because two spellings of a citation are one study and a
    DOI is the identifier that says so. The source text otherwise, which is weaker and is
    why a corpus mixing the two can overcount studies rather than undercount them -
    erring towards CROSS_STUDY, which is the direction that grants MORE scope, so
    `scope_is_understated_risk` reports it rather than leaving it silent.
    """
    doi = getattr(record, "doi", None)
    if doi:
        return f"doi:{str(doi).strip().casefold()}"
    source = getattr(record, "source", None)
    return f"source:{str(source or '').strip().casefold()}"


@dataclass(frozen=True)
class ScopeStamp:
    """The provenance behind a figure. Scope is derived from it and is not a field.

    THERE IS DELIBERATELY NO `scope` FIELD. `scope` is a property of `studies`, so
    widening the claim means naming another study - which is data, not a flag. A frozen
    dataclass with a scope field would still have a `scope=` constructor argument, and a
    guard that can be satisfied by a keyword argument is not a guard.

    Build these with `stamp_over`, never by hand: the classmethod derives every field
    from records so that none of them is an assertion.
    """

    studies: tuple[str, ...]
    instruments: tuple[str, ...]
    platforms: tuple[str, ...]
    records_behind_it: int
    """DISTINCT MEASUREMENTS, counted once each however many comparisons used them.

    This is the number the caveat sentence quotes, and it is the conservative one on
    purpose. Until 20 September 2026 the model-level stamp reported 1402 here, arrived at
    by summing eighteen per-stratum counts - so a measurement paired into three strata was
    counted three times and a corpus of 517 measurements was published, to every caller,
    as 1402. Per-stratum the count was always right; only the combination was wrong.
    """
    pairings: int = 0
    """CROSS-PLATFORM PAIRINGS behind the figure, where the builder knows it. 0 = not stated.

    The second number, kept because it answers a different question and the two were being
    conflated. `records_behind_it` is how many measurements exist; this is how many times
    one was compared against another. One matched ion compared in three strata is one
    pairing in each - three pairings over two measurements - and neither number is a
    substitute for the other.
    """

    @property
    def scope(self) -> ComparisonScope:
        """DERIVED. One study is within-study; more than one is cross-study."""
        return ComparisonScope.CROSS_STUDY if len(self.studies) > 1 else ComparisonScope.WITHIN_STUDY

    @property
    def scope_is_understated_risk(self) -> bool:
        """Whether a study here is identified by a source string rather than a DOI.

        Two spellings of one citation would count as two studies and grant CROSS_STUDY
        wrongly, so this is reported beside the scope rather than assumed away.
        """
        return any(not study.startswith("doi:") for study in self.studies)

    @property
    def supports(self) -> frozenset[Claim]:
        """Every claim this provenance supports. Never includes reproducibility."""
        allowed = {Claim.PLATFORM_DIFFERENCE_WITHIN_A_STUDY}
        if self.scope is ComparisonScope.CROSS_STUDY:
            allowed.add(Claim.PLATFORM_DIFFERENCE_ACROSS_STUDIES)
        return frozenset(allowed)

    def caveat(self) -> str:
        """The scope in words, for a report. Travels WITH the number, never instead of it."""
        named = ", ".join(self.studies)
        if self.scope is ComparisonScope.WITHIN_STUDY:
            return (
                f"WITHIN ONE STUDY ({named}), between {len(self.platforms)} platform(s) over"
                f" {self._counted()}. This is a difference between the instruments that"
                f" study used. It is NOT interlaboratory reproducibility."
            )
        return (
            f"ACROSS {len(self.studies)} STUDIES ({named}), between {len(self.platforms)} platform(s) over"
            f" {self._counted()}. Still NOT interlaboratory reproducibility, which needs"
            f" a distribution over laboratories rather than a comparison between platforms."
        )

    def _counted(self) -> str:
        """The size of the evidence, with each number labelled as what it counts.

        The headline is always the DISTINCT measurement count, because the caveat is the
        conservative sentence and the conservative sentence must not overstate. The pairing
        count follows it, named, where the builder knew it - two numbers that were being
        conflated are less misleading stated together than either is alone.
        """
        measurements = f"{self.records_behind_it} distinct measurement(s)"
        if not self.pairings:
            return measurements
        return f"{measurements} entering {self.pairings} cross-platform pairing(s)"


def stamp_over(records: Iterable[object], *, pairings: int = 0) -> ScopeStamp:
    """Derive a stamp from records. The only way one should be built.

    Refuses an empty set: a stamp over nothing would report the narrowest scope and pass
    the narrowest check, which is the wrong direction to fail in.

    `records` MUST already be distinct. This counts what it is handed and cannot tell a
    repeated record from two records, so de-duplication belongs to the caller who knows
    what a duplicate is. `pairings` is passed through because it cannot be derived from
    records at all - it is a property of how they were compared.
    """
    studies: set[str] = set()
    instruments: set[str] = set()
    platforms: set[str] = set()
    count = 0
    for record in records:
        count += 1
        studies.add(provenance_atom(record))
        instrument = getattr(record, "instrument", None)
        if instrument:
            instruments.add(str(instrument))
        platform = getattr(record, "platform", None)
        if platform is None:
            ims = getattr(record, "ims_type", None)
            method = getattr(record, "dtims_method", None)
            if ims is not None:
                platform = f"{ims}/{method}" if method else str(ims)
        if platform:
            platforms.add(str(platform))
    if not count:
        raise ScopeExceededError(NO_PROVENANCE)
    return ScopeStamp(
        studies=tuple(sorted(studies)),
        instruments=tuple(sorted(instruments)),
        platforms=tuple(sorted(platforms)),
        records_behind_it=count,
        pairings=pairings,
    )


def stamp_over_stamps(stamps: Sequence[ScopeStamp], *, records_behind_it: int) -> ScopeStamp:
    """Combine stamps without re-reading records. Unions the provenance; never narrows it.

    `records_behind_it` IS REQUIRED AND IS NOT DERIVED, which is the whole point of the
    keyword. This function cannot see records, so it cannot tell whether two stamps describe
    overlapping sets - and it used to answer anyway, by summing. A measurement paired into
    three strata was therefore counted three times, and the model published 1402 records
    behind a corpus of 517. Summing is right only when the stamps are disjoint, a fact
    nothing here can check and nothing then checked.

    So the count now has to come from whoever holds the records and knows what a duplicate
    is. A caller with no better answer can pass the sum explicitly - but they have to write
    it down, and a number written down is a number somebody chose.
    """
    if not stamps:
        raise ScopeExceededError(NO_PROVENANCE)
    return ScopeStamp(
        studies=tuple(sorted({s for stamp in stamps for s in stamp.studies})),
        instruments=tuple(sorted({i for stamp in stamps for i in stamp.instruments})),
        platforms=tuple(sorted({p for stamp in stamps for p in stamp.platforms})),
        records_behind_it=records_behind_it,
        pairings=sum(stamp.pairings for stamp in stamps),
    )


def assert_may_be_quoted_as(stamp: ScopeStamp, claim: Claim) -> None:
    """Raise unless the provenance supports the claim. Called where a figure is published.

    `INTERLABORATORY_REPRODUCIBILITY` raises for EVERY stamp, whatever it holds, because
    no comparison between platforms is a reproducibility figure for either of them. That
    branch is unreachable-by-data on purpose: it exists so the request has an answer
    rather than a silent success.
    """
    if claim in stamp.supports:
        return
    why = NEVER_REPRODUCIBILITY if claim is Claim.INTERLABORATORY_REPRODUCIBILITY else WITHIN_ONE_STUDY
    raise ScopeExceededError(
        SCOPE_TOO_NARROW.format(
            studies=len(stamp.studies),
            plural="y" if len(stamp.studies) == 1 else "ies",
            named=", ".join(stamp.studies),
            claim=claim.value,
            why=why,
        )
    )
