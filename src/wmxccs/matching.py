"""Matched-ion construction: which measurements are of the same ion, and which are not.

THE MODULE EVERYTHING DOWNSTREAM DEPENDS ON. Every statistic, every bias figure
and every harmonized value in later milestones is computed over what this module
decides is a matched set. A pairing that is too loose invents agreement between
things that were never the same ion; one that is too tight reports a clean run
over pairs it failed to find. Both are silent.

BUILT AGAINST SYNTHETIC FIXTURES, ON PURPOSE.
There is no real cross-platform data in this repository and there may not be for
some time: the only identified source is the steroid interplatform study, whose
licence has not been read. Waiting for that answer before writing the pairing
logic would put the deadline at risk for no benefit, because the logic does not
depend on which dataset fills it. So the rules below are exercised by synthetic
records that declare themselves synthetic, and there are guards - tested ones -
that make it impossible to quote a synthetic result as a real one.

WHAT A MATCHED SET IS
---------------------
Two or more measurements that share a matched-ion key and come from MORE THAN ONE
PLATFORM. Three-way and wider sets are first-class, not pairs with extras: the
likely benchmark carries DTIMS, TWIMS and TIMS for one compound, and flattening
that into three pairs would count one compound three times and weight it as
though three laboratories had agreed independently.

FIVE WAYS TWO RECORDS CAN LOOK LIKE A MATCH AND NOT BE ONE
----------------------------------------------------------
1. THE SAME MEASUREMENT, TWICE. A value republished in a review, or transcribed
   twice, is one measurement. It is collapsed before matching, and every DOI it
   appeared under is kept. Collapsing is on the value and its provenance, never
   on the compound name: two laboratories measuring one compound to the same
   number are a matched pair, and treating them as a duplicate would delete the
   only kind of evidence this platform exists to find.
2. TWO CONFORMERS. One structure giving two arrival-time peaks under one set of
   conditions is a legitimate pair from M0 - not a duplicate, and not a
   cross-platform match either. They are two values of one ion from one
   instrument.
3. TWO IONS WITH AN UNSTATED CHARGE CARRIER. Neither is known to be the other.
   The key already makes this structural; this module asserts it survives.
4. ONE PLATFORM, TWICE. Two DTIMS values of one ion are replicates. Useful, and
   not a cross-platform comparison.
5. A SET NOBODY MAY USE. A matched set is only as usable as its least usable
   member, and the refusal has to name WHICH member, because the remedy is to go
   and read that source's licence.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

from .identity import MatchedIonKey
from .licensing import declares_synthetic
from .reuse import ReuseStatus, as_reuse_status, can_train_commercial

# Why a matched set may not enter a fit.
SYNTHETIC_SET = (
    "this set is synthetic: {synthetic} of {held} measurements are fixtures declared in code. No number"
    " from it describes real data, and it may not be quoted as a result"
)
UNUSABLE_MEMBER = (
    "the measurement from {source!r} ({doi}) has reuse status '{status}', which is not cleared for use;"
    " a matched set is only as usable as its least usable member"
)
CONFORMER_INDICES_UNCONFIRMED = (
    "this set pairs conformer {index} across platforms, and conformer numbering is source-local: that"
    " conformer {index} in one source is conformer {index} in another has not been established by"
    " anything. A person has to confirm the two peaks are the same conformer before this pair is used"
)


def _platform_of(record: object) -> str:
    """What counts as a platform for the purpose of pairing.

    STEPPED-FIELD AND SINGLE-FIELD DTIMS COUNT AS TWO DIFFERENT PLATFORMS, and
    that is deliberate rather than a side effect of reusing the calibration
    group's label. One ion measured both ways is a matched set here.

    The reason is that the two are not the same measurement in the way that
    matters to this platform. Stepped-field DTIMS is primary: the cross section
    follows from first principles and no reference value enters it. Single-field
    DTIMS is calibrated, like TWIMS and TIMS, and inherits whatever bias its
    reference set carries. Comparing them measures exactly the thing this
    pipeline exists to measure, and the literature already treats them apart -
    the interlaboratory reproducibility figures quoted in CONTEXT.md are
    separate numbers for the two, 0.29 per cent RSD stepped-field against 0.54
    per cent mean absolute bias single-field.

    Folding them together would lose a real comparison and would also make the
    primary-versus-derived distinction invisible at exactly the point where a
    harmonization model most needs it: a fit whose "DTIMS" arm silently mixed
    primary and calibrated values would be regressing against a moving anchor.

    This is a wider reading of "two platforms" than the DTIMS/TWIMS/TIMS framing
    in the module docstring, so it is stated here, tested, and carries a mutation
    of its own rather than being left for somebody to discover in M3.
    """
    ims = getattr(record, "ims_type", None)
    method = getattr(record, "dtims_method", None)
    return str(ims) if method is None else f"{ims}/{method}"


def _doi_of(record: object) -> str | None:
    doi = getattr(record, "doi", None)
    return str(doi).strip().casefold() if doi else None


def _source_of(record: object) -> str:
    return str(getattr(record, "source", "") or "(no source recorded)")


def _status_of(record: object) -> ReuseStatus | None:
    try:
        return as_reuse_status(getattr(record, "reuse_status", None))
    except (TypeError, ValueError):
        return None


@dataclass(frozen=True)
class Measurement:
    """One measurement after deduplication, with every source it was published under.

    `record` is the first record seen for this measurement; `also_published_as`
    holds the provenance of the rest. Nothing is merged and no value is changed:
    a duplicate is dropped from the COUNT, and its citation is kept, because the
    fact that two papers carry the number is worth knowing and is not evidence
    that two measurements were made.
    """

    record: object
    also_published_as: tuple[tuple[str, str | None], ...] = ()  # (source, doi)

    @property
    def platform(self) -> str:
        return _platform_of(self.record)

    @property
    def sources(self) -> tuple[str, ...]:
        return (_source_of(self.record), *(source for source, _doi in self.also_published_as))

    @property
    def dois(self) -> tuple[str, ...]:
        found = [_doi_of(self.record), *(doi for _source, doi in self.also_published_as)]
        return tuple(doi for doi in found if doi)

    @property
    def republished(self) -> bool:
        return bool(self.also_published_as)


def _duplicate_key(record: object) -> tuple:
    """What makes two records the SAME MEASUREMENT rather than two of them.

    The matched-ion key, the platform, the value, and the conformer index.

    - the VALUE, because two laboratories do not independently produce the same
      CCS to the last digit. Same key, same platform, same number is one
      measurement quoted twice far more often than it is a coincidence.
    - the PLATFORM, because a DTIMS and a TWIMS value that happen to agree
      exactly are the matched pair this platform exists to find, and collapsing
      them would delete it.
    - the CONFORMER INDEX, because two conformers of one ion carry different
      values anyway, and including it makes it impossible for a future change to
      the value comparison to collapse them by accident.

    The DOI is deliberately NOT here. Two records of one measurement under two
    DOIs is exactly the republication case, and keying on the DOI would keep them
    apart - which is the failure this check exists to prevent. The compound name
    is not here either, and never is.
    """
    return (
        getattr(record, "matched_ion_key", None),
        _platform_of(record),
        getattr(record, "ccs", None),
        getattr(record, "conformer", None),
    )


def deduplicate(records: Iterable[object]) -> tuple[tuple[Measurement, ...], int]:
    """Collapse republications. Returns the measurements and how many records were collapsed."""
    seen: dict[tuple, list[object]] = {}
    order: list[tuple] = []
    for record in records:
        key = _duplicate_key(record)
        if key not in seen:
            seen[key] = []
            order.append(key)
        seen[key].append(record)

    measurements: list[Measurement] = []
    collapsed = 0
    for key in order:
        group = seen[key]
        first, rest = group[0], group[1:]
        collapsed += len(rest)
        measurements.append(
            Measurement(
                record=first,
                also_published_as=tuple((_source_of(other), _doi_of(other)) for other in rest),
            )
        )
    return tuple(measurements), collapsed


@dataclass(frozen=True)
class MatchedIon:
    """Measurements of one ion, and everything needed to decide whether they may be used."""

    key: MatchedIonKey
    conformer: int | None
    measurements: tuple[Measurement, ...]

    @property
    def platforms(self) -> tuple[str, ...]:
        return tuple(sorted({measurement.platform for measurement in self.measurements}))

    @property
    def is_cross_platform(self) -> bool:
        """Two or more PLATFORMS, not two or more measurements. Replicates are not a match."""
        return len(self.platforms) >= 2

    @property
    def platform_pairs(self) -> tuple[str, ...]:
        found = self.platforms
        return tuple(f"{a} vs {b}" for i, a in enumerate(found) for b in found[i + 1 :])

    @property
    def sources(self) -> tuple[str, ...]:
        found: list[str] = []
        for measurement in self.measurements:
            for source in measurement.sources:
                if source not in found:
                    found.append(source)
        return tuple(found)

    @property
    def dois(self) -> tuple[str, ...]:
        found: list[str] = []
        for measurement in self.measurements:
            for doi in measurement.dois:
                if doi not in found:
                    found.append(doi)
        return tuple(found)

    @property
    def reuse_statuses(self) -> Mapping[str, int]:
        return dict(Counter(str(_status_of(m.record) or "unreadable") for m in self.measurements))

    @property
    def synthetic_measurements(self) -> int:
        """Measurements declaring themselves invented, asked of the whole record."""
        return sum(1 for m in self.measurements if declares_synthetic(m.record))

    @property
    def declares_synthetic(self) -> bool:
        return self.synthetic_measurements > 0

    @property
    def blockers(self) -> tuple[str, ...]:
        """Every reason this set may not enter a fit. Each names the member at fault."""
        problems: list[str] = []
        if self.declares_synthetic:
            problems.append(
                SYNTHETIC_SET.format(synthetic=self.synthetic_measurements, held=len(self.measurements))
            )
        for measurement in self.measurements:
            status = _status_of(measurement.record)
            if status is None:
                problems.append(
                    f"the measurement from {_source_of(measurement.record)!r} has an unreadable reuse"
                    " status, so it is refused rather than assumed"
                )
            elif not can_train_commercial(status):
                problems.append(
                    UNUSABLE_MEMBER.format(
                        source=_source_of(measurement.record),
                        doi=_doi_of(measurement.record) or "no DOI",
                        status=status.value,
                    )
                )
        if self.conformer is not None and self.is_cross_platform:
            problems.append(CONFORMER_INDICES_UNCONFIRMED.format(index=self.conformer))
        return tuple(problems)

    @property
    def usable(self) -> bool:
        return not self.blockers

    def __str__(self) -> str:
        where = ", ".join(self.platforms)
        conformer = "" if self.conformer is None else f" conformer {self.conformer}"
        return f"{self.key}{conformer}  [{len(self.measurements)} measurements across {where}]"


@dataclass(frozen=True)
class MatchingReport:
    """What pairing found, what it refused to pair, and why.

    Conservation holds and is tested:
        records_in == measurements + duplicates_collapsed
        measurements == sum of the members of every group, matched and unmatched
    """

    records_in: int
    duplicates_collapsed: int
    matched: tuple[MatchedIon, ...] = ()  # cross-platform
    single_platform: tuple[MatchedIon, ...] = ()  # one platform only
    unmatchable: tuple[MatchedIon, ...] = ()  # charge carrier unstated: can never pair

    @property
    def measurements(self) -> int:
        return sum(len(group.measurements) for group in self.all_groups)

    @property
    def all_groups(self) -> tuple[MatchedIon, ...]:
        return self.matched + self.single_platform + self.unmatchable

    @property
    def usable(self) -> tuple[MatchedIon, ...]:
        """Matched sets that may enter a fit. Synthetic and licence-blocked sets are not here."""
        return tuple(group for group in self.matched if group.usable)

    @property
    def synthetic_groups(self) -> tuple[MatchedIon, ...]:
        return tuple(group for group in self.all_groups if group.declares_synthetic)

    @property
    def conformer_groups(self) -> tuple[MatchedIon, ...]:
        return tuple(group for group in self.all_groups if group.conformer is not None)

    @property
    def platform_pair_counts(self) -> Mapping[str, int]:
        counted: Counter[str] = Counter()
        for group in self.matched:
            for pair in group.platform_pairs:
                counted[pair] += 1
        return dict(counted)

    @property
    def widest_match(self) -> int:
        """The largest number of platforms any one ion was measured on."""
        return max((len(group.platforms) for group in self.matched), default=0)

    @property
    def quotable(self) -> bool:
        """Whether any number from this report may be presented as describing real data.

        False if ANY group anywhere in the report is synthetic, not merely if the
        matched ones are. A report mixing real and invented records cannot be
        quoted without saying which is which, and the safe answer to "may I quote
        this" is no.
        """
        return not self.synthetic_groups

    def refusal(self) -> str | None:
        """Why this report may not be quoted as a result, or None if it may."""
        if self.synthetic_groups:
            synthetic = sum(group.synthetic_measurements for group in self.synthetic_groups)
            return (
                f"this report covers {synthetic} synthetic measurement(s) in {len(self.synthetic_groups)}"
                " group(s). They are fixtures declared in code, and no number computed over them describes"
                " real data. A report holding one may be read, and may not be quoted as a result"
            )
        return None

    def summary(self) -> str:
        lines = [
            "Matched-ion construction",
            f"  records in                {self.records_in}",
            f"  duplicates collapsed      {self.duplicates_collapsed}"
            "  (the same measurement republished, not a second measurement)",
            f"  measurements              {self.measurements}",
            f"  distinct ions             {len(self.all_groups)}",
            f"  MATCHED (>1 platform)     {len(self.matched)}",
            f"    usable                  {len(self.usable)}",
            f"    widest match            {self.widest_match} platforms",
            f"  one platform only         {len(self.single_platform)}",
            f"  unmatchable               {len(self.unmatchable)}"
            "  (charge carrier unstated; can never pair, not even with each other)",
        ]
        if self.platform_pair_counts:
            lines.append("  by platform pair:")
            for pair, count in sorted(self.platform_pair_counts.items()):
                lines.append(f"    {pair}: {count}")
        if self.conformer_groups:
            lines.append(
                f"  conformer-resolved ions   {len(self.conformer_groups)}"
                "  (held apart by conformer index; never merged, never deduplicated)"
            )
        if self.synthetic_groups:
            lines.append("")
            lines.append(f"  SYNTHETIC: {self.refusal()}")
        blocked = [group for group in self.matched if not group.usable]
        if blocked:
            lines.append("")
            lines.append(f"  {len(blocked)} matched set(s) may not enter a fit:")
            for group in blocked:
                lines.append(f"    {group}")
                for blocker in group.blockers:
                    lines.append(f"      - {blocker}")
        return "\n".join(lines)


def build_matched_ions(records: Sequence[object]) -> MatchingReport:
    """Group measurements into matched ions.

    The grouping key is the matched-ion key TOGETHER WITH the conformer index.
    Two conformers of one ion therefore never merge into one set, and never
    cancel each other out as duplicates: they are two values of one ion under one
    set of conditions, which M0 established is legitimate and not a collision.

    A set whose key is unmatchable - an ion whose charge carrier the source never
    named - is reported separately and can never be cross-platform, because such
    a key is unique to its own record by construction.
    """
    measurements, collapsed = deduplicate(records)

    grouped: dict[tuple, list[Measurement]] = {}
    order: list[tuple] = []
    for measurement in measurements:
        key = getattr(measurement.record, "matched_ion_key", None)
        conformer = getattr(measurement.record, "conformer", None)
        group_key = (key, conformer)
        if group_key not in grouped:
            grouped[group_key] = []
            order.append(group_key)
        grouped[group_key].append(measurement)

    matched: list[MatchedIon] = []
    single: list[MatchedIon] = []
    unmatchable: list[MatchedIon] = []
    for group_key in order:
        key, conformer = group_key
        ion = MatchedIon(key=key, conformer=conformer, measurements=tuple(grouped[group_key]))
        if key is not None and not getattr(key, "matchable", True):
            unmatchable.append(ion)
        elif ion.is_cross_platform:
            matched.append(ion)
        else:
            single.append(ion)

    return MatchingReport(
        records_in=len(records),
        duplicates_collapsed=collapsed,
        matched=tuple(matched),
        single_platform=tuple(single),
        unmatchable=tuple(unmatchable),
    )


class NotQuotableError(Exception):
    """A report holding synthetic measurements was asked to be quoted as a result.

    Deliberately not a ValueError, for the same reason TrainingGateError is not:
    a caller catching ValueError around a numeric routine must not swallow this
    along with a bad float.
    """


def assert_quotable(report: MatchingReport) -> None:
    """Raise unless every number in `report` describes real data.

    The one call that stands between a synthetic fixture and a figure in a
    document. Whatever reports a result in a later milestone calls this first.
    """
    refusal = report.refusal()
    if refusal is not None:
        raise NotQuotableError(refusal)
