"""Confidence grading: rules, not a fitted model.

WHAT A GRADE IS FOR
-------------------
A harmonized CCS is a corrected number, and a correction can be wrong in ways the
number itself cannot show. The grade is what travels beside it saying how much
weight it will bear. It is not a probability, it is not calibrated against
anything, and it does not come out of a model: it is a set of checks, each of
which can be evaluated today, each of which names a specific reason the
correction might not apply to this particular ion.

EVERY RULE HERE IS CHECKABLE NOW, AND NONE NEEDS A TRAINED MODEL. That is the
point of writing them as rules. When a model arrives it will produce a number and
an interval; these rules will still be what decides whether that number should be
trusted, and they will not need re-deriving.

A GRADE ONLY EVER FALLS
-----------------------
Every rule is a DEMOTION. There is no rule that raises a grade and no way to
offset one concern against another, because they are not commensurable: an ion
outside the calibration range is not made safer by its stratum being well
populated. The final grade is the WORST demotion found, not a score, not an
average, and not a count. Two mild concerns do not add up to a severe one, and a
single severe one is not diluted by everything else being fine.

WHAT IS DELIBERATELY NOT HERE
-----------------------------
No number is invented. Where a rule cannot be evaluated - no stratum, no fitted
correction - it returns nothing rather than assuming the best, and the caller is
left to decide. A rule that silently passed when it could not be checked would be
worse than no rule, because it would look like evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from statistics import median
from typing import Sequence

from .reuse import ReuseStatus, as_reuse_status, can_use
from .statistics import (
    MIN_POINTS_FOR_DEMING,
    MIN_POINTS_FOR_LIMITS_OF_AGREEMENT,
    StratumStatistics,
    deming_slope_intercept,
)

# --- thresholds -----------------------------------------------------------------------------

# POLICY. A correction fitted over a range of cross sections is interpolation
# inside that range and extrapolation outside it. A linear correction carried a
# little way past its data is arguable; carried a long way it is an assertion
# about ions nothing in the fit resembles. Ten per cent is where this draws the
# line, and it is a choice rather than a derivation.
EXTRAPOLATION_MARGIN = 0.10

# POLICY, anchored to something real. If dropping a stratum's outliers would move
# the correction it implies by more than this, at a typical ion, then the
# correction is substantially the outliers' doing rather than the bulk's. Half a
# per cent is a little over the published interlaboratory reproducibility of
# stepped-field DTIMS (0.29 per cent RSD), so a correction that moves by more than
# this on the strength of a few ions is moving by more than the reference method's
# own spread.
LEVERAGE_LIMIT_PERCENT = 0.5

# The population at which a stratum stops being worth a caveat of its own. Shared
# with readiness.py rather than redefined, so the two cannot drift apart.
from .readiness import TARGET_MATCHED_IONS  # noqa: E402  (imported here so the comment above it reads)


class ConfidenceGrade(StrEnum):
    """How much weight a harmonized value will bear. Ordered, worst wins.

    The words are deliberately not letters. A reader seeing "B" has to look up
    what B means and will guess; a reader seeing "qualified" already knows they
    are being told to read the caveat.
    """

    # Nothing counted against it. NOT a promise that it is right - it is a
    # statement that none of the checks below found a reason to doubt it.
    SUPPORTED = "supported"
    # Usable with its caveat attached. The caveat is not decoration.
    QUALIFIED = "qualified"
    # Usable only if the caveat is as prominent as the number.
    WEAK = "weak"
    # Do not use this number. Something is known to be wrong with applying this
    # correction to this ion, or the data behind it may not be used at all.
    UNSUPPORTED = "unsupported"


# Worst first. The order is explicit rather than relying on the enum's, so that
# adding a grade in the middle cannot silently reorder the comparison.
_SEVERITY = {
    ConfidenceGrade.SUPPORTED: 0,
    ConfidenceGrade.QUALIFIED: 1,
    ConfidenceGrade.WEAK: 2,
    ConfidenceGrade.UNSUPPORTED: 3,
}


@dataclass(frozen=True)
class Demotion:
    """One reason a grade fell, and how far."""

    rule: str
    grade: ConfidenceGrade
    detail: str

    def __str__(self) -> str:
        return f"[{self.grade}] {self.rule}: {self.detail}"


@dataclass(frozen=True)
class Confidence:
    """A grade, and every reason it is not higher."""

    grade: ConfidenceGrade
    demotions: tuple[Demotion, ...] = ()
    # Checks that could not be run at all, with why. Kept apart from demotions:
    # "this was not checked" is a different statement from "this was checked and
    # it was fine", and collapsing them would let an unevaluable rule read as a
    # passed one.
    not_checked: tuple[str, ...] = ()

    @property
    def usable(self) -> bool:
        return self.grade is not ConfidenceGrade.UNSUPPORTED

    def summary(self) -> str:
        lines = [f"confidence: {self.grade}"]
        for demotion in self.demotions:
            lines.append(f"  {demotion}")
        for skipped in self.not_checked:
            lines.append(f"  [not checked] {skipped}")
        if not self.demotions and not self.not_checked:
            lines.append("  no check found a reason to doubt this correction")
        return "\n".join(lines)


# --- the rules ------------------------------------------------------------------------------


def outside_calibration_range(ccs: float, stratum: StratumStatistics) -> Demotion | None:
    """The ion sits outside the range of cross sections the correction was fitted over.

    Inside the range the correction interpolates. A little outside it
    extrapolates, which is a caveat. A long way outside it asserts a relationship
    for ions that nothing in the fit resembles, which is not a caveat but a
    refusal.
    """
    if not stratum.points:
        return None
    values = [point.reference_ccs for point in stratum.points]
    low, high = min(values), max(values)
    if low <= ccs <= high:
        return None
    span = high - low
    margin = span * EXTRAPOLATION_MARGIN if span > 0 else abs(high) * EXTRAPOLATION_MARGIN
    if low - margin <= ccs <= high + margin:
        return Demotion(
            rule="outside the calibration range",
            grade=ConfidenceGrade.WEAK,
            detail=(
                f"{ccs:g} is outside the {low:g} to {high:g} the correction was fitted over, though within"
                f" {EXTRAPOLATION_MARGIN:.0%} of it. The correction is being extrapolated"
            ),
        )
    return Demotion(
        rule="far outside the calibration range",
        grade=ConfidenceGrade.UNSUPPORTED,
        detail=(
            f"{ccs:g} is more than {EXTRAPOLATION_MARGIN:.0%} beyond the {low:g} to {high:g} the correction"
            " was fitted over. Nothing in the fit resembles this ion, so the correction is an assertion"
            " rather than an estimate"
        ),
    )


def thinly_populated(stratum: StratumStatistics) -> Demotion | None:
    """The correction rests on too few matched ions to be worth much."""
    n = stratum.n
    if n < MIN_POINTS_FOR_DEMING:
        return Demotion(
            rule="no correction can be fitted",
            grade=ConfidenceGrade.UNSUPPORTED,
            detail=f"{n} matched ion(s): below the {MIN_POINTS_FOR_DEMING} a slope needs, there is no correction",
        )
    if n < MIN_POINTS_FOR_LIMITS_OF_AGREEMENT:
        return Demotion(
            rule="thinly populated calibration group",
            grade=ConfidenceGrade.WEAK,
            detail=(
                f"{n} matched ion(s): a correction from this few has limits of agreement that cannot be"
                f" quoted at all, which needs {MIN_POINTS_FOR_LIMITS_OF_AGREEMENT}"
            ),
        )
    if n < TARGET_MATCHED_IONS:
        return Demotion(
            rule="thinly populated calibration group",
            grade=ConfidenceGrade.QUALIFIED,
            detail=(
                f"{n} matched ion(s): below the {TARGET_MATCHED_IONS} at which a platform-pair figure is"
                " worth quoting without a caveat"
            ),
        )
    return None


def flagged_as_an_outlier(matched_ion_key: object, stratum: StratumStatistics) -> Demotion | None:
    """This ion is one the stratum already reported as not transferring.

    The most direct reason a grade can fall. Every other rule says the correction
    may not apply here; this one says the data already shows it does not. The
    published work finds a small fraction of ions disagreeing by as much as seven
    per cent while the bulk sit inside two, and those are exactly the ions a
    harmonized number would be confidently wrong about.
    """
    for point in stratum.outliers:
        if point.ion.key == matched_ion_key:
            return Demotion(
                rule="flagged as an outlier in its own stratum",
                grade=ConfidenceGrade.UNSUPPORTED,
                detail=(
                    f"this ion differs by {point.difference_percent:+.2f}% where its stratum's centre is"
                    f" {stratum.centre_percent:+.2f}%. The correction is fitted to the bulk and this ion is"
                    " not in it, so a harmonized value here would be confidently wrong"
                ),
            )
    return None


def the_outlier_rule_can_run_for(matched_ion_key: object, stratum: StratumStatistics) -> bool:
    """Whether `flagged_as_an_outlier` is capable of answering for THIS ion.

    IT IS A CORPUS LOOKUP. It asks whether this ion's key is among the ions the stratum
    already reported as not transferring, and a key that is not in the stratum at all
    cannot be among them. So for any ion NOT IN THE SEED CORPUS the rule does not return
    "not an outlier" - it returns nothing, having been unable to look.

    That distinction was invisible until 20 September 2026: the rule returned None either
    way, `grade_correction` treated None as "checked and fine", and every novel ion was
    graded as though the check had passed. The population this platform exists to serve is
    precisely the ions that are NOT already in its corpus, so the only ion-specific rule in
    the scheme was inert for all of them while the response asserted it had run.
    """
    return any(point.ion.key == matched_ion_key for point in stratum.points)


OUTLIER_RULE_UNEVALUABLE = (
    "outlier check: this ion is not in the stratum's own data, so whether it is one of the ions that"
    " already fail to transfer COULD NOT BE LOOKED UP. The rule is a corpus lookup keyed on identity and"
    " it can only ever answer for ions already measured on both platforms here. For a new ion it is not"
    " 'passed', it is unanswered - and this rule is the only one in the scheme that speaks about THIS ion"
    " rather than about the stratum around it. If you are bringing chemistry this platform has not"
    " measured, expect this note on EVERY response rather than occasionally: for that use the rule is"
    " inert always"
)
UNEVALUABLE_DEMOTION = (
    "{count} rule(s) of the scheme could not be evaluated for this ion, so this grade rests on fewer"
    " checks than the scheme advertises and is not the same grade as one that passed them all: {names}"
)


def unverified_source_behind_it(stratum: StratumStatistics) -> Demotion | None:
    """Some measurement the correction rests on may not be used.

    The licence gate already keeps such records out of a fit, so in a well-formed
    pipeline this should never fire. It is here because a grade that trusted the
    gate to have run would be trusting something it cannot see, and because the
    remedy is specific: read that source's terms, do not relabel the record.
    """
    offending: list[str] = []
    for point in stratum.points:
        for record in (point.reference, point.other):
            try:
                status = as_reuse_status(getattr(record, "reuse_status", None))
            except (TypeError, ValueError):
                offending.append(f"{getattr(record, 'source', '?')!r} (unreadable reuse status)")
                continue
            # can_use, not can_train_commercial: see the same correction in
            # matching.blockers. A commercial platform reading this repository
            # gets the stricter answer by setting reuse.PLATFORM_USE_CONTEXT, not
            # by changing this line.
            if not can_use(status):
                offending.append(f"{getattr(record, 'source', '?')!r} ({status.value})")
    if not offending:
        return None
    named = ", ".join(sorted(set(offending))[:3])
    return Demotion(
        rule="a source behind the correction may not be used",
        grade=ConfidenceGrade.UNSUPPORTED,
        detail=(
            f"{len(set(offending))} measurement(s) the correction rests on are not cleared for use: {named}."
            " The remedy is to read those terms and record them, never to relabel the record"
        ),
    )


def slope_leverage_percent(stratum: StratumStatistics) -> float | None:
    """How far the stratum's outliers move the correction, at a typical ion.

    A DIAGNOSTIC, AND ONLY A DIAGNOSTIC. This refits the slope with the outliers
    excluded and compares, purely to measure how much of the correction they are
    responsible for. NOTHING IS REMOVED from the reported statistics, from the
    correction, or from any count: the refit is thrown away and only the
    difference between the two is kept. Dropping outliers to get a nicer fit is
    the one thing this repository refuses hardest, and measuring their influence
    is the opposite of doing it.

    Returns None where the comparison cannot be made: no outliers to exclude, or
    too few points left to fit once they are.
    """
    outlier_keys = {id(point) for point in stratum.outliers}
    kept = [point for point in stratum.points if id(point) not in outlier_keys]
    if not outlier_keys or len(kept) < MIN_POINTS_FOR_DEMING:
        return None
    if stratum.agreement.deming_slope is None or stratum.agreement.deming_intercept is None:
        return None

    without = deming_slope_intercept(
        [point.reference_ccs for point in kept],
        [point.other_ccs for point in kept],
        stratum.agreement.lambda_used,
    )
    if without is None:
        return None
    slope_without, intercept_without = without

    # Compared at the middle of the range rather than at zero, because the
    # intercept alone says nothing about the correction an actual ion receives.
    middle = median([point.reference_ccs for point in stratum.points])
    if middle == 0:
        return None
    with_outliers = stratum.agreement.deming_slope * middle + stratum.agreement.deming_intercept
    without_outliers = slope_without * middle + intercept_without
    return 100.0 * abs(with_outliers - without_outliers) / middle


def correction_driven_by_outliers(stratum: StratumStatistics) -> Demotion | None:
    """The slope is substantially the work of ions that do not transfer.

    This is the M3 finding made into a rule. In the benchmark corpus three ions
    out of twenty-four carried a seven per cent bias and moved the fitted slope
    from the injected 1.020 to 1.031 - so the correction applied to every
    well-behaved ion was being set, in part, by the three that were not.
    """
    leverage = slope_leverage_percent(stratum)
    if leverage is None or leverage <= LEVERAGE_LIMIT_PERCENT:
        return None
    return Demotion(
        rule="the correction is driven by ions that do not transfer",
        grade=ConfidenceGrade.WEAK,
        detail=(
            f"excluding this stratum's {len(stratum.outliers)} outlier(s) would move the correction by"
            f" {leverage:.2f}% of CCS at a typical ion, past the {LEVERAGE_LIMIT_PERCENT:g}% this treats as"
            " the limit. The correction applied to well-behaved ions is substantially the work of ions that"
            " are not. The outliers are NOT removed; this is a measurement of their influence"
        ),
    )


def grade_correction(ccs: float, matched_ion_key: object, stratum: StratumStatistics) -> Confidence:
    """Every rule, applied to one ion against the stratum that would correct it.

    The grade is the worst demotion found. There is no scoring and no offsetting: a
    well-populated stratum does not make an extrapolation safe, and a licence that may not
    be used is not redeemed by anything at all.

    A RULE THAT CANNOT BE EVALUATED IS NEVER SILENTLY SKIPPED. It is named in
    `not_checked` with why it could not run, AND it demotes the grade by one notch. The
    second half is the part that was missing: an empty `not_checked` is used throughout
    this API as a positive claim that every rule ran - it is populated on 133 of 417
    estimates from the seed corpus - so a silently skipped rule was an active assertion
    that a check had passed when it had not been made.

    THE CORPUS FIGURES HERE UNDERSTATE THE EFFECT IN KIND, NOT IN DEGREE. They are measured
    by replaying the seed corpus against itself, and every record in that replay is by
    definition already in the corpus, so the one identity-keyed rule can be evaluated for
    most of them. For the use this platform exists for - an ion it has not measured - that
    rule is inert ALWAYS. 100 per cent, not 30.
    """
    # WHICH RULES COULD RUN AT ALL is decided before any of them is applied, because a rule
    # that returns None because it found nothing and a rule that returns None because it
    # could not look are different answers and were indistinguishable here until
    # 20 September 2026.
    not_checked: list[str] = []

    outlier_evaluable = the_outlier_rule_can_run_for(matched_ion_key, stratum)
    if not outlier_evaluable:
        not_checked.append(OUTLIER_RULE_UNEVALUABLE)
    if slope_leverage_percent(stratum) is None:
        not_checked.append(
            "leverage: this stratum has no outliers to exclude, or too few points left once they are, so"
            " how much the correction depends on them could not be measured"
        )
    if not stratum.points:
        not_checked.append("calibration range: the stratum holds no points, so there is no range to be inside")

    found = [
        rule
        for rule in (
            outside_calibration_range(ccs, stratum),
            thinly_populated(stratum),
            flagged_as_an_outlier(matched_ion_key, stratum) if outlier_evaluable else None,
            unverified_source_behind_it(stratum),
            correction_driven_by_outliers(stratum),
        )
        if rule is not None
    ]

    grade = ConfidenceGrade.SUPPORTED
    for demotion in found:
        if _SEVERITY[demotion.grade] > _SEVERITY[grade]:
            grade = demotion.grade

    # AN UNEVALUABLE RULE DEMOTES ONE NOTCH. Not to unsupported - the ion may be perfectly
    # fine - but a grade computed from fewer rules than the scheme advertises is not the
    # same grade as one computed from all of them, and nothing else in the response
    # distinguishes them. Clamped at unsupported, which has no notch below it.
    if not_checked:
        # One notch is a step in RANK, found by position in the ordered list. Indexing with
        # the severity number itself would give the same answer today only because those
        # numbers happen to run 0,1,2,3 with no gaps - a grade inserted at 5 would silently
        # skip rungs. Nothing would have caught that.
        order = sorted(_SEVERITY, key=lambda g: _SEVERITY[g])
        demoted = order[min(order.index(grade) + 1, len(order) - 1)]
        found.append(
            Demotion(
                rule="a rule of the scheme could not be evaluated for this ion",
                grade=demoted,
                detail=UNEVALUABLE_DEMOTION.format(
                    count=len(not_checked),
                    names="; ".join(note.split(":")[0] for note in not_checked),
                ),
            )
        )
        grade = demoted

    return Confidence(
        grade=grade,
        demotions=tuple(sorted(found, key=lambda d: -_SEVERITY[d.grade])),
        not_checked=tuple(not_checked),
    )


def grade_rules() -> tuple[dict, ...]:
    """The scheme itself, as data, so it can be published and argued with.

    Served by the API. A grading scheme nobody can read is a grading scheme
    nobody can challenge, and this one is rules rather than a fitted model
    precisely so that it CAN be challenged.
    """
    return (
        {
            "rule": "outside the calibration range",
            "falls_to": [ConfidenceGrade.WEAK.value, ConfidenceGrade.UNSUPPORTED.value],
            "why": "a correction interpolates inside the range it was fitted over and extrapolates outside it",
            "threshold": f"within {EXTRAPOLATION_MARGIN:.0%} beyond the range is weak, further is unsupported",
            "basis": "policy",
            "applies_to": "any ion: it compares the submitted value against the fitted range",
        },
        {
            "rule": "thinly populated calibration group",
            "falls_to": [
                ConfidenceGrade.QUALIFIED.value,
                ConfidenceGrade.WEAK.value,
                ConfidenceGrade.UNSUPPORTED.value,
            ],
            "why": "a correction from few matched ions carries limits nobody can quote",
            "threshold": (
                f"under {MIN_POINTS_FOR_DEMING} there is no slope at all; under"
                f" {MIN_POINTS_FOR_LIMITS_OF_AGREEMENT} no limits of agreement; under {TARGET_MATCHED_IONS}"
                " a caveat"
            ),
            "basis": "derived, except the last which is policy",
            "applies_to": "any ion: it is a property of the stratum, not of the ion",
        },
        {
            "rule": "flagged as an outlier in its own stratum",
            "falls_to": [ConfidenceGrade.UNSUPPORTED.value],
            "why": "the data already shows the correction does not fit this ion",
            "threshold": "more than the outlier margin from the stratum's median difference",
            "basis": "the published interplatform envelope",
            # A SCOPE LIMIT ON THE PRODUCT, published rather than left in a per-response
            # field. This is the only rule in the scheme that speaks about the SUBMITTED
            # ION rather than about the stratum around it, and it is the one that cannot
            # run for most callers.
            "applies_to": (
                "ONLY IONS ALREADY IN THIS CORPUS. This rule is a lookup: it asks whether the submitted"
                " ion is one the corpus already recorded as failing to transfer. A new ion is not in the"
                " corpus, so the rule cannot answer for it - and a new ion is the ordinary case for anyone"
                " using this service. FOR THAT CASE THE RULE IS INERT ALWAYS, NOT SOMETIMES: the rate is"
                " 100 per cent, not the 30 per cent that replaying this corpus against itself reports,"
                " because every record in that replay is by definition already present. The corpus figure"
                " is not a smaller version of the real one; it measures a different population. So a"
                " caller bringing new chemistry should read this rule as one the scheme lists and cannot"
                " apply to them at all. Where it cannot run it is named in the response's not_checked and"
                " the grade is demoted one notch; it is never treated as passed."
            ),
        },
        {
            "rule": "a source behind the correction may not be used",
            "falls_to": [ConfidenceGrade.UNSUPPORTED.value],
            "why": "an unverified or blocked licence is not a question of quality but one of permission",
            "threshold": "any measurement in the stratum whose reuse status is not cleared",
            "basis": "the licence gate, default deny",
            "applies_to": "any ion: it is a property of the records behind the correction, not of the submitted one",
        },
        {
            "rule": "the correction is driven by ions that do not transfer",
            "falls_to": [ConfidenceGrade.WEAK.value],
            "why": "a slope levered by a few outliers is applied to every well-behaved ion",
            "threshold": (
                f"excluding the outliers would move the correction by more than {LEVERAGE_LIMIT_PERCENT:g}%"
                " of CCS at a typical ion"
            ),
            "basis": "policy, anchored to stepped-field DTIMS interlaboratory reproducibility",
            "applies_to": "any ion, where leverage is measurable at all; where it is not, the response says so in not_checked",
        },
        {
            # A RULE ABOUT THE SCHEME ITSELF, and it is published for the same reason as the
            # other five: a demotion a caller can receive but cannot look up is a demotion
            # they cannot argue with. Its absence here was caught by the scheme-completeness
            # test the moment the demotion was added, which is what that test is for.
            "rule": "a rule of the scheme could not be evaluated for this ion",
            "falls_to": [
                ConfidenceGrade.QUALIFIED.value,
                ConfidenceGrade.WEAK.value,
                ConfidenceGrade.UNSUPPORTED.value,
            ],
            "why": (
                "a grade computed from fewer rules than the scheme advertises is not the same grade as one"
                " that passed them all, and an empty not_checked is used throughout this API as a positive"
                " claim that every rule was evaluated"
            ),
            "threshold": (
                "any rule that could not run for this ion. The grade falls ONE NOTCH from wherever the"
                " rules that did run left it - not to unsupported, because the ion may be perfectly fine."
                " Which rule, and why it could not run, is named in the response's not_checked"
            ),
            "basis": "policy",
            "applies_to": "any ion: it is about the scheme's own coverage, not about the data",
        },
    )
