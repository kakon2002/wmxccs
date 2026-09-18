"""Cross-platform statistics: association held apart from agreement, and both refused where they cannot be earned.

Every number this module produces is a claim about two instruments, so the tests
here are written against the CLAIM rather than against the arithmetic that
happens to make it. Each rule is asserted in both directions, because a module
that refused everything and a module that pooled everything each pass a suite
built from one half.

Four things are load-bearing and easy to lose quietly:

- ASSOCIATION AND AGREEMENT ANSWER DIFFERENT QUESTIONS. The benchmark corpus
  correlates at r = 0.998 and runs two per cent high on every ion. If Lin's
  concordance ever stops falling below Pearson r on that corpus, the one figure
  that distinguishes the two questions has stopped distinguishing them, and a
  reader is being told a biased platform agrees.
- THE ANSWER IS KNOWN IN ADVANCE. fixtures declares the biases it injects, so
  these tests compare a recovered figure against a DECLARED one rather than
  against whatever the code produced. The median recovers the injected bias; the
  mean and the Deming slope are dragged by the three ions that do not transfer,
  and that is asserted too, because it is the reason the centre is a median.
- OUTLIERS ARE FLAGGED AND KEPT. The three ions carrying the seven per cent bias
  are still in stratum.points, and every n counts them. A flag that removed them
  would improve every figure in the report and describe an instrument pair that
  does not exist.
- THE SYNTHETIC GUARD NEEDS A CONTROL. The corpora below are invented, and the
  guard that stops them being quoted is only meaningful beside records that are
  NOT invented. The house records here are internal data with no DOI, which
  clear the gate, so a module that refused everything would turn this file red.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from wmxccs import fixtures, statistics
from wmxccs.fixtures import OUTLIER_BIAS_PERCENT, TIMS_BIAS_PERCENT, TWIMS_BIAS_PERCENT
from wmxccs.identity import SmallMoleculeAnalyte
from wmxccs.licensing import declares_synthetic
from wmxccs.loader import load_measurements_file
from wmxccs.matching import build_matched_ions
from wmxccs.models import UncertaintyType
from wmxccs.reuse import ReuseStatus
from wmxccs.statistics import (
    MIN_POINTS_FOR_CORRELATION,
    MIN_POINTS_FOR_DEMING,
    MIN_POINTS_FOR_LIMITS_OF_AGREEMENT,
    OUTLIER_MARGIN_PERCENT,
    Agreement,
    Association,
    NotQuotableError,
    assert_quotable,
    choose_reference,
    compare_platforms,
    deming_slope_intercept,
    error_variance_ratio,
    lins_concordance,
    pearson_r,
)

REPO = Path(__file__).resolve().parents[1]
SEED_FILE_2016 = REPO / "data" / "seed" / "struwe2016_chemcommun.csv"
SEED_FILE_2015 = REPO / "data" / "seed" / "struwe2015_analyst.csv"

# Platform labels as statistics.py spells them. Stepped-field and single-field
# DTIMS are two platforms here; see matching._platform_of.
DTIMS = "DTIMS/stepped_field"
DTIMS_SINGLE = "DTIMS/single_field"
TWIMS = "TWIMS"
TIMS = "TIMS"
CYCLIC = "CYCLIC"

# Eight ions with real scatter on both axes. Near-collinear points make every
# regression agree with every other, which would hide exactly the differences
# between Deming and least squares that this file exists to assert.
SCATTERED_REFERENCE = (100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0)
SCATTERED_OTHER = (110.0, 115.0, 160.0, 150.0, 205.0, 190.0, 245.0, 230.0)

# Words that would name a point being thrown away. Nothing in the module may
# carry one, as a name or as a parameter.
REMOVAL_WORDS = ("drop", "remove", "exclude", "trim", "reject", "winsor", "discard", "prune", "clip")


# --- building corpora ----------------------------------------------------------------------


def comparison_of(records):
    """The comparison report for a corpus, built the way the pipeline builds it."""
    return compare_platforms(build_matched_ions(records))


def paired_corpus(values, *, tag, other_platform="twims", **overrides):
    """One synthetic ion per (reference, other) value pair: primary DTIMS against another platform.

    Every ion gets its own molecule. Two ions sharing one would share a
    matched-ion key, merge into a single set, and quietly shrink the comparison.
    """
    records = []
    for index, (reference_ccs, other_ccs) in enumerate(values):
        analyte = fixtures.small_molecule(fixtures.synthetic_inchikey(f"{tag}{index:02d}"))
        records.append(
            fixtures.measurement(
                analyte=analyte,
                platform="dtims",
                ccs=reference_ccs,
                source="fixture reference platform",
                doi=fixtures.DOI_A,
                **overrides,
            )
        )
        records.append(
            fixtures.measurement(
                analyte=analyte,
                platform=other_platform,
                ccs=other_ccs,
                source="fixture other platform",
                doi=fixtures.DOI_B,
                **overrides,
            )
        )
    return tuple(records)


def sloping_values(count: int, *, bias_percent: float = 2.0) -> list[tuple[float, float]]:
    """`count` ion pairs biased by `bias_percent`, with a little scatter so nothing is collinear."""
    values = []
    for index in range(count):
        reference = 150.0 + 10.0 * index
        values.append((reference, round(reference * (1 + bias_percent / 100.0) + (index % 3) * 0.1, 4)))
    return values


def sole_pair(report):
    assert len(report.pairs) == 1, f"expected one platform pair, got {len(report.pairs)}"
    return report.pairs[0]


def pair_of(report, reference_platform: str, other_platform: str):
    found = [
        pair
        for pair in report.pairs
        if (pair.reference_platform, pair.other_platform) == (reference_platform, other_platform)
    ]
    assert len(found) == 1, f"expected one {other_platform} against {reference_platform} pair, got {found}"
    return found[0]


def sole_stratum(pair):
    assert len(pair.strata) == 1, f"expected one calibration-group stratum, got {len(pair.strata)}"
    return pair.strata[0]


def only_stratum(report):
    return sole_stratum(sole_pair(report))


def refusal_mentioning(block, word: str) -> str:
    found = [refusal for refusal in block.refusals if word in refusal]
    assert len(found) == 1, f"expected one refusal mentioning {word!r}, got {block.refusals}"
    return found[0]


def with_uncertainty(spread, kind, **overrides):
    """One measurement whose spread is of a stated kind. The kind is the point."""
    return fixtures.measurement(ccs_uncertainty=spread, uncertainty_type=kind, **overrides)


def ols_slope(xs, ys) -> float:
    """Ordinary least squares of y on x, computed HERE because the module offers none.

    It exists in this file only so that the two properties OLS lacks - symmetry,
    and a slope that does not depend on which axis is called x - can be shown
    against the Deming fit rather than merely asserted about it.
    """
    mean_x, mean_y = sum(xs) / len(xs), sum(ys) / len(ys)
    sxx = sum((x - mean_x) ** 2 for x in xs)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    return sxy / sxx


# --- records that are NOT synthetic ---------------------------------------------------------
#
# Internal data with no DOI: trainable, needing no licence record, and declaring
# nothing synthetic in the record or in the analyte inside it. Without these the
# synthetic tests below would be satisfied by a module that refused everything.

HOUSE = ReuseStatus.INTERNAL_PROPRIETARY
HOUSE_SOURCE = "an in-house run"


def house_record(tag: str, platform: str, ccs: float):
    analyte = SmallMoleculeAnalyte(
        inchikey=fixtures.synthetic_inchikey(tag), source=HOUSE_SOURCE, reuse_status=HOUSE
    )
    return fixtures.measurement(
        analyte=analyte, platform=platform, ccs=ccs, source=HOUSE_SOURCE, doi=None, reuse_status=HOUSE
    )


def house_corpus(count: int = 4):
    """A genuine cross-platform comparison of our own, with nothing invented in it."""
    records = []
    for index in range(count):
        reference = 150.0 + 10.0 * index
        records.append(house_record(f"house{index:02d}", "dtims", reference))
        records.append(house_record(f"house{index:02d}", "twims", round(reference * 1.02 + 0.1 * (index % 2), 4)))
    return tuple(records)


# --- the corpora the assertions below rest on -----------------------------------------------


@pytest.fixture(scope="module")
def benchmark():
    return comparison_of(fixtures.benchmark_corpus())


@pytest.fixture(scope="module")
def twims_against_dtims(benchmark):
    return sole_stratum(pair_of(benchmark, DTIMS, TWIMS))


@pytest.fixture(scope="module")
def tims_against_dtims(benchmark):
    return sole_stratum(pair_of(benchmark, DTIMS, TIMS))


@pytest.fixture(scope="module")
def real_corpus():
    """Every record the two seed files hold, cleared or not. 117 of them."""
    held = []
    for path in (SEED_FILE_2016, SEED_FILE_2015):
        held.extend(load_measurements_file(path).records)
    return tuple(held)


def test_the_benchmark_corpus_is_the_corpus_these_tests_believe_it_is(benchmark):
    # Everything below reads figures off this corpus, so its shape is asserted
    # once here. A fixture that quietly changed size or platforms would make the
    # recovered biases below meaningless rather than wrong.
    assert benchmark.ions_considered == fixtures.BENCHMARK_IONS == 24
    assert len(benchmark.pairs) == 3
    assert {(pair.reference_platform, pair.other_platform) for pair in benchmark.pairs} == {
        (DTIMS, TWIMS),
        (DTIMS, TIMS),
        (TIMS, TWIMS),
    }
    assert benchmark.n_points == 72
    assert benchmark.ions_skipped_for_replicates == ()
    assert benchmark.unusable_sets_skipped == ()


# --- 1. association and agreement are different questions, and are never merged ---------------


def test_lins_concordance_is_lower_than_pearson_r_because_a_platform_can_correlate_and_still_run_high(
    twims_against_dtims,
):
    """The whole reason Lin's CCC is here.

    TWIMS is injected two per cent high on every ion in this corpus. Correlation
    cannot see that - it asks whether the points lie on A line - and concordance
    can, because it asks whether they lie on THE line, y = x.
    """
    association = twims_against_dtims.association
    agreement = twims_against_dtims.agreement
    assert association.pearson_r == pytest.approx(1.0, abs=0.01)  # 0.998: nearly perfect association
    assert agreement.mean_difference_percent == pytest.approx(TWIMS_BIAS_PERCENT, abs=1.0)
    assert agreement.median_difference_percent == pytest.approx(TWIMS_BIAS_PERCENT, abs=0.05)
    assert agreement.lins_ccc is not None
    assert agreement.lins_ccc < association.pearson_r
    # and lower by an amount a reader would notice, not by a rounding difference
    assert association.pearson_r - agreement.lins_ccc > 0.001


def test_a_perfectly_correlated_platform_running_two_per_cent_high_scores_one_on_r_and_less_on_concordance():
    # The same point with no scatter at all, so nothing but the bias can explain
    # the gap: r is exactly 1 and concordance is not.
    reference = [150.0, 200.0, 250.0, 300.0]
    biased = [value * 1.02 for value in reference]
    assert pearson_r(reference, biased) == pytest.approx(1.0, abs=1e-12)
    concordance = lins_concordance(reference, biased)
    assert concordance is not None
    assert concordance < 0.999


def test_association_and_agreement_stay_in_separate_blocks_with_no_combined_score_anywhere(benchmark):
    stratum = sole_stratum(pair_of(benchmark, DTIMS, TWIMS))
    assert isinstance(stratum.association, Association)
    assert isinstance(stratum.agreement, Agreement)
    # Averaging the two would produce a number answering neither question, so
    # nothing anywhere may be named one.
    subjects = [benchmark, benchmark.pairs[0], stratum, stratum.association, stratum.agreement, stratum.points[0]]
    for subject in subjects:
        offenders = [name for name in dir(subject) if "score" in name.lower()]
        assert offenders == [], f"{type(subject).__name__} offers a combined figure: {offenders}"
    assert [name for name in vars(statistics) if "score" in name.lower()] == []


def test_the_two_blocks_are_reported_under_headings_that_say_which_question_each_answers(twims_against_dtims):
    summary = "\n".join(twims_against_dtims.summary())
    assert "do they move together?" in summary
    assert "do they give the same number?" in summary
    assert "concordance, not correlation" in summary


# --- 2. the injected bias is recovered --------------------------------------------------------


def test_the_injected_twims_bias_is_recovered_by_the_median_difference_and_by_the_deming_slope(
    twims_against_dtims,
):
    agreement = twims_against_dtims.agreement
    assert agreement.median_difference_percent == pytest.approx(TWIMS_BIAS_PERCENT, abs=0.05)
    # The slope carries the same bias as a multiplier. The tolerance is looser
    # than the median's on purpose: three ions carry the seven per cent bias and
    # they are spread across the range, so they lever the slope upwards. That
    # they do is asserted below rather than tidied away.
    assert agreement.deming_slope == pytest.approx(1 + TWIMS_BIAS_PERCENT / 100.0, abs=0.015)


def test_the_smaller_injected_tims_bias_is_recovered_too_so_no_single_number_can_be_hard_coded(
    tims_against_dtims,
):
    agreement = tims_against_dtims.agreement
    assert agreement.median_difference_percent == pytest.approx(TIMS_BIAS_PERCENT, abs=0.05)
    assert agreement.mean_difference_percent == pytest.approx(TIMS_BIAS_PERCENT, abs=0.05)
    # No outlier was injected into the TIMS half, so the slope recovers its bias
    # much more tightly than the TWIMS slope does.
    assert agreement.deming_slope == pytest.approx(1 + TIMS_BIAS_PERCENT / 100.0, abs=0.002)
    assert TIMS_BIAS_PERCENT < TWIMS_BIAS_PERCENT  # the two are genuinely different numbers


def test_the_mean_difference_is_dragged_above_the_median_by_the_ions_that_do_not_transfer(
    twims_against_dtims,
):
    # Why the stratum centre is a median. Three of twenty-four ions carry a seven
    # per cent bias, and they move the mean by well over half a per cent.
    agreement = twims_against_dtims.agreement
    assert agreement.mean_difference_percent > agreement.median_difference_percent + 0.5
    assert agreement.median_difference_percent == pytest.approx(TWIMS_BIAS_PERCENT, abs=0.05)


# --- 3. no ordinary least squares -------------------------------------------------------------


def test_the_module_offers_no_ordinary_least_squares_under_any_name():
    forbidden = ("ols", "least_squares", "leastsquares", "linregress", "regress")
    offenders = sorted(name for name in vars(statistics) if any(word in name.lower() for word in forbidden))
    assert offenders == [], f"an OLS fit has appeared in the module: {offenders}"
    assert not hasattr(statistics, "ols")


def test_the_deming_fit_is_symmetric_at_lambda_one_where_ordinary_least_squares_is_not():
    """The property that makes Deming the right tool and OLS the wrong one.

    Both axes are measurements, so which platform is called x must not change the
    relationship. Fitting each way round gives reciprocal slopes here. The same
    two OLS fits do not: their slopes multiply to r squared, which is less than
    one whenever the points are not perfectly collinear.
    """
    forward, _intercept = deming_slope_intercept(SCATTERED_REFERENCE, SCATTERED_OTHER, 1.0)
    backward, _other_intercept = deming_slope_intercept(SCATTERED_OTHER, SCATTERED_REFERENCE, 1.0)
    assert forward * backward == pytest.approx(1.0, rel=1e-9)
    assert backward == pytest.approx(1.0 / forward, rel=1e-9)

    r = pearson_r(SCATTERED_REFERENCE, SCATTERED_OTHER)
    ols_product = ols_slope(SCATTERED_REFERENCE, SCATTERED_OTHER) * ols_slope(SCATTERED_OTHER, SCATTERED_REFERENCE)
    assert ols_product == pytest.approx(r**2, rel=1e-9)
    assert ols_product < 0.95, "these points are too collinear to tell the two fits apart"


def test_the_deming_slope_moves_with_the_error_variance_ratio_towards_the_least_squares_line():
    # Lambda is var(error in y) / var(error in x). Push all the error into y and
    # the fit has to approach OLS of y on x; at lambda 1 it does not, which is
    # the whole reason the ratio is carried and reported.
    at_one, _ = deming_slope_intercept(SCATTERED_REFERENCE, SCATTERED_OTHER, 1.0)
    at_four, _ = deming_slope_intercept(SCATTERED_REFERENCE, SCATTERED_OTHER, 4.0)
    nearly_all_error_in_y, _ = deming_slope_intercept(SCATTERED_REFERENCE, SCATTERED_OTHER, 1e6)
    least_squares = ols_slope(SCATTERED_REFERENCE, SCATTERED_OTHER)
    assert at_four < at_one
    assert least_squares < at_four
    assert nearly_all_error_in_y == pytest.approx(least_squares, rel=1e-4)
    assert at_one != pytest.approx(least_squares, rel=1e-3)


def test_a_deming_fit_refuses_an_error_variance_ratio_of_zero_or_less():
    # Lambda zero is OLS wearing a Deming name, and a negative ratio is not a
    # ratio of variances at all.
    for lam in (0.0, -1.0):
        with pytest.raises(ValueError) as caught:
            deming_slope_intercept(SCATTERED_REFERENCE, SCATTERED_OTHER, lam)
        assert "positive" in str(caught.value)


def test_a_deming_fit_returns_nothing_where_the_points_have_no_covariance_to_fit():
    # The negative direction of the fit itself: a vertical cloud has no direction.
    flat = [200.0, 200.0, 200.0, 200.0]
    assert deming_slope_intercept(flat, [180.0, 190.0, 200.0, 210.0], 1.0) is None
    assert pearson_r(flat, [180.0, 190.0, 200.0, 210.0]) is None


# --- 4. nothing is pooled ----------------------------------------------------------------------


def test_a_platform_pair_holding_two_calibration_strata_gets_no_pooled_figure():
    report = comparison_of(fixtures.benchmark_corpus_with_two_calibrants())
    pair = pair_of(report, DTIMS, TWIMS)
    assert len(pair.strata) == 2
    assert pair.pooled is None
    assert {stratum.n for stratum in pair.strata} == {12}
    assert pair.n == 24  # every ion is still counted, it is just not averaged


def test_the_pooling_refusal_names_how_many_strata_there_are_and_why_they_do_not_average():
    report = comparison_of(fixtures.benchmark_corpus_with_two_calibrants())
    pair = pair_of(report, DTIMS, TWIMS)
    refusal = pair.pooling_refusal
    assert refusal is not None
    assert "2 calibration-group strata" in refusal
    assert "not the same comparison" in refusal
    assert "NO POOLED FIGURE" in "\n".join(pair.summary())


def test_the_two_strata_differ_only_in_the_calibrant_which_is_the_difference_that_must_not_be_averaged():
    report = comparison_of(fixtures.benchmark_corpus_with_two_calibrants())
    pair = pair_of(report, DTIMS, TWIMS)
    other_groups = sorted(stratum.other_group for stratum in pair.strata)
    assert len(set(other_groups)) == 2
    assert any("dextran" in group for group in other_groups)
    assert any("polyalanine" in group for group in other_groups)
    assert len({stratum.reference_group for stratum in pair.strata}) == 1


def test_a_platform_pair_holding_one_stratum_does_get_a_pooled_figure(benchmark):
    # The positive direction. Without it a module that always refused to pool
    # would pass every test above.
    pair = pair_of(benchmark, DTIMS, TWIMS)
    assert len(pair.strata) == 1
    assert pair.pooled is not None
    assert pair.pooled is pair.strata[0]
    assert pair.pooling_refusal is None
    assert pair.pooled.agreement.deming_slope is not None
    assert "NO POOLED FIGURE" not in "\n".join(pair.summary())


# --- 5. outliers are reported, never removed, and are relative to their own stratum --------------


def expected_outlier_count() -> int:
    """How many benchmark ions carry the seven per cent bias, read off the fixture's own rule."""
    return sum(1 for index in range(fixtures.BENCHMARK_IONS) if index % fixtures.OUTLIER_EVERY == 0)


def test_exactly_the_ions_carrying_the_injected_seven_per_cent_bias_are_flagged(twims_against_dtims):
    assert expected_outlier_count() == 3
    outliers = twims_against_dtims.outliers
    assert len(outliers) == 3
    for point in outliers:
        assert point.difference_percent == pytest.approx(OUTLIER_BIAS_PERCENT, abs=0.05)


def test_the_stratum_centre_is_the_median_so_badly_transferring_ions_cannot_hide_behind_themselves(
    twims_against_dtims,
):
    # A mean centre would sit at 2.6 per cent here, dragged towards the very ions
    # the flag exists to find, and would report the pair's systematic offset as
    # larger than the bias that was injected.
    assert twims_against_dtims.centre_percent == pytest.approx(TWIMS_BIAS_PERCENT, abs=0.05)
    assert twims_against_dtims.centre_percent < twims_against_dtims.agreement.mean_difference_percent


def test_outliers_are_measured_from_their_own_stratums_centre_and_not_from_zero(twims_against_dtims):
    # Every ordinary ion in this stratum sits at about two per cent, which is the
    # margin itself. Measured from zero, half of them would trip the flag and the
    # three that genuinely did not transfer would be buried among them.
    ordinary = [
        point
        for point in twims_against_dtims.points
        if point.difference_percent == pytest.approx(TWIMS_BIAS_PERCENT, abs=0.1)
    ]
    assert len(ordinary) == 21
    from_zero = [point for point in twims_against_dtims.points if abs(point.difference_percent) > OUTLIER_MARGIN_PERCENT]
    assert len(from_zero) > len(twims_against_dtims.outliers)
    for point in ordinary:
        assert point not in twims_against_dtims.outliers
        assert not point.is_outlier_against(twims_against_dtims.centre_percent)


def test_a_flagged_ion_is_still_in_the_points_and_still_counted_by_every_statistic(twims_against_dtims):
    kept = [id(point) for point in twims_against_dtims.points]
    assert twims_against_dtims.outliers, "nothing was flagged, so this proves nothing"
    for point in twims_against_dtims.outliers:
        assert id(point) in kept
    assert twims_against_dtims.n == fixtures.BENCHMARK_IONS == 24
    assert twims_against_dtims.association.n == 24
    assert twims_against_dtims.agreement.n == 24
    assert twims_against_dtims.residuals.n == 24
    assert "reported and kept" in "\n".join(twims_against_dtims.summary())


def test_the_tims_pair_has_no_injected_outliers_and_flags_none(tims_against_dtims):
    # The negative direction. A flag that fired on every stratum would be as
    # useless as one that never fired.
    assert tims_against_dtims.n == 24
    assert tims_against_dtims.outliers == ()
    assert tims_against_dtims.centre_percent == pytest.approx(TIMS_BIAS_PERCENT, abs=0.05)
    assert "OUTLIERS" not in "\n".join(tims_against_dtims.summary())


def test_nothing_in_this_module_takes_a_parameter_that_would_remove_a_point():
    """There is no drop, no trim and no exclusion anywhere, and adding one is a bug.

    Asserted against the names rather than the behaviour because the behaviour is
    the absence of a feature: the guard has to fail when somebody ADDS the
    parameter, which is the only moment it could ever be wrong.
    """
    offenders: list[str] = []
    for name, member in vars(statistics).items():
        if getattr(member, "__module__", None) != statistics.__name__:
            continue
        candidates: list[tuple[str, object]] = [(name, member)]
        if inspect.isclass(member):
            candidates += [(f"{name}.{attribute}", value) for attribute, value in vars(member).items()]
        for label, value in candidates:
            if isinstance(value, property):
                value = value.fget
            if not callable(value):
                continue
            if any(word in label.lower() for word in REMOVAL_WORDS):
                offenders.append(label)
            try:
                signature = inspect.signature(value)
            except (TypeError, ValueError):  # a builtin or a slot wrapper
                continue
            offenders += [
                f"{label}({parameter})"
                for parameter in signature.parameters
                if any(word in parameter.lower() for word in REMOVAL_WORDS)
            ]
    assert offenders == [], f"something here can remove a point: {offenders}"


def test_the_report_counts_the_flagged_ions_without_taking_them_out_of_the_comparison(benchmark):
    assert len(benchmark.outliers) == 6  # three in TWIMS against DTIMS, three in TWIMS against TIMS
    assert benchmark.n_points == 72
    summary = benchmark.summary()
    assert "reported, never removed" in summary


# --- 6. every statistic carries its n, and refuses below its floor --------------------------------


def test_the_three_floors_are_pinned_and_say_which_are_derived_and_which_is_policy():
    # DERIVED: with two points Pearson r is exactly +1 or -1 whatever the points
    # are, so the third point is the first that can disagree with a line.
    assert MIN_POINTS_FOR_CORRELATION == 3
    # DERIVED: a Deming fit spends two degrees of freedom on a slope and an
    # intercept, so two points fit perfectly by construction.
    assert MIN_POINTS_FOR_DEMING == 3
    # POLICY, and the loosest figure in the module: the conventional floor below
    # which an SD of differences is too unstable to quote as a limit. It is not
    # derived from anything in this repository, and somebody with a real corpus
    # is meant to be able to argue with it.
    assert MIN_POINTS_FOR_LIMITS_OF_AGREEMENT == 10


def test_every_statistic_block_carries_the_n_it_was_computed_from(benchmark):
    for pair in benchmark.pairs:
        for stratum in pair.strata:
            assert stratum.n == len(stratum.points)
            assert stratum.association.n == stratum.n
            assert stratum.agreement.n == stratum.n
            assert stratum.residuals.n == stratum.n
            assert f"n = {stratum.n}" in "\n".join(stratum.association.summary())
            assert f"n = {stratum.n}" in "\n".join(stratum.agreement.summary())


def test_a_correlation_from_two_points_is_refused_because_r_is_plus_or_minus_one_by_construction():
    reference, other = [180.0, 260.0], [183.6, 265.2]
    # what would have been reported: two points always lie on a line
    assert pearson_r(reference, other) == pytest.approx(1.0, abs=1e-12)
    assert pearson_r(reference, [265.2, 183.6]) == pytest.approx(-1.0, abs=1e-12)

    stratum = only_stratum(comparison_of(paired_corpus(list(zip(reference, other)), tag="twopt")))
    assert stratum.n == 2
    assert stratum.association.n == 2
    assert stratum.association.pearson_r is None
    assert stratum.association.r_squared is None
    refusal = refusal_mentioning(stratum.association, "correlation")
    assert "+1 or -1" in refusal
    assert f"at least {MIN_POINTS_FOR_CORRELATION}" in refusal
    assert "2 paired ion(s)" in refusal


def test_a_deming_slope_from_two_points_is_refused_with_the_reason_it_cannot_be_earned():
    stratum = only_stratum(comparison_of(paired_corpus(sloping_values(2), tag="demtwo")))
    assert stratum.agreement.n == 2
    assert stratum.agreement.deming_slope is None
    assert stratum.agreement.deming_intercept is None
    refusal = refusal_mentioning(stratum.agreement, "Deming slope")
    assert f"at least {MIN_POINTS_FOR_DEMING}" in refusal
    assert "no residual left to disagree with it" in refusal


def test_a_correlation_and_a_slope_are_produced_at_the_floor_itself():
    # The positive direction of both floors. A module that refused at three as
    # well would satisfy every refusal test above.
    stratum = only_stratum(comparison_of(paired_corpus(sloping_values(3), tag="floor")))
    assert stratum.n == MIN_POINTS_FOR_CORRELATION == MIN_POINTS_FOR_DEMING == 3
    assert stratum.association.pearson_r is not None
    assert stratum.association.r_squared is not None
    assert stratum.agreement.deming_slope is not None
    assert stratum.agreement.lins_ccc is not None
    assert stratum.agreement.refusals  # the limits of agreement still refuse at three


def test_limits_of_agreement_are_refused_below_ten_points_and_produced_at_ten():
    nine = only_stratum(comparison_of(paired_corpus(sloping_values(9), tag="loanine")))
    assert nine.agreement.n == 9
    assert nine.agreement.limits_of_agreement is None
    refusal = refusal_mentioning(nine.agreement, "limits of agreement")
    assert "standard deviation of the differences" in refusal
    assert str(MIN_POINTS_FOR_LIMITS_OF_AGREEMENT) in refusal

    ten = only_stratum(comparison_of(paired_corpus(sloping_values(10), tag="loaten")))
    assert ten.agreement.n == MIN_POINTS_FOR_LIMITS_OF_AGREEMENT == 10
    assert ten.agreement.limits_of_agreement is not None
    low, high = ten.agreement.limits_of_agreement
    assert low < ten.agreement.bias < high
    assert high - low == pytest.approx(2 * statistics.LOA_MULTIPLIER * ten.agreement.bias_sd, rel=1e-12)
    assert not [refusal for refusal in ten.agreement.refusals if "limits of agreement" in refusal]


def test_a_small_comparison_warns_about_its_own_weakness_rather_than_refusing_everything(benchmark):
    # Constraint 8: refuse on inevaluability, warn on weakness.
    weak = only_stratum(comparison_of(paired_corpus(sloping_values(4), tag="weak")))
    assert weak.agreement.bias is not None
    assert weak.agreement.deming_slope is not None
    assert any("small comparison" in warning for warning in weak.agreement.warnings)
    assert "weak:" in "\n".join(weak.agreement.summary())
    # and a comparison that is not small does not warn
    assert sole_stratum(pair_of(benchmark, DTIMS, TWIMS)).agreement.warnings == ()


def test_a_single_pair_still_reports_its_difference_and_refuses_everything_that_needs_a_spread():
    # One disagreement is a real observation. A limit, a slope and a correlation
    # are not, and each refuses separately.
    stratum = only_stratum(comparison_of(paired_corpus([(200.0, 204.0)], tag="lonely")))
    assert stratum.n == 1
    assert stratum.agreement.bias == pytest.approx(4.0, rel=1e-12)
    assert stratum.agreement.mean_difference_percent == pytest.approx(2.0, rel=1e-12)
    assert stratum.agreement.bias_sd is None
    assert stratum.agreement.limits_of_agreement is None
    assert stratum.agreement.deming_slope is None
    assert stratum.association.pearson_r is None
    assert stratum.residuals.correlation_with_reference_ccs is None


def test_r_squared_is_the_square_of_the_correlation_and_not_a_regression_fit(twims_against_dtims):
    # A reader who sees "R squared" will think it describes how well a line fits.
    # It does not: it is the shared variation between two platforms, and it says
    # nothing about whether they give the same number.
    association = Association(n=5, pearson_r=-0.6)
    assert association.r_squared == pytest.approx(0.36, rel=1e-12)
    assert Association(n=4).r_squared is None

    real = twims_against_dtims.association
    assert real.r_squared == pytest.approx(real.pearson_r**2, rel=1e-12)
    assert real.r_squared < abs(real.pearson_r)
    assert "NOT a regression fit" in "\n".join(real.summary())


def test_a_platform_with_no_variation_at_all_is_refused_by_name():
    flat = only_stratum(
        comparison_of(paired_corpus([(200.0, 204.0), (200.0, 208.0), (200.0, 212.0)], tag="flat"))
    )
    assert flat.association.pearson_r is None
    refusal = refusal_mentioning(flat.association, "no variation")
    assert "reference" in refusal


# --- 7. uncertainty conversion, which is where uncertainty_type earns its place ------------------


def test_a_standard_deviation_converts_to_itself_so_the_ratio_is_the_ratio_of_the_variances():
    ratio, basis = error_variance_ratio(
        [with_uncertainty(1.0, UncertaintyType.SD)], [with_uncertainty(2.0, UncertaintyType.SD)]
    )
    assert ratio == pytest.approx(4.0, rel=1e-12)  # (2 / 1) squared
    assert "measured from the reported uncertainties" in basis
    assert "ASSUMED" not in basis


def test_a_two_standard_deviation_spread_is_halved_so_it_weighs_the_same_as_one_standard_deviation():
    reference = [with_uncertainty(1.0, UncertaintyType.SD)]
    as_two_sd, two_sd_basis = error_variance_ratio(reference, [with_uncertainty(2.0, UncertaintyType.TWO_SD)])
    as_one_sd, one_sd_basis = error_variance_ratio(reference, [with_uncertainty(1.0, UncertaintyType.SD)])
    assert as_two_sd == pytest.approx(as_one_sd, rel=1e-12)
    assert as_two_sd == pytest.approx(1.0, rel=1e-12)
    assert two_sd_basis == one_sd_basis
    assert "measured from the reported uncertainties" in two_sd_basis


def test_a_standard_error_converts_only_where_the_replicate_count_is_recorded():
    reference = [with_uncertainty(1.0, UncertaintyType.SD)]
    # SD = SEM * sqrt(replicates), so a standard error of 1.0 over four
    # replicates is a standard deviation of 2.0, and a variance ratio of 4.
    with_count, measured = error_variance_ratio(
        reference, [with_uncertainty(1.0, UncertaintyType.SEM, replicates=4)]
    )
    assert with_count == pytest.approx(4.0, rel=1e-12)
    assert "measured from the reported uncertainties" in measured

    without_count, assumed = error_variance_ratio(reference, [with_uncertainty(1.0, UncertaintyType.SEM)])
    assert without_count == pytest.approx(1.0, rel=1e-12)
    assert "ASSUMED EQUAL" in assumed


def test_a_ninety_five_per_cent_interval_is_refused_outright_because_its_width_is_never_fixed():
    # A CI95 may be a half-width or a full width and this repository has never
    # said which. Dividing by 1.96 on a guess would silently halve or double
    # every weight built on it.
    ratio, basis = error_variance_ratio(
        [with_uncertainty(1.0, UncertaintyType.SD)], [with_uncertainty(1.96, UncertaintyType.CI95)]
    )
    assert ratio == pytest.approx(1.0, rel=1e-12)
    assert "ASSUMED EQUAL" in basis


def test_an_unknown_uncertainty_type_converts_to_nothing_which_is_the_point_of_the_value():
    ratio, basis = error_variance_ratio(
        [with_uncertainty(1.0, UncertaintyType.SD)], [with_uncertainty(2.0, UncertaintyType.UNKNOWN)]
    )
    assert ratio == pytest.approx(1.0, rel=1e-12)
    assert "ASSUMED EQUAL" in basis


def test_an_assumed_lambda_says_it_was_assumed_rather_than_measured():
    """The honesty of the basis string is the guard.

    Falling back to one is not wrong. Reporting the fallback as though the ratio
    had been measured from the data is, because a reader has no other way to tell
    an orthogonal fit from a weighted one.
    """
    bare = fixtures.measurement(ccs_uncertainty=None, uncertainty_type=None)
    ratio, basis = error_variance_ratio([bare], [bare])
    assert ratio == pytest.approx(1.0, rel=1e-12)
    assert "ASSUMED EQUAL" in basis
    assert "not every paired record reports a convertible uncertainty" in basis
    assert basis != "measured from the reported uncertainties of both platforms"


def test_the_assumed_lambda_travels_into_the_agreement_block_and_into_its_summary():
    stratum = only_stratum(
        comparison_of(
            paired_corpus(sloping_values(4), tag="nolam", ccs_uncertainty=None, uncertainty_type=None)
        )
    )
    assert stratum.agreement.lambda_used == pytest.approx(1.0, rel=1e-12)
    assert "ASSUMED EQUAL" in stratum.agreement.lambda_basis
    assert "lambda ASSUMED EQUAL" in "\n".join(stratum.agreement.summary())


def test_a_measured_lambda_travels_into_the_agreement_block_when_both_sides_report_a_spread(benchmark):
    # The positive direction: every benchmark record carries an SD, so the ratio
    # is measured rather than assumed.
    agreement = sole_stratum(pair_of(benchmark, DTIMS, TWIMS)).agreement
    assert agreement.lambda_used == pytest.approx(1.0, rel=1e-12)
    assert "measured from the reported uncertainties" in agreement.lambda_basis
    assert "ASSUMED" not in agreement.lambda_basis


# --- 8. the synthetic guard, and a control that is not synthetic ------------------------------------


def test_a_report_over_synthetic_fixtures_is_not_quotable_and_says_why(benchmark):
    assert benchmark.synthetic is True
    assert benchmark.quotable is False
    refusal = benchmark.refusal()
    assert refusal is not None
    assert "synthetic" in refusal
    assert "may not be quoted" in refusal
    assert "SYNTHETIC" in benchmark.summary()


def test_assert_quotable_refuses_a_synthetic_statistics_report(benchmark):
    with pytest.raises(NotQuotableError) as caught:
        assert_quotable(benchmark)
    assert "synthetic" in str(caught.value)


def test_not_quotable_error_is_not_a_value_error():
    # A caller catching ValueError around a numeric routine must not swallow this
    # along with a bad float.
    assert issubclass(NotQuotableError, Exception)
    assert not issubclass(NotQuotableError, ValueError)


def test_the_control_records_this_file_relies_on_really_are_not_synthetic():
    # If these ever started declaring themselves synthetic, every "is quotable"
    # test below would pass for the wrong reason and prove nothing.
    for record in house_corpus():
        assert not declares_synthetic(record)


def test_a_report_over_real_records_is_quotable_and_assert_quotable_returns_cleanly():
    report = comparison_of(house_corpus())
    assert report.synthetic is False
    assert report.quotable is True
    assert report.refusal() is None
    assert assert_quotable(report) is None
    # and it really is a comparison, not an empty report that refused everything
    stratum = only_stratum(report)
    assert stratum.n == 4
    assert stratum.agreement.deming_slope is not None
    assert "SYNTHETIC" not in report.summary()


# --- 9. the real corpus refuses, and says why --------------------------------------------------------


def test_the_real_corpus_supports_no_cross_platform_comparison_at_all(real_corpus):
    assert len(real_corpus) == 117
    report = comparison_of(real_corpus)
    assert report.pairs == ()
    assert report.n_points == 0
    assert report.outliers == ()
    assert report.ions_considered == 0


def test_the_real_corpus_names_the_single_platform_reason_and_says_more_of_it_will_not_help(real_corpus):
    """The most important thing this report tells a reader.

    A corpus that is too small and a corpus that can never work however much of
    it arrives call for completely different decisions, and "nothing here" does
    not tell them apart. These 117 records are all TWIMS: more TWIMS values will
    not produce a single cross-platform pair.
    """
    report = comparison_of(real_corpus)
    reason = report.no_comparison_reason
    assert reason is not None
    assert "single platform" in reason
    assert "TWIMS" in reason
    assert "More values on the same platform will not change that" in reason
    summary = report.summary()
    assert "NO CROSS-PLATFORM COMPARISON IS POSSIBLE FROM THIS CORPUS." in summary
    assert reason in summary


def test_the_real_corpus_is_quotable_even_though_it_has_nothing_to_say(real_corpus):
    # The synthetic guard must not be refusing everything: these records are real.
    report = comparison_of(real_corpus)
    assert report.synthetic is False
    assert report.quotable is True
    assert assert_quotable(report) is None


def test_a_corpus_whose_only_matched_set_is_unusable_says_so_rather_than_saying_nothing():
    report = comparison_of(fixtures.a_match_with_an_unusable_member())
    assert report.pairs == ()
    reason = report.no_comparison_reason
    assert reason is not None
    assert "1 matched set(s) exist and none may be used" in reason
    assert "least usable member" in reason


def test_a_corpus_with_nothing_matched_at_all_says_that_instead():
    report = comparison_of(())
    assert report.pairs == ()
    assert report.no_comparison_reason == statistics.NOTHING_MATCHED


# --- 10. the decisions behind each number --------------------------------------------------------


def test_the_difference_percent_is_taken_relative_to_the_reference_and_not_to_the_pair_mean():
    # 100 to 110 is ten per cent of the reference and 9.52 per cent of the pair
    # mean. The published figures this is compared against - within 1 per cent
    # for TIMS and 2 per cent for TWIMS RELATIVE TO DTIMS - are the former.
    stratum = only_stratum(comparison_of(paired_corpus([(100.0, 110.0)], tag="relref")))
    (point,) = stratum.points
    assert point.reference_ccs == pytest.approx(100.0, rel=1e-12)
    assert point.other_ccs == pytest.approx(110.0, rel=1e-12)
    assert point.difference == pytest.approx(10.0, rel=1e-12)
    assert point.difference_percent == pytest.approx(10.0, rel=1e-12)
    assert point.difference_percent != pytest.approx(100.0 * 10.0 / 105.0, rel=1e-3)
    assert stratum.agreement.mean_difference_percent == pytest.approx(10.0, rel=1e-12)
    assert stratum.agreement.mape == pytest.approx(10.0, rel=1e-12)


def test_the_difference_keeps_its_sign_so_a_platform_running_low_is_not_reported_as_running_high():
    stratum = only_stratum(comparison_of(paired_corpus([(200.0, 190.0)], tag="signed")))
    (point,) = stratum.points
    assert point.difference == pytest.approx(-10.0, rel=1e-12)
    assert point.difference_percent == pytest.approx(-5.0, rel=1e-12)
    assert stratum.agreement.bias == pytest.approx(-10.0, rel=1e-12)
    assert stratum.agreement.mae == pytest.approx(10.0, rel=1e-12)  # the unsigned figure is separate


def test_a_primary_stepped_field_value_is_the_reference_over_a_calibrated_platform():
    """A primary value is the only CCS in a comparison that rests on nobody else's.

    Both cases here are ones where the alphabetical fallback would have chosen
    the other platform, so a fit that ignored the primary rule would anchor the
    comparison on a calibrated value and measure nothing.
    """
    against_cyclic = sole_pair(comparison_of(fixtures.a_cyclic_pair_at_different_pass_counts()))
    assert against_cyclic.reference_platform == DTIMS
    assert against_cyclic.other_platform == CYCLIC
    assert CYCLIC < DTIMS  # alphabetically the cyclic value would have won

    against_single_field = sole_pair(comparison_of(fixtures.a_primary_and_a_calibrated_dtims_pair()))
    assert against_single_field.reference_platform == DTIMS
    assert against_single_field.other_platform == DTIMS_SINGLE
    assert DTIMS_SINGLE < DTIMS


def test_with_no_primary_side_the_reference_is_the_alphabetically_first_platform():
    # The negative direction, and it is honest about being arbitrary: comparing
    # two calibrated platforms anchors nothing, and the Deming slope is symmetric
    # so the choice changes no conclusion.
    records = {
        TWIMS: fixtures.measurement(platform="twims"),
        TIMS: fixtures.measurement(platform="tims"),
    }
    assert choose_reference(TWIMS, TIMS, records) == (TIMS, TWIMS)
    assert choose_reference(TIMS, TWIMS, records) == (TIMS, TWIMS)


def test_an_ion_with_replicates_on_one_side_is_held_rather_than_averaged_or_picked_from():
    # Averaging would merge two originals, which this platform does not do
    # anywhere, and picking the first would make the result depend on file order.
    twin = fixtures.small_molecule(fixtures.synthetic_inchikey("reptwin"))
    replicated = (
        fixtures.measurement(analyte=twin, platform="dtims", ccs=200.0, doi=fixtures.DOI_A),
        fixtures.measurement(analyte=twin, platform="twims", ccs=204.0, doi=fixtures.DOI_B),
        fixtures.measurement(
            analyte=twin, platform="twims", ccs=206.0, source="a second run", doi=fixtures.DOI_B
        ),
    )
    report = comparison_of(paired_corpus(sloping_values(3), tag="repclean") + replicated)

    assert report.ions_considered == 4
    stratum = only_stratum(report)
    assert stratum.n == 3
    assert len(report.ions_skipped_for_replicates) == 1
    assert "REPTWIN" in report.ions_skipped_for_replicates[0]
    other_values = {point.other_ccs for point in stratum.points}
    assert 204.0 not in other_values and 206.0 not in other_values  # neither replicate was picked
    assert 205.0 not in other_values  # and they were not averaged either
    assert "ions held for replicates" in report.summary()


def test_a_matched_set_with_a_member_nobody_may_use_is_skipped_and_counted_rather_than_compared():
    report = comparison_of(
        paired_corpus(sloping_values(3), tag="usable") + fixtures.a_match_with_an_unusable_member()
    )
    assert report.ions_considered == 4
    assert len(report.unusable_sets_skipped) == 1
    assert "UNUSABLE" in report.unusable_sets_skipped[0]
    stratum = only_stratum(report)
    assert stratum.n == 3
    assert report.n_points == 3
    assert "matched sets not used" in report.summary()


def test_a_synthetic_set_is_still_compared_because_the_arithmetic_has_to_be_exercisable(benchmark):
    # The one blocker carried through rather than skipped. Everything in the
    # benchmark corpus is synthetic, and if that blocked the comparison the
    # arithmetic could never be exercised at all. `quotable` is what stops the
    # result being presented as real.
    assert benchmark.unusable_sets_skipped == ()
    assert benchmark.n_points == 72
    assert benchmark.quotable is False


def test_coverage_counts_the_ions_inside_the_band_and_not_the_ions_outside_it(
    twims_against_dtims, tims_against_dtims
):
    twims = twims_against_dtims.agreement.coverage
    # Every TWIMS ion is about two per cent high, so none of them is within one.
    assert twims[1.0] == pytest.approx(0.0, abs=1e-9)
    inside_two = sum(1 for point in twims_against_dtims.points if abs(point.difference_percent) <= 2.0)
    assert 0 < inside_two < twims_against_dtims.n, "the band does not discriminate here"
    assert twims[2.0] == pytest.approx(100.0 * inside_two / twims_against_dtims.n, rel=1e-12)

    # The TIMS half is one per cent out, so all of it is inside the two per cent band.
    assert tims_against_dtims.agreement.coverage[2.0] == pytest.approx(100.0, abs=1e-9)
    assert "within 1%" in "\n".join(tims_against_dtims.agreement.summary())


def test_coverage_is_measured_against_zero_while_outlier_flagging_is_measured_against_the_centre(
    twims_against_dtims,
):
    # Two different questions: how close is this platform to the reference, and
    # which ions do not behave like their neighbours. The same stratum answers
    # them differently, which is the evidence they are not the same question.
    assert twims_against_dtims.agreement.coverage[2.0] < 100.0  # measured from zero
    assert len(twims_against_dtims.outliers) == 3  # measured from the stratum centre
    assert twims_against_dtims.centre_percent == pytest.approx(TWIMS_BIAS_PERCENT, abs=0.05)


def test_the_residual_trend_is_reported_against_reference_ccs_and_says_mass_is_not_held(
    twims_against_dtims,
):
    # CCS is the size proxy because no record in this repository carries a mass.
    assert twims_against_dtims.residuals.n == 24
    assert twims_against_dtims.residuals.correlation_with_reference_ccs is not None
    assert "mass is not held by any record" in "\n".join(twims_against_dtims.residuals.summary())
