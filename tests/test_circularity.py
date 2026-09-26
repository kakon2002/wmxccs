"""A calibration must not be compared against its own reference set.

THE RISK THE SCHEMA WAS BUILT FOR AND NOTHING ACTED ON. `CalibrationReference` has recorded,
since M0, which reference set a calibrated value came from, which publication those reference
values were in, and which platform they were themselves measured on. Its docstring names the
consequence exactly: "if the platform then compares that TWIMS value against the same DTIMS
values it calibrated from, the agreement it reports is partly an artefact of the calibration
rather than a fact about the instruments."

Until 26 September 2026 nothing read those fields at comparison time. The lineage was
recorded and the circularity was undetected - a guard that looks present and is not, which is
the class LIMITATIONS 4.5 collects. `statistics.circularity_between` is the detection and
these are its tests.

WHY IT IS A REFUSAL AND NOT A FLAG. The case in hand is the Bush Lab MicroSource records:
1,437 TWIMS values whose reference values are their own laboratory's nitrogen DTIMS
measurements. A flag on a figure is read by whoever goes looking for it. If Bush Lab DTIMS
data arrives later, the pair must not form a stratum at all - so the refusal happens before a
stratum exists, and is counted on the report rather than dropped, because a refusal nobody can
count cannot be told apart from an absence of data.
"""

from __future__ import annotations

import pytest

from wmxccs import fixtures
from wmxccs.matching import build_matched_ions
from wmxccs.models import (
    CalibrationReference,
    CCSMeasurement,
    DriftGas,
    DTIMSMethod,
    IMSType,
    Polarity,
)
from wmxccs.statistics import circularity_between, compare_platforms

REFERENCE_PAPER = "10.1021/acs.analchem.7b01709"
OTHER_PAPER = "10.1021/jasms.2c00196"

LINEAGE = CalibrationReference(
    reference_set="the authors' own nitrogen DTIMS values, measured for this paper",
    doi=REFERENCE_PAPER,
    platform=IMSType.DTIMS,
    method=DTIMSMethod.STEPPED_FIELD,
)


def calibrated(analyte, **overrides) -> CCSMeasurement:
    """A TWIMS value calibrated against DTIMS reference values from REFERENCE_PAPER."""
    fields = dict(
        analyte=analyte,
        adduct="[M+H]+",
        charge=1,
        polarity=Polarity.POSITIVE,
        ims_type=IMSType.TWIMS,
        drift_gas=DriftGas.N2,
        calibrant="PolyAla and nine drug standards",
        calibration_reference=LINEAGE,
        ccs=157.1,
        source="Bush Lab CCS database",
        source_locator=f"MicroSource sheet; published in doi:{REFERENCE_PAPER}",
        reuse_status=fixtures.FIXTURE,
    )
    fields.update(overrides)
    return CCSMeasurement(**fields)


def the_reference_platform(analyte, **overrides) -> CCSMeasurement:
    """THE SYNTHETIC DTIMS RECORD CARRYING THE SAME LINEAGE, which is what must be refused.

    Stepped-field DTIMS, so it is a PRIMARY value and carries no calibrant of its own - it is
    the thing the TWIMS value above was calibrated against.
    """
    fields = dict(
        analyte=analyte,
        adduct="[M+H]+",
        charge=1,
        polarity=Polarity.POSITIVE,
        ims_type=IMSType.DTIMS,
        dtims_method=DTIMSMethod.STEPPED_FIELD,
        drift_gas=DriftGas.N2,
        ccs=155.0,
        source="Bush Lab CCS database",
        source_locator=f"DTIM reference values; published in doi:{REFERENCE_PAPER}",
        reuse_status=fixtures.FIXTURE,
    )
    fields.update(overrides)
    return CCSMeasurement(**fields)


# --- the rule, on a pair -------------------------------------------------------------------


def test_a_calibrated_value_against_its_own_reference_paper_is_refused():
    analyte = fixtures.small_molecule(fixtures.synthetic_inchikey("circ0001"))
    why = circularity_between(calibrated(analyte), the_reference_platform(analyte))
    assert why is not None
    assert "refused as circular" in why
    assert REFERENCE_PAPER in why
    assert "measures how well the calibration reproduces its own reference set" in why


def test_the_refusal_is_symmetric_in_the_order_of_the_pair():
    """Which record is handed in first is an accident of iteration, not a fact about them."""
    analyte = fixtures.small_molecule(fixtures.synthetic_inchikey("circ0002"))
    one, other = calibrated(analyte), the_reference_platform(analyte)
    assert circularity_between(one, other) is not None
    assert circularity_between(other, one) is not None


def test_a_different_paper_on_the_same_platform_is_NOT_refused():
    """The other direction, and the one that keeps this from refusing everything.

    A DTIMS value from a different study is an independent comparison, which is the entire
    purpose of this platform. Refusing it would be worse than the bug.
    """
    analyte = fixtures.small_molecule(fixtures.synthetic_inchikey("circ0003"))
    independent = the_reference_platform(
        analyte,
        doi=OTHER_PAPER,
        source="another laboratory",
        source_locator="Table 2 of a different study",
    )
    assert circularity_between(calibrated(analyte), independent) is None


def test_a_value_with_no_recorded_lineage_is_not_refused():
    """Null lineage is "nobody recorded it", not "it is circular"."""
    analyte = fixtures.small_molecule(fixtures.synthetic_inchikey("circ0004"))
    no_lineage = calibrated(analyte, calibration_reference=None)
    assert circularity_between(no_lineage, the_reference_platform(analyte)) is None


def test_a_reference_set_measured_on_another_platform_is_not_refused():
    """The lineage must name the platform in the pair, or it is a different chain."""
    analyte = fixtures.small_molecule(fixtures.synthetic_inchikey("circ0005"))
    twims_referenced = calibrated(
        analyte,
        calibration_reference=CalibrationReference(
            reference_set="somebody else's TWIMS ladder",
            doi=REFERENCE_PAPER,
            platform=IMSType.TWIMS,
        ),
    )
    assert circularity_between(twims_referenced, the_reference_platform(analyte)) is None


def test_the_dtims_method_must_match_the_lineage_when_the_lineage_states_one():
    """A stepped-field reference set is not the single-field values from the same paper.

    Those are different numbers measured differently, and only one of them is what the
    calibration used.
    """
    analyte = fixtures.small_molecule(fixtures.synthetic_inchikey("circ0006"))
    single_field = the_reference_platform(
        analyte, dtims_method=DTIMSMethod.SINGLE_FIELD, calibrant="a ladder"
    )
    assert circularity_between(calibrated(analyte), single_field) is None


def test_the_publication_is_matched_in_the_locator_as_well_as_in_the_doi():
    """The records this was written for keep the paper in the locator, not in `doi`.

    Permission for the Bush Lab values comes from the database page, not from the ACS paper,
    so the licence gate requires `doi` to be null and the publication is named in
    `source_locator`. A test that only exercised the `doi` route would pass while the guard
    never fired on the one source it exists for.
    """
    analyte = fixtures.small_molecule(fixtures.synthetic_inchikey("circ0007"))
    reference_by_locator = the_reference_platform(analyte, doi=None)
    assert reference_by_locator.doi is None
    assert REFERENCE_PAPER in (reference_by_locator.source_locator or "")
    assert circularity_between(calibrated(analyte), reference_by_locator) is not None


# --- the rule, through compare_platforms ---------------------------------------------------


def paired_corpus(count: int, *, circular: bool) -> list[CCSMeasurement]:
    records: list[CCSMeasurement] = []
    for index in range(count):
        analyte = fixtures.small_molecule(fixtures.synthetic_inchikey(f"circpair{index:03d}"))
        reference = the_reference_platform(analyte, ccs=150.0 + 10.0 * index)
        if not circular:
            reference = the_reference_platform(
                analyte,
                ccs=150.0 + 10.0 * index,
                doi=OTHER_PAPER,
                source="another laboratory",
                source_locator="Table 2 of a different study",
            )
        records.append(reference)
        records.append(calibrated(analyte, ccs=(150.0 + 10.0 * index) * 1.02))
    return records


def test_a_circular_pair_never_reaches_a_stratum_and_is_counted():
    """THE ASSERTION THAT MATTERS. Not flagged on a figure: absent from the comparison."""
    matched = build_matched_ions(paired_corpus(12, circular=True))
    assert len(matched.matched) == 12, "the ions must pair, or this proves nothing"
    comparison = compare_platforms(matched)
    assert comparison.pairs == (), "a circular pair reached a stratum"
    assert comparison.n_points == 0
    assert len(comparison.pairs_refused_as_circular) == 12
    assert all("refused as circular" in why for why in comparison.pairs_refused_as_circular)


def test_an_independent_pair_of_the_same_shape_does_reach_a_stratum():
    """The other direction. A comparison that refused every calibrated value would be useless.

    Same platforms, same analytes, same calibration lineage recorded - the only difference is
    that the DTIMS values come from a different publication.
    """
    matched = build_matched_ions(paired_corpus(12, circular=False))
    assert len(matched.matched) == 12
    comparison = compare_platforms(matched)
    assert comparison.pairs_refused_as_circular == ()
    assert comparison.n_points == 12
    assert sum(len(pair.strata) for pair in comparison.pairs) == 1


def test_the_refusal_is_reported_rather_than_silently_thinning_the_corpus():
    """A refusal nobody can count is indistinguishable from an absence of data."""
    comparison = compare_platforms(build_matched_ions(paired_corpus(3, circular=True)))
    assert comparison.ions_considered == 3
    assert len(comparison.pairs_refused_as_circular) == 3
    assert comparison.no_comparison_reason is not None, (
        "no pairs were built, so the report must say why rather than reading as an empty corpus"
    )
