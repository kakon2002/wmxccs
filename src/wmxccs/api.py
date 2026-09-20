"""HTTP API: health, the confidence scheme, and harmonization against a fitted model.

`/harmonize` ANSWERS NOW, AND STILL REFUSES WHAT IT CANNOT DO
------------------------------------------------------------
A model exists as of M4, so this endpoint returns 200 with harmonized values for
measurements it covers. What it does NOT do is extrapolate to cover the rest. Three
outcomes, and which one a caller gets is the useful part of the answer:

- A HARMONIZED VALUE, where the measurement's platform and calibration group match a
  stratum the model applies. It carries the interval, BOTH coverage figures, the grade,
  the scope and the provenance.
- NO VALUE, WITH THE REASON, per measurement, where the model does not cover it: a
  platform it never compared against the primary method, a calibration group it has no
  stratum for, a value already on the primary platform, or a grade of `unsupported`
  meaning the correction would be an assertion rather than an interpolation. The
  original comes back untouched and `not_harmonized_because` says which.
- 501 FOR THE WHOLE REQUEST, where NOTHING in it could be harmonized. A response full of
  absent values is a refusal, and it should carry the status code of one rather than a
  200 that looks like success to anything reading the status alone.

THE 501 IS STILL THE CONTRACT, NOT AN ERROR PATH
------------------------------------------------
It hands back every measurement exactly as sent, the provenance and licence terms of
each, and the maturity stamp. What it never returns is a number: not a zero, not a null
standing in for one, not an interval of infinite width, not a grade computed against
nothing. An API that returned a plausible-looking value with a caveat buried in a field
would be used and the caveat would not be read.

The request is validated BEFORE any of this, so malformed input still gets a 422.

WHAT THE 200 DOES NOT MEAN
--------------------------
The maturity stamp is `provisional` on every response and cannot be otherwise while the
corpus is one study, and every estimate carries a scope saying so. A 200 means "this is
what the model says", not "this is validated". See LIMITATIONS 7F.
"""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI, Request, Response

from . import __version__
from .contracts import (
    ConfidenceReason,
    ConfidenceReport,
    ConfidenceRule,
    ConfidenceRulesResponse,
    HarmonizationUnavailable,
    HarmonizedEstimate,
    HarmonizedMeasurement,
    HarmonizeRequest,
    HarmonizeResponse,
    HealthResponse,
    IntervalKind,
    ModelVersion,
    ScopeReport,
    SourceProvenance,
)
from .grading import ConfidenceGrade, grade_rules
from .harmonization import HarmonizationModel, harmonize as harmonize_one
from .readiness import DataMaturity, MaturityStamp
from .reuse import ReuseStatus, as_reuse_status
from .sources import licence_for

NO_MODEL_LOADED = (
    "No harmonization model is loaded in this process. A model exists and can be fitted from the seed"
    " corpus - see build_default_model - but this application was started without one, so there is nothing"
    " to harmonize against. Your measurements come back untouched with their provenance, and no placeholder"
    " value is returned."
)
NOTHING_IN_THIS_REQUEST_COULD_BE_HARMONIZED = (
    "A model is loaded and it covers none of the {count} measurement(s) in this request. Each one carries"
    " its own reason under `not_harmonized_because`, and the originals come back untouched. This is a 501"
    " rather than a 200 full of absent values, because a response where nothing could be answered is a"
    " refusal and should read as one to anything checking the status code alone."
)

# Left at its old name so nothing importing it breaks, and pointed at the honest text.
# The previous wording said there was not one cross-platform matched ion in this
# repository, which stopped being true on 19 September 2026 and was shipping inside a
# response body.
HARMONIZATION_UNAVAILABLE = NO_MODEL_LOADED

# "Every one can be evaluated today" stood here until 20 September 2026 and was false in the
# way that mattered most: the one rule that cannot be evaluated for a new ion is the only rule
# in the scheme that is ABOUT the submitted ion, and a new ion is what this service is for. A
# published note asserting the opposite is worse than no note.
CONFIDENCE_NOTE = (
    "These are rules, not a fitted model: none needs training data, and a grade only ever falls - the"
    " final grade is the worst demotion found, never a score and never an average. Two mild concerns do"
    " not add up to a severe one, and a severe one is not offset by everything else being fine."
    " NOT EVERY RULE CAN ANSWER FOR EVERY ION. Of the five substantive rules, four are properties of the"
    " calibration group - how many ions it holds, how far your value sits outside its range, whether a"
    " source behind it may be used, how much its slope is levered - and those answer for any ion. The"
    " fifth, the outlier check, is the only one that is about YOUR ion: it is a lookup against ions"
    " already in this corpus, so for an ion this platform has not measured it cannot run AT ALL - not"
    " rarely, not usually, but every time. If you are bringing chemistry we have not measured, you are"
    " graded by four rules, not five, on every request. (The sixth rule published here is about that"
    " situation itself: it is the demotion applied when any rule could not be evaluated.) Each rule's"
    " applies_to says which ions it can answer for; a rule that could not run is named in the response's"
    " not_checked and demotes the grade one notch, and is never counted as passed."
)


def provenance_of(record) -> SourceProvenance:
    """Where one measurement came from and on what terms.

    The licence is read from the REGISTRY by DOI, never from the record. A record may
    claim any reuse status it likes; what the registry holds is what somebody read.

    BOTH ARE REPORTED, AND SO IS WHETHER THEY AGREE. That last part was the gap: this
    docstring already promised the two were "reported together precisely so a caller can
    see when they disagree", and only the registry's LICENCE TEXT was returned, never its
    status - so a response could show a claimed `excluded` beside an open-access licence
    string and a reader had to notice the contradiction unaided. `claim_backed_by_registry`
    now states it.
    """
    entry = licence_for(getattr(record, "doi", None))
    claimed = as_reuse_status(record.reuse_status)
    return SourceProvenance(
        source=str(getattr(record, "source", "")),
        doi=getattr(record, "doi", None),
        reuse_status_claimed=claimed,
        reuse_status_in_registry=entry.reuse_status if entry else None,
        # None where there is no registry entry to agree with. False is the case worth
        # seeing: the record asserts terms nobody recorded.
        claim_backed_by_registry=None if entry is None else entry.reuse_status is claimed,
        licence=entry.licence if entry else None,
        licence_reported_by=entry.reported_by if entry else None,
        licence_reported_on=entry.reported_on if entry else None,
        registered=entry is not None,
    )


SEED_DIRECTORY = Path(__file__).resolve().parents[2] / "data" / "seed"


def build_default_model(directory: Path | None = None) -> HarmonizationModel | None:
    """Fit a model from every seed file on disk, or None if none of them yields one.

    THE DEPLOYMENT WIRING, deliberately at the edge. harmonization.py takes a comparison
    report and knows nothing about files; this reads the files. Keeping them apart is what
    lets the model be fitted in a test from three points without touching a disk.

    It loads every CSV in data/seed, pools the CLEARED records only - the licence gate runs
    first, as it must - builds matched ions across all of them, and fits. Returns None
    rather than an empty model where there is nothing to fit, so that /harmonize answers
    501 rather than 200-with-nothing.
    """
    from .loader import load_measurements_file
    from .matching import build_matched_ions
    from .statistics import compare_platforms
    from .harmonization import fit_harmonization

    directory = directory or SEED_DIRECTORY
    if not directory.is_dir():
        return None
    cleared: list[object] = []
    for path in sorted(directory.glob("*.csv")):
        cleared.extend(load_measurements_file(path).cleared)
    if not cleared:
        return None
    comparison = compare_platforms(build_matched_ions(cleared))
    model = fit_harmonization(comparison)
    return model if model.applied else None


def _estimate_of(result, fingerprint) -> HarmonizedEstimate:
    """Convert a harmonization result into the contract's estimate.

    Both coverage figures travel: `interval_coverage` is the nominal level the quantiles
    were taken at, `guaranteed_coverage` is what jackknife+ actually proves. The interval
    kind says which method, so the second cannot be mistaken for the first.
    """
    correction = result.correction
    band = result.interval
    return HarmonizedEstimate(
        ccs=result.harmonized_ccs,
        basis=result.basis,
        slope_derived_ccs=result.slope_derived_ccs,
        median_derived_ccs=result.median_derived_ccs,
        interval_low=band.low,
        interval_high=band.high,
        interval_coverage=band.nominal_coverage,
        interval_kind=IntervalKind.JACKKNIFE_PLUS,
        reference_platform=correction.reference_platform,
        matched_ions_behind_it=correction.n,
        scope=ScopeReport.of(result.scope),
        interval_is_informative=band.interval_is_informative,
        guaranteed_coverage=band.guaranteed_coverage,
        model_version=ModelVersion.of(fingerprint),
    )


def _confidence_of(result) -> ConfidenceReport | None:
    if result.confidence is None:
        return None
    return ConfidenceReport(
        grade=result.confidence.grade,
        reasons=tuple(
            ConfidenceReason(rule=demotion.rule, falls_to=demotion.grade, detail=demotion.detail)
            for demotion in result.confidence.demotions
        ),
        not_checked=tuple(result.confidence.not_checked),
    )


# Worded from the loader's own refusal of the same status, because it is the same objection:
# that status is a declaration the record is not a measurement, and this service will not
# derive a number from something declared not to be data. The loader says "a row in a file is
# a real record"; the equivalent here is that a request is a real request.
SYNTHETIC_IN_A_REQUEST = (
    "a request may not carry synthetic_fixture: that status is a test's declaration of a record built in"
    " code, and a measurement submitted to this endpoint is a real measurement. This is refused for the"
    " integrity of the answer rather than for a licence - a harmonized value derived from a record that"
    " declares itself invented would be indistinguishable, once returned, from one derived from data."
    " The loader refuses the identical claim in a file, and this is the same refusal at the other door."
)

WITHHELD_UNSUPPORTED = (
    "the correction for this measurement grades 'unsupported', which means do not use the number - so it is"
    " not returned. {why}"
)


def _harmonized_measurement(model: HarmonizationModel | None, record) -> HarmonizedMeasurement:
    """One input, harmonized if the model covers it, with the reason if not.

    NOTHING IS EXTRAPOLATED. A grade of `unsupported` withholds the value rather than
    shipping it with a warning, because the grade's own definition is "do not use this
    number" and handing over a number we have just told the caller not to use is the kind
    of contradiction a caller resolves in favour of the number.
    """
    provenance = provenance_of(record)
    # BEFORE the model is consulted: a record declaring itself invented gets no number
    # whether or not a model is loaded, because the objection is to the record rather than
    # to the state of the service.
    if provenance.reuse_status_claimed is ReuseStatus.SYNTHETIC_FIXTURE:
        return HarmonizedMeasurement(
            original=record, provenance=provenance, not_harmonized_because=SYNTHETIC_IN_A_REQUEST
        )
    if model is None:
        return HarmonizedMeasurement(
            original=record, provenance=provenance, not_harmonized_because=NO_MODEL_LOADED
        )
    result = harmonize_one(model, record)
    confidence = _confidence_of(result)
    if not result.was_corrected:
        return HarmonizedMeasurement(
            original=record,
            provenance=provenance,
            confidence=confidence,
            not_harmonized_because=" ".join(result.refusals) or "the model does not cover this measurement",
        )
    if confidence is not None and confidence.grade is ConfidenceGrade.UNSUPPORTED:
        return HarmonizedMeasurement(
            original=record,
            provenance=provenance,
            confidence=confidence,
            not_harmonized_because=WITHHELD_UNSUPPORTED.format(
                why="; ".join(reason.detail for reason in confidence.reasons)
            ),
        )
    return HarmonizedMeasurement(
        original=record,
        provenance=provenance,
        harmonized=_estimate_of(result, model.fingerprint),
        confidence=confidence,
    )


def create_app(model: HarmonizationModel | None = None) -> FastAPI:
    app = FastAPI(
        title="wmxccs",
        version=__version__,
        description=(
            "Cross-platform CCS harmonization. Stores DTIMS, TWIMS, TIMS and cyclic measurements without"
            " merging them, pairs the same ion across platforms, and returns a harmonized CCS with an"
            " interval and a confidence grade ALONGSIDE the originals, never in place of them."
        ),
    )
    # The fitted harmonization model, passed in rather than built here so that a test can
    # serve a three-point model and a deployment can serve the seed corpus, with no
    # difference in the code between them.
    app.state.model = model

    def maturity(request: Request) -> MaturityStamp:
        """Read from the MODEL, not from a mutable flag beside it.

        This used to read `app.state.model_validated`, which meant a response could say
        'validated' while every estimate in it said within-study - and setting that flag
        anywhere would have done it silently. The model derives its own maturity from its
        own scope, and VALIDATED is not reachable while the corpus is one study.
        """
        loaded = request.app.state.model
        if loaded is None:
            return MaturityStamp(data_maturity=DataMaturity.PROVISIONAL, matched_ion_count=0)
        return loaded.maturity

    @app.get("/health", response_model=HealthResponse)
    def health(request: Request) -> HealthResponse:
        loaded = request.app.state.model
        return HealthResponse(
            version=__version__,
            model_loaded=loaded is not None,
            matched_ions_available=0 if loaded is None else loaded.maturity.matched_ion_count,
            model_version=None if loaded is None else ModelVersion.of(loaded.fingerprint),
        )

    @app.get("/confidence/rules", response_model=ConfidenceRulesResponse)
    def confidence_rules() -> ConfidenceRulesResponse:
        """The grading scheme, published as data so it can be challenged."""
        return ConfidenceRulesResponse(
            grades=tuple(grade.value for grade in ConfidenceGrade),
            # BUILT FROM THE DICT, NOT FIELD BY FIELD. Naming the fields here meant a key
            # added to the scheme was silently not served, which is how `applies_to` - the
            # field carrying the outlier rule's scope limit - reached no caller at all.
            # ConfidenceRule forbids extras, so the next such key raises instead.
            rules=tuple(ConfidenceRule(**rule) for rule in grade_rules()),
            note=CONFIDENCE_NOTE,
        )

    @app.post(
        "/harmonize",
        responses={
            200: {"model": HarmonizeResponse, "description": "At least one measurement was harmonized."},
            501: {
                "model": HarmonizationUnavailable,
                "description": "No model is loaded, or the model covers nothing in this request.",
            },
        },
    )
    def harmonize(body: HarmonizeRequest, request: Request, response: Response):
        """Harmonize what the model covers, refuse what it does not, extrapolate nothing.

        TWO RESPONSE SHAPES AND TWO STATUS CODES, declared rather than inferred. They are
        different shapes, and collapsing them would let a caller write code against a
        harmonized value that is sometimes absent without the status saying so.
        """
        loaded = request.app.state.model
        measurements = tuple(_harmonized_measurement(loaded, record) for record in body.measurements)
        if any(measurement.harmonized is not None for measurement in measurements):
            return HarmonizeResponse(measurements=measurements, maturity=maturity(request))
        response.status_code = 501
        detail = (
            NO_MODEL_LOADED
            if loaded is None
            else NOTHING_IN_THIS_REQUEST_COULD_BE_HARMONIZED.format(count=len(measurements))
        )
        return HarmonizationUnavailable(
            detail=detail, measurements=measurements, maturity=maturity(request)
        )

    return app


# The served application. Built with a model fitted from the seed corpus, so that
# `uvicorn wmxccs.api:app` is a working deployment rather than a skeleton. If the seed
# directory is missing or yields no applicable stratum, the model is None and /harmonize
# answers 501 - which is the honest behaviour for an installation with no data.
app = create_app(build_default_model())

__all__ = [
    "HarmonizationUnavailable",
    "HarmonizeResponse",
    "app",
    "build_default_model",
    "create_app",
    "provenance_of",
]
