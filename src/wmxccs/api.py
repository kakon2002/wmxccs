"""HTTP API: health, the confidence scheme, and a harmonize endpoint that refuses.

THREE ENDPOINTS, AND ONLY ONE OF THEM REFUSES
---------------------------------------------
`/health` and `/confidence/rules` answer today, because neither needs a model.
`/harmonize` answers 501 and will keep answering 501 until a harmonization model
has been fitted on real cross-platform matched ions, of which this repository
holds none.

THE 501 IS NOT AN ERROR PATH, IT IS THE CONTRACT
------------------------------------------------
A caller posting measurements gets back a complete, useful response: every
measurement exactly as they sent it, the provenance and licence terms of each, the
maturity stamp, and a statement of why there is no harmonized value. What they do
not get is a number. Not a zero, not a null standing in for one, not an interval
of infinite width, not a grade computed against nothing.

That distinction is the whole reason this endpoint exists in this state. An API
that returned a plausible-looking value with a caveat buried in a field would be
used, and the caveat would not be read. An API that returns 501 cannot be used by
accident.

The request is validated BEFORE the refusal, so malformed input still gets a 422.
A 501 that swallowed a bad request would tell a caller their measurement was fine
when it was not, and they would find out when the model arrived.
"""

from __future__ import annotations

from fastapi import FastAPI, Request

from . import __version__
from .contracts import (
    ConfidenceRule,
    ConfidenceRulesResponse,
    HarmonizationUnavailable,
    HarmonizedMeasurement,
    HarmonizeRequest,
    HarmonizeResponse,
    HealthResponse,
    SourceProvenance,
)
from .grading import ConfidenceGrade, grade_rules
from .readiness import DataMaturity, MaturityStamp
from .reuse import as_reuse_status
from .sources import licence_for

HARMONIZATION_UNAVAILABLE = (
    "No harmonization model has been fitted, because there is not one cross-platform matched ion in this"
    " repository to fit it on: every measurement held is travelling-wave, from one laboratory. This"
    " endpoint answers 501 until a model fitted on licence-cleared matched ions, with a calibrated"
    " interval, is loaded. It never returns a placeholder value, and it returns your measurements"
    " untouched so that nothing is lost by the refusal."
)

CONFIDENCE_NOTE = (
    "These are rules, not a fitted model. Every one can be evaluated today, none needs training data, and"
    " a grade only ever falls: the final grade is the worst demotion found, never a score and never an"
    " average. Two mild concerns do not add up to a severe one, and a severe one is not offset by"
    " everything else being fine."
)


def provenance_of(record) -> SourceProvenance:
    """Where one measurement came from and on what terms.

    The licence is read from the REGISTRY by DOI, never from the record. A record
    may claim any reuse status it likes; what the registry holds is what somebody
    read, and the two are reported together precisely so a caller can see when
    they disagree.
    """
    entry = licence_for(getattr(record, "doi", None))
    return SourceProvenance(
        source=str(getattr(record, "source", "")),
        doi=getattr(record, "doi", None),
        reuse_status=as_reuse_status(record.reuse_status),
        licence=entry.licence if entry else None,
        licence_reported_by=entry.reported_by if entry else None,
        licence_reported_on=entry.reported_on if entry else None,
        registered=entry is not None,
    )


def create_app() -> FastAPI:
    app = FastAPI(
        title="wmxccs",
        version=__version__,
        description=(
            "Cross-platform CCS harmonization. Stores DTIMS, TWIMS, TIMS and cyclic measurements without"
            " merging them, pairs the same ion across platforms, and returns a harmonized CCS with an"
            " interval and a confidence grade ALONGSIDE the originals, never in place of them."
        ),
    )
    # The fitted harmonization model. None until one exists, and nothing sets it.
    app.state.model = None
    # Cross-platform matched ions the loaded model rests on, and whether it has been
    # checked against held-out ions. Both stay at these values until a licence-cleared
    # matched-ion set exists, so every output is stamped provisional.
    app.state.matched_ion_count = 0
    app.state.model_validated = False

    def maturity(request: Request) -> MaturityStamp:
        return MaturityStamp(
            data_maturity=(
                DataMaturity.VALIDATED if request.app.state.model_validated else DataMaturity.PROVISIONAL
            ),
            matched_ion_count=request.app.state.matched_ion_count,
        )

    @app.get("/health", response_model=HealthResponse)
    def health(request: Request) -> HealthResponse:
        return HealthResponse(
            version=__version__,
            model_loaded=request.app.state.model is not None,
            matched_ions_available=request.app.state.matched_ion_count,
        )

    @app.get("/confidence/rules", response_model=ConfidenceRulesResponse)
    def confidence_rules() -> ConfidenceRulesResponse:
        """The grading scheme, published as data so it can be challenged."""
        return ConfidenceRulesResponse(
            grades=tuple(grade.value for grade in ConfidenceGrade),
            rules=tuple(
                ConfidenceRule(
                    rule=rule["rule"],
                    falls_to=tuple(rule["falls_to"]),
                    why=rule["why"],
                    threshold=rule["threshold"],
                    basis=rule["basis"],
                )
                for rule in grade_rules()
            ),
            note=CONFIDENCE_NOTE,
        )

    @app.post("/harmonize", status_code=501, response_model=HarmonizationUnavailable)
    def harmonize(body: HarmonizeRequest, request: Request) -> HarmonizationUnavailable:
        """Refuses, and hands back everything it was given.

        The response model is HarmonizationUnavailable rather than HarmonizeResponse
        because the two are different shapes and pretending otherwise would let a
        caller write code against a harmonized value that is never there. The 200
        shape is specified in contracts.HarmonizeResponse for callers building
        ahead of the model.
        """
        return HarmonizationUnavailable(
            detail=HARMONIZATION_UNAVAILABLE,
            measurements=tuple(
                # harmonized and confidence are left absent, not empty. There is no
                # number to put in them and no grade to compute against nothing.
                HarmonizedMeasurement(original=record, provenance=provenance_of(record))
                for record in body.measurements
            ),
            maturity=maturity(request),
        )

    return app


app = create_app()

# Named so a reader of contracts.py can find where the 200 shape is used, which is
# nowhere yet, deliberately.
__all__ = ["HarmonizeResponse", "app", "create_app", "provenance_of"]
