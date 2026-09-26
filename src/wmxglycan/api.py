"""The six endpoints the spec names. Real services over a frozen store, not a demonstration.

    POST /v1/predictions                    create from a composition plus ion metadata
    GET  /v1/predictions/{id}               the frozen prediction, candidates, confidence, provenance
    POST /v1/predictions/{id}/validation    attach experimental measurements
    GET  /v1/predictions/{id}/comparison    delta CCS, interval coverage, status
    GET  /v1/runs                           search history
    GET  /v1/models/current                 version, fingerprint, data snapshot, domain

THREE RULES THAT SHAPE EVERY ONE OF THEM

1. NOTHING IS EVER OVERWRITTEN. A prediction is frozen on creation and the store has no UPDATE
   statement. A second create under the same client reference answers 409 and names the
   prediction that stands; it does not merge, and it does not quietly make a second one the
   caller will believe is the same. Attaching measurements APPENDS a validation record and the
   response carries the prediction's payload digest read before and after, so "the prediction
   was not touched" is a measurement in the response rather than a promise in a docstring.

2. EVERY RESPONSE CARRIES THE STAMP - package version, both pipeline digests, and the data
   snapshot - so a prediction read back later can be checked against the pipeline that made it.
   The frozen payload is served BYTE FOR BYTE from storage rather than re-rendered from live
   objects, which is what makes a past prediction reconstructible rather than merely re-derived
   under whatever is installed now.

3. THE CCS CORE IS ON THE OTHER SIDE OF A WALL THIS MODULE DOES NOT CROSS. `create_app` takes a
   `CCSEvidenceLookup` by injection and imports nothing from `wmxccs`. The default `app` below is
   built with NO evidence source, so it honestly reports NOT_CONSULTED rather than claiming to
   have looked; `tools/glycan_service.py` is the composition root that wires the two packages
   together for a real deployment, and it is the only thing that imports both.

WHY THE DEFAULT APP HAS NO CCS EVIDENCE SOURCE

Because the alternative is worse. Wiring the adapter in here would mean this module imports
`wmxccs`, which is the one thing the two-package structure exists to prevent; and defaulting to
something that pretends to look would put a scientific claim - "nothing found in the literature" -
behind a deployment nobody configured. NOT_CONSULTED is the honest default and it is
distinguishable from an absence, which is the whole point of the five-state evidence enum.
"""

from __future__ import annotations

import dataclasses
from collections import defaultdict
from pathlib import Path
from typing import Mapping

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse

from . import __version__
from .attestation import AttestationIndex, default_attestation_index
from .ccs_evidence import CCSEvidence, CCSEvidenceLookup, CCSEvidenceState
from .composition import Composition, CompositionError, parse_composition
from .enumeration import Enumerator
from .fingerprint import PipelineFingerprint, default_fingerprint
from .prediction import (
    AGREEMENT_LIMIT_PERCENT,
    INTERVAL_UNEVALUABLE,
    NO_PREDICTED_VALUE,
    AttachedSpread,
    BandOut,
    CandidateOut,
    ClassOut,
    ComparisonResponse,
    ComparisonState,
    ConfidenceOut,
    CoverageOut,
    CurrentModelResponse,
    DecisionOut,
    DecisionRuleOut,
    DomainOut,
    MeasurementComparison,
    ModelStamp,
    PredictionRequest,
    PredictionResponse,
    ReasonOut,
    RuleAccountingOut,
    RunOut,
    RunsResponse,
    ValidationRequest,
    ValidationResponse,
)
from .prediction import MeasurementIn
from .ranking import (
    Calibration,
    Decision,
    RankedSet,
    decide,
    decision_reachability,
    rank,
)
from .store import AlreadyFrozen, Store, StoredPrediction, new_id, now

DEFAULT_DATABASE = Path(__file__).resolve().parents[2] / "data" / "glycan_service.sqlite3"
# The dashboard ships inside the package so that an installed deployment serves the same
# bytes a checkout does. It is one self-contained file - no build step, no bundler, no CDN -
# and it is served from THIS service so it is same-origin: no CORS configuration to get
# wrong, and no second process to keep in step with the API it is a client of.
DASHBOARD = Path(__file__).resolve().parent / "static" / "dashboard.html"

COMPOSITION_REFUSED = (
    "the composition {given!r} could not be read: {problem}. Nothing is partly parsed here - a"
    " composition that does not parse completely is refused, because a partly read composition"
    " would silently become a different molecule"
)
NO_SUCH_PREDICTION = (
    "no prediction {prediction_id!r} exists. Predictions are never deleted, so an id that is not"
    " here was never issued by this service"
)
NOTHING_ATTACHED = (
    "no experimental measurements have been attached to prediction {prediction_id!r}, so there is"
    " nothing to compare. POST to /v1/predictions/{prediction_id}/validation first"
)

# Served beside the number, because a reader who adds the candidate shares up needs the answer
# where they are looking rather than in a document they would have to know to open.
SHARES_DO_NOT_SUM_TO_ONE = (
    "They cannot, and a version of this that did would be the bug. The shares are over a"
    " HYPOTHESIS SPACE, not over the candidate list: it holds every candidate shown, plus the"
    " reference structures of this composition that the enumerator did not propose, plus one"
    " catch-all for a structure neither proposed here nor deposited anywhere. Only the first group"
    " is returned to you, so `candidate_shares_sum` is below 1 by exactly the mass sitting on the"
    " other two - and those three DO account to 1, which `shares_account_to_one` reports as a"
    " checked identity. The catch-all is never zero, so this can never close. An earlier version"
    " reported shares that summed to exactly 1.0 for a composition with one deposited structure,"
    " which asserted the answer was in the set on the strength of that one deposition."
    " `mass_on_unattested_classes` is NOT a fourth term: it is already inside"
    " `candidate_shares_sum` and counts candidates that are shown but attested by nothing."
)


def _band_shares(result: RankedSet) -> tuple[float, ...]:
    return tuple(
        band.band_total_share for band in result.bands if band.band_total_share is not None
    )


def _candidate_shares_sum(result: RankedSet) -> float | None:
    """Every share served over the candidates, added up, so nobody has to add them up."""
    shares = _band_shares(result)
    return float(sum(shares)) if shares else None


def _mass_off_the_candidates(result: RankedSet) -> float | None:
    total = _candidate_shares_sum(result)
    return None if total is None else 1.0 - total


def _shares_account_to_one(result: RankedSet) -> bool | None:
    """The identity, CHECKED against this result rather than asserted about the design.

    A promise in a docstring would stay true-looking after the arithmetic stopped holding. This
    returns False in that case, and a test asserts it is True for every composition it ranks - so
    a change that breaks the accounting is a failing test rather than a wrong number in a field.
    """
    total = _candidate_shares_sum(result)
    if total is None or result.share_not_enumerated is None or result.share_not_proposed is None:
        return None
    return abs(total + result.share_not_enumerated + result.share_not_proposed - 1.0) < 1e-9


# --- rendering a ranked set into the served shapes -------------------------------------------------


def _reasons(candidate) -> tuple[ReasonOut, ...]:
    return tuple(
        ReasonOut(
            kind=reason.kind,
            text=reason.text,
            enzymes=tuple(reason.enzymes),
            reference=reason.reference,
        )
        for reason in candidate.rationale
    )


def _render(
    result: RankedSet, fingerprint: PipelineFingerprint, index: AttestationIndex
) -> dict:
    """The body of a prediction, as a plain dict, ready to be frozen and served.

    A DICT AND NOT A MODEL, deliberately. This is what gets stored, and it is stored as the JSON
    that was served rather than as objects that a later schema change could re-render
    differently. The response model validates it on the way out, so the shape is still checked.
    """
    class_of = {key: one.class_id for one in result.classes.values() for key in one.members}
    return {
        "candidates_total": result.candidates_total,
        "classes_total": result.classes_total,
        "tied_candidates": result.tied_candidates,
        "largest_indistinguishable_class": result.largest_tie,
        "enumerator_built": result.enumerator_built,
        "enumerator_rejected": result.enumerator_rejected,
        "enumerator_truncated": result.enumerator_capped,
        "is_a_ranking": result.is_a_ranking,
        "refusal": result.refusal,
        "weaknesses": list(result.weaknesses),
        "bands": [
            BandOut(
                rank=band.rank,
                classes=tuple(sorted(band.classes)),
                candidates=band.candidates,
                attested_structures=band.attested_structures,
                share_per_class=band.share_per_class,
                band_total_share=band.band_total_share,
                share_under_priors=dict(band.share_under_priors),
            ).model_dump()
            for band in result.bands
        ],
        "classes": [
            ClassOut(
                class_id=one.class_id,
                members=one.rendered_members,
                is_a_tie=one.is_a_tie,
                attested_structures=one.attested_structures,
                attesting_accessions=one.attesting_accessions,
                reference_rows=one.reference_rows,
                evidence_share=one.evidence_share,
                evidence_share_under_priors=dict(one.share_under_priors),
            ).model_dump()
            for one in sorted(result.classes.values(), key=lambda c: c.class_id)
        ],
        "candidates": [
            CandidateOut(
                canonical_key=key,
                iupac_condensed=candidate.iupac_condensed,
                class_id=class_of[key],
                # PER CANDIDATE, from the index, and NOT read off the class. The first version
                # reported `class.attested_structures > 0`, which marks every member of a partly
                # attested class as attested - so a class of 10 holding 4 attested structures
                # would have reported 10 attested candidates. The class-level count stays on the
                # class, where it means what it says.
                attested=index.attests(key),
                attesting_accessions=index.accessions_for(key),
                reference_rows=index.rows_for(key),
                rationale=_reasons(candidate),
                rules_applicable=result.rule_accounting.rules_applicable.get(key, 0),
                ordering_caveats=result.rule_accounting.ordering_caveats.get(key, 0),
            ).model_dump()
            for key, candidate in sorted(result.candidates.items())
        ],
        "confidence": ConfidenceOut(
            prior=result.prior,
            prior_pseudocount=result.prior_pseudocount,
            priors=tuple(sorted(result.mass_under_priors)),
            calibration=result.calibration,
            what_would_calibrate_it=result.what_would_calibrate_it,
            what_the_number_means=result.share_means,
            candidate_shares_sum=_candidate_shares_sum(result),
            mass_not_on_any_candidate=_mass_off_the_candidates(result),
            mass_on_reference_structures_not_enumerated=result.share_not_enumerated,
            mass_on_a_structure_nobody_proposed=result.share_not_proposed,
            mass_on_unattested_classes=result.share_on_unattested_classes,
            shares_account_to_one=_shares_account_to_one(result),
            why_the_shares_do_not_sum_to_one=SHARES_DO_NOT_SUM_TO_ONE,
            mass_under_priors={
                name: dict(masses) for name, masses in result.mass_under_priors.items()
            },
        ).model_dump(),
        "coverage": CoverageOut(
            composition=result.coverage.composition,
            state=result.coverage.state,
            completeness=result.coverage.completeness,
            truth_may_not_be_in_the_candidate_set=(
                result.coverage.truth_may_not_be_in_the_candidate_set
            ),
            reference_rows=result.coverage.reference_rows,
            reference_structures=result.coverage.reference_structures,
            attested_by_a_candidate=result.coverage.attested_by_a_candidate,
            not_enumerated=result.coverage.not_enumerated,
            not_enumerated_examples=result.coverage.not_enumerated_examples,
            summary=result.coverage.summary(),
        ).model_dump(),
        "rules": RuleAccountingOut(
            rules_in_scheme=result.rule_accounting.rules_in_scheme,
            rules_violated=result.rule_accounting.rules_violated,
            verified_by_running_the_check=(
                result.rule_accounting.rules_violated_verified_by_running_the_check
            ),
            rules_applicable_range=tuple(
                sorted(set(result.rule_accounting.rules_applicable.values()))
            ),
            summary=result.rule_accounting.summary(),
        ).model_dump(),
        "decision": DecisionOut(
            decision=result.decision,
            rules=tuple(
                DecisionRuleOut(
                    name=rule.name,
                    fires=rule.fires,
                    applies_to=rule.applies_to,
                    because=rule.because,
                    decision_if_it_fires=rule.decision_if_it_fires,
                )
                for rule in result.decision_rules
            ),
        ).model_dump(),
        "ccs_evidence": _evidence_body(result.ccs_evidence),
    }


def _evidence_body(evidence: CCSEvidence) -> dict:
    """The CCS evidence, flattened for the response, with its own summary sentence."""
    body = evidence.model_dump(mode="json")
    body["summary"] = evidence.summary()
    return body


def _stamp(fingerprint: PipelineFingerprint) -> ModelStamp:
    return ModelStamp.of(fingerprint)


def _response_from(row: StoredPrediction, fingerprint: PipelineFingerprint) -> PredictionResponse:
    """Build the response from the STORED BYTES, never from a fresh computation.

    This is what makes a frozen prediction frozen in practice rather than in principle: reading
    it back does not re-run the enumerator, so a change to the tables or the corpus cannot alter
    what a past prediction says. The stamp carried in the body is the one it was made under.
    """
    body = row.body
    return PredictionResponse(
        prediction_id=row.prediction_id,
        created_at=row.created_at,
        composition=row.composition,
        adduct=row.adduct,
        charge=row.charge,
        client_reference=row.client_reference,
        stamp=ModelStamp(**body["stamp"]),
        **{key: value for key, value in body.items() if key != "stamp"},
    )


# --- the comparison ---------------------------------------------------------------------------------


def _compare(
    row: StoredPrediction,
    measurements: tuple[MeasurementIn, ...],
    evidence: CCSEvidence,
    validations: int,
    fingerprint: PipelineFingerprint,
) -> ComparisonResponse:
    """Delta CCS, interval coverage and status, each with the reason it is what it is.

    THREE AXES, kept apart because collapsing them would lose the only one that works:

      against the prediction   UNDEFINED, always. There is no predicted cross section, so a delta
                               is not unmeasured, it does not exist.
      against a reference      a real delta where the platform holds a releasable reference for
                               this ion, and the evidence state's own reason where it does not.
      among the attached       the caller's own measurements against each other. This needs no
                               model and no reference, so it always works, and it is the caller's
                               reproducibility - which is worth reporting rather than omitting
                               because the other two axes are unevaluable.
    """
    releasable = (
        evidence.reference
        if evidence.state is CCSEvidenceState.MEASURED_REFERENCE and evidence.reference is not None
        else None
    )
    state_reason = {
        CCSEvidenceState.HELD_NOT_RELEASABLE: ComparisonState.REFERENCE_HELD_NOT_RELEASABLE,
        CCSEvidenceState.NONE_IN_CORPUS_SEARCHED: ComparisonState.NO_REFERENCE_IN_CORPUS,
        CCSEvidenceState.NOT_CONSULTED: ComparisonState.REFERENCE_NOT_CONSULTED,
        CCSEvidenceState.LOOKUP_FAILED: ComparisonState.REFERENCE_NOT_CONSULTED,
    }.get(evidence.state, ComparisonState.NO_REFERENCE_IN_CORPUS)

    per_measurement: list[MeasurementComparison] = []
    for measured in measurements:
        if releasable is not None:
            delta = measured.ccs - releasable.ccs
            percent = 100.0 * delta / releasable.ccs if releasable.ccs else None
            per_measurement.append(
                MeasurementComparison(
                    ccs=measured.ccs,
                    uncertainty=measured.uncertainty,
                    uncertainty_type=measured.uncertainty_type,
                    adduct=measured.adduct,
                    charge=measured.charge,
                    ims_type=measured.ims_type,
                    drift_gas=measured.drift_gas,
                    source=measured.source,
                    against_reference=ComparisonState.COMPARED,
                    reference_ccs=releasable.ccs,
                    delta_ccs=delta,
                    delta_percent=percent,
                    within_agreement_limit=(
                        None if percent is None else abs(percent) <= AGREEMENT_LIMIT_PERCENT
                    ),
                )
            )
            continue
        per_measurement.append(
            MeasurementComparison(
                ccs=measured.ccs,
                uncertainty=measured.uncertainty,
                uncertainty_type=measured.uncertainty_type,
                adduct=measured.adduct,
                charge=measured.charge,
                ims_type=measured.ims_type,
                drift_gas=measured.drift_gas,
                source=measured.source,
                against_reference=state_reason,
                unevaluable_because=evidence.summary(),
            )
        )

    # The caller's own measurements against each other, grouped by the ion and the conditions
    # that make two values comparable at all. Grouping on the gas as well as the platform is
    # hard constraint 4: values measured against different gases are not interchangeable.
    grouped: Mapping[tuple, list[MeasurementIn]] = defaultdict(list)
    for measured in measurements:
        grouped[
            (measured.adduct, measured.charge, measured.ims_type, measured.drift_gas)
        ].append(measured)
    spreads: list[AttachedSpread] = []
    for (adduct, charge, ims_type, gas), group in sorted(grouped.items(), key=lambda kv: kv[0]):
        if len(group) < 2:
            continue
        values = [one.ccs for one in group]
        lowest, highest = min(values), max(values)
        spread = highest - lowest
        midpoint = (highest + lowest) / 2
        percent = 100.0 * spread / midpoint if midpoint else 0.0
        spreads.append(
            AttachedSpread(
                adduct=adduct,
                charge=charge,
                ims_type=ims_type,
                drift_gas=gas,
                measurements=len(group),
                lowest=lowest,
                highest=highest,
                spread=spread,
                spread_percent=percent,
                within_agreement_limit=percent <= AGREEMENT_LIMIT_PERCENT,
            )
        )

    compared = [one for one in per_measurement if one.against_reference is ComparisonState.COMPARED]
    return ComparisonResponse(
        prediction_id=row.prediction_id,
        composition=row.composition,
        adduct=row.adduct,
        charge=row.charge,
        validations=validations,
        measurements=len(measurements),
        status=ComparisonState.COMPARED if compared else ComparisonState.NO_PREDICTED_VALUE,
        against_prediction=ComparisonState.NO_PREDICTED_VALUE,
        against_prediction_because=NO_PREDICTED_VALUE,
        interval_coverage=None,
        interval_coverage_because=INTERVAL_UNEVALUABLE,
        per_measurement=tuple(per_measurement),
        among_attached=tuple(spreads),
        stamp=_stamp(fingerprint),
    )


# --- the application ---------------------------------------------------------------------------------


def create_app(
    *,
    store: Store | None = None,
    evidence: CCSEvidenceLookup | None = None,
    enumerator: Enumerator | None = None,
    index: AttestationIndex | None = None,
    fingerprint: PipelineFingerprint | None = None,
) -> FastAPI:
    """Build the service. Everything is injected so a test serves the same code a deployment does.

    `evidence` is the only door to the CCS core and it is a Protocol, so this module imports
    nothing from `wmxccs`. Passing None is not a broken configuration: it means no evidence
    source was wired, and every response says NOT_CONSULTED rather than claiming an absence.
    """
    app = FastAPI(
        title="wmxglycan",
        version=__version__,
        description=(
            "Glycan and isomer prediction. Ranks candidate N-glycan structures for a composition"
            " on the evidence that attests them, refuses to order a set nothing separates, and"
            " reports cross sections only as measured values - never as predictions, because"
            " there is no glycan CCS model."
        ),
    )
    app.state.store = Store() if store is None else store
    app.state.evidence = evidence
    app.state.enumerator = Enumerator() if enumerator is None else enumerator
    app.state.index = default_attestation_index() if index is None else index
    app.state.fingerprint = default_fingerprint() if fingerprint is None else fingerprint

    def stamp(request: Request) -> ModelStamp:
        return _stamp(request.app.state.fingerprint)

    def look_up(request: Request, composition: Composition, adduct: str, charge: int) -> CCSEvidence:
        """Ask the wired evidence source, and never let a failure become an absence."""
        source = request.app.state.evidence
        if source is None:
            return CCSEvidence(
                state=CCSEvidenceState.NOT_CONSULTED, composition=composition.canonical
            )
        try:
            return source.evidence_for(composition, adduct, charge)
        except Exception as failure:  # noqa: BLE001 - any failure, and it must not read as absence
            return CCSEvidence(
                state=CCSEvidenceState.LOOKUP_FAILED,
                composition=composition.canonical,
                adduct=adduct,
                charge=charge,
                failure=f"{type(failure).__name__}: {failure}",
                key_attempted=f"{composition.canonical} {adduct} {charge:+d}",
            )

    # --- 1. create ---------------------------------------------------------------------------------

    @app.post(
        "/v1/predictions",
        status_code=201,
        response_model=PredictionResponse,
        responses={
            409: {"description": "A prediction already exists under this client_reference."},
            422: {"description": "The composition or the ion could not be read."},
        },
    )
    def create_prediction(body: PredictionRequest, request: Request) -> PredictionResponse:
        """Rank the candidates for a composition and ion, freeze the answer, and return it.

        FROZEN ON CREATION rather than after a draft stage, because the answer is a pure function
        of the request and the pipeline: there is no state in which it is half-decided.
        """
        try:
            composition = parse_composition(body.composition)
        except (CompositionError, TypeError) as problem:
            raise HTTPException(
                status_code=422,
                detail=COMPOSITION_REFUSED.format(given=body.composition, problem=problem),
            ) from problem

        store: Store = request.app.state.store
        if body.client_reference is not None:
            standing = store.by_client_reference(body.client_reference)
            if standing is not None:
                # 409 AND THE ID THAT STANDS. Not a merge, and not a silent second prediction the
                # caller would believe is the same one.
                raise HTTPException(
                    status_code=409,
                    detail={
                        "error": "already_frozen",
                        "prediction_id": standing.prediction_id,
                        "detail": (
                            f"a prediction is already frozen under client_reference"
                            f" {body.client_reference!r}. It was NOT modified and no second"
                            " prediction was created. Fetch it, or create with a different"
                            " reference"
                        ),
                    },
                )

        fingerprint: PipelineFingerprint = request.app.state.fingerprint
        evidence = look_up(request, composition, body.adduct, body.charge)
        result = rank(
            request.app.state.enumerator.enumerate(composition),
            index=request.app.state.index,
            enumerator=request.app.state.enumerator,
            ccs=None,
            adduct=None,
            charge=None,
        )
        # The evidence is looked up once, here, and carried into the frozen body: a prediction
        # must not change because a later read consults a corpus that has moved.
        result = _with_evidence(result, evidence)

        body_dict = _render(result, fingerprint, request.app.state.index)
        body_dict["stamp"] = stamp(request).model_dump()
        prediction_id = new_id("pred")
        created_at = now()
        try:
            row = store.freeze(
                prediction_id=prediction_id,
                created_at=created_at,
                composition=composition.canonical,
                adduct=body.adduct,
                charge=body.charge,
                pipeline_version=fingerprint.version,
                fingerprint=fingerprint.short,
                snapshot_digest=fingerprint.snapshot.digest,
                payload=body_dict,
                client_reference=body.client_reference,
            )
        except AlreadyFrozen as clash:
            raise HTTPException(status_code=409, detail=str(clash)) from clash
        store.record_run(
            run_id=new_id("run"),
            created_at=created_at,
            kind="prediction.created",
            prediction_id=row.prediction_id,
            composition=row.composition,
            fingerprint=fingerprint.short,
            snapshot_digest=fingerprint.snapshot.digest,
            detail=(
                f"{result.candidates_total} candidates in {result.classes_total} classes,"
                f" {len(result.bands)} band(s), {result.decision.value}"
                + ("" if result.is_a_ranking else ", RANKING REFUSED")
            ),
        )
        return _response_from(row, fingerprint)

    # --- 2. read -----------------------------------------------------------------------------------

    @app.get(
        "/v1/predictions/{prediction_id}",
        response_model=PredictionResponse,
        responses={404: {"description": "No such prediction."}},
    )
    def get_prediction(prediction_id: str, request: Request) -> PredictionResponse:
        """The frozen prediction, served from the stored bytes and never recomputed."""
        row = request.app.state.store.get(prediction_id)
        if row is None:
            raise HTTPException(
                status_code=404, detail=NO_SUCH_PREDICTION.format(prediction_id=prediction_id)
            )
        return _response_from(row, request.app.state.fingerprint)

    # --- 3. attach ---------------------------------------------------------------------------------

    @app.post(
        "/v1/predictions/{prediction_id}/validation",
        status_code=201,
        response_model=ValidationResponse,
        responses={404: {"description": "No such prediction."}},
    )
    def attach_validation(
        prediction_id: str, body: ValidationRequest, request: Request
    ) -> ValidationResponse:
        """Attach experimental measurements. APPEND ONLY; the prediction is not touched.

        The response carries the prediction's payload digest read before and after the write, so
        a caller does not have to trust that nothing was merged - the evidence is in the reply.
        """
        store: Store = request.app.state.store
        row = store.get(prediction_id)
        if row is None:
            raise HTTPException(
                status_code=404, detail=NO_SUCH_PREDICTION.format(prediction_id=prediction_id)
            )
        fingerprint: PipelineFingerprint = request.app.state.fingerprint
        before = row.payload_digest
        created_at = now()
        record = store.attach(
            prediction_id=prediction_id,
            validation_id=new_id("val"),
            created_at=created_at,
            fingerprint=fingerprint.short,
            snapshot_digest=fingerprint.snapshot.digest,
            payload={
                "measurements": [one.model_dump(mode="json") for one in body.measurements],
                "note": body.note,
                "stamp": stamp(request).model_dump(),
            },
        )
        after = store.get(prediction_id)
        assert after is not None
        store.record_run(
            run_id=new_id("run"),
            created_at=created_at,
            kind="validation.attached",
            prediction_id=prediction_id,
            composition=row.composition,
            fingerprint=fingerprint.short,
            snapshot_digest=fingerprint.snapshot.digest,
            detail=f"{len(body.measurements)} measurement(s) attached",
        )
        return ValidationResponse(
            validation_id=record.validation_id,
            prediction_id=prediction_id,
            created_at=record.created_at,
            measurements_attached=len(body.measurements),
            validations_on_this_prediction=len(store.validations(prediction_id)),
            prediction_digest_before=before,
            prediction_digest_after=after.payload_digest,
            prediction_unchanged=before == after.payload_digest,
            stamp=stamp(request),
        )

    # --- 4. compare --------------------------------------------------------------------------------

    @app.get(
        "/v1/predictions/{prediction_id}/comparison",
        response_model=ComparisonResponse,
        responses={
            404: {"description": "No such prediction, or nothing has been attached to it."},
        },
    )
    def get_comparison(prediction_id: str, request: Request) -> ComparisonResponse:
        """Delta CCS, interval coverage and status - each carrying the reason it is what it is."""
        store: Store = request.app.state.store
        row = store.get(prediction_id)
        if row is None:
            raise HTTPException(
                status_code=404, detail=NO_SUCH_PREDICTION.format(prediction_id=prediction_id)
            )
        records = store.validations(prediction_id)
        if not records:
            raise HTTPException(
                status_code=404, detail=NOTHING_ATTACHED.format(prediction_id=prediction_id)
            )
        measurements = tuple(
            MeasurementIn(**one)
            for record in records
            for one in record.body.get("measurements", ())
        )
        composition = parse_composition(row.composition)
        evidence = look_up(request, composition, row.adduct, row.charge)
        return _compare(
            row, measurements, evidence, len(records), request.app.state.fingerprint
        )

    # --- 5. history --------------------------------------------------------------------------------

    @app.get("/v1/runs", response_model=RunsResponse)
    def search_runs(
        request: Request,
        kind: str | None = Query(default=None, description="prediction.created or validation.attached"),
        composition: str | None = Query(default=None),
        prediction_id: str | None = Query(default=None),
        since: str | None = Query(default=None, description="ISO-8601; runs at or after this."),
        limit: int = Query(default=50, ge=1, le=500),
        offset: int = Query(default=0, ge=0),
    ) -> RunsResponse:
        """Search the history. Every filter is optional and they compose with AND.

        `total` is the count BEFORE paging, so a caller can page without guessing whether there
        is more. A response that returned only `returned` would make the last page
        indistinguishable from a full one.
        """
        store: Store = request.app.state.store
        # The composition filter is CANONICALISED where it parses, so a caller searching for
        # "HexNAc2Hex5" finds runs stored as "Hex5HexNAc2". Where it does not parse it is passed
        # through unchanged rather than refused: a run log is searched with whatever a caller
        # remembers, and answering "no runs" is more useful than a validation error on a filter.
        filters = {
            "kind": kind,
            "composition": _canonical_or_as_given(composition),
            "prediction_id": prediction_id,
            "since": since,
        }
        applied = {key: value for key, value in filters.items() if value is not None}
        rows = store.runs(limit=limit, offset=offset, **applied)
        return RunsResponse(
            total=store.count_runs(**applied),
            returned=len(rows),
            limit=limit,
            offset=offset,
            filters=applied,
            runs=tuple(
                RunOut(
                    run_id=one.run_id,
                    created_at=one.created_at,
                    kind=one.kind,
                    prediction_id=one.prediction_id,
                    composition=one.composition,
                    fingerprint=one.fingerprint,
                    data_snapshot=one.snapshot_digest,
                    detail=one.detail,
                )
                for one in rows
            ),
            stamp=stamp(request),
        )

    # --- the dashboard, served from the same origin as the endpoints it calls -----------------------

    @app.get("/", response_class=HTMLResponse, include_in_schema=False)
    def dashboard() -> HTMLResponse:
        """The six-page dashboard. A client of the six endpoints below, not a separate rendering.

        Served here rather than from a static file server so that it is same-origin: the page
        calls `/v1/...` with no base URL and no CORS, and there is no second process that could
        drift out of step with the API it is a client of. Read from disk on each request, which
        costs nothing at this size and means editing the file does not need a restart.
        """
        if not DASHBOARD.is_file():  # pragma: no cover - only if the package data is missing
            raise HTTPException(
                status_code=500,
                detail=(
                    f"the dashboard is not installed at {DASHBOARD}. It ships as package data;"
                    " a wheel built without it will serve the endpoints and not the page"
                ),
            )
        return HTMLResponse(DASHBOARD.read_text(encoding="utf-8"))

    # --- 6. what is serving ------------------------------------------------------------------------

    @app.get("/v1/models/current", response_model=CurrentModelResponse)
    def current_model(request: Request) -> CurrentModelResponse:
        """The pipeline identity, the data snapshot, and the applicability domain."""
        enumerator: Enumerator = request.app.state.enumerator
        store: Store = request.app.state.store
        from .composition import Residue

        return CurrentModelResponse(
            stamp=stamp(request),
            calibration=Calibration.NEVER_CALIBRATED,
            what_would_calibrate_it=(
                "a set of compositions whose true structure is independently known, enumerated"
                " and scored blind, with the realised frequency of the true structure landing in"
                " each band compared against that band's share. This repository holds no such"
                " set and nothing in it can construct one"
            ),
            decision_values=tuple(value.value for value in Decision),
            # ACTUALLY DERIVED NOW. This said "DERIVED FROM THE RULES, not listed" above a
            # hand-written one-element tuple, which is the shape LIMITATIONS 4.5 collects: a
            # comment asserting a property the code next to it does not have. `ranking`
            # reads the two gates and answers, so a gate that opens cannot leave this behind.
            decision_reachable_today=tuple(
                one.decision.value for one in decision_reachability() if one.reachable_today
            ),
            decision_unreachable_today={
                one.decision.value: one.what_would_reach_it or ""
                for one in decision_reachability()
                if not one.reachable_today
            },
            domain=DomainOut(
                residues_supported=tuple(residue.value for residue in Residue),
                biosynthetic_rules=len(enumerator.constraints),
            ),
            predictions_frozen=len(store.predictions(limit=1_000_000)),
        )

    return app


def _canonical_or_as_given(text: str | None) -> str | None:
    """A composition filter in canonical form where it parses, unchanged where it does not."""
    if text is None:
        return None
    try:
        return parse_composition(text).canonical
    except Exception:  # noqa: BLE001 - any unparseable filter is passed through, not refused
        return text


def _with_evidence(result: RankedSet, evidence: CCSEvidence) -> RankedSet:
    """The ranked set with its CCS evidence replaced, and the decision re-derived from it.

    `rank()` is called without an evidence source so that the lookup happens once here and the
    answer it produced is frozen into the body. The decision has to be recomputed rather than
    kept, because one of its published rules reads the evidence - keeping the old one would serve
    rules that disagree with the evidence printed beside them.
    """
    decision, rules = decide(
        classes=result.classes,
        bands=result.bands,
        coverage=result.coverage,
        evidence=evidence,
        refused=result.refusal is not None,
    )
    return dataclasses.replace(
        result, ccs_evidence=evidence, decision=decision, decision_rules=rules
    )


# The served application. NO CCS EVIDENCE SOURCE, deliberately: wiring one here would mean this
# module imports wmxccs, and defaulting to something that pretends to look would put a claim
# about the literature behind a deployment nobody configured. `tools/glycan_service.py` is the
# composition root that wires both packages for a real deployment.
app = create_app(store=Store(DEFAULT_DATABASE))

__all__ = ["app", "create_app", "DASHBOARD", "DEFAULT_DATABASE"]
