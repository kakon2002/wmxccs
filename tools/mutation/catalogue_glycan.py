"""Mutations for src/wmxglycan. One entry per behaviour in the glycan layer that is
supposed to have a test behind it.

SEPARATE FROM catalogue.py, AND NOT BECAUSE THE GLYCAN LAYER IS SECOND-CLASS.
The wmxccs catalogue is closed at 282 and its floor test is written against that
number; mixing two packages into one list would make that floor a statement about
two things at once, and a glycan module arriving would look like CCS coverage
growing. runner.all_mutations() concatenates the two, derived rather than
hand-listed, and tests/test_glycan_mutation_catalogue.py asserts the concatenation
holds every entry of both.

WHAT IS DELIBERATELY NOT HERE: THE 1,220 PORTED TESTS.
Project2's catalogue held 154 mutations, and roughly seventy of them anchored into
the modules this repository ported on 26 September 2026 - the enumerator, the
glycan graph, the featuriser, splits.py, evaluation.py, contracts.py. They are NOT
retrofitted here, on the owner's instruction: those modules came across
byte-identical and were mutation-verified where they were written, and re-covering
them costs days the deadline does not have.

That is a real gap and it is recorded in docs/GLYCAN_LIMITATIONS.md rather than
implied by an empty section. The honest reading: the ported modules are covered by
evidence gathered elsewhere against an identical copy, and the code written HERE is
covered by evidence gathered here.

THE ONE PORTED MODULE THAT DOES APPEAR IS enumeration.py, and only at the line this
repository CHANGED. The ordering-note context gate is new code in a ported file -
before the fix, 114 of 167 Hex5HexNAc4Fuc1 candidates published a cited sentence
that was false about 108 of them - so it is covered like any other new guard.

MOST OF THE [N] ENTRIES BELOW GUARD A DEFECT THAT WAS REALLY IN THIS CODE.
The ranker was reviewed adversarially after it was built and before it was
committed, and the review found the hypothesis space closing the world, the rule
check measured and then ignored, a refused composition asserting the corpus held
nothing of it, a self-asserted flag trusted where the level should have been read,
and a licence gate that required its field to exist without requiring it to permit
anything. Each of those is a mutation here, so the fix cannot quietly come undone.

Sections:
  [N] ranking.py      - bands, ties, shares, coverage and the decision
  [Y] attestation.py  - what may count as evidence
  [V] ccs_evidence.py - the provenance guard and the five states
  [O] enumeration.py  - the ordering-note context gate, the one line changed here
"""

from __future__ import annotations

from . import Mutation

GLYCAN = "wmxglycan"

GLYCAN_MUTATIONS: tuple[Mutation, ...] = (
    # --- [N] ranking.py: the unit of ranking, and refusing to order what cannot be ordered -----
    Mutation(
        label="[N] candidates are ranked individually, so isomers no feature separates get an order",
        file="ranking.py",
        package=GLYCAN,
        find="        grouped.setdefault(class_key_for(candidate), []).append(candidate)",
        replace="        grouped.setdefault(candidate.canonical_key, []).append(candidate)",
    ),
    Mutation(
        label="[N] a class key holding a NaN is accepted, so every class silently becomes a singleton",
        file="ranking.py",
        package=GLYCAN,
        find="    if any(isinstance(value, float) and value != value for _name, value in row):",
        replace="    if False:",
    ),
    Mutation(
        label="[N] a class id depends on the order candidates arrived in",
        file="ranking.py",
        package=GLYCAN,
        find="        class_id = min(keys)  # derived from the members, so it cannot depend on enumerator order",
        replace="        class_id = next(iter(keys))",
    ),
    Mutation(
        label="[N] attestation counts reference ROWS, so a re-spelled structure outranks the rest",
        file="ranking.py",
        package=GLYCAN,
        find="            attested_structures=len(attested),",
        replace="            attested_structures=sum(index.rows_for(key) for key in attested),",
    ),
    Mutation(
        # RE-ANCHORED AND RELABELLED 26 September 2026, not retired. It used to read
        # "the hypothesis space is the candidate set", which was accurate when the space was
        # classes + missing; the catch-all arm was added by the implementation review, so
        # dropping `missing` alone no longer closes the world and the label would have
        # overstated what the mutation does.
        label="[N] the mass reserved for reference structures nobody enumerated is dropped",
        file="ranking.py",
        package=GLYCAN,
        find="    space = len(counts) + missing + 1",
        replace="    space = len(counts) + 1",
    ),
    Mutation(
        label="[N] the hypothesis space omits the catch-all, so the shares can close the world",
        file="ranking.py",
        package=GLYCAN,
        find="    observations = sum(counts) + missing",
        replace="    space = len(counts) + missing\n    observations = sum(counts) + missing",
    ),
    Mutation(
        label="[N] the mass on structures nobody enumerated is reported as zero",
        file="ranking.py",
        package=GLYCAN,
        find='        "not_enumerated": missing * (1.0 + a) / total,',
        replace='        "not_enumerated": 0.0,',
    ),
    Mutation(
        label="[N] the mass on a structure nobody proposed is reported as zero",
        file="ranking.py",
        package=GLYCAN,
        find='        "not_proposed": a / total,',
        replace='        "not_proposed": 0.0,',
    ),
    Mutation(
        label="[N] an undefined share is published as zero instead of None",
        file="ranking.py",
        package=GLYCAN,
        find='            "by_count": {count: None for count in set(counts) | {1}},',
        replace='            "by_count": {count: 0.0 for count in set(counts) | {1}},',
    ),
    Mutation(
        label="[N] completeness is claimed whenever the corpus happens to miss nothing",
        file="ranking.py",
        package=GLYCAN,
        # RE-ANCHORED 27 September 2026. The property returned a bare `True`; it now reads
        # CORPUS_CAN_ESTABLISH_COMPLETENESS so that this fact and the reachability published by
        # `decision_reachability()` come from one place. The INTENT is unchanged - claim
        # completeness whenever the corpus happens to miss nothing - and the anchor moved with the
        # code rather than being deleted, because a mutation that stops applying is a failure here.
        find="        return not CORPUS_CAN_ESTABLISH_COMPLETENESS",
        replace="        return self.not_enumerated > 0",
    ),
    Mutation(
        label="[N] a candidate set that breaks a curated rule is ranked rather than refused",
        file="ranking.py",
        package=GLYCAN,
        find="    if accounting.rules_violated:",
        replace="    if False:",
    ),
    Mutation(
        label="[N] a single candidate is refused with the wording for a set nothing separates",
        file="ranking.py",
        package=GLYCAN,
        find="    elif len(classes) == 1 and not any(one.is_a_tie for one in classes.values()):",
        replace="    elif False:",
    ),
    Mutation(
        label="[N] a set whose classes all carry the same evidence is served as a ranking",
        file="ranking.py",
        package=GLYCAN,
        find="    elif len(bands) == 1:",
        replace="    elif False:",
    ),
    Mutation(
        label="[N] rules violated is assumed to be zero instead of measured by running the check",
        file="ranking.py",
        package=GLYCAN,
        find="        violated += len(working.broken_rules(graph))",
        replace="        violated += 0",
    ),
    Mutation(
        label="[N] a refused composition claims the corpus holds nothing of it",
        file="ranking.py",
        package=GLYCAN,
        find="    reference_keys = index.structures_of_composition(composition_text)\n    coverage = Coverage(",
        replace="    reference_keys = frozenset()\n    coverage = Coverage(",
    ),
    Mutation(
        label="[N] the evidence level is not read, so a shared measurement discriminates",
        file="ranking.py",
        package=GLYCAN,
        find="        and evidence.level is EvidenceLevel.STRUCTURE\n        and evidence.discriminates_between_candidates",
        replace="        and evidence.discriminates_between_candidates",
    ),
    Mutation(
        label="[N] a band reports its per-class share as its total",
        file="ranking.py",
        package=GLYCAN,
        find='                else policy["by_count"][count] * len(group)',
        replace='                else policy["by_count"][count]',
    ),
    Mutation(
        label="[N] the candidates mapping keeps the enumerator's depth-first order",
        file="ranking.py",
        package=GLYCAN,
        find="            for c in sorted(result.candidates, key=lambda one: one.canonical_key)",
        replace="            for c in result.candidates",
    ),
    Mutation(
        label="[N] a truncated candidate set is published without saying so",
        file="ranking.py",
        package=GLYCAN,
        find="        enumerator_capped=result.capped,",
        replace="        enumerator_capped=False,",
    ),
    Mutation(
        label="[N] a lookup that raises is reported as an absence of evidence",
        file="ranking.py",
        package=GLYCAN,
        find="            state=CCSEvidenceState.LOOKUP_FAILED,",
        replace="            state=CCSEvidenceState.NONE_IN_CORPUS_SEARCHED,",
    ),
    Mutation(
        label="[N] a ranker with no evidence source claims it looked and found nothing",
        file="ranking.py",
        package=GLYCAN,
        find="            state=CCSEvidenceState.NOT_CONSULTED, composition=composition.canonical",
        replace=(
            "            state=CCSEvidenceState.NONE_IN_CORPUS_SEARCHED,"
            " composition=composition.canonical,\n"
            '            corpus_searched="nothing was searched", records_consulted=0'
        ),
    ),
    Mutation(
        label="[N] a decision rule fires and the decision does not follow it",
        file="ranking.py",
        package=GLYCAN,
        find="    if any(rule.fires for rule in rules):",
        replace="    if False:",
    ),
    Mutation(
        label="[N] the platform claims a validated model exists, opening AI_ONLY",
        file="ranking.py",
        package=GLYCAN,
        find="    return VALIDATED_MODEL is not None",
        replace="    return True",
    ),
    # --- [Y] attestation.py: what may count as evidence ----------------------------------------
    Mutation(
        label="[Y] an index that loaded nothing is returned, so every candidate scores zero",
        file="attestation.py",
        package=GLYCAN,
        find="        if self.structures_indexed <= 0 or not self.rows:",
        replace="        if False:",
    ),
    Mutation(
        label="[Y] one structure spelled twice attests twice",
        file="attestation.py",
        package=GLYCAN,
        find="        return 1 if canonical_key in self.rows else 0",
        replace="        return self.rows.get(canonical_key, 0)",
    ),
    Mutation(
        label="[Y] a record that declares itself resolved while its string is not is believed",
        file="attestation.py",
        package=GLYCAN,
        find="        if (record.has_unresolved_linkage, record.has_unresolved_anomericity) != (",
        replace="        if False and (record.has_unresolved_linkage, record.has_unresolved_anomericity) != (",
    ),
    Mutation(
        label="[Y] an unresolved reference structure is allowed to attest a linkage assignment",
        file="attestation.py",
        package=GLYCAN,
        find="        if parsed_linkage or parsed_anomer:",
        replace="        if False:",
    ),
    # --- [V] ccs_evidence.py: the provenance guard and the five states -------------------------
    Mutation(
        label="[V] a reference value with no stated uncertainty type is shown",
        file="ccs_evidence.py",
        package=GLYCAN,
        find='        if self.uncertainty_type.strip().casefold() in {"", "unknown", "none"}:',
        replace="        if False:",
    ),
    Mutation(
        label="[V] a reference value whose gas the source never stated is shown",
        file="ccs_evidence.py",
        package=GLYCAN,
        find='        if self.drift_gas.strip().casefold() in {"", "unstated", "unknown"}:',
        replace="        if False:",
    ),
    Mutation(
        label="[V] a reference under a status with no permitted use is shown",
        file="ccs_evidence.py",
        package=GLYCAN,
        find="        if not permitted:",
        replace="        if False:",
    ),
    Mutation(
        label="[V] composition-level evidence may claim to discriminate between candidates",
        file="ccs_evidence.py",
        package=GLYCAN,
        find="        if self.discriminates_between_candidates and self.level is not EvidenceLevel.STRUCTURE:",
        replace="        if False:",
    ),
    Mutation(
        label="[V] composition-level evidence is attached to an individual candidate",
        file="ccs_evidence.py",
        package=GLYCAN,
        find="        if self.level is not EvidenceLevel.STRUCTURE:",
        replace="        if False:",
    ),
    Mutation(
        label="[V] a searched absence is reported without naming the corpus or the count",
        file="ccs_evidence.py",
        package=GLYCAN,
        # RE-ANCHORED: the condition grew two clauses when the zero-denominator hole was
        # closed, so the old anchor went STALE and the sweep said so.
        find="            or self.records_consulted <= 0",
        replace="            or False",
    ),
    Mutation(
        label="[V] a state that is not a measured reference carries a reference value anyway",
        file="ccs_evidence.py",
        package=GLYCAN,
        find="        if state is not CCSEvidenceState.MEASURED_REFERENCE and self.reference is not None:",
        replace="        if False:",
    ),
    # --- [S] store.py: frozen means frozen, and appending never merges --------------------------
    Mutation(
        label="[S] a second write to a frozen prediction is merged instead of refused",
        file="store.py",
        package=GLYCAN,
        find="            raise AlreadyFrozen(ALREADY_FROZEN.format(prediction_id=standing)) from clash",
        replace="            return row",
    ),
    Mutation(
        label="[S] the payload-moved guard stops noticing that a prediction was touched",
        file="store.py",
        package=GLYCAN,
        find="        if after.payload_digest != before.payload_digest:",
        replace="        if False:",
    ),
    Mutation(
        label="[S] a validation forgets which prediction digest it was attached to",
        file="store.py",
        package=GLYCAN,
        find="            predicted_digest=before.payload_digest,",
        replace='            predicted_digest="",',
    ),
    Mutation(
        label="[S] foreign keys stay off, so a validation may reference a prediction that is gone",
        file="store.py",
        package=GLYCAN,
        find='        self._db.execute("PRAGMA foreign_keys = ON")',
        replace='        self._db.execute("PRAGMA foreign_keys = OFF")',
    ),
    Mutation(
        label="[S] a run filter is interpolated into the query instead of parameterised",
        file="store.py",
        package=GLYCAN,
        find='                where.append(f"{column} = ?")\n                values.append(value)',
        replace='                where.append(f"{column} = \'{value}\'")',
    ),
    Mutation(
        label="[S] the run total is the size of the page rather than of the result",
        file="store.py",
        package=GLYCAN,
        find="        return len(self.runs(limit=1_000_000, offset=0, **filters))",
        replace="        return len(self.runs(limit=50, offset=0, **filters))",
    ),
    # --- [G] fingerprint.py: what identifies a pipeline -----------------------------------------
    Mutation(
        label="[G] the corpus digest ignores the structures, so two corpora look identical",
        file="fingerprint.py",
        package=GLYCAN,
        find="            keys_digest=_digest(sorted(index.rows)),",
        replace='            keys_digest="",',
    ),
    Mutation(
        label="[G] the prior policy is left out of the pipeline identity",
        file="fingerprint.py",
        package=GLYCAN,
        find='        "prior_pseudocount": PRIOR_PSEUDOCOUNT,',
        replace='        "prior_pseudocount": 0,',
    ),
    Mutation(
        label="[G] the feature columns are left out of the pipeline identity",
        file="fingerprint.py",
        package=GLYCAN,
        find='        "features": list(FEATURE_NAMES),',
        replace='        "features": [],',
    ),
    Mutation(
        label="[G] the curated rules are left out of the pipeline identity",
        file="fingerprint.py",
        package=GLYCAN,
        find="            for rule in enumerator.constraints\n        ),",
        replace="            for rule in ()\n        ),",
    ),
    # --- [A] api.py: the endpoints ---------------------------------------------------------------
    Mutation(
        label="[A] a repeated client reference creates a second prediction instead of 409",
        file="api.py",
        package=GLYCAN,
        find="            if standing is not None:",
        replace="            if False:",
    ),
    Mutation(
        label="[A] a candidate's attestation is read off its class, so partly attested reads as all",
        file="api.py",
        package=GLYCAN,
        find="                attested=index.attests(key),",
        replace="                attested=result.classes[class_of[key]].attested_structures > 0,",
    ),
    Mutation(
        label="[A] the attach response claims the prediction was unchanged without checking",
        file="api.py",
        package=GLYCAN,
        find="            prediction_unchanged=before == after.payload_digest,",
        replace="            prediction_unchanged=True,",
    ),
    Mutation(
        label="[A] a single measurement is reported as a spread of zero",
        file="api.py",
        package=GLYCAN,
        find="        if len(group) < 2:",
        replace="        if False:",
    ),
    Mutation(
        label="[A] reading a prediction recomputes it instead of serving the frozen bytes",
        file="api.py",
        package=GLYCAN,
        find="        return _response_from(row, request.app.state.fingerprint)\n\n    # --- 3. attach",
        replace=(
            "        return _response_from(\n"
            "            dataclasses.replace(\n"
            "                row, payload=row.payload.replace('\"candidates_total\":', '\"x\":', 1)\n"
            "            ),\n"
            "            request.app.state.fingerprint,\n"
            "        )\n\n    # --- 3. attach"
        ),
    ),
    Mutation(
        label="[A] a failing evidence lookup is reported as an absence rather than a failure",
        file="api.py",
        package=GLYCAN,
        find="                state=CCSEvidenceState.LOOKUP_FAILED,",
        replace="                state=CCSEvidenceState.NOT_CONSULTED,",
    ),
    Mutation(
        label="[A] the decision is not re-derived from the evidence that was frozen with it",
        file="api.py",
        package=GLYCAN,
        find="        result = _with_evidence(result, evidence)",
        replace="        result = dataclasses.replace(result, ccs_evidence=evidence)",
    ),
    Mutation(
        label="[A] a held reference value reaches the comparison response",
        file="api.py",
        package=GLYCAN,
        find="        if evidence.state is CCSEvidenceState.MEASURED_REFERENCE and evidence.reference is not None",
        replace="        if evidence.reference is not None",
    ),
    Mutation(
        label="[A] a missing dashboard is served as an empty page rather than an error",
        file="api.py",
        package=GLYCAN,
        find="        if not DASHBOARD.is_file():  # pragma: no cover - only if the package data is missing",
        replace="        if False:",
    ),
    # --- [P] prediction.py: what a caller may attach ---------------------------------------------
    Mutation(
        label="[P] a measurement with no uncertainty type is accepted",
        file="prediction.py",
        package=GLYCAN,
        find='        if self.uncertainty_type.strip().casefold() in {"", "unknown", "none", "n/a"}:',
        replace="        if False:",
    ),
    Mutation(
        label="[P] a measurement against an unstated gas is accepted",
        file="prediction.py",
        package=GLYCAN,
        find='        if self.drift_gas.strip().casefold() in {"", "unstated", "unknown", "n/a"}:',
        replace="        if False:",
    ),
    Mutation(
        label="[P] a request for charge zero is accepted as an ion",
        file="prediction.py",
        package=GLYCAN,
        find='            raise ValueError("charge 0 is not an ion; give a signed charge such as 1 or -1")',
        replace="            pass",
    ),
    # --- [O] enumeration.py: the one line this repository changed -------------------------------
    Mutation(
        label="[O] an ordering caveat is published for a candidate that lacks the context it names",
        file="enumeration.py",
        package=GLYCAN,
        find="            if is_order_constraint(rule) and any(",
        replace="            if is_order_constraint(rule) and True or any(",
    ),
    # --- [R] the four claims moved out of documents and into the product, 27 September 2026 -----
    #
    # Each of these is a fact a reader would otherwise have to be told by a document: the shares
    # do not close, two decisions are unreachable, and no measured cross section reaches any
    # candidate. A fact served by the API is only as good as the arithmetic behind it, so each is
    # broken here and a test has to notice.
    Mutation(
        label="[R] the shares accounting claims to close when it does not",
        file="api.py",
        package=GLYCAN,
        find="    return abs(total + result.share_not_enumerated + result.share_not_proposed - 1.0) < 1e-9",
        replace="    return True",
    ),
    Mutation(
        label="[R] the candidate share sum silently drops a band",
        file="api.py",
        package=GLYCAN,
        find="    return float(sum(shares)) if shares else None",
        replace="    return float(sum(shares[1:])) if shares else None",
    ),
    Mutation(
        label="[R] the mass off the candidates is reported as the mass on them",
        file="api.py",
        package=GLYCAN,
        find="    return None if total is None else 1.0 - total",
        replace="    return None if total is None else total",
    ),
    Mutation(
        label="[R] the explanation of why the shares fall short is dropped from the response",
        file="api.py",
        package=GLYCAN,
        find="            why_the_shares_do_not_sum_to_one=SHARES_DO_NOT_SUM_TO_ONE,",
        replace='            why_the_shares_do_not_sum_to_one="",',
    ),
    Mutation(
        label="[R] AI_ONLY is published as reachable today",
        file="ranking.py",
        package=GLYCAN,
        find="            reachable_today=not ai_only_blockers,",
        replace="            reachable_today=True,",
    ),
    Mutation(
        label="[R] an unreachable decision is published with nothing that blocks it",
        file="ranking.py",
        package=GLYCAN,
        find="            blocked_by=ai_only_blockers,",
        replace="            blocked_by=(),",
    ),
    # RE-ANCHORED AND RELABELLED 27 September 2026, same day it was written. It targeted
    # `recommended_blockers = (no_completeness,) if not completeness else ()`, a line that existed
    # because the derivation wrongly held that the completeness gate alone blocked
    # IM_VALIDATION_RECOMMENDED. An adversarial read showed no gate state reaches that value at all,
    # the line went away with the correction, and the mutation follows the intent to where the claim
    # now lives: the distinction between "a gate is shut" and "no input can reach this".
    Mutation(
        label="[R] an unreachable-under-any-gate decision is published as merely gated",
        file="ranking.py",
        package=GLYCAN,
        find="            unreachable_under_any_gate=True,",
        replace="            unreachable_under_any_gate=False,",
    ),
    Mutation(
        label="[R] the fall-through is blamed on a data gate rather than on the rule set",
        file="ranking.py",
        package=GLYCAN,
        find="            blocked_by=(the_rule_set,),",
        replace="            blocked_by=(no_completeness,),",
    ),
    Mutation(
        label="[R] the served reachability goes back to a hand-written list",
        file="api.py",
        package=GLYCAN,
        find="            decision_reachable_today=tuple(\n                one.decision.value for one in decision_reachability() if one.reachable_today\n            ),",
        replace="            decision_reachable_today=(Decision.IM_VALIDATION_REQUIRED.value,),",
    ),
    Mutation(
        label="[R] the routes out of an unreachable decision are dropped from the response",
        file="api.py",
        package=GLYCAN,
        find="                one.decision.value: one.what_would_reach_it or \"\"",
        replace='                one.decision.value: ""',
    ),
    Mutation(
        label="[R] the service claims to hold a measured cross section for a candidate",
        file="prediction.py",
        package=GLYCAN,
        find="    holds_a_measured_cross_section_for_any_candidate: bool = False",
        replace="    holds_a_measured_cross_section_for_any_candidate: bool = True",
    ),
)
