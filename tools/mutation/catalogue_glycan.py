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
        find='        constant, and `completeness` carries the distinction that actually varies.\n        """\n        return True',
        replace='        constant, and `completeness` carries the distinction that actually varies.\n        """\n        return self.not_enumerated > 0',
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
    # --- [O] enumeration.py: the one line this repository changed -------------------------------
    Mutation(
        label="[O] an ordering caveat is published for a candidate that lacks the context it names",
        file="enumeration.py",
        package=GLYCAN,
        find="            if is_order_constraint(rule) and any(",
        replace="            if is_order_constraint(rule) and True or any(",
    ),
)
