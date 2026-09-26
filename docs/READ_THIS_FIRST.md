# Read this first

**This repository holds two things at two very different maturities, and the second one does not
predict cross sections.**

`src/wmxccs` is a cross-platform CCS harmonization pipeline with a model fitted on 142
cross-platform matched ions. `src/wmxglycan` is a glycan and isomer layer that enumerates the
candidate N-glycan structures for a composition and bands them by the evidence that attests them.
It has **no fitted model of any kind** and V1 will not have one. Every cross section it reports is
a measured value or an explicit statement that none is held.

Two packages, one repository, no shared code. They are at `0.8.0` and `0.1.0` and the difference
between those numbers is the point.

---

This document exists because the seven things below are what a reader is most likely to get wrong
in their first hour, ranked by how expensive the mistake is. Each one is a place where the honest
answer is less than the obvious reading. None of them is a bug, and four of them have been
mistaken for bugs already.

Where a point is enforced rather than described, the enforcing file is named. Prefer it to this
document: a paragraph can go stale and a test cannot.

## 1. The two maturities, and the wall between them

The repository is called `wmxccs`, the glycan layer lives inside it, and the specification this
was built against asks for cross-section prediction. The natural reading is "this predicts glycan
cross sections". **It does not.**

What the CCS core actually is: 142 matched ions, all from one study
(`DOI 10.1021/jasms.2c00196`). Eighteen strata, of which only the nine anchored on stepped-field
DTIMS are ever applied. Its maturity is `provisional` and cannot be anything else while the corpus
is one study — the maturity is derived from the data's provenance, not declared. The prediction
interval is guaranteed at 80 per cent and not 90, and measured coverage is arithmetically pinned
so it is not evidence of anything.

What the glycan layer is: an enumerator, an attestation index over 3,640 distinct reference
structures, a ranker that refuses to order what nothing separates, six endpoints and a dashboard.
`training.py` is a refuse-to-fit guard in which every path to a fit raises, and it raises on an
absence of records to fit on rather than on unfinished code.

The wall is enforced by `tests/test_glycan_boundary.py`, not by this paragraph. The two packages
never import each other in either direction. `tools/glycan_ccs_evidence.py` is the only production
file that imports both, and it is in neither package — which is what keeps the wall up, because the
traffic is carried from outside rather than by either side reaching across. Tests that check the
relationship between the two packages load both as well, necessarily — you cannot check a wall
without loading what is on either side of it — so what
`test_exactly_one_production_file_imports_both_packages` counts is the crossings outside the suite.

(README said `tools/glycan_service.py` was that one file until 27 September 2026. It was wrong:
that file imports `wmxglycan` and the adapter, and nothing from `wmxccs`. Corrected in place.)

## 2. A band is not a ranking, and "Band 1" is not the answer

For `Hex5HexNAc4Fuc1` the platform returns 167 candidates in 61 indistinguishable classes and
three bands. **152 of those 167 candidates are tied with at least one other candidate**, and the
largest tied class holds ten. Band 1 is one class of four candidates resting on four of the 34
reference structures known for that composition.

A band is an evidence tier. It is not a rank, and the candidates inside one are not ordered
because nothing available can order them: linkage-position isomers on an identical shape have
byte-identical feature rows. **No candidate is numbered anywhere in this product**, and that is
deliberate — a list from 1 to 167 would present an arbitrary order as a ranking, which is the
failure this project was built against.

The thing to be careful about is the word "band" doing work a hurried reader will not do. "Band 1
of 3" reads like a podium. It means "the tier with the most attestation", where the most is four
deposited structures out of 34.

Also: the curated rules **cannot order a candidate set at all.** All 15 govern MGAT branching and
bisecting, and not one of them constrains what the candidates for a single composition differ in.
Attestation in the reference corpus carries the whole discriminating load, and where the corpus
attests nothing, a tie is the honest output.

## 3. The shares do not sum to 100%, and a version that did would be the bug

Add up the per-class shares for `Hex5HexNAc4Fuc1` and you get about **69%**. Nothing is missing.

The shares are over a hypothesis space, not over the candidate list. That space holds every
candidate shown, **plus** the reference structures of the composition that the enumerator did not
propose, **plus** one catch-all for a structure neither proposed here nor deposited anywhere. Only
the first group comes back to you, so the candidate shares are below 1 by exactly the mass sitting
on the other two. Those three do account to 1:

| | Hex5HexNAc4Fuc1 |
|---|---|
| on the candidates shown | 69.03% |
| on reference structures not enumerated | 30.09% |
| on a structure nobody proposed | 0.88% |

The catch-all is never zero, so this can never close. An earlier version reported shares summing
to exactly 1.0 for a composition with a single deposited structure, which amounted to asserting
that the answer was in the candidate set on the strength of that one deposition.

**You do not have to know any of that to find it out.** The API response carries
`candidate_shares_sum`, `mass_not_on_any_candidate`, a `shares_account_to_one` identity checked
against the actual numbers, and `why_the_shares_do_not_sum_to_one` — all in the same `confidence`
object as the figure itself. The dashboard's Confidence card says the same thing next to the same
numbers.

One trap inside the trap: `mass_on_unattested_classes` (41.59% for G2F) is a **subset** of the
mass on candidates shown, not a fourth term. It counts candidates that are listed but attested by
nothing. Adding it to the three above exceeds 100%, and both the field description and the
dashboard say so where the number appears.

Last thing about the number: it is never called a probability. The served field is `confidence`,
the prior is a declared **policy** (`uniform_1`, pseudocount 1.0), all four standard
non-informative priors travel in the response, and the choice of prior is the single largest lever
in the output. Calibration is `never_calibrated`, and what would calibrate it is stated in the
response rather than left as an exercise.

## 4. The decision field is constant, and the two values you cannot reach are not stubs

Every response says `IM_VALIDATION_REQUIRED`. `AI_ONLY` and `IM_VALIDATION_RECOMMENDED` never
appear, and someone who greps the enum, finds two values nothing produces, and concludes the field
is half-written would be making the most understandable mistake in this repository.

It is constant because two of the published decision rules fire on **every** input: no model here
has been validated against independently known structures, and the completeness of a candidate set
is not establishable from a corpus that records depositions rather than existence. A platform with
no validated CCS model asking for instrument validation on every answer is the correct output, and
it is precisely what the specification means by identifying predictions that require experimental
validation.

Do not fix it by relaxing a threshold. What would actually reach the other two values is written
at the enum in `src/wmxglycan/ranking.py`, and it is **derived** from the two gates rather than
described beside them — `decision_reachability()` reads the same predicates the decision function
reads, so a gate that opens changes the published reachability in the same commit. The API serves
that derivation as `decision_reachable_today` and `decision_unreachable_today`, the latter naming
what would reach each one; the dashboard's Model Monitor shows both.

Both branches are live code, not dead enum members: `tests/test_glycan_ranking.py` patches each
gate and asserts the other values then appear. That matters because "AI_ONLY never appears" is
equally satisfied by dead code, a misspelled comparison, or an enum member nothing references.

## 5. The measured-reference path is unreachable. Not thin — unreachable

Twenty-four glycan cross sections clear the licence gate. That reads like 24 usable values. **It
is zero.**

All 24 are milk oligosaccharides and not one of them records a composition, so nothing can key a
reference value to a candidate set. The state `MEASURED_REFERENCE` is therefore unreachable for
every candidate set this platform can produce, in this release, without exception. The 89 Struwe
2015 records that *are* N-glycans — including eight for `Hex5HexNAc2` as `[M-H]-` — are held on a
drift gas and an uncertainty type genuinely absent from the source. Reading another paper for the
gas does not release them, and the gate is not being loosened to reach a number.

So the states a real deployment returns are `HELD_NOT_RELEASABLE`, `NONE_IN_CORPUS_SEARCHED` and
`LOOKUP_FAILED`, and every response carries a **stated absence** of CCS evidence rather than
omitting the field. The absence is the output, not a gap in the wiring.

None of that is asserted in prose where it matters: the adapter's `reachable_states()` **counts**
it, reporting `0 of 541 cleared record(s) carry a composition`, and a test fails the moment that
count changes. The claim cannot quietly become false. It is stated at the enum in
`src/wmxglycan/ccs_evidence.py`, served in the API's `domain` block as
`holds_a_measured_cross_section_for_any_candidate: false`, and recorded in
`docs/GLYCAN_LIMITATIONS.md`.

There is a second, larger step behind the first. Even a cleared record carrying a composition
would give a value keyed on composition and ion — shared by every isomer of that composition, and
therefore unable to confirm one of them. Separating candidates needs a measurement that resolves
to one **structure**. Nothing in this repository has one.

## 6. Not everything counted here is mutation-verified here

The suite is 4,105 tests and the mutation catalogue is 339 anchors. Both numbers are real and they
do not mean the same thing, and the second is smaller than the first in a way that matters.

The 339 anchors are 282 over the CCS core, verified at `v0.7.0-mvp`, plus 57 over the glycan code
written in this repository. **The 1,771 ported glycan tests are not covered by this repository's
sweep.** They were mutation-verified in the repository they were ported from, and re-covering them
here was declined on cost with the deadline three days out. That is a decision, not an oversight,
and it is the claim most likely to be repeated without its qualifier.

What a green suite and a green sweep each mean: the suite says the tests ran, the sweep says they
would have caught something. A mutation that survives is a behaviour with no test behind it and is
treated here as a failure rather than as a note. The sweep of `c145eaf` ran 120 mutations and six
survived, all six in code written that week; six tests were written and each was confirmed to kill
its mutation before `v0.8.0-mvp` was tagged. `LIMITATIONS.md` section 4.5 records it as instance
Thirteen, and it is the entry to read if you only read one — a defect was found, fixed, annotated
with a comment at the site explaining exactly what it had been, and no test was written.

The sweep refuses a working tree with uncommitted changes. It copies `src/` once at the start and
reads `tests/` live, so an edit made mid-run changes the suite between one mutation and the next.
`--dirty` overrides it; a figure measured that way is not a sweep result.

## 7. `Project2` is not a dependency and you cannot fetch it

`docs/GLYCAN_PORT.md` references a repository called `Project2` throughout, against a digest for
every ported file. Fifteen of the glycan modules are byte-identical to their source there.

It is **not** a dependency. Code was copied and adapted, never imported; it is not on this
repository's path; and a checkout anywhere else will not have it. Nothing in the install or the
test suite reaches for it. The digests are there so that a port can be **audited** — so you can
confirm that a file arrived unchanged — not so that anything can be resolved. If a test fails
looking for it, that is a bug in the test and not a missing package.

---

## Where to go next

| | |
|---|---|
| `README.md` | what it is, how to install and run it, and how to check it |
| `docs/HANDOVER.md` | picking the whole thing up cold |
| `LIMITATIONS.md` | the long, authoritative account of what the CCS core does not know. Section 4.5 is the recurring failure class |
| `docs/GLYCAN_LIMITATIONS.md` | the same for the glycan layer, including what the dashboard refuses to show |
| `CLAUDE.md` | the hard constraints and the process rules, each learned the expensive way |
| `CONTEXT.md` | everything established before this repository existed |

`docs/CCS_Harmonization_Plan.md` and `docs/KICKOFF_PROMPT.md` are historical. They record what was
asked in September, they predate the glycan layer, and they are not a description of what is here.
