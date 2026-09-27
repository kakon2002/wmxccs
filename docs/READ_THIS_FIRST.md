# Read this first

**This repository holds two things at two very different maturities, and the second one does not
predict cross sections.**

`src/wmxccs` is a cross-platform CCS harmonization pipeline with a model fitted on 142
cross-platform matched ions. `src/wmxglycan` is a glycan and isomer layer that enumerates the
candidate N-glycan structures for a composition and bands them by the evidence that attests them.
It has **no fitted model of any kind** and V1 will not have one. Every cross section it reports is
a measured value or an explicit statement that none is held.

Two packages, one repository, no shared code. They are at `0.8.1` and `0.1.0` and the difference
between those numbers is the point.

---

This document exists because the seven things below are what a reader is most likely to get wrong
in their first hour, ranked by how expensive the mistake is. Each one is a place where the honest
answer is less than the obvious reading. None of them is a bug, and four of them have been
mistaken for bugs already.

Where a point is enforced rather than described, the enforcing file is named. Prefer it to this
document: a paragraph can go stale and a test cannot.

## What is built

Before the seven: the glycan layer is not a scaffold waiting for a model. These exist, they are
tested, and they are what the platform does.

**The enumerator** (`enumeration.py`) takes a composition and returns the N-glycan structures the
curated biosynthetic rules permit on a complete Man3GlcNAc2 core, **up to a stated tree budget of
60,000 arrangements**. For `Hex5HexNAc4Fuc1` it builds 1,012 trees and returns 167 candidates, well
inside the budget. It reports what it rejected and why.

**The budget is a scope limit and the response always says whether it was hit** —
`enumerator_built`, `enumerator_truncated`, and a refusal that names the budget rather than blaming
the curated rules. Both halves of this paragraph were wrong until 28 September 2026: it claimed
*every* permitted structure and claimed the enumerator refuses rather than truncating. Measured, a
tri-antennary tri-sialylated composition exhausts 60,000 and yields nothing, while
`Hex6HexNAc5NeuAc3` completes at 54,204 trees and returns 11,672 candidates in five bands — so no
budget makes the space complete, and the old default of 5,000 was turning answerable compositions
into refusals that blamed glycobiology.

**The class algebra** (`ranking.py`) groups candidates that no feature the platform computes can
separate. A class is a `frozenset` of candidate keys and a band is a `frozenset` of class ids, so an
order inside either is unrepresentable rather than merely undocumented. Bands are a tuple, because
between bands there is a real ordering. The key is built from the None-bearing feature mapping, not
from `as_row()`, because `as_row()` maps absent values to a NaN that compares equal to itself only
by identity.

**The attestation index** (`attestation.py`) reads SugarBase v12 through glycowork and indexes 4,001
fully-resolved structures as 3,640 distinct ones, over 687 compositions. Resolvedness is derived
from parsing each structure, not read off the record's own flags. It counts distinct structures and
never rows, because 352 keys are spelled more than once in the source.

**The store** (`store.py`) freezes a prediction on creation and appends to it afterwards. It
contains no `UPDATE` and no `DELETE` statement, asserted by a test that reads the module, and no
mutating HTTP verb is exposed anywhere. A second write to a frozen prediction returns 409 naming the
record that stands. Attaching measurements returns the payload digest read before and after the
write, so a caller does not have to take "nothing was merged" on trust.

**The six endpoints** (`api.py`) create a prediction, read it back, attach measurements, compare,
list the run history, and report the current model. Every response carries the pipeline fingerprint
and the data snapshot it was made against, so a prediction read in a year can be checked against the
pipeline that produced it.

**The dashboard** (`static/dashboard.html`) is a client of those six endpoints, served by the same
process at `/`. Six pages, no build step, no CDN. It shows bands, tied groups and how many candidates
share each position, and it numbers no candidate anywhere.

Behind them: a test suite, and a mutation harness that breaks each guard one at a time and requires a
test to notice. Both counts are in `README.md` § Checking it, where they are derived rather than
written down. **This file states no live TEST, MUTATION or ANCHOR count** — that is the rule
`tests/test_repository_counts.py` enforces, and it is item 6. It does carry other measured figures,
candidate counts and corpus sizes among them, and those are checked by
`tests/test_glycan_served_claims.py` and the ranking suite rather than by that guard.

**The refusals are the deliverable, not an unfinished state.** A platform that says these four
candidates are indistinguishable and will not guess between them is doing the job it was built for;
the alternative is a ranked list whose order is arbitrary, which is the failure this project was
built against. Read the seven below as the boundaries of a working thing, not as a list of gaps.

## 1. The two maturities, and the wall between them

The repository is called `wmxccs`, the glycan layer lives inside it, and the specification this
was built against asks for cross-section prediction. The natural reading is "this predicts glycan
cross sections". **It does not.**

What the CCS core actually is: 142 matched ions, all from one study
(`DOI 10.1021/jasms.2c00196`). Eighteen strata, of which only the nine anchored on stepped-field
DTIMS are ever applied. Its maturity is `provisional` and cannot be anything else while the corpus
is one study. The *scope* stamp is derived from the data's provenance and re-derived by a
validator, so an interlaboratory claim is unrepresentable; the maturity LEVEL is a declared constant
(`harmonization.py` returns `DataMaturity.PROVISIONAL` as a literal, under a docstring saying
"Always provisional while the corpus is one study"). This said the maturity was derived rather than
declared, which credited the wrong field — a second study would not move it. The prediction
interval is guaranteed at 80 per cent and not 90, and measured coverage is arithmetically pinned
so it is not evidence of anything.

What the glycan layer is: an enumerator, an attestation index over 3,640 distinct reference
structures, a ranker that refuses to order what nothing separates, six endpoints and a dashboard.
`training.py` is a refuse-to-fit guard in which every path to a fit raises. With today's corpus it
raises on an absence of records to fit on — not on a licence problem, and the error class says which
— though past the readiness gate the terminal raise is a `NotImplementedError`, because no model has
ever been trained here.

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

Also: the curated rules **cannot order a candidate set at all.** Eleven of the fifteen reach
N-glycan enumeration — nine MGAT rules on branching order and bisecting interference, FUT8 on core
fucosylation, one class-agnostic blood group rule; the other four are O-glycan core rules — and
**not one of the eleven constrains galactosylation type, fucose position, chain extension or
LacdiNAc**, which is what the candidates for one composition differ in.

That is the served field's wording (`domain.rules_govern`) and it is deliberate. This said "not one
of the eleven constrains what the candidates for a single composition differ in", which is stronger
and false: the number of curated rules bearing on a candidate varies — measured over the 167
`Hex5HexNAc4Fuc1` candidates, 1 rule on 41 of them, 2 on 31, 3 on 50, 4 on 39, 5 on 6 — and FUT8's
ordering caveat separates 6 candidates from the other 161. Neither quantity is used as a score, and
`ranking.py` says why in the code: *"a field documented as constant is a field a reader is told not
to check."* The document was telling a reader not to check a quantity the code publishes per
candidate.
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

Do not fix it by relaxing a threshold. What would reach the other two values is written at the enum
in `src/wmxglycan/ranking.py`. `decision_reachability()` reads the same two gates the decision
function reads, so a gate that opens changes the published reachability in the same commit, and the
API serves that as `decision_reachable_today` and `decision_unreachable_today` — the latter naming
what would reach each blocked value. The dashboard's Model Monitor shows both.

**The two are unreachable in two different ways, and the difference matters more than the fact.**

`AI_ONLY` is live code behind two gates. Open both — a validated model, and a way to establish
completeness — supply structure-level evidence that discriminates between candidates, and the value
appears. `tests/test_glycan_ranking.py` does exactly that, which is what shows the branch is real;
"AI_ONLY never appears" is otherwise equally satisfied by dead code, a misspelled comparison, or an
enum member nothing references.

`IM_VALIDATION_RECOMMENDED` is different: **no gate state and no input reaches it.** Its `return` is
a fall-through, and reaching a fall-through means silencing every rule — but the rule "no cross
section is held for this structure" fires whenever the evidence is not structure-level, which is the
same condition the `AI_ONLY` branch above it tests. So whenever every rule is silent, `AI_ONLY` is
returned first. Only a change to the published rules would reach it.

This paragraph said something else until it was checked. It claimed both branches had positive
controls and that the completeness gate alone blocked the second, and the code and its test agreed
with each other because I wrote both from the same wrong reading of the decision function. A probe
over both gates × three evidence shapes × one and two classes × refused-or-not settled it: 48
combinations, 47 `IM_VALIDATION_REQUIRED`, 1 `AI_ONLY`, **0 `IM_VALIDATION_RECOMMENDED`**. That
probe is now the test, so the claim is checked against the decision function rather than against
this document. It is the clearest example in the repository of the failure `LIMITATIONS.md` § 4.5
collects, and it was found by an adversarial read rather than by the suite.

## 5. The measured-reference path is unreachable. Not thin — unreachable

Twenty-four glycan cross sections clear every gate the loader applies. That reads like 24 usable
values. **It is zero.**

(And it is not the licence gate that does the work, though the shorthand suggests it: all 117
glycan rows are `open_attribution` and the licence refuses none of them. 4 are dropped as a
suspected shared peak, 89 for an unstated gas and an absent uncertainty type. Permission is not
what is missing here — which is worth knowing, because it is the one blocker that could in
principle be resolved by correspondence.)

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
therefore unable to confirm one of them. Separating candidates needs a measurement that resolves to
one **structure** and is keyed to a candidate this platform produces.

Be precise about which half is missing, because the obvious phrasing is wrong: **the repository does
hold 24 structure-resolved measurements.** The milk oligosaccharides are fully linkage- and
anomericity-specified, and each carries an IUPAC string. What it holds none of is structure-level
evidence keyed to a candidate — those 24 are outside the N-glycan candidate space and record no
composition, so no candidate set can reach one. `EvidenceLevel.STRUCTURE` is the term for the thing
that is absent, and it is absent because of the keying, not because nothing was resolved.

## 6. Not everything counted here is mutation-verified here

The suite and the mutation catalogue both have counts, they are both real, and they do not mean the
same thing. **No CURRENT count is repeated here** — a figure tied to a named commit is allowed and
one is used below, which is what the guard enforces. They are in `README.md` § Checking it, where
`tests/test_repository_counts.py` derives them from pytest and the catalogue and fails if the file
disagrees. Three of the numbers that used to sit in this paragraph went stale within one commit of
being written, which is the whole argument for not writing them here.

The relationship, which does not change when the counts do: the anchors split into those over the
CCS core, verified at `v0.7.0-mvp`, and those over the glycan code written in this repository.
**The ported glycan tests are not covered by this repository's sweep at all.** They were
mutation-verified in the repository they came from, and re-covering them here was declined on cost
with the deadline three days out. `docs/GLYCAN_PORT.md` § 2 lists exactly which fifteen MODULES
came across, with a digest each. **There is no list of the ported TEST files** — four are named in
passing and the count is stated, but nothing enumerates them, and the prefix `test_glycan_` does not
identify them either, because most files carrying it were written here. So a
sentence like "every test in this repository is mutation-verified" would be false, and it is the
claim most likely to be repeated without its qualifier.

What a green suite and a green sweep each mean: the suite says the tests ran, the sweep says they
would have caught something. A mutation that survives is a behaviour with no test behind it and is
treated here as a failure rather than as a note. The sweep of `c145eaf` ran 120 mutations and six
survived, all six in code written that week; six tests were written and each was confirmed to kill
its mutation before `v0.8.0-mvp` was tagged. `LIMITATIONS.md` section 4.5 records it as instance
Thirteen: a defect was found, fixed, annotated with a comment at the site explaining exactly what it
had been, and no test was written.

**And neither gate catches everything.** Instance Fourteen in the same section is the one to read
first, because it is the only entry there that no amount of testing would have found — the test and
the code it guarded were written in the same hour from the same wrong reading, so they agreed with
each other and not with the function. *A suite checks that the code does what its author believed.
It does not check the belief.* What found it was an independent reader, which is why that is now
treated here as a category of tool rather than a nicety.

The sweep refuses a working tree with uncommitted changes. It copies `src/` once at the start and
reads `tests/` live, so an edit made mid-run changes the suite between one mutation and the next.
`--dirty` overrides it; a figure measured that way is not a sweep result.

## 7. `Project2` is not a dependency and you cannot fetch it

`docs/GLYCAN_PORT.md` references a repository called `Project2` throughout, against a digest for
every ported module. **Fourteen of the fifteen** are byte-identical to their source there;
`enumeration.py` was deliberately adapted, and GLYCAN_PORT.md records the original's digest so the
edit can be told from a bad copy. (This document said "fifteen" until it was checked against
GLYCAN_PORT.md, which says fourteen. The `v0.8.0-mvp` tag message says fifteen too — see the note
at the end of this file.) The digests cover the ported MODULES; the ported test files carry none.

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

---

## Three figures in the `v0.8.0-mvp` tag message are wrong

The tag was written and pushed before the claims in this document were checked against the code, and
the check found errors the tag repeats. A pushed annotated tag is not rewritten here without being
asked, so the corrections live in the repository instead. If you read the tag, read this:

| the tag says | the truth |
|---|---|
| "15 of them are byte-identical to the repository they came from" | **Fourteen of the fifteen.** `enumeration.py` was deliberately adapted, and `docs/GLYCAN_PORT.md` § 2 has always said fourteen |
| "All 15 govern MGAT branching and bisecting" | **Eleven reach N-glycan enumeration** — nine MGAT, plus FUT8 on core fucosylation, plus one class-agnostic blood group rule. Four are O-glycan core rules. The conclusion the sentence draws is unaffected |
| "`tests/test_glycan_ranking.py` patches both gates and asserts the other two values then appear" | Only `AI_ONLY` appears. `IM_VALIDATION_RECOMMENDED` is unreachable under any gate state — see item 4 |

None of the three changes a served number, a fingerprint, or a test result. All three are
overstatements — of how much the curated rules cover, of how clean the port was, and of how much the
suite proves — and they are recorded here rather than quietly fixed, because a tag is meant to be a
fixed record of what was claimed at a release.

**How they were found** is the part worth carrying forward. Seven read-only agents were pointed at
this document with one instruction: check every factual claim against the repository and report only
what is false, unsupported, or stale. They checked 127 claims and found several real errors,
including the `IM_VALIDATION_RECOMMENDED` one — which the full suite did not catch and could not have
caught, because the test asserting it was written from the same wrong reading of the code as the code
itself. A suite checks that the code does what its author believed. It does not check the belief.
