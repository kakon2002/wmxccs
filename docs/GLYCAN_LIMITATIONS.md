# Glycan layer: known limitations

For `src/wmxglycan`. The CCS core has its own, far longer, register in
[`LIMITATIONS.md`](../LIMITATIONS.md); this file does not repeat it and the two are not
interchangeable. What ported and what did not is in
[`GLYCAN_PORT.md`](GLYCAN_PORT.md).

Started 26 September 2026, at the port. Numbers here were measured in **this** repository
against the installed `glycowork 1.10.0`; where a figure is carried over from `Project2`
without being re-measured it says so.

---

## 1. The guard class travelled with the code

**The single most useful thing the port produced, and it is a finding about the class rather
than about the bug.** Recorded in full as instance **Twelve** of
[`LIMITATIONS.md` § 4.5](../LIMITATIONS.md), "THE RECURRING CLASS: a guard that looks tested
and is not", and named here because it is the glycan layer's first entry in that register and
because of where it came from.

Fourteen test files were copied verbatim from `Project2` and run against the copy before
anything was adapted. **1,164 passed and one failed** — a guard that importing the package
loads no estimator library, defeated by a hand-written list of eleven module names that
included `wmxglycan.api`, which this repository deliberately does not carry.

Three things make it worth this much space:

- **It is the fourth occurrence of one shape.** Instances Ten, Eleven and Twelve are all a
  hand-maintained list of what a guard covers, going stale. The shape survives being written
  by careful people to a written standard.
- **The direction it failed in was the harmless one.** A module went missing, so it was loud.
  The silent direction is a module being *added* and the list not growing with it: green
  suite, guard covering less than it appears to, nobody looking.
- **The fix needed its own floor.** The list is now derived from the package directory, and
  deriving it moves the failure mode rather than removing it — an empty glob imports nothing
  and leaves the guard green. `test_the_derived_module_list_actually_covers_the_package`
  asserts the derivation is not empty.

**The lesson for this layer specifically: expect this class in ported code and go looking.**
This one surfaced by luck. Nothing about a green copied suite says its guards still cover
what they covered where they were written.

## 2. No glycan CCS model exists, and V1 will not have one

**The owner's decision of 26 September 2026, and it shapes everything downstream.** Recorded
as a decision rather than a defect.

What the repository actually holds for glycan CCS, measured:

| seed file | read | cleared | held | what they are |
|---|---|---|---|---|
| `struwe2016_chemcommun.csv` | 28 | **24** | 4 | milk oligosaccharides (LNH, LNnH and relatives) |
| `struwe2015_analyst.csv` | 89 | **0** | 89 | N-glycans, high mannose: Hex3–Hex10HexNAc2 |

**The 24 cleared values cannot serve any N-glycan prediction.** They are human milk
oligosaccharides on a lactose core, they carry `composition = None` (only an IUPAC string),
and the enumerator builds N-glycans on a complete branched Man3GlcNAc2 core. There is no
composition on which the two meet.

**The 89 held values are N-glycans and are held for a reason that reading a paper does not
fix.** Every one is blocked as *"not defined well enough to train: the analyte, the ion, the
reported uncertainty, or the record itself"*. In the seed rows the drift gas is `UNSTATED`
and `uncertainty_type` is `unknown` with no uncertainty value — and **hard constraint 7 refuses
a number without an uncertainty type.** The uncertainty is genuinely absent from the source,
so reading Hofmann 2014 for the drift gas would unblock the gas half and leave the
uncertainty half exactly where it is. **The gate is not loosened to reach a number.**

### Say the consequence plainly: `MEASURED_REFERENCE` is UNREACHABLE, not thin

Put the two halves above together and the result is stronger than either, so it is worth one
sentence with no hedging in it: **the measured-reference path is unreachable for every candidate
set this platform can produce in this release.** Not rare, not thin, not "sparse coverage" —
unreachable, without exception.

"24 cleared glycan records" reads like 24 usable values. The usable number is **zero**: 24 with no
composition to key them by, and 89 that have one but are held. So every response carries a stated
ABSENCE of CCS evidence, and the states a real deployment returns are `HELD_NOT_RELEASABLE`,
`NONE_IN_CORPUS_SEARCHED` and `LOOKUP_FAILED`.

**It is counted, not claimed.** The adapter's `reachable_states()` reports
`NOT REACHABLE ... 0 of 541 cleared record(s) carry a composition`, derived by counting the corpus,
and `test_the_unreachability_claim_fails_the_moment_a_cleared_record_carries_a_composition` fails
the moment that count changes. The same fact is stated at `CCSEvidenceState` in
`src/wmxglycan/ccs_evidence.py` and served in the API's `domain` block as
`holds_a_measured_cross_section_for_any_candidate: false`.

**And there is a second step behind the first.** A cleared record carrying a composition would
make the state reachable, and would still not separate candidates: a value keyed on composition and
ion is shared by every isomer of that composition. Discriminating between candidates needs a
measurement that resolves to one STRUCTURE, which is a different and much larger ask. The decision
rule `no cross section is held for this structure` reads the evidence LEVEL for exactly this
reason, and `EvidenceLevel.STRUCTURE` is what it is looking for.

The consequence for the ranker is in § 3 below: CCS is evidence where we hold it and a stated
absence where we do not, never a prediction.

## 3. What the ranker can and cannot rank on

### The curated rules cannot order a candidate set. This is arithmetic, not an omission.

**11 of the 15 curated biosynthetic rules reach N-glycan enumeration at all**: nine MGAT rules
governing branching order and bisecting interference, FUT8 on core fucosylation, and one
class-agnostic blood group rule. The other four are O-glycan core rules the enumerator never
applies. (This said "All 15 ... govern MGAT branching order and bisecting interference" until
27 September 2026. Wrong three ways, and in four files; `tests/test_glycan_constraints.py` pins the
split at 11 for N-linked and 5 for O-linked.) The conclusion is unaffected and is the point:
**not one of the eleven constrains galactosylation type, fucose position, chain extension or
LacdiNAc** —
which is precisely what the 167 candidates for `Hex5HexNAc4Fuc1` differ in. The enumerator
rejects every candidate that breaks a rule, so **every survivor satisfies every rule
identically**, and a score built from rule compliance is a constant across the set.

**AND "EVERY SURVIVOR SATISFIES EVERY RULE IDENTICALLY" IS ITSELF TOO STRONG, so it is split
into three named quantities.** What is constant is the number of rules VIOLATED, which is zero
— and the ranker verifies that by RUNNING the check rather than trusting that the enumerator
rejected the violators, because "zero" is also what an unrun check returns. "Rules satisfied" is
a different quantity and it measurably varies: the number of curated rules that actually bear on
a candidate is 1 to 5 across the 167 `Hex5HexNAc4Fuc1` candidates, and 0 for all six
`Hex5HexNAc2` candidates — so "15 rules satisfied" would have been false for every Man5
candidate.

A varying quantity published as a constant is worse than an omission, because the field's own
label tells the reader not to check it. So `rules_in_scheme`, `rules_violated` and
`rules_applicable` are three separate fields, the varying one is shown per candidate, and none
of them enters the number. The reason it does not: the count of applicable rules is largely a
count of how branched a candidate is, and no curated rule says a more branched structure is more
likely, so folding it in would be an uncited branching prior.

### Attestation discriminates, and it covers a tenth of the set

Measured here, against SugarBase v12 via glycowork 1.10.0: **12,664** reference records build,
**4,001** are fully resolved (8,663 carry an unresolved linkage or anomericity and cannot be
ground truth for a task about linkages), over **687** distinct compositions among the resolved
ones. Resolvedness is DERIVED FROM THE PARSE and never read off the record's own flags; a record
whose declared flags disagree with its own structure string is refused by name, because
`GlycanStructure` validates only that some structure identifier exists and never compares the two.

| composition | candidates | reference rows | distinct reference structures | candidates attested | reference structures NOT enumerated |
|---|---|---|---|---|---|
| `Hex3HexNAc4Fuc1` | 10 | 26 | 22 | **6** | 16 |
| `Hex4HexNAc4Fuc1` | 63 | 27 | 26 | **16** | 10 |
| `Hex5HexNAc4Fuc1` | 167 | 38 | 34 | **17** | 17 |
| `Hex5HexNAc2` | 6 | 30 | 28 | **6** | 22 |

Pooled to the indistinguishable class — which is the unit the ranker bands, because members of a
class cannot be separated by anything the platform computes — the same corpus gives:

| composition | candidates | classes | largest class | bands | pooled attestation |
|---|---|---|---|---|---|
| `Hex3HexNAc4Fuc1` | 10 | 9 | 2 | 2 | 6 classes at 1, 3 at 0 |
| `Hex4HexNAc4Fuc1` | 63 | 32 | 4 | 3 | 2 at 2, 12 at 1, 18 at 0 |
| `Hex5HexNAc4Fuc1` | 167 | 61 | 10 | 3 | **1 at 4**, 13 at 1, 47 at 0 |
| `Hex5HexNAc2` | 6 | 6 | 1 | **1** | 6 classes at 1 — REFUSED |

Two limits fall straight out of the first table and both are load-bearing:

1. **150 of the 167 G2F candidates have no attestation at all.** Absence in a 12,664-record
   literature database is not evidence of absence, so a score that assigns them zero asserts
   something the data does not support. See § 4 on the weight they get and why it is a declared
   policy rather than a derived fact.
2. **The enumerator misses between a third and four fifths of the attested structures of a
   composition, and it is DIAGNOSED rather than open.** For `Hex5HexNAc4Fuc1`, 17 of the 34
   fully-resolved reference structures of that composition are not among its 167 candidates; for
   `Hex5HexNAc2` it is 22 of 28. **This was new in this repository and was not measured in
   `Project2`.**

   Diagnosed on the smallest case, `Hex3HexNAc4Fuc1` — 10 candidates against 22 resolved
   reference structures, 16 not enumerated — by comparing every bond in the missing set against
   every bond the candidates carry. **It is not a key-matching bug**, which was the first thing
   to rule out: 6 of the 10 candidates ARE matched, so the canonicalisation works on both sides.
   The 16 decompose into three causes:

   | cause | example | verdict |
   |---|---|---|
   | **the deliberate mammalian scope** | `Fuc(a1-3)` on the reducing GlcNAc, 3 structures. Insect and plant N-glycans carry a core alpha1,3-fucose; the enumerator places core fucose only at a1-6. | CORRECT. Out of scope by design, and `n_glycan_implausibility_reasons` already records the mammalian assumption. |
   | **structures the curated rules deliberately reject** | an antenna on the a1-6 arm with none on the a1-3 arm, which MGAT2's REQUIRES rule forbids. | CORRECT, and already recorded: these are among the 280 observed structures that contradict the curated mammalian rules, which are deliberately not loosened. |
   | **an elaboration-vocabulary gap** | five monolinks no `SITES` entry places: `GlcNAc(a1-2)Man`, `GlcNAc(b1-6)Man` on the a1-3 arm, `Man(b1-3)Man`, `Man(b1-6)Man`, and `Fuc(a1-4)GlcNAc` (Lewis a). | A REAL GAP. Any narrowing or widening of the vocabulary must come from a rule with a citation, never from choosing a vocabulary that lands on a target. |

   **"The true structure is among the candidates" remains an unsafe assumption**, and the ranker
   no longer makes it: every share is computed over a hypothesis space that includes the
   reference structures no candidate matches, and `share_not_enumerated` reports the mass sitting
   on them — **30.1% for `Hex5HexNAc4Fuc1` and 77.2% for `Hex5HexNAc2`**, with a further 0.9% and 1.8% on the catch-all hypothesis described in § 4.

### A count of 2 was not two reports. It was one structure spelled twice.

**This section said the wrong thing on 26 September 2026, and the correction is the most
valuable thing an adversarial pass produced.** It said "a count of 2 against a count of 1 is two
papers against one paper". That is not what a count of 2 was.

Measured here: the 4,001 fully-resolved reference records index to **3,640 distinct canonical
keys**. 352 keys carry more than one row — 343 twice, 9 three times, **361 duplicate rows in
total** — and not one of those collisions is two identical IUPAC strings. Every one is the same
structure written with its branches in a different order, and several carry the **same GlyTouCan
accession on both rows**, which is the source's own statement that they are one structure and
not two reports.

The consequence for a ranker was total. Counting rows put the top of every ranking in the hands
of how SugarBase happens to spell things: for `Hex5HexNAc4Fuc1` the three highest-scoring
candidates were the three the database spells twice, and for `Hex5HexNAc2` — Man5, the most
studied N-glycan there is — one candidate would have "won" on the strength of a single structure
entered once with accession G83351GR and once with none.

**The asymmetry was visible inside one package and nobody had looked.** `enumeration.py` already
deduplicates the CANDIDATE side on `graph.canonical_key()`; the first version of the attestation
index applied no such deduplication to the REFERENCE side, using the very same function. So the
unit of evidence is now the distinct structure: `structures_for()` returns 0 or 1 and is the only
thing a score may use, and `rows_for()` exists so the discrepancy stays visible rather than being
corrected away silently.

What survives of the original claim, and it is still load-bearing: SugarBase carries no
abundance, no tissue and for most rows no species, so **a count here is a count of deposited
structures and supports no statement about how common one is in a sample.**

### A tie is not one thing, and conflating the two would hide the worse one

- **Indistinguishable by features.** Two candidates whose 44-column feature rows are
  byte-identical. Established in `Project2`: two isomers differing only in a linkage position
  on an otherwise identical shape match on all 44 columns. The platform cannot order these
  even in principle.
- **Indistinguishable by evidence.** Two candidates with different feature rows and the same
  score — which, with 150 of 167 unattested, is most of the set.

Both must be returned as unordered groups. An arbitrary order rendered as a ranking is the
failure this project has been built against throughout, and the second kind is the common one.

## 4. The confidence number: what it is normalised over, why the prior is policy, and why it
has never been calibrated

The number the ranker returns is **a share of evidence weight, not a probability of being
correct**, and it is graded as such. The distinction is the same one
[`LIMITATIONS.md` § 7C](../LIMITATIONS.md) draws for the CCS confidence grades: "the
confidence grades are rules and have never been calibrated".

**It is not normalised over the candidate set.** Normalising over the candidates alone asserts
the answer is among them, and § 3 measures that false between a third and four fifths of the
time. The hypothesis space is therefore every candidate class PLUS every fully-resolved
reference structure of the composition that no candidate matches, and the share sitting on that
second group is reported. For `Hex5HexNAc2` it is **77.2%**, which is the honest headline for
that composition: over three quarters of the evidence is on structures the platform did not
propose.

**AND THAT WAS STILL NOT ENOUGH.** Reserving mass only for reference structures the corpus knows
about closed the world whenever it knew of none that were missed. `Hex6HexNAc3Fuc3` has exactly
ONE fully-resolved reference structure, that structure IS among its candidates, so
`not_enumerated` was 0 and the class shares summed to exactly **1.0** — the response asserting
the answer was in the set on the strength of one deposition, and the published decision rule
"the candidate set may not contain the answer" reporting `fires: False`.

So the space carries a third kind of hypothesis, always: **a structure neither proposed by the
enumerator nor deposited in the corpus.** It has no observations, so its weight is the prior
pseudocount alone, and its share is published as `share_not_proposed`. The corpus can show a
candidate set is INCOMPLETE; it can never show one is COMPLETE, which is why `SetCompleteness`
has no `COMPLETE` member and `truth_may_not_be_in_the_candidate_set` is a constant `True` with
its reason written into the property.

**The event the number describes is a literature deposition, not a molecule in a sample.** That
sentence ships in the response as `share_means`, not only here, because the consumer reads the
response.

**RULING, 27 September 2026: the number is never labelled a probability, and that is settled.**
The owner's reading of the spec is that it asks for a probability-LIKE confidence rather than a
probability, so the served field is `evidence_share` and no field anywhere in the API is named
`probability`, `p_correct` or `likelihood`. A test walks the response models and fails on any
such name. The declared POLICY prior and all four standard alternatives travel in every
response that carries the number, so a consumer can see the lever without reading this file.

**THE PRIOR IS A CHOSEN CONSTANT AND IS PUBLISHED AS ONE.** An earlier draft of this section
claimed there was "no free parameter" and that the uniform Dirichlet prior was *forced* by the
enumerator having already applied every rule it has. **That was wrong, and it was the most
seductive error in the design**, because it dressed a choice as an absence of choice.

Rule-equivalence forces the prior to be SYMMETRIC. It says nothing about its CONCENTRATION, and
uniformity on the simplex is a choice of coordinates rather than an absence of information:
Jeffreys (1/2), Perks (1/K) and Haldane (0) are each standard and each has a published
justification. The choice is also the single largest lever in the output: for
`Hex5HexNAc4Fuc1` the top class's share ranges from about 4.5% to about 11.5% across the four,
and the mass on structures nobody enumerated ranges from 30% to 50%.

So `PRIOR_PSEUDOCOUNT` is declared POLICY, chosen rather than measured, in the same idiom this
repository already uses for `MIN_GROUP_RECORDS` and `LEVERAGE_LIMIT_PERCENT` — and **every
response publishes the share under all four priors.** A field that ships four numbers for one
quantity cannot be read as "the probability".

**What would calibrate it:** a set of compositions whose true structure is independently known,
scored blind, and the realised frequency of the true structure appearing in each score band
compared against the band. Nothing in this repository can do that today, and the ranker says so
on every response rather than in this file only.

**The ordering constraints are reported and not scored.** Four of the nine FORBIDS rules
describe assembly order rather than coexistence, so they travel with a candidate as an
`ordering` note instead of excluding it. Turning "this needs a particular assembly order" into
a number requires a weight, and no curated source provides one. Reported with enzyme,
rationale and reference; absent from the score.

## 4A. Three defects in the ported code, found by attacking the ranker's design

None of these were in the new code. All three were found by adversarially reviewing a DESIGN
that was going to consume them, which is worth noting on its own: the ported suite is green and
was mutation-verified where it was written, and **none of these failed a test.**

### FIXED: 108 of 167 candidates published a cited claim that was false about them

`Enumerator._candidate` attached an ordering caveat whenever a curated rule's PRODUCT was
present, and the caveat's own words are *"This candidate carries {product} alongside that
context"*. Measured on `Hex5HexNAc4Fuc1` before the fix: **114 of 167 candidates carried an
ordering note and only 18 contained the bisecting GlcNAc that all four order rules name as their
context.** So 108 candidates published a sentence that was false about them, carrying an enzyme
and a *Schachter 1986* citation — the form a reader trusts most.

Fixed by gating the caveat on the context being present. After the fix: **6 notes on
`Hex5HexNAc4Fuc1`, all FUT8, all on candidates that really are bisected and core-fucosylated,
and zero false claims.** Candidate counts are unmoved at 10 / 63 / 167 / 6, because the caveat
never excluded anything — only the rationale changed.

### NOT FIXED, and recorded: `_context_holds` treats two different GlcNAc as one residue

The first attempt at the fix above used `_context_holds`, the same predicate `_broken_rules`
uses. That suppressed FUT8's caveat entirely — including on the bisected, core-fucosylated
structures it exists to explain, of which there are 321 in the reference set. The cause:
`_context_holds` anchors two fragments to one residue when they end at the same
(residue, linkage) pair, and FUT8's product `Fuc(a1-6)GlcNAc` and the bisecting context
`GlcNAc(b1-4)Man(b1-4)GlcNAc` both end at `('GlcNAc', None)` — **the reducing GlcNAc and the
chitobiose GlcNAc, which are different residues.**

So an order rule's context is searched anywhere in the molecule instead, which is also the right
reading of what an order rule says: *"bisecting GlcNAc blocks subsequent MGAT2 action"* is a
claim about the state of the molecule, not about two residues sharing a position. **The anchoring
defect itself is left alone**, because `_context_holds` also gates `_broken_rules` for the
non-order rules, where changing it could change which candidates are refused. That is a
rule-semantics decision with a candidate-count consequence, and it is not made as a side effect
of building a ranker.

### NOT FIXED, worked around: `FeatureVector.as_row()` equality holds by a CPython accident

`as_row()` maps an absent value to the module-level `math.nan` singleton.
`(1.0, math.nan) == (1.0, math.nan)` is `True` — but **only because tuple comparison
short-circuits on element identity.** `(1.0, float("nan")) == (1.0, float("nan"))` is `False`,
and so is a numpy NaN or a JSON round trip. Verified in this environment and pinned by a test.

Keyed on `as_row()`, the ranker's indistinguishable classes would silently become singletons the
moment anyone wrote `float("nan")`, routed a row through numpy for an estimator — which
`splits.py` already feeds toward — or serialised a row. The platform would then begin ordering
152 of 167 candidates it cannot distinguish, **with no test failing.** So `ranking.class_key_for`
is built from the None-bearing mapping, where `None == None` is identity-true, and it raises if a
NaN reaches it.

**This is a latent defect in the ported featuriser, not only in the ranker's use of it.** Any
other consumer of `as_row()` that compares or hashes rows has the same hole, and `splits.py`'s
leakage guard is the first place to look.

## 4B. What an adversarial pass over the FINISHED ranker found

The design was reviewed before the code was written, which caught three fatal flaws in the
design. The code was then reviewed again before it was committed, and that caught **eleven more
in the implementation** — which is the more useful fact of the two: a design that has survived
one adversarial pass reads as finished, and this one was not.

**Fatal, and it was the design's own central promise.** For some real compositions the shares
summed to exactly 1.0 over the candidate set — the closed-world claim the module docstring said
three separate things prevented. Fixed with the catch-all hypothesis above.

**The rule check was measured and then inert.** `rules_violated` was computed by running the
check, published, and acted on by nothing: a candidate set breaking curated rules was banded and
returned with no refusal and no weakness. A non-zero count means the enumerator and the rule
check disagree, so the ranker now refuses.

**A refused composition fabricated an absence.** `_refused()` never consulted the index it was
handed, so a composition the enumerator declines was reported as one the corpus knows nothing
about. It holds 15 fully-resolved structures of `Hex2HexNAc2Fuc2`, and the claim was attached to
a refusal — where a reader takes it as the reason nothing could be ranked. It also published
`rules_in_scheme=0` with the verified flag left at `True`.

**A self-asserted flag was trusted where the level should have been read.** `_decide` consulted
`discriminates_between_candidates` without checking `EvidenceLevel.STRUCTURE`, and the validator
allowed a COMPOSITION-level finding to set it — so a measurement shared by all 167 candidates
could report that the rule "no cross section is held for this structure" did not fire. That is
the same shape as the `is_measured_reference` flag this module had already removed for being
self-asserted. Both ends are guarded now.

**The licence gate required its field to exist without requiring it to permit anything.**
`CCSReference` accepted every `ReuseStatus`, including the `unverified` DEFAULT and `excluded`.
Three statuses have no permitted use at all and are now refused, derived from the predicates
rather than listed.

**An absence accepted a denominator of zero.** "We searched a corpus of 0 records and found
nothing" was a valid finding, and it is what an unsearched corpus looks like too. `HeldValues.records`
next door was already pinned `gt=0`; that was the asymmetry.

**A wiring bug was delivered as the honest default.** A supplied lookup with no adduct or charge
returned `NOT_CONSULTED` — "nothing looked" — indistinguishable from a deployment with no
adapter. It raises now.

**Four tests could not fail**, and three of them are kept with the tautology quoted in the
docstring because the shape is more useful than the fix:

| test | why it could not fail |
|---|---|
| every candidate in a class shares one band | `ranks = {band_of[...]}` is a one-element set literal |
| candidates in a class are different structures | `len(frozenset) == len(set(frozenset))` holds for any code |
| ordering caveats are kept out of the number | the stated property was a comment; the dict was built and discarded |
| AI_ONLY is unreachable across the input space | two rules fire on every input, so the sweep never reached the branch it asserted about — it would have passed with the branch deleted |

The last one is the instructive one. **An input sweep cannot prove a branch is live**, and the
fix is not a better sweep: the sweep now asserts what it really establishes (no constructible
input yields AI_ONLY, and the two universal rules are why) and a separate positive control
patches both gates and reaches the branch. Neither alone is a proof.

**And one test file committed the failure class in the act of citing it.**
`test_every_module_the_ranker_introduced_has_at_least_one_mutation` hand-wrote the set
`{"ranking.py", "attestation.py", "ccs_evidence.py"}` in the same breath as a comment saying
"a hand-written list of what is covered is instance Twelve". It is now derived from the port
record in `GLYCAN_PORT.md` — a module is written here exactly when it is absent from the digest
table — with a floor under the derivation.

**A crash nothing had reached.** `summary()` raised `TypeError` on `Hex5HexNAc2` because the
unattested-class mass was `None` where the honest value is a measured `0.0`. No test had
rendered a refused set. A method nothing calls is a method nothing protects.

### The decision field is CONSTANT at IM_VALIDATION_REQUIRED. This is not a broken field.

**RULING, 27 September 2026: it stays constant, and this section exists so that nobody later
reads a single-valued field as a defect and "fixes" it.** Everything landing on
`IM_VALIDATION_REQUIRED` is the true state of a platform with no CCS model, and it is exactly
what the spec means by identifying predictions that require experimental validation. A field
that said otherwise for some inputs would be the defect.

The brief asked for three values and expected "almost everything" to land on REQUIRED. The
measured answer is stronger: **everything does**, for two independent reasons that hold for
every possible input.

| value | reachable? | what is holding it shut |
|---|---|---|
| `AI_ONLY` | no | `no validated model exists` fires on every input — nothing here has been validated against known truth. AND `the candidate set may not contain the answer` fires on every input, because completeness is not establishable from a corpus of depositions. **Two independent gates**, and `tests/test_glycan_ranking.py` opens both to prove the branch is live code rather than an enum member nothing references. |
| `IM_VALIDATION_RECOMMENDED` | no | the same two rules. Any firing rule yields REQUIRED. |
| `IM_VALIDATION_REQUIRED` | yes, always | |

**IF YOU ARE READING THIS BECAUSE THE FIELD LOOKS BROKEN, IT IS NOT.** What would move it:

1. a model validated against compositions whose true structure is independently known; and
2. a cross section that resolves to ONE structure rather than to a composition and an ion —
   a composition-level measurement is shared by every isomer and cannot confirm one of them.

This repository has neither, and the fields that would carry them exist and are refused rather
than absent. **Nothing was tuned to make some cases look better**, which was an explicit
instruction; the rules are published per response with `applies_to` on each, so a caller can
read why each one fires instead of trusting this table.

### Deferred, with the reason

- **8 of Man5's 22 un-enumerated reference structures have feature rows identical to a class the
  ranker returned.** By the module's own pooling rule an observation of an indistinguishable
  structure is an observation of the class, so pooling those would give Man5 two bands instead of
  one and it would no longer refuse. That is a real refinement and it changes a headline result,
  so it needs its own review rather than being taken while tidying up.
- **`dataclasses.asdict()` drops every honesty-bearing field**, because they are properties, and
  `members` is a `frozenset` that `json.dumps` refuses. A serialisation contract belongs with the
  API layer that will consume this, not ahead of it.

### The sweep found two guards that no test protected, after two adversarial passes

The first glycan sweep ran 35 mutations and killed 33. **Both survivors were guards this
repository had just written, in code two adversarial passes had already been over.** That is the
argument for the harness in one sentence: a green suite and two reviews still left two behaviours
with nothing behind them.

**One: a validator made a redundant guard untestable.** `_decide` reads
`evidence.level is EvidenceLevel.STRUCTURE` as well as the `discriminates_between_candidates`
flag — defence in depth, added by the implementation review. But the `CCSEvidence` validator now
*refuses* that combination outright, so no ordinary construction can tell the two versions of
`_decide` apart, and the mutation removing the level check survived.

The two fixes were each correct and together they hid one of themselves. The test now builds the
illegal object with `model_construct`, which skips every validator, exactly as the CCS core's
suite does for a record that reached it without validation — and that is the real threat model
rather than a contrivance: the validator guards the constructor, and a producer using
`model_construct` or `model_copy(update=...)` reaches `_decide` without ever passing it.

**Two: a clause added and never exercised.** The `records_consulted <= 0` arm was added to close
the zero-denominator hole, and the test beside it covered `None` and a missing corpus name and
not zero. The guard was present and unprotected for as long as it existed.

Both are killed now, and both are the same shape as instance Twelve: **a guard whose coverage
nobody checked.** Neither was found by writing more tests. They were found by breaking the code
on purpose and noticing that nothing complained.

### A process failure worth recording: a review agent edited the live source

The implementation review ran against the working tree while it was being edited, and one
verifier reproduced a finding by writing `RuleAccounting(rules_in_scheme=999, rules_violated=7)`
into `ranking.py` and did not revert it. It was caught within ten minutes, by an anchor mismatch
rather than by a failing test, and only because the file was **untracked** — so `git diff` showed
nothing and the edit was invisible to the obvious check. Every other new file was audited by
modification time and content and was intact.

**The lesson is about sequencing, not about the agents.** A reviewer that can write to the tree
it is reviewing needs the tree committed first, so that any edit it makes shows up in a diff.
Reviewing uncommitted work with write-capable agents means the only record of what changed is
the file itself.

## 4C. The six endpoints, and what each of them cannot do

Built 27 September 2026. Real services over a SQLite store that survives a restart, not a
demonstration UI.

| endpoint | what it does | what it cannot do |
|---|---|---|
| `POST /v1/predictions` | ranks the candidates for a composition and ion, freezes the answer, returns 201 | it cannot predict a cross section, and it does not try |
| `GET /v1/predictions/{id}` | the frozen prediction, served from the stored bytes | it never recomputes, so a corpus that has moved cannot alter a past answer |
| `POST /v1/predictions/{id}/validation` | appends experimental measurements, returns 201 | it cannot touch the prediction; the response carries its digest read before and after |
| `GET /v1/predictions/{id}/comparison` | delta CCS, interval coverage, status | no delta against a PREDICTION exists; see below |
| `GET /v1/runs` | history, filterable by kind, composition, prediction and time | |
| `GET /v1/models/current` | pipeline fingerprint, data snapshot, applicability domain | it reports no fitted CCS model, because there is none |

### Immutability is enforced by the schema, not by the caller

`store.py` contains **no UPDATE statement and no DELETE**, and a test asserts that by reading
the module. A second write to a frozen prediction raises `AlreadyFrozen` and writes nothing; a
second create under the same `client_reference` answers **409 and names the prediction that
stands**, so a caller fetches it rather than guessing. Attaching is append-only, and the attach
response carries `prediction_digest_before` and `prediction_digest_after` read from storage, so
"the prediction was not touched" is a measurement in the reply rather than a promise in a
docstring. All of it still holds after a restart, which is the only version of immutability
worth having.

### The comparison endpoint: three axes, and only one of them works today

Collapsing them into a single status would have lost the one that does work.

| axis | state today | why |
|---|---|---|
| against the **prediction** | `no_predicted_value`, always | there is no predicted cross section, so a delta is not unmeasured - it is undefined. An interval coverage figure would describe an interval that was never produced. |
| against a **reference** the platform holds | `no_reference_in_corpus` for the fucosylated compositions; `reference_held_not_releasable` for `Hex5HexNAc2`, where 8 records exist and are blocked | the delta arithmetic is live and is exercised by a synthetic releasable reference, because every path against the real corpus is an unevaluable one and an unexercised arithmetic path is an unguarded one |
| among the **caller's own** measurements | **computed, always** | it needs no model and no reference. It is the caller's own reproducibility, grouped by adduct, charge, platform AND drift gas - hard constraint 4, so values against different gases are never spread against each other |

`agreement_limit_percent` is **2.0 and declared POLICY**, chosen rather than measured, in the
same idiom as the CCS core's limits-of-agreement threshold. It is published in every comparison
response so the limit a verdict rests on travels with the verdict.

### The pipeline fingerprint, and the one thing it deliberately does not catch

There is no fitted model, so `GET /v1/models/current` reports the fingerprint of a
**deterministic pipeline**: `ac288c6b0f83/95e925999909` at the time of writing, over the curated
rules, the enzyme table, the placement vocabulary, the 44 feature columns, the prior policy, and
a digest of the 3,640 canonical keys in the corpus.

**It does not digest module source, and the cost of that is stated rather than hidden: a logic
change that touches no table, no corpus and no policy constant will not move it.** A fingerprint
that moved on a comment would train its readers to ignore it, and the CCS core's has held across
five releases precisely because it digests inputs rather than code. The mutation catalogue is
what guards logic; this digest guards the inputs.

## 4D. The dashboard, and the two places it refuses to look finished

Wired 27 September 2026 from the static prototype: six pages, the prototype's layout and design
tokens kept, and every number on it now from a response. It is served by the service itself at
`/`, so it is same-origin - no CORS to configure, no build step, no second process that could
drift out of step with the API it is a client of. One self-contained file, as the prototype was.

### The candidate list: bands are ordered, nothing inside one is

This is where false certainty would have entered. The prototype rendered `#1 Candidate
GLY-ISO-001` with a score of `0.54`. For `Hex5HexNAc4Fuc1` that shape would number 167
candidates, of which **152 sit in a class the platform's own 24 structure columns cannot
separate** - an arbitrary order presented as a ranking.

So no candidate is numbered anywhere. Rendered from a real response and checked by executing the
page's own JavaScript:

| what the page emits | count |
|---|---|
| positions given to anything | **3** - "Band 1 of 3", "Band 2 of 3", "Band 3 of 3" |
| tied groups shown, each with its size | **46**, sizes 10 / 8 / 6 / 5 / 4 / 3 / 2 |
| candidates inside those groups | **152**, listed unordered |
| "alone in its class" badges | **15**, the singleton classes |
| candidates numbered | **0** |

Each tied group says, in the page: *"N candidates sharing this position"*, *"tied · not
separable"*, and *"The platform cannot tell you which of them it is, and this dashboard will not
guess."* A refused set - Man5 - is headed **THIS SET IS NOT ORDERED** and still lists its
candidates: what is withheld is an order, not an answer.

### The comparison: an unevaluable axis shows its reason, never a blank

A dash in a delta column reads as "about zero" to anyone skimming, so there are none. All three
axes render with their own state and the service's own reason text:

- **Axis 1, against the AI prediction** - `undefined`, always, with the full reason. Not empty.
- **Axis 2, against a reference** - a real delta where one is releasable; otherwise
  *"No delta, and not because it is zero:"* followed by the evidence state's own sentence.
- **Axis 3, the caller's own measurements against each other** - computed, always.

Interval coverage is never rendered as YES or NO. The prototype showed `Interval coverage: YES`;
there is no interval, so that was false. It now shows `not evaluable` and why.

### What the page does NOT pretend to send

The prototype's prediction form offered drift gas, instrument, antibody class, glycosylation
site, retention time and MS/MS evidence. The service takes composition, adduct and charge and
nothing else. A field that looks submitted and is not is the same class of false certainty as a
fake rank, so they are kept in the layout under the heading **"Not sent, and not silently
dropped"**, with the note that drift gas, instrument and calibrant *are* real on the validation
page, where they describe a measurement.

### Two guards on the client that a browser-less test can still hold

- **Every field the JavaScript reads off a response is a real field of a real response model**,
  derived from the pydantic models on one side and from the shipped script on the other, with
  snake_case as the discriminator (the API's fields are snake_case; the browser's own APIs are
  camelCase). 79 field reads, 0 unknown. A mistyped field renders `undefined` in a browser and
  fails nothing, which is the worst available silence. A floor test proves the check would
  notice a typo.
- **The page's own rendering functions are executed over real responses** via node, cut out of
  `dashboard.html` rather than re-implemented, so what is checked is what ships. That is how the
  figures in the table above were obtained.

## 5. Mutation coverage: both packages are shadowed, and the ported modules are deliberately bare

`tools/mutation` was scoped to `src/wmxccs`. It now shadows **both** packages on every run —
on every run, even a filtered one, so a filtered run and a full sweep measure the same code —
and a mutation carries the package its file sits in. That field is load-bearing rather than
tidy: both packages have a `models.py`, a `sources.py`, a `reuse.py`, a `licensing.py` and a
`contracts.py`, so a mutation identified by filename alone would name two files and the sweep
would pick one of them by accident. That is this tool's own failure mode, so it is not left to a
convention.

The catalogue is in two modules and the concatenation is DERIVED: `catalogue.py` holds the 282
wmxccs entries and is closed, `catalogue_glycan.py` holds 22 glycan entries, and
`runner.all_mutations()` joins them. A split is a hand-maintained duality, which is instance
Twelve's shape, so `tests/test_glycan_mutation_catalogue.py` asserts the derivation holds every
entry of both and that no label collides across them. It is a FUNCTION and not a module
constant, because a constant computed at import time silently ignores a test that patches
`MUTATIONS` — four existing tests do exactly that, and they caught it the first time it was
written as a constant.

**The ported modules are deliberately NOT covered**, on the owner's instruction: they were
mutation-verified in `Project2` and re-covering them costs days the deadline does not have.
There is one exception and it is asserted rather than left to drift — `enumeration.py` carries
exactly one entry, for the ordering-note context gate, because that is the one line this
repository changed.

The honest reading of the ported 1,220 tests is that they are green and were mutation-verified
elsewhere against a copy that was byte-identical **at the moment the port landed**. Fourteen of
the fifteen modules still are; `enumeration.py` is not, because the ordering-note context gate was
fixed here. That one region carries a catalogue entry of its own; the rest of `enumeration.py`
does not. Good evidence, and not the same as evidence gathered here.
