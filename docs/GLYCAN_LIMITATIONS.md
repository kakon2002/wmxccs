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

The consequence for the ranker is stated in § 3 below: CCS is evidence where we hold it and a
stated absence where we do not, never a prediction.

## 3. What the ranker can and cannot rank on

### The curated rules cannot order a candidate set. This is arithmetic, not an omission.

All 15 curated biosynthetic rules govern MGAT branching order and bisecting interference.
**Not one constrains galactosylation type, fucose position, chain extension or LacdiNAc** —
which is precisely what the 167 candidates for `Hex5HexNAc4Fuc1` differ in. The enumerator
rejects every candidate that breaks a rule, so **every survivor satisfies every rule
identically**, and a score built from rule compliance is a constant across the set.

The ranker therefore reports the rule component as the constant it is, with the count, rather
than folding it into a number where it would look like a contribution. A component that
cannot discriminate and is presented as if it might is how a score stops meaning anything.

### Attestation discriminates, and it covers a tenth of the set

Measured here, against SugarBase v12 via glycowork 1.10.0: **12,664** reference records
build, **4,001** are fully resolved (8,663 carry an unresolved linkage or anomericity and
cannot be ground truth for a task about linkages), over **1,837** distinct compositions.

| composition | candidates | resolved reference structures of that composition | candidates attested | reference structures NOT enumerated |
|---|---|---|---|---|
| `Hex3HexNAc4Fuc1` | 10 | 22 | **6** | 16 |
| `Hex4HexNAc4Fuc1` | 63 | 26 | **16** | 10 |
| `Hex5HexNAc4Fuc1` | 167 | 34 | **17** | 17 |
| `Hex5HexNAc2` | 6 | 28 | **6** | 22 |

Two limits fall straight out of that table and both are load-bearing:

1. **150 of the 167 G2F candidates have no attestation at all.** Absence in a 12,664-record
   literature database is not evidence of absence, so a score that assigns them zero asserts
   something the data does not support. See § 4 on why the weight they get is derived rather
   than chosen.
2. **The enumerator misses half the attested structures it should cover.** For
   `Hex5HexNAc4Fuc1`, 17 of the 34 fully-resolved reference structures of that composition are
   not among its 167 candidates. For `Hex5HexNAc2` it is 22 of 28. **This is new in this
   repository and was not measured in `Project2`.** Some of it is the known complete-core
   limit (421 of 4,001 resolved reference structures have no complete branched core, carried
   over, not re-measured here); the rest is elaboration vocabulary the enumerator does not
   place. **Until it is diagnosed, no coverage claim should be made for the candidate set, and
   "the true structure is among the candidates" is not a safe assumption.**

### Attestation counts are literature reports, not abundances

SugarBase records that a structure has been reported. It carries no abundance, no tissue and,
for most rows, no species. **A count of 2 against a count of 1 is two papers against one
paper, not twice as common.** This is the single largest reason the ranker's number is not a
probability of occurrence, and the reason is in the data rather than in the implementation.

### A tie is not one thing, and conflating the two would hide the worse one

- **Indistinguishable by features.** Two candidates whose 44-column feature rows are
  byte-identical. Established in `Project2`: two isomers differing only in a linkage position
  on an otherwise identical shape match on all 44 columns. The platform cannot order these
  even in principle.
- **Indistinguishable by evidence.** Two candidates with different feature rows and the same
  score — which, with 150 of 167 unattested, is most of the set.

Both must be returned as unordered groups. An arbitrary order rendered as a ranking is the
failure this project has been built against throughout, and the second kind is the common one.

## 4. The confidence number is normalised, derived, and has never been calibrated

The number the ranker returns is **a share of evidence weight, not a probability of being
correct**, and it is graded as such. The distinction is the same one
[`LIMITATIONS.md` § 7C](../LIMITATIONS.md) draws for the CCS confidence grades: "the
confidence grades are rules and have never been calibrated".

**There is no free parameter, and that is a deliberate constraint rather than an achievement.**
The weight of an unattested candidate is one observation-equivalent, which is the posterior
mean of a multinomial under a uniform Dirichlet prior — the prior that says no candidate is
favoured before the corpus is consulted, which is exactly the position the enumerator leaves us
in, since it has already applied every rule it has and the survivors are rule-equivalent. Any
other value for that weight would be a number chosen to move the answer.

**What would calibrate it:** a set of compositions whose true structure is independently known,
scored blind, and the realised frequency of the true structure appearing in each score band
compared against the band. Nothing in this repository can do that today, and the ranker says so
on every response rather than in this file only.

**The ordering constraints are reported and not scored.** Four of the nine FORBIDS rules
describe assembly order rather than coexistence, so they travel with a candidate as an
`ordering` note instead of excluding it. Turning "this needs a particular assembly order" into
a number requires a weight, and no curated source provides one. Reported with enzyme,
rationale and reference; absent from the score.

## 5. The glycan layer has no mutation coverage for the ported code

`tools/mutation` and its catalogue were scoped to `src/wmxccs`. The port did not extend it to
the 15 ported modules, on the owner's instruction: they were mutation-verified in `Project2`
and re-covering them costs days the deadline does not have. **New code written here is
covered.** The honest reading of the ported 1,220 tests is therefore that they are green and
were mutation-verified elsewhere, against a copy that is byte-identical — which is good
evidence and is not the same as evidence gathered here.
