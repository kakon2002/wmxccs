# Limitations

## READ THIS FIRST: what these numbers are, in plain terms

This platform takes a molecule's "collision cross section" - a measured size, in square
angstroms - as reported by four different kinds of instrument, and reports how much those
instruments disagree. It then offers a corrected value that puts one instrument's number
onto another's scale, **alongside the original, never instead of it**.

**Where the numbers come from.** One published study of 87 steroids
(DOI 10.1021/jasms.2c00196). That is the whole corpus. 142 ions are measured on more than
one kind of instrument; 93 of them on all four.

**What the corrections say.** Between one instrument and another, in that one study, the
values differ by roughly 0.1 to 1.1 per cent depending on the instrument pair and the ion.
The corrections move a value by about that much.

### Seven things it would be easy to read into these numbers that are not there

**1. "Drift tube versus trapped ion" is one laboratory comparing its own two instruments.**
Those measurements were made by the study's authors on their own equipment. If a different
laboratory ran the same ions, it might not get the same difference.

**2. "Travelling wave" is not that laboratory's measurement at all.** Those values are
taken from an earlier paper (DOI 10.1021/acs.analchem.9b05247) and are already an average
over four instruments in several laboratories. So a travelling-wave comparison is one
laboratory's number against somebody else's average - which is a different thing again from
the first case, and neither is wrong, but they are not the same kind of comparison.

**3. NEITHER IS A REPRODUCIBILITY FIGURE.** "How much does this instrument type vary
between laboratories" is a different question, answered by measuring the same ions
independently in many laboratories and looking at the spread. That has been done elsewhere -
the published figure for one drift-tube method is about 0.29 per cent - and it is not what
this platform produces. The software will not let a number from here be labelled that way:
the claim cannot be expressed in the code at all.

**4. The uncertainty ranges are guaranteed at 80 per cent, not 90.** Each corrected value
comes with a range. The range is calculated at a 90 per cent setting, but the method used
(jackknife+) mathematically guarantees only 80 per cent - that is, about one in five values
may fall outside its stated range rather than one in ten. Both numbers are reported with
every range, and the smaller one is the one to rely on.

**5. One of the internal safety checks currently does nothing.** A check exists to warn when
a correction is being driven by a handful of unusual ions rather than by the data as a whole.
On this corpus it never triggers, and on eight of the eighteen comparisons it cannot even be
calculated. Its threshold has deliberately NOT been lowered to make it trigger, because
adjusting a safety check until it fires is fitting the check to the data.

**6. A second safety check cannot run at all for a molecule we have not already measured.**
One of the checks asks "is this particular ion a known bad case - one where the correction
is already known not to work?" It answers that by looking the ion up in our own data. A
molecule that is not in our data cannot be looked up, so for a new molecule the check does
not say "fine", it says nothing. **A new molecule is the ordinary case for anyone using
this service.** Until 20 September 2026 the response reported that silence as though the
check had passed. It now says which checks could not be run, and lowers the confidence
grade by one step when any of them could not. That step is not a claim the value is wrong -
it means less was verified than the grading scheme describes.

If you are bringing a molecule this platform has not measured, **expect that check to be
unavailable every single time, not now and then.** Counting it across our own stored data
makes it look like an occasional gap, roughly one time in three; that count is misleading
because our own stored data is, by definition, data we already have. For anything new, the
check simply does not apply.

**7. The confidence grades are capped by how much data we have, not by how good your
measurement is.** Every corrected value today comes back `qualified` or `weak`, and none
comes back `supported`. That reads like a verdict on the measurement and it is not one. The
grading scheme lowers its confidence for any comparison group built from fewer than 100
matched ions, and **every group this platform applies holds between 23 and 46** - so that
rule fires on every single record, whatever the measurement is like. `qualified` is
therefore the best grade available at present, and `weak` is usually just that one step
lower because a second check could not be run at all.

**What would lift it** is more overlapping data, not a change to the rules: one comparison
group reaching 100 matched ions makes `supported` reachable for ions already in our corpus.
The threshold has deliberately not been lowered to make the top grade appear, because
adjusting a check until it passes is fitting the check to the data. For a molecule we have
never measured, `qualified` stays the ceiling however large the corpus becomes, because the
check described in point 6 can never run for it.

### What to do with a number from this platform

Use it with the range and the grade that came with it, for comparing instruments within this
kind of study. Do not quote it as a reproducibility figure, do not quote the range without
its 80 per cent figure, and do not assume it transfers to molecules unlike steroids - nothing
here has been tested on proteins, peptides or sugars.

Everything below this section is the detailed version, written for somebody working on the
code.

---

What this repository does not do, does not know, and must not be read as claiming.
Written at M0 and updated through each milestone since; sections 7A to 7D were added
after it. Everything here is a live limitation unless it says it is closed.

---

## 1. The thing that matters most

**CLOSED as of 19 September 2026, and this section previously said the opposite.**
Until real data arrived it read "There is not one cross-platform matched ion in this
repository", which was true of the seed corpus and is no longer true of the
repository. The correction is recorded rather than quietly made, because that
sentence was the headline of this whole document.

There are now **142 cross-platform matched ions**, all from one source: the steroid
interplatform study, DOI 10.1021/jasms.2c00196. Section 7D describes what they are and
what is wrong with them.

Two histograms, and it matters which one is quoted:

| over | two platforms | three | four |
| --- | --- | --- | --- |
| all 521 records the file holds | 2 | 43 | **97** |
| the 517 that clear the gate, the only set that may train | 2 | 47 | **93** |

The difference is four records the shared-peak check holds, each of which drops its ion
from four platforms to three. **M4 is fitted on the 93.**

What has NOT changed is the shape of the limitation, only its size:

- **one source.** Every matched ion in the repository comes from a single paper by a
  single group. Cross-platform bias measured within one study is bias between that
  study's instruments, which is not the same quantity as bias between platforms in
  general, and nothing here separates the two. CCSbase is registered and not ingested;
  the Bush Lab database is registered and PARTLY ingested - 1,437 MicroSource records as
  a reference library that contributes no matched ion and cannot, so this limitation is
  untouched by it.
- **one compound class.** 87 steroids. Nothing about a glycan, a peptide, a protein
  or an antibody has a matched ion, so the biopharmaceutical identity layer built in
  M1 has still never seen real data.
- **no pooling across sources is possible yet**, even once the other two are
  ingested, because these compounds are identified by this dataset's own names
  rather than by structure. See 7D.
- **the two Struwe seed files still pair nothing, for the three independent reasons
  below**, and that has not changed either. They remain the only glycan data here.

### The seed corpus, unchanged

The 117 measurements in `data/seed/struwe*.csv` yield zero matched ions, and it is
over-determined - each of these alone is sufficient:

1. both seed files are TWIMS, so there is only one platform;
2. the analyte sets are disjoint — the 2015 file is high-mannose N-glycans
   (Man3 to Man9Glc), the 2016 file is milk oligosaccharides (LNH, LNnH, LNT,
   LNnT). No molecule appears in both;
3. the 2015 rows carry `drift_gas = UNSTATED` and the 2016 rows carry `He`, and
   gas is part of the matched-ion key, so even a shared analyte could not pair.

Reason 3 became stronger on 19 September 2026. A record whose drift gas the source
never stated now keys UNIQUELY TO ITSELF, so the 2015 file holds no key capable of
pairing with anything - not 46 keys shared among its 89 rows, but 89 unmatchable
ones. `matched_ion_keys_held` is 0 for that file and `unmatchable_held` is 89. It
does not merely fail to pair with the 2016 file; it cannot pair at all until
somebody reads Hofmann 2014 and records the gas. See 4.6.

`readiness.py` reports that as a blocker, in those terms, rather than as a small
number.

## 2. No number in this repository describes model performance

No harmonization model has been fitted. No CCS value has been corrected,
harmonized or predicted. Every readiness report carries a `MaturityStamp` of
`provisional`, and the count beside it is the number of cross-platform matched ions
behind it - which was zero when this was written and is 142 for the steroid source.
There is no code path that can produce a `validated` stamp, because the milestone
that would earn one does not exist.

The cross-platform statistics in section 7D ARE computed over real measurements and
are `quotable`. They describe published values of real instruments, which is what
they say they describe. None of them is a model, a correction, or a prediction.

## 3. Licences: what is settled and what is not

**The commercial-versus-research question is ANSWERED.** On 19 September 2026 the
CEO stated that this is academic and research use, not commercial. That is recorded
as `reuse.PLATFORM_USE_CONTEXT = UseContext.ACADEMIC_RESEARCH`, a named constant and
deliberately not a parameter, and it is what reopened three sources: CCSbase, the
Bush Lab CCS database, and the steroid interplatform study.

Those three are registered as `academic_only`, **not** as `open_attribution`. The
distinction is load-bearing rather than pedantic: a reuse status says what a
SOURCE'S TERMS permit, and the context says what WE are. Setting
`PLATFORM_USE_CONTEXT` back to `COMMERCIAL` refuses all three again immediately,
without anybody re-reading a licence page. Each of the three also carries a
`context_basis` field naming who decided the platform may use it and when -
required by `sources.py` for exactly the entries whose usability rests on a decision
rather than on their own terms, because otherwise they read like entries whose terms
permit anybody.

Five of the eleven registered sources are usable: two `open_attribution` and the
three `academic_only`. Two remain **unverified** (METLIN-CCS, and the 2026 Bayesian
paper) and the default-deny gate blocks both. Two are non-commercial-no-derivatives
and would be blocked whatever the context. Two are excluded outright. Each unverified entry carries a `what_to_check` field naming the
single thing that would settle it.

The steroid study's licence is recorded CONSERVATIVELY. Its PMC record carries an
ACS AuthorChoice banner and Europe PMC reports "cc by" with `isOpenAccess Y`; both
were read by this code and NEITHER is accepted as evidence, because CONTEXT.md's own
standard lists a Europe PMC open-access flag among the things that do not count, and
a banner is not licence text. If the article really is CC BY then `academic_only`
understates it and `open_attribution` is correct - understating a licence costs
nothing while this platform is academic, and overstating one cannot be undone.

Nothing in this repository has read a licence. Every entry in `sources.py` records
a report by a named person on a stated date. The code does not verify those reports
and does not claim to.

## 4. Known defects and weaknesses, carried over deliberately

### 4.1 The lone-conformer check has a hole where no locator is given

`loader._lone_conformers` groups sibling conformers by `(calibration group,
analyte, declared total, source_locator)`. Where a file gives no locator, rows
from different samples pool into one key and the check weakens to requiring that
the indices present cover the set and appear equally often. Two samples of one
molecule that have each lost a *different* sibling satisfy that count and are not
caught.

This is inherited from the glycan platform, where it was documented as a
weakening but described as still catching a missing sibling. It does not, in that
one case. It is recorded here rather than fixed because fixing it means deciding
what identifies a sample, which is a schema question for M1.

**It currently bites nothing**, because the 2015 seed conversion folds sample
origin into the locator exactly as the original adapter did. That folding is
load-bearing: the eleven conformer pairs in that file are four compounds measured
from four sample origins, and the locator is the only thing keeping one origin's
pair from pooling with another's. There is a test pinning it.

### 4.2 The conformal calibration floor is derived, tested, and not yet on a fit path

`smallest_calibration_set`, `calibration_refusal` and `calibration_warning` are
correct and tested: 9 calibration points for 90 per cent coverage, 19 for 95.
`interval_readiness` wires them into one reachable call.

Nothing calls it, because nothing computes a prediction interval yet. In the
glycan platform the equivalent functions existed, were correct, and were called
by nothing but their own tests — the floor reached a human only as a line in a
summary. **Whatever computes a prediction interval in M4 must call
`interval_readiness` first and refuse on a refusal.** That obligation is written
down here because an uncalled guard is the failure shape this project keeps
meeting.

### 4.3 The mutation CLI is unhelpful on an empty catalogue

With no mutations, `python -m tools.mutation` prints "no mutation label contains
any of []" and exits 1, which describes a filter that matched nothing rather than
an empty catalogue. Latent only: the catalogue holds 282, so nothing reaches it.
Left alone because it is cosmetic and unreachable, not because the harness is
untouchable: see the entry in section 4A for a case where the opposite call was
made.

**Do not "fix" this by falling back to running everything.** That change was made
during M0, by an agent writing tests, and it is worth recording because of how it
failed. Replacing the refusal with `chosen = list(MUTATIONS)` made one test pass
and turned a filter typo into a full sweep - which is precisely what the comment
three lines above it warns against ("Selecting on the remaining arguments and
ignoring the rest would quietly start a full sweep in answer to a typo - this
tool's own failure mode, turned on itself"). The consequence was a test suite that
ran for two hours and thirty-eight minutes instead of three seconds, because a test
calling the CLI with a non-matching filter now kicked off a complete mutation sweep,
each mutation of which ran the whole suite again. It was found by timing, not by a
failing test. The line is restored. `tools/mutation` is otherwise the original,
byte for byte apart from the package name and the one deliberate change recorded
in section 4A.

### 4.4 Two agents mutating one working tree can make a mutation permanent

Also from M0, also worth recording. Two agents ran mutation checks against the same
checkout at the same time; one read a file its sibling had already mutated, took
that as the pristine text, and wrote the mutation back. It was caught because the
anchor check independently proves every `find` string is present and every `replace`
string absent, and because a second agent noticed the source did not match what it
had read earlier. There was no git history to recover from, which is the part that
made it dangerous rather than merely annoying.

The harness itself is not at fault: it copies the package to a temporary directory
and never writes the repository. The hazard is in doing mutation work by hand over a
shared tree. Serialise such runs, or give each worker its own worktree, and commit
before starting. It recurred during M1 and M2, without damage: one agent solved it
correctly on its own by shadowing the package into its scratchpad rather than
editing in place, which is what the harness itself does and what anyone doing this
by hand should copy.

### 4.5 THE RECURRING CLASS: a guard that looks tested and is not

Twelve separate instances in this repository so far, in twelve different shapes. They
are collected here rather than filed apart, because the shape is the point: in every
one, the suite was green, the coverage looked complete, and a behaviour nobody was
actually protecting could have been deleted without a single test going red.

**One: a guard removed to make a test pass.** `runner.py`'s refusal on a filter that
matches nothing was replaced with a fallback that runs everything. One test went
green, a typo became a full sweep - exactly what the comment three lines above
warns against - and the suite went from three seconds to two hours thirty-eight.
Found by timing, not by a failing test. (Section 4.3.)

**Two: guards with tests and no mutation.** Nine biopharmaceutical identity guards
were named in the brief, implemented, and tested, and the catalogue anchored none of
them. Every test could have been deleted and the sweep would still have reported
that every mutation did what it was expected to do. Found by an agent listing what
it had tested against what the catalogue held, and noticing the two lists differed.

**Three: a mutation that could never be killed.** The catalogue entry for "the
refusal stops returning the measurements it was given" replaced only the opening
of a function argument and left the generator body dangling, so the mutated file
did not compile. An uncompilable mutation is UNRUNNABLE, and this harness
deliberately does not count UNRUNNABLE as a kill - so it was a mutation no test
could ever kill. This one is the LEAST bad of the four, and the reason is worth
noting: the sweep would have reported it as a problem forever rather than quietly
passing, so the harness's own design caught it by construction. It was still
wrong, and it was found by an agent reasoning about why a mutation would not die
rather than by the sweep.

**Four: a fixture that made a guard untestable.** The conformer-pair fixture gave
its two peaks different CCS values, so the VALUE alone kept the two records apart
and the conformer index in the deduplication key never did any work. The mutation
that removes the index from that key was in the catalogue, it was being killed, and
it was being killed by the value rather than by the thing it was written to test.
Found by an agent reasoning about why a mutation died rather than being satisfied
that it did.

**Five: a fixture helper that silently merged two molecules.** `synthetic_inchikey`
first stripped digits, so twenty-four benchmark ions named `bench00` to `bench23`
became one molecule and the whole corpus collapsed to a single matched ion. Fixed,
and then it did it again: it truncated to fourteen letters, so any two tags
agreeing in their first fourteen collapsed the same way, quietly shrinking a
stratum to one point. Both times the symptom was a count that was obviously wrong
if anybody looked - "1 matched ion considered" where twenty-four were expected -
and nothing failed. Long tags are now folded rather than cut.

**Six: two guards whose every test used data where the wrong answer is the right
one.** When the platform's use context arrived - the CEO's answer that this is
academic and not commercial work - `can_train_commercial` was split from `can_use`,
and every call site was meant to be reviewed. Two were missed: `matching.blockers`
and the unusable-source rule in `grading.py`. Both went on asking whether a
COMMERCIAL platform could use a record, in code that meant "may THIS platform use
it". The suite stayed green at 1,784 tests, and the reason is the whole of this
entry: every matching and grading test used `synthetic_fixture`, `open_attribution`
or `unverified`, and all three give the SAME answer to both questions. The one
branch where the two differ had no test at all.

It was found by loading real data. The steroid interplatform study is
`academic_only`, so the licence gate cleared most of its 521 records and matching
then reported all 142 matched ions as blocked by a licence - two modules
contradicting each other about the same records in the same run. Had the first real
dataset been open-licensed instead, this would have sat undetected until the first
academic-only source arrived, which is to say until the numbers mattered.

The durable fix is not the two-line correction. It is
`test_no_module_outside_the_three_that_should_asks_the_commercial_licence_question`,
which asserts that only `reuse.py` and `sources.py` CALL the commercial predicate at
all. A new call site now fails a test, so choosing between the two questions has to
be a decision somebody writes down rather than a default nobody notices.

**Seven: a rule applied to one unknown and not to its twin.** An unstated charge
carrier made a matched-ion key unique to its own record, so two such records could
never match. An unstated DRIFT GAS did the opposite and matched every other
gas-unstated record, across platforms. Both are "the source did not say", both sit
in the same key, and the argument for keying the first uniquely applies to the
second word for word. There was no test either way for the gas, because nobody had
noticed there was a question. Section 4.6 has the detail. Found, again, by real data:
CCSbase states no gas for 25,020 records across three platforms, so the rule was one
ingest away from manufacturing cross-platform pairs in bulk.

**Eight: a guard that existed only in its own documentation.**
`SmallMoleculeAnalyte.dataset_compound_id` carries two paragraphs explaining that the
namespace is the whole safety property - that within a dataset the id is a real
identity and across datasets it matches NOTHING because the namespace differs - and
it had no validator. A BARE COMPOUND NAME was accepted, and two bare names written by
two different adapters matched each other. That is compound-name matching, which this
package refuses everywhere else, arriving through the one field built to prevent it.

Nothing would have failed. The key is well formed, the record validates, the pair
looks like a cross-source matched ion, and the resulting figure would have been
reported as inter-platform bias. The reasoning was written down, reviewed, and
believed, and belief is not enforcement.

It was found by looking for the hole deliberately, after the owner said "do not bridge
on compound names" - that is, by asking what would actually stop it rather than
checking whether the intention was recorded. That is the only one of the eight found
by reading a guard's documentation and then testing the guard against it.

**A NINTH, in the same session and the same shape, worth counting separately because
the mechanism differs.** `calibration_warning` fired only at a calibration set of
exactly nine, the conformal floor. But the quantile index is ceil((n+1)(1-alpha)), and
that index equals n - making the quantile the largest observed score and the interval
the full observed range - for every n from nine to eighteen. The repository was silent
for ten through eighteen, which is precisely where M4's strata land, and a test
asserted that silence at ten. So the wrong belief was pinned rather than merely held:
correcting the code turned three tests red. Section 4.7 has the detail.

**M4 PRODUCED THREE MORE AND THE HARNESS CAUGHT ALL THREE BEFORE THEY SHIPPED**, which is
the first time that has happened and is worth recording as the process working rather
than as three more failures. The sweep reported them as survivors:

- a test for the Passing-Bablok shift overflow that branched on the result and asserted
  something true in BOTH branches, so a clamped answer satisfied it;
- leave-one-COMPOUND-out and leave-one-ION-out being indistinguishable on this corpus,
  because every compound contributes exactly one ion per stratum - the same shape as
  instance four, a guard whose case does not occur in the data it was written against;
- a correction matched on the platform alone rather than on the calibration group, which
  every existing test passed because all three adducts' corrections are small enough to
  look plausible when swapped.

All three now have tests that fail when the guard is removed, and all three needed a
FIXTURE rather than the real corpus, because the real corpus cannot distinguish the
correct behaviour from the broken one. That is the lesson repeated: the data you have
does not exercise the guard you wrote.

**M5 PRODUCED THREE MORE AND THE HARNESS CAUGHT ALL THREE AGAIN.** Two were contract
validators - a measurement with no harmonized value and no reason, and one carrying both a
value and a reason - which every endpoint test passed because those tests build VALID
responses and check what they carry. Nothing constructed the invalid state, so nothing
covered the validator that forbids it.

The third is worth its own sentence because the test meant to cover it existed and did not.
`build_default_model` returns None both when no records load AND when records load but yield
no applicable correction. The test passed an EMPTY directory, which returns None one branch
earlier - so the branch being tested was never reached, and the mutation that breaks it
survived. Reaching it needs data that loads, clears the gate, and still pairs nothing, which
is what the Struwe seed files are. **A test that exercises an earlier return is not a test of
a later one**, and the only thing that distinguished them was the sweep.

Running tally: nine instances found after the fact, six caught by the sweep before shipping
(three in M4, three in M5). The harness is now catching them faster than they are written,
which is the only acceptable direction.

WHAT THE NINE HAVE IN COMMON, and what to do about it. Green is not evidence. The
question that catches all nine is not "do the tests pass" but "what would have to
break for this to go red, and is that the thing I think I am protecting". Note that
three were found by somebody asking that question about a mutation that had apparently
behaved correctly, two by real data contradicting the code, two by reading a guard's
own documentation and then testing the guard against it, and NONE by a test going red.
Concretely:

- when a mutation is killed, check WHICH test killed it and whether that test is
  about the guard. A mutation killed incidentally is a mutation with no cover.
- when a guard is added, add its mutation in the same change. A test without a
  mutation is a test that can be deleted silently.
- when a fixture is built for a guard, make sure the guard is the ONLY thing keeping
  the case from passing. If anything else in the fixture would also separate the
  records, the fixture tests that other thing.
- never fix a failing test by changing the code it is testing, without first
  establishing which of the two is wrong.
- when a mutation is written, check the mutated file still COMPILES. An
  uncompilable mutation is UNRUNNABLE, which this harness does not count as a
  kill, so it is a mutation nothing can ever satisfy.
- when a fixture builds several cases, check the cases are actually distinct. Two
  of these six were one helper quietly giving two different things the same
  identity, and in both the symptom was a count that was obviously wrong to anyone
  who looked at it.
- when a predicate is SPLIT in two, find the inputs where the two disagree and test
  the call site on those. A test whose data answers both questions identically
  cannot tell which question the code is asking - and every existing test is, by
  construction, data that passed before the split.
- when a rule is written for one "the source did not say" case, ask which OTHER
  fields have the same case, and whether the same rule holds for them. Two of these
  were one rule that should have been two, and in both the second place had no
  test because nobody had noticed there was a question to answer.
- **when a docstring explains why something is safe, go and find the code that makes
  it safe.** If there is none, the docstring is the guard, and a docstring stops
  nothing. Two of these nine were exactly that: a namespace described as load-bearing
  with no validator behind it, and a degeneracy described correctly at one value and
  enforced only there.
- when a guard fires on a THRESHOLD, check the whole range the threshold was meant to
  cover rather than the single value in the test. A boundary condition tested at one
  point is a boundary condition believed, not established.

**Ten: a rule that could not answer, reported as a rule that passed.** CLOSED 20 September
2026, and the worst of the ten because the system was making a positive claim rather than
merely failing to check. `flagged_as_an_outlier` is a CORPUS LOOKUP: it asks whether the
submitted ion is among the ions this stratum already recorded as not transferring. For an
ion that is not in the corpus there is nothing to look in, and it returned `None` - the same
value it returns for an ion that IS in the corpus and is fine. `grade_correction` read
`None` as "checked and clean". `not_checked` came back empty, and an empty `not_checked` is
used throughout this API as a positive assertion that every rule was evaluated. It is not
vestigial: it populates on 133 of 417 graded responses from the seed corpus, which is what
makes the empty case a claim rather than an absence. (417 is the number of measurements that
get a correction computed and graded; 402 of those are served a value, the other 15 grading
`unsupported`. The two are different counts and neither stands in for the other.)

So the only rule in the scheme that speaks about the SUBMITTED ION rather than about the
stratum around it was inert for every ion not already in the corpus - which is the entire
population this platform exists to serve - while the response asserted it had been applied.
Found by renaming a corpus ion and resubmitting it: the service returned a value it had just
refused for the same measurement under its real name, with a 90 per cent interval that
excluded that ion's own true DTIMS values.

The fix is in three parts, because the first two alone would have left the same hole one
level up. Evaluability is decided BEFORE the rules run, so "found nothing" and "could not
look" are different answers. An unevaluable rule is named in `not_checked` and demotes the
grade one notch - not to `unsupported`, because the ion may be perfectly fine, but a grade
resting on fewer checks than the scheme advertises is not the same grade. And the demotion
itself is published in `/confidence/rules`, because a demotion a caller can receive but
cannot look up is a demotion they cannot argue with.

The third part caught a fresh instance of the same class on the way in. `applies_to` - the
field carrying the scope limit - was added to the scheme and served to nobody, because the
endpoint hand-copied five fields by name and the test that checked the endpoint against the
scheme hand-listed the same four. Both passed, because both only checked the fields they
already knew about. The endpoint now builds the response model straight from the scheme
dict, `ConfidenceRule` sets `extra="forbid"` so the next added key raises instead of being
dropped, and the test compares the two key sets rather than a list of names.

Effect on the seed corpus, measured rather than estimated: grades moved from 395 qualified /
15 weak / 7 unsupported to **186 qualified / 216 weak / 15 unsupported**, and values served
fell from 410 to 402.

**Do not carry the corpus figure forward as the real one.** Replaying the seed corpus reports
the outlier rule as unevaluable on 124 of 417 graded responses - about 30 per cent - and that
number flatters itself, because every record in the replay is BY DEFINITION already in the
corpus it is being looked up in. For the use this platform exists for, an ion it has not
measured, the rule is inert on **100 per cent** of requests. The corpus figure is not a
smaller version of the true one; it is a measurement of a different population, and 30 per
cent reads like a minor gap where the truth is that the scheme's only ion-specific rule never
applies to the intended user at all.

- when a predicate returns None, ask whether it has TWO reasons to do so - "I looked and
  found nothing" and "I could not look" - and whether the call site can tell them apart. If
  it cannot, the caller is being told the stronger of the two. This is the split-predicate
  lesson above, arrived at from the other end: there the split was visible in the source and
  the call site was untested; here the split was invisible because one function was quietly
  answering two questions with one value.
- when a response field means "nothing to report", check whether anything ever populates it.
  If it does, an empty one is a CLAIM, and it needs the same evidence as any other claim the
  response makes.
- when a rule is keyed on identity, ask which population it can answer for, and whether that
  population is the one the software is for. A check that works perfectly on the data you
  have and not at all on the data you expect is worse than no check, because it reports.

**Eleven: a hand-copying constructor that silently omits, and the `extra="forbid"` that
did not stop it.** CLOSED 21 September 2026, found in an adversarial pass over the fixes for
the other ten, and it is the THIRD occurrence of one shape.

`ScopeReport.of` built the served scope by copying five fields off `ScopeStamp` by name. When
`pairings` was added to the stamp on 20 September - as part of correcting the measurement
count from 1402 to 517 - the constructor did not copy it, `ScopeReport` had no such field,
and no structured field carried the number. The `caveat` STRING inside the same object went
on saying "31 cross-platform pairing(s)". So a caller parsing JSON got strictly less than the
prose sitting beside it, and had to scrape a sentence for a count the object could simply
have carried.

**The three occurrences are one shape.** Instance Ten was `applies_to`: added to the grading
scheme and dropped by an endpoint that hand-copied five fields, past a test that hand-listed
four. Then the same endpoint's fix - build from the source dict, set `extra="forbid"` - was
applied, and this happened anyway, one module away.

**WHY THE FIX FOR TEN DID NOT PREVENT ELEVEN, and this is the part worth remembering:
`extra="forbid"` REJECTS UNKNOWN KEYS, NOT MISSING ONES.** It guards the direction where
something offers the model a key the model does not know. It cannot guard the direction where
the SOURCE grows a field and nothing ever offers it, because there is no input to reject -
the key simply never appears. Forbid looked like protection for both directions and is
protection for one. That asymmetry is invisible when you write it, because the failure it
does catch is loud and the failure it does not catch is silent.

The fix is in two parts, and the second is the one that generalises:

- `of` builds from `dataclasses.fields(stamp)` rather than naming fields, so a field added to
  the stamp IS offered - which turns a silent omission into the one thing forbid is good at,
  an unknown key, loudly refused;
- two tests compare the field sets IN BOTH DIRECTIONS. Every stamp field must be a report
  field or be named in `OMITTED_FROM_THE_STAMP`; every report field must be a stamp field or
  be named in `DERIVED_HERE`. The first is the direction nothing was ever guarding.

- **when a field set is copied from one type to another, assert the two SETS, in both
  directions, and name the deliberate omissions in code.** A constructor that lists fields is
  a constructor that will be out of date, and the test that checks it will be out of date in
  the same way and at the same moment, because both were written by someone looking at the
  same list.
- **ask what a guard is asymmetric about.** `extra="forbid"`, a minimum without a maximum, a
  validator on one spelling of a field and not the other: each of these has caught something
  here, and each has a blind direction that reads as covered.

**Twelve: THE CLASS TRAVELS WITH THE CODE. A hand-written module list, ported in from another
repository, failed on the straight copy.** CLOSED 26 September 2026, during the glycan port.
It is the FOURTH occurrence of the hand-maintained-list shape, after Ten and Eleven, and the
first that was not written here.

`test_training.py` in `Project2` guarded something worth guarding: that importing the package
loads no estimator library, so an accidental fit anywhere in the suite turns the suite red. It
did it by importing every module in a list called `PACKAGE_MODULES` - eleven names, written
out by hand. The glycan port copied fourteen test files verbatim and ran them against the
copy, and exactly one test failed:

```
FAILED test_training.py::test_no_estimator_library_is_imported_by_the_package
  ModuleNotFoundError: No module named 'wmxglycan.api'
```

`api.py` is deliberately not ported, so a guard about estimator libraries failed for a reason
with nothing to do with estimator libraries.

**Why this is an instance and not a porting inconvenience.** The failure was loud here only
because a module went MISSING. The dangerous direction is the other one, and it is silent: the
day anyone ADDS a module to that package and does not think about this list, the guard stops
covering it, the suite stays green, and an estimator import in the new module goes unnoticed.
That is the same asymmetry as instance Eleven, in a different mechanism - a hand-written list
of what a guard covers can only ever be checked against the world by someone looking at both,
and the person who adds a module is not looking at a test about estimators.

The list is now derived from the package directory. **And the derivation needed its own floor,
which is the part worth carrying forward:** a glob that returns nothing would import nothing,
leave the guard green, and be the identical defect one level up - so
`test_the_derived_module_list_actually_covers_the_package` asserts the derived list has at
least fifteen entries and names the six modules the port was commissioned for. Replacing a
hand-written list with a derived one moves the failure mode rather than removing it, unless
something asserts the derivation is not empty.

- **a guard's coverage is data, and data goes stale.** Wherever a test enumerates what it
  protects - modules, fields, call sites, rules - derive the enumeration from the thing itself
  and then assert the derivation is non-empty. Two of the four occurrences of this shape were
  found by accident and one was found by a port.
- **the failure class is not a property of this codebase's habits.** It arrived intact from a
  separate repository written to the same standards. Expect it in anything ported in, and look
  for it deliberately rather than hoping a copied suite surfaces it - this one surfaced only
  because a module was absent, which was luck.

### 4.6 An unstated drift gas keyed records together; an unstated carrier did not

CLOSED 19 September 2026, and recorded because the asymmetry stood for four
milestones without anybody deciding it.

A record whose CHARGE CARRIER the source never named - `[M+24?]24+` - has always
keyed uniquely to itself, so it can never match another record, including another
unstated-carrier one. That was built in M0 on an explicit instruction, and the
reason is that two papers both reporting "the 24+ ion" of one protein have not
reported the same ion: one may be 24 protons and the other 24 ammonium adducts,
408 Da apart.

A record whose DRIFT GAS the source never stated did the opposite. `DriftGas.UNSTATED`
is one enum value, so every gas-unstated record matched every other one - **including
across platforms.** The argument against that is the carrier argument word for word: a
value referenced to helium and a value referenced to nitrogen are not values of the
same ion, which is precisely why the gas is one of the five key components. "Neither
source said" is the absence of evidence, not evidence of agreement.

**It was found by real data, and it was about to do damage at scale.** CCSbase holds
25,020 measurements across DTIMS, TWIMS and TIMS and states no gas per record.
Ingesting it under the old rule would have manufactured cross-platform matched ions
in bulk - a helium drift-tube value sitting against a nitrogen travelling-wave value
as though they were one ion - and every one of those pairs would have fed a bias
figure that looked like a finding.

Now both unknowns make a key unmatchable, the key names which one is missing (so a
curator knows whether resolving the gas is enough), and a record missing both says
both. Nothing narrows for records that DO state their gas: two laboratories on two
platforms both stating nitrogen still pair, which is what the key is for. The steroid
corpus states N2 throughout and its 142 matched ions are unaffected.

What it cost: the 2015 seed file went from 46 shared keys to 89 unmatchable ones,
which is a documented count that changed. The M0 conclusion is strengthened rather
than weakened by it - see section 1.

### 4.7 The conformal floor was enforced at one point of a ten-wide range

CLOSED 19 September 2026, found while sizing M4's strata and before any M4 code existed.

Split conformal takes the k-th smallest conformity score as its quantile, where
k = ceil((n+1)(1-alpha)). There are THREE regimes, and the repository modelled two:

| | |
| --- | --- |
| k > n | no such score exists, so no finite interval: **refused** |
| k == n | the quantile IS the maximum, so the interval is the **full observed range** - valid, correctly derived, and excluding nothing that was seen |
| k < n | informative |

At 90 per cent coverage `k == n` holds for every n from 9 to 18. `calibration_warning`
warned only at 9, so nothing was said for 10 through 18 - and
`test_only_a_calibration_set_exactly_at_the_floor_carries_the_degenerate_warning`
asserted no warning at 10, so the wrong belief was pinned. Correcting the derivation
turned three tests red, which is the signature of a belief rather than an oversight.

**It would have mattered immediately.** M4's strata are 23, 29 and 41 ions grouped by
66 compounds, so a grouped split lands calibration sets at roughly 11, 14 and 20 - two
of the three inside the silent range. M4 would have quoted full-observed-range
intervals as 90 per cent prediction intervals, with no warning attached.

Two derived sizes now exist, both computed from alpha rather than tabulated:

- `smallest_calibration_set` - an interval EXISTS: 9 at 90 per cent, 19 at 95.
- `smallest_informative_calibration_set` - it is narrower than the observed data:
  19 at 90 per cent, 39 at 95. Solving ceil((n+1)(1-a)) <= n-1 gives n >= (2-a)/a.

The warning now tests the quantile index itself rather than comparing against a
remembered number, so the floor, the refusal and the warning cannot drift apart.

### 4.8 THE MUTATION SWEEP MEASURES A SUITE WITH NO SERVED MODEL

OPEN, found 26 September 2026 while extending the harness to a second package, and it has been
true of **every mutation sweep ever run in this repository** - including the 282-of-282 at
`v0.7.0-mvp`. It is not caused by the glycan work; the glycan work is what noticed it.

The harness copies `src/wmxccs` into a temporary directory and puts that first on `PYTHONPATH`.
`api.py` locates the corpus relative to its own file:

```python
SEED_DIRECTORY = Path(__file__).resolve().parents[2] / "data" / "seed"
```

Inside the shadow, `__file__` is `<temp>/shadow/wmxccs/api.py`, so `parents[2]` is `<temp>` and
`<temp>/data/seed` does not exist. **`build_default_model()` therefore returns `None` during
every sweep.** Measured, not inferred:

```
SEED_DIRECTORY : C:\Users\...\Temp\probe2-eqsrmwmz\data\seed
exists         : False
build_default_model() -> None
```

WHAT IS ACTUALLY LOST, and it is much less than that sounds

One test. The shadow baseline reads `3841 passed, 1 skipped` where a direct run reads `3842
passed`, and the skip is
`test_entry_point.py::test_the_banner_carries_the_model_version_an_answer_can_be_reproduced_from`,
which skips itself when the banner says `NO MODEL LOADED`. So during a sweep nothing checks that
the startup banner carries the model version, the scope line, or
"NOT interlaboratory reproducibility".

The 26 other call sites of `build_default_model` are safe, and one of them is safe on purpose:
`test_api_harmonize.py` passes the seed directory **explicitly**, computed from the TEST file's
location rather than the package's, and then asserts

```python
assert model is not None, "the seed corpus must yield a model or this whole file is vacuous"
```

Somebody had already met this trap and guarded it. That assertion is the reason twenty
harmonization tests measure a real fit inside the shadow instead of quietly measuring the
501 branch.

IT CAUGHT A NEW TEST THE SAME DAY. `tests/test_glycan_boundary.py` asserts that fitting the
served model loads no `numpy`, `networkx`, `glycowork` or `sklearn`. Written with a
zero-argument `build_default_model()`, it passed inside the shadow because **nothing was
imported at all** rather than because the right things were not - green in both places,
meaningless in one. It now passes the seed directory explicitly and asserts the model is not
None, exactly as `test_api_harmonize.py` does.

WHY IT IS NOT FIXED HERE

The fix is small - expose `data/` inside the shadow - and it changes the environment the
282-of-282 figure was measured in, so taking it would require a fresh full sweep of roughly six
hours to re-establish that number. The owner asked for the ranker and to stop after it. Two
things make deferring defensible rather than convenient:

- the direction of the error is towards LESS coverage, never towards a false kill. A test that
  cannot see a model either skips or exercises the no-model branch; neither can make a mutation
  look dead when it is alive.
- every mutation in the catalogue was killed at `v0.7.0-mvp`, so no mutation depends on the
  skipped test as its only killer. The loss is bounded to that test's own assertions.

**The general lesson is the one worth carrying:** a test harness that relocates the package
relocates every path the package derives from `__file__`, and anything the package finds that
way silently disappears. Any guard written against a resource located relative to the package
rather than to the repository is a guard that does not exist during a sweep.

## 4A. Defects found during M0 and fixed

Recorded because each was a real hole, and because the test that found it is the
test that keeps it shut. None of these is outstanding.

| What was wrong | Consequence had it shipped |
|---|---|
| `AntibodyIdentity` was listed in `component_records()` and carries no reuse status | The gate refused EVERY intact-antibody and ADC measurement whatever its licence, closing the biopharmaceutical layer entirely - the platform's stated differentiator |
| `GlycanAnalyte.identity_key` fell through to `self.composition.canonical` on a record whose composition is null | A glycan stating only an unverified-form accession, or a structure naming no linkage, loaded, cleared the gate, sat in the corpus, and raised `AttributeError` the first time anything asked for its key - in M2, over real data, naming neither the record nor its file |
| `SourceLicence` coerced `evidence` with `tuple()` | A bare string became one piece of evidence per character. "RSCarticlepage" passed every check as seven pieces of evidence, leaving the entry reading as thoroughly evidenced while holding nothing |
| The loader caught bare `Exception` around `_build_analyte` | A row naming an analyte kind the union does not hold was counted as failed validation rather than under its own reason, making a wrong column indistinguishable from bad data |
| `MeasurementLoadReport` had no `uncertainty_types_held` | A file where nothing clears said nothing at all about its spreads. For the 2015 seed an uncertainty type of `unknown` is one of the two things blocking all 89 rows, so the report was silent about half its own headline |
| `runner.sweep` and the `--check` path read every target file with no error handling | A mutation naming a module that had gone away raised `FileNotFoundError`, taking the whole run down before a single mutation was judged and losing every other verdict with it |

**The runner fix is the one place `tools/mutation` deliberately departs from the
original.** The port was byte-for-byte apart from the package name, and that
faithfulness was the right default: it is what made the removed guard in section
4.3 detectable by a one-line diff. It was the wrong default here. A module can
vanish exactly as a line can - renamed, split or dropped between milestones - and
both are the same failure the harness exists to make loud: a mutation that has
silently stopped applying. Crashing named neither the mutation nor the remedy and
cost every other result in the run, which is loss of coverage in its purest form,
inside the tool built to catch loss of coverage. It now reports STALE and says the
module went away rather than that a line moved. Two tests hold it, one of which
asserts that a missing file costs exactly one result and not the other 109.

## 5. Where this repository departs from CONTEXT.md, and why

CONTEXT.md has since been corrected in place for the factual errors listed in
section 6. These are the remaining places where its instructions could not be
followed exactly. Each was a judgement call and each is reversible.

| CONTEXT.md says | What was done | Why |
|---|---|---|
| Port `licensing.py`, "no change beyond the package name" | `reuse.py` ported as well | `ReuseStatus` and the tier predicates live in `reuse.py`, not `licensing.py`, which only re-exports them. Folding them together re-creates the import cycle `reuse.py` exists to break |
| Port `sources.py` | Ported without its SugarBase dataset entry | `sources.py` imported `sugarbase.py`, which is on the do-not-port list, and called it at module import. The `DATASETS` mechanism is kept and seeded with the two seed transcriptions, so the dataset-backed branch of the gate stays exercised instead of becoming unreachable |
| Glycan analyte: "reuse the existing representation" and do not port `composition.py` | The canonicalisation half of `composition.py` was ported into `identity.py`; the N-glycan plausibility rules and the Man3GlcNAc2 core check were not | Those two instructions conflict. Canonical spelling is identity machinery — without it `Hex5HexNAc4` and `HexNAc4Hex5` are two matched ions instead of one. The biosynthetic rules are glycan biology and do not come across |
| `readiness.py` "ported, rethresholded" | Rewritten around matched-ion and platform counts | `Readiness` lives inside `training.py`, whose thresholds come from `splits.py` and `evaluation.py`, neither of which is in the port list and both of which import do-not-port modules. The thresholds that carried over are marked as ported; the platform-count blocker is new, and is the one this platform actually needs |
| `loader.py` "adapted to new columns" | One canonical format, plus a separate conversion step in `tools/seed_struwe.py` | The glycan platform had per-source adapter modules. A second entrance to the record model is a second way to skip the gate, which its own loader docstring warns against. Each source is now converted once into a reviewable file and then loaded by the one loader |
| `reducing_end_label` and `derivatisation` on the measurement | Moved onto the glycan analyte | A steroid has no reducing end. A field that can only ever be null on five of six analyte kinds is a glycan assumption leaking into every record |

Two further decisions that CONTEXT.md does not cover:

- **Composition is optional on a glycan** when a structure identifier is present.
  The 2016 transcription has no composition column; the glycan platform derived
  one from the IUPAC string using a monosaccharide-class table, which is glycan
  chemistry and does not port. Rather than derive or invent one, those rows record
  no composition and are identified by their structure, which is finer anyway.
- **`[M+24?]24+`** is a new adduct form meaning "twenty-four charges, carrier not
  stated". Native-MS papers routinely report a charge state without saying whether
  the carriers are protons, sodium or ammonium, and the Bush Lab protein data is
  expected to be full of them, so this is the common case rather than the exotic
  one. Without this the adduct is mandatory and its charge must match, so such a
  record could not be built at all; writing `[M+24H]24+` instead would assert
  protons, which the paper does not say. It blocks training, following the
  `DriftGas.UNSTATED` precedent exactly.

  It also does one thing that precedent does not. **Its matched-ion key is made
  unique to its own record**, so two unstated-carrier records can never match each
  other. Two papers both reporting "the 24+ ion" of one protein have not reported
  the same ion: one may be twenty-four protons and the other twenty-four ammonium
  adducts, which differ by 408 Da and do not have the same cross section. Letting
  the two land on one key would pair them, and the pipeline would report the
  difference between two different ions as inter-platform bias. The protection is
  structural rather than a rule downstream code has to remember: the key carries
  the record's own provenance, so it is equal to nothing but itself. The loader
  counts such records apart from its matched-ion figures, under `unmatchable_held`.

## 6. Corrections to CONTEXT.md

Stated plainly because CONTEXT.md is the document everyone reads first.

- **"1,796 tests" is wrong.** The glycan repository collects **1,868** test cases
  from 705 test functions at the commit read. The string "1796" appears nowhere in
  that repository.
- **"136 curated mutations" is wrong.** The catalogue holds **154**.
- **"Roughly 40 percent of its code carries over" understates it as an import
  fact.** The transitive closure of the named port set is 15 of 20 modules, four
  of which are on the explicit do-not-port list. The glycan core is not a
  separable layer; it is load-bearing under the registry, the record model and the
  loader, and cutting it is most of the porting work.
- **"89 stored records stay blocked until the drift gas is known" is wrong.**
  Resolving the gas alone unblocks nothing. All 89 rows of the 2015 file carry
  `uncertainty_type = unknown` with no spread, which is an independent blocker.
  Both must be resolved. There is a test pinning both so the day somebody reads
  Hofmann 2014 the suite tells them there is a second one.
- **The four held rows in the 2016 file are not marked in the file.** CONTEXT.md
  describes them as held pending curation. Nothing in the CSV says so: the hold is
  *derived* by the loader's shared-peak detector from exact CCS equality within a
  calibration group across differing analyte identity. Port the detector away and
  all 28 rows clear with nothing to show it happened.
- **"~6,500 rows" for the Bush Lab database is wrong, and four of its eight per-sheet
  figures were artefacts.** The workbook holds **1,804 data rows**, not roughly 6,500 -
  an overstatement of 3.6x on the only biopharmaceutical dataset this project has
  identified. CONTEXT.md and LIMITATIONS 7E both gave 1,000 for native-like protein
  cations, denatured protein cations, anionic homopolymers and other peptides; the true
  figures are 72, 27, 35 and 12. A round 1,000 repeated four times is a reader default or
  a pre-allocated range, not a count. The four small figures in the old table (33, 24,
  1,441 and 989) were real but counted the header row, which is the tell: accurate counts
  and artefacts sat side by side and looked alike. Corrected in both documents on
  23 September 2026, with the per-sheet header-row counts shown so the arithmetic can be
  checked. Counting non-empty rows instead of data rows gives 1,815, and treating every
  sheet as having a single header row gives 1,805; three sheets carry a merged banner
  above their column names, so 1,804 is the figure.
- **Calibration reference lineage was not missing from the glycan repository.**
  The objective document has no field for it, but `ccs_is_calibrated`,
  `calibrant_reference` and the primary-versus-calibrated distinction all existed
  and ported. What was genuinely new is the structured reference set: which
  publication the reference values came from and what platform they were measured
  on, which is what makes circularity detectable at all.

## 7. What is not known about the seed data

- The drift gas the 2015 values refer to. The ESI says a dextran ladder of known
  DTCCS and never says which gas. Recorded as `UNSTATED`, not inferred from the
  same group's 2016 paper. Resolving it means reading Hofmann 2014, Anal. Chem.
  86, 10789.

- The platform the calibration references were measured on, for both files.
  "DTCCS" in the 2015 cell plainly suggests a drift tube, and that is a reading of
  an abbreviation rather than of a paper, so `calref_platform` is null on both and
  `traces_to_primary` returns `None` for every seeded record. That is the honest
  answer.
- Whether either transcription is faithful to its ESI. Neither has been checked
  against the source within this repository. Every count reported is a count of
  the transcription as delivered.
- The seven high-mannose structures for the 2015 file. They are in Figure 1 of the
  paper, which nobody here has read. The composition column is taken as
  transcribed rather than derived from a structure, so a composition typo in that
  file cannot be caught the way a structure would catch one.

### 7.1 The 2015 uncertainty type is probably not resolvable at all

All 89 rows of the Struwe 2015 file are blocked by **two** independent things:
`drift_gas = UNSTATED` and `uncertainty_type = unknown`. Reading Hofmann 2014 and
settling the gas unblocks **zero** records on its own.

The second one is very likely permanent. **The owner reports, on 19 September
2026, that the Analyst ESI carries no uncertainty column at all.** If that is so,
`unknown` is not a transcription gap waiting on a more careful reading; it is a
correct and complete record of what the source states, and no amount of going back
to that paper will change it.

**So do not treat these 89 records as pending curation.** They are not a task on
anybody's list. Unblocking them needs one of:

- a different source reporting a spread for the same ions, in which case those are
  different records with their own provenance, not a repair of these; or
- a decision that a CCS value with no reported spread may be used for some purpose
  that does not need one. That is a modelling decision for whoever builds the
  harmonization model, and it must be made explicitly, in the open, and never by
  quietly defaulting the uncertainty type to something usable.

The uncertainty type is not droppable, for the reason the field exists: a spread
of unknown kind is read as whatever the reader assumes, and a two-standard-
deviation figure loaded into a field read as one standard deviation halves every
interval built on it with nothing downstream able to detect it.

This corrects the reference documents, which listed the gas as the single blocker
on these records. Both `CONTEXT.md` and `PROJECT_REFERENCE.md` have been amended.

## 7A. M1 and M2: what is built, and what it has never seen

### The cyclic layer has no data behind it

`CyclicSettings` accepts the fields the objective document asks for and enforces
the rules that follow from them, and **not one cyclic measurement exists in this
repository**. Public cyclic CCS data is very thin, which the objective document
notes in its own risk table, and no internal instrument access has been confirmed.

Two consequences worth stating plainly:

- `arrival_time_correction` is FREE TEXT, where every other controlled vocabulary
  here is an enum. It is not an enum because nobody has read a real cyclic methods
  section, and inventing the members from what is usual would be exactly the kind
  of plausible detail that has burned this project before. When a real dataset
  arrives, read its methods and make it an enum from what is actually there.
- the wrap-around rules are reasoned, not measured. That a multipass value with no
  wrap-around statement is unusable follows from what wrap-around does to arrival
  time; it has not been checked against a paper that reports one.

### M2 is built entirely against synthetic fixtures

`matching.py` has never seen a real cross-platform pair, because there is not one
to see. Everything it does is exercised by `fixtures.py`, whose records all
declare `SYNTHETIC_FIXTURE`, and `assert_quotable` refuses to let any report
covering one be quoted as a result.

That is a deliberate trade: the pairing rules do not depend on which dataset fills
them, and waiting for the steroid study's licence to be read would have put the
deadline at risk for no benefit. But it means **the first real dataset is also the
first test of whether these rules fit real data**, and two of them are judgement
calls that real data may overturn:

- **deduplication keys on the exact CCS value.** Two records agreeing to the last
  digit on one platform are treated as one measurement republished. That is right
  far more often than it is wrong, but a genuine independent replicate that lands
  on the same number would be collapsed. The provenance keeps both citations, so
  the evidence survives; the count does not.
- **conformer indices are never paired across sources.** A matched set carrying a
  conformer index is built and then refused, because conformer numbering is
  source-local and nothing establishes that conformer 1 in one paper is conformer 1
  in another. If a real benchmark numbers its conformers consistently, this refuses
  pairs it should have found - which is the safe direction, and still a cost.

### Two pairing decisions that real data may argue with

Both were deliberate, both are tested, and both are the sort of thing that only a
real dataset can settle.

**Stepped-field and single-field DTIMS count as two platforms.** One ion measured
both ways is a matched set here. The case for it is that a primary value and a
calibrated one are genuinely different measurements - the published
reproducibility figures are separate numbers for the two - and that folding them
together would hide the primary-versus-derived distinction from any later fit,
whose DTIMS arm would then be regressing against a moving anchor. The case
against it is that a benchmark reporting both would produce matched sets a reader
expecting "DTIMS versus TWIMS" did not ask for. If that reading is wrong, the fix
is one line in `matching._platform_of`, and a mutation will tell you the moment
somebody changes it by accident.

**Conformer indices are never paired across sources.** A matched set carrying a
conformer index is built and then refused, because conformer numbering is
source-local and nothing establishes that conformer 1 in one paper is conformer 1
in another.

THIS IS A NAMING PROBLEM, NOT PHYSICS, and the distinction matters enough to state
plainly so that nobody later files it as a bug. Two arrival-time peaks of one ion
are a real, physical thing: one structure folded two ways, confirmed by MS/MS in the
Struwe 2015 data. What is arbitrary is the NUMBER each laboratory hangs on each
peak. One group may number by arrival time, another by abundance, another by the
order the peaks came out of their fitting software. Nothing in a published table
says which. So pairing on the index would not be resolving a hard physical question
slightly conservatively; it would be manufacturing matches out of two independent
labelling conventions that happen to use the same integers.

The consequence is a cost, and it is the right cost: if a real benchmark does number
its conformers consistently, this refuses pairs it could have found, and a person
can confirm the correspondence and lift it. The reverse mistake is not recoverable -
a fabricated pair enters the statistics as evidence and nothing downstream can tell
it from a real one.

### A protein's sequence atom carries its subunit label

`ProteinAnalyte.identity_atoms()` spells its sequence atom
`sequence:<sequence>|<subunit>`, so two records with one sequence and different
subunit labels do not share an atom and will not be merged by matched-ion
construction. That matches `identity_key()`, which holds them apart.

There is a real argument the other way: an identical sequence IS the same
molecule, and a subunit label is only naming. If a real dataset shows that
laboratories label the same chain inconsistently, this refuses merges it should
have made, and the fix is one line. It is the conservative direction - it costs a
merge a person can still make by hand, rather than making one nobody asked for -
and it was chosen because two identity surfaces that can contradict each other
are a latent bug whose only symptom is a pair nobody can explain.

### The fixture corpus is not a benchmark

Eighteen cases, chosen to exercise one rule each. It says nothing about how often
these situations occur, and no number computed over it describes anything. Its
only claim is that each rule fires when it should and stays quiet when it should
not.

## 7B. M3: statistics computed over nothing real

Every figure `statistics.py` can produce has been produced from records declared
in code. The module has never seen a measurement made by an instrument.

### Mass is out of scope by decision, not by oversight

The brief for M3 asks for residuals against mass as well as against CCS, so that
size-dependent bias is visible. **No record in this repository carries a mass**,
so the mass half is not computed, and that is a DELIBERATE M0 SCOPING CHOICE that
has been reviewed and confirmed rather than an omission.

What was decided in M0: `composition.py` came across as its canonicalisation half
only. The residue-mass table, `RESIDUE_MASS`, `WATER_MASS` and the
`monoisotopic_mass` property were left behind with the N-glycan plausibility rules,
because a mass table is chemistry and this platform takes any ion. Carrying a
glycan mass table into a repository that also holds steroids, peptides, antibodies
and ADCs would have given exactly one of the seven analyte kinds a mass and left
the rest without one, which is worse than none having it.

What is computed instead is the residual against the reference CCS. CCS is a
reasonable size proxy for the purpose intended and it is not the same thing: two
ions of one cross section can differ in mass, and a bias that tracks mass rather
than size would not show.

**What adding it would require**, if a later milestone wants it:

1. A FIELD on `CCSMeasurement`. Almost certainly two, not one: a neutral
   monoisotopic mass and an m/z are different quantities, and an m/z without its
   charge is not a mass. The charge is already on the record, so `mz` plus the
   existing `charge` would give a derived mass - but see 3.
2. A PROVENANCE DECISION, which is the hard part. A mass read off a published
   table is a transcribed observation and belongs in the loader with every other
   transcribed value. A mass computed from a formula or a sequence is a DERIVED
   number, and this platform has no field that distinguishes the two. Storing a
   derived mass in a field a reader assumes was transcribed is precisely the
   failure class recorded in section 4.5. Whichever is chosen, the record has to
   say which it is.
3. A PER-KIND SOURCE for the derived case, and this is where it gets expensive:
   a peptide mass follows from its sequence, a glycan mass from a residue table,
   a small-molecule mass from its formula, and an antibody mass from neither. Six
   of the seven analyte kinds would need their own derivation, each with its own
   reference data and its own licence question.
4. A LOADER COLUMN and a validator refusing a mass that disagrees with the
   adduct's charge by more than rounding.

None of that is hard. All of it is a schema and provenance decision that should
be made when something actually needs mass-based residuals, rather than
speculatively now.

### The outlier rule was wrong the first time, and the second version is still a choice

Outliers are flagged more than 2 per cent from **their own stratum's median**
difference, not from zero. The first version used an absolute threshold and was
demonstrably wrong: over a corpus with an injected 2 per cent bias it flagged 14
of 24 ions, burying the three that genuinely did not transfer among the ordinary
ones.

The relative version is defensible - the systematic offset is what harmonization
corrects, so what matters is which ions are still wrong after it is removed - but
the margin itself is borrowed from the published TWIMS envelope and is not derived
from anything measured here. A platform pair that genuinely agrees to 0.2 per cent
would have its real outliers hidden by a 2 per cent margin. When a real corpus
exists, the margin should be reconsidered against that corpus's own spread.

### Lambda is assumed equal more often than it is measured

The Deming fit weights by the ratio of the two platforms' error variances. That
ratio is measured only when EVERY paired record on both sides reports an
uncertainty this repository can convert to a standard deviation. In practice that
will often fail, because:

- `CI95` is refused outright. A 95 per cent interval may be reported as a half
  width or a full width and this repository has never fixed which. Dividing by
  1.96 on the assumption it is a half width would silently halve or double every
  weight built on it.
- `SEM` needs the replicate count, which most published tables do not give.
- `unknown` converts to nothing, which is the point of the value.

When the ratio cannot be measured the fit falls back to orthogonal regression and
the report SAYS it assumed rather than measured. That string is the guard, and it
has a mutation of its own, because a fit reporting an assumed weighting as a
measured one is exactly the kind of quiet overclaim this project exists to refuse.

### Limits of agreement rest on the one threshold that is pure policy

`MIN_POINTS_FOR_LIMITS_OF_AGREEMENT = 10` is not derived from anything. The two
other floors are: a correlation needs three points because r is exactly plus or
minus one at two, and a Deming slope needs three because two parameters need one
residual left over. Ten is the conventional floor for quoting an SD of
differences, it is marked as policy in the source, and somebody with a real corpus
should argue with it.

### What M4 must do with the leverage finding

Recorded here because it is a conclusion from M3 that constrains a milestone that
does not exist yet, and it would otherwise live only in a conversation.

The benchmark showed three outlying ions in twenty-four moving a fitted Deming
slope from the injected 1.020 to 1.031. The correction applied to every
well-behaved ion was therefore, in part, the work of the three that were not. So
the harmonization model, when it is built:

- **reports BOTH the slope-derived and the median-derived correction.** The
  contract already has fields for both, and for which one the headline value came
  from. They agree when the fit is well behaved and diverge exactly when it is
  not, so a caller comparing them learns something a single number hides.
- **prefers a robust fit.** Three ions should not drive the correction applied to
  everything else. Robust does NOT mean dropping the outliers: they stay in the
  corpus, in the statistics and in the reported counts, and `grading.py` already
  flags any ion the fit does not cover. It means fitting so that a few points
  cannot dominate.
- **keeps `grading.correction_driven_by_outliers` on the path.** It measures the
  leverage directly and demotes the grade when it exceeds half a per cent of CCS
  at a typical ion. That check is the M3 finding made operational, and a robust
  fit should make it fire less often rather than make it unnecessary.

### The benchmark corpus is constructed, not observed

`fixtures.benchmark_corpus()` injects a 2 per cent TWIMS bias, a 1 per cent TIMS
bias and a 7 per cent tail on three ions. Those numbers were chosen to be
recoverable and to be roughly the shape of the published envelope, so that the
arithmetic is exercised on a case resembling the real problem. **The resemblance
is a convenience and is not evidence.** No figure computed over that corpus
describes any instrument, `quotable` is False for every report built from it, and
`assert_quotable` refuses.

## 7C. The API and the confidence scheme, built ahead of the model

Both were built while waiting for a licence answer, because neither needs data.
Both stop short of anything that needs a fitted model, and the stopping is
structural rather than a convention.

### /harmonize will answer 501 for as long as there is no model

Not an error path: a caller posting measurements gets back every measurement
exactly as sent, the provenance and licence terms of each, the maturity stamp, and
a statement of why there is no harmonized value. What they do not get is a number.
The 200 shape is fully specified in `contracts.HarmonizeResponse` so callers can
build against it, and no endpoint returns it.

The response model for the endpoint is `HarmonizationUnavailable`, NOT
`HarmonizeResponse`. That is deliberate: one shape has a harmonized estimate and
the other does not, and declaring the richer one would let a caller write code
against a field that is never populated.

### The confidence grades are rules and have never been calibrated

A grade is not a probability. It is not calibrated against anything, because
calibrating it would need exactly the held-out matched ions this repository does
not have. It is a set of checks, each of which names a specific reason a
correction might not apply to a particular ion.

Two of the five thresholds are POLICY rather than derived, and are marked as such
in the source and in the published rules: the ten per cent extrapolation margin,
and the half a per cent leverage limit. The second is at least anchored to the
published stepped-field DTIMS interlaboratory reproducibility, which is a real
number from a real paper, quoted in CONTEXT.md and not read from the paper here.

**The grading has never graded a real harmonized value**, because none exists. It
has been exercised against synthetic strata only.

### One of the five rules cannot answer for an ion that is not already in the corpus

Stated here as a scope limit on the product rather than only as a per-response field,
because a caller decides whether to trust this service BEFORE they send anything to it.

Four of the five rules are properties of the CALIBRATION GROUP - how many ions it holds, how
far the submitted value sits outside its range, whether a source behind it may be used, how
much its slope is levered by ions that do not transfer. Those answer for any ion, because
they are not about the ion.

`flagged as an outlier in its own stratum` is the exception and the only rule that is about
the submitted ion. It is a lookup against ions already measured on both platforms here, so
it cannot answer for an ion the corpus has not seen - and an ion the corpus has not seen is
the ordinary case for anyone using this service.

**For that ordinary case the rule is inert always, not occasionally.** The figure to carry is
100 per cent, not the 30 per cent that replaying this corpus against itself produces: that
replay consists entirely of ions already present, so it measures a population nobody using
this service belongs to. Stated this way round because 124 of 417 reads like a minor gap, and
the true figure for the intended use is different in kind - the scheme publishes five
substantive rules and a caller bringing new chemistry is graded by four of them.

Where it cannot run:

- it is named in the response's `not_checked`, with why;
- the grade is demoted one notch, and the demotion names it;
- it is never counted as passed.

The scope limit is served on `/confidence/rules` in each rule's `applies_to` field. It is
the only rule there whose `applies_to` is not "any ion", and a test asserts that it stays
the only one - a caveat carried by every rule marks nothing.

### No record grades `supported`, and the reason is corpus size rather than data quality

OPEN, and deliberately not closed by moving a threshold. Measured 20 September 2026 over
the whole seed corpus: **186 qualified, 216 weak, 15 unsupported, 0 supported**, of 417
corrected records.

`thinly_populated` demotes any stratum holding fewer than `TARGET_MATCHED_IONS` = 100
matched ions. The nine strata this model applies hold **23, 23, 23, 29, 31, 31, 41, 46 and
46**. So the population rule fires on 417 of 417 records - every one - and the top grade is
unreachable for reasons that have nothing to do with the ion being graded.

That matters for how a grade READS. `qualified` is not a middling result here; it is the
CEILING. And `weak`, which 216 records carry, is usually that ceiling minus one notch: 209
of the 216 are `qualified` demoted once because a rule could not be evaluated, not because
anything about the measurement is wrong. Only 7 of the 216 are weak for a reason about the
submitted value itself - it sat outside the range the correction was fitted over.

**The top grade is live, not dead code.** Built on a synthetic stratum of 100 matched ions
with one ion that does not transfer, so that leverage is measurable, `grade_correction`
returns `supported` for an ion in that stratum. What stands between a real caller and it is
the corpus, not the code:

| applied stratum | ion already in the corpus | ion the corpus has never seen |
|---|---|---|
| 23 to 46 ions, as today | `qualified` | `weak` |
| 100 ions or more | **`supported`** | `qualified` |

**What would lift the ceiling**, in order:

1. **An applied stratum reaching 100 matched ions.** That is the only change needed for
   `supported` to become reachable, and it is data rather than code. Today's largest applied
   stratum holds 46, so it needs 54 more ions measured on BOTH platforms of one pair, within
   one calibration group and one adduct - not 54 more measurements, which would spread
   across strata and lift none of them.
2. **Leverage measurable in that stratum**, which needs at least one flagged outlier and
   three points remaining once they are excluded. Where it is not, the response says so in
   `not_checked` and the grade falls a notch for that instead.
3. **For the intended user, nothing lifts it past `qualified`.** An ion this corpus has
   never measured cannot have the outlier rule evaluated for it, and an unevaluable rule
   costs one notch - so a genuinely new ion tops out at `qualified` at ANY corpus size under
   this scheme. That is a consequence of a rule keyed on corpus membership, recorded here
   rather than left for a caller to discover when their grade never improves.

**The threshold has NOT been lowered** and will not be, to make the top grade appear.
Adjusting a check until it passes is fitting the check to the data - the same objection this
file already records against touching the leverage limit. 100 is a policy figure, stated as
policy in `/confidence/rules`, and it is arguable; but it must be argued on what a
platform-pair figure needs before being quoted, not on what would make today's corpus look
better.

### Three mismatches between what we show and what we serve

OPEN, all three, found in the adversarial pass of 21 September 2026 and left as they are
because each is a disclosure question rather than a defect in a number.

**1. A demotion a caller can receive is not in the published scheme.** A response can carry
`rule: "far outside the calibration range"` - it does, on 2 of the 417 corrected seed records
- and `/confidence/rules` publishes six rules, none of them by that name. The nearest is
`outside the calibration range`, which is what it collapses to: one published rule that fires
at two severities under two different names. A caller who receives the severe one and goes to
look it up finds nothing.

This is the principle stated for the meta-rule in Group 1a - *a demotion a caller can receive
but cannot look up is one they cannot argue with* - left unapplied to this instance. The test
suite does not catch it because `PUBLISHED_BY_DEMOTION` in `tests/test_grading.py` maps the
two names deliberately, so the suite is satisfied by a mapping the caller cannot see.
`no correction can be fitted` has the same shape and collapses to `thinly populated
calibration group`, though no seed record produces it. *To close it:* publish both names, or
give the rule one name and carry the severity separately. Either is a contract change.

**2. A correction is computed and shown but never served.** `robust_derived_ccs` is on the
`Harmonized` object and `tools/demo_end_to_end.py` prints it beside the slope- and
median-derived values, so the demonstration shows THREE corrections where the API returns
two. `HarmonizedEstimate` deliberately carries only the slope- and median-derived values -
the contract docstring says why - so the demo is showing a number a caller cannot obtain.
The demo now labels it `<- computed, and NOT served by the API`, which is a disclosure rather
than a fix. *To close it:* serve it, or stop printing it. Serving it is the better answer if
anybody wants it, since it is already computed and already digested into `parameters_sha256`.

**3. The rounding mode is unstated.** Served values use Python's `round`, which is
round-half-to-even, so 0.125 rounds to 0.12 while 166.65 rounds to 166.7. The displacement is
at most half a unit in the last served place and the interval is at least ten units wide by
construction, so this cannot move a value outside its own interval or change any comparison
that matters - but nothing anywhere says which rounding is in use, and a reader who notices
the asymmetry has no way to find out that it is deliberate.

### Relabelling an ion still extracts a value the service refuses under its real name

OPEN, and it cannot be closed by grading. Measured on 20 September 2026, after the fix above.

Of the 15 seed records the service refuses outright under their own identity, **5 return a
value when the same measurement is resubmitted with a new compound name**. Worked example:
androstanedione's travelling-wave record is refused as itself - the corpus has it flagged as
an ion that does not transfer, and that rule falls straight to `unsupported`. Renamed to a
compound the corpus has never seen, the same numbers come back as 176.774 with an interval
of 173.358 to 180.054, graded `weak`.

The mechanism is not a bug in the fix; it is the fix working as ruled. The outlier rule is
the only rule keyed on the analyte's identity, identity is supplied by the caller and cannot
be verified here, and the ruling of 20 September 2026 is explicit that an unevaluable rule
demotes ONE NOTCH rather than refusing - because a genuinely new ion may be perfectly fine,
and refusing every new ion would refuse the platform's entire purpose. So a relabelled
outlier is indistinguishable from a new ion, and gets a new ion's treatment.

**What changed is what the response says**, and that is the whole of the improvement. Before,
the renamed record came back with `not_checked: []` - a positive claim that every rule had
been evaluated, including the one that would have refused it. Now it comes back a grade
lower, with the outlier check named as unevaluable and a demotion saying the grade rests on
fewer checks than the scheme advertises. The number is still served; the claim that it was
fully checked is not.

This is the identity analogue of a limit already recorded for licences in
`SourceProvenance`: a field the caller supplies cannot gate anything, because a caller
refused on its value edits it and resubmits. It obstructs only honest callers. The same
holds here, and the honest disclosure is that **this service cannot tell a new ion from a
renamed one, and does not claim to.**

Closing it needs corpus growth or structure-based identity, not a stricter rule: an ion
resubmitted under a new name but with the same InChIKey already matches, and it is only the
dataset-local compound name that can be freely rewritten.

**What would remove the limit** is not a code change. The rule can only answer for ions
measured on both platforms of a pair, so it becomes generally useful exactly as the corpus
grows toward covering the chemistry callers actually send. Lowering it to a distance check
against the stratum median would make it fire for everyone, but that is a different rule
answering a different question, and it would no longer mean "this ion is known not to
transfer".

## 7D. The first real data: the steroid interplatform study

DOI 10.1021/jasms.2c00196, supporting file SI_3, sheet `S2_Interplatform CCS
Database`. Converted by `tools/ingest_steroid.py`, which also writes the sheet out
verbatim to `data/seed/as_delivered/js2c00196_si_003_S2.csv` so the conversion is
diffable against its input from a clone with no `data/raw` directory.

**142 cross-platform matched ions.** That is the first result in this project that
describes real instruments, and it is what every remaining milestone was waiting
for. 521 records, every row of the sheet accounted for, nothing transcribed by
hand, no row lost. Over ALL 521 records, two ions pair across two platforms, 43 across
three and 97 across all four. Over the 517 that clear the gate - the only set that may
train, and the one every M4 figure comes from - it is 2, 47 and **93**: the four records
the shared-peak check holds each drop their ion from four platforms to three. This
document quoted only the first triple until 19 September 2026, which read as though 97
ions were trainable when 93 are.

### The coverage is uneven, and that is the real problem arriving early

| platform | values |
| --- | --- |
| TWIMS, cross-laboratory average (n=12, four instruments) | 142 |
| TIMS, timsTOF pro | 142 |
| DTIMS single-field | 135 |
| DTIMS stepped-field | 102 |

**Forty ions have no stepped-field value at all**, and stepped-field is the only
primary method here - the only one whose CCS is not read off a calibration curve.
So for 40 of 142 ions there is no first-principles anchor in this dataset, and any
harmonization of them rests on calibrated values alone. That is not a defect in
the data; it is the condition the platform exists to handle, and it showed up in
the very first source.

### The trapped-ion calibrant was unstated for one day, and the article resolved it

**CLOSED 19 September 2026.** The history is kept because the first reading was
right and was still not the end of it.

The supporting information names the travelling-wave calibrant (Waters Major Mix)
and the single-field drift-tube calibrant (Agilent ESI-L G1969-85000). For trapped
ion mobility it states only the MASS calibration - "10 mM sodium formate and a 7th
order high-performance calibration" - and then says: *"In addition to the external
calibration, each sample was automatically post-run calibrated by injecting a 1:1
mixture of both calibrants."*

**"Both calibrants", and only one named.** On the SI alone the honest record was
`UNSTATED_CALIBRANT`, and all 142 trapped-ion values were held. The obvious guess -
the Agilent mix, which the SI background calls typical for TIM-MS with an "e.g.",
and which this same laboratory used for its drift-tube work - was deliberately NOT
written in, because an inference recorded as a stated method is indistinguishable
downstream from a fact.

**The article states it.** The Experimental section, read from the Europe PMC open
full text for PMC9545150:

> "Prior to analysis, the instrument was mass calibrated with sodium formate
> clusters (10 mM in 50:50 2-propanol/water) and TIM CCS N2 was calibrated using
> ions from **Agilent ESI-L Tune Mix** via a linear function."

So "both calibrants" is sodium formate for mass and the Agilent mix for mobility.
The guess was right, and that it was right is not what makes the value usable - the
sentence is. **142 values moved from held to trainable and the corpus gained its
third technology:** 517 of 521 records now clear, six platform pairs instead of
three, 701 paired points instead of 326.

**The general lesson, which is why this is written at length: the supporting
information is not the source.** It is one document of the source, and it omitted a
method the article states plainly. Any future `UNSTATED_CALIBRANT` should be checked
against the article's own Experimental section before anybody emails an author.

It also settles which platforms share a calibrant, which is the paper's own central
finding: *"DT CCS N2 and TIM CCS N2 are routinely calibrated with the same
commercially available compound mixture (i.e., reference ions and reference values)
established by Stow et al., while TW CCS N2 systems were calibrated using a
different commercial calibrant mix."* This corpus reproduces that directly - the two
platforms sharing a calibrant agree to within a third of a per cent (bias -0.32%,
CCC 0.998), and the one that does not is three times further out (-1.07%). The
drift-tube single-field and trapped-ion records therefore carry the SAME calibrant
string here, on the paper's authority rather than because two spellings looked
alike.

One methodological difference from the paper, which is why our figures will not
match theirs exactly: they *"excluded"* ions whose residuals fell outside the
whiskers from their linear models. This repository reports outliers and never
removes them, and uses Deming rather than OLS.

### There is no pooled figure for any platform pair, and the adducts are why

Each platform pair splits into three calibration-group strata, one per adduct,
because the adduct is part of the calibration group. `statistics.py` refuses to
pool them. That refusal was built against synthetic data and it earns itself here:
the bias really does differ by adduct. Travelling wave against single-field DTIMS
is -1.07 per cent for [M+H]+, -0.14 per cent for [M+Na]+ and -0.20 per cent for
[M-H]-. A pooled number would have been an average of three different things, and
would have looked like a finding.

Negative mode is by some distance the worst case in the set: Lin's concordance
0.963 against stepped-field DTIMS, limits of agreement spanning roughly 16 Å², and
the largest lambda measured from reported uncertainties anywhere in this corpus
(459.8, meaning the travelling-wave standard deviations are far larger than the
drift tube's). Twenty-five ions is a small stratum, and nothing here separates "a
genuinely harder measurement" from "a smaller n".

### Five outliers, kept, and three compounds between them

The outlier rule flagged five ion-stratum pairs, which are three compounds
recurring across strata. The two worst are the largest and most flexible ions in
the set: an undecylenate ester (+2.3 per cent between platforms) and a
diglucuronide (+6.8 to +7.3 per cent). Every one stays in the data, in the
statistics and in the reported counts.

**A 7 per cent cross-platform disagreement on a real published ion is a result, not
noise.** It is precisely where harmonization is hard, and removing it would erase
the finding while improving every summary statistic. Whether it is a conformational
difference, a misassigned peak or something else is a question for somebody who
knows steroids; this repository's job was to surface it, and it did.

### The shared-peak check fired on real data, four times

Two different steroids - 4-androstene-17α-methyl-17β-ol-3-one and
4-estrene-17α-ethinyl-17β-ol-3-one - are reported with an identical CCS of 177.76
in one calibration group. Four records are held for review. Two isomers can
genuinely agree to two decimal places, so this is not an accusation of a
transcription fault in the paper. It is exactly what the check is for: hold it,
name it, let a person decide. The check was written against synthetic collisions
and had never seen a real one.

### What was NOT ingested, and why

The sheet carries a FIFTH platform column: a single-laboratory travelling-wave
database from one Synapt G2-S, alongside the cross-laboratory average of four
instruments. **It is not ingested.** Both columns are travelling wave in nitrogen
against the same calibrant, so they land in one calibration group on one platform,
which makes them a replicate pair - and matched-ion construction HOLDS an ion with
replicates on one side of a comparison rather than averaging them or picking one.
Ingesting both would have held all 142 ions and the dataset would have yielded
nothing.

That is a real limitation and not a tidy-up: 142 published values are in the
repository's reach and are not in it. Using them needs a decision about what a
same-platform replicate pair means for a cross-platform comparison, which is an M4
question and is not answered by a conversion script. The count and the reason are
printed by the adapter on every run rather than left in this file.

### The compounds are identified by this dataset's names and nothing else

The sheet gives a compound name, a commercial name, a formula and an m/z. Nothing
identifies a structure - no InChIKey, no SMILES. None was invented and none was
looked up: resolving 87 names against a structure database fails silently and
wrongly on exactly the compounds that matter here, which are isomers with similar
names, and a wrong structure assignment would be invisible downstream.

So every compound is recorded as `steroid_jasms2022:<name>`. Within this dataset
that is a real identity on the paper's own authority, which is what lets the 142
ions pair at all. Across datasets it matches NOTHING, which is correct rather than
unfortunate - and it means **this source cannot yet be pooled with CCSbase or the
Bush Lab data even after those are ingested.** Resolving these names to structures
is a curation act needing a provenance trail, and it is the prerequisite for any
cross-source pairing.

Two compound names carry an asterisk, which the sheet's own footnote explains as
multiple conformations observed by drift tube or trapped ion mobility. The asterisk
is KEPT in the recorded identity rather than stripped, because stripping it would
merge a starred compound with an unstarred one of the same name.

### The adapter reads columns by position, and nothing can mutation-test it

`tools/ingest_steroid.py` reads platform columns by index, so it checks that every
column it reads carries the header text it expects and refuses the sheet outright
otherwise. Without that, a revised supporting file with one column inserted would
be ingested silently, recording travelling-wave values as trapped-ion ones and
producing a bias figure that looked like a result.

That guard has a test. It does NOT have a mutation, and it cannot: the harness
shadows `src/wmxccs` and only mutates files inside the package, by design - that
design is what lets it guarantee it never writes the repository. So the adapter, and
everything else under `tools/`, is tested but not mutation-verified. Given section
4.5, that gap is named here rather than assumed harmless. Closing it means giving
the harness a second shadow root, which touches the invariant that makes it safe,
and was not worth doing inside this change.

## 7E. CCSbase and the Bush Lab database: retrieved, characterised, not ingested

Both were downloaded on 19 September 2026 and neither has been converted. What
follows is what they turn out to hold, recorded because in both cases it differs from
what CONTEXT.md assumed.

### CCSbase: 25,020 records, three platforms, and no stated drift gas

**There is no bulk download.** The site's `/download` button returns a 182-byte
`batch_query.csv` - the TEMPLATE for the batch-query upload feature, not an export.
The paper behind it (Ross, Cho & Xu 2020, Anal. Chem.) is not open access and has no
supplementary data in Europe PMC. The table is paginated at ten rows over 2,502 pages.

What worked is the site's own search form, which POSTs to `/results`: one broad query
returned all 25,020 rows in a single 19.5 MB response. That is the site used as
intended rather than crawled - one request, not 2,502 - and the terms restrict the
PURPOSE of use (academic, non-commercial) rather than the method. There is no
robots.txt. The response is kept in `data/raw/ccsbase/` with its sha256.

The table carries twelve columns, and they are better than expected:

| column | what it gives |
| --- | --- |
| ID | a stable per-record id, `CCSBASE_A4F2E9AA6E` |
| Name, Adduct, m/z, CCS, Z | the measurement |
| SMI | **SMILES**, so a structural identity is reachable |
| Type | compound class |
| Ref | a primary paper, linked to its DOI - 36 listed, 35 with records, blank on 3,572 |
| CCS Type | **the platform**: DT 10,967, TW 8,529, TIMS 5,524 |
| CCS method | 23 distinct methods, most naming a calibrant |

So CCSbase is genuinely platform-aware, as CONTEXT.md said, and it carries the
calibration method too. Three things stand between it and ingestion, and none is a
code problem:

1. **THE DRIFT GAS IS NOT A COLUMN.** Not for any of the 25,020 records. One method
   string mentions "helium and nitrogen drift gas" (120 records), which is proof that
   the database is not uniformly nitrogen and therefore that assuming nitrogen would
   be writing in a method detail. Under the keying rule of section 4.6 every one of
   these records is unmatchable until its gas is resolved, so **CCSbase currently
   contributes zero matched ions.** The gas is stated in the primary papers, so this is
   a bounded curation task rather than an open one, and each paper resolved unlocks its
   own records.

   **THE BOUND DOES NOT COVER EVERYTHING, and that is the part to carry.** The
   bibliography lists 36 papers; only **35 contribute any record** to this export (ref
   30, ToxBase, contributes none); and **3,572 of the 25,020 records - 14 per cent -
   carry no reference at all.** So the bounded task resolves at most 21,448 records, and
   a seventh of the database has no stated origin to resolve the gas from. Those records
   need CCSbase itself to say where they came from, which is a question for the
   maintainers rather than a curation task.
2. **8,388 records name no calibrant.** `CCS method` reads "single field, calibrated"
   for 5,233 DT and 2,950 TIMS records, and "?" for 205 more. Those take
   `UNSTATED_CALIBRANT` and are held - exactly as the steroid trapped-ion values were until
   the article resolved them, which is the first thing to try here too.
3. **PROVENANCE IS TWO-LAYERED.** The values belong to the 36 listed primary papers,
   where one is named at all; CCSbase is the compilation. Its terms govern the compilation and were the terms read on
   11 September 2026, so a record's `source` should be CCSbase - that is where the
   value was obtained and whose terms permit the use - with the primary paper's
   reference and DOI recorded in `source_locator`. That keeps the licence claim
   backed by the entry that actually covers it, and keeps the real origin visible.
   It does NOT establish that each of the 36 papers permits reuse, and nothing here
   should be read as claiming it does.

Also newly read, verbatim, from the site's own terms and not previously recorded:
*"Any derivative works (e.g. softwares, websites) must reproduce the above copyright
notice"*, and use *"must be for academic, non-commercial purposes only"*, with
commercial users directed to Dr Libin Xu and UW CoMotion. The notice-reproduction
clause is an obligation this repository would take on by ingesting, and it is
stronger than the plain citation requirement the registry recorded.

### The Bush Lab database: the only real biopharmaceutical data identified

Downloaded from the documented published-to-web xlsx URL, which needs no
workaround. 213 KB, eight sheets, **1,804 data rows**:

| sheet | header rows | data rows |
| --- | --- | --- |
| Native-Like Protein Cations | 1 | 72 |
| Native-Like Protein Cations and (complexes) | 2 | 164 |
| Denatured Protein Cations | 1 | 27 |
| Polyalanine Cations | 2 | 31 |
| Anionic Homopolymers | 2 | 35 |
| Other Peptides | 1 | 12 |
| Small Molecular Ions | 1 | 23 |
| MicroSource Collection | 1 | 1,440 |
| **total** | **11** | **1,804** |

CORRECTED 23 September 2026. This table previously read "roughly 6,500 rows" and gave
1,000 for four of the eight sheets. Four figures were an artefact rather than a count -
a reader default or a pre-allocated range - and the true total is 1,804, so the only
biopharmaceutical dataset this project has identified was overstated by a factor of 3.6.
The four small figures in the old table were accurate but counted the header, which is
the tell: real counts and artefacts sat side by side.

**Count data rows, not non-empty rows, and mind the two-row headers.** Three sheets put
a merged banner above the column names - Cations/Anions on the complexes sheet,
`z = 1 | z = 2 | z = 3` on Polyalanine, and a grouping row on Anionic Homopolymers - so
non-empty rows total 1,815 and data rows total 1,804. A count that treats every sheet as
having one header row gives 1,805 and is wrong by one, on Polyalanine.

This is the data the M1 biopharmaceutical identity layer was built for and has never
seen: native-like and denatured protein cations, protein complexes, peptides. Five
things are visible from opening the sheets, and each one changes what a converter has
to do.

**1. The gas IS stated, in the column headers.** `Ω(He) / nm^2` and `Ω(N2) / nm^2`,
exactly as the steroid sheet states its gas in `TWCCSN2`. So unlike CCSbase, these
records are MATCHABLE the moment they are converted - the blocker of section 4.6 does
not apply. Most ions carry both a helium and a nitrogen value, which become two
records on two keys rather than one: a helium value and a nitrogen value are not
values of the same ion, which is the whole reason the gas is in the key.

**2. THE UNITS DIFFER BETWEEN SHEETS, by a factor of 100.** The protein and peptide
sheets report nm², the polyalanine, homopolymer, small-molecule and MicroSource sheets
report Å². Getting this wrong would put a protein cross section out by two orders of
magnitude. It is the single most dangerous detail in the file and the easiest to miss,
because both spellings look like a unit and neither looks like a mistake.

**3. THE PROTEIN SHEETS CARRY NO ADDUCT, only a charge.** A row says z = 3 and never
says what the three charges are. This is precisely the case CONTEXT.md predicted would
be common rather than rare in this source, and the `[M+24?]24+` form was built in M0
for it. **Prediction confirmed.** Every such record keys uniquely to itself and can
never pair, so the native and denatured protein data - roughly 3,000 rows - arrives
unmatchable unless the carrier can be established from the cited papers.

**4. The MicroSource Collection is the clean part.** 1,440 CCS VALUES over 1,424 distinct
COMPOUND NAMES - see the figures below for which number counts what - with a
stated adduct (`[M+H]+`), charge, nitrogen CCS, a per-row standard deviation in
`s / Å^2`, a formula, a CAS number and a reference. Adduct, charge, gas and
uncertainty all stated: this is the most immediately usable table found in any source
so far, including the steroid study.

**5. 41 cells in one sheet are corrupt**, holding date serials far outside any valid
range (21955915 in C99, 9677161 in C108, and 39 more) which openpyxl refuses to read.
A converter must count and report them, never let them arrive as nulls.

One discrepancy to settle before ingesting: the registry records these values as
DTIMS, and the database page describes the ions as *"primarily from traveling-wave ion
mobility spectrometry"*. Those are different platforms and the difference is the
entire subject of this repository. The per-row `Ref` column names the paper for each
value, so this is answerable per row rather than per file, and it must be answered
rather than assumed.

Its terms ask only that users cite the appropriate publications, which is the
lightest obligation of the three sources, and the registry records it as
`academic_only` on the conservative reading.

### Bush Lab: converted, and every value blocked on a platform the file never states

OPEN as of 26 September 2026, and the blocker is one column.

`tools/ingest_bushlab.py` reads all eight sheets, resolves every structural trap, and writes
the flat conversion to `data/seed/as_delivered/bushlab_ccs_database_converted.csv` - one row
per CCS value, carrying the native unit, the native value, the conversion factor and the
converted value, so the arithmetic is checkable rather than trusted. **1,804 data rows yield
2,045 CCS values**, because a row can hold a helium and a nitrogen column and the polymer
sheets hold one column pair per charge state.

**IT WRITES NO SEED FILE, and that is the result rather than a gap.** The workbook states NO
instrument, method or platform anywhere - not on any sheet, not in any cell, not in any
header. It states the gas, the charge, the cross section and a per-row `Ref` naming one of
nine papers. A `CCSMeasurement` requires an `ims_type`, `IMSType` has no UNSTATED member, and
DTIMS additionally requires `dtims_method`. So not one of the 2,045 values can become a
record until the nine papers are read. The database page's "primarily traveling-wave" is a
sentence about the collection, and `primarily` is not a value a row can carry.

| blocker | values | what resolves it |
| --- | --- | --- |
| no platform stated | **2,045** (all) | the nine cited papers, each stating its own |
| no charge stated | 46 | `Small Molecular Ions` has no charge column at all |
| adduct the schema refuses | 3 | `M+` twice, `[M+3H]+3` once; NOT repaired, see below |

The three adducts are left exactly as written. `[M+3H]+3` is almost certainly `[M+3H]3+`, and
that is the reason it is not corrected: "almost certainly" is a guess about somebody else's
data, in the column that decides which ions match each other.

**MATCHED IONS UNDER THE STRICT RULE: ZERO, AND ZERO EVEN ONCE THE PLATFORMS ARE KNOWN.**
This is the finding worth carrying, and it is not the platform blocker. Of the 2,045 values,
559 carry an UNSTATED CHARGE CARRIER - the protein, peptide and polymer sheets give a charge
and never say what the charges are - and such a key is unique per record by design, so those
values can never match anything, including each other. 46 more have no charge and 3 have an
unusable adduct. That leaves **1,437 values that could ever participate, and every one of
them is on the MicroSource sheet, from a single paper**: 1,425 distinct ion keys, **none of them
appearing in more than one cited paper**. One paper is one platform, and a matched ion needs
two. So the ceiling is zero, and resolving the nine papers does not raise it.

**CROSS-GAS PAIRS: 165**, over 347 values, being one analyte at one charge measured in both
helium and nitrogen. 146 of the 165 have both gases on one sheet - usually the two columns of
a single row. By sheet: Native-Like Protein Cations 96, Anionic Homopolymers 88, Polyalanine
86, Denatured Protein Cations 36, Other Peptides 22, complexes 19.

### Hines 2017 IS READ, and it resolves the platform for 70 per cent of the source

READ 26 September 2026, open access, and the one paper worth reading of the nine. DOI
10.1021/acs.analchem.7b01709, PMC5616088. Hines, Ross, Davidson, **Bush** and Xu,
"Large-Scale Structural Characterization of Drug and Drug-Like Compounds by High-Throughput
Ion Mobility-Mass Spectrometry", Anal. Chem. 89, 9023-9030 (2017). Read in full text from
Europe PMC; the ACS page refuses automated access and PubMed answers with a CAPTCHA.

It governs the MicroSource sheet, which is 1,440 of the 2,045 CCS values in this source - 70
per cent - and the only sheet whose values could ever participate in a matched ion.

| what the ingest needed | what the paper states, quoted |
| --- | --- |
| platform | **TWIMS**: "IM-MS analysis was performed on a Waters Synapt G2-Si HDMS" |
| drift gas | **nitrogen**: "using nitrogen as the drift gas" |
| calibrated? | **yes**: "Using a combination of small molecule and polypeptide CCS calibrants" |
| charge carrier | **stated**: "masses corresponding to the protonated, sodiated, and water-loss ions were extracted for each data file" |
| values / compounds | "A total of 1440 CCS values representing 1425 unique compounds (71% coverage) were obtained" |

So the sheet's own adduct column is the paper's carriers, and the three adducts it holds in
bulk - 1,279 [M+H]+, 89 [M+Na]+, 40 [M+H-H2O]+ - are exactly protonated, sodiated and
water-loss. THE CARRIER IS STATED FOR THESE RECORDS, which is what separates MicroSource from
the 559 protein, peptide and polymer values that give a charge and never say what it is.

**THE CALIBRATION IS THE FINDING, AND IT IS A CIRCULARITY RISK.** The calibrants are PolyAla,
acetaminophen, betaine hydrochloride, alprenolol hydrochloride, clozapine N-oxide,
erythromycin, ondansetron hydrochloride, reserpine, vancomycin hydrochloride, verapamil
hydrochloride, and a peptide Ac-ETDYYRKG-NH2. Their reference values are **the authors' own
nitrogen DTIMS measurements on a modified Waters Synapt G2 HDMS** - an RF-confining drift
cell in the same laboratory.

That matters to this repository specifically. These TWIMS values are DERIVED from Bush Lab
DTIMS values, so a comparison of Bush Lab TWIMS against Bush Lab DTIMS would be partly
circular: it would be measuring how well a calibration reproduces the thing it was calibrated
against. The `calref_reference_set`, `calref_doi`, `calref_platform` and `calref_method`
fields exist to make exactly this visible, and they must be filled for these records rather
than left blank. It does not block ingestion and it does block one comparison.

**Coverage is 71 per cent of the plate, not of the drug universe.** "The CCS values for the
rest of the 560 drugs ... were not successfully determined due to either low peak intensity
(<1 x 10^3 counts) or the peak being too wide (>=25 bins)." So roughly 1,985 compounds were
attempted. The absent 560 are absent for measurement reasons, which is a selection this
corpus inherits and should not describe as coverage of anything wider.

**What it does NOT change: the matched-ion count is still zero.** MicroSource is one paper on
one platform, and a matched ion needs two platforms. Resolving Hines 2017 turns 1,440 blocked
values into 1,440 records that can be BUILT and held as a reference library; it creates no
pair. That was established before the paper was read and the paper does not alter it.

### MicroSource is INGESTED as a reference library: 1,437 records, none trainable

INGESTED 26 September 2026, by `tools/ingest_bushlab.py --platform-map
data/seed/bushlab_platforms.json`, into `data/seed/bushlab_microsource.csv`.

| | |
| --- | --- |
| values offered by the sheet | 1,440 |
| records BUILT | **1,437** (the 3 unusable adducts are not built) |
| records CLEARED to train | **0** |
| records HELD | **1,437** |
| matched ions contributed | **0**, unchanged at 142 for the corpus |
| model fingerprint | `064eb9fba603/0d69f799f6f1`, UNMOVED |

Zero cleared is the correct outcome and not a failure of the ingest. These are a reference
library by instruction and by structure: MicroSource is one paper on one platform, so it can
form no matched ion, so nothing was ever going to enter a fit. The gate refuses for four
reasons and each is worth reading.

**552 refused because permission does not reach the publication.** The row claims
`academic_only`, and the gate requires a published source to name a DOI with a licence record.
The licence that was read is the DATABASE PAGE's citation request; the ACS paper's licence has
not been read - Europe PMC flags it open access and states no licence, and an open-access flag
is not a licence. So `doi` is left null, following the decision already taken for CCSbase, and
the paper is named in `source_locator` instead. The gate then refuses because a null DOI cannot
back a published claim, which is correct: nobody has read terms that cover these values as a
publication. *To clear them:* read the licence on the ACS page and register it, or get the
written grant from the laboratory that this entry has recommended since 19 September.

**868 refused as suspected shared peaks, and THIS HEURISTIC DOES NOT TRANSFER.** The detector
refuses a record whose CCS is identical to another analyte's in the same calibration group. On
the steroid study that caught 4 rows in 521 and meant something: an unresolved peak shared
between two compounds. Here all 1,437 records sit in ONE calibration group and the sheet
reports **one decimal place** across a range of 108.8 to 355.8 square angstroms - **2,470
available 0.1 slots for 1,437 records**. Collisions are arithmetically inevitable: 833 distinct
CCS values, 375 of them shared, 979 records sharing one. Eight different compounds report
159.4.

That is the pigeonhole principle, not evidence of a shared peak. The heuristic is sound on a
few hundred records across several calibration groups and produces mostly false positives on a
single-platform library of this size and precision. **It has NOT been changed**: it is doing
what it was written to do, and loosening a detector because a new source trips it is how a
guard gets quietly weakened. What it needs is a decision - whether a reference library should
be gated by a rule written for a harmonization corpus - and that is a decision, not a code
change. Recorded here so the 868 is not read as 868 suspicious measurements.

**12 held for curation review**, being the six conformer pairs the paper does not name - see
below. **5 held as incomplete conformer sets**, which is the holds cascading: when one member
of a pair is refused for a shared peak, its sibling claims a sibling that is no longer in the
cleared set, and the loader holds it rather than letting it stand as though it were the only
peak. That interaction is correct and worth knowing about.

### The calibration lineage is recorded, and a comparison against it is now REFUSED

Condition 1 of the ingest, and the reason for the rest. Every one of the 1,437 records carries:

    calref_reference_set  the authors' own nitrogen DTIMS values for PolyAla n=2-21 and nine
                          drug standards, measured for this paper on a modified Waters Synapt
                          G2 HDMS with an RF-confining drift cell
    calref_doi            10.1021/acs.analchem.7b01709
    calref_platform       DTIMS
    calref_method         stepped_field

**Recording it was not enough, so the consequence is now encoded.**
`statistics.circularity_between` refuses any pair where one member's calibration reference
traces to the other member's publication and platform, and the refusal happens BEFORE a
stratum exists, so a circular pair cannot enter a fit even as one point among many. It is
counted on `ComparisonReport.pairs_refused_as_circular`, because a refusal nobody can count
cannot be told apart from an absence of data.

`CalibrationReference` has carried this lineage since M0 and its docstring named the risk
exactly - and until 26 September 2026 NOTHING READ IT at comparison time. The circularity was
detectable and undetected, which is the LIMITATIONS 4.5 class again. Ten tests in
`tests/test_circularity.py` pin it, including the case the guard exists for: a synthetic
stepped-field DTIMS record carrying the same lineage is refused, and an otherwise identical
record from a DIFFERENT publication is not. Four mutations, all killed.

### THE CIRCULARITY GUARD DOES NOT CATCH A SAME-LABORATORY CHAIN ACROSS TWO PAPERS

OPEN, and the most important thing to know before adding any more Bush Lab data. Read this
before ingesting Bush 2010, Bush 2012, Allen 2012, Allen 2013, Allen 2016, Salbo 2012,
Campuzano 2012 or Forsythe 2015.

**What the guard catches.** A pair where one member's `calibration_reference` names the
publication AND platform of the other member. That is provable circularity from the records
themselves: these MicroSource TWIMS values were calibrated against nitrogen DTIMS values
measured in the same paper, so a DTIMS record from `doi:10.1021/acs.analchem.7b01709` is
refused against them, and the refusal is counted.

**What it does not catch, and will not warn you about.** A chain that runs through the same
laboratory and the same instrument but a DIFFERENT publication. The Hines 2017 calibration used
an RF-confining drift cell on a modified Waters Synapt G2 at the University of Washington. Bush
2010, Bush 2012, Allen 2012, Allen 2013 and Allen 2016 are that laboratory's own DTIMS work on
that class of instrument. **If any of them is ingested, the guard will compare it against these
TWIMS records without objection**, because no field in either record says the two are related:
the DOIs differ, and the schema has no laboratory or instrument-identity field to match on.

**Why it was built this narrow rather than broader.** Whether Bush 2010's values are
independent of a calibration built in 2017 is a JUDGEMENT ABOUT TWO PAPERS - did the later
calibration use the earlier values, or re-measure the same standards? - and it is not a fact
recorded in either record. A guard that refused every same-laboratory pair would refuse real
comparisons, and one that claimed to catch this and did not would be worse than one that says
plainly that it does not. Widening it on a guess is the failure class LIMITATIONS 4.5 collects:
a guard that reads broader than it is.

*To close it,* in ascending order of effort:

1. **Read the papers before pairing them.** Each of the five states whether it measured its own
   standards or took published values. That answer, recorded as the `calref_doi` of whichever
   record is the calibrated one, makes the existing guard fire with no code change. This is the
   cheap answer and it is the right one.
2. **Give the schema a laboratory or instrument identity** and match on it. That is a real
   field on `CCSMeasurement` and a real decision about what counts as one laboratory across
   twelve years and an instrument modification, which is why it is not done here.
3. **Refuse same-`source` cross-platform pairs by default** and require an explicit
   declaration of independence to compare them. Blunt, safe, and it would refuse a legitimate
   comparison the day this corpus holds one - so it needs the CEO's decision, not a commit.

Until one of those happens, the honest statement is: **this platform can prove circularity
inside one publication and cannot see it across two papers from one laboratory.** The
MicroSource records are safe against the specific paper they were calibrated from, and unguarded
against the rest of the Bush Lab corpus.

### The sixteen doubled compounds: twelve conformers, four separate ions, six flagged

Condition 3. Hines 2017: the sixteen "display two peaks, had two major adducts, or were
mixtures for which we have reported a CCS value of each component" - three situations the sheet
does not separate. Reading the adducts separates them part of the way.

- **4 differ by adduct** and are ordinary separate ions, not conformers: glycocholic acid and
  methyldopa ([M+H]+ and [M+Na]+), methoxamine hydrochloride and podofilox ([M+H]+ and
  [M+H-H2O]+).
- **12 share an adduct**, so their two values are two peaks of one ion or two components of a
  mixture. All 24 rows are kept, numbered `conformer` 1 and 2 of `conformers_total` 2, ordered
  by cross section so the numbering is reproducible. Neither averaged nor dropped.
- **6 of those 12 are flagged for curation**, because the paper NAMES only the fluoroquinolone
  protomers (ciprofloxacin, norfloxacin, enoxacin, pefloxacine mesylate, sarafloxacin) and its
  new cephalosporin finding (cefpodoxime proxetil). The other six - antimycin A, bacampicillin,
  irigenin 7-benzyl ether, methimazole, montelukast sodium, temefos - are recorded as
  conformers AND carry a flag saying the sheet does not distinguish a second peak from a
  mixture component. **A mixture's two components are two analytes**, and calling them
  conformers of one ion would merge two compounds. Antimycin A is the clearest worry: the sheet
  labels it "antimycin a (a1 shown)", and antimycin A is a mixture of A1 to A4. Resolving these
  six means reading the paper's figures.

### The selection this corpus inherits: 1,440 of about 1,985 attempted

Condition 4. Hines 2017: "The CCS values for the rest of the 560 drugs ... were not
successfully determined due to either low peak intensity (<1 x 10^3 counts) or the peak being
too wide (>=25 bins)." So roughly 1,985 compounds were attempted and 1,425 succeeded, a stated
71 per cent.

**The 560 absences are not random.** They are compounds that ionised poorly or gave broad
arrival-time distributions, and broad peaks are what conformationally flexible or multi-protomer
species give. So this library is biased toward compounds that behave well in TWIMS, and any
statement about coverage of drug chemistry inherits that. It is a reference library of what
could be measured, not of what exists.

### The other eight papers are unread, deliberately, and none can change the matched-ion count

The economics are lopsided and were checked before choosing. Hines 2017 alone governs 1,440 of
2,045 values. The remaining eight govern 605, and **559 of those 605 are the unstated-carrier
values** - protein, peptide and polymer rows that give a charge and never say what the charges
are. An unstated carrier keys uniquely per record by design, so those values can never match
anything whatever platform they turn out to have been measured on. Reading their papers would
resolve `ims_type` and would not make one pair.

| unread paper | values it governs | what reading it would add |
| --- | --- | --- |
| Bush 2010 | 174 | reference-library records only |
| Allen 2013 | 102 | reference-library records only |
| Allen 2016 | 92 | reference-library records only |
| Bush 2012 | 86 | reference-library records only |
| Allen 2012 | 76 | reference-library records only |
| Campuzano 2012 | 46 | nothing: that sheet also states no charge |
| Forsythe 2015 | 16 | reference-library records only |
| Salbo 2012 | 13 | reference-library records only |

Campuzano 2012 is worth its own line: it governs `Small Molecular Ions`, which has no charge
column at all, so its 46 values stay unbuildable even with a platform. It is the one paper of
the eight that would resolve nothing on its own.

**The matched-ion count for this source is structurally zero and no reading changes it.** Not
zero-for-now, and not zero-pending-curation: zero because of what the file contains. 559
values carry a key that cannot match by construction, 46 have no charge, 3 have an adduct the
schema refuses, and the 1,437 that remain are one paper on one platform. Reading all nine
papers would produce a reference library of up to 2,042 records and zero matched ions. That is
the honest description of this source and it should be what the coverage matrix says.

### The three MicroSource figures, and which counts what

Stated once here because they differ by one and by sixteen, and a document that uses the
wrong one is wrong in a way nobody notices.

| figure | what it counts | where it comes from |
| --- | --- | --- |
| **1,440** | CCS VALUES, one per row | counted from the sheet, and the paper's own figure |
| **1,424** | DISTINCT COMPOUND NAMES in the sheet | counted from the sheet |
| **1,425** | the paper's stated count of UNIQUE COMPOUNDS | Hines 2017, abstract and results |
| 1,425 | also, coincidentally, the distinct (name, adduct, charge, gas) ION KEYS | counted from the sheet |

**1,425 is not a count of anything we hold.** It is the figure Hines 2017 reports - "A total
of 1440 CCS values representing 1425 unique compounds (71% coverage) were obtained" - and the
sheet does not reproduce it: 1,440 values over exactly 16 doubled names gives 1,424, which is
what the file contains. One of the two is off by one and this repository cannot settle which,
so both are recorded and neither is presented as the compound count without saying whose it
is. That the distinct ION-KEY count also lands on 1,425 is a coincidence and must not be used
to reconcile them.

CORRECTED 26 September 2026. An earlier report of this to the owner said 1,425 was the
ion-key count and therefore the origin of the figure in the coverage matrix. That was wrong:
it is the paper's compound count, and the key count matching it is chance.

### What the sixteen doubled compounds are, and why they are not one ion each

Hines 2017 says it plainly: "16 of these drugs display two peaks, had two major adducts, or
were mixtures for which we have reported a CCS value of each component". Three different
situations, and the sheet does not say which applies to which compound. Reading the adducts
separates them partly:

- **12 share one adduct**, so their two values are two PEAKS or a mixture, not two ions:
  ciprofloxacin, norfloxacin, enoxacin, pefloxacine mesylate and sarafloxacin - all
  fluoroquinolones, matching the paper's "known fluoroquinolone protomers" - plus cefpodoxime
  proxetil, which is the paper's "new finding of cephalosporin protomers", and antimycin A,
  bacampicillin, irigenin 7-benzyl ether, methimazole, montelukast sodium and temefos.
- **4 differ by adduct** and are therefore legitimately different ions: glycocholic acid and
  methyldopa ([M+H]+ and [M+Na]+), methoxamine and podofilox ([M+H]+ and [M+H-H2O]+).

The twelve matter more than their number suggests. A protomer is a CONFORMER of one ion, not a
second measurement of it, and the schema carries `conformer` and `conformers_total` for
exactly this. Ingesting the pair as one ion would average two structures; ingesting it as two
identical keys would make one ion look like a disagreeing duplicate. Which of the twelve are
protomers and which are mixtures is in the paper's figures, not in the sheet.

### Cross-gas comparison is a different kind of comparison, and is NOT a matched ion

Asked and answered rather than assumed, 26 September 2026. **Under our own rule it is not a
matched ion, and the rule is right.** The gas is part of the matched-ion key because a helium
cross section and a nitrogen cross section of one ion are not two measurements of one
quantity - they are measurements against different collision partners, and the nitrogen value
is larger by a margin that varies with the ion. Pooling them would be averaging two different
physical quantities, which is the failure the key exists to prevent.

So Bush Lab's headline strength - the same ions in both gases - is real and is NOT
harmonization data. It is **gas-dependence** data: same laboratory, same instrument, same ion,
two collision gases. Our comparison is the opposite shape: two platforms, one gas.

*What building it would take,* if it is wanted, and it is NOT built:

1. **Its own pairing unit.** A cross-gas pair keys on everything EXCEPT gas - analyte, adduct,
   charge, state - which is a different key from the matched-ion key and must not be confused
   with it or reuse its code path.
2. **Its own claim, and a guard.** `Claim` would need a member such as
   `GAS_DEPENDENCE_WITHIN_A_STUDY`, and `assert_may_be_quoted_as` must refuse to let a
   gas-dependence figure be quoted as a platform difference or as reproducibility, exactly as
   it refuses the interlaboratory claim today.
3. **Its own statistics.** The helium-to-nitrogen relationship is not a bias to be corrected
   out; it is a conversion to be fitted, and it is strongly size-dependent, so a single ratio
   would be wrong at both ends. A regression, not an offset.
4. **It must never feed the harmonization model.** A gas relationship is not a platform
   correction. Letting one into a stratum is precisely the forcing-of-incompatible-conditions
   the brief forbids.
5. **AND IT HAS THE SAME BLOCKER.** A gas relationship cannot be fitted without knowing
   whether both values came from one instrument, which is the platform question again. The
   nine papers gate this too.

### The bridge between sources is compound identity, and it is 30 compounds wide

Worth stating on its own, because it is the highest-value thing found while
characterising these two sources and it is not obvious from either.

Every matched ion in this repository is INSIDE one study. Bias measured that way is
bias between that study's own instruments, which is not the quantity this platform
exists to report. Pairing ACROSS sources would give bias between laboratories, and
what stands in the way is not the platforms or the gases - it is that nobody has
established which compound is which.

Measured, not estimated:

- **2 of 87** steroid compounds share their SYSTEMATIC name with a CCSbase compound.
- **30 of 87** share their COMMERCIAL name. Betamethasone, dexamethasone, cortisone,
  danazol, androstenedione, 17-hydroxyprogesterone and 24 more.

So a bridge exists and it is 30 compounds wide. **It must not be crossed on the
names.** The trap is visible in that very list: "androstenedione" names both
4-androstene-3,17-dione and 5-androstene-3,17-dione, which are different molecules
with different cross sections, and betamethasone and dexamethasone are C16 epimers of
one another. A name-based join would pair some of these correctly and some
incorrectly, and the incorrect ones would appear as inter-laboratory bias of a few per
cent - indistinguishable from the real thing, which is the entire measurement.

### The evidence for that rule, measured

`reports/steroid_name_resolution.md` is the attempt, made 20 September 2026 against
PubChem PUG REST and NIH CIR, to resolve all 87 names by machine. **It is the evidence
anybody questioning this rule should read first.** Its headline:

> 87 names, 6 resolved clean, 67 ambiguous, 14 failed - manual chemist review needed for 81.

**67 of 87 resolved only through the COMMERCIAL name, never the systematic one.** These
names are not machine-readable as published, and the commercial name is the weaker
identifier: eight of them name two compounds at once with a slash.

The sharpest case is the pair the demonstration happens to use:

| name | result |
| --- | --- |
| `4,9,11-estratiene-17β-ol-3-one` | resolved, **commercial name only** ("trenbolone") |
| `4,9,11-estratiene-17α-ol-3-one` | **FAILED** - its commercial name, "α-trenbolone/epitrenbolone", names two compounds |

Neither systematic name resolved: both carry the malformed "estratiene" parent described
below, and it defeats both resolvers. The β form was rescued only by its commercial name;
the α form was not rescued at all. **A pipeline that accepted the first and skipped the
second would keep one epimer and silently drop its partner** - and the two have different
cross sections, so the survivor would then be paired against whatever else carried that
name.

Two of the six clean resolutions are marked NEEDS CONFIRMATION in the report for the
converse reason: they resolved correctly AND share connectivity with another name here,
so a wrong answer in those two would have looked exactly as right as the correct one.

Nothing from that report has entered the corpus. It is a report, written outside `data/`,
and it exists to size the manual job rather than to do it.

The resolution is asymmetric between the sources, and that asymmetry decides the order
of work:

- **CCSbase can be resolved by machine.** All 25,020 records carry SMILES, with no
  exceptions. SMILES to InChIKey is a deterministic computation, so every CCSbase
  compound has a structural identity available.
- **The Bush Lab MicroSource sheet can be resolved too**, carrying a molecular formula
  and a CAS number for each of its 1,440 rows, which are 1,424 distinct compounds.
- **The steroid study cannot be resolved at all without a person.** It carries no
  structure, no SMILES, no CAS, no InChIKey. Its systematic names do fully specify the
  structures to a reader who knows steroid nomenclature -
  "4-androstene-17alpha-methyl-17beta-ol-3-one" is a complete description - so this is
  a bounded task of 87 compounds for somebody competent to do it, not an impossible
  one. It is simply not a task code can do silently.

#### A worked example, from the ion the demonstration happens to use

The point above is easy to read as fastidiousness. It is not, and one row of the
published sheet shows why. Verbatim, as the authors wrote it:

| compound name | commercial name | formula | ion |
| --- | --- | --- | --- |
| `4,9,11-estratiene-17β-ol-3-one` | trenbolone | C18H22O2 | [M+H]+ |
| `4,9,11-estratiene-17α-ol-3-one` | α-trenbolone/epitrenbolone | C18H22O2 | [M+H]+ |

**The first name is malformed.** Three locants - 4, 9 and 11 - describe three double
bonds, so the parent must be a TRIENE and the name required is
*estra-4,9,11-trien-17β-ol-3-one*. "Estratiene" names a diene. The compound is
trenbolone, which the sheet's own commercial-name column says and the systematic name
does not quite.

**And the row beneath it is the same molecule at one stereocentre.** Same formula, same
m/z of 271.1693, differing only in `17α` against `17β`. They are epimers with different
cross sections, and this repository holds both.

So a machine resolution of these names has to do three things at once that nothing does
reliably: recognise a malformed parent and repair it rather than fail; refuse to fall back
on the formula or the mass, which are identical for the pair; and carry the stereo
descriptor through, which is the ONLY thing distinguishing them. A lookup that gets
"trenbolone" from the commercial-name column and stops has silently merged an epimer pair.
A lookup that fails on the malformed name and falls back to C18H22O2 has done the same.

This is also the reason the corpus is NOT corrected here. The name is the source's, it is
what the supporting information published, and rewriting a published identifier inside the
data would create a record the source does not contain. The malformation is recorded; the
data is left as it was found.

None of this is implemented, and the dataset-scoped identity in `identity.py` is what
keeps it honest in the meantime: these compounds pair inside their own source and
nowhere else, which is the correct answer until the work is done rather than a
limitation to route around.

## 7F. M4: the harmonization model

Built 19 September 2026 on the 93 four-platform matched ions that clear the gate, from
one study. `robust.py`, `scope.py` and `harmonization.py`.

### What it produces

Eighteen strata - six platform pairs by three adducts, never pooled - each fitted three
ways. **Only nine are ever applied**: the ones whose reference is stepped-field DTIMS, the
only primary method in the corpus. A correction between two CALIBRATED platforms anchors
nothing, since both sides already rest on somebody else's reference values, so those nine
are published as diagnostics and never applied to a measurement. There is no transitive
composition: this package will not refer TWIMS to TIMS and then TIMS to the drift tube.

Of 517 cleared records, **417 are corrected and 402 of those are served a value**; the 15
between the two grade `unsupported` and are withheld with their reason. 100 records are
already on the primary platform,
which is not a failure but the correct answer for a value with nothing to be referred to.

### The headline correction is chosen by a rule, not by preference

    basis = ROBUST_SLOPE   if the Passing-Bablok rank interval on the slope EXCLUDES 1
            MEDIAN         otherwise

Measured: ROBUST_SLOPE in 10 of 18 strata, MEDIAN in 8. Where the slope cannot be
distinguished from 1, a slope-derived correction and a constant offset describe the same
data and the offset is the one that does not extrapolate. All three corrections - Deming,
Passing-Bablok and the median offset - are computed and reported for every stratum
whichever is chosen, because they diverge exactly when a few ions are levering the fit.

**Why Passing-Bablok, in one measurement.** On a clean line of eleven points, wrecking one
point moves the Deming slope from 1.04 to over 6 and does not move the robust slope at
all. On the real corpus the same effect is concentrated in negative mode:

| stratum | Deming | Passing-Bablok |
| --- | --- | --- |
| DTIMS/stepped vs TWIMS [M-H]- | 1.19024 | 1.05052 |
| DTIMS/stepped vs TIMS [M-H]- | 1.14698 | 1.00160 |
| DTIMS/single vs TIMS [M-H]- | 1.10959 | 1.01609 |

Those Deming slopes carry intercepts of -31 to -40 square angstrom to compensate. A slope
of 1.147 between two platforms whose median difference is -0.6 per cent is two leveraged
ions and an intercept, and Lin's concordance above 0.97 in the same strata flags none of it.

### The interval, and why there is no train/calibrate split

Section 4.7 establishes that a 90 per cent split-conformal interval is the FULL OBSERVED
RANGE for any calibration set below 19. The strata hold 23, 29 and 41 ions, so a split
would either starve the fit or produce an interval that excludes nothing.

Grouped leave-one-compound-out **jackknife+** avoids the split: every ion serves as both
fit and calibration, so n_cal = n = 23/29/41 and every interval is informative. The price
is stated rather than hidden - jackknife+ proves coverage of 1-2*alpha, so a nominally 90
per cent interval is **guaranteed at 80**, and both numbers travel with every interval.

Widths at a typical ion run from **0.55 per cent of CCS** (single-field DTIMS, [M+Na]+) to
**3.73 per cent** (TIMS, [M+H]+).

An interval is a property of the QUERY, not of the stratum. The first implementation
reduced it to one band per stratum and produced intervals of 170 to 217 square angstrom -
the entire range of steroid cross sections - which is what a spread looks like when it is
mislabelled as an interval.

### COVERAGE IS NOT EVIDENCE, and this is the most important limitation here

Nested leave-one-compound-out coverage on the nine applied strata: 27/29, 37/41, 21/23,
28/31, 42/46, 21/23, 28/31, 42/46, 21/23. Every one of those is **exactly**
ceil(n(1-alpha))/n. The rate is pinned by the quantile index and could not have come out
otherwise, so it says nothing whatever about whether the interval is well calibrated.
`CoverageCheck.is_evidence_of_calibration` returns False, always, as a property rather
than as a sentence in a docstring.

It is not useless: it is falsifiable against GROSS error. Scaling every residual by 0.5
takes coverage to 0.759, by 0.25 to 0.448, by 0.1 to 0.069 - all below the 0.80 guarantee.
So the check catches an interval that is plainly too narrow and catches nothing subtler.
**The informative figures are the WIDTH and the TAIL RATIO**, and those are what the
reports carry.

Three of the four independent designs proposed coverage as validation and one gated
`DataMaturity.VALIDATED` on it.

### The scope caveat is structural

`scope.ScopeStamp` has NO `scope` field. `scope` is a computed property of `studies`, so
there is no constructor argument to set and widening the claim requires naming a second
study - which is data somebody has to produce. `ComparisonScope` has two members and
neither is "interlaboratory reproducibility": the claim is not representable anywhere in
this package, so it cannot be recorded, serialised or returned.
`assert_may_be_quoted_as` raises for that claim on EVERY stamp however many studies it
holds, because a comparison between platforms is not a reproducibility figure for either
of them.

`contracts.ScopeReport` is REQUIRED on every `HarmonizedEstimate` with no default, and
re-derives the scope from the studies so a hand-built report cannot widen it.
`HarmonizeResponse` refuses a validated maturity beside a within-study scope.

**The scope is not uniform across the platform pairs, and this is a refinement of the
obvious statement.** The drift-tube and trapped-ion values are the authors' own; the
travelling-wave values are republished from doi:10.1021/acs.analchem.9b05247 and are
themselves an interlaboratory average over four instruments. So a DTIMS-versus-TIMS
correction is within one laboratory, and a TWIMS-versus-anything correction is one
laboratory against an aggregate from another paper. Neither is an interlaboratory
reproducibility figure. See 7D.

### What M4 is not

- **Not validated, and validation is not reachable from this corpus.** The maturity stamp
  is PROVISIONAL and cannot be otherwise while the scope is within-study, because
  validation means checked against data the model was not fitted on and one study has
  none by definition.
- **Not tested on anything but steroids.** 66 compounds, one chemical class, three
  adducts, all small molecules. Nothing here has seen a glycan, a peptide or a protein.
- **The corrections are small and the intervals are not.** A typical correction moves a
  value by 0.1 to 1.0 per cent; a typical interval is 0.5 to 3.7 per cent wide. For most
  ions the interval contains the uncorrected value, which is an honest statement of how
  much this corpus supports and not a flaw in the arithmetic.
- **The rank interval on the slope is a normal approximation** (the Kendall tau variance),
  and at n=23 it is indicative rather than exact. It is also the BASIS SELECTOR, so a
  stratum near the boundary could pick the other basis on a slightly different corpus.
- **Nothing measures whether a LINEAR correction is the right model.** A slope and an
  intercept fitted robustly are still a slope and an intercept.
- **The leverage check is silent on ten strata and INERT on eight.** Measured leverage
  runs 0.044 to 0.340 per cent against a 0.5 per cent limit, so
  `correction_driven_by_outliers` fires nowhere - and on 8 of the 18 strata
  `slope_leverage_percent` returns None outright, because they have no outliers to exclude
  or too few points left once they are. "It never fires" and "it cannot fire here" are
  different statements and the second is the one that matters: a check returning None on
  44 per cent of strata is not covering them. The threshold is NOT lowered to make it
  fire, because it is anchored to published stepped-field DTIMS reproducibility and
  tuning a guard until it triggers is fitting the guard to the data.

## 7G. M5: the deployable MVP, and what "deployable" does not cover

Built 19 September 2026. `/harmonize` is wired to the fitted model, there is a
demonstration that runs one real ion through every stage, and the package installs and
serves from a clean virtual environment - verified, not assumed.

### What the API does and does not do

Three outcomes, and which one a caller gets is the useful part:

- **200 with a harmonized value** where the measurement's platform and calibration group
  match one of the nine applied strata. It carries both coverage figures, the interval, the
  grade, the scope and the provenance.
- **No value, with a per-measurement reason**, where the model does not cover it. The field
  is `not_harmonized_because` and it is REQUIRED whenever the value is absent - a validator
  on the contract enforces it, because an absent value with no reason is indistinguishable
  from an oversight.
- **501 for the whole request** where nothing in it could be harmonized. A response full of
  absent values is a refusal and should read as one to anything checking the status alone.

**Nothing is extrapolated.** A correction that grades `unsupported` has its VALUE WITHHELD
rather than returned with a warning, because that grade's definition is "do not use this
number" and handing over a number while saying not to use it is a contradiction a caller
resolves in favour of the number. The grade and its reasons are still returned, so the
caller learns why.

### What a submitted record's licence status does, and does not, decide

**Found by adversarial testing on 20 September 2026, not by a test.** `/harmonize` applied
no licence check to submitted measurements at all: a record claiming `unverified`,
`excluded`, `non_commercial_no_derivatives` or `synthetic_fixture` was corrected and a
number returned. `assert_trainable` and `claim_problem` were called nowhere in `api.py`.

The two halves of that turned out to need opposite answers.

**`synthetic_fixture` IS NOW REFUSED.** It is not a licence claim - it declares that the
record was built in code and is not a measurement. The loader has always refused the
identical claim in a file, with the reasoning that "a row in a file is a real record"; the
API accepted it over HTTP. **The same declaration was fatal at one entrance and ignored at
the other**, and a number derived from a record that declares itself invented is
indistinguishable, once returned, from one derived from data. Refused now whether or not a
model is loaded, because the objection is to the record rather than to the state of the
service.

**THE LICENCE STATUSES DELIBERATELY DO NOT GATE THE ANSWER**, and this is a decision
recorded rather than an oversight left standing:

- The licence gate governs what may enter a FIT. The model is already fitted, on records
  that passed the gate. A submitted measurement enters no fit and changes no parameter.
- **The status is self-asserted and unverifiable.** A caller refused for `excluded` edits
  the field to `open_attribution` and resubmits. A gate on a field the caller supplies
  protects nothing and obstructs only honest callers. Shipping it would be theatre.
- The caller's relationship to their own data is not knowable here. They may be its author.

**Responsibility for the input's terms therefore remains the caller's**, and that is now a
stated disclosure in `contracts.SourceProvenance` rather than an inference. It matters most
for `non_commercial_no_derivatives`, where the no-derivatives clause bears directly on the
fact that a harmonized value IS a derivative of the input: a caller holding data under
those terms is the party making a derivative of it, and this service does not and cannot
check that. Recorded as a reservation rather than treated as settled.

### A response could contradict itself about a licence, and now cannot

Also found by the same testing. `SourceProvenance.reuse_status` was named as though it were
established fact and carried the RECORD'S CLAIM, while the registry contributed only its
licence TEXT - so a response could show `reuse_status: "excluded"` beside
`licence: "ACS AuthorChoice open access"` and a reader had to notice the contradiction
unaided. CLAUDE.md constraint 2 says a claim is valid only when a registry entry backs it,
and nothing said whether one did.

Worse, `provenance_of`'s own docstring already promised the two were "reported together
precisely so a caller can see when they disagree". **The intention was written down and half
implemented** - the ninth instance of the pattern in section 4.5, and the second where a
docstring stood in for the code.

Three fields now, and the names carry the distinction:

| field | what it is |
| --- | --- |
| `reuse_status_claimed` | what the submitted record says. A claim, self-asserted. |
| `reuse_status_in_registry` | what the registry holds for that DOI, or None if unregistered. |
| `claim_backed_by_registry` | whether they agree. **None** where unregistered, so there is nothing to agree with; **False** is the interesting case. |

### Deployment is from a checkout, and a wheel is not enough

The seed CSVs are not declared as package data, so a built wheel carries no measurements
and serves a model-less API. **This is a licensing decision, not an oversight.** The steroid
data is `academic_only` and carries an attribution obligation; bundling it into a
redistributable artefact is a decision nobody has made, and making it silently as a
packaging convenience is exactly the kind of thing this repository exists not to do.

An editable install from a clone resolves `data/seed` in the source tree and serves a real
model. That is the supported deployment and it is what was verified.

### Model versioning: DONE, and what it does and does not promise

Every response carries a `model_version` with two sha256 digests, and `/health` carries the
same pair, so an answer can be reproduced or told apart from another.

| | |
| --- | --- |
| `corpus_sha256` | over the records BEHIND THE FIT |
| `parameters_sha256` | over everything that decides an answer: which stratum applies, on what basis, the slopes, intercepts and offsets, and the leave-one-out residuals that set every interval |

Two rather than one, because the pair says WHAT changed: a different corpus with the same
parameters means the data moved without moving the fit; the same corpus with different
parameters means the code did. Measured: changing one cross section by 0.001 square angstrom
moves both.

**What "corpus" does NOT mean, because the name suggests otherwise.** It covers the records
that feed a correction, not every file on disk. The Struwe seed files pair nothing, so
adding or removing them leaves the digest unchanged - and leaves every answer unchanged too,
which is the point. The digest moves exactly when an answer could. Hashing every file would
be a weaker promise: it would make two identical answers look like they came from different
models whenever unrelated data moved.

The model is still refitted at every startup. That stays, and it is now safe: the same data
gives the same digest, so two deployments that agree can be shown to agree.

### Six deliberate omissions, each with what closing it would take

Recorded as decisions rather than as a to-do list. Each needs somebody to DECIDE, not
somebody to code, and that is why none of them is half-built. Four of the six wait on the
same unanswered question - who is allowed to call this service - and building any of them
before that answer means building it wrongly.

**1. No authentication and no authorisation.** Anybody who can reach the port can call every
endpoint. *To close it:* a decision from the CEO about who may call this - internal only,
named collaborators, or public - because that decision determines the mechanism. An internal
service needs a network boundary and nothing else; named collaborators need API keys and a
way to revoke them; public access needs accounts, and accounts mean a user store, a password
or token flow, and somebody who resets them. Building any of the three before the decision
means building two of them wrongly. The licence terms bear on this before the engineering
does: the seed data is `academic_only`, so an interface serving derived values to unknown
callers is a licence question first.

**2. Loopback binding is the only access control there is.** `python -m wmxccs` defaults to
`--host 127.0.0.1`, so out of the box the service is reachable only from the machine it runs
on. This is a real decision and it is load-bearing precisely BECAUSE of omission 1: it is not
defence in depth, it is the whole of the defence. `--host 0.0.0.0` exists and works, and
using it publishes an unauthenticated service to whatever network the host is on. A test
asserts the default so it cannot drift. *To close it:* omission 1. Until then, anyone
deploying this behind a reverse proxy or on a shared host is responsible for the boundary,
and should know that nothing inside the application will stop a request.

**3. No CORS policy.** Nothing sets `Access-Control-Allow-Origin`, so no browser page on
another origin can call this. That is the safe default and it is deliberate. *To close it:*
the same decision as omission 1, plus an explicit list of origins. A permissive `*` would be
the wrong answer for academic-only data regardless of who is asking.

**4. No rate limiting.** One caller can issue requests as fast as the process will serve
them. Request size is bounded - 1000 measurements per request, and every text field capped -
but that is not a rate limit and is not offered as one. *To close it:* it follows omission 1,
and needs one more fact that nobody has: what a legitimate caller's peak volume looks like,
which is unknowable while there are no callers. Guessing a limit now would either throttle a
real user or protect nothing.

**5. No seed data in a built wheel.** The seed CSVs are deliberately not declared as package
data, so `pip install wmxccs` from a wheel gives a working API whose `/harmonize` answers 501
for everything. Deployment is from a checkout. This is a licensing decision rather than an
oversight: the steroid data is `academic_only` with an attribution obligation, and bundling
it into a redistributable artefact is a decision nobody has made. *To close it:* somebody
with the authority to redistribute that data says so in writing, and the attribution
obligation is satisfied inside the artefact - or the wheel ships with a loader that fetches
the data from a location the licence does cover.

**6. A request never grows the corpus** - a design choice rather than an omission, and it
should stay one. A caller submitting two platforms' values for one ion gets each corrected
against the stored model; the pair they sent is added to nothing. Accepting submissions into
the corpus would mean taking unverified data through the licence gate on a stranger's
assertion, refitting the model between requests so that two callers get different answers
from one version, and losing the provenance chain that makes any of these numbers quotable.
*To close it:* a curation queue where submissions are held, licence-checked by a named
person, and merged deliberately - which is a product, not a feature.

## 8. Scope of the test suite

The tests assert the constraints in CLAUDE.md, not only the happy path, and the
mutation catalogue is what demonstrates that they bite. But:

- the wmxccs catalogue holds 282 mutations against that package's modules, and is closed at that number; tools/mutation/catalogue_glycan.py holds 35 more against src/wmxglycan, and runner.all_mutations() joins the two so a sweep covers 317. It began smaller than
  the glycan platform's 154 because 46 of those anchored into modules that do not come
  across and 26 into modules not in this milestone; it has since passed it. The floor in
  the catalogue test goes up, never quietly down;
- a mutation that survives is a behaviour with no test behind it. There are no
  documented expected survivors in this catalogue, so any survivor is a finding;
- numerical behaviour IS tested as of M4 - the robust slope, the leave-one-out
  intervals, the conformal floor, the digests and the served precision all have tests
  that fail when the arithmetic changes. What is still NOT tested is accuracy against an
  independent reference, because no such reference exists for these ions: every figure
  here is checked against the corpus it was fitted on, or against a synthetic fixture
  whose answer is known by construction. A test that the arithmetic is self-consistent
  is not a test that the answer is right, and nothing in this repository claims it is.
