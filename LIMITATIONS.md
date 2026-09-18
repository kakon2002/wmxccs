# Limitations

What this repository does not do, does not know, and must not be read as claiming.
Written at M0. Everything here is a live limitation unless it says it is closed.

---

## 1. The thing that matters most

**There is not one cross-platform matched ion in this repository.** The corpus is
117 measurements, all travelling-wave, all from one laboratory, and the platform's
entire premise is comparing the same ion measured on two different platforms.
Nothing downstream of matched-ion construction can be built, validated or
demonstrated on the data now held.

This is not a shortage that more of the same data would fix. It is
over-determined, and each of these alone is sufficient:

1. both seed files are TWIMS, so there is only one platform;
2. the analyte sets are disjoint — the 2015 file is high-mannose N-glycans
   (Man3 to Man9Glc), the 2016 file is milk oligosaccharides (LNH, LNnH, LNT,
   LNnT). No molecule appears in both;
3. the 2015 rows carry `drift_gas = UNSTATED` and the 2016 rows carry `He`, and
   gas is part of the matched-ion key, so even a shared analyte could not pair.

`readiness.py` reports this as a blocker, in those terms, rather than as a small
number. That is the M0 deliverable.

## 2. No number in this repository describes model performance

No harmonization model has been fitted. No CCS value has been corrected,
harmonized or predicted. Every readiness report carries a `MaturityStamp` of
`provisional`, and the count beside it is the number of cross-platform matched
ions behind it, which is zero. There is no code path that can produce a
`validated` stamp, because the milestone that would earn one does not exist.

## 3. Licences: what is settled and what is not

Four of the eleven registered sources are **unverified**, and the default-deny
gate means none of them may be ingested. Each carries a `what_to_check` field
naming the single thing that would settle it.

The one that matters is the **steroid interplatform study**
(DOI 10.1021/jasms.2c00196): 87 steroids, 142 values, DTIMS + TWIMS + TIMS. It is
the only identified source that would supply cross-platform matched ions
directly. ACS publishes open access under both CC-BY and CC-BY-NC-ND, and which
one applies has not been read. One page decides whether a benchmark set exists.

Nothing in this repository has read a licence. Every entry in `sources.py` records
a report by a named person on a stated date. The code does not verify those
reports and does not claim to.

**The commercial-versus-research question is still open with the CEO.** Until it
is answered the platform is treated as commercial, which is what blocks CCSBase.
If the answer is internal research, CCSBase becomes usable and the data problem
largely changes shape.

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
an empty catalogue. Latent only: the catalogue holds 110, so nothing reaches it.
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

Five separate instances in this repository so far, in five different shapes. They
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

WHAT THE FIVE HAVE IN COMMON, and what to do about it. Green is not evidence. The
question that catches all five is not "do the tests pass" but "what would have to
break for this to go red, and is that the thing I think I am protecting". Note that
three of the five were found by somebody asking that question about a mutation that
had apparently behaved correctly, and none by a test going red. Concretely:

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
  of these five were one helper quietly giving two different things the same
  identity, and in both the symptom was a count that was obviously wrong to anyone
  who looked at it.

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

## 8. Scope of the test suite

The tests assert the constraints in CLAUDE.md, not only the happy path, and the
mutation catalogue is what demonstrates that they bite. But:

- the catalogue holds 195 mutations against twelve modules. It is smaller than the
  glycan platform's 154 because 46 of those anchored into modules that do not come
  across and 26 into modules not in this milestone. The floor in the catalogue test
  is the real current count and goes up, never quietly down;
- a mutation that survives is a behaviour with no test behind it. There are no
  documented expected survivors in this catalogue, so any survivor is a finding;
- nothing here tests numerical accuracy of anything, because nothing numerical is
  computed yet.
