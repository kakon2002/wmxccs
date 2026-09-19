# CONTEXT.md

Everything established before this repo existed. Read once, then trust it.

## Where this came from

Wellmatix spent 12 to 17 September building a glycan isomer and CCS prediction
platform (`Project2`, package `wmxglycan`, 1,868 tests). On 17 September the CEO
replaced the objective with a cross-platform CCS harmonization pipeline, to be built
as a separate platform. This repo is that platform.

The glycan work is not wasted. Its data layer is the data layer this project needs.
Its glycan-specific layer is not needed and does not come across.

> **Counts corrected 18 September, from the repository rather than from memory.**
> The test figure was written here as 1,796. `pytest` collects **1,868** cases from
> 705 test functions, and the string "1796" appears nowhere in that repository. The
> mutation catalogue holds **154** entries, not the 136 quoted in the plan document.
> "Roughly 40 percent of its code carries over" is also low as an import fact: the
> transitive closure of the port set below is 15 of the 20 package modules, four of
> which are on the do-not-port list. The glycan core is not a separable layer that
> can be left behind; cutting it out is most of the porting work.

## What ports from Project2, and what does not

Port these, adapting names and removing glycan assumptions:

| Module | What it does | Change on port |
|---|---|---|
| **`reuse.py`** | **The ReuseStatus enum, DEFAULT_REUSE_STATUS, SYNTHETIC_FIXTURE and the three tier predicates** | **None beyond package name. MANDATORY, and it was missing from this list** |
| `licensing.py` | The default-deny training gate. It re-exports `reuse.py`; it does not define the statuses | None beyond package name, once `reuse.py` is there |
| `sources.py` / registry | Licence records with reader and date, refuses unbacked claims | Add the sources listed below, **and cut its `sugarbase` import first** |
| `models.py`, CCSMeasurement only | Platform, gas split, calibrant, adduct, conformer, uncertainty type, provenance | Generalise the analyte: not a glycan, any ion. See "identity" below |
| Loader hardening (in `measurements.py`) | Strict CSV, placeholder rejection, date format, row conservation | The hardening itself needs none. The module around it imports `features.py` and `splits.py` and has to be cut free of both |
| Adduct normalisation (in `models.py`) | Sorted components, charge consistency | None |
| Conformer sibling check (in `measurements.py`, not `models.py`) | Two values, one structure, one adduct | None |
| Readiness and refuse-to-fit (`Readiness` is inside `training.py`) | Blockers vs warnings, thresholds as named constants | Rethreshold for matched-ion counts. Its thresholds come from `splits.py` and `evaluation.py`, which are not on this list |
| `tools/mutation` | Catalogue-driven mutation harness with anchor tests | Port whole, rewrite catalogue for new modules. The only coupling is one line naming the package directory |

Do not port: `composition.py`, `enumeration.py`, `glycan_graph.py`, `enzymes.py`,
the SugarBase loader, the 44-column glycan featuriser, any biosynthetic rule.

**Three things this table got wrong, found while porting on 18 September.** Each
cost real time, so they are recorded here rather than left to be rediscovered.

1. **`reuse.py` was missing and is mandatory.** `ReuseStatus`, the default status
   and the tier predicates all live there, not in `licensing.py`, which only
   re-exports them. `licensing.py`, `sources.py` and `models.py` all import it.
   Copying only the files this table named gives an `ImportError` in three places.
   Do NOT respond by pasting the enum into `licensing.py`: `licensing.py` imports
   `sources.py`, which would then have to import `licensing.py`, and `reuse.py`
   exists precisely to break that cycle. Its own docstring says so.

2. **`sources.py` does not import cleanly.** It does
   `from .sugarbase import SUGARBASE_ATTRIBUTION, SUGARBASE_LICENCE, dataset_version`
   and then CALLS `dataset_version()` at module scope while building its SugarBase
   dataset entry, which reads the installed `glycowork` version. `sugarbase.py` is
   on the do-not-port list. So `import sources` fails outright, and because
   `licensing.py` imports `sources.py`, so does every test that touches the gate.
   This is the single edge to cut first, before anything else can even be run.
   Keep the `DatasetLicence` mechanism and register something real in it, or the
   dataset-backed branch of the gate becomes unreachable code that nothing notices
   is dead.

3. **The port set is larger than the named modules.** `measurements.py` imports
   `features.py` (do-not-port) and three names from `splits.py`, one of them
   private. `training.py` imports from `splits.py`, `evaluation.py` and
   `contracts.py`. None of those four is named anywhere in this document, and the
   `Readiness` figures cannot be produced without deciding what to do about each.

The port list should also have said that **`struwe2015.py` and `struwe2016.py`**
exist and are the only things that can read the two transcriptions. They are not
named here, and loading the seed files is impossible without either porting them
or replacing them with a conversion step.

## Identity: the one schema change that matters

`CCSMeasurement.glycan` becomes `CCSMeasurement.analyte`, and an analyte is one of:

- small molecule: InChIKey required, SMILES optional
- peptide: sequence, modifications
- glycan: composition OR a structure identifier, and the structure where resolved
- protein or subunit: accession or sequence, native or denatured state
- intact antibody: identity, glycoform, native or denatured, CIU state
- ADC: antibody identity, DAR or conjugation state, linker payload class

The matched-ion key is analyte identity + adduct + charge + gas + structural state.
For proteins and antibodies, charge state is part of identity: a native 24+ ion
and a denatured 40+ ion of the same antibody are different matched-ion keys.

**Corrections to this section, 18 September.**

- This said "glycan: composition, structure if resolved (reuse the existing
  representation)", which conflicts with `composition.py` being on the
  do-not-port list. What was actually done: the canonicalisation half of that
  module carries over (residue parsing and canonical spelling, which is identity
  machinery, because without it `Hex5HexNAc4` and `HexNAc4Hex5` are two matched
  ions instead of one); the N-glycan plausibility rules and the Man3GlcNAc2 core
  check do not, because they are biosynthetic rules.
- **Composition is OPTIONAL**, not required. This platform takes any ion, and a
  steroid or a protein has no composition to give. It is also the only honest
  option for the 2016 transcription, which carries an IUPAC structure and no
  composition column: the glycan platform DERIVED the composition from the
  structure using a monosaccharide-class table, which is glycan chemistry and
  does not come across. Those rows record no composition and are identified by
  their structure, which is finer anyway. A glycan must still state at least one
  of the four identifiers, or it has no key at all.
- **The glycan identity key is the finest identifier stated, NOT the composition.**
  Worth being explicit, because the glycan platform's rule was the opposite and
  copying it here would be a serious error. That rule ("split on composition, not
  structure") was about TRAINING FOLDS, and its reasoning does not transfer: a
  prediction request carries no structure, so a structure-derived split key is
  finer than the model's own identity. Matching is a different question. LNH and
  LNnH share the composition Hex4HexNAc2 and differ by 11.4 per cent as [M-H]-,
  so keying them alike would pair two different molecules and report the
  difference between them as inter-platform bias. Grouping for a fit and matching
  for a comparison are two different keys and one name must not serve both.
- The reducing-end label and the derivatisation moved from the MEASUREMENT, where
  the glycan platform kept them, onto the GLYCAN ANALYTE. A steroid has no
  reducing end, and a field that can only ever be null on five of six analyte
  kinds is a glycan assumption leaking into every record.
- **A new adduct form: `[M+24?]24+`**, meaning "twenty-four charges, carrier not
  stated". Native-MS papers routinely report a charge state without saying whether
  the carriers are protons, sodium or ammonium, and the Bush Lab protein data is
  expected to be full of them, so this is the common case rather than the exotic
  one. Without it such a record cannot be built at all, because the adduct is
  mandatory and its charge must match; writing `[M+24H]24+` instead would assert
  protons, which the paper does not say. It blocks training like the other
  unknowns, and its key is made unique to its own record, so two sources both
  reporting "the 24+ ion" can never be paired as though they were the same ion.

## New fields required by the objective

cIMS, stored separately from single-pass TWIMS:
pass number, effective path length, TW velocity, TW height, arrival-time correction
method, wrap-around flag, single-pass vs multipass CCS.

All platforms: calibration reference lineage. Whether a value is primary
(stepped-field DTIMS) or calibration-derived, and if derived, what reference set.
This matters because TWIMS calibrations often derive from DTIMS values, and the
objective names reference circularity as a risk.

**Correction, 18 September: most of this already existed and ported.** The
objective document has no field for it, but the glycan repository did:
`ccs_is_calibrated` derives primary-versus-derived from the platform and the DTIMS
method, `calibrant_reference` records what the calibrant was referenced against,
and a validator already refuses a calibrant reference on a primary value. What was
genuinely new is the STRUCTURED reference set: which publication the reference
values came from, and what platform they were themselves measured on. Those two are
what make circularity detectable at all, because `calibrant` alone does not say -
two laboratories both writing "dextran" may be using different ladders with
different assumed values.

Both seeded files leave the reference platform null, deliberately. "DTCCS" in the
2015 cell plainly suggests a drift tube, and that is a reading of an abbreviation
rather than of a paper, so `traces_to_primary` returns "not known" for every seeded
record. Resolving it means reading Hofmann 2014, Anal. Chem. 86, 10789.

## Data sources and their status

Every entry here was checked by a named person on a stated date, or is marked
unverified. Do not upgrade any of these without a registry entry.

### Usable now

| Source | Licence | Evidence | What it holds |
|---|---|---|---|
| Struwe et al. 2016, Chem Commun 52:12353, DOI 10.1039/c6cc06247d | CC BY 3.0 | RSC article page, issue listing, PDF header. Read by Shawon Chakrabarty Kakon, 15 Sep 2026 | 28 TWIMS values, 4 milk oligosaccharides, 4 adducts. Already transcribed |
| Struwe et al. 2015, Analyst 140:6799, DOI 10.1039/c5an01092f | CC BY 3.0 | PDF footer every page. Read 17 Sep 2026 | 89 TWIMS values, 7 high-mannose N-glycans, 6 adducts. Drift gas UNSTATED. Already transcribed |

Both are single-platform TWIMS from one lab. They give zero cross-platform matched
pairs. They are useful for the glycan class layer and for nothing else in the MVP.

**Confirmed by loading them, 18 September, and stronger than stated.** 117 rows in,
117 held, none lost. Zero matched pairs is **over-determined**: any one of these
alone is sufficient, so fixing one changes nothing.

1. both files are TWIMS, so there is one platform;
2. the analyte sets are disjoint. The 2015 file is high-mannose N-glycans, Man3 to
   Man9Glc; the 2016 file is milk oligosaccharides, LNH, LNnH, LNT and LNnT. No
   molecule appears in both, and no identity atom is shared between them;
3. the 2015 rows carry `drift_gas = UNSTATED` and the 2016 rows carry He, and gas
   is part of the key, so even a shared molecule could not pair.

Two further things this table did not say, both of which matter:

- **24 of the 28 rows in the 2016 file clear the gate, not 28.** Four are held as a
  suspected shared peak: LNH and LNnH report identical values for [M+H]+ (228.9) and
  [M+Cl]- (245.0), almost certainly one unresolved peak transcribed against both
  compounds. That hold is **derived by the loader**, not marked anywhere in the CSV.
  Port the shared-peak detector away and all 28 clear, with nothing to show it
  happened.
- **NONE of the 89 rows in the 2015 file clears the gate, and the drift gas is only
  half the reason.** Every one of them also carries `uncertainty_type = unknown`
  with no spread, which is an independent blocker. Resolving the gas alone, by
  reading Hofmann 2014, unblocks nothing at all.
- **And the second half is probably permanent.** The owner reports, 19 September
  2026, that the Analyst ESI carries no uncertainty column at all. If so, `unknown`
  is a correct and complete record of what the source states, not a transcription
  gap, and **these 89 rows are not pending curation**. Do not put them on a list.
  Unblocking them needs a different source, or an explicit decision that a value
  with no reported spread is usable for some purpose that does not need one - never
  a default. See LIMITATIONS.md section 7.1.

### Named in the objective: usable since the CEO answered, or still unverified

| Source | Status | Detail |
|---|---|---|
| Steroid interplatform study, JASMS 2022, DOI 10.1021/jasms.2c00196 | **`academic_only`, INGESTED** | 87 steroids, 142 values, DTIMS + TWIMS + TIMS. Converted by `tools/ingest_steroid.py`: 521 records, **142 cross-platform matched ions** across up to four platforms, 517 clear to train. The benchmark set for RQ1, and the first real data in the repository. See LIMITATIONS.md 7D |
| CCSBase, ccsbase.net | **`academic_only`, RETRIEVED, blocked on drift gas** | 25,020 records read in full 19 Sep 2026: DT 10,967, TW 8,529, TIMS 5,524, 23 calibration methods, 36 primary references, with SMILES and a per-record id. **The drift gas is not a column**, so every record is unmatchable under the keying rule and CCSbase yields zero matched ions until the gas is resolved from the 36 papers. Not a downloadable file: `/download` is the batch-query template. See LIMITATIONS 7E |
| Bush Lab CCS Database, biophysicalms.org/ccsdatabase | **`academic_only`, DOWNLOADED, not yet ingested** | 8 sheets, ~6,500 rows: native-like protein cations (1,000), complexes (989), denatured (1,000), polyalanine (33), anionic homopolymers (1,000), other peptides (1,000), small molecules (24), MicroSource collection (1,441). The only identified source that would give the biopharma layer real data. 41 corrupt date-serial cells a converter must account for. See LIMITATIONS 7E |
| METLIN-CCS, Nature Methods 2023 | Unverified | Paper says "freely available" at METLIN, XCMS Online and PanoramaWeb. That is access, not a licence. 185,589 TIMS values, 27,633 standards. Largest TIMS source |
| Bayesian harmonization study, Anal Chem 2026, DOI 10.1021/acs.analchem.5c06667 | Unverified, ACS | 840 measurements, 347 compounds, three platforms. If its data is released and licensed, it is a ready-made matched-ion set |

All three usable entries are recorded as `academic_only`, **not** `open_attribution`.
The distinction is enforced: it is what makes the gate refuse them again the moment
`reuse.PLATFORM_USE_CONTEXT` is set back to `COMMERCIAL`.

### Previously excluded, still excluded

AllCCS2 (no licence posted). GlycoMob (no terms, possibly offline). Manabe 2022
(CC BY-NC-ND). Sastre Toraño 2025 Nat Commun (CC BY-NC-ND, both clauses block).

## The question that decided the data plan: ANSWERED 19 September 2026

Was this platform internal research or a product? **The CEO answered: academic and
research use, not commercial.**

That is recorded in code as `reuse.PLATFORM_USE_CONTEXT = UseContext.ACADEMIC_RESEARCH`
- a named constant, deliberately not a parameter and not a configuration value, so
that changing it is a visible edit to the package rather than a deployment setting.
It reopened the three academic-only sources above.

Two things about how it is recorded, both deliberate:

1. **The sources keep their own terms.** A reuse status says what a SOURCE'S TERMS
   permit; the use context says what WE are. CCSbase's terms still restrict it to
   academic non-commercial use and the registry still says so. Nothing was relabelled
   to make it loadable. `can_train_commercial` still exists and still answers the
   commercial question, so "would this be usable if we were a product" remains a
   question the code can answer at any time.
2. **Every entry that is usable only because of the context records who decided
   that.** `sources.py` REQUIRES a `context_basis` for exactly those entries and
   refuses to construct one without it. An entry usable on a decision rather than on
   its own terms would otherwise read like an entry whose terms permit anybody.

If the platform is ever commercialised, one constant changes and the gate refuses all
three again without anybody re-reading a licence page.

## M4 is built, and one piece of the guidance for it turned out to be wrong

Built 19 September 2026 on the 93 four-platform matched ions that clear the gate.
`robust.py`, `scope.py`, `harmonization.py`. See LIMITATIONS 7F.

The guidance was: report both the slope-derived and median-derived correction, prefer a
robust fit, keep the leverage check. Two of those held. The third needs correcting, and
the correction came from the data rather than from an argument:

**"Prefer a robust fit" is right; "prefer the robust SLOPE" is not.** Passing-Bablok
does exactly what was wanted - on one clean line, wrecking a single point moves the Deming
slope from 1.04 to over 6 and does not move the robust slope at all, and in negative mode
it flattens Deming slopes of 1.11-1.19 back to 1.00-1.07. But having flattened them, the
slope is then indistinguishable from 1 in 8 of 18 strata, and an affine correction fitted
through a slope that is indistinguishable from 1 is a constant offset with extra
extrapolation risk. So the headline basis is DERIVED:

    ROBUST_SLOPE   where the Passing-Bablok rank interval on the slope excludes 1
    MEDIAN         otherwise

All three corrections are still computed and reported every time, exactly as instructed.

**The leverage check is kept, it never fires on this corpus, and on half the strata it
cannot even be computed.** Measured: leverage runs 0.044 to 0.340 per cent against a limit
of 0.5, so `correction_driven_by_outliers` is silent on all 18 strata - and in 8 of the 18
`slope_leverage_percent` returns None, because those strata have no outliers to exclude or
too few points left once they are. So the check is silent on ten and inert on eight.

That is recorded rather than tuned. The threshold is anchored to published stepped-field
DTIMS interlaboratory reproducibility, and lowering it until something fired would be
fitting the guard to the data. But "it never fires" and "it cannot fire here" are
different statements, and the second one is the one that should worry a reader: a check
that returns None on 44 per cent of strata is not covering them.

**One number in this document was a pre-gate count presented as trainable.** "97 across
all four" is the histogram over all 521 records; over the 517 that clear the gate it is
93, the difference being four records the shared-peak check holds. Both are true and they
answer different questions. Corrected in place here and in LIMITATIONS 7D.

## M5: the deployable MVP

Built 19 September 2026. `/harmonize` is wired to the fitted model; `python -m wmxccs` or
`wmxccs-serve` starts it; `tools/demo_end_to_end.py` runs one real ion through every stage.
Installed and served from a clean virtual environment and verified there, not assumed.

Three things worth carrying forward:

**Nothing is extrapolated.** A correction grading `unsupported` has its value WITHHELD, not
returned with a warning. The grade means "do not use this number", and returning one while
saying that is a contradiction a caller resolves in favour of the number. The grade and its
reasons still come back.

**An absent value must say why.** `not_harmonized_because` is required by a contract
validator whenever there is no harmonized value. A platform the model does not cover, an ion
outside the fitted range, and a value already on the primary platform are three different
situations needing three different actions from the caller.

**DEPLOY FROM A CHECKOUT, NOT A WHEEL.** The seed CSVs are deliberately not package data, so
a built wheel carries no measurements and its `/harmonize` answers 501. That is a licensing
decision: the steroid data is academic_only with an attribution obligation, and bundling it
into a redistributable artefact is a decision nobody has made. Making it silently, as a
packaging convenience, is the sort of thing this project exists not to do. If a wheel
deployment is ever wanted, that decision has to be taken explicitly and recorded in the
registry.

Not done, and each needs a decision rather than code: authentication, rate limiting and CORS
(the server binds 127.0.0.1, which is the only concession); model versioning, since the model
is refitted from the seed files at every startup and a changed seed file changes the answers
with nothing recording that it did; and growing the corpus from submitted measurements, which
the API deliberately does not do.

## Domain facts already established

- Interlab DTIMS reproducibility: 0.29% RSD stepped field, 0.54% bias single field
  (Stow et al. 2017). TIMS within about 1% of DTIMS, TWIMS within about 2%.
- Adduct changes separability. For LNH/LNnH, [M+Na]+ differ 3.8%, [M-H]- differ
  11.4%, [M+Cl]- differ 0.08%. Adduct is part of the matched-ion key, not metadata.
- TWIMS values calibrated against helium references while measured in nitrogen
  are helium-referenced values. The gas the value refers to and the gas in the
  cell are separate fields. A stepped-field DTIMS value cannot claim a gas it was
  not measured in.
- Deprotonated ions can show multiple conformers of one structure (Struwe 2015,
  Man5/Man6/Man9). Two CCS values for one matched-ion key is legitimate when
  the conformer index differs. Do not average, do not drop.
- Single-field DTIMS is calibrated. Only stepped-field DTIMS is primary.

## What the MVP can and cannot claim in seven days

The objective document is a research programme. The MVP list has eight items. Items
1 to 5 (ingestion, matched-ion table, pairwise plots, statistics, first harmonization
model with cross-validation and uncertainty) are buildable in seven days if a
matched-ion dataset exists. Items 6 to 8 (cIMS layer, biopharma library, API and
scorecard connector) get scaffolds and stubs, not full implementations.

A hierarchical Bayesian model with proper validation is not a seven-day deliverable.
The MVP harmonization model is Deming regression per platform pair, optionally
stratified by molecular class, with grouped cross-validation and prediction
intervals. The Bayesian model is the v2 roadmap item and the code leaves a seam for it.

The objective document itself says generic harmonization is no longer novel, citing
the 2026 study. The novelty it wants is biopharma and cIMS, and public data for
both is sparse. The MVP therefore demonstrates the pipeline on whatever matched-ion
data is licence-clean, and the biopharma layer is where Wellmatix's own
measurements go.
