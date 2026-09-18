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
  reading Hofmann 2014, unblocks nothing at all. Both must be resolved.

### Named in the objective, licence unverified or blocking

| Source | Status | Detail |
|---|---|---|
| CCSBase, ccsbase.net | **Blocked for commercial use** | Terms restrict use to academic, non-commercial purposes. Commercial users directed to UW CoMotion. Read 11 Sep 2026 |
| METLIN-CCS, Nature Methods 2023 | Unverified | Paper says "freely available" at METLIN, XCMS Online and PanoramaWeb. That is access, not a licence. 185,589 TIMS values, 27,633 standards. Largest TIMS source |
| Bush Lab CCS Database, biophysicalms.org/ccsdatabase | Unverified | Page says only "please cite the appropriate publication(s)". Native proteins, denatured proteins, peptides, small molecules. Most measured in both He and N2 by DTIMS. Downloadable Google Sheet. Ideal for the biopharma layer if usable |
| Steroid interplatform study, JASMS 2022, DOI 10.1021/jasms.2c00196 | Unverified, ACS | 87 steroids, 142 values, DTIMS + TWIMS + TIMS. The obvious benchmark set for RQ1 |
| Bayesian harmonization study, Anal Chem 2026, DOI 10.1021/acs.analchem.5c06667 | Unverified, ACS | 840 measurements, 347 compounds, three platforms. If its data is released and licensed, it is a ready-made matched-ion set |

### Previously excluded, still excluded

AllCCS2 (no licence posted). GlycoMob (no terms, possibly offline). Manabe 2022
(CC BY-NC-ND). Sastre Toraño 2025 Nat Commun (CC BY-NC-ND, both clauses block).

## The open question that decides the data plan

Is this platform internal research or a product? The objective document is titled
"research concept document" and its MVP exposes results through an API. If the
platform is internal research only, CCSBase and academic-only sources may be usable
and the data plan is comfortable. If it is a product, or feeds one, they are not,
and the matched-ion benchmark must come from licence-clean papers, of which none
are confirmed yet.

This question has been put to the CEO. Until answered, treat the platform as
commercial and the academic-only sources as blocked.

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
