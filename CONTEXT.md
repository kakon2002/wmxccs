# CONTEXT.md

Everything established before this repo existed. Read once, then trust it.

## Where this came from

Wellmatix spent 12 to 17 September building a glycan isomer and CCS prediction
platform (`Project2`, package `wmxglycan`, 1,796 tests). On 17 September the CEO
replaced the objective with a cross-platform CCS harmonization pipeline, to be built
as a separate platform. This repo is that platform.

The glycan work is not wasted. Its data layer is the data layer this project needs.
Its glycan-specific layer is not needed and does not come across.

## What ports from Project2, and what does not

Port these, adapting names and removing glycan assumptions:

| Module | What it does | Change on port |
|---|---|---|
| `licensing.py` | ReuseStatus enum, default-deny training gate, SYNTHETIC_FIXTURE | None beyond package name |
| `sources.py` / registry | Licence records with reader and date, refuses unbacked claims | Add the sources listed below |
| `models.py`, CCSMeasurement only | Platform, gas split, calibrant, adduct, conformer, uncertainty type, provenance | Generalise the analyte: not a glycan, any ion. See "identity" below |
| Loader hardening | Strict CSV, placeholder rejection, date format, row conservation | None |
| Adduct normalisation | Sorted components, charge consistency | None |
| Conformer sibling check | Two values, one structure, one adduct | None |
| Readiness and refuse-to-fit | Blockers vs warnings, thresholds as named constants | Rethreshold for matched-ion counts |
| `tools/mutation` | Catalogue-driven mutation harness with anchor tests | Port whole, rewrite catalogue for new modules |

Do not port: `composition.py`, `enumeration.py`, `glycan_graph.py`, `enzymes.py`,
the SugarBase loader, the 44-column glycan featuriser, any biosynthetic rule.

## Identity: the one schema change that matters

`CCSMeasurement.glycan` becomes `CCSMeasurement.analyte`, and an analyte is one of:

- small molecule: InChIKey required, SMILES optional
- peptide: sequence, modifications
- glycan: composition, structure if resolved (reuse the existing representation)
- protein or subunit: accession or sequence, native or denatured state
- intact antibody: identity, glycoform, native or denatured, CIU state
- ADC: antibody identity, DAR or conjugation state, linker payload class

The matched-ion key is analyte identity + adduct + charge + gas + structural state.
For proteins and antibodies, charge state is part of identity: a native 24+ ion
and a denatured 40+ ion of the same antibody are different matched-ion keys.

## New fields required by the objective

cIMS, stored separately from single-pass TWIMS:
pass number, effective path length, TW velocity, TW height, arrival-time correction
method, wrap-around flag, single-pass vs multipass CCS.

All platforms: calibration reference lineage. Whether a value is primary
(stepped-field DTIMS) or calibration-derived, and if derived, what reference set.
This matters because TWIMS calibrations often derive from DTIMS values, and the
objective names reference circularity as a risk.

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
