> **HISTORICAL DOCUMENT.** This records what was asked in September 2026 and predates
> the glycan layer added on 26 September. It is kept as a record and is NOT a description
> of what this repository contains now; for that, read `README.md`.

# Cross-Platform CCS Harmonization Platform: Plan

**Prepared for:** James Kang
**Date:** 18 September 2026
**Target:** deployable version 25 September, 27 September at the latest
**Status:** new repository, separate from the glycan platform

---

## 1. What this platform does

It stores CCS measurements from DTIMS, TWIMS, TIMS and cIMS without merging them,
pairs the same ion measured on different platforms, quantifies the bias and
agreement between those platforms, and returns a harmonized CCS with an
uncertainty interval and a confidence grade. The original measurements are always
returned alongside the harmonized value and are never overwritten.

The unit that makes this work is the matched ion: molecule, adduct, charge state,
gas, and structural state together. Compound name alone never matches two records.

---

## 2. What carries over from the previous work

The glycan platform is not being extended, but its data layer is close to what
section 5 of the objective document specifies. These parts port directly and save
roughly two days.

| Already built | Used here for |
|---|---|
| Measurement record with platform, gas, calibrant, adduct, charge, conformer, uncertainty type, provenance | The database schema in section 5.2 |
| Separate fields for the gas a value refers to and the gas in the cell | Preventing helium-referenced and nitrogen values being pooled |
| Licence gate with a registry naming who read the terms and when | Section 17, data provenance and reuse |
| Loader hardening: strict parsing, placeholder rejection, row conservation | Section 12, dataset audit |
| Conformer handling: two values for one ion is legitimate | Section 8 and the cIMS conformer index |
| Refuse to fit when evaluation is impossible, warn when weak | Section 10, honest confidence |
| Grouped splitting so one analyte never spans train and test | Section 12, cross-validation |
| Mutation testing harness, 154 curated mutations | Verifying the guards actually bite |

What does not carry over: everything glycan-specific. Isomer enumeration,
biosynthetic rules, the structure database, glycan graph features.

## 3. What is new

Matched-ion construction across platforms. The cross-platform statistics module.
The harmonization model. cIMS metadata including pass number and path length.
Biopharmaceutical identity fields for antibodies, subunits, ADCs and glycopeptides.
Calibration reference lineage. Confidence grading. The API.

One field the objective document does not list but its own risk table requires:
**calibration reference lineage**. Whether a value is primary or calibration-derived,
and if derived, which reference set it came from. Section 17 names reference-value
circularity as a risk, where TWIMS calibrations trace back to DTIMS values. That
risk cannot be detected without this field.

---

## 4. The main technical risk

This platform needs the same ion measured on two or more platforms. That is the
entire premise. Nothing in the statistics or the model works without it.

Our existing 117 CCS measurements are all travelling wave, from one laboratory.
They contribute zero matched pairs. They are useful as a molecular class layer and
nothing more.

So the benchmark must come from published interplatform work. Status of each
source named in the objective document:

| Source | Content | Licence status |
|---|---|---|
| Steroid interplatform study, JASMS 2022 | 87 steroids, 142 values, DTIMS + TWIMS + TIMS. The ideal benchmark | Open access, supporting information free of charge. ACS uses two open licences and which one applies is **being confirmed** |
| METLIN-CCS | 185,589 values, 27,633 standards, TIMS | Described as freely available. That is access, not a licence. **Unverified** |
| Bush Lab CCS database | Native and denatured proteins, peptides. Most ions measured in both helium and nitrogen | Page asks only for citation. **Unverified** |
| CCSBase | Platform-aware experimental database | Terms restrict use to academic, non-commercial purposes. **Blocked for commercial use** |
| Bayesian harmonization study, Anal Chem 2026 | 840 measurements, 347 compounds, three platforms | **Unverified** |

**The decision needed:** is this platform internal research, or a product, or does
it feed one? If internal research, most of these sources are usable and the data
problem largely disappears. If it is commercial, only clearly licensed sources
qualify and the benchmark set is much smaller.

Everything else in the plan proceeds either way. Only the data question is blocked.

---

## 5. Validation targets, taken from published results

The steroid study gives concrete numbers to validate the harmonization model
against. It reported that 95 percent of ions fell within 1 percent bias for TIMS
and 2 percent for TWIMS relative to DTIMS, while under 1.5 percent of ions showed
biases as large as 7 percent.

That last figure is the important one, and it is what the confidence grading is
for. A model that reduces average disagreement while hiding the few ions that
genuinely do not transfer is worse than no model. The pipeline must flag those
rather than smooth them.

Correlation figures from the same study: 0.9949 for TWIMS against DTIMS, 0.9953
for TIMS against DTIMS, 0.9989 for TWIMS against TIMS. High correlation with real
bias underneath is exactly why the pipeline reports agreement as well as
association.

---

## 6. Schedule

| Day | Date | Work | Done when |
|---|---|---|---|
| 1 | Thu 18 | Repo, data layer ported, licence registry seeded | Existing records load, counts reported |
| 2 | Fri 19 | Identity model, cIMS fields, calibration lineage | Matched-ion key tested across analyte types |
| 3 | Sat 20 | Matched-ion construction, deduplication, provenance | Pairs found in the benchmark set |
| 4 | Sun 21 | Statistics: Pearson, Deming, Bland-Altman, bias, concordance, outliers | Reproduces the published steroid figures |
| 5 | Mon 22 | Harmonization model, grouped cross-validation, prediction intervals | Validated against published bias figures |
| 6 | Tue 23 | Confidence grading, API, dashboard | End to end run |
| 7 | Wed 24 | Packaging, deployment, documentation | Deployable |
| 8-9 | Thu 25 to Fri 26 | Buffer | |

Day 4 is the checkpoint that matters. If the statistics reproduce the published
steroid results, the pipeline is correct. If they do not, something upstream is
wrong and there is still time to fix it.

---

## 7. What version 1 will and will not include

**Included:** ingestion with provenance and licence control, matched-ion
construction, the full statistical comparison set, a harmonization model with
cross-validation and prediction intervals, confidence grading, and an API.

**Scaffolded but not populated:** the cIMS layer and the biopharmaceutical
reference library. Both have the correct fields and accept data, but public CCS
coverage for intact antibodies and ADCs is very thin, which the objective document
notes in its own risk table. These fill as Wellmatix generates measurements.

**Deferred to version 2:** the hierarchical Bayesian model. A Deming regression per
platform pair, with a class-stratified variant, is the seven-day deliverable. The
code leaves a defined seam for the Bayesian model to drop into.

This matters for positioning. The objective document states that generic
harmonization alone is no longer novel, citing the 2026 study. The differentiation
it identifies is biopharmaceutical and cIMS coverage, and that is exactly the part
that depends on internal measurements rather than public data.

---

## 8. Immediate decisions needed

1. **Internal research or commercial product.** Decides which data sources are
   usable. Blocking for the benchmark set.
2. **Instrument access for cIMS.** The cIMS layer is a stated differentiator and
   public cIMS data is minimal. Whether any internal cIMS measurements exist or
   are planned changes what that layer can contain by the 25th.

Everything else proceeds without input.
