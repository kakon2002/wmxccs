> **HISTORICAL DOCUMENT.** This records what was asked in September 2026 and predates
> the glycan layer added on 26 September. It is kept as a record and is NOT a description
> of what this repository contains now; for that, read `README.md`.

# Kickoff prompt

Paste this as the first message in a fresh Claude Code session opened in a new
empty folder. CLAUDE.md and CONTEXT.md should already be in that folder.

---

Read CLAUDE.md and CONTEXT.md in full before doing anything.

This is a new project: a cross-platform CCS harmonization pipeline. It is a separate
platform from the glycan work in C:\Users\User\Project2. That repo stays untouched
and is never imported. You port code from it by copying and adapting, as CONTEXT.md
specifies, and nothing glycan-specific comes across.

Deadline is a deployable version by 25 September, 27 at the latest. Today is the 18th.
Work in the milestones below, stop after each, report with numbers.

## M0. Repo and ported data layer (today)

1. Create the repo. Package name `wmxccs`. Layout:

   ```
   src/wmxccs/
     licensing.py      ported from Project2
     sources.py        ported, registry seeded from CONTEXT.md
     models.py         ported CCSMeasurement, analyte generalised per CONTEXT.md
     identity.py       new: Analyte union, matched-ion key
     loader.py         ported hardening, adapted to new columns
     readiness.py      ported, rethresholded for matched-ion counts
   tools/mutation/     ported whole, catalogue emptied, to be refilled per module
   tests/
   data/raw/           gitignored
   data/seed/          versioned, curated transcriptions
   CLAUDE.md, CONTEXT.md, LIMITATIONS.md, README.md, pyproject.toml
   ```

2. Port exactly what CONTEXT.md lists. For each ported module, run its original
   tests against the port and confirm they pass before adapting. Then adapt.

3. Generalise the analyte. `CCSMeasurement.analyte` becomes a tagged union:
   small molecule (InChIKey required), peptide, glycan (reuse the existing
   representation), protein or subunit, intact antibody, ADC. Define the
   matched-ion key exactly as CONTEXT.md states. Write tests that assert a native
   24+ antibody ion and a denatured 40+ ion of the same antibody produce different
   keys, and that compound name alone never produces a match.

4. Seed the registry with every source in CONTEXT.md, at the status stated there,
   with the reader and date given. Unverified sources get entries too, marked
   unverified, so their status is visible and nobody re-checks them.

5. Load the two transcribed Struwe files from Project2/data/raw into data/seed.
   They should load and hold. Report the counts. Expect zero matched-ion pairs
   across platforms, since both are TWIMS from one lab.

6. Port the mutation harness. Refill the catalogue for the ported modules only.
   Run it. Report killed, survived, stale.

Stop after M0. Report: tests passing, records held, matched pairs (expect zero),
mutation result, and anything in Project2 that resisted porting.

## What comes after M0, for planning only

- M1: cIMS and calibration-lineage fields, the pass-number and reference-source
  metadata from CONTEXT.md.
- M2: matched-ion construction across platforms. Pairing, dedup by DOI, provenance
  graph. This is the module everything downstream depends on.
- M3: cross-platform statistics. Pearson, R², Deming regression, slope and
  intercept, mean ΔCCS%, MAE/MAPE/RMSE, Bland-Altman, concordance correlation,
  residual analysis. All on matched pairs only.
- M4: harmonization model. Deming per platform pair, class-stratified variant,
  grouped cross-validation, prediction intervals. Never overwrites originals.
  Leave a seam for the Bayesian model.
- M5: confidence grading, API, and a minimal dashboard.
- M6: packaging and deployment, report.

Do not start M1 until M0 is reported. The data question in CONTEXT.md is open and
may change M2 entirely.
