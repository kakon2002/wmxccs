# CLAUDE.md

Project instructions for Claude Code. Read this before every task.

## What this is

A cross-platform CCS harmonization pipeline. It stores collision cross section
measurements from DTIMS, TWIMS, TIMS and cIMS without merging them, pairs the same
ion across platforms, quantifies inter-platform bias and agreement, and returns a
harmonized CCS estimate with an uncertainty interval and a confidence grade.

It never overwrites an original measurement with a corrected one. Both are returned.

Deadline: deployable version by 25 September 2026, 27 September at the latest.
Daily report to the owner every evening.

## Two packages, and the wall between them

`src/wmxccs` is the CCS core described above. `src/wmxglycan` is the glycan and
isomer layer, ported beside it on 26 September 2026 on the owner's brief. One
repository, two packages, and they share no code.

**This section said the opposite until 26 September 2026** and the old wording is
worth keeping, because it was correct for every milestone up to v0.7.0: "Not the
glycan platform. That repo (`Project2`, package `wmxglycan`) is a separate project
and stays separate... nothing glycan-specific comes across." The brief of
26 September supersedes it. What survives from it, unchanged and load-bearing:

- **`Project2` is still not a dependency.** Code was COPIED AND ADAPTED from it,
  never imported. It is not on this repository's path and a checkout elsewhere
  will not have it. Every ported file is listed in `docs/GLYCAN_PORT.md` against
  the digest of the file it came from.
- **The two packages never import each other.** Not in either direction. They meet
  through an explicit schema and nowhere else. This is enforced by
  `tests/test_glycan_boundary.py`, not by this paragraph.
- **The CCS core does not change because the glycan layer exists.** It is the
  owner's P0 item. `networkx`, `glycowork` and `scikit-learn` are declared under
  the `glycan` extra, so a CCS-only install is what it always was, and the model
  fingerprint `064eb9fba603/0d69f799f6f1` is the evidence: record it before
  touching anything and confirm it at every stop. If it moves, something was
  touched that should not have been.

The glycan layer holds its own `IMSType`, `DriftGas`, `Polarity` and
`CCSMeasurement`. That is deliberate duplication, not an oversight: one shared
enum would make every change to the CCS core a silent change to the glycan
layer's accepted inputs, where two definitions with a schema between them turn
the same disagreement into a validation error somebody can read.

## Hard constraints

These are not style. Breaking any of them makes the output unusable.

1. **Never fabricate.** No invented CCS values, citations, licence statuses, record
   counts or benchmark numbers. Unknown is written as UNKNOWN with what would
   resolve it. A plausible placeholder is worse than a blank.

2. **Licence gate, default deny.** Every record carries a reuse status defaulting
   to unverified. Unverified never enters a fit. The gate raises, it does not drop
   rows silently. A licence claim is only valid when backed by a registry entry
   that names who read the terms and when.

3. **Original values are never replaced.** A harmonized CCS is a separate output
   alongside the original. Any function that mutates an original measurement is
   a bug.

4. **No pooling across platforms or gases without a model.** Values from different
   platforms, or measured against different reference gases, never enter one
   average or one regression as if interchangeable. They are paired through the
   matched-ion key and related through a fitted model, never assumed equal.

5. **Matched-ion key is the comparison unit.** Molecule + adduct + charge + gas +
   structural state. Compound name alone never matches two records.

6. **Grouped splits before any error number.** The same analyte must not appear
   in both train and test. Report the split strategy beside every metric.

7. **Every uncertainty carries its type.** SD, 2SD, SEM, CI95 or unknown. A
   number without a type is refused.

8. **Refuse on inevaluability, warn on weakness.** A fit on a small set proceeds
   and reports its own weakness. A fit that cannot be evaluated is refused.

9. **Tests with every module.** Nothing is done until tests pass, and the tests
   must assert the constraints above, not only the happy path. Break each guard
   deliberately and confirm the suite notices.

## Data reality

Read CONTEXT.md before touching data. Three of the sources named in the objective
document have unverified or blocking licences. Do not ingest from any source whose
licence is not recorded in the registry with a reader and a date.

## Reporting

After each milestone: what was built, the numbers, what surprised you, what the
next step depends on. Flag any point where evidence contradicts an assumption in
CONTEXT.md rather than working around it.
