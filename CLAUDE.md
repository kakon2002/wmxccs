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

## What this is not

Not the glycan platform. That repo (`Project2`, package `wmxglycan`) is a separate
project and stays separate. Code is ported from it where noted in CONTEXT.md,
never imported from it, and nothing glycan-specific comes across.

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
