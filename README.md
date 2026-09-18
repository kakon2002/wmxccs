# wmxccs

Cross-platform CCS harmonization pipeline. Stores DTIMS, TWIMS, TIMS and cIMS
measurements without merging them, pairs the same ion across platforms,
quantifies inter-platform bias and agreement, and returns a harmonized CCS with
uncertainty and a confidence grade alongside the original values.

It never overwrites an original measurement with a corrected one. Both are
returned.

Separate from the glycan platform. See `CLAUDE.md` for the constraints,
`CONTEXT.md` for everything established before this repository existed, and
`LIMITATIONS.md` for what this does not do and does not know.

**Deadline: deployable by 25 September 2026, 27 at the latest.**

## Status: M0, M1 and M2 complete

The data layer is in place and the seed corpus loads. **No harmonization model
exists, and there is not one cross-platform matched ion in the repository** — the
117 seeded measurements are all travelling-wave from one laboratory. That is the
expected M0 result, not a failure, and `readiness.py` reports it as a blocker with
the reason rather than as a small number. See `LIMITATIONS.md` section 1.

## Layout

```
src/wmxccs/
  reuse.py        reuse statuses and their tiers; sits below everything else
  licensing.py    the default-deny training gate
  sources.py      the licence registry: who read which terms, when, and where
  identity.py     the analyte union, and the matched-ion key
  models.py       CCSMeasurement, its validators, calibration lineage, cyclic IMS
  loader.py       strict CSV loading that loses nothing and invents nothing
  readiness.py    what refuses, what merely warns, and the honest zero
  matching.py     matched-ion construction: which measurements are the same ion
  fixtures.py     a synthetic corpus that cannot be quoted as a result
tools/
  mutation/       the mutation harness: break a guard, require a test to notice
  seed_struwe.py  one-off conversion of the two transcriptions into seed format
tests/
data/
  raw/            gitignored
  seed/           versioned: the converted seed files
  seed/as_delivered/   the transcriptions verbatim, filenames included
```

## The matched-ion key

The unit everything rests on. Two measurements are of the same ion when all five
of these agree:

**analyte identity + adduct + signed charge + the gas the value refers to +
structural state**

It deliberately carries **no platform**. Excluding it is the whole point: a DTIMS
value and a TWIMS value of one ion must land on the same key, or there is nothing
to compare and no bias to measure.

For a protein or an antibody, charge state and native-or-denatured condition are
part of identity rather than metadata: a native 24+ antibody ion and a denatured
40+ ion of the same antibody are two different ions with two different cross
sections.

A compound name is never an identity. Every analyte identifies itself by a
declared identifier — an InChIKey, a sequence, an accession, an INN, a structure,
a canonical composition. `display_name` is in no key.

Do not confuse the matched-ion key with the **calibration group**, which is close
to its complement: the calibration group carries the platform, the method, the
calibrant and the cyclic pass count, and carries no analyte, because it answers
"were these produced the same way", not "are these the same ion". Using either in
place of the other produces a check that can never fire.

The cyclic pass count shows the distinction working. A six-pass and a single-pass
cyclic value of one ion share a matched-ion key, so they can be compared, and sit
in different calibration groups, so they are never averaged together.

## Matched ions

A **matched set** is two or more measurements sharing a matched-ion key from two
or more different platforms. Three-way sets are first class: the likely benchmark
carries DTIMS, TWIMS and TIMS for one compound, and flattening that into three
pairs would count one compound three times.

Five things look like a match and are not one, and each has a rule and a test:
the same measurement republished in two papers; two conformers of one ion; two
ions whose charge carrier was never stated; two values from one platform; and a
set one of whose members nobody may use, which is refused with the member named.

`matching.py` is built and exercised entirely against `fixtures.py`, a synthetic
corpus, because there is no real cross-platform data yet. Everything in it
declares `SYNTHETIC_FIXTURE` and `assert_quotable` refuses to let a report
covering one be presented as a result.

## Running it

```bash
python -m venv .venv
.venv/Scripts/python.exe -m pip install -e ".[dev]"

.venv/Scripts/python.exe -m pytest -q            # the suite
.venv/Scripts/python.exe -m tools.mutation --check   # anchors only, about a second
.venv/Scripts/python.exe -m tools.mutation           # the full sweep
```

The mutation sweep breaks each guard on purpose and requires a test to notice. A
green suite says the tests ran; the sweep says they would have caught something. It
never writes the repository: the package is copied to a temporary directory, the
mutation is applied to the copy, and the copy is put first on `PYTHONPATH`.

To rebuild the seed files from the transcriptions:

```bash
.venv/Scripts/python.exe tools/seed_struwe.py
```

## What licences allow

Every record carries a reuse status defaulting to `unverified`, and unverified
never enters a fit. The gate **raises rather than dropping rows**, because a
loader that filters on a boolean loses rows silently and the problem only surfaces
when somebody audits counts.

A licence claim in a data cell is not a licence. Claims are checked against
`sources.py`, which records who read the terms, on what date, and at what URL. Four
of the eleven registered sources are unverified and each says what would settle it.
