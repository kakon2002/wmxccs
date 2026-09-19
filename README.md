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

## Status: M0 to M3 complete, the API contract and confidence scheme done, and the first real data ingested

**142 cross-platform matched ions**, from the steroid interplatform study
(DOI 10.1021/jasms.2c00196). 521 records, 517 of them clear to train, 2 ions paired
across two platforms, 43 across three and 97 across all four. Converted by
`tools/ingest_steroid.py`; nothing transcribed by hand.

That is the first result here that describes real instruments, and it is what M4 was
waiting for. **No harmonization model exists yet**, and every readiness report still
stamps `provisional`. What the data does and does not support - uneven platform
coverage, 142 trapped-ion values held by a calibrant the paper never names, compounds
identified only within their own dataset - is in `LIMITATIONS.md` section 7D.

The two Struwe seed files still pair nothing, for three independent reasons, and
`readiness.py` reports that as a blocker with the reason rather than as a small
number. See `LIMITATIONS.md` section 1.

The other two licence-clear sources are retrieved and characterised, not ingested:

- **CCSbase** - 25,020 records read in full, DT 10,967 / TW 8,529 / TIMS 5,524, with
  SMILES and a per-record calibration method. **Its drift gas is not a column**, so
  every record is unmatchable and it yields no matched ions until the gas is resolved
  from the 36 primary papers it aggregates.
- **Bush Lab** - eight sheets, about 6,500 rows of protein, peptide and drug-like
  ions: the only real data the biopharmaceutical layer has been offered. Its gas IS
  stated, so it is matchable. Its protein sheets state a charge with no adduct, which
  is the unstated-charge-carrier case CONTEXT.md predicted for exactly this source.

Both are in `LIMITATIONS.md` section 7E, with what each one still needs.

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
  statistics.py   association and agreement, kept apart and never pooled
  grading.py      confidence as rules, not a fitted model
  contracts.py    the request and response bodies
  api.py          health, the confidence scheme, and a harmonize that refuses
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

## Association is not agreement

`statistics.py` reports the two separately and never combines them, because they
answer different questions and only one of them is the question this platform
exists to answer.

**Association** (Pearson r, r squared) is whether two platforms move together.
**Agreement** (Deming regression, Bland-Altman bias and limits, Lin's concordance,
delta CCS per cent, MAE, MAPE, RMSE) is whether they give the same number. A
platform can correlate at 0.998 and run two per cent high on every single ion.

There is no ordinary least squares anywhere, not behind a flag. Both axes carry
measurement error, so an OLS slope is biased towards zero by however noisy the x
platform happens to be, and it gives two different answers depending which way
round the platforms are put. Deming is symmetric and is the only regression here.

Nothing is pooled across gases, calibrants or cyclic pass counts: figures are per
platform pair **and** per calibration-group pair within it, and a pair holding more
than one stratum gets no pooled figure at all.

Outliers are flagged relative to their own stratum's offset and **never removed**.
There is no parameter that would remove them. Those ions are the result: they are
where a harmonization model will be confidently wrong.

## The API, and what it refuses

Three endpoints. Two answer today; one will answer 501 until a harmonization model
has been fitted. Real cross-platform matched ions now exist - 142 of them - so what
stands between `/harmonize` and a number is the fitting itself, no longer the data.

| | |
|---|---|
| `GET /health` | version, whether a model is loaded, how many matched ions it rests on |
| `GET /confidence/rules` | the grading scheme as data, so it can be argued with |
| `POST /harmonize` | **501**, with your measurements returned untouched |

The 501 is the contract, not an error path. A caller gets back every measurement
exactly as sent, the provenance and licence terms of each, the maturity stamp, and
a statement of why there is no harmonized value. What they do not get is a number:
no zero, no null standing in for one, no interval of infinite width, no grade
computed against nothing. An API returning a plausible value with a caveat in a
field would be used and the caveat would not be read. One returning 501 cannot be
used by accident.

A malformed body still gets 422. The request is validated before the refusal, so a
501 never tells a caller their measurement was fine when it was not.

The eventual 200 shape is fully specified in `contracts.HarmonizeResponse` so
callers can build against it now. No endpoint returns it, and `/harmonize` declares
`HarmonizationUnavailable` instead, so nobody can write code against a harmonized
value that is never there.

## Confidence is graded by rules, not by a model

A grade is not a probability and has never been calibrated, because calibrating it
would need the held-out matched ions this repository does not have. It is a set of
checks, each evaluable today, each naming a specific reason a correction might not
apply to a particular ion.

A grade **only ever falls**, and the final grade is the worst demotion found: never
a score, never an average, never a count. Two mild concerns do not add up to a
severe one, and a severe one is not offset by everything else being fine. A grade
falls when the ion sits outside the calibration range, its calibration group is
thinly populated, it was flagged as an outlier in its own stratum, a source behind
the correction may not be used, or the correction is substantially driven by ions
that do not transfer.

That last rule is the M3 leverage finding made operational: it refits the slope
without the stratum's outliers purely to **measure** how much they move the
correction, and demotes when that exceeds half a per cent of CCS at a typical ion.
Nothing is removed from the statistics, the correction or any count. Measuring the
influence of outliers is the opposite of dropping them, and a test holds the two
apart.

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
