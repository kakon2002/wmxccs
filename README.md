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

## Status: M0 to M5 complete. The model is fitted and served over HTTP

**142 cross-platform matched ions**, from the steroid interplatform study
(DOI 10.1021/jasms.2c00196). 521 records, 517 of them clear to train, 2 ions paired
across two platforms, 43 across three and 97 across all four - and 93 of those 97 span
four platforms in records that may TRAIN, the other four being held by the shared-peak
check. Converted by
`tools/ingest_steroid.py`; nothing transcribed by hand.

**M4 is built on those ions.** Eighteen strata - six platform pairs by three adducts,
never pooled - each fitted three ways, of which only the nine anchored on stepped-field
DTIMS are ever applied. Of 517 cleared records, **417 are corrected and 402 of those are
served a value**; the other 100 are already on the primary platform, and the 15 between 417
and 402 grade `unsupported`, which means the correction would be an assertion rather than an
interpolation, so the number is withheld and the reason is returned instead. The original
measurement is returned untouched beside every corrected one.

Three things about it are worth knowing before reading any number it produces:

- **The headline correction is chosen, not fixed.** Passing-Bablok where its rank interval
  on the slope excludes 1, the median offset otherwise - 10 strata and 8. All three
  corrections are always reported, because they diverge exactly when a few ions are
  levering the fit. On one clean line, wrecking a single point moves the Deming slope from
  1.04 to over 6 and does not move the robust slope at all.
- **The interval is guaranteed at 80 per cent, not 90.** Grouped leave-one-compound-out
  jackknife+ needs no train/calibrate split, which matters because a split-conformal
  interval on these stratum sizes would be the full observed range. The method proves
  1-2*alpha, and both numbers travel with every interval.
- **Coverage is not evidence.** Measured leave-one-out coverage equals the arithmetically
  pinned value in all nine applied strata, so it could not have come out otherwise. What
  is reported instead is the interval width and the tail ratio.

Every M4 output carries a scope that is DERIVED from the provenance of its data rather
than declared: `ComparisonScope` has no member for interlaboratory reproducibility, so
that claim is not representable, and a harmonized cross section cannot be serialised
without its scope. The maturity stays `provisional` and cannot be otherwise while the
corpus is one study. See `LIMITATIONS.md` sections 7D and 7F.

The two Struwe seed files still pair nothing, for three independent reasons, and
`readiness.py` reports that as a blocker with the reason rather than as a small
number. See `LIMITATIONS.md` section 1.

The other two licence-clear sources are retrieved and characterised, not ingested:

- **CCSbase** - 25,020 records read in full, DT 10,967 / TW 8,529 / TIMS 5,524, with
  SMILES and a per-record calibration method. **Its drift gas is not a column**, so
  every record is unmatchable and it yields no matched ions until the gas is resolved
  from the primary papers it aggregates. It lists **36**, of which **35 contribute any
  records**, and **14 per cent of records (3,572 of 25,020) cite no reference at all** -
  so resolving all 36 papers still leaves a seventh of the database without a stated
  origin to resolve it from.
- **Bush Lab** - eight sheets, 1,804 data rows of protein, peptide and drug-like
  ions: the only real data the biopharmaceutical layer has been offered. Its gas IS
  stated, so it is matchable. Its protein sheets state a charge with no adduct, which
  is the unstated-charge-carrier case CONTEXT.md predicted for exactly this source.

Both are in `LIMITATIONS.md` section 7E, with what each one still needs.

## Running it

Python 3.11 or newer; developed and verified on 3.14. The LIBRARY is platform-neutral;
the commands are not, because a virtual environment puts its interpreter in
`.venv/Scripts` on Windows and `.venv/bin` everywhere else. Every command below is given
for both, and every one names the environment's own interpreter rather than a bare
`python` - a bare `python` after this install runs whichever interpreter is on PATH, which
is usually not the one the package was just installed into.

Windows:

```
python -m venv .venv
.venv/Scripts/python -m pip install -e ".[dev,serve]"
.venv/Scripts/python -m wmxccs                  # serves on http://127.0.0.1:8000
```

macOS and Linux:

```
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[dev,serve]"
.venv/bin/python -m wmxccs                      # serves on http://127.0.0.1:8000
```

Activating the environment first (`.venv\Scripts\Activate.ps1`, or
`source .venv/bin/activate`) lets you write `python` for the rest of the session
instead.

`-m wmxccs` prints what it is serving before it starts - how many corrections the
model holds, how many matched ions are behind it, and the scope caveat - so an operator can
see whether a model was found. Interactive documentation is at `/docs`.

Three endpoints:

| | |
|---|---|
| `GET /health` | version, whether a model is loaded, how many matched ions it rests on, and the model's two sha256 digests |
| `GET /confidence/rules` | the grading scheme as data, so it can be challenged |
| `POST /harmonize` | 200 with harmonized values, or 501 where the model covers nothing |

The model is fitted at startup from the CSV files in `data/seed`, which are committed. An
installation whose seed directory is empty serves a working API whose `/harmonize` answers
501 - which is correct rather than broken, and `/health` says so.

`pip install -e .` alone installs the library without a web server; `[serve]` adds uvicorn
and `[ingest]` adds openpyxl, which only the source-workbook adapters in `tools/` need.

**Deploy from a checkout, not from a built wheel.** The seed CSVs are deliberately not
declared as package data, so a wheel carries no measurements and its `/harmonize` answers
501. That is a licensing decision rather than an oversight: the steroid data is
`academic_only` with an attribution obligation, and bundling it into a redistributable
artefact is a decision nobody has made. An editable install from a clone - which is what
the commands above do - resolves `data/seed` in the source tree and serves a real model.

### Seeing it work end to end

```
.venv/Scripts/python tools/demo_end_to_end.py   # Windows
.venv/bin/python tools/demo_end_to_end.py       # macOS, Linux
```

One real ion from the steroid corpus through every stage in order - ingestion, identity,
matched-ion construction, statistics, the three corrections, the interval, the grade, and
what the result may and may not be called. It asserts nothing; the assertions are in
`tests/`. It exists so the pipeline is legible without reading the code.

### Checking it

Run these from the repository root. Substitute `.venv/bin/python` on macOS and Linux.

```
.venv/Scripts/python -m pytest -q                    # the suite: 3766 tests (2317 wmxccs, 1449 wmxglycan)
.venv/Scripts/python -m tools.mutation --check       # anchors only, about a second
.venv/Scripts/python -m tools.mutation               # the full sweep: 304 mutations (282 wmxccs, 22 wmxglycan)
```

The sweep breaks each guard in the package one at a time and requires a test to notice. A
green suite says the tests ran; the sweep says they would have caught something. A mutation
that survives is a behaviour with no test behind it, and is treated as a failure rather
than as a note.

It never writes the repository: the package is copied to a temporary directory, the
mutation is applied to the copy, and the copy is put first on `PYTHONPATH`.

To rebuild the seed files from the transcriptions:

```
.venv/Scripts/python tools/seed_struwe.py
```

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

`matching.py`'s own tests run against `fixtures.py`, a synthetic corpus, so that each
pairing rule can be exercised on a case built to isolate it - real data does not come with
one clean example of each. It is also exercised against the real steroid corpus through
the ingest, grading and harmonization tests. Everything in `fixtures.py` declares
`SYNTHETIC_FIXTURE` and `assert_quotable` refuses to let a report covering one be
presented as a result.

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

Three endpoints, all three answering. `/harmonize` returned 501 unconditionally until
M4 fitted a model; it now returns values for the measurements the model covers, and
refuses the rest individually rather than as a whole.

| | |
|---|---|
| `GET /health` | version, whether a model is loaded, how many matched ions it rests on, and the model's two sha256 digests |
| `GET /confidence/rules` | the grading scheme as data, so it can be argued with |
| `POST /harmonize` | 200 with harmonized values, or 501 where the model covers nothing in the request |

**The refusals are still the contract, not an error path**, and there are three outcomes
rather than two:

- A HARMONIZED VALUE, where the measurement's platform and calibration group match a
  stratum the model applies. It carries the interval, both coverage figures, the grade,
  the scope caveat and the provenance.
- NO VALUE, WITH THE REASON, per measurement: a platform never compared against the
  primary method, a calibration group with no stratum, a value already on the primary
  platform, or a grade of `unsupported` - meaning the correction would be an assertion
  rather than an interpolation. The original comes back untouched and
  `not_harmonized_because` says which.
- 501 FOR THE WHOLE REQUEST, where nothing in it could be harmonized. A response full of
  absent values is a refusal and should carry the status code of one.

What a caller never gets is a placeholder: no zero, no null standing in for a value, no
interval of infinite width, no grade computed against nothing. An API returning a
plausible value with a caveat in a field would be used and the caveat would not be read.

A malformed body still gets 422, validated before any of this, so a refusal never tells a
caller their measurement was fine when it was not.

A 200 does NOT mean validated. The maturity stamp is `provisional` on every response and
cannot be otherwise while the corpus is one study, and every estimate carries a scope
saying so.

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

## What licences allow

Every record carries a reuse status defaulting to `unverified`, and unverified
never enters a fit. The gate **raises rather than dropping rows**, because a
loader that filters on a boolean loses rows silently and the problem only surfaces
when somebody audits counts.

A licence claim in a data cell is not a licence. Claims are checked against
`sources.py`, which records who read the terms, on what date, and at what URL. Two
of the eleven registered sources are unverified and each says what would settle it: the
METLIN CCS database and the Bayesian harmonization paper. The other nine are read and
recorded - two `open_attribution`, three `academic_only`, two
`non_commercial_no_derivatives` and two `excluded`.
