# Handover

> **Before anything else: this repository holds TWO packages at TWO maturities, and the glycan one
> does not predict cross sections.** `src/wmxccs` is the CCS harmonization pipeline, `0.8.1`, with
> a model fitted on 142 cross-platform matched ions. `src/wmxglycan` ranks candidate glycan
> structures on the evidence that attests them, `0.1.0`, with **no fitted model of any kind** —
> there is no glycan CCS model and V1 will not have one. Nothing in the CCS sections below carries
> across that wall.
>
> **[READ_THIS_FIRST.md](READ_THIS_FIRST.md) is the seven things a reader gets wrong in the first
> hour.** Read it before this document. Four of the seven have been mistaken for bugs already, and
> two of them would change a number you reported to somebody.

For an engineer picking up `wmxccs` cold. Written 21 September 2026 against `v0.6.0-mvp`;
figures and tag reference corrected 23 September 2026 for `v0.6.2-mvp`. The banner above and the
document it points at were written on 27 September 2026 and first ship in **`v0.8.1-mvp`** —
`v0.8.0-mvp` predates them by five commits, which is why that tag was superseded rather than
reused.

This is an orientation document. It does not repeat `LIMITATIONS.md`, which is the long
and authoritative account of what this platform does not do and does not know. What this
does is tell you what the thing is, how to run it, what is actually in it, the two limits
that shape every number it returns, and the order to read things in.

---

## 1. What the platform does

It harmonizes collision cross section (CCS) measurements across ion mobility platforms.

Four platform types measure CCS: drift tube (DTIMS), travelling wave (TWIMS), trapped ion
(TIMS) and cyclic (cIMS). The same ion measured on two of them does not give the same
number, and the difference is systematic rather than random. Published values are therefore
not directly comparable, and people compare them anyway.

This platform:

- **stores measurements without merging them.** A DTIMS value and a TWIMS value of the same
  ion stay two measurements. There is no averaging step anywhere;
- **pairs the same ion across platforms** through a matched-ion key;
- **quantifies the bias** between each pair of platforms, per calibration group, with a
  regression a few bad ions cannot move;
- **returns a harmonized CCS with an uncertainty interval and a confidence grade, BESIDE
  the original** — never instead of it.

The last point is the hard rule of the codebase. Any function that mutates an original
measurement is a bug, and the API echoes every submitted measurement back unchanged and
unrounded, whatever else it does.

### What it is not

**It is not a reproducibility measurement**, and the code will not let a figure from it be
labelled as one. The corpus is one study, so a comparison between platforms here is a
comparison between the instruments that study used — not a statement about how much a
platform type varies between laboratories. That claim is not representable: `ComparisonScope`
has no member for it, and `assert_may_be_quoted_as` raises for every stamp when asked.

**CORRECTED 27 September 2026.** This paragraph said "It is also **not the glycan platform.** That
is a separate repository (`Project2`, package `wmxglycan`)", which was true when it was written and
contradicts the banner at the top of this file. `src/wmxglycan` is in THIS repository, ported beside
the CCS core on 26 September 2026. What survives from the old wording, unchanged: **`Project2` is
not a dependency.** Code was copied and adapted, never imported; it is not on this repository's
path; and `docs/GLYCAN_PORT.md` lists every ported module against the digest of its source. The two
packages never import each other, and `tests/test_glycan_boundary.py` enforces that.

The sentence above this one is still true and is about the CCS core alone: it is not a
reproducibility measurement.

---

## 2. How to run it

Python 3.11 or newer; developed on 3.14. The library is platform-neutral, the commands are
not — a virtual environment puts its interpreter in `.venv/Scripts` on Windows and
`.venv/bin` elsewhere. Always name the environment's interpreter rather than a bare
`python`, which after installing runs whatever is on PATH and generally not this.

**Windows**

```
python -m venv .venv
.venv/Scripts/python -m pip install -e ".[dev,serve]"
.venv/Scripts/python -m wmxccs                   # http://127.0.0.1:8000, docs at /docs
```

**macOS and Linux**

```
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[dev,serve]"
.venv/bin/python -m wmxccs
```

Then, from the repository root:

```
.venv/Scripts/python tools/demo_end_to_end.py    # one ion through every stage
.venv/Scripts/python -m pytest -q                # counts: README.md, section "Checking it"
.venv/Scripts/python -m tools.mutation --check   # anchors only, about a second
.venv/Scripts/python -m tools.mutation           # the full sweep, ~7 h. REFUSES a dirty tree; commit first

# the two services
.venv/Scripts/python -m uvicorn wmxccs.api:app --port 8000          # the CCS core
.venv/Scripts/python -m uvicorn tools.glycan_service:app --port 8010  # the glycan service, WIRED
#   then open http://127.0.0.1:8010/ for the six-page dashboard, served by the service itself
.venv/Scripts/python tools/glycan_service.py                        # what is wired, serves nothing
```

**Deploy from a checkout, not from a built wheel.** The seed CSVs are deliberately not
package data, so a wheel carries no measurements and its `/harmonize` answers 501 for
everything. That is a licensing decision — see §5 below.

### The three endpoints

| | |
|---|---|
| `GET /health` | version, whether a model is loaded, how many matched ions it rests on, and the model's two digests |
| `GET /confidence/rules` | the grading scheme as data, so it can be argued with |
| `POST /harmonize` | 200 with harmonized values, or 501 where the model covers nothing in the request |

`/harmonize` has three outcomes and the distinction matters: a value; no value *with the
reason*, per measurement; or 501 for the whole request when nothing in it could be
harmonized. What a caller never gets is a placeholder — no zero, no null standing in for a
value, no interval of infinite width.

### The model fingerprint

Every estimate and `/health` carry two sha256 digests:

```
corpus      064eb9fba6039147b92c623fafbff342854becbdff60421ea15d2c46ceb83109
parameters  0d69f799f6f1e91cbfce67d9e090b0e7a6656ae852bd0dbd81b26987abb5f3f9
short       064eb9fba603/0d69f799f6f1
```

The model is refitted from the seed files at every startup. Without these, two deployments
could give different numbers for one input and neither response would say so. **Quote them
with any number that leaves this system.**

They are also your regression check. `corpus` covers the records behind the fit; `parameters`
covers everything that decides an answer. If you change something that is supposed to be
presentational and a digest moves, you changed an answer. Verify with `python -m wmxccs`, or:

```python
from wmxccs.api import build_default_model
print(build_default_model().fingerprint.short)
```

---

## 3. The corpus it holds

One study. Everything below comes from
**DOI 10.1021/jasms.2c00196**, a steroid interplatform study, converted by
`tools/ingest_steroid.py`. Nothing was transcribed by hand.

```
rows read          2,075   across four seed files, 0 failed to parse
records built      2,075
records cleared      541   UNCHANGED by the Bush Lab ingest: see below
records held       1,534   1,437 Bush Lab, 89 struwe2015, 4 struwe2016, 4 steroid

held, by reason      876   suspected shared peak
                     552   reuse claim not backed by a licence record
                      89   gas AND uncertainty both unstated (all of struwe2015)
                      12   held for curation review
                       5   conformer set incomplete

matched ions         142   on two or more platforms
                       2   span two platforms
                      47   span three
                      93   span all four

strata                18   fitted
                       9   applied, anchored on stepped-field DTIMS
                           sizes 23, 23, 23, 29, 31, 31, 41, 46, 46
                       9   reported as diagnostics only

corrections graded   417
values served        402   the other 15 grade `unsupported` and are withheld with their reason
not corrected        124   each with its reason returned

grades                 0   supported
                     186   qualified
                     216   weak
                      15   unsupported

scope                517   distinct measurements, entering 701 cross-platform pairings
```

THE CORPUS GREW AND THE FIT DID NOT. 1,437 Bush Lab MicroSource records were ingested on
26 September 2026 as a REFERENCE LIBRARY: stored, complete, and none of them trainable. They
contribute no matched ion and structurally cannot - one paper on one platform, where a matched
ion needs two - so `records cleared`, `matched ions`, every grade figure and the fingerprint are
all unchanged by their arrival. Every one of the 1,437 also fails the licence check
independently: permission for those values comes from the database page's citation request,
and the gate requires a published source to name a DOI with a licence record. See LIMITATIONS 7E.

Two of the four seed files contribute nothing to the fit. `struwe2015_analyst.csv` is held
entirely on licence — 89 records, none cleared. That is the licence gate working, not a bug.

Note the two counts in the last line. **517 is distinct measurements; 701 is pairings.**
A measurement used in three strata is one measurement and three pairings. These were
conflated until 20 September 2026, and the model published 1402 — 517 counted twice over,
summed across strata — in the caveat sentence that travels with every figure.

### Reading the data

`data/seed/*.csv` is the corpus, one row per measured value, committed.
`data/seed/as_delivered/` holds the source sheet written out verbatim, so the conversion is
diffable against its input from a clone that has no `data/raw/`.

---

## 4. The two ceilings

Both are limits of the DATA, not of the code. Neither is fixed by editing a threshold, and
both are the things a new reader most often misreads as defects.

### Ceiling one: confidence stops at `qualified`

**No record in this corpus grades `supported`.** 186 qualified, 216 weak, 15 unsupported, of
417.

That reads like a verdict on the measurements and it is not one. The grading scheme demotes
any calibration group holding fewer than `TARGET_MATCHED_IONS` = 100 matched ions, and every
stratum this model applies holds between **23 and 46**. So the population rule fires on 417
of 417 records — every one — and the top grade is out of reach for reasons that have nothing
to do with the ion being graded.

Two consequences worth holding on to:

- **`qualified` is the ceiling, not a middling result.**
- **`weak` is usually that ceiling minus one notch.** 209 of the 216 weak grades are
  `qualified` demoted once because a rule could not be evaluated. Only 7 are weak for a
  reason about the submitted value itself.

The grade is live code, not dead. On a synthetic stratum of 100 matched ions with one ion
that does not transfer, `grade_correction` returns `supported`. What stands between a real
caller and it is the corpus:

| applied stratum | ion already in the corpus | ion the corpus has never seen |
|---|---|---|
| 23 to 46 ions, as today | `qualified` | `weak` |
| 100 ions or more | **`supported`** | `qualified` |

**What lifts it:** one applied stratum reaching 100 matched ions. That means 54 more ions
measured on *both* platforms of *one* pair, within *one* calibration group and *one* adduct
— not 54 more measurements, which would spread across strata and lift none of them.

**What does not lift it, ever, for the ordinary user:** the right-hand column. An ion this
corpus has never measured cannot have the outlier rule evaluated for it, and an unevaluable
rule always costs one notch — so a genuinely new ion tops out at `qualified` at *any* corpus
size under this scheme. If you are wondering why a caller's grade never improves, this is
why.

**Do not lower the threshold.** It is policy, it is published as policy in
`/confidence/rules`, and it is arguable — but argue it on what a platform-pair figure needs
before being quoted, not on what would make today's corpus look better. A test fails if 100
drifts down to a value the corpus already reaches.

### Ceiling two: scope stops at within-study

Every figure this platform produces is **WITHIN ONE STUDY**, and the maturity stamp is
`provisional` on every response. Neither can be otherwise while the corpus is one study,
because validation means checking against data the model was not fitted on, and one study has
none by definition.

This is enforced structurally rather than by prose. `ScopeStamp` has **no `scope` field** —
scope is a computed property of which studies are behind the figure, so widening the claim
means naming another study, which is data rather than a flag. `assert_may_be_quoted_as`
raises for `INTERLABORATORY_REPRODUCIBILITY` against every stamp, whatever it holds.

Two specifics a reader will otherwise get wrong:

- drift tube against trapped ion is **one laboratory comparing its own two instruments**;
- travelling wave is not that laboratory's measurement at all. It is republished from
  DOI 10.1021/acs.analchem.9b05247 and is itself a four-instrument interlaboratory average.
  So a TWIMS comparison is one laboratory's number against somebody else's average, which is
  a third kind of comparison again.

**What lifts it:** a second study in the corpus, ingested through the licence gate. Section
7E of `LIMITATIONS.md` describes two candidate sources already retrieved and characterised
but not ingested, with what each still needs.

### A third thing, which is a method limit rather than a data one

Intervals are **guaranteed at 80 per cent, not the nominal 90.** Jackknife+ proves
1−2·alpha. Both figures travel with every estimate, adjacent and with the guaranteed one
first, and the smaller is the one to rely on. This does not improve with more data; it is
what the method proves.

---

## 5. Where to start reading

### Read these two first

Before the source, before the README, read these two passages of `LIMITATIONS.md`. They are
where the reasoning of this project is densest, and both describe traps you can otherwise
fall into within a day.

**1. §4.5, instance Ten — "a rule that could not answer, reported as a rule that passed."**

`flagged_as_an_outlier` is a corpus lookup: it asks whether the submitted ion is one the
corpus already recorded as not transferring. For an ion not in the corpus there is nothing to
look in, and it returned `None` — the same value it returns for an ion that *is* present and
is fine. The caller was told a check had passed that had never been made, on every novel ion,
which is the entire population this platform exists to serve.

Read it for the fix, but read it above all for the failure class it heads: **a guard that
looks tested and is not.** Ten instances are catalogued, with the lessons drawn from each.
The one that recurs most: a predicate returning `None` for two different reasons is telling
its caller the stronger one, and a response field meaning "nothing to report" is a *claim*
whenever anything ever populates it.

Note that the tenth instance includes one committed *while fixing the ninth*. That is not
irony, it is the point — the class is easy to re-enter.

**2. §7C — the API and the confidence scheme.**

Everything in ceiling one above, in full, plus the two findings that are recorded as OPEN
rather than fixed:

- the outlier rule cannot answer for an ion the corpus has not measured — **100 per cent of
  the time for the intended user**, not the 30 per cent that replaying our own corpus
  reports, because that replay consists entirely of ions already present;
- **relabelling an ion still extracts a value the service refuses under its real name**,
  5 of 15 such records. Identity is caller-supplied and unverifiable here, and the ruling
  that an unevaluable rule *demotes* rather than refuses is what leaves the gap — refusing
  every new ion would refuse the platform's purpose. What changed is that the response no
  longer claims the check was made.

### Then, in this order

| | |
|---|---|
| `README.md` | what it is, how to run it, the layout |
| `LIMITATIONS.md` §"READ THIS FIRST" | seven things it is easy to read into these numbers that are not there. Written for a non-specialist; read it anyway |
| `CLAUDE.md` | the nine hard constraints. Breaking any of them makes the output unusable |
| `CONTEXT.md` | everything established before this repository existed |
| `PROJECT_REFERENCE.md` | everything established since. Read before re-researching anything |
| `LIMITATIONS.md` §7G | six deliberate omissions, each with what closing it would take |

### Then the source, in dependency order

Each module's docstring says what it is for and, usually, what it refuses to do.

```
identity.py       what was measured, and the key that decides when two measurements are one ion
models.py         the measurement record and its validators
reuse.py          reuse statuses and their tiers
licensing.py      the gate
sources.py        where each source came from, on what terms, and who recorded reading them
loader.py         a file of measurements, read without losing or inventing a row
matching.py       which measurements are of the same ion, and which are not
statistics.py     how two platforms agree, and how they merely correlate
robust.py         Passing-Bablok: a slope a few bad ions cannot move
harmonization.py  the model. A corrected value beside the original, never instead of it
grading.py        confidence: rules, not a fitted model
scope.py          what a figure is entitled to be quoted as
readiness.py      what must be true before anything may be fitted
contracts.py      the request and response bodies
api.py            the three endpoints
```

The two concepts to get right before anything else, both in `identity.py`:

- **the matched-ion key** — analyte identity, adduct, signed charge, the gas the value refers
  to, and structural state. It carries **no platform**, deliberately: that is what lets two
  platforms' values be recognised as the same ion. A compound name alone never matches two
  records;
- **the calibration group** — platform, method, gas, calibrant, adduct, state, and **no
  analyte**. It is not the matched-ion key and confusing the two will produce a pooling bug
  that looks like a statistics bug.

Also: every "the source did not say" unknown keys **uniquely**. An unstated drift gas and an
unstated charge carrier each key to themselves, so two such records never match each other.
This is deliberate and is tested in both directions.

---

## 6. How work is verified here

Two gates, and the second is the one that means something. **The commands are in section 1
and are not repeated here**; the counts are in `README.md` § Checking it and are not repeated
anywhere.

The first gate is the suite. The second is the mutation sweep, which breaks each guard one at a
time and requires a test to notice: a green suite says the tests ran, the sweep says they would
have caught something. A mutation that survives is a behaviour with no test behind it and is
treated as a failure rather than as a note.

**WHY THERE IS NO SECOND COMMAND BLOCK HERE.** There was one, carrying the same two counts as
section 1, and **they drifted three times**: 2,269 against 2,317, then 3,766 against 3,842, then
3,842 against 4,088 - each time correct in one block and wrong in the other, and each time found
after the fact. On 27 September 2026 the duplicate was deleted rather than refreshed a fourth
time, on the owner's instruction: a figure maintained in two places is a figure that will be
wrong in one of them.

The counts now live in one file and `tests/test_repository_counts.py` derives the real numbers
and asserts that file matches them, so the figure cannot go stale without the suite failing. It
is no longer maintained by hand at all, which is the only arrangement that has ever survived a
busy week here.

## 7. Status

M0 to M5 complete and tagged `v0.6.2-mvp`. Deadline: deployable by 25 September 2026, 27 at
the latest — so the remaining time is for whatever the owner decides, not for a rewrite.

Six things are recorded as decided-not-to-build in `LIMITATIONS.md` §7G, each with what
closing it would require. They wait on decisions rather than code: authentication, loopback
binding as the only access control there is, CORS, rate limiting, seed data in a built wheel,
and corpus growth from submitted measurements. The `academic_only` licence bears on several
of them before the engineering does.

If you change anything, **check the fingerprint before and after.** Anything presentational
must leave `064eb9fba603/0d69f799f6f1` untouched. If it moves, you changed an answer.
