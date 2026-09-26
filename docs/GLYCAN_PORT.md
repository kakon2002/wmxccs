# The glycan port

What came across from `Project2` into `src/wmxglycan`, what did not, and what had to
change. Written on 26 September 2026 against commit `ccf0db1c` (tag `v0.7.0-mvp`).

Read this beside `CLAUDE.md` § "Two packages, and the wall between them", which states
the rules this document reports compliance with.

---

## 1. The P0 item first: the CCS core did not move

The owner's first instruction was to protect what exists, record the model fingerprint,
and confirm it at every stop.

```
recorded before any glycan work, and again after the whole port:

  corpus     064eb9fba6039147b92c623fafbff342854becbdff60421ea15d2c46ceb83109
  parameters 0d69f799f6f1e91cbfce67d9e090b0e7a6656ae852bd0dbd81b26987abb5f3f9
  short      064eb9fba603/0d69f799f6f1        UNMOVED
  matched ions 142   strata 18 fitted, 9 applied
```

Unmoved since `v0.5.0-provisional`, and unmoved across this port. Three independent
checks back that up, because a fingerprint is a hash of inputs and a hash agreeing is
weaker evidence than it looks:

1. **`git status src/wmxccs/` reports nothing.** Not one byte of the CCS package
   changed. The port added a sibling directory and never reached into this one.
2. **The CCS suite still collects exactly 2,317 tests and they all pass.** Same
   number as at `v0.7.0-mvp`, measured, not carried forward.
3. **`tests/test_glycan_boundary.py` proves the separation rather than asserting it.**
   Importing `wmxccs.api` and fitting the served model loads no `wmxglycan`, and no
   `networkx`, `glycowork`, `sklearn`, `numpy` or `pandas` either — which is a stronger
   claim than the port needed, and is the one the `pyproject` makes in prose about the
   CCS statistics being pure Python on purpose.

The only pre-existing files this port touched are `CLAUDE.md` and `pyproject.toml`, both
because they asserted the opposite of the new brief. Both corrections keep the old
wording and say when it stopped being true.

## 2. What ported

Fifteen modules. The five the owner named need ten more to import at all, so the whole
dependency closure came across; nothing outside that closure did.

**Fourteen of the fifteen are still byte-identical to the file they came from.** One is not,
and that is what this table is for: the digest is of the ORIGINAL in `Project2`, so a
deliberate edit here can be told from a bad copy. Updated 26 September 2026, when building the
ranker turned up a defect in `enumeration.py`; what changed and why is in § 4A of
[`GLYCAN_LIMITATIONS.md`](GLYCAN_LIMITATIONS.md) and summarised below the table.

| module | lines | sha256 (first 16) of the original | state |
|---|---|---|---|
| `composition.py` | 310 | `80cf439cbfcd5434` | identical |
| `constraints.py` | 148 | `f739d4f5a384b643` | identical |
| `contracts.py` | 92 | `37116dc88cc1886e` | identical |
| `enumeration.py` | 630 | `291e61875a0fa0d7` | **ADAPTED, see below** |
| `enzymes.py` | 156 | `25b56780dc437d9d` | identical |
| `evaluation.py` | 344 | `5c6a676b61513f6a` | identical |
| `features.py` | 529 | `69b529363353098e` | identical |
| `glycan_graph.py` | 246 | `2ecfa140bc2332af` | identical |
| `licensing.py` | 262 | `974cec4c1d27e0e1` | identical |
| `models.py` | 785 | `a9e63cc8abc01ac9` | identical |
| `reuse.py` | 118 | `01f56c3e52667594` | identical |
| `sources.py` | 320 | `fb8bd551c929d1b5` | identical |
| `splits.py` | 709 | `131c7d8e63afc864` | identical |
| `sugarbase.py` | 243 | `f696929dc7a9a715` | identical |
| `training.py` | 412 | `8569792ac44b15e8` | identical |

`__init__.py` is new. `attestation.py`, `ranking.py` and `ccs_evidence.py` were written here
and are not ported at all; `tools/glycan_ccs_evidence.py` is the adapter between the two
packages and is the only file in the repository that imports both.

### The one adapted module: `enumeration.py`

Two changes, both confined to the RATIONALE a candidate carries. **No candidate is included or
excluded by either, and the counts are unmoved at 10 / 63 / 167 / 6**, which is asserted in the
tests rather than claimed here.

1. **The ordering caveat is gated on its context.** It used to attach whenever a curated rule's
   product was present, while its own wording claims the candidate "carries {product} alongside
   that context". Measured on `Hex5HexNAc4Fuc1`: 114 of 167 candidates carried the note and only
   18 contained the bisecting GlcNAc the note names, so **108 published a cited sentence that was
   false about them.** After the fix there are 6 notes, all true.
2. **`Enumerator.broken_rules()` was exposed**, a one-line public wrapper, so the ranker can
   verify that a candidate breaks no rule by running the check instead of trusting that the
   enumerator rejected the violators. Those are different claims, and "zero" is also what an
   unrun check returns.

Mapped against the five items in the brief:

| the brief asked for | module | what carries it |
|---|---|---|
| 1. composition parsing and canonicalisation | `composition.py` | `parse_composition`, `Composition.canonical`. Residue order is the enum's definition order; `dHex` reads as `Fuc`; a string that does not parse completely is refused rather than partly parsed. |
| 2. the candidate generator | `enumeration.py` + `constraints.py` + `enzymes.py` + `glycan_graph.py` + `sugarbase.py` | `Enumerator.enumerate()`. Curated rules from the glycowork tables, deduplicated on `canonical_key`. |
| 3. the feature extractor | `features.py` | 44 columns in three blocks: 9 analyte, 11 condition, 24 structure. |
| 4. grouped splits keyed on composition | `splits.py` | `analyte_key` defaults to `AnalyteLevel.COMPOSITION` and reads the composition and nothing else. |
| 5. the evaluation harness and the refuse-to-fit guard | `evaluation.py` + `training.py` | `m3_report`, the maturity stamp, `TrainingSet.readiness`, and `fit_ccs_baseline`, which refuses. |

### Unknown linkage and anomericity stay unknown

The brief singled this out. `glycan_graph.py` keeps what an IUPAC-condensed string left
open as open, `features.py` counts it (`anomer_unknown`, `position_unknown`,
`sia_position_ambiguous` are three of the 44 columns), and nothing completes it. The
tests for it ported unchanged and pass: `test_glycan_glycan_graph.py::test_what_the_string_leaves_open_stays_open`
and `test_glycan_sugarbase.py::test_the_shipped_records_keep_what_the_strings_left_open`.

## 3. What did not port, and why

- **The ranker. It does not exist.** The brief said not to port it; there was nothing to
  refrain from. `Project2`'s M5 ("Ranking. Score candidates using CCS interval
  consistency plus...") was never built, and no function in its package ranks or scores
  anything — checked, not assumed. This is recorded because "the ranker was not ported"
  and "there was no ranker to port" are different facts and only the second one is true.

  **A ranker was subsequently BUILT here**, on 26 September 2026, on the owner's decision that
  V1 has no glycan CCS model: `src/wmxglycan/ranking.py`, with `attestation.py` and
  `ccs_evidence.py` beside it. It is new code and not a port, and it shares no design with
  anything in `Project2`, because there was nothing there to share a design with.
- **`api.py`** — the glycan HTTP surface. Deliberately left. The brief puts the glycan
  layer *on top of the existing CCS service*, so what the endpoints are and what crosses
  between the two packages is the boundary decision, not a copy. Copying `Project2`'s app
  would have decided it by default.
- **`measurements.py`** — the CCS measurement loader. `wmxccs.loader` already holds that
  role for this repository, over a corpus of 2,075 records. A second loader for the same
  kind of data is how two record definitions drift apart.
- **`struwe2015.py`, `struwe2016.py`** — source adapters for the two Struwe papers.
  `wmxccs` already ingests both (`tools/seed_struwe.py`, and the 89 held struwe2015
  records in the corpus). Porting these would give the repository two readers for one
  paper.
- **`Project2`'s mutation catalogue and runner.** Out of scope for this task. See the
  gap in § 6.

## 4. What resisted

**One test failed on the straight copy, out of 1,165.** That is the honest headline:
the code itself resisted nothing, and no source module needed a single edit.

```
1 failed, 1164 passed in 40.85s
FAILED test_training.py::test_no_estimator_library_is_imported_by_the_package
  ModuleNotFoundError: No module named 'wmxglycan.api'
```

`PACKAGE_MODULES` was a **hand-written list of eleven module names**, and it named
`wmxglycan.api`, which this repository deliberately does not carry. So a guard about
estimator libraries failed for a reason with nothing to do with estimator libraries.

This is the recurring failure class from `LIMITATIONS.md` § 4.5 — a guard that looks
tested and is not — arriving from the other repository intact. The list would have gone
stale the first time anyone added a module and forgot it, and the guard would then have
covered less than it appeared to while staying green. **Adapted rather than patched:**
the list is now derived from the package directory, with a floor under it
(`test_the_derived_module_list_actually_covers_the_package`) so that a glob returning
nothing cannot make the guard vacuous — which would be the same defect one level up.

Two test files were adapted, and only these two things changed:

| file | change | why |
|---|---|---|
| `tests/test_glycan_training.py` | `PACKAGE_MODULES` derived from the package, plus a test that the derived list is not empty | as above |
| `tests/test_glycan_splits.py` | a subprocess's `from tests.test_splits import many` follows the renamed file | the 14 ported test files are prefixed `test_glycan_` so no basename collides with a CCS test; three otherwise would (`test_models`, `test_licensing`, `test_sources`) |

Worth recording as a non-event: the copy ran under **networkx 3.7 and pandas 3.0.6**,
where `Project2` runs 3.6.1 and 3.0.5. No test failed because of the drift.

## 5. The numbers this port had to reproduce

Measured from the copy, in this repository:

| check | expected | measured |
|---|---|---|
| `Hex5HexNAc4Fuc1` candidates | 167 | **167** |
| glycowork rows in the shipped SugarBase release | 50,461 | **50,461** |
| feature columns | 44 | **44** (9 + 11 + 24) |

Candidate counts for the four compositions asked for, with the enumerator's own
accounting — built, rejected by rule, duplicates removed:

| composition | candidates | trees built | rejected | duplicates |
|---|---|---|---|---|
| `Hex3HexNAc4Fuc1` | 10 | 34 | 24 | 0 |
| `Hex4HexNAc4Fuc1` | 63 | 266 | 203 | 0 |
| `Hex5HexNAc4Fuc1` | **167** | 1,012 | 845 | 0 |
| `Hex5HexNAc2` | 6 | 6 | 0 | 0 |

`duplicates removed: 0` in all four is not the deduplication failing to run. Candidates
are accumulated into a dict keyed on `canonical_key`, so a duplicate never becomes a
second entry; the count confirms the generator does not build the same tree twice by two
routes, and the distinct-key count equals the candidate count in every case.

`Hex5HexNAc2` is Man5, and 6 is the number of distinct arrangements of five mannoses on
the core — no rule rejects any of them, which is why its rejected column is empty.

## 6. Where the glycan layer stands, stated plainly

- **Nothing is fitted and nothing is predicted.** `fit_ccs_baseline` refuses on every
  path: not a training set, readiness refused, unknown baseline, synthetic records in the
  fit — and then `NotImplementedError("no model has ever been trained")`. The
  `from sklearn` import sits below every one of those guards, and there is a test that it
  stays there.
- **The ranker exists, and it ranks on evidence rather than on a prediction.** Built
  26 September 2026 in `src/wmxglycan/ranking.py`, on the owner's decision that V1 has no
  glycan CCS model. It bands indistinguishable classes of candidates by how many distinct
  fully-resolved reference structures attest them, refuses to order a set nothing separates,
  and reports CCS as measured evidence or a stated absence and never as a prediction. See
  `docs/GLYCAN_LIMITATIONS.md` §§ 3-5.
- **The boundary schema now exists and carries traffic.** It did not when this document was
  first written, and the sentence here said so. `wmxglycan.ccs_evidence` declares the schema —
  a `CCSEvidence` model with five states and a `CCSEvidenceLookup` Protocol — and
  `tools/glycan_ccs_evidence.py` is the adapter, in neither package, the only file that imports
  both. `tests/test_glycan_boundary.py` still holds the wall up in both directions, and a test
  asserts the Protocol is satisfiable by a stub defined inside the glycan layer's own tests, so
  the declaration cannot quietly come to depend on the CCS core.
- **The code written here has mutation coverage; the ported code deliberately does not.**
  `tools/mutation` now shadows BOTH packages on every run and a mutation names the package its
  file sits in. The catalogue is two modules joined by a derived function: 282 wmxccs entries,
  closed, and 22 glycan entries covering `ranking.py`, `attestation.py`, `ccs_evidence.py` and
  the one changed line in `enumeration.py`. The 1,220 ported tests are NOT re-covered, on the
  owner's instruction, and § 4A of `GLYCAN_LIMITATIONS.md` is what that costs: three defects in
  ported code, none of which failed a test.
- **`glycowork` is pinned exactly, at 1.10.0.** The 50,461-structure SugarBase v12
  release and the curated enzyme and constraint tables ship inside the package, so the
  version is part of the dataset version string (`SugarBase v12 via glycowork 1.10.0`)
  that the licence record quotes. A version range here would let the data under a cited
  figure change without the citation changing.
