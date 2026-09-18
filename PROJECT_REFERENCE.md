# Project reference

Everything established on this project. Written 18 September 2026. Read this before
re-researching anything.

---

## 1. Project history

**12 to 17 September 2026.** Built a glycan isomer and CCS prediction platform for
Wellmatix. Repo `Project2`, package `wmxglycan`, 1,868 tests, milestones M0 to M3
complete. Candidate generation and ranking worked on real data. The CCS model was
never trained because no licence-clean training data existed at sufficient scale.

**17 September 2026.** The CEO replaced the objective with a cross-platform CCS
harmonization pipeline, to be built as a **separate platform**, explicitly not
integrated into the glycan work. Deadline 25 September, 27 at the latest.

**The glycan repo must not be deleted** until the port is complete. Roughly 40
percent of its code carries over.

> Corrected 18 September, from the repository: 1,868 is the collected test count
> (705 test functions); "1,796" appears nowhere in that repo. The 40 percent figure
> is low as an import fact - the transitive closure of the port set is 15 of the 20
> package modules, four of them on the do-not-port list, so cutting the glycan core
> out is most of the porting work rather than a tidy-up afterwards.

---

## 2. Licence status of every data source checked

This table is the most valuable thing in this document. Each entry cost real time
to establish. Do not re-check a settled one.

### Usable

| Source | Licence | Evidence | Content |
|---|---|---|---|
| glycowork / SugarBase | MIT, Copyright 2021 Daniel Bojar | Package LICENSE file, downloaded and read 17 Sep 2026 | 50,461 glycan structures with GlyTouCan accessions, 15,712 N-linked, 15 curated biosynthetic rules with citations, 355 glycoenzymes |
| Struwe et al. 2016, Chem Commun 52:12353, DOI 10.1039/c6cc06247d | CC BY 3.0 Unported | RSC article page, issue listing and PDF header. Read by Shawon, 15 Sep 2026 | 28 TWIMS CCS values, 4 milk oligosaccharides, 4 adducts. Transcribed to CSV |
| Struwe et al. 2015, Analyst 140:6799, DOI 10.1039/c5an01092f | CC BY 3.0 Unported | PDF footer on every page. Read 17 Sep 2026 | 89 TWIMS CCS values, 7 high-mannose N-glycans, 6 adducts, 4 sample origins. Transcribed to CSV |

### Blocked

| Source | Status | Detail |
|---|---|---|
| CCSbase, ccsbase.net | Academic, non-commercial only | About page states use must be for academic non-commercial purposes; commercial users directed to Libin Xu and UW CoMotion. Read 11 Sep 2026. 25,020 entries, mostly lipids and metabolites |
| Manabe et al. 2022, JASMS | CC BY-NC-ND 4.0 | Both clauses block a company independently. 71 PA-labelled N-glycans |
| Sastre Toraño et al. 2025, Nat Commun, DOI 10.1038/s41467-025-67069-w | CC BY-NC-ND 4.0 | Verbatim: "You do not have permission under this licence to share adapted material derived from this article or parts of it." Read 17 Sep 2026. Held 164 intact-ion CCS values and 394 SNFG structure images |
| AllCCS2, allccs.zhulab.cn | No licence posted | Copyright notice only, which defaults to all rights reserved. Built for small molecules |
| GlycoMob | No terms stated | Over 900 glycan CCS values. Site blocks automated access, a recent review reports it became inaccessible |

### Unverified, and blocking the new project

| Source | Why it matters | What to check |
|---|---|---|
| Steroid interplatform study, JASMS 2022, DOI 10.1021/jasms.2c00196 | **The benchmark dataset.** 87 steroids, 142 CCS values, DTIMS + TWIMS + TIMS. Open access with supporting information free of charge. PMC9545150 | ACS publishes open access under both CC-BY and CC-BY-NC-ND. Read the licence badge on the article page |
| METLIN-CCS, Nature Methods 2023 | 185,589 CCS values, 27,633 standards, TIMS on timsTOF Pro. Largest non-DTIMS resource | Paper says "freely available" at METLIN, XCMS Online and PanoramaWeb. That is access, not a licence. Check each host's terms |
| Bush Lab CCS Database, biophysicalms.org/ccsdatabase | Native and denatured proteins, peptides, small molecules. Most ions measured in **both helium and nitrogen** by similar methods. Downloadable Google Sheet | Page asks only to "cite the appropriate publication(s)", which is a norm not terms |
| Bayesian harmonization study, Anal Chem 2026, DOI 10.1021/acs.analchem.5c06667 | 840 measurements, 347 compounds, three platforms. A ready-made matched-ion set if released | Licence and data availability both unchecked |

### The standard for a licence claim

Verbatim licence text from the publisher's own article page, with the URL and the
date read, recorded by a named person. None of these count: an open access badge,
the journal's usual policy, a PDF copy hosted on ResearchGate or similar, a
Europe PMC open-access flag (observed to be stale for RSC titles), or a download
button with no terms beside it.

---

## 3. Domain facts established

**Platform reproducibility.** Interlaboratory DTIMS: 0.29 percent RSD stepped
field, 0.54 percent mean absolute bias single field (Stow et al. 2017). The steroid
interplatform study found 95 percent of ions within 1 percent bias for TIMS and 2
percent for TWIMS relative to DTIMS, with under 1.5 percent of ions showing biases
up to 7 percent. Correlations from that study: 0.9949 TWIMS-DTIMS, 0.9953
TIMS-DTIMS, 0.9989 TWIMS-TIMS.

**High correlation does not mean agreement.** A platform can correlate at r near 1
while being systematically 2 percent high. Report association and agreement
separately. This is why Deming regression and Bland-Altman are in the plan, not
just Pearson.

**Primary versus calibrated CCS.** Only stepped-field DTIMS gives primary,
first-principles CCS. Single-field DTIMS, TWIMS, TIMS and cIMS all calibrate
against reference values. The calibrant identity changes the answer, and for
glycans the calibrant class matters specifically.

**The gas a value refers to is not always the gas in the cell.** TWIMS measured in
nitrogen but calibrated against helium reference values produces a helium-referenced
CCS. These must be separate fields. Pooling them is an invisible error.

**Adduct choice determines whether isomer information exists at all.** For the
LNH/LNnH pair: [M-H]- differ by 11.4 percent, [M+Na]+ by 3.8 percent, [M+Cl]- by
0.08 percent. Chloride is negative mode and separates nothing. The effect is
specific to deprotonation, which the source paper attributes to charge mobility.

**Deprotonated ions can show multiple conformers of one structure.** Struwe 2015
found two arrival-time peaks for Man5, Man6, Man9 and Man9Glc as [M-H]-, confirmed
by MS/MS to be the same structure folded differently. Two CCS values for one
matched-ion key is legitimate when the conformer index differs. Do not average,
do not drop one. MS/MS is needed alongside [M-H]- to distinguish real isomers from
conformers.

**Isomer CCS differences are often under 1 percent**, which is at or below
measurement reproducibility. Where the difference is smaller than the error, no
model resolves it.

---

## 4. Design decisions and why

**Split on composition, not structure.** Established on the glycan project and it
overturned the original instruction. A prediction request carries no structure, so
two isomers of one composition are identical to the model. Any structure-derived
split key is finer than the model's own identity and therefore manufactures the
leak it appears to prevent.

**Licence status is a required field defaulting to unverified, and the gate raises
rather than drops.** A silent drop hides a licence problem until someone audits
row counts. The gate also re-validates the record, because `model_copy(update=...)`
bypasses validators.

**A licence claim in a data cell is not a licence.** Claims are checked against a
registry that records who read the terms and when. Unbacked claims are refused.

**Uncertainty carries its type.** SD, 2SD, SEM, CI95 or unknown. Loading 2SD into
an SD field corrupts every interval with nothing able to detect it. `unknown` may
stand alone without a number, since it asserts ignorance rather than describing one.

**Placeholder strings are rejected.** Blank-cell fillers like `n.d.`, `not reported`,
`#N/A`, `<NA>`, `NaT` and bare dashes were being accepted as real calibrant names.
Enum values must not be spelt like something a loader would write, which is why
`none` became `underivatised`.

**Refuse on inevaluability, warn on weakness.** A fit on 40 records proceeds and
reports its own weakness. A fit that cannot be evaluated at all is refused. The
binding constraint on conformal intervals is the calibration set size, not the
training set: a 90 percent interval needs at least 9 calibration points, 95 percent
needs 19.

**Any flag making a biological claim must derive from graph structure, never string
matching on a linkage name.** Three separate bugs of this shape reached review,
including an antennary fucose scoring as core fucose, and a chitinase credited with
forming a bond it degrades.

**Results carry their evaluation mode and cannot be misreported.** A diagnostic run
must be structurally incapable of being quoted as deployment performance.

---

## 5. What ports from the glycan repo

Port by copying and adapting, never importing.

Port: `licensing.py`, the source registry, `CCSMeasurement` with the analyte
generalised, loader hardening, adduct normalisation, conformer sibling checks,
readiness thresholds, and `tools/mutation` with its catalogue rewritten.

Do not port: `composition.py`, `enumeration.py`, `glycan_graph.py`, `enzymes.py`,
the SugarBase loader, the 44-column glycan featuriser, biosynthetic rules.

**The one schema change that matters:** the analyte becomes a tagged union covering
small molecules, peptides, glycans, proteins and subunits, intact antibodies and
ADCs. For proteins and antibodies, charge state and native/denatured condition are
part of identity, not metadata. A native 24+ ion and a denatured 40+ ion of the
same antibody are different matched ions.

---

## 6. Open questions

| Question | Blocks | Who answers |
|---|---|---|
| Is this platform internal research or commercial? | Which data sources are usable, and therefore whether a benchmark set exists | Kang |
| Which licence does the steroid study carry? | The benchmark dataset | Shawon, one page read |
| Is there internal cIMS instrument access? | The cIMS layer, a stated differentiator with almost no public data | Kang |
| Which gas were the dextran references in, for Struwe 2015? | One of TWO blockers on those 89 records; resolving it alone unblocks nothing. Check Hofmann 2014, Anal Chem 86:10789 | Shawon |

---

## 7. Carried over from the glycan project, still unresolved

These do not block the new platform but were never closed:

- Figure 1 of Analyst 2015 holds the seven high-mannose structures needed to
  unblock those 89 records.
- **Corrected 18 September: the drift gas is only HALF of what blocks those 89
  rows.** Every one of them also carries `uncertainty_type = unknown` with no
  spread, which is an independent training blocker. Reading Hofmann 2014 and
  resolving the gas would unblock zero records on its own. Both must be resolved,
  and there is now a test pinning both so that whoever resolves the gas is told
  about the second one rather than discovering it afterwards.
- Four rows in the 2016 dataset are held pending curation: human milk LNH and LNnH
  share identical values for [M+H]+ (228.9) and [M+Cl]- (245.0), almost certainly
  one unresolved peak reported against both compounds.
- **Corrected 18 September: those four rows are not marked as held anywhere in the
  CSV.** The hold is DERIVED by the loader's shared-peak detector, from exact CCS
  equality within one calibration group across differing analyte identity. It is a
  result, not data. A port that dropped the detector would clear all 28 rows and
  report a clean run.
- A curated glycosyltransferase list with citations was needed before the ranking
  layer shipped. glycowork carries no EC numbers or CAZy references, so it cannot
  be derived from shipped data.

---

## 8. Working notes

Claude Code has repeatedly caught real defects that review alone missed, several
of them in instructions given to it. When it pushes back, check before overriding.
Examples: the composition split key, the conformer detector limitation, the
`declares_synthetic` hole, and refusing to reconstruct data from a description.

The mutation harness is the reason those guards are known to bite rather than
merely pass. It lives in `tools/mutation` and must be ported.

Two failure classes have recurred and are worth watching for: **something that
reports success while doing less than it claims** (a swallowed CSV row, a
miscounted source, a stale mutation anchor), and **a plausible biological or
numerical claim that was never actually read from a source**.
