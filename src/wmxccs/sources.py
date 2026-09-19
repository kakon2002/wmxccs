"""Where each source came from, on what terms, and who reported checking.

A row in a measurements file may claim any reuse status it likes. The claim is
only as good as the evidence behind it, and this registry is that evidence. Two
kinds of source are recorded, because a value and the analyte identity it was
measured on can come from different places:

- a PUBLICATION OR RESOURCE, keyed by DOI where it has one and by URL where it
  does not, for values transcribed from a paper or drawn from a database;
- a DATASET, keyed by the exact provenance string its own loader writes into an
  analyte's source field, for identities drawn from a shipped or seeded corpus.

The licence gate checks every trainable claim against both, so "open_attribution"
on a record traces back to something a person read rather than something a
person typed, whether the record came from a file or was built in code. The one
trainable status that is not a claim about a source is synthetic_fixture, a
test's declaration that the record is invented; it needs no record, and a record
citing a DOI that IS on record is that paper's record and may not call itself
synthetic.

Nothing here is inferred by this code. A publication's licence is recorded after
a person has read it on the publisher's page and reported what they read; the
record carries that report, attributed, and does not claim to quote the page
verbatim unless it says so. A third party's metadata flag is not evidence either
way. The registries are read-only mappings, and an entry is added by editing this
module, with its evidence, and by no other route.

UNVERIFIED SOURCES GET ENTRIES TOO
----------------------------------
That is the point of the `what_to_check` field, and an unverified entry is
refused without one. An unverified source with no record is indistinguishable
from a source nobody thought of, so somebody re-reads the same page next week
and spends the same day on it. An entry saying "nobody has read this, and here
is the one page that would settle it" is the thing that stops that happening.
It also makes the registry satisfy the platform's own first rule: what is
unknown is written as unknown, with what would resolve it.

PORTED, AND WHAT CHANGED
------------------------
The record shape, the claim checks and their messages are ported from the glycan
platform. Four changes:

- the SugarBase dataset entry is gone, with the glycan structure database it
  described. The DATASETS mechanism stays and is seeded with the two seed
  transcriptions, so the dataset-backed branch of the gate keeps being exercised
  rather than becoming unreachable code that nothing notices is dead.
- `url` is new. The project's own standard for a licence claim is the verbatim
  text from the publisher's page WITH THE URL AND THE DATE READ, and the ported
  record had nowhere to put the URL. Three of the sources below have no DOI at
  all, so a URL is the only key they have.
- `what_to_check` is new, and is required for an unverified entry.
- `evidence` and `notes` are coerced to tuples rather than merely checked, so a
  list cannot leave a frozen record holding a mutable field.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
from types import MappingProxyType
from typing import Mapping

from .reuse import (  # not from licensing: the gate imports this module
    ReuseStatus,
    can_train_commercial,
    can_use,
)


@dataclass(frozen=True)
class SourceLicence:
    """One source's licence, and the report of who read it where.

    Keyed by DOI where the source has one, and by URL where it does not. At
    least one of the two is required: a source nobody can navigate to cannot be
    re-checked, and an entry that cannot be re-checked is not evidence.
    """

    citation: str
    title: str
    licence: str
    reuse_status: ReuseStatus
    evidence: tuple[str, ...]  # what the reporter said they read, and where, in their words
    reported_by: str  # who reported the reading; not a claim that this code verified anything
    reported_on: date  # when the report was received
    attribution: str  # what the reporter said the licence asks for in return
    doi: str | None = None
    url: str | None = None
    # Required when the status is UNVERIFIED: the one thing that would settle it.
    # An unverified entry without this is just a shrug with a citation attached.
    what_to_check: str = ""
    # Required when this platform can use the source ONLY because of what this
    # platform is - academic_only and plain non_commercial terms. It records WHO
    # decided the platform is that, and when.
    #
    # Enforced rather than encouraged, because these are the entries whose
    # usability rests on a decision rather than on a licence, and an entry that
    # did not say so would read exactly like one whose terms permit anybody. If
    # the platform is ever commercialised these are the rows to revisit, and this
    # field is how somebody finds them.
    context_basis: str = ""
    # What the source holds, so the registry can be read as a data inventory
    # without opening five papers. Free text, and never treated as a count.
    content: str = ""
    notes: tuple[str, ...] = ()  # discrepancies met and how they were settled; what is still awaited

    def __post_init__(self) -> None:
        # A bare string is refused rather than coerced. tuple("RSC article page")
        # is seventeen one-character pieces of evidence, every one of which passes
        # the non-blank check below, and the entry then reads as thoroughly
        # evidenced while holding nothing. A string with spaces fails instead, for
        # the wrong reason ("records no evidence"). Both are silent corruption of
        # the single field this registry exists to hold, so the type is checked
        # before anything is coerced.
        for name in ("evidence", "notes"):
            if isinstance(getattr(self, name), str):
                raise ValueError(
                    f"{name} is a sequence of separate statements, not one string; wrap it as a tuple,"
                    f' as in {name}=("...",). A bare string would be split into one entry per character'
                )
        # Coerce before checking, so a list cannot leave a frozen record holding
        # a mutable, unhashable field that every later `in` test then trips on.
        object.__setattr__(self, "evidence", tuple(self.evidence))
        object.__setattr__(self, "notes", tuple(self.notes))
        if not (self.doi or self.url):
            raise ValueError(f"the licence for {self.citation!r} has neither a DOI nor a URL, so nobody can re-check it")
        if not self.licence.strip():
            raise ValueError(f"the licence for {self.key} names no licence")
        if not self.evidence or not all(item.strip() for item in self.evidence):
            raise ValueError(f"the licence for {self.key} records no evidence, so it cannot be relied on")
        if not self.reported_by.strip():
            raise ValueError(f"the licence for {self.key} names nobody who reported it")
        if not isinstance(self.reported_on, date):
            raise ValueError(f"the licence for {self.key} has no date for its report")
        if not isinstance(self.reuse_status, ReuseStatus):
            raise ValueError(f"the licence for {self.key} has a reuse status that is not a ReuseStatus")
        if self.reuse_status is ReuseStatus.UNVERIFIED and not self.what_to_check.strip():
            raise ValueError(
                f"the entry for {self.key} is unverified and does not say what would settle it. An unverified"
                " source recorded without that is indistinguishable from one nobody thought of, and the next"
                " person spends the same day on the same page"
            )
        # An entry the platform may use only because of what the platform is has
        # to say who decided that. See context_basis.
        if (
            not can_train_commercial(self.reuse_status)
            and can_use(self.reuse_status)
            and not self.context_basis.strip()
        ):
            raise ValueError(
                f"the entry for {self.key} has status '{self.reuse_status.value}', which this platform may"
                " use only because it is an academic and research platform rather than a commercial one."
                " Record context_basis: who decided that, and when. Without it this row is"
                " indistinguishable from one whose terms permit anybody"
            )
        # The one thing a permissive licence actually asks for in return. An
        # entry that clears training and names no attribution is a licence
        # obligation nobody can discharge.
        if can_use(self.reuse_status) and not self.attribution.strip():
            raise ValueError(
                f"the licence for {self.key} permits training but names no attribution; record what the"
                " licence asks for in return"
            )
        if self.doi:
            object.__setattr__(self, "doi", self.doi.strip().casefold())
        if self.url:
            object.__setattr__(self, "url", self.url.strip())

    @property
    def key(self) -> str:
        """How this entry is found: its DOI, or its URL where it has no DOI."""
        return self.doi or self.url or ""

    @property
    def permits_training(self) -> bool:
        return can_use(self.reuse_status)


@dataclass(frozen=True)
class DatasetLicence:
    """A shipped or seeded dataset's licence, keyed by the provenance its loader writes."""

    provenance: str  # exact; an analyte whose source is this string came from this dataset
    licence: str
    reuse_status: ReuseStatus
    attribution: str
    derived_from_doi: str | None = None  # the publication the dataset was transcribed from, where there is one

    def __post_init__(self) -> None:
        if not self.provenance.strip() or not self.licence.strip():
            raise ValueError("a dataset licence needs a provenance string and a licence")
        if not isinstance(self.reuse_status, ReuseStatus):
            raise ValueError(f"the dataset {self.provenance!r} has a reuse status that is not a ReuseStatus")


# =========================================================================================
# USABLE NOW
# =========================================================================================

STRUWE_2016 = SourceLicence(
    doi="10.1039/c6cc06247d",
    url="https://doi.org/10.1039/c6cc06247d",
    # Authors and volume/page as the owner gave them; the journal abbreviation, the
    # title and the page range are as the RSC listing shows them.
    citation="Struwe, Baldauf, Hofmann, Rudd & Pagel, Chem. Commun. 2016, 52, 12353-12356",
    title="Ion mobility separation of deprotonated oligosaccharide isomers - evidence for gas-phase charge migration",
    licence="Creative Commons Attribution 3.0 Unported (CC BY 3.0)",
    reuse_status=ReuseStatus.OPEN_ATTRIBUTION,
    content="28 TWIMS CCS values, 4 milk oligosaccharides, 4 adducts. Transcribed to CSV.",
    # The owner's report, in the owner's words. Nothing here was read by this code.
    evidence=(
        "RSC article page: states the article is licensed under a Creative Commons Attribution 3.0 Unported"
        " Licence, with material reusable without further permission given correct acknowledgement",
        "The issue listing: says Open Access, Creative Commons BY",
        "PDF header: marked Open Access Article",
    ),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 15),
    attribution="Correct acknowledgement. The owner reports that commercial use and derivatives are both permitted.",
    notes=(
        "Discrepancy met and settled: Europe PMC's record for this DOI was marked isOpenAccess N,"
        " 'subscription required'. Asked to confirm, the owner read the publisher's own page and reported"
        " the CC BY statement there; the publisher's statement governs and the Europe PMC flag is taken as stale.",
        "No CCS values from this paper were reconstructed from memory at any point. The values were awaited"
        " from a person with access to the article and its supplementary information.",
        "The values: the transcription in data/seed was delivered by the owner on 2026-09-16, reported as"
        " obtained from the paper's ESI, 28 rows. Nobody has independently checked it against the ESI within"
        " this repository; the loader's counts are of the transcription as delivered. Its ccs_2sd column is"
        " loaded as ccs_uncertainty with uncertainty_type two_sd, on the owner's instruction that the column"
        " reports two standard deviations.",
        "Four rows are held pending curation: LNH and LNnH share identical values for [M+H]+ (228.9) and"
        " [M+Cl]- (245.0), almost certainly one unresolved peak reported against both compounds. The hold is"
        " DERIVED by the loader's shared-peak check, not marked in the file.",
    ),
)

STRUWE_2015 = SourceLicence(
    doi="10.1039/c5an01092f",
    url="https://doi.org/10.1039/c5an01092f",
    citation="Struwe, Benesch, Harvey & Pagel, Analyst 2015, 140, 6799",
    title="(not recorded: the owner reported the DOI, authors, journal, volume and page, but not the title)",
    licence="Creative Commons Attribution 3.0 Unported (CC BY 3.0)",
    reuse_status=ReuseStatus.OPEN_ATTRIBUTION,
    content="89 TWIMS CCS values, 7 high-mannose N-glycans, 6 adducts, 4 sample origins. Transcribed to CSV.",
    evidence=(
        "PDF footer, on every page: states the article is licensed under a Creative Commons Attribution 3.0"
        " Unported Licence",
    ),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 17),
    attribution="Correct acknowledgement, as CC BY 3.0 asks.",
    notes=(
        "The values: 89 rows, delivered by the owner on 2026-09-17 and reported as transcribed from ESI"
        " Table S1. Nobody has checked the transcription against the ESI within this repository.",
        "The drift gas is UNSTATED, and deliberately so. The ESI says the travelling-wave values were"
        " calibrated against a dextran ladder of known DTCCS, but never says which gas those reference values"
        " were in. It is NOT inferred from the same group's 2016 paper. Resolving it means reading the"
        " calibration reference the owner names: Hofmann 2014, Anal. Chem. 86, 10789.",
        "Resolving the gas alone unblocks NOTHING. Every one of the 89 rows also carries uncertainty_type"
        " 'unknown' with no spread, which is an independent training blocker. Both would have to be"
        " resolved before any of these rows could train.",
        "And the second one probably cannot be resolved from this paper. The owner reports, on"
        " 2026-09-19, that the Analyst ESI carries no uncertainty column at all. If so, 'unknown' is not"
        " a transcription gap awaiting a more careful reading: it is a correct and complete record of"
        " what the source states, and these 89 rows are not pending curation. Unblocking them needs"
        " either a different source reporting a spread for the same ions, which would be different"
        " records with their own provenance, or an explicit decision that a value with no reported"
        " spread may be used for some purpose that does not need one. Never by defaulting the"
        " uncertainty type to something usable.",
        "No structure strings were transcribed. The table names compounds (Man5 and so on) and the structures"
        " sit in the paper's Figure 1, which nobody has read here; the owner did not fill them in from memory."
        " The composition column is therefore taken as transcribed rather than derived from a structure.",
        "Twenty-two of the values are conformer pairs: Man5, Man6, Man9 and Man9Glc each give two"
        " arrival-time peaks as [M-H]-, which MS/MS confirmed to be one structure folded two ways rather than"
        " two isomers. They are stored as conformers of one structure, never averaged.",
    ),
)

# =========================================================================================
# NAMED IN THE OBJECTIVE, LICENCE UNVERIFIED
#
# Every one of these is recorded so that nobody re-checks it from scratch, and
# every one names the single thing that would settle it. None may be ingested
# while it stands here: the gate is default-deny and unverified is not trainable.
# =========================================================================================

STEROID_INTERPLATFORM_2022 = SourceLicence(
    doi="10.1021/jasms.2c00196",
    url="https://doi.org/10.1021/jasms.2c00196",
    citation="Journal of the American Society for Mass Spectrometry, 2022",
    title="(not recorded: an interplatform steroid CCS comparison; the title has not been read here)",
    licence="ACS AuthorChoice open access; the specific open licence has not been read from the"
    " publisher's own page",
    reuse_status=ReuseStatus.ACADEMIC_ONLY,
    content="87 steroids, 142 CCS values, DTIMS + TWIMS + TIMS. The obvious benchmark set: the only named"
    " source that would give cross-platform matched ions directly. Open access with supporting information"
    " free of charge. PMC9545150.",
    evidence=(
        "Reported by Shawon Chakrabarty Kakon on 19 September 2026 as ACS AuthorChoice open access, with"
        " supporting information free of charge.",
        "The PMC record for PMC9545150 carries an ACS AuthorChoice banner, and the Europe PMC record for"
        " the DOI reports license 'cc by' with isOpenAccess Y. BOTH WERE READ BY THIS CODE ON 19 SEPTEMBER"
        " 2026 AND NEITHER IS ACCEPTED AS EVIDENCE OF THE LICENCE. CONTEXT.md's own standard lists a Europe"
        " PMC open-access flag among the things that do not count, having been observed stale for RSC"
        " titles, and a banner is not licence text.",
        "So the status here is the CONSERVATIVE reading, not the reported one. If the article really is"
        " CC BY then academic_only understates it and open_attribution would be correct; understating a"
        " licence costs nothing while this platform is academic, and overstating one cannot be undone.",
    ),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 19),
    attribution="Cite Feuerstein et al., J. Am. Soc. Mass Spectrom. 2022, as the paper asks. If the licence"
    " is confirmed as CC BY, correct attribution is all it requires.",
    context_basis="The CEO answered on 19 September 2026 that this platform is academic and research use, not commercial. That answer is the question CONTEXT.md records as the one deciding the data plan, and it is recorded in reuse.PLATFORM_USE_CONTEXT. The source's own terms are unchanged and are still recorded as what they are, so the gate refuses this row again on its own if the platform is ever commercialised.",
    what_to_check="ONE thing now. LICENCE: read the badge on the ACS article page itself and record the"
    " verbatim text with the URL and the date. That would settle whether this is open_attribution rather"
    " than academic_only, and it is the one reading that would let the data be used if the platform were"
    " ever commercialised. (The trapped-ion calibrant question that stood here is RESOLVED - see notes.)",
    notes=(
        "THE TRAPPED-ION CALIBRANT IS RESOLVED, 19 September 2026, and the history is worth keeping because"
        " the first reading was right and was still not the end of it. The supporting information gives only"
        " the mass calibration for TIM-MS ('10 mM sodium formate and a 7th order high-performance"
        " calibration') and then says 'In addition to the external calibration, each sample was"
        " automatically post-run calibrated by injecting a 1:1 mixture of BOTH calibrants' - naming one"
        " calibrant and presupposing a second. On the SI alone the honest record was UNSTATED_CALIBRANT and"
        " 142 values were held. THE ARTICLE STATES IT: 'Prior to analysis, the instrument was mass"
        " calibrated with sodium formate clusters (10 mM in 50:50 2-propanol/water) and TIM CCS N2 was"
        " calibrated using ions from Agilent ESI-L Tune Mix via a linear function.' Read from the Europe PMC"
        " open full text for PMC9545150, kept with its sha256 under data/raw/steroid_jasms2022/. All 142"
        " trapped-ion values now train and this corpus has three technologies. THE GENERAL LESSON: the"
        " supporting information is NOT the source, it is one document of it, and it omitted a method the"
        " article states plainly. Check the article before recording any calibrant as unstated.",
        "THE TRAVELLING-WAVE VALUES IN THIS FILE ARE NOT THIS STUDY'S OWN MEASUREMENTS, which matters for"
        " what a correction fitted on them means. The article: 'TWIM-MS data sets were reported in two of our"
        " previous publications and publicly available data was used for all comparisons.' The column"
        " ingested here is the INTERLABORATORY library - an average over four Waters instruments across"
        " several laboratories - from Hernandez-Mesa et al., Anal. Chem. 2020, 92, 5013-5022,"
        " doi:10.1021/acs.analchem.9b05247. The single-laboratory library, which is NOT ingested, is"
        " doi:10.1021/acs.analchem.7b05117. The drift-tube and trapped-ion values ARE the authors' own, new"
        " for this paper. So a DTIMS-versus-TIMS correction is within one laboratory and a TWIMS-versus-"
        " anything correction is one laboratory against an interlaboratory aggregate - which is still not an"
        " interlaboratory reproducibility figure, because there is an aggregate on one side and a single"
        " laboratory on the other, but it is not within-laboratory either.",
        "NO LICENCE LAUNDERING IS IMPLIED BY THAT. The values were obtained from this study's own supporting"
        " information, which is the artefact whose terms are recorded here, and Hernandez-Mesa, Dervilly and"
        " Le Bizec are authors of both papers: this is the same group republishing its own earlier data. The"
        " origin is recorded because provenance should name it, not because the licence chain needs it.",
        "THE PAPER ALSO SETTLES WHICH PLATFORMS SHARE A CALIBRANT, which is its own central finding:"
        " 'DT CCS N2 and TIM CCS N2 are routinely calibrated with the same commercially available compound"
        " mixture (i.e., reference ions and reference values) established by Stow et al., while TW CCS N2"
        " systems were calibrated using a different commercial calibrant mix.' So the drift-tube"
        " single-field and trapped-ion records carry the SAME calibrant string here, on the paper's"
        " authority rather than because two spellings looked alike.",
        "THIS IS THE ONE THAT MATTERS MOST. It is the only unverified source that would supply matched ions"
        " across three platforms, which is the premise of the whole pipeline. One page read decides whether"
        " the benchmark set exists.",
        "Its published figures are the validation targets for M3: 95 per cent of ions within 1 per cent bias"
        " for TIMS and 2 per cent for TWIMS relative to DTIMS, under 1.5 per cent of ions showing biases up to"
        " 7 per cent, and correlations of 0.9949 TWIMS-DTIMS, 0.9953 TIMS-DTIMS, 0.9989 TWIMS-TIMS. Those"
        " numbers are reported in the reference document; they have not been read from the paper here.",
        "The data itself: supporting file js2c00196_si_003.xlsx, sheet 'S2_Interplatform CCS Database',"
        " 142 data rows over five platform columns - travelling wave single- and cross-laboratory, trapped"
        " ion, single-field drift tube and stepped-field drift tube. Downloaded 19 September 2026 from"
        " Europe PMC's supplementaryFiles REST endpoint for PMC9545150, which is a documented open API."
        " The publisher's own copy answers 403 to an automated request and PMC serves a proof-of-work"
        " challenge for binary downloads; neither was worked around.",
        "It carries an m/z column, which makes it the first data in this repository from which a mass-based"
        " residual could be computed. See LIMITATIONS on why mass is still out of scope.",
    ),
)

METLIN_CCS = SourceLicence(
    url="https://metlin.scripps.edu/",
    citation="METLIN-CCS, Nature Methods, 2023",
    title="(not recorded)",
    licence="(UNVERIFIED: no licence text has been read)",
    reuse_status=ReuseStatus.UNVERIFIED,
    content="185,589 CCS values, 27,633 standards, TIMS on timsTOF Pro. The largest non-DTIMS resource.",
    evidence=(
        'The paper says the resource is "freely available" at METLIN, XCMS Online and PanoramaWeb. That is a'
        " statement about ACCESS, not a licence, and it was read as a description rather than as terms.",
    ),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 18),
    attribution="(unknown until the licence is read)",
    what_to_check="Read the terms of use on each of the three hosts named: METLIN, XCMS Online and"
    " PanoramaWeb. They may differ, and the one the data is actually taken from is the one that governs."
    " A download button with no terms beside it is not a licence.",
    notes=(
        "Single-platform (TIMS). Large, but on its own it yields no cross-platform pairs; its value is in"
        " pairing against a DTIMS source, which makes it useful only if a DTIMS source is also clear.",
    ),
)

BUSH_LAB_CCS = SourceLicence(
    url="https://biophysicalms.org/ccsdatabase",
    citation="Bush Lab CCS Database, biophysicalms.org",
    title="(not recorded)",
    licence="No terms of use posted; the page asks only that users cite the appropriate publications",
    reuse_status=ReuseStatus.ACADEMIC_ONLY,
    content="DOWNLOADED 19 September 2026 and counted from the file rather than reported: 213 KB, eight"
    " sheets, roughly 6,500 rows. Native-Like Protein Cations 1,000; Native-Like Protein Cations and"
    " complexes 989; Denatured Protein Cations 1,000; Polyalanine Cations 33; Anionic Homopolymers 1,000;"
    " Other Peptides 1,000; Small Molecular Ions 24; MicroSource Collection 1,441. Retrieved from the"
    " documented published-to-web xlsx URL, which needs no workaround. Most ions are reported to be measured"
    " in BOTH helium and nitrogen by similar DTIMS methods, which has not yet been checked against the file.",
    evidence=(
        'The page asks only that users "cite the appropriate publication(s)". Read by Shawon Chakrabarty'
        " Kakon; reported again 19 September 2026.",
        "NO TERMS OF USE ARE POSTED. That is the honest state of this entry and it is the weakest of the"
        " three reopened on 19 September 2026: the other two state their terms, and this one states a"
        " citation request and nothing else.",
    ),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 19),
    attribution="Cite the appropriate publications, as the page asks. That is the only thing the source"
    " asks for in return, and it is the only thing it says.",
    context_basis="The CEO answered on 19 September 2026 that this platform is academic and research use, not commercial. That answer is the question CONTEXT.md records as the one deciding the data plan, and it is recorded in reuse.PLATFORM_USE_CONTEXT. The source's own terms are unchanged and are still recorded as what they are, so the gate refuses this row again on its own if the platform is ever commercialised.",
    notes=(
        "The single best fit for the biopharmaceutical layer, which the objective document names as the"
        " differentiator: native AND denatured protein ions, which is exactly the folding-state distinction"
        " the matched-ion key was built to carry.",
        "Also the only named source measuring the same ions in both helium and nitrogen, which would let the"
        " gas-reference question be tested rather than assumed.",
        "WHAT THE FILE SHOWS BEFORE ANY CONVERSION, from opening the sheets on 19 September 2026. (1) THE"
        " GAS IS STATED, in the column headers - 'Omega(He)' and 'Omega(N2)' - so unlike CCSbase these"
        " records are matchable as soon as they are converted. Most ions carry both, which becomes two"
        " records on two keys rather than one. (2) THE UNITS DIFFER BETWEEN SHEETS BY A FACTOR OF 100: the"
        " protein and peptide sheets report nm^2, the polyalanine, homopolymer, small-molecule and"
        " MicroSource sheets report A^2. Getting this wrong puts a protein cross section out by two orders"
        " of magnitude, and both spellings look like a unit. (3) THE PROTEIN SHEETS CARRY NO ADDUCT, only a"
        " charge: a row says z = 3 and never says what the three charges are. CONTEXT.md predicted exactly"
        " this for exactly this source and the [M+24?]24+ form was built for it in M0 - PREDICTION"
        " CONFIRMED - so roughly 3,000 protein rows arrive unmatchable. (4) The MicroSource Collection is"
        " the clean part: 1,441 drug-like molecules with adduct, charge, nitrogen CCS, a per-row standard"
        " deviation, a formula, a CAS number and a reference. The most immediately usable table found in"
        " any source so far. (5) 41 cells in one sheet are corrupt, holding date serials far outside any"
        " valid range - 21955915 in C99, 9677161 in C108, and 39 more - which openpyxl refuses to read."
        " They must be counted and reported, never allowed to arrive as nulls.",
        "A DISCREPANCY TO SETTLE BEFORE INGESTING. This entry records the values as DTIMS; the database page"
        " describes the ions as 'primarily from traveling-wave ion mobility spectrometry'. Those are"
        " different platforms and the difference is the entire subject of this repository. The per-row Ref"
        " column names the paper for each value, so it is answerable per row rather than per file - and it"
        " must be answered rather than assumed.",
        "NOT INGESTED. The file is in hand and the adapter is not written; it is third in the stated order,"
        " after the steroid study and CCSbase.",
        "THE ONE RESERVATION ON THIS ENTRY, recorded because it is the only one of the three that rests on"
        " inference rather than on stated terms. A source posting NO terms defaults to all rights reserved,"
        " which is how AllCCS2 was treated and why AllCCS2 is excluded. This entry reads the citation"
        " request as an implicit grant for scholarly use, which is defensible for a public database"
        " published to a research community and is still a reading rather than a quotation. If the data"
        " matters enough to build on, the safe move is a written grant from the laboratory, and that is"
        " cheap to ask for.",
    ),
)

BAYESIAN_HARMONIZATION_2026 = SourceLicence(
    doi="10.1021/acs.analchem.5c06667",
    url="https://doi.org/10.1021/acs.analchem.5c06667",
    citation="Analytical Chemistry, 2026",
    title="(not recorded: a Bayesian CCS harmonization study)",
    licence="(UNVERIFIED: no licence text has been read)",
    reuse_status=ReuseStatus.UNVERIFIED,
    content="840 measurements, 347 compounds, three platforms. A ready-made matched-ion set if its data is"
    " released and licensed.",
    evidence=("Neither the licence nor the data availability statement has been read.",),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 18),
    attribution="(unknown until the licence is read)",
    what_to_check="Read both the licence badge and the data availability statement on the ACS article page."
    " The licence governing the ARTICLE may not govern the supporting data, and it is the data that is wanted.",
    notes=(
        "This is also the paper the objective document cites when it says generic harmonization is no longer"
        " novel. It is simultaneously the best available benchmark and the reason the platform's novelty has"
        " to come from the biopharmaceutical and cyclic-IMS layers instead.",
    ),
)

# =========================================================================================
# BLOCKED, AND RECORDED SO THAT NOBODY RE-CHECKS THEM
# =========================================================================================

CCSBASE = SourceLicence(
    url="https://ccsbase.net/",
    citation="CCSbase, ccsbase.net",
    title="(not recorded)",
    licence="Academic, non-commercial use only",
    reuse_status=ReuseStatus.ACADEMIC_ONLY,
    content="25,020 records, READ IN FULL on 19 September 2026 and counted from the response rather than"
    " reported: DT 10,967, TW 8,529, TIMS 5,524, across 23 named calibration methods and 36 primary"
    " references. Twelve columns: a stable per-record id, name, adduct, m/z, CCS, charge, SMILES, compound"
    " class, the reference, the platform and the calibration method. THE DRIFT GAS IS NOT ONE OF THEM, for"
    " any record - see notes, because that one absence decides whether any of this can pair.",
    evidence=(
        "The About page states that use must be for academic non-commercial purposes, and directs commercial"
        " users to Libin Xu and UW CoMotion.",
        "READ AGAIN ON 19 SEPTEMBER 2026 at https://ccsbase.net/about, and it says more than was recorded"
        " the first time. Verbatim: use 'must be for academic, non-commercial purposes only'; 'Any derivative"
        " works (e.g. softwares, websites) must reproduce the above copyright notice'; and published work"
        " must cite Ross, Cho & Xu 2020, Anal. Chem., doi:10.1021/acs.analchem.9b05772. The"
        " notice-reproduction clause is an OBLIGATION THIS PLATFORM WOULD TAKE ON by ingesting, and it is"
        " stronger than the plain citation requirement recorded before.",
    ),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 11),
    attribution="TWO obligations, both from the About page read 19 September 2026. (1) Cite Ross, Cho & Xu"
    " 2020, Anal. Chem., doi:10.1021/acs.analchem.9b05772. (2) Reproduce the site's copyright notice in any"
    " derivative work, which includes software and websites - so a deployed platform carrying these values"
    " has to carry the notice, not merely a citation. Each record also names one of 36 primary papers, and"
    " those are the origin of the values; citing the compilation does not cite them.",
    what_to_check="THE DRIFT GAS, per reference. It is not a column, so every record is unmatchable under"
    " the keying rule in LIMITATIONS 4.6 and CCSbase currently contributes ZERO matched ions. This is a"
    " bounded task rather than an open one: each record names one of 36 primary papers, those papers state"
    " the gas, and each paper resolved unlocks its own records. Do NOT default the database to nitrogen -"
    " one of its own method strings reads 'stepped-field DT-IMS, helium and nitrogen drift gas' for 120"
    " records, which is proof that the database is not uniformly nitrogen.",
    context_basis="The CEO answered on 19 September 2026 that this platform is academic and research use, not commercial. That answer is the question CONTEXT.md records as the one deciding the data plan, and it is recorded in reuse.PLATFORM_USE_CONTEXT. The source's own terms are unchanged and are still recorded as what they are, so the gate refuses this row again on its own if the platform is ever commercialised.",
    notes=(
        "REOPENED 19 September 2026. This was the source most affected by the open question put to the CEO,"
        " and the answer was academic and research use, so it is usable. Nothing about the source changed:"
        " its terms still restrict use to academic non-commercial purposes and the entry still says so.",
        "The terms are the clearest of the three reopened sources: the About page states the restriction"
        " explicitly and names where commercial users should go instead. This is what academic_only is for.",
        "Platform-aware, and that is why it matters here: it records which instrument a value came from,"
        " so it can supply the same compound on more than one platform. CONFIRMED against the data on"
        " 19 September 2026: the CCS Type column holds DT, TW or TIMS for every one of the 25,020 records,"
        " and a CCS method column names the calibrant for most of them.",
        "HOW IT WAS RETRIEVED, because it is not a downloadable file and CONTEXT.md implied it was. The"
        " site's /download button returns a 182-byte batch_query.csv, which is the TEMPLATE for the"
        " batch-query upload feature and not an export. The paper behind the site is not open access and has"
        " no supplementary data in Europe PMC. The table is paginated at ten rows over 2,502 pages. What"
        " worked was the site's own search form, which POSTs to /results: one broad query returned all"
        " 25,020 rows in a single response. One request rather than 2,502, the site used as intended rather"
        " than crawled, no robots.txt present, and the terms restrict the PURPOSE of use rather than the"
        " method. The response is kept with its sha256 under data/raw/ccsbase/.",
        "NOT INGESTED, and the reason is the drift gas rather than the licence. See what_to_check.",
        "PROVENANCE IS TWO-LAYERED and this entry covers only the outer layer. The values belong to 36"
        " primary papers; CCSbase is the compilation, and these are the compilation's terms. A converted"
        " record should therefore name CCSbase as its source - that is where the value was obtained and"
        " whose terms permit the use - and carry the primary paper's reference and DOI in source_locator."
        " NOTHING HERE ESTABLISHES that each of those 36 papers permits reuse, and this entry must not be"
        " read as claiming it does.",
    ),
)

SASTRE_TORANO_2025 = SourceLicence(
    doi="10.1038/s41467-025-67069-w",
    url="https://doi.org/10.1038/s41467-025-67069-w",
    citation="Sastre Torano et al., Nature Communications, 2025",
    title="De novo sequencing of glycans by ion mobility-mass spectrometry using a self-expanding database",
    licence="Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)",
    reuse_status=ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES,
    content="164 intact-ion CCS values and 394 SNFG structure images.",
    evidence=(
        "Rights and permissions section of the article page, read 15 September 2026, quoted verbatim by the"
        ' owner: "You do not have permission under this licence to share adapted material derived from this'
        ' article or parts of it."',
    ),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 15),
    attribution="None available to this platform: the licence permits neither commercial use nor derivatives.",
    notes=(
        "Recorded so that nobody rechecks it. Both clauses block this platform independently; the"
        " no-derivatives clause is the stricter, which is why this uses non_commercial_no_derivatives rather"
        " than non_commercial: a plain non-commercial record is inference-only, and this one has no permitted"
        " use here at all.",
    ),
)

MANABE_2022 = SourceLicence(
    url="https://pubs.acs.org/journal/jamsef",
    citation="Manabe et al., Journal of the American Society for Mass Spectrometry, 2022",
    title="(not recorded)",
    licence="Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)",
    reuse_status=ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES,
    content="71 PA-labelled N-glycans.",
    evidence=("Reported as CC BY-NC-ND 4.0; both clauses block a company independently.",),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 18),
    attribution="None available to this platform: the licence permits neither commercial use nor derivatives.",
    notes=(
        "No DOI was recorded for this source in the reference document, so it is keyed by the journal URL"
        " rather than by a DOI. That is a weaker key than the others and should be replaced with the DOI when"
        " somebody has it in front of them. The licence status itself is not in doubt.",
        "Glycan-specific, and so of limited use to this platform in any case; recorded to close it out.",
    ),
)

ALLCCS2 = SourceLicence(
    url="https://allccs.zhulab.cn/",
    citation="AllCCS2, allccs.zhulab.cn",
    title="(not recorded)",
    licence="No licence posted; a copyright notice only, which defaults to all rights reserved",
    reuse_status=ReuseStatus.EXCLUDED,
    content="Built for small molecules.",
    evidence=("The site carries a copyright notice and no licence or terms of use.",),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 18),
    attribution="None: no licence is offered.",
    notes=("No terms means all rights reserved. Excluded rather than unverified: the absence was checked.",),
)

GLYCOMOB = SourceLicence(
    url="http://www.glycomob.org/",
    citation="GlycoMob, glycomob.org",
    title="(not recorded)",
    licence="No terms stated",
    reuse_status=ReuseStatus.EXCLUDED,
    content="Over 900 glycan CCS values.",
    evidence=(
        "No terms are stated on the site. The site blocks automated access, and a recent review reports that"
        " it became inaccessible.",
    ),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 18),
    attribution="None: no licence is offered.",
    notes=("Excluded on terms, and very likely moot on availability.",),
)


# --- the seeded transcriptions, as datasets -----------------------------------------------
#
# These keep the dataset-backed branch of the licence gate alive. The glycan
# platform's only dataset was the structure database that does not come across;
# with DATASETS empty, dataset_for would always return None, the dataset branch
# of the claim check would become unreachable, and nothing would fail - it would
# simply stop being exercised, and the first real dataset in a later milestone
# would land on an untested path. These two entries are real: an analyte identity
# in either seed file was assigned by the transcription of that paper's ESI, and
# the paper's own CC BY 3.0 licence is what governs it.

SEED_STRUWE_2016 = DatasetLicence(
    provenance="wmxccs seed transcription: Struwe 2016 Chem Commun ESI Table S1",
    licence="Creative Commons Attribution 3.0 Unported (CC BY 3.0)",
    reuse_status=ReuseStatus.OPEN_ATTRIBUTION,
    attribution="Correct acknowledgement of Struwe et al., Chem. Commun. 2016, 52, 12353-12356.",
    derived_from_doi="10.1039/c6cc06247d",
)

SEED_STRUWE_2015 = DatasetLicence(
    provenance="wmxccs seed transcription: Struwe 2015 Analyst ESI Table S1",
    licence="Creative Commons Attribution 3.0 Unported (CC BY 3.0)",
    reuse_status=ReuseStatus.OPEN_ATTRIBUTION,
    attribution="Correct acknowledgement of Struwe et al., Analyst 2015, 140, 6799.",
    derived_from_doi="10.1039/c5an01092f",
)


_ENTRIES = (
    STRUWE_2016,
    STRUWE_2015,
    STEROID_INTERPLATFORM_2022,
    METLIN_CCS,
    BUSH_LAB_CCS,
    BAYESIAN_HARMONIZATION_2026,
    CCSBASE,
    SASTRE_TORANO_2025,
    MANABE_2022,
    ALLCCS2,
    GLYCOMOB,
)

# Keyed by DOI where there is one and by URL where there is not. A MappingProxyType
# and not a dict: the registry is evidence, and evidence that a caller can append
# to at runtime is not evidence.
REGISTRY: Mapping[str, SourceLicence] = MappingProxyType({entry.key: entry for entry in _ENTRIES})
DATASETS: Mapping[str, DatasetLicence] = MappingProxyType(
    {entry.provenance: entry for entry in (SEED_STRUWE_2016, SEED_STRUWE_2015)}
)


def licence_for(doi: str | None) -> SourceLicence | None:
    """The registry entry for a DOI, spelled any way, or None if there is none."""
    if not doi:
        return None
    return REGISTRY.get(doi.strip().casefold())


def source_for(key: str | None) -> SourceLicence | None:
    """The registry entry for a DOI or a URL, or None if there is none."""
    if not key:
        return None
    return REGISTRY.get(key.strip().casefold()) or REGISTRY.get(key.strip())


def dataset_for(provenance: str | None) -> DatasetLicence | None:
    """The dataset an analyte came from, by the exact provenance its loader wrote. No fuzzy matching."""
    if not provenance:
        return None
    return DATASETS.get(provenance.strip())


def unverified_sources() -> tuple[SourceLicence, ...]:
    """Every entry nobody has read the terms for, with what would settle each.

    The work list. Reported rather than buried, because an unverified source is
    a decision waiting on one page being read, and the platform's data plan is
    blocked on exactly these.
    """
    return tuple(entry for entry in _ENTRIES if entry.reuse_status is ReuseStatus.UNVERIFIED)


NO_DOI_FOR_PUBLISHED_CLAIM = (
    "the row claims {claimed} with no DOI, so the claim cannot be checked; a published source needs its DOI"
    " and a licence record, and only an internal measurement may be trainable without one"
)
NO_RECORD_FOR_DOI = (
    "the row claims {claimed} for DOI {doi}, which has no licence record; a trainable claim about a"
    " publication needs the licence read on the publisher's page and recorded first"
)
CLAIM_EXCEEDS_RECORD = (
    "the row claims {claimed} for DOI {doi}, but the licence record for it says {recorded}"
    " ({licence}); the record governs"
)
COMPONENT_CLAIM_UNBACKED = (
    "the analyte claims {claimed} but neither the row's DOI nor a registered dataset backs it"
    " (analyte source: {source!r}); an analyte's licence needs the paper it was read from, named as the"
    " measurement's own source so the row's DOI applies to it, or the dataset it was drawn from"
)
COMPONENT_CLAIM_EXCEEDS_DATASET = (
    "the analyte claims {claimed}, but it came from {source!r}, whose licence record says {recorded}"
    " ({licence}); the record governs"
)


def claim_problem(doi: str | None, claimed: ReuseStatus) -> str | None:
    """Why a row's reuse-status claim is not backed by the registry, or None if it is.

    Only a claim that would admit the row to training is checked: a row that
    says it may not train is telling the truth in the safe direction, and the
    gate refuses it regardless. An internal measurement has no DOI and needs no
    record; a claim of open licensing does, because that claim is about a
    publication and a publication has a page somebody can read.
    """
    if not isinstance(claimed, ReuseStatus):
        return f"the row's reuse status {claimed!r} is not a ReuseStatus, so its claim cannot be checked"
    if not can_use(claimed):
        return None
    entry = licence_for(doi)
    if claimed is ReuseStatus.SYNTHETIC_FIXTURE and entry is None:
        # A test's declaration that the record is invented, DOI included. Only a
        # DOI that is on record contradicts it, and then the record governs.
        return None
    if claimed is ReuseStatus.INTERNAL_PROPRIETARY and not doi:
        return None
    if not doi:
        return NO_DOI_FOR_PUBLISHED_CLAIM.format(claimed=claimed)
    if entry is None:
        return NO_RECORD_FOR_DOI.format(claimed=claimed, doi=doi)
    if entry.reuse_status is not claimed:
        return CLAIM_EXCEEDS_RECORD.format(claimed=claimed, doi=doi, recorded=entry.reuse_status, licence=entry.licence)
    return None


def component_claim_problem(
    doi: str | None,
    component_source: str | None,
    claimed: ReuseStatus,
    *,
    measurement_source: str | None = None,
) -> str | None:
    """Why a component's reuse-status claim is not backed, or None if it is.

    A component is anything a measurement is built from that carries a licence of
    its own and no DOI field: in this platform, the analyte identity record. In
    the glycan platform this rule was named for the glycan structure, which was
    the only such component; the rule itself never depended on what the component
    was, and it does not now.

    What backs a component is decided by what its source field NAMES, never by
    the claim alone:

    - a registered dataset, matched on its exact provenance string: the dataset's
      record governs, and a claim it does not match is refused, an internal claim
      included;
    - the row's own paper, which the component names by giving the same source as
      the measurement: the row's DOI and its licence record govern, exactly as
      for the measurement's claim. A row's DOI does not back a component whose
      source names something else, or the DOI would launder any identity at all;
    - anything else: an internal claim is an in-house assignment and needs no
      record, a synthetic claim is a test's declaration and needs none, and any
      other trainable claim is unbacked. That last case is the one the check
      exists for: "open_attribution" typed against an identity of unknown origin.
    """
    if not isinstance(claimed, ReuseStatus):
        return f"the analyte's reuse status {claimed!r} is not a ReuseStatus, so its claim cannot be checked"
    if not can_use(claimed):
        return None
    dataset = dataset_for(component_source)
    if dataset is not None:
        if dataset.reuse_status is claimed:
            return None
        return COMPONENT_CLAIM_EXCEEDS_DATASET.format(
            claimed=claimed, source=component_source, recorded=dataset.reuse_status, licence=dataset.licence
        )
    names_the_paper = (
        bool(doi)
        and component_source is not None
        and measurement_source is not None
        and component_source.strip() == measurement_source.strip()
    )
    if names_the_paper:
        return claim_problem(doi, claimed)
    if claimed in (ReuseStatus.INTERNAL_PROPRIETARY, ReuseStatus.SYNTHETIC_FIXTURE):
        return None
    return COMPONENT_CLAIM_UNBACKED.format(claimed=claimed, source=component_source)
