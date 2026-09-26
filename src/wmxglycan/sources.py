"""Where each source came from, on what terms, and who reported checking.

A row in a measurements file may claim any reuse status it likes. The claim is
only as good as the evidence behind it, and this registry is that evidence. Two
kinds of source are recorded, because a value and the structure it was measured
on can come from different places:

- a PUBLICATION, keyed by DOI, for values transcribed from a paper;
- a DATASET, keyed by the exact provenance string its own loader writes into a
  structure's source field, for structures drawn from a shipped database.

The licence gate checks every trainable claim against both, so "open_attribution"
on a record traces back to something a person read rather than something a
person typed, whether the record came from a file or was built in code. The one
trainable status that is not a claim about a source is synthetic_fixture, a
test's declaration that the record is invented; it needs no record, and a
record citing a DOI that IS on record is that paper's record and may not call
itself synthetic.

Nothing here is inferred by this code. A publication's licence is recorded
after a person has read it on the publisher's page and reported what they read;
the record carries that report, attributed, and does not claim to quote the page
verbatim unless it says so. A third party's metadata flag is not evidence either
way. A dataset's licence is the one its own loader module records from the
package it ships in, and the dataset entry here points at that module rather
than restating a reading of its own.

The registries are read-only mappings. An entry is added by editing this
module, with its evidence, and by no other route.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
from types import MappingProxyType
from typing import Mapping

from .reuse import ReuseStatus, can_train_commercial  # not from licensing: the gate imports this module
from .sugarbase import SUGARBASE_ATTRIBUTION, SUGARBASE_LICENCE, dataset_version


@dataclass(frozen=True)
class SourceLicence:
    """One publication's licence, and the report of who read it where."""

    doi: str
    citation: str
    title: str
    licence: str
    reuse_status: ReuseStatus
    evidence: tuple[str, ...]  # what the reporter said they read, and where, in their words
    reported_by: str  # who reported the reading; not a claim that this code verified anything
    reported_on: date  # when the report was received
    attribution: str  # what the reporter said the licence asks for in return
    notes: tuple[str, ...] = ()  # discrepancies met and how they were settled; what is still awaited

    def __post_init__(self) -> None:
        if not self.doi.strip():
            raise ValueError("a source licence needs a DOI")
        if not self.licence.strip():
            raise ValueError(f"the licence for {self.doi} names no licence")
        if not self.evidence or not all(item.strip() for item in self.evidence):
            raise ValueError(f"the licence for {self.doi} records no evidence, so it cannot be relied on")
        if not self.reported_by.strip():
            raise ValueError(f"the licence for {self.doi} names nobody who reported it")
        if not isinstance(self.reported_on, date):
            raise ValueError(f"the licence for {self.doi} has no date for its report")
        if not isinstance(self.reuse_status, ReuseStatus):
            raise ValueError(f"the licence for {self.doi} has a reuse status that is not a ReuseStatus")
        object.__setattr__(self, "doi", self.doi.strip().casefold())

    @property
    def permits_training(self) -> bool:
        return can_train_commercial(self.reuse_status)


@dataclass(frozen=True)
class DatasetLicence:
    """A shipped dataset's licence, keyed by the provenance its loader writes."""

    provenance: str  # exact; a structure whose source is this string came from this dataset
    licence: str
    reuse_status: ReuseStatus
    attribution: str

    def __post_init__(self) -> None:
        if not self.provenance.strip() or not self.licence.strip():
            raise ValueError("a dataset licence needs a provenance string and a licence")
        if not isinstance(self.reuse_status, ReuseStatus):
            raise ValueError(f"the dataset {self.provenance!r} has a reuse status that is not a ReuseStatus")


STRUWE_2016 = SourceLicence(
    doi="10.1039/c6cc06247d",
    # Authors and volume/page as the owner gave them; the journal abbreviation, the
    # title and the page range are as the RSC listing shows them.
    citation="Struwe, Baldauf, Hofmann, Rudd & Pagel, Chem. Commun. 2016, 52, 12353-12356",
    title="Ion mobility separation of deprotonated oligosaccharide isomers - evidence for gas-phase charge migration",
    licence="Creative Commons Attribution 3.0 Unported (CC BY 3.0)",
    reuse_status=ReuseStatus.OPEN_ATTRIBUTION,
    # The owner's report, in the owner's words. Nothing here was read by this code.
    evidence=(
        "RSC article page: states the article is licensed under a Creative Commons Attribution 3.0 Unported"
        " Licence, with material reusable without further permission given correct acknowledgement",
        "The issue listing: says Open Access, Creative Commons BY",
        "PDF header: marked Open Access Article",
    ),
    # Stated directly by the owner, in the session of 2026-09-16 by the session clock:
    # "reported_by is Shawon Chakrabarty Kakon. I read the RSC page."
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 15),  # the session in which the licence was first reported, by the session clock
    attribution="Correct acknowledgement. The owner reports that commercial use and derivatives are both permitted.",
    notes=(
        "Discrepancy met and settled: Claude found Europe PMC's record for this DOI marked isOpenAccess N,"
        " 'subscription required'. Asked to confirm, the owner read the publisher's own page and reported"
        " the CC BY statement there; the publisher's statement governs and the Europe PMC flag is taken as stale.",
        "Process note by Claude: no CCS values from this paper were reconstructed from memory at any point."
        " The automated routes to the full text that Claude tried in the session were closed, and the values"
        " were awaited from a person with access to the article and its supplementary information.",
        "Scientific finding, as the owner summarised it pending the paper itself: LNH and LNnH show comparable"
        " CCS as [M+Na]+ and differ by over 7 per cent as [M-H]-. The owner corrected an earlier, wider"
        " summary as not from any source they had read, and said nothing else about that pair should be"
        " treated as established until the paper is in hand.",
        "The values: the transcription in data/raw was delivered by the owner on 2026-09-16 (session clock),"
        " reported as obtained from the paper's ESI, 28 rows. Nobody has independently checked it against the"
        " ESI within this repository; the loader's counts are of the transcription as delivered. Its ccs_2sd"
        " column is loaded as ccs_uncertainty with uncertainty_type two_sd, on the owner's instruction that"
        " the column reports two standard deviations.",
    ),
)

STRUWE_2015 = SourceLicence(
    doi="10.1039/c5an01092f",
    # Authors, journal, year, volume and page exactly as the owner gave them.
    citation="Struwe, Benesch, Harvey & Pagel, Analyst 2015, 140, 6799",
    title="(not recorded: the owner reported the DOI, authors, journal, volume and page, but not the title)",
    licence="Creative Commons Attribution 3.0 Unported (CC BY 3.0)",
    reuse_status=ReuseStatus.OPEN_ATTRIBUTION,
    evidence=(
        "PDF footer, on every page: states the article is licensed under a Creative Commons Attribution 3.0"
        " Unported Licence",
    ),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 17),  # the session in which the licence was reported, by the session clock
    attribution="Correct acknowledgement, as CC BY 3.0 asks.",
    notes=(
        "The values: data/raw/struwe2015_analyst_ccs.csv, 89 rows, delivered by the owner on 2026-09-17"
        " (session clock) and reported as transcribed from ESI Table S1. Nobody has checked the transcription"
        " against the ESI within this repository; the loader's counts are of the file as delivered.",
        "The drift gas is UNSTATED, and deliberately so. The owner reports that the ESI says the"
        " travelling-wave values were calibrated against a dextran ladder of known DTCCS, but never says which"
        " gas those reference values were in. It is NOT inferred from the same group's 2016 paper. Resolving"
        " it means reading the calibration reference the owner names: Hofmann 2014, Anal. Chem. 86, 10789.",
        "No structure strings were transcribed. The table names compounds (Man5 and so on) and the structures"
        " sit in the paper's Figure 1, which nobody has read here; the owner did not fill them in from memory."
        " The composition column is therefore taken as transcribed rather than derived from a structure, so a"
        " composition typo in this file cannot be caught the way the 2016 adapter catches one.",
        "Twenty-two of the values are conformer pairs: Man5, Man6, Man9 and Man9Glc each give two"
        " arrival-time peaks as [M-H]-, which the owner reports MS/MS confirmed to be one structure folded two"
        " ways rather than two isomers. They are stored as conformers of one structure, never averaged.",
    ),
)

SASTRE_TORANO_2025 = SourceLicence(
    doi="10.1038/s41467-025-67069-w",
    citation="Sastre Toraño et al., Nature Communications, 2025",
    title="De novo sequencing of glycans by ion mobility-mass spectrometry using a self-expanding database",
    licence="Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)",
    reuse_status=ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES,
    evidence=(
        "Rights and permissions section of the article page, read 15 September 2026, quoted verbatim by the"
        ' owner: "You do not have permission under this licence to share adapted material derived from this'
        ' article or parts of it."',
    ),
    reported_by="Shawon Chakrabarty Kakon",
    reported_on=date(2026, 9, 15),
    attribution="None available to this platform: the licence permits neither commercial use nor derivatives.",
    notes=(
        "Recorded so that nobody rechecks it. This is the dataset LIMITATIONS.md carried as the last pending"
        " CCS candidate: roughly 87 usable intact species, with structures embedded as 394 PNG images that"
        " would need manual annotation. The licence closes it as BLOCKED rather than pending, and the"
        " annotation cost is moot.",
        "Both clauses block this platform independently, as the owner put it. The no-derivatives clause is the"
        " stricter, and is why this entry uses non_commercial_no_derivatives rather than non_commercial: a"
        " plain non-commercial record is inference-only, and this one has no permitted use here at all.",
    ),
)

# The structure database this repository already loads. Its licence is the one
# sugarbase.py records from the glycowork package it ships in, not a reading of
# this module's own; this entry lets a measurement whose structure came from it
# be backed without guessing from how a source is spelt.
SUGARBASE = DatasetLicence(
    provenance=dataset_version().provenance,
    licence=SUGARBASE_LICENCE,
    reuse_status=ReuseStatus.OPEN_ATTRIBUTION,
    attribution=SUGARBASE_ATTRIBUTION,
)

REGISTRY: Mapping[str, SourceLicence] = MappingProxyType(
    {entry.doi: entry for entry in (STRUWE_2016, STRUWE_2015, SASTRE_TORANO_2025)}
)
DATASETS: Mapping[str, DatasetLicence] = MappingProxyType({entry.provenance: entry for entry in (SUGARBASE,)})


def licence_for(doi: str | None) -> SourceLicence | None:
    """The registry entry for a DOI, spelled any way, or None if there is none."""
    if not doi:
        return None
    return REGISTRY.get(doi.strip().casefold())


def dataset_for(provenance: str | None) -> DatasetLicence | None:
    """The dataset a structure came from, by the exact provenance its loader wrote. No fuzzy matching."""
    if not provenance:
        return None
    return DATASETS.get(provenance.strip())


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
STRUCTURE_CLAIM_UNBACKED = (
    "the structure claims {claimed} but neither the row's DOI nor a registered dataset backs it"
    " (structure source: {source!r}); a structure's licence needs the paper it was read from, named as the"
    " measurement's own source so the row's DOI applies to it, or the database it was drawn from"
)
STRUCTURE_CLAIM_EXCEEDS_DATASET = (
    "the structure claims {claimed}, but it came from {source!r}, whose licence record says {recorded}"
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
    if not can_train_commercial(claimed):
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


def structure_claim_problem(
    doi: str | None,
    structure_source: str | None,
    claimed: ReuseStatus,
    *,
    measurement_source: str | None = None,
) -> str | None:
    """Why a structure's reuse-status claim is not backed, or None if it is.

    What backs a structure is decided by what its source field NAMES, never by
    the claim alone:

    - a registered dataset, matched on its exact provenance string: the
      dataset's record governs, and a claim it does not match is refused, an
      internal claim included;
    - the row's own paper, which the structure names by giving the same source
      as the measurement: the row's DOI and its licence record govern, exactly
      as for the measurement's claim. A row's DOI does not back a structure
      whose source names something else, or the DOI would launder any
      structure at all;
    - anything else: an internal claim is an in-house assignment and needs no
      record, a synthetic claim is a test's declaration and needs none, and
      any other trainable claim is unbacked. That last case is the one the
      check exists for: "open_attribution" typed against a structure of
      unknown origin.
    """
    if not isinstance(claimed, ReuseStatus):
        return f"the structure's reuse status {claimed!r} is not a ReuseStatus, so its claim cannot be checked"
    if not can_train_commercial(claimed):
        return None
    dataset = dataset_for(structure_source)
    if dataset is not None:
        if dataset.reuse_status is claimed:
            return None
        return STRUCTURE_CLAIM_EXCEEDS_DATASET.format(
            claimed=claimed, source=structure_source, recorded=dataset.reuse_status, licence=dataset.licence
        )
    names_the_paper = (
        bool(doi)
        and structure_source is not None
        and measurement_source is not None
        and structure_source.strip() == measurement_source.strip()
    )
    if names_the_paper:
        return claim_problem(doi, claimed)
    if claimed in (ReuseStatus.INTERNAL_PROPRIETARY, ReuseStatus.SYNTHETIC_FIXTURE):
        return None
    return STRUCTURE_CLAIM_UNBACKED.format(claimed=claimed, source=structure_source)
