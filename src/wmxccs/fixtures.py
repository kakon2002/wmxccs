"""A synthetic corpus that exercises matched-ion construction, and cannot be mistaken for data.

WHY THIS EXISTS. There is no real cross-platform CCS data in this repository, and
there may not be until one ACS licence page has been read. The pairing rules do
not depend on which dataset fills them, so they are built and exercised here
instead of waiting.

WHY IT IS SAFE. Every record below carries ReuseStatus.SYNTHETIC_FIXTURE, on the
measurement and on the analyte inside it. That status is the one trainable status
that is not a claim about any source: it needs no licence record, the loader
refuses it from any file so it can never arrive as data, and
`matching.assert_quotable` refuses to let any report covering one be presented as
a result. The corpus is reachable from the package, so it is easy to run - and
impossible to quote.

The DOIs below are invented, and are invented on purpose. A synthetic record
citing a DOI that IS in the licence registry would be refused as synthetic by the
registry check, which is the behaviour that keeps a fixture from borrowing a real
paper's licence. These use the 10.9999 registrant, which is not allocated.

WHAT IT COVERS. One case per rule in matching.py, and for every rule a case that
must NOT match as well as one that must. The pairs that must not match are the
more important half: a matcher that pairs everything passes every test built only
from things that should pair.
"""

from __future__ import annotations

from .identity import (
    ADCAnalyte,
    AntibodyIdentity,
    DriftGas,
    FoldingState,
    GlycanAnalyte,
    GlycopeptideAnalyte,
    PeptideAnalyte,
    Polarity,
    ProteinAnalyte,
    SmallMoleculeAnalyte,
)
from .models import CCSMeasurement, CyclicSettings, DTIMSMethod, IMSType, PassMode, UncertaintyType
from .reuse import ReuseStatus

FIXTURE = ReuseStatus.SYNTHETIC_FIXTURE

# Not allocated to any registrant. A fixture must not cite a real paper.
DOI_A = "10.9999/fixture.a"
DOI_B = "10.9999/fixture.b"
DOI_REVIEW = "10.9999/fixture.review"

CAFFEINE = "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
GLUCOSE = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"


def synthetic_inchikey(tag: str) -> str:
    """A well-formed InChIKey that is obviously not a real one, distinct per tag.

    EVERY CASE BELOW GETS ITS OWN MOLECULE, and that is not decoration. The cases
    live in one corpus, and two cases sharing a molecule share a matched-ion key:
    the two-way match and the three-way match would merge into one five-way set,
    the republication case would collapse into them, and every count in the
    report would describe something nobody wrote. The first version of this file
    did exactly that.

    The middle block is FIXTUREAAA so that a real InChIKey can never be confused
    with one of these, in a log line or in a bug report.
    """
    letters = "".join(character for character in tag.upper() if character.isalpha())
    body = (letters + "XXXXXXXXXXXXXX")[:14]
    return f"{body}-FIXTUREAAA-N"

_PLATFORMS = {
    "dtims": dict(ims_type=IMSType.DTIMS, dtims_method=DTIMSMethod.STEPPED_FIELD, calibrant=None),
    "dtims_single": dict(ims_type=IMSType.DTIMS, dtims_method=DTIMSMethod.SINGLE_FIELD, calibrant="polyalanine"),
    "twims": dict(ims_type=IMSType.TWIMS, calibrant="polyalanine"),
    "tims": dict(ims_type=IMSType.TIMS, calibrant="polyalanine"),
    "cyclic": dict(
        ims_type=IMSType.CYCLIC,
        calibrant="polyalanine",
        cyclic=CyclicSettings(passes=1, pass_mode=PassMode.SINGLE_PASS, wrap_around=False),
    ),
}


def small_molecule(inchikey: str = CAFFEINE) -> SmallMoleculeAnalyte:
    return SmallMoleculeAnalyte(inchikey=inchikey, source="a fixture", reuse_status=FIXTURE)


def measurement(
    *,
    analyte=None,
    platform: str = "twims",
    ccs: float = 180.0,
    adduct: str = "[M+H]+",
    charge: int = 1,
    polarity: Polarity = Polarity.POSITIVE,
    drift_gas: DriftGas = DriftGas.N2,
    source: str = "fixture laboratory A",
    doi: str | None = DOI_A,
    conformer: int | None = None,
    conformers_total: int | None = None,
    reuse_status: ReuseStatus = FIXTURE,
    **overrides,
) -> CCSMeasurement:
    """One synthetic measurement. Everything a test wants to vary is a keyword."""
    fields = dict(
        analyte=analyte if analyte is not None else small_molecule(),
        adduct=adduct,
        charge=charge,
        polarity=polarity,
        drift_gas=drift_gas,
        ccs=ccs,
        ccs_uncertainty=1.0,
        uncertainty_type=UncertaintyType.SD,
        source=source,
        doi=doi,
        conformer=conformer,
        conformers_total=conformers_total,
        reuse_status=reuse_status,
        **_PLATFORMS[platform],
    )
    fields.update(overrides)
    return CCSMeasurement(**fields)


# --- the cases ----------------------------------------------------------------------------
#
# Each builder returns records for ONE rule, and its name says what should happen.


def a_three_way_match() -> tuple[CCSMeasurement, ...]:
    """One ion on DTIMS, TWIMS and TIMS. ONE matched set of three, not three pairs."""
    return tuple(
        measurement(analyte=small_molecule(synthetic_inchikey("threeway")), platform=platform, ccs=ccs, source=f"fixture laboratory {letter}", doi=doi)
        for platform, ccs, letter, doi in (
            ("dtims", 180.0, "A", DOI_A),
            ("twims", 183.6, "B", DOI_B),
            ("tims", 181.2, "C", "10.9999/fixture.c"),
        )
    )


def a_two_way_match() -> tuple[CCSMeasurement, ...]:
    """One ion on two platforms. MATCHES."""
    return (
        measurement(analyte=small_molecule(synthetic_inchikey("twoway")), platform="dtims", ccs=180.0),
        measurement(analyte=small_molecule(synthetic_inchikey("twoway")), platform="twims", ccs=183.6, source="fixture laboratory B", doi=DOI_B),
    )


def replicates_on_one_platform() -> tuple[CCSMeasurement, ...]:
    """Two DTIMS values of one ion, differing slightly. DOES NOT MATCH: one platform."""
    return (
        measurement(analyte=small_molecule(synthetic_inchikey("replicates")), platform="dtims", ccs=180.0),
        measurement(analyte=small_molecule(synthetic_inchikey("replicates")), platform="dtims", ccs=180.4, source="fixture laboratory B", doi=DOI_B),
    )


def one_measurement_republished() -> tuple[CCSMeasurement, ...]:
    """The same DTIMS value in an original paper and in a review. ONE measurement, not a pair.

    Same ion, same platform, same value, different DOI. This is the case the
    deduplication rule exists for, and it must not become a matched pair.
    """
    return (
        measurement(analyte=small_molecule(synthetic_inchikey("republished")), platform="dtims", ccs=180.0, source="fixture laboratory A", doi=DOI_A),
        measurement(analyte=small_molecule(synthetic_inchikey("republished")), platform="dtims", ccs=180.0, source="a review reprinting it", doi=DOI_REVIEW),
    )


def two_laboratories_agreeing_exactly() -> tuple[CCSMeasurement, ...]:
    """Two PLATFORMS reporting the identical value. MATCHES, and is not a duplicate.

    The counterweight to the republication case. Collapsing on the value alone
    would delete exactly the evidence this platform exists to find.
    """
    return (
        measurement(analyte=small_molecule(synthetic_inchikey("agreeing")), platform="dtims", ccs=180.0, doi=DOI_A),
        measurement(analyte=small_molecule(synthetic_inchikey("agreeing")), platform="twims", ccs=180.0, source="fixture laboratory B", doi=DOI_B),
    )


def two_conformers_on_one_platform() -> tuple[CCSMeasurement, ...]:
    """One ion, two arrival-time peaks, one instrument. NOT duplicates, NOT a match."""
    return (
        measurement(analyte=small_molecule(synthetic_inchikey("conformers")), platform="dtims", ccs=180.0, conformer=1, conformers_total=2),
        measurement(analyte=small_molecule(synthetic_inchikey("conformers")), platform="dtims", ccs=192.0, conformer=2, conformers_total=2),
    )


def the_same_conformer_on_two_platforms() -> tuple[CCSMeasurement, ...]:
    """Conformer 1 on two platforms. Groups, and carries a blocker.

    Conformer numbering is source-local. That conformer 1 in one paper is
    conformer 1 in another has not been established by anything, so the set is
    built and refused rather than silently used.
    """
    return (
        measurement(analyte=small_molecule(synthetic_inchikey("sameconformer")), platform="dtims", ccs=180.0, conformer=1, conformers_total=2),
        measurement(analyte=small_molecule(synthetic_inchikey("sameconformer")), 
            platform="twims", ccs=183.0, conformer=1, conformers_total=2, source="fixture laboratory B", doi=DOI_B
        ),
    )


def two_unstated_charge_carriers() -> tuple[CCSMeasurement, ...]:
    """Two 24+ protein ions, neither naming its carrier. NEVER MATCH, not even each other."""
    native = ProteinAnalyte(
        accession="P01857", folding_state=FoldingState.NATIVE, source="a fixture", reuse_status=FIXTURE
    )
    return (
        measurement(
            analyte=native, platform="dtims", adduct="[M+24?]24+", charge=24, ccs=7000.0, doi=DOI_A
        ),
        measurement(
            analyte=native,
            platform="twims",
            adduct="[M+24?]24+",
            charge=24,
            ccs=7100.0,
            source="fixture laboratory B",
            doi=DOI_B,
        ),
    )


def a_native_and_a_denatured_antibody() -> tuple[CCSMeasurement, ...]:
    """Same antibody, different conformation and charge. DO NOT MATCH."""
    identity = AntibodyIdentity(inn="trastuzumab")
    from .identity import IntactAntibodyAnalyte

    native = IntactAntibodyAnalyte(
        antibody=identity, folding_state=FoldingState.NATIVE, source="a fixture", reuse_status=FIXTURE
    )
    denatured = IntactAntibodyAnalyte(
        antibody=identity, folding_state=FoldingState.DENATURED, source="a fixture", reuse_status=FIXTURE
    )
    return (
        measurement(analyte=native, platform="dtims", adduct="[M+24H]24+", charge=24, ccs=7000.0),
        measurement(
            analyte=denatured,
            platform="twims",
            adduct="[M+40H]40+",
            charge=40,
            ccs=11000.0,
            source="fixture laboratory B",
            doi=DOI_B,
        ),
    )


def two_adcs_differing_only_in_dar() -> tuple[CCSMeasurement, ...]:
    """DAR 2 and DAR 4 of one antibody, on two platforms. DO NOT MATCH: different molecules."""

    def adc(dar: int) -> ADCAnalyte:
        return ADCAnalyte(
            antibody=AntibodyIdentity(inn="trastuzumab"),
            linker_payload_class="vc-MMAE",
            dar=dar,
            folding_state=FoldingState.NATIVE,
            source="a fixture",
            reuse_status=FIXTURE,
        )

    return (
        measurement(analyte=adc(2), platform="dtims", adduct="[M+24H]24+", charge=24, ccs=7200.0),
        measurement(
            analyte=adc(4),
            platform="twims",
            adduct="[M+24H]24+",
            charge=24,
            ccs=7400.0,
            source="fixture laboratory B",
            doi=DOI_B,
        ),
    )


def a_glycopeptide_and_its_backbone() -> tuple[CCSMeasurement, ...]:
    """A glycopeptide and its bare peptide, on two platforms. DO NOT MATCH."""
    glycopeptide = GlycopeptideAnalyte(
        sequence="NLTK",
        glycan_composition="Hex5HexNAc2",
        attachment_site="N297",
        source="a fixture",
        reuse_status=FIXTURE,
    )
    backbone = PeptideAnalyte(sequence="NLTK", source="a fixture", reuse_status=FIXTURE)
    return (
        measurement(analyte=glycopeptide, platform="dtims", ccs=420.0),
        measurement(
            analyte=backbone, platform="twims", ccs=260.0, source="fixture laboratory B", doi=DOI_B
        ),
    )


def two_gases_for_one_ion() -> tuple[CCSMeasurement, ...]:
    """One ion, one platform pair, two reference gases. DO NOT MATCH."""
    return (
        measurement(analyte=small_molecule(synthetic_inchikey("twogases")), platform="dtims", ccs=180.0, drift_gas=DriftGas.N2),
        measurement(analyte=small_molecule(synthetic_inchikey("twogases")), 
            platform="twims", ccs=125.0, drift_gas=DriftGas.HE, source="fixture laboratory B", doi=DOI_B
        ),
    )


def a_cyclic_pair_at_different_pass_counts() -> tuple[CCSMeasurement, ...]:
    """One ion on DTIMS and on cyclic at six passes. MATCHES: passes do not change the ion."""
    six = CyclicSettings(passes=6, pass_mode=PassMode.MULTIPASS, wrap_around=False)
    return (
        measurement(analyte=small_molecule(synthetic_inchikey("cyclicpair")), platform="dtims", ccs=180.0),
        measurement(analyte=small_molecule(synthetic_inchikey("cyclicpair")), 
            platform="cyclic",
            ccs=181.5,
            cyclic=six,
            source="fixture laboratory B",
            doi=DOI_B,
        ),
    )


def a_primary_and_a_calibrated_dtims_pair() -> tuple[CCSMeasurement, ...]:
    """One ion, stepped-field and single-field DTIMS. MATCHES: they are two platforms.

    Stepped-field is primary and single-field is calibrated, so comparing them
    measures the thing this pipeline exists to measure. See matching._platform_of.
    """
    return (
        measurement(analyte=small_molecule(synthetic_inchikey("primarycal")), platform="dtims", ccs=180.0),
        measurement(
            analyte=small_molecule(synthetic_inchikey("primarycal")),
            platform="dtims_single",
            ccs=181.9,
            source="fixture laboratory B",
            doi=DOI_B,
        ),
    )


def two_conformers_sharing_one_value() -> tuple[CCSMeasurement, ...]:
    """Two conformers of one ion reported at the SAME value, on one platform.

    Contrived, and it has to exist. The ordinary conformer case gives its two
    peaks different values, so the VALUE alone keeps those records apart and the
    conformer index in the deduplication key is never actually load-bearing
    there. Here it is the only thing standing between two legitimate rows and
    being collapsed into one, which is what a conformer pair must never be.
    """
    return (
        measurement(
            analyte=small_molecule(synthetic_inchikey("sharedvalue")),
            platform="dtims",
            ccs=180.0,
            conformer=1,
            conformers_total=2,
        ),
        measurement(
            analyte=small_molecule(synthetic_inchikey("sharedvalue")),
            platform="dtims",
            ccs=180.0,
            conformer=2,
            conformers_total=2,
        ),
    )


def a_match_with_an_unusable_member() -> tuple[CCSMeasurement, ...]:
    """A genuine cross-platform pair, one member of which nobody may use.

    The set is built - it IS the same ion - and it is refused, with the refusal
    naming which member and why. Dropping it silently would report a smaller
    corpus with nothing to say a source is waiting on a licence check.
    """
    return (
        measurement(analyte=small_molecule(synthetic_inchikey("unusable")), platform="dtims", ccs=180.0),
        measurement(analyte=small_molecule(synthetic_inchikey("unusable")), 
            platform="twims",
            ccs=183.6,
            source="a source nobody has checked",
            doi=DOI_B,
            reuse_status=ReuseStatus.UNVERIFIED,
        ),
    )


def two_different_molecules() -> tuple[CCSMeasurement, ...]:
    """Two InChIKeys, two platforms. DO NOT MATCH."""
    return (
        measurement(analyte=small_molecule(CAFFEINE), platform="dtims", ccs=180.0),
        measurement(
            analyte=small_molecule(GLUCOSE),
            platform="twims",
            ccs=181.0,
            source="fixture laboratory B",
            doi=DOI_B,
        ),
    )


def two_glycans_of_one_composition() -> tuple[CCSMeasurement, ...]:
    """Two resolved isomers sharing a composition, on two platforms. DO NOT MATCH.

    The LNH/LNnH case, which is the reason glycan identity keys on the finest
    structure rather than on the composition.
    """
    lnh = GlycanAnalyte(
        iupac_condensed="Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)GlcNAc(b1-3)]Gal(b1-4)Glc",
        composition="Hex4HexNAc2",
        source="a fixture",
        reuse_status=FIXTURE,
    )
    lnnh = GlycanAnalyte(
        iupac_condensed="Gal(b1-4)GlcNAc(b1-3)[Gal(b1-4)GlcNAc(b1-6)]Gal(b1-4)Glc",
        composition="Hex4HexNAc2",
        source="a fixture",
        reuse_status=FIXTURE,
    )
    return (
        measurement(analyte=lnh, platform="dtims", ccs=228.9),
        measurement(
            analyte=lnnh, platform="twims", ccs=245.0, source="fixture laboratory B", doi=DOI_B
        ),
    )


# Every case, by what it demonstrates. The corpus is the union.
CASES = {
    "a three-way match across DTIMS, TWIMS and TIMS": a_three_way_match,
    "a two-way match": a_two_way_match,
    "replicates on one platform": replicates_on_one_platform,
    "one measurement republished in a review": one_measurement_republished,
    "two laboratories agreeing exactly": two_laboratories_agreeing_exactly,
    "two conformers on one platform": two_conformers_on_one_platform,
    "the same conformer index on two platforms": the_same_conformer_on_two_platforms,
    "two unstated charge carriers": two_unstated_charge_carriers,
    "a native and a denatured antibody": a_native_and_a_denatured_antibody,
    "two ADCs differing only in DAR": two_adcs_differing_only_in_dar,
    "a glycopeptide and its bare backbone": a_glycopeptide_and_its_backbone,
    "one ion referred to two gases": two_gases_for_one_ion,
    "a cyclic pair at different pass counts": a_cyclic_pair_at_different_pass_counts,
    "a primary and a calibrated DTIMS pair": a_primary_and_a_calibrated_dtims_pair,
    "two conformers sharing one value": two_conformers_sharing_one_value,
    "a match with an unusable member": a_match_with_an_unusable_member,
    "two different molecules": two_different_molecules,
    "two glycan isomers of one composition": two_glycans_of_one_composition,
}


def corpus() -> tuple[CCSMeasurement, ...]:
    """Every case at once. What matched-ion construction is exercised against."""
    return tuple(record for build in CASES.values() for record in build())
