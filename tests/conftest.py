"""Shared fixtures and builders.

The glycan platform had no conftest at all: every test file redeclared its own
helpers, and the minimal valid record was written out by hand in twenty places.
That is why a change to one required field there touched a dozen files. One place
here, so that adding a required field breaks one builder rather than the suite.

EVERY BUILDER HERE PRODUCES A SYNTHETIC FIXTURE. The reuse status says so, in the
record and in the analyte inside it. That is the status's whole purpose: a record
built in code is distinguishable from a real one at a glance, it needs no licence
record, and the loader refuses it from any file, so it cannot leak into data.
"""

from __future__ import annotations

import pytest

from wmxccs.identity import (
    ADCAnalyte,
    AntibodyIdentity,
    DriftGas,
    FoldingState,
    GlycanAnalyte,
    IntactAntibodyAnalyte,
    PeptideAnalyte,
    Polarity,
    ProteinAnalyte,
    SmallMoleculeAnalyte,
)
from wmxccs.models import CCSMeasurement, CyclicSettings, IMSType, PassMode
from wmxccs.reuse import ReuseStatus

FIXTURE = ReuseStatus.SYNTHETIC_FIXTURE
SOURCE = "a test"

# A real InChIKey shape; the molecule it names is not asserted anywhere.
INCHIKEY = "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
OTHER_INCHIKEY = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

# The two milk-oligosaccharide structures that share a composition and do not
# share a cross section. Used wherever a test needs two analytes that a
# composition-keyed platform would wrongly call one.
LNH_IUPAC = "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)GlcNAc(b1-3)]Gal(b1-4)Glc"
LNNH_IUPAC = "Gal(b1-4)GlcNAc(b1-3)[Gal(b1-4)GlcNAc(b1-6)]Gal(b1-4)Glc"


def small_molecule(**overrides) -> SmallMoleculeAnalyte:
    fields = dict(inchikey=INCHIKEY, source=SOURCE, reuse_status=FIXTURE)
    fields.update(overrides)
    return SmallMoleculeAnalyte(**fields)


def peptide(**overrides) -> PeptideAnalyte:
    fields = dict(sequence="PEPTIDE", source=SOURCE, reuse_status=FIXTURE)
    fields.update(overrides)
    return PeptideAnalyte(**fields)


def glycan(**overrides) -> GlycanAnalyte:
    fields = dict(composition="Hex5HexNAc2", source=SOURCE, reuse_status=FIXTURE)
    fields.update(overrides)
    return GlycanAnalyte(**fields)


def protein(**overrides) -> ProteinAnalyte:
    fields = dict(accession="P01857", source=SOURCE, reuse_status=FIXTURE)
    fields.update(overrides)
    return ProteinAnalyte(**fields)


def antibody_identity(**overrides) -> AntibodyIdentity:
    fields = dict(inn="trastuzumab")
    fields.update(overrides)
    return AntibodyIdentity(**fields)


def antibody(**overrides) -> IntactAntibodyAnalyte:
    fields = dict(antibody=antibody_identity(), source=SOURCE, reuse_status=FIXTURE)
    fields.update(overrides)
    return IntactAntibodyAnalyte(**fields)


def adc(**overrides) -> ADCAnalyte:
    fields = dict(
        antibody=antibody_identity(),
        linker_payload_class="vc-MMAE",
        dar=2,
        source=SOURCE,
        reuse_status=FIXTURE,
    )
    fields.update(overrides)
    return ADCAnalyte(**fields)


def cyclic_settings(**overrides) -> CyclicSettings:
    """A single-pass cyclic run with wrap-around explicitly ruled out.

    Single pass and wrap_around=False on purpose: both are the states that do NOT
    block training, so a test that wants a blocker has to ask for one rather than
    getting it by accident from the default.
    """
    fields = dict(passes=1, pass_mode=PassMode.SINGLE_PASS, wrap_around=False)
    fields.update(overrides)
    return CyclicSettings(**fields)


def cyclic_measurement(**overrides) -> CCSMeasurement:
    """A cyclic measurement. Cyclic records MUST state their cyclic settings."""
    fields = dict(ims_type=IMSType.CYCLIC, calibrant="dextran", cyclic=cyclic_settings())
    fields.update(overrides)
    return measurement(**fields)


def measurement(**overrides) -> CCSMeasurement:
    """The minimal valid measurement: a TWIMS glycan value, everything else defaulted.

    TWIMS rather than DTIMS because TWIMS is calibrated, which exercises the
    calibrant rule in the direction a real record usually takes. A stepped-field
    DTIMS record has to drop the calibrant, so tests that want a primary value
    ask for one explicitly.
    """
    fields = dict(
        analyte=glycan(),
        adduct="[M+H]+",
        charge=1,
        polarity=Polarity.POSITIVE,
        ims_type=IMSType.TWIMS,
        drift_gas=DriftGas.N2,
        calibrant="dextran",
        ccs=300.0,
        source=SOURCE,
        reuse_status=FIXTURE,
    )
    fields.update(overrides)
    return CCSMeasurement(**fields)


def primary(**overrides) -> CCSMeasurement:
    """A stepped-field DTIMS measurement: the only primary kind, and it takes no calibrant."""
    from wmxccs.models import DTIMSMethod

    fields = dict(
        ims_type=IMSType.DTIMS,
        dtims_method=DTIMSMethod.STEPPED_FIELD,
        calibrant=None,
    )
    fields.update(overrides)
    return measurement(**fields)


@pytest.fixture
def glycan_measurement() -> CCSMeasurement:
    return measurement()


@pytest.fixture
def antibody_measurement() -> CCSMeasurement:
    """A native 24+ antibody ion. The charge is carried by protons, stated.

    An antibody whose charge carrier the source did not name is written
    "[M+24?]24+" instead; there are tests for both.
    """
    return measurement(
        analyte=antibody(folding_state=FoldingState.NATIVE),
        adduct="[M+24H]24+",
        charge=24,
        ccs=7000.0,
    )
