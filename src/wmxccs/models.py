"""The measurement record: one CCS value, the ion it was measured on, and where it came from.

Records are frozen, so a reuse status, calibrant or charge cannot change after
validation. They forbid unknown fields, so a misspelt field name is an error
rather than an input that is silently ignored. A copy made with
model_copy(update=...), or a record made with model_construct, skips validation;
the training gate validates such records again before trusting them.

Ported from the glycan platform's CCSMeasurement. What changed:

- `glycan: GlycanStructure` became `analyte: Analyte`, the tagged union in
  identity.py. Everything that assumed the analyte was a glycan went with it.
- the calibration group lost `reducing_end_label` and `derivatisation` as fields
  of its own and takes the analyte's structural state instead, which is those
  two for a glycan and the folding state for a protein.
- calibration reference lineage is NEW. See CalibrationReference below.

WHAT THIS MODULE DOES NOT DO
----------------------------
It never merges, averages or corrects a value. A harmonized CCS is a separate
output produced later, alongside the original and never in place of it. There is
no method here that returns a changed measurement, and adding one would be a bug
rather than a feature.
"""

from __future__ import annotations

import re
from datetime import date, datetime
from enum import StrEnum
from typing import Annotated, NamedTuple, Self

from pydantic import AfterValidator, BeforeValidator, Field, model_validator

from .identity import (
    UNDEFINED_GASES,
    Adduct,
    Analyte,
    DriftGas,
    Polarity,
    _Record,
    _Text,
    adduct_carrier_is_unstated,
    is_one_of,
    parse_adduct,
)
from .reuse import DEFAULT_REUSE_STATUS, ReuseStatus

# A DOI is a registered prefix "10." plus a registrant code, then a suffix.
# This shape is specified rather than inferred from examples, so a value that
# does not match is rejected rather than warned about.
_DOI = re.compile(r"10\.\d{4,9}/\S+")
_ISO_DATE = re.compile(r"[0-9]{4}-[0-9]{2}-[0-9]{2}")


def _check_doi(value: str) -> str:
    if not _DOI.fullmatch(value):
        raise ValueError(
            f"{value!r} is not a DOI: the form is '10.' then a registrant code, a slash and a suffix,"
            " as in 10.1038/s41467-021-00001-2. Record the bare DOI, not a URL or a citation"
        )
    return value


DOIField = Annotated[_Text, AfterValidator(_check_doi)]


def _date_as_written(value: object) -> object:
    """A date is a date object or the text YYYY-MM-DD, and nothing else.

    Left to pydantic, a run of digits - "1462060800", or a cell that held a
    count - is read as a Unix timestamp and becomes a plausible date with
    nothing to show it happened. A datetime is refused too: a time of day is a
    different claim from a date, and cutting it off would be a silent change.
    """
    if isinstance(value, datetime):
        raise ValueError("a measurement date is a date, not a datetime; give the date alone")
    if isinstance(value, date):
        return value
    if isinstance(value, str) and _ISO_DATE.fullmatch(value.strip()):
        return date.fromisoformat(value.strip())
    raise ValueError("a measurement date is written YYYY-MM-DD; a number or any other spelling is not a date")


_DateAsWritten = Annotated[date, BeforeValidator(_date_as_written)]


# --- the platform ----------------------------------------------------------------------


class IMSType(StrEnum):
    """The ion mobility technique the value was measured with.

    Deliberately NOT part of the matched-ion key. Pairing the same ion ACROSS
    these is the entire purpose of the platform.
    """

    DTIMS = "DTIMS"
    TWIMS = "TWIMS"
    TIMS = "TIMS"
    CYCLIC = "CYCLIC"


class DTIMSMethod(StrEnum):
    """How a drift-tube CCS was obtained."""

    # Drift times at several drift voltages; CCS follows from first
    # principles. Primary: no calibrant enters the value.
    STEPPED_FIELD = "stepped_field"
    # One drift voltage, calibrated against reference ions of known CCS.
    SINGLE_FIELD = "single_field"


class UncertaintyType(StrEnum):
    """What a reported uncertainty IS. It travels with the number, always.

    Papers report different things under the same column heading: one standard
    deviation, two, a standard error, a 95 per cent interval. A two-standard-
    deviation spread loaded into a field read as one standard deviation halves
    every interval built on it, and nothing downstream can detect that it
    happened. So a value may not be stored without its type, and "unknown" can
    be stored but blocks training, like the other unknowns.
    """

    SD = "sd"
    TWO_SD = "two_sd"
    SEM = "sem"
    CI95 = "ci95"
    UNKNOWN = "unknown"


TRAINABLE_UNCERTAINTY_TYPES = frozenset(UncertaintyType) - {UncertaintyType.UNKNOWN}


class CalibrationLineage(StrEnum):
    """Whether a CCS value stands on first principles or on somebody else's numbers."""

    # Stepped-field DTIMS. CCS follows from the Mason-Schamp relation and the
    # measured drift times; no reference value enters it.
    PRIMARY = "primary"
    # Everything else: TWIMS, TIMS, cyclic and single-field DTIMS all read their
    # CCS off a calibration curve built from reference ions of assumed CCS.
    DERIVED = "derived"


def _platform_label(ims_type: IMSType, dtims_method: DTIMSMethod | None) -> str:
    # str(), not .value, so a raw string left by model_copy still prints.
    return str(ims_type) if dtims_method is None else f"{ims_type}/{dtims_method}"


# --- cyclic ion mobility ----------------------------------------------------------------


class PassMode(StrEnum):
    """Whether the reported cyclic value came from one pass of the array or several."""

    SINGLE_PASS = "single_pass"
    MULTIPASS = "multipass"
    # The source reports a cyclic measurement without saying which. Spelt like
    # DriftGas.UNSTATED and for the same reason: a positive record of what the
    # paper does not say, not a filler. It blocks training, because a value that
    # does not say how far the ion travelled cannot be compared with one that does.
    UNSTATED = "UNSTATED"


UNDEFINED_PASS_MODES = frozenset({PassMode.UNSTATED})


class CyclicSettings(_Record):
    """What a cyclic IMS measurement did to the ion, held apart from single-pass TWIMS.

    A cyclic device sends ions round a closed path any number of times. More
    passes means more separation, and it also means a longer flight in which the
    fastest ions can catch and lap the slowest. So a cyclic value is not simply a
    TWIMS value from a fancier instrument: it carries its own conditions, and the
    ones below decide whether two cyclic values were produced the same way.

    PASS NUMBER IS PART OF THE CALIBRATION GROUP, NOT METADATA. A six-pass value
    and a single-pass value of the same ion are not interchangeable: the
    calibration that converts arrival time to CCS is specific to the path the ion
    took, and pooling them would average two quantities that were measured
    against different effective lengths. They are the SAME ION, so they share a
    matched-ion key - that is what lets them be compared - but they are not the
    same measurement, so they never pool.

    Nothing here is required beyond the pass count and the mode. The wave
    settings and the path length are recorded where a source states them and left
    null where it does not; guessing a travelling-wave velocity from what is
    usual would be inventing an instrument setting.
    """

    passes: int | None = Field(
        default=None,
        strict=True,
        ge=1,
        description="How many times round the array the reported value is for, counting from one."
        " Null where the source reports a multipass value without saying how many.",
    )
    pass_mode: PassMode = Field(
        default=PassMode.UNSTATED,
        description="Single-pass or multipass. Stated separately from the count because a source may report"
        " one without the other, and because it is what decides pooling when the count is null.",
    )
    effective_path_length_m: float | None = Field(
        default=None,
        strict=True,
        gt=0,
        allow_inf_nan=False,
        description="The path the ion actually travelled, in metres. Units are in the field name on purpose:"
        " a length recorded in centimetres and read as metres is the same class of silent error as a"
        " two-standard-deviation spread read as one.",
    )
    tw_velocity_m_per_s: float | None = Field(
        default=None, strict=True, gt=0, allow_inf_nan=False, description="Travelling-wave velocity, m/s."
    )
    tw_height_v: float | None = Field(
        default=None, strict=True, gt=0, allow_inf_nan=False, description="Travelling-wave height, volts."
    )
    arrival_time_correction: _Text | None = Field(
        default=None,
        description="How the arrival time was corrected for the time ions spend outside the separation"
        " region, as the source describes it. FREE TEXT DELIBERATELY, and it is not in the calibration"
        " group. It should be an enum, like every other controlled vocabulary here, and it is not one"
        " because no real cyclic dataset has been read in this repository yet. Inventing the members"
        " from what is usual would be exactly the kind of plausible detail this project refuses. When a"
        " real dataset arrives, read the methods and make this an enum from what is actually there.",
    )
    wrap_around: bool | None = Field(
        default=None,
        strict=True,
        description="Whether faster ions lapped slower ones, so that arrival time no longer orders mobility."
        " Null where the source does not say. On a multipass value that silence is a gap, not a no: the"
        " longer the flight, the likelier the lap.",
    )

    @model_validator(mode="after")
    def the_count_and_the_mode_agree(self) -> Self:
        if self.passes is None:
            return self
        expected = PassMode.SINGLE_PASS if self.passes == 1 else PassMode.MULTIPASS
        if self.pass_mode is not PassMode.UNSTATED and self.pass_mode is not expected:
            raise ValueError(
                f"{self.passes} pass(es) is {expected.value}, but pass_mode says {self.pass_mode.value};"
                " one of the two is wrong and guessing which would change what the value means"
            )
        return self

    @property
    def is_multipass(self) -> bool:
        """True where the value is known to be multipass, by count or by mode."""
        if self.passes is not None:
            return self.passes > 1
        return self.pass_mode is PassMode.MULTIPASS

    def training_blockers(self) -> list[str]:
        blockers = super().training_blockers()
        if is_one_of(getattr(self, "pass_mode", None), UNDEFINED_PASS_MODES) and self.passes is None:
            blockers.append(
                "a cyclic value that states neither a pass count nor a pass mode does not say how far the"
                " ion travelled, so it cannot be pooled with any other cyclic value"
            )
        if self.is_multipass and self.wrap_around is None:
            blockers.append(
                "a multipass cyclic value does not say whether wrap-around occurred; if faster ions lapped"
                " slower ones the arrival time no longer orders mobility, and the value may not be of the"
                " species it is assigned to. Read the source's methods rather than assuming it did not happen"
            )
        if self.wrap_around is True and self.arrival_time_correction is None:
            blockers.append(
                "wrap-around is reported and no arrival-time correction is recorded, so the reported value"
                " rests on an ordering that had already broken down"
            )
        return blockers


# --- calibration reference lineage ------------------------------------------------------


class CalibrationReference(_Record):
    """Where a derived CCS value's calibration curve came from.

    THIS IS THE FIELD THE OBJECTIVE DOCUMENT DOES NOT HAVE, AND ITS OWN RISK
    TABLE REQUIRES. The risk is reference-value circularity: a TWIMS value is
    calibrated against reference CCS values, those reference values were
    themselves usually measured by DTIMS, and if the platform then compares that
    TWIMS value against the same DTIMS values it calibrated from, the agreement
    it reports is partly an artefact of the calibration rather than a fact about
    the instruments. A pipeline cannot detect that without knowing WHICH
    reference set a value came from, and `calibrant` alone does not say: two
    laboratories both writing "dextran" may be using different ladders with
    different assumed values.

    So the reference set is recorded by name, by the publication its values were
    published in where there is one, and by the platform those values were
    themselves measured on. `traces_to_primary` then answers, per record,
    whether the chain reaches first principles or stops at another calibration.

    Nothing here is inferred. A record whose source does not say which reference
    set was used records the name it gave and leaves the DOI null; it does not
    guess at the usual one.
    """

    reference_set: _Text = Field(
        description="What the calibration curve was built from, as the source names it, e.g."
        " 'dextran DTCCS ladder' or 'Bush denatured protein He CCS values'.",
    )
    doi: DOIField | None = Field(
        default=None,
        description="The publication those reference values were taken from, where the source names one."
        " This is what makes circularity detectable: a value calibrated from a paper that is also in the"
        " comparison set is not independent of it.",
    )
    platform: IMSType | None = Field(
        default=None,
        description="The technique the REFERENCE values were themselves measured on. Null where the source"
        " does not say. A reference set measured by TWIMS is itself calibrated, so a value calibrated"
        " against it does not trace to first principles however many steps are in between.",
    )
    method: DTIMSMethod | None = Field(
        default=None,
        description="Where the reference set is DTIMS, whether it was stepped-field (primary) or single-field"
        " (itself calibrated). Allowed only for DTIMS, exactly as on a measurement.",
    )

    @model_validator(mode="after")
    def method_belongs_to_dtims(self) -> Self:
        if self.method is not None and self.platform is not IMSType.DTIMS:
            raise ValueError(
                f"dtims_method applies only to DTIMS, and the reference set is {self.platform or 'unstated'}"
            )
        return self

    @property
    def traces_to_primary(self) -> bool | None:
        """Whether this reference set is itself primary. None where the source does not say enough.

        True only for stepped-field DTIMS. False for any other stated platform,
        including single-field DTIMS, which is calibrated like the rest. None
        where the platform is not stated or a DTIMS reference does not say which
        method - which is a gap to be reported, not an assumption to be made.
        """
        if self.platform is None:
            return None
        if self.platform is not IMSType.DTIMS:
            return False
        if self.method is None:
            return None
        return self.method is DTIMSMethod.STEPPED_FIELD


# --- the conditions that fix what a value means -----------------------------------------


class CalibrationGroup(NamedTuple):
    """Conditions that must match before CCS values are pooled WITHIN one platform.

    NOT the matched-ion key, and close to its complement. This carries the
    platform, the method and the calibrant and carries no analyte identity; the
    matched-ion key carries the analyte and no platform. They exist for opposite
    purposes: this one says "these values were produced the same way", the
    matched-ion key says "these values are of the same ion". Using either in
    place of the other produces a check that can never fire, and a check that
    can never fire passes every test.

    Necessary, not sufficient: one calibrant name can still hide different
    reference values or calibration procedures, which is what
    CalibrationReference exists to record.
    """

    ims_type: IMSType
    dtims_method: DTIMSMethod | None
    drift_gas: DriftGas
    calibrant: str | None
    adduct: str
    structural_state: tuple
    # How far the ion travelled, for a cyclic measurement, and None for every
    # other platform. IN THE GROUP, not metadata: a six-pass value and a
    # single-pass value of one ion were calibrated against different effective
    # path lengths and do not pool. They still share a matched-ion KEY, which is
    # what lets them be compared; this is what stops them being averaged.
    cyclic_passes: tuple | None = None

    def __str__(self) -> str:
        # "-" cannot collide with a real calibrant: it is rejected as a placeholder.
        state = ",".join("" if part is None else str(part) for part in self.structural_state)
        parts = [
            _platform_label(self.ims_type, self.dtims_method),
            self.drift_gas,
            self.calibrant or "-",
            self.adduct,
            state,
        ]
        if self.cyclic_passes is not None:
            passes, mode = self.cyclic_passes
            parts.append(f"{'?' if passes is None else passes} pass/{mode}")
        return "|".join(str(part) for part in parts)


class MeasurementConditions(_Record):
    """The conditions that fix what a CCS value means: ion, platform, reference gas and calibrant.

    Shared by measurements and, later, by harmonization queries: a query about
    an ion has to be typed by the same conditions as a stored record, or the two
    cannot be compared. A calibrant is required exactly when one enters the CCS
    value: TWIMS, TIMS, cyclic and single-field DTIMS. Any platform added later
    counts as calibrated until listed otherwise.
    """

    adduct: Adduct
    charge: int = Field(strict=True, description="Signed charge state, e.g. 2 or -1.")
    polarity: Polarity
    ims_type: IMSType
    dtims_method: DTIMSMethod | None = Field(
        default=None, description="Required for DTIMS, and allowed only for DTIMS."
    )
    drift_gas: DriftGas = Field(
        description="The gas the reported CCS value refers to, not necessarily the gas in the cell. A TWIMS"
        " value measured in an N2 cell but calibrated against He reference values is a He value, and does not"
        " pool with N2 values. A measurement records the gas in the cell separately, as cell_gas."
    )
    calibrant: _Text | None = Field(
        default=None,
        description="Required when a calibrant enters the CCS value (TWIMS, TIMS, cyclic, single-field DTIMS);"
        " not allowed for stepped-field DTIMS.",
    )

    @property
    def ccs_is_calibrated(self) -> bool:
        """True if a calibrant enters the CCS value."""
        if self.ims_type is IMSType.DTIMS:
            return self.dtims_method is not DTIMSMethod.STEPPED_FIELD
        return True

    @property
    def calibration_lineage(self) -> CalibrationLineage:
        """Whether this value stands on first principles or on a calibration curve."""
        return CalibrationLineage.DERIVED if self.ccs_is_calibrated else CalibrationLineage.PRIMARY

    @model_validator(mode="after")
    def check_conditions(self) -> Self:
        if self.ims_type is IMSType.DTIMS and self.dtims_method is None:
            raise ValueError(
                "DTIMS records must state dtims_method (stepped_field or single_field):"
                " single-field CCS is calibrated, stepped-field CCS is primary"
            )
        if self.ims_type is not IMSType.DTIMS and self.dtims_method is not None:
            raise ValueError(f"dtims_method applies only to DTIMS, not {self.ims_type.value}")
        platform = _platform_label(self.ims_type, self.dtims_method)
        if self.ccs_is_calibrated and self.calibrant is None:
            raise ValueError(
                f"{platform} gives calibrated CCS, so the calibrant is part of the measurement and changes"
                " the value; a record without one is rejected"
            )
        if not self.ccs_is_calibrated and self.calibrant is not None:
            raise ValueError(
                f"{platform} gives primary CCS: no calibrant enters the value, so none may be recorded"
                " (a standard run only as a check is not a calibrant)"
            )
        if self.charge == 0:
            raise ValueError("charge cannot be 0: ion mobility measures ions")
        if (self.charge > 0) != (self.polarity is Polarity.POSITIVE):
            raise ValueError(f"charge {self.charge:+d} does not match {self.polarity.value} polarity")
        _, adduct_charge = parse_adduct(self.adduct)
        if adduct_charge != self.charge:
            raise ValueError(f"adduct {self.adduct} carries charge {adduct_charge:+d}, but charge is {self.charge:+d}")
        return self


class CCSMeasurement(MeasurementConditions):
    """One experimental CCS value, the analyte it was measured on, and where it came from."""

    analyte: Analyte
    ccs: float = Field(strict=True, gt=0, allow_inf_nan=False, description="Collision cross section in square angstrom.")
    ccs_uncertainty: float | None = Field(
        default=None,
        strict=True,
        gt=0,
        allow_inf_nan=False,
        description="The reported spread of the CCS value, in square angstrom. Meaningless without uncertainty_type,"
        " which must accompany it. Null where the source reports none.",
    )
    uncertainty_type: UncertaintyType | None = Field(
        default=None,
        description="What ccs_uncertainty is: one standard deviation, two, a standard error, or a 95 per cent"
        " interval. Required whenever ccs_uncertainty is given. 'unknown' can be stored but blocks training.",
    )
    calibration_reference: CalibrationReference | None = Field(
        default=None,
        description="Where the calibration curve came from, for a derived value. Not allowed where no calibrant"
        " enters the value. Null where the source does not say which reference set was used, which is a gap"
        " worth seeing rather than a value worth guessing.",
    )
    conformer: int | None = Field(
        default=None,
        strict=True,
        ge=1,
        description="Which resolved conformer this value is, counting from one, where one structure gives more"
        " than one arrival-time peak under a single set of conditions. Null where the source reports one peak."
        " Two rows differing only in this are a legitimate pair, not a duplicate and not a collision.",
    )
    conformers_total: int | None = Field(
        default=None,
        strict=True,
        ge=1,
        description="How many conformers the source resolved for this structure under these conditions. More"
        " than one is a claim that the other rows exist: the loader holds a row whose siblings are missing,"
        " rather than letting a value stand as though it were the only peak.",
    )
    cell_gas: DriftGas | None = Field(
        default=None,
        description="The gas that was in the ion mobility cell, recorded for provenance only; drift_gas, the gas"
        " the CCS value refers to, decides pooling. Null if not reported. For a primary (stepped-field DTIMS)"
        " value the two must be the same gas.",
    )
    cyclic: CyclicSettings | None = Field(
        default=None,
        description="Cyclic-specific conditions. Required for a CYCLIC record and refused on every other"
        " platform, so a pass count can never be recorded against an instrument that has no passes.",
    )
    instrument: _Text | None = None
    source: _Text
    # Provenance per value. `source` is free text and cannot be verified to be a
    # citation, so these four carry the parts a later milestone has to be able to
    # act on: a DOI to group by study without scraping prose, a locator to find
    # the value again in the paper, a replicate count to tell a single shot from
    # a mean of five, and a date, without which no chronological split is
    # possible at all. Null where a source did not report them, never a filler.
    doi: DOIField | None = Field(
        default=None,
        description="The DOI of the publication the value came from, bare and without a URL prefix."
        " Used to group measurements by study.",
    )
    source_locator: _Text | None = Field(
        default=None,
        description="Where in that source the value sits: a table, figure or supplementary file, as in"
        " 'Table S2' or 'Supplementary Data 3'. A DOI alone does not let anyone find the number again.",
    )
    replicates: int | None = Field(
        default=None,
        strict=True,
        ge=1,
        description="How many independent measurements this value is drawn from. Null where the source did"
        " not say, which is not the same as one.",
    )
    measured_on: _DateAsWritten | None = Field(
        default=None,
        description="The date the measurement was made, where the source reports it, written YYYY-MM-DD. Required"
        " before any chronological split can be built; a publication date is not a measurement date.",
    )
    reuse_status: ReuseStatus = DEFAULT_REUSE_STATUS

    @model_validator(mode="after")
    def check_measurement(self) -> Self:
        if self.cell_gas is not None and not self.ccs_is_calibrated and self.cell_gas is not self.drift_gas:
            raise ValueError(
                f"a primary CCS refers to the gas it was measured in, so cell_gas {self.cell_gas.value} must"
                f" equal drift_gas {self.drift_gas.value}"
            )
        # The type travels with the number, in both directions: a spread with no
        # stated type would be read as whatever the reader assumed, and a type with
        # no spread describes nothing.
        if self.ccs_uncertainty is not None and self.uncertainty_type is None:
            raise ValueError(
                "ccs_uncertainty needs uncertainty_type: a spread of unknown kind is read as whatever the reader"
                " assumes, so say what it is, or 'unknown' if the source does not"
            )
        # "unknown" is the one type that may stand alone. It does not describe a
        # number; it records that the source states no spread, which is a fact
        # worth keeping and which blocks training on its own. Every other type
        # describes a number, and describing a number that is not there says
        # nothing.
        if (
            self.uncertainty_type is not None
            and self.uncertainty_type is not UncertaintyType.UNKNOWN
            and self.ccs_uncertainty is None
        ):
            raise ValueError(
                f"uncertainty_type {self.uncertainty_type.value!r} describes ccs_uncertainty, which is not"
                " given; only 'unknown' may stand alone"
            )
        if self.conformer is not None and self.conformers_total is None:
            raise ValueError(
                "a conformer index says nothing without conformers_total: state how many the source resolved"
            )
        if self.conformers_total is not None and self.conformers_total > 1 and self.conformer is None:
            raise ValueError(
                f"conformers_total is {self.conformers_total}, so this row must say which of them it is"
            )
        if self.conformer is not None and self.conformers_total is not None and self.conformer > self.conformers_total:
            raise ValueError(f"conformer {self.conformer} of {self.conformers_total}: the index exceeds the total")
        if self.ims_type is IMSType.CYCLIC and self.cyclic is None:
            raise ValueError(
                "a cyclic record must state its cyclic settings: how far the ion travelled is part of what"
                " the value means, and a cyclic value without a pass count cannot be told apart from a"
                " single-pass travelling-wave value"
            )
        if self.ims_type is not IMSType.CYCLIC and self.cyclic is not None:
            raise ValueError(f"cyclic settings apply only to CYCLIC records, not {self.ims_type.value}")
        if self.calibration_reference is not None and not self.ccs_is_calibrated:
            raise ValueError(
                "a primary CCS has no calibrant, so it cannot have a calibration reference"
                " (a standard run only as a check is not a calibrant)"
            )
        return self

    @property
    def matched_ion_key(self):
        """The key that decides whether this value may be compared with another.

        Built from the analyte's identity and structural state, the canonical
        adduct, the signed charge and the gas the value REFERS TO. It carries no
        platform, no calibrant, no instrument and no DOI: a DTIMS value and a
        TWIMS value of one ion land on one key, which is the whole point.
        """
        from .identity import MatchedIonKey

        # An ion whose charge carrier the source never named cannot be matched
        # with anything, including another record of the same shape, so its key
        # carries its own provenance and is unique to it. See MatchedIonKey.
        # Derived from the record rather than generated, so the key is stable
        # across runs and a report is reproducible.
        unmatchable = None
        if adduct_carrier_is_unstated(self.adduct):
            unmatchable = (
                "charge carrier not stated",
                self.source,
                self.doi,
                self.source_locator,
                self.ccs,
                self.conformer,
            )

        return MatchedIonKey(
            analyte=self.analyte.identity_key(),
            adduct=self.adduct,
            charge=self.charge,
            drift_gas=self.drift_gas,
            state=self.analyte.structural_state(),
            unmatchable=unmatchable,
        )

    @property
    def calibration_group(self) -> CalibrationGroup:
        return CalibrationGroup(
            self.ims_type,
            self.dtims_method,
            self.drift_gas,
            self.calibrant,
            self.adduct,
            self.analyte.structural_state(),
            None if self.cyclic is None else (self.cyclic.passes, str(self.cyclic.pass_mode)),
        )

    @property
    def traces_to_primary(self) -> bool | None:
        """Whether this value's calibration chain reaches first principles.

        True for a primary value itself, and for a derived value whose reference
        set is stepped-field DTIMS. False where the reference set is itself
        calibrated. None where the value is derived and the reference set is
        unrecorded or does not say what it was measured on - the honest answer
        being that nobody here knows, which is the state most published TWIMS
        values are actually in.
        """
        if not self.ccs_is_calibrated:
            return True
        if self.calibration_reference is None:
            return None
        return self.calibration_reference.traces_to_primary

    def component_records(self) -> tuple:
        """Records this one is built from. The training gate checks them as well.

        The cyclic settings are NOT here, deliberately. They carry no licence of
        their own - they are conditions of this measurement, not a separately
        sourced record - and the gate refuses anything in this tuple that cannot
        state a reuse status. That mistake closed the whole biopharmaceutical
        layer once already, with AntibodyIdentity. Their training blockers are
        folded into this record's own instead.
        """
        parts: list = [self.analyte]
        nested = getattr(self.analyte, "component_records", None)
        if nested is not None:
            parts.extend(nested())
        return tuple(parts)

    def training_blockers(self) -> list[str]:
        """Reasons other than licence that this record may not enter training."""
        blockers = super().training_blockers()
        # Read defensively: a record built with model_construct may have no
        # drift_gas at all, and an AttributeError here would be caught upstream
        # and reported as "blockers could not be checked", hiding the specific
        # reasons this method exists to give.
        if is_one_of(getattr(self, "drift_gas", None), UNDEFINED_GASES):
            blockers.append(
                "drift_gas is 'UNSTATED', so the gas this CCS refers to is not defined: it cannot be pooled"
                " with helium or with nitrogen values, and resolving it means reading the calibration"
                " reference the source cites, not inferring it from the same group's other work"
            )
        adduct = getattr(self, "adduct", None)
        if isinstance(adduct, str):
            try:
                unstated_carrier = adduct_carrier_is_unstated(adduct)
            except ValueError:  # an unvalidated copy; revalidation reports it
                unstated_carrier = False
            if unstated_carrier:
                blockers.append(
                    f"adduct {adduct} does not name the charge carrier, so the ion is not defined: its mass"
                    " depends on whether those charges are protons or sodium, and so does its cross section."
                    " Resolving it means reading the source's methods, never assuming protons"
                )
        cyclic = getattr(self, "cyclic", None)
        if cyclic is not None:
            try:
                blockers.extend(f"cyclic settings: {blocker}" for blocker in cyclic.training_blockers())
            except Exception as exc:  # fail closed, as everywhere else in this method
                blockers.append(f"cyclic settings could not be checked ({type(exc).__name__}: {exc})")
        kind = self.uncertainty_type
        if kind is not None and not is_one_of(kind, TRAINABLE_UNCERTAINTY_TYPES):
            blockers.append(
                f"uncertainty_type is {str(kind)!r}, so the reported spread cannot be interpreted; an interval built"
                " on it would be wrong by an unknown factor"
            )
        return blockers

    @property
    def validation_warnings(self) -> list[str]:
        """Problems that do not block the record but should be reviewed, including its analyte's."""
        return [f"analyte: {warning}" for warning in self.analyte.validation_warnings]
