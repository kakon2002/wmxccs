"""Pydantic records for glycan structures and CCS measurements.

Records are frozen, so a reuse status, calibrant or charge cannot change after
validation. They forbid unknown fields, so a misspelt field name is an error
rather than an input that is silently ignored. A copy made with
model_copy(update=...), or a record made with model_construct, skips
validation; the training gate validates such records again before trusting them.
"""

import re
import warnings
from collections import Counter
from datetime import date, datetime
from enum import StrEnum
from typing import Annotated, NamedTuple, Self

from pydantic import (
    AfterValidator,
    BaseModel,
    BeforeValidator,
    ConfigDict,
    Field,
    PlainSerializer,
    PlainValidator,
    StrictBool,
    StringConstraints,
    ValidationError,
    WithJsonSchema,
    model_validator,
)

from .composition import Composition, parse_composition
from .reuse import DEFAULT_REUSE_STATUS, ReuseStatus  # from below the gate: the gate imports the registry, which imports this


class UnverifiedFormatWarning(UserWarning):
    """A value does not match a format rule that has not itself been verified yet.

    A UserWarning, so Python shows it by default.
    """


def _to_composition(value: object) -> Composition:
    if isinstance(value, Composition):
        return value
    if isinstance(value, str):
        return parse_composition(value)  # CompositionError is a ValueError, so pydantic reports it
    raise ValueError(f"composition must be a string such as 'Hex5HexNAc4Fuc1', not {type(value).__name__}")


# A Composition in Python, its canonical string in JSON.
CompositionField = Annotated[
    Composition,
    PlainValidator(_to_composition),
    PlainSerializer(lambda comp: comp.canonical, return_type=str),
    WithJsonSchema(
        {
            "type": "string",
            "description": "Glycan composition as <Residue><count> tokens in any order; dHex is read as Fuc.",
            "examples": ["Hex5HexNAc4Fuc1"],
        }
    ),
]

# Strings that stand in for a missing value, as found in published tables and
# written by spreadsheet and dataframe tools. They are rejected wherever text
# is required, and so is any string with no letter or digit in it (dashes of
# any kind, "?", "..."): a genuinely unknown value is null, in fields that
# allow null. Only these spellings are caught, so a loader should still map
# its own fillers to null. No enum value may be spelt like a placeholder unless
# it blocks training, so a filler can never land on a value that trains;
# "unknown" is the one such spelling, and it blocks training.
_PLACEHOLDERS = frozenset(
    {
        "-",
        "?",
        "#n/a",
        "<na>",
        "missing",
        "n.a.",
        "n.d.",
        "n/a",
        "na",
        "nan",
        "nat",
        "nd",
        "none",
        "not available",
        "not determined",
        "not reported",
        "null",
        "tbd",
        "unknown",
    }
)


def _is_placeholder(value: str) -> bool:
    """True for a listed placeholder spelling, or for a string with no letter or digit in it."""
    return value.strip().casefold() in _PLACEHOLDERS or not any(ch.isalnum() for ch in value)


def _not_placeholder(value: str) -> str:
    if _is_placeholder(value):
        raise ValueError(f"{value!r} is a placeholder, not a value; use null where the field allows it")
    return value


_Text = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1), AfterValidator(_not_placeholder)]

# The form G + five digits + two capital letters comes from accessions seen in
# use, not from GlyTouCan's documentation. Until it is verified, a mismatch
# warns instead of rejecting.
_GLYTOUCAN_AC = re.compile(r"G[0-9]{5}[A-Z]{2}")

# An IUPAC-condensed structure states its linkages in parentheses, as in
# "Man(a1-3)Man". A string with no linkage in it says nothing about how any
# residues are joined, so it cannot be the evidence behind a claim that the
# linkages are resolved, however structure-shaped it looks.
_IUPAC_LINKAGE = re.compile(r"\([ab?]?[0-9?]+-[0-9?]+(?:/[0-9?]+)*\)")


def _glytoucan_format_problem(value: str) -> str | None:
    if _GLYTOUCAN_AC.fullmatch(value):
        return None
    return (
        f"GlyTouCan accession {value!r} does not match the expected form (G, five digits, two capital"
        " letters); kept, because that form is not yet verified"
    )


def _warn_on_glytoucan_format(value: str) -> str:
    problem = _glytoucan_format_problem(value)
    if problem is not None:
        warnings.warn(problem, UnverifiedFormatWarning)
    return value


# A DOI is a registered prefix "10." plus a registrant code, then a suffix.
# Unlike the GlyTouCan form, this shape is specified rather than inferred from
# examples, so a value that does not match is rejected rather than warned about.
_DOI = re.compile(r"10\.\d{4,9}/\S+")

_ISO_DATE = re.compile(r"[0-9]{4}-[0-9]{2}-[0-9]{2}")


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


def _check_doi(value: str) -> str:
    if not _DOI.fullmatch(value):
        raise ValueError(
            f"{value!r} is not a DOI: the form is '10.' then a registrant code, a slash and a suffix,"
            " as in 10.1038/s41467-021-00001-2. Record the bare DOI, not a URL or a citation"
        )
    return value


def _check_wurcs(value: str) -> str:
    if not value.startswith("WURCS=") or any(ch.isspace() for ch in value):
        raise ValueError("WURCS must start with 'WURCS=' and contain no whitespace")
    return value


# Bracket notation with the charge after the bracket: [M+H]+, [M+2H]2+, [M-H]-, [2M+Na]+.
_ADDUCT = re.compile(r"\[(?P<mult>[1-9][0-9]*)?M(?P<parts>[^\[\]\s]*)\](?P<n>[1-9][0-9]*)?(?P<sign>[+-])")
# One component of the ion: a sign, an optional count and a species such as H, Na, NH4, HCOO or H2O.
_ADDUCT_PART = re.compile(r"([+-])([1-9][0-9]*)?([A-Z][A-Za-z0-9]*)")
_MAX_ADDUCT_LENGTH = 100  # bound for strings arriving through the API; real adducts are far shorter


def _parse_adduct(adduct: str) -> tuple[str, int]:
    """Canonical spelling and signed charge of an adduct.

    Components are sorted, additions before losses and then by species, so
    "[M+Na+H]2+" and "[M+H+Na]2+" give the same key. Counts, multipliers and
    charges of 1 are dropped: "[M+1H]1+" -> ("[M+H]+", 1). A species named
    twice is rejected rather than merged.
    """
    match = _ADDUCT.fullmatch(adduct)
    if match is None:
        raise ValueError(f"adduct {adduct!r} is not in bracket notation such as [M+H]+, [M+2H]2+ or [M-H]-")
    body = match["parts"]
    parts: list[tuple[str, str, str]] = []
    pos = 0
    while pos < len(body):
        part = _ADDUCT_PART.match(body, pos)
        if part is None:
            raise ValueError(
                f"adduct {adduct!r}: cannot read {body[pos:]!r}; expected components such as +H, +2Na or -H2O"
            )
        sign, count, species = part.groups()
        parts.append((sign, "" if count in (None, "1") else count, species))
        pos = part.end()
    named = Counter(species for _, _, species in parts)
    repeated = sorted(name for name, times in named.items() if times > 1)
    if repeated:
        raise ValueError(
            f"adduct {adduct!r} names {', '.join(repeated)} more than once; write each species once, with its count"
        )
    parts.sort(key=lambda part: (part[0] != "+", part[2]))
    multiplier = "" if match["mult"] in (None, "1") else match["mult"]
    ion = multiplier + "M" + "".join(sign + count + species for sign, count, species in parts)
    magnitude = int(match["n"] or 1)
    canonical = f"[{ion}]{magnitude if magnitude > 1 else ''}{match['sign']}"
    return canonical, magnitude if match["sign"] == "+" else -magnitude


def adduct_components(adduct: str) -> tuple[tuple[str, int, str], ...]:
    """The (sign, count, species) parts of an adduct: "[M+Na+H]2+" -> (("+", 1, "H"), ("+", 1, "Na")).

    _parse_adduct keeps only the canonical spelling and the signed charge, which
    is all a record needs in order to be well formed. Anything reading the ion
    itself needs the parts: a sodiated ion is not a protonated one of the same
    charge, and the two differ in mass and in collision cross section. Parts come
    back in canonical order, additions before losses and then by species, so two
    spellings of one ion give one sequence.

    Raises ValueError for an adduct that is not in bracket notation.
    """
    canonical, _ = _parse_adduct(adduct)  # validates first; a malformed adduct raises here
    match = _ADDUCT.fullmatch(canonical)
    assert match is not None  # _parse_adduct built this string, so it matches by construction
    body = match["parts"]
    parts: list[tuple[str, int, str]] = []
    pos = 0
    while pos < len(body):
        part = _ADDUCT_PART.match(body, pos)
        assert part is not None  # same reason: canonical text, already parsed once
        sign, count, species = part.groups()
        parts.append((sign, int(count) if count else 1, species))
        pos = part.end()
    return tuple(parts)


_Adduct = Annotated[
    str,
    StringConstraints(strip_whitespace=True, max_length=_MAX_ADDUCT_LENGTH),
    AfterValidator(lambda value: _parse_adduct(value)[0]),
]


class Polarity(StrEnum):
    POSITIVE = "positive"
    NEGATIVE = "negative"


class ReducingEndLabel(StrEnum):
    """The label chemistry at the reducing end of the measured glycan.

    A labelled glycan is a different analyte from the native one, with its own
    mass and CCS, so the label is part of the calibration group. "other" and
    "unknown" can be stored but block training: an open bucket would pool
    unrelated labels in one conformal bin. A real label missing from this list
    is added to it, not recorded as "other".
    """

    NATIVE = "native"
    REDUCED = "reduced"
    PA = "PA"
    TWO_AB = "2-AB"
    TWO_AA = "2-AA"
    PROCAINAMIDE = "procainamide"
    RAPIFLUOR_MS = "RapiFluor-MS"
    APTS = "APTS"
    OTHER = "other"
    UNKNOWN = "unknown"


class Derivatisation(StrEnum):
    """Whole-molecule modification of the measured glycan, independent of the reducing-end label.

    Part of the calibration group. The value describes the measured analyte,
    not the sample protocol: "underivatised" means nothing was modified apart
    from the reducing-end label, and a derivatisation that was not reported is
    "unknown", which can be stored but blocks training.

    The two sialic-acid values are open buckets, like an "other" label: each
    covers chemistries of different mass, and a linkage-specific protocol that
    treats alpha2,3- and alpha2,6-linked sialic acids differently has no
    correct value. They can be stored but block training, and they need a
    composition with NeuAc or NeuGc. A glycan with no sialic acid from such a
    sample is "underivatised" only if nothing else on it reacted: the carboxyl
    group of a 2-AA label can react with the same reagents, so a neutral 2-AA
    glycan from such a sample is "unknown". A specific value is added here
    when a real dataset needs one.
    """

    UNDERIVATISED = "underivatised"
    PERMETHYLATION = "permethylation"
    SIALIC_ACID_AMIDATION = "sialic_acid_amidation"
    SIALIC_ACID_ESTERIFICATION = "sialic_acid_esterification"
    UNKNOWN = "unknown"


# What defines the measured analyte well enough to train on. Membership, not
# identity, so a raw string that skipped validation is judged by its value.
_TRAINABLE_LABELS = frozenset(ReducingEndLabel) - {ReducingEndLabel.OTHER, ReducingEndLabel.UNKNOWN}
_OPEN_LABELS = frozenset({ReducingEndLabel.OTHER})
_SIALIC_ACID_DERIVATISATIONS = frozenset(
    {Derivatisation.SIALIC_ACID_AMIDATION, Derivatisation.SIALIC_ACID_ESTERIFICATION}
)
_TRAINABLE_DERIVATISATIONS = frozenset(Derivatisation) - {Derivatisation.UNKNOWN} - _SIALIC_ACID_DERIVATISATIONS


def _is_one_of(value: object, allowed: frozenset) -> bool:
    try:
        return value in allowed
    except (TypeError, ValueError):  # unhashable, or a comparison with no truth value, after skipped validation
        return False


def check_derivatisation_fits(composition: Composition, derivatisation: Derivatisation) -> None:
    """Raise ValueError if a sialic-acid derivatisation is claimed for a composition with no sialic acid."""
    if derivatisation in _SIALIC_ACID_DERIVATISATIONS and composition.neuac + composition.neugc == 0:
        raise ValueError(
            f"derivatisation {str(derivatisation)!r} needs a sialic acid to act on, and {composition.canonical}"
            " has no NeuAc or NeuGc. Record what was done to this analyte: 'underivatised' if nothing on it was"
            " modified, 'unknown' if something else may have reacted, such as the carboxyl of a 2-AA label"
        )


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


_TRAINABLE_UNCERTAINTY_TYPES = frozenset(UncertaintyType) - {UncertaintyType.UNKNOWN}


class IMSType(StrEnum):
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


class DriftGas(StrEnum):
    """A gas used for ion mobility.

    drift_gas is the gas a CCS value refers to; a measurement's cell_gas is
    the gas that was in the cell. They differ when, for example, a TWIMS value
    from an N2 cell is calibrated against He reference values.
    """

    N2 = "N2"
    HE = "He"
    # The source did not say which gas its values refer to. A positive fact, not
    # a filler: it records what the paper does not state, so the gap is visible
    # rather than guessed at, and it is deliberately not spelt "unknown", which
    # is a placeholder this package refuses everywhere else. It blocks training,
    # because a CCS whose gas is undefined cannot be pooled with helium or with
    # nitrogen values, and it keeps a calibration group of its own for the same
    # reason. Resolving it means reading the calibration reference a paper
    # cites, never inferring from what the same group used elsewhere.
    UNSTATED = "UNSTATED"


# A drift gas that does not define the quantity the value refers to.
_UNDEFINED_GASES = frozenset({DriftGas.UNSTATED})


def _platform_label(ims_type: IMSType, dtims_method: DTIMSMethod | None) -> str:
    # str(), not .value, so a raw string left by model_copy still prints.
    return str(ims_type) if dtims_method is None else f"{ims_type}/{dtims_method}"


class CalibrationGroup(NamedTuple):
    """Conditions that must match before CCS values are pooled. Reused later as the conformal bin key.

    Necessary, not sufficient: one calibrant name can still hide different
    reference values or calibration procedures.
    """

    ims_type: IMSType
    dtims_method: DTIMSMethod | None
    drift_gas: DriftGas
    calibrant: str | None
    adduct: str
    reducing_end_label: ReducingEndLabel
    derivatisation: Derivatisation

    def __str__(self) -> str:
        # "-" cannot collide with a real calibrant: it is rejected as a placeholder.
        parts = (
            _platform_label(self.ims_type, self.dtims_method),
            self.drift_gas,
            self.calibrant or "-",
            self.adduct,
            self.reducing_end_label,
            self.derivatisation,
        )
        return "|".join(str(part) for part in parts)


def _revalidation_problem(record: BaseModel) -> str | None:
    """Why `record` cannot be trusted as validated, or None if it can.

    It must pass validation, and equal its validated form: a copy with an
    unsorted adduct or unstripped text passes validation but is not what
    validation would have stored.
    """
    advice = "build it through the model rather than model_construct or model_copy"
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # validation_warnings carries these; repeating them on every check is noise
            dumped = record.model_dump()
            validated = type(record).model_validate(dumped)
            unchanged = validated.model_dump() == dumped
    except ValidationError as exc:
        first = exc.errors()[0]
        where = ".".join(str(part) for part in first["loc"]) or "record"
        return f"it does not pass validation ({where}: {first['msg']}); {advice}"
    except Exception as exc:  # contents that cannot even be serialised for checking
        return f"it cannot be re-validated ({type(exc).__name__}: {exc}); {advice}"
    if not unchanged:
        return f"it differs from its validated form, so it was built or changed without validation; {advice}"
    return None


class _Record(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    def training_blockers(self) -> list[str]:
        """Reasons other than licence that this record may not enter training.

        The gate does not trust validation it cannot see, so a record built or
        changed without validation is validated again here.
        """
        problem = _revalidation_problem(self)
        return [] if problem is None else [problem]


class GlycanStructure(_Record):
    """A glycan, as far as it is actually known.

    Linkage and anomericity are unresolved unless the record says otherwise,
    and only a record with a structure identifier can say otherwise: a WURCS,
    a GlyTouCan accession in the expected form, or the IUPAC-condensed string
    the structure was read from.
    """

    composition: CompositionField
    wurcs: Annotated[str, StringConstraints(strip_whitespace=True), AfterValidator(_check_wurcs)] | None = None
    glytoucan_ac: Annotated[_Text, AfterValidator(_warn_on_glytoucan_format)] | None = None
    iupac_condensed: Annotated[_Text, StringConstraints(pattern=r"^\S+$")] | None = Field(
        default=None,
        description="The structure as an IUPAC-condensed string, e.g. Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc."
        " It states the linkages, so it counts as a structure identifier.",
    )
    has_unresolved_linkage: StrictBool = True
    has_unresolved_anomericity: StrictBool = True
    source: _Text
    reuse_status: ReuseStatus = DEFAULT_REUSE_STATUS

    @model_validator(mode="after")
    def composition_only_stays_unresolved(self) -> Self:
        resolved = not (self.has_unresolved_linkage and self.has_unresolved_anomericity)
        well_formed_accession = self.glytoucan_ac is not None and _glytoucan_format_problem(self.glytoucan_ac) is None
        states_linkages = self.iupac_condensed is not None and _IUPAC_LINKAGE.search(self.iupac_condensed) is not None
        identified = self.wurcs is not None or well_formed_accession or states_linkages
        if resolved and not identified:
            raise ValueError(
                "a record with only a composition cannot claim resolved linkage or anomericity: it has no WURCS,"
                " no IUPAC-condensed structure stating its linkages, and no GlyTouCan accession in the expected"
                " form (neither an accession that does not match the form, nor a structure string with no"
                " linkage in it, counts as a structure identifier)"
            )
        return self

    @property
    def validation_warnings(self) -> list[str]:
        """Problems that do not block the record but should be reviewed."""
        if self.glytoucan_ac is None:
            return []
        problem = _glytoucan_format_problem(self.glytoucan_ac)
        return [] if problem is None else [problem]


class MeasurementConditions(_Record):
    """The conditions that fix what a CCS value means: ion, analyte chemistry, platform, reference gas and calibrant.

    Shared by measurements and prediction requests. A calibrant is required
    exactly when one enters the CCS value: TWIMS, TIMS, cyclic and single-field
    DTIMS. Any platform added later counts as calibrated until listed otherwise.
    """

    adduct: _Adduct
    charge: int = Field(strict=True, description="Signed charge state, e.g. 2 or -1.")
    polarity: Polarity
    reducing_end_label: ReducingEndLabel = Field(
        default=ReducingEndLabel.UNKNOWN,
        description="Label chemistry at the reducing end. Part of the calibration group; 'other' and 'unknown'"
        " block training.",
    )
    derivatisation: Derivatisation = Field(
        default=Derivatisation.UNKNOWN,
        description="Whole-molecule modification of the measured analyte, independent of the reducing-end label."
        " Part of the calibration group. 'underivatised' means nothing apart from the label was modified; a value"
        " that was not reported is 'unknown', which blocks training. The sialic-acid values block training and"
        " need a composition with NeuAc or NeuGc.",
    )
    ims_type: IMSType
    dtims_method: DTIMSMethod | None = Field(
        default=None, description="Required for DTIMS, and allowed only for DTIMS."
    )
    drift_gas: DriftGas = Field(
        description="The gas the reported CCS value refers to, not necessarily the gas in the cell. A TWIMS value"
        " measured in an N2 cell but calibrated against He reference values is a He value, and does not pool"
        " with N2 values. A measurement records the gas in the cell separately, as cell_gas."
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
        _, adduct_charge = _parse_adduct(self.adduct)
        if adduct_charge != self.charge:
            raise ValueError(f"adduct {self.adduct} carries charge {adduct_charge:+d}, but charge is {self.charge:+d}")
        return self

    @property
    def calibration_group(self) -> CalibrationGroup:
        return CalibrationGroup(
            self.ims_type,
            self.dtims_method,
            self.drift_gas,
            self.calibrant,
            self.adduct,
            self.reducing_end_label,
            self.derivatisation,
        )


class CCSMeasurement(MeasurementConditions):
    """One experimental CCS value, the glycan it was measured on, and where it came from."""

    glycan: GlycanStructure
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
    calibrant_reference: _Text | None = Field(
        default=None,
        description="What the calibrant was referenced against, as in 'helium CCS reference values': where the"
        " calibration curve came from. Provenance, kept out of the calibration group; drift_gas already says"
        " which gas the value refers to. Not allowed where no calibrant enters the value.",
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
    instrument: _Text | None = None
    source: _Text
    # Provenance per value. `source` is free text and cannot be verified to be a
    # citation, so these four carry the parts a later milestone has to be able to
    # act on: a DOI to group by study without scraping prose, a locator to find
    # the value again in the paper, a replicate count to tell a single shot from
    # a mean of five, and a date, without which no chronological split is
    # possible at all. Null where a source did not report them, never a filler.
    doi: Annotated[_Text, AfterValidator(_check_doi)] | None = Field(
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
        # Pydantic does not re-validate a glycan built with model_construct, so read it defensively.
        composition = getattr(self.glycan, "composition", None)
        if not isinstance(composition, Composition):
            raise ValueError(
                "glycan has no parsed composition, so it was built without validation; build it through"
                " GlycanStructure"
            )
        check_derivatisation_fits(composition, self.derivatisation)
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
        if self.calibrant_reference is not None and not self.ccs_is_calibrated:
            raise ValueError(
                "a primary CCS has no calibrant, so it cannot have a calibrant reference"
                " (a standard run only as a check is not a calibrant)"
            )
        return self

    def component_records(self) -> tuple[GlycanStructure, ...]:
        """Records this one is built from. The training gate checks them as well."""
        return (self.glycan,)

    def training_blockers(self) -> list[str]:
        """Reasons other than licence that this record may not enter training."""
        blockers = super().training_blockers()
        label = self.reducing_end_label
        if _is_one_of(label, _OPEN_LABELS):
            blockers.append(
                "reducing-end label is 'other', an open bucket that would pool unrelated labels in one"
                " calibration group; add the label to ReducingEndLabel instead"
            )
        elif not _is_one_of(label, _TRAINABLE_LABELS):
            blockers.append(f"reducing-end label is {str(label)!r}, so the measured analyte is not defined")
        derivatisation = self.derivatisation
        if _is_one_of(derivatisation, _SIALIC_ACID_DERIVATISATIONS):
            blockers.append(
                f"derivatisation is {str(derivatisation)!r}, an open bucket: it covers chemistries of different"
                " mass, and linkage-specific protocols have no correct value; add a specific value to"
                " Derivatisation instead"
            )
        elif not _is_one_of(derivatisation, _TRAINABLE_DERIVATISATIONS):
            blockers.append(f"derivatisation is {str(derivatisation)!r}, so the measured analyte is not defined")
        # Read defensively: a record built with model_construct may have no
        # drift_gas at all, and an AttributeError here would be caught upstream
        # and reported as "blockers could not be checked", hiding the specific
        # reasons this method exists to give. An absent gas is revalidation's
        # problem, not this check's.
        if _is_one_of(getattr(self, "drift_gas", None), _UNDEFINED_GASES):
            blockers.append(
                "drift_gas is 'UNSTATED', so the gas this CCS refers to is not defined: it cannot be pooled"
                " with helium or with nitrogen values, and resolving it means reading the calibration"
                " reference the source cites, not inferring it from the same group's other work"
            )
        kind = self.uncertainty_type
        if kind is not None and not _is_one_of(kind, _TRAINABLE_UNCERTAINTY_TYPES):
            blockers.append(
                f"uncertainty_type is {str(kind)!r}, so the reported spread cannot be interpreted; an interval built"
                " on it would be wrong by an unknown factor"
            )
        return blockers

    @property
    def validation_warnings(self) -> list[str]:
        """Problems that do not block the record but should be reviewed, including its glycan's."""
        return [f"glycan: {warning}" for warning in self.glycan.validation_warnings]
