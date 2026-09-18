"""What was measured, and the key that decides when two measurements are the same ion.

THE MATCHED ION IS THE UNIT THIS PLATFORM IS BUILT ON
-----------------------------------------------------
Every number the platform produces compares one ion measured on one platform
with the SAME ion measured on another. "The same ion" is not "the same
compound": it is the analyte's identity, the adduct, the signed charge, the gas
the value refers to, and the structural state, all five together. Two rows
naming one compound and differing in any of the five are different ions, and
averaging them would pool quantities that are not the same quantity.

So the key deliberately does NOT contain the platform. Excluding it is the
whole point: a DTIMS value and a TWIMS value of one ion must land on the same
key, or there is nothing to compare and no bias to measure.

WHY A NAME IS NEVER AN IDENTITY
-------------------------------
Every variant below identifies itself by a declared identifier - an InChIKey, a
sequence, an accession, an INN, a canonical composition - and never by the free
text a paper happens to print in its compound column. `display_name` is carried
on every analyte and is in no key: it is what a report prints, not what a match
is made on. Two rows both saying "trastuzumab" match only when their adduct,
charge, gas and structural state agree as well.

WHERE THIS MODULE SITS
----------------------
Below models.py, for the same reason reuse.py sits below licensing.py: models.py
imports the analyte and the key from here, so nothing here may import models.py.
The record primitives both need - placeholder rejection, the frozen-record base,
the revalidation check - therefore live here too rather than in the module that
reads more naturally, because a cycle is a worse problem than a surprising home.

PORTED, AND WHAT CHANGED
------------------------
The adduct grammar, the placeholder rejection, the revalidation check and the
composition representation are ported from the glycan platform unchanged in
behaviour. Two things moved on purpose:

- `reducing_end_label` and `derivatisation` were fields of the measurement
  there, where every record was a glycan. Here they are fields of the GLYCAN
  ANALYTE, because a steroid has no reducing end and a field that can only ever
  be null on five of six analyte kinds is not a field, it is a leak.
- the composition representation is ported WITHOUT the N-glycan plausibility
  rules and the Man3GlcNAc2 core check that surrounded it in the glycan
  platform. Those are biosynthetic rules, which do not come across. Canonical
  spelling does, because without it "Hex5HexNAc4" and "HexNAc4Hex5" are two
  matched ions instead of one, which is a missed pair reported as a clean run.
"""

from __future__ import annotations

import re
import warnings
from collections import Counter
from dataclasses import dataclass, fields
from enum import StrEnum
from typing import Annotated, Literal, NamedTuple, Self, Union

from pydantic import (
    AfterValidator,
    BaseModel,
    ConfigDict,
    Field,
    PlainSerializer,
    PlainValidator,
    StringConstraints,
    ValidationError,
    WithJsonSchema,
    model_validator,
)

from .reuse import DEFAULT_REUSE_STATUS, ReuseStatus


class UnverifiedFormatWarning(UserWarning):
    """A value does not match a format rule that has not itself been verified yet.

    A UserWarning, so Python shows it by default.
    """


# --- text that is actually a value ----------------------------------------------------

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


# --- the frozen record base ------------------------------------------------------------


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


# --- the ion: adduct and charge --------------------------------------------------------

# Bracket notation with the charge after the bracket: [M+H]+, [M+2H]2+, [M-H]-, [2M+Na]+.
_ADDUCT = re.compile(r"\[(?P<mult>[1-9][0-9]*)?M(?P<parts>[^\[\]\s]*)\](?P<n>[1-9][0-9]*)?(?P<sign>[+-])")

# The species standing for a charge carrier the source did not name. Written
# "[M+24?]24+": twenty-four charges, carrier unstated.
#
# It exists because of native mass spectrometry, and because the alternative is
# to invent chemistry. Papers reporting intact proteins and antibodies routinely
# give a charge state and nothing else - "the 24+ ion" - without saying whether
# those charges are protons, sodium or ammonium. Without this form such a record
# cannot be built at all, because the adduct is mandatory and its charge must
# match the charge field; and writing "[M+24H]24+" instead would state that the
# carriers are protons, which the paper does not say. That is fabrication, and
# it is the one thing this platform refuses above all others.
#
# So it follows the DriftGas.UNSTATED precedent exactly: a positive record of
# what the source does not state, keying APART from any named carrier so the two
# can never be pooled, and blocking training because the ion is not defined.
UNSTATED_CARRIER = "?"

# One component of the ion: a sign, an optional count and a species such as H,
# Na, NH4, HCOO or H2O - or "?" for a carrier the source did not name.
_ADDUCT_PART = re.compile(r"([+-])([1-9][0-9]*)?([A-Z][A-Za-z0-9]*|\?)")
_MAX_ADDUCT_LENGTH = 100  # bound for strings arriving through the API; real adducts are far shorter


def parse_adduct(adduct: str) -> tuple[str, int]:
    """Canonical spelling and signed charge of an adduct.

    Components are sorted, additions before losses and then by species, so
    "[M+Na+H]2+" and "[M+H+Na]2+" give the same key. Counts, multipliers and
    charges of 1 are dropped: "[M+1H]1+" -> ("[M+H]+", 1). A species named
    twice is rejected rather than merged.

    This is what makes two spellings of one ion one matched ion. Without it the
    pairing step would report a clean run over half the pairs it should have
    found, which is this platform's worst failure shape.
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

    parse_adduct keeps only the canonical spelling and the signed charge, which
    is all a record needs in order to be well formed. Anything reading the ion
    itself needs the parts: a sodiated ion is not a protonated one of the same
    charge, and the two differ in mass and in collision cross section. Parts come
    back in canonical order, additions before losses and then by species, so two
    spellings of one ion give one sequence.

    Raises ValueError for an adduct that is not in bracket notation.
    """
    canonical, _ = parse_adduct(adduct)  # validates first; a malformed adduct raises here
    match = _ADDUCT.fullmatch(canonical)
    assert match is not None  # parse_adduct built this string, so it matches by construction
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


def adduct_carrier_is_unstated(adduct: str) -> bool:
    """True if any component of the ion is a charge carrier the source did not name.

    Such an adduct is well formed and keys apart from every named carrier, but
    the ion it describes is not defined: the mass depends on whether those
    charges are protons or sodium, and so does the cross section. It blocks
    training, and resolving it means reading the paper's methods, never assuming
    protons because protons are usual.
    """
    return any(species == UNSTATED_CARRIER for _sign, _count, species in adduct_components(adduct))


Adduct = Annotated[
    str,
    StringConstraints(strip_whitespace=True, max_length=_MAX_ADDUCT_LENGTH),
    AfterValidator(lambda value: parse_adduct(value)[0]),
]


class Polarity(StrEnum):
    POSITIVE = "positive"
    NEGATIVE = "negative"


class DriftGas(StrEnum):
    """A gas used for ion mobility.

    drift_gas is the gas a CCS value REFERS TO; a measurement's cell_gas is the
    gas that was in the cell. They differ when, for example, a TWIMS value from
    an N2 cell is calibrated against He reference values. drift_gas is in the
    matched-ion key and cell_gas is not, because what a value refers to is what
    decides whether it may be compared with another.
    """

    N2 = "N2"
    HE = "He"
    # The source did not say which gas its values refer to. A positive fact, not
    # a filler: it records what the paper does not state, so the gap is visible
    # rather than guessed at, and it is deliberately not spelt "unknown", which
    # is a placeholder this package refuses everywhere else. It blocks training,
    # because a CCS whose gas is undefined cannot be pooled with helium or with
    # nitrogen values, and it keeps a matched-ion key of its own for the same
    # reason. Resolving it means reading the calibration reference a paper
    # cites, never inferring from what the same group used elsewhere.
    UNSTATED = "UNSTATED"


# A drift gas that does not define the quantity the value refers to.
UNDEFINED_GASES = frozenset({DriftGas.UNSTATED})


# --- glycan composition: the canonicalisation half of the ported representation ---------


class CompositionError(ValueError):
    """A composition that cannot be accepted as written."""


class Residue(StrEnum):
    """Supported residues. Definition order is the canonical output order."""

    HEX = "Hex"
    HEXNAC = "HexNAc"
    FUC = "Fuc"
    NEUAC = "NeuAc"
    NEUGC = "NeuGc"


# Accepted spellings, case-sensitive.
_SPELLINGS: dict[str, Residue] = {residue.value: residue for residue in Residue} | {"dHex": Residue.FUC}
# Composition field names are the lower-cased residue names.
_FIELD: dict[Residue, str] = {residue: residue.value.lower() for residue in Residue}

_TOKEN = re.compile(r"([A-Za-z]+)([0-9]+)")  # [0-9], not \d, which also matches non-ASCII digits
_LETTERS = re.compile(r"[A-Za-z]+")
_MAX_COMPOSITION_LENGTH = 256  # bound for strings arriving through the API; real compositions are far shorter


@dataclass(frozen=True)
class Composition:
    """Monosaccharide counts. Immutable; two spellings of one composition compare equal."""

    hex: int = 0
    hexnac: int = 0
    fuc: int = 0
    neuac: int = 0
    neugc: int = 0

    def __post_init__(self) -> None:
        for field in fields(self):
            value = getattr(self, field.name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise CompositionError(f"{field.name} count must be an int, not {value!r}")
            if value < 0:
                raise CompositionError(f"{field.name} count cannot be negative, got {value}")
        if not any(getattr(self, field.name) for field in fields(self)):
            raise CompositionError("a composition needs at least one residue")

    @classmethod
    def parse(cls, text: str) -> "Composition":
        return parse_composition(text)

    def count(self, residue: Residue) -> int:
        return getattr(self, _FIELD[residue])

    def counts(self) -> dict[str, int]:
        """Non-zero counts keyed by residue name, in canonical order."""
        return {residue.value: n for residue in Residue if (n := self.count(residue))}

    @property
    def canonical(self) -> str:
        """Canonical spelling: fixed residue order, Fuc rather than dHex, zero counts omitted."""
        return "".join(f"{name}{n}" for name, n in self.counts().items())

    def __str__(self) -> str:
        return self.canonical


def _what_went_wrong(s: str, pos: int) -> str:
    rest = s[pos:]
    letters = _LETTERS.match(rest)
    if letters is not None:
        return f"{letters.group()!r} is not followed by a count"
    return f"{rest[:1]!r} is not a residue name"


def _hint(name: str) -> str:
    known = ", ".join(sorted(_SPELLINGS))
    for spelling in _SPELLINGS:
        if spelling.casefold() == name.casefold():
            return f"residue names are case-sensitive; did you mean {spelling!r}? Known residues: {known}"
    return f"known residues: {known}"


def parse_composition(text: str) -> Composition:
    """Parse a composition such as "Hex5HexNAc4Fuc1", residues in any order.

    Raises CompositionError unless the whole string is <Residue><count> tokens:
    unknown or mis-cased residue names, missing counts, stray characters,
    repeated residues (dHex and Fuc are the same residue), counts with leading
    zeros and compositions with no residues are all rejected.
    """
    if not isinstance(text, str):
        raise TypeError(f"composition must be a str, not {type(text).__name__}")
    s = text.strip()
    if not s:
        raise CompositionError("empty composition")
    if len(s) > _MAX_COMPOSITION_LENGTH:
        raise CompositionError(f"composition longer than {_MAX_COMPOSITION_LENGTH} characters")

    counts: dict[Residue, int] = {}
    pos = 0
    while pos < len(s):
        token = _TOKEN.match(s, pos)
        if token is None:
            raise CompositionError(f"cannot parse {text!r} at position {pos}: {_what_went_wrong(s, pos)}")
        name, digits = token.groups()
        residue = _SPELLINGS.get(name)
        if residue is None:
            raise CompositionError(f"unknown residue {name!r} in {text!r}; {_hint(name)}")
        if residue in counts:
            note = " (dHex is read as Fuc)" if residue is Residue.FUC else ""
            raise CompositionError(f"{residue.value} given more than once in {text!r}{note}")
        if len(digits) > 1 and digits[0] == "0":
            raise CompositionError(f"count {digits!r} for {name} has a leading zero in {text!r}")
        counts[residue] = int(digits)
        pos = token.end()

    try:
        return Composition(**{_FIELD[residue]: n for residue, n in counts.items()})
    except CompositionError as exc:
        raise CompositionError(f"{exc}: {text!r}") from None


def canonical_composition(text: str) -> str:
    """Canonical spelling of `text`, e.g. "HexNAc4dHex1Hex5" -> "Hex5HexNAc4Fuc1"."""
    return parse_composition(text).canonical


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


# --- structural state ------------------------------------------------------------------


class FoldingState(StrEnum):
    """Whether the ion was sprayed from folded or unfolded conditions.

    Part of IDENTITY for a protein, a subunit, an antibody or an ADC, not
    metadata about one: a native 24+ antibody ion and a denatured 40+ ion of the
    same antibody are two different ions with two different cross sections, and
    a platform that pooled them would report a bias that is really the
    difference between a folded and an unfolded molecule.
    """

    NATIVE = "native"
    DENATURED = "denatured"
    # The source did not say. Spelt like DriftGas.UNSTATED and for the same
    # reason: it is a positive record of what the paper does not state, not a
    # filler, and it is deliberately not "unknown", which this package refuses
    # as a placeholder everywhere else. It blocks training and keys apart from
    # both real states, because an ion of unrecorded conformation cannot be
    # matched with one that states its conformation.
    UNSTATED = "UNSTATED"


UNDEFINED_FOLDING_STATES = frozenset({FoldingState.UNSTATED})


class ReducingEndLabel(StrEnum):
    """The label chemistry at the reducing end of a measured glycan.

    A labelled glycan is a different analyte from the native one, with its own
    mass and CCS, so the label is part of the structural state. "other" and
    "unknown" can be stored but block training: an open bucket would pool
    unrelated labels under one key. A real label missing from this list is
    added to it, not recorded as "other".
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
    """Whole-molecule modification of a measured glycan, independent of the reducing-end label.

    Part of the structural state. The value describes the measured analyte, not
    the sample protocol: "underivatised" means nothing was modified apart from
    the reducing-end label, and a derivatisation that was not reported is
    "unknown", which can be stored but blocks training.

    The two sialic-acid values are open buckets, like an "other" label: each
    covers chemistries of different mass, and a linkage-specific protocol that
    treats alpha2,3- and alpha2,6-linked sialic acids differently has no
    correct value. They can be stored but block training, and they need a
    composition with NeuAc or NeuGc. A glycan with no sialic acid from such a
    sample is "underivatised" only if nothing else on it reacted: the carboxyl
    group of a 2-AA label can react with the same reagents, so a neutral 2-AA
    glycan from such a sample is "unknown". A specific value is added here when
    a real dataset needs one.

    Spelt "underivatised" and never "none", because "none" is a placeholder
    spelling this package refuses, and an enum value spelt like a filler is a
    value a loader can write by accident.
    """

    UNDERIVATISED = "underivatised"
    PERMETHYLATION = "permethylation"
    SIALIC_ACID_AMIDATION = "sialic_acid_amidation"
    SIALIC_ACID_ESTERIFICATION = "sialic_acid_esterification"
    UNKNOWN = "unknown"


# What defines the measured analyte well enough to train on. Membership, not
# identity, so a raw string that skipped validation is judged by its value.
TRAINABLE_LABELS = frozenset(ReducingEndLabel) - {ReducingEndLabel.OTHER, ReducingEndLabel.UNKNOWN}
OPEN_LABELS = frozenset({ReducingEndLabel.OTHER})
SIALIC_ACID_DERIVATISATIONS = frozenset(
    {Derivatisation.SIALIC_ACID_AMIDATION, Derivatisation.SIALIC_ACID_ESTERIFICATION}
)
TRAINABLE_DERIVATISATIONS = frozenset(Derivatisation) - {Derivatisation.UNKNOWN} - SIALIC_ACID_DERIVATISATIONS


def is_one_of(value: object, allowed: frozenset) -> bool:
    try:
        return value in allowed
    except (TypeError, ValueError):  # unhashable, or a comparison with no truth value, after skipped validation
        return False


class AnalyteKind(StrEnum):
    """What sort of thing was measured. The tag of the analyte union."""

    SMALL_MOLECULE = "small_molecule"
    PEPTIDE = "peptide"
    GLYCAN = "glycan"
    PROTEIN = "protein"
    INTACT_ANTIBODY = "intact_antibody"
    ADC = "adc"


# --- identifier formats ----------------------------------------------------------------

# Specified by IUPAC: 14 characters of skeleton, 10 of stereo and version, one
# of protonation. A specified shape, so a mismatch is rejected rather than
# warned about - unlike the GlyTouCan form below, which was inferred from
# accessions seen in use.
_INCHIKEY = re.compile(r"[A-Z]{14}-[A-Z]{10}-[A-Z]")

# The form G + five digits + two capital letters comes from accessions seen in
# use, not from GlyTouCan's documentation. Until it is verified, a mismatch
# warns instead of rejecting.
_GLYTOUCAN_AC = re.compile(r"G[0-9]{5}[A-Z]{2}")

# An IUPAC-condensed structure states its linkages in parentheses, as in
# "Man(a1-3)Man". A string with no linkage in it says nothing about how any
# residues are joined, so it cannot be the evidence behind a claim that the
# linkages are resolved, however structure-shaped it looks.
_IUPAC_LINKAGE = re.compile(r"\([ab?]?[0-9?]+-[0-9?]+(?:/[0-9?]+)*\)")

# The twenty standard amino acids, plus selenocysteine and pyrrolysine, which
# are real residues with real masses. The ambiguity codes are deliberately
# absent: B is Asp-or-Asn, Z is Glu-or-Gln, J is Leu-or-Ile and X is anything
# at all. A sequence holding one does not define a molecule, so it cannot
# define an ion, and a key built on it would match things that are not the same.
_PEPTIDE_SEQUENCE = re.compile(r"[ACDEFGHIKLMNPQRSTVWYUO]+")


def _check_inchikey(value: str) -> str:
    if not _INCHIKEY.fullmatch(value):
        raise ValueError(
            f"{value!r} is not an InChIKey: the form is 14 capitals, a hyphen, 10 capitals, a hyphen and one"
            " capital, as in RYYVLZVUVIJVGH-UHFFFAOYSA-N"
        )
    return value


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


def _check_wurcs(value: str) -> str:
    if not value.startswith("WURCS=") or any(ch.isspace() for ch in value):
        raise ValueError("WURCS must start with 'WURCS=' and contain no whitespace")
    return value


def _check_peptide_sequence(value: str) -> str:
    sequence = value.strip().upper()
    if not _PEPTIDE_SEQUENCE.fullmatch(sequence):
        offending = sorted({ch for ch in sequence if not _PEPTIDE_SEQUENCE.fullmatch(ch)})
        raise ValueError(
            f"peptide sequence holds {offending or ['no residues']} which are not single-letter amino acids."
            " B, Z, J and X are ambiguity codes and are refused: a sequence holding one does not define a"
            " molecule, so it cannot define an ion"
        )
    return sequence


_Sequence = Annotated[str, StringConstraints(strip_whitespace=True), AfterValidator(_check_peptide_sequence)]


def _sorted_unique(values: tuple[str, ...]) -> tuple[str, ...]:
    """Modifications in a canonical order, so two spellings of one set give one key.

    Sorted rather than left as written, because "Phospho@S5, Acetyl@K2" and
    "Acetyl@K2, Phospho@S5" are the same molecule. A modification named twice is
    refused rather than merged: it is either a transcription fault or a claim
    about two sites, and guessing which would be inventing data.
    """
    repeated = sorted(name for name, times in Counter(values).items() if times > 1)
    if repeated:
        raise ValueError(f"modification(s) {repeated} named more than once; write each site once")
    return tuple(sorted(values))


_Modifications = Annotated[tuple[_Text, ...], AfterValidator(_sorted_unique)]


# --- the analyte union -----------------------------------------------------------------


class _Analyte(_Record):
    """What was measured. Every variant carries where its ASSIGNMENT came from, and on what terms.

    `source` and `reuse_status` are the analyte's own, not the measurement's:
    where a structure or an accession was assigned is not where the CCS value
    was read, and a value from an open paper does not clear an identity taken
    from a restricted one. The licence gate walks into these as component
    records, which is why every variant must carry both - a component with no
    reuse_status refuses its whole parent, by design.

    `display_name` is in NO key. It is what a report prints - "LNH", "Man5",
    "trastuzumab" - and a paper's compound column is free text that two groups
    spell two ways. Identity is always a declared identifier.
    """

    display_name: _Text | None = Field(
        default=None,
        description="What to print for this analyte. Never part of any key: a compound name is free text,"
        " and two records matching on a name alone would not be the same ion.",
    )
    source: _Text
    reuse_status: ReuseStatus = DEFAULT_REUSE_STATUS

    @property
    def kind(self) -> AnalyteKind:
        return AnalyteKind(self.kind_tag)  # type: ignore[attr-defined]

    def identity_key(self) -> tuple:
        """The identity half of the matched-ion key: the finest identifier this analyte states.

        Deterministic and single-valued, so two records either key alike or do
        not. Where one molecule is named through two different identifiers - a
        WURCS here, an IUPAC string there - the two do not key alike, and
        merging them is the job of matched-ion construction in M2, which has
        `identity_atoms` to do it with. Merging them here would mean choosing
        silently between two identifiers that may not describe one molecule.
        """
        raise NotImplementedError

    def identity_atoms(self) -> frozenset[str]:
        """Every identifier this analyte states, as namespaced atoms.

        The seam for M2. Two analytes sharing any atom are the same molecule and
        may be merged into one matched ion; this is the input that pairing needs,
        and the reason a record keeps every identifier it was given rather than
        only the finest one.
        """
        raise NotImplementedError

    def structural_state(self) -> tuple:
        """The structural-state half of the matched-ion key. Empty where the kind has none."""
        return ()

    @property
    def validation_warnings(self) -> list[str]:
        """Problems that do not block the record but should be reviewed."""
        return []


class SmallMoleculeAnalyte(_Analyte):
    """A small molecule, identified by its InChIKey.

    The InChIKey is REQUIRED and is the whole identity. It is a hash of the
    structure, so two groups that measured the same molecule produce the same
    key without ever agreeing on what to call it, which is exactly what a
    compound name cannot do. SMILES is optional and is not the key: one molecule
    has many valid SMILES strings and one InChIKey.
    """

    kind_tag: Literal[AnalyteKind.SMALL_MOLECULE] = AnalyteKind.SMALL_MOLECULE
    inchikey: Annotated[_Text, AfterValidator(_check_inchikey)]
    smiles: _Text | None = Field(
        default=None,
        description="Optional, and never the key: one molecule has many valid SMILES and one InChIKey.",
    )

    def identity_key(self) -> tuple:
        return (AnalyteKind.SMALL_MOLECULE.value, "inchikey", self.inchikey)

    def identity_atoms(self) -> frozenset[str]:
        return frozenset({f"inchikey:{self.inchikey}"})


class PeptideAnalyte(_Analyte):
    """A peptide, identified by its sequence and its modifications.

    Modifications are part of the identity, not metadata: a phosphopeptide is
    eighty daltons heavier than its unmodified form and has its own cross
    section. They are held in canonical order, so two spellings of one set give
    one key.
    """

    kind_tag: Literal[AnalyteKind.PEPTIDE] = AnalyteKind.PEPTIDE
    sequence: _Sequence
    modifications: _Modifications = Field(
        default=(),
        description="Modifications as the source names them, e.g. 'Phospho@S5'. Held in a canonical order."
        " Part of the identity: a modified peptide is a different molecule.",
    )

    def identity_key(self) -> tuple:
        return (AnalyteKind.PEPTIDE.value, "sequence", self.sequence, self.modifications)

    def identity_atoms(self) -> frozenset[str]:
        return frozenset({f"peptide:{self.sequence}|{'+'.join(self.modifications)}"})


class GlycanAnalyte(_Analyte):
    """A glycan, as far as it is actually known, and what was done to it before measuring.

    Linkage and anomericity are unresolved unless the record says otherwise, and
    only a record with a structure identifier can say otherwise: a WURCS, a
    GlyTouCan accession in the expected form, or the IUPAC-condensed string the
    structure was read from.

    IDENTITY IS THE FINEST IDENTIFIER STATED, NOT THE COMPOSITION. The glycan
    platform split its TRAINING folds on composition, deliberately, because a
    prediction request carries no structure and a structure-derived split key is
    finer than the model's own identity. That rule is about splitting and it
    does not come across to matching. LNH and LNnH share the composition
    Hex4HexNAc2 and differ by 11.4 per cent as [M-H]-; keying them alike would
    pair two different molecules and report the difference between them as
    inter-platform bias. Grouping for a fit and matching for a comparison are
    two different keys, and one name must not serve both.

    The reducing-end label and the derivatisation live here, on the analyte,
    rather than on the measurement where the glycan platform kept them. A
    labelled glycan is a different molecule from the native one, and a steroid
    has no reducing end at all.
    """

    kind_tag: Literal[AnalyteKind.GLYCAN] = AnalyteKind.GLYCAN
    composition: CompositionField | None = Field(
        default=None,
        description="Monosaccharide counts, where the source gives them. OPTIONAL, unlike in the glycan"
        " platform, where every record was a glycan and a composition was always derived from the structure."
        " Here it is derived from nothing: deriving a composition from an IUPAC string needs a table mapping"
        " monosaccharide names onto residue classes, which is glycan chemistry and does not come across."
        " A source that gives a structure and no composition therefore records null, and is identified by its"
        " structure, which is finer anyway. A glycan must still state SOMETHING - see states_an_identifier.",
    )
    wurcs: Annotated[str, StringConstraints(strip_whitespace=True), AfterValidator(_check_wurcs)] | None = None
    glytoucan_ac: Annotated[_Text, AfterValidator(_warn_on_glytoucan_format)] | None = None
    iupac_condensed: Annotated[_Text, StringConstraints(pattern=r"^\S+$")] | None = Field(
        default=None,
        description="The structure as an IUPAC-condensed string, e.g. Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc."
        " It states the linkages, so it counts as a structure identifier.",
    )
    has_unresolved_linkage: bool = Field(default=True, strict=True)
    has_unresolved_anomericity: bool = Field(default=True, strict=True)
    reducing_end_label: ReducingEndLabel = Field(
        default=ReducingEndLabel.UNKNOWN,
        description="Label chemistry at the reducing end. Part of the structural state; 'other' and 'unknown'"
        " block training.",
    )
    derivatisation: Derivatisation = Field(
        default=Derivatisation.UNKNOWN,
        description="Whole-molecule modification, independent of the reducing-end label. Part of the structural"
        " state. 'underivatised' means nothing apart from the label was modified; a value that was not reported"
        " is 'unknown', which blocks training. The sialic-acid values block training and need a composition"
        " with NeuAc or NeuGc.",
    )

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

    @model_validator(mode="after")
    def states_an_identifier(self) -> Self:
        """A glycan identified by nothing at all has no key, so it cannot be a matched ion."""
        if self.composition is None and self.wurcs is None and self.glytoucan_ac is None and self.iupac_condensed is None:
            raise ValueError(
                "a glycan needs a composition, a WURCS, a GlyTouCan accession or an IUPAC-condensed structure:"
                " with none of them it has no identity, and a measurement of it could never be matched to"
                " another measurement of the same molecule"
            )
        return self

    @model_validator(mode="after")
    def derivatisation_fits_the_composition(self) -> Self:
        """A sialic-acid derivatisation needs a sialic acid to have acted on.

        Where no composition is recorded the claim cannot be checked, so it is
        refused rather than waved through: a derivatisation that names a residue
        must be stated against a record that says whether the residue is there.
        """
        if self.derivatisation in SIALIC_ACID_DERIVATISATIONS and self.composition is None:
            raise ValueError(
                f"derivatisation {str(self.derivatisation)!r} names a sialic acid, but this record gives no"
                " composition, so whether there is one to act on cannot be checked. Record the composition,"
                " or record what was done to this analyte in terms that do not name a residue"
            )
        if (
            self.derivatisation in SIALIC_ACID_DERIVATISATIONS
            and self.composition is not None
            and self.composition.neuac + self.composition.neugc == 0
        ):
            raise ValueError(
                f"derivatisation {str(self.derivatisation)!r} needs a sialic acid to act on, and"
                f" {self.composition.canonical} has no NeuAc or NeuGc. Record what was done to this analyte:"
                " 'underivatised' if nothing on it was modified, 'unknown' if something else may have reacted,"
                " such as the carboxyl of a 2-AA label"
            )
        return self

    def identity_key(self) -> tuple:
        """The finest identifier stated, in a fixed preference order."""
        if self.wurcs is not None:
            return (AnalyteKind.GLYCAN.value, "wurcs", self.wurcs)
        if self.glytoucan_ac is not None and _glytoucan_format_problem(self.glytoucan_ac) is None:
            return (AnalyteKind.GLYCAN.value, "glytoucan", self.glytoucan_ac)
        if self.iupac_condensed is not None and _IUPAC_LINKAGE.search(self.iupac_condensed) is not None:
            return (AnalyteKind.GLYCAN.value, "iupac", self.iupac_condensed)
        return (AnalyteKind.GLYCAN.value, "composition", self.composition.canonical)

    def identity_atoms(self) -> frozenset[str]:
        atoms = set()
        if self.composition is not None:
            atoms.add(f"composition:{self.composition.canonical}")
        if self.wurcs is not None:
            atoms.add(f"wurcs:{self.wurcs}")
        if self.glytoucan_ac is not None:
            atoms.add(f"glytoucan:{self.glytoucan_ac}")
        if self.iupac_condensed is not None:
            atoms.add(f"iupac:{self.iupac_condensed}")
        return frozenset(atoms)

    def structural_state(self) -> tuple:
        return (str(self.reducing_end_label), str(self.derivatisation))

    def training_blockers(self) -> list[str]:
        blockers = super().training_blockers()
        label = self.reducing_end_label
        if is_one_of(label, OPEN_LABELS):
            blockers.append(
                "reducing-end label is 'other', an open bucket that would pool unrelated labels under one"
                " matched-ion key; add the label to ReducingEndLabel instead"
            )
        elif not is_one_of(label, TRAINABLE_LABELS):
            blockers.append(f"reducing-end label is {str(label)!r}, so the measured analyte is not defined")
        derivatisation = self.derivatisation
        if is_one_of(derivatisation, SIALIC_ACID_DERIVATISATIONS):
            blockers.append(
                f"derivatisation is {str(derivatisation)!r}, an open bucket: it covers chemistries of different"
                " mass, and linkage-specific protocols have no correct value; add a specific value to"
                " Derivatisation instead"
            )
        elif not is_one_of(derivatisation, TRAINABLE_DERIVATISATIONS):
            blockers.append(f"derivatisation is {str(derivatisation)!r}, so the measured analyte is not defined")
        return blockers

    @property
    def validation_warnings(self) -> list[str]:
        if self.glytoucan_ac is None:
            return []
        problem = _glytoucan_format_problem(self.glytoucan_ac)
        return [] if problem is None else [problem]


class _FoldedAnalyte(_Analyte):
    """Shared by the three kinds whose conformation is part of what was measured.

    A protein, an antibody and an ADC are sprayed either from a solution that
    keeps them folded or from one that does not, and the two give different
    cross sections for the same molecule. So folding_state is in the key.
    """

    folding_state: FoldingState = Field(
        default=FoldingState.UNSTATED,
        description="Native or denatured, as the source states it. Part of the IDENTITY of the ion, not"
        " metadata about it. 'UNSTATED' records that the source does not say, keys apart from both, and"
        " blocks training.",
    )

    def training_blockers(self) -> list[str]:
        blockers = super().training_blockers()
        if is_one_of(getattr(self, "folding_state", None), UNDEFINED_FOLDING_STATES):
            blockers.append(
                "folding_state is 'UNSTATED', so it is not known whether this ion was folded or unfolded."
                " A native and a denatured ion of one molecule have different cross sections, so a value"
                " whose conformation is unrecorded cannot be compared with one whose is"
            )
        return blockers


class ProteinAnalyte(_FoldedAnalyte):
    """A protein or a subunit, identified by an accession or a sequence.

    At least one of the two is required. A name is not enough: "the light chain"
    identifies nothing across two laboratories, and a subunit in particular is
    whatever the digest produced.
    """

    kind_tag: Literal[AnalyteKind.PROTEIN] = AnalyteKind.PROTEIN
    accession: _Text | None = Field(
        default=None,
        description="A database accession, e.g. a UniProt one. No format is enforced: the accession may come"
        " from any of several databases, and refusing an unrecognised shape would refuse real records.",
    )
    sequence: _Sequence | None = None
    subunit: _Text | None = Field(
        default=None,
        description="Which part of a larger molecule this is, where it is a part, e.g. 'light chain'."
        " Part of the identity, because a light chain is not a heavy chain.",
    )

    @model_validator(mode="after")
    def states_an_identifier(self) -> Self:
        if self.accession is None and self.sequence is None:
            raise ValueError(
                "a protein or subunit needs an accession or a sequence: a display name identifies nothing"
                " across two laboratories, and matching on one would pair ions that are not the same ion"
            )
        return self

    def identity_key(self) -> tuple:
        which = ("accession", self.accession) if self.accession is not None else ("sequence", self.sequence)
        return (AnalyteKind.PROTEIN.value, *which, self.subunit)

    def identity_atoms(self) -> frozenset[str]:
        atoms = set()
        if self.accession is not None:
            atoms.add(f"accession:{self.accession}|{self.subunit or ''}")
        if self.sequence is not None:
            atoms.add(f"sequence:{self.sequence}")
        return frozenset(atoms)

    def structural_state(self) -> tuple:
        return (str(self.folding_state),)


class AntibodyIdentity(_Record):
    """Which antibody, by a declared identifier rather than by what a paper calls it.

    At least one of an INN, an accession or a sequence. The INN - "trastuzumab",
    "rituximab" - is a controlled name assigned by the WHO, not free text, and
    for an intact monoclonal it is very often the only identifier a paper gives.
    It is held in its own field, and lower-cased, so that it can be an identity
    without `display_name` ever becoming one.
    """

    inn: _Text | None = Field(
        default=None,
        description="International nonproprietary name, e.g. 'trastuzumab'. A controlled name, held"
        " lower-cased so two spellings give one key.",
    )
    accession: _Text | None = None
    sequence: _Sequence | None = None

    @model_validator(mode="after")
    def states_an_identifier(self) -> Self:
        if self.inn is None and self.accession is None and self.sequence is None:
            raise ValueError(
                "an antibody needs an INN, an accession or a sequence; a display name is not an identity"
            )
        return self

    @model_validator(mode="after")
    def normalise_inn(self) -> Self:
        if self.inn is not None and self.inn != self.inn.casefold():
            object.__setattr__(self, "inn", self.inn.casefold())
        return self

    def key(self) -> tuple:
        if self.inn is not None:
            return ("inn", self.inn)
        if self.accession is not None:
            return ("accession", self.accession)
        return ("sequence", self.sequence)

    def atoms(self) -> frozenset[str]:
        atoms = set()
        if self.inn is not None:
            atoms.add(f"inn:{self.inn}")
        if self.accession is not None:
            atoms.add(f"accession:{self.accession}")
        if self.sequence is not None:
            atoms.add(f"sequence:{self.sequence}")
        return frozenset(atoms)


class IntactAntibodyAnalyte(_FoldedAnalyte):
    """A whole antibody, measured intact.

    Glycoform is part of the identity because it is part of the mass: a G0F/G0F
    antibody and a G2F/G2F one differ by four hexoses. It is optional and not
    required to be stated, because public intact-antibody CCS measurements
    almost never resolve it and requiring it would refuse the entire biopharma
    layer this platform exists to build. An unstated glycoform keys apart from
    a stated one, so the two are never pooled.

    CIU state is the collision-induced-unfolding step the value was taken at.
    An antibody stepped through an activation ramp gives several cross sections
    for one charge state, and they are not replicates of one number.
    """

    kind_tag: Literal[AnalyteKind.INTACT_ANTIBODY] = AnalyteKind.INTACT_ANTIBODY
    antibody: AntibodyIdentity
    glycoform: _Text | None = Field(
        default=None,
        description="The glycoform as the source states it, e.g. 'G0F/G0F'. Part of the identity: it is part"
        " of the mass. Null where the source does not resolve it, which keys apart from any stated glycoform.",
    )
    ciu_state: _Text | None = Field(
        default=None,
        description="Which collision-induced-unfolding state this value was taken at, as the source labels it."
        " Null for a value not taken from a CIU ramp. Part of the structural state: an unfolded intermediate"
        " is not the same ion population as the compact form.",
    )

    def identity_key(self) -> tuple:
        return (AnalyteKind.INTACT_ANTIBODY.value, *self.antibody.key(), self.glycoform)

    def identity_atoms(self) -> frozenset[str]:
        return frozenset(f"{atom}|glycoform={self.glycoform or ''}" for atom in self.antibody.atoms())

    def structural_state(self) -> tuple:
        return (str(self.folding_state), self.ciu_state)

    def component_records(self) -> tuple[_Record, ...]:
        return (self.antibody,)


class ADCAnalyte(_FoldedAnalyte):
    """An antibody-drug conjugate: an antibody, a payload class, and how much of it is attached.

    DAR - the drug-to-antibody ratio - is part of the identity, because a DAR 2
    species and a DAR 4 species are different molecules with different masses
    and different cross sections. A record states a resolved DAR, or a named
    conjugation state where the source describes the species some other way,
    and at least one of the two is required: an ADC whose loading is unrecorded
    is not a defined molecule.
    """

    kind_tag: Literal[AnalyteKind.ADC] = AnalyteKind.ADC
    antibody: AntibodyIdentity
    linker_payload_class: _Text = Field(
        description="The class of linker and payload, e.g. 'vc-MMAE'. Required: two conjugates of one"
        " antibody at one DAR are different molecules if the payload differs.",
    )
    dar: int | None = Field(
        default=None,
        strict=True,
        ge=0,
        description="Drug-to-antibody ratio for a RESOLVED species, counting from zero. Not an average over"
        " a distribution: an average describes a mixture, and a mixture is not an ion.",
    )
    conjugation_state: _Text | None = Field(
        default=None,
        description="How the source names the species where it does not give a resolved DAR. Used only when"
        " dar is null.",
    )
    glycoform: _Text | None = None
    ciu_state: _Text | None = None

    @model_validator(mode="after")
    def states_its_loading(self) -> Self:
        if self.dar is None and self.conjugation_state is None:
            raise ValueError(
                "an ADC needs a resolved DAR or a named conjugation state: a conjugate whose loading is"
                " unrecorded is not a defined molecule, and its mass and cross section are not defined either"
            )
        return self

    def identity_key(self) -> tuple:
        loading = ("dar", self.dar) if self.dar is not None else ("conjugation", self.conjugation_state)
        return (
            AnalyteKind.ADC.value,
            *self.antibody.key(),
            self.linker_payload_class,
            *loading,
            self.glycoform,
        )

    def identity_atoms(self) -> frozenset[str]:
        loading = f"dar={self.dar}" if self.dar is not None else f"conjugation={self.conjugation_state}"
        return frozenset(
            f"{atom}|{self.linker_payload_class}|{loading}|glycoform={self.glycoform or ''}"
            for atom in self.antibody.atoms()
        )

    def structural_state(self) -> tuple:
        return (str(self.folding_state), self.ciu_state)

    def component_records(self) -> tuple[_Record, ...]:
        return (self.antibody,)


# The tagged union. Discriminated EXPLICITLY on kind_tag, never left to
# pydantic's smart union. _revalidation_problem checks a record by dumping it,
# re-validating the dump and comparing the two dumps; under a smart union a
# dumped protein can re-validate as a different arm while the dumps still
# compare equal, and the guard would pass on a record whose analyte kind had
# silently changed underneath it. A discriminator makes that unrepresentable.
Analyte = Annotated[
    Union[
        SmallMoleculeAnalyte,
        PeptideAnalyte,
        GlycanAnalyte,
        ProteinAnalyte,
        IntactAntibodyAnalyte,
        ADCAnalyte,
    ],
    Field(discriminator="kind_tag"),
]

ANALYTE_TYPES: dict[AnalyteKind, type[_Analyte]] = {
    AnalyteKind.SMALL_MOLECULE: SmallMoleculeAnalyte,
    AnalyteKind.PEPTIDE: PeptideAnalyte,
    AnalyteKind.GLYCAN: GlycanAnalyte,
    AnalyteKind.PROTEIN: ProteinAnalyte,
    AnalyteKind.INTACT_ANTIBODY: IntactAntibodyAnalyte,
    AnalyteKind.ADC: ADCAnalyte,
}


# --- the matched-ion key ---------------------------------------------------------------


class MatchedIonKey(NamedTuple):
    """The unit of comparison: when two measurements are of the same ion.

    Five parts, and every one of them is necessary:

    - `analyte`   which molecule, by a declared identifier and never by a name;
    - `adduct`    canonical, so two spellings of one ion give one key;
    - `charge`    signed, and part of identity for a protein or an antibody,
                  where a 24+ and a 40+ ion of one molecule are two ions;
    - `drift_gas` the gas the VALUE REFERS TO, not the gas in the cell. A TWIMS
                  value measured in nitrogen but calibrated against helium
                  references is a helium value and pools with helium values;
    - `state`     the structural state: the label and derivatisation of a
                  glycan, the folding state of a protein, the folding and CIU
                  state of an antibody or an ADC.

    WHAT IS DELIBERATELY ABSENT, AND WHY IT MATTERS MORE THAN WHAT IS PRESENT:

    - the PLATFORM. A DTIMS value and a TWIMS value of one ion must land on the
      same key or there is nothing to compare and no bias to measure. This is
      the whole reason the key exists, and it is what makes this key different
      from the calibration group, which is conditions-only and DOES carry the
      platform and the calibrant.
    - the CONFORMER INDEX. One structure can give more than one arrival-time
      peak under a single set of conditions. Two CCS values under one matched-
      ion key differing only in conformer index are a legitimate pair, not a
      duplicate and not a collision: they are not averaged and neither is
      dropped. Putting the index in the key would make them two ions and hide
      that they are one.
    - the CELL GAS, the instrument, the calibrant, the laboratory and the DOI.
      All provenance. If any of them entered the key, no two laboratories could
      ever produce a matched pair, and the platform would report zero pairs
      over a corpus full of them.
    """

    analyte: tuple
    adduct: str
    charge: int
    drift_gas: DriftGas
    state: tuple

    def __str__(self) -> str:
        analyte = ":".join("" if part is None else str(part) for part in self.analyte)
        state = ",".join("" if part is None else str(part) for part in self.state)
        return f"{analyte} | {self.adduct} | {self.charge:+d} | {self.drift_gas} | {state}"
