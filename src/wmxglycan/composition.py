"""Glycan compositions: parsing, canonical form, monoisotopic mass, the Man3GlcNAc2
core check and mammalian N-glycan rules.

A composition string is a run of <Residue><count> tokens in any order, such as
"Hex5HexNAc4Fuc1" or "HexNAc4Fuc1Hex5". dHex is read as Fuc, the only
deoxyhexose in mammalian N-glycans. A string that does not parse completely is
rejected; nothing is ever partially parsed.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, fields
from enum import StrEnum
from typing import Iterable


class CompositionError(ValueError):
    """A composition that cannot be accepted as written."""


class Residue(StrEnum):
    """Supported residues. Definition order is the canonical output order."""

    HEX = "Hex"
    HEXNAC = "HexNAc"
    FUC = "Fuc"
    NEUAC = "NeuAc"
    NEUGC = "NeuGc"


# Monoisotopic residue masses in Da: each monosaccharide minus one H2O, as it
# sits inside a chain.
RESIDUE_MASS: dict[Residue, float] = {
    Residue.HEX: 162.052824,
    Residue.HEXNAC: 203.079373,
    Residue.FUC: 146.057909,
    Residue.NEUAC: 291.095417,
    Residue.NEUGC: 307.090331,
}

# The residues of a free glycan add up to one water short: the H and OH at the
# two ends of the chain.
WATER_MASS = 18.010565

# Accepted spellings, case-sensitive.
_SPELLINGS: dict[str, Residue] = {residue.value: residue for residue in Residue} | {"dHex": Residue.FUC}

# Composition field names are the lower-cased residue names.
_FIELD: dict[Residue, str] = {residue: residue.value.lower() for residue in Residue}

_TOKEN = re.compile(r"([A-Za-z]+)([0-9]+)")  # [0-9], not \d, which also matches non-ASCII digits
_LETTERS = re.compile(r"[A-Za-z]+")
_MAX_LENGTH = 256  # bound for strings arriving through the API; real compositions are far shorter

# The Man3GlcNAc2 core.
_CORE_HEX = 3
_CORE_HEXNAC = 2


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
    def parse(cls, text: str) -> Composition:
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

    @property
    def residue_mass(self) -> float:
        """Sum of monoisotopic residue masses, without the free-glycan water."""
        return sum(RESIDUE_MASS[residue] * self.count(residue) for residue in Residue)

    @property
    def monoisotopic_mass(self) -> float:
        """Monoisotopic neutral mass of the free glycan: residue masses plus one water."""
        return self.residue_mass + WATER_MASS

    def __str__(self) -> str:
        return self.canonical


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
    if len(s) > _MAX_LENGTH:
        raise CompositionError(f"composition longer than {_MAX_LENGTH} characters")

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


# Monosaccharide names, as written in IUPAC-condensed structures, that map onto
# the five residue classes. Everything else is left out on one of two grounds.
# Mass: a pentose, a uronic acid, Kdn, and sulfated, phosphorylated or
# methylated residues do not weigh what the class weighs, so calling one Hex
# would state a mass the molecule does not have. Identity: Galf, Rha, D-Fuc and
# the other furanose forms and epimers weigh exactly what their class weighs,
# but they are not the same sugar, and both ring form and epimer shift CCS,
# which is the value this platform exists to predict.
MONOSACCHARIDE_CLASSES: dict[str, Residue] = {
    "Hex": Residue.HEX,
    "Man": Residue.HEX,
    "Gal": Residue.HEX,
    "Glc": Residue.HEX,
    "HexNAc": Residue.HEXNAC,
    "GlcNAc": Residue.HEXNAC,
    "GalNAc": Residue.HEXNAC,
    "Fuc": Residue.FUC,
    "dHex": Residue.FUC,
    "Neu5Ac": Residue.NEUAC,
    "Neu5Gc": Residue.NEUGC,
}


class UnsupportedResidueError(CompositionError):
    """A structure holds residues that a composition of the five classes cannot represent."""

    def __init__(self, names: Iterable[str]) -> None:
        self.names = tuple(sorted(set(names)))
        super().__init__(
            f"cannot express {', '.join(self.names)} as Hex, HexNAc, Fuc, NeuAc or NeuGc: the class would"
            " state either the wrong mass or the wrong sugar"
        )


def composition_from_residue_names(names: Iterable[str]) -> Composition:
    """Composition of a structure given its monosaccharide names, e.g. ("Man", "Man", "GlcNAc").

    Raises UnsupportedResidueError naming every residue that has no class here,
    so a caller can report what it skipped rather than quietly rounding it off.
    """
    names = list(names)
    unsupported = [name for name in names if name not in MONOSACCHARIDE_CLASSES]
    if unsupported:
        raise UnsupportedResidueError(unsupported)
    counts: dict[str, int] = {}
    for name in names:
        field = _FIELD[MONOSACCHARIDE_CLASSES[name]]
        counts[field] = counts.get(field, 0) + 1
    return Composition(**counts)


def man3glcnac2_core_warnings(composition: Composition | str) -> list[str]:
    """Warnings for each residue class short of the full Man3GlcNAc2 core (3 Hex, 2 HexNAc).

    A missing core is a warning, not a rejection, because truncated N-glycans
    exist. Paucimannosidic N-glycans are Man1-3GlcNAc2, with or without core
    fucose, so those with one or two mannoses lack the full core. A glycan
    released by an endoglycosidase such as Endo H keeps only one core GlcNAc.

    An empty list means the counts cover the core. Counts alone cannot show
    that the residues are arranged as one.
    """
    comp = _as_composition(composition)
    found: list[str] = []
    if comp.hex < _CORE_HEX:
        note = ""
        if comp.hex >= 1 and comp.hexnac == _CORE_HEXNAC:
            note = "; with HexNAc2 this fits a truncated paucimannosidic N-glycan (Man1-2GlcNAc2)"
        found.append(f"Hex{comp.hex}: short of the 3 Hex in the full Man3GlcNAc2 core{note}")
    if comp.hexnac == 1:
        found.append(
            "HexNAc1: short of the 2 HexNAc in the full Man3GlcNAc2 core; fits a glycan released by an"
            " endoglycosidase such as Endo H, which cuts between the two core GlcNAc and leaves one on the protein"
        )
    elif comp.hexnac == 0:
        found.append("HexNAc0: none of the 2 HexNAc in the full Man3GlcNAc2 core")
    return found


def has_full_man3glcnac2_core(composition: Composition | str) -> bool:
    """True if the counts cover the full Man3GlcNAc2 core: at least 3 Hex and 2 HexNAc.

    Counts cannot show that the residues are arranged as that core.
    """
    return not man3glcnac2_core_warnings(composition)


def n_glycan_implausibility_reasons(composition: Composition | str) -> list[str]:
    """Readable reasons `composition` cannot be a mammalian N-glycan.

    Rules:
      1. At least one HexNAc. Every N-glycan keeps at least one core GlcNAc:
         both after release with the full chitobiose core, one after release
         by an endoglycosidase. This holds for any N-glycan, mammalian or not.
      2. With no HexNAc beyond two there is no antenna, so the glycan can carry
         a. no sialic acid, which needs a Gal or GalNAc on an antennary GlcNAc;
         b. at most one Fuc: in mammals the only fucose site left is the core
            alpha1,6 position.

    Rule 2 assumes mammalian biosynthesis. Rule b does not hold outside
    mammals: insect N-glycans, for example, can carry both an alpha1,3- and an
    alpha1,6-fucose on the Asn-linked core GlcNAc. Rule 2 also assumes that at
    HexNAc2 both HexNAc are the core GlcNAc, as in a glycan released with the
    full chitobiose core (by PNGase F, for example). A glycan released by an
    endoglycosidase keeps only one core GlcNAc, so its HexNAc2 can include an
    antennary GlcNAc, and rule 2 would misjudge it.

    A missing core is not a reason here; see man3glcnac2_core_warnings. An
    empty list means no rule is broken, which is plausibility, not proof.
    """
    comp = _as_composition(composition)
    reasons: list[str] = []
    if comp.hexnac == 0:
        reasons.append(
            "HexNAc0: every N-glycan keeps at least one core GlcNAc, so a composition with no HexNAc is not"
            " an N-glycan"
        )
    elif comp.hexnac <= _CORE_HEXNAC:
        sialic = [f"{r.value}{comp.count(r)}" for r in (Residue.NEUAC, Residue.NEUGC) if comp.count(r)]
        if sialic:
            reasons.append(
                f"{'+'.join(sialic)} with HexNAc{comp.hexnac}: sialic acid needs a Gal or GalNAc on an"
                " antenna, and with no HexNAc beyond the core there is no antenna"
            )
        if comp.fuc > 1:
            reasons.append(
                f"Fuc{comp.fuc} with HexNAc{comp.hexnac}: with no antenna, the core alpha1,6 position is"
                " the only mammalian fucose site, so at most Fuc1"
            )
    return reasons


def is_plausible_n_glycan(composition: Composition | str) -> bool:
    """True if `composition` breaks none of the rules in n_glycan_implausibility_reasons."""
    return not n_glycan_implausibility_reasons(composition)


def _as_composition(composition: Composition | str) -> Composition:
    return composition if isinstance(composition, Composition) else parse_composition(composition)


def _what_went_wrong(s: str, pos: int) -> str:
    letters = _LETTERS.match(s, pos)
    if letters is None:
        return f"expected a residue name, found {s[pos]!r}"
    end = letters.end()
    found = repr(s[end]) if end < len(s) else "end of string"
    return f"expected a count after {letters.group()!r}, found {found}"


def _hint(name: str) -> str:
    for spelling in _SPELLINGS:
        if spelling.casefold() == name.casefold():
            return f"did you mean {spelling!r}?"
    return f"supported: {', '.join(_SPELLINGS)}"
