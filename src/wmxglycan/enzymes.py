"""The glycoenzymes that build mammalian glycans, read from the installed glycowork package.

glycowork ships network/monolink_to_enzyme.csv: 355 rows, all human. A monolink
is a residue and the linkage it arrives on, such as "Gal(b1-4)", so the table
says which residue-and-linkage combinations an enzyme can produce at all. That
makes it a vocabulary rather than a set of rules: an enumerator can use it to
refuse a linkage no known enzyme makes. Where each linkage may sit is a
separate question, which the curated constraints in constraints.py answer part of.

Not every row states a monolink. Of the 355:

- 161 name a residue and a linkage, and those are the vocabulary;
- 161 record an enzyme with no monolink at all, such as the glycosidases and
  chaperones (AGA, AGL, ALG14);
- 33 record a modification rather than an attachment, spelling six distinct
  marks (3S, 4S, 6S, OS, GlcA2S, GlcNAc3S) written by the sulfotransferases.

The two skipped groups are counted and kept visible rather than dropped
quietly, because a shrinking vocabulary would silently narrow what any
enumerator built on it can propose. Sulfation is a real modification this
vocabulary cannot express, so it is reported rather than ignored.

glycowork is MIT licensed (Copyright (c) 2021 Daniel Bojar), so the table is
usable with attribution.
"""

from __future__ import annotations

import csv
from importlib import resources
from typing import Iterable, Iterator

ENZYME_FILE = "network/monolink_to_enzyme.csv"


class Monolink:
    """One residue and the linkage it arrives on, with the enzymes that make it."""

    __slots__ = ("residue", "linkage", "enzymes")

    def __init__(self, residue: str, linkage: str, enzymes: Iterable[str]) -> None:
        self.residue = residue
        self.linkage = linkage
        self.enzymes = tuple(sorted(set(enzymes)))

    @property
    def text(self) -> str:
        """As written in the table, e.g. "Gal(b1-4)"."""
        return f"{self.residue}({self.linkage})"

    def __repr__(self) -> str:
        return f"Monolink({self.text!r}, {len(self.enzymes)} enzymes)"


class EnzymeCatalogue:
    """Which residue-and-linkage combinations a human glycoenzyme is known to make."""

    def __init__(
        self,
        monolinks: Iterable[Monolink],
        *,
        rows_with_monolink: int = 0,
        rows_without_monolink: int = 0,
        modification_rows: Iterable[str] = (),
    ) -> None:
        self._by_text = {monolink.text: monolink for monolink in monolinks}
        self.rows_with_monolink = rows_with_monolink  # rows that named a residue and a linkage
        self.rows_without_monolink = rows_without_monolink
        written = list(modification_rows)
        self.rows_with_modification = len(written)
        self.modifications = tuple(sorted(set(written)))  # the distinct marks, e.g. ("3S", "6S")
        if not self._by_text:
            raise ValueError("no glycoenzymes were loaded")

    @property
    def rows_read(self) -> int:
        """Every row of the file, each counted in exactly one bucket."""
        return self.rows_with_monolink + self.rows_without_monolink + self.rows_with_modification

    @classmethod
    def from_glycowork(cls) -> EnzymeCatalogue:
        text = resources.files("glycowork").joinpath(ENZYME_FILE).read_text(encoding="utf-8-sig")
        return cls.from_csv(text)

    @classmethod
    def from_csv(cls, text: str) -> EnzymeCatalogue:
        reader = csv.DictReader(text.lstrip("﻿").splitlines())
        expected = {"glycoenzyme", "monolink", "species"}
        # The header decides which columns exist; a file with a good header and no
        # rows is empty, not malformed, and must not be blamed for missing columns.
        missing = expected - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"the glycoenzyme file is missing the columns {sorted(missing)}")
        rows = list(reader)
        makers: dict[tuple[str, str], set[str]] = {}
        with_monolink = 0
        without_monolink = 0
        modification_rows: list[str] = []
        for line, row in enumerate(rows, start=2):
            written = (row["monolink"] or "").strip()
            enzyme = (row["glycoenzyme"] or "").strip()
            if not enzyme:
                raise ValueError(f"glycoenzyme row on line {line} names no enzyme")
            if not written:
                without_monolink += 1  # an enzyme with no attachment recorded, such as a glycosidase
                continue
            if not (written.endswith(")") and "(" in written):
                modification_rows.append(written)  # a mark such as "6S", not a residue on a linkage
                continue
            residue, _, linkage = written[:-1].partition("(")
            if not residue or not linkage:
                raise ValueError(f"cannot read the monolink {written!r} on line {line} as residue(linkage)")
            makers.setdefault((residue, linkage), set()).add(enzyme)
            with_monolink += 1
        return cls(
            (Monolink(residue, linkage, enzymes) for (residue, linkage), enzymes in makers.items()),
            rows_with_monolink=with_monolink,
            rows_without_monolink=without_monolink,
            modification_rows=modification_rows,
        )

    def makes(self, residue: str, linkage: str) -> bool:
        """True if some human enzyme is known to attach that residue on that linkage."""
        return f"{residue}({linkage})" in self._by_text

    def enzymes_for(self, residue: str, linkage: str) -> tuple[str, ...]:
        """The enzymes known to make that monolink, or an empty tuple if none is."""
        monolink = self._by_text.get(f"{residue}({linkage})")
        return monolink.enzymes if monolink else ()

    def monolinks_for(self, residue: str) -> tuple[Monolink, ...]:
        """Every linkage that residue is known to arrive on."""
        return tuple(sorted((m for m in self._by_text.values() if m.residue == residue), key=lambda m: m.linkage))

    @property
    def residues(self) -> tuple[str, ...]:
        return tuple(sorted({monolink.residue for monolink in self._by_text.values()}))

    @property
    def enzyme_count(self) -> int:
        return len({enzyme for monolink in self._by_text.values() for enzyme in monolink.enzymes})

    def summary(self) -> str:
        return (
            f"{self.rows_with_monolink} of {self.rows_read} rows name a residue and a linkage,"
            f" giving {len(self._by_text)} monolinks made by {self.enzyme_count} enzymes;"
            f" {self.rows_without_monolink} rows record an enzyme with no monolink"
            f" and {self.rows_with_modification} record a modification rather than an attachment"
            f" ({', '.join(self.modifications)})"
        )

    def __len__(self) -> int:
        return len(self._by_text)

    def __iter__(self) -> Iterator[Monolink]:
        return iter(sorted(self._by_text.values(), key=lambda m: m.text))
