"""Biosynthetic constraints, read from the installed glycowork package.

Each rule says that a product is forbidden or required in a given context, and
carries the enzyme, the rationale and the literature reference behind it. Those
three travel with the rule because a ranked candidate has to be able to show
why it was ranked that way; a rule that arrives without them, or with a
placeholder in their place, is an error rather than a row to load with the
reasoning dropped.

glycowork is MIT licensed (Copyright (c) 2021 Daniel Bojar), so the rules are
trainable and redistributable with attribution.
"""

from __future__ import annotations

import csv
from enum import StrEnum
from importlib import resources
from typing import Annotated, Iterable, Iterator, Self

from pydantic import AfterValidator, BaseModel, ConfigDict, StringConstraints, model_validator

from .models import _not_placeholder  # the guard the records use: "N/A" is not a citation

CONSTRAINTS_FILE = "network/biosynthetic_constraints.csv"

_Text = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1), AfterValidator(_not_placeholder)]


class ConstraintKind(StrEnum):
    FORBIDS = "forbids"
    REQUIRES = "requires"


class GlycanClass(StrEnum):
    N_LINKED = "N"
    O_LINKED = "O"


class BiosyntheticConstraint(BaseModel):
    """One curated rule about what an enzyme can and cannot build.

    The context is the precondition, in both directions: a "forbids" rule says
    the product cannot form while the context is present, and a "requires" rule
    says the product can only form once the context is already there. Neither
    says anything about what must happen wherever the context appears.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    product: _Text
    contexts: tuple[_Text, ...]
    kind: ConstraintKind
    glycan_class: GlycanClass | None  # None where the rule holds for any class
    enzyme: _Text
    rationale: _Text
    reference: _Text

    @model_validator(mode="after")
    def check_rule(self) -> Self:
        if not self.contexts:
            raise ValueError(f"the rule for {self.product!r} has no context to apply in")
        return self

    def applies_to(self, glycan_class: GlycanClass | str | None) -> bool:
        """True if this rule governs that class of glycan. A rule with no class governs all of them."""
        return self.glycan_class is None or glycan_class is None or self.glycan_class == glycan_class

    def explain(self) -> str:
        """The sentence to show beside a ranked candidate.

        The product stays the subject and the context stays the precondition,
        so neither kind of rule is read backwards.
        """
        contexts = " or ".join(self.contexts)
        credit = f"({self.enzyme}; {self.reference})"
        if self.kind is ConstraintKind.FORBIDS:
            return f"{self.product} cannot form where {contexts} is present: {self.rationale} {credit}"
        return f"{self.product} can only form where {contexts} is already present: {self.rationale} {credit}"


class ConstraintSet:
    """The curated rules, queryable by class, enzyme or product."""

    def __init__(self, constraints: Iterable[BiosyntheticConstraint]) -> None:
        self.constraints = tuple(constraints)
        if not self.constraints:
            raise ValueError("no biosynthetic constraints were loaded")

    @classmethod
    def from_glycowork(cls) -> ConstraintSet:
        """Read the rules that ship inside the installed glycowork package."""
        text = resources.files("glycowork").joinpath(CONSTRAINTS_FILE).read_text(encoding="utf-8-sig")
        return cls.from_csv(text)

    @classmethod
    def from_csv(cls, text: str) -> ConstraintSet:
        # The shipped file carries a byte order mark; a caller who read it as plain
        # UTF-8 hands us that mark on the front of the first column name.
        rows = list(csv.DictReader(text.lstrip("﻿").splitlines()))
        expected = {"product", "context", "kind", "glycan_class", "enzyme", "rationale", "reference"}
        missing = expected - set(rows[0] if rows else {})
        if missing:
            raise ValueError(f"the constraints file is missing the columns {sorted(missing)}")
        return cls(_constraint(row, line) for line, row in enumerate(rows, start=2))

    def for_glycan_class(self, glycan_class: GlycanClass | str) -> tuple[BiosyntheticConstraint, ...]:
        """Rules for that class, including the ones that hold for any class."""
        return tuple(rule for rule in self.constraints if rule.applies_to(glycan_class))

    def of_kind(self, kind: ConstraintKind | str) -> tuple[BiosyntheticConstraint, ...]:
        return tuple(rule for rule in self.constraints if rule.kind == kind)

    def for_enzyme(self, enzyme: str) -> tuple[BiosyntheticConstraint, ...]:
        return tuple(rule for rule in self.constraints if rule.enzyme == enzyme)

    def for_product(self, product: str) -> tuple[BiosyntheticConstraint, ...]:
        return tuple(rule for rule in self.constraints if rule.product == product)

    @property
    def enzymes(self) -> tuple[str, ...]:
        """The enzyme labels as curated; one of them names two enzymes together."""
        return tuple(sorted({rule.enzyme for rule in self.constraints}))

    def __len__(self) -> int:
        return len(self.constraints)

    def __iter__(self) -> Iterator[BiosyntheticConstraint]:
        return iter(self.constraints)


def _constraint(row: dict[str, str], line: int) -> BiosyntheticConstraint:
    glycan_class = (row.get("glycan_class") or "").strip()
    try:
        return BiosyntheticConstraint(
            product=row["product"],
            contexts=tuple(part.strip() for part in (row["context"] or "").split("|") if part.strip()),
            kind=(row["kind"] or "").strip(),
            glycan_class=glycan_class or None,
            enzyme=row["enzyme"],
            rationale=row["rationale"],
            reference=row["reference"],
        )
    except Exception as exc:
        # A rule that arrives without its enzyme, rationale or reference, or with a
        # placeholder standing in for one, cannot be shown beside a candidate. It stops
        # the load rather than loading half-explained.
        raise ValueError(f"constraint on line {line} cannot be loaded: {exc}") from exc
