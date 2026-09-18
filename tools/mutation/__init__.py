"""Mutation testing: break the code on purpose, and require a test to notice.

A green suite says the tests ran. It does not say they would have caught
anything. This tool answers the other question. For each behaviour worth
protecting, change the source so that behaviour is wrong, run the suite, and
require it to fail. A mutation that survives is a behaviour with no test behind
it, however many tests appear to cover the file.

The catalogue is hand-written, not generated. Every entry was authored against
one specific guard, almost always after a real defect or a review finding, so
its label says what would go WRONG rather than which operator was flipped. That
is what makes a survivor readable: it names a promise the code no longer keeps.

WHY AN UNAPPLIED MUTATION IS A FAILURE, NOT A NOTE
--------------------------------------------------
An earlier version of this harness printed "SKIPPED (pattern not found)" and
still exited zero. That is the same failure this repository keeps meeting in
other clothes: a report of success that covers less than it claims, alongside a
CSV row swallowed by a stray quote and a source count that silenced its own
warning. Twice, an ordinary edit moved the line a mutation anchored on, the
mutation quietly stopped being applied, and the sweep went on reporting that
everything was killed.

So there are three ways to not-apply, and all three are failures of the run:

- STALE      the anchor matches nothing. The guard moved, or went away.
- AMBIGUOUS  the anchor matches more than once, so which site it would hit is
             a coin toss. The old harness took the first silently.
- INERT      it applied and changed nothing, which tests nothing.

A mutation that goes stale should be RE-ANCHORED, not retired. Retiring one is
a deliberate act: delete it, and say in a comment what now covers its intent.

EXPECTATIONS, SO THAT SURVIVAL CAN ALSO BE A FAILURE
-----------------------------------------------------
One guard in this repository is known to be unreachable and is kept anyway, so
its mutation survives every run by design. If survival were always a failure,
that one would fail forever and the signal would be trained away. So each
mutation declares what it expects, and an exemption must give a reason.

That buys a property worth having in both directions: an exemption that stops
being true also fails. If a test starts killing a documented survivor, the
documentation has rotted, and the run says so rather than quietly agreeing.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from enum import StrEnum
from typing import Iterable


class Expect(StrEnum):
    """What a mutation is expected to do to the suite."""

    KILLED = "killed"
    SURVIVES = "survives"


class Anchor(StrEnum):
    """Whether a mutation can be applied to the source at all.

    Everything but OK is a failure of the tool run. See the module docstring:
    an unapplied mutation reporting success is the whole reason this exists.
    """

    OK = "ok"
    STALE = "stale"
    AMBIGUOUS = "ambiguous"
    INERT = "inert"


class Verdict(StrEnum):
    """What the suite did once the mutation was applied."""

    KILLED = "killed"
    SURVIVED = "survived"
    # The suite could not run at all: a collection error, a usage error, an
    # internal pytest failure. Deliberately NOT a kill. A mutation that stops
    # the suite from starting has not been noticed by anything - it has merely
    # broken the runner - and counting it as killed would be a success the run
    # did not earn, which is the failure this tool exists to refuse.
    UNRUNNABLE = "unrunnable"


@dataclass(frozen=True)
class Mutation:
    """One deliberate defect, and what is expected to happen when it is made.

    `find` must identify exactly one site in `file`. It is matched literally,
    not as a pattern, because a regex that drifts is harder to notice than a
    string that stops matching.
    """

    label: str
    file: str
    find: str
    replace: str
    expect: Expect = Expect.KILLED
    reason: str = ""  # required when the mutation is expected to survive
    # Further edits, applied in order after the first. Needed for a mutation
    # that MOVES something - a guard placed after the check it was meant to
    # precede, say - which is two coordinated changes and cannot be expressed
    # as one. Each is anchored against the text as it stands when its turn
    # comes, so an earlier edit may legitimately create the site a later one
    # names.
    also: tuple[tuple[str, str], ...] = ()

    def __post_init__(self) -> None:
        if not self.label.strip():
            raise ValueError("a mutation needs a label saying what would go wrong")
        if not self.file.strip():
            raise ValueError(f"{self.label}: names no file to mutate")
        for find, _replace in self.edits:
            if not find:
                raise ValueError(f"{self.label}: an empty anchor would match anywhere")
        if all(find == replace for find, replace in self.edits):
            raise ValueError(f"{self.label}: every edit is identical, so it would change nothing")
        object.__setattr__(self, "expect", Expect(self.expect))
        if self.expect is Expect.SURVIVES and not self.reason.strip():
            raise ValueError(
                f"{self.label}: a mutation expected to survive must say why no test can reach it."
                " Without a reason, an exemption cannot be told apart from a gap in the tests."
            )

    @property
    def edits(self) -> tuple[tuple[str, str], ...]:
        """Every edit this mutation makes, in the order it makes them."""
        return ((self.find, self.replace), *self.also)

    def anchor_in(self, source: str) -> Anchor:
        """Whether every edit identifies exactly one site, and the whole changes something."""
        text = source
        for find, replace in self.edits:
            found = text.count(find)
            if found == 0:
                return Anchor.STALE
            if found > 1:
                return Anchor.AMBIGUOUS
            text = text.replace(find, replace, 1)
        return Anchor.INERT if text == source else Anchor.OK

    def apply_to(self, source: str) -> str:
        text = source
        for find, replace in self.edits:
            text = text.replace(find, replace, 1)
        return text

    @property
    def expected_verdict(self) -> Verdict:
        return Verdict.KILLED if self.expect is Expect.KILLED else Verdict.SURVIVED


def repeated_labels(mutations: Iterable[Mutation]) -> list[str]:
    """Labels used more than once. A duplicate makes a filtered run ambiguous."""
    counted = Counter(mutation.label for mutation in mutations)
    return sorted(label for label, count in counted.items() if count > 1)
