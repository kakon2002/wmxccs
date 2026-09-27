"""Candidate N-glycan structures for a composition, under mammalian biosynthetic constraints.

Candidates are built forwards from the Man3GlcNAc2 core by attaching residues
only where a human glycoenzyme is known to attach them. Two data sources decide
what is possible, and both ship with glycowork:

- the glycoenzyme table (enzymes.py) says which residue-and-linkage
  combinations an enzyme can make at all, which is what keeps the space from
  including linkages no enzyme produces;
- the curated constraints (constraints.py) say which combinations are ruled out
  or required by the order enzymes act in, and each candidate keeps the rules
  that bear on it so a ranking layer can show its reasoning.

What this enumerator does not do, stated plainly:

- It builds on the complete branched Man3GlcNAc2 core and cannot represent a
  truncated one, so TRUNCATED paucimannosidic species (Man1-2GlcNAc2) and
  degradation species are outside its reach by construction rather than by rule.
  Man3GlcNAc2 and Man3GlcNAc2Fuc1 DO enumerate - one candidate each - and both
  are paucimannosidic under the usual Man1-3GlcNAc2 definition, so the limit is
  the truncated core and not the word. This is a stated scope limit for an
  antibody platform, where paucimannose is rare on the Fc: 421 of the 4,001
  fully resolved reference structures lack a complete core, so no coverage
  figure above 89.5 per cent is reachable without changing the design.

- Of the fifteen curated rules, ELEVEN reach N-glycan enumeration: nine MGAT
  rules governing branching order and bisecting interference, FUT8 on core
  fucosylation, and one class-agnostic blood group rule. The other four are
  O-glycan core rules and are filtered out before enumeration. (This paragraph
  said "the curated rules govern branching order and bisecting interference
  only" until 27 September 2026 - the sixth copy of that sentence, and the one a
  search for "15 curated rules" could not find, because it states the claim with
  no number in it.) Nothing in the eleven constrains galactosylation type, fucose
  position, chain extension or LacdiNAc, so the antenna space is unconstrained. The 167
  candidates for Hex5HexNAc4Fuc1 are what an unbounded elaboration vocabulary
  produces; the figure is not tuned to any expectation. Excluding type-1 LacNAc
  and poly-LacNAc would bring it to 39, but no curated rule licenses that
  exclusion, so it is not applied. Any future narrowing of this number must come
  from a rule with a citation, never from choosing a vocabulary to hit a target.
- Every placement it makes is at a stated position, so every candidate comes out
  fully resolved. Ambiguity lives in the candidate *set*, not inside a candidate:
  several candidates for one composition are the statement that the structure is
  not determined. Nothing here completes an unresolved input, because the input
  is a composition, which carries no linkages to leave open.

Candidates are generated, not observed. They carry reuse status
internal_proprietary: they are our own output, derived from MIT-licensed rules.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable, Iterator, Sequence

from pydantic import BaseModel, ConfigDict

from .composition import Composition, Residue, n_glycan_implausibility_reasons
from .constraints import ConstraintKind, ConstraintSet, GlycanClass
from .enzymes import EnzymeCatalogue
from .glycan_graph import GlycanGraph
from .licensing import ReuseStatus
from .models import GlycanStructure

# THE TREE BUDGET. A STATED SCOPE LIMIT, not a tuning knob, and named here because it was an
# undisclosed default argument until 28 September 2026 - at which point a composition that exhausted
# it was told "every one of the N arrangements broke a biosynthetic rule", which is a claim about
# glycobiology and was false.
#
# RAISED FROM 5,000, and the reason is that 5,000 was suppressing answers rather than protecting
# anything. Measured completion, trees built and wall clock:
#
#     Hex5HexNAc4Fuc1     1,012   0.4s  complete
#     Hex6HexNAc5Fuc1    12,844   5.1s  complete  ->  1,729 candidates in 4 bands, RANKABLE
#     Hex6HexNAc5NeuAc2  24,052   8.7s  complete  ->  4,620 candidates in 4 bands, RANKABLE
#     Hex6HexNAc5NeuAc3  54,204  26.0s  complete  -> 11,672 candidates in 5 bands, RANKABLE
#
# At 5,000 the last three were refusals. They are ordinary human N-glycans and the platform can
# answer for them, so the budget was the defect.
#
# AND NO BUDGET FIXES THE TAIL. Hex7HexNAc6NeuAc3 and Hex8HexNAc7NeuAc4 exhaust 60,000 as well
# (28.9s and 37.9s, zero candidates), because the arrangement count grows combinatorially with the
# antennae. So this is a scope limit that is reported, not a number that can be raised until the
# problem goes away - which is why every response carries it and why exhaustion now says so.
#
# The cost is stated rather than hidden: a composition that cannot complete burns about thirty
# seconds before refusing, where at 5,000 it burned two and said something untrue.
DEFAULT_TREE_BUDGET = 60_000

CANDIDATE_SOURCE = "enumerated by wmxglycan from a composition under mammalian biosynthetic constraints"

# Which composition class a residue belongs to.
RESIDUE_CLASS = {
    "Man": Residue.HEX,
    "Gal": Residue.HEX,
    "GlcNAc": Residue.HEXNAC,
    "GalNAc": Residue.HEXNAC,
    "Fuc": Residue.FUC,
    "Neu5Ac": Residue.NEUAC,
    "Neu5Gc": Residue.NEUGC,
}


@dataclass(frozen=True)
class Site:
    """One place a residue can be attached, and what attaching it means."""

    linkage: str  # as written, e.g. "b1-2"
    residue: str  # what is attached
    kind: str  # the kind the attached residue becomes, which decides its own sites
    role: str  # what this placement is, in words, for the report and the rationale

    @property
    def monolink(self) -> str:
        return f"{self.residue}({self.linkage})"


# What can be attached to what, in a mammalian N-glycan. Each entry is checked
# against the glycoenzyme table before it is used, so a site whose monolink no
# human enzyme makes is dropped rather than quietly enumerated.
SITES: dict[str, tuple[Site, ...]] = {
    "reducing_glcnac": (Site("a1-6", "Fuc", "fuc", "core fucose"),),
    "chitobiose_glcnac": (),
    "beta_man": (Site("b1-4", "GlcNAc", "bisecting_glcnac", "bisecting GlcNAc"),),
    "alpha3_man": (
        Site("b1-2", "GlcNAc", "antenna_glcnac", "antenna on the alpha1-3 arm"),
        Site("b1-4", "GlcNAc", "antenna_glcnac", "second antenna on the alpha1-3 arm"),
        Site("a1-2", "Man", "arm_man", "alpha1-2 mannose"),
    ),
    "alpha6_man": (
        Site("b1-2", "GlcNAc", "antenna_glcnac", "antenna on the alpha1-6 arm"),
        Site("b1-6", "GlcNAc", "antenna_glcnac", "second antenna on the alpha1-6 arm"),
        Site("a1-3", "Man", "arm_man", "alpha1-3 mannose"),
        Site("a1-6", "Man", "arm_man", "alpha1-6 mannose"),
    ),
    "arm_man": (Site("a1-2", "Man", "arm_man", "alpha1-2 mannose cap"),),
    "bisecting_glcnac": (),  # the bisecting GlcNAc is not elongated
    "antenna_glcnac": (
        Site("b1-4", "Gal", "gal", "type-2 LacNAc"),
        Site("b1-3", "Gal", "gal", "type-1 LacNAc"),
        Site("b1-4", "GalNAc", "galnac", "LacdiNAc"),
        Site("a1-3", "Fuc", "fuc", "Lewis fucose"),
    ),
    "gal": (
        Site("a2-3", "Neu5Ac", "neuac", "alpha2-3 sialic acid"),
        Site("a2-6", "Neu5Ac", "neuac", "alpha2-6 sialic acid"),
        Site("a1-2", "Fuc", "fuc", "H antigen fucose"),
        Site("a1-3", "Gal", "agal", "alpha-Gal"),
        Site("b1-3", "GlcNAc", "antenna_glcnac", "poly-LacNAc extension"),
    ),
    "galnac": (Site("a2-6", "Neu5Ac", "neuac", "alpha2-6 sialic acid on LacdiNAc"),),
    "agal": (),
    "neuac": (Site("a2-8", "Neu5Ac", "neuac", "alpha2-8 sialic acid"),),
    "fuc": (),
}

_FRAGMENT = re.compile(r"([A-Za-z0-9]+)(?:\(([^)]+)\))?")


class EnumerationError(ValueError):
    """A fragment or rule this module cannot read."""


class Node:
    """A residue in a candidate under construction."""

    __slots__ = ("name", "kind", "children")

    def __init__(self, name: str, kind: str) -> None:
        self.name = name
        self.kind = kind
        self.children: list[tuple[str, Node]] = []  # (linkage, child)

    def render(self) -> str:
        """The candidate as an IUPAC-condensed string."""
        if not self.children:
            return self.name
        written = sorted(
            (f"{child.render()}({linkage})" for linkage, child in self.children),
            key=lambda text: (-len(text), text),
        )
        return written[0] + "".join(f"[{branch}]" for branch in written[1:]) + self.name

    def walk(self) -> Iterator["Node"]:
        yield self
        for _, child in self.children:
            yield from child.walk()


def _core() -> Node:
    """Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc, the core every N-glycan is built on."""
    reducing = Node("GlcNAc", "reducing_glcnac")
    chitobiose = Node("GlcNAc", "chitobiose_glcnac")
    beta = Node("Man", "beta_man")
    reducing.children.append(("b1-4", chitobiose))
    chitobiose.children.append(("b1-4", beta))
    beta.children.append(("a1-3", Node("Man", "alpha3_man")))
    beta.children.append(("a1-6", Node("Man", "alpha6_man")))
    return reducing


class CandidateReason(BaseModel):
    """Why a candidate is plausible: a placement an enzyme makes, or a rule that bears on it.

    For an enzyme reason, `enzymes` are those known to make that *monolink* --
    the residue on that linkage -- somewhere in human glycobiology. The shipped
    table records no acceptor for all but six of its rows, so this is not
    attribution of the individual bond: a monolink's list can include enzymes
    that make it in a different context, and even hydrolases that only cleave
    it. The list is never truncated, so nothing is dropped without being seen.
    """

    model_config = ConfigDict(frozen=True)

    kind: str  # "enzyme", "rule" or "ordering"
    text: str
    enzymes: tuple[str, ...] = ()
    reference: str | None = None


class GlycanCandidate(BaseModel):
    """One enumerated structure, with the reasoning that supports it."""

    model_config = ConfigDict(frozen=True)

    structure: GlycanStructure
    canonical_key: str
    rationale: tuple[CandidateReason, ...]

    @property
    def iupac_condensed(self) -> str:
        return self.structure.iupac_condensed or ""


@dataclass(frozen=True)
class EnumerationResult:
    """The candidates for one composition, and what shaped the set."""

    composition: Composition
    candidates: tuple[GlycanCandidate, ...]
    built: int = 0  # trees generated before the rules and deduplication
    rejected: int = 0  # trees a rule threw out; each may have broken several rules
    rejected_by_rule: dict[str, int] = field(default_factory=dict)
    duplicates: int = 0
    capped: bool = False
    refusal: str | None = None  # why no candidate could be built at all

    def __len__(self) -> int:
        return len(self.candidates)

    def summary(self) -> str:
        if self.refusal:
            return f"{self.composition}: no candidates ({self.refusal})"
        lines = [
            f"{self.composition}: {len(self.candidates)} candidates"
            f" (built {self.built}, {self.rejected} rejected,"
            f" {self.duplicates} duplicates removed{', capped' if self.capped else ''})"
        ]
        for rule, count in sorted(self.rejected_by_rule.items(), key=lambda item: -item[1]):
            lines.append(f"    rejected by {rule}: {count}")
        return "\n".join(lines)


def is_order_constraint(rule) -> bool:
    """True if a rule describes the order enzymes act in rather than what cannot coexist.

    The curated rationales draw the line themselves. Four of the nine FORBIDS
    rules read "blocks subsequent ...", and FUT8's finishes "the reverse order is
    permitted and common". All four have the bisecting GlcNAc as their context:
    bisecting stops an enzyme acting *later*, it does not unmake what already
    happened. This enumerator generates finished structures, not assembly
    histories, so a structure the permitted order can reach is a real structure.
    Such a rule therefore travels with the candidate as an ordering note instead
    of excluding it.

    The word is the curated evidence, so it is what is tested against; a reworded
    table fails the test rather than quietly changing what is enumerated.
    """
    return rule.kind is ConstraintKind.FORBIDS and "subsequent" in rule.rationale.lower()


class Enumerator:
    """Builds candidate N-glycans for a composition from the two glycowork sources."""

    def __init__(
        self,
        catalogue: EnzymeCatalogue | None = None,
        constraints: ConstraintSet | None = None,
        sites: dict[str, tuple[Site, ...]] | None = None,
    ) -> None:
        self.catalogue = catalogue if catalogue is not None else EnzymeCatalogue.from_glycowork()
        self.constraints = constraints if constraints is not None else ConstraintSet.from_glycowork()
        self.dropped_sites: tuple[Site, ...] = ()
        self.sites = self._usable(sites if sites is not None else SITES)
        self._rules = tuple(
            (rule, _chain(rule.product), tuple(_chain(context) for context in rule.contexts))
            for rule in self.constraints.for_glycan_class(GlycanClass.N_LINKED)
        )
        self.order_constraints = tuple(rule for rule, _, _ in self._rules if is_order_constraint(rule))

    def _usable(self, sites: dict[str, tuple[Site, ...]]) -> dict[str, tuple[Site, ...]]:
        """Keep only the placements a human glycoenzyme is known to make."""
        kept: dict[str, tuple[Site, ...]] = {}
        dropped: list[Site] = []
        for kind, options in sites.items():
            usable = []
            for site in options:
                if self.catalogue.makes(site.residue, site.linkage):
                    usable.append(site)
                else:
                    dropped.append(site)
            kept[kind] = tuple(usable)
        self.dropped_sites = tuple(dropped)
        return kept

    @property
    def placeable_residues(self) -> frozenset[Residue]:
        """The residues this enumerator can actually attach, read off the site vocabulary.

        THE APPLICABILITY DOMAIN, and the single source for it. `enumerate` refuses a composition
        needing anything outside this set, and `api.py` serves it as `residues_supported`, so the
        advertised domain and the refusal are the same fact. Deriving it also means that adding a
        site for a residue makes it supported in both places at once rather than in one of them.

        NeuGc is in the composition parser's alphabet and is NOT in here: no site in the kept
        vocabulary attaches it, because no enzyme in the curated table licenses one. Until
        27 September 2026 `residues_supported` was filled from the parser's alphabet, so the
        response advertised NeuGc while `enumerate` refused every composition containing it.
        """
        return frozenset(
            residue
            for residue in Residue
            if any(
                RESIDUE_CLASS.get(site.residue) is residue
                for options in self.sites.values()
                for site in options
            )
        )

    def enumerate(
        self, composition: Composition | str, limit: int = DEFAULT_TREE_BUDGET
    ) -> EnumerationResult:
        """Every candidate for `composition`, deduplicated on canonical structure identity."""
        composition = composition if isinstance(composition, Composition) else Composition.parse(composition)
        if composition.hex < 3 or composition.hexnac < 2:
            return EnumerationResult(
                composition, (), refusal="fewer residues than the Man3GlcNAc2 core needs"
            )
        budget = {
            Residue.HEX: composition.hex - 3,
            Residue.HEXNAC: composition.hexnac - 2,
            Residue.FUC: composition.fuc,
            Residue.NEUAC: composition.neuac,
            Residue.NEUGC: composition.neugc,
        }
        # DERIVED FROM `placeable_residues`, which is also what the API serves as the
        # applicability domain. Until 27 September 2026 this expression was inline here and the
        # served field was filled from the composition parser's alphabet instead, so the response
        # advertised NeuGc while this check refused every composition containing it.
        placeable = self.placeable_residues
        unplaceable = [
            residue.value for residue, count in budget.items() if count and residue not in placeable
        ]
        if unplaceable:
            return EnumerationResult(
                composition, (), refusal=f"no human enzyme in the table attaches {', '.join(unplaceable)}"
            )

        core = _core()
        trees: list[str] = []
        state = {"capped": False}
        slots = [(node, site) for node in core.walk() for site in self.sites.get(node.kind, ())]
        self._place(core, slots, 0, budget, trees, limit, state)

        seen: dict[str, GlycanCandidate] = {}
        rejected_by_rule: dict[str, int] = {}
        rejected = duplicates = 0
        for rendered in trees:
            graph = GlycanGraph.from_iupac_condensed(rendered)
            broken = self._broken_rules(graph)
            if broken:
                rejected += 1
                for reason in broken:  # a tree may break several; each is credited
                    rejected_by_rule[reason] = rejected_by_rule.get(reason, 0) + 1
                continue
            key = graph.canonical_key()
            if key in seen:
                duplicates += 1
                continue
            seen[key] = self._candidate(graph, rendered, key)

        refusal = None
        if not seen:
            reasons = tuple(n_glycan_implausibility_reasons(composition))
            if reasons:
                refusal = "; ".join(reasons)
            elif state["capped"]:
                # THE BUDGET, NEVER THE RULES. This branch used to fall through to "every one of
                # the N arrangements broke a biosynthetic rule", which asserts a fact about
                # glycobiology on the strength of a search that stopped early. Measured:
                # Hex6HexNAc5NeuAc2 gave that refusal at 5,000 trees and 4,620 candidates at
                # 60,000. A budget exhaustion and a biological refusal are different findings and
                # must never be reported as each other.
                refusal = (
                    f"the tree budget of {limit:,} arrangements was exhausted before any candidate"
                    f" passed the curated rules. {rejected:,} of the {len(trees):,} arrangements"
                    " examined broke a rule, and the rest of the space was NOT examined. THIS IS"
                    " NOT A STATEMENT THAT NO CANDIDATE EXISTS: a larger budget may find some, and"
                    " for some compositions it does. The arrangement count grows combinatorially"
                    " with the antennae, so this is a stated scope limit rather than a number that"
                    " can be raised until it goes away"
                )
            elif rejected:
                refusal = (
                    f"every one of the {rejected} arrangements broke a biosynthetic rule, and the"
                    " space was searched EXHAUSTIVELY - the tree budget was not reached, so this is"
                    " a finding about the composition and not about the search"
                )
            elif not trees:
                refusal = "no arrangement of these residues fits the mammalian sites"
        return EnumerationResult(
            composition=composition,
            candidates=tuple(seen.values()),
            built=len(trees),
            rejected=rejected,
            rejected_by_rule=rejected_by_rule,
            duplicates=duplicates,
            capped=state["capped"],
            refusal=refusal,
        )

    def _place(
        self,
        core: Node,
        slots: list[tuple[Node, Site]],
        index: int,
        budget: dict[Residue, int],
        out: list[str],
        limit: int,
        state: dict[str, bool],
    ) -> None:
        """Decide each slot once, in order, so every distinct tree is built exactly once."""
        if len(out) >= limit:
            state["capped"] = True  # the walk really did stop early, rather than land on the limit
            return
        if index == len(slots):
            if not any(budget.values()):
                out.append(core.render())
            return
        node, site = slots[index]
        self._place(core, slots, index + 1, budget, out, limit, state)  # leave it empty
        residue = RESIDUE_CLASS[site.residue]
        if budget.get(residue, 0) <= 0:
            return
        if any(linkage == site.linkage for linkage, _ in node.children):
            return  # the position is taken
        child = Node(site.residue, site.kind)
        node.children.append((site.linkage, child))
        budget[residue] -= 1
        added = [(child, option) for option in self.sites.get(child.kind, ())]
        slots.extend(added)
        self._place(core, slots, index + 1, budget, out, limit, state)
        del slots[len(slots) - len(added):]
        budget[residue] += 1
        node.children.pop()

    def broken_rules(self, graph: GlycanGraph) -> tuple[str, ...]:
        """Public: every curated rule this structure breaks, in words.

        Exists so the ranker can VERIFY that a candidate breaks no rule rather than trusting
        that the enumerator rejected the violators. Those are different claims, and the second
        one is the sort that stays true right up until it does not.
        """
        return self._broken_rules(graph)

    @property
    def rules_the_check_evaluates(self) -> int:
        """How many curated rules `_broken_rules` actually consults.

        NOT the scheme size, and not the N-linked subset either: order rules cannot exclude a
        finished structure, so `_broken_rules` skips them and they are reported as caveats instead.
        Derived from the same list that method walks, and by the same predicate, so a response
        saying "N evaluated" cannot disagree with what was evaluated.

        The response used to assert "0 violated, verified" beside the SCHEME size, which overstated
        the check by eight rules - four O-glycan rows that never reach N-glycan enumeration and four
        order rules that are never evaluated as violations at all.
        """
        return sum(1 for rule, _product, _contexts in self._rules if not is_order_constraint(rule))

    @property
    def rules_skipped_as_order(self) -> int:
        """Curated rules reported as ordering caveats rather than evaluated as violations."""
        return sum(1 for rule, _product, _contexts in self._rules if is_order_constraint(rule))

    def _broken_rules(self, graph: GlycanGraph) -> tuple[str, ...]:
        """Every curated rule this structure breaks, in words. Empty if it breaks none.

        A rule about the order enzymes act in cannot exclude a finished
        structure, so it is not consulted here; see is_order_constraint.
        """
        broken: list[str] = []
        for rule, product, contexts in self._rules:
            if is_order_constraint(rule):
                continue
            anchors = _nodes_matching(graph, product)
            if not anchors:
                continue
            present = any(_context_holds(graph, product, context, anchors) for context in contexts)
            if rule.kind is ConstraintKind.FORBIDS and present:
                broken.append(f"{rule.enzyme}: {rule.product} forbidden in {' or '.join(rule.contexts)}")
            elif rule.kind is ConstraintKind.REQUIRES and not present:
                broken.append(f"{rule.enzyme}: {rule.product} requires {' or '.join(rule.contexts)}")
        return tuple(broken)

    def _candidate(self, graph: GlycanGraph, rendered: str, key: str) -> GlycanCandidate:
        structure = GlycanStructure(
            composition=graph.composition(),
            iupac_condensed=rendered,
            has_unresolved_linkage=graph.has_unresolved_linkage,
            has_unresolved_anomericity=graph.has_unresolved_anomericity,
            source=CANDIDATE_SOURCE,
            reuse_status=ReuseStatus.INTERNAL_PROPRIETARY,
        )
        reasons: list[CandidateReason] = []
        for parent, child, data in graph.graph.edges(data=True):
            residue = graph.graph.nodes[child]["name"]
            linkage = str(data["linkage"])
            enzymes = self.catalogue.enzymes_for(residue, linkage)
            if enzymes:
                reasons.append(
                    CandidateReason(
                        kind="enzyme",
                        text=(
                            f"the monolink {residue}({linkage}) is one human enzymes are known to make;"
                            f" this candidate places it on {graph.graph.nodes[parent]['name']}"
                        ),
                        enzymes=enzymes,
                    )
                )
        # THE ORDERING NOTE IS GATED ON THE CONTEXT, NOT ONLY ON THE PRODUCT.
        #
        # This loop used to attach the note whenever the rule's PRODUCT was present, and the
        # note's own words are "This candidate carries {product} alongside that context". For a
        # candidate that carries the product and NOT the context, that sentence is false - and
        # it is false while carrying an enzyme and a literature citation, which is the form a
        # reader trusts most. Measured before the fix, on Hex5HexNAc4Fuc1: 114 of 167 candidates
        # carried an ordering note and only 18 contained the bisecting GlcNAc that every one of
        # the four order rules names as its context, so 108 candidates published a cited claim
        # about themselves that was not true of them.
        #
        # THE CONTEXT IS SEARCHED ANYWHERE IN THE MOLECULE, not anchored to the product's own
        # residue, and the difference matters. `_context_holds` anchors the two fragments to one
        # residue when they end at the same residue-and-linkage pair, which is right for MGAT4 -
        # both its fragments name the SAME alpha1-3 mannose. It is wrong for all four ORDER
        # rules, whose context is the bisecting GlcNAc: FUT8's product `Fuc(a1-6)GlcNAc` ends at
        # ('GlcNAc', None) and the bisecting context `GlcNAc(b1-4)Man(b1-4)GlcNAc` also ends at
        # ('GlcNAc', None), so the anchoring treats the REDUCING GlcNAc and the CHITOBIOSE
        # GlcNAc as one residue, the intersection is empty, and FUT8's note can never fire - on
        # exactly the bisected, core-fucosylated structures the note exists to explain. There
        # are 321 of them in the reference set.
        #
        # Searched anywhere is also the right reading of what an order rule says. "Bisecting
        # GlcNAc blocks SUBSEQUENT MGAT2 action" is a claim about the state of the molecule, not
        # about two residues sharing a position: once the bisecting GlcNAc is there, that enzyme
        # can no longer act, wherever it would have acted.
        #
        # This is deliberately confined to the rationale and never reaches `_broken_rules`,
        # which skips order rules entirely. So no candidate is included or excluded by this
        # change - only what a candidate says about itself. The candidate counts are unmoved:
        # 10 / 63 / 167 / 6 for the four reported compositions, asserted in the tests.
        #
        # The ANCHORING DEFECT ITSELF IS NOT FIXED HERE. `_context_holds` still identifies two
        # different GlcNAc as one residue, and it still gates `_broken_rules` for the non-order
        # rules, where changing it could change which candidates are refused. That is a
        # rule-semantics decision with a candidate-count consequence, and it is recorded in
        # docs/GLYCAN_LIMITATIONS.md rather than made as a side effect of building a ranker.
        #
        # Where the product is present and the context is not, the rule is not an ordering
        # caveat at all - it is a FORBIDS rule this candidate satisfies non-vacuously, because
        # the thing that would forbid the product is absent. So it is reported as an ordinary
        # rule reason, which is how every other satisfied FORBIDS rule is already reported.
        for rule, product, contexts in self._rules:
            anchors = _nodes_matching(graph, product)
            if not anchors:
                continue
            if is_order_constraint(rule) and any(
                _nodes_matching(graph, context) for context in contexts
            ):
                reasons.append(
                    CandidateReason(
                        kind="ordering",
                        text=(
                            f"{rule.explain()} This candidate carries {rule.product} alongside that context,"
                            " which the permitted order reaches."
                        ),
                        enzymes=(rule.enzyme,),
                        reference=rule.reference,
                    )
                )
            else:
                reasons.append(
                    CandidateReason(
                        kind="rule", text=rule.explain(), enzymes=(rule.enzyme,), reference=rule.reference
                    )
                )
        return GlycanCandidate(structure=structure, canonical_key=key, rationale=tuple(reasons))


def _chain(fragment: str) -> tuple[tuple[str, str | None], ...]:
    """An IUPAC fragment as residues from the outside in, e.g. "GlcNAc(b1-2)Man(a1-3)".

    A branched fragment is refused rather than flattened: dropping the brackets
    would turn "GlcNAc(b1-2)[GlcNAc(b1-4)]Man(a1-3)", two residues on one
    mannose, into a three-residue stack and match the wrong structures both ways.
    No rule in the shipped table is branched; one that appears must be handled
    rather than silently misread.
    """
    if "[" in fragment or "]" in fragment:
        raise EnumerationError(f"cannot read the branched fragment {fragment!r}; only a single chain is supported")
    parts = [(match.group(1), match.group(2)) for match in _FRAGMENT.finditer(fragment) if match.group(1)]
    return tuple(parts)


def _nodes_matching(graph: GlycanGraph, chain: Sequence[tuple[str, str | None]]) -> set[int]:
    """Every node that the chain's innermost residue can sit on, the chain holding outwards."""
    if not chain:
        return set()
    deepest, arrival = chain[-1]
    found: set[int] = set()
    for node, data in graph.graph.nodes(data=True):
        if data["name"] != deepest:
            continue
        # A fragment ending in a linkage, such as "Man(a1-3)Man(a1-6)", says how its
        # innermost residue attaches to whatever carries it. Skipping that check lets
        # the core beta-mannose, which arrives on b1-4, stand in for an arm mannose.
        if arrival is not None and not _arrives_on(graph, node, arrival):
            continue
        if _matches(graph, node, chain, len(chain) - 1):
            found.add(node)
    return found


def _contains(graph: GlycanGraph, chain: Sequence[tuple[str, str | None]]) -> bool:
    """True if the structure holds that chain of residues, joined by those linkages."""
    return bool(_nodes_matching(graph, chain))


def _context_holds(
    graph: GlycanGraph,
    product: Sequence[tuple[str, str | None]],
    context: Sequence[tuple[str, str | None]],
    anchors: set[int],
) -> bool:
    """Whether the context is present, on the same residue as the product where that is what is meant.

    MGAT4 requires GlcNAc(b1-4)Man(a1-3) to sit on a mannose that already carries
    GlcNAc(b1-2)Man(a1-3): both fragments end at the same residue arriving the
    same way, so they speak about one mannose and the match must be anchored to
    it. MGAT2 requires GlcNAc(b1-2)Man(a1-6) given GlcNAc(b1-2)Man(a1-3): those
    end at different arms and are two residues by construction, so the context is
    searched for anywhere.
    """
    found = _nodes_matching(graph, context)
    if not found:
        return False
    if product and context and product[-1] == context[-1]:
        return bool(found & anchors)
    return True


def _matches(graph: GlycanGraph, node: int, chain: Sequence[tuple[str, str | None]], index: int) -> bool:
    if index == 0:
        return True
    residue, linkage = chain[index - 1]
    for child in graph.graph.successors(node):
        if graph.graph.nodes[child]["name"] != residue:
            continue
        if linkage is not None and str(graph.graph.edges[node, child]["linkage"]) != linkage:
            continue
        if _matches(graph, child, chain, index - 1):
            return True
    return False


def _arrives_on(graph: GlycanGraph, node: int, linkage: str) -> bool:
    """True if `node` is attached to its parent by that linkage.

    A residue written in braces states a linkage but not what carries it, so
    where it sits is unknown. Unknown is not a match: accepting it would let a
    fragment of undetermined position discharge a rule's precondition, which is
    an unresolved structure being quietly treated as a resolved one.
    """
    parents = list(graph.graph.predecessors(node))
    if not parents:
        return False
    return str(graph.graph.edges[parents[0], node]["linkage"]) == linkage


def contains_motif(graph: GlycanGraph, fragment: str) -> bool:
    """True if the structure holds that IUPAC fragment, such as "Fuc(a1-6)GlcNAc".

    Read what this does and does not say. The fragment is matched from the
    outside in, and the INNERMOST residue's own arrival linkage is checked only
    where the fragment states one. "Fuc(a1-6)GlcNAc" therefore means "some
    GlcNAc, anywhere, carrying an a1-6 fucose": it is true of an antennary
    GlcNAc, and true of a GlcNAc inside a braced fragment whose position is
    unknown. It is NOT a test for core fucose, which is a claim about the
    reducing-end residue specifically, and this function cannot express that
    because it has no notion of where the chain ends.

    A bare residue name is weaker still: "Fuc" is satisfied by a fucose whose
    position is entirely unknown. So a motif flag built on a bare name claims a
    placement the structure never stated, and a flag built on a fragment claims
    only the chain it spells. A caller wanting a positional claim, such as core
    fucose or a bisecting GlcNAc, must anchor it against the structure itself
    rather than trust a substring: a string pattern is not a biological claim.

    Raises EnumerationError for a branched fragment, which cannot be read as a
    single chain.
    """
    return _contains(graph, _chain(fragment))


def motif_anchors(graph: GlycanGraph, fragment: str) -> tuple[int, ...]:
    """Every node the fragment's innermost residue sits on, in ascending order.

    Empty when the motif is absent, so it answers "where" where contains_motif
    answers "whether".
    """
    return tuple(sorted(_nodes_matching(graph, _chain(fragment))))
