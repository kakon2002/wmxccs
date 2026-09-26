"""Glycan structures as graphs: monosaccharides as nodes, linkages as edges.

The strings are parsed by glycowork, which writes each linkage as a node of its
own; here those nodes are collapsed into the edges between the residues they
join. Edges run from the residue nearer the reducing end to the residue it
carries.

What a string leaves open stays open. An ambiguous position such as "a2-3/6"
keeps both positions rather than picking one, "?" stays unknown, and a residue
written in braces, whose attachment point is undetermined, becomes a node with
no edge at all.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from typing import Iterator

import networkx as nx

from .composition import Composition, composition_from_residue_names

# A linkage as glycowork writes it: "a1-3", "b1-4", "a2-3/6", "?1-?", and the
# less common "b6-2", "a2-11" and "1-5", where the anomer is left out entirely.
# Residue names always carry letters, so no residue can match this.
_LINKAGE = re.compile(r"(?P<anomer>[ab?])?(?P<donor>[0-9?]+)-(?P<acceptors>[0-9?]+(?:/[0-9?]+)*)")


class GlycanGraphError(ValueError):
    """A structure string that cannot be turned into a graph."""


def is_linkage_label(label: str) -> bool:
    """True if a glycowork node label is a linkage rather than a monosaccharide."""
    return _LINKAGE.fullmatch(label) is not None


@dataclass(frozen=True)
class Linkage:
    """A glycosidic linkage as written, with whatever it leaves open left open."""

    raw: str
    anomer: str | None  # "a", "b", or None where the string says "?"
    donor: int | None  # the anomeric carbon of the residue being carried
    acceptors: tuple[int, ...]  # positions on the residue carrying it; several when ambiguous, none when "?"
    acceptor_unknown: bool = False  # the string offered "?" as one of the alternatives, as in "a1-3/?"

    @classmethod
    def parse(cls, text: str) -> Linkage:
        match = _LINKAGE.fullmatch(text)
        if match is None:
            raise GlycanGraphError(f"cannot read {text!r} as a linkage such as a1-3, a2-3/6 or ?1-?")
        alternatives = match["acceptors"].split("/")
        return cls(
            raw=text,
            anomer=match["anomer"] if match["anomer"] != "?" else None,
            donor=int(match["donor"]) if match["donor"] != "?" else None,
            acceptors=tuple(int(p) for p in alternatives if p != "?"),
            acceptor_unknown=any(p == "?" for p in alternatives),
        )

    @property
    def canonical(self) -> str:
        """One spelling per linkage, so two writings of the same ambiguity agree.

        "a2-3/6" and "a2-6/3" say the same thing; keyed on the raw text they
        would look like different structures and a duplicate would escape.
        """
        anomer = self.anomer or "?"
        donor = self.donor if self.donor is not None else "?"
        acceptors = [str(position) for position in sorted(self.acceptors)]
        if self.acceptor_unknown:
            acceptors.append("?")
        return f"{anomer}{donor}-{'/'.join(acceptors) if acceptors else '?'}"

    @property
    def anomer_known(self) -> bool:
        return self.anomer is not None

    @property
    def position_known(self) -> bool:
        """True only when both ends are pinned down: one donor carbon and one acceptor position.

        An alternative of "?" leaves the position open even beside a stated one,
        so "a1-3/?" is not resolved.
        """
        return self.donor is not None and len(self.acceptors) == 1 and not self.acceptor_unknown

    def __str__(self) -> str:
        return self.raw


class GlycanGraph:
    """A parsed glycan: residues as nodes, linkages as edges, rooted at the reducing end."""

    def __init__(self, graph: nx.DiGraph, iupac_condensed: str) -> None:
        self.graph = graph
        self.iupac_condensed = iupac_condensed

    @classmethod
    def from_iupac_condensed(cls, text: str) -> GlycanGraph:
        """Parse an IUPAC-condensed string such as "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc"."""
        if not isinstance(text, str) or not text.strip():
            raise GlycanGraphError("empty glycan string")
        parsed = _parse(text)
        labels = {index: data["string_labels"] for index, data in parsed.nodes(data=True)}
        graph = nx.DiGraph()
        for index, label in labels.items():
            if not is_linkage_label(label):
                graph.add_node(index, name=label, attached=True, linkage=None)
        dangling: dict[int, Linkage] = {}  # residue -> the linkage it was written with, its other end unknown
        for index, label in labels.items():
            if not is_linkage_label(label):
                continue
            linkage = Linkage.parse(label)
            carried = [n for n in parsed.successors(index) if n in graph]
            carriers = [n for n in parsed.predecessors(index) if n in graph]
            if len(carried) != 1 or len(carriers) > 1:
                raise GlycanGraphError(
                    f"linkage {label!r} joins {len(carriers)} residues to {len(carried)} in {text!r}"
                )
            if not carriers:  # written in braces: the linkage is given, the residue it attaches to is not
                dangling[carried[0]] = linkage
                continue
            graph.add_edge(carriers[0], carried[0], linkage=linkage)
        if not graph:
            raise GlycanGraphError(f"no monosaccharides in {text!r}")

        # Braces can hold a whole fragment, not just one residue. A fragment hangs off
        # the structure at an unknown point, so every residue in it is unattached.
        main: list[set[int]] = []
        for component in nx.weakly_connected_components(graph):
            roots = [node for node in component if graph.in_degree(node) == 0]
            if len(roots) != 1:
                raise GlycanGraphError(f"a fragment of {text!r} has {len(roots)} roots, not one")
            root = roots[0]
            if root in dangling:
                for node in component:
                    graph.nodes[node]["attached"] = False
                graph.nodes[root]["linkage"] = dangling[root]
            else:
                main.append(component)
        if len(main) != 1:
            raise GlycanGraphError(f"{text!r} holds {len(main)} attached structures, not one")
        if not nx.is_tree(graph.subgraph(main[0])):
            raise GlycanGraphError(f"the attached residues of {text!r} do not form a tree")
        return cls(graph, text)

    @property
    def residue_names(self) -> tuple[str, ...]:
        """Every monosaccharide, floating ones included, in the order glycowork numbered them."""
        return tuple(data["name"] for _, data in sorted(self.graph.nodes(data=True)))

    @property
    def unattached_residues(self) -> tuple[str, ...]:
        """Residues written in braces: known to be there, not known where."""
        return tuple(data["name"] for _, data in sorted(self.graph.nodes(data=True)) if not data["attached"])

    @property
    def root(self) -> int | None:
        """The reducing-end residue, or None if every residue is floating."""
        roots = [
            node
            for node, data in self.graph.nodes(data=True)
            if data["attached"] and self.graph.in_degree(node) == 0
        ]
        return roots[0] if len(roots) == 1 else None

    def linkages(self) -> Iterator[Linkage]:
        """Every linkage, including those of floating residues, whose other end is unknown."""
        for _, _, data in self.graph.edges(data=True):
            yield data["linkage"]
        for _, data in self.graph.nodes(data=True):
            if data["linkage"] is not None:
                yield data["linkage"]

    @property
    def has_unresolved_linkage(self) -> bool:
        """True if any position is ambiguous or unknown, or any residue floats."""
        if self.unattached_residues:
            return True
        return any(not linkage.position_known for linkage in self.linkages())

    @property
    def has_unresolved_anomericity(self) -> bool:
        return any(not linkage.anomer_known for linkage in self.linkages())

    def composition(self) -> Composition:
        """Composition of the whole structure, floating residues included.

        Raises UnsupportedResidueError if any residue has no place among the
        five classes.
        """
        return composition_from_residue_names(self.residue_names)

    def canonical_key(self) -> str:
        """A form that two graphs share exactly when they are the same structure.

        Branches are sorted, so the same tree written two ways gives one key.
        Floating fragments are kept apart from the attached structure and
        sorted among themselves, so a residue of unknown attachment is never
        confused with one that is attached.
        """

        def subtree(node: int) -> str:
            branches = sorted(
                f"{self.graph.edges[node, child]['linkage'].canonical}:{subtree(child)}"
                for child in self.graph.successors(node)
            )
            name = self.graph.nodes[node]["name"]
            return f"{name}({','.join(branches)})" if branches else name

        roots = [node for node in self.graph if self.graph.in_degree(node) == 0]
        attached = [node for node in roots if self.graph.nodes[node]["attached"]]
        floating = sorted(
            f"{self.graph.nodes[node]['linkage'].canonical}:{subtree(node)}"
            for node in roots
            if not self.graph.nodes[node]["attached"]
        )
        key = subtree(attached[0]) if attached else ""
        return key + ("{" + ";".join(floating) + "}" if floating else "")

    def __len__(self) -> int:
        return self.graph.number_of_nodes()

    def __repr__(self) -> str:
        return f"GlycanGraph({self.iupac_condensed!r})"


@lru_cache(maxsize=1)
def _glycowork_parser():
    """glycowork's IUPAC-condensed parser, imported on first use: it pulls in pandas."""
    from glycowork.motif.graph import glycan_to_nxGraph

    return glycan_to_nxGraph


def _parse(text: str) -> nx.DiGraph:
    try:
        return _glycowork_parser()(text)
    except GlycanGraphError:
        raise
    except Exception as exc:  # glycowork raises whatever the string happens to break
        raise GlycanGraphError(f"glycowork could not parse {text!r} ({type(exc).__name__}: {exc})") from exc
