"""What a model would see: one record turned into a fixed row of numbers.

Three declared blocks, in a fixed width and a fixed order:

- the analyte block, from the record's own validated composition;
- the condition block, from the measurement conditions;
- the structure block, from the IUPAC-condensed string where there is one.

Only the analyte and condition blocks are fitted by default, because those are
exactly what a prediction request can supply. The structure block is extracted,
reported and available for diagnostics, but fitting a column that is always
absent at inference while populated in training is a defect rather than a
missing value, so turning it on is a deliberate act with a contract change
behind it.

Nothing here is fitted, learned or shared between records. A feature vector for
one record is the same whether it is extracted alone or inside a batch, which is
what makes it safe to featurise before splitting: a statistic taken over the
whole dataset would carry test information into training.

Missing is None, and at the wire it becomes NaN. It is never 0. Zero means the
extractor looked and found none; None means it could not look, and the two must
not be the same number. An ambiguity is likewise never resolved into one of its
alternatives: a sialic acid written a2-3/6 counts as ambiguous and as neither
position.

None of these features has been checked against a single CCS measurement. The
whole set is an untested physical hypothesis: collision cross section is a
rotationally averaged three-dimensional projected area, and everything here is
two-dimensional graph and count data.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import networkx as nx

from .composition import Composition, has_full_man3glcnac2_core
from .glycan_graph import GlycanGraph, GlycanGraphError
from .models import Derivatisation, DriftGas, IMSType, Polarity, ReducingEndLabel, adduct_components

FITTED_ANALYTE_FEATURES = (
    "hex",
    "hexnac",
    "fuc",
    "neuac",
    "neugc",
    "residue_total",
    "monoisotopic_mass",
    "antenna_hexnac",
    "has_full_core",
)

FITTED_CONDITION_FEATURES = (
    "charge_magnitude",
    "polarity_positive",
    "adduct_protons",
    "adduct_sodium",
    "adduct_other_additions",
    "adduct_losses",
    "ccs_is_calibrated",
    "drift_gas_he",
    "ims_type",
    "reducing_end_label",
    "derivatisation",
)

STRUCTURE_FEATURES = (
    "residues_attached",
    "unplaced_residues",
    "bonds",
    "max_depth",
    "mean_depth",
    "branch_points",
    "max_out_degree",
    "terminal_residues",
    "wiener_index",
    "mean_path_length",
    "arm_a3_size",
    "arm_a6_size",
    "arm_asymmetry",
    "core_fucose",
    "bisecting_glcnac",
    "linkages_total",
    "alpha_linkages",
    "beta_linkages",
    "anomer_unknown",
    "position_unknown",
    "flexible_1_6_linkages",
    "sia_a2_3",
    "sia_a2_6",
    "sia_position_ambiguous",
)

FEATURE_NAMES = FITTED_ANALYTE_FEATURES + FITTED_CONDITION_FEATURES + STRUCTURE_FEATURES

# All three blocks are fitted. The model scores an enumerated CANDIDATE, not the
# request the user sent, so a structure is available at inference and the
# structure block carries the only signal in the whole extractor that can
# separate two isomers of one composition. Without it a ranking layer has
# nothing to rank on.
FITTED_FEATURES = FEATURE_NAMES

# What a bare prediction request can supply by itself. A request is no longer a
# complete model input: it becomes one once a candidate has been enumerated for
# it. Kept named because the distinction decides which split mode is honest.
DEPLOYMENT_SUPPLIED_FEATURES = FITTED_ANALYTE_FEATURES + FITTED_CONDITION_FEATURES

REQUIRES_STRUCTURE = frozenset(STRUCTURE_FEATURES)
CATEGORICAL_FEATURES = frozenset({"ims_type", "reducing_end_label", "derivatisation"})

# Declared here rather than taken from enum declaration order, so that adding a
# member upstream cannot silently renumber an existing code and change what a
# stored model means. A test asserts these cover every member.
_IMS_ORDER = (IMSType.DTIMS, IMSType.TWIMS, IMSType.TIMS, IMSType.CYCLIC)
_LABEL_ORDER = (
    ReducingEndLabel.NATIVE,
    ReducingEndLabel.REDUCED,
    ReducingEndLabel.PA,
    ReducingEndLabel.TWO_AB,
    ReducingEndLabel.TWO_AA,
    ReducingEndLabel.PROCAINAMIDE,
    ReducingEndLabel.RAPIFLUOR_MS,
    ReducingEndLabel.APTS,
    ReducingEndLabel.OTHER,
    ReducingEndLabel.UNKNOWN,
)
_DERIVATISATION_ORDER = (
    Derivatisation.UNDERIVATISED,
    Derivatisation.PERMETHYLATION,
    Derivatisation.SIALIC_ACID_AMIDATION,
    Derivatisation.SIALIC_ACID_ESTERIFICATION,
    Derivatisation.UNKNOWN,
)

# Fields deliberately kept out of the feature set. The reason ships as code, so
# it cannot drift away from the decision the way a comment can.
EXCLUDED: Mapping[str, str] = {
    "ccs": "the target; reading it would be leakage of the most direct kind",
    "source": "a provenance string that stands in for the laboratory, so a model would learn the lab",
    "instrument": "names a model of machine, which on a small corpus is very nearly a laboratory identifier",
    "cell_gas": "provenance only; the gas the value refers to is already carried as drift_gas",
    "reuse_status": "a licensing fact about the record, with no bearing on the ion",
    "glytoucan_ac": "an accession identifies the structure, so fitting it is memorisation",
    "wurcs": "as above: an identifier, not a property of the ion",
    "iupac_condensed": "the structure string is read to build the graph, never fitted as text",
    "calibrant": "a genuine physical determinant of a calibrated value, and also very nearly a laboratory name",
}

# Failure reasons, as module constants so a caller matches a name rather than prose.
STRUCTURE_PARSE_FAILED = "structure did not parse"
UNSUPPORTED_RESIDUES = "structure holds residues outside the five composition classes"
COMPOSITION_DISAGREES = "the structure and the recorded composition disagree"


@dataclass(frozen=True)
class FeatureVector:
    """One record as a fixed-width row, with what could not be determined left as None."""

    values: tuple[float | None, ...]
    resolution_gaps: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if len(self.values) != len(FEATURE_NAMES):
            raise ValueError(f"a feature vector holds {len(FEATURE_NAMES)} values, not {len(self.values)}")

    def as_mapping(self) -> dict[str, float | None]:
        return dict(zip(FEATURE_NAMES, self.values))

    def as_row(self, names: Sequence[str] = FITTED_FEATURES) -> tuple[float, ...]:
        """The named features as plain floats, with None becoming NaN at the wire."""
        mapping = self.as_mapping()
        return tuple(math.nan if mapping[name] is None else float(mapping[name]) for name in names)

    def analyte_row(self) -> tuple[float, ...]:
        """The analyte half of the fitted row: what a model can see of the molecule itself.

        This is the row the leakage law is checked over. Two measurements of one
        glycan under different conditions share it, which is the point: they are
        one analyte and must not be separated by a split.
        """
        return self.as_row(FITTED_ANALYTE_FEATURES)

    def known_mask(self, names: Sequence[str] = FEATURE_NAMES) -> tuple[bool, ...]:
        mapping = self.as_mapping()
        return tuple(mapping[name] is not None for name in names)


def _code(value: object, order: Sequence[object]) -> float | None:
    for index, member in enumerate(order):
        if value == member:
            return float(index)
    return None


def extract_composition(composition: Composition) -> dict[str, float | None]:
    """The analyte block, from a validated composition.

    Always taken from the record's own composition, never from the graph: a
    structure carrying a residue outside the five classes still parses into a
    graph, but asking that graph for a composition raises.
    """
    total = composition.hex + composition.hexnac + composition.fuc + composition.neuac + composition.neugc
    return {
        "hex": float(composition.hex),
        "hexnac": float(composition.hexnac),
        "fuc": float(composition.fuc),
        "neuac": float(composition.neuac),
        "neugc": float(composition.neugc),
        "residue_total": float(total),
        "monoisotopic_mass": composition.monoisotopic_mass,
        # HexNAc beyond the two core residues: an antenna count only where a core exists.
        "antenna_hexnac": float(max(composition.hexnac - 2, 0)),
        "has_full_core": float(has_full_man3glcnac2_core(composition)),
    }


def extract_conditions(conditions: object) -> dict[str, float | None]:
    """The condition block, from anything carrying measurement conditions.

    Read by attribute rather than by type, so a prediction request and a
    measurement go through one path.
    """
    parts = adduct_components(getattr(conditions, "adduct"))
    protons = sum(count for sign, count, species in parts if sign == "+" and species == "H")
    sodium = sum(count for sign, count, species in parts if sign == "+" and species == "Na")
    other = sum(count for sign, count, species in parts if sign == "+" and species not in {"H", "Na"})
    losses = sum(count for sign, count, _ in parts if sign == "-")
    charge = getattr(conditions, "charge")
    return {
        "charge_magnitude": float(abs(charge)),
        "polarity_positive": float(getattr(conditions, "polarity") == Polarity.POSITIVE),
        "adduct_protons": float(protons),
        "adduct_sodium": float(sodium),
        "adduct_other_additions": float(other),
        "adduct_losses": float(losses),
        # Whether a calibrant entered the value at all; the calibrant's identity is excluded.
        "ccs_is_calibrated": float(getattr(conditions, "ccs_is_calibrated")),
        "drift_gas_he": float(getattr(conditions, "drift_gas") == DriftGas.HE),
        "ims_type": _code(getattr(conditions, "ims_type"), _IMS_ORDER),
        "reducing_end_label": _code(getattr(conditions, "reducing_end_label"), _LABEL_ORDER),
        "derivatisation": _code(getattr(conditions, "derivatisation"), _DERIVATISATION_ORDER),
    }


def _attached(graph: GlycanGraph) -> nx.DiGraph:
    """The part of the structure whose placement is known.

    A residue written in braces is known to be present but not where, so it
    counts toward the composition and the mass and toward unplaced_residues, and
    toward no shape feature at all. This is not cosmetic: such a graph is not a
    tree, and a walk from the root silently skips the floating residue anyway.
    """
    return graph.graph.subgraph([node for node, data in graph.graph.nodes(data=True) if data["attached"]])


def _child_on(attached: nx.DiGraph, node: int | None, name: str, linkage: str) -> bool:
    """True if `node` carries a residue of that name on exactly that linkage.

    Anchored to a named node, so it states WHERE. A fragment match cannot: the
    innermost residue of "Fuc(a1-6)GlcNAc" has no arrival linkage, so the
    fragment is satisfied by any GlcNAc anywhere, an antennary one and a
    floating one included.
    """
    if node is None:
        return False
    return any(
        attached.nodes[child]["name"] == name and str(attached.edges[node, child]["linkage"]) == linkage
        for child in attached.successors(node)
    )


def _reducing_end(attached: nx.DiGraph) -> int | None:
    """The residue everything else hangs from, or None if that is not one node."""
    roots = [node for node in attached if attached.in_degree(node) == 0]
    return roots[0] if len(roots) == 1 else None


def _beta_mannose(attached: nx.DiGraph) -> int | None:
    """The core beta-mannose: b1-4 on the chitobiose, carrying both arms. None if absent."""
    for node, data in attached.nodes(data=True):
        if data["name"] != "Man":
            continue
        parents = list(attached.predecessors(node))
        if not parents or attached.nodes[parents[0]]["name"] != "GlcNAc":
            continue
        if str(attached.edges[parents[0], node]["linkage"]) != "b1-4":
            continue
        grand = list(attached.predecessors(parents[0]))
        if not grand or attached.nodes[grand[0]]["name"] != "GlcNAc":
            continue
        if str(attached.edges[grand[0], parents[0]]["linkage"]) != "b1-4":
            continue
        arms = {str(attached.edges[node, child]["linkage"]) for child in attached.successors(node)}
        if "a1-3" in arms and "a1-6" in arms:
            return node
    return None


def _arm_sizes(attached: nx.DiGraph, beta: int | None) -> tuple[float | None, float | None]:
    """Residues carried on the alpha1-3 and alpha1-6 arms.

    The only features in the whole extractor that can see the arm-isomer axis:
    two glycans of one composition differing only in which arm carries what have
    the same everything else and differ here.
    """
    if beta is None:
        return None, None
    sizes: dict[str, float] = {}
    for child in attached.successors(beta):
        linkage = str(attached.edges[beta, child]["linkage"])
        if linkage in {"a1-3", "a1-6"}:
            sizes[linkage] = float(len(nx.descendants(attached, child)) + 1)
    return sizes.get("a1-3"), sizes.get("a1-6")


def extract_structure(iupac_condensed: str) -> dict[str, float | None]:
    """The structure block. Raises GlycanGraphError if the string does not parse."""
    graph = GlycanGraph.from_iupac_condensed(iupac_condensed)
    attached = _attached(graph)
    unplaced = graph.graph.number_of_nodes() - attached.number_of_nodes()

    root = _reducing_end(attached)
    depths = nx.single_source_shortest_path_length(attached, root) if root is not None else {}
    undirected = attached.to_undirected()
    paths = [
        length
        for _, lengths in nx.all_pairs_shortest_path_length(undirected)
        for length in lengths.values()
    ]
    pairs = [length for length in paths if length > 0]

    alpha = beta = anomer_unknown = position_unknown = flexible = 0
    sia_a2_3 = sia_a2_6 = sia_ambiguous = 0
    for parent, child, data in attached.edges(data=True):
        linkage = data["linkage"]
        if linkage.anomer == "a":
            alpha += 1
        elif linkage.anomer == "b":
            beta += 1
        else:
            anomer_unknown += 1
        if not linkage.position_known:
            position_unknown += 1
        if linkage.position_known and linkage.donor == 1 and linkage.acceptors == (6,):
            # A 1-6 bond has an extra rotatable torsion, so it matters to shape.
            # An unknown position is not counted: "a1-6/?" says it might be here.
            flexible += 1
        if attached.nodes[child]["name"] in {"Neu5Ac", "Neu5Gc"}:
            # Neither the position nor the anomer may be assumed. An ambiguous
            # a2-3/6 is neither a2-3 nor a2-6, and a stated beta or an unstated
            # anomer must not be counted in a column named for alpha.
            if not linkage.position_known or linkage.anomer != "a":
                sia_ambiguous += 1
            elif linkage.acceptors == (3,):
                sia_a2_3 += 1
            elif linkage.acceptors == (6,):
                sia_a2_6 += 1

    beta_man = _beta_mannose(attached)
    arm_a3, arm_a6 = _arm_sizes(attached, beta_man)
    out_degrees = [attached.out_degree(node) for node in attached]

    return {
        "residues_attached": float(attached.number_of_nodes()),
        "unplaced_residues": float(unplaced),
        "bonds": float(attached.number_of_edges()),
        "max_depth": float(max(depths.values())) if depths else None,
        "mean_depth": (sum(depths.values()) / len(depths)) if depths else None,
        "branch_points": float(sum(1 for degree in out_degrees if degree >= 2)),
        "max_out_degree": float(max(out_degrees)) if out_degrees else None,
        "terminal_residues": float(sum(1 for degree in out_degrees if degree == 0)),
        "wiener_index": float(sum(pairs) / 2) if pairs else 0.0,
        "mean_path_length": (sum(pairs) / len(pairs)) if pairs else None,
        "arm_a3_size": arm_a3,
        "arm_a6_size": arm_a6,
        "arm_asymmetry": None if arm_a3 is None or arm_a6 is None else abs(arm_a3 - arm_a6),
        # Both are POSITIONAL claims and are anchored to the node they are about,
        # on the attached subgraph like every other shape feature. A fragment
        # match cannot express either: "Fuc(a1-6)GlcNAc" is satisfied by an
        # antennary GlcNAc and by one inside a braced fragment of unknown
        # placement, so it would call antennary fucose core fucose.
        "core_fucose": float(_child_on(attached, root, "Fuc", "a1-6")),
        "bisecting_glcnac": float(_child_on(attached, beta_man, "GlcNAc", "b1-4")),
        "linkages_total": float(attached.number_of_edges()),
        "alpha_linkages": float(alpha),
        "beta_linkages": float(beta),
        "anomer_unknown": float(anomer_unknown),
        "position_unknown": float(position_unknown),
        "flexible_1_6_linkages": float(flexible),
        "sia_a2_3": float(sia_a2_3),
        "sia_a2_6": float(sia_a2_6),
        "sia_position_ambiguous": float(sia_ambiguous),
    }


def has_structure(record: object) -> bool:
    """Whether this record can teach a structure feature anything.

    A measurement whose glycan carries only a composition has nothing to say
    about isomers: every structure column comes out None, and fitting on it
    would teach the model that those columns are absent rather than what they
    mean. Such a record is refused from a structure-feature fit, not silently
    filled in.
    """
    glycan = getattr(record, "glycan", None)
    return bool(getattr(glycan, "iupac_condensed", None))


def resolution_gaps(graph: GlycanGraph) -> list[str]:
    """What this structure leaves open, in words."""
    gaps: list[str] = []
    if graph.unattached_residues:
        gaps.append(f"{len(graph.unattached_residues)} residue(s) of unknown attachment")
    if any(not linkage.position_known for linkage in graph.linkages()):
        gaps.append("at least one linkage position is ambiguous or unknown")
    if any(not linkage.anomer_known for linkage in graph.linkages()):
        gaps.append("at least one anomeric configuration is unknown")
    return gaps


def is_fully_resolved(graph: GlycanGraph) -> bool:
    return not resolution_gaps(graph)


def extract(record: object) -> FeatureVector:
    """One record as a feature vector. Structure features are None where there is no structure."""
    glycan = getattr(record, "glycan", None)
    composition = glycan.composition if glycan is not None else getattr(record, "composition")
    values: dict[str, float | None] = {}
    values.update(extract_composition(composition))
    values.update(extract_conditions(record))

    gaps: tuple[str, ...] = ()
    notes: list[str] = []
    text = getattr(glycan, "iupac_condensed", None) if glycan is not None else None
    if text is None:
        values.update({name: None for name in STRUCTURE_FEATURES})
    else:
        values.update(extract_structure(text))
        graph = GlycanGraph.from_iupac_condensed(text)
        gaps = tuple(resolution_gaps(graph))
        try:
            if graph.composition() != composition:
                notes.append(COMPOSITION_DISAGREES)
        except Exception:
            notes.append(UNSUPPORTED_RESIDUES)  # the graph is fine; its residues have no composition class
    return FeatureVector(
        values=tuple(values[name] for name in FEATURE_NAMES),
        resolution_gaps=gaps,
        notes=tuple(notes),
    )


@dataclass(frozen=True)
class ExtractionReport:
    """What a batch produced, and what it could not. Every row lands in exactly one bucket."""

    vectors: tuple[FeatureVector, ...]
    rows_in: int
    records_with_structure: int
    structures_parsed: int
    vectors_built: int
    failure_counts: Mapping[str, int] = field(default_factory=dict)
    failure_examples: Mapping[str, str] = field(default_factory=dict)

    @property
    def rows_failed(self) -> int:
        return sum(self.failure_counts.values())

    def summary(self) -> str:
        lines = [
            f"  rows in                  {self.rows_in}",
            f"  records with a structure {self.records_with_structure}",
            f"  structures parsed        {self.structures_parsed}",
            f"  vectors built            {self.vectors_built}",
        ]
        for reason, count in sorted(self.failure_counts.items(), key=lambda item: -item[1]):
            lines.append(f"    {reason}: {count}  e.g. {self.failure_examples.get(reason, '')}")
        return "\n".join(lines)


def extract_all(records: Iterable[object]) -> ExtractionReport:
    """Featurise a batch, counting every row that did not make it.

    Raises ValueError on an empty batch: nothing to featurise is a caller
    mistake, not a data gap, and it must not be mistaken for a clean run.
    """
    records = list(records)
    if not records:
        raise ValueError("no records were given to featurise")
    vectors: list[FeatureVector] = []
    counts: dict[str, int] = {}
    examples: dict[str, str] = {}
    with_structure = parsed = 0

    def note(reason: str, subject: str) -> None:
        counts[reason] = counts.get(reason, 0) + 1
        examples.setdefault(reason, subject)

    for record in records:
        glycan = getattr(record, "glycan", None)
        text = getattr(glycan, "iupac_condensed", None) if glycan is not None else None
        if text is not None:
            with_structure += 1
        try:
            vector = extract(record)
        except GlycanGraphError:
            note(STRUCTURE_PARSE_FAILED, str(text))
            continue
        if text is not None:
            parsed += 1
        for reason in vector.notes:
            note(reason, str(text))
        vectors.append(vector)

    return ExtractionReport(
        vectors=tuple(vectors),
        rows_in=len(records),
        records_with_structure=with_structure,
        structures_parsed=parsed,
        vectors_built=len(vectors),
        failure_counts=counts,
        failure_examples=examples,
    )
