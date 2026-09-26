"""What a prediction was made with, so that a past prediction can be reconstructed.

THERE IS NO FITTED GLYCAN MODEL, so this is not a model fingerprint in the sense the CCS
core uses one. It is a fingerprint of a DETERMINISTIC PIPELINE, and that is the more
useful thing here: the enumerator and the ranker are pure functions of a composition, two
curated tables, one reference corpus and one declared policy constant. Pin those and the
output is reproducible; change any of them and every prediction made before the change is
a prediction about a different pipeline.

TWO DIGESTS, mirroring `wmxccs`'s corpus/parameters split, because they answer different
questions and a single digest would collapse them:

  parameters   everything the platform DECIDES: the curated rules, the enzyme table, the
               placement vocabulary, the 44 feature column names, the prior policy. A change
               here means the same corpus now gives a different answer.
  corpus       the reference evidence: which release, how many structures, and a digest over
               the canonical keys themselves. A change here means the same code now gives a
               different answer.

WHAT IS DELIBERATELY NOT DIGESTED: THE SOURCE TEXT. Digesting module source would move the
fingerprint on a comment, and a fingerprint that moves for reasons that cannot change an
answer trains its readers to ignore it. `wmxccs` learned this the other way round - its
fingerprint has been unmoved across five releases precisely because it digests inputs and
parameters rather than code - and the same discipline is used here. The consequence is
stated plainly rather than hidden: A LOGIC CHANGE THAT ALTERS NO TABLE, NO CORPUS AND NO
POLICY CONSTANT WILL NOT MOVE THIS FINGERPRINT. The mutation catalogue is what guards
against that, not this digest.

The version string beside the digests is what moves for a logic change, and it is the
package version. Both travel on every response.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache

from . import __version__
from .attestation import AttestationIndex, default_attestation_index
from .constraints import ConstraintKind
from .enumeration import SITES, Enumerator
from .features import FEATURE_NAMES
from .ranking import POLICY_PRIOR, PRIOR_PSEUDOCOUNT, PRIORS

_SHORT = 12


def _digest(payload: object) -> str:
    """A stable sha256 over a JSON rendering with sorted keys.

    `sort_keys` and `separators` are load-bearing: without them the digest would depend on
    dict insertion order and on whitespace, so it would move for reasons that change nothing.
    """
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class DataSnapshot:
    """The reference evidence a prediction was made against, identified well enough to refetch.

    `keys_digest` is over the sorted canonical keys, so two snapshots with the same counts but
    different structures are distinguishable. Counts alone would not be: a release that swapped
    one structure for another would look identical.
    """

    dataset: str
    licence: str
    attribution: str
    structures_considered: int
    structures_indexed: int
    distinct_structures: int
    duplicate_rows: int
    compositions: int
    keys_digest: str

    @classmethod
    def of(cls, index: AttestationIndex) -> DataSnapshot:
        return cls(
            dataset=str(index.dataset_version),
            licence=index.dataset_version.licence,
            attribution=index.dataset_version.attribution,
            structures_considered=index.structures_considered,
            structures_indexed=index.structures_indexed,
            distinct_structures=index.distinct_structures,
            duplicate_rows=index.duplicate_rows,
            compositions=len(index.keys_by_composition),
            keys_digest=_digest(sorted(index.rows)),
        )

    @property
    def digest(self) -> str:
        return _digest(
            {
                "dataset": self.dataset,
                "considered": self.structures_considered,
                "indexed": self.structures_indexed,
                "distinct": self.distinct_structures,
                "compositions": self.compositions,
                "keys": self.keys_digest,
            }
        )

    @property
    def short(self) -> str:
        return self.digest[:_SHORT]


@dataclass(frozen=True)
class PipelineFingerprint:
    """The two digests, and the package version that moves when logic does."""

    version: str
    parameters: str
    corpus: str
    snapshot: DataSnapshot

    @property
    def short(self) -> str:
        """`<parameters>/<corpus>`, the form the CCS core prints, so the two read alike."""
        return f"{self.parameters[:_SHORT]}/{self.corpus[:_SHORT]}"

    def summary(self) -> str:
        return (
            f"wmxglycan {self.version}  pipeline {self.short}\n"
            f"  parameters {self.parameters}\n"
            f"  corpus     {self.corpus}\n"
            f"  data       {self.snapshot.dataset}"
            f"  [{self.snapshot.licence}, {self.snapshot.attribution}]\n"
            f"  structures {self.snapshot.structures_indexed} resolved,"
            f" {self.snapshot.distinct_structures} distinct,"
            f" {self.snapshot.compositions} compositions"
        )


def parameters_of(enumerator: Enumerator) -> dict:
    """Everything the platform decides with, in a form a digest can be taken over.

    Returned rather than digested directly so that a test can diff two of these and say WHICH
    parameter moved, instead of only that something did.
    """
    return {
        # The curated rules, by their content and not by their row order in the shipped CSV.
        # Rendered to JSON strings so the list is sortable: a dict is not.
        "constraints": sorted(
            json.dumps(
                {
                    "product": rule.product,
                    "contexts": sorted(rule.contexts),
                    "kind": rule.kind.value,
                    "enzyme": rule.enzyme,
                    "glycan_class": (
                        None if rule.glycan_class is None else rule.glycan_class.value
                    ),
                    # The rationale decides whether a rule is read as an ORDER constraint, via
                    # the word "subsequent", so a reworded table is a different pipeline. The
                    # rationale text itself is not digested - only the reading it produces.
                    "orders_only": rule.kind is ConstraintKind.FORBIDS
                    and "subsequent" in rule.rationale.lower(),
                },
                sort_keys=True,
            )
            for rule in enumerator.constraints
        ),
        # The enzyme table, by the monolinks it licenses rather than by its row count: two
        # tables with 355 rows each can license different bonds.
        "monolinks": sorted(monolink.text for monolink in enumerator.catalogue),
        "enzyme_count": enumerator.catalogue.enzyme_count,
        # The placement vocabulary, including the sites DROPPED because no enzyme makes them.
        "sites": sorted(
            f"{kind}:{site.linkage}:{site.residue}:{site.kind}"
            for kind, options in SITES.items()
            for site in options
        ),
        "sites_dropped": sorted(site.monolink for site in enumerator.dropped_sites),
        # The featuriser's columns, in order: the class key is built from a subset of these.
        "features": list(FEATURE_NAMES),
        # The ranker's declared policy. The prior is the one free parameter in the module and
        # it moves the answer more than the data does, so it is part of the identity.
        "prior_policy": POLICY_PRIOR,
        "prior_pseudocount": PRIOR_PSEUDOCOUNT,
        "priors": sorted(PRIORS),
    }


def fingerprint_of(
    enumerator: Enumerator | None = None, index: AttestationIndex | None = None
) -> PipelineFingerprint:
    """The fingerprint for a given enumerator and reference corpus."""
    enumerator = Enumerator() if enumerator is None else enumerator
    index = default_attestation_index() if index is None else index
    snapshot = DataSnapshot.of(index)
    return PipelineFingerprint(
        version=__version__,
        parameters=_digest(parameters_of(enumerator)),
        corpus=snapshot.digest,
        snapshot=snapshot,
    )


@lru_cache(maxsize=1)
def default_fingerprint() -> PipelineFingerprint:
    """The fingerprint of the installed tables and the installed corpus, built once."""
    return fingerprint_of()
