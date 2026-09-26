"""Reuse statuses and their tiers: what a record's source terms allow.

This module sits below every other module that carries a status, so that the
licence gate in licensing.py can import the licence registry in sources.py -
which imports the dataset loader, which imports the record models, which
import this - without a cycle. licensing.py re-exports everything here; the
gate and the statuses are one interface, and callers import from there.

Each status sits in exactly one tier: trainable, inference-only, or no use at
all (UNVERIFIED, EXCLUDED). Redistribution is decided separately, and never
for a status with no permitted use. The tests pin every status to its tier, so
a new status cannot slip in unclassified.
"""

from __future__ import annotations

from enum import StrEnum


class ReuseStatus(StrEnum):
    """What a record's source terms allow us to do with it."""

    # Commercial reuse allowed with attribution, e.g. CC BY 4.0.
    OPEN_ATTRIBUTION = "open_attribution"
    # Commercial reuse allowed; copies and adaptations that are shared must
    # carry the same licence, e.g. CC BY-SA 4.0. Not trainable in V1: whether
    # a model trained on such data is an adaptation of it is unsettled.
    OPEN_SHARE_ALIKE = "open_share_alike"
    # Reuse for non-commercial purposes only, e.g. CC BY-NC 4.0.
    NON_COMMERCIAL = "non_commercial"
    # Non-commercial AND no derivatives, e.g. CC BY-NC-ND 4.0. Deliberately not
    # folded into NON_COMMERCIAL: the no-derivatives clause forbids sharing
    # adapted material, so unlike a plain non-commercial record this one has no
    # permitted use here AT ALL - not training, not redistribution, and not the
    # inference-side reference use that NON_COMMERCIAL allows. It therefore sits
    # in no usable tier, beside UNVERIFIED and EXCLUDED. A CC BY-ND licence,
    # permitting commercial use but forbidding derivatives, is a different
    # combination and would need its own member; none has been met yet.
    NON_COMMERCIAL_NO_DERIVATIVES = "non_commercial_no_derivatives"
    # Terms restrict use to academic or research purposes.
    ACADEMIC_ONLY = "academic_only"
    # Our own data: trainable, and confidential.
    INTERNAL_PROPRIETARY = "internal_proprietary"
    # A record built in code by a test, and declared to be one. Trainable, so a
    # test can build a training set, and needing no licence record, because it
    # is not a claim about any source: a synthetic record's DOI, if it carries
    # one, is invented with the rest of it. A real record never carries this
    # status. The measurement loader refuses it from any file, so it cannot
    # arrive as data; a record citing a paper whose licence is on record is
    # that paper's record and is refused as synthetic; and a fit refuses a set
    # holding one unless told, explicitly, that it is a test run. It exists so
    # that a synthetic record is distinguishable from a real one at a glance,
    # and so that no invented DOI or licence has to enter the registry to let a
    # test past the gate.
    SYNTHETIC_FIXTURE = "synthetic_fixture"
    # Nobody has checked the terms yet. The default for every record.
    UNVERIFIED = "unverified"
    # Checked, and deliberately kept out of the platform.
    EXCLUDED = "excluded"


DEFAULT_REUSE_STATUS = ReuseStatus.UNVERIFIED

_TRAINABLE = frozenset(
    {ReuseStatus.OPEN_ATTRIBUTION, ReuseStatus.INTERNAL_PROPRIETARY, ReuseStatus.SYNTHETIC_FIXTURE}
)
_INFERENCE_ONLY = frozenset(
    {ReuseStatus.OPEN_SHARE_ALIKE, ReuseStatus.NON_COMMERCIAL, ReuseStatus.ACADEMIC_ONLY}
)
_REDISTRIBUTABLE = frozenset({ReuseStatus.OPEN_ATTRIBUTION, ReuseStatus.OPEN_SHARE_ALIKE})


def as_reuse_status(status: ReuseStatus | str) -> ReuseStatus:
    """The status as a ReuseStatus; ValueError for an unknown string, TypeError for anything else."""
    if isinstance(status, ReuseStatus):
        return status
    if isinstance(status, str):
        try:
            return ReuseStatus(status)
        except ValueError:
            known = ", ".join(s.value for s in ReuseStatus)
            raise ValueError(f"unknown reuse status {status!r}; expected one of: {known}") from None
    raise TypeError(f"reuse status must be a ReuseStatus or its string value, not {type(status).__name__}")


def can_train_commercial(status: ReuseStatus | str) -> bool:
    """True if a model used commercially may be trained on the record.

    For a trainable status this is the tier only. Whether the CLAIM behind the
    status is backed by a licence record is the gate's question, in
    licensing.assert_trainable.
    """
    return as_reuse_status(status) in _TRAINABLE


def can_redistribute(status: ReuseStatus | str) -> bool:
    """True if the record itself may be passed to third parties.

    Share-alike records only under their original licence. Internal data is
    ours but confidential: releasing it is a business decision made outside
    this gate, so the gate says no. A synthetic fixture is nobody's to pass on.
    """
    return as_reuse_status(status) in _REDISTRIBUTABLE


def is_inference_only(status: ReuseStatus | str) -> bool:
    """True if the record may be used, but never trained on.

    That covers terms that forbid commercial training (non-commercial,
    academic-only) and share-alike, which V1 keeps out of training by policy.
    Such a record may be held for reference on the inference side, for example
    shown with its citation next to a prediction. Whether it may also be passed
    on is can_redistribute's question, and whether a particular inference-side
    use is allowed is a question for the source's terms; this tier only
    guarantees the record cannot reach training. UNVERIFIED and EXCLUDED
    records are not inference-only: they have no permitted use at all.
    """
    return as_reuse_status(status) in _INFERENCE_ONLY
