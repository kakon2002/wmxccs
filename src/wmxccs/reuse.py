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


class UseContext(StrEnum):
    """What THIS platform is. Not a property of any source.

    A reuse status says what a source's TERMS allow. Whether those terms permit
    what we are actually doing depends on what we are: the same CC BY-NC dataset
    is usable by a university group and not by a company, and the licence has not
    changed between the two. Permission is therefore a function of BOTH, and
    keeping them apart is what lets the second be changed without rewriting the
    first.

    Before this existed the tiers assumed commercial use, so academic-only terms
    were simply untrainable. That was correct for a commercial platform and wrong
    the moment the platform stopped being one.
    """

    # The strictest reading, and the DEFAULT. Only terms that permit commercial
    # use clear the gate.
    COMMERCIAL = "commercial"
    # Academic and research use. Terms restricting use to academic or
    # non-commercial purposes are satisfied, and are still recorded as what they
    # are, so the distinction survives if this ever changes back.
    ACADEMIC_RESEARCH = "academic_research"


# Which statuses a model may be trained on, per context. Read as: what do these
# terms permit a platform of this kind to do.
_TRAINABLE_BY_CONTEXT: dict[UseContext, frozenset[ReuseStatus]] = {
    UseContext.COMMERCIAL: frozenset(
        {ReuseStatus.OPEN_ATTRIBUTION, ReuseStatus.INTERNAL_PROPRIETARY, ReuseStatus.SYNTHETIC_FIXTURE}
    ),
    UseContext.ACADEMIC_RESEARCH: frozenset(
        {
            ReuseStatus.OPEN_ATTRIBUTION,
            ReuseStatus.INTERNAL_PROPRIETARY,
            ReuseStatus.SYNTHETIC_FIXTURE,
            # Terms restricting use to academic or research purposes are satisfied
            # by an academic platform. That is the whole point of the context.
            ReuseStatus.ACADEMIC_ONLY,
            # A plain non-commercial licence likewise. No registered source
            # currently carries this, so it changes nothing today; it is here
            # because leaving it out would be the wrong reading of the licence,
            # and a tier that is wrong in a direction nobody exercises is still
            # wrong.
            ReuseStatus.NON_COMMERCIAL,
        }
    ),
}

# THE PLATFORM'S CONTEXT, AND WHAT IT RESTS ON.
#
# Set to academic research on 19 September 2026 on the CEO's instruction, in
# answer to the question CONTEXT.md records as the one that decides the data plan:
# "Is this platform internal research or a product?" The answer was academic and
# research use, not commercial.
#
# It is a single named constant and not a parameter, for the same reason the
# floors in readiness.py are not parameters: a context that can be passed in is a
# context that can be passed COMMERCIAL by a caller who does not know what that
# implies, or ACADEMIC_RESEARCH by one who is not entitled to decide it. Changing
# what this platform is should be an edit to this line, with a reason beside it,
# in a commit somebody signed.
#
# If the platform is ever commercialised, set this back to COMMERCIAL. Nothing
# else needs to change: every academic-only source is still recorded as
# academic_only, so the gate will refuse them again on its own, and the licence
# registry will not have to be re-read.
PLATFORM_USE_CONTEXT = UseContext.ACADEMIC_RESEARCH

# Statuses that permit SOME use short of training, in any context. Share-alike is
# here by V1 policy rather than by its licence; the two non-commercial-flavoured
# ones are here because a reference use is not a commercial one.
_INFERENCE_ONLY_EVER = frozenset(
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


def can_use(status: ReuseStatus | str, context: UseContext | None = None) -> bool:
    """True if a platform of this kind may train on a record carrying this status.

    THE GATE'S PREDICATE. `context` defaults to the platform's own,
    PLATFORM_USE_CONTEXT, which is where the CEO's answer is recorded.

    For a permitted status this is the tier only. Whether the CLAIM behind the
    status is backed by a licence record is a separate question, and the gate asks
    both: see licensing.assert_trainable.
    """
    permitted = _TRAINABLE_BY_CONTEXT[context if context is not None else PLATFORM_USE_CONTEXT]
    return as_reuse_status(status) in permitted


def can_train_commercial(status: ReuseStatus | str) -> bool:
    """True if a model used COMMERCIALLY may be trained on the record.

    Deliberately unchanged in meaning now that the platform is academic, and
    deliberately still here. It answers a question about the SOURCE's terms that
    does not depend on what this platform currently is, so it is what a report
    should use to say "this data would not be usable if we commercialised", and it
    is what keeps that distinction checkable rather than theoretical.

    The gate itself uses can_use, not this.
    """
    return as_reuse_status(status) in _TRAINABLE_BY_CONTEXT[UseContext.COMMERCIAL]


def can_redistribute(status: ReuseStatus | str) -> bool:
    """True if the record itself may be passed to third parties.

    Share-alike records only under their original licence. Internal data is
    ours but confidential: releasing it is a business decision made outside
    this gate, so the gate says no. A synthetic fixture is nobody's to pass on.
    """
    return as_reuse_status(status) in _REDISTRIBUTABLE


def is_inference_only(status: ReuseStatus | str, context: UseContext | None = None) -> bool:
    """True if the record may be used, but never trained on.

    That covers terms that forbid commercial training (non-commercial,
    academic-only) and share-alike, which V1 keeps out of training by policy.
    Such a record may be held for reference on the inference side, for example
    shown with its citation next to a prediction. Whether it may also be passed
    on is can_redistribute's question, and whether a particular inference-side
    use is allowed is a question for the source's terms; this tier only
    guarantees the record cannot reach training. UNVERIFIED and EXCLUDED
    records are not inference-only: they have no permitted use at all.

    Context-aware, because the tier is. Under academic research an academic-only
    record is trainable and therefore NOT inference-only; under commercial use it
    is. The answer changes with what the platform is, which is correct, and it is
    why the context is named in the signature rather than assumed.
    """
    coerced = as_reuse_status(status)
    if can_use(coerced, context):
        return False
    return coerced in _INFERENCE_ONLY_EVER
