"""Licence gate for data records.

Every record carries a ReuseStatus and starts as UNVERIFIED. The gate is
default-deny: a status is trainable only if it is listed as trainable in
reuse.py, so a status nobody classified can never reach training.

A trainable status is a CLAIM, and the gate checks the claim against the
licence registry in sources.py. "open_attribution" on a record is backed by
the record's DOI and a licence record for that DOI; on a component it is
built from - in this platform, its analyte identity - by the enclosing
measurement's DOI when the component names the measurement's own source, or by
the registered dataset the component was drawn from. A claim
nothing backs is refused as UnbackedClaimError, which is a LicenceGateError.
The check lives here and not only in the file loader because a TrainingSet
can be built from records constructed directly, and a convention that every
record goes through the loader fails silently the first time it is not
followed.

Two statuses need no record. internal_proprietary without a DOI is our own
data, which has no publication to point at. synthetic_fixture is a test's
declaration that the record was built in code and is not a claim about any
source; the loader refuses it from any file, so it cannot arrive as data.

assert_trainable raises instead of returning False. A loader that filters on a
boolean drops refused rows silently, and the problem only surfaces when
someone audits row counts. An exception stops the load at the offending row.
"""

from __future__ import annotations

from .reuse import (  # noqa: F401  re-exported: the gate and the statuses are one interface
    DEFAULT_REUSE_STATUS,
    PLATFORM_USE_CONTEXT,
    ReuseStatus,
    UseContext,
    as_reuse_status,
    can_redistribute,
    can_train_commercial,
    can_use,
    is_inference_only,
)
from .sources import claim_problem, component_claim_problem


class TrainingGateError(Exception):
    """A record may not enter training.

    Deliberately not a ValueError: a loader that catches ValueError to skip
    malformed rows must not swallow gate refusals along with them.
    """


class LicenceGateError(TrainingGateError):
    """A training refusal in which at least one reuse status is at fault."""


class UnbackedClaimError(LicenceGateError):
    """A trainable reuse status that no licence record backs.

    The status itself is fine; the record's word for it is all there is. The
    remedy is to read the licence on the publisher's page and record it in
    sources.py, never to relabel the record.
    """


def assert_trainable(subject: object) -> None:
    """Raise unless `subject` may enter commercial training.

    `subject` is a ReuseStatus, its string value, or a record with a
    `reuse_status` attribute. For a record the gate checks

    - every reuse status that governs it: its own, and those of the records
      listed by its `component_records()`. A CCS value from an open paper does
      not clear an analyte identity taken from a restricted one;
    - that every trainable status is backed by the licence registry: a record
      on its own DOI, a component on the enclosing record's DOI when it names
      the enclosing record's source, or on a registered dataset;
    - anything else a record reports through `training_blockers()`, such as an
      undefined drift gas, or an adduct that does not name its charge carrier.

    A bare status has no claim to check: nothing can be trained from a status
    alone, and its tier is all there is to say about it.

    Raises LicenceGateError if any reuse status refuses, else UnbackedClaimError
    (itself a LicenceGateError) if any claim is unbacked, else TrainingGateError
    if any other blocker applies; all are TrainingGateError. A record hook that
    raises counts as a refusal, so no other exception leaves the gate in place
    of one. Written with explicit raises, not `assert`, so it still runs under
    `python -O`.
    """
    licence: list[str] = []
    unbacked: list[str] = []
    other: list[str] = []
    _collect(subject, licence, unbacked, other)
    if licence or unbacked or other:
        error = LicenceGateError if licence else UnbackedClaimError if unbacked else TrainingGateError
        raise error(
            "training gate refused a record: "
            + "; ".join(licence + unbacked + other)
            + ". Correct a record only where its source supports the correction,"
            " and do not filter it out quietly."
        )


def _collect(
    subject: object, licence: list[str], unbacked: list[str], other: list[str], enclosing: object = None
) -> None:
    """Add every refusal for `subject`, and for each record it is built from, to the three lists.

    The record hooks are called defensively. A raw ValueError or TypeError
    escaping here could be swallowed by a loader that skips malformed rows,
    taking a licence refusal with it, so a hook that raises becomes a refusal.
    """
    if subject is None or isinstance(subject, str):  # a bare status; ReuseStatus is a str
        if (refusal := _refusal("record", subject)) is not None:
            licence.append(refusal)
        return
    # Reading the record is itself defensive. getattr's default absorbs only
    # AttributeError, so an accessor raising anything else would leave the gate
    # as a raw exception - past every `except TrainingGateError` a caller has,
    # and past the loader's buckets, aborting a whole load instead of counting
    # one row.
    try:
        status = subject.reuse_status
    except AttributeError:
        raise LicenceGateError(
            f"cannot tell the reuse status of a {type(subject).__name__}; the gate accepts"
            " a ReuseStatus or a record with a reuse_status"
        ) from None
    except Exception as exc:
        raise LicenceGateError(
            f"the reuse status of a {type(subject).__name__} could not be read"
            f" ({type(exc).__name__}: {exc}), so it is refused rather than assumed"
        ) from exc
    try:
        where = getattr(subject, "source", "(no source recorded)")
    except Exception as exc:
        where = f"(source unreadable: {type(exc).__name__})"
    what = f"{type(subject).__name__} from {where!r}"
    if (refusal := _refusal(what, status)) is not None:
        licence.append(refusal)
    elif (problem := _claim_refusal(subject, enclosing)) is not None:
        unbacked.append(f"{what}: {problem}")

    blockers = getattr(subject, "training_blockers", None)
    if blockers is not None:
        try:
            found = [str(blocker) for blocker in blockers()]
        except Exception as exc:  # fail closed: an unreadable record is refused, not waved through
            found = [f"its training blockers could not be checked ({type(exc).__name__}: {exc})"]
        other.extend(f"{what}: {blocker}" for blocker in found)

    components = getattr(subject, "component_records", None)
    if components is not None:
        try:
            parts = list(components())
        except Exception as exc:  # fail closed: licences that cannot be read are not cleared
            licence.append(
                f"{what}: the records it is built from could not be listed, so their licences are unchecked"
                f" ({type(exc).__name__}: {exc})"
            )
            return
        for component in parts:
            _collect(component, licence, unbacked, other, enclosing=subject)


def _claim_refusal(subject: object, enclosing: object) -> str | None:
    """Why the subject's trainable status is not backed by the registry, or None if it is.

    A top-level record with a `doi` of its own is checked on it. A component
    that names a DOI of its own is checked on THAT, never on the enclosing
    record's: otherwise a matching source string would let one record's
    citation back a component that cites something else, which is the
    laundering the component rule exists to prevent, reached from the other
    side. A component naming no DOI - an analyte identity record, which has no
    such field - is checked as a component: backed by a registered dataset named
    in its `source`, or by the enclosing record's DOI when it names the
    enclosing record's source. Wrapped like the record hooks: an exception here becomes a
    refusal, so a claim that cannot be checked is not waved through.
    """
    try:
        status = as_reuse_status(subject.reuse_status)  # type: ignore[attr-defined]
        if enclosing is None and hasattr(subject, "doi"):
            return claim_problem(subject.doi, status)  # type: ignore[attr-defined]
        if (own_doi := getattr(subject, "doi", None)) is not None:
            return claim_problem(own_doi, status)
        return component_claim_problem(
            getattr(enclosing, "doi", None),
            getattr(subject, "source", None),
            status,
            measurement_source=getattr(enclosing, "source", None),
        )
    except Exception as exc:  # fail closed
        return f"its reuse-status claim could not be checked ({type(exc).__name__}: {exc})"


def declares_synthetic(subject: object) -> bool:
    """True if this record, or anything it is built from, declares itself invented.

    The declaration can sit on a COMPONENT: an analyte identity built by a test
    inside a measurement claiming a real status. A count taken from the
    top-level status alone would present such a set as real data, which is the
    one thing the status exists to prevent, so the whole record is asked, not
    just its outermost layer.

    Fail-closed in the direction that matters here: a record whose components
    cannot be listed counts as synthetic, because the alternative is to present
    something unreadable as real. The gate refuses such a record anyway, so in
    practice this only ever labels, never admits.
    """
    seen: set[int] = set()

    def walk(record: object) -> bool:
        if id(record) in seen:  # a cycle is not a reason to loop forever
            return False
        seen.add(id(record))
        if isinstance(record, str):  # a bare status; ReuseStatus is a str
            try:
                return as_reuse_status(record) is ReuseStatus.SYNTHETIC_FIXTURE
            except (TypeError, ValueError):
                return False
        try:
            if as_reuse_status(record.reuse_status) is ReuseStatus.SYNTHETIC_FIXTURE:  # type: ignore[attr-defined]
                return True
        except Exception:
            pass  # unreadable status: the gate's business, not this count's
        components = getattr(record, "component_records", None)
        if components is None:
            return False
        try:
            parts = list(components())
        except Exception:
            return True  # label it rather than present it as data
        return any(walk(part) for part in parts)

    return walk(subject)


# Why each untrainable status is refused. Share-alike is refused by V1 policy,
# not by its licence, and the message must not suggest otherwise: blaming the
# licence would invite relabelling the record.
_REFUSAL_REASONS = {
    ReuseStatus.OPEN_SHARE_ALIKE: (
        "which V1 keeps out of training by policy: the licence allows commercial use, but whether a model"
        " trained on the data is an adaptation of it is unsettled, so do not relabel the record"
    ),
    ReuseStatus.NON_COMMERCIAL: "whose terms do not permit commercial use",
    ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES: (
        "whose terms forbid commercial use AND derivative works; each clause blocks training on its own, so"
        " no reading of the licence reaches it and relabelling the record would be a false claim"
    ),
    ReuseStatus.ACADEMIC_ONLY: "whose terms restrict use to academic or research purposes",
    ReuseStatus.UNVERIFIED: "meaning nobody has checked its terms yet",
    ReuseStatus.EXCLUDED: "meaning it was checked and deliberately kept out of the platform",
}


def _refusal(what: str, status: object) -> str | None:
    try:
        coerced = as_reuse_status(status)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return f"{what} has unrecognised reuse status {status!r}"
    # can_use, not can_train_commercial: what these terms permit depends on what
    # this platform is, and PLATFORM_USE_CONTEXT is where that is recorded.
    if can_use(coerced):
        return None
    reason = _REFUSAL_REASONS.get(coerced, "which is not cleared for training")
    return f"{what} has reuse status '{coerced.value}', {reason}"
