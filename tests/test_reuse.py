"""What each reuse status permits, pinned so that no status is ever unclassified.

reuse.py sorts every status into exactly one tier: trainable, inference-only, or
no permitted use at all. The tiers are frozensets, so a status added to the enum
and forgotten by both would land in the third one silently. Silently is the
problem: the safe direction is still a policy decision nobody made. The table
below IS the policy, one row per status, and the tests only read the module back
against it. A new status with no row fails the first test in the file.

Nothing here touches the licence registry or the gate. Whether the CLAIM behind
a status is backed by evidence is licensing.assert_trainable's question, and
test_licensing.py asks it. This file asks only what a status would permit if it
were true.
"""

from __future__ import annotations

import pytest

from wmxccs.reuse import (
    PLATFORM_USE_CONTEXT,
    UseContext,
    can_use,
    DEFAULT_REUSE_STATUS,
    ReuseStatus,
    as_reuse_status,
    can_redistribute,
    can_train_commercial,
    is_inference_only,
)

# status -> (may train commercially, inference-only, may be redistributed).
#
# Written out rather than derived from the module, because a table computed from
# the thing it checks agrees with it by construction. Redistribution is a
# separate question from use, which is why it is a third column and not implied
# by the first two: share-alike may be passed on under its own licence but may
# not be trained on, and internal data is trainable but confidential.
# FOUR columns now, not three, because permission stopped being a property of the
# status alone the day the CEO answered. Read each row as:
#
#   (usable IF COMMERCIAL, usable BY THIS PLATFORM, inference-only here, redistributable)
#
# The first and second columns differ for exactly the statuses the answer
# reopened, and that gap IS the decision: an academic-only source is usable now
# and would not be if this platform were commercialised. Keeping both columns in
# one table is what makes that visible to a reader rather than buried in a tier.
POLICY: dict[ReuseStatus, tuple[bool, bool, bool, bool]] = {
    ReuseStatus.OPEN_ATTRIBUTION: (True, True, False, True),
    ReuseStatus.OPEN_SHARE_ALIKE: (False, False, True, True),
    # Reopened by the academic context. Usable now, so no longer inference-only.
    ReuseStatus.NON_COMMERCIAL: (False, True, False, False),
    # The no-derivatives clause is not about commerce, so the context does not
    # touch it: still no permitted use of any kind.
    ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES: (False, False, False, False),
    # Reopened by the academic context. This is the row CCSbase sits on.
    ReuseStatus.ACADEMIC_ONLY: (False, True, False, False),
    ReuseStatus.INTERNAL_PROPRIETARY: (True, True, False, False),
    ReuseStatus.SYNTHETIC_FIXTURE: (True, True, False, False),
    ReuseStatus.UNVERIFIED: (False, False, False, False),
    ReuseStatus.EXCLUDED: (False, False, False, False),
}

# The statuses with no permitted use of any kind. Held separately from the table
# so that the claim "no usable tier" is asserted in its own right and not only as
# two False columns that a reader has to notice are both False.
NO_PERMITTED_USE = (
    ReuseStatus.UNVERIFIED,
    ReuseStatus.EXCLUDED,
    ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES,
)

NOT_A_STATUS_STRING = (None, 1, 1.0, True, b"unverified", ["unverified"], ReuseStatus, object())


def test_every_reuse_status_has_a_row_in_the_policy_table() -> None:
    """The guard against a status being added to the enum and classified by nobody.

    Both directions: a member with no row is unclassified, and a row naming a
    member that no longer exists is a policy about nothing.
    """
    assert set(POLICY) == set(ReuseStatus), (
        "every ReuseStatus needs a row in POLICY saying what it permits."
        f" Missing: {sorted(str(s) for s in set(ReuseStatus) - set(POLICY))}."
        f" Stale: {sorted(str(s) for s in set(POLICY) - set(ReuseStatus))}"
    )


@pytest.mark.parametrize(("status", "policy"), sorted(POLICY.items(), key=lambda item: str(item[0])))
def test_each_status_permits_exactly_what_the_policy_table_says(
    status: ReuseStatus, policy: tuple[bool, bool, bool, bool]
) -> None:
    commercial, usable_now, inference_only, redistributable = policy
    assert can_train_commercial(status) is commercial
    assert can_use(status) is usable_now
    assert is_inference_only(status) is inference_only
    assert can_redistribute(status) is redistributable


@pytest.mark.parametrize(("status", "policy"), sorted(POLICY.items(), key=lambda item: str(item[0])))
def test_the_commercial_answer_does_not_move_when_the_platform_is_academic(
    status: ReuseStatus, policy: tuple[bool, bool, bool, bool]
) -> None:
    """can_train_commercial answers about the SOURCE, not about us.

    The whole point of keeping it beside can_use is that it does not follow the
    platform. If it did, the question "would this data still be usable if we
    commercialised" would have no way of being asked, and the answer would
    silently become yes.
    """
    commercial = policy[0]
    assert can_train_commercial(status) is commercial
    assert can_use(status, UseContext.COMMERCIAL) is commercial


@pytest.mark.parametrize("status", sorted(ReuseStatus, key=str))
def test_no_status_is_usable_commercially_without_also_being_usable_academically(
    status: ReuseStatus,
) -> None:
    """The academic tier is a superset. Anything a company may use, a university may.

    Asserted as an invariant rather than read off the table, because a future edit
    that permitted something commercially and not academically would be a
    contradiction nobody would think to look for.
    """
    if can_use(status, UseContext.COMMERCIAL):
        assert can_use(status, UseContext.ACADEMIC_RESEARCH)


def test_the_platform_context_is_academic_research_and_is_not_a_parameter() -> None:
    """Where the CEO's answer is recorded, and that it is a constant.

    A context that can be passed in is one that can be passed COMMERCIAL by a
    caller who does not know what that implies, or ACADEMIC_RESEARCH by one not
    entitled to decide it. Changing what this platform is should be an edit to
    that line with a reason beside it.
    """
    import inspect

    assert PLATFORM_USE_CONTEXT is UseContext.ACADEMIC_RESEARCH
    # can_use takes a context so a caller can ASK about another one; it defaults
    # to the platform's, which is the only one the gate ever uses.
    signature = inspect.signature(can_use)
    assert signature.parameters["context"].default is None


@pytest.mark.parametrize("status", sorted(ReuseStatus, key=str))
def test_no_status_is_both_trainable_and_inference_only(status: ReuseStatus) -> None:
    """The tiers are exclusive: inference-only means it cannot reach training."""
    assert not (can_train_commercial(status) and is_inference_only(status))


@pytest.mark.parametrize("status", NO_PERMITTED_USE)
def test_a_status_with_no_permitted_use_is_in_no_usable_tier_at_all(status: ReuseStatus) -> None:
    assert can_train_commercial(status) is False
    assert is_inference_only(status) is False
    assert can_redistribute(status) is False


@pytest.mark.parametrize("status", NO_PERMITTED_USE)
@pytest.mark.parametrize("context", sorted(UseContext, key=str))
def test_a_status_with_no_permitted_use_stays_unusable_in_every_context(
    status: ReuseStatus, context: UseContext
) -> None:
    """The context reopens academic terms. It must not reopen these.

    Unverified means nobody read the terms, excluded means somebody read them and
    said no, and no-derivatives forbids sharing adapted material whoever is doing
    the sharing. None of the three is a question about commerce, so becoming an
    academic platform answers none of them.
    """
    assert can_use(status, context) is False


def test_no_derivatives_has_no_permitted_use_though_plain_non_commercial_is_inference_only() -> None:
    """The pair that must not be folded together.

    A plain non-commercial record may be held for reference beside a prediction.
    A no-derivatives one may not, because sharing adapted material is what the
    extra clause forbids, so it has no use here at all. Collapsing the two
    members would be invisible except in this direction.
    """
    # Under COMMERCIAL use the pair separates the way it always did: plain
    # non-commercial may be held for reference, no-derivatives may not.
    assert is_inference_only(ReuseStatus.NON_COMMERCIAL, UseContext.COMMERCIAL) is True
    assert is_inference_only(ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES, UseContext.COMMERCIAL) is False
    assert can_use(ReuseStatus.NON_COMMERCIAL, UseContext.COMMERCIAL) is False
    assert can_use(ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES, UseContext.COMMERCIAL) is False
    # Under this platform's own context the pair separates harder, which is the
    # point: plain non-commercial becomes trainable and no-derivatives does not
    # move at all, because the clause that blocks it is not about commerce.
    assert can_use(ReuseStatus.NON_COMMERCIAL) is True
    assert can_use(ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES) is False
    # And neither is ever usable commercially, whatever this platform is.
    assert can_train_commercial(ReuseStatus.NON_COMMERCIAL) is False
    assert can_train_commercial(ReuseStatus.NON_COMMERCIAL_NO_DERIVATIVES) is False


def test_the_default_status_is_unverified_and_the_default_is_not_trainable() -> None:
    assert DEFAULT_REUSE_STATUS is ReuseStatus.UNVERIFIED
    assert can_train_commercial(DEFAULT_REUSE_STATUS) is False
    assert is_inference_only(DEFAULT_REUSE_STATUS) is False


@pytest.mark.parametrize("status", sorted(ReuseStatus, key=str))
def test_as_reuse_status_accepts_every_status_and_its_own_string_spelling(status: ReuseStatus) -> None:
    assert as_reuse_status(status) is status
    assert as_reuse_status(status.value) is status


def test_as_reuse_status_refuses_an_unknown_string_and_says_what_is_known() -> None:
    """An unknown string must not be read as anything, least of all as permissive.

    The message lists the statuses that do exist, because the caller met this
    with a spelling in front of them and the next thing they need is the list.
    """
    with pytest.raises(ValueError, match="unknown reuse status"):
        as_reuse_status("open_sesame")
    with pytest.raises(ValueError) as refusal:
        as_reuse_status("open_sesame")
    assert "open_sesame" in str(refusal.value)
    for status in ReuseStatus:
        assert status.value in str(refusal.value)


@pytest.mark.parametrize("value", NOT_A_STATUS_STRING)
def test_as_reuse_status_refuses_anything_that_is_not_a_status_or_a_string(value: object) -> None:
    with pytest.raises(TypeError, match="must be a ReuseStatus or its string value"):
        as_reuse_status(value)


@pytest.mark.parametrize("ask", (can_train_commercial, is_inference_only, can_redistribute))
def test_an_unrecognised_status_string_raises_rather_than_answering_no(ask) -> None:
    """Fail loud, not closed-and-quiet.

    Every tier question returns False for a status with no permitted use, so an
    unrecognised spelling that merely answered False would look identical to a
    correct refusal and the typo would never be found. It raises instead.
    """
    with pytest.raises(ValueError, match="unknown reuse status"):
        ask("open_sesame")


@pytest.mark.parametrize("ask", (can_train_commercial, is_inference_only, can_redistribute))
def test_the_tier_questions_answer_on_a_string_spelling_as_on_the_status(ask) -> None:
    for status in ReuseStatus:
        assert ask(status.value) is ask(status)
