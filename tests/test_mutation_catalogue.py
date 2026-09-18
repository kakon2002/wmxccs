"""Every mutation must still identify exactly one site in the file it names.

This is the guard the old harness could not give itself. A full sweep takes
about three quarters of an hour and runs when somebody remembers; this runs in
under a second on every ordinary pytest invocation, so an anchor that goes
stale is caught within seconds of the edit that moved it, rather than at the
next sweep - by which time the sweep would have spent three quarters of an hour
reporting success over less than it claimed.

That is not hypothetical. It happened twice on the platform this harness was
ported from: an ordinary edit moved a line, the mutation quietly stopped being
applied, and the harness printed "SKIPPED (pattern not found)" and exited zero.
Same failure class as a CSV row swallowed by a stray quote, and as a source
count that silenced its own warning: a report of success that covers less than
it says.

A stale anchor is fixed by RE-ANCHORING it. Retiring one is a deliberate act:
delete it from the catalogue and leave a comment saying what now covers its
intent.
"""

from functools import lru_cache

import pytest

from tools.mutation import Anchor, Expect, Mutation, repeated_labels
from tools.mutation.catalogue import MUTATIONS
from tools.mutation.runner import SRC, pytest_command

# The catalogue this one replaces held 154 mutations. Forty-six of those
# anchored into modules that do not come across at all, and another twenty-six
# into modules that are not in this milestone. The floor below is what THIS
# catalogue actually holds, so that a botched edit which drops part of it fails
# here rather than passing quietly with less coverage.
CATALOGUE_FLOOR = 103


@lru_cache(maxsize=None)
def source_of(name: str) -> str:
    return (SRC / name).read_text(encoding="utf-8")


def describe(mutation: Mutation) -> str:
    return mutation.label


@pytest.mark.parametrize("mutation", MUTATIONS, ids=describe)
def test_every_anchor_identifies_exactly_one_site(mutation):
    source = source_of(mutation.file)
    anchor = mutation.anchor_in(source)
    occurrences = source.count(mutation.find)
    assert anchor is Anchor.OK, (
        f"{mutation.label}: anchor is {anchor.value} in {mutation.file}"
        f" ({occurrences} match(es)). Re-anchor it so it names exactly one site,"
        " or delete it and record what now covers its intent."
    )


@pytest.mark.parametrize("mutation", MUTATIONS, ids=describe)
def test_every_mutation_actually_changes_its_file(mutation):
    source = source_of(mutation.file)
    assert mutation.apply_to(source) != source, f"{mutation.label}: applying it leaves the file unchanged"


@pytest.mark.parametrize("mutation", MUTATIONS, ids=describe)
def test_every_target_file_is_a_package_module(mutation):
    assert (SRC / mutation.file).is_file(), f"{mutation.label}: {mutation.file} is not in the package"


def test_no_label_is_used_twice():
    # A duplicate makes a filtered run ambiguous and a report unreadable.
    assert repeated_labels(MUTATIONS) == []


def test_every_expected_survivor_says_why():
    # An exemption without a reason cannot be told apart from a gap in the tests.
    # Vacuous while the catalogue has no survivors, and deliberately kept: the
    # day one is added, this is the test that asks it to justify itself.
    for mutation in MUTATIONS:
        if mutation.expect is Expect.SURVIVES:
            assert mutation.reason.strip(), mutation.label


def test_this_catalogue_has_no_expected_survivors_at_all():
    # The catalogue this was ported from had exactly one, a documented
    # unreachable guard in splits.py, and splits.py is not ported. Nothing here
    # is exempt yet, so every mutation must be killed by something. A survivor
    # appearing without discussion is worth failing over: it would otherwise
    # read as a permission to leave one behaviour untested.
    survivors = [m.label for m in MUTATIONS if m.expect is Expect.SURVIVES]
    assert survivors == []


def test_the_catalogue_did_not_shrink_by_accident():
    # Not a target to grow towards; a FLOOR. It goes UP as modules arrive and
    # never quietly down: lowering it to match a shrunken catalogue is how a
    # loss of coverage gets ratified instead of noticed. A floor far below the
    # real count lets a great deal vanish unseen, which is the thing this file
    # exists to stop, so raise it whenever the catalogue grows.
    assert len(MUTATIONS) >= CATALOGUE_FLOOR


def test_this_file_no_longer_needs_excluding_from_the_suite_the_runner_starts():
    # It used to. While the harness mutated the repository in place, this module
    # asserted that every anchor still matched, which is false BY DESIGN mid-run
    # once a neighbouring line has moved - so it had to be ignored, or every
    # mutation would be reported killed by this file rather than by a test that
    # noticed anything. The runner now mutates a COPY and never writes the
    # repository, so the anchors checked here are never the ones under mutation
    # and the special case is gone rather than managed.
    command = pytest_command()
    assert "--ignore" not in command
    assert "-x" in command  # a killed mutation stops at the first test that notices


def test_the_runner_runs_the_whole_suite_rather_than_a_list_of_files():
    # The first harness passed a hand-written list of test files and warned in
    # its own comment that a missing one would score a mutation SURVIVED for the
    # wrong reason. That list is exactly the thing that goes stale quietly.
    command = pytest_command()
    assert not [part for part in command if part.startswith("tests/")]


# --- the vocabulary refuses a mutation that could not test anything ---------------


def test_a_mutation_that_changes_nothing_is_refused():
    with pytest.raises(ValueError, match="identical"):
        Mutation(label="x", file="models.py", find="same", replace="same")


def test_a_multi_edit_mutation_whose_every_edit_is_identical_is_refused():
    # An inert mutation is not always a single no-op edit. A pair of edits that
    # each change nothing is the same failure spelt longer, and it must be
    # refused at construction rather than discovered halfway through a sweep.
    with pytest.raises(ValueError, match="identical"):
        Mutation(label="x", file="models.py", find="same", replace="same", also=(("also", "also"),))


def test_an_empty_anchor_is_refused():
    with pytest.raises(ValueError, match="match anywhere"):
        Mutation(label="x", file="models.py", find="", replace="y")


def test_an_empty_anchor_in_a_later_edit_is_refused_too():
    # The later edits go through the same check as the first, so a mutation
    # cannot smuggle a match-anywhere anchor in behind a valid one.
    with pytest.raises(ValueError, match="match anywhere"):
        Mutation(label="x", file="models.py", find="a", replace="b", also=(("", "c"),))


def test_an_expected_survivor_without_a_reason_is_refused():
    with pytest.raises(ValueError, match="why no test can reach it"):
        Mutation(label="x", file="models.py", find="a", replace="b", expect=Expect.SURVIVES)


def test_an_expected_survivor_with_a_reason_is_accepted():
    # The negative direction of the test above: the guard must refuse a bare
    # exemption without also refusing a documented one.
    documented = Mutation(
        label="x",
        file="models.py",
        find="a",
        replace="b",
        expect=Expect.SURVIVES,
        reason="no path reaches it",
    )
    assert documented.expect is Expect.SURVIVES


def test_a_mutation_needs_a_label_and_a_file():
    with pytest.raises(ValueError, match="label"):
        Mutation(label="  ", file="models.py", find="a", replace="b")
    with pytest.raises(ValueError, match="names no file"):
        Mutation(label="x", file="", find="a", replace="b")


@pytest.mark.parametrize(
    "source, expected",
    [
        ("one hit here", Anchor.OK),
        ("no match at all", Anchor.STALE),
        ("one hit here and one hit here", Anchor.AMBIGUOUS),
    ],
)
def test_the_anchor_classes(source, expected):
    mutation = Mutation(label="x", file="models.py", find="one hit here", replace="changed")
    assert mutation.anchor_in(source) is expected


def test_an_edit_pair_that_cancels_out_is_inert_rather_than_ok():
    # Refused at construction only when every edit is textually identical. A
    # pair that undoes itself passes that check and has to be caught when it is
    # applied, or it would be swept as though it tested something.
    circular = Mutation(
        label="there and back", file="models.py", find="after", replace="later", also=(("later", "after"),)
    )
    assert circular.anchor_in("before after\n") is Anchor.INERT
