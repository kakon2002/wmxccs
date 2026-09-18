"""The mutation engine's own tests.

The engine exists to notice when a guard silently stops working, so an engine
with no tests of its own would be the least defensible thing in the repository.
`sweep` takes the suite runner as a callable precisely so this can run in
milliseconds against a fake one, rather than starting pytest a hundred times.

The `src` fixture is a real minimal package called `wmxccs`, the same name as
the package this repository pip-installs editable. That is deliberate and is
the whole point of the fixture: the shadow verification runs for real here and
has to BEAT AN INSTALLED PACKAGE OF THAT NAME, rather than resolving in an
empty field. Renaming the fixture to something nothing else provides would keep
every test green and prove nothing.
"""

from pathlib import Path

import pytest

from tools.mutation import Anchor, Expect, Mutation, Verdict
from tools.mutation.runner import (
    Outcome,
    ShadowFailed,
    main,
    pytest_command,
    report,
    select,
    shadow_of,
    shadowing_worked,
    sweep,
    verified_shadow,
    where_package_resolves,
)

PACKAGE = "wmxccs"
BODY = "before\nif guard:\n    stop()\nafter\n"
GUARD = Mutation(label="the guard is removed", file="thing.py", find="if guard:", replace="if False:")


def killed(_shadow):
    return 1, ["FAILED tests/test_thing.py::test_it - AssertionError"]


def survived(_shadow):
    return 0, ["12 passed"]


def never(_shadow):
    raise AssertionError("the suite must not run for a mutation that cannot be applied")


@pytest.fixture
def src(tmp_path):
    # Named wmxccs, the name this repository installs editable, so the shadow
    # really does have to shadow something importable.
    package = tmp_path / PACKAGE
    package.mkdir()
    (package / "__init__.py").write_text("", encoding="utf-8")
    (package / "thing.py").write_text(BODY, encoding="utf-8")
    return package


def quiet(_message):
    pass


# --- the repository is never written ------------------------------------------------


def test_the_source_is_never_modified(src):
    before = (src / "thing.py").read_text(encoding="utf-8")
    sweep([GUARD], src=src, run=killed, log=quiet)
    assert (src / "thing.py").read_text(encoding="utf-8") == before


def test_the_source_is_untouched_even_when_the_suite_runner_raises(src):
    before = (src / "thing.py").read_text(encoding="utf-8")

    def explode(_shadow):
        raise RuntimeError("pytest could not start")

    with pytest.raises(RuntimeError):
        sweep([GUARD], src=src, run=explode, log=quiet)
    assert (src / "thing.py").read_text(encoding="utf-8") == before


def test_the_mutation_is_applied_to_the_copy_and_only_to_the_copy(src):
    seen = []

    def look(shadow):
        seen.append((shadow / PACKAGE / "thing.py").read_text(encoding="utf-8"))
        # The real tree, meanwhile, still says what it always said.
        assert (src / "thing.py").read_text(encoding="utf-8") == BODY
        return killed(shadow)

    sweep([GUARD], src=src, run=look, log=quiet)
    assert "if False:" in seen[0] and "if guard:" not in seen[0]


def test_the_copy_is_put_back_between_mutations(src):
    other = Mutation(label="after is removed", file="thing.py", find="after", replace="gone")
    seen = []

    def look(shadow):
        seen.append((shadow / PACKAGE / "thing.py").read_text(encoding="utf-8"))
        return killed(shadow)

    sweep([GUARD, other], src=src, run=look, log=quiet)
    # The second run must not still be carrying the first mutation.
    assert "if False:" in seen[0] and "gone" not in seen[0]
    assert "gone" in seen[1] and "if False:" not in seen[1]


def test_the_copy_is_put_back_even_when_the_suite_runner_raises_mid_sweep(src, tmp_path, monkeypatch):
    # The restore is in a finally block, and this is the only way to see that it
    # ran: the temporary directory is normally removed on the way out, so a copy
    # left holding a deliberate defect would vanish before anything could look
    # at it. The directory is kept here instead, so the copy can be read after
    # the crash. Without the finally, a mutated copy outlives the run.
    kept = tmp_path / "kept"
    kept.mkdir()

    class KeptDirectory:
        def __init__(self, *_args, **_kwargs):
            pass

        def __enter__(self):
            return str(kept)

        def __exit__(self, *_exc):
            return False

    monkeypatch.setattr("tools.mutation.runner.tempfile.TemporaryDirectory", KeptDirectory)

    seen = []

    def look_then_explode(shadow):
        seen.append((shadow / PACKAGE / "thing.py").read_text(encoding="utf-8"))
        raise RuntimeError("pytest could not start")

    with pytest.raises(RuntimeError):
        sweep([GUARD], src=src, run=look_then_explode, log=quiet)
    assert "if False:" in seen[0]  # it really was applied before the crash
    assert (kept / "shadow" / PACKAGE / "thing.py").read_text(encoding="utf-8") == BODY


# --- a shadow that does not take is refused, never reported --------------------------


def test_the_shadow_really_shadows(src, tmp_path):
    # The package name here is installed editable in this environment, so this
    # asserts the real property: PYTHONPATH wins over the install.
    shadow = verified_shadow(src, tmp_path / "into")
    resolved = where_package_resolves(shadow, PACKAGE)
    assert Path(resolved).resolve().is_relative_to((shadow / PACKAGE).resolve())


def test_a_shadow_that_cannot_be_imported_at_all_is_refused(tmp_path):
    # A directory whose name is not a Python identifier: the import cannot even
    # be spelt, so nothing resolves and the run must refuse rather than proceed
    # against unmutated code.
    package = tmp_path / "not-importable"
    package.mkdir()
    (package / "__init__.py").write_text("", encoding="utf-8")
    with pytest.raises(ShadowFailed, match="would look like a survivor"):
        verified_shadow(package, tmp_path / "into")


def test_shadowing_is_judged_by_where_the_import_landed(tmp_path):
    # The failure that matters cannot be staged honestly - PYTHONPATH winning is
    # the thing being relied on - so the judgement drawn from each answer is
    # checked directly instead.
    package = tmp_path / "shadow" / PACKAGE
    package.mkdir(parents=True)
    inside = package / "__init__.py"
    inside.write_text("", encoding="utf-8")
    installed = tmp_path / "site-packages" / PACKAGE / "__init__.py"
    installed.parent.mkdir(parents=True)
    installed.write_text("", encoding="utf-8")

    assert shadowing_worked(str(inside), package)
    # Imported fine, but found the installed package: every mutation would then
    # be measured against unmutated code and every one would look like a survivor.
    assert not shadowing_worked(str(installed), package)
    assert not shadowing_worked("", package)  # nothing resolved at all


def test_a_shadow_that_did_not_take_is_refused_rather_than_swept(src, tmp_path, monkeypatch):
    # The whole sweep must stop, not merely note it. A sweep against unmutated
    # code reports every mutation as a survivor, which reads as a collapse of
    # the test suite rather than as a broken tool.
    monkeypatch.setattr("tools.mutation.runner.where_package_resolves", lambda _shadow, _name: "")
    with pytest.raises(ShadowFailed):
        sweep([GUARD], src=src, run=never, log=quiet)


def test_the_shadow_is_a_copy_not_a_link(src, tmp_path):
    shadow = shadow_of(src, tmp_path / "into")
    (shadow / PACKAGE / "thing.py").write_text("changed in the copy\n", encoding="utf-8")
    assert (src / "thing.py").read_text(encoding="utf-8") == BODY


# --- an unapplied mutation is a failure, never a note --------------------------------


def test_a_stale_anchor_is_a_failure_and_the_suite_is_never_run(src):
    stale = Mutation(label="moved away", file="thing.py", find="if vanished:", replace="if False:")
    (outcome,) = sweep([stale], src=src, run=never, log=quiet)
    assert outcome.anchor is Anchor.STALE
    assert not outcome.ok
    assert "matches nothing" in outcome.problem and "Re-anchor" in outcome.problem


def test_an_ambiguous_anchor_is_a_failure(src):
    (src / "thing.py").write_text("if guard:\n    a()\nif guard:\n    b()\n", encoding="utf-8")
    (outcome,) = sweep([GUARD], src=src, run=killed, log=quiet)
    assert outcome.anchor is Anchor.AMBIGUOUS
    assert not outcome.ok
    assert "more than once" in outcome.problem


def test_the_old_behaviour_of_taking_the_first_match_is_gone(src):
    # The first harness used .replace(find, replace, 1) and mutated whichever
    # site came first, reporting a kill for a guard it had never touched.
    (src / "thing.py").write_text("if guard:\n    a()\nif guard:\n    b()\n", encoding="utf-8")
    before = (src / "thing.py").read_text(encoding="utf-8")
    sweep([GUARD], src=src, run=killed, log=quiet)
    assert (src / "thing.py").read_text(encoding="utf-8") == before


@pytest.mark.xfail(
    strict=True,
    reason=(
        "sweep reads every target file with no try/except while building `pristine`, so a mutation"
        " naming a module that is not there raises FileNotFoundError and takes the whole run down"
        " before a single mutation is judged. The catalogue is supposed to fail LOUDLY and"
        " READABLY on a target that went away, the same as it does on an anchor that went away:"
        " the anchor tests call that STALE and say re-anchor it or say what covers its intent."
        " A traceback says neither, and it loses the other 102 results with it."
    ),
)
def test_a_mutation_naming_a_file_that_is_not_there_is_stale_rather_than_a_crash(src):
    # A module can vanish the same way a line can: renamed, split, or dropped
    # between milestones. That is the case this harness exists to make loud.
    ghost = Mutation(label="its module went away", file="vanished.py", find="if guard:", replace="if False:")
    (outcome,) = sweep([ghost], src=src, run=never, log=quiet)
    assert outcome.anchor is Anchor.STALE
    assert not outcome.ok
    assert "matches nothing" in outcome.problem


# --- a mutation that moves something needs more than one edit -------------------------


def test_a_mutation_can_make_several_coordinated_edits(src):
    moved = Mutation(
        label="the guard is moved below the thing it guards",
        file="thing.py",
        find="if guard:\n    stop()\n",
        replace="",
        also=(("after\n", "after\nif guard:\n    stop()\n"),),
    )
    seen = []

    def look(shadow):
        seen.append((shadow / PACKAGE / "thing.py").read_text(encoding="utf-8"))
        return killed(shadow)

    (outcome,) = sweep([moved], src=src, run=look, log=quiet)
    assert outcome.anchor is Anchor.OK
    assert seen[0] == "before\nafter\nif guard:\n    stop()\n"


def test_a_later_edit_is_anchored_against_the_text_the_earlier_one_produced():
    staged = Mutation(
        label="two steps",
        file="thing.py",
        find="one",
        replace="two",
        also=(("two", "three"),),
    )
    assert staged.anchor_in("one\n") is Anchor.OK
    assert staged.apply_to("one\n") == "three\n"


def test_edits_that_cancel_out_are_inert(src):
    circular = Mutation(
        label="there and back", file="thing.py", find="after", replace="later", also=(("later", "after"),)
    )
    (outcome,) = sweep([circular], src=src, run=killed, log=quiet)
    assert outcome.anchor is Anchor.INERT
    assert "changed nothing" in outcome.problem


def test_an_inert_mutation_never_starts_the_suite(src):
    # It applied and changed nothing, so running the suite could only produce a
    # SURVIVED that means nothing at all.
    circular = Mutation(
        label="there and back", file="thing.py", find="after", replace="later", also=(("later", "after"),)
    )
    (outcome,) = sweep([circular], src=src, run=never, log=quiet)
    assert outcome.anchor is Anchor.INERT


# --- expectations, in both directions ------------------------------------------------


def test_a_killed_mutation_that_was_expected_to_be_killed_is_fine(src):
    (outcome,) = sweep([GUARD], src=src, run=killed, log=quiet)
    assert outcome.verdict is Verdict.KILLED and outcome.ok and outcome.problem is None
    assert outcome.killer.startswith("FAILED")


def test_a_survivor_that_was_expected_to_be_killed_fails_the_run(src):
    (outcome,) = sweep([GUARD], src=src, run=survived, log=quiet)
    assert outcome.verdict is Verdict.SURVIVED and not outcome.ok
    assert "NO TEST CAUGHT THIS" in outcome.problem


def documented(label="unreachable guard"):
    return Mutation(
        label=label,
        file="thing.py",
        find="if guard:",
        replace="if False:",
        expect=Expect.SURVIVES,
        reason="no path reaches it",
    )


def test_a_documented_survivor_that_survives_is_fine(src):
    (outcome,) = sweep([documented()], src=src, run=survived, log=quiet)
    assert outcome.ok and outcome.problem is None


def test_a_documented_survivor_that_gets_killed_also_fails_the_run(src):
    # The exemption has rotted: a test now covers it, so the documentation is
    # wrong and should be deleted. An exemption that cannot expire would train
    # the signal away.
    (outcome,) = sweep([documented()], src=src, run=killed, log=quiet)
    assert not outcome.ok
    assert "exemption has rotted" in outcome.problem and "no path reaches it" in outcome.problem


@pytest.mark.parametrize("code", [2, 3, 4, 5])
def test_a_suite_that_could_not_run_is_not_a_kill(src, code):
    # pytest exits 1 when a test failed and 2 or more when it could not run at
    # all. Treating every non-zero code as a kill would let a mutation that
    # breaks collection report a success nothing earned.
    (outcome,) = sweep([GUARD], src=src, run=lambda _s: (code, ["ERROR collecting"]), log=quiet)
    assert outcome.verdict is Verdict.UNRUNNABLE
    assert not outcome.ok
    assert "could not run" in outcome.problem and f"exit {code}" in outcome.problem
    assert outcome.killer == ""  # nothing noticed it, so nothing is named


def test_only_exit_one_counts_as_a_kill(src):
    (one,) = sweep([GUARD], src=src, run=killed, log=quiet)
    (zero,) = sweep([GUARD], src=src, run=survived, log=quiet)
    assert one.verdict is Verdict.KILLED and one.exit_code == 1
    assert zero.verdict is Verdict.SURVIVED and zero.exit_code == 0


def test_a_documented_survivor_is_not_excused_from_being_unrunnable(src):
    # An exemption says "no test can reach this", not "any outcome is fine". A
    # mutation that stops the suite starting has tested nothing either way.
    (outcome,) = sweep([documented()], src=src, run=lambda _s: (3, ["INTERNALERROR"]), log=quiet)
    assert outcome.verdict is Verdict.UNRUNNABLE
    assert not outcome.ok


# --- the suite the runner starts -------------------------------------------------------


def test_the_runner_runs_the_whole_suite_and_excludes_nothing():
    # The first harness passed a hand-written list of test files and warned in
    # its own comment that a missing one would score a mutation SURVIVED for the
    # wrong reason. The second had to exclude the catalogue's anchor tests,
    # because mutating the repository in place made their assertions false
    # mid-run. Mutating a copy removes both problems rather than managing them.
    command = pytest_command()
    assert "-x" in command  # a killed mutation stops at the first test that notices
    assert "--ignore" not in command
    assert not [part for part in command if part.startswith("tests/")]


# --- the summary rounds nothing up ------------------------------------------------------


def test_the_run_fails_when_any_mutation_did_not_do_what_was_expected():
    outcomes = [
        Outcome(mutation=GUARD, anchor=Anchor.OK, verdict=Verdict.KILLED),
        Outcome(mutation=GUARD, anchor=Anchor.STALE),
    ]
    assert report(outcomes, selected=2, total=2, log=quiet) == 1


def test_the_run_passes_only_when_everything_did():
    outcomes = [Outcome(mutation=GUARD, anchor=Anchor.OK, verdict=Verdict.KILLED)]
    assert report(outcomes, selected=1, total=1, log=quiet) == 0


def test_the_summary_names_every_problem_rather_than_only_counting_them():
    # A count alone would say the run failed without saying which promise is no
    # longer kept, and the label is the whole readable part of a survivor.
    lines = []
    report(
        [Outcome(mutation=GUARD, anchor=Anchor.OK, verdict=Verdict.SURVIVED)],
        selected=1,
        total=1,
        log=lines.append,
    )
    printed = "\n".join(lines)
    assert "the guard is removed" in printed
    assert "NO TEST CAUGHT THIS" in printed


def test_a_filtered_run_says_so():
    lines = []
    report([Outcome(mutation=GUARD, anchor=Anchor.OK, verdict=Verdict.KILLED)], selected=1, total=140, log=lines.append)
    assert any("FILTERED RUN" in line and "not a full sweep" in line for line in lines)


def test_a_full_run_does_not_claim_to_be_filtered():
    lines = []
    report([Outcome(mutation=GUARD, anchor=Anchor.OK, verdict=Verdict.KILLED)], selected=1, total=1, log=lines.append)
    assert not any("FILTERED RUN" in line for line in lines)


# --- the command line refuses what it does not understand ---------------------------------


def test_an_unrecognised_option_is_refused_rather_than_dropped(capsys):
    # Dropping it would leave nothing selected, which means EVERYTHING selected,
    # which means a typo starts a three-quarter-hour sweep. Same failure as the
    # rest of this tool guards against, pointed at itself.
    assert main(["--dry-run"]) == 1
    printed = capsys.readouterr().out
    assert "unrecognised option(s): --dry-run" in printed
    assert "python -m tools.mutation" in printed  # and it says what it does understand


@pytest.mark.parametrize("flag", ["--check", "--list", "--help"])
def test_the_flags_the_tool_does_understand_are_not_refused(flag, capsys):
    # The negative direction of the test above. A refusal list that refused
    # everything would also pass that test. Only the refusal is asserted here,
    # not the exit code: --check exits one when an anchor has gone stale, and
    # that is the tool working, not the flag being rejected.
    main([flag])
    assert "unrecognised option" not in capsys.readouterr().out


def test_help_prints_the_usage_and_runs_nothing(capsys):
    assert main(["--help"]) == 0
    printed = capsys.readouterr().out
    assert "the full sweep" in printed
    assert "never written" in printed  # the property worth knowing before running it


def test_listing_the_catalogue_runs_no_suite(capsys):
    assert main(["--list"]) == 0
    printed = capsys.readouterr().out
    assert "mutations" in printed and "[models.py]" in printed


def test_a_filter_matching_no_label_is_refused(capsys):
    assert main(["no-such-mutation-anywhere"]) == 1
    assert "no mutation label contains" in capsys.readouterr().out


def test_a_duplicate_label_stops_the_run_before_anything_is_swept(capsys, monkeypatch):
    # A filtered run selects by label substring, so two mutations sharing a
    # label make "run just that one" ambiguous and the report unreadable.
    twins = (
        Mutation(label="same name", file="models.py", find="a", replace="b"),
        Mutation(label="same name", file="models.py", find="c", replace="d"),
    )
    monkeypatch.setattr("tools.mutation.runner.MUTATIONS", twins)
    assert main(["--list"]) == 1
    assert "label more than once" in capsys.readouterr().out


def test_the_check_path_reports_a_stale_anchor_without_running_the_suite(capsys, monkeypatch):
    ghost = Mutation(label="its guard moved", file="models.py", find="no such line anywhere", replace="x")
    monkeypatch.setattr("tools.mutation.runner.MUTATIONS", (ghost,))
    assert main(["--check"]) == 1
    printed = capsys.readouterr().out
    assert "STALE" in printed and "its guard moved" in printed


@pytest.mark.xfail(
    strict=True,
    reason=(
        "the --check path reads (SRC / mutation.file) with no try/except, so a mutation naming a"
        " module that is not there raises FileNotFoundError instead of reporting STALE. --check is"
        " the one-second guard that is supposed to be the readable answer to 'has anything stopped"
        " applying', and a traceback from it names neither the mutation nor what to do, and hides"
        " every anchor after it."
    ),
)
def test_the_check_path_reports_stale_rather_than_crashing_when_a_target_file_is_gone(capsys, monkeypatch):
    ghost = Mutation(label="its module went away", file="vanished.py", find="a", replace="b")
    monkeypatch.setattr("tools.mutation.runner.MUTATIONS", (ghost,))
    assert main(["--check"]) == 1
    printed = capsys.readouterr().out
    assert "STALE" in printed and "its module went away" in printed


# --- selection ----------------------------------------------------------------------------


def test_selecting_nothing_selects_everything():
    a = Mutation(label="alpha", file="f.py", find="a", replace="b")
    b = Mutation(label="beta", file="f.py", find="c", replace="d")
    assert select([a, b], []) == [a, b]


def test_selection_is_by_label_substring():
    a = Mutation(label="[S] alpha", file="f.py", find="a", replace="b")
    b = Mutation(label="[R] beta", file="f.py", find="c", replace="d")
    assert select([a, b], ["[S] "]) == [a]
    assert select([a, b], ["alpha", "beta"]) == [a, b]


def test_selection_keeps_the_catalogue_order_rather_than_the_order_asked_for():
    # The report reads in catalogue order, which is section order. A filter that
    # reordered would make two runs of the same sweep hard to compare.
    a = Mutation(label="[I] alpha", file="f.py", find="a", replace="b")
    b = Mutation(label="[M] beta", file="f.py", find="c", replace="d")
    assert select([a, b], ["beta", "alpha"]) == [a, b]
