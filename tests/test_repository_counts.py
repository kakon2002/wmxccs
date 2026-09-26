"""The counts in README are DERIVED-CHECKED, so they cannot go stale without this failing.

WHY THIS FILE EXISTS. The test count and the mutation count were written by hand into two
separate blocks of docs/HANDOVER.md and drifted three times: 2,269 against 2,317, then 3,766
against 3,842, then 3,842 against 4,088 - each time correct in one block and wrong in the other.
Refreshing them a fourth time would have been the fourth application of a method that had already
failed three times. So the duplicate was deleted, the figure lives in exactly one file, and this
derives the real numbers and compares.

THE SPLIT IS THREE-WAY, AND THAT IS NOT PEDANTRY. "2,317 wmxccs tests" is the owner's P0 evidence
that the glycan work changed nothing in the CCS core, and it is defined by which files the tests
live in. A repository-level test file - this one - belongs to NEITHER package, and counting it as
a wmxccs test would have quietly moved that number from 2,317 to 2,326 and made the P0 claim
false. A dry run of this file against a clean clone caught exactly that, before it shipped.

THE FLOOR UNDER IT. Every assertion here would pass vacuously if its regex stopped matching, so
each first asserts that it found something to check. That is the recurring failure this
repository records in LIMITATIONS 4.5: a guard that looks tested and is not. The count of
instances lives in that one section, for the same reason these counts live in one file.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

from tools.mutation.runner import all_mutations

REPO = Path(__file__).resolve().parents[1]
README = REPO / "README.md"
HANDOVER = REPO / "docs" / "HANDOVER.md"
# This file. Named here so the three-way split is derived from one place rather than repeated.
REPOSITORY_LEVEL = "tests/test_repository_counts.py"

# The CCS suite at v0.7.0-mvp, and it must not move: it is the owner's P0 item expressed as a
# number. If the glycan work ever changes it, this fails and says so in those words.
CCS_AT_V070 = 2317


@pytest.fixture(scope="module")
def readme() -> str:
    return README.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def collected() -> dict[str, int]:
    """The real counts, from pytest itself, in subprocesses so this run is not re-entered."""

    def count(*extra: str) -> int:
        result = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", "--collect-only", *extra],
            cwd=REPO, capture_output=True, text=True, timeout=900,
        )
        assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-2000:]
        found = re.search(r"(\d+) tests collected", result.stdout)
        assert found, f"could not read a collected count from:\n{result.stdout[-2000:]}"
        return int(found.group(1))

    # IGNORE-GLOBS RATHER THAN A PATH PREFIX. `pytest tests/test_glycan_` is not a path and
    # pytest refuses it; a dry run against a clean clone caught that before it shipped.
    total = count()
    ccs = count("--ignore-glob=tests/test_glycan_*", f"--ignore={REPOSITORY_LEVEL}")
    repository = count(REPOSITORY_LEVEL)
    return {
        "total": total,
        "ccs": ccs,
        "repository": repository,
        "glycan": total - ccs - repository,
    }


def test_the_three_way_split_accounts_for_every_test(collected):
    assert collected["total"] > 4000, collected
    assert collected["glycan"] > 1700, collected
    assert collected["repository"] > 5, collected
    assert collected["ccs"] + collected["glycan"] + collected["repository"] == collected["total"]


def test_the_ccs_suite_is_still_exactly_what_it_was_at_v0_7_0(collected):
    """The owner's P0 item, as a number. The glycan work must not have changed the CCS suite."""
    assert collected["ccs"] == CCS_AT_V070, (
        f"the CCS suite is now {collected['ccs']} and was {CCS_AT_V070} at v0.7.0-mvp. Either a"
        " CCS test was added or removed, or a new test file is being counted on the wrong side of"
        " the split - check REPOSITORY_LEVEL above before changing this number."
    )


# --- the README figure is the real figure -----------------------------------------------------------


def test_the_documented_test_counts_are_the_real_ones(readme, collected):
    found = re.search(
        r"# (\d[\d,]*) tests \((\d[\d,]*) wmxccs, (\d[\d,]*) wmxglycan, (\d[\d,]*) repository-level\)",
        readme,
    )
    assert found, "README no longer states a test count in the form this test reads"
    read = [int(group.replace(",", "")) for group in found.groups()]
    documented = dict(zip(("total", "ccs", "glycan", "repository"), read))
    assert documented == collected, (
        f"README says {documented} and pytest collects {collected}. Correct the README rather"
        " than this test: the figure is derived here precisely so it cannot drift."
    )
    assert documented["ccs"] + documented["glycan"] + documented["repository"] == documented["total"]


def test_the_documented_mutation_count_is_the_real_one(readme):
    found = re.search(r"# the full sweep: (\d[\d,]*) mutations", readme)
    assert found, "README no longer states a mutation count in the form this test reads"
    assert int(found.group(1).replace(",", "")) == len(all_mutations())


def test_the_documented_anchor_count_is_the_real_one(readme):
    found = re.search(r"# (\d[\d,]*) anchors", readme)
    assert found, "README no longer states an anchor count in the form this test reads"
    assert int(found.group(1).replace(",", "")) == len(all_mutations())


def test_this_check_would_notice_a_wrong_figure(collected):
    """The floor. If the regex could not see a bad number, every test above would pass on
    anything, which is the shape LIMITATIONS 4.5 collects."""
    pattern = (
        r"# (\d[\d,]*) tests \((\d[\d,]*) wmxccs, (\d[\d,]*) wmxglycan, (\d[\d,]*) repository-level\)"
    )
    broken = (
        f"# {collected['total'] + 1} tests ({collected['ccs']} wmxccs,"
        f" {collected['glycan']} wmxglycan, {collected['repository']} repository-level)"
    )
    found = re.search(pattern, broken)
    assert found, "the pattern cannot even read a well-formed line, so it checks nothing"
    assert int(found.group(1)) != collected["total"]


# --- the figure lives in exactly one file ------------------------------------------------------------


def test_handover_states_no_count_and_has_one_command_block():
    """It carried them twice and drifted three times. One block, and no figure in it."""
    handover = HANDOVER.read_text(encoding="utf-8")
    offenders = re.findall(r"#\s*[\d,]{3,7}\s+(?:tests|mutations|anchors)\b", handover)
    assert offenders == [], (
        f"docs/HANDOVER.md states a count again: {offenders}. It points at README on purpose;"
        " see the note in that file about why the duplicate was deleted rather than refreshed."
    )
    assert handover.count(".venv/Scripts/python -m pytest -q") == 1, (
        "the duplicate command block is back; one block, or the figure drifts again"
    )
    assert 'README.md`' in handover or "README.md," in handover


def test_the_counts_appear_in_readme_and_no_other_maintained_document():
    """LIMITATIONS and CONTEXT carry HISTORICAL counts - "the suite stayed green at 1,784 tests" -
    which are statements about a past moment and correctly frozen. What must not exist is a second
    place stating a CURRENT count, which is what a command-block comment is."""
    for name in ("LIMITATIONS.md", "CONTEXT.md", "docs/GLYCAN_LIMITATIONS.md", "docs/GLYCAN_PORT.md"):
        text = (REPO / name).read_text(encoding="utf-8")
        offenders = re.findall(r"#\s*[\d,]{3,7}\s+(?:tests|mutations|anchors)\b", text)
        assert offenders == [], f"{name} states a current count in a command block: {offenders}"


# --- the version, in one place, and the tag agrees with it -------------------------------------------


def test_the_distribution_version_is_written_in_exactly_one_place():
    pyproject = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    assert 'version = { attr = "wmxccs.__version__" }' in pyproject
    assert 'dynamic = ["version"]' in pyproject
    source = (REPO / "src" / "wmxccs" / "__init__.py").read_text(encoding="utf-8")
    assignments = re.findall(r"^__version__\s*=", source, re.M)
    assert len(assignments) == 1, (
        f"{len(assignments)} version assignments in wmxccs/__init__.py; the last one silently wins"
    )


def test_both_packages_state_a_version_and_neither_is_a_dev_marker():
    import wmxccs
    import wmxglycan

    for package in (wmxccs, wmxglycan):
        assert re.fullmatch(r"\d+\.\d+\.\d+", package.__version__), package.__version__
        assert "dev" not in package.__version__, (
            f"{package.__name__} is deployed, so its version should not carry a dev marker"
        )
    assert wmxccs.__version__ != wmxglycan.__version__, (
        "two packages with two histories: if these ever coincide, say so on purpose"
    )


def test_the_glycan_version_travels_on_every_glycan_response():
    """Because the stamp is what makes a past prediction reconstructible."""
    import wmxglycan
    from wmxglycan.fingerprint import default_fingerprint

    assert default_fingerprint().version == wmxglycan.__version__
