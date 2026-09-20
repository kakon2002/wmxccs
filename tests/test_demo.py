"""The end-to-end demonstration, asserted to run and to say what it claims to say.

WHY THIS FILE EXISTS. Until 20 September 2026 nothing in this repository imported,
executed or read `tools/demo_end_to_end.py`. It is the first thing a new reader is pointed
at - the README sends them there - and it was the only substantial file with no test
behind it at all. A demonstration that crashes, or that silently stops explaining itself,
is worse than none: it is read as a statement about the pipeline rather than about itself.

This is deliberately NOT a test of the pipeline. Every number the demo prints is asserted
properly elsewhere, against the objects rather than against stdout. What is asserted here
is the demo's own job: that it runs, that it does not truncate its explanations mid-word,
and that the sentences it exists to deliver are actually delivered.

The claims below are the ones a reader would be misled by their absence. `weak` is the
ordinary grade on this corpus - 216 of 417 corrected records - so a reader meeting one and
assuming breakage would be wrong, and the line preventing that reading is the deliverable.
"""

from __future__ import annotations

import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]
DEMO = REPO / "tools" / "demo_end_to_end.py"
RULE_WIDTH = 92


@pytest.fixture(scope="module")
def flowed(output) -> str:
    """`output` with every run of whitespace collapsed to one space.

    The demo WRAPS its prose, which is the fix for the mid-word truncation, so a sentence
    it prints is a sentence broken across lines. Asserting an exact substring against the
    raw text would therefore fail for the very formatting that was added, and matching a
    short fragment instead would assert almost nothing. Flowing the text first lets these
    tests name the whole sentence.
    """
    return " ".join(output.split())


@pytest.fixture(scope="module")
def output() -> str:
    finished = subprocess.run(
        [sys.executable, str(DEMO)],
        capture_output=True,
        text=True,
        encoding="utf-8",
        timeout=300,
    )
    assert finished.returncode == 0, (
        f"the demonstration the README points new readers at exits {finished.returncode}\n"
        f"{finished.stderr[-2000:]}"
    )
    return finished.stdout


def test_the_demonstration_runs_to_the_end(output):
    assert "1. INGESTION" in output
    assert "7. SCOPE" in output
    assert "model version" in output


def test_the_grade_reasons_are_not_truncated_mid_word(output, flowed):
    """They were cut at 96 characters with no ellipsis.

    The first reason ended at "though within 10% of i" and the second at "without a cave".
    A reader cannot tell a sentence that was cut from one that was never finished, so the
    tool reads as broken at exactly the point it is explaining itself.
    """
    assert "The correction is being extrapolated" in flowed
    assert "worth quoting without a caveat" in flowed
    assert "though within 10% of i\n" not in output
    assert "without a cave\n" not in output


def test_the_answer_and_scope_screens_stay_inside_the_rule_they_draw(output):
    """The demo rules its own output at 92 columns; sections 6 and 7 are prose and must fit.

    The earlier sections print identifiers and measurement tables that are wider than the
    rule and are left alone - breaking a matched-ion key across lines would make it less
    readable, not more. Prose has no such excuse.
    """
    section = output.split("6. THE ANSWER")[1]
    too_wide = [line for line in section.splitlines() if len(line) > RULE_WIDTH]
    assert not too_wide, "lines past the rule on the answer and scope screens:\n" + "\n".join(
        f"{len(line)}: {line}" for line in too_wide
    )


def test_the_demonstration_shows_a_weak_grade_and_says_why_that_is_expected(output, flowed):
    """216 of 417 corrected records grade weak, so this is the grade a reader will meet.

    Someone who sees `weak` and concludes something failed has been misled by us, and the
    two structural reasons are the whole answer: every applied stratum is below the size
    the scheme asks for, and a second rule then takes it down one further notch.
    """
    assert "grade               weak" in output
    assert "why weak, and why it is the ordinary answer here rather than a failure" in flowed
    assert "fires on EVERY corrected record in this corpus" in flowed
    assert "`supported` is out of reach by construction" in flowed
    assert "not a bad measurement, and not a broken pipeline" in flowed


def test_the_demonstration_explains_the_tail_ratio_rather_than_printing_it_bare(output, flowed):
    assert "tail ratio" in output
    assert "largest leave-one-out residual over the one that sets the interval width" in flowed


def test_the_demonstration_says_why_a_correction_smaller_than_its_interval_is_worth_applying(
    flowed,
):
    """The obvious objection to this whole platform, answered on the screen that provokes it.

    The shift is under half a per cent and the interval is four per cent wide. A reader who
    concludes the correction is therefore noise has drawn the wrong inference from a true
    observation, and the answer is that the two measure different things.
    """
    assert "the correction is smaller than the uncertainty around it" in flowed
    assert "SYSTEMATIC offset between two platforms" in flowed
    assert "A bias you know about does not average out over many measurements" in flowed
    assert "Scatter does." in flowed


def test_the_demonstration_does_not_invent_precision_the_interval_cannot_support(output, flowed):
    """It printed three decimals for every value while the API serves one or two.

    A demonstration of the service that shows numbers the service would not return is
    describing something other than the service.
    """
    assert "The model holds more digits than these" in flowed
    assert "it is what the API serves" in flowed
    assert "The ORIGINAL above is never rounded" in flowed
