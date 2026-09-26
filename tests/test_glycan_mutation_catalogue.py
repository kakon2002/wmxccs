"""The glycan half of the mutation catalogue, and the guard on the split itself.

The catalogue lives in two modules: `catalogue.py` holds the 282 wmxccs entries and is closed,
`catalogue_glycan.py` holds the glycan ones. A split is a hand-maintained duality, which is the
shape of instance Twelve in LIMITATIONS 4.5 - so the concatenation is DERIVED by
`runner.all_mutations()` and this file asserts the derivation holds every entry of both. Without
that, a glycan mutation could sit in a list nothing sweeps and the run would report success over
a catalogue smaller than the one on disk.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tools.mutation import Anchor, Expect, repeated_labels
from tools.mutation.catalogue import MUTATIONS
from tools.mutation.catalogue_glycan import GLYCAN_MUTATIONS
from tools.mutation.runner import GLYCAN_SRC, PACKAGE_SOURCES, SRC, SRC_ROOT, all_mutations

# A floor, not a target. It goes up as guards arrive and never quietly down.
GLYCAN_FLOOR = 57
CCS_CATALOGUE = 282


def source_of(mutation) -> str:
    return (SRC_ROOT / mutation.package / mutation.file).read_text(encoding="utf-8")


# --- the split cannot rot -------------------------------------------------------------------------


def test_the_derived_catalogue_holds_every_entry_of_both_halves():
    combined = all_mutations()
    assert len(combined) == len(MUTATIONS) + len(GLYCAN_MUTATIONS)
    assert set(id(m) for m in MUTATIONS) <= set(id(m) for m in combined)
    assert set(id(m) for m in GLYCAN_MUTATIONS) <= set(id(m) for m in combined)


def test_the_catalogue_is_assembled_on_every_call_rather_than_at_import(monkeypatch):
    # A module-level constant would silently ignore a test that patches MUTATIONS, and four
    # existing tests in test_mutation_runner.py do exactly that. They caught this the first time
    # it was written as a constant, and this pins the lesson where the split lives.
    monkeypatch.setattr("tools.mutation.runner.GLYCAN_MUTATIONS", ())
    assert len(all_mutations()) == len(MUTATIONS)


def test_no_label_is_used_twice_across_the_two_halves():
    # A duplicate label makes a filtered run ambiguous, and across two files it is exactly the
    # kind of collision nobody looks for.
    assert repeated_labels(all_mutations()) == []


def test_the_ccs_catalogue_did_not_change_when_the_glycan_one_arrived():
    assert len(MUTATIONS) == CCS_CATALOGUE
    assert all(mutation.package == "wmxccs" for mutation in MUTATIONS)


def test_the_glycan_catalogue_did_not_shrink_by_accident():
    assert len(GLYCAN_MUTATIONS) >= GLYCAN_FLOOR


# --- every glycan entry is well formed and lands in the glycan package -----------------------------


@pytest.mark.parametrize("mutation", GLYCAN_MUTATIONS, ids=lambda m: m.label)
def test_every_glycan_mutation_targets_the_glycan_package(mutation):
    assert mutation.package == "wmxglycan"
    assert (GLYCAN_SRC / mutation.file).is_file(), f"{mutation.file} is not in the glycan package"


@pytest.mark.parametrize("mutation", GLYCAN_MUTATIONS, ids=lambda m: m.label)
def test_every_glycan_anchor_identifies_exactly_one_site(mutation):
    anchor = mutation.anchor_in(source_of(mutation))
    assert anchor is Anchor.OK, f"{mutation.label}: anchor is {anchor.value} in {mutation.target}"


@pytest.mark.parametrize("mutation", GLYCAN_MUTATIONS, ids=lambda m: m.label)
def test_every_glycan_mutation_actually_changes_its_file(mutation):
    source = source_of(mutation)
    assert mutation.apply_to(source) != source


def test_no_glycan_mutation_is_an_expected_survivor():
    # Every one of these is a guard this repository wrote, so every one must be killable. An
    # expected survivor here would mean a guard nothing reaches, which is worth arguing about
    # in a comment rather than accepting silently.
    survivors = [m.label for m in GLYCAN_MUTATIONS if m.expect is Expect.SURVIVES]
    assert survivors == []


# --- the target is derived from the package rather than hand-kept ----------------------------------


# The modules this repository wrote, as opposed to the ones it ported. DERIVED, and the
# derivation is the point: the first version of the test below hand-wrote this set in the same
# breath as citing instance Twelve, which is the failure class committed inside the test that
# names it. A module is "written here" exactly when it is absent from the port record, so the
# port record is what decides it.
PORT_RECORD = Path(__file__).resolve().parents[1] / "docs" / "GLYCAN_PORT.md"


def written_here() -> set[str]:
    ported = PORT_RECORD.read_text(encoding="utf-8")
    on_disk = {path.name for path in GLYCAN_SRC.glob("*.py")} - {"__init__.py"}
    # A ported module appears in the digest table as `| `name.py` |`.
    return {name for name in on_disk if f"| `{name}` |" not in ported}


PORTED_MODULES = 15  # the port record's digest table. A fact, not a target.


def test_the_port_record_and_the_package_agree_on_what_was_written_here():
    """The floor under the derivation, as a PARTITION rather than a list of names.

    An empty `written_here()` would make the coverage test below vacuous, which is the same
    defect one level up - exactly what happened when PACKAGE_MODULES was replaced by a glob
    without a floor. An earlier version of this test hand-wrote the three modules that existed
    at the time, which meant adding a module made THIS test fail rather than the coverage test
    it exists to protect, and the natural repair was to extend the hand-written list: instance
    Twelve again, one level up again.

    So what is asserted is the shape: every module is in exactly one half, the ported half is
    the size the port record says, and neither half is empty.
    """
    mine = written_here()
    on_disk = {path.name for path in GLYCAN_SRC.glob("*.py")} - {"__init__.py"}
    ported = on_disk - mine
    assert mine, "the derivation returned nothing, so the coverage test would be vacuous"
    assert ported, "no module reads as ported, so the port record is not being read at all"
    assert mine | ported == on_disk and not (mine & ported), "the two halves must partition"
    assert len(ported) == PORTED_MODULES, sorted(ported)


def test_every_module_the_ranker_introduced_has_at_least_one_mutation():
    mine = written_here()
    assert mine, "the derivation returned nothing, so this test would pass vacuously"
    covered = {mutation.file for mutation in GLYCAN_MUTATIONS}
    assert mine <= covered, f"uncovered: {sorted(mine - covered)}"


def test_the_ported_modules_are_deliberately_not_covered_except_where_this_repo_changed_one():
    # The owner's instruction, recorded as an assertion so it is a decision rather than a drift.
    # enumeration.py is the exception and the reason is one changed line: the ordering-note
    # context gate. If coverage spreads further into ported code, this fails and someone has to
    # say why.
    covered = {mutation.file for mutation in GLYCAN_MUTATIONS}
    ported_and_covered = covered - written_here()
    assert ported_and_covered == {"enumeration.py"}, ported_and_covered
    gate = [m for m in GLYCAN_MUTATIONS if m.file == "enumeration.py"]
    assert len(gate) == 1
    assert "ordering caveat" in gate[0].label


# --- the shadow covers both packages ----------------------------------------------------------------


def test_both_packages_are_shadowed_on_every_run():
    # Shadowing only the packages a filtered run names would mean a filtered run and a full
    # sweep measure different code, and the difference would not appear in the report.
    assert set(PACKAGE_SOURCES) == {SRC, GLYCAN_SRC}
    for src in PACKAGE_SOURCES:
        assert src.is_dir()
        assert (src / "__init__.py").is_file()


def test_a_mutation_naming_an_unshadowed_package_is_a_fatal_error_not_a_stale_anchor():
    # STALE says "the guard moved", which is a fact about the source. This is a fact about the
    # tool's configuration, and reporting it as stale would send someone to re-anchor a mutation
    # that is perfectly well anchored in a file nobody copied.
    from tools.mutation import Mutation
    from tools.mutation.runner import UnshadowedPackage, sweep

    stray = Mutation(
        label="names a package nobody shadows",
        file="thing.py",
        package="not_a_package",
        find="a",
        replace="b",
    )
    with pytest.raises(UnshadowedPackage, match="not being shadowed"):
        sweep([stray], srcs=[SRC], run=lambda _shadow: (1, []), log=lambda _m: None, verify=False)


def test_the_two_packages_both_have_a_models_module_so_the_package_field_is_load_bearing():
    # The reason Mutation carries a package at all: a mutation identified by filename alone
    # would name two files and the sweep would pick one of them by accident.
    shared = {path.name for path in SRC.glob("*.py")} & {path.name for path in GLYCAN_SRC.glob("*.py")}
    assert "models.py" in shared
    assert len(shared) >= 5, sorted(shared)


def test_a_glycan_mutation_resolves_to_a_different_file_than_a_ccs_one_of_the_same_name():
    ccs = SRC_ROOT / "wmxccs" / "models.py"
    glycan = SRC_ROOT / "wmxglycan" / "models.py"
    assert ccs.is_file() and glycan.is_file()
    assert ccs.read_text(encoding="utf-8") != glycan.read_text(encoding="utf-8")
