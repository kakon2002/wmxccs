"""The wall between the two packages, enforced rather than described.

One repository, two packages. `wmxccs` is the CCS core; `wmxglycan` is the glycan
and isomer layer ported beside it. The rule is that they talk through an explicit
schema and never by reaching into each other's internals, and the only way that
rule survives contact with a deadline is if breaking it turns the suite red.

RIGHT NOW NOTHING CROSSES AT ALL, and that is worth stating rather than implying:
there is no boundary schema yet, because there is no traffic yet to carry. What
these tests protect is the precondition for one - that neither package can start
importing the other by accident, so that when the schema arrives it is the only
road between them and not a convenience beside three shortcuts.

The consequence, which is the point of the rule: this package carries its own
IMSType, DriftGas, Polarity and CCSMeasurement, so the repository holds two
independent definitions of a measurement condition. That looks like duplication
and is not. A single shared enum would make every change to the CCS core a change
to the glycan layer's accepted inputs, silently; two definitions with a schema
between them turn the same disagreement into a validation error at the boundary,
where someone can read it.
"""

from __future__ import annotations

import ast
import re
import subprocess
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parent.parent / "src"
CCS = SRC / "wmxccs"
GLYCAN = SRC / "wmxglycan"


def _modules(package: Path) -> list[Path]:
    found = sorted(package.glob("*.py"))
    # A glob that comes back empty would make every test below vacuously green,
    # which is the failure class this project keeps finding.
    assert len(found) >= 15, f"{package.name}: found {len(found)} modules"
    return found


def _imported_names(path: Path) -> set[str]:
    """Top-level package names this module imports, by AST rather than by substring.

    A substring scan for "wmxccs" hits every docstring that mentions the other
    package by name, and this port's docstrings mention it repeatedly and on
    purpose. Only real import statements count.
    """
    names: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Import):
            names.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            # level > 0 is a relative import, which cannot leave its own package.
            if node.level == 0 and node.module:
                names.add(node.module.split(".")[0])
    return names


@pytest.mark.parametrize("path", _modules(GLYCAN), ids=lambda p: p.name)
def test_the_glycan_layer_never_imports_the_ccs_core(path):
    assert "wmxccs" not in _imported_names(path)


@pytest.mark.parametrize("path", _modules(CCS), ids=lambda p: p.name)
def test_the_ccs_core_never_imports_the_glycan_layer(path):
    # This direction matters more than the other one. The CCS core is the P0 item:
    # it is fitted, fingerprinted and served, and an import edge pointing out of it
    # would put the glycan layer's dependencies - networkx, glycowork, scikit-learn -
    # inside the thing that is deliberately pure Python.
    assert "wmxglycan" not in _imported_names(path)


@pytest.mark.parametrize("path", _modules(GLYCAN), ids=lambda p: p.name)
def test_the_glycan_layer_never_imports_from_the_repository_it_was_ported_from(path):
    # Project2 is where this code came from. It is not a dependency, it is not on
    # this repository's path, and a checkout on another machine will not have it.
    # Copying and adapting was the instruction; an import would look identical in
    # a diff and work only on the machine the port was done on.
    imported = _imported_names(path)
    assert not {name for name in imported if "project2" in name.lower()}


def test_importing_the_ccs_core_does_not_load_the_glycan_layer():
    # The AST tests above cannot see an import inside a function body reached at
    # runtime, or one performed by importlib. This can.
    code = (
        "import sys\n"
        "import wmxccs.api\n"
        "wmxccs.api.build_default_model()\n"
        "print(sorted(m for m in sys.modules if m.split('.')[0] in"
        " ('wmxglycan', 'networkx', 'glycowork', 'sklearn', 'numpy', 'pandas')))\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    # Not merely "no wmxglycan": fitting the served model must not pull in a
    # numerical stack at all. The CCS statistics are pure Python on purpose, and
    # that is a claim the pyproject makes in prose, so something has to check it.
    assert result.stdout.strip() == "[]", result.stdout


def test_importing_the_glycan_layer_does_not_load_the_ccs_core():
    code = (
        "import sys\n"
        "import wmxglycan.enumeration, wmxglycan.features, wmxglycan.splits\n"
        "import wmxglycan.evaluation, wmxglycan.training, wmxglycan.composition\n"
        "print(sorted(m for m in sys.modules if m.split('.')[0] == 'wmxccs'))\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "[]", result.stdout


def test_the_two_packages_define_their_own_measurement_vocabulary():
    # The positive form of the rule. If someone "removes the duplication" by making
    # one of these an alias for the other, this fails - which is the point, because
    # the aliased version would look tidier and would couple the CCS core's accepted
    # values to the glycan layer's without anybody deciding to.
    from wmxccs.models import IMSType as CCSIMSType
    from wmxglycan.models import IMSType as GlycanIMSType

    assert CCSIMSType is not GlycanIMSType
    assert CCSIMSType.__module__ == "wmxccs.models"
    assert GlycanIMSType.__module__ == "wmxglycan.models"


def test_the_glycan_layer_is_not_on_the_ccs_service_dependency_list():
    # wmxccs must remain installable and servable without networkx, glycowork or
    # scikit-learn. They are declared under the `glycan` extra for exactly that
    # reason, and a dependency moved into the base list would go unnoticed.
    #
    # PARSED, NOT SCANNED. The first version of this test split the file on the
    # extras header and looked for the package names as substrings, and it failed
    # on its first run - because the comment above the dependency list NAMES all
    # three while explaining that they are not dependencies. A substring scan
    # cannot tell a requirement from a sentence about a requirement, which is the
    # mistake the other tests in this file avoid by walking the AST.
    import tomllib

    config = tomllib.loads((SRC.parent / "pyproject.toml").read_text(encoding="utf-8"))
    base = config["project"]["dependencies"]
    extras = config["project"]["optional-dependencies"]
    declared = {name for requirements in extras.values() for name in requirements}

    def named(requirements):
        # "networkx>=3.2" -> "networkx"; enough for a name check, and it does not
        # need a requirements parser to be correct about it.
        return {
            re.split(r"[<>=!~\[; ]", requirement, maxsplit=1)[0].strip() for requirement in requirements
        }

    for package in ("networkx", "glycowork", "scikit-learn"):
        assert package not in named(base), f"{package} must not be a base dependency of wmxccs"
        assert package in named(declared), f"{package} should be declared under an extra"
    # And the extra is the one the prose names, so the install instruction is real.
    assert named(extras["glycan"]) == {"networkx", "glycowork", "scikit-learn"}
