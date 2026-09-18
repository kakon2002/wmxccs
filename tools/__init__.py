"""Development tooling. Not part of the wmxccs package and never shipped.

`[tool.setuptools.packages.find]` looks only under `src`, so nothing here is
built into a distribution. It is on the path during tests only, via
`pythonpath` in pyproject.
"""
