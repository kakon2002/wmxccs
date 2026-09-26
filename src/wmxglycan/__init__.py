"""Wellmatix glycan and isomer layer, ported beside the CCS core.

PORTED FROM Project2 (package wmxglycan), COPIED AND ADAPTED, NEVER IMPORTED.
Project2 is not a dependency of this repository and is not on its path. Each
module here began as a byte-identical copy; every subsequent edit is recorded in
docs/GLYCAN_PORT.md against the digest of the file it came from.

THIS PACKAGE NEVER IMPORTS wmxccs, AND wmxccs NEVER IMPORTS THIS PACKAGE.
They meet through an explicit schema and nowhere else. The one-line consequence,
which is the whole reason the rule is worth keeping: this package carries its own
IMSType, DriftGas, Polarity and CCSMeasurement, so the two packages hold two
independent definitions of a measurement condition and any disagreement between
them shows up at the boundary as a validation error rather than silently inside
one of them. See docs/GLYCAN_PORT.md, "the boundary".
"""

__version__ = "0.1.0.dev0"
