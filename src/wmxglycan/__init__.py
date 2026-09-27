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

# 0.1.0, and the `.dev0` marker is dropped: as of 27 September 2026 this package is deployed.
# It serves six endpoints and hosts a dashboard, and a version that says "dev" on something a
# CEO is asked to open is a version that misdescribes it.
#
# A SEPARATE NUMBER FROM wmxccs ON PURPOSE. Two packages with two histories: the CCS core is at
# 0.8.1 with a model fitted on 142 matched ions, this layer is at 0.1.0 with no fitted model at
# all, and one number across both would make the second look like the first. The distribution
# version is wmxccs's, because that is what pyproject reads; this one travels on every glycan
# response inside the pipeline stamp.
__version__ = "0.1.0"
