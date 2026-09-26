"""Cross-platform CCS harmonization."""

# 0.8.0 because THE DISTRIBUTION GAINED A SECOND PACKAGE AND TWO SERVICES. src/wmxglycan was
# ported in beside this one on 26 September 2026, its six endpoints and a dashboard were wired on
# 27 September, and the mutation harness now covers both packages. None of that is a change to
# this package - the CCS model is untouched and its fingerprint is unmoved at
# 064eb9fba603/0d69f799f6f1 - but it is unquestionably a change to the thing a tag names.
#
# THIS LINE IS THE ONLY CHANGE TO src/wmxccs SINCE v0.7.0-mvp, and it is here because the tag and
# the package have to agree. That was settled at v0.6.2, when a tag was redirected onto a later
# commit for exactly this reason and the alternative was called "the worst way round". pyproject
# reads the distribution version from this attribute, so a v0.8.0-mvp tag with 0.7.0 written here
# would be the same defect in the other direction.
#
# The fingerprint does not read this string. It digests the corpus and the fitted parameters, so
# it is the evidence that bumping a version changed no arithmetic.
__version__ = "0.8.0"

# 0.7.0 because the CORPUS GREW: 1,437 Bush Lab MicroSource records arrived as a reference
# library, which is more than a patch even though not one of them is trainable and the fitted
# model is untouched. The model fingerprint is unmoved at 064eb9fba603/0d69f799f6f1, which is
# the evidence: the corpus this version adds cannot reach a fit.
#
# 0.6.2 was DOCUMENTATION ONLY in substance: four stale or wrong figures in the prose, none
# of them a served number and none of them touching the model. The one code change is this
# string, because pyproject reads the version from here and a tag has to match something.
# The model fingerprint is unmoved at 064eb9fba603/0d69f799f6f1.
#
# 0.6.1 carried two fixes found by an adversarial pass over 0.6.0 itself: a demonstration
# whose interval line contradicted itself by subtraction, and a response that dropped a field
# the object it was built from carried. The second changes the served shape - `pairings` is
# now a field of the scope report and not only a number inside its prose - which is additive
# for a consumer and is why this is a patch rather than a minor.
#
# Still the milestone number rather than a claim about stability: the model is fitted on one
# study, every response says so, and the maturity stamp is `provisional` on all of them.
#
# pyproject.toml reads this attribute, so it is the single place the version is written.
