"""Cross-platform CCS harmonization."""

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
__version__ = "0.7.0"
