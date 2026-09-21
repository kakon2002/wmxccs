"""Cross-platform CCS harmonization."""

# 0.6.1 carries two fixes found by an adversarial pass over 0.6.0 itself: a demonstration
# whose interval line contradicted itself by subtraction, and a response that dropped a field
# the object it was built from carried. The second changes the served shape - `pairings` is
# now a field of the scope report and not only a number inside its prose - which is additive
# for a consumer and is why this is a patch rather than a minor.
#
# Still the milestone number rather than a claim about stability: the model is fitted on one
# study, every response says so, and the maturity stamp is `provisional` on all of them.
#
# pyproject.toml reads this attribute, so it is the single place the version is written.
__version__ = "0.6.1"
