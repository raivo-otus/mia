# Shared fixtures, sourced by testthat before any test file (both under
# R CMD check and devtools::test()).
#
# gp_small is a seeded 300-feature subset of GlobalPatterns with a pruned
# 300-tip rowTree. It keeps all 26 samples (and their SampleType groups), a few
# dozen phyla and genera, NA-containing ranks and a consistent tree, which makes
# it a drop-in replacement for the full 19216-feature object in tests whose
# assertions do not depend on the exact size or content of the data.
#
# Tests must treat it as read-only: copy it to a local object before modifying.
gp_small <- local({
    data(GlobalPatterns, package = "mia", envir = environment())
    set.seed(42)
    x <- GlobalPatterns[sort(sample(nrow(GlobalPatterns), 300)), ]
    rowTree(x) <- ape::keep.tip(rowTree(x), rowLinks(x)$nodeLab)
    x
})
