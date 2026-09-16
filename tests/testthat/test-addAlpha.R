test_that("Estimate Alpha Diversity Indices with Rarefaction", {
    # Rarefaction repeats the subsampling 'niter' times per index, so the data
    # must be small. Family-level GlobalPatterns (341 features at full
    # sequencing depth) keeps the relations between rarefied and non-rarefied
    # estimates that the sanity checks below rely on; a random feature subset
    # does not.
    data(GlobalPatterns, package = "mia")
    tse <- agglomerateByRank(GlobalPatterns, rank = "Family")
    depth <- min(colSums(assay(tse, "counts")), na.rm = TRUE)
    ## Testing diversity
    # Calculate the default Shannon index with no rarefaction with 3 different
    # ways: default, niter=NULL, niter=0
    tse <- addAlpha(tse, assay.type = "counts", index = "shannon")
    tse <- addAlpha(
        tse, assay.type = "counts", index = "shannon_diversity", niter = NULL)
    tse <- addAlpha(
        tse, assay.type = "counts", index = "shannon", name = "shannon2",
        niter = 0)
    # Check that index was calculated
    expect_true(any(grepl("shannon", colnames(colData(tse)))))
    expect_true(any(grepl("shannon_diversity", colnames(colData(tse)))))
    expect_true(any(grepl("shannon2", colnames(colData(tse)))))
    # They should be equal
    expect_equal(tse$shannon, tse$shannon_diversity)
    expect_equal(tse$shannon, tse$shannon2)

    # Calculate same index with 3 rarefaction rounds
    tse <- addAlpha(
        tse, assay.type = "counts", index = "shannon", sample = depth,
        niter = 3, name = "shannon_3")
    # Check that index was calculated
    expect_true(any(grepl("shannon_3", colnames(colData(tse)))))
    # They should differ little bit
    expect_false( all(tse$shannon_diversity == tse$shannon_3) )
    # However, they should be the same with some tolerance
    expect_equal(tse$shannon_diversity, tse$shannon_3, tolerance = 1e-2)
    expect_true( cor(tse$shannon_diversity, tse$shannon_3) > 0.9 )

    ## Testing dominance
    # Calculate the default gini_dominance index with no rarefaction
    tse <- addAlpha(tse, assay.type = "counts", index = "gini_dominance")
    # Calculate same index with 3 rarefaction rounds
    tse <- addAlpha(
        tse, assay.type = "counts", index = "gini_dominance", sample = depth,
        niter = 3, name = "gini_dominance_3")
    # Check that index was calculated
    expect_true( any(grepl("gini_dominance", colnames(colData(tse)))) )
    expect_true(any(grepl("gini_dominance_3", colnames(colData(tse)))))
    # They should differ little bit
    expect_false(all(tse$gini_dominance == tse$gini_dominance_3))
    # However, they should be the same with some tolerance
    expect_equal(tse$gini_dominance, tse$gini_dominance_3, tolerance = 1e-2)

    ## Testing evenness
    # Calculate the default pielou index with no rarefaction
    tse <- addAlpha(tse, assay.type = "counts", index = "pielou")
    # Calculate same index with 3 rarefaction rounds
    tse <- addAlpha(
        tse, assay.type = "counts", index = "pielou", sample = depth,
        niter = 3, name = "pielou_3")
    # Check that index was calculated
    expect_true(any(grepl("pielou", colnames(colData(tse)))))
    expect_true(any(grepl("pielou_3", colnames(colData(tse)))))
    # They should differ little bit
    expect_false(all(tse$pielou == tse$pielou_3))
    # However, they should be the same with some tolerance
    expect_equal(tse$pielou, tse$pielou_3, tolerance = 2e-1)

    ## Testing richness
    # Calculate the default chao1 index with no rarefaction
    tse <- addAlpha(tse, assay.type = "counts", index = "chao1")
    # Calculate same index with 3 rarefaction rounds (default depth)
    tse <- addAlpha(
        tse, assay.type = "counts", index = "chao1", niter = 3,
        name = "chao1_3")
    # Check that index was calculated
    expect_true(any(grepl("chao1", colnames(colData(tse)))))
    expect_true(any(grepl("chao1_3", colnames(colData(tse)))))
    # They should differ. The difference should be same with some tolerance
    expect_false(all(tse$chao1 == tse$chao1_3))
    expect_equal(tse$chao1, tse$chao1_3, tolerance = mean(tse$chao1))
    expect_true( cor(tse$chao1, tse$chao1_3) > 0.6 )

    # test non existing index
    expect_error(addAlpha(tse, assay.type = "counts", index = "test"))

    # comparing 3 iter with 5 iters estimates
    tse <- addAlpha(
        tse, assay.type = "counts", index = "shannon", sample = depth,
        niter = 5, name = "shannon_5")
    # They should differ little bit
    expect_false(all(tse$shannon_5 == tse$shannon_3))
    # However, they should be the same with some tolerance
    expect_equal(tse$shannon_3, tse$shannon_5, tolerance = 2e-2)

    # Testing with multiple indices
    tse <- addAlpha(
        tse, assay.type = "counts",
        index = c("coverage", "absolute", "camargo", "ace"))
    # Check that indices were calculated
    expect_true(any(grepl("coverage", colnames(colData(tse)))))
    expect_true(any(grepl("absolute", colnames(colData(tse)))))
    expect_true(any(grepl("camargo", colnames(colData(tse)))))
    expect_true(any(grepl("ace", colnames(colData(tse)))))

    # Testing with multiple indices with rarefaction
    tse <- addAlpha(
        tse, assay.type = "counts", sample = depth, niter = 3,
        index = c("coverage", "absolute", "camargo", "ace"),
        name = c("coverage_3", "absolute_3", "camargo_3", "ace_3"))
    # Check that indices were calculated
    expect_true(any(grepl("coverage_3", colnames(colData(tse)))))
    expect_true(any(grepl("absolute_3", colnames(colData(tse)))))
    expect_true(any(grepl("camargo_3", colnames(colData(tse)))))
    expect_true(any(grepl("ace_3", colnames(colData(tse)))))
    # Check that values differ little bit
    expect_false(all(tse$coverage == tse$coverage_3))
    expect_false(all(tse$absolute == tse$absolute_3))
    expect_false(all(tse$camargo == tse$camargo_3))
    expect_false(all(tse$ace == tse$ace_3))
    # However, they should be the same with some tolerance. The correlation
    # thresholds are sanity checks that rarefied estimates track the
    # non-rarefied ones sample by sample; on family-level data camargo and
    # absolute correlate less strongly than on the full 19216-feature object.
    expect_equal(tse$coverage, tse$coverage_3, tolerance = 0.05)
    expect_true( cor(tse$camargo, tse$camargo_3) > 0.5 )
    expect_true( cor(tse$ace, tse$ace_3) > 0.6 )
    expect_true( cor(tse$absolute, tse$absolute_3) > 0.8 )

    # Check that we get error if 'sample' is too high and all samples were
    # dropped
    expect_error(
    tse <- addAlpha(
        tse, assay.type = "counts",
        sample = 1e10, niter = 1,
        index = "absolute", name = "absolute_fail")
    )

    # Check with random rarefaction depth and check that correct samples are
    # returned. Also user should get warning about missing values.
    sample <- sort(colSums(assay(tse, "counts")), decreasing = TRUE)
    cols <- names(sample)[1:2]
    sample <- sample[[2]]
    expect_warning(
    tse <- addAlpha(
        tse, assay.type = "counts",
        sample = sample, niter = 1,
        index = "absolute", name = "absolute_missing")
    )
    res <- tse$absolute_missing
    expect_true( all(names(res)[!is.na(res)] %in% cols) )

    # Test that results of getAlpha equals to addAlpha
    expect_error(getAlpha(tse, index = 1))
    expect_error(getAlpha(tse, index = "shannon", name = TRUE))
    expect_error(getAlpha(tse, index = "shannon", name = c("test", "test2")))
    index <- c("shannon", "observed", "ace")
    name <- c("test", "test2", "random_name")
    res <- getAlpha(tse, index = index, name = name)
    colData(tse) <- NULL
    tse <- addAlpha(tse, index = index, name = name)
    res2 <- colData(tse)
    expect_equal(res, res2)
})
