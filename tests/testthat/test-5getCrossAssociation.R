
context("getCrossAssociation")

test_that("getCrossAssociation", {
    # Bundled MultiAssayExperiment instead of a network fetch from
    # ExperimentHub. Experiment 1 (microbiota) is reduced to the 50 most
    # abundant genera so that the significance tests stay cheap; experiment 2
    # (metabolites, assay "nmr") has 38 features.
    data(HintikkaXOData, package = "mia")
    mae <- HintikkaXOData
    mae[[1]] <- agglomerateByRank(mae[[1]], rank = "Genus")
    mae[[1]] <- mae[[1]][rowSds(assay(mae[[1]], "counts")) > 0, ]
    mae[[1]] <- mae[[1]][
        order(rowMeans(assay(mae[[1]], "counts")), decreasing = TRUE)[1:50], ]

    ############################### Test input ###############################
    # One valid set of arguments; each entry of 'wrong' overrides some of them
    # with an invalid value and must produce an error.
    base <- list(
        x = mae, experiment1 = 1, experiment2 = 2,
        assay.type1 = "counts", assay.type2 = "nmr",
        method = "spearman", mode = "table", p.adj.method = "fdr",
        p.adj.threshold = 0.05, cor.threshold = NULL, sort = FALSE,
        filter.self.cor = FALSE, verbose = FALSE, show.warnings = FALSE)
    expect_s3_class(do.call(getCrossAssociation, base), "data.frame")
    wrong <- list(
        list(experiment1 = 4),
        list(altexp1 = 1),
        list(altexp1 = FALSE),
        list(altexp2 = "test"),
        list(experiment2 = TRUE),
        list(assay.type1 = "test"),
        list(method = 1),
        list(method = FALSE),
        list(mode = TRUE),
        list(p.adj.method = 1),
        list(p.adj.threshold = 2),
        list(p.adj.threshold = TRUE),
        list(cor.threshold = 2),
        list(cor.threshold = TRUE),
        list(sort = 1),
        list(filter.self.cor = 1),
        list(verbose = 1),
        list(x = mae[[1]], experiment2 = assay(mae[[2]], "nmr")),
        list(x = mae[[1]], experiment2 = NULL),
        # Only one of assay.type, col.var and dimred may be given, so the
        # assay.type is cleared to reach the col.var checks themselves
        list(assay.type1 = NULL, col.var1 = FALSE),
        list(assay.type2 = NULL, col.var2 = 1),
        list(assay.type1 = NULL, col.var1 = "test"),
        list(experiment2 = 1, assay.type2 = "counts", test.signif = TRUE,
             symmetric = "TRUE"),
        list(experiment2 = 1, assay.type2 = "counts", test.signif = TRUE,
             symmetric = 1),
        list(experiment2 = 1, assay.type2 = "counts", test.signif = TRUE,
             symmetric = NULL),
        list(experiment2 = 1, assay.type2 = "counts", test.signif = TRUE,
             symmetric = c(TRUE, TRUE))
    )
    for( override in wrong ){
        args <- base
        args[names(override)] <- override
        expect_error(do.call(getCrossAssociation, args))
    }
    ############################# Test input end #############################

    # Test that association is calculated correctly with numeric data: compare
    # against stats::cor.test and p.adjust applied to the same pairs
    cor <- getCrossAssociation(
        mae, method = "pearson", assay.type1 = "counts", assay.type2 = "nmr",
        p.adj.threshold = NULL, show.warnings = FALSE, test.signif = TRUE,
        verbose = FALSE)
    d1 <- t(assay(mae[[1]], "counts"))
    d2 <- t(assay(mae[[2]], "nmr"))
    ref <- lapply(seq_len(nrow(cor)), function(i){
        test <- cor.test(
            d1[, as.character(cor$Var1[i])], d2[, as.character(cor$Var2[i])],
            method = "pearson")
        c(cor = unname(test$estimate), pval = test$p.value)
    })
    ref <- do.call(rbind, ref)
    expect_equal(cor$cor, unname(ref[, "cor"]))
    expect_equal(cor$pval, unname(ref[, "pval"]))
    expect_equal(cor$p_adj, p.adjust(ref[, "pval"], method = "fdr"))

    # Test that association is calculated correctly with factor data
    # Create a dummy data
    assay1 <- matrix(rep(c("A", "B", "B"), 20*30/3),
                     nrow = 20, ncol = 30)
    assay2 <- matrix(rep(c("A", "B", "A", "A", "B", "B"), 20*30/6),
                     nrow = 20, ncol = 30)
    # Reference
    # ref <- c()
    # for(i in 1:20){
    #     ref <- c(ref, GoodmanKruskal::GKtau(assay1[i, ], assay2[i, ])$tauxy)
    # }
    ref <- c(0.25, 1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 1.00, 0.25,
             1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 1.00, 0.25, 1.00)
    # Calculate values for 20 feature-pairs
    result <- c()
    for(i in 1:20){
        result <- c(result, .calculate_gktau(assay1[i, ], assay2[i, ])$estimate)
    }
    # Values should be the same
    expect_equal(round(result, 4), round(ref, 4))

    mae_sub <- mae[1:10, 1:10]
    # Test that output is in correct type
    expect_true( is.data.frame(
        getCrossAssociation(mae_sub, assay.type1 = "counts", assay.type2 = "nmr", p.adj.threshold = NULL, show.warnings = FALSE, test.signif = TRUE)) )
    expect_true( is.data.frame(getCrossAssociation(mae_sub, assay.type1 = "counts", assay.type2 = "nmr", show.warnings = FALSE)) )
    # There should not be any p-values that are under 0
    expect_true( is.null(
        getCrossAssociation(mae_sub, assay.type1 = "counts", assay.type2 = "nmr", p.adj.threshold = 0, show.warnings = FALSE, test.signif = TRUE)) )
    # Test that output is in correct type
    expect_true( is.list(
        getCrossAssociation(mae_sub, assay.type1 = "counts", assay.type2 = "nmr", mode = "matrix", p.adj.threshold = NULL, show.warnings = FALSE, test.signif = TRUE)) )

    expect_true( is.matrix(getCrossAssociation(mae_sub, assay.type1 = "counts", assay.type2 = "nmr", mode = "matrix", show.warnings = FALSE)) )

    # There should not be any p-values that are under 0
    expect_true( is.null(
        getCrossAssociation(mae_sub, assay.type1 = "counts", assay.type2 = "nmr", p.adj.threshold = 0, mode = "matrix", show.warnings = FALSE, test.signif = TRUE)) )

    # When correlation between same assay is calculated, calculation is made
    # faster by not calculating duplicates
    cor <- getCrossAssociation(mae, experiment1 = 1, experiment2 = 1, assay.type1 = "counts", assay.type2 = "counts", show.warnings = FALSE, symmetric = TRUE, test.signif = TRUE)
    cor2 <- getCrossAssociation(mae, experiment1 = 1, experiment2 = 1, assay.type1 = "counts", assay.type2 = "counts", show.warnings = FALSE, test.signif = TRUE)
    # Get random variables and test that their duplicates are equal
    set.seed(374)
    for(i in 1:10 ){
        random_var1 <- sample(cor$Var1, 1)
        random_var2 <- sample(cor$Var1, 1)
        expect_equal(as.numeric(cor[cor$Var1 == random_var1 & cor$Var2 == random_var2, c("cor", "pval", "p_adj")]),
                     as.numeric(cor[cor$Var1 == random_var2 & cor$Var2 == random_var1, c("cor", "pval", "p_adj")]))
    }
    expect_equal(cor, cor2)
    # Test that paired samples work correctly
    tse1 <- mae[[1]]
    tse2 <- mae[[1]]
    # Convert assay to have random values
    set.seed(6478)
    mat <- matrix(sample(0:100, nrow(tse2)*ncol(tse2), replace = TRUE),
                  nrow = nrow(tse2), ncol = ncol(tse2))
    colnames(mat) <- colnames(tse2)
    rownames(mat) <- rownames(tse2)
    assay(tse2) <- mat
    # Calculate with paired samples
    cor_paired <- getCrossAssociation(
        tse1, experiment2 = tse2, assay.type1 = "counts",
        assay.type2 = "counts", paired = TRUE, by = 2, show.warnings = FALSE,
        test.signif = TRUE)
    # Calculate all pairs
    cor <- getCrossAssociation(
        tse1, experiment2 = tse2, assay.type1 = "counts",
        assay.type2 = "counts", by = 2, show.warnings = FALSE,
        test.signif = TRUE)
    # Take only pairs that are paired
    cor <- cor[cor$Var1 == cor$Var2, ]
    rownames(cor) <- NULL

    # Should be equal
    expect_equal(cor[, c("cor", "pval")], cor_paired[, c("cor", "pval")])

    # Test that result does not depend on names (if there are equal names).
    # 'cor2' above is the same calculation on the original names.
    tse <- mae[[1]]
    rownames(tse)[1:10] <- rep("Unknown", 10)
    cor_table <- getCrossAssociation(tse, assay.type1 = "counts", assay.type2 = "counts", show.warnings = FALSE, test.signif = TRUE)
    expect_equal(cor_table[ , 3:5], cor2[ , 3:5])
    mat <- getCrossAssociation(tse, assay.type1 = "counts", assay.type2 = "counts", mode = "matrix", show.warnings = FALSE)
    expect_true( is.matrix(mat) )
    expect_true(nrow(mat) == nrow(tse) && ncol(mat) == nrow(tse))
    mat <- getCrossAssociation(tse, assay.type1 = "counts", assay.type2 = "counts", mode = "matrix", show.warnings = FALSE, cor.threshold = 0.8, filter.self.cor = TRUE)
    expect_true(nrow(mat) < nrow(tse) && ncol(mat) < nrow(tse))


    # Test user's own function
    expect_true( is.data.frame(getCrossAssociation(tse, assay.type1 = "counts", assay.type2 = "counts", method = "canberra", mode = "table", show.warnings = TRUE, association.fun = stats::dist) ) )

    expect_true( is.matrix( getCrossAssociation(tse,  assay.type1 = "counts", assay.type2 = "counts",method = "bray", show.warnings = FALSE, mode = "matrix", association.fun = vegan::vegdist, test.signif = TRUE) ) )
    expect_error( getCrossAssociation(tse, assay.type1 = "counts", assay.type2 = "counts", method = "bray", show.warnings = FALSE, mode = "matrix", association.fun = DelayedMatrixStats::rowSums2, test.signif = TRUE) )

    # Test that output has right columns
    tab1 <- getCrossAssociation(tse, assay.type1 = "counts", assay.type2 = "counts", show.warnings = FALSE)
    expect_equal(colnames(tab1), c("Var1", "Var2", "cor"))
    expect_equal(colnames(cor_table), c("Var1", "Var2", "cor", "pval", "p_adj"))

    # Test that the table have same information with different levels
    tab1_levels1 <- levels(tab1$Var1)
    tab1_levels2 <- levels(tab1$Var2)
    tab1$Var1 <- as.character(tab1$Var1)
    tab1$Var2 <- as.character(tab1$Var2)
    tab2 <- getCrossAssociation(tse, assay.type1 = "counts", assay.type2 = "counts", show.warnings = FALSE, sort = TRUE)
    tab2_levels1 <- levels(tab2$Var1)
    tab2_levels2 <- levels(tab2$Var2)
    tab2$Var1 <- as.character(tab2$Var1)
    tab2$Var2 <- as.character(tab2$Var2)
    expect_equal(tab1, tab2)
    expect_true( !all(tab1_levels1 == tab2_levels1) )
    expect_true( !all(tab1_levels2 == tab2_levels2) )

    # Test altexps
    # (tse is already agglomerated, so splitByRanks warns about overwriting
    # the agglomeration metadata)
    altExps(tse) <- suppressWarnings(
        splitByRanks(tse, ranks = c("Phylum", "Family")))
    # Test that output has right columns
    expect_equal(getCrossAssociation(tse, tse, assay.type1 = "counts", assay.type2 = "counts", show.warnings = FALSE, altexp1 = 1, altexp2 = "Phylum"),
                 getCrossAssociation(altExps(tse)[[1]], altExp(tse, "Phylum"), assay.type1 = "counts", assay.type2 = "counts", show.warnings = FALSE))
    expect_equal(getCrossAssociation(tse, tse, assay.type1 = "counts", assay.type2 = "counts", show.warnings = FALSE, altexp1 = "Family", altexp2 = NULL),
                 getCrossAssociation(altExp(tse, "Family"), tse, assay.type1 = "counts", assay.type2 = "counts", show.warnings = FALSE))

    # Test colData_variable
    # Check that all the correct names are included
    indices <- c("shannon", "gini_simpson")
    tse <- addAlpha(tse, index = indices)
    res <- getCrossAssociation(tse, tse, assay.type1 = "counts", col.var2 = indices)
    unique_var1 <- unfactor(unique(res$Var1))
    unique_var2 <- unfactor(unique(res$Var2))
    rownames <- rownames(tse)

    expect_true( all(rownames %in% unique_var1) && all(unique_var1 %in% rownames) &&
        all(indices %in% unique_var2) && all(unique_var2 %in% indices) )
    # Check that assay.type is disabled
    res2 <- getCrossAssociation(tse, assay.type1 = "counts", col.var2 = indices)
    expect_equal(res, res2)

    colData(tse)[, "test"] <- rep("a")
    expect_error(getCrossAssociation(tse, assay.type1 = "counts", col.var2 = c("shannon", "test")))

    # Check that dimred works
    tse <- runMDS(tse, assay.type = "counts")
    # Either col.var or dimred must be specified, not both
    expect_error(getCrossAssociation(tse, col.var1 = "shannon", assay.type2 = "counts", dimred1 = "MDS"))
    expect_error(getCrossAssociation(tse, dimred1 = c("test", "test"), assay.type2 = "counts",))
    expect_error(getCrossAssociation(tse, dimred1 = TRUE, assay.type2 = "counts",))
    #
    test <- getCrossAssociation(tse, col.var1 = "shannon", dimred2 = "MDS")
    test2 <- getCrossAssociation(tse, col.var1 = "shannon", dimred2 = 1)
    expect_equal(test, test2)
    #
    tse <- transformAssay(tse, method = "relabundance")
    test <- getCrossAssociation(tse, dimred1 = "MDS", assay.type2 = "relabundance", mode = "matrix", method = "pearson")
    test2 <- cor(reducedDim(tse, "MDS"), t(assay(tse, "relabundance")))
    expect_equal(test, test2, check.attributes = FALSE)
})
