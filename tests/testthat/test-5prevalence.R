context("prevalence")

test_that("getPrevalence", {

    data(GlobalPatterns, package="mia")
    expect_error(getPrevalence(GlobalPatterns, detection="test"),
                 "'detection' must be a single numeric value or coercible to one")
    expect_error(getPrevalence(GlobalPatterns, include.lowest="test"),
                 "'include.lowest' must be TRUE or FALSE")
    expect_error(getPrevalence(GlobalPatterns, sort="test"),
                 "'sort' must be TRUE or FALSE")
    expect_error(getPrevalence(GlobalPatterns, as.relative="test"),
                 "'as.relative' must be TRUE or FALSE")
    expect_error(getPrevalence(GlobalPatterns, assay.type="test"),
                 "'assay.type' must be a valid name")
    # Output should be always a frequency between 0 to 1
    pr <- getPrevalence(GlobalPatterns, detection=0.1/100, as.relative=TRUE)
    expect_true(min(pr) >= 0 && max(pr) <= 1)
    pr <- getPrevalence(GlobalPatterns, detection=0.1/100, as.relative=FALSE)
    expect_true(min(pr) >= 0 && max(pr) <= 1)

    # Same prevalences should be returned for as.relative T/F in certain cases.
    pr1 <- getPrevalence(GlobalPatterns, detection=1, include.lowest=TRUE, as.relative=FALSE)
    pr2 <- getPrevalence(GlobalPatterns, detection=0/100, include.lowest=FALSE, as.relative=TRUE)
    expect_true(all(pr1 == pr2))
    pr2 <- getPrevalence(GlobalPatterns, detection=0, include.lowest=FALSE, as.relative=FALSE)
    expect_true(all(pr1 == pr2))

    # Different ways to use relative abundance should yield the same output
    pr2 <- getPrevalence(GlobalPatterns, as.relative=TRUE, assay.type = "counts")
    GlobalPatterns <- transformAssay(GlobalPatterns, method="relabundance")
    pr1 <- getPrevalence(GlobalPatterns, as.relative=FALSE, assay.type = "relabundance")
    expect_true(all(pr1 == pr2))

    # Sorting should put the top values first
    pr <- getPrevalence(GlobalPatterns, sort=TRUE, detection = 0.1/100)
    expect_equal(as.vector(which.max(pr)), 1)
    pr <- names(head(getPrevalence(GlobalPatterns, sort=TRUE,  include.lowest = TRUE), 5L))
    actual <- getTop(GlobalPatterns,
                         method="prevalence",
                         top=5,
                         assay.type="counts")
    expect_equal(pr, actual)
    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL

    pr1 <- unname(getPrevalence(GlobalPatterns, detection=0.004, as.relative=TRUE))
    # If there are no rownames, getPrevalence returns names "featurex" where
    # x denotes the index.
    pr2 <- unname(getPrevalence(gp_null, detection=0.004, as.relative=TRUE))
    expect_equal(pr1, pr2)

    pr1 <- getPrevalence(GlobalPatterns, detection=0.004, as.relative=TRUE, rank = "Family")
    pr2 <- getPrevalence(gp_null, detection=0.004, as.relative=TRUE, rank = "Family")
    expect_equal(pr1, pr2)

    # Check that na.rm works correctly
    tse <- GlobalPatterns
    # Get reference value
    ref <- getPrevalence(tse, assay.type = "counts")
    # Add NA values to matrix
    remove <- c(1, 3, 10)
    assay(tse, "counts")[remove, ] <- NA
    # There should be 3 NA values if na.rm = FALSE. Otherwise there should be 0
    expect_warning(
        res <- getPrevalence(tse, assay.type = "counts", na.rm = FALSE) )
    expect_true( sum(is.na(res)) == 3)
    expect_warning(
        res <- getPrevalence(tse, assay.type = "counts", na.rm = TRUE) )
    expect_true( sum(is.na(res)) == 0)
    # Expect that other than features with NA values are the same as in reference
    removed <- rownames(tse)[remove]
    res <- res[ !names(res) %in% removed ]
    ref <- ref[ !names(ref) %in% removed ]
    expect_equal( res[ names(ref) ], ref )

    # Now test that the number of samples where feature was detected is correct
    tse <- GlobalPatterns
    ref <- getPrevalence(tse, assay.type = "counts")
    # Add NA values to specific feature that has non-zero value
    feature <- rownames(tse)[[7]]
    assay(tse, "counts")[feature, 1] <- NA
    expect_warning(
        res <- getPrevalence(tse, assay.type = "counts", na.rm = TRUE))
    # Get the feature values and check that they have correct number of samples
    res <- res[ feature ]
    ref <- ref[ feature ]
    expect_true(res*ncol(tse) == 2)
    expect_true(ref*ncol(tse) == 3)

    #
    tse <- gp_small
    rank <- "Genus"
    # Add NA values to two features that have no genus, so that they end up in
    # the group that 'empty.rm' controls
    remove <- which(is.na(rowData(tse)[[rank]]))[c(1, 2)]
    assay(tse, "counts")[remove, ] <- NA
    # Check that agglomeration works
    tse_agg <- agglomerateByRank(tse, ignore.taxonomy = FALSE, empty.rm = FALSE, rank = rank)
    expect_warning(ref <- getPrevalence(tse_agg, na.rm = FALSE))
    expect_warning(res <- getPrevalence(tse, rank = "Genus", empty.rm = FALSE))
    expect_true( all(res == ref, na.rm = TRUE) )
    #
    tse_agg <- agglomerateByRank(
        tse, ignore.taxonomy = FALSE, empty.rm = TRUE, rank = rank)
    ref <- getPrevalence(tse_agg, na.rm = TRUE)
    res <- getPrevalence(
        tse, na.rm = TRUE, rank = "Genus", empty.rm = TRUE)
    expect_true( all(res == ref, na.rm = TRUE) )

    # Test that results of addPrevalence equals to getPrevalence
    # sort should be disabled
    tse <- gp_small
    tse2 <- addPrevalence(tse, sort = TRUE)
    prev <- getPrevalence(tse, sort = FALSE)
    expect_equal(rowData(tse2)[["prevalence"]], prev)
    tse2 <- addPrevalence(tse, rank = "Family")
    prev <- getPrevalence(tse, rank = "Family")
    expect_equal(rowData(tse2)[["prevalence"]], prev)
    expect_error(addPrevalence(tse, name = 1))
    expect_error(addPrevalence(tse, name = NULL))
    expect_error(addPrevalence(tse, name = c("asd", "test")))
})

test_that("getPrevalent", {

    data(GlobalPatterns, package="mia")
    expect_error(getPrevalent(GlobalPatterns, prevalence="test"),
                 "'prevalence' must be a single numeric value or coercible to one")
    # Results compatible with getPrevalence; the same sorting for top taxa
    # obtained in different ways
    pr1 <- getPrevalent(GlobalPatterns, detection=0.1/100, as.relative=TRUE, sort=TRUE)
    pr2 <- names(getPrevalence(GlobalPatterns, rank = "Kingdom", detection=0.1/100, as.relative=TRUE, sort=TRUE))
    expect_true(all(pr1 == pr2))

    # Retrieved taxa are the same for counts and relative abundances
    pr1 <- getPrevalent(GlobalPatterns, prevalence=0.1/100, as.relative=TRUE)
    pr2 <- getPrevalent(GlobalPatterns, prevalence=0.1/100, as.relative=FALSE)
    expect_true(all(pr1 == pr2))

    # Prevalence and detection threshold at 0 has the same impact on counts and relative abundances
    pr1 <- getPrevalent(GlobalPatterns, detection=0, prevalence=0, as.relative=TRUE)
    pr2 <- getPrevalent(GlobalPatterns, detection=0, prevalence=0, as.relative=FALSE)
    expect_true(all(pr1 == pr2))

    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL

    pr1 <- getPrevalent(GlobalPatterns, detection=0.0045, prevalence = 0.25, as.relative=TRUE)
    pr1 <- which(rownames(GlobalPatterns) %in% pr1)
    pr2 <- getPrevalent(gp_null, detection=0.0045, prevalence = 0.25, as.relative=TRUE)
    expect_equal(pr1, pr2)

    pr1 <- getPrevalent(gp_small, detection=0.004, prevalence = 0.1,
                            as.relative=TRUE, rank = "Family")
    gp_null <- gp_small
    rownames(gp_null) <- NULL
    pr2 <- getPrevalent(gp_null, detection=0.004, prevalence = 0.1,
                            as.relative=TRUE, rank = "Family")
    expect_equal(pr1, pr2)

})

test_that("getPrevalentAbundance", {
    tse <- gp_small
    tse <- addPrevalentAbundance(
        tse, rank = "Class", name = "res", prevalence = 0.35, detection = 2)
    res <- getPrevalentAbundance(
        tse, rank = "Class", prevalence = 0.35, detection = 2)
    expect_equal(tse[["res"]], res)
})

test_that("getRare", {

    data(GlobalPatterns, package="mia")
    expect_error(getRare(GlobalPatterns, prevalence="test"),
                 "'prevalence' must be a single numeric value or coercible to one")

    ############# Test that output type is correct #############
    expect_type(taxa <- getRare(GlobalPatterns,
                                    detection = 130,
                                    prevalence = 90/100), "character")

    ##### Test that getPrevalent and getRare has all the taxa ####

    # Gets rownames for all the taxa
    all_taxa <- rownames(GlobalPatterns)

    # Gets prevalent taxa
    prevalent_taxa <- getPrevalent(GlobalPatterns,
                                       detection = 0,
                                       prevalence = 90/100,
                                       include.lowest = FALSE)
    # Gets rare taxa
    rare_taxa <- getRare(GlobalPatterns,
                             detection = 0,
                             prevalence = 90/100,
                             include.lowest = FALSE)

    # Concatenates prevalent and rare taxa
    prevalent_and_rare_taxa <- c(prevalent_taxa, rare_taxa)

    # If all the elements are in another vector, vector of TRUEs
    # Opposite --> negative of FALSEs
    # If every element is FALSE, any is FALSE --> opposite --> TRUE
    expect_true( !any( !(all_taxa %in% prevalent_and_rare_taxa) ) &&
                     !any( !( prevalent_and_rare_taxa %in% all_taxa) ) )

    ##### Test that getPrevalent and getRare has all the taxa, ####
    ##### but now with detection limit and with rank #####

    # The partition into prevalent and rare taxa does not depend on the rank,
    # so two ranks are enough to check the 'rank' argument.
    for( rank in c("Phylum", "Genus") ){

        # Agglomerates data by rank
        se <- agglomerateByRank(gp_small, rank = rank)

        # Gets rownames for all the taxa
        all_taxa <- rownames(se)

        # Gets prevalent taxa
        prevalent_taxa <- getPrevalent(gp_small,
                                           prevalence = 0.05,
                                           detection = 0.1,
                                           rank = rank,
                                           include.lowest = TRUE, as.relative = TRUE)
        # Gets rare taxa
        rare_taxa <- getRare(gp_small,
                                 prevalence = 0.05,
                                 detection = 0.1,
                                 rank = rank,
                                 include.lowest = TRUE, as.relative = TRUE)

        # Concatenates prevalent and rare taxa
        prevalent_and_rare_taxa <- c(prevalent_taxa, rare_taxa)

        # If all the elements are in another vector, vector of TRUEs
        # Opposite --> negative of FALSEs
        # If every element is FALSE, any is FALSE --> opposite --> TRUE
        expect_true( !any( !(all_taxa %in% prevalent_and_rare_taxa) ) &&
                         !any( !( prevalent_and_rare_taxa %in% all_taxa) ) )

        # Expect that there are no duplicates
        expect_true(!anyDuplicated(prevalent_taxa))
        expect_true(!anyDuplicated(rare_taxa))

    }

    # Check that works also when rownames is NULL
    gp_null <- GlobalPatterns
    rownames(gp_null) <- NULL

    pr1 <- getRare(GlobalPatterns, detection=0.0045, prevalence = 0.25, as.relative=TRUE)
    pr1 <- which(rownames(GlobalPatterns) %in% pr1)
    pr2 <- getRare(gp_null, detection=0.0045, prevalence = 0.25, as.relative=TRUE)
    expect_equal(pr1, pr2)

    pr1 <- getRare(gp_small, detection=0.004, prevalence = 0.1,
                       as.relative=TRUE, rank = "Family")
    gp_null <- gp_small
    rownames(gp_null) <- NULL
    pr2 <- getRare(gp_null, detection=0.004, prevalence = 0.1,
                       as.relative=TRUE, rank = "Family")
    expect_equal(pr1, pr2)

})

# subsetByPrevalent and subsetByRare share the same interface and are checked
# with the same assertions against getPrevalent and getRare respectively.
for( fun_name in c("subsetByPrevalent", "subsetByRare") ){
test_that(fun_name, {
    subset_fun <- get(fun_name)
    get_fun <- get(if( fun_name == "subsetByPrevalent" ) "getPrevalent" else "getRare")
    tse <- gp_small
    expect_error(subset_fun(tse, prevalence="test"),
                 "'prevalence' must be a single numeric value or coercible to one")
    # Expect TSE object
    expect_equal(class(subset_fun(tse)), class(tse))

    # Results compatible with getPrevalent / getRare
    pr1 <- rownames(subset_fun(
        tse, rank = "Class", detection=0.1/100, as.relative=TRUE, sort=TRUE))
    pr2 <- get_fun(tse, rank = "Class", detection=0.1/100,
                            as.relative=TRUE, sort=TRUE)
    expect_true(all(pr1 == pr2))

    # Retrieved taxa are the same for counts and relative abundances
    pr1 <- assay(subset_fun(tse, prevalence=0.1/100, as.relative=TRUE), "counts")
    pr2 <- assay(subset_fun(tse, prevalence=0.1/100, as.relative=FALSE), "counts")
    expect_true(all(pr1 == pr2))

    # Prevalence and detection threshold at 0 has the same impact on counts and relative abundances
    pr1 <- rownames(subset_fun(tse, detection=0, prevalence=0, as.relative=TRUE))
    pr2 <- rownames(subset_fun(tse, detection=0, prevalence=0, as.relative=FALSE))
    expect_true(all(pr1 == pr2))

    # Check that works also when rownames is NULL
    gp_null <- tse
    rownames(gp_null) <- NULL

    pr1 <- subset_fun(tse, detection=12, prevalence = 0.33)
    pr1 <- unname(assay(pr1, "counts"))
    pr2 <- subset_fun(gp_null, detection=12, prevalence = 0.33)
    pr2 <- unname(assay(pr2, "counts"))
    expect_equal(pr1, pr2)

    pr1 <- subset_fun(tse, detection=5, prevalence = 0.33, rank = "Phylum")
    pr1 <- unname(assay(pr1, "counts"))
    pr2 <- subset_fun(gp_null, detection=5, prevalence = 0.33, rank = "Phylum")
    pr2 <- unname(assay(pr2, "counts"))
    expect_equal(pr1, pr2)

    # Check that tree subsetting works
    expect_error(subset_fun(tse, update.tree = 1))
    expect_error(subset_fun(tse, update.tree = NULL))
    expect_error(subset_fun(tse, update.tree = c(TRUE, FALSE)))
    tse_sub <- subset_fun(tse, prevalence = 0.4, rank = "Genus", update.tree = TRUE)
    expect_equal(length(rowTree(tse_sub)$tip.label), nrow(tse_sub))
})
}

test_that("subsetByRare + subsetByPrevalent include all the taxa", {
    tse <- gp_small
    set.seed(2358)
    d <- runif(1, 0.0001, 0.1)
    p <- runif(1, 0.0001, 0.5)
    rare <- rownames(subsetByRare(tse, detection=d, prevalence=p,
                                      as.relative=TRUE))
    prevalent <- rownames(subsetByPrevalent(tse, detection=d, prevalence=p,
                                                as.relative=TRUE))
    all_taxa <- c(rare, prevalent)
    expect_true( all(all_taxa %in% rownames(tse)) )
    expect_setequal(all_taxa, rownames(tse))
})

test_that("agglomerateByPrevalence", {

    data(GlobalPatterns, package="mia")
    expect_error(agglomerateByPrevalence(GlobalPatterns, other.name=TRUE),
                 "'other.name' must be a single character value")
    actual <- agglomerateByPrevalence(GlobalPatterns, rank = "Kingdom")
    expect_s4_class(actual,class(GlobalPatterns))
    expect_equal(dim(actual),c(2,26))

    actual <- agglomerateByPrevalence(
        GlobalPatterns,
        rank = "Phylum",
        detection = 1/100,
        prevalence = 50/100,
        as.relative = TRUE,
        other.name = "test",
        update.tree = FALSE)
    expect_s4_class(actual,class(GlobalPatterns))
    expect_equal(dim(actual),c(6,26))
    expect_equal(rowData(actual)[6,"Phylum"],"test")
    expect_false(length(rowTree(actual)$tip.label) == length(rownames(actual)))

    actual <- agglomerateByPrevalence(GlobalPatterns,
                                      rank = NULL,
                                      detection = 0.0001,
                                      prevalence = 50/100,
                                      as.relative = TRUE,
                                      other.name = "test",
                                      update.tree = TRUE)
    expect_s4_class(actual,class(GlobalPatterns))
    expect_equal(dim(actual),c(6,26))
    expect_true(all(is.na(rowData(actual)[6,])))
    expect_equal(length(rowTree(actual)$tip.label), length(rownames(actual)))
    actual <- agglomerateByPrevalence(GlobalPatterns,
                                      rank = "Phylum",
                                      detection = 1/100,
                                      prevalence = 50/100,
                                      as_relative = TRUE,
                                      other_label = "test",
                                      update.tree = TRUE)
    expect_equal(length(rowTree(actual)$tip.label), length(rownames(actual)))

    # Load data from miaTime package
    skip_if_not_installed("miaTime")
    data("SilvermanAGutData", package = "miaTime")
    se <- SilvermanAGutData

    # checking reference consensus sequence generation
    actual <- agglomerateByPrevalence(se,"Genus", update.refseq = FALSE)
    # There should be only one exact match for each sequence
    seqs_test <- as.character( referenceSeq(actual) )
    seqs_ref <- as.character( referenceSeq(se) )
    expect_true(all(vapply(
    seqs_test, function(seq) sum(seqs_ref %in% seq) == 1,
    FUN.VALUE = logical(1) )) )
    # checking reference consensus sequence generation using 'Genus:Alistipes'
    reference <- se[rowData(se)[["Genus"]] %in% "Alistipes", ]
    reference <- reference[1, ]
    expect_equal(
        as.character(referenceSeq(actual)[["Alistipes"]]),
        as.character(referenceSeq(reference)[["seq_1"]])
    )

    # Merging creates consensus sequences.
    th <- 0.1
    actual <- agglomerateByPrevalence(
      se, "Genus", update.refseq = TRUE, threshold = th)
    seqs_test <- referenceSeq(actual)
    # Get single taxon as reference. Merge those sequences and test that it
    # equals to one that is output of agglomerateByPrevalence
    seqs_ref <- referenceSeq(se)
    set.seed(9536)
    feature <- sample(na.omit(rowData(actual)[["Genus"]]), 1)
    seqs_ref <- seqs_ref[ rowData(se)[["Genus"]] %in% feature ]
    seqs_ref <- .merge_refseq(
      seqs_ref, factor(rep(feature, length(seqs_ref))),
      threshold = th)
    seqs_test <- seqs_test[ names(seqs_test) %in% feature ]
    expect_equal(seqs_test, seqs_ref)
    # The same for 'Genus:Alistipes' against DECIPHER directly
    reference <- se[rowData(se)[["Genus"]] %in% "Alistipes", ]
    reference <- ConsensusSequence(referenceSeq(reference), threshold = th)
    expect_equal(
        as.character(referenceSeq(actual)[["Alistipes"]]),
        as.character(reference)
    )
})
