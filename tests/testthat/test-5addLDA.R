context("addLDA")
test_that("addLDA", {
  skip_if_not_installed("topicmodels")
  data(GlobalPatterns, package="mia")
  # Fit the models on genus-level data instead of 19216 features. On this data
  # the fit depends on the seed, and topicmodels seeds from the clock by
  # default, so both fits below get the same fixed seed.
  tse <- agglomerateByRank(GlobalPatterns, rank = "Genus")
  tse <- addLDA(tse, control = list(seed = 123))
  expect_named(reducedDims(tse),"LDA")
  expect_true(is.matrix(reducedDim(tse,"LDA")))
  expect_equal(dim(reducedDim(tse,"LDA")),c(26,2))
  red <- reducedDim(tse,"LDA")
  expect_equal(names(attributes(red)),
               c("dim","dimnames","loadings", "model", "eval_metrics"))
  expect_equal(dim(attr(red,"loadings")),c(nrow(tse),2))
  # Check if ordination matrix returned by topicmodels::LDA is the same as
  # the addLDA one (addLDA stores the getLDA result)
  df <- as.data.frame(t(assay(tse, "counts")))
  lda_model <- topicmodels::LDA(df, 2, control = list(seed = 123))
  posteriors <- topicmodels::posterior(lda_model, df)
  scores1 <- t(as.data.frame(posteriors$topics))
  loadings <- t(as.data.frame(posteriors$terms))
  # Compare topicmodels::LDA and addLDA
  expect_equal(loadings, attr(red, "loadings"), tolerance = 10**-3)
  # ERRORs
  expect_error(
    addLDA(GlobalPatterns, k = "test", assay.type = "counts", name = "LDA")
  )
  expect_error(
    addLDA(GlobalPatterns, k = 1.5, assay.type = "counts", name = "LDA")
  )
  expect_error(
    addLDA(GlobalPatterns, k = TRUE, assay.type = "counts", name = "LDA")
  )
  expect_error(
    addLDA(GlobalPatterns, k = 2, assay.type = "test", name = "LDA")
  )
  expect_error(
    addLDA(GlobalPatterns, k = 2, assay.type = 1, name = "LDA")
  )
  expect_error(
    addLDA(GlobalPatterns, k = 2, assay.type = TRUE, name = "LDA")
  )
  expect_error(
    addLDA(GlobalPatterns, k = 2, assay.type = "counts", name = 1)
  )
  expect_error(
    addLDA(GlobalPatterns, k = 2, assay.type = "counts", name = TRUE)
  )
  # Check that perplexity is calculated correctly, from the model fitted above
  lda <- reducedDim(tse, "LDA")
  ref <- topicmodels::perplexity(attr(lda, "model"))
  test <- attr(lda, "eval_metrics")[["perplexity"]]
  expect_equal(test, ref)
  # Check coherence
  skip_if_not_installed("topicdoc")
  ref <- mean( topicdoc::topic_coherence(attr(lda, "model"), df) )
  test <- attr(lda, "eval_metrics")[["coherence"]]
  expect_equal(test, ref)
})
