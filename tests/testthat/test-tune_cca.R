test_that("prep_cv_folds creates valid structure", {
  skip_if_not_installed("testthat")

  prism_test <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
    readRDS()

  # Create mock predictor by modifying prism
  coarse <- prism_test %>%
    mutate(tmean = tmean * 0.8 + units::set_units(rnorm(length(tmean), 0, 0.5), "°C"))

  cv <- prep_cv_folds(coarse, prism_test,
                      kfolds = 3,
                      max_k_pred = 5,
                      max_k_resp = 5,
                      weight = FALSE)

  # Check class

  expect_s3_class(cv, "cv_folds")

  # Check structure
  expect_equal(cv$kfolds, 3)
  expect_equal(cv$max_k_pred, 5)
  expect_equal(cv$max_k_resp, 5)
  expect_equal(length(cv$folds), 3)

  # Check fold contents
  fold1 <- cv$folds[[1]]
  expect_s3_class(fold1$train_pred_patterns, "patterns")
  expect_s3_class(fold1$train_resp_patterns, "patterns")
  expect_s3_class(fold1$test_pred_data, "stars")
  expect_s3_class(fold1$test_resp_data, "stars")
  expect_equal(fold1$fold_id, 1L)

  # Check that patterns have max_k modes
  expect_equal(fold1$train_pred_patterns$k, 5)
  expect_equal(fold1$train_resp_patterns$k, 5)
})

test_that("truncation gives same CCA as direct computation", {
  skip_if_not_installed("testthat")

  prism_test <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
    readRDS()

  # Use full data (don't filter to subset that's too small)
  # Compute patterns at k=5, then truncate to k=3
  patterns_5 <- patterns(prism_test, k = 5, weight = FALSE)
  patterns_3_truncated <- patterns_5[1:3]

  # Compute directly at k=3
  patterns_3_direct <- patterns(prism_test, k = 3, weight = FALSE)

  # The eigenvalues should match
  expect_equal(
    patterns_3_truncated$eigenvalues$eigenvalues[1:3],
    patterns_3_direct$eigenvalues$eigenvalues[1:3]
  )

  # Amplitudes should match (up to sign)
  for (i in 1:3) {
    pc_name <- paste0("PC", i)
    truncated_amps <- patterns_3_truncated$amplitudes[[pc_name]]
    direct_amps <- patterns_3_direct$amplitudes[[pc_name]]

    # Check correlation is ~1 (allowing for sign flip)
    correlation <- cor(truncated_amps, direct_amps)
    expect_true(abs(correlation) > 0.999)
  }
})

test_that("k_cca is automatically set to min(k_pred, k_resp)", {
  skip_if_not_installed("testthat")

  prism_test <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
    readRDS()

  coarse <- prism_test %>%
    mutate(tmean = tmean * 0.8)

  cv <- prep_cv_folds(coarse, prism_test,
                      kfolds = 2,
                      max_k_pred = 3,
                      max_k_resp = 3,
                      weight = FALSE)

  # k_cca should be automatically set to min(k_pred, k_resp)
  results <- tune_cca(cv,
                      k_pred = 2:3,
                      k_resp = 2:3,
                      metrics = "rmse")

  # Check k_cca = min(k_pred, k_resp) for each row
  expect_true(all(results$k_cca == pmin(results$k_pred, results$k_resp)))
})

test_that("metrics are computed correctly", {
  skip_if_not_installed("testthat")

  prism_test <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
    readRDS()

  coarse <- prism_test %>%
    mutate(tmean = tmean * 0.8)

  cv <- prep_cv_folds(coarse, prism_test,
                      kfolds = 2,
                      max_k_pred = 3,
                      max_k_resp = 3,
                      weight = FALSE)

  # Test all metrics (k_cca will be auto-set to min(2,2)=2)
  results <- tune_cca(cv,
                      k_pred = 2,
                      k_resp = 2,
                      metrics = c("rmse", "cor_spatial", "cor_temporal"))

  # Check all metric columns present
  expect_true("rmse" %in% names(results))
  expect_true("cor_spatial" %in% names(results))
  expect_true("cor_temporal" %in% names(results))

  # RMSE should be positive
  expect_true(all(results$rmse > 0))

  # Correlations should be between -1 and 1
  expect_true(all(results$cor_spatial >= -1 & results$cor_spatial <= 1))
  expect_true(all(results$cor_temporal >= -1 & results$cor_temporal <= 1))
})

test_that("full workflow runs end-to-end", {
  skip_if_not_installed("testthat")

  prism_test <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
    readRDS()

  # Create coarse predictor
  coarse <- prism_test %>%
    mutate(tmean = tmean * 0.9 + units::set_units(rnorm(length(tmean), 0, 0.3), "°C"))

  # Step 1: Prepare folds
  cv <- prep_cv_folds(coarse, prism_test,
                      kfolds = 3,
                      max_k_pred = 5,
                      max_k_resp = 5,
                      weight = FALSE)

  expect_s3_class(cv, "cv_folds")

  # Step 2: Grid search (k_cca auto-set to min(k_pred, k_resp))
  results <- tune_cca(cv,
                      k_pred = 1:3,
                      k_resp = 1:3,
                      metrics = c("rmse", "cor_spatial"))

  expect_s3_class(results, "tbl_df")
  expect_true(nrow(results) > 0)

  # Step 3: Summarize
  summary <- summarize_cv(results, metric = "rmse", minimize = TRUE)

  expect_s3_class(summary, "tbl_df")
  expect_true("rmse_mean" %in% names(summary))
  expect_true("rmse_sd" %in% names(summary))

  # Best params should be attached
  best <- attr(summary, "best_params")
  expect_true(!is.null(best))
  expect_true("k_pred" %in% names(best))
  expect_true("k_resp" %in% names(best))
  expect_true("k_cca" %in% names(best))

  # Check RMSE values are reasonable (not NA, not infinite)
  expect_true(all(!is.na(summary$rmse_mean)))
  expect_true(all(is.finite(summary$rmse_mean)))
})

test_that("summarize_cv sorts correctly", {
  skip_if_not_installed("testthat")

  prism_test <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
    readRDS()

  coarse <- prism_test %>%
    mutate(tmean = tmean * 0.9)

  cv <- prep_cv_folds(coarse, prism_test,
                      kfolds = 2,
                      max_k_pred = 3,
                      max_k_resp = 3,
                      weight = FALSE)

  results <- tune_cca(cv,
                      k_pred = 1:2,
                      k_resp = 1:2,
                      metrics = c("rmse", "cor_spatial"))

  # Test minimize = TRUE (RMSE)
  summary_rmse <- summarize_cv(results, metric = "rmse", minimize = TRUE)
  expect_equal(summary_rmse$rmse_mean[1], min(summary_rmse$rmse_mean))

  # Test minimize = FALSE (correlation)
  summary_cor <- summarize_cv(results, metric = "cor_spatial", minimize = FALSE)
  expect_equal(summary_cor$cor_spatial_mean[1], max(summary_cor$cor_spatial_mean))
})

test_that("prep_cv_folds validates inputs", {
  skip_if_not_installed("testthat")

  prism_test <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
    readRDS()

  # Test with non-overlapping times
  pred_early <- filter(prism_test, time <= 1995)
  resp_late <- filter(prism_test, time >= 2005)

  expect_error(
    prep_cv_folds(pred_early, resp_late, kfolds = 3),
    "No common time steps"
  )
})

test_that("tune_cca validates k values", {
  skip_if_not_installed("testthat")

  prism_test <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
    readRDS()

  coarse <- prism_test %>%
    mutate(tmean = tmean * 0.9)

  cv <- prep_cv_folds(coarse, prism_test,
                      kfolds = 2,
                      max_k_pred = 3,
                      max_k_resp = 3,
                      weight = FALSE)

  # k_pred exceeds max_k_pred
  expect_error(
    tune_cca(cv, k_pred = 1:5, k_resp = 1:3),
    "exceeds max_k_pred"
  )

  # k_resp exceeds max_k_resp
  expect_error(
    tune_cca(cv, k_pred = 1:3, k_resp = 1:5),
    "exceeds max_k_resp"
  )
})

test_that("print.cv_folds works", {
  skip_if_not_installed("testthat")

  prism_test <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
    readRDS()

  coarse <- prism_test %>%
    mutate(tmean = tmean * 0.9)

  cv <- prep_cv_folds(coarse, prism_test,
                      kfolds = 2,
                      max_k_pred = 3,
                      max_k_resp = 3,
                      weight = FALSE)

  # Should not error
  expect_no_error(print(cv))

  # Check class and structure instead of output (cli output is hard to capture)
  expect_s3_class(cv, "cv_folds")
  expect_equal(cv$kfolds, 2)
  expect_equal(cv$max_k_pred, 3)
})
