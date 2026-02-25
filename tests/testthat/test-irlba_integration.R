test_that("perform_pca_smart uses base prcomp for small data", {
  set.seed(123)
  small_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)

  # Below threshold: no IRLBA message, just runs quietly
  result <- perform_pca_smart(small_matrix, k = 3, center = FALSE, size_threshold = 200)

  expect_true(is.list(result))
  expect_true("rotation" %in% names(result))
  expect_true("x" %in% names(result))
  # base prcomp returns all components, truncation happens downstream in get_eofs()
  expect_equal(ncol(result$rotation), ncol(small_matrix))
})

test_that("perform_pca_smart triggers IRLBA for large data", {
  set.seed(123)
  large_matrix <- matrix(rnorm(1000), nrow = 50, ncol = 20)

  if (requireNamespace("irlba", quietly = TRUE)) {
    expect_message(
      result <- perform_pca_smart(large_matrix, k = 3, center = FALSE, size_threshold = 500),
      "Using IRLBA"
    )
    expect_equal(ncol(result$rotation), 3)
  } else {
    expect_message(
      result <- perform_pca_smart(large_matrix, k = 3, center = FALSE, size_threshold = 500),
      "irlba.*unavailable"
    )
  }

  expect_true(is.list(result))
  expect_true("rotation" %in% names(result))
  expect_true("x" %in% names(result))
})

test_that("patterns accepts irlba_threshold parameter", {
  args <- formals(patterns)
  expect_true("irlba_threshold" %in% names(args))
  expect_equal(args$irlba_threshold, 500000)
})

test_that("perform_pca_smart requires size_threshold argument", {
  set.seed(123)
  test_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)

  expect_error(
    perform_pca_smart(test_matrix, k = 3, center = FALSE),
    "size_threshold"
  )
})

test_that("high threshold always uses base prcomp", {
  test_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)

  # Very high threshold — data size is well below, so no IRLBA message
  expect_no_message(
    result <- perform_pca_smart(test_matrix, k = 3, center = FALSE, size_threshold = 1e10)
  )
  expect_true(is.list(result))

  # Very low threshold should try IRLBA if available
  if (requireNamespace("irlba", quietly = TRUE)) {
    expect_message(
      perform_pca_smart(test_matrix, k = 3, center = FALSE, size_threshold = 1),
      "IRLBA"
    )
  }
})
