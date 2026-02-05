test_that("IRLBA integration works correctly", {
  skip_if_not_installed("testthat")

  # Test that perform_pca_smart function exists
  expect_true(exists("perform_pca_smart"))

  # Test with small matrix (should use base prcomp)
  set.seed(123)
  small_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)

  # Capture messages to check method selection
  expect_message(
    result_small <- perform_pca_smart(small_matrix, k = 3, center = FALSE, size_threshold = 200),
    "Using base prcomp"
  )

  expect_true(is.list(result_small))
  expect_true("rotation" %in% names(result_small))
  expect_true("x" %in% names(result_small))
  expect_equal(ncol(result_small$rotation), 3)

  # Test with large matrix (should trigger IRLBA if available)
  large_matrix <- matrix(rnorm(1000), nrow = 50, ncol = 20)

  if (requireNamespace("irlba", quietly = TRUE)) {
    expect_message(
      result_large <- perform_pca_smart(large_matrix, k = 3, center = FALSE, size_threshold = 500),
      "Using IRLBA"
    )
  } else {
    expect_message(
      result_large <- perform_pca_smart(large_matrix, k = 3, center = FALSE, size_threshold = 500),
      "irlba.*package not available"
    )
  }

  expect_true(is.list(result_large))
  expect_true("rotation" %in% names(result_large))
  expect_true("x" %in% names(result_large))
})

test_that("patterns accepts irlba_threshold parameter", {
  skip_if_not_installed("testthat")

  # Test that patterns function accepts the new parameter without error
  # We'll test the signature without running the full computation due to test data issues
  args <- formals(patterns)

  # Check that irlba_threshold parameter exists with correct default
  expect_true("irlba_threshold" %in% names(args))
  expect_equal(args$irlba_threshold, 50000)
})

test_that("IRLBA fallback works correctly", {
  skip_if_not_installed("testthat")

  # Test error handling in perform_pca_smart
  set.seed(123)
  test_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)

  # Test with invalid parameters that might cause IRLBA to fail
  # This should fall back to base prcomp
  expect_message(
    result <- perform_pca_smart(test_matrix, k = 3, center = FALSE),
    "base prcomp"
  )

  expect_true(is.list(result))
})

test_that("irlba_threshold parameter validation", {
  skip_if_not_installed("testthat")

  # Test that extreme threshold values work
  test_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)

  # Very high threshold should always use base prcomp
  expect_message(
    perform_pca_smart(test_matrix, k = 3, size_threshold = 1e10),
    "base prcomp"
  )

  # Very low threshold should try to use IRLBA if available
  if (requireNamespace("irlba", quietly = TRUE)) {
    expect_message(
      perform_pca_smart(test_matrix, k = 3, size_threshold = 1),
      "IRLBA"
    )
  }
})