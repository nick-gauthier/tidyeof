# Test that project_patterns returns consistent results with stored amplitudes
# These tests catch bugs like sign flips from flip_patterns not being applied
# to pca$rotation

test_that("project_patterns returns same amplitudes as stored for original data", {
  pat <- patterns(prism, k = 5)
  stored <- pat$amplitudes
  projected <- project_patterns(pat, prism)

  # Should be identical (same data projected onto same patterns)
  expect_equal(projected, stored, tolerance = 1e-10)
})

test_that("project_patterns consistent after subsetting patterns", {
  pat <- patterns(prism, k = 5)

  for (i in 1:5) {
    pat_sub <- pat[1:i]
    stored_sub <- pat_sub$amplitudes
    projected_sub <- project_patterns(pat_sub, prism)

    expect_equal(
      projected_sub, stored_sub,
      tolerance = 1e-10,
      label = paste("k =", i)
    )
  }
})

test_that("project_patterns consistent with scale = TRUE", {
  pat <- patterns(prism, k = 5, scale = TRUE)
  stored <- pat$amplitudes
  projected <- project_patterns(pat, prism)

  expect_equal(projected, stored, tolerance = 1e-10)
})

test_that("project_patterns consistent with weight = FALSE", {
  pat <- patterns(prism, k = 5, weight = FALSE)
  stored <- pat$amplitudes
  projected <- project_patterns(pat, prism)

  expect_equal(projected, stored, tolerance = 1e-10)
})

test_that("project_patterns consistent with scale = TRUE and weight = FALSE", {
  pat <- patterns(prism, k = 5, scale = TRUE, weight = FALSE)
  stored <- pat$amplitudes
  projected <- project_patterns(pat, prism)

  expect_equal(projected, stored, tolerance = 1e-10)
})

test_that("reconstruct with amplitudes matches reconstruct with stars input", {
  pat <- patterns(prism, k = 5, scale = TRUE)

  # Reconstruct using stored amplitudes
  rec_stored <- reconstruct(pat)

  # Reconstruct by projecting original data

  rec_projected <- reconstruct(pat, prism)

  # Should be identical

  expect_equal(rec_projected, rec_stored, tolerance = 1e-10)
})

test_that("reconstruct consistent after subsetting", {
  pat <- patterns(prism, k = 5, scale = TRUE)

  for (i in 1:5) {
    pat_sub <- pat[1:i]

    rec_stored <- reconstruct(pat_sub)
    rec_projected <- reconstruct(pat_sub, prism)

    # Need to align dimensions for comparison
    rec_projected <- rec_projected %>%
      `st_dimensions<-`(st_dimensions(rec_stored))

    expect_equal(
      rec_projected, rec_stored,
      tolerance = 1e-10,
      label = paste("k =", i)
    )
  }
})
