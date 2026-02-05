
test_that("patterns() %>% reconstruct() give consistent outputs", {
  expect_equal(patterns(filter(prism), rotate = TRUE, k = 4) %>% reconstruct(nonneg = FALSE),
               patterns(filter(prism), rotate = FALSE, k = 4) %>% reconstruct(nonneg = FALSE))
})
