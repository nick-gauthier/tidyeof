prism <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
  readRDS()

test_that("patterns() %>% reconstruct() give consistent outputs", {
  expect_equal(patterns(filter(prism), rotate = TRUE, k = 4) %>% reconstruct(),
               patterns(filter(prism), rotate = FALSE, k = 4) %>% reconstruct())
})
