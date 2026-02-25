library(testthat)
library(tidyeof)
prism <- system.file('testdata/prism_test.RDS', package = 'tidyeof') %>%
  readRDS()

test_check("tidyeof")
