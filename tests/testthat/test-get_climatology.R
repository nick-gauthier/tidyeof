prism <- system.file('testdata/prism_test.RDS', package = 'tidyEOF') %>%
  readRDS()

test_that("get_climatology() returns correct format", {
  # Basic format checks
  clim <- get_climatology(prism)
  clim_nounits <- get_climatology(units::drop_units(prism))
  clim_mon <- get_climatology(prism, cycle_length = 12)
  clim_mon_nounits <- get_climatology(units::drop_units(prism), cycle_length = 12)

  # Class checks
  expect_true(inherits(prism, "stars"))
  expect_true(inherits(clim, "stars"))
  expect_true(inherits(clim_nounits, "stars"))
  expect_true(inherits(clim_mon, "stars"))
  expect_true(inherits(clim_mon_nounits, "stars"))

  # Name checks
  expect_equal(names(prism), "tmean")
  expect_equal(names(clim), "tmean")
  expect_equal(names(clim_nounits), "tmean")
  expect_equal(names(clim_mon), "tmean")
  expect_equal(names(clim_mon_nounits), "tmean")

  # Units checks
  expect_true(inherits(prism[[1]], "units"))
  expect_true(inherits(clim[[1]], "units"))
  expect_true(inherits(clim_mon[[1]], "units"))
  expect_false(inherits(clim_nounits[[1]], "units"))
  expect_false(inherits(clim_mon_nounits[[1]], "units"))

  # Dimension checks
  expect_equal(length(dim(prism)), 3)
  expect_equal(length(dim(clim)), 3)
  expect_equal(length(dim(clim_nounits)), 3)
  expect_equal(length(dim(clim_mon)), 4)
  expect_equal(length(dim(clim_mon_nounits)), 4)

  # Monthly names check
  expect_equal(st_get_dimension_values(clim_mon, "group"), month.name)
  expect_equal(st_get_dimension_values(clim_mon_nounits, "group"), month.name)
})

test_that("get_climatology handles units correctly for multiple attributes", {
  # Create multi-attribute test data
  prism_multi <- c(temp = prism, precip = prism)

  # Set different units for precip while preserving values
  prec_vals <- units::drop_units(prism_multi$precip)
  prism_multi$precip <- prec_vals * units::as_units("mm")

  # Test annual climatology
  clim <- get_climatology(prism_multi)
  expect_equal(units(clim$temp), units(prism_multi$temp))
  expect_equal(units(clim$precip), units(prism_multi$precip))

  # Test monthly climatology
  clim_mon <- get_climatology(prism_multi, cycle_length = 12)
  expect_equal(units(clim_mon$temp), units(prism_multi$temp))
  expect_equal(units(clim_mon$precip), units(prism_multi$precip))

  # Test mixed units case
  prism_mixed <- c(temp = prism, precip = units::drop_units(prism))
  clim_mixed <- get_climatology(prism_mixed)
  expect_true(inherits(clim_mixed$temp, "units"))
  expect_false(inherits(clim_mixed$precip, "units"))

  # Test no units case
  prism_none <- c(temp = units::drop_units(prism),
                  precip = units::drop_units(prism))
  clim_none <- get_climatology(prism_none)
  expect_false(inherits(clim_none$temp, "units"))
  expect_false(inherits(clim_none$precip, "units"))
})

test_that("get_climatology() calculates correct values", {
  clim_mon <- get_climatology(prism, cycle_length = 12)

  # Test February averages
  feb_mean <- st_apply(prism[,,,c(2, 14, 26)], 1:2, mean, rename = FALSE)
  feb_sd <- st_apply(prism[,,,c(2, 14, 26)], 1:2, sd, rename = FALSE)

  expect_equal(
    units::drop_units(feb_mean),
    units::drop_units(abind::adrop(clim_mon[,,,2,1]))
  )
  expect_equal(
    units::drop_units(feb_sd),
    units::drop_units(abind::adrop(clim_mon[,,,2,2]))
  )
})