#' Calculate teleconnections between a patterns object and another layer
#'
#' Computes pixel-wise correlations between a spatiotemporal field and each
#' PC amplitude time series. Checks for overlapping time steps between the
#' raster field and the PC amplitudes.
#'
#' @param dat A stars object with a time dimension
#' @param patterns A patterns object from patterns()
#' @param amplitudes Optional amplitudes tibble (defaults to patterns$amplitudes)
#'
#' @return A stars object with correlation values for each PC
#' @export
get_correlation <- function(dat, patterns, amplitudes = NULL) {
  if(is.null(amplitudes)) amplitudes <- patterns$amplitudes

  times_amps <- amplitudes$time
  times_dat <- st_get_dimension_values(dat, 'time')
  # Use subsetting to preserve Date/POSIXct class (intersect can strip it)
  times_cor <- times_amps[times_amps %in% times_dat]

  if(length(times_cor) <= 0) cli::cli_abort("Need at least two time steps in common.")
  if(!(identical(times_cor, times_amps) & identical(times_cor, times_dat))) cli::cli_inform("Using the time period {first(times_cor)} to {last(times_cor)}.")

  amps <- filter(amplitudes, time %in% times_cor) %>%
    select(-time)

  suppressWarnings( # suppress warnings that sd is zero
  filter(dat, time %in% times_cor) %>%
  st_apply(get_spatial_dimensions(.), function(x) cor(x, amps), .fname = 'PC') %>%
    st_set_dimensions(., 'PC', values = paste0('PC', st_get_dimension_values(., 'PC'))) %>%
    aperm(c(2,3,1))
  )
}

#' Calculate FDR-corrected significance contours for teleconnections
#'
#' Computes pixel-wise correlations with FDR correction and returns significance
#' contour lines as sf polygons.
#'
#' @param dat A stars object with a time dimension
#' @param patterns A patterns object from patterns()
#' @param fdr False discovery rate threshold (default 0.1)
#' @param amplitudes Optional amplitudes tibble (defaults to patterns$amplitudes)
#'
#' @return An sf object with significance contour polygons for each PC
#' @export
get_fdr <- function(dat, patterns, fdr = 0.1, amplitudes = NULL) {
  if(is.null(amplitudes)) amplitudes <- patterns$amplitudes

  times_amps <- amplitudes$time
  times_dat <- st_get_dimension_values(dat, 'time')
  # Use subsetting to preserve Date/POSIXct class (intersect can strip it)
  times_cor <- times_amps[times_amps %in% times_dat]

  if(length(times_cor) <= 0) cli::cli_abort("Need at least two time steps in common.")
  if(!(identical(times_cor, times_amps) & identical(times_cor, times_dat))) cli::cli_inform("Using the time period {first(times_cor)} to {last(times_cor)}.")

  amps <- filter(amplitudes, time %in% times_cor) %>%
    select(-time)

  suppressWarnings( # suppress warnings that sd is zero
    fdr_rast <- filter(dat, time %in% times_cor) %>%
      st_apply(get_spatial_dimensions(.), fdr_fun, amps = amps, .fname = 'PC') %>%
      aperm(c(2,3,1)) %>%
      st_apply('PC', adjust) %>%
      setNames('FDR')
  )

  fdr_rast %>%
    st_get_dimension_values('PC') %>%
    seq_along() %>%
    map(~slice(fdr_rast, 'PC', .x) %>%
          st_contour(contour_lines = TRUE, breaks = fdr) %>%
          transmute(PC = paste0('PC', .x))) %>%
    do.call(rbind, .)
}

fdr_fun <- function(x, amps) {
  if(!any(is.na(x))) apply(amps, 2, function(y) cor.test(x, y)$p.value) else rep(NA, ncol(amps))
}

adjust <- function(x) {
  p.adjust(x, method = 'fdr', n = sum(!is.na(x)))
}

