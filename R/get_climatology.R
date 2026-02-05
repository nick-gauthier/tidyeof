#' Calculate climatological mean and standard deviation for spatial data
#'
#' Computes climatological statistics (mean and standard deviation) for a spatial field,
#' either annually or monthly. Preserves spatial dimensions and units from the input data.
#'
#' @param dat A stars object containing a spatial field with dimensions (x, y, time)
#' @param monthly Logical. If TRUE, computes monthly climatology. If FALSE (default),
#'   computes statistics over the entire period.
#'
#' @return A stars object with:
#'   - Original spatial dimensions (x, y)
#'   - For monthly climatology: Additional month dimension with values from month.name
#'   - Additional statistic dimension containing 'mean' and 'sd' statistics
#'   - Original units preserved from input data
#'
#' @examples
#' # Create sample data
#' library(stars)
#' times <- seq(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "month")
#' x <- seq(0, 1, length.out = 10)
#' y <- seq(0, 1, length.out = 10)
#' dat <- stars::st_as_stars(array(rnorm(10*10*36), c(10, 10, 36))) %>%
#'   st_set_dimensions(1, values = x, name = "x") %>%
#'   st_set_dimensions(2, values = y, name = "y") %>%
#'   st_set_dimensions(3, values = times, name = "time")
#'
#' # Calculate annual climatology
#' clim <- get_climatology(dat)
#'
#' # Calculate monthly climatology
#' monthly_clim <- get_climatology(dat, monthly = TRUE)
#'
#' @export
get_climatology <- function(dat, monthly = FALSE) {
  # Basic validation
  if (!inherits(dat, "stars")) {
    rlang::abort("Input must be a stars object", class = "tidyeof_invalid_input")
  }

  if (monthly) {
    # Check for complete years
    times <- st_get_dimension_values(dat, 'time')
    n_times <- length(times)
    if (n_times %% 12 != 0) {
      warning("Data does not contain complete years")
    }

    # Monthly climatology calculation
    # NOTE: stars::aggregate has a bug where factor-returning grouping functions
    # cause the new dimension to be named 'geometry' instead of 'time'. This happens
    # because `by` gets reassigned from function -> character vector (factor levels)
    # before the dimension naming logic runs, so is.function(by) returns FALSE.
    # See: https://github.com/r-spatial/stars/blob/main/R/aggregate.R lines 151-156, 220-224
    #
    # We work around this by renaming 'geometry' -> 'month' after aggregation.
    #
    # LIMITATION: This approach will BREAK for sf-based stars objects (with a real
    # 'geometry' dimension) because aggregate() will collide with the existing
    # geometry dimension. When we add support for sf-based stars, this will need
    # to be rewritten — likely using a manual grouping approach instead of aggregate().
    mean_result <- aggregate(dat, by_months, FUN = mean) %>%
      aperm(c(2, 3, 1)) %>%
      st_set_dimensions('geometry', names = 'month')

    sd_result <- aggregate(dat, by_months, FUN = sd) %>%
      aperm(c(2, 3, 1)) %>%
      st_set_dimensions('geometry', names = 'month')
  } else {
    # Annual climatology calculation
    spatial_dims <- get_spatial_dimensions(dat)
    mean_result <- st_apply(dat, spatial_dims, mean, na.rm = TRUE, rename = FALSE)
    sd_result <- st_apply(dat, spatial_dims, sd, na.rm = TRUE, rename = FALSE)
  }

  # Combine and restore units
  result <- c(mean = mean_result, sd = sd_result, along = 'statistic') %>%
    restore_units(dat)

  return(result)
}

#' Calculate anomalies from a climatological mean
#'
#' @param dat A stars object with dimensions (x, y, time)
#' @param clim Optional climatology from get_climatology(). If NULL, computed internally
#' @param scale Logical. If TRUE, divide by standard deviation
#' @param monthly Logical. If TRUE, compute monthly anomalies
#' @return A stars object with anomalies
#' @export
get_anomalies <- function(dat, clim = NULL, scale = FALSE, monthly = FALSE) {
  # Basic validation
  if (!inherits(dat, "stars")) {
    rlang::abort("Input must be a stars object", class = "tidyeof_invalid_input")
  }

  # Get or validate climatology
  if (is.null(clim)) {
    clim <- get_climatology(dat, monthly = monthly)
  } else if (!inherits(clim, "stars")) {
    rlang::abort("climatology must be a stars object", class = "tidyeof_invalid_input")
  }

  # Extract statistics
  mn <- slice(clim, 'statistic', 1)
  stdev <- slice(clim, 'statistic', 2)

  # Compute anomalies
  if (monthly) {
    # Get time values and dimensions
    times <- st_get_dimension_values(dat, 'time')
    dims <- dim(dat)
    n_times <- length(times)

    if (n_times %% 12 != 0) {
      rlang::abort("Data does not contain complete years. Need 12-month intervals for monthly climatology.", class = "tidyeof_invalid_time")
    }

    # Get the actual months in the data
    month_names <- as.character(unique(lubridate::month(times, label = TRUE, abbr = FALSE)))

    # Redimension to monthly structure and calculate anomalies
    # NOTE: Hardcoded x, y dimension names assume raster-based stars.
    # Will need revision for sf-based stars with geometry dimension.
    dat <- st_redimension(dat,
                                c(x = dims[[1]],
                                  y = dims[[2]],
                                  month = 12,
                                  year = dims[[3]] / 12)) |>
      st_set_dimensions('month', values = month_names)
  }

  # calculate anomalies
    out <- dat - mn
    if (scale) {
      out <- out / stdev
    }

    if(monthly) {
      # Restore original time dimension
      out <- st_redimension(out,
                            c(x = dims[[1]],
                              y = dims[[2]],
                              time = dims[[3]])) |>
        st_set_dimensions('time', values = times)
    }

  return(out)
}


#' @export
restore_climatology <- function(anomalies, clim, scale = FALSE, monthly = FALSE) {
  # Basic validation
  if (!inherits(anomalies, "stars")) {
    rlang::abort("Anomalies must be a stars object", class = "tidyeof_invalid_input")
  }
  if (!inherits(clim, "stars")) {
    rlang::abort("Climatology must be a stars object", class = "tidyeof_invalid_input")
  }

  # Extract statistics and handle units properly
  target_mean <- slice(clim, 'statistic', 1)
  target_sd <- slice(clim, 'statistic', 2)

  # Check if anomalies and climatology have matching units
  anomalies_has_units <- any(sapply(names(anomalies), function(var) !is.null(tryCatch(units(anomalies[[var]]), error = function(e) NULL))))
  clim_has_units <- any(sapply(names(target_mean), function(var) !is.null(tryCatch(units(target_mean[[var]]), error = function(e) NULL))))

  # If units mismatch, drop units from climatology to match anomalies
  if(!anomalies_has_units && clim_has_units) {
    target_mean <- units::drop_units(target_mean)
    target_sd <- units::drop_units(target_sd)
  }

  if (monthly) {
    # Get time values and dimensions
    times <- st_get_dimension_values(anomalies, 'time')
    dims <- dim(anomalies)
    n_times <- length(times)

    # Check for complete years
    if (n_times %% 12 != 0) {
      warning("Data does not contain complete years")
    }

    # Get the actual months in the data
    month_names <- as.character(unique(lubridate::month(times, label = TRUE, abbr = FALSE)))

    # Redimension to monthly structure
    # NOTE: Hardcoded x, y dimension names assume raster-based stars.
    # Will need revision for sf-based stars with geometry dimension.
    anomalies <- st_redimension(anomalies,
                                 c(x = dims[[1]],
                                   y = dims[[2]],
                                   month = 12,
                                   year = dims[[3]] / 12)) |>
      st_set_dimensions('month', values = month_names)
  }
    # Restore climatology
    if (scale) {
      anomalies <- anomalies * target_sd
    }
    out <- anomalies + target_mean

    if(monthly) {
      # Restore original time dimension
      # NOTE: Hardcoded x, y dimension names assume raster-based stars.
      # Will need revision for sf-based stars with geometry dimension.
      out <- st_redimension(out,
                            c(x = dims[[1]],
                              y = dims[[2]],
                              time = dims[[3]])) |>
        st_set_dimensions('time', values = times)
    }

  # Restore units
  out <- restore_units(out, clim)

  return(out)
}

# convenience function for monthly aggregation, based on example in aggregate.stars
# In get_climatology:
by_months <- function(x) {
  mon <- lubridate::month(x, label = TRUE, abbr = FALSE)
  factor(mon, levels = unique(mon))  # Preserve order from data
}


# Restore units from reference stars object to new stars object
restore_units <- function(new, ref) {
  for (var in names(new)) {
    ref_units <- tryCatch(units(ref[[var]]), error = function(e) NULL)
    if (!is.null(ref_units)) {
      new[[var]] <- units::set_units(new[[var]], ref_units, mode = 'standard')
    }
  }
  new
}

# Previous implementation (had potential bug if attribute order differed):
# restore_units <- function(new, ref) {
#   old_units <- purrr::map(ref, purrr::possibly(units))
#   apply_units <- function(x, y) purrr::modify_in(x, names(old_units)[y], ~units::set_units(.x, old_units[[y]], mode = 'standard'))
#   seq_along(new) %>%
#     purrr::reduce(apply_units, .init = new)
# }
