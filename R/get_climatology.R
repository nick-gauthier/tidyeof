#' Calculate climatological mean and standard deviation for spatial data
#'
#' Computes climatological statistics (mean and standard deviation) for a spatial field,
#' either annually or monthly. Preserves spatial dimensions and units from the input data.
#'
#' @param dat A stars object containing a spatial field with dimensions (x, y, time)
#' @param monthly Logical. If TRUE, computes monthly climatology. If FALSE (default),
#'   computes statistics over the entire period.
#'
#' @return A list with two stars objects:
#'   \item{mean}{Climatological mean with original spatial dimensions and units}
#'   \item{sd}{Climatological standard deviation with same structure}
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
    # FIXME: This warns, but get_anomalies() aborts on the same condition.
    # Pick one behavior and use it consistently.
    times <- st_get_dimension_values(dat, 'time')
    n_times <- length(times)
    if (n_times %% 12 != 0) {
      warning("Data does not contain complete years")
    }

    # Monthly climatology calculation
    if (has_geometry_dimension(dat)) {
      # For geometry stars: direct array manipulation avoids aggregate() collision
      # with the real geometry dimension name
      months <- by_months(times)
      month_levels <- levels(months)
      n_months <- length(month_levels)
      mat <- dat[[1]]  # geometry × time matrix
      n_geom <- nrow(mat)

      mean_arr <- array(NA_real_, dim = c(n_geom, n_months))
      sd_arr   <- array(NA_real_, dim = c(n_geom, n_months))
      for (m in seq_along(month_levels)) {
        idx <- which(months == month_levels[m])
        mean_arr[, m] <- rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
        sd_arr[, m]   <- apply(mat[, idx, drop = FALSE], 1, sd, na.rm = TRUE)
      }

      # Build result from template: take first n_months time slices, rename dim to month
      template <- dat[, , seq_len(n_months), drop = FALSE]
      mean_result <- template
      mean_result[[1]] <- mean_arr
      mean_result <- st_set_dimensions(mean_result, 'time',
                                       values = month_levels, names = 'month')
      sd_result <- template
      sd_result[[1]] <- sd_arr
      sd_result <- st_set_dimensions(sd_result, 'time',
                                     values = month_levels, names = 'month')
    } else {
      # For raster: use aggregate() with workaround for stars naming bug
      # NOTE: stars::aggregate with factor-returning grouping functions causes
      # the new dimension to be named 'geometry' instead of 'time'. We rename after.
      mean_result <- aggregate(dat, by_months, FUN = mean) %>%
        aperm(c(2, 3, 1)) %>%
        st_set_dimensions('geometry', names = 'month')

      sd_result <- aggregate(dat, by_months, FUN = sd) %>%
        aperm(c(2, 3, 1)) %>%
        st_set_dimensions('geometry', names = 'month')
    }
  } else {
    # Annual climatology calculation
    spatial_dims <- get_spatial_dimensions(dat)
    mean_result <- st_apply(dat, spatial_dims, mean, na.rm = TRUE, rename = FALSE)
    sd_result <- st_apply(dat, spatial_dims, sd, na.rm = TRUE, rename = FALSE)
  }

  list(
    mean = restore_units(mean_result, dat),
    sd = restore_units(sd_result, dat)
  )
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
  } else if (!is.list(clim) || !all(c("mean", "sd") %in% names(clim))) {
    rlang::abort("climatology must be a list with 'mean' and 'sd' stars objects (from get_climatology())",
                 class = "tidyeof_invalid_input")
  }

  mn <- clim$mean
  stdev <- clim$sd

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
    month_names <- unique(format(times, "%B"))

    if (has_geometry_dimension(dat)) {
      # Geometry path: direct array manipulation
      # NOTE (pre-existing): month_idx assumes clim month ordering matches
      # month_names from dat. If clim was computed from data with a different
      # start month (e.g., July-start clim applied to January-start dat),
      # the column indices will be wrong. The raster path has the same
      # assumption via st_redimension. Fix would be to index by month name.
      mat <- dat[[1]]  # geometry × time
      mn_mat <- mn[[1]]  # geometry × month
      sd_mat <- stdev[[1]]  # geometry × month
      month_idx <- match(format(times, "%B"), month_names)

      for (t in seq_len(n_times)) {
        m <- month_idx[t]
        mat[, t] <- mat[, t] - mn_mat[, m]
        if (scale) mat[, t] <- mat[, t] / sd_mat[, m]
      }

      out <- dat
      out[[1]] <- mat
      return(out)
    }

    # Raster path: redimension to monthly structure
    # NOTE (pre-existing, cosmetic): st_redimension sets names() on the
    # underlying array's dim attribute, so a roundtrip through redimension
    # can change dim names on dat[[1]] (e.g., X1/X2/X3 -> x/y/time).
    # Values are unaffected; only array metadata differs.
    spatial_dims <- get_spatial_dimensions(dat)
    redim_spec <- stats::setNames(
      c(dims[spatial_dims], month = 12L, year = as.integer(dims[["time"]] / 12)),
      c(spatial_dims, "month", "year")
    )
    dat <- st_redimension(dat, redim_spec) |>
      st_set_dimensions('month', values = month_names)
  }

  # calculate anomalies
    out <- dat - mn
    if (scale) {
      out <- out / stdev
    }

    if(monthly) {
      # Restore original time dimension (spatial_dims captured before redimension)
      redim_spec <- stats::setNames(
        c(dims[spatial_dims], time = dims[["time"]]),
        c(spatial_dims, "time")
      )
      out <- st_redimension(out, redim_spec) |>
        st_set_dimensions('time', values = times)
    }

  return(out)
}


#' Restore original field from anomalies and climatology
#'
#' Reverses the operation of \code{get_anomalies()}, adding the climatological
#' mean (and optionally multiplying by standard deviation) back to anomaly fields.
#'
#' @param anomalies A stars object containing anomalies (from \code{get_anomalies()})
#' @param clim A climatology list with \code{mean} and \code{sd} stars objects
#'   (from \code{get_climatology()})
#' @param scale Logical. If TRUE, multiply by standard deviation before adding mean
#'   (use when anomalies were standardized)
#' @param monthly Logical. If TRUE, restore using monthly climatology
#'
#' @return A stars object with the original field restored
#'
#' @examples
#' \dontrun{
#' clim <- get_climatology(dat)
#' anom <- get_anomalies(dat, clim)
#' restored <- restore_climatology(anom, clim)
#' }
#'
#' @export
restore_climatology <- function(anomalies, clim, scale = FALSE, monthly = FALSE) {
  # Basic validation
  if (!inherits(anomalies, "stars")) {
    rlang::abort("Anomalies must be a stars object", class = "tidyeof_invalid_input")
  }
  if (!is.list(clim) || !all(c("mean", "sd") %in% names(clim))) {
    rlang::abort("Climatology must be a list with 'mean' and 'sd' stars objects (from get_climatology())",
                 class = "tidyeof_invalid_input")
  }

  target_mean <- clim$mean
  target_sd <- clim$sd

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
    month_names <- unique(format(times, "%B"))

    if (has_geometry_dimension(anomalies)) {
      # Geometry path: direct array manipulation (reverse of get_anomalies)
      # NOTE (pre-existing): same month-ordering assumption as get_anomalies —
      # clim must have been computed from data with the same start month.
      mat <- anomalies[[1]]  # geometry × time
      mn_mat <- target_mean[[1]]  # geometry × month
      sd_mat <- target_sd[[1]]  # geometry × month
      month_idx <- match(format(times, "%B"), month_names)

      for (t in seq_len(n_times)) {
        m <- month_idx[t]
        if (scale) mat[, t] <- mat[, t] * sd_mat[, m]
        mat[, t] <- mat[, t] + mn_mat[, m]
      }

      out <- anomalies
      out[[1]] <- mat
      # Restore units and return early
      out <- restore_units(out, clim$mean)
      return(out)
    }

    # Raster path: redimension to monthly structure
    spatial_dims <- get_spatial_dimensions(anomalies)
    redim_spec <- stats::setNames(
      c(dims[spatial_dims], month = 12L, year = as.integer(dims[["time"]] / 12)),
      c(spatial_dims, "month", "year")
    )
    anomalies <- st_redimension(anomalies, redim_spec) |>
      st_set_dimensions('month', values = month_names)
  }
    # Restore climatology
    if (scale) {
      anomalies <- anomalies * target_sd
    }
    out <- anomalies + target_mean

    if(monthly) {
      # Restore original time dimension (spatial_dims captured before redimension)
      redim_spec <- stats::setNames(
        c(dims[spatial_dims], time = dims[["time"]]),
        c(spatial_dims, "time")
      )
      out <- st_redimension(out, redim_spec) |>
        st_set_dimensions('time', values = times)
    }

  # Restore units
  out <- restore_units(out, clim$mean)

  return(out)
}

# convenience function for monthly aggregation, based on example in aggregate.stars
# In get_climatology:
by_months <- function(x) {
  mon <- format(x, "%B")
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
