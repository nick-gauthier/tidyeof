#' Compute climatological statistics for a spatial field
#'
#' @param dat A stars object with two spatial dimensions and a temporal dimension
#' @param cycle_length Integer length of the temporal cycle (e.g., 12 for monthly data).
#'   If NULL, computes statistics over entire period (default).
#'
#' @return A stars object with:
#'   - Original spatial dimensions (x, y)
#'   - Group dimension (cycle positions or month names for monthly data)
#'   - Statistic dimension (mean, sd)
#'   Original units are preserved.
#'
#' @examples
#' # Compute annual climatology
#' clim <- get_climatology(temperature)
#'
#' # Compute monthly climatology
#' monthly_clim <- get_climatology(temperature, cycle_length = 12)
#'
#' @export
#'
get_climatology <- function(dat, cycle_length = NULL) {
  # Basic input validation
  if (!inherits(dat, "stars")) {
    stop("Input must be a stars object with at least 3 dimensions")
  }

  # Store original units
  orig_units <- lapply(dat, function(x) {
    if (inherits(x, "units")) units(x)
  })

  # Handle entire period climatology
  if (is.null(cycle_length)) {
    mean_result <- st_apply(dat, 1:2, mean, na.rm = TRUE)
    sd_result <- st_apply(dat, 1:2, sd, na.rm = TRUE)

    final <- c(mean = mean_result, sd = sd_result, along = "statistic") %>%
      setNames(names(dat))
  }
  # Handle cyclic climatology
  else {
    # Validate cycle_length
    n_periods <- dim(dat)[[3]]
    if (cycle_length < 1 || cycle_length > n_periods || n_periods %% cycle_length != 0) {
      stop("Invalid cycle_length. Must be positive and divide total periods evenly.")
    }

    # Split into cycles and compute statistics
    positions <- seq_len(n_periods)
    groups <- split(positions, (positions - 1) %% cycle_length + 1)

    means <- lapply(groups, function(idx) {
      st_apply(dat[,,,idx], 1:2, mean, na.rm = TRUE)
    })
    sds <- lapply(groups, function(idx) {
      st_apply(dat[,,,idx], 1:2, sd, na.rm = TRUE)
    })

    # Set group names
    group_names <- seq_len(cycle_length)
    if (cycle_length == 12 &&
        inherits(st_get_dimension_values(dat, "time"), c("Date", "POSIXct"))) {
      group_names <- month.name
    }

    mean_result <- do.call(c, c(means, list(along = "group"))) %>%
      st_set_dimensions("group", values = group_names)
    sd_result <- do.call(c, c(sds, list(along = "group"))) %>%
      st_set_dimensions("group", values = group_names)

    final <- c(mean = mean_result, sd = sd_result, along = "statistic") %>%
      setNames(names(dat))
  }

  # Restore units
  for (attr_name in names(dat)) {
    if (!is.null(orig_units[[attr_name]])) {
      final[[attr_name]] <- units::set_units(final[[attr_name]],
                                             orig_units[[attr_name]],
                                             mode = "standard")
    }
  }

  return(final)
}

#' @export
get_anomalies <- function(dat, clim = NULL, scale = FALSE, monthly = FALSE) {

  if(is.null(clim)) clim <- get_climatology(dat, monthly = monthly)

  mn <- slice(clim, 'var', 1) # using slice instead of filter here because filter resets offset
  stdev <- slice(clim, 'var', 2)

  if(monthly) { # what if the months don't start with January or are uneven? does this still work?
    out <- sweep_months(dat, mn, '-')
    if(scale) {
      out <- sweep_months(out, stdev, '/')
    }

  } else {
    out <- dat - mn  # center the field
    if(scale) out <- out / stdev  ## scale the field (optional)
  }
  return(out)
}

#' @export
restore_climatology <- function(anomalies, clim, scale = FALSE, monthly = FALSE) {
  target_mean <- slice(clim, 'var', 1) %>%
    units::drop_units() # is this necessary
  target_sd <- slice(clim, 'var', 2) %>%
    units::drop_units()

  if(monthly) {
    if(scale) anomalies <- sweep_months(anomalies, target_sd, '*')

    final <- sweep_months(anomalies, target_mean, '+')
  } else {
    if(scale) anomalies <- anomalies * target_sd

    final <- anomalies + target_mean
  }
# should restore units here?
  return(final)
}

# convenience function for monthly aggregation, based on example in aggregate.stars
by_months = function(x) {
  lubridate::month(x, label = TRUE, abbr = FALSE)
}

# convenience functions for sweeping monthly summary statistics.
# should check that time is posix?
sweep_months <- function(e1, e2, FUN) {
  FUN <- match.fun(FUN)

  purrr::map(1:12, ~ FUN(filter(e1, lubridate::month(time) == .x), abind::adrop(filter(e2, month == month.name[.x])))) %>%
    do.call(c, .) %>%
    slice(., 'time', order(time(.))) %>% # reshuffle so months/years in right order
  st_set_dimensions(., 'time', values = time(.)) # the lubridate command above results in interval times not dates
}

# convenience function to add back units . . . there has to be a better way!
restore_units <- function(new, ref) {
  old_units <- purrr::map(ref, purrr::possibly(units))
  apply_units <- function(x, y) purrr::modify_in(x, names(old_units)[y], ~units::set_units(.x, old_units[[y]], mode = 'standard'))

  seq_along(new) %>%
    purrr::reduce(apply_units, .init = new)
}