#' Title
#'
#' @param target_patterns
#' @param amplitudes
#'
#' @return
#' @export
#'
#' @examples
#' Detect if variable should be non-negative
#' @param patterns A patterns object
#' @return Logical indicating if non-negative constraint should apply
#' @keywords internal
should_be_nonnegative <- function(patterns) {
  # First check units - precipitation is typically in mm, cm, or inches
  if(!is.null(patterns$units)) {
    for(unit in patterns$units) {
      if(!is.null(unit)) {
        unit_str <- as.character(unit$numerator)
        # Check for length units that typically indicate precipitation/depth
        if(any(unit_str %in% c("mm", "cm", "m", "in", "inches"))) {
          return(TRUE)
        }
      }
    }
  }

  # Fall back to checking variable names
  precip_patterns <- c("precip", "prec", "prcp", "rain", "snow", "ppt", "rr")
  var_names <- tolower(patterns$names)
  return(any(sapply(precip_patterns, function(p) grepl(p, var_names))))
}

reconstruct_field <- function(target_patterns, amplitudes = NULL, nonneg = "auto") {
  if(is.null(amplitudes)) amplitudes <- target_patterns$amplitudes
  if(is(amplitudes, 'stars')) amplitudes <- project_patterns(target_patterns, amplitudes)

  # Auto-detect whether to apply non-negative constraint
  if(nonneg == "auto") {
    nonneg <- should_be_nonnegative(target_patterns)
  }
  # check (ncol(amplitudes) - 1) == number of PCs in eofs?
  # check margin 3 is time?

  # this should work on new -raw- data and get the pc projection like cca

  anomalies <- amplitudes %>%
    rowwise() %>%
    mutate(PCs = list(c_across(-time)), .keep = 'unused') %>%
    ungroup() %>%
    tibble::deframe() %>%
    purrr::map(~sweep(target_patterns$eofs, MARGIN = 3, STATS = .x, FUN = "*")) %>%
    do.call('c', .) %>%
    stars::st_apply(1:2, sum) %>%
    merge(name = 'time') %>%
    stars::st_set_dimensions('time', values = amplitudes$time) %>%
    setNames(target_patterns$names)

  final <- restore_climatology(anomalies,
                               clim = target_patterns$climatology,
                               scale = target_patterns$scaled,
                               monthly = target_patterns$monthly)

  if(target_patterns$weight) {
    # Handle units properly for division by area weights
    weights <- area_weights(final)
    final <- units::drop_units(final) / weights
    # Units will be restored later in the function
  }
  if(nonneg) final <- mutate(final, across(everything(), ~pmax(.x, 0 * .x)))

  # Restore units for each variable
  if(!is.null(target_patterns$units) && any(!sapply(target_patterns$units, is.null))) {
    for(var_name in names(target_patterns$units)) {
      if(!is.null(target_patterns$units[[var_name]])) {
        final[[var_name]] <- units::set_units(final[[var_name]], target_patterns$units[[var_name]], mode = 'standard')
      }
    }
  }

  final
}