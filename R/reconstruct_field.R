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
  validate_patterns(target_patterns)

  if(is.null(amplitudes)) amplitudes <- target_patterns$amplitudes
  if(is(amplitudes, 'stars')) amplitudes <- project_patterns(target_patterns, amplitudes)

  # Auto-detect whether to apply non-negative constraint
  if(nonneg == "auto") {
    nonneg <- should_be_nonnegative(target_patterns)
  }
  # check (ncol(amplitudes) - 1) == number of PCs in eofs?
  # check margin 3 is time?

  # With Hannachi-style handling we store EOFs in physical units and PCs as
  # unscaled scores, so multiplying amplitudes by EOFs yields anomalies directly

  amps_matrix <- amplitudes %>%
    dplyr::select(-time) %>%
    as.matrix()

  if (ncol(amps_matrix) != target_patterns$k) {
    cli::cli_abort(
      "Amplitude matrix columns ({ncol(amps_matrix)}) must match {.field k} ({target_patterns$k})."
    )
  }

  eof_array <- target_patterns$eofs[[1]]
  spatial_sizes <- dim(eof_array)[-length(dim(eof_array))]
  eof_matrix <- matrix(eof_array,
                       nrow = prod(spatial_sizes),
                       ncol = target_patterns$k)

  valid_pixels <- target_patterns$valid_pixels
  eof_valid <- eof_matrix[valid_pixels, , drop = FALSE]

  anomalies_valid <- amps_matrix %*% t(eof_valid)

  anomalies <- matrix_to_spacetime(
    anomalies_valid,
    template_eofs = target_patterns$eofs,
    spatial_template = target_patterns$climatology,
    valid_pixels = valid_pixels,
    times = amplitudes$time,
    var_names = target_patterns$names
  )

  final <- restore_climatology(anomalies,
                               clim = target_patterns$climatology,
                               scale = target_patterns$scaled,
                               monthly = target_patterns$monthly)

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
