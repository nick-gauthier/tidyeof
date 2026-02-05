#' Reconstruct spatial field from EOF amplitudes
#'
#' Converts PC amplitudes back into a full spatial-temporal field by multiplying
#' by the EOF patterns and adding back the climatology.
#'
#' @param target_patterns A patterns object containing EOFs and climatology
#' @param amplitudes Amplitudes to use for reconstruction. Can be:
#'   - NULL (default): uses original amplitudes from target_patterns
#'   - tibble: with time column and PC columns
#'   - stars object: will be projected onto patterns first
#'
#' @details
#' For bounded variables like precipitation, the reconstructed field may contain
#' small negative values due to EOF truncation. To clamp these, use
#' `mutate(result, across(everything(), ~pmax(.x, 0 * .x)))` (the `0 * .x`
#' trick preserves units).
#'
#' @return A stars object with reconstructed spatial-temporal data
#' @export
reconstruct <- function(target_patterns, amplitudes = NULL) {
  validate_patterns(target_patterns)

  if(is.null(amplitudes)) amplitudes <- target_patterns$amplitudes
  if(inherits(amplitudes, 'stars')) amplitudes <- project_patterns(target_patterns, amplitudes)

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
    spatial_template = target_patterns$climatology$mean,
    valid_pixels = valid_pixels,
    times = amplitudes$time,
    var_names = target_patterns$names
  )

  final <- restore_climatology(anomalies,
                               clim = target_patterns$climatology,
                               scale = target_patterns$scaled,
                               monthly = target_patterns$monthly)

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
