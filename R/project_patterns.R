#' Project New Data onto Existing Patterns
#'
#' This function projects new spatial-temporal data onto existing EOF patterns,
#' returning the corresponding principal component time series. This is a core
#' function used in pattern-based downscaling and reconstruction.
#'
#' @param patterns A patterns object containing EOFs, climatology, and other metadata
#' @param newdata A stars object with new spatial-temporal data to project
#'
#' @return A tibble with time column and PC amplitude columns
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get patterns from training data
#' patterns <- get_patterns(training_data, k = 5)
#'
#' # Project new data onto these patterns
#' new_amplitudes <- project_patterns(patterns, new_data)
#' }
project_patterns <- function(patterns, newdata) {

  validate_patterns(patterns)

  new_times <- st_get_dimension_values(newdata, 'time')

  anomalies <- newdata %>%
    get_anomalies(patterns$climatology,
                  scale = patterns$scaled,
                  monthly = patterns$monthly)

  # Apply area weighting if specified in patterns
  if (patterns$weight) {
    anomalies <- anomalies * area_weights(newdata)
  }

  anomalies <- units::drop_units(anomalies)

  flattened <- flatten_time_space(anomalies[1])
  data_matrix <- flattened$matrix

  valid_pixels <- if (!is.null(patterns$valid_pixels)) {
    patterns$valid_pixels
  } else {
    seq_len(ncol(data_matrix))
  }
  data_matrix <- data_matrix[, valid_pixels, drop = FALSE]

  complete_rows <- stats::complete.cases(data_matrix)
  if (!all(complete_rows)) {
    data_matrix <- data_matrix[complete_rows, , drop = FALSE]
    times_complete <- new_times[complete_rows]
  } else {
    times_complete <- new_times
  }

  if (nrow(data_matrix) == 0) {
    cli::cli_abort("No complete time steps available after removing missing values in `newdata`.")
  }

  rotation_names <- rownames(patterns$pca$rotation)
  if (!is.null(rotation_names)) {
    colnames(data_matrix) <- rotation_names
  }

  projected <- predict(patterns$pca, data_matrix) %>%
    .[, 1:patterns$k, drop = FALSE]

  # Apply rotation if present
  if (is.matrix(patterns$rotation)) {
    projected <- projected %*% patterns$rotation
  }

  # Convert to tibble with proper names
  result <- projected %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(names0(patterns$k, 'PC')) %>%
    mutate(time = times_complete, .before = 1)

  return(result)
}

#' Validate Spatial Compatibility Between Patterns and New Data
#'
#' Checks that new data is spatially compatible with existing patterns
#'
#' @param patterns A patterns object
#' @param newdata A stars object with new data
#'
#' @return Invisible(TRUE) if compatible, stops with error if not
#'
#' @keywords internal
validate_spatial_compatibility <- function(patterns, newdata) {

  # Check that newdata is a stars object
  if (!inherits(newdata, "stars")) {
    cli::cli_abort("`newdata` must be a {.cls stars} object.")
  }

  # Check that spatial dimensions match
  if (!has_geometry_dimension(newdata)) {
    cli::cli_abort("`newdata` must include spatial dimensions.")
  }

  # Check that spatial dimensions are compatible
  # Use existing get_spatial_dimensions function from get_patterns.R
  pattern_spatial_dims <- get_spatial_dimensions(patterns$climatology)
  newdata_spatial_dims <- get_spatial_dimensions(newdata)

  # Basic compatibility check - both should have spatial dimensions
  if (length(pattern_spatial_dims) == 0) {
    cli::cli_abort("`patterns` object has no spatial dimensions.")
  }

  if (length(newdata_spatial_dims) == 0) {
    cli::cli_abort("`newdata` has no spatial dimensions.")
  }

  # For more detailed validation, we'd need to check dimension sizes
  # but this is a basic compatibility check

  # Check variable names are compatible
  pattern_vars <- names(patterns$climatology)
  newdata_vars <- names(newdata)

  if (!all(pattern_vars %in% newdata_vars)) {
    missing_vars <- setdiff(pattern_vars, newdata_vars)
    cli::cli_abort(
      "`newdata` is missing variables present in patterns: {paste(missing_vars, collapse = ', ')}."
    )
  }

  return(invisible(TRUE))
}

# Removed get_spatial_dimensions_patterns to avoid conflicts with existing function

# Removed duplicate has_geometry_dimension function - using the one from get_patterns.R

#' Project Multiple Datasets onto Patterns
#'
#' Convenience function to project multiple datasets onto the same patterns
#'
#' @param patterns A patterns object
#' @param data_list List of stars objects to project
#' @param names Optional names for the datasets
#'
#' @return List of projected amplitude tibbles
#'
#' @export
project_patterns_multiple <- function(patterns, data_list, names = NULL) {

  if (!is.list(data_list)) {
    cli::cli_abort("`data_list` must be a list of {.cls stars} objects.")
  }

  # Apply projection to each dataset
  results <- purrr::map(data_list, ~ project_patterns(patterns, .x))

  # Add names if provided
  if (!is.null(names)) {
    if (length(names) != length(results)) {
      cli::cli_abort("Length of `names` must match length of `data_list`.")
    }
    names(results) <- names
  }

  return(results)
}
