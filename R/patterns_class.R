#' Create a patterns object with validation
#'
#' Low-level constructor for patterns objects. Validates that all required
#' components are present and have the correct types.
#'
#' @param eofs A stars object containing EOF spatial patterns
#' @param amplitudes A data.frame/tibble with time column and PC amplitudes
#' @param eigenvalues A data.frame/tibble with eigenvalue statistics
#' @param k Number of components retained
#' @param proj_matrix Projection matrix for new data
#' @param center Centering values from PCA (or FALSE)
#' @param scale Scaling values from PCA (or FALSE)
#' @param rotation Rotation matrix if varimax was applied (or NULL)
#' @param climatology Stars object with climatology
#' @param units List of original data units
#' @param names Variable names
#' @param scaled Logical, whether data was scaled
#' @param monthly Logical, whether monthly climatology was used
#' @param rotate Logical, whether rotation was applied
#' @param weight Logical, whether area weighting was applied
#' @param valid_pixels Indices of valid (non-NA) pixels
#'
#' @return A patterns object
#' @keywords internal
new_patterns <- function(eofs, amplitudes, eigenvalues, k, proj_matrix,
                         center = FALSE, scale = FALSE, rotation = NULL,
                         climatology = NULL, units = NULL, names = NULL,
                         scaled = FALSE, monthly = FALSE, rotate = FALSE,
                         weight = TRUE, valid_pixels = NULL) {
  stopifnot(inherits(eofs, "stars"))
  stopifnot(is.data.frame(amplitudes))
  stopifnot(is.data.frame(eigenvalues))
  stopifnot(is.numeric(k), length(k) == 1)
  stopifnot(is.matrix(proj_matrix))

  structure(
    list(
      eofs = eofs,
      amplitudes = amplitudes,
      climatology = climatology,
      eigenvalues = eigenvalues,
      rotation = rotation,
      units = units,
      names = names,
      scaled = scaled,
      monthly = monthly,
      rotate = rotate,
      k = k,
      weight = weight,
      valid_pixels = valid_pixels,
      proj_matrix = proj_matrix,
      center = center,
      scale = scale
    ),
    class = "patterns"
  )
}

#' Print method for patterns objects
#'
#' @param x A patterns object
#' @param ... Additional arguments passed to print
#' @export
print.patterns <- function(x, ...) {
  cli::cli_h1("Patterns Object")
  cli::cli_text("Modes: {.field {x$k}}")
  cli::cli_text("Time steps: {.field {nrow(x$amplitudes)}}")

  cli::cli_h2("Processing Options")
  cli::cli_text("Scale: {.field {x$scaled}}")
  cli::cli_text("Monthly: {.field {x$monthly}}")
  cli::cli_text("Rotated: {.field {x$rotate}}")
  cli::cli_text("Area weighted: {.field {x$weight}}")

  cli::cli_h2("Eigenvalues (% variance)")
  top_eigs <- head(x$eigenvalues, x$k)
  variance_pct <- round(top_eigs$eigenvalues / sum(x$eigenvalues$eigenvalues) * 100, 1)
  cli::cli_text("{.val {variance_pct}}")

  invisible(x)
}

#' Extract subset of PCs from a patterns object
#'
#' @param x A patterns object
#' @param i Index or names of PCs to keep
#' @param ... Additional arguments (unused)
#' @return A patterns object with only selected PCs
#' @export
`[.patterns` <- function(x, i, ...) {
  # Handle different types of indexing
  if(is.character(i)) {
    pc_names <- paste0("PC", seq_len(x$k))
    i <- match(i, pc_names)
    if(any(is.na(i))) rlang::abort(glue("Invalid PC names. Available PCs: {glue_collapse(pc_names, sep = ', ')}"), class = "tidyeof_invalid_subset")
  }

  if(any(i > x$k)) rlang::abort(glue("Index out of bounds. Requested: {glue_collapse(i[i > x$k], sep = ', ')}, but only {x$k} components available."), class = "tidyeof_invalid_subset")

  # Subset components while preserving scaling
  x$eofs <- dplyr::slice(x$eofs, i, along = "PC", drop = FALSE)
  x$amplitudes <- x$amplitudes[, c(1, i + 1)] # +1 because time is first column
  x$k <- length(i)

  if(!is.null(x$rotation)) {
    x$rotation <- x$rotation[, i, drop = FALSE]
  }

  if (!is.null(x$proj_matrix)) {
    x$proj_matrix <- x$proj_matrix[, i, drop = FALSE]
  }

  # Eigenvalues stay the same but we update k
  x
}

