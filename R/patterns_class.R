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
    if(any(is.na(i))) rlang::abort(glue("Invalid PC names. Available PCs: {glue_collapse(pc_names, sep = ', ')}"), class = "tidyEOF_invalid_subset")
  }

  if(any(i > x$k)) rlang::abort(glue("Index out of bounds. Requested: {glue_collapse(i[i > x$k], sep = ', ')}, but only {x$k} components available."), class = "tidyEOF_invalid_subset")

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

