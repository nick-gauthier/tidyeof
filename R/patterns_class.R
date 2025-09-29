#' Print method for patterns objects
#'
#' @param x A patterns object
#' @param ... Additional arguments passed to print
#' @export
print.patterns <- function(x, ...) {
  cli::cli_h1("Patterns Object")
  cli::cli_text("Modes: {.field {x$k}}")

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

#' Summary method for patterns objects
#'
#' @param object A patterns object
#' @param ... Additional arguments passed to summary
#' @export
summary.patterns <- function(object, ...) {
  # Calculate total variance explained by retained PCs
  var_explained <- sum(head(object$eigenvalues$eigenvalues, object$k)) /
    sum(object$eigenvalues$eigenvalues) * 100

  # Get pattern correlations if rotated
  pattern_cors <- if(object$rotate && object$k > 1) {
    cors <- cor(object$amplitudes[,-1])  # exclude time column
    round(cors[upper.tri(cors)], 3)
  } else {
    NULL
  }

  # Create summary object
  structure(
    list(
      k = object$k,
      n_times = nrow(object$amplitudes),
      var_explained = var_explained,
      pattern_cors = pattern_cors,
      scaled = object$scaled,
      monthly = object$monthly,
      rotate = object$rotate,
      weight = object$weight,
      units = object$units
    ),
    class = "summary.patterns"
  )
}

#' Print method for summary.patterns objects
#'
#' @param x A summary.patterns object
#' @param ... Additional arguments passed to print
#' @export
print.summary.patterns <- function(x, ...) {
  cli::cli_h1("Patterns Summary")
  cli::cli_text("Number of modes: {.field {x$k}}")
  cli::cli_text("Time steps: {.field {x$n_times}}")
  cli::cli_text("Total variance explained: {.field {x$var_explained:.1f}%}")

  if(!is.null(x$pattern_cors)) {
    cat("\nPattern correlations:\n")
    print(x$pattern_cors)
  }

  cli::cli_h2("Processing Options")
  cli::cli_text("Scale: {.field {x$scaled}}")
  cli::cli_text("Monthly: {.field {x$monthly}}")
  cli::cli_text("Rotated: {.field {x$rotate}}")
  cli::cli_text("Area weighted: {.field {x$weight}}")

  invisible(x)
}

#' Format patterns as a data frame
#'
#' @param x A patterns object
#' @param ... Additional arguments
#' @export
as.data.frame.patterns <- function(x, ...) {
  # Convert EOF spatial patterns to long format
  eofs_df <- x$eofs %>%
    as_tibble() %>%
    tidyr::pivot_longer(
      -c(x, y),
      names_to = "PC",
      values_to = "weight"
    )

  # Convert amplitudes to long format
  amps_df <- x$amplitudes %>%
    tidyr::pivot_longer(
      -time,
      names_to = "PC",
      values_to = "amplitude"
    )

  list(
    eofs = eofs_df,
    amplitudes = amps_df,
    eigenvalues = x$eigenvalues
  )
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
  x$eofs <- x$eofs[, , i, drop = FALSE]
  x$amplitudes <- x$amplitudes[, c(1, i + 1)] # +1 because time is first column
  x$k <- length(i)

  if(!is.null(x$rotation)) {
    x$rotation <- x$rotation[, i, drop = FALSE]
  }

  # Eigenvalues stay the same but we update k
  x
}

#' Extract components from patterns object
#'
#' @param x A patterns object
#' @param name Name of component to extract
#' @return The requested component
#' @export
`$.patterns` <- function(x, name) {
  # Handle different types of requests
  if(grepl("^PC", name)) {
    # Extract PC amplitude time series
    pc_num <- as.numeric(sub("^PC", "", name))
    if(pc_num > x$k) rlang::abort(glue("Invalid PC number {pc_num}. Only {x$k} components available."), class = "tidyEOF_invalid_access")
    return(x$amplitudes[[paste0("PC", pc_num)]])
  }

  if(grepl("^EOF", name)) {
    # Extract EOF spatial pattern
    eof_num <- as.numeric(sub("^EOF", "", name))
    if(eof_num > x$k) rlang::abort(glue("Invalid EOF number {eof_num}. Only {x$k} components available."), class = "tidyEOF_invalid_access")
    return(x$eofs[, , eof_num])
  }

  # Default to standard list behavior
  x[[name]]
}
