#' Compute Common EOF Patterns from Multiple Sources
#'
#' Performs joint PCA on multiple datasets that share a spatial domain. Each
#' dataset is anomalized with its own climatology, then concatenated along time
#' for joint PCA. The resulting shared spatial patterns get source-specific
#' amplitudes and climatologies, enabling CCA coupling and cross-source prediction.
#'
#' @param datasets Named list of stars objects sharing the same spatial grid.
#'   Names become the source identifiers used for extraction.
#' @param k Number of EOF modes to retain
#' @param scale Logical, whether to scale anomalies by standard deviation (default TRUE)
#' @param rotate Logical, whether to apply varimax rotation (default FALSE)
#' @param monthly Logical, whether to use monthly climatology (default FALSE)
#' @param weight Logical, whether to apply area weighting (default TRUE)
#' @param irlba_threshold Minimum data elements to trigger IRLBA (default 500000)
#'
#' @return A `common_patterns` S3 object. Source-specific patterns are extracted
#'   with `$` or `[[` using source names (e.g., `cpat$era`). Each extracted
#'   element is a standard `patterns` object with shared EOFs but source-specific
#'   climatology and amplitudes.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cpat <- common_patterns(
#'   list(era = era_coarse, phyda = phyda_coarse),
#'   k = 11, scale = TRUE
#' )
#'
#' # Extract source-specific patterns
#' cpat$era     # patterns object with ERA climatology + amplitudes
#' cpat$phyda   # patterns object with PHYDA climatology + amplitudes
#'
#' # Use with existing couple/predict workflow
#' fine_pat <- patterns(era_fine, k = 13)
#' coupled <- couple(cpat$era, fine_pat, k = 9)
#' predict(coupled, era_new)
#'
#' # Cross-source prediction
#' predict(coupled, phyda_new, predictor_patterns = cpat$phyda)
#' }
common_patterns <- function(datasets, k = 4, scale = TRUE, rotate = FALSE,
                            monthly = FALSE, weight = TRUE,
                            irlba_threshold = 500000) {

  # --- Input validation ---
  if (!is.list(datasets) || length(datasets) == 0) {
    cli::cli_abort(
      "{.arg datasets} must be a non-empty list of {.cls stars} objects.",
      class = "tidyeof_invalid_input"
    )
  }

  if (is.null(names(datasets)) || any(names(datasets) == "")) {
    cli::cli_abort(
      "{.arg datasets} must be a named list. Each element needs a source name.",
      class = "tidyeof_invalid_input"
    )
  }

  for (nm in names(datasets)) {
    if (!inherits(datasets[[nm]], "stars")) {
      cli::cli_abort(
        "Element {.val {nm}} of {.arg datasets} must be a {.cls stars} object, not {.cls {class(datasets[[nm]])}}.",
        class = "tidyeof_invalid_input"
      )
    }
  }

  if (isTRUE(rotate) && k <= 1) {
    cli::cli_abort(
      "Rotation requires k > 1.",
      class = "tidyeof_invalid_option"
    )
  }

  source_names <- names(datasets)
  n_sources <- length(datasets)

  # --- Per-source: compute climatology and anomalies ---
  source_climatologies <- purrr::map(datasets, get_climatology, monthly = monthly)
  source_anomalies <- purrr::imap(datasets, function(dat, nm) {
    get_anomalies(dat, clim = source_climatologies[[nm]], scale = scale, monthly = monthly)
  })

  # Capture units from first dataset
  first_dat <- datasets[[1]]
  original_units <- setNames(
    purrr::map(names(first_dat), ~tryCatch(units(first_dat[[.x]]), error = function(e) NULL)),
    names(first_dat)
  )

  # --- Compute spatial weights once from first dataset ---
  weights <- if (weight) area_weights(first_dat) else NULL

  # --- Per-source: flatten, find valid pixels, weight ---
  flattened_sources <- purrr::map(source_anomalies, function(anom) {
    flat <- flatten_time_space(units::drop_units(anom)[1])
    valid <- which(!apply(flat$matrix, 2, anyNA))
    list(flat = flat, valid = valid)
  })

  # --- Compute valid_pixels intersection ---
  all_valid <- purrr::map(flattened_sources, "valid")
  shared_valid <- Reduce(intersect, all_valid)

  if (length(shared_valid) == 0) {
    cli::cli_abort(
      "No spatial pixels are valid across all sources. Check for non-overlapping NA patterns.",
      class = "tidyeof_no_valid_pixels"
    )
  }

  n_pixels_total <- ncol(flattened_sources[[1]]$flat$matrix)

  # Compute weights for shared valid pixels
  if (!is.null(weights)) {
    if (length(weights) != n_pixels_total) {
      cli::cli_abort(
        "Length of weights ({length(weights)}) must match number of spatial points ({n_pixels_total}).",
        class = "tidyeof_weight_mismatch"
      )
    }
    weights_valid <- weights[shared_valid]
  } else {
    weights_valid <- rep(1, length(shared_valid))
  }

  # --- Per-source: extract valid pixels, apply weights, track row ranges ---
  weighted_matrices <- list()
  row_ranges <- list()
  source_times <- list()
  current_row <- 1L

  for (nm in source_names) {
    mat <- flattened_sources[[nm]]$flat$matrix[, shared_valid, drop = FALSE]
    mat_weighted <- sweep(mat, 2, weights_valid, `*`)
    n_rows <- nrow(mat_weighted)

    weighted_matrices[[nm]] <- mat_weighted
    row_ranges[[nm]] <- c(start = current_row, end = current_row + n_rows - 1L)
    source_times[[nm]] <- stars::st_get_dimension_values(datasets[[nm]], "time")
    current_row <- current_row + n_rows
  }

  # --- Concatenate all weighted anomaly matrices ---
  concat_matrix <- do.call(rbind, weighted_matrices)

  # --- Validate k ---
  max_k <- min(nrow(concat_matrix) - 1, length(shared_valid))
  check_k_valid(k, max_k)

  # --- PCA on concatenated matrix ---
  pca_result <- perform_pca_smart(
    concat_matrix,
    k = k,
    center = FALSE,
    size_threshold = irlba_threshold
  )

  # Keep matrices for synchronized operations
  loadings_matrix <- pca_result$rotation[, 1:k, drop = FALSE]
  scores_matrix <- pca_result$x[, 1:k, drop = FALSE]
  sdev_vector <- pca_result$sdev[1:k]

  rotation_matrix <- NULL

  # --- Optional rotation ---
  if (rotate && k > 1) {
    rotation_result <- rotate_pca_components(loadings_matrix, scores_matrix, sdev_vector)
    loadings_weighted <- rotation_result$loadings
    all_scores <- rotation_result$scores
    rotation_matrix <- rotation_result$rotation_matrix
    component_sdev <- rotation_result$sdev
    component_variance <- rotation_result$eigenvalues
  } else {
    loadings_weighted <- loadings_matrix
    all_scores <- scores_matrix
    component_sdev <- sdev_vector
    component_variance <- sdev_vector^2
  }

  # --- Remove weights from loadings to get physical-unit EOFs ---
  loadings <- sweep(loadings_weighted, 1, weights_valid, `/`)

  # --- Build shared EOF stars object ---
  pc_names <- names0(k, 'PC')
  dims <- dim(first_dat)

  full_patterns <- array(NA, dim = c(k, n_pixels_total))
  full_patterns[, shared_valid] <- t(loadings)

  if (has_geometry_dimension(first_dat)) {
    template <- first_dat[,,1:k, drop = FALSE]
    spatial_patterns <- template %>%
      stars::st_set_dimensions('time', values = pc_names, names = 'PC') %>%
      setNames("weight")
    spatial_patterns[[1]] <- t(full_patterns)
  } else {
    template <- first_dat[,,,1:k, drop = FALSE]
    spatial_patterns <- template %>%
      stars::st_set_dimensions('time', values = pc_names, names = 'PC') %>%
      setNames("weight")
    pattern_array <- array(full_patterns, dim = c(k, dims[[1]], dims[[2]]))
    reordered_patterns <- aperm(pattern_array, c(2, 3, 1))
    spatial_patterns[[1]] <- reordered_patterns
  }

  # --- Build eigenvalues tibble ---
  n_total <- length(pca_result$sdev)
  eigenvalues <- pca_result |>
    broom::tidy(matrix = 'pcs') |>
    mutate(eigenvalues = std.dev ^ 2,
           percent = percent * 100,
           cumulative = cumulative * 100,
           error = sqrt(2 / n_total),
           low = eigenvalues * (1 - error) * 100 / sum(eigenvalues),
           hi = eigenvalues * (1 + error) * 100 / sum(eigenvalues),
           cumvar_line = hi + 0.02 * max(hi))

  if (rotate && k > 1) {
    total_variance <- sum(pca_result$sdev^2)
    rotated_percent <- component_variance / total_variance * 100
    rotated_cumulative <- cumsum(rotated_percent)

    eigenvalues <- eigenvalues %>%
      mutate(
        std.dev = if_else(PC <= k, component_sdev[PC], std.dev),
        eigenvalues = if_else(PC <= k, component_variance[PC], eigenvalues),
        percent = if_else(PC <= k, rotated_percent[PC], percent),
        cumulative = if_else(PC <= k, rotated_cumulative[PC], cumulative)
      )
  }

  # --- Compute shared sign vector from shared EOFs ---
  # Build a temporary patterns-like object just to compute signs
  signs <- compute_eof_signs(spatial_patterns)

  # --- Per-source: create patterns objects ---
  source_patterns <- list()

  for (nm in source_names) {
    rr <- row_ranges[[nm]]
    source_scores <- all_scores[rr["start"]:rr["end"], , drop = FALSE]
    colnames(source_scores) <- pc_names

    source_amps <- source_scores %>%
      as_tibble() %>%
      mutate(time = source_times[[nm]], .before = 1)

    # Capture source-specific units
    src_dat <- datasets[[nm]]
    src_units <- setNames(
      purrr::map(names(src_dat), ~tryCatch(units(src_dat[[.x]]), error = function(e) NULL)),
      names(src_dat)
    )

    pat <- new_patterns(
      eofs = spatial_patterns,
      amplitudes = source_amps,
      eigenvalues = eigenvalues,
      k = k,
      proj_matrix = loadings_weighted,
      center = pca_result$center,
      scale = pca_result$scale,
      rotation = rotation_matrix,
      climatology = source_climatologies[[nm]],
      units = src_units,
      names = names(src_dat),
      scaled = scale,
      monthly = monthly,
      rotate = rotate,
      weight = weight,
      valid_pixels = shared_valid
    )

    # Apply consistent signs across all sources
    pat <- apply_sign_flips(pat, signs)

    source_patterns[[nm]] <- pat
  }

  # --- Build common_patterns S3 object ---
  structure(
    list(
      source_patterns = source_patterns,
      sources = source_names,
      k = k,
      shared_eigenvalues = eigenvalues,
      pattern_opts = list(
        scale = scale,
        rotate = rotate,
        monthly = monthly,
        weight = weight
      )
    ),
    class = "common_patterns"
  )
}

#' @export
`$.common_patterns` <- function(x, name) {
  sources <- .subset2(x, "sources")
  if (name %in% sources) {
    .subset2(.subset2(x, "source_patterns"), name)
  } else {
    .subset2(x, name)
  }
}

#' @export
`[[.common_patterns` <- function(x, i, ...) {
  sources <- .subset2(x, "sources")
  if (is.character(i) && length(i) == 1 && i %in% sources) {
    .subset2(.subset2(x, "source_patterns"), i)
  } else {
    .subset2(x, i)
  }
}

#' @export
print.common_patterns <- function(x, ...) {
  cli::cli_h1("Common Patterns Object")
  cli::cli_text("Sources: {.val {x$sources}}")
  cli::cli_text("Shared EOF modes: {.field {x$k}}")

  cli::cli_h2("Time steps per source")
  for (nm in x$sources) {
    n_times <- nrow(x$source_patterns[[nm]]$amplitudes)
    cli::cli_text("{.field {nm}}: {n_times}")
  }

  cli::cli_h2("Processing Options")
  cli::cli_text("Scale: {.field {x$pattern_opts$scale}}")
  cli::cli_text("Rotate: {.field {x$pattern_opts$rotate}}")
  cli::cli_text("Monthly: {.field {x$pattern_opts$monthly}}")
  cli::cli_text("Weight: {.field {x$pattern_opts$weight}}")

  invisible(x)
}
