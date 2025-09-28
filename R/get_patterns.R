# Input validation helpers ----

#' Check if input is a stars object
#' @param x Object to check
#' @param arg Argument name for error messages
#' @param call Calling environment for error messages
#' @keywords internal
check_stars_object <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!inherits(x, "stars")) {
    abort(
      glue("Argument {.arg {arg}} must be a stars object, not {.cls {class(x)}}."),
      class = "tidyEOF_invalid_input",
      call = call
    )
  }
}

#' Check if k is valid
#' @param k Number of components
#' @param max_k Maximum allowed components
#' @param arg Argument name for error messages
#' @param call Calling environment for error messages
#' @keywords internal
check_k_valid <- function(k, max_k, arg = rlang::caller_arg(k), call = rlang::caller_env()) {
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k > max_k) {
    abort(
      glue("Argument {.arg {arg}} must be a single integer between 1 and {max_k}."),
      class = "tidyEOF_invalid_k",
      call = call
    )
  }
}

#' Get EOFs and PCs from spatiotemporal data
#'
#' @param dat A `stars` object containing spatial and temporal dimensions
#' @param k The number of PC/EOF modes to retain
#' @param scale Logical, whether to scale before PCA
#' @param rotate Logical, whether to apply Varimax rotation
#' @param monthly Logical, whether to use monthly climatology
#' @param weight Logical, whether to apply latitude weighting
#'
#' @return A `patterns` object containing EOFs, amplitudes, and metadata
#' @export
get_patterns <- function(dat, k = 4, scale = FALSE, rotate = FALSE, monthly = FALSE, weight = TRUE){

  # Input validation
  check_stars_object(dat)

  # Capture units from original data before any modifications
  original_units <- setNames(purrr::map(names(dat), ~tryCatch(units(dat[[.x]]), error = function(e) NULL)), names(dat))

  # Apply latitude weighting if requested
  if(weight) dat <- dat * area_weights(dat) # do we do this before or after scaling?

  # Get climatology and anomalies
  climatology <- get_climatology(dat, monthly = monthly)
  anomalies <- get_anomalies(dat, clim = climatology, scale = scale, monthly = monthly)

  # Get EOFs, which includes PCA/rotation handling internally
  eofs <- get_eofs(anomalies, k = k, rotate = rotate)

  # Create output pattern object with all components
  patterns <- list(
    eofs = eofs$spatial_patterns,
    amplitudes = eofs$amplitudes,
    climatology = climatology,
    pca = eofs$pca,
    eigenvalues = eofs$eigenvalues,
    rotation = eofs$rotation_matrix,
    units = original_units,
    names = names(dat),
    scaled = scale,
    monthly = monthly,
    rotate = rotate,
    k = k,
    weight = weight
  )

  # Align patterns so EOFs have roughly similar dominant signs
  patterns <- flip_patterns(patterns)

  class(patterns) <- "patterns"
  return(patterns)
}

#' Matrix-centric rotation helper that keeps all components synchronized
#' @keywords internal
rotate_pca_components <- function(loadings_matrix, scores_matrix, sdev_vector, k) {
  # Varimax rotation with scaling for robust rotation
  rot <- varimax(loadings_matrix %*% diag(sdev_vector))

  # Extract rotated components
  rotated_loadings <- unclass(rot$loadings)
  rotation_matrix <- rot$rotmat

  # Calculate rotated eigenvalues for reordering
  rotated_eigenvals <- diag(t(rotated_loadings) %*% rotated_loadings)

  # Reorder by explained variance
  ev_order <- order(rotated_eigenvals, decreasing = TRUE)

  # Apply reordering to ALL components simultaneously (this prevents sync issues)
  final_loadings <- rotated_loadings[, ev_order, drop = FALSE]
  final_rotation_matrix <- rotation_matrix[, ev_order, drop = FALSE]
  final_scores <- scores_matrix %*% final_rotation_matrix
  final_sdev <- sqrt(rotated_eigenvals[ev_order])

  list(
    loadings = final_loadings,
    scores = final_scores,
    sdev = final_sdev,
    rotation_matrix = final_rotation_matrix
  )
}

#' Internal function to calculate EOFs and related components
#' @keywords internal
get_eofs <- function(dat, k, rotate = FALSE) {
  # Save original dimension info
  orig_dims <- stars::st_dimensions(dat)
  orig_crs <- sf::st_crs(dat)
  times <- stars::st_get_dimension_values(dat, "time")
  dims <- dim(dat)

  pc_names <- names0(k, 'PC')

  # Handle spatial dimensions differently for geometry vs raster
  if (has_geometry_dimension(dat)) {
    # For sf geometry, spatial dimension is already a single dimension
    n_pixels <- dims[[1]]  # geometry dimension size
    space_dat <- dat %>%
      aperm(c("time", "geometry")) # Move spatial dimension to last
  } else {
    # For raster, combine x,y dimensions into a single spatial dimension
    n_pixels <- dims[[1]] * dims[[2]]
    space_dat <- dat %>%
      stars::st_redimension(c(location = n_pixels, time = length(times))) %>%
      aperm(c("time", "location")) # Move spatial dimension to last
  }

  # Get valid pixel indices (non-NA)
  valid_pixels <- which(!apply(space_dat[[1]], 2, anyNA)) # only works on first attribute

  # Validate k value
  max_k <- min(length(times) - 1, length(valid_pixels))
  check_k_valid(k, max_k)

  # Perform PCA on valid pixels, without centering since data are anomalies
  pca_result <- prcomp(space_dat[[1]][, valid_pixels], center = FALSE)

  # Keep everything as matrices initially for synchronized operations
  loadings_matrix <- pca_result$rotation[, 1:k, drop = FALSE]
  scores_matrix <- pca_result$x[, 1:k, drop = FALSE]
  sdev_vector <- pca_result$sdev[1:k]

  rotation_matrix <- NULL

  # Handle rotation with synchronized matrix operations
  if (rotate && k > 1) {
    rotation_result <- rotate_pca_components(loadings_matrix, scores_matrix, sdev_vector, k)
    loadings <- rotation_result$loadings
    amplitudes <- rotation_result$scores
    final_sdev <- rotation_result$sdev
    rotation_matrix <- rotation_result$rotation_matrix
  } else {
    loadings <- loadings_matrix
    amplitudes <- scores_matrix
    final_sdev <- sdev_vector
  }

  # Create EOF spatial patterns
  full_patterns <- array(NA, dim = c(k, n_pixels))
  full_patterns[, valid_pixels] <- t(loadings)

  if (has_geometry_dimension(dat)) {
    # For sf geometry, create patterns directly
    spatial_patterns <- dat[,,1:k] |>
      stars::st_set_dimensions('time', values = pc_names, names = 'PC') %>%
      setNames("weight")

    # Fill in the pattern values (geometry x PC)
    spatial_patterns[[1]] <- t(full_patterns)  # transpose to get geometry x PC
  } else {
    # For raster, reshape back to x,y,PC structure
    spatial_patterns <-  dat[,,,1:k] |>
      stars::st_set_dimensions('time', values = pc_names, names = 'PC') %>%
      setNames("weight")

    spatial_patterns[1] <- full_patterns |>
      array(dim = c(k, dims[[1]], dims[[2]])) %>%
      aperm(c(2, 3, 1))
  }

  # Format amplitudes (already rotated and reordered if needed)
  # Note: amplitudes from prcomp$x are already scaled by sdev internally
  # We should NOT divide by sdev again
  colnames(amplitudes) <- pc_names
  amplitudes <- amplitudes %>%
    as_tibble() %>%
    mutate(time = times, .before = 1)

  # Calculate eigenvalues with proper reordering for rotated cases
  n <- length(pca_result$sdev)
  if (rotate && k > 1) {
    # For rotated case: reorder the first k sdevs, keep rest unchanged
    reordered_sdev <- c(final_sdev, pca_result$sdev[(k+1):n])
  } else {
    reordered_sdev <- pca_result$sdev
  }

  eigenvalues <- tibble(
    PC = seq_len(n),
    std.dev = reordered_sdev,
    eigenvalues = std.dev^2,
    percent = eigenvalues / sum(eigenvalues),
    cumulative = cumsum(percent),
    error = sqrt(2/n),
    low = eigenvalues * (1 - error) * 100 / sum(eigenvalues),
    hi = eigenvalues * (1 + error) * 100 / sum(eigenvalues),
    cumvar_line = hi + 0.02 * max(hi)
  )

  list(
    spatial_patterns = spatial_patterns,
    amplitudes = amplitudes,
    pca = pca_result,
    eigenvalues = eigenvalues,
    rotation_matrix = rotation_matrix,
    valid_pixels = valid_pixels
  )
}

#' Check if stars object has sf geometry
#' @param dat A stars object
#' @return Logical indicating if object has geometry dimension
#' @keywords internal
has_geometry_dimension <- function(dat) {
  "geometry" %in% names(stars::st_dimensions(dat))
}

#' Get spatial dimension names for stars object
#' @param dat A stars object
#' @return Character vector of spatial dimension names
#' @keywords internal
get_spatial_dimensions <- function(dat) {
  if (has_geometry_dimension(dat)) {
    return("geometry")
  } else {
    # For raster, look for common spatial dimension name patterns
    dims <- stars::st_dimensions(dat)
    dim_names <- names(dims)

    # Common x-y dimension name patterns
    x_patterns <- c("x", "lon", "longitude")
    y_patterns <- c("y", "lat", "latitude")

    x_dim <- dim_names[tolower(dim_names) %in% x_patterns][1]
    y_dim <- dim_names[tolower(dim_names) %in% y_patterns][1]

    if (!is.na(x_dim) && !is.na(y_dim)) {
      return(c(x_dim, y_dim))
    } else {
      # Fallback: look for dimensions with spatial reference systems
      spatial_dims <- sapply(dims, function(d) !all(is.na(d$refsys)))
      spatial_names <- names(dims)[spatial_dims]

      # Filter out time and other non-spatial dimensions
      spatial_names <- spatial_names[!spatial_names %in% c("time", "band", "PC", "statistic")]

      if (length(spatial_names) >= 2) {
        return(spatial_names[1:2])
      } else if (length(spatial_names) == 1) {
        return(spatial_names)
      } else {
        # Final fallback: first two dimensions that aren't time
        non_time_dims <- dim_names[!dim_names %in% c("time", "band", "PC", "statistic")]
        return(non_time_dims[1:min(2, length(non_time_dims))])
      }
    }
  }
}

#' Compute area-based weights for spatial data
#'
#' Calculates area weights for spatial data using st_area(). Works uniformly
#' for raster grids, irregular geometries, and different coordinate systems.
#'
#' @param dat A stars object with spatial dimensions
#' @return Numeric weights (sqrt of normalized areas)
#' @export
area_weights <- function(dat) {
  if (has_geometry_dimension(dat)) {
    # For sf geometry, extract the geometry and compute areas directly
    geom <- sf::st_geometry(dat)
    areas <- sf::st_area(geom)
    area_values <- as.numeric(areas)
  } else {
    # For raster grids, use st_area on the stars object
    areas <- sf::st_area(dat)
    area_values <- as.numeric(areas[[1]])
  }

  # Normalize by mean area and take sqrt
  # sqrt so the covariance matrix is weighted by area
  sqrt(area_values / mean(area_values, na.rm = TRUE))
}

# from tidymodels/recipes
names0 <- function(num, prefix = "PC") {
  if (num < 1) {
    rlang::abort("`k` should be > 0.")
  }
  ind <- format(seq_len(num))
  ind <- gsub(" ", "0", ind)
  paste0(prefix, ind)
}

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
  cli::cli_text("Latitude weighted: {.field {x$weight}}")

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
  cli::cli_text("Latitude weighted: {.field {x$weight}}")

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

#' Plot method for patterns objects
#'
#' @param x A patterns object
#' @param type Type of plot: "combined" (default), "eofs", or "amplitudes"
#' @param scaled For EOFs: show correlations (TRUE) or raw loadings (FALSE)
#' @param rawdata Optional raw data for correlation calculation when scaled = TRUE
#' @param scale For amplitudes: scaling method ("standardized", "variance", "raw")
#' @param scale_y For amplitudes: y-axis scaling ("fixed" or "free")
#' @param events For amplitudes: optional dates to mark with vertical lines
#' @param ... Additional arguments (currently unused)
#' @return A ggplot2 object or patchwork object for combined plots
#' @export
plot.patterns <- function(x,
                          type = "combined",
                          scaled = FALSE,
                          rawdata = NULL,
                          scale = c("standardized", "variance", "raw"),
                          scale_y = c("fixed", "free"),
                          events = NULL,
                          ...) {

  type <- match.arg(type, c("combined", "eofs", "amplitudes"))
  scale <- match.arg(scale)
  scale_y <- match.arg(scale_y)

  if(type == "combined") {
    # Smart combined plot: EOFs above, PCs below with KISS layout logic
    if(!requireNamespace("patchwork", quietly = TRUE)) {
      warning("patchwork package needed for combined plots. Install with: install.packages('patchwork')\nShowing EOFs only.")
      return(.plot_eofs_internal(x, scaled = scaled, rawdata = rawdata))
    } else {
      # Simple smart layout - covers 80% of use cases perfectly
      layout <- .simple_smart_layout(x$k)

      # Create EOF plot with smart layout
      p_eofs <- .plot_eofs_internal(x, scaled = scaled, rawdata = rawdata, layout = layout$eof)

      # Create amplitude plot with matching layout
      p_amps <- .plot_amplitudes_internal(x, scale = scale, scale_y = scale_y, events = events, layout = layout$pc)

      # Combine with patchwork using top-over-bottom layout and alignment tricks
      return(p_eofs / p_amps +
             patchwork::plot_layout(
               heights = layout$heights,
               guides = "collect",
               axis_titles = "collect"
             ) +
             patchwork::plot_annotation(
               title = glue("EOF Analysis: {glue_collapse(x$names, sep = ', ')}"),
               subtitle = glue("k = {x$k} modes | {ifelse(x$scaled, 'scaled', 'unscaled')} | {ifelse(x$rotate, 'rotated', 'unrotated')}"),
               theme = theme(plot.title = element_text(size = 14, face = "bold"))
             ) &
             theme(plot.margin = margin(5, 5, 5, 5)))
    }
  } else if(type == "eofs") {
    return(.plot_eofs_internal(x, scaled = scaled, rawdata = rawdata))
  } else if(type == "amplitudes") {
    return(.plot_amplitudes_internal(x, scale = scale, scale_y = scale_y, events = events))
  }
}

#' Simple smart layout function - KISS principle
#' @keywords internal
.simple_smart_layout <- function(k) {
  if(k <= 3) {
    # Single row: works for most cases, clean alignment
    list(eof = list(nrow = 1), pc = list(nrow = 1), heights = c(5, 1))
  } else {
    # k >= 4: Use square-ish grid
    ncol <- ceiling(sqrt(k))
    list(eof = list(ncol = ncol), pc = list(ncol = ncol), heights = c(5, 1))
  }
}

#' Internal function for EOF plotting
#' @keywords internal
.plot_eofs_internal <- function(x, scaled = FALSE, rawdata = NULL, layout = NULL) {
  # Use provided layout or default to old behavior
  if(is.null(layout)) {
    facet_args <- list(nrow = 1)
  } else {
    facet_args <- layout
  }

  if(scaled) {
    if(is.null(rawdata)) {
      rlang::abort("rawdata must be provided when scaled = TRUE for correlation calculation", class = "tidyEOF_missing_rawdata")
    }
    ggplot() +
      geom_stars(data = get_correlation(rawdata, x)) +
      do.call(facet_wrap, c(list(~PC), facet_args)) +
      scale_fill_distiller(palette = 'RdBu', na.value = NA, limits = c(-1, 1)) +
      coord_sf() +
      theme_void() +
      theme(legend.position = "right") +
      labs(fill = "Correlation")
  } else {
    ggplot() +
      geom_stars(data = x$eofs) +
      do.call(facet_wrap, c(list(~PC), facet_args)) +
      scico::scale_fill_scico(palette = 'vik', midpoint = 0, na.value = NA) +
      coord_sf() +
      theme_void() +
      theme(legend.position = "right") +
      labs(fill = "Loading")
  }
}

#' Internal function for amplitude plotting
#' @keywords internal
.plot_amplitudes_internal <- function(x, scale = "standardized", scale_y = "fixed", events = NULL, layout = NULL) {

  # Use provided layout or default to old behavior
  if(is.null(layout)) {
    facet_args <- list(nrow = 1)
  } else {
    facet_args <- layout
  }

  # Calculate scaling factors based on method
  if(scale == "variance") {
    # Multiply by sqrt(eigenvalue) to show variance contribution
    eigs <- x$eigenvalues %>%
      select(PC, std.dev) %>%
      mutate(PC = paste0("PC", PC))

    amps <- x$amplitudes %>%
      pivot_longer(-time, names_to = 'PC', values_to = 'amplitude') %>%
      left_join(eigs, by = 'PC') %>%
      mutate(amplitude = amplitude * std.dev)
  } else if(scale == "raw") {
    # Use EOF loadings to convert back to original units
    eigs <- split(x$eofs) %>%
      as_tibble() %>%
      # get the sqrt of sum of squared loadings (works even if rotated)
      summarise(across(starts_with('PC'), ~sqrt(sum(.x^2, na.rm = TRUE)))) %>%
      pivot_longer(everything(), names_to = 'PC', values_to = 'scaling_factor')

    amps <- x$amplitudes %>%
      pivot_longer(-time, names_to = 'PC', values_to = 'amplitude') %>%
      left_join(eigs, by = 'PC') %>%
      mutate(amplitude = amplitude * scaling_factor)
  } else {
    # Default: standardized (already in x$amplitudes)
    amps <- x$amplitudes %>%
      pivot_longer(-time, names_to = 'PC', values_to = 'amplitude')
  }

  p <- ggplot(amps, aes(time, amplitude)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
    labs(x = "Time", y = 'Amplitude') +
    theme_bw()

  # Add faceting with smart layout and appropriate scales
  if(scale_y == "free") {
    p <- p + do.call(facet_wrap, c(list(~PC, scales = "free_y"), facet_args))
  } else {
    p <- p + do.call(facet_wrap, c(list(~PC), facet_args))
  }

  # Add event markers if provided
  if(!is.null(events)) {
    p <- p + geom_vline(xintercept = events, color = 'red', linetype = 2, alpha = 0.7)
  }

  # Force panel heights to match EOF panels
  p
}


#' Compare two patterns objects
#'
#' @param x First patterns object
#' @param y Second patterns object
#' @param method Comparison method
#' @return Comparison metrics
#' @export
compare_patterns <- function(x, y,
                             method = c("congruence", "procrustes", "correlation")) {

  method <- match.arg(method)

  # Ensure patterns have same spatial domain
  validate_spatial_compatibility(x, y)

  # Handle scaling based on method
  x_patterns <- prepare_patterns_for_comparison(x, method)
  y_patterns <- prepare_patterns_for_comparison(y, method)

  switch(method,
         "congruence" = {
           # Factor congruence coefficients
           # Use unrotated patterns for comparison
           pattern_congruence(x_patterns, y_patterns)
         },
         "procrustes" = {
           # Procrustes rotation to align patterns
           # Returns rotation matrix and congruence
           procrustes_align(x_patterns, y_patterns)
         },
         "correlation" = {
           # Temporal correlation of amplitudes
           # Handle potential time mismatches
           temporal_correlation(x_patterns, y_patterns)
         }
  )
}

# Helper functions for standardization and comparison

#' Standardize patterns for comparison or plotting
#' @keywords internal
standardize_patterns <- function(patterns, method) {
  switch(method,
         "unit_variance" = {
           # Scale to unit variance while preserving relative signs
           scale_patterns_unit_variance(patterns)
         },
         "max_amplitude" = {
           # Scale by maximum amplitude
           scale_patterns_max(patterns)
         },
         "none" = patterns
  )
}

#' Calculate symmetric limits for plotting
#' @keywords internal
symmetric_limits <- function(x) {
  max_abs <- max(abs(x$weight), na.rm = TRUE)
  c(-max_abs, max_abs)
}

#' Validate spatial compatibility of patterns
#' @keywords internal
validate_spatial_compatibility <- function(x, y) {
  # Check CRS, resolution, extent
  if(!identical(stars::st_dimensions(x$eofs)[c("x", "y")],
                stars::st_dimensions(y$eofs)[c("x", "y")])) {
    rlang::abort("Patterns must have identical spatial dimensions for comparison", class = "tidyEOF_incompatible_patterns")
  }
}

#' Scale patterns to unit variance while preserving relative signs
#'
#' @param patterns A stars object containing EOF patterns
#' @return A stars object with scaled patterns
#' @keywords internal
scale_patterns_unit_variance <- function(patterns) {
  # Convert to array for calculations
  pat_array <- stars::st_as_stars(patterns)

  # Calculate variance of each pattern (ignoring NAs)
  pattern_vars <- apply(pat_array, 3, function(x) sqrt(mean(x^2, na.rm = TRUE)))

  # Scale each pattern by its standard deviation while preserving sign
  scaled <- sweep(pat_array, 3, pattern_vars, "/")

  # Preserve attributes from original
  stars::st_dimensions(scaled) <- stars::st_dimensions(patterns)
  names(scaled) <- names(patterns)

  scaled
}

#' Scale patterns by their maximum absolute amplitude
#'
#' @param patterns A stars object containing EOF patterns
#' @return A stars object with scaled patterns
#' @keywords internal
scale_patterns_max <- function(patterns) {
  # Convert to array for calculations
  pat_array <- stars::st_as_stars(patterns)

  # Find maximum absolute value for each pattern
  pattern_max <- apply(pat_array, 3, function(x) max(abs(x), na.rm = TRUE))

  # Scale each pattern by its maximum
  scaled <- sweep(pat_array, 3, pattern_max, "/")

  # Preserve attributes from original
  stars::st_dimensions(scaled) <- stars::st_dimensions(patterns)
  names(scaled) <- names(patterns)

  scaled
}

#' Scale a complete patterns object
#'
#' @param patterns A patterns object
#' @param method Scaling method ("unit_variance" or "max_amplitude")
#' @param scale_pcs Logical, whether to also scale the PC time series (default TRUE)
#' @return A patterns object with scaled components
#' @export
scale_patterns <- function(patterns,
                           method = c("unit_variance", "max_amplitude"),
                           scale_pcs = TRUE) {
  method <- match.arg(method)

  # Get scaling factors based on method
  if (method == "unit_variance") {
    scaled_eofs <- scale_patterns_unit_variance(patterns$eofs)
    scaling_factors <- apply(stars::st_as_stars(patterns$eofs),
                             3,
                             function(x) sqrt(mean(x^2, na.rm = TRUE)))
  } else {
    scaled_eofs <- scale_patterns_max(patterns$eofs)
    scaling_factors <- apply(stars::st_as_stars(patterns$eofs),
                             3,
                             function(x) max(abs(x), na.rm = TRUE))
  }

  # Scale PCs inversely if requested to preserve reconstruction
  if (scale_pcs) {
    pc_cols <- grep("^PC", names(patterns$amplitudes), value = TRUE)
    scaled_amps <- patterns$amplitudes
    scaled_amps[pc_cols] <- sweep(scaled_amps[pc_cols], 2, scaling_factors, "*")
  } else {
    scaled_amps <- patterns$amplitudes
  }

  # Create new patterns object with scaled components
  patterns$eofs <- scaled_eofs
  patterns$amplitudes <- scaled_amps
  patterns$scaling_method <- method
  patterns$scaling_factors <- scaling_factors

  patterns
}

#' Generate symmetric limits for pattern plotting
#'
#' @param patterns A patterns object or stars object
#' @param buffer Fraction to expand limits by (default 0.05)
#' @return A numeric vector of length 2 with symmetric limits
#' @keywords internal
symmetric_limits <- function(patterns, buffer = 0.05) {
  if (inherits(patterns, "patterns")) {
    patterns <- patterns$eofs
  }

  # Find maximum absolute value
  max_abs <- max(abs(stars::st_as_stars(patterns)), na.rm = TRUE)

  # Add buffer and make symmetric
  lim <- max_abs * (1 + buffer)
  c(-lim, lim)
}

#' Scale patterns to unit variance while preserving relative signs
#'
#' @param patterns A stars object containing EOF patterns
#' @return A stars object with scaled patterns
#' @keywords internal
scale_patterns_unit_variance <- function(patterns) {
  # Get the array from the stars object
  pat_array <- patterns[[1]]

  # Calculate variance of each pattern (ignoring NAs)
  pattern_vars <- apply(pat_array, 3, function(x) sqrt(mean(x^2, na.rm = TRUE)))

  # Scale each pattern by its standard deviation while preserving sign
  scaled <- sweep(pat_array, 3, pattern_vars, "/")

  # Create new stars object with scaled values
  patterns[[1]] <- scaled
  patterns
}

#' Scale patterns by their maximum absolute amplitude
#'
#' @param patterns A stars object containing EOF patterns
#' @return A stars object with scaled patterns
#' @keywords internal
scale_patterns_max <- function(patterns) {
  # Get the array from the stars object
  pat_array <- patterns[[1]]

  # Find maximum absolute value for each pattern
  pattern_max <- apply(pat_array, 3, function(x) max(abs(x), na.rm = TRUE))

  # Scale each pattern by its maximum
  scaled <- sweep(pat_array, 3, pattern_max, "/")

  # Create new stars object with scaled values
  patterns[[1]] <- scaled
  patterns
}

#' Generate symmetric limits for pattern plotting
#'
#' @param patterns A patterns object or stars object
#' @param buffer Fraction to expand limits by (default 0.05)
#' @return A numeric vector of length 2 with symmetric limits
#' @keywords internal
symmetric_limits <- function(patterns, buffer = 0.05) {
  if (inherits(patterns, "patterns")) {
    patterns <- patterns$eofs
  }

  # Get array and find maximum absolute value
  max_abs <- max(abs(patterns[[1]]), na.rm = TRUE)

  # Add buffer and make symmetric
  lim <- max_abs * (1 + buffer)
  c(-lim, lim)
}