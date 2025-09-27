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

  # Apply latitude weighting if requested
  if(weight) dat <- dat * lat_weights(dat) # do we do this before or after scaling?

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
    units = purrr::map(dat, purrr::possibly(units)),
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

#' Internal function to calculate EOFs and related components
#' @keywords internal
get_eofs <- function(dat, k, rotate = FALSE) {
  # Save original dimension info
  orig_dims <- stars::st_dimensions(dat)
  orig_crs <- sf::st_crs(dat)
  times <- stars::st_get_dimension_values(dat, "time")
  dims <- dim(dat)
  n_pixels <- dims[[1]] * dims[[2]]

  pc_names <- names0(k, 'PC')

  # Combine x,y dimensions into a single spatial dimension
  space_dat <- dat %>%
    stars::st_redimension(c(location = n_pixels, time = length(times))) %>%
    aperm(c("time", "location")) # Move spatial dimension to last

  # Get valid pixel indices (non-NA)
  valid_pixels <- which(!apply(space_dat[[1]], 2, anyNA)) # only works on first attribute

  # Perform PCA on valid pixels, without centering since data are anomalies
  pca_result <- prcomp(space_dat[[1]][, valid_pixels], center = FALSE)

  # Handle rotation if requested
  rotation_matrix <- NULL
  loadings <- pca_result$rotation[, 1:k, drop = FALSE]

  if (rotate && k > 1) {
    rot <- varimax(loadings %*% diag(pca_result$sdev, k, k)) # scale by sdev for more robust rotation
    loadings <- unclass(rot$loadings) # what happens if we don't unclass?
    rotation_matrix <- rot$rotmat # save the rotation matrix to use later on the amplitudes

    # Reorder by explained variance, after psych::principal()
    ev.rotated <- diag(t(loadings) %*% loadings)
    ev.order <- order(ev.rotated, decreasing = TRUE)
    loadings <- loadings[, ev.order]
    rotation_matrix <- rotation_matrix[, ev.order]
  }

  # Create EOF spatial patterns
  full_patterns <- array(NA, dim = c(k, n_pixels))
  full_patterns[, valid_pixels] <- t(loadings)

  spatial_patterns <-  dat[,,,1:k] |>
    stars::st_set_dimensions('time', values = pc_names, names = 'PC') %>%
    setNames("weight")

  spatial_patterns[1] <- full_patterns |>
    array(dim = c(k, dims[[1]], dims[[2]])) %>%
    aperm(c(2, 3, 1))

  # Calculate amplitudes
  amplitudes <- pca_result$x[, 1:k]
  if (rotate && k > 1) amplitudes <- amplitudes %*% rotation_matrix # do these need to be reordered in the case of reofs or is it preserved by new rotation matrix?
  colnames(amplitudes) <- pc_names
  amplitudes <- amplitudes %>%
    sweep(2, pca_result$sdev[1:k], "/") %>% # does this happen before or after rotation? is this doing anything?
    as_tibble() %>%
    mutate(time = times, .before = 1)

  # Calculate eigenvalues
  n <- length(pca_result$sdev)
  eigenvalues <- tibble(
    PC = seq_len(n),
    std.dev = pca_result$sdev,
    eigenvalues = std.dev^2,
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

#' @export
lat_weights <- function(dat) {
  lats <- stars::st_dim_to_attr(dat, which = 2) # get the y coordinates -- this is brittle if not in x, y, time order
  # convert to radians then apply cosine weighting
  # sqrt so the covariance matrix is weighted by cosine latitude
  sqrt(cos(lats * pi / 180))
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
  cat(sprintf("A patterns object with %d modes\n", x$k))
  cat(sprintf("  Processing options:\n"))
  cat(sprintf("    - Scale: %s\n", x$scaled))
  cat(sprintf("    - Monthly: %s\n", x$monthly))
  cat(sprintf("    - Rotated: %s\n", x$rotate))
  cat(sprintf("    - Latitude weighted: %s\n", x$weight))

  cat("\nEigenvalues:\n")
  top_eigs <- head(x$eigenvalues, x$k)
  print(round(top_eigs$eigenvalues / sum(x$eigenvalues$eigenvalues) * 100, 1))

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
  cat("Summary of patterns object:\n\n")
  cat(sprintf("Number of modes: %d\n", x$k))
  cat(sprintf("Time steps: %d\n", x$n_times))
  cat(sprintf("Total variance explained: %.1f%%\n", x$var_explained))

  if(!is.null(x$pattern_cors)) {
    cat("\nPattern correlations:\n")
    print(x$pattern_cors)
  }

  cat("\nProcessing options:\n")
  cat(sprintf("  Scale: %s\n", x$scaled))
  cat(sprintf("  Monthly: %s\n", x$monthly))
  cat(sprintf("  Rotated: %s\n", x$rotate))
  cat(sprintf("  Latitude weighted: %s\n", x$weight))

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
    if(any(is.na(i))) stop("Invalid PC names")
  }

  if(any(i > x$k)) stop("Index out of bounds")

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
    if(pc_num > x$k) stop("Invalid PC number")
    return(x$amplitudes[[paste0("PC", pc_num)]])
  }

  if(grepl("^EOF", name)) {
    # Extract EOF spatial pattern
    eof_num <- as.numeric(sub("^EOF", "", name))
    if(eof_num > x$k) stop("Invalid EOF number")
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
               title = paste("EOF Analysis:", x$names, collapse = ", "),
               subtitle = paste("k =", x$k, "modes |",
                              if(x$scaled) "scaled" else "unscaled", "|",
                              if(x$rotate) "rotated" else "unrotated"),
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
      stop("rawdata must be provided when scaled = TRUE")
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
  p# + theme(
    #panel.spacing = unit(0.1, "lines"),
    #strip.text = element_text(size = 10),
    #aspect.ratio = 0.6  # Force consistent aspect ratio for alignment
  #)
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
    stop("Patterns must have identical spatial dimensions")
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