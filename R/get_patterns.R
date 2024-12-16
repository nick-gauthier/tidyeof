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
  orig_dims <- st_dimensions(dat)
  orig_crs <- st_crs(dat)
  times <- st_get_dimension_values(dat, "time")
  dims <- dim(dat)
  n_pixels <- dims[[1]] * dims[[2]]

  pc_names <- names0(k, 'PC')

  # Combine x,y dimensions into a single spatial dimension
  space_dat <- dat %>%
    st_redimension(c(location = n_pixels, time = length(times))) %>%
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
    st_set_dimensions('time', values = pc_names, names = 'PC') %>%
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
  lats <- st_dim_to_attr(dat, which = 2) # get the y coordinates -- this is brittle if not in x, y, time order
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

#' Plot method for patterns objects
#'
#' @param x A patterns object
#' @param type Character string indicating type of plot ("eofs", "amplitudes", or "scree")
#' @param ... Additional arguments passed to plotting functions
#' @export
plot.patterns <- function(x, type = "eofs", ...) {
  switch(
    type,
    "eofs" = plot_eofs(x, ...),
    "amplitudes" = plot_amps(x, ...),
    "scree" = plot_scree(x$eigenvalues, k = x$k, ...),
    stop("Invalid plot type. Must be one of 'eofs', 'amplitudes', or 'scree'")
  )
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