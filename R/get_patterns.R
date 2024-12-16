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
  dims <- dim(anom)
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
  amplitudes <- if (rotate && k > 1) {
    pca_result$x[, 1:k] %*% rotation_matrix # do these need to be reordered in the case of reofs or is it preserved by new rotation matrix?
  } else {
    pca_result$x[, 1:k]
  } %>%
    sweep(2, pca_result$sdev[1:k], "/") %>% # does this happen before or after rotation?
    as_tibble() %>%
    setNames(pc_names) %>%
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
    eigenvalues = eigenvalues
    #rotation_matrix = rotation_matrix,
    #valid_pixels = valid_pixels
  )
}

#' @export
lat_weights <- function(dat) {
  lats <- st_dim_to_attr(dat, which = 2) # get the y coordinates -- this is brittle if not in x, y, time order
  # convert to radians then apply cosine weighting
  # sqrt so the covariance matrix is weighted by cosine latitude
  sqrt(cos(lats * pi / 180))
}

#' @export
print.patterns <- function(obj) {
  print(paste0('A `pattern` object with k = ', obj$k, ', scale = ', obj$scale,
               ', monthly = ', obj$monthly, ', and rotate = ', obj$rotate))
  print(obj$eofs)
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
