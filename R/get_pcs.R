#' Calculate principal components from spatiotemporal data
#'
#' @param dat A stars object containing spatiotemporal data
#' @param scale Logical, whether to scale the data before computing anomalies
#' @param clim Optional climatology object from get_climatology()
#' @param monthly Logical, whether to use monthly climatology
#'
#' @return A prcomp object containing the PCA results
#' @export
get_pcs <- function(dat, scale = FALSE, clim = NULL, monthly = FALSE, weight = TRUE) {
  # weight by sqrt cosine latitude, in radians
  if(weight) dat <- dat * lat_weights(dat)

  # Get anomalies (centers and optionally scales the data)
  anomalies <- get_anomalies(dat,
                             clim = clim,
                             scale = scale,
                             monthly = monthly)

  # Convert to matrix format for PCA
  data_matrix <- anomalies %>%
    split('time') %>% # split along the time dimension
    setNames(st_get_dimension_values(dat, 'time')) %>%
    as_tibble() %>%
    dplyr::select(-c(x,y)) %>%
    na.omit() %>%
    t() # transpose to space-as-columns format

  # Perform PCA without additional centering since data are anomalies
  prcomp(data_matrix, center = FALSE)
}

#' @export
lat_weights <- function(dat) {
  lats <- st_dim_to_attr(dat, which = 2) # get the y coordinates -- this is brittle if not in x, y, time order
    # convert to radians then apply cosine weighting
    # sqrt so the covariance matrix is weighted by cosine latitude
  sqrt(cos(lats * pi / 180))
}

#' @export
get_eigenvalues <- function(pca){
  n <- length(pca$sdev)

  pca %>%
    broom::tidy(matrix = 'pcs') %>%
    mutate(eigenvalues = std.dev ^ 2,
           error = sqrt(2 / n),
           low =  eigenvalues * (1 - error) * 100 / sum(eigenvalues),
           hi = eigenvalues * (1 + error) * 100 / sum(eigenvalues),
           cumvar_line = hi + 0.02 * max(hi))
}
