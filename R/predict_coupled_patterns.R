#' Predict Method for Coupled Patterns
#'
#' This function makes predictions using a coupled_patterns object created by couple().
#' It applies the learned CCA relationship to new predictor data to predict response patterns.
#'
#' @param object A coupled_patterns object from couple()
#' @param newdata New predictor data (stars object) for making predictions
#' @param k Number of CCA modes to use for prediction. If NULL, uses all available modes
#' @param reconstruct Logical, whether to reconstruct the full spatial field (default: TRUE)
#' @param nonneg How to handle non-negative constraints: "auto", TRUE, FALSE, or "inherit"
#' @param ... Additional arguments (currently unused)
#'
#' @return If reconstruct=TRUE, returns a stars object with reconstructed spatial fields.
#'         If reconstruct=FALSE, returns a tibble with predicted amplitudes.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create coupled patterns
#' coupled <- couple(pred_patterns, resp_patterns, k = 3)
#'
#' # Make predictions on new data
#' predictions <- predict(coupled, new_predictor_data)
#'
#' # Just get predicted amplitudes without spatial reconstruction
#' amplitudes <- predict(coupled, new_predictor_data, reconstruct = FALSE)
#' }
predict.coupled_patterns <- function(object, newdata, k = NULL, reconstruct = TRUE,
                                   nonneg = "auto", ...) {


  # Validate inputs
  if (!inherits(object, "coupled_patterns")) {
    cli::cli_abort("object must be a coupled_patterns object from couple()",
                   class = "tidyEOF_invalid_input")
  }

  if (object$method != "cca") {
    cli::cli_abort("Only CCA-based coupled_patterns objects are supported",
                   class = "tidyEOF_unsupported_method")
  }

  # Use object's k if not specified
  if (is.null(k)) {
    k <- object$k
  }

  # Validate k
  if (k > object$k) {
    warning("Requested k (", k, ") exceeds available modes (", object$k, "). Using k = ", object$k)
    k <- object$k
  }

  # Project new data onto predictor patterns to get PC amplitudes
  new_amplitudes <- project_patterns(object$predictor_patterns, newdata)

  # Apply CCA transformation
  predicted_amplitudes <- apply_cca_prediction(
    new_amplitudes = new_amplitudes,
    cca_result = object$cca,
    k = k
  )

  if (!reconstruct) {
    return(predicted_amplitudes)
  }

  # Reconstruct spatial field
  reconstructed <- reconstruct(
    target_patterns = object$response_patterns,
    amplitudes = predicted_amplitudes,
    nonneg = nonneg
  )

  return(reconstructed)
}

#' Apply CCA Prediction Transform
#'
#' Internal function that applies the CCA transformation to predict response amplitudes
#' from predictor amplitudes.
#'
#' @param new_amplitudes Tibble with time and predictor PC amplitudes
#' @param cca_result CCA result object from cancor()
#' @param k Number of CCA modes to use
#'
#' @return Tibble with time and predicted response PC amplitudes
#'
#' @keywords internal
apply_cca_prediction <- function(new_amplitudes, cca_result, k) {

  # Extract times
  new_times <- new_amplitudes$time

  # Convert to matrix for CCA transformation
  pred_matrix <- new_amplitudes %>%
    dplyr::select(-time) %>%
    as.matrix()

  # Apply training centering if CCA was fit with centering
  if (!identical(cca_result$xcenter, FALSE)) {
    xcenter <- cca_result$xcenter
    if (!is.null(names(xcenter)) && !is.null(colnames(pred_matrix))) {
      xcenter <- xcenter[colnames(pred_matrix)]
    }
    pred_matrix <- sweep(pred_matrix, 2, xcenter, "-")
  }

  # Apply CCA transformation:
  # 1. Transform predictors to canonical variables
  # 2. Apply canonical correlations
  # 3. Transform back to response space
  canonical_predictors <- pred_matrix %*% cca_result$xcoef[, 1:k, drop = FALSE]
  canonical_responses <- canonical_predictors %*% diag(cca_result$cor[1:k], nrow = k)

  # Transform canonical responses back to PC space
  # Use pseudo-inverse of TRUNCATED ycoef, not subset of full inverse
  # (full inverse has cross-contributions from unused modes)
  ycoef_k <- cca_result$ycoef[, 1:k, drop = FALSE]
  response_amplitudes <- canonical_responses %*% MASS::ginv(ycoef_k)

  # Add back response centering if used during training
  if (!identical(cca_result$ycenter, FALSE)) {
    ycenter <- cca_result$ycenter
    if (!is.null(names(ycenter)) && !is.null(colnames(response_amplitudes))) {
      ycenter <- ycenter[colnames(response_amplitudes)]
    }
    response_amplitudes <- sweep(response_amplitudes, 2, ycenter, "+")
  }

  # Convert back to tibble with proper column names
  n_response_pcs <- ncol(response_amplitudes)
  pc_names <- paste0("PC", 1:n_response_pcs)

  result <- response_amplitudes %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(pc_names) %>%
    mutate(time = new_times, .before = 1)

  return(result)
}

#' Get Canonical Variables from Coupled Patterns
#'
#' Extract canonical variables from either predictor or response patterns
#'
#' @param object A coupled_patterns object
#' @param data Original data (patterns object or amplitudes tibble)
#' @param type Either "predictor" or "response"
#' @param k Number of canonical modes to extract
#'
#' @return Tibble with canonical variables
#'
#' @export
get_canonical_variables <- function(object, data, type = c("predictor", "response"), k = NULL) {

  type <- match.arg(type)

  if (is.null(k)) {
    k <- object$k
  }

  # Extract amplitudes
  if ("patterns" %in% class(data)) {
    amplitudes <- data$amplitudes
  } else {
    amplitudes <- data
  }

  # Get transformation coefficients
  if (type == "predictor") {
    coef_matrix <- object$cca$xcoef[, 1:k, drop = FALSE]
  } else {
    coef_matrix <- object$cca$ycoef[, 1:k, drop = FALSE]
  }

  # Apply transformation
  times <- amplitudes$time
  amp_matrix <- amplitudes %>%
    dplyr::select(-time) %>%
    as.matrix()

  # Apply training centering so canonical variables match cancor() inputs
  if (type == "predictor" && !identical(object$cca$xcenter, FALSE)) {
    xcenter <- object$cca$xcenter
    if (!is.null(names(xcenter)) && !is.null(colnames(amp_matrix))) {
      xcenter <- xcenter[colnames(amp_matrix)]
    }
    amp_matrix <- sweep(amp_matrix, 2, xcenter, "-")
  }
  if (type == "response" && !identical(object$cca$ycenter, FALSE)) {
    ycenter <- object$cca$ycenter
    if (!is.null(names(ycenter)) && !is.null(colnames(amp_matrix))) {
      ycenter <- ycenter[colnames(amp_matrix)]
    }
    amp_matrix <- sweep(amp_matrix, 2, ycenter, "-")
  }

  canonical_vars <- amp_matrix %*% coef_matrix

  # Return as tibble
  canonical_names <- paste0("CV", 1:k)
  result <- canonical_vars %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(canonical_names) %>%
    mutate(time = times, .before = 1)

  return(result)
}

#' Get Canonical Spatial Patterns from Coupled Patterns
#'
#' Computes the spatial patterns corresponding to each canonical mode by
#' taking linear combinations of the original EOFs weighted by CCA coefficients.
#' These are the spatial patterns that, when projected onto the data, yield
#' the canonical variates.
#'
#' @param object A coupled_patterns object
#' @param type Either "predictor" or "response"
#' @param k Number of canonical modes to extract (default: all available)
#'
#' @return A stars object with canonical spatial patterns (dimension "CV" instead of "PC")
#'
#' @export
#'
#' @examples
#' \dontrun{
#' coupled <- couple(pred_patterns, resp_patterns, k = 3)
#'
#' # Get canonical patterns for response side
#' resp_canonical <- get_canonical_patterns(coupled, type = "response")
#' plot(resp_canonical)
#'
#' # Compare to original EOFs
#' plot(coupled$response_patterns$eofs)
#' }
get_canonical_patterns <- function(object, type = c("predictor", "response"), k = NULL) {

  type <- match.arg(type)

  if (is.null(k)) {
    k <- object$k
  }

  # Get the patterns object and CCA coefficients

  if (type == "predictor") {
    patterns <- object$predictor_patterns
    coef_matrix <- object$cca$xcoef[, 1:k, drop = FALSE]
  } else {
    patterns <- object$response_patterns
    coef_matrix <- object$cca$ycoef[, 1:k, drop = FALSE]
  }

  # Extract EOF array (spatial dims + PC)
  eof_stars <- patterns$eofs
  eof_dims <- stars::st_dimensions(eof_stars)

  # Find spatial dimensions (everything except PC)
  spatial_dim_names <- setdiff(names(eof_dims), "PC")
  n_pcs <- length(stars::st_get_dimension_values(eof_stars, "PC"))

  # Get the underlying array
  eof_array <- eof_stars[[1]]

  # Determine array dimension order
  dim_names <- names(dim(eof_array))
  if (is.null(dim_names)) {
    # Assume standard order from stars: spatial dims first, PC last
    dim_names <- c(spatial_dim_names, "PC")
  }
  pc_dim <- which(dim_names == "PC")

  # Reshape to matrix: (spatial pixels) x (PCs)
  # Move PC dimension to last if not already
  if (pc_dim != length(dim(eof_array))) {
    perm_order <- c(setdiff(seq_along(dim(eof_array)), pc_dim), pc_dim)
    eof_array <- aperm(eof_array, perm_order)
  }

  spatial_shape <- dim(eof_array)[-length(dim(eof_array))]
  n_spatial <- prod(spatial_shape)
  eof_matrix <- matrix(eof_array, nrow = n_spatial, ncol = n_pcs)

  # Compute canonical patterns: (spatial) x (canonical modes)
  # canonical_pattern[i] = Σ_j EOF[j] × coef[j,i]
  canonical_matrix <- eof_matrix %*% coef_matrix

  # Reshape back to spatial array with CV dimension
  canonical_array <- array(canonical_matrix, dim = c(spatial_shape, k))

  # Build new stars object with CV dimension instead of PC
  new_dims <- eof_dims[spatial_dim_names]

  cv_dim <- list(
    from = 1L,
    to = k,
    offset = NA_real_,
    delta = NA_real_,
    refsys = NA_character_,
    point = FALSE,
    values = paste0("CV", 1:k)
  )
  class(cv_dim) <- "dimension"
  new_dims$CV <- cv_dim
  class(new_dims) <- "dimensions"

  result <- stars::st_as_stars(canonical_array, dimensions = new_dims)
  names(result) <- names(eof_stars)

  result
}

#' Get Canonical Correlations from Coupled Patterns
#'
#' Extract the canonical correlations and related statistics
#'
#' @param object A coupled_patterns object
#' @param k Number of modes to return (default: all available)
#'
#' @return Data frame with canonical correlation statistics
#'
#' @export
get_canonical_correlations <- function(object, k = NULL) {

  if (is.null(k)) {
    k <- object$k
  }

  correlations <- object$cca$cor[1:k]

  # TODO: variance_explained calculation is misleading. Squared canonical

  # correlations are not directly interpretable as variance explained like PCA
  # eigenvalues. For proper variance decomposition in CCA, need to compute
  # redundancy indices or use Stewart-Love indices. Fix or rename this column.
  result <- data.frame(
    mode = 1:k,
    correlation = correlations,
    correlation_squared = correlations^2,
    variance_explained = correlations^2 / sum(correlations^2) * 100
  )

  return(result)
}
