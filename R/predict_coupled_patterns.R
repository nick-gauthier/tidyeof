#' Predict Method for Coupled Patterns
#'
#' This function makes predictions using a coupled_patterns object created by couple_patterns().
#' It applies the learned CCA relationship to new predictor data to predict response patterns.
#'
#' @param object A coupled_patterns object from couple_patterns()
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
#' coupled <- couple_patterns(pred_patterns, resp_patterns, k = 3)
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
    cli::cli_abort("object must be a coupled_patterns object from couple_patterns()",
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
  reconstructed <- reconstruct_field(
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


  # Apply CCA transformation:
  # 1. Transform predictors to canonical variables
  # 2. Apply canonical correlations
  # 3. Transform back to response space
  canonical_predictors <- pred_matrix %*% cca_result$xcoef[, 1:k, drop = FALSE]
  canonical_responses <- canonical_predictors %*% diag(cca_result$cor[1:k], nrow = k)

  # Transform canonical responses back to PC space
  # Use generalized inverse (ginv) to handle case where k < number of response PCs
  ycoef_subset <- cca_result$ycoef[, 1:k, drop = FALSE]
  response_amplitudes <- canonical_responses %*% MASS::ginv(ycoef_subset)

  # Convert back to tibble with proper column names
  n_response_pcs <- ncol(response_amplitudes)
  pc_names <- paste0("PC", 1:n_response_pcs)

  result <- response_amplitudes %>%
    as_tibble() %>%
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

  canonical_vars <- amp_matrix %*% coef_matrix

  # Return as tibble
  canonical_names <- paste0("CV", 1:k)
  result <- canonical_vars %>%
    as_tibble() %>%
    setNames(canonical_names) %>%
    mutate(time = times, .before = 1)

  return(result)
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

  result <- data.frame(
    mode = 1:k,
    correlation = correlations,
    correlation_squared = correlations^2,
    variance_explained = correlations^2 / sum(correlations^2) * 100
  )

  return(result)
}