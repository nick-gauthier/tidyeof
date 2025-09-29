#' Couple Pattern Relationships Using CCA
#'
#' This function couples predictor and response patterns using Canonical Correlation Analysis (CCA)
#' as the primary method. CCA finds linear combinations of predictor and response patterns that
#' maximize correlation between them.
#'
#' @param predictor_patterns A patterns object containing predictor patterns (e.g., from get_patterns())
#' @param response_patterns A patterns object containing response patterns (e.g., from get_patterns())
#' @param k Number of CCA modes to retain. If NULL, uses min(ncol(predictor), ncol(response))
#' @param method Coupling method. Currently only "cca" is supported in the refactored version
#' @param center Logical, whether to center the data before CCA (default: FALSE)
#' @param validate Logical, whether to validate input patterns compatibility
#'
#' @return A coupled_patterns object containing:
#'   \item{cca}{The CCA results from cancor()}
#'   \item{predictor_patterns}{The original predictor patterns}
#'   \item{response_patterns}{The original response patterns}
#'   \item{k}{Number of CCA modes retained}
#'   \item{method}{Coupling method used}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get patterns from your data
#' pred_patterns <- get_patterns(predictor_data, k = 5)
#' resp_patterns <- get_patterns(response_data, k = 5)
#'
#' # Couple the patterns
#' coupled <- couple_patterns(pred_patterns, resp_patterns, k = 3)
#'
#' # Make predictions
#' predictions <- predict(coupled, new_predictor_data)
#' }
couple_patterns <- function(predictor_patterns, response_patterns, k = NULL,
                           method = "cca", center = FALSE, validate = TRUE) {

  # Validate inputs
  if (validate) {
    validate_patterns_compatibility(predictor_patterns, response_patterns)
  }

  if (method != "cca") {
    stop("Only 'cca' method is supported in the refactored version. ",
         "For legacy methods, see vignettes/legacy_pattern_coupling.qmd")
  }

  # Extract amplitude matrices
  pred_amps <- extract_amplitudes_matrix(predictor_patterns)
  resp_amps <- extract_amplitudes_matrix(response_patterns)

  # Determine k if not specified
  if (is.null(k)) {
    k <- min(ncol(pred_amps), ncol(resp_amps))
  }

  # Validate k
  max_k <- min(ncol(pred_amps), ncol(resp_amps))
  if (k > max_k) {
    warning("k = ", k, " exceeds maximum possible (", max_k, "). Setting k = ", max_k)
    k <- max_k
  }

  # Perform CCA
  cca_result <- cancor(pred_amps, resp_amps, xcenter = center, ycenter = center)

  # Create coupled_patterns object
  coupled <- list(
    cca = cca_result,
    predictor_patterns = predictor_patterns,
    response_patterns = response_patterns,
    k = k,
    method = method,
    center = center
  )

  class(coupled) <- "coupled_patterns"
  return(coupled)
}

# Helper function to extract amplitude matrices from patterns objects
extract_amplitudes_matrix <- function(patterns) {
  if ('patterns' %in% class(patterns)) {
    amps <- patterns$amplitudes %>%
      dplyr::select(-time) %>%
      as.matrix()
  } else {
    # Assume it's already in the right format
    amps <- patterns %>%
      dplyr::select(-time) %>%
      as.matrix()
  }
  return(amps)
}

# Helper function to validate pattern compatibility
validate_patterns_compatibility <- function(predictor_patterns, response_patterns) {
  # Check that both are patterns objects or compatible
  if (!('patterns' %in% class(predictor_patterns)) && !is.data.frame(predictor_patterns)) {
    stop("predictor_patterns must be a patterns object or data frame with amplitudes")
  }

  if (!('patterns' %in% class(response_patterns)) && !is.data.frame(response_patterns)) {
    stop("response_patterns must be a patterns object or data frame with amplitudes")
  }

  # Check time dimensions match
  pred_times <- if ('patterns' %in% class(predictor_patterns)) {
    predictor_patterns$amplitudes$time
  } else {
    predictor_patterns$time
  }

  resp_times <- if ('patterns' %in% class(response_patterns)) {
    response_patterns$amplitudes$time
  } else {
    response_patterns$time
  }

  common_times <- intersect(pred_times, resp_times)
  if (length(common_times) == 0) {
    stop("No common time steps found between predictor and response patterns")
  }

  if (length(common_times) < 10) {
    warning("Only ", length(common_times), " common time steps found. Consider using more data.")
  }
}

# Legacy GAM fitting function - moved to legacy support
# This function is preserved for backward compatibility but not used in the main CCA workflow
fit_model_legacy <- function(data_in) {
  .Deprecated("couple_patterns",
              msg = "GAM-based coupling is deprecated. Use couple_patterns() with method='cca' instead.")

  gam_formula <- data_in %>%
    dplyr::select(-PC, -amplitude, -time) %>%
    names() %>%
    purrr::map( ~ paste0("s(", ., ", bs = 'cr', k = 3)")) %>%
    paste(collapse = ' + ') %>%
    paste('amplitude ~ ', .) %>%
    as.formula()

  model <- data_in %>%
    dplyr::group_by(PC) %>%
    tidyr::nest() %>%
    dplyr::mutate(mod = purrr::map(data, ~ mgcv::gam(gam_formula, data = ., method = 'REML', select = TRUE))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(r2 = purrr::map_dbl(mod, ~ summary(.)$r.sq)) %>%
    dplyr::mutate(
      data = purrr::map2(data, mod, modelr::add_predictions),
      data = purrr::map2(data, mod, modelr::add_residuals)
    )

  return(model)
}

# Legacy PCR fitting function - moved to legacy support
fit_pcr_legacy <- function(data_in) {
  .Deprecated("couple_patterns",
              msg = "PCR-based coupling is deprecated. Use couple_patterns() with method='cca' instead.")

  lm_formula <- data_in %>%
    dplyr::select(-PC, -amplitude, -time) %>%
    names() %>%
    paste(collapse = ' + ') %>%
    paste('amplitude ~ ', .) %>%
    as.formula

  data_in %>%
    dplyr::group_by(PC) %>%
    tidyr::nest() %>%
    dplyr::mutate(mod = purrr::map(data, ~ eval(getCall(MuMIn::dredge(lm(lm_formula, data = ., na.action = "na.fail")), 1)))) %>%
    dplyr::ungroup()
}

#' Print method for coupled_patterns
#' @param x A coupled_patterns object
#' @param ... Additional arguments (ignored)
#' @export
print.coupled_patterns <- function(x, ...) {
  cat("Coupled Patterns Object\n")
  cat("Method:", x$method, "\n")
  cat("CCA modes retained:", x$k, "\n")
  cat("Canonical correlations:", round(x$cca$cor[1:x$k], 3), "\n")
  cat("Predictor patterns:", ncol(extract_amplitudes_matrix(x$predictor_patterns)), "PCs\n")
  cat("Response patterns:", ncol(extract_amplitudes_matrix(x$response_patterns)), "PCs\n")
}

#' Summary method for coupled_patterns
#' @param object A coupled_patterns object
#' @param ... Additional arguments (ignored)
#' @export
summary.coupled_patterns <- function(object, ...) {
  cat("Coupled Patterns Summary\n")
  cat("======================\n")
  print(object)
  cat("\nCanonical Correlation Details:\n")
  for (i in 1:object$k) {
    cat("Mode", i, "- Correlation:", round(object$cca$cor[i], 4), "\n")
  }
}