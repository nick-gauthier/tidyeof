#' Couple Pattern Relationships Using CCA
#'
#' This function couples predictor and response patterns using Canonical Correlation Analysis (CCA)
#' as the primary method. CCA finds linear combinations of predictor and response patterns that
#' maximize correlation between them.
#'
#' @param predictor_patterns A patterns object containing predictor patterns (e.g., from get_patterns())
#' @param response_patterns A patterns object containing response patterns (e.g., from get_patterns())
#' @param k Number of CCA modes to retain. If NULL, uses min(ncol(predictor), ncol(response))
#' @param method Coupling method. Currently only "cca" is supported
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

  # Validate inputs and get common times
  common_times <- if (validate) {
    validate_patterns_compatibility(predictor_patterns, response_patterns)
  } else {
    # If not validating, still need to find common times
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
    # Preserve Date/POSIXct class
    pred_times[pred_times %in% resp_times]
  }

  if (method != "cca") {
    cli::cli_abort("Only 'cca' method is supported. For legacy methods, see vignettes/legacy_pattern_coupling.qmd",
                   class = "tidyEOF_unsupported_method")
  }

  # Extract amplitude matrices filtered to common times
  pred_amps <- extract_amplitudes_matrix(predictor_patterns, common_times)
  resp_amps <- extract_amplitudes_matrix(response_patterns, common_times)

  # Notify user if filtering occurred

  pred_n <- if ('patterns' %in% class(predictor_patterns)) nrow(predictor_patterns$amplitudes) else nrow(predictor_patterns)
  resp_n <- if ('patterns' %in% class(response_patterns)) nrow(response_patterns$amplitudes) else nrow(response_patterns)

  if (length(common_times) < pred_n || length(common_times) < resp_n) {
    cli::cli_inform("Filtered to {length(common_times)} common time steps (predictor had {pred_n}, response had {resp_n}).")
  }

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
extract_amplitudes_matrix <- function(patterns, times = NULL) {
  if ('patterns' %in% class(patterns)) {
    amps <- patterns$amplitudes
  } else {
    amps <- patterns
  }

  # Filter to specified times if provided
  if (!is.null(times)) {
    # Use semi_join for robust time matching (handles Date/POSIXct better than %in%)
    times_df <- tibble::tibble(time = times)
    amps <- dplyr::semi_join(amps, times_df, by = "time")
  }

  amps %>%
    dplyr::select(-time) %>%
    as.matrix()
}

# Helper function to validate pattern compatibility and return common times
validate_patterns_compatibility <- function(predictor_patterns, response_patterns) {
  # Check that both are patterns objects or compatible
  if (!('patterns' %in% class(predictor_patterns)) && !is.data.frame(predictor_patterns)) {
    cli::cli_abort("predictor_patterns must be a patterns object or data frame with amplitudes",
                   class = "tidyEOF_invalid_input")
  }

  if (!('patterns' %in% class(response_patterns)) && !is.data.frame(response_patterns)) {
    cli::cli_abort("response_patterns must be a patterns object or data frame with amplitudes",
                   class = "tidyEOF_invalid_input")
  }

  # Get time vectors

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

  # Use match-based intersection to preserve Date/POSIXct class
  common_times <- pred_times[pred_times %in% resp_times]

  if (length(common_times) == 0) {
    cli::cli_abort("No common time steps found between predictor and response patterns",
                   class = "tidyEOF_no_common_times")
  }

  if (length(common_times) < 10) {
    cli::cli_warn("Only {length(common_times)} common time steps found. Consider using more data.")
  }

  # Return common times for filtering
  common_times
}

# Legacy GAM/PCR fitting functions - preserved for reference but no longer used
# See vignettes/legacy_pattern_coupling.qmd for historical context
#
# fit_model_legacy <- function(data_in) {
#   gam_formula <- data_in %>%
#     dplyr::select(-PC, -amplitude, -time) %>%
#     names() %>%
#     purrr::map( ~ paste0("s(", ., ", bs = 'cr', k = 3)")) %>%
#     paste(collapse = ' + ') %>%
#     paste('amplitude ~ ', .) %>%
#     as.formula()
#
#   model <- data_in %>%
#     dplyr::group_by(PC) %>%
#     tidyr::nest() %>%
#     dplyr::mutate(mod = purrr::map(data, ~ mgcv::gam(gam_formula, data = ., method = 'REML', select = TRUE))) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(r2 = purrr::map_dbl(mod, ~ summary(.)$r.sq)) %>%
#     dplyr::mutate(
#       data = purrr::map2(data, mod, modelr::add_predictions),
#       data = purrr::map2(data, mod, modelr::add_residuals)
#     )
#   return(model)
# }
#
# fit_pcr_legacy <- function(data_in) {
#   lm_formula <- data_in %>%
#     dplyr::select(-PC, -amplitude, -time) %>%
#     names() %>%
#     paste(collapse = ' + ') %>%
#     paste('amplitude ~ ', .) %>%
#     as.formula
#
#   data_in %>%
#     dplyr::group_by(PC) %>%
#     tidyr::nest() %>%
#     dplyr::mutate(mod = purrr::map(data, ~ eval(getCall(MuMIn::dredge(lm(lm_formula, data = ., na.action = "na.fail")), 1)))) %>%
#     dplyr::ungroup()
# }

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