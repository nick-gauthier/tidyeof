# =============================================================================
# LEGACY WRAPPERS - COMMENTED OUT
# =============================================================================
# This file contains deprecated wrapper functions for backward compatibility.
# Kept for reference when replicating old analyses that used previous versions.
# None of this code is currently active in the package.
#
# If you need to use any of this, copy the relevant functions to your own script
# or uncomment selectively.
# =============================================================================

# #' Legacy Wrapper Functions for Backward Compatibility
# #'
# #' These functions provide backward compatibility for the original prediction methods.
# #' They are deprecated and will be removed in a future version. Use couple_patterns()
# #' with the new CCA-based architecture instead.
#
# #' Legacy CCA Prediction (Deprecated)
# #'
# #' @param preds Predictor patterns object
# #' @param obs Response patterns object
# #' @param newdata New data for prediction
# #' @param k Number of CCA modes
# #'
# #' @return Reconstructed field
# #'
# #' @export
# predict_cca <- function(preds, obs, newdata, k) {
#   .Deprecated("couple_patterns",
#               msg = paste("predict_cca() is deprecated.",
#                          "Use couple_patterns() followed by predict() instead.",
#                          "See vignettes/legacy_pattern_coupling.qmd for migration guide."))
#
#   # Use new architecture under the hood
#   coupled <- couple_patterns(preds, obs, k = k)
#   predict(coupled, newdata, k = k)
# }
#
# #' Legacy GAM Prediction (Deprecated)
# #'
# #' @param preds Predictor patterns object
# #' @param obs Response patterns object
# #' @param newdata New data for prediction
# #' @param k Number of modes (ignored in GAM, kept for compatibility)
# #'
# #' @return Reconstructed field
# #'
# #' @export
# predict_gam <- function(preds, obs, newdata, k) {
#   .Deprecated("couple_patterns",
#               msg = paste("predict_gam() is deprecated.",
#                          "GAM-based coupling is no longer supported in the main workflow.",
#                          "Use couple_patterns() with CCA instead, or see legacy_pattern_coupling.qmd"))
#
#   stop("GAM-based prediction is no longer supported in the refactored version. ",
#        "Use couple_patterns() with method='cca' instead, or refer to ",
#        "vignettes/legacy_pattern_coupling.qmd for the original implementation.")
# }
#
# #' Legacy PCR Prediction (Deprecated)
# #'
# #' @param preds Predictor patterns object
# #' @param obs Response patterns object
# #' @param newdata New data for prediction
# #' @param k Number of modes (ignored in PCR, kept for compatibility)
# #'
# #' @return Reconstructed field
# #'
# #' @export
# predict_pcr <- function(preds, obs, newdata, k) {
#   .Deprecated("couple_patterns",
#               msg = paste("predict_pcr() is deprecated.",
#                          "PCR-based coupling is no longer supported in the main workflow.",
#                          "Use couple_patterns() with CCA instead, or see legacy_pattern_coupling.qmd"))
#
#   stop("PCR-based prediction is no longer supported in the refactored version. ",
#        "Use couple_patterns() with method='cca' instead, or refer to ",
#        "vignettes/legacy_pattern_coupling.qmd for the original implementation.")
# }
#
# #' Legacy EOT Prediction (Deprecated)
# #'
# #' @param cv Cross-validation object from prep_eot
# #' @param k Number of modes
# #'
# #' @return Reconstructed field
# #'
# #' @export
# predict_eot <- function(cv, k) {
#   .Deprecated("couple_patterns",
#               msg = paste("predict_eot() is deprecated.",
#                          "EOT-based coupling is no longer supported in the main workflow.",
#                          "Use couple_patterns() with CCA instead, or see legacy_pattern_coupling.qmd"))
#
#   stop("EOT-based prediction is no longer supported in the refactored version. ",
#        "Use couple_patterns() with method='cca' instead, or refer to ",
#        "vignettes/legacy_pattern_coupling.qmd for the original implementation.")
# }
#
# #' Legacy Data Preparation for GAM Coupling (Deprecated)
# #'
# #' @param patterns.x Predictor patterns
# #' @param patterns.y Response patterns
# #'
# #' @return Prepared data for GAM fitting
# #'
# #' @export
# prep_data <- function(patterns.x, patterns.y) {
#   .Deprecated("couple_patterns",
#               msg = paste("prep_data() is deprecated.",
#                          "Use couple_patterns() instead which handles data preparation internally."))
#
#   amps.x <- if('patterns' %in% class(patterns.x)) patterns.x$amplitudes else patterns.x
#   amps.y <- if('patterns' %in% class(patterns.y)) patterns.y$amplitudes else patterns.y
#
#   predictors <- amps.x
#   predictands <- tidyr::pivot_longer(amps.y, -time, names_to = 'PC', values_to = 'amplitude')
#
#   dplyr::inner_join(predictands, predictors, by = 'time')
# }
#
# #' Legacy GAM Model Fitting (Deprecated)
# #'
# #' @param data_in Prepared data from prep_data
# #'
# #' @return Fitted GAM models
# #'
# #' @export
# fit_model <- function(data_in) {
#   .Deprecated("couple_patterns",
#               msg = paste("fit_model() is deprecated.",
#                          "GAM-based coupling is no longer supported in the main workflow.",
#                          "Use couple_patterns() instead."))
#
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
#
#   return(model)
# }
#
# #' Legacy PCR Model Fitting (Deprecated)
# #'
# #' @param data_in Prepared data from prep_data
# #'
# #' @return Fitted PCR models
# #'
# #' @export
# fit_pcr <- function(data_in) {
#   .Deprecated("couple_patterns",
#               msg = paste("fit_pcr() is deprecated.",
#                          "PCR-based coupling is no longer supported in the main workflow.",
#                          "Use couple_patterns() instead."))
#
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
#
# #' Convert New Coupled Patterns to Legacy Format
# #'
# #' @param coupled_patterns A coupled_patterns object from the new couple_patterns()
# #'
# #' @return A list mimicking the old coupled_patterns structure
# #'
# #' @export
# as_legacy_coupled_patterns <- function(coupled_patterns) {
#   .Deprecated(msg = "Legacy coupled_patterns format is deprecated. Use the new coupled_patterns object directly.")
#
#   if (!inherits(coupled_patterns, "coupled_patterns")) {
#     stop("Input must be a coupled_patterns object")
#   }
#
#   canonical_corrs <- get_canonical_correlations(coupled_patterns)
#
#   model_data <- data.frame(
#     PC = paste0("PC", 1:coupled_patterns$k),
#     r2 = canonical_corrs$correlation_squared,
#     stringsAsFactors = FALSE
#   )
#
#   legacy_format <- list(
#     model = model_data,
#     data = NULL,
#     method = "cca_converted",
#     note = "Converted from new CCA-based coupled_patterns object"
#   )
#
#   class(legacy_format) <- c("coupled_patterns_legacy", "list")
#
#   return(legacy_format)
# }
#
# #' Migration Helper: Show How to Update Code
# #'
# #' @param legacy_function Name of the legacy function you were using
# #'
# #' @export
# show_migration_guide <- function(legacy_function = c("predict_cca", "predict_gam", "predict_pcr", "predict_eot", "couple_patterns")) {
#
#   legacy_function <- match.arg(legacy_function)
#
#   cat("Migration Guide for", legacy_function, "\n")
#   cat("=" , rep("=", nchar(legacy_function) + 15), "\n", sep = "")
#
#   switch(legacy_function,
#     "predict_cca" = {
#       cat("OLD CODE:\n")
#       cat("result <- predict_cca(preds, obs, newdata, k = 5)\n\n")
#       cat("NEW CODE:\n")
#       cat("coupled <- couple_patterns(preds, obs, k = 5)\n")
#       cat("result <- predict(coupled, newdata)\n\n")
#     },
#     "predict_gam" = {
#       cat("OLD CODE:\n")
#       cat("result <- predict_gam(preds, obs, newdata, k = 5)\n\n")
#       cat("NEW CODE:\n")
#       cat("# GAM coupling is deprecated. Use CCA instead:\n")
#       cat("coupled <- couple_patterns(preds, obs, k = 5, method = 'cca')\n")
#       cat("result <- predict(coupled, newdata)\n\n")
#       cat("# For original GAM implementation, see:\n")
#       cat("# vignettes/legacy_pattern_coupling.qmd\n\n")
#     },
#     "predict_pcr" = {
#       cat("OLD CODE:\n")
#       cat("result <- predict_pcr(preds, obs, newdata, k = 5)\n\n")
#       cat("NEW CODE:\n")
#       cat("# PCR coupling is deprecated. Use CCA instead:\n")
#       cat("coupled <- couple_patterns(preds, obs, k = 5, method = 'cca')\n")
#       cat("result <- predict(coupled, newdata)\n\n")
#       cat("# For original PCR implementation, see:\n")
#       cat("# vignettes/legacy_pattern_coupling.qmd\n\n")
#     },
#     "predict_eot" = {
#       cat("OLD CODE:\n")
#       cat("cv_data <- prep_eot(preds, obs, k_preds, k_obs, k_cca)\n")
#       cat("result <- predict_eot(cv_data, k = 5)\n\n")
#       cat("NEW CODE:\n")
#       cat("# EOT coupling is deprecated. Use CCA instead:\n")
#       cat("coupled <- couple_patterns(preds, obs, k = 5, method = 'cca')\n")
#       cat("result <- predict(coupled, newdata)\n\n")
#       cat("# For original EOT implementation, see:\n")
#       cat("# vignettes/legacy_pattern_coupling.qmd\n\n")
#     },
#     "couple_patterns" = {
#       cat("OLD CODE (GAM-based):\n")
#       cat("coupled <- couple_patterns(patterns.x, patterns.y)\n\n")
#       cat("NEW CODE (CCA-based):\n")
#       cat("coupled <- couple_patterns(patterns.x, patterns.y, k = 5, method = 'cca')\n")
#       cat("# Then use: predict(coupled, newdata)\n\n")
#     }
#   )
#
#   cat("For complete examples and original code, see:\n")
#   cat("vignette('legacy_pattern_coupling', package = 'tidyEOF')\n")
# }
