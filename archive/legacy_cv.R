# Legacy CV functions — moved from R/prep_cv.R
# Superseded by prep_cv_folds() + tune_cca()

#' Prepare CCA cross-validation data (deprecated)
#'
#' Use [prep_cv_folds()] instead for modern cross-validation workflows.
#'
#' @param preds Predictor stars object
#' @param obs Observed/response stars object
#' @param k_preds Number of predictor PCs to retain
#' @param k_obs Number of response PCs to retain
#' @param kfolds Number of cross-validation folds (default 5)
#' @param scale Passed to patterns()
#' @param rotate Passed to patterns()
#' @param monthly Passed to patterns()
#' @param weight Passed to patterns()
#'
#' @return Tibble with train_obs, train_preds, and test columns for each fold
prep_cca <- function(preds, obs, k_preds, k_obs, kfolds = 5, scale = FALSE, rotate = FALSE, monthly = FALSE, weight = TRUE) {
  .Deprecated("prep_cv_folds")
  pred_times <- st_get_dimension_values(preds, 'time')
  obs_times <- st_get_dimension_values(obs, 'time')
  time_steps <- pred_times[pred_times %in% obs_times]
  preds <- filter(preds, time %in% time_steps)
  obs <- filter(obs, time %in% time_steps)
  folds <- prep_folds(time_steps, kfolds = kfolds)

  train_obs <- purrr::map(folds, ~ filter(obs, !(time %in% .))  %>%
                            patterns(k = k_obs, scale = scale, rotate = rotate, monthly = monthly, weight = weight))

  train_preds <- purrr::map(folds, ~ filter(preds, !(time %in% .)) %>%
                              patterns(k = k_preds, scale = scale, rotate = rotate, monthly = monthly, weight = weight))

  test <- purrr::map(folds, ~ filter(preds, time %in% .))

  tibble(train_obs, train_preds, test)
}

#' Prepare delta method cross-validation data
#'
#' @param preds Predictor stars object
#' @param obs Observed/response stars object
#'
#' @return Tibble with train_preds, train_obs, and test columns for each fold
prep_delta <- function(preds, obs) {
  pred_times <- st_get_dimension_values(preds, 'time')
  obs_times <- st_get_dimension_values(obs, 'time')
  time_steps <- pred_times[pred_times %in% obs_times]
  preds <- filter(preds, time %in% time_steps)
  obs <- filter(obs, time %in% time_steps)
  folds <- prep_folds(time_steps)

  train_obs <- purrr::map(folds, ~ filter(obs, !(time %in% .)))
  train_preds <- purrr::map(folds, ~ filter(preds, !(time %in% .)))
  test <- purrr::map(folds, ~ filter(preds, time %in% .))

  tibble(train_preds, train_obs, test)
}
