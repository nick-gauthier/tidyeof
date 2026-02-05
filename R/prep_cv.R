#' Prepare contiguous cross-validation folds
#'
#' Divides time steps into k contiguous folds for temporal cross-validation.
#' Folds are kept contiguous to preserve temporal structure.
#'
#' @param times Vector of time values to split
#' @param kfolds Number of folds (default 5)
#'
#' @return List of k vectors, each containing the time values for that fold
#'
#' @seealso [prep_cv_folds()] for preparing complete CV folds with pre-computed patterns
#' @seealso [tune_cca()] for hyperparameter grid search
#' @export
prep_folds <- function(times, kfolds = 5){
  # divide years into kfolds contiguous folds
  n <- length(times)
  r <- n %% kfolds
  fold_times <- rep(n %/% kfolds, kfolds)
  if(r > 0) fold_times[1:r] <- fold_times[1:r] + 1
  tibble(time = times, fold = rep(1:kfolds, times = fold_times)) %>%
    group_nest(fold) %>%
    pull(data) %>%
    purrr::map(pull)
}

#' Prepare CCA cross-validation data (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Use [prep_cv_folds()] instead for modern cross-validation workflows.
#'
#' Sets up k-fold cross-validation for CCA-based pattern coupling. Extracts
#' patterns from training folds and prepares test data for each fold.
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
#' @export
prep_cca <- function(preds, obs, k_preds, k_obs, kfolds = 5, scale = FALSE, rotate = FALSE, monthly = FALSE, weight = TRUE) {
  .Deprecated("prep_cv_folds")
  # find the years of overlap between the response and predictor fields
  pred_times <- st_get_dimension_values(preds, 'time')
  obs_times <- st_get_dimension_values(obs, 'time')
  # Use subsetting to preserve Date/POSIXct class (intersect can strip it)
  time_steps <- pred_times[pred_times %in% obs_times]
  preds <- filter(preds, time %in% time_steps)
  obs <- filter(obs, time %in% time_steps)
  folds <- prep_folds(time_steps, kfolds = kfolds)

  # preprocess the training data for each fold
  train_obs <- purrr::map(folds, ~ filter(obs, !(time %in% .))  %>%
                            patterns(k = k_obs, scale = scale, rotate = rotate, monthly = monthly, weight = weight))

  train_preds <- purrr::map(folds, ~ filter(preds, !(time %in% .)) %>%
                              patterns(k = k_preds, scale = scale, rotate = rotate, monthly = monthly, weight = weight))

  # preprocess test data for each fold
  test <- purrr::map(folds, ~ filter(preds, time %in% .))

  tibble(train_obs, train_preds, test)
}

#' Prepare delta method cross-validation data
#'
#' Sets up k-fold cross-validation for delta-based downscaling.
#'
#' @param preds Predictor stars object
#' @param obs Observed/response stars object
#'
#' @return Tibble with train_preds, train_obs, and test columns for each fold
#' @export
prep_delta <- function(preds, obs) {
  # find the years of overlap between the response and predictor fields
  pred_times <- st_get_dimension_values(preds, 'time')
  obs_times <- st_get_dimension_values(obs, 'time')
  # Use subsetting to preserve Date/POSIXct class (intersect can strip it)
  time_steps <- pred_times[pred_times %in% obs_times]
  preds <- filter(preds, time %in% time_steps)
  obs <- filter(obs, time %in% time_steps)
  folds <- prep_folds(time_steps) # this does 5 fold by default, but could change

  # preprocess the training data for each fold
  train_obs <- purrr::map(folds, ~ filter(obs, !(time %in% .)))

  train_preds <- purrr::map(folds, ~ filter(preds, !(time %in% .)))

  # preprocess test data for each fold
  test <- purrr::map(folds, ~ filter(preds, time %in% .))

  tibble(train_preds, train_obs, test)
}