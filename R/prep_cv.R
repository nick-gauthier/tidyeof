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

#' Prepare CCA cross-validation data
#'
#' Sets up k-fold cross-validation for CCA-based pattern coupling. Extracts
#' patterns from training folds and prepares test data for each fold.
#'
#' @param preds Predictor stars object
#' @param obs Observed/response stars object
#' @param k_preds Number of predictor PCs to retain
#' @param k_obs Number of response PCs to retain
#' @param kfolds Number of cross-validation folds (default 5)
#' @param scale Passed to get_patterns()
#' @param rotate Passed to get_patterns()
#' @param monthly Passed to get_patterns()
#' @param weight Passed to get_patterns()
#'
#' @return Tibble with train_obs, train_preds, and test columns for each fold
#' @export
prep_cca <- function(preds, obs, k_preds, k_obs, kfolds = 5, scale = FALSE, rotate = FALSE, monthly = FALSE, weight = TRUE) {
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
                            get_patterns(k = k_obs, scale = scale, rotate = rotate, monthly = monthly, weight = weight))

  train_preds <- purrr::map(folds, ~ filter(preds, !(time %in% .)) %>%
                              get_patterns(k = k_preds, scale = scale, rotate = rotate, monthly = monthly, weight = weight))

  # preprocess test data for each fold
  test <- purrr::map(folds, ~ filter(preds, time %in% .))

  tibble(train_obs, train_preds, test)
}


# =============================================================================
# LEGACY: prep_eot - EOT-based cross-validation using remote package
# Moved to R_archive/ for reference. Uses anomalize(), denoise(), eot() from
# remote package and as('Raster') conversion. Not compatible with current
# stars-based architecture.
# =============================================================================
# #' @export
# prep_eot <- function(preds, obs, k_preds, k_obs, k_cca){
#   # find the years of overlap between the response and predictor fields
#   time_steps <- intersect(st_get_dimension_values(preds, 'time'),
#                           st_get_dimension_values(obs, 'time'))
#   preds <- filter(preds, time %in% time_steps)
#   obs <- filter(obs, time %in% time_steps)
#   folds <- prep_folds(time_steps) # this does 5 fold by default, but could change
#
#   # preprocess the training data for each fold
#   obs_train_rast <- purrr::map(folds, ~ filter(obs, !(time %in% .)) %>%
#                                  as('Raster'))
#
#   pred_train_rast <-  purrr::map(folds, ~ filter(preds, !(time %in% .)) %>%
#                                    as('Raster'))
#
#   obs_train <- purrr::map(obs_train_rast, ~ anomalize(.) %>%
#                             denoise(k = k_obs, weighted = TRUE, verbose = FALSE))
#
#   pred_train <-  purrr::map(pred_train_rast, ~ anomalize(.) %>%
#                               denoise(k = k_preds, weighted = TRUE, verbose = FALSE))
#
#   # preprocess test data for each fold
#   test <- purrr::map(folds, ~ filter(preds, time %in% .) %>%
#                        as('Raster')) %>% # test years
#     map2(pred_train_rast, ~.x - mean(.y)) # subtract the training predictor mean from the test predictors
#
#   means <- purrr::map(obs_train_rast, raster::mean)
#
#   # fit model to training data
#   eots <- map2(pred_train, obs_train, ~eot(.x, .y, n = k_cca, standardised = FALSE, verbose = FALSE))
#
#   list(eot = eots, # this hard codes 10 patterns, but could be changed to min(k_preds, k_obs)
#        test = test,
#        mean = means)
# }

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