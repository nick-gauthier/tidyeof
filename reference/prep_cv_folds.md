# Prepare cross-validation folds for CCA tuning

Computes EOF patterns at maximum k once per fold, enabling cheap
truncation during grid search. This is the expensive step - call it
once, then use \`tune_cca()\` for fast hyperparameter exploration.

## Usage

``` r
prep_cv_folds(
  predictor,
  response,
  kfolds = 5,
  max_k_pred = 10,
  max_k_resp = 10,
  scale = FALSE,
  scale_pred = NULL,
  scale_resp = NULL,
  rotate = FALSE,
  monthly = FALSE,
  weight = TRUE,
  common_with = NULL
)
```

## Arguments

- predictor:

  A stars object with predictor data

- response:

  A stars object with response data

- kfolds:

  Number of cross-validation folds (default 5)

- max_k_pred:

  Maximum number of predictor EOFs to compute (default 10)

- max_k_resp:

  Maximum number of response EOFs to compute (default 10)

- scale:

  Logical, whether to scale data before EOF extraction (default FALSE).
  Sets the default for both predictor and response; use \`scale_pred\`
  and/or \`scale_resp\` to override individually.

- scale_pred:

  Logical, whether to scale predictor data before EOF extraction.
  Overrides \`scale\` for the predictor when not NULL (default NULL).

- scale_resp:

  Logical, whether to scale response data before EOF extraction.
  Overrides \`scale\` for the response when not NULL (default NULL).

- rotate:

  Logical, whether to apply varimax rotation (default FALSE)

- monthly:

  Logical, whether to compute monthly climatology (default FALSE)

- weight:

  Logical, whether to apply area weighting (default TRUE)

- common_with:

  Optional named list of additional stars objects to include in common
  EOF computation via \[common_patterns()\]. When provided, predictor
  patterns are computed jointly with these datasets. The primary
  predictor is included under the name \`.primary\`. Test times are
  excluded from all datasets to prevent leakage.

## Value

A cv_folds S3 object containing:

- folds:

  List of fold data, each containing train patterns and test data

- max_k_pred:

  Maximum predictor k used

- max_k_resp:

  Maximum response k used

- kfolds:

  Number of folds

- common_times:

  Vector of overlapping time steps

- pattern_opts:

  List of pattern extraction options

- common_with_sources:

  Names of common EOF sources (if used)

## Examples

``` r
if (FALSE) { # \dontrun{
# Standard folds
cv <- prep_cv_folds(coarse_data, fine_data,
                    kfolds = 5, max_k_pred = 10, max_k_resp = 10)

# With common EOFs from additional sources
cv <- prep_cv_folds(coarse_data, fine_data,
                    common_with = list(phyda = phyda_coarse),
                    kfolds = 5, max_k_pred = 10, max_k_resp = 10)

# Then tune over k_pred and k_resp (unchanged)
results <- tune_cca(cv, k_pred = 1:10, k_resp = 1:10)
} # }
```
