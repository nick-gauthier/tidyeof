# Evaluate a single fold with given parameters

Evaluate a single fold with given parameters

## Usage

``` r
evaluate_fold(fold, k_pred, k_resp, k_cca, metrics)
```

## Arguments

- fold:

  A single fold list from cv_folds\$folds

- k_pred:

  Number of predictor EOFs

- k_resp:

  Number of response EOFs

- k_cca:

  Number of CCA modes

- metrics:

  Metrics to compute

## Value

Tibble with fold_id and metric values
