# Evaluate EOF reconstruction for a single fold

Evaluate EOF reconstruction for a single fold

## Usage

``` r
evaluate_eof_fold(fold, k, metrics)
```

## Arguments

- fold:

  A fold list containing train_patterns and test_data

- k:

  Number of EOFs to use

- metrics:

  Metrics to compute

## Value

Tibble with fold_id and metric values
