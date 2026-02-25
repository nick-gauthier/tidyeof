# Summarize EOF cross-validation results

Aggregates results across folds and identifies optimal k.

## Usage

``` r
summarize_eof_cv(cv_results, metric = "rmse", minimize = TRUE)
```

## Arguments

- cv_results:

  A tibble from \`tune_eof()\`

- metric:

  Which metric to optimize (default "rmse")

- minimize:

  Logical, whether to minimize (TRUE for RMSE) or maximize (FALSE for
  correlations). Default TRUE.

## Value

A tibble with mean and sd of each metric per k value, sorted by the
target metric. Best k is attached as attribute "best_k".
