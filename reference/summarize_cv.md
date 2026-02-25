# Summarize cross-validation results

Aggregates results across folds and identifies best hyperparameters.

## Usage

``` r
summarize_cv(cv_results, metric = "rmse", minimize = TRUE)
```

## Arguments

- cv_results:

  A tibble from \`tune_cca()\`

- metric:

  Which metric to optimize (default "rmse")

- minimize:

  Logical, whether to minimize (TRUE for RMSE) or maximize (FALSE for
  correlations). Default TRUE.

## Value

A tibble with mean and sd of each metric per parameter combination,
sorted by the target metric. Best parameters are attached as attribute
"best_params".

## Examples

``` r
if (FALSE) { # \dontrun{
results <- tune_cca(cv, k_pred = 1:5, k_resp = 1:5)
summary <- summarize_cv(results, metric = "rmse", minimize = TRUE)

# Get best parameters (k_pred, k_resp, k_cca)
attr(summary, "best_params")
} # }
```
