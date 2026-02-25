# Cross-validate EOF truncation for a single field

Evaluates reconstruction skill for different numbers of EOFs using
k-fold cross-validation. For each fold, EOFs are fit on training data,
test data is projected onto those EOFs, and the reconstruction is
compared to the original test data.

## Usage

``` r
tune_eof(
  data,
  k = 1:10,
  kfolds = 5,
  max_k = max(k),
  metrics = c("rmse", "cor_spatial", "cor_temporal"),
  scale = FALSE,
  monthly = FALSE,
  weight = TRUE
)
```

## Arguments

- data:

  A stars object with spatial-temporal data

- k:

  Vector of EOF counts to evaluate (default 1:10)

- kfolds:

  Number of cross-validation folds (default 5)

- max_k:

  Maximum EOFs to compute per fold (default max(k))

- metrics:

  Character vector of metrics to compute. Options: "rmse",
  "cor_spatial", "cor_temporal" (default: all three)

- scale:

  Logical, whether to scale data before EOF extraction (default FALSE)

- monthly:

  Logical, whether to compute monthly climatology (default FALSE)

- weight:

  Logical, whether to apply area weighting (default TRUE)

## Value

A tibble with columns: k, fold, and one column per metric.

## Examples

``` r
if (FALSE) { # \dontrun{
# Find optimal k for precipitation field
results <- tune_eof(precip_data, k = 1:15, kfolds = 5)
summary <- summarize_eof_cv(results, metric = "rmse")

# Plot reconstruction skill vs k
library(ggplot2)
results %>%
  group_by(k) %>%
  summarize(rmse = mean(rmse)) %>%
  ggplot(aes(k, rmse)) + geom_line() + geom_point()
} # }
```
