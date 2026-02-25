# Compute spatial error metrics between predicted and observed fields

Compute spatial error metrics between predicted and observed fields

## Usage

``` r
compute_spatial_metrics(
  predicted,
  observed,
  metrics = c("rmse", "cor_spatial", "cor_temporal")
)
```

## Arguments

- predicted:

  A stars object with predicted values

- observed:

  A stars object with observed values

- metrics:

  Character vector of metrics to compute. Options: "rmse",
  "cor_spatial", "cor_temporal"

## Value

Named list of computed metrics
