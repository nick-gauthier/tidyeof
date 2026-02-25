# Calculate temporal correlation (averaged over space)

For each spatial location, compute correlation across time, then average
across all locations.

## Usage

``` r
calc_cor_temporal(pred_matrix, obs_matrix)
```

## Arguments

- pred_matrix:

  Predicted values matrix (time x space)

- obs_matrix:

  Observed values matrix (time x space)

## Value

Mean temporal correlation across spatial locations
