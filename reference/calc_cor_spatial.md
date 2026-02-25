# Calculate spatial correlation (averaged over time)

For each time step, compute correlation across spatial locations, then
average across all time steps.

## Usage

``` r
calc_cor_spatial(pred_matrix, obs_matrix)
```

## Arguments

- pred_matrix:

  Predicted values matrix (time x space)

- obs_matrix:

  Observed values matrix (time x space)

## Value

Mean spatial correlation across time steps
