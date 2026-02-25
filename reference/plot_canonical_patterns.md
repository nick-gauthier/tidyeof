# Plot canonical spatial patterns

Internal helper to plot canonical patterns computed from CCA.

## Usage

``` r
plot_canonical_patterns(x, side = "both", k = NULL, scaled = FALSE)
```

## Arguments

- x:

  A coupled_patterns object

- side:

  "predictor", "response", or "both"

- k:

  Number of modes

- scaled:

  Whether to scale the patterns

## Value

A ggplot or patchwork object
