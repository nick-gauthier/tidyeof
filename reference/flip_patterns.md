# Flip EOF patterns to have consistent sign

Ensures all EOF patterns have positive-dominant loadings by flipping the
sign of both the spatial pattern and corresponding amplitude time series
when needed. This makes plotting and interpretation more consistent
across analyses.

## Usage

``` r
flip_patterns(patterns)
```

## Arguments

- patterns:

  A patterns object from patterns()

## Value

The patterns object with signs adjusted for consistency
