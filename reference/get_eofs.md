# Internal function to calculate EOFs and related components

Internal function to calculate EOFs and related components

## Usage

``` r
get_eofs(dat, k, rotate = FALSE, irlba_threshold, weights = NULL)
```

## Arguments

- weights:

  Optional numeric vector of spatial weights (one per spatial location)
  to be applied to the anomalies before PCA
