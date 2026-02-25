# Calculate anomalies from a climatological mean

Calculate anomalies from a climatological mean

## Usage

``` r
get_anomalies(dat, clim = NULL, scale = FALSE, monthly = FALSE)
```

## Arguments

- dat:

  A stars object with dimensions (x, y, time)

- clim:

  Optional climatology from get_climatology(). If NULL, computed
  internally

- scale:

  Logical. If TRUE, divide by standard deviation

- monthly:

  Logical. If TRUE, compute monthly anomalies

## Value

A stars object with anomalies
