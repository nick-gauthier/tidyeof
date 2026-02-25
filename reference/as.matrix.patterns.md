# Extract amplitude matrix from a patterns object or data frame

Provides a unified way to get the numeric amplitude matrix (without the
time column) from either a patterns object or a bare amplitudes data
frame. This is the single dispatch point used throughout the package.

## Usage

``` r
# S3 method for class 'patterns'
as.matrix(x, ...)
```

## Arguments

- x:

  A patterns object or a data.frame/tibble with a time column

- ...:

  Additional arguments (currently unused)

## Value

A numeric matrix with rows = time steps, columns = PCs
