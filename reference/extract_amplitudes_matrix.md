# Extract amplitude matrix from a patterns object or data frame

Like \[as.matrix.patterns()\] but also works on bare data frames.
Optionally filters to specified times before extracting.

## Usage

``` r
extract_amplitudes_matrix(x, times = NULL)
```

## Arguments

- x:

  A patterns object or a data.frame/tibble with a time column

- times:

  Optional time vector to filter to before extraction

## Value

A numeric matrix
