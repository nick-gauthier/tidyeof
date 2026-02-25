# Flatten a time-indexed stars object into a matrix

Internal helper to turn a single-attribute stars object with a \`time\`
dimension into a matrix of shape time x space along with the metadata
needed to reconstruct the original spatial layout.

## Usage

``` r
flatten_time_space(dat)
```

## Arguments

- dat:

  A stars object with a \`time\` dimension

## Value

A list containing the flattened matrix, spatial dimension names, their
lengths, and coordinate values
