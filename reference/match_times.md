# Match time steps between amplitudes and a stars object

Match time steps between amplitudes and a stars object

## Usage

``` r
match_times(amplitudes, dat)
```

## Arguments

- amplitudes:

  A tibble with a time column

- dat:

  A stars object with a time dimension

## Value

A list with \`times\` (common time steps) and \`amps\` (filtered
amplitude matrix)
