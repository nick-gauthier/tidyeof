# Compute sign vector for consistent EOF orientation

Determines the sign (+1 or -1) needed for each EOF component so that the
dominant loading is positive. Works for both raster and sf geometry
stars objects.

## Usage

``` r
compute_eof_signs(eofs)
```

## Arguments

- eofs:

  A stars object with EOF spatial patterns (must have a PC dimension)

## Value

Named numeric vector of +1/-1 values, one per PC
