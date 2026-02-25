# Compute area-based weights for spatial data

Calculates area weights for spatial data using st_area(). Works
uniformly for raster grids, irregular geometries, and different
coordinate systems.

## Usage

``` r
area_weights(dat)
```

## Arguments

- dat:

  A stars object with spatial dimensions

## Value

Numeric weights (sqrt of normalized areas)
