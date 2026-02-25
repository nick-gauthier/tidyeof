# Plot method for common_patterns objects

Shows shared EOFs on top and overlaid amplitude time series (colored by
source) on the bottom.

## Usage

``` r
# S3 method for class 'common_patterns'
plot(
  x,
  scale = c("standardized", "variance", "raw"),
  scale_y = c("fixed", "free"),
  overlay = NULL,
  overlay_color = "grey30",
  overlay_fill = NA,
  ...
)
```

## Arguments

- x:

  A common_patterns object

- scale:

  Amplitude scaling: "standardized" (default), "variance", or "raw"

- scale_y:

  Y-axis scaling: "fixed" (default) or "free"

- overlay:

  Optional sf object to overlay on EOF maps

- overlay_color:

  Color for overlay geometry (default "grey30")

- overlay_fill:

  Fill for overlay geometry (default NA)

- ...:

  Additional arguments (currently unused)

## Value

A patchwork object (EOFs + amplitudes)
