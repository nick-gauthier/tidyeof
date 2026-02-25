# Plot method for patterns objects

Plot method for patterns objects

## Usage

``` r
# S3 method for class 'patterns'
plot(
  x,
  type = "combined",
  scaled = FALSE,
  rawdata = NULL,
  scale = c("standardized", "variance", "raw"),
  scale_y = c("fixed", "free"),
  events = NULL,
  overlay = NULL,
  overlay_color = "grey30",
  overlay_fill = NA,
  ...
)
```

## Arguments

- x:

  A patterns object

- type:

  Type of plot: "combined" (default), "eofs", or "amplitudes"

- scaled:

  For EOFs: show correlations (TRUE) or raw loadings (FALSE)

- rawdata:

  Optional raw data for correlation calculation when scaled = TRUE

- scale:

  For amplitudes: scaling method ("standardized", "variance", "raw")

- scale_y:

  For amplitudes: y-axis scaling ("fixed" or "free")

- events:

  For amplitudes: optional dates to mark with vertical lines

- overlay:

  Optional sf object to overlay on EOF maps (e.g., coastlines,
  boundaries)

- overlay_color:

  Color for overlay geometry (default "grey30")

- overlay_fill:

  Fill for overlay geometry (default NA for no fill)

- ...:

  Additional arguments (currently unused)

## Value

A ggplot2 object or patchwork object for combined plots
