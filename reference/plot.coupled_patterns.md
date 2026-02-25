# Plot coupled patterns diagnostics

Provides visualization helpers for \`coupled_patterns\` objects. The
default \`type = "combined"\` mirrors the patterns plotting workflow by
displaying the predictor and response spatial patterns alongside their
canonical variate time series. Additional types include canonical
correlation bars, standalone canonical variate panels, canonical spatial
patterns, or direct access to the underlying predictor / response
pattern plots.

## Usage

``` r
# S3 method for class 'coupled_patterns'
plot(
  x,
  type = c("combined", "correlations", "canonical", "canonical_patterns", "predictor",
    "response"),
  side = c("predictor", "response", "both"),
  data = NULL,
  k = NULL,
  scaled = FALSE,
  ...
)
```

## Arguments

- x:

  A \`coupled_patterns\` object.

- type:

  Plot type: one of \`"combined"\`, \`"correlations"\`, \`"canonical"\`,
  \`"canonical_patterns"\`, \`"predictor"\`, or \`"response"\`.

- side:

  When \`type = "canonical"\` or \`type = "canonical_patterns"\`, choose
  from the predictor side, response side, or both (values:
  \`"predictor"\`, \`"response"\`, \`"both"\`). Ignored for other plot
  types.

- data:

  Optional amplitudes or patterns for the canonical variate plots. For
  \`side = "both"\`, a list with elements \`predictor\` and \`response\`
  may be supplied. Defaults to the training patterns stored in \`x\`
  when omitted.

- k:

  Number of canonical modes to display (defaults to all available).

- scaled:

  Logical, passed to the underlying pattern plots when relevant
  (defaults to \`FALSE\`).

- ...:

  Additional arguments forwarded to \`plot.patterns()\` for \`type =
  "predictor"\` or \`type = "response"\` calls.

## Value

A ggplot object (or a patchwork object when \`type = "combined"\` or
\`type = "canonical_patterns"\` with \`side = "both"\`).
