# Get Canonical Variables from Coupled Patterns

Extract canonical variables from either predictor or response patterns

## Usage

``` r
get_canonical_variables(
  object,
  data,
  type = c("predictor", "response"),
  k = NULL
)
```

## Arguments

- object:

  A coupled_patterns object

- data:

  Original data (patterns object or amplitudes tibble)

- type:

  Either "predictor" or "response"

- k:

  Number of canonical modes to extract

## Value

Tibble with canonical variables
