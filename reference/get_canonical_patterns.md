# Get Canonical Spatial Patterns from Coupled Patterns

Computes the spatial patterns corresponding to each canonical mode by
taking linear combinations of the original EOFs weighted by CCA
coefficients. These are the spatial patterns that, when projected onto
the data, yield the canonical variates.

## Usage

``` r
get_canonical_patterns(object, type = c("predictor", "response"), k = NULL)
```

## Arguments

- object:

  A coupled_patterns object

- type:

  Either "predictor" or "response"

- k:

  Number of canonical modes to extract (default: all available)

## Value

A stars object with canonical spatial patterns (dimension "CV" instead
of "PC")

## Examples

``` r
if (FALSE) { # \dontrun{
coupled <- couple(pred_patterns, resp_patterns, k = 3)

# Get canonical patterns for response side
resp_canonical <- get_canonical_patterns(coupled, type = "response")
plot(resp_canonical)

# Compare to original EOFs
plot(coupled$response_patterns$eofs)
} # }
```
