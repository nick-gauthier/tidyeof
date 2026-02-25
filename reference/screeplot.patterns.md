# Scree plot for EOF patterns

Scree plot for EOF patterns

## Usage

``` r
# S3 method for class 'patterns'
screeplot(x, k = NULL, kmax = 10, rule_n = FALSE, ...)
```

## Arguments

- x:

  A patterns object from patterns()

- k:

  Optional number of components to highlight with vertical line

- kmax:

  Maximum number of components to show (default 10)

- rule_n:

  Logical, whether to show the modified Rule N significance cutoff as a
  dashed blue line (default FALSE)

- ...:

  Additional arguments (currently unused)

## Value

A ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
pat <- patterns(data, k = 5)
screeplot(pat)
screeplot(pat, k = 3, kmax = 8)
screeplot(pat, rule_n = TRUE)  # show significance cutoff
} # }
```
