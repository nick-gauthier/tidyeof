# Find Rule N significance cutoff for a patterns object

Tests each eigenvalue in sequence and returns the index of the last
significant mode. This is the maximum k supported by the data according
to the modified Rule N test.

## Usage

``` r
rule_n_cutoff(x, p = 0.05)
```

## Arguments

- x:

  A patterns object

- p:

  Significance level (default 0.05)

## Value

Integer: index of last significant eigenvalue, or 0 if none are
significant
