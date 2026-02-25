# Tidy PCA standard deviations into a tibble

Replaces \`broom::tidy(pca, matrix = "pcs")\` with a dependency-free
version. Returns a tibble with columns PC, std.dev, percent, cumulative
— identical to the broom output.

## Usage

``` r
tidy_pca_sdev(pca_result)
```

## Arguments

- pca_result:

  A prcomp (or prcomp_irlba) result

## Value

Tibble with PC, std.dev, percent, cumulative
