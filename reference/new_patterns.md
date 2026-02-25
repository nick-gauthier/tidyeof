# Create a patterns object with validation

Low-level constructor for patterns objects. Validates that all required
components are present and have the correct types.

## Usage

``` r
new_patterns(
  eofs,
  amplitudes,
  eigenvalues,
  k,
  proj_matrix,
  rotation = NULL,
  climatology = NULL,
  units = NULL,
  names = NULL,
  scaled = FALSE,
  monthly = FALSE,
  rotate = FALSE,
  weight = TRUE,
  valid_pixels = NULL
)
```

## Arguments

- eofs:

  A stars object containing EOF spatial patterns

- amplitudes:

  A data.frame/tibble with time column and PC amplitudes

- eigenvalues:

  A data.frame/tibble with eigenvalue statistics

- k:

  Number of components retained

- proj_matrix:

  Projection matrix for new data

- rotation:

  Rotation matrix if varimax was applied (or NULL)

- climatology:

  List with \`mean\` and \`sd\` stars objects (from get_climatology)

- units:

  List of original data units

- names:

  Variable names

- scaled:

  Logical, whether data was scaled

- monthly:

  Logical, whether monthly climatology was used

- rotate:

  Logical, whether rotation was applied

- weight:

  Logical, whether area weighting was applied

- valid_pixels:

  Indices of valid (non-NA) pixels

## Value

A patterns object
