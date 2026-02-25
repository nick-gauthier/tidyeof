# Convert a time-by-space matrix back to a stars object

Convert a time-by-space matrix back to a stars object

## Usage

``` r
matrix_to_spacetime(
  mat,
  template_eofs,
  spatial_template,
  valid_pixels,
  times,
  var_names
)
```

## Arguments

- mat:

  Matrix with rows = time, columns = flattened spatial cells

- template_eofs:

  EOF stars object providing spatial metadata

- spatial_template:

  Stars object supplying spatial dimension metadata (e.g., the
  climatology)

- valid_pixels:

  Integer indices of spatial cells with valid data

- times:

  Vector of time values

- var_names:

  Character vector of attribute names for the result
