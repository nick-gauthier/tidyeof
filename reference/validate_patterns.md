# Validate a patterns object for internal consistency

Ensures required components are present and aligned so downstream
functions can rely on a predictable structure. This is chiefly intended
for guarding public entry points like project_patterns() and
reconstruct().

## Usage

``` r
validate_patterns(x, arg = rlang::caller_arg(x), call = rlang::caller_env())
```

## Arguments

- x:

  Object to validate

- arg:

  Argument name for informative error messages

- call:

  Calling environment for error localisation
