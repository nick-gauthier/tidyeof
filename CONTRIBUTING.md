# Contributing to tidyeof

Thank you for considering contributing to tidyeof!

## Reporting Issues

If you find a bug or have a feature request, please [open an
issue](https://github.com/nick-gauthier/tidyEOF/issues) on GitHub. When
reporting a bug, please include:

- A minimal reproducible example
- The output of
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)
- What you expected to happen vs what actually happened

## Development Setup

1.  Clone the repository:

    ``` bash
    git clone https://github.com/nick-gauthier/tidyEOF.git
    ```

2.  Install development dependencies:

    ``` r
    install.packages(c("devtools", "testthat", "roxygen2", "patchwork"))
    ```

3.  Load the package for development:

    ``` r
    devtools::load_all()
    ```

4.  Run tests:

    ``` r
    devtools::test()
    ```

## Pull Requests

1.  Fork the repository and create a new branch for your changes
2.  Make your changes and add tests if applicable
3.  Run `devtools::check()` to ensure no new issues are introduced
4.  Submit a pull request with a clear description of the changes

## Code Style

- Follow tidyverse style conventions
- Use `cli` for user-facing messages and errors
- Use [`rlang::abort()`](https://rlang.r-lib.org/reference/abort.html)
  for error conditions
- Document functions with roxygen2
