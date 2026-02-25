#' Validate a patterns object for internal consistency
#'
#' Ensures required components are present and aligned so downstream functions
#' can rely on a predictable structure. This is chiefly intended for
#' guarding public entry points like project_patterns() and reconstruct().
#'
#' @param x Object to validate
#' @param arg Argument name for informative error messages
#' @param call Calling environment for error localisation
#' @keywords internal
validate_patterns <- function(x,
                              arg = rlang::caller_arg(x),
                              call = rlang::caller_env()) {
  if (!inherits(x, "patterns")) {
    cli::cli_abort(
      "Argument {.arg {arg}} must be a {.cls patterns} object, not {.cls {class(x)}}.",
      call = call
    )
  }

  required <- c("eofs", "amplitudes", "climatology", "k")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    cli::cli_abort(
      "Patterns object is missing required component{?s}: {missing*}.",
      call = call
    )
  }

  if (!inherits(x$eofs, "stars")) {
    cli::cli_abort(
      "Patterns object must store EOFs as a {.cls stars} object.",
      call = call
    )
  }

  eof_dims <- stars::st_dimensions(x$eofs)
  if (!"PC" %in% names(eof_dims)) {
    cli::cli_abort(
      "EOF object must contain a {.field PC} dimension.",
      call = call
    )
  }

  pc_values <- stars::st_get_dimension_values(x$eofs, "PC")
  n_pc <- length(pc_values)
  if (!isTRUE(all.equal(n_pc, x$k))) {
    cli::cli_abort(
      "Number of EOF modes ({n_pc}) does not match {.field k} ({x$k}).",
      call = call
    )
  }

  if (!"time" %in% names(x$amplitudes)) {
    cli::cli_abort(
      "Amplitude tibble must include a {.field time} column.",
      call = call
    )
  }

  amplitude_names <- setdiff(names(x$amplitudes), "time")
  if (!identical(amplitude_names, pc_values)) {
    cli::cli_abort(
      "Amplitude columns {amplitude_names*} do not match expected {pc_values*}.",
      call = call
    )
  }

  if (is.matrix(x$rotation) && ncol(x$rotation) != x$k) {
    cli::cli_abort(
      "Rotation matrix must have {.field k} columns (found {ncol(x$rotation)}).",
      call = call
    )
  }

  if (is.null(x$valid_pixels)) {
    cli::cli_abort(
      "Patterns object must carry {.field valid_pixels} indices for spatial alignment.",
      call = call
    )
  }

  invisible(x)
}
