#' Plot EOF spatial patterns (deprecated)
#'
#' @param patterns A patterns object
#' @param scaled Logical, whether to show correlations (TRUE) or raw loadings (FALSE)
#' @param rawdata Optional raw data for correlation calculation when scaled = TRUE
#' @return A ggplot object
#' @export
#' @seealso \code{\link{plot.patterns}} for the recommended plotting interface
#'
#' @details
#' \strong{This function is deprecated.} Use \code{plot(patterns)} or
#' \code{plot(patterns, type = "eofs")} instead.
plot_eofs <- function(patterns, scaled = FALSE, rawdata = NULL){
  .Deprecated("plot.patterns", package = "tidyEOF",
              msg = "plot_eofs() is deprecated. Use plot(patterns) or plot(patterns, type = 'eofs') instead.")

  # Call the new method
  plot.patterns(patterns, type = "eofs", scaled = scaled, rawdata = rawdata)
}

#' Plot PC amplitude time series (deprecated)
#'
#' @param patterns A patterns object
#' @param scale How to scale amplitudes: "standardized" (default), "variance", or "raw"
#' @param scale_y Y-axis scaling: "fixed" (default) or "free"
#' @param events Optional vector of dates to mark with vertical lines
#' @return A ggplot object
#' @export
#' @seealso \code{\link{plot.patterns}} for the recommended plotting interface
#'
#' @details
#' \strong{This function is deprecated.} Use \code{plot(patterns, type = "amplitudes")}
#' with additional parameters instead.
plot_amps <- function(patterns,
                     scale = c("standardized", "variance", "raw"),
                     scale_y = c("fixed", "free"),
                     events = NULL) {

  .Deprecated("plot.patterns", package = "tidyEOF",
              msg = "plot_amps() is deprecated. Use plot(patterns, type = 'amplitudes', scale = '...', scale_y = '...') instead.")

  scale <- match.arg(scale)
  scale_y <- match.arg(scale_y)

  # Call the new method
  plot.patterns(patterns, type = "amplitudes", scale = scale, scale_y = scale_y, events = events)
}

#' Scree plot for EOF patterns
#'
#' @param x A patterns object from patterns()
#' @param k Optional number of components to highlight with vertical line
#' @param kmax Maximum number of components to show (default 10)
#' @param ... Additional arguments (currently unused)
#' @return A ggplot object
#' @export
#' @examples
#' \dontrun{
#' pat <- patterns(data, k = 5)
#' screeplot(pat)
#' screeplot(pat, k = 3, kmax = 8)
#' }
screeplot.patterns <- function(x, k = NULL, kmax = 10, ...) {
  x$eigenvalues %>%
    dplyr::mutate(separated = if_else(is.na(lag(low)), TRUE, hi < lag(low)),
           multiplet = as.factor(cumsum(separated))) %>%
    filter(PC <= kmax) %>%
    ggplot2::ggplot(aes(x = PC, y = percent)) +
    ggplot2::geom_linerange(aes(x = PC, ymin = low, ymax = hi)) +
    ggplot2::geom_point(size = 2, aes(color = multiplet)) +
    ggplot2::geom_text(aes(x = PC, y = cumvar_line, label = glue::glue("{round(cumulative, 0)}%")), size = 2.5, vjust = 0) +
    ggplot2::labs(x = "Principal Component", y = "Normalized Eigenvalue") +
    ggplot2::geom_vline(xintercept = k + .5, linetype = 2, color = 'red', alpha = .7) +
    ggplot2::theme_bw() +
    ggplot2::guides(color = 'none') +
    ggplot2::scale_x_continuous(breaks = seq(0, 25, 5)) +
    ggplot2::scale_color_brewer(palette = 'Spectral')
}

#' Scree plot (deprecated)
#'
#' @param dat Data object
#' @param k Number of components
#' @param kmax Maximum components to show
#' @param scale Scaling option
#' @param monthly Monthly option
#' @param weight Weighting option
#' @export
#' @details
#' \strong{This function is deprecated.} Use \code{screeplot(pat)} instead
#' where pat is the result of \code{patterns()}.
plot_scree <- function(dat, k = NULL, kmax = 10, scale = FALSE, monthly = FALSE, weight = TRUE) {
  .Deprecated("screeplot.patterns", package = "tidyEOF",
              msg = "plot_scree() is deprecated. Use screeplot(patterns(dat, ...)) instead.")

  # For backward compatibility, try to create patterns and call new method
  tryCatch({
    pat <- patterns(dat, k = if(is.null(k)) 5 else k, scale = scale, monthly = monthly, weight = weight)
    screeplot.patterns(pat, k = k, kmax = kmax)
  }, error = function(e) {
    stop("plot_scree() is deprecated and could not convert to new interface. Use screeplot(patterns(dat, ...)) instead.")
  })
}

#' Calculate significance of EOF patterns using North's Rule of Thumb
#'
#' @param patterns A patterns object
#' @param confidence Confidence level (default 0.95)
#' @return A stars object with significance regions
#' @keywords internal
calculate_significance <- function(patterns, confidence = 0.95) {
  # Get key dimensions
  n_times <- nrow(patterns$amplitudes)
  eigenvals <- patterns$eigenvalues$eigenvalues

  # Calculate sampling error for eigenvalues using North's rule
  # Delta lambda = lambda * sqrt(2/n)
  sampling_errors <- eigenvals * sqrt(2/n_times)

  # Test if eigenvalues are distinct
  # If error bars overlap, patterns may be mixed
  is_distinct <- logical(length(eigenvals))
  for(i in seq_along(eigenvals)) {
    if(i < length(eigenvals)) {
      # Check if error bars overlap with next eigenvalue
      delta_i <- sampling_errors[i]
      delta_j <- sampling_errors[i + 1]
      separation <- eigenvals[i] - eigenvals[i + 1]
      is_distinct[i] <- separation > (delta_i + delta_j)
    } else {
      # Last eigenvalue is compared with zero
      is_distinct[i] <- eigenvals[i] > sampling_errors[i]
    }
  }

  # Create significance field for each EOF
  # Use same dimensions as original EOFs
  signif_array <- patterns$eofs[[1]] * 0  # Initialize with zeros

  for(i in seq_len(patterns$k)) {
    if(is_distinct[i]) {
      # If eigenvalue is distinct, mark areas where EOF loading
      # is above threshold as significant
      pattern <- abs(patterns$eofs[[1]][,,i])
      # Use standard error to set threshold
      threshold <- sqrt(1/n_times)  # Simple approximation for significance
      signif_array[,,i] <- pattern > threshold
    }
  }

  # Convert to stars object
  result <- patterns$eofs
  result[[1]] <- signif_array
  names(result) <- "significance"

  # Store significance info in attributes
  attr(result, "is_distinct") <- is_distinct
  attr(result, "sampling_errors") <- sampling_errors
  attr(result, "n_effective") <- n_times

  result
}

#' Print significance test results
#'
#' @param signif_obj Significance object from calculate_significance
#' @export
print_significance <- function(signif_obj) {
  is_distinct <- attr(signif_obj, "is_distinct")
  sampling_errors <- attr(signif_obj, "sampling_errors")
  n_eff <- attr(signif_obj, "n_effective")

  cat("\nEOF Significance Testing (North's Rule of Thumb)\n")
  cat("Number of effective samples:", n_eff, "\n\n")

  for(i in seq_along(is_distinct)) {
    cat(glue::glue("EOF {i}: {ifelse(is_distinct[i], 'Distinct', 'May be degenerate')} (±{round(sampling_errors[i], 3)})\n"))
  }
}
#' Plot method for patterns objects
#'
#' @param x A patterns object
#' @param type Type of plot: "combined" (default), "eofs", or "amplitudes"
#' @param scaled For EOFs: show correlations (TRUE) or raw loadings (FALSE)
#' @param rawdata Optional raw data for correlation calculation when scaled = TRUE
#' @param scale For amplitudes: scaling method ("standardized", "variance", "raw")
#' @param scale_y For amplitudes: y-axis scaling ("fixed" or "free")
#' @param events For amplitudes: optional dates to mark with vertical lines
#' @param overlay Optional sf object to overlay on EOF maps (e.g., coastlines, boundaries)
#' @param overlay_color Color for overlay geometry (default "grey30")
#' @param overlay_fill Fill for overlay geometry (default NA for no fill)
#' @param ... Additional arguments (currently unused)
#' @return A ggplot2 object or patchwork object for combined plots
#' @export
plot.patterns <- function(x,
                          type = "combined",
                          scaled = FALSE,
                          rawdata = NULL,
                          scale = c("standardized", "variance", "raw"),
                          scale_y = c("fixed", "free"),
                          events = NULL,
                          overlay = NULL,
                          overlay_color = "grey30",
                          overlay_fill = NA,
                          ...) {

  type <- match.arg(type, c("combined", "eofs", "amplitudes"))
  scale <- match.arg(scale)
  scale_y <- match.arg(scale_y)

  if(type == "combined") {
    # Smart combined plot: EOFs above, PCs below with KISS layout logic
    if(!requireNamespace("patchwork", quietly = TRUE)) {
      warning("patchwork package needed for combined plots. Install with: install.packages('patchwork')\nShowing EOFs only.")
      return(.plot_eofs_internal(x, scaled = scaled, rawdata = rawdata, overlay = overlay, overlay_color = overlay_color, overlay_fill = overlay_fill))
    } else {
      # Simple smart layout - covers 80% of use cases perfectly
      layout <- .simple_smart_layout(x$k)

      # Create EOF plot with smart layout
      p_eofs <- .plot_eofs_internal(x, scaled = scaled, rawdata = rawdata, layout = layout$eof, overlay = overlay, overlay_color = overlay_color, overlay_fill = overlay_fill)

      # Create amplitude plot with matching layout
      p_amps <- .plot_amplitudes_internal(x, scale = scale, scale_y = scale_y, events = events, layout = layout$pc)

      # Combine with patchwork using top-over-bottom layout and alignment tricks
      return(p_eofs /
             p_amps +
             patchwork::plot_layout(
               heights = layout$heights,
               guides = "collect",
               axis_titles = "collect"
             ) +
             patchwork::plot_annotation(
               title = glue::glue("EOF Analysis: {glue::glue_collapse(x$names, sep = ', ')}"),
               subtitle = glue::glue("k = {x$k} modes | {ifelse(x$scaled, 'scaled', 'unscaled')} | {ifelse(x$rotate, 'rotated', 'unrotated')}")
             ) &
             ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5)))
    }
  } else if(type == "eofs") {
    return(.plot_eofs_internal(x, scaled = scaled, rawdata = rawdata, overlay = overlay, overlay_color = overlay_color, overlay_fill = overlay_fill))
  } else if(type == "amplitudes") {
    return(.plot_amplitudes_internal(x, scale = scale, scale_y = scale_y, events = events))
  }
}

#' Simple smart layout function - KISS principle
#' @keywords internal
.simple_smart_layout <- function(k) {
  if(k <= 3) {
    # Single row: works for most cases, clean alignment
    list(eof = list(nrow = 1), pc = list(nrow = 1), heights = c(5, 1))
  } else {
    # k >= 4: Use square-ish grid
    ncol <- ceiling(sqrt(k))
    list(eof = list(ncol = ncol), pc = list(ncol = ncol), heights = c(5, 1))
  }
}

#' Internal function for EOF plotting
#' @keywords internal
.plot_eofs_internal <- function(x, scaled = FALSE, rawdata = NULL, layout = NULL,
                                 overlay = NULL, overlay_color = "grey30", overlay_fill = NA) {
  # Use provided layout or let facet_wrap choose defaults
  facet_args <- if(!is.null(layout)) layout else list()

  # Build overlay layer if provided
  overlay_layer <- if(!is.null(overlay)) {
    ggplot2::geom_sf(data = overlay, fill = overlay_fill, color = overlay_color, inherit.aes = FALSE)
  } else {
    NULL
  }

  if(scaled) {
    if(is.null(rawdata)) {
      rlang::abort("rawdata must be provided when scaled = TRUE for correlation calculation", class = "tidyEOF_missing_rawdata")
    }
    ggplot2::ggplot() +
      stars::geom_stars(data = get_correlation(rawdata, x)) +
      overlay_layer +
      do.call(ggplot2::facet_wrap, c(list(~PC), facet_args)) +
      ggplot2::scale_fill_distiller(palette = 'RdBu', na.value = NA, limits = c(-1, 1)) +
      ggplot2::coord_sf() +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "right") +
      ggplot2::labs(fill = "Correlation")
  } else {
    ggplot2::ggplot() +
      stars::geom_stars(data = x$eofs) +
      overlay_layer +
      do.call(ggplot2::facet_wrap, c(list(~PC), facet_args)) +
      scico::scale_fill_scico(palette = 'vik', midpoint = 0, na.value = NA) +
      ggplot2::coord_sf() +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "right") +
      ggplot2::labs(fill = "Loading")
  }
}

#' Internal function for amplitude plotting
#' @keywords internal
.plot_amplitudes_internal <- function(x, scale = "standardized", scale_y = "fixed", events = NULL, layout = NULL) {

  # Use provided layout or let facet_wrap choose defaults
  facet_args <- if(!is.null(layout)) layout else list()


  # Get PC names directly from amplitudes to ensure matching

  pc_names <- names(x$amplitudes)[-1]  # exclude 'time' column

  # Calculate scaling factors based on method
  if(scale == "variance") {
    # Multiply by sqrt(eigenvalue) to show variance contribution
    eigs <- x$eigenvalues %>%
      dplyr::filter(PC <= x$k) %>%
      dplyr::select(PC, std.dev) %>%
      dplyr::mutate(PC = pc_names[PC])

    amps <- x$amplitudes %>%
      tidyr::pivot_longer(-time, names_to = 'PC', values_to = 'amplitude') %>%
      dplyr::left_join(eigs, by = 'PC') %>%
      dplyr::mutate(amplitude = amplitude * std.dev)
  } else if(scale == "raw") {
    # Use EOF loadings to convert back to original units
    eigs <- split(x$eofs) %>%
      as_tibble() %>%
      dplyr::summarise(dplyr::across(dplyr::starts_with('PC'), ~sqrt(sum(.x^2, na.rm = TRUE)))) %>%
      tidyr::pivot_longer(dplyr::everything(), names_to = 'PC', values_to = 'scale_factor')

    amps <- x$amplitudes %>%
      tidyr::pivot_longer(-time, names_to = 'PC', values_to = 'amplitude') %>%
      dplyr::left_join(eigs, by = 'PC') %>%
      dplyr::mutate(amplitude = amplitude * scale_factor)
  } else {
    # Default: standardized (divide by sdev to get unit variance)
    eigs <- x$eigenvalues %>%
      dplyr::filter(PC <= x$k) %>%
      dplyr::select(PC, std.dev) %>%
      dplyr::mutate(PC = pc_names[PC])

    amps <- x$amplitudes %>%
      tidyr::pivot_longer(-time, names_to = 'PC', values_to = 'amplitude') %>%
      dplyr::left_join(eigs, by = 'PC') %>%
      dplyr::mutate(amplitude = amplitude / std.dev)
  }

  p <- ggplot2::ggplot(amps, ggplot2::aes(time, amplitude)) +
    ggplot2::geom_line() +
    do.call(ggplot2::facet_wrap, c(list(~PC), facet_args)) +
    ggplot2::labs(y = "Amplitude", x = NULL)

  if(scale_y == "fixed") {
    p <- p + ggplot2::coord_cartesian(ylim = range(amps$amplitude, na.rm = TRUE))
  }

  if(!is.null(events)) {
    p <- p + ggplot2::geom_vline(xintercept = events, linetype = 2, alpha = 0.5)
  }

  p + ggplot2::theme_bw()
}
