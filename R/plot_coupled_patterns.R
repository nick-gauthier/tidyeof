#' Plot coupled patterns diagnostics
#'
#' Provides visualization helpers for `coupled_patterns` objects. The default
#' `type = "combined"` mirrors the patterns plotting workflow by displaying the
#' predictor and response spatial patterns alongside their canonical variate time
#' series. Additional types include canonical correlation bars, standalone
#' canonical variate panels, or direct access to the underlying predictor / response
#' pattern plots.
#'
#' @param x A `coupled_patterns` object.
#' @param type Plot type: one of `"combined"`, `"correlations"`,
#'   `"canonical"`, `"predictor"`, or `"response"`.
#' @param side When `type = "canonical"`, choose canonical variates from the
#'   predictor side, response side, or both (values: `"predictor"`,
#'   `"response"`, `"both"`). Ignored for other plot types.
#' @param data Optional amplitudes or patterns for the canonical variate plots.
#'   For `side = "both"`, a list with elements `predictor` and `response` may be
#'   supplied. Defaults to the training patterns stored in `x` when omitted.
#' @param k Number of canonical modes to display (defaults to all available).
#' @param scaled Logical, passed to the underlying pattern plots when relevant
#'   (defaults to `FALSE`).
#' @param ... Additional arguments forwarded to `plot.patterns()` for
#'   `type = "predictor"` or `type = "response"` calls.
#'
#' @return A ggplot object (or a patchwork object when `type = "combined"`).
#' @export
plot.coupled_patterns <- function(x,
                                  type = c("combined", "correlations", "canonical",
                                           "predictor", "response"),
                                  side = c("predictor", "response", "both"),
                                  data = NULL,
                                  k = NULL,
                                  scaled = FALSE,
                                  ...) {
  type <- match.arg(type)

  if (is.null(k)) {
    k <- x$k
  } else {
    k <- min(k, x$k)
  }

  if (type == "correlations") {
    corr_df <- get_canonical_correlations(x, k)

    return(
      ggplot2::ggplot(corr_df, ggplot2::aes(x = factor(mode), y = correlation)) +
        ggplot2::geom_col(fill = "#3B8BC4") +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", correlation)),
                           vjust = -0.4, size = 3) +
        ggplot2::labs(
          x = "Canonical mode",
          y = "Correlation",
          title = "Canonical correlations"
        ) +
        ggplot2::ylim(0, max(corr_df$correlation) * 1.1) +
        ggplot2::theme_bw()
    )
  }

  if (type %in% c("predictor", "response")) {
    target <- if (type == "predictor") x$predictor_patterns else x$response_patterns
    return(plot(target, scaled = scaled, ...))
  }

  side <- match.arg(side)

  make_canonical_plot <- function(side, data, k) {
    if (side == "both") {
      predictor_data <- if (is.list(data) && "predictor" %in% names(data)) data$predictor else x$predictor_patterns
      response_data  <- if (is.list(data) && "response"  %in% names(data)) data$response  else x$response_patterns

      cv_pred <- get_canonical_variables(x, data = predictor_data, type = "predictor", k = k) %>%
        dplyr::mutate(side = "Predictor")
      cv_resp <- get_canonical_variables(x, data = response_data, type = "response", k = k) %>%
        dplyr::mutate(side = "Response")

      cv_long <- dplyr::bind_rows(cv_pred, cv_resp) %>%
        tidyr::pivot_longer(-c(time, side), names_to = "CV", values_to = "value")

      ggplot2::ggplot(cv_long, ggplot2::aes(x = time, y = value, color = side)) +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(~CV, scales = "free_y") +
        ggplot2::scale_color_manual(values = c(Predictor = "#4C9F38", Response = "#E17D30")) +
        ggplot2::labs(
          x = NULL,
          y = "Canonical variate",
          color = NULL,
          title = glue::glue("Canonical variates ({k} mode{ifelse(k == 1, '', 's')})")
        ) +
        ggplot2::theme_bw()
    } else {
      data_side <- if (is.null(data)) {
        if (side == "predictor") x$predictor_patterns else x$response_patterns
      } else {
        data
      }

      cv_tbl <- get_canonical_variables(x, data = data_side, type = side, k = k)
      cv_long <- cv_tbl %>%
        tidyr::pivot_longer(-time, names_to = "CV", values_to = "value")

      colour <- if (side == "predictor") "#4C9F38" else "#E17D30"

      ggplot2::ggplot(cv_long, ggplot2::aes(x = time, y = value)) +
        ggplot2::geom_line(color = colour) +
        ggplot2::facet_wrap(~CV, scales = "free_y") +
        ggplot2::labs(
          x = NULL,
          y = "Canonical variate",
          title = glue::glue("{stringr::str_to_title(side)} canonical variates")
        ) +
        ggplot2::theme_bw()
    }
  }

  if (type == "canonical") {
    return(make_canonical_plot(side, data, k))
  }

  # Combined view: predictor and response spatial patterns plus canonical series
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("patchwork package needed for combined coupled plots. Install with: install.packages('patchwork')")
    return(list(
      predictor = plot(x$predictor_patterns, type = "combined", scaled = scaled, ...),
      response  = plot(x$response_patterns,  type = "combined", scaled = scaled, ...),
      canonical = make_canonical_plot("both", data, k)
    ))
  }

  predictor_eofs <- plot(x$predictor_patterns, type = "eofs", scaled = scaled, ...)
  response_eofs  <- plot(x$response_patterns,  type = "eofs", scaled = scaled, ...)
  canonical_plot <- make_canonical_plot("both", data, k)

  (predictor_eofs /
     response_eofs /
     canonical_plot) +
    patchwork::plot_layout(heights = c(3, 3, 2), guides = "collect") +
    patchwork::plot_annotation(
      title = glue::glue("Coupled patterns overview (k = {k})"),
      subtitle = "Top: predictor patterns | middle: response patterns | bottom: canonical variates"
    ) &
    ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5))
}
