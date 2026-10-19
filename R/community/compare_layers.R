# The three community products against each other, cell by cell (#28): the
# community mean of all species on the log scale, in cells where both
# products have one. Modelled cells flagged `extrapolated` are left out.

community_layer_pair_labels <- c(
  presence_modelled = "Presence (x) vs modelled (y)",
  local_presence = "Local records (x) vs presence (y)",
  local_modelled = "Local records (x) vs modelled (y)"
)

#' Cell-level pairs of log community means.
#'
#' @param local,presence,modelled the three products (`community_local`,
#'   `community_presence`, `community_modelled`).
#' @param group woodiness group compared; the modelled product has `all` only.
#' @return tibble(pair, trait_name, cell_x, cell_y, x, y, log_x, log_y).
community_layer_pairs <- function(local, presence, modelled, group = "all") {
  pick <- function(d) {
    d |>
      dplyr::filter(woodiness_group == group, !is.na(community_mean)) |>
      dplyr::select(trait_name, cell_x, cell_y, x, y, log_mean = community_mean) |>
      dplyr::mutate(log_mean = log(log_mean))
  }
  m <- dplyr::filter(modelled, !extrapolated)
  pair <- function(a, b, name) {
    dplyr::inner_join(pick(a), dplyr::select(pick(b), -x, -y),
                      by = c("trait_name", "cell_x", "cell_y"),
                      suffix = c("_x", "_y")) |>
      dplyr::transmute(pair = name, trait_name, cell_x, cell_y, x, y,
                       log_x = log_mean_x, log_y = log_mean_y)
  }
  dplyr::bind_rows(
    pair(presence, m, "presence_modelled"),
    pair(local, presence, "local_presence"),
    pair(local, m, "local_modelled")
  )
}

#' Agreement per pair and trait: cells, correlation of the log community
#' means, median and spread of their difference (y minus x, log units).
summarise_layer_pairs <- function(pairs) {
  pairs |>
    dplyr::group_by(pair, trait_name) |>
    dplyr::summarise(.groups = "drop",
                     n_cells = dplyr::n(),
                     cor_log = stats::cor(log_x, log_y),
                     median_diff_log = stats::median(log_y - log_x),
                     sd_diff_log = stats::sd(log_y - log_x),
                     sd_log_x = stats::sd(log_x),
                     sd_log_y = stats::sd(log_y)) |>
    dplyr::mutate(pair = factor(pair, names(community_layer_pair_labels))) |>
    dplyr::arrange(pair, trait_name) |>
    dplyr::mutate(pair = as.character(pair))
}

#' One pair: a scatter per trait, with the 1:1 line, cells and correlation.
plot_layer_pair <- function(pairs, summary, which, units) {
  d <- dplyr::filter(pairs, pair == which) |>
    dplyr::add_count(trait_name) |>
    dplyr::mutate(alpha = pmin(0.8, pmax(0.08, 2000 / n)))
  s <- dplyr::filter(summary, pair == which)
  lab <- function(tr) {
    u <- units$unit[units$trait_name == tr]
    if (length(u) == 0 || is.na(u)) tr else paste0(tr, " (", u, ")")
  }
  d$trait_label <- vapply(d$trait_name, lab, character(1))
  s$trait_label <- vapply(s$trait_name, lab, character(1))
  s$text <- sprintf("r = %.2f\nn = %s", s$cor_log, scales::label_comma()(s$n_cells))
  # label at each panel's top left, in data units (log axes cannot place Inf)
  corner <- d |>
    dplyr::group_by(trait_name) |>
    dplyr::summarise(lx = exp(min(log_x)), ly = exp(max(log_y)))
  s <- dplyr::left_join(s, corner, by = "trait_name")
  ggplot2::ggplot(d, ggplot2::aes(exp(log_x), exp(log_y))) +
    ggplot2::geom_abline(linetype = "dashed", colour = "#6b6a66", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(alpha = alpha), size = 0.8, stroke = 0,
                        colour = "#1c5cab") +
    ggplot2::geom_text(data = s, ggplot2::aes(x = lx, y = ly, label = text),
                       inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.6,
                       colour = "#2b2b28") +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~trait_label, scales = "free", ncol = 3) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0("Community mean, all species: ",
                                 community_layer_pair_labels[[which]]),
                  subtitle = paste0("One point per 10 km cell with a value in both; dashed: 1:1",
                                    if (grepl("modelled", which))
                                      "; modelled cells outside the fitted range left out")) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0),
                   plot.subtitle = ggplot2::element_text(colour = "#6b6a66"))
}

#' Map of modelled minus presence (log units) per trait, all species.
plot_layer_difference <- function(pairs, limit = 1.5) {
  d <- dplyr::filter(pairs, pair == "presence_modelled") |>
    dplyr::mutate(difference = log_y - log_x)
  diverging <- c("#7a2c10", "#eb6834", "#f6a57f", "#f2f1ee", "#86b6ef", "#3987e5", "#0d366b")
  map_cells(d, "difference", diverging,
            paste0("Modelled minus presence community mean, all species (log units; ",
                   "blue: modelled higher; squished at ±", limit, ")"),
            limits = c(-limit, limit)) +
    ggplot2::facet_wrap(~trait_name, ncol = 3) +
    ggplot2::theme(legend.key.height = grid::unit(1.5, "cm"))
}

#' Write the comparison: summary CSV, one scatter figure per pair, and the
#' difference maps. Returns the file paths.
write_layer_comparison <- function(pairs, units, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  summary <- summarise_layer_pairs(pairs)
  csv <- file.path(dir, "community_layer_agreement.csv")
  readr::write_csv(dplyr::mutate(summary, dplyr::across(where(is.double), \(v) round(v, 3))), csv)
  figs <- c(
    unlist(lapply(names(community_layer_pair_labels), function(p) {
      save_figure(plot_layer_pair(pairs, summary, p, units),
                  file.path(dir, "figures", paste0("pair_", p)), width = 9, height = 9)
    }), use.names = FALSE),
    save_figure(plot_layer_difference(pairs),
                file.path(dir, "figures", "difference_modelled_minus_presence"),
                width = 10, height = 10)
  )
  c(csv, unname(figs))
}
