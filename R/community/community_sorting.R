# Sorting factor (roadmap.md, 25 November 2026): species with the same
# climate centre sort among locations (soil, topography), so the species
# that grow together are more alike than the location model's species SD
# implies; the location effect only takes up shifts shared by species
# measured in the same cells. The factor k scales the species-spread term of
# the community variance (community_statistics_climate()).
#
# k is estimated from the species effects u-hat of the location model and
# each cell's GBIF species: in cells with at least `min_species` measured
# species of a growth form, the variance of u-hat among them, against what
# the model expects for independent species, the mean over them of
# sd_u(s)^2 less the estimate's own error variance (u-hat is shrunk). The
# same ratio over all species of the form is a check and should be about 1.
# k is capped at 1 (co-occurring species are not expected to be less alike
# than random ones; above 1 is noise) and set to 1 with too few cells.

#' @param fits list of fit_location_target() results (with centres).
#' @param species output of `cell_species()` (GBIF names x cells).
#' @param name_map output of `gbif_to_species()`.
#' @param min_species measured species a cell needs.
#' @param min_cells cells a growth form needs for k to be estimated.
#' @return tibble(trait, woodiness_group, n_species, ratio_all, cells,
#'   ratio_within, k).
estimate_sorting_factors <- function(fits, species, name_map, min_species = 5, min_cells = 20) {
  cell_sp <- species |>
    dplyr::inner_join(name_map, by = "gbif_taxon_name") |>
    dplyr::distinct(cell_x, cell_y, taxon_name, woodiness_group)
  purrr::map(purrr::compact(fits), function(m) {
    ce <- m$centres |>
      dplyr::mutate(woodiness_group = as.character(woodiness_group)) |>
      dplyr::left_join(m$species_variance, by = "taxon_name") |>
      dplyr::transmute(taxon_name, woodiness_group, u = g_0,
                       expected = pmax(species_variance - log_se^2, 0))
    overall <- ce |>
      dplyr::group_by(woodiness_group) |>
      dplyr::summarise(n_species = dplyr::n(), ratio_all = stats::var(u) / mean(expected))
    within <- cell_sp |>
      dplyr::inner_join(ce, by = c("taxon_name", "woodiness_group")) |>
      dplyr::group_by(cell_x, cell_y, woodiness_group) |>
      dplyr::filter(dplyr::n() >= min_species) |>
      dplyr::summarise(.groups = "drop", v = stats::var(u), e = mean(expected)) |>
      dplyr::group_by(woodiness_group) |>
      dplyr::summarise(cells = dplyr::n(), ratio_within = mean(v) / mean(e))
    dplyr::left_join(overall, within, by = "woodiness_group") |>
      dplyr::mutate(trait = m$trait, .before = 1,
                    cells = dplyr::coalesce(cells, 0L),
                    k = dplyr::if_else(cells >= min_cells, pmin(ratio_within, 1), 1))
  }) |> purrr::list_rbind()
}

# The sorting factor of one trait and growth form (1 without one).
.sorting_k <- function(sorting, trait, group) {
  if (is.null(sorting)) return(1)
  k <- sorting$k[sorting$trait == trait & sorting$woodiness_group == group]
  if (length(k) == 1 && is.finite(k)) k else 1
}
