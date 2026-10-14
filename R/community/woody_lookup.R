# Provisional woody / non-woody lookup for the community layers (#22).
#
# PROVISIONAL, pending the decision in #19 (applied in #27). The rule:
#
#   1. The modal `woodiness_detailed` value: woody -> woody; herbaceous ->
#      non-woody; the intermediate categories (semi_woody, woody_base,
#      woody_like_inflorescence, woody_like_stem, woody_root) -> non-woody,
#      flagged `intermediate`. A tie between values of different classes is a
#      conflict and falls through to step 2.
#   2. Otherwise the growth-form keywords of #19's candidate rule: a taxon with
#      any growth form containing herb, fern, tussock, geophyte, lycophyte,
#      basal_large, graminoid or hummock (so also herbaceous climbers and any
#      mixed form) -> non-woody; any other recorded growth form -> woody.
#   3. Otherwise unassigned.
#
# Inputs are the categorical rows of the geometric observed means, so taxa
# are already rolled up to species.

woodiness_detailed_class <- c(
  woody = "woody",
  herbaceous = "non_woody",
  semi_woody = "non_woody",
  woody_base = "non_woody",
  woody_like_inflorescence = "non_woody",
  woody_like_stem = "non_woody",
  woody_root = "non_woody"
)

non_woody_growth_form_keywords <- c("herb", "fern", "tussock", "geophyte",
                                    "lycophyte", "basal_large", "graminoid",
                                    "hummock")

#' Provisional woody / non-woody class per taxon.
#'
#' @param observed_means output of `estimate_observed_means_geometric()`.
#' @return one row per taxon with any woodiness or growth-form record:
#'   `woodiness_group` (woody, non_woody, unassigned), `source`
#'   (woodiness_detailed, plant_growth_form, none), `intermediate` and
#'   `conflict` flags, the inputs used, and `provisional = TRUE`.
build_woody_lookup <- function(observed_means) {
  cat_rows <- dplyr::filter(observed_means, type == "categorical")
  values_of <- function(trait, column) {
    cat_rows |>
      dplyr::filter(trait_name == trait) |>
      dplyr::select(taxon_name, value_mode, value_count) |>
      stats::setNames(c("taxon_name", column, paste0(column, "_counts")))
  }

  wd <- values_of("woodiness_detailed", "woodiness_detailed_mode")
  gf <- values_of("plant_growth_form", "growth_form_mode")
  wo <- values_of("woodiness", "woodiness_mode")

  classes_in <- function(values) {
    v <- strsplit(values, " ", fixed = TRUE)
    lapply(v, \(x) unique(stats::na.omit(woodiness_detailed_class[x])))
  }
  counted_values <- function(counts) {
    lapply(strsplit(dplyr::coalesce(counts, ""), "; ", fixed = TRUE),
           \(x) sub(" \\([0-9]+\\)$", "", x))
  }

  taxa <- dplyr::bind_rows(
    dplyr::select(wd, taxon_name), dplyr::select(gf, taxon_name),
    dplyr::select(wo, taxon_name)
  ) |>
    dplyr::distinct() |>
    dplyr::left_join(wd, by = "taxon_name") |>
    dplyr::left_join(gf, by = "taxon_name") |>
    dplyr::left_join(wo, by = "taxon_name")

  mode_classes <- classes_in(dplyr::coalesce(taxa$woodiness_detailed_mode, ""))
  all_classes <- lapply(counted_values(taxa$woodiness_detailed_mode_counts),
                        \(x) unique(stats::na.omit(woodiness_detailed_class[x])))
  forms <- counted_values(taxa$growth_form_mode_counts)
  non_woody_form <- vapply(forms, \(x) any(vapply(
    non_woody_growth_form_keywords, \(k) any(grepl(k, x, fixed = TRUE)),
    logical(1))), logical(1))
  any_form <- vapply(forms, \(x) any(nzchar(x)), logical(1))

  from_wd <- vapply(mode_classes, \(x) if (length(x) == 1) x else NA_character_,
                    character(1))
  intermediate_mode <- vapply(
    strsplit(dplyr::coalesce(taxa$woodiness_detailed_mode, ""), " ", fixed = TRUE),
    \(x) any(x %in% setdiff(names(woodiness_detailed_class), c("woody", "herbaceous"))),
    logical(1))
  from_gf <- dplyr::case_when(non_woody_form ~ "non_woody",
                              any_form ~ "woody",
                              TRUE ~ NA_character_)

  taxa |>
    dplyr::mutate(
      woodiness_group = dplyr::coalesce(from_wd, from_gf, "unassigned"),
      source = dplyr::case_when(!is.na(from_wd) ~ "woodiness_detailed",
                                !is.na(from_gf) ~ "plant_growth_form",
                                TRUE ~ "none"),
      intermediate = !is.na(from_wd) & intermediate_mode,
      # records of the taxon carry both woody and non-woody woodiness_detailed
      # values (whether or not the mode resolves it)
      conflict = vapply(all_classes, length, integer(1)) > 1,
      provisional = TRUE,
      .after = taxon_name
    ) |>
    dplyr::arrange(taxon_name)
}

#' Counts per class and source, and the number of conflicting taxa.
summarise_woody_lookup <- function(lookup) {
  dplyr::bind_rows(
    dplyr::count(lookup, woodiness_group, source, name = "n_taxa"),
    tibble::tibble(woodiness_group = "(all)", source = "intermediate (flagged)",
                   n_taxa = sum(lookup$intermediate)),
    tibble::tibble(woodiness_group = "(all)", source = "conflicting values",
                   n_taxa = sum(lookup$conflict)),
    tibble::tibble(woodiness_group = "unassigned", source = "of which have a woodiness value",
                   n_taxa = sum(lookup$woodiness_group == "unassigned" &
                                  !is.na(lookup$woodiness_mode)))
  )
}

#' Write the lookup and its summary to `dir`.
#' @return the file paths.
write_woody_lookup <- function(lookup, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, c("woody_lookup_provisional.csv",
                            "woody_lookup_provisional_summary.csv"))
  readr::write_csv(lookup, paths[1], na = "")
  readr::write_csv(summarise_woody_lookup(lookup), paths[2])
  paths
}
