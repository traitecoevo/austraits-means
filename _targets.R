# ============================================================================
# targets pipeline for austraits-means
# ============================================================================
# Orchestration lives here; all logic stays in R/ functions (loaded via
# R/init.R). The two Quarto reports in reports/ are *reporting layers* rendered
# at the end of the pipeline (tar_quarto): they tar_load() results rather than
# compute.
#
# Run:   targets::tar_make()
# Graph: targets::tar_visnetwork()  /  targets::tar_outdated()
#
# TIERING ----------------------------------------------------------------------
# The heavy upstream artefacts (raw GBIF cleaning ~1.2GB, the 22k-species O(n^2)
# medoid build) are expensive and rarely rebuilt, and their builders live
# outside the pipeline (R/data/build_sources.R, run by R/scripts/build_gbif.R
# and R/scripts/build_medoids.R). They are tracked here as *file inputs* — targets
# watches their hashes and invalidates downstream when they change, but never
# recomputes them as part of tar_make(). Regenerate them deliberately with the
# builders in R/data/build_sources.R, then re-run tar_make().
#
# Everything from `unified` forward (the modelling pipeline) is fully computed
# and incrementally rebuilt by targets — this is where the manual `update`
# toggles and `if (file.exists)` guards used to live.
# ============================================================================

library(targets)
library(tarchetypes)

# Render/parallel backends used implicitly (tar_quarto needs `quarto`; the
# optional crew controller below needs `crew`). Referenced here so renv's
# dependency scanner tracks them in renv.lock.
requireNamespace("quarto", quietly = TRUE)
if (FALSE) requireNamespace("crew", quietly = TRUE)

source("R/init.R")

tar_option_set(
  packages = c(
    "tidyverse", "arrow", "terra", "sf",
    "glmmTMB", "broom.mixed", "marginaleffects",
    "austraits", "APCalign", "furrr"
  ),
  format = "rds"  # built-in; switch to "qs" if the qs2 package is installed
  # Cross-target parallelism (optional): uncomment to fan traits across workers.
  # TMB is fast and brms manages its own chains, so this is off by default for
  # a robust first run. Enable once the pipeline is established.
  # , controller = crew::crew_controller_local(
  #     workers = max(1, parallelly::availableCores(methods = "system") - 2))
)

# ---- configuration ---------------------------------------------------------
out_dir <- "outputs/traits"

# austraits.build development export the trait records come from (height
# included, via the AVH_2026 dataset). How it was built:
# data-raw/austraits/README.md.
austraits_dev_rds <- "data/austraits/austraits-7.0.0.9000-d535de9.rds"

# The geometric observed means that austraits.build's `post_export()` wrote for
# that same build (austraits.build#899). This repository reads them rather
# than computing its own; how to make the file: data-raw/austraits/README.md.
austraits_dev_observed_means <-
  "data/austraits/austraits-7.0.0.9000-d535de9-observed-means.parquet"

# APC/APNI release that trait-record names and GBIF names are both aligned to
# for the crosswalk (an APCalign release tag; #13). Changing it rebuilds the
# crosswalk and everything downstream.
apc_version <- "2026-07-30"

# Static trait list for per-trait branching (must be known at pipeline-build
# time); also the traits read from the austraits.build export.
trait_names <- c(
  "plant_height", "leaf_mass_per_area", "leaf_area", "leaf_N_per_dry_mass",
  "leaf_N_per_area", "leaf_P_per_area", "leaf_P_per_dry_mass",
  "seed_dry_mass", "wood_density"
)

# Opt-in BRMS backend (full posterior; hours per large trait). Default off.
run_brms <- FALSE

# Per-trait fan-out (TMB): one branch per trait — fit + export centres, plus a
# one-row summary and a figure. Only the changed trait refits. Captured here so
# the tar_combine() steps below can reference the mapped branch targets.
mapped <- tar_map(
  values = list(trait = trait_names),
  names = "trait",

  # The earlier pooled model, kept as a diagnostic; it writes no files.
  tar_target(tmb_fit,
             fit_trait_target(unified, trait, backend = "tmb",
                              n_min = 5, out_dir = NULL)),

  # The location model (R/model/spread_rtmb.R; roadmap.md, 25 November
  # 2026): the one model for the modelled means and, when built, the
  # community statistics. Free fit, then shrunk, with each species' modelled
  # mean. Height takes about 25 minutes.
  tar_target(location_fit, fit_location_target(unified_cells, trait, woody_lookup)),
  # its species means as the trait's working export,
  # <trait>_climate_response.csv.gz and _meta.csv
  tar_target(trait_export, export_location_centres(location_fit, out_dir),
             format = "file"),

  # Woody and non-woody coefficients, locations pooled: for the comparison
  # below only.
  tar_target(woody_fit, fit_woodiness_model(unified, trait, woody_lookup)),

  # The pooled models' one role: what the location model fixes, step by
  # step (pooled -> locations pooled -> location model;
  # R/process/compare_pooled_location.R).
  tar_target(pooled_vs_location,
             compare_pooled_location(location_fit, tmb_fit, woody_fit)),

  # Per-trait figures, from the location model, woody and non-woody apart:
  # the species' modelled means against climate, over the model's
  # distribution of species means (line and 50% / 90% bands) ...
  tar_target(trait_fig, {
    if (is.null(location_fit)) character(0) else {
      p <- plot_trait_centres(location_fit$centres, trait,
                              x = "precip", log_y = TRUE, fit = location_fit)
      save_figure(p, file.path(out_dir, "figures",
                               clean_species_filename(trait)),
                  width = 7, height = 5.5)
    }
  }, format = "file"),

  # ... and the within- vs between-species lines, as in the earlier height
  # report, one row per growth form.
  tar_target(trait_wb_fig, {
    if (is.null(location_fit)) character(0) else {
      save_figure(plot_within_between(location_fit, unit = location_fit$data$unit[1]),
                  file.path(out_dir, "figures",
                            paste0(clean_species_filename(trait), "_within_between")),
                  width = 13, height = 6 * length(location_fit$groups))
    }
  }, format = "file")
)

list(

  # ===== Tier B: tracked file inputs (raw + seeded heavy artefacts) =========
  # Raw sources (tracked for provenance; consumed by the builders in
  # R/data/build_sources.R when regenerating the seeded artefacts below).
  tar_target(climate_australia_tif, "data/australia.tif", format = "file"),
  tar_target(climate_bio1_tif, "data/climate/wc2.1_2.5m_bio_1.tif",
             format = "file"),
  tar_target(climate_bio12_tif, "data/climate/wc2.1_2.5m_bio_12.tif",
             format = "file"),

  # Seeded heavy artefacts — produced upstream, reused as inputs (not rebuilt).
  tar_target(gbif_cleaned_file,
             "data-processed/data_gbif_cleaned_0027765-260928105237408.parquet",
             format = "file"),
  tar_target(medoids_file,
             "data-processed/climate_centres_0027765-260928105237408.parquet",
             format = "file"),
  tar_target(austraits_dev_file, austraits_dev_rds, format = "file"),
  tar_target(austraits_dev_observed_means_file, austraits_dev_observed_means,
             format = "file"),

  # ===== Tier A: computed modelling pipeline ================================
  # Version string and commit of the austraits.build export, carried into the
  # product metadata.
  tar_target(austraits_dev_info, austraits_build_info(austraits_dev_file)),

  # Field records of the analysis traits with climate attached, then the
  # value-type rule (per-dataset exceptions: dataset_value_type_rules()).
  tar_target(austraits_climate,
             build_austraits_climate(
               austraits_dev_file, traits = trait_names,
               rasters = load_climate_rasters(
                 vars = basename(c(climate_bio1_tif, climate_bio12_tif)),
                 climate_dir = dirname(climate_bio1_tif),
                 australia_tif = climate_australia_tif))),
  # Infraspecific taxa are rolled up to species (roadmap.md).
  tar_target(records_all,
             roll_up_to_species(load_austraits_records(austraits_climate),
                                dev_taxa)),

  # Name crosswalk: every taxon name in the build and every GBIF name, both
  # aligned to the pinned APC release.
  tar_target(dev_taxa,
             readRDS(austraits_dev_file)$taxa |>
               dplyr::select(taxon_name, taxon_rank, binomial)),
  tar_target(gbif_names, gbif_taxon_names(gbif_cleaned_file)),
  tar_target(name_crosswalk,
             build_name_crosswalk(
               taxon_names = dev_taxa$taxon_name,
               gbif_names = gbif_names, apc_version = apc_version)),

  # Per-species medoid lookup (crosswalk-bridged).
  tar_target(medoid_lookup,
             build_medoid_lookup(medoids = medoids_file,
                                 crosswalk = name_crosswalk)),
  tar_target(name_matching,
             summarise_name_matching(records_all, name_crosswalk, gbif_names,
                                     medoid_lookup)),

  # The linchpin: unified record-level table, medoid-centred.
  tar_target(unified, build_dataset(records_all, medoid_lookup)),
  # The same records keyed to the location model's cells.
  tar_target(unified_cells,
             add_location_key(unified, "cell", cell_km = location_cell_km())),

  # Medoid reasonableness check (vs the GBIF occurrence cloud). Reads the large
  # GBIF parquet, so depends on gbif_cleaned_file.
  tar_target(medoid_gbif_check,
             check_medoid_vs_gbif(medoid_lookup = medoid_lookup,
                                  gbif_path = gbif_cleaned_file,
                                  n_sample = 300, seed = 1)),
  tar_target(medoid_check_fig, {
    p <- plot_medoid_check(medoid_gbif_check)  # build only; save below for paths
    save_figure(p, "outputs/checks/medoid_gbif_check", width = 7, height = 6)
  }, format = "file"),

  # ---- per-trait fan-out (TMB) --------------------------------------------
  # The branch targets (tmb_fit_<trait>, location_fit_<trait>, woody_fit_<trait>,
  # pooled_vs_location_<trait>, trait_fig_<trait>, ...),
  # defined as `mapped` above.
  mapped,

  # Lists of fits across traits.
  tar_combine(tmb_fits, mapped[["tmb_fit"]], command = list(!!!.x)),
  tar_combine(location_fits, mapped[["location_fit"]], command = list(!!!.x)),

  # Fixed-effects summary of the location model (one row per trait x growth
  # form).
  tar_target(trait_summary, summarise_location_fits(location_fits)),
  tar_target(trait_summary_file, {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    p <- file.path(out_dir, "trait_model_summary.csv")
    readr::write_csv(trait_summary, p)
    p
  }, format = "file"),

  # Pooled models against the location model, across traits.
  tar_combine(pooled_location_comparison, mapped[["pooled_vs_location"]],
              command = list(!!!.x)),
  tar_target(pooled_location_comparison_files,
             write_pooled_location_comparison(pooled_location_comparison),
             format = "file"),

  # Data basis of the fits (build version, APC release, GBIF extract, medoid
  # build), written next to the exports.
  tar_target(data_basis_info,
             data_basis(austraits_dev_info, apc_version, gbif_cleaned_file,
                        medoids_file)),
  tar_target(data_basis_file, {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    p <- file.path(out_dir, "data_basis.json")
    jsonlite::write_json(data_basis_info, p, auto_unbox = TRUE, pretty = TRUE)
    p
  }, format = "file"),

  # Shifts against the June 2026 outputs, snapshotted before the refit (#14):
  # the product then (pooled model) against the product now (location model).
  tar_target(baseline_dir, "docs/baselines/2026-06", format = "file"),
  tar_target(baseline_comparison,
             compare_with_baseline(location_fits, unified, baseline_dir)),
  tar_target(baseline_comparison_files,
             write_baseline_comparison(baseline_comparison), format = "file"),

  # ---- Modelled means product (all traits, product schema; #17) ----------
  # From the location model since 25 November 2026; the pooled TMB fits
  # above remain as diagnostics (fixed-effect summary, figures, baseline).
  tar_target(modelled_means, modelled_means_table(location_fits)),
  tar_target(modelled_means_files,
             write_species_means(
               modelled_means, modelled_means_dictionary(),
               dir = file.path("outputs/species_means",
                               austraits_dev_info$austraits_version),
               stem = "modelled_means",
               metadata = modelled_means_metadata(data_basis_info)),
             format = "file"),

  # ===== Observed means, portal method: species means 7.0.0 ================
  # Reproduces the AusTraits data portal's species averages for release 7.0.0
  # (R/means/observed_means_portal.R, docs/observed-means-portal-method.md).
  tar_target(austraits_release_file, "data/austraits/austraits-7.0.0.rds",
             format = "file"),
  tar_target(portal_core_traits_file,
             "data-raw/portal/trait_groups_for_portal.csv", format = "file"),
  tar_target(portal_records,
             prepare_portal_records(
               read_austraits_release(austraits_release_file, "7.0.0"),
               read_portal_core_traits(portal_core_traits_file))),
  tar_target(observed_means_portal,
             estimate_observed_means_portal(portal_records)),
  tar_target(observed_means_portal_files,
             write_species_means(
               label_estimates(
                 add_portal_reliability(observed_means_portal,
                                        portal_reliability(portal_records)),
                 "observed_mean", "portal_arithmetic"),
               dplyr::bind_rows(observed_means_portal_dictionary(),
                                reliability_dictionary()),
               dir = "outputs/species_means/7.0.0",
               stem = "observed_means_portal_arithmetic",
               metadata = observed_means_portal_metadata()),
             format = "file"),

  # Reproduction check against the parquet the portal deployed (built 9 Feb
  # 2026; copied from the portal's git-ignored data build). Also asserted by
  # tests/testthat/test-observed-means-portal.R. Optional: without that file
  # these targets error and the rest of the pipeline carries on.
  tar_target(portal_deployed_file,
             "data/portal/austraits-7.0.0-full/austraits-species-averages.parquet",
             format = "file", error = "continue"),
  tar_target(portal_reproduction,
             compare_portal_reproduction(observed_means_portal,
                                         arrow::read_parquet(portal_deployed_file)),
             error = "continue"),
  tar_target(portal_reproduction_files,
             write_portal_reproduction(portal_reproduction), format = "file",
             error = "continue"),

  # Draft Zenodo deposit for species means 7.0.0 (#10): the files, README and
  # deposit metadata. Uploading is a separate, deliberate step
  # (R/scripts/create_zenodo_draft.R); nothing publishes. Retired 24 November
  # 2026 with the portal-method product (roadmap.md); kept, not maintained.
  tar_target(portal_method_doc, "docs/observed-means-portal-method.md",
             format = "file"),
  tar_target(zenodo_package_portal,
             prepare_zenodo_package(observed_means_portal_files,
                                    portal_method_doc,
                                    dir = "outputs/zenodo/7.0.0",
                                    version = "7.0.0",
                                    austraits_doi = "10.5281/zenodo.15718081"),
             format = "file"),

  # ===== Observed means, geometric method (development build; #15) =========
  # Published by austraits.build with the build, not computed here
  # (R/means/observed_means_release.R, docs/observed-means-geometric-method.md).
  # The woody lookup needs the woodiness and growth-form traits too.
  tar_target(observed_means_geometric,
             read_observed_means(austraits_dev_observed_means_file,
                                 traits = c(trait_names, "woodiness_detailed",
                                            "plant_growth_form"))),

  # ===== Community statistics (roadmap.md, 25 November 2026) ===============
  # For each cell of the location model's grid (50 km) and growth form, the
  # community mean and community range of the species GBIF records there,
  # from the location model (R/community/community_statistics.R). Replaces
  # the three 10 km community layers and their comparisons.
  # Woody / non-woody per taxon, by the rule of #19.
  tar_target(woody_lookup, build_woody_lookup(observed_means_geometric)),
  tar_target(community_dir,
             file.path("outputs/community", austraits_dev_info$austraits_version)),
  tar_target(woody_lookup_files, write_woody_lookup(woody_lookup, community_dir),
             format = "file"),
  tar_target(community_units, trait_units(observed_means_geometric, trait_names)),

  # Combined species means (modelled where available, else observed), a
  # species means export in its own right.
  tar_target(species_means_combined,
             combine_species_means(location_fits, observed_means_geometric)),
  tar_target(species_means_combined_files,
             write_species_means(
               species_means_combined, species_means_combined_dictionary(),
               dir = file.path("outputs/species_means",
                               austraits_dev_info$austraits_version),
               stem = "species_means_combined",
               metadata = species_means_combined_metadata(data_basis_info)),
             format = "file"),

  # The species of each cell (GBIF), with growth form and climate centre, and
  # the cell's climate.
  tar_target(gbif_presence, gbif_presence_cells(gbif_cleaned_file)),
  tar_target(gbif_species_map,
             gbif_to_species(gbif_names, name_crosswalk, dev_taxa,
                             woody_lookup$taxon_name)),
  tar_target(community_species,
             cell_species(gbif_presence, medoids_file, gbif_species_map, woody_lookup)),
  tar_target(cell_climate_grid,
             cell_climate(load_climate_rasters(
               vars = basename(c(climate_bio1_tif, climate_bio12_tif)),
               climate_dir = dirname(climate_bio1_tif),
               australia_tif = climate_australia_tif))),

  # The product: from climate alone. The co-occurring species' centres are
  # learned as smooth surfaces of climate from the GBIF cells, so the
  # statistics need only a location's climate (R/community/community_climate.R).
  tar_target(centre_summaries, cell_centre_summaries(community_species, cell_climate_grid)),
  tar_target(centre_surfaces, fit_centre_surfaces(centre_summaries)),
  tar_target(community_stats,
             community_statistics_climate(
               location_fits, predict_centre_surfaces(centre_surfaces, cell_climate_grid))),
  tar_target(community_stats_files,
             write_community_layers(community_stats, community_dir,
                                    community_statistics_metadata(data_basis_info),
                                    statistics = c("community_mean", "community_sd",
                                                   "range_lower", "range_upper")),
             format = "file"),
  tar_target(community_stats_maps,
             write_community_maps(community_stats, community_units, community_dir),
             format = "file"),

  # The check: the same statistics from each cell's own GBIF species
  # (community_statistics()), against the climate-alone ones.
  tar_target(community_stats_gbif,
             community_statistics(location_fits, cell_climate_grid, community_species)),
  tar_target(community_climate_check,
             compare_community_climate_gbif(community_stats, community_stats_gbif)),
  tar_target(community_climate_check_files,
             write_community_climate_check(community_climate_check,
                                           file.path(community_dir, "check")),
             format = "file"),

  # ===== Observed vs modelled means (#24) ===================================
  tar_target(observed_modelled,
             compare_observed_modelled(location_fits, observed_means_geometric, unified)),
  tar_target(portal_geometric,
             compare_portal_geometric(observed_means_portal,
                                      observed_means_geometric)),
  tar_target(means_comparison_files,
             write_means_comparison(observed_modelled, portal_geometric),
             format = "file"),

  # ===== Figure set (#25) ====================================================
  tar_target(figure_set,
             write_figure_set(observed_means_portal, observed_means_geometric,
                              modelled_means, location_fits, observed_modelled,
                              community_stats, means_comparison_files,
                              c(community_stats_maps, pooled_location_comparison_files)),
             format = "file"),

  # ---- reporting layers (rendered notebooks) ------------------------------
  # rendered next to the .qmd; code runs from the project root
  tar_quarto(report_dataset, "reports/dataset_construction.qmd"),
  tar_quarto(report_analysis, "reports/analysis.qmd")

) |>
  # ===== opt-in BRMS backend (full posterior) ==============================
  # Appended only when run_brms is TRUE. Height stays on TMB (too large for
  # NUTS). Fits are cached on disk under fits_dir (targets stores the summary
  # objects; the Stan fits live in data-processed/brms_fits/).
  (\(pipeline) if (!run_brms) pipeline else c(pipeline, list(
    tar_target(brms_traits, setdiff(sort(unique(unified$trait_name)),
                                    "plant_height")),
    tar_target(brms_fits,
               run_trait_models(unified, traits = brms_traits,
                                 backend = "brms", n_min = 5, out_dir = NULL,
                                 fits_dir = "data-processed/brms_fits",
                                 threads = max(1, floor(
                                   parallelly::availableCores(
                                     methods = "system") / 4)),
                                 iter = 2000, warmup = 1000,
                                 refresh = 250, silent = 1)),
    tar_target(brms_diag, summarise_trait_diagnostics(brms_fits)),
    tar_target(brms_diag_file, {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      p <- file.path(out_dir, "brms_diagnostics.csv")
      readr::write_csv(brms_diag, p); p
    }, format = "file")
  )))()
