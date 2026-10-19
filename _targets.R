# ============================================================================
# targets pipeline for austraits-means
# ============================================================================
# Orchestration lives here; all logic stays in R/ functions (loaded via
# R/init.R). The two Quarto notebooks are now *reporting layers* rendered at the
# end of the pipeline (tar_quarto): they tar_load() results rather than compute.
#
# Run:   targets::tar_make()
# Graph: targets::tar_visnetwork()  /  targets::tar_outdated()
#
# TIERING ----------------------------------------------------------------------
# The heavy upstream artefacts (raw GBIF cleaning ~1.2GB, the 22k-species O(n^2)
# medoid build) are expensive and rarely rebuilt, and their producers live
# outside the two notebooks
# (producer notebooks climate_centres.qmd, Data_height_AVH.qmd at the root; builders in
# R/data/build_sources.R). They are tracked here as *file inputs* — targets
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

  tar_target(tmb_fit,
             fit_trait_target(unified, trait, backend = "tmb",
                              n_min = 5, out_dir = out_dir)),

  tar_target(tmb_row,
             if (is.null(tmb_fit)) NULL
             else summarise_trait_fits(stats::setNames(list(tmb_fit), trait))),

  tar_target(tmb_fig, {
    if (is.null(tmb_fit)) character(0) else {
      p <- plot_trait_centres(tmb_fit$centres, trait,
                              x = "precip", log_y = TRUE)
      save_figure(p, file.path(out_dir, "figures",
                               clean_species_filename(trait)),
                  width = 7, height = 5)
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
  # The branch targets (tmb_fit_<trait>, tmb_row_<trait>, tmb_fig_<trait>),
  # defined as `mapped` above.
  mapped,

  # Fixed-effects summary across traits (one row per trait).
  tar_combine(trait_summary, mapped[["tmb_row"]],
              command = dplyr::bind_rows(!!!.x)),
  tar_target(trait_summary_file, {
    p <- file.path(out_dir, "trait_model_summary.csv")
    readr::write_csv(trait_summary, p)
    p
  }, format = "file"),

  # List of fits — used to gate the combined product on all exports existing.
  tar_combine(tmb_fits, mapped[["tmb_fit"]], command = list(!!!.x)),

  # Data basis of the fits (build version, APC release, GBIF extract, medoid
  # build), written next to the exports.
  tar_target(data_basis_info,
             data_basis(austraits_dev_info, apc_version, gbif_cleaned_file,
                        medoids_file)),
  tar_target(data_basis_file, {
    p <- file.path(out_dir, "data_basis.json")
    jsonlite::write_json(data_basis_info, p, auto_unbox = TRUE, pretty = TRUE)
    p
  }, format = "file"),

  # Shifts against the June 2026 outputs, snapshotted before the refit (#14).
  tar_target(baseline_dir, "docs/baselines/2026-06", format = "file"),
  tar_target(baseline_comparison,
             compare_with_baseline(tmb_fits, unified, baseline_dir)),
  tar_target(baseline_comparison_files,
             write_baseline_comparison(baseline_comparison), format = "file"),

  # ---- Modelled means product (all traits, product schema; #17) ----------
  tar_target(modelled_means, modelled_means_table(tmb_fits)),
  tar_target(modelled_means_files,
             write_species_means(
               modelled_means, modelled_means_dictionary(),
               dir = file.path("outputs/species_means",
                               austraits_dev_info$austraits_version),
               stem = "modelled_means_tmb",
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
  # (R/scripts/create_zenodo_draft.R); nothing publishes.
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
  # Successor to the portal method: same recipe on the log scale
  # (R/means/observed_means_geometric.R, docs/observed-means-geometric-method.md).
  tar_target(geometric_traits,
             geometric_mean_traits(read_portal_core_traits(portal_core_traits_file))),
  tar_target(dev_records, prepare_records(austraits_dev_file, geometric_traits)),
  tar_target(observed_means_geometric,
             estimate_observed_means_geometric(dev_records)),
  tar_target(observed_means_geometric_files,
             write_species_means(
               label_estimates(observed_means_geometric,
                               "observed_mean", "geometric"),
               observed_means_geometric_dictionary(),
               dir = file.path("outputs/species_means",
                               austraits_dev_info$austraits_version),
               stem = "observed_means_geometric",
               metadata = observed_means_geometric_metadata(austraits_dev_info,
                                                            geometric_traits)),
             format = "file"),

  # ===== Community layers (#27, #28) ========================================
  # Three products on a 10 km equal-area grid, each for all species and for
  # woody and non-woody (docs/community-layers-design.md).
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
             combine_species_means(tmb_fits, observed_means_geometric)),
  tar_target(species_means_combined_files,
             write_species_means(
               species_means_combined, species_means_combined_dictionary(),
               dir = file.path("outputs/species_means",
                               austraits_dev_info$austraits_version),
               stem = "species_means_combined",
               metadata = species_means_combined_metadata(data_basis_info)),
             format = "file"),

  # Observed, local records: each species' records in a cell are averaged,
  # then the species. No model, no occurrence data.
  tar_target(community_local, community_layers_local(records_all, woody_lookup)),
  tar_target(community_local_files,
             write_community_layers(
               community_local, file.path(community_dir, "observed_local"),
               community_local_metadata(data_basis_info),
               statistics = c("community_mean", "community_variance", "n_species")),
             format = "file"),
  tar_target(community_local_maps,
             write_community_maps(community_local, community_units,
                                  file.path(community_dir, "observed_local"),
                                  "observed, local records", "observed_local"),
             format = "file"),

  # Observed, GBIF presence: the species recorded in a cell, each given its
  # geometric observed mean.
  tar_target(community_species_values,
             observed_species_values(observed_means_geometric, trait_names)),
  tar_target(gbif_presence, gbif_presence_cells(gbif_cleaned_file)),
  tar_target(gbif_species_map,
             gbif_to_species(gbif_names, name_crosswalk, dev_taxa,
                             union(community_species_values$taxon_name,
                                   woody_lookup$taxon_name))),
  tar_target(community_presence,
             community_layers_presence(gbif_presence, gbif_species_map,
                                       woody_lookup, community_species_values)),
  tar_target(community_presence_files,
             write_community_layers(
               community_presence, file.path(community_dir, "observed_presence"),
               community_presence_metadata(data_basis_info),
               statistics = c("community_mean", "community_variance", "n_species",
                              "coverage")),
             format = "file"),
  tar_target(community_presence_maps,
             write_community_maps(community_presence, community_units,
                                  file.path(community_dir, "observed_presence"),
                                  "observed, GBIF presence x species means",
                                  "observed_presence"),
             format = "file"),

  # Modelled: each trait model's net average for the cell, from the cell's
  # climate and the mean climate centre of the species present.
  tar_target(cell_climate_grid,
             cell_climate(load_climate_rasters(
               vars = basename(c(climate_bio1_tif, climate_bio12_tif)),
               climate_dir = dirname(climate_bio1_tif),
               australia_tif = climate_australia_tif))),
  tar_target(cell_centres,
             cell_species_centres(gbif_presence, medoids_file, gbif_species_map,
                                  woody_lookup)),
  tar_target(community_modelled,
             community_layers_modelled(tmb_fits, cell_climate_grid, cell_centres)),
  tar_target(community_modelled_files,
             write_community_layers(
               community_modelled, file.path(community_dir, "modelled"),
               community_modelled_metadata(data_basis_info),
               statistics = c("community_mean", "lower", "upper", "extrapolated")),
             format = "file"),
  tar_target(community_modelled_maps,
             write_modelled_maps(community_modelled, community_presence,
                                 community_units,
                                 file.path(community_dir, "modelled")),
             format = "file"),

  # The three products against each other, cell by cell.
  tar_target(community_pairs,
             community_layer_pairs(community_local, community_presence,
                                   community_modelled)),
  tar_target(community_comparison_files,
             write_layer_comparison(community_pairs, community_units,
                                    file.path(community_dir, "comparison")),
             format = "file"),

  # Native species only, for the presence and modelled products: a test of
  # how much introduced species move them (native anywhere in Australia per
  # the APC). The models are still fitted on all species.
  tar_target(gbif_native, gbif_native_status(gbif_names, apc_version)),
  tar_target(gbif_presence_native, native_presence(gbif_presence, gbif_native)),
  tar_target(community_presence_native,
             community_layers_presence(gbif_presence_native, gbif_species_map,
                                       woody_lookup, community_species_values)),
  tar_target(community_presence_native_files,
             write_community_layers(
               community_presence_native,
               file.path(community_dir, "observed_presence_native"),
               community_presence_metadata(data_basis_info, native_only = TRUE),
               statistics = c("community_mean", "community_variance", "n_species",
                              "coverage")),
             format = "file"),
  tar_target(community_presence_native_maps,
             write_community_maps(community_presence_native, community_units,
                                  file.path(community_dir, "observed_presence_native"),
                                  "observed, GBIF presence x species means, native species",
                                  "observed_presence_native"),
             format = "file"),
  tar_target(cell_centres_native,
             cell_species_centres(gbif_presence_native, medoids_file,
                                  gbif_species_map, woody_lookup)),
  tar_target(community_modelled_native,
             community_layers_modelled(tmb_fits, cell_climate_grid,
                                       cell_centres_native)),
  tar_target(community_modelled_native_files,
             write_community_layers(
               community_modelled_native, file.path(community_dir, "modelled_native"),
               community_modelled_metadata(data_basis_info, native_only = TRUE),
               statistics = c("community_mean", "lower", "upper", "extrapolated")),
             format = "file"),
  tar_target(community_modelled_native_maps,
             write_modelled_maps(community_modelled_native, community_presence_native,
                                 community_units,
                                 file.path(community_dir, "modelled_native"),
                                 species_label = "native species"),
             format = "file"),
  tar_target(community_comparison_native_files,
             write_layer_comparison(
               community_layer_pairs(community_local, community_presence_native,
                                     community_modelled_native),
               community_units, file.path(community_dir, "comparison_native"),
               species_label = "native species"),
             format = "file"),
  tar_target(community_native_shift_files,
             write_layer_comparison(
               native_layer_pairs(community_presence, community_presence_native,
                                  community_modelled, community_modelled_native),
               community_units, file.path(community_dir, "comparison_native_vs_all"),
               differences = c("presence_all_native", "modelled_all_native"),
               species_label = ""),
             format = "file"),

  # ===== Observed vs modelled means (#24) ===================================
  tar_target(observed_modelled,
             compare_observed_modelled(tmb_fits, observed_means_geometric, unified)),
  tar_target(portal_geometric,
             compare_portal_geometric(observed_means_portal,
                                      observed_means_geometric)),
  tar_target(means_comparison_files,
             write_means_comparison(observed_modelled, portal_geometric),
             format = "file"),

  # ===== Figure set (#25) ====================================================
  tar_target(figure_set,
             write_figure_set(observed_means_portal, observed_means_geometric,
                              modelled_means, tmb_fits, observed_modelled,
                              community_presence, means_comparison_files,
                              c(community_local_maps, community_presence_maps,
                                community_modelled_maps,
                                community_comparison_files)),
             format = "file"),

  # ---- reporting layers (rendered notebooks) ------------------------------
  tar_quarto(report_dataset, "dataset_construction.qmd"),
  tar_quarto(report_analysis, "analysis.qmd")

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
                                 backend = "brms", n_min = 5, out_dir = out_dir,
                                 fits_dir = "data-processed/brms_fits",
                                 threads = max(1, floor(
                                   parallelly::availableCores(
                                     methods = "system") / 4)),
                                 iter = 2000, warmup = 1000,
                                 refresh = 250, silent = 1)),
    tar_target(brms_diag, summarise_trait_diagnostics(brms_fits)),
    tar_target(brms_diag_file, {
      p <- file.path(out_dir, "brms_diagnostics.csv")
      readr::write_csv(brms_diag, p); p
    }, format = "file")
  )))()
