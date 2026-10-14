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
austraits_dev_rds <- "data/austraits/austraits-7.0.0.9000-cd26bb2.rds"

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

# Community layers per species-means input (#20 is open, so both are built):
# combined (modelled where available, else observed) and modelled means
# only (at each species' climate centre).
community_variants <- tar_map(
  values = list(input = c("combined", "modelled_centre")),
  names = "input",

  tar_target(community_layers,
             compute_community_layers(gbif_presence, gbif_species_map,
                                      woody_lookup,
                                      select_species_means(species_means_combined,
                                                           input))),
  tar_target(community_layers_files,
             write_community_layers(community_layers,
                                    file.path(community_dir, input),
                                    community_layers_metadata(data_basis_info,
                                                              input = input)),
             format = "file"),
  tar_target(community_maps,
             write_community_maps(community_layers, species_means_combined,
                                  file.path(community_dir, input), input = input),
             format = "file")
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
  tar_target(gbif_cleaned_file, "data-processed/data_gbif_cleaned.parquet",
             format = "file"),
  tar_target(medoids_file, "data-processed/climate_centres.parquet",
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

  # ---- Round 2: simple means + combined best-estimate product -------------
  tar_target(simple_means, estimate_simple_means(records_all)),
  tar_target(combined, {
    tmb_fits  # dependency: ensures per-trait centre CSVs are exported first
    centre_estimates <- load_centre_estimates(
      sort(unique(simple_means$trait_name)), dir = out_dir)
    combine_estimates(simple_means, centre_estimates)
  }),
  tar_target(combined_file, {
    p <- file.path(out_dir, "species_trait_means_combined.csv.gz")
    readr::write_csv(combined, p); p
  }, format = "file"),

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
               label_estimates(observed_means_portal,
                               "observed_mean", "portal_arithmetic"),
               observed_means_portal_dictionary(),
               dir = "outputs/species_means/7.0.0",
               stem = "observed_means_portal_arithmetic",
               metadata = observed_means_portal_metadata()),
             format = "file"),

  # Reproduction check against the parquet the portal deployed (built 9 Feb
  # 2026; copied from the portal's git-ignored data build). Also asserted by
  # tests/testthat/test-observed-means-portal.R.
  tar_target(portal_deployed_file,
             "data/portal/austraits-7.0.0-full/austraits-species-averages.parquet",
             format = "file"),
  tar_target(portal_reproduction,
             compare_portal_reproduction(observed_means_portal,
                                         arrow::read_parquet(portal_deployed_file))),
  tar_target(portal_reproduction_files,
             write_portal_reproduction(portal_reproduction), format = "file"),

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

  # ===== Community layers, first draft (provisional choices; #22, #23) =====
  # Woody / non-woody per taxon, provisional pending #19.
  tar_target(woody_lookup, build_woody_lookup(observed_means_geometric)),
  tar_target(woody_lookup_files,
             write_woody_lookup(woody_lookup,
                                file.path("outputs/community",
                                          austraits_dev_info$austraits_version)),
             format = "file"),

  # Combined species means (modelled where available, else observed), GBIF
  # presence per 0.25 degree cell, and the community statistics (#23).
  tar_target(community_dir,
             file.path("outputs/community", austraits_dev_info$austraits_version)),
  tar_target(species_means_combined,
             combine_species_means(tmb_fits, observed_means_geometric)),
  tar_target(gbif_presence,
             gbif_presence_cells(gbif_cleaned_file,
                                 community_layer_settings()$resolution)),
  tar_target(gbif_species_map,
             gbif_to_species(gbif_names, name_crosswalk, dev_taxa,
                             union(species_means_combined$taxon_name,
                                   woody_lookup$taxon_name))),
  # The layer variants, one per species-means input (#20 decides).
  community_variants,

  # Model surface: each trait model's population-level prediction at the
  # cell's climate (species averaged out), all species together.
  tar_target(cell_climate_grid,
             cell_climate(load_climate_rasters(
               vars = basename(c(climate_bio1_tif, climate_bio12_tif)),
               climate_dir = dirname(climate_bio1_tif),
               australia_tif = climate_australia_tif))),
  tar_target(model_surface, predict_model_surface(tmb_fits, cell_climate_grid)),
  tar_target(model_surface_files,
             write_model_surface(model_surface,
                                 file.path(community_dir, "model_surface"),
                                 model_surface_metadata(data_basis_info)),
             format = "file"),
  tar_target(model_surface_maps,
             write_model_surface_maps(model_surface, community_layers_combined,
                                      species_means_combined,
                                      file.path(community_dir, "model_surface")),
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
