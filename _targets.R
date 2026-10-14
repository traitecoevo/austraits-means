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
# medoid build, the AVH<->GBIF crosswalk, the AusTraits flatten) are expensive
# and rarely rebuilt, and their producers live outside the two notebooks
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

# Static trait list for per-trait branching (must be known at pipeline-build
# time). Matches the sources in build_dataset(): height (AVH) + AusTraits.
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
  tar_target(gbif_cleaned_file, "data-processed/data_gbif_cleaned.parquet",
             format = "file"),
  tar_target(medoids_file, "data-processed/climate_centres.parquet",
             format = "file"),
  tar_target(crosswalk_file,
             "data-processed/gbif_avh_name_crosswalk.parquet",
             format = "file"),
  tar_target(height_climate_file,
             "data-processed/data_height_climate.parquet", format = "file"),
  tar_target(austraits_climate_file,
             "data-processed/data_traits_climate.parquet", format = "file"),

  # ===== Tier A: computed modelling pipeline ================================
  # Per-species medoid lookup (crosswalk-bridged), rebuilt if either input
  # parquet changes.
  tar_target(medoid_lookup,
             build_medoid_lookup(medoids = medoids_file,
                                 crosswalk = crosswalk_file)),

  # The linchpin: unified record-level table (AVH height + AusTraits), medoid-
  # centred. Was data-processed/data_traits_climate_unified.parquet.
  tar_target(unified,
             build_dataset(height_climate = height_climate_file,
                           austraits_climate = austraits_climate_file,
                           medoid_lookup = medoid_lookup)),

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

  # ---- Round 2: simple means + combined best-estimate product -------------
  tar_target(records_all,
             dplyr::bind_rows(load_height_records(height_climate_file),
                              load_austraits_records(austraits_climate_file))),
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
