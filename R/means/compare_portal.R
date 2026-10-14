# Reproduction check for the portal-method observed means: compare this
# pipeline's output with the species averages the portal deployed for
# AusTraits 7.0.0.
#
# Numeric rows are compared at key level (taxon_name × trait_name × unit):
# numbers within a relative tolerance, strings exactly. Categorical rows are not
# unique per key (the portal keeps one row per distinct value's sources), so
# they are compared as multisets of whole rows. `row_id` and row order are
# excluded because they depend on package versions the portal did not pin.

#' Compare two portal-method observed means tables.
#'
#' @param ours output of `estimate_observed_means_portal()` (labels allowed).
#' @param deployed the deployed `austraits-species-averages.parquet`.
#' @param tol relative tolerance for numeric columns.
#' @param exclude columns left out of the comparison.
#' @return list: `pass` (logical), `summary` (one row per check), `columns`
#'   (mismatches per column) and `examples` (up to 20 mismatching rows per
#'   column, as key + both values).
compare_portal_reproduction <- function(ours, deployed, tol = 1e-8,
                                        exclude = c("row_id", "estimate_type",
                                                    "method")) {
  ours <- dplyr::select(ours, -dplyr::any_of(exclude))
  deployed <- dplyr::select(deployed, -dplyr::any_of(exclude))
  key <- c("taxon_name", "trait_name", "unit")
  cols <- intersect(names(ours), names(deployed))
  value_cols <- setdiff(cols, key)

  # ---- numeric rows: one row per key --------------------------------------
  num_o <- dplyr::filter(ours, type == "numerical")
  num_d <- dplyr::filter(deployed, type == "numerical")
  dup_o <- nrow(num_o) - nrow(dplyr::distinct(num_o, dplyr::across(dplyr::all_of(key))))
  dup_d <- nrow(num_d) - nrow(dplyr::distinct(num_d, dplyr::across(dplyr::all_of(key))))

  joined <- dplyr::full_join(
    dplyr::mutate(num_o[, cols], .in_ours = TRUE),
    dplyr::mutate(num_d[, cols], .in_deployed = TRUE),
    by = key, suffix = c(".ours", ".deployed"), na_matches = "na"
  )
  both <- dplyr::filter(joined, !is.na(.in_ours), !is.na(.in_deployed))

  col_checks <- purrr::map(value_cols, function(cc) {
    a <- both[[paste0(cc, ".ours")]]
    b <- both[[paste0(cc, ".deployed")]]
    bad <- !.same_values(a, b, tol)
    list(
      row = tibble::tibble(part = "numerical", column = cc, n_compared = length(a),
                           n_mismatch = sum(bad)),
      examples = if (any(bad)) {
        both[bad, key] |>
          dplyr::mutate(column = cc, ours = as.character(a[bad]),
                        deployed = as.character(b[bad])) |>
          utils::head(20)
      }
    )
  })

  # ---- categorical rows: multisets of whole rows -----------------------------
  cat_o <- dplyr::filter(ours, type == "categorical")[, cols]
  cat_d <- dplyr::filter(deployed, type == "categorical")[, cols]
  counted <- dplyr::full_join(
    dplyr::count(cat_o, dplyr::across(dplyr::everything()), name = "n_ours"),
    dplyr::count(cat_d, dplyr::across(dplyr::everything()), name = "n_deployed"),
    by = cols, na_matches = "na"
  ) |>
    dplyr::mutate(dplyr::across(c(n_ours, n_deployed), \(x) dplyr::coalesce(x, 0L)))
  cat_bad <- dplyr::filter(counted, n_ours != n_deployed)

  # For categorical mismatches, attribute them to columns: for each affected
  # taxon × trait, which columns' value sets differ.
  cat_keys <- dplyr::distinct(cat_bad, taxon_name, trait_name)
  cat_col_rows <- purrr::map(setdiff(value_cols, "type"), function(cc) {
    if (!nrow(cat_keys)) {
      return(tibble::tibble(part = "categorical", column = cc, n_compared = 0L,
                            n_mismatch = 0L))
    }
    set_of <- function(df) {
      df |>
        dplyr::semi_join(cat_keys, by = c("taxon_name", "trait_name")) |>
        dplyr::group_by(taxon_name, trait_name) |>
        dplyr::summarise(.groups = "drop",
                         v = paste(sort(unique(as.character(.data[[cc]])), na.last = TRUE),
                                   collapse = "|"))
    }
    s <- dplyr::full_join(set_of(cat_o), set_of(cat_d),
                          by = c("taxon_name", "trait_name"), suffix = c("_o", "_d"))
    tibble::tibble(part = "categorical", column = cc, n_compared = nrow(s),
                   n_mismatch = sum(!.same_values(s$v_o, s$v_d, tol)))
  })

  summary <- tibble::tribble(
    ~part, ~check, ~n,
    "all", "rows in ours", nrow(ours),
    "all", "rows deployed", nrow(deployed),
    "all", "columns only in ours", length(setdiff(names(ours), names(deployed))),
    "all", "columns only deployed", length(setdiff(names(deployed), names(ours))),
    "numerical", "rows in ours", nrow(num_o),
    "numerical", "rows deployed", nrow(num_d),
    "numerical", "duplicate keys in ours", dup_o,
    "numerical", "duplicate keys deployed", dup_d,
    "numerical", "keys only in ours", sum(is.na(joined$.in_deployed)),
    "numerical", "keys only deployed", sum(is.na(joined$.in_ours)),
    "numerical", "keys with any value mismatch",
      sum(Reduce(`|`, purrr::map(value_cols, \(cc) !.same_values(
        both[[paste0(cc, ".ours")]], both[[paste0(cc, ".deployed")]], tol)))),
    "categorical", "rows in ours", nrow(cat_o),
    "categorical", "rows deployed", nrow(cat_d),
    "categorical", "rows only in ours", sum(pmax(cat_bad$n_ours - cat_bad$n_deployed, 0)),
    "categorical", "rows only deployed", sum(pmax(cat_bad$n_deployed - cat_bad$n_ours, 0)),
    "categorical", "taxon x trait with any mismatch", nrow(cat_keys)
  )

  mismatch_checks <- c("columns only in ours", "columns only deployed",
                       "duplicate keys in ours", "keys only in ours",
                       "keys only deployed", "keys with any value mismatch",
                       "rows only in ours", "rows only deployed")
  list(
    pass = all(summary$n[summary$check %in% mismatch_checks] == 0) &&
      nrow(ours) == nrow(deployed),
    summary = summary,
    columns = dplyr::bind_rows(purrr::map(col_checks, "row"), cat_col_rows),
    examples = dplyr::bind_rows(purrr::map(col_checks, "examples"))
  )
}

# Elementwise equality: numbers within relative tolerance, everything else
# exactly; NA equals NA.
.same_values <- function(a, b, tol) {
  both_na <- is.na(a) & is.na(b)
  if (is.numeric(a) && is.numeric(b)) {
    close <- a == b | abs(a - b) <= tol * pmax(abs(a), abs(b))
  } else {
    close <- as.character(a) == as.character(b)
  }
  both_na | (!is.na(close) & close)
}

#' Write the reproduction check to `dir` as three CSVs.
#' @return the file paths.
write_portal_reproduction <- function(check, dir = "outputs/checks") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, paste0("portal_reproduction_",
                                 c("summary", "columns", "examples"), ".csv"))
  readr::write_csv(check$summary, paths[1])
  readr::write_csv(check$columns, paths[2])
  readr::write_csv(check$examples, paths[3])
  paths
}
