# Zenodo deposit "AusTraits species means": one version per AusTraits release
# (or development build), each linked to the release it derives from (#10;
# publication is #30).
#
# `prepare_zenodo_package()` assembles one version's files, README and
# deposit metadata into a directory. `create_zenodo_draft()` uploads that
# directory as a new *draft* deposition. Nothing here publishes: publication
# needs the maintainers' go-ahead and is done by hand on Zenodo.

#' Creators of the deposit: as on the AusTraits 7.0.0 Zenodo record.
zenodo_creators <- function() {
  list(
    list(name = "Falster, Daniel", affiliation = "University of New South Wales",
         orcid = "0000-0002-9814-092X"),
    list(name = "Gallagher, Rachael", affiliation = "Western Sydney University",
         orcid = "0000-0002-4680-8115"),
    list(name = "Wenk, Elizabeth", affiliation = "University of New South Wales",
         orcid = "0000-0001-5640-5910"),
    list(name = "Sauquet, Hervé",
         affiliation = "National Herbarium of New South Wales, Royal Botanic Gardens and Domain Trust Sydney",
         orcid = "0000-0001-8305-3236")
  )
}

#' README for one version of the species means deposit (markdown).
#'
#' @param version product version, e.g. "7.0.0".
#' @param austraits_doi DOI of the AusTraits release it derives from.
#' @param n_rows,n_taxa,n_traits size of the table.
#' @param commit this repository's commit that built it.
species_means_readme <- function(version, austraits_doi, n_rows, n_taxa,
                                 n_traits, commit) {
  fmt <- function(x) format(x, big.mark = ",")
  paste0(
"# AusTraits species means ", version, "

Species-level summaries of the plant trait records in AusTraits ", version, " (doi:", austraits_doi, "), so that the species means shown by the AusTraits data portal can be archived and cited.

This version holds the **observed means by the portal method** (`estimate_type = observed_mean`, `method = portal_arithmetic`): the species means the AusTraits data portal computes for release ", version, ", reproduced exactly. For numeric traits the value is an **arithmetic** mean of per-location means and per-observation flora midpoints; for categorical traits, counts of each value. The table has ", fmt(n_rows), " rows covering ", fmt(n_taxa), " taxa and ", n_traits, " traits.

Read the limitations in `observed-means-portal-method.md` before use. In particular, the column `value_geom_mean` is not a geometric mean (it equals `value_mean`, or is NaN where that is negative), records are not filtered by basis of record or life stage, and infraspecific names are not rolled up to species.

## Files

- `observed_means_portal_arithmetic.parquet`: the table (Apache Parquet).
- `observed_means_portal_arithmetic.csv.gz`: the same table as gzipped CSV; missing values are empty.
- `observed_means_portal_arithmetic_dictionary.csv`: definition and type of every column.
- `observed_means_portal_arithmetic_metadata.json`: provenance of the table.
- `observed-means-portal-method.md`: the workflow document: inputs, method, limitations, column dictionary and the reproduction test.
- `README.md`: this file.

## Versions

Each version carries the same version tag as the AusTraits release it derives from, so species means ", version, " derives from AusTraits ", version, ". A re-issue for the same release appends a fourth component (for example ", version, ".1). Later versions are planned to add the geometric observed means and the modelled means, the model-estimated mean at the centre of each species' climate distribution.

## How it was made

Built by the `targets` pipeline in https://github.com/traitecoevo/austraits-means (commit ", commit, ") from the AusTraits ", version, " release, using the method of the AusTraits data portal at commit 13479b9 (https://github.com/traitecoevo/austraits.portal). The output was tested against the species means file the portal deployed and matches it in every row and column.

## Licence and citation

CC BY 4.0. Please cite this deposit together with the AusTraits release it derives from (doi:", austraits_doi, ") and the AusTraits data paper (Falster et al. 2021, Scientific Data 8: 254, doi:10.1038/s41597-021-01006-6).
")
}

#' Deposit metadata for one version (Zenodo deposition API, legacy format).
#'
#' @param readme_md the README, whose text becomes the description.
zenodo_species_means_metadata <- function(version, austraits_doi, readme_md) {
  description <- readme_md |>
    sub(pattern = "^# [^\n]*\n", replacement = "") |>
    commonmark::markdown_html()
  list(
    upload_type = "dataset",
    title = "AusTraits species means",
    version = version,
    creators = zenodo_creators(),
    description = description,
    access_right = "open",
    license = "cc-by-4.0",
    keywords = list("plant traits", "Australia", "species means", "AusTraits",
                    "trait database", "functional ecology"),
    related_identifiers = list(
      list(identifier = austraits_doi, relation = "isSupplementTo",
           resource_type = "dataset"),
      list(identifier = "https://github.com/traitecoevo/austraits-means",
           relation = "isCompiledBy", resource_type = "software")
    ),
    communities = list(list(identifier = "austraits"))
  )
}

#' Assemble one version of the deposit into a directory.
#'
#' @param product_files the product's files (`write_species_means()` output).
#' @param doc_file the workflow document.
#' @param dir output directory (emptied of earlier package files first).
#' @return the paths written, including `zenodo_metadata.json` (not uploaded
#'   as a file; it is the deposit's metadata).
prepare_zenodo_package <- function(product_files, doc_file, dir, version,
                                   austraits_doi) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(dir, full.names = TRUE))
  data <- arrow::read_parquet(product_files[grepl("\\.parquet$", product_files)])
  meta <- jsonlite::read_json(product_files[grepl("_metadata\\.json$", product_files)])

  readme <- species_means_readme(
    version, austraits_doi, n_rows = nrow(data),
    n_taxa = dplyr::n_distinct(data$taxon_name),
    n_traits = dplyr::n_distinct(data$trait_name), commit = meta$commit)

  file.copy(c(product_files, doc_file), dir, overwrite = TRUE)
  writeLines(readme, file.path(dir, "README.md"))
  jsonlite::write_json(
    zenodo_species_means_metadata(version, austraits_doi, readme),
    file.path(dir, "zenodo_metadata.json"), auto_unbox = TRUE, pretty = TRUE)
  list.files(dir, full.names = TRUE)
}

#' Upload a prepared package as a draft deposition. Never publishes.
#'
#' Creates a new draft, or with `id` resumes an existing one. The metadata is
#' set first and files go up smallest first, so an interrupted upload leaves a
#' draft that lacks only the largest files; files already in the draft with a
#' matching MD5 are skipped, so the upload can simply be rerun. Requests are
#' retried on transient failures.
#'
#' @param dir a directory from `prepare_zenodo_package()`.
#' @param token a Zenodo personal access token (scope `deposit:write`).
#' @param sandbox TRUE for sandbox.zenodo.org (for trying the upload), FALSE
#'   for zenodo.org.
#' @param id an existing draft's id to resume, or NULL for a new draft.
#' @return list(id, url): the draft's id and the page to review it.
create_zenodo_draft <- function(dir, token = Sys.getenv("ZENODO_TOKEN"),
                                sandbox = TRUE, id = NULL) {
  if (!nzchar(token)) stop("Set ZENODO_TOKEN (a personal access token).")
  base <- if (sandbox) "https://sandbox.zenodo.org/api" else "https://zenodo.org/api"
  authed <- function(url) {
    httr2::request(url) |>
      httr2::req_auth_bearer_token(token) |>
      httr2::req_retry(max_tries = 5, retry_on_failure = TRUE,
                       is_transient = \(resp) httr2::resp_status(resp) %in%
                         c(429, 502, 503, 504))
  }
  api <- function(path) authed(base) |> httr2::req_url_path_append(path)

  metadata <- jsonlite::read_json(file.path(dir, "zenodo_metadata.json"))
  files <- setdiff(list.files(dir, full.names = TRUE),
                   file.path(dir, "zenodo_metadata.json"))
  files <- files[order(file.size(files))]

  draft <- if (is.null(id)) {
    api("deposit/depositions") |>
      httr2::req_body_raw("{}", type = "application/json") |>
      httr2::req_perform() |>
      httr2::resp_body_json()
  } else {
    api(paste0("deposit/depositions/", id)) |>
      httr2::req_perform() |>
      httr2::resp_body_json()
  }
  if (!identical(draft$submitted, FALSE)) stop("Deposition ", draft$id, " is not a draft.")

  api(paste0("deposit/depositions/", draft$id)) |>
    httr2::req_method("PUT") |>
    httr2::req_body_json(list(metadata = metadata)) |>
    httr2::req_perform()

  present <- authed(draft$links$bucket) |>
    httr2::req_perform() |>
    httr2::resp_body_json()
  present_md5 <- vapply(present$contents, \(x) sub("^md5:", "", x$checksum %||% ""),
                        character(1))
  names(present_md5) <- vapply(present$contents, \(x) x$key, character(1))

  for (f in files) {
    if (identical(unname(present_md5[basename(f)]), unname(tools::md5sum(f)))) next
    message(format(Sys.time(), "%H:%M:%S"), " uploading ", basename(f), " (",
            round(file.size(f) / 1e6, 1), " MB)")
    # No in-request retry here: a retried request cannot re-read the file
    # stream. Rerun the whole upload instead; finished files are skipped.
    httr2::request(draft$links$bucket) |>
      httr2::req_auth_bearer_token(token) |>
      httr2::req_url_path_append(basename(f)) |>
      httr2::req_method("PUT") |>
      httr2::req_body_file(f, type = "application/octet-stream") |>
      httr2::req_perform()
  }

  list(id = draft$id, url = draft$links$html)
}
