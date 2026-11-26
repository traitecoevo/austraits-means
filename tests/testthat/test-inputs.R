test_that("check_file_md5 passes a matching file and stops on any other", {
  f <- tempfile(fileext = ".txt")
  writeLines("austraits", f)
  md5 <- unname(tools::md5sum(f))

  expect_identical(check_file_md5(f, md5), f)
  # a development build has no published checksum
  expect_identical(check_file_md5(f, NULL), f)

  writeLines("austraits, edited", f)
  expect_error(check_file_md5(f, md5), "not the published")
  expect_error(check_file_md5(tempfile(), md5), "is missing")
})

test_that("austraits_build_info carries the version, commit and DOI", {
  f <- tempfile(fileext = ".rds")
  saveRDS(list(build_info = list(version = "8.0.0", git_SHA = "abc123")), f)

  expect_identical(
    austraits_build_info(f, "10.5281/zenodo.1"),
    list(austraits_version = "8.0.0", austraits_build_commit = "abc123",
         austraits_doi = "10.5281/zenodo.1"))
  # no DOI for a development build
  expect_identical(austraits_build_info(f)$austraits_doi, NA_character_)
})
