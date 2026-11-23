# roll_up_to_species() (R/data/build_dataset.R), mirrored from
# austraits.build's observed_means_records() (austraits.build#899).

test_that("infraspecific taxa roll up to the name before their rank marker", {
  # Names and binomials as the development build's taxa table has them: for a
  # hybrid, `binomial` is the first two words, which is not the species
  taxa <- tibble::tibble(
    taxon_name = c(
      "Acacia a", "Acacia a subsp. b", "Acacia c var. d", "Acacia",
      "Mentha x piperita var. citrata",
      "Acacia trachycarpa x Acacia tumida var. pilbarensis",
      "Myoporum acuminatum - Myoporum boninense subsp. australe"
    ),
    taxon_rank = c("species", "subspecies", "variety", "genus",
                   "variety", "variety", "subspecies"),
    binomial = c("Acacia a", "Acacia a", "Acacia c", NA,
                 "Mentha x", "Acacia trachycarpa", "Myoporum acuminatum")
  )
  records <- tibble::tibble(taxon_name = taxa$taxon_name, value = seq_along(taxa$taxon_name))

  rolled <- roll_up_to_species(records, taxa)

  expect_equal(rolled$taxon_name, c(
    "Acacia a", "Acacia a", "Acacia c", "Acacia", "Mentha x piperita",
    # hybrid formulae and intergrades belong to neither parent
    "Acacia trachycarpa x Acacia tumida var. pilbarensis",
    "Myoporum acuminatum - Myoporum boninense subsp. australe"
  ))
  expect_equal(rolled$taxon_name_aligned, taxa$taxon_name)
})
