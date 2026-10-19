# Combine all individual species panel PDFs into a single appendix document
# with a title page. Panels are ordered alphabetically by species name.
#
# Usage:  Rscript R/combine_species_panels.R

library(qpdf)

out_dir_species <- "outputs/height/species"
jobs_csv        <- "outputs/height/species_panel_jobs.csv"
combined_pdf    <- "outputs/height/species_panels_appendix.pdf"

# Inclusion criterion used when generating panels (see the retired height
# report Data_height_AVH.qmd, in git history:
# jobs_to_run filters n_records > 50).
min_records <- 50

# Page size of the individual panels (see ggsave() in species_panel_functions.R)
page_w <- 14
page_h <- 15

# --- Determine panel order ---------------------------------------------------
# Order alphabetically by species name. Filenames are cleaned species names, so
# sorting the file paths gives alphabetical species order.

if (file.exists(jobs_csv)) {
  jobs <- read.csv(jobs_csv, stringsAsFactors = FALSE)
  jobs <- jobs[order(jobs$taxon_name), ]
  pdfs <- jobs$output_file[file.exists(jobs$output_file)]
} else {
  pdfs <- sort(list.files(out_dir_species,
                          pattern = "\\.pdf$", full.names = TRUE))
}

pdfs <- pdfs[basename(pdfs) != basename(combined_pdf)]

# Drop empty or damaged PDFs (e.g. a species whose generation failed) so one
# bad file doesn't halt the merge. Report any that are skipped.

is_valid <- function(f) {
  file.info(f)$size > 0 &&
    !inherits(try(qpdf::pdf_length(f), silent = TRUE), "try-error")
}
ok <- vapply(pdfs, is_valid, logical(1))
if (any(!ok)) {
  message("Skipping ", sum(!ok), " empty/damaged PDF(s):")
  message(paste0("  ", basename(pdfs[!ok]), collapse = "\n"))
}
pdfs <- pdfs[ok]

if (length(pdfs) == 0) {
  stop("No valid species panel PDFs found in ", out_dir_species)
}

n_species <- length(pdfs)
message(n_species, " species panels to combine")

# --- Build a title page ------------------------------------------------------
# Drawn to match the front matter of the retired report Data_height_AVH.qmd
# (in git history). This document is
# an appendix to that report.

# Simple word-wrap to a target number of characters per line.
wrap <- function(txt, width = 90) {
  strwrap(txt, width = width)
}

body <- paste(
  "This document is an appendix to the report “Analysis of height variation",
  "in Australia: Report for NVIS team”. It collects the individual",
  "species-level diagnostic panels referenced in that report.",
  "",
  "The work forms part of a project for the Australian Government Department",
  "of Climate Change, Energy, the Environment and Water (DCCEEW) to strengthen",
  "the plant trait data underpinning the National Vegetation Information",
  "System (NVIS)."
)

# One entry per subplot in each panel (see build_species_plot() in
# species_panel_functions.R; the layout is a, b, d, e).
panel_legend <- c(
  "a)  Species distribution — mapped occurrence records (AVH height records and GBIF range records).",
  "b)  Environmental space — the species' position in precipitation–temperature space, against all Australian taxa.",
  "d)  Difference by state — variation in measured height across Australian states.",
  "e)  Rainfall response — the fitted height–rainfall relationship for the species."
)

inclusion <- paste0(
  "Panels are included for the ", n_species, " species with more than ",
  min_records, " georeferenced height records in the Australasian Virtual ",
  "Herbarium (AVH) — below which there are too few data to characterise a ",
  "response. Panels are ordered alphabetically by species name."
)

slope_note <- paste0(
  "Fitted per-species parameters — including the within-species precipitation ",
  "slope — are tabulated in the results CSV accompanying the report ",
  "(species_height_climate_response.csv.gz)."
)

title_pdf <- tempfile(fileext = ".pdf")
pdf(title_pdf, width = page_w, height = page_h)
op <- par(mar = c(0, 0, 0, 0))
plot.new()

# Title block
text(0.5, 0.93, "Appendix",
     cex = 1.6, col = "grey40", adj = 0.5)
text(0.5, 0.88, "Species height–climate response panels",
     cex = 3.0, font = 2, adj = 0.5)
text(0.5, 0.83, "Daniel Falster & Elizabeth Wenk",
     cex = 1.5, adj = 0.5)

# Body text, left-aligned within a centred column
y <- 0.73
for (line in unlist(lapply(body, wrap, width = 95))) {
  text(0.12, y, line, cex = 1.15, adj = 0)
  y <- y - 0.021
}

# Panel legend
y <- y - 0.025
text(0.12, y, "What each panel shows", cex = 1.3, font = 2, adj = 0)
y <- y - 0.030
for (entry in panel_legend) {
  for (line in wrap(entry, width = 92)) {
    text(0.12, y, line, cex = 1.1, adj = 0)
    y <- y - 0.021
  }
  y <- y - 0.006
}

# Inclusion criteria block
y <- y - 0.020
for (line in wrap(inclusion, width = 95)) {
  text(0.12, y, line, cex = 1.1, font = 3, adj = 0)
  y <- y - 0.021
}

# Slope / results CSV note
y <- y - 0.020
for (line in wrap(slope_note, width = 95)) {
  text(0.12, y, line, cex = 1.1, font = 3, adj = 0)
  y <- y - 0.021
}

text(0.5, 0.04, format(Sys.Date(), "%d %B %Y"),
     cex = 1.2, col = "grey40", adj = 0.5)
par(op)
dev.off()

# --- Merge -------------------------------------------------------------------

qpdf::pdf_combine(c(title_pdf, pdfs), output = combined_pdf)
unlink(title_pdf)
message("Wrote ", combined_pdf, " (", n_species + 1, " pages, ",
        round(file.info(combined_pdf)$size / 1e6), " MB)")
