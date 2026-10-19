# Save a figure as both PNG (for embedding/preview) and PDF (vector output).

#' Save a ggplot/grob to both <stem>.png and <stem>.pdf.
#'
#' @param plot a ggplot or grob.
#' @param path output path; any .png/.pdf extension is stripped and both are written.
#' @param width,height device size in inches.
#' @param dpi raster resolution for the PNG.
#' @param ... passed to ggplot2::ggsave().
#' @return (invisibly) a named character vector c(png=..., pdf=...).
save_figure <- function(plot, path, width = 8, height = 6, dpi = 300, ...) {
  stem <- sub("\\.(png|pdf)$", "", path, ignore.case = TRUE)
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)

  out <- c(png = paste0(stem, ".png"), pdf = paste0(stem, ".pdf"))
  ggplot2::ggsave(out[["png"]], plot, width = width, height = height, dpi = dpi, ...)
  ggplot2::ggsave(out[["pdf"]], plot, width = width, height = height, ...)

  invisible(out)
}
