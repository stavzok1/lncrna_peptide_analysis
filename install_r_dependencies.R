#!/usr/bin/env Rscript
# Install R packages required by tr_limma_de.R (limma + jsonlite).
# Safe to re-run; skips already-installed packages.

pkgs_cran <- c("jsonlite", "BiocManager")
for (p in pkgs_cran) {
  if (!suppressWarnings(requireNamespace(p, quietly = TRUE))) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
if (!suppressWarnings(requireNamespace("limma", quietly = TRUE))) {
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}
message("R deps OK: limma, jsonlite")
