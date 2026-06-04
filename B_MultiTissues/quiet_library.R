# ─────────────────────────────────────────────────────────────────────────────
#  Library loading
#  Install if missing (CRAN, Bioconductor, GitHub), load silently, report version
# ─────────────────────────────────────────────────────────────────────────────

for (pkg in c("BiocManager", "remotes")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

cran_packages <- c(
  # Data formatting / handling
  "dplyr", "data.table", "matrixStats", "tidyr", "tibble",
  "purrr", "forcats", "stringr",
  # HDF5 / parallel
  "parallel", "rhdf5",
  # Graphics
  "ggplot2", "ggrastr", "ggrepel", "ggExtra", "ggVennDiagram", "patchwork",
  "scales", "viridis", "cowplot", "gridGraphics", "grid", "Cairo",
  "UpSetR",
  # Stats
  "boot", "emmeans",
  # Reporting
  "gt", "progress"
)

bioc_packages <- c(
  # Methylation / genomics
  "rtracklayer", "genomation",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "GenomicRanges", "IRanges",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "GSEABase",
  "methylKit", 
  # Enrichment
  "rGREAT", "simplifyEnrichment", "org.Hs.eg.db", "BioMartGOGeneSets",
  "AnnotationHub"
)

# ─────────────────────────────────────────────────────────────────────────────
#  Installer / loader
# ─────────────────────────────────────────────────────────────────────────────

install_and_load <- function(packages, source = "cran") {
  not_loaded <- character()
  
  for (pkg in packages) {
    pkg_name <- basename(pkg)  # handles "user/repo" GitHub format
    
    # Install if missing
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      message("Installing: ", pkg)
      tryCatch(
        switch(source,
               cran   = install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE),
               bioc   = BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE),
               github = remotes::install_github(pkg, quiet = TRUE)
        ),
        error = function(e) message("  Failed to install ", pkg, ": ", e$message)
      )
    }
    
    # Load silently
    loaded <- suppressPackageStartupMessages(
      suppressMessages(suppressWarnings(
        require(pkg_name, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      ))
    )
    
    if (!loaded) {
      not_loaded <- c(not_loaded, pkg_name)
    } else {
      ver <- as.character(utils::packageVersion(pkg_name))
      cat(sprintf("  %-55s %s\n", pkg_name, ver))
    }
  }
  
  if (length(not_loaded) > 0)
    message("  Could not load: ", paste(not_loaded, collapse = ", "))
}

# ─────────────────────────────────────────────────────────────────────────────
#  Load
# ─────────────────────────────────────────────────────────────────────────────

message("── CRAN packages ──────────────────────────────────────────────")
install_and_load(cran_packages, source = "cran")

message("── Bioconductor packages ───────────────────────────────────────")
install_and_load(bioc_packages, source = "bioc")

# ─────────────────────────────────────────────────────────────────────────────
#  Cleanup
# ─────────────────────────────────────────────────────────────────────────────

rm(install_and_load, cran_packages, bioc_packages)
libLoaded <- TRUE