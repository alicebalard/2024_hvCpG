quiet_library <- function(pkg) {
  # Check if installed
  installed <- requireNamespace(pkg, quietly = TRUE)
  
  # Install if missing
  if (!installed) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      suppressMessages(suppressWarnings(
        install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
      ))
    }
    
    cran_pkgs <- suppressMessages(available.packages(repos = "https://cloud.r-project.org"))
    if (pkg %in% rownames(cran_pkgs)) {
      suppressMessages(suppressWarnings(
        install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
      ))
    } else {
      suppressMessages(suppressWarnings(
        BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
      ))
    }
  }
  
  # Load silently
  suppressPackageStartupMessages(
    suppressMessages(
      suppressWarnings(
        library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      )
    )
  )
  
  # Get version after loading
  version <- as.character(utils::packageVersion(pkg))
  
  # Print only our message
  cat(sprintf("Load package %s v%s\n", pkg, version))
}

quiet_library_all <- function(pkgs) {
  invisible(lapply(pkgs, quiet_library))
}

quiet_library_all(
  c("dplyr", "data.table", "matrixStats", "reshape2","tidyr", "tibble", ## data formatting
    "parallel", "rhdf5",  "stringr", ## data handling
    "ggplot2", "progress", "ggrastr", "ggrepel", "scales", 
    "UpSetR", "gridGraphics", "grid", "cowplot","ggExtra", ## graphical
    "boot", "emmeans", ## stats
    "rtracklayer", "IlluminaHumanMethylation450kanno.ilmn12.hg19", 
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "GenomicRanges", "rGREAT", ## methylation
    "testthat" # grammar & packaging
  ))
## NB: not all libraries are necessary; to clean when packaging

rm(quiet_library_all, quiet_library)

## load only in local, error in HPC CS
#ggVennDiagram
#methylKit
#Cairo

libLoaded = TRUE
