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

quiet_library_all(c("dplyr", "data.table", "matrixStats", "ggplot2", "reshape2", "ggrepel",
                    "parallel", "rhdf5", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "tidyr", "dplyr",
                    "rhdf5", "ggplot2", "data.table", "GenomicRanges", "rtracklayer", "stringr"
))
## NB: not all libraries are necessary; to clean when packaging