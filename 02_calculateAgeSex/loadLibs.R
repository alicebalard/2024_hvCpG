## load all libraries needed for the project
list.of.packages <- c(
  "tidyverse")

###################################################################
message("Install CRAN packages if missing, and load CRAN packages...")

## install from CRAN and require all libraries from CRAN and github
install_if_missing <- function(packages, dependencies = TRUE) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    message("Installing missing packages: ", paste(new_packages, collapse = ", "))
    install.packages(new_packages, dependencies = dependencies)
  } else {
    message("All CRAN packages are already installed.\n")
  }
}

install_if_missing(list.of.packages)

message("Loading CRAN packages...")
load_packages <- function(package_list) {
  not_loaded <- character()
  for (pkg in package_list) {
    if (suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
      # Package loaded successfully, do nothing
    } else {
      not_loaded <- c(not_loaded, pkg)
    }
  }  
  if (length(not_loaded) > 0) {
    message("The following packages could not be loaded: ", paste(not_loaded, collapse = ", "))
  } else {
    message("All packages loaded successfully.\n")
  }
}

load_packages(list.of.packages)

#####################################################
## install from biocmanager and require all libraries
## Biocmanager packages 
bioc_packages <- c("plyranges", 
                   "GenomicRanges") 

install_and_load_bioc_packages <- function(package_list) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  not_loaded <- character()
  
  for (pkg in package_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing", pkg, "..."))
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
    
    if (suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))) {
      # Package loaded successfully, do nothing
    } else {
      not_loaded <- c(not_loaded, pkg)
    }
  }
  
  if (length(not_loaded) > 0) {
    message("The following packages could not be loaded: ", paste(not_loaded, collapse = ", "))
  } else {
    message("All packages loaded successfully.\n")
  }
}

message("Loading bioconductor packages...")
install_and_load_bioc_packages(bioc_packages)


## MammalMethylClock installed in R 4.2.1
library(MammalMethylClock)
