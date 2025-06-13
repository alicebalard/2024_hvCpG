## Arguments:
which="MariaArrays"

cat(paste0("Prepare ", which, " dataset and source functions for optimisation:\n"))
source("hvCpG_algorithm_detection.R")
cat("Functions sourced.\n")

cistrans_GoDMC_hvCpG_matched_control <- read.table("/home/alice/2024_hvCpG/03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt", header = T)

## Check
table(cistrans_GoDMC_hvCpG_matched_control$hvCpG_name %in% MariasCpGs$CpG)

cpgnames[cpgnames %in% cistrans_GoDMC_hvCpG_matched_control$hvCpG_name] %>% length
cpgnames[cpgnames %in% cistrans_GoDMC_hvCpG_matched_control$controlCpG_name]  %>% length

## Let's take 1500 hvCpG and matched mQTL controls
set.seed(123)

sub1500 <- cistrans_GoDMC_hvCpG_matched_control[sample(1:nrow(cistrans_GoDMC_hvCpG_matched_control), 1500), ]

test3000CpGsvec <- c(sub1500$hvCpG_name, sub1500$controlCpG_name)

## Nelder-Mead or L-BFGS-B?
result <- getAllOptimAlpha_parallel(cpgvec = test3000CpGsvec, optimMeth="Nelder-Mead", NCORES=8, p0=0.95)
save(result, test3000CpGsvec, MariasCpGs, file = "/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/results_NM_test3000_0.95.RDA")

result <- getAllOptimAlpha_parallel(cpgvec = test3000CpGsvec, optimMeth="L-BFGS-B", NCORES=8, p0=0.95)
save(result, test3000CpGsvec, MariasCpGs, file = "/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/results_LB_test3000_0.95.RDA")
## We choose Nelder-Mead, as L-BFGS-B produces weird zero values (see Rmd plots).

## Vary p0 true negative:
result <- getAllOptimAlpha_parallel(cpgvec = test3000CpGsvec, optimMeth="Nelder-Mead", NCORES=8, p0=0.80)
save(result, test3000CpGsvec, MariasCpGs, file = "/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/results_NM_test3000_0.80.RDA")

result <- getAllOptimAlpha_parallel(cpgvec = test3000CpGsvec, optimMeth="Nelder-Mead", NCORES=8, p0=0.99)
save(result, test3000CpGsvec, MariasCpGs, file = "/home/alice/2024_hvCpG/05_hvCpGalgorithm/resultsDir/results_NM_test3000_0.99.RDA")