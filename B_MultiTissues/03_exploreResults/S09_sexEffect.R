#####################################################################
## Prepare
library(here)
## Load libraries
source(here("B_MultiTissues", "quiet_library.R"))

## Load functions
if (!exists("functionsLoaded")) {
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))}
#####################################################################

################
## Sex effect ##
################

### only 6 groups, 3 people

## 1/ male
if (!file.exists(file.path(here::here("B_MultiTissues/dataOut/figures/correlations/correlation_Atlas_4_vs_6_maleEffect_chrX.pdf")))){
  X = readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_04_maleOnly.rds"))
  Y = readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_06_bothsexes6gp.rds"))
  X_autosomes = X[X$chr %in% 1:22,]
  Y_autosomes = Y[Y$chr %in% 1:22,]
  X_chrX = X[X$chr %in% "X",]
  Y_chrX = Y[Y$chr %in% "X",]

  makeCompPlot(
    X = X_autosomes,
    Y = Y_autosomes,
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_4_vs_6_maleEffect_autosomes",
    xlab = "Pr(hv) on WGBS atlas with 6 datasets of only males",
    ylab = "Pr(hv) on WGBS atlas with 6 datasets of mix males/females")
  
  makeCompPlot(
    X = X_chrX,
    Y = Y_chrX,
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_4_vs_6_maleEffect_chrX",
    xlab = "Pr(hv) on WGBS atlas with 6 datasets of only males",
    ylab = "Pr(hv) on WGBS atlas with 6 datasets of mix males/females")
}

## 2/ female
if (!file.exists(file.path(here::here("B_MultiTissues/dataOut/figures/correlations/correlation_Atlas_5_vs_6_femaleEffect_chrX.pdf")))){
  X = readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_05_femaleOnly6gp.rds"))
  Y = readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_06_bothsexes6gp.rds"))
  X_autosomes = X[X$chr %in% 1:22,]
  Y_autosomes = Y[Y$chr %in% 1:22,]
  X_chrX = X[X$chr %in% "X",]
  Y_chrX = Y[Y$chr %in% "X",]
  
  makeCompPlot(
    X = X_autosomes,
    Y = Y_autosomes,
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_5_vs_6_femaleEffect_autosomes",
    xlab = "Pr(hv) on WGBS atlas with 6 datasets of only females",
    ylab = "Pr(hv) on WGBS atlas with 6 datasets of mix males/females")
  
  makeCompPlot(
    X = X_chrX,
    Y = Y_chrX,
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_5_vs_6_femaleEffect_chrX",
    xlab = "Pr(hv) on WGBS atlas with 6 datasets of only females",
    ylab = "Pr(hv) on WGBS atlas with 6 datasets of mix males/females")
}


### only 2 people, 22 groups (still running on cluster...)

## 1/ male
if (!file.exists(file.path(here::here("B_MultiTissues/dataOut/figures/correlations/correlation_Atlas_maleEffect_pairs_autosomes.pdf")))){
  X = readRDS(here("gitignore/fullres_10X_15_pairs_MM"))
  X_autosomes = X[X$chr %in% 1:22,]
  Y = readRDS(here("gitignore/fullres_10X_17_pairs_MF"))
  Y_autosomes = Y[Y$chr %in% 1:22,]
  
  makeCompPlot(
    X = X_autosomes,
    Y = Y_autosomes,
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_maleEffect_pairs (autosomes)",
    xlab = "Pr(hv) on WGBS atlas with 22 datasets of 2 males",
    ylab = "Pr(hv) on WGBS atlas with 22 datasets of 1 male / 1 female")
}

if (!file.exists(file.path(here::here("B_MultiTissues/dataOut/figures/correlations/correlation_Atlas_maleEffect_pairs_chrX.pdf")))){
  X = readRDS(here("gitignore/fullres_10X_15_pairs_MM"))
  X_chrX = X[X$chr %in% "X",]
  Y = readRDS(here("gitignore/fullres_10X_17_pairs_MF"))
  Y_chrX = Y[Y$chr %in% "X",]
  
  makeCompPlot(
    X = X_chrX,
    Y = Y_chrX,
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_maleEffect_pairs (chrX)",
    xlab = "Pr(hv) on WGBS atlas with 22 datasets of 2 males",
    ylab = "Pr(hv) on WGBS atlas with 22 datasets of 1 male / 1 female")
}

## 2/ female
if (!file.exists(file.path(here::here("B_MultiTissues/dataOut/figures/correlations/correlation_Atlas_femaleEffect_pairs_autosomes.pdf")))){
  X = readRDS(here("gitignore/fullres_10X_16_pairs_FF"))
  X_autosomes = X[X$chr %in% 1:22,]
  Y = readRDS(here("gitignore/fullres_10X_17_pairs_MF"))
  Y_autosomes = Y[Y$chr %in% 1:22,]
  
  makeCompPlot(
    X = X_autosomes,
    Y = Y_autosomes,
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_femaleEffect_pairs (autosomes)",
    xlab = "Pr(hv) on WGBS atlas with 22 datasets of 2 females",
    ylab = "Pr(hv) on WGBS atlas with 22 datasets of 1 male / 1 female")
}

if (!file.exists(file.path(here::here("B_MultiTissues/dataOut/figures/correlations/correlation_Atlas_femaleEffect_pairs_chrX.pdf")))){
  X = readRDS(here("gitignore/fullres_10X_16_pairs_FF"))
  X_chrX = X[X$chr %in% "X",]
  Y = readRDS(here("gitignore/fullres_10X_17_pairs_MF"))
  Y_chrX = Y[Y$chr %in% "X",]
  
  makeCompPlot(
    X = X_chrX,
    Y = Y_chrX,
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_femaleEffect_pairs (chrX)",
    xlab = "Pr(hv) on WGBS atlas with 22 datasets of 2 females",
    ylab = "Pr(hv) on WGBS atlas with 22 datasets of 1 male / 1 female")
}
