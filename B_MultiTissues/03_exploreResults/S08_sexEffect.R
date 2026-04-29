if (!exists("libLoaded")) {
  source(here("B_MultiTissues", "quiet_library.R"))}

################
## Sex effect ##
################

### only 6 groups, 3 people

## 1/ male
if (!file.exists(file.path(here::here("B_MultiTissues/dataOut/figures/correlations/correlation_Atlas_4_vs_6_maleEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here("gitignore/fullres_Atlas10X_4_maleOnly")),
    Y = readRDS(here("gitignore/fullres_Atlas10X_6_bothsexes6gp")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_4_vs_6_maleEffect",
    xlab = "Pr(hv) on WGBS atlas with 6 datasets of only males",
    ylab = "Pr(hv) on WGBS atlas with 6 datasets of mix males/females")
}

## 2/ female
if (!file.exists(file.path(here::here("B_MultiTissues/dataOut/figures/correlations/correlation_Atlas_5_vs_6_femaleEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here("gitignore/fullres_Atlas10X_5_femaleOnly6gp")),
    Y = readRDS(here("gitignore/fullres_Atlas10X_6_bothsexes6gp")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_5_vs_6_femaleEffect",
    xlab = "Pr(hv) on WGBS atlas with 6 datasets of only females",
    ylab = "Pr(hv) on WGBS atlas with 6 datasets of mix males/females")
}

### only 2 people, 22 groups

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
