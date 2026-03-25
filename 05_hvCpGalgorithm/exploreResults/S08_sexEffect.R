
################
## Sex effect ##
################

## 1/ male
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_4_vs_6_maleEffect.pdf")))){
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
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_5_vs_6_femaleEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here("gitignore/fullres_Atlas10X_5_femaleOnly6gp")),
    Y = readRDS(here("gitignore/fullres_Atlas10X_6_bothsexes6gp")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_5_vs_6_femaleEffect",
    xlab = "Pr(hv) on WGBS atlas with 6 datasets of only females",
    ylab = "Pr(hv) on WGBS atlas with 6 datasets of mix males/females")
}
