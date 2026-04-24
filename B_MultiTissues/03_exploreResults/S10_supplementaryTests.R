
#########################
## correspMariaTissues ## --> very similar than full atlas
#########################

if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_0_vs_3_correspMariaTissues.pdf")))){
  makeCompPlot(
    X = readRDS(here("gitignore/fullres_Atlas10X")),
    Y = readRDS(here("gitignore/fullres_Atlas10X_3_correspMariaTissues")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_0_vs_3_correspMariaTissues",
    xlab = "Pr(hv) calculated on WGBS atlas datasets",
    ylab = "Pr(hv) calculated on WGBS atlas datasets with only cells found in array")
}

##############
## byTissue ##
##############

if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_0_vs_8_byTissue.pdf")))){
  ## Cut by tissue rather than by cell type. Is is closer to array data?
  makeCompPlot(
    X = readRDS(here::here("gitignore/fullres_Atlas10X")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_8_byTissue")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_0_vs_8_byTissue",
    xlab = "Pr(hv) calculated on WGBS atlas datasets (cut by cell types)",
    ylab = "Pr(hv) calculated on WGBS atlas datasets cut by tissues")
}

if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Array_vs_8_byTissue.pdf")))){
  makeCompPlot(
    X = resCompArray,
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_8_byTissue")),
    whichAlphaX = "alpha_array_all",
    whichAlphaY = "alpha",          
    title = "Array_vs_8_byTissue",
    xlab = "Pr(hv) calculated on array datasets",
    ylab = "Pr(hv) calculated on WGBS atlas datasets cut by tissues")
}


###################
## rmMultSamples ## 
###################

## Some individuals have multiple cells sampled. Does that affect our results? NOPE
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_0_vs_2_rmMultSamples.pdf")))){
  makeCompPlot(
    X = readRDS(here::here("gitignore/fullres_Atlas10X")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_2_rmMultSamples")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_0_vs_2_rmMultSamples",
    xlab = "Pr(hv) calculated on WGBS atlas datasets",
    ylab = "Pr(hv) calculated on WGBS atlas datasets keeping one sample/individual only")
}