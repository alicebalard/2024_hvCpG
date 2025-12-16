###################################################################
## Plot secondary analyses results of algorithm ran on atlas data ##
####################################################################
library(here)

if (!exists("libLoaded")) {
  source(here("05_hvCpGalgorithm", "quiet_library.R"))
}

if (!exists("functionsLoaded")) {
  source(here("05_hvCpGalgorithm/exploreResults", "functions.R"))
}

if (!exists("resCompArray")) {
  source(here("05_hvCpGalgorithm/exploreResults/S02_analyseResultsArray_local.R"))  
}

##################################
## Save all data in RDS objects ##
##################################

for (file in list.files(here("05_hvCpGalgorithm/resultsDir/Atlas/"))){
  if (!file.exists(here(paste0("gitignore/fullres_", file)))){
    ## Add previous MEs including Maria's results if not sourced yet
    if (!exists("KesslerSIV_GRanges_hg38")) {
      source(here("05_hvCpGalgorithm/exploreResults/prepPreviousSIV.R"))
    }
    
    system.time(Atlas_dt <- prepAtlasdt(file))
    
    print("Number of CpG tested:")
    print(nrow(Atlas_dt))
    
    print(paste0("Saving results for ", file, "in ", here(paste0("gitignore/fullres_", file))))
    saveRDS(Atlas_dt, file = here(paste0("gitignore/fullres_", file)))
    print("Saved")
  }
}

##################
## 1_byDevLayer ## --> not very useful
##################

if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_0_vs_1_byDevLayer.pdf")))){
  makeCompPlot(
    X = readRDS(here::here("gitignore/fullres_Atlas10X")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_1_byDevLayer")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_0_vs_1_byDevLayer",
    xlab = "Pr(hv) calculated on WGBS atlas datasets",
    ylab = "Pr(hv) calculated on WGBS atlas datasets cut in 3 groups (germ layers)")
}

#####################
## 2_rmMultSamples ## 
#####################

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

###########################
## 3_correspMariaTissues ## --> very similar than full atlas
###########################

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

################
## Sex effect ##
################

## 1/ male
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_4_vs_5_maleEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here("gitignore/fullres_Atlas10X_4_maleOnly")),
    Y = readRDS(here("gitignore/fullres_Atlas10X_5_allButMaleOnly")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_4_vs_5_maleEffect",
    xlab = "Pr(hv) calculated on WGBS atlas males only datasets",
    ylab = "Pr(hv) calculated on WGBS atlas all but males only datasets")
}

## 2/ female
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_6_vs_7_femaleEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here("gitignore/fullres_Atlas10X_6_femaleOnly")),
    Y = readRDS(here("gitignore/fullres_Atlas10X_7_allButfemaleOnly")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_6_vs_7_femaleEffect",
    xlab = "Pr(hv) calculated on WGBS atlas females only datasets",
    ylab = "Pr(hv) calculated on WGBS atlas all but females only datasets")
}

################
## 8_byTissue ##
################

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

#########################
## 9_immune cells only ##
#########################

if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_0_vs_9_immuneOnly.pdf")))){
  makeCompPlot(
    X = readRDS(here::here("gitignore/fullres_Atlas10X")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_0_vs_9_immuneOnly",
    xlab = "Pr(hv) calculated on WGBS atlas datasets",
    ylab = "Pr(hv) calculated on WGBS atlas datasets, immune cells only")
}

#############################
## 10_immune cells removed ##
#############################

## Hypothesis testing:

## Prediction 1:
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_10_vs_9_immuneEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here::here("gitignore/fullres_Atlas10X_10_noImmune")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_10_vs_9_immuneEffect",
    xlab = "Pr(hv) calculated on WGBS atlas datasets, immune cells removed",
    ylab = "Pr(hv) calculated on WGBS atlas datasets, only immune cells")
}

Z_inner_immvsnoimm <- makeZ_inner(
  X = readRDS(here::here("gitignore/fullres_Atlas10X_10_noImmune")),
  Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
  whichAlphaX = "alpha",
  whichAlphaY = "alpha")

stable <- Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X < 0.5 & Z_inner_immvsnoimm$alpha_Y < 0.5]
cellUniversal <- Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X > 0.5 & Z_inner_immvsnoimm$alpha_Y > 0.5]
immune <- Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X < 0.5 & Z_inner_immvsnoimm$alpha_Y > 0.5]
rest <-  Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X > 0.5 & Z_inner_immvsnoimm$alpha_Y < 0.5]

message(paste0("We found ", 
               length(stable), " (", round(length(stable)/length(Z_inner_immvsnoimm$name)*100), "%) stable CpGs, ",
               length(cellUniversal), " (", round(length(cellUniversal)/length(Z_inner_immvsnoimm$name)*100), "%) cell universal hvCpGs, ", 
               length(immune), " (", round(length(immune)/length(Z_inner_immvsnoimm$name)*100), "%) immune cell hvCpGs and ",
               length(rest), " (", round(length(rest)/length(Z_inner_immvsnoimm$name)*100), "%) hvCpGs undetected in immune cells, out of a total of ",
               length(Z_inner_immvsnoimm$name), " CpGs investigated"))

## We found 19672473 (86%) stable CpGs, 1047081 (5%) cell universal hvCpGs and 1464480 (6%) immune cell hvCpGs, out of a total of 22805115 CpGs investigated

## Test enrichment in the 3 categories

## Create GRange objects
makeGR_CpGset <- function(ids){
  chr = sub("_.+$", "", ids, perl = TRUE)
  pos = as.integer(sub("^.*_", "", ids, perl = TRUE))
  
  gr = GRanges(
    seqnames = Rle(chr),
    ranges   = IRanges(start = pos, width = 1),
    genome   = "hg38")
  
  gr = keepStandardChromosomes(gr, pruning.mode = "coarse")
  seqlevelsStyle(gr) <- "UCSC"   # Keep standard chroms
  return(gr)
}

gr_bgr <- makeGR_CpGset(Z_inner_immvsnoimm$name)

for (x in c("cellUniversal", "immune", "rest", "stable")){
  if (!file.exists(here(paste0("05_hvCpGalgorithm/exploreResults/annotations_rGREAT/", x, ".RDS")))){
    gr <- makeGR_CpGset(get(x))
    res <- great(gr, "C5", "hg38", background = gr_bgr)
    saveRDS(res, file = here(paste0("05_hvCpGalgorithm/exploreResults/annotations_rGREAT/", x, ".RDS")))
  }
}


## Info:
# https://www.youtube.com/watch?v=_ycOp3P4AG0

tb <- getEnrichmentTable(res_immune)

head(tb)

head(tb[["GO Biological Process"]])
## Use "binom" no "hyper"

plotVolcano(res_immune)

## To get info for a given region
getRegionGeneAssociations()

## Simplify GO terms
simplifyGO



## Prediction 2:
Z_inner_all <- makeZ_inner(
  X = resCompArray,
  Y = readRDS(here::here("gitignore/fullres_Atlas10X")),
  whichAlphaX = "alpha_array_all",
  whichAlphaY = "alpha")

Z_inner_noimm <- makeZ_inner(
  X = resCompArray,
  Y = readRDS(here::here("gitignore/fullres_Atlas10X_10_noImmune")),
  whichAlphaX = "alpha_array_all",
  whichAlphaY = "alpha")

Z_inner_onlyimm <- makeZ_inner(
  X = resCompArray,
  Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
  whichAlphaX = "alpha_array_all",
  whichAlphaY = "alpha")

c_noimm <- cor.test(Z_inner_noimm$alpha_X, Z_inner_noimm$alpha_Y)
c_all <- cor.test(Z_inner_all$alpha_X, Z_inner_all$alpha_Y)
c_onlyimm <- cor.test(Z_inner_onlyimm$alpha_X, Z_inner_onlyimm$alpha_Y)

slope_noimm <- lm(data = Z_inner_noimm, alpha_Y ~ alpha_X)$coefficients[["alpha_X"]]
slope_all <- lm(data = Z_inner_all,  alpha_Y ~ alpha_X)$coefficients[["alpha_X"]]
slope_onlyimm <- lm(data = Z_inner_onlyimm,  alpha_Y ~ alpha_X)$coefficients[["alpha_X"]]

p <- data.frame(x=seq(0,1,0.01)) %>%
  dplyr::mutate(noImmun = slope_noimm * x,
                all = slope_all * x,
                onlyimm = slope_onlyimm * x) %>%
  pivot_longer(cols = c(noImmun, all, onlyimm)) %>%
  ggplot(aes(x = x, y = value, group = name, col = name)) +
  geom_abline(slope = 1, linetype = 3) +
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  theme_minimal(base_size = 14) +
  labs(x = "Pr(hv) array datasets", y = "Pr(hv) WGBS datasets", col = "WGBS datasets") +
  ylim(c(0,1))

ggplot2::ggsave(
  filename = here::here("05_hvCpGalgorithm/figures/correlations/correlation_prediction2.pdf"),
  plot = p, width = 9, height = 7
)
