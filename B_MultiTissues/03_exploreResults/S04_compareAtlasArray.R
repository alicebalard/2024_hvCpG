#############################################
## Overlap plot: Atlas (x) vs Array (y)    ##
#############################################

#####################################################################
## Prepare
library(here)
## Load libraries
source(here("B_MultiTissues", "quiet_library.R"))

## Load functions
if (!exists("functionsLoaded")) {
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))}

## Add previous MEs including Maria's results
## Load the set of previously tested MEs & vmeQTL
if (!exists("previousSIVprepared")) {
  source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))}
#####################################################################

## Prepare Array data 
source(here("B_MultiTissues/03_exploreResults/S02_analyseResultsArray_local.R"))

makeCompPlot(
  X = resCompArray_allvs3,
  Y = readRDS(here::here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_atlas_general.rds")),
  whichAlphaX = "alpha_array_all",
  whichAlphaY = "alpha",          
  title = "Array_vs_Atlas",
  xlab = "Pr(hv) calculated on array datasets",
  ylab = "Pr(hv) calculated on WGBS atlas datasets")

################################################################
## What are alpha array for different cutoffs of alpha atlas? ##
################################################################
Z_inner <- makeZ_inner(
  X = resCompArray_allvs3,
  Y = readRDS(here::here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_atlas_general.rds")),
  whichAlphaX = "alpha_array_all",
  whichAlphaY = "alpha")

df_plot <- bind_rows(
  Z_inner %>% filter(alpha_Y > 0.8) %>% mutate(threshold = "> 0.8"),
  Z_inner %>% filter(alpha_Y > 0.7) %>% mutate(threshold = "> 0.7"),
  Z_inner %>% filter(alpha_Y > 0.6) %>% mutate(threshold = "> 0.6")
)

pdf(here("B_MultiTissues/dataOut/figures/histCutoff-arrayAtlas.pdf"), width = 8, height = 5)
ggplot(df_plot, aes(x = alpha_X, fill=threshold)) +
  geom_histogram(bins = 30,color = "white") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("yellow", "orange", "red"))+
  xlab("p(hv) in array") +
  ylab("Count")
dev.off()

table(df_plot$threshold)

########################################################################
## What are the positions with a high alpha in array but low in WGBS? ##
########################################################################

saveRDS(object = Z_inner[!is.na(Z_inner$alpha_X) & !is.na(Z_inner$alpha_Y) &
                           Z_inner$alpha_X > 0.95 & Z_inner$alpha_Y < 0.5,] %>%
          dplyr::select(name) %>% dplyr::rename(chrpos = name) %>% as.vector(), 
        here("B_MultiTissues/dataOut/CpGArray95moreAtlas50less.RDS"))

saveRDS(object = Z_inner[!is.na(Z_inner$alpha_X) & !is.na(Z_inner$alpha_Y) &
                           Z_inner$alpha_X > 0.95 & Z_inner$alpha_Y > 0.5,] %>%
          dplyr::select(name) %>% dplyr::rename(chrpos = name) %>% as.vector(),
        here("B_MultiTissues/dataOut/CpGArray95moreAtlas50more.RDS"))

saveRDS(object = Z_inner[!is.na(Z_inner$alpha_X) & !is.na(Z_inner$alpha_Y) &
                           Z_inner$alpha_X > 0.5 & Z_inner$alpha_Y > 0.5,] %>%
          dplyr::select(name) %>% dplyr::rename(chrpos = name) %>% as.vector(), 
        here("B_MultiTissues/dataOut/CpGArray50moreAtlas50more.RDS"))
