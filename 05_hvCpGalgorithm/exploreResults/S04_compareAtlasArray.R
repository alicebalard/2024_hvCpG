#############################################
## Overlap plot: Atlas (x) vs Array (y)    ##
#############################################
library(here)
source(here("05_hvCpGalgorithm", "quiet_library.R"))
source(here("05_hvCpGalgorithm/exploreResults", "functions.R"))

## Prepare Array data 
source(here("05_hvCpGalgorithm/exploreResults/S02_analyseResultsArray_local.R"))

makeCompPlot(
  X = resCompArray,
  Y = readRDS(here::here("gitignore/fullres_Atlas10X")),
  whichAlphaX = "alpha_array_all",
  whichAlphaY = "alpha",          
  title = "Array_vs_Atlas",
  xlab = "Pr(hv) calculated on array datasets",
  ylab = "Pr(hv) calculated on WGBS atlas datasets")

#################################################################
## Is the slope different if I take array reduced to 3 ind/ds? ##
#################################################################

makeCompPlot(
  X = resCompArray,
  Y = readRDS(here::here("gitignore/fullres_Atlas10X")),
  whichAlphaX = "alpha_array_3ind",
  whichAlphaY = "alpha",          
  title = "Array_vs_Atlas",
  xlab = "Pr(hv) calculated on array datasets with only 3 ind/ds",
  ylab = "Pr(hv) calculated on WGBS atlas datasets")

################################################################
## What are alpha array for different cutoffs of alpha atlas? ##
################################################################
Z_inner <- makeZ_inner(
  X = resCompArray,
  Y = readRDS(here::here("gitignore/fullres_Atlas10X")),
  whichAlphaX = "alpha_array_all",
  whichAlphaY = "alpha")

df_plot <- bind_rows(
  Z_inner %>% filter(alpha_Y > 0.8) %>% mutate(threshold = "> 0.8"),
  Z_inner %>% filter(alpha_Y > 0.7) %>% mutate(threshold = "> 0.7"),
  Z_inner %>% filter(alpha_Y > 0.6) %>% mutate(threshold = "> 0.6")
)

pdf(here("05_hvCpGalgorithm/figures/histCutoff-arrayAtlas.pdf"), width = 8, height = 5)
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
        here("05_hvCpGalgorithm/dataOut/CpGArray95moreAtlas50less.RDS"))

saveRDS(object = Z_inner[!is.na(Z_inner$alpha_X) & !is.na(Z_inner$alpha_Y) &
                           Z_inner$alpha_X > 0.95 & Z_inner$alpha_Y > 0.5,] %>%
          dplyr::select(name) %>% dplyr::rename(chrpos = name) %>% as.vector(),
        here("05_hvCpGalgorithm/dataOut/CpGArray95moreAtlas50more.RDS"))

saveRDS(object = Z_inner[!is.na(Z_inner$alpha_X) & !is.na(Z_inner$alpha_Y) &
                           Z_inner$alpha_X > 0.5 & Z_inner$alpha_Y > 0.5,] %>%
          dplyr::select(name) %>% dplyr::rename(chrpos = name) %>% as.vector(), 
        here("05_hvCpGalgorithm/dataOut/CpGArray50moreAtlas50more.RDS"))
