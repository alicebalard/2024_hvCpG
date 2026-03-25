library(here)

if (!exists("functionsLoaded")) {
  source(here("05_hvCpGalgorithm/exploreResults", "functions.R"))
}

library(rGREAT)
library(stringr)

## Read vectors saved in R06
overlapLayers <- readRDS(here("gitignore/overlapLayers.RDS"))
topIntersect90 <- readRDS(here("gitignore/topIntersect90.RDS"))

background <- makeGRfromMyCpGPos(overlapLayers, "background")
foreground <- makeGRfromMyCpGPos(topIntersect90, "topIntersect90")

system.time(res <- great(gr = foreground, gene_sets = "GO:BP", biomart_dataset = "hg38", background = background, cores = 10))

saveRDS(res, file = here(paste0("05_hvCpGalgorithm/exploreResults/annotations/topIntersect90_rGREAT.RDS")))
