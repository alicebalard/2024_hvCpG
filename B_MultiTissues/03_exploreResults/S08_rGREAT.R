library(here)
if (!exists("libLoaded")) {
  source(here("B_MultiTissues", "quiet_library.R"))}

library(rGREAT)
library(stringr)

overlapLayers <- readRDS(here("gitignore/overlapLayers.RDS"))
topIntersect90 <- readRDS(here("gitignore/topIntersect90.RDS"))

background <- makeGRfromMyCpGPos(overlapLayers, "background")
foreground <- makeGRfromMyCpGPos(topIntersect90, "topIntersect90")

system.time(res <- great(gr = foreground, gene_sets = "GO:BP", biomart_dataset = "hg38", background = background, cores = 10))
  
saveRDS(res, here("B_MultiTissues/03_exploreResults/annotations/topIntersect90_rGREAT.RDS"))
print("SUCCESS: rGREAT analysis complete!")
