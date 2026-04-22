library(here)

source(here("05_hvCpGalgorithm", "quiet_library.R"))
if (!exists("functionsLoaded")) {
  source(here("05_hvCpGalgorithm/exploreResults", "functions.R"))
}

library(rGREAT)

overlapLayers <- readRDS(here("gitignore/overlapLayers.RDS"))
topIntersect90 <- readRDS(here("gitignore/topIntersect90.RDS"))

background <- makeGRfromMyCpGPos(overlapLayers, "background")
foreground <- makeGRfromMyCpGPos(topIntersect90, "topIntersect90")

# ✅ ONLINE GREAT - NO COMPILATION ISSUES
system.time({
  job <- submitGreatJob(
    regions = foreground,
    species = "hg38",
    background = background,
    request_id = "topIntersect90_hg38"
  )
  
  # Wait for completion & get GO BP results
  res <- getEnrichmentTables(job, gene_sets = "GO Biological Process")
})

saveRDS(res, here("05_hvCpGalgorithm/exploreResults/annotations/topIntersect90_rGREAT.RDS"))
print("SUCCESS: rGREAT analysis complete!")
