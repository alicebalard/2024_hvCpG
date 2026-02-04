###############################################################
## Explore the top 90% alpha overlapping between germ layers ##
###############################################################
library(here)

if (!exists("libLoaded")) {
  source(here("05_hvCpGalgorithm", "quiet_library.R"))
}
if (!exists("functionsLoaded")) {
  source(here("05_hvCpGalgorithm/exploreResults", "functions.R"))
}

## Read vectors saved in R06
overlapLayers <- readRDS(here("gitignore/overlapLayers.RDS"))
topIntersect90 <- readRDS(here("gitignore/topIntersect90.RDS"))

#################################
## Test enrichment of features ##
#################################

# Create GRanges
overlapLayers_GR <- makeGRfromMyCpGPos(overlapLayers, "overlapLayers")
topIntersect90_GR <- makeGRfromMyCpGPos(topIntersect90, "topIntersect90")

# Import bed file
bed_features <- genomation::readTranscriptFeatures(here("gitignore/hg38_GENCODE_V47.bed"))

# Annotate CpGs with features
topAnno <- genomation::annotateWithGeneParts(
  target  = topIntersect90_GR,
  feature = bed_features
)
topIntersect90_GR$featureType <- ifelse(topAnno@members[, "prom"] == 1, "promoter",
                                        ifelse(topAnno@members[, "exon"] == 1, "exon",
                                               ifelse(topAnno@members[, "intron"] == 1, "intron",
                                                      "intergenic")))
allAnno <- genomation::annotateWithGeneParts(
  target  = overlapLayers_GR,
  feature = bed_features
)
overlapLayers_GR$featureType <- ifelse(allAnno@members[, "prom"] == 1, "promoter",
                                       ifelse(allAnno@members[, "exon"] == 1, "exon",
                                              ifelse(allAnno@members[, "intron"] == 1, "intron",
                                                     "intergenic")))

## 1. Build counts for subset vs background 

# Feature levels in fixed order
feat_levels <- c("promoter", "exon", "intron", "intergenic")

# Count CpGs per feature
bg_counts   <- table(factor(overlapLayers_GR$featureType,   levels = feat_levels))
sub_counts  <- table(factor(topIntersect90_GR$featureType, levels = feat_levels))

# Combine into a 2x4 contingency table
cont_tab <- rbind(
  subset    = as.numeric(sub_counts),
  background = as.numeric(bg_counts)
)
colnames(cont_tab) <- feat_levels
cont_tab

enrich_list <- lapply(feat_levels, function(f) {
  # 2x2 table for feature f vs not‑f
  a <- sub_counts[f]                      # subset in feature f
  b <- sum(sub_counts)  - a              # subset not in f
  c <- bg_counts[f] - a                  # background in f but not in subset
  d <- sum(bg_counts) - bg_counts[f] - b # background not in f and not in subset
  
  mat <- matrix(c(a, b, c, d), nrow = 2,
                dimnames = list(
                  set      = c("subset", "background"),
                  inFeat   = c("yes", "no")
                ))
  
  ft <- fisher.test(mat, alternative = "greater")
  
  data.frame(
    feature      = f,
    subset_n     = as.numeric(a),
    bg_n         = as.numeric(bg_counts[f]),
    subset_prop  = as.numeric(a) / sum(sub_counts),
    bg_prop      = as.numeric(bg_counts[f]) / sum(bg_counts),
    odds_ratio   = unname(ft$estimate),
    p_value      = ft$p.value
  )
})

enrich_df <- bind_rows(enrich_list) |>
  mutate(p_adj = p.adjust(p_value, method = "BH")) |>
  mutate(log2_or = log2(odds_ratio))

ggplot(enrich_df, aes(x = reorder(feature, log2_or), y = log2_or)) +
  geom_col(aes(fill = log2_or > 0)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "indianred3"), 
                    guide = "none") +
  labs(x = "Feature", y = "log₂(odds ratio)") +
  coord_flip()+
  theme_minimal(base_size = 14)
