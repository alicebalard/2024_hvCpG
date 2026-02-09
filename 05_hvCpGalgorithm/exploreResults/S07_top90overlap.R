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

retest = FALSE # takes some time to run

if (retest == TRUE){
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
}

###############################
## Test GO of top candidates ##
###############################

## 1. Keep CpGs in regions where at least 2 CpGs are in 50bp distance to each other
## 2. annotate with associated genes (in gene body or +/- 10kb from TSS)
## 3. run GO term enrichment with clusterProfiler::enrichGO

minimum_CpG_per_cluster = 2

## Create universe
universe <- annotateCpGs_txdb(
  clusterCpGs(overlapLayers, max_gap = 50, min_size = minimum_CpG_per_cluster),
  tss_window = 10000
)

print(paste0("Gene universe contains ", length(universe), " genes"))
## "Gene universe contains 25159 genes"

## Annotate 
resAnnot_topIntersect90 <- CpG_GO_pipeline(
  topIntersect90, universe = universe, 
  max_gap = 50, min_size = minimum_CpG_per_cluster, tss_window = 10000)
## Reduced from 174494 to 6906 clustered CpGs

df_all <- purrr::imap_dfr(resAnnot_topIntersect90, function(er, ont_name) {
  if (is.null(er) || nrow(er@result) == 0) return(tibble())
  as_tibble(er@result) |> 
    mutate(group_raw = "topIntersect90", ontology = ont_name)
}) |> bind_rows()

df_all <- df_all |>
  mutate(
    group = "topIntersect90",
    group = factor(group),
    ontology = factor(ontology, levels = c("BP", "MF", "CC"))
  )

# Filter significant terms
df_sig <- df_all |>
  filter(!is.na(p.adjust) & p.adjust < 0.05) |> 
  filter(Count > 10 & FoldEnrichment > 2) 

# Reorder by enrichment strength
df_sig <- df_sig |>
  group_by(ontology) |>
  mutate(Description = fct_reorder(Description, FoldEnrichment, .desc = TRUE)) |>
  ungroup()

# Plot
p <- ggplot(df_sig, aes(x = group, y = Description)) +
  geom_point(aes(size = FoldEnrichment, color = p.adjust), alpha = 0.9) +
  scale_size_continuous(name = "Fold Enrichment", 
                        range = c(1.5, 8), breaks = c(2, 2.5, 3, 3.5)) +  
  scale_color_viridis_c(name = "FDR", option = "plasma", direction = -1) +
  facet_wrap(ontology ~ ., scales = "free", space = "free_x") +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = NULL, 
       title = "GO Enrichment: topIntersect90 (FDR < 0.05)") +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold")
  )

print(p)
