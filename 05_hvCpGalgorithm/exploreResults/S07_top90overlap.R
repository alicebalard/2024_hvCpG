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

## Add previous MEs including Maria's results
source(here("05_hvCpGalgorithm/exploreResults/prepPreviousSIV.R"))

## Read vectors saved in R06
overlapLayers <- readRDS(here("gitignore/overlapLayers.RDS"))
overlap.1or2 <- readRDS(here("gitignore/overlap.1or2.RDS"))
topIntersect90 <- readRDS(here("gitignore/topIntersect90.RDS"))

###################################################################################################
## Extract high alpha for test in 05_hvCpGalgorithm/exploreResults/fetalSIV/testFetalSIV_ingp5.R ##
###################################################################################################
length(topIntersect90) # 174494

## Map on arrays
matches <- match(x = topIntersect90, table = dico$chrpos_hg38)

Pos <- dico[na.omit(matches), ]

table(Pos$array)
# 450k    450k and EPIC      EPIC 
# 49          592          514 
# 49+592 = 641 on the 450k array
# 592+514 = 1106 on the EPIC array

saveRDS(Pos, here("05_hvCpGalgorithm/exploreResults/fetalSIV/topIntersect90_pos.RDS"))

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

# Method 1. ClusterProfiler

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

# Method 2. rgreat
# GREAT analysis was conducted using the rGREAT R package (Gu and Hübschmann 2023).
# GREAT performs hypergeometric tests using a foreground and background set and returns 
# annotations of gene sets near the foreground CpGs. 
# The foreground set consists of topIntersect90
# The background set consists of all the sequenced CpGs 
# We reported ontologies after an FDR-adjusted hypergeometric p value<0.05 and an unadjusted hypergeometric p value<0.001. 
## The GREAT settings used were hg38 for the species assembly
background <- makeGRfromMyCpGPos(overlapLayers, "background")
foreground <- makeGRfromMyCpGPos(topIntersect90, "topIntersect90")

system.time(res <- great(gr = foreground, gene_sets = "GO:BP", biomart_dataset = "hg38", background = background, cores = 10))
saveRDS(res, file = here(paste0("05_hvCpGalgorithm/exploreResults/annotations/topIntersect90_rGREAT.RDS")))

#######################
## Enrichement in TE ##
#######################

# UCSC RepeatMasker annotations (Oct2022) for Human (hg38) from AnnotationHub
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("UCSC", "RepeatMasker", "Homo sapiens"))

# Retrieve the desired resource, UCSC RepeatMasker annotations for hg38:
rmskhg38 <- ah[["AH111333"]]

# Filter for ERV1 & ERVK
te_regions <- rmskhg38[mcols(rmskhg38)$repFamily %in% c("ERV1", "ERVK")]

# View summary
table(mcols(te_regions)$repFamily)
length(te_regions)  # Total TE regions

topIntersect90_GR <- makeGRfromMyCpGPos(topIntersect90, "topIntersect90")
overlap.1or2_GR <- makeGRfromMyCpGPos(overlap.1or2, "overlap.1or2")
overlapLayers_GR <- makeGRfromMyCpGPos(overlapLayers, "overlapLayers")

getERVenrichTable <- function(te_regions, mysubtitle, target, nameTarget = "topIntersect90"){
  # Find overlaps
  te_intersectTarget <- te_regions[overlapsAny(te_regions, target, ignore.strand=TRUE)]
  te_overlapLayers <- te_regions[overlapsAny(te_regions, overlapLayers_GR, ignore.strand=TRUE)]
  
  # Function for Fisher's test by repFamily
  fisher_test_erv <- function(family, foreground, background) {
    hits <- sum(mcols(foreground)$repFamily == family)
    bg <- sum(mcols(background)$repFamily == family)
    
    contingency <- matrix(c(
      hits,
      length(target) - hits,
      bg,
      length(overlapLayers) - bg
    ), nrow=2, byrow=TRUE)
    
    rownames(contingency) <- c(nameTarget, "background")
    colnames(contingency) <- c(family, "other")
    
    test <- fisher.test(contingency)
    list(
      family = family,
      contingency = contingency,
      pvalue = test$p.value,
      odds_ratio = test$estimate
    )
  }
  
  # Test ERV1 and ERVK separately
  erv1_test <- fisher_test_erv("ERV1", te_intersectTarget, te_overlapLayers)
  ervk_test <- fisher_test_erv("ERVK", te_intersectTarget, te_overlapLayers)
  
  # Results
  summary_df <- data.frame(
    ERV = c("ERV1", "ERVK"),
    contingency = I(list(erv1_test$contingency, ervk_test$contingency)),
    n_hits = c(erv1_test$contingency[1,1], ervk_test$contingency[1,1]),
    n_total = c(sum(erv1_test$contingency[1,]), sum(ervk_test$contingency[1,])),
    n_bg_hits = c(erv1_test$contingency[2,1], ervk_test$contingency[2,1]),
    n_bg_total = c(sum(erv1_test$contingency[2,]), sum(ervk_test$contingency[2,])),
    p_value = c(erv1_test$pvalue, ervk_test$pvalue),  # Keep raw numbers
    odds_ratio = c(erv1_test$odds_ratio, ervk_test$odds_ratio)
  ) %>%
    mutate(
      pct_hits = n_hits / n_total,
      pct_bg = n_bg_hits / n_bg_total,
      p_value_fmt = case_when(
        p_value < 2.2e-16 ~ "<2.2e-16",
        p_value < 0.001 ~ sprintf("%.2e", p_value),
        TRUE ~ sprintf("%.3f", p_value)
      ))
  
  # Pretty GT table
  summary_df %>%
    dplyr::select(ERV, n_hits, n_total, pct_hits, n_bg_hits, n_bg_total, pct_bg, p_value_fmt, odds_ratio) %>%
    gt() %>%
    fmt_number(columns = c(n_hits, n_total, n_bg_hits, n_bg_total), decimals = 0) %>%
    fmt_percent(columns = c(pct_hits, pct_bg), decimals = 1) %>%
    # NO fmt_scientific() - use pre-formatted column
    fmt_number(columns = odds_ratio, decimals = 3) %>%
    cols_label(
      ERV = "ERV Type",
      n_hits = "Hits",
      n_total = "Total", 
      pct_hits = "%",
      n_bg_hits = "BG Hits",
      n_bg_total = "Total background",
      pct_bg = "% in Background",
      p_value_fmt = "P-value",  # Renamed column
      odds_ratio = "Odds Ratio"
    ) %>%
    tab_style(style = list(cell_fill(color = "lightblue")),
              locations = cells_column_labels()) %>%
    tab_header(title = "ERV Enrichment Analysis (Fisher's Exact Test)",
               subtitle = paste0(mysubtitle, " in ", nameTarget)) %>%
    tab_options(table.font.size = 12, data_row.padding = px(4)) %>% print()
}

getERVenrichTable(te_regions = te_regions, mysubtitle = "exact TE region", target = topIntersect90_GR, nameTarget = "topIntersect90")

## Try for proximal of TE (<= 3kb or <= 10 kb)
getERVenrichTable(te_regions = promoters(te_regions, upstream = 3000, downstream = 3000), 
                  mysubtitle = "proximal TE region (<= 3kb)", target = topIntersect90_GR, nameTarget = "topIntersect90")
getERVenrichTable(te_regions = promoters(te_regions, upstream = 10000, downstream = 10000), 
                  mysubtitle = "proximal TE region (<= 10kb)", target = topIntersect90_GR, nameTarget = "topIntersect90")

## And for CpG p90% in 1 or 2 layers only?
getERVenrichTable(te_regions = te_regions, mysubtitle = "exact TE region", target = overlap.1or2_GR, nameTarget = "overlap.1or2")

## Try for proximal of TE (<= 3kb or <= 10 kb)
getERVenrichTable(te_regions = promoters(te_regions, upstream = 3000, downstream = 3000), 
                  mysubtitle = "proximal TE region (<= 3kb)", target = overlap.1or2_GR, nameTarget = "overlap.1or2")
getERVenrichTable(te_regions = promoters(te_regions, upstream = 10000, downstream = 10000), 
                  mysubtitle = "proximal TE region (<= 10kb)", target = overlap.1or2_GR, nameTarget = "overlap.1or2")

################################
## Genomic positions of top90 ##
################################
## Reload if needed
# system.time(Atlas_dt <- prepAtlasdt("Atlas10X"))
# nrow(Atlas_dt) # 23036026

topIntersect90_dt <- Atlas_dt[name %in% topIntersect90]
head(topIntersect90_dt)

# Compute chromosome centers for x-axis labeling
df2 <- topIntersect90_dt[, .(center = mean(range(pos2, na.rm = TRUE))), by = chr]
df2 <- merge(data.frame(chr = factor(c(1:22, "X", "Y", "M"), levels=as.character(c(1:22, "X", "Y", "M")))),
             df2, by = "chr", all.x = TRUE, sort = TRUE)
df2 <- na.omit(df2)

# Compute chromosome boundaries
df_bounds <- topIntersect90_dt[, .(min_pos = min(pos2, na.rm = TRUE), 
                                   max_pos = max(pos2, na.rm = TRUE)), by = chr]

# Midpoints between chromosomes = where to draw dotted lines
df_bounds[, next_start := data.table::shift(min_pos, n = 1, type = "lead")]
vlines <- df_bounds[!is.na(next_start), .(xintercept = (max_pos + next_start)/2)]


global_dens <- density(topIntersect90_dt$pos2, na.rm = TRUE)
global_df <- data.frame(x = global_dens$x, y = global_dens$y)

ggplot() +
  geom_line(data = global_df,
            aes(x = x, y = y),
            color = "darkred", linewidth = 0.7) +
  geom_vline(data = vlines,
             aes(xintercept = xintercept),
             linetype = 3, color = "grey60") +
  scale_x_continuous(breaks = df2$center,
                     labels = as.character(df2$chr),
                     expand = c(0, 0)) +
  labs(x = "Chromosome", y = "CpG density") +
  theme_minimal(base_size = 14)

######################################
## Test proximity to specific genes ##
######################################
# topIntersect90_GR <- makeGRfromMyCpGPos(vec = topIntersect90, setname = "topIntersect90")
# bed_features <- genomation::readTranscriptFeatures(here("gitignore/hg38_GENCODE_V47.bed"))
# 
# ens = "ENSG00000103126"
# 
# # Axin (human) as a GRanges
# gr_axin <- c(bed_features$promoters[grep(ens, bed_features$promoters$name)],
#              bed_features$exons[grep(ens, bed_features$exons$name)],
#              bed_features$introns[grep(ens, bed_features$introns$name)],
#              bed_features$TSSes[grep(ens, bed_features$TSSes$name)])
# "ENST00000429538.8"
# 287440..352723, complement)
# gr_totest = topIntersect90_GR
# gr_gene = gr_PAX8
# 
# # Find overlaps between CpGs and the gene region
# hits <- findOverlaps(gr_totest, gr_gene)
# 
# # Extract matching rows from original df
# df_hits <- gr_totest[queryHits(hits), ]
# 
# # Add all annotations
# df_hits <- Atlas_dt[match(paste0(df_hits@seqnames, "_", df_hits@ranges), Atlas_dt$name)]
# 
# # Determine limits and breaks
# x_min <- floor(min(df_hits$pos) / 5000) * 5000
# x_max <- ceiling(max(df_hits$pos) / 5000) * 5000
# breaks_seq <- seq(x_min, x_max, by = 5000)
# 
# ggplot(df_hits, aes(x = pos, y = alpha)) +
#   geom_smooth(col = "black") +
#   geom_point(aes(fill = region_type), pch = 21) +
#   theme_minimal(base_size = 14) +
#   ylab("p(hv)") +
#   scale_x_continuous(
#     breaks = breaks_seq,
#     labels = function(x) paste0(formatC(x / 1000, format = "f", digits = 0), "k")
#   ) +
#   xlab("Genomic position (chr2)") +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )