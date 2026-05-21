#####################################
## One run per developmental layer ##
#####################################

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

## To avoid re-running everything:
load(here("gitignore/fullTable3layers.Rda"))

############################################

endo <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_12_endo.rds"))
meso <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_meso.rds"))
ecto <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_14_ecto.rds"))
Atlas_dt <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_atlas_general.rds"))

## With only 6 groups (to see if power)
endo6gp <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_12_2_endo6gp.rds"))
meso6gp <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_2_meso6gp.rds"))

endoGR <- GRanges(seqnames = endo$chr,
                  ranges = IRanges(start = endo$pos, end = endo$pos),
                  alpha_endo = endo$alpha)

ectoGR <- GRanges(seqnames = ecto$chr,
                  ranges = IRanges(start = ecto$pos, end = ecto$pos),
                  alpha_ecto = ecto$alpha)

mesoGR <- GRanges(seqnames = meso$chr,
                  ranges = IRanges(start = meso$pos, end = meso$pos),
                  alpha_meso = meso$alpha)

allLayersGR <- GRanges(seqnames = Atlas_dt$chr,
                       ranges = IRanges(start = Atlas_dt$pos, end = Atlas_dt$pos),
                       alpha_allLayers = Atlas_dt$alpha)

################################
## Test: are 6 groups enough? ##
################################

## Endoderm
if (!file.exists(here::here(
  "B_MultiTissues/dataOut/figures/correlations/correlation_endoFullvsReduced6gp.pdf"))){
  ## Use data table to handle large data
  setDT(endo6gp)
  setDT(endo)
  
  x <- endo6gp[, .(name, alpha_6gp = alpha)]; setkey(x, name)
  y <- endo[, .(name, alpha_endo = alpha)]; setkey(y, name)
  
  m <- x[y, nomatch = 0]   # keeps matched names only
  mycor <- cor(m$alpha_6gp, m$alpha_endo, use = "complete.obs")
  set.seed(1234)
  p <- ggplot(m[sample(nrow(m), 100000),], aes(x = alpha_6gp, y = alpha_endo))+
    geom_point(pch = 21, alpha = 0.1) +
    theme_minimal(base_size = 14) +
    ylim(c(0,1)) +
    annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
    labs(title = "Probability of being hypervariable in WGBS atlas endoderm cell types",
         subtitle = "(100k random CpG plotted)",
         x = "Pr(hv) calculated on a subset of cell types (N=6)",
         y = "Pr(hv) calculated on all cell types (N=21)")
  
  ggplot2::ggsave(
    filename = here::here(paste0("B_MultiTissues/dataOut/figures/correlations/correlation_endoFullvsReduced6gp.pdf")),
    plot = p, width = 8, height = 8
  )
}

## Mesoderm
if (!file.exists(here::here(
  "B_MultiTissues/dataOut/figures/correlations/correlation_mesoFullvsReduced6gp.pdf"))){
  setDT(meso6gp)
  setDT(meso)
  
  x <- meso6gp[, .(name, alpha_6gp = alpha)]; setkey(x, name)
  y <- meso[, .(name, alpha_meso = alpha)]; setkey(y, name)
  
  m <- x[y, nomatch = 0]   # keeps matched names only
  mycor <- cor(m$alpha_6gp, m$alpha_meso, use = "complete.obs")
  
  set.seed(1234)
  p <- ggplot(m[sample(nrow(m), 100000),], aes(x = alpha_6gp, y = alpha_meso))+
    geom_point(pch = 21, alpha = 0.1) +
    theme_minimal(base_size = 14) +
    ylim(c(0,1)) +
    annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
    labs(title = "Probability of being hypervariable in WGBS atlas mesoderm cell types",
         subtitle = "(100k random CpG plotted)",
         x = "Pr(hv) calculated on a subset of cell types (N=6)", 
         y = "Pr(hv) calculated on all cell types (N=19)")
  
  ggplot2::ggsave(
    filename = here::here(paste0("B_MultiTissues/dataOut/figures/correlations/correlation_mesoFullvsReduced6gp.pdf")),
    plot = p, width = 8, height = 8
  )
}

####################################################################
## Create a table with all CpG sites & pr(hv) for each germ layer ##
####################################################################

## 1. Create union of all unique CpG positions
table3layers <- union(union(allLayersGR, union(ectoGR, mesoGR)), endoGR)

## 2. Use findOverlaps to map alpha values back
# we want endoGR[i] -> table3layers[endoHits[[i]]] for each i, etc.

endoHits <- findOverlaps(endoGR, table3layers, select = "first")
ectoHits <- findOverlaps(ectoGR, table3layers, select = "first")
mesoHits <- findOverlaps(mesoGR, table3layers, select = "first")
allLayersHits <- findOverlaps(allLayersGR, table3layers, select = "first")

# initialize columns with NA
mcols(table3layers)$alpha_endo <- NA_real_
mcols(table3layers)$alpha_ecto <- NA_real_
mcols(table3layers)$alpha_meso <- NA_real_
mcols(table3layers)$alpha_allLayers <- NA_real_

# copy only the hits
mcols(table3layers)$alpha_endo[endoHits]   <- mcols(endoGR)$alpha_endo
mcols(table3layers)$alpha_ecto[ectoHits]   <- mcols(ectoGR)$alpha_ecto
mcols(table3layers)$alpha_meso[mesoHits]   <- mcols(mesoGR)$alpha_meso
mcols(table3layers)$alpha_allLayers[allLayersHits]   <- mcols(allLayersGR)$alpha_allLayers

## Add chr_pos column to identify positions
table3layers$chr_pos <- paste0("chr", table3layers@seqnames, "_", table3layers@ranges@start)

## add a geometric mean between the 3 layers
table3layers$alpha_geomean <- exp(rowMeans(
  log(cbind(table3layers$alpha_endo, 
            table3layers$alpha_ecto, 
            table3layers$alpha_meso)),
  na.rm = FALSE))

### SAVED ###
save(table3layers, file = here("gitignore/fullTable3layers.Rda"))

df <- as.data.frame(table3layers)

p <- ggplot(df[sample(nrow(df), 100000),],
            aes(x = alpha_geomean, y = alpha_allLayers)) +
  geom_point(pch = 21, alpha = 0.1) +
  geom_abline(slope = 1) +
  theme_minimal(base_size = 14) +
  labs(title = "Pr(hv) using all layers vs geometric mean of each 3 layers",
       x = "Geometric mean Pr(hv) on the three layers",
       y = "Pr(hv) on all cell types",
       subtitle = "(100k random CpG plotted)")

table(is.na(df$alpha_geomean))
# FALSE     TRUE 
# 21522541  3027679 
table(is.na(df$alpha_meso))
table(is.na(df$alpha_endo))
table(is.na(df$alpha_ecto))

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/geomMean_vs_all.png"),
  plot = p, width = 7, height = 7,
  dpi = 300, bg = "white")

############################################
## Plot Manhattan of geometric mean alpha ##
############################################

table3layersdt <- as.data.table(table3layers)
names(table3layersdt)[names(table3layersdt) %in% "alpha_geomean"] <- "alpha"
names(table3layersdt)[names(table3layersdt) %in% "seqnames"] <- "chr"
names(table3layersdt)[names(table3layersdt) %in% "start"] <- "pos"

# Compute cumulative position offsets for Manhattan plot
setorder(table3layersdt, chr, pos)

offsets <- table3layersdt[, .(max_pos = max(pos, na.rm = TRUE)), by = chr]
offsets[, cum_offset := c(0, head(cumsum(as.numeric(max_pos)), -1))]

table3layersdt <- merge(table3layersdt,
                        offsets[, .(chr, cum_offset)], 
                        by = "chr", all.x = TRUE, sort = FALSE)

## Mark group membership in dt
table3layersdt[, group := NA_character_]
table3layersdt[chr_pos %in% DerakhshanhvCpGs_hg38, group := "hvCpG_Derakhshan"]
table3layersdt[chr_pos %in% mQTLcontrols_hg38, group := "mQTLcontrols"]

# Convert to integer/numeric if not already
table3layersdt[, cum_offset := as.numeric(cum_offset)]
table3layersdt[, pos2 := pos + cum_offset]

table(is.na(table3layersdt$alpha))
# FALSE     TRUE 
# 21522541  3027679 

plotManhattan3 <- plotManhattanFromdt(table3layersdt, plotDerakhshan = FALSE)
ggplot2::ggsave(
  filename = here::here(
    "B_MultiTissues/dataOut/figures/Manhattan/ManhattanAlphaPlot_atlas_geomMean.png"),
  plot = plotManhattan3, width = 14, height = 4,
  dpi = 300, bg = "white")

#######################################################
## Check alpha for the different MEs: are they high? ##
#######################################################

# Build GRanges from geometric mean
geomMeanGR <- GRanges(seqnames = table3layers@seqnames,
                      ranges = IRanges(start = table3layers@ranges@start, 
                                       end = table3layers@ranges@start),
                      alpha_geomean = table3layers$alpha_geomean)

geomMeanGR <- geomMeanGR[!is.na(geomMeanGR$alpha_geomean)]

# Fix chromosome names in geomMeanGR (1 -> chr1)
seqlevels(geomMeanGR) <- paste0("chr", seqlevels(geomMeanGR))

sets <- list(
  mQTLcontrols = makeGRfromMyCpGPos(vec = mQTLcontrols_hg38, setname = "mQTLcontrols"),
  HarrisSIV = HarrisSIV_hg38_GR,
  VanBaakESS = VanBaakESS_hg38_GR,
  KesslerSIV = KesslerSIV_GRanges_hg38,
  CoRSIV = corSIV_GRanges_hg38,
  hvCpG = DerakhshanhvCpGs_hg38_GR
)

# Now do overlap join for each set
MEsetdt <- rbindlist(lapply(names(sets), function(nm) {
  hits <- findOverlaps(sets[[nm]], geomMeanGR)
  data.table(
    alpha_geomean = geomMeanGR$alpha_geomean[subjectHits(hits)],
    ME    = nm
  )
}))

MEsetdt <- na.omit(MEsetdt) ## 63078

# Set controls as baseline
MEsetdt[, ME := relevel(factor(ME), ref = "mQTLcontrols")]

p1 <- ggplot(MEsetdt, aes(x = ME, y = alpha_geomean)) +
  geom_jitter(data = MEsetdt,
              aes(fill=ME), pch=21, size = 3, alpha = .05)+
  geom_violin(aes(col=ME))+
  geom_boxplot(aes(col=ME), width = .1) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylab("Pr(hv) (geometric mean)")

## Statistical comparisons of alpha between MEs

# Fit the model with controls as baseline
fit <- lm(alpha_geomean ~ ME, data = MEsetdt)

# Get estimated marginal means and contrasts vs baseline
emm <- emmeans(fit, ~ ME)
contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = "mQTLcontrols", adjust = "sidak") %>%
  as.data.frame()

emm
# ME           emmean      SE    df lower.CL upper.CL
# mQTLcontrols  0.211 0.00635 69709    0.199    0.224
# CoRSIV        0.315 0.00141 69709    0.312    0.318
# HarrisSIV     0.361 0.00974 69709    0.342    0.380
# hvCpG         0.526 0.00656 69709    0.513    0.539
# KesslerSIV    0.404 0.00684 69709    0.391    0.418
# VanBaakESS    0.474 0.01070 69709    0.453    0.495
# 
# Confidence level used: 0.95 

contrasts <- contrasts %>%
  mutate(ME = contrast,  # rename for clarity
         lower = estimate - 1.96*SE,
         upper = estimate + 1.96*SE)

p2 <- ggplot(contrasts, aes(x = ME, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    y = "Difference in Pr(hv) (geometric mean) vs mQTLcontrols",
    x = "",
    title = "Comparison of ME groups to mQTLcontrols",
    subtitle = "lm with multiple comparison correction (Sidak)"
  ) +
  theme_minimal()

pdf(here("B_MultiTissues/dataOut/figures/alphaComparisonBetweenMEtypes.pdf"), 
    width = 14, height = 4)
cowplot::plot_grid(p1,p2, rel_widths = c(1, .8))
dev.off()

##############################
## Save the top alpha > 90% ##
##############################

totalSiteswGeomMean <- table3layers[!is.na(table3layers$alpha_geomean), ]$chr_pos

top90SNPrm <- table3layers[!is.na(table3layers$alpha_geomean) & 
                             (table3layers$alpha_geomean >= .9), ]$chr_pos

message(paste0("Total CpG sites with non NA geometric mean: ", length(totalSiteswGeomMean)))
message(paste0("Total top90 CpG sites: ", length(top90SNPrm), " (", 
               round(length(top90SNPrm)/length(totalSiteswGeomMean)*100,2), "% of total)"))

saveRDS(top90SNPrm, file = here("gitignore/top90SNPrm.RDS"))
saveRDS(totalSiteswGeomMean, file = here("gitignore/alButtop90SNPrm.RDS"))

####################################################################################
## Test enrichement of the most likely germ layer-universal hvCpG in previous MEs ##
####################################################################################
if (!exists("top90SNPrm")){
  top90SNPrm <- readRDS(overlap, file = here("gitignore/top90SNPrm.RDS"))}

listGR <- list(top90 = makeGRfromMyCpGPos(vec = top90SNPrm, setname = "top90CpGs"),
               allButTop90 = makeGRfromMyCpGPos(
                 setdiff(totalSiteswGeomMean, top90SNPrm), "allButTop90"))

# ---- Run it (ME sets in putativeME_GR$set will be tested separately)
res_quadrants <- test_enrichment_quadrants(listGR, putativeME_GR, me_col = "set")

# Order quadrants within each facet by log2OR
res_plot2 <- res_quadrants %>%
  mutate(
    log2OR = log2(odds_ratio),
    signif  = p_adj_BH < 0.05
  ) %>%
  dplyr::group_by(CpG_set) %>%
  mutate(quadrant_ord = reorder(quadrant, log2OR)) %>%
  ungroup()

plot <- ggplot(res_plot2, aes(x = quadrant_ord, y = log2OR, fill = signif)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("grey", "black")) +
  labs(
    x = NULL,
    y = expression(log[2]~"(odds ratio)"),
    title = "ME enrichment by group (vs other group)",
    subtitle = "2x2 Fisher's exact test "
  ) +
  facet_wrap(~ CpG_set, scales = "free_x", nrow = 1) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none" , ## if all significant
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

pdf(here("B_MultiTissues/dataOut/figures/top90CpGsEnrichME.pdf"), width = 8, height = 6)
plot
dev.off()

############################################################
## How many of each putative ME is actually in the top90? ##
############################################################

# Overlaps of ME ranges with CpG sites (top90 and allButTop90)
hits_top <- findOverlaps(putativeME_GR, listGR$top90, ignore.strand = TRUE)
hits_all <- findOverlaps(putativeME_GR, listGR$allButTop90, ignore.strand = TRUE)

# Logical flags per ME range (query index)
overlap_top <- logical(length(putativeME_GR))
overlap_top[unique(queryHits(hits_top))] <- TRUE

overlap_all <- logical(length(putativeME_GR))
overlap_all[unique(queryHits(hits_all))] <- TRUE

# Build a small data.frame with one row per ME range
df_me <- as.data.frame(mcols(putativeME_GR)) |>
  mutate(
    in_top90       = overlap_top,
    in_allButTop90 = overlap_all
  )

# Summaries per set
summary_df <- df_me |>
  group_by(set) |>
  summarise(
    n_total = n(),
    n_in_top90 = sum(in_top90),
    pc_in_top90 = n_in_top90/n_total*100,
    n_in_allButTop90 = sum(in_allButTop90),
    pc_in_allButTop90 = n_in_allButTop90/n_total*100,
    n_in_both = sum(in_top90 & in_allButTop90),   # non 0 if region rather than CpG
    pc_in_both = n_in_both/n_total*100,
    n_in_neither = n_total - n_in_top90 - n_in_allButTop90 + n_in_both,
    pc_in_neither = n_in_neither/n_total*100)

## Format pretty
summary_df %>%
  mutate(across(starts_with("pc_"), ~ scales::percent(.x / 100, accuracy = 0.1))) %>%
  gt() %>%
  fmt_number(columns = starts_with("n_"), decimals = 0) %>%
  cols_label(
    set = "Putative ME set",
    n_total = "Total CpG sites or regions",
    n_in_top90 = "N in top 90% CpGs",
    pc_in_top90 = "%",
    n_in_allButTop90 = "N in all but top 90% CpGs",
    pc_in_allButTop90 = "%",
    n_in_both = "N in both groups",
    pc_in_both = "%",
    n_in_neither = "N CpGs not covered",
    pc_in_neither = "%"
  ) %>% 
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_column_labels()
  ) %>%
  tab_options(
    table.font.size = 13,
    data_row.padding = px(3)
  )

## Screenshot saved in figures/topCpGsEnrichME_table.png

######################################################################
## Check enrichement of telomeres and centromeres for high geomMean ##
######################################################################

# Centromeres (from cytoBand - acen bands)
# wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz \
# | zcat | awk '$5=="acen"' > centromeres_hg38.bed

# 1. Parse coordinates from name vectors 
parse_cpg_names <- function(names_vec) {
  dt <- data.table(name = names_vec)
  dt[, chr := sub("^(chr[^_]+)_.*", "\\1", name)]        # "chr1"
  dt[, pos := as.integer(sub("^chr[^_]+_(\\d+)$", "\\1", name))]
  dt[, pos_end := pos]
  dt
}

getEnrichCentroTelo <- function(threshold = 0.90){
  top <- table3layers[!is.na(table3layers$alpha_geomean) & 
                        (table3layers$alpha_geomean >= threshold), ]$chr_pos
  hv_dt <- parse_cpg_names(top)
  total_dt <- parse_cpg_names(totalSiteswGeomMean)  # all background sites
  
  # 2. Load regions
  centro <- fread(here("gitignore/centromeres_hg38.bed"),
                  col.names = c("chr", "start", "end", "band", "stain"))
  centro[, region := "centromere"]
  
  # Add 1Mb subtelomeric buffer using chrom sizes
  chrom_sizes <- fread("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz",
                       col.names = c("chr", "size", "file"))
  chrom_sizes  <- chrom_sizes[chr %in% paste0("chr", c(1:22, "X", "Y"))]
  SUBTELO_DIST <- 1e6
  
  subtelo <- rbind(
    chrom_sizes[, .(chr, start = 0L, end = as.integer(SUBTELO_DIST), region = "subtelomere")],
    chrom_sizes[, .(chr, start = as.integer(size - SUBTELO_DIST), end = size, region = "subtelomere")]
  )
  
  regions <- rbind(
    centro[, .(chr, start, end, region)],
    subtelo[, .(chr, start, end, region)]
  )
  setkey(regions, chr, start, end)
  
  # 3. Overlap function 
  get_region_hits <- function(dt, regions) {
    setkey(dt, chr, pos, pos_end)
    hits <- foverlaps(dt, regions,
                      by.x = c("chr", "pos", "pos_end"),
                      by.y = c("chr", "start", "end"),
                      type = "within", nomatch = NULL)
    unique(hits$name)  # CpG names overlapping any region
  }
  
  hv_in_centro    <- get_region_hits(copy(hv_dt),    regions[region == "centromere"])
  hv_in_subtelo   <- get_region_hits(copy(hv_dt),    regions[region == "subtelomere"])
  bg_in_centro    <- get_region_hits(copy(total_dt), regions[region == "centromere"])
  bg_in_subtelo   <- get_region_hits(copy(total_dt), regions[region == "subtelomere"])
  
  # 4. Contingency tables + Fisher test 
  enrich_test <- function(hv_in, bg_in, hv_all, bg_all, label) {
    a <- length(hv_in)                        # hvCpG in region
    b <- length(hv_all) - a                   # hvCpG outside
    c <- length(bg_in)                        # background in region
    d <- length(bg_all) - c                   # background outside
    
    mat <- matrix(c(a, b, c, d), nrow = 2,
                  dimnames = list(c("in_region", "outside"),
                                  c("hvCpG", "background")))
    
    ft  <- fisher.test(mat, alternative = "greater")
    pct_hv <- round(100 * a / length(hv_all), 2)
    pct_bg <- round(100 * c / length(bg_all), 2)
    fold   <- round(pct_hv / pct_bg, 2)
    
    cat("\n──", label, "──\n")
    cat("  hvCpGs in region:     ", a, "/", length(hv_all),
        paste0("(", pct_hv, "%)"), "\n")
    cat("  Background in region: ", c, "/", length(bg_all),
        paste0("(", pct_bg, "%)"), "\n")
    cat("  Fold enrichment:      ", fold, "\n")
    cat("  Fisher p (one-sided): ", ft$p.value, "\n")
    cat("  Odds ratio:           ", round(ft$estimate, 2), "\n")
  }
  
  enrich_test(hv_in_centro,  bg_in_centro,  top, totalSiteswGeomMean, "Centromere")
  enrich_test(hv_in_subtelo, bg_in_subtelo, top, totalSiteswGeomMean, "Subtelomere (1Mb)")
}

getEnrichCentroTelo(0.9)
## For 90%, no enrichment
# ── Centromere ──
# hvCpGs in region:      1871 / 196333 (0.95%) 
# Background in region:  190918 / 21522541 (0.89%) 
# Fold enrichment:       1.07 
# Fisher p (one-sided):  0.001123841 
# Odds ratio:            1.08 
# 
# ── Subtelomere (1Mb) ──
# hvCpGs in region:      4040 / 196333 (2.06%) 
# Background in region:  482013 / 21522541 (2.24%) 
# Fold enrichment:       0.92 
# Fisher p (one-sided):  1 
# Odds ratio:            0.92 

getEnrichCentroTelo(0.8)
## enriched in centromeres 
# ── Centromere ──
# hvCpGs in region:      3701 / 282564 (1.31%) 
# Background in region:  190918 / 21522541 (0.89%) 
# Fold enrichment:       1.47 
# Fisher p (one-sided):  1.526111e-109 
# Odds ratio:            1.48 
# 
# ── Subtelomere (1Mb) ──
# hvCpGs in region:      6261 / 282564 (2.22%) 
# Background in region:  482013 / 21522541 (2.24%) 
# Fold enrichment:       0.99 
# Fisher p (one-sided):  0.8037475 
# Odds ratio:            0.99 

getEnrichCentroTelo(0.70)
## enriched 
# ── Centromere ──
# hvCpGs in region:      8246 / 415760 (1.98%) 
# Background in region:  190918 / 21522541 (0.89%) 
# Fold enrichment:       2.22 
# Fisher p (one-sided):  0 
# Odds ratio:            2.26 
# 
# ── Subtelomere (1Mb) ──
# hvCpGs in region:      9858 / 415760 (2.37%) 
# Background in region:  482013 / 21522541 (2.24%) 
# Fold enrichment:       1.06 
# Fisher p (one-sided):  9.630318e-09 
# Odds ratio:            1.06 

############################################################################
## Test in B_MultiTissues/03_exploreResults/fetalSIV/testFetalSIV_ingp5.R ##
############################################################################
length(top90SNPrm) # 196.333

## Map on arrays
matches <- match(x = top90SNPrm, table = dico$chrpos_hg38)

Pos <- dico[na.omit(matches), ]

table(Pos$array)
# 450k 450k and EPIC          EPIC 
# 80          1026           718

# 80+1026 = 1106 on the 450k array
# 1026+718 = 1744 on the EPIC array

#################################
## Test enrichment of features ##
#################################

retest = TRUE # takes some time to run

if (retest == TRUE){
  # Import bed file
  bed_features <- genomation::readTranscriptFeatures(here("gitignore/hg38_GENCODE_V47.bed"))
  
  # Annotate CpGs with features
  topAnno <- genomation::annotateWithGeneParts(
    target  = listGR$top90,
    feature = bed_features
  )
  listGR$top90$featureType <- ifelse(topAnno@members[, "prom"] == 1, "promoter",
                                     ifelse(topAnno@members[, "exon"] == 1, "exon",
                                            ifelse(topAnno@members[, "intron"] == 1, "intron",
                                                   "intergenic")))
  allButTop90Anno <- genomation::annotateWithGeneParts(
    target  = listGR$allButTop90,
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
  sub_counts  <- table(factor(top90SNPrm_GR$featureType, levels = feat_levels))
  
  # Combine into a 2x4 contingency table
  cont_tab <- rbind(
    subset    = as.numeric(sub_counts),
    background = as.numeric(bg_counts)
  )
  colnames(cont_tab) <- feat_levels
  cont_tab
  
  enrich_list <- lapply(feat_levels, function(f) {
    # 2x2 table for feature f vs not‑f
    a <- sub_counts[f]                     # subset in feature f
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
  
  p <- ggplot(enrich_df, aes(x = reorder(feature, log2_or), y = log2_or)) +
    geom_col(aes(fill = log2_or > 0)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "indianred3"), 
                      guide = "none") +
    labs(x = "Feature", y = "log₂(odds ratio)") +
    coord_flip()+
    theme_minimal(base_size = 14)
  
  p
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
saveRDS(res, file = here(paste0("B_MultiTissues/03_exploreResults/annotations/topIntersect90_rGREAT.RDS")))

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

top90SNPrm_GR <- makeGRfromMyCpGPos(topIntersect90, "topIntersect90")
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



##############################
## new candidate locus Matt ##
##############################

## Matt's data are in hg19
dataMatt <- readxl::read_xlsx(here("gitignore/DEGCAGS_intersect_repeats_Alice.xlsx"))

# --- Download and import hg19 → hg38 chain file ---
chain_dir <- here("B_MultiTissues/dataIn")
chain_gz <- file.path(chain_dir, "hg19ToHg38.over.chain.gz")
chain_file <- file.path(chain_dir, "hg19ToHg38.over.chain")

if (!file.exists(chain_file)) {
  message("⬇️  Downloading UCSC liftOver chain file...")
  dir.create(chain_dir, showWarnings = FALSE, recursive = TRUE)
  download.file(
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    destfile = chain_gz,
    quiet = TRUE
  )
  R.utils::gunzip(chain_gz, destname = chain_file, remove = FALSE)
}
chain <- import.chain(chain_file)

# --- Liftover (hg19 → hg38) ---
dataMatt_gr <- GRanges(
  seqnames = dataMatt$chromosome,
  ranges = IRanges(start = dataMatt$start, end = dataMatt$end),
  alpha_endo = NA,  alpha_meso = NA,  alpha_ecto = NA,  alpha_all = NA,
  hg19_chr = dataMatt$chromosome, hg19_start = dataMatt$start, hg19_end = dataMatt$end,
  `%change` = dataMatt$`%change`, padj = dataMatt$padj,
  TE_chromosome = dataMatt$TE_chromosome, TE_start = dataMatt$TE_start,
  TE_end = dataMatt$TE_end, TE_family = dataMatt$TE_family,TE_type = dataMatt$TE_type
)

mapped <- liftOver(dataMatt_gr, chain)

# Keep one-to-one mappings only
keep <- lengths(mapped) == 1
dataMatt_hg38_gr <- unlist(mapped[keep])

# Initialise the column with NA first
dataMatt_hg38_gr$alpha_geomean <- NA_real_
# Find overlapping ranges
overlaps <- findOverlaps(dataMatt_hg38_gr, geomMeanGR)
dataMatt_hg38_gr[queryHits(overlaps),]$alpha_geomean <- 
  geomMeanGR[subjectHits(overlaps),]$alpha_geomean

ggplot(as.data.frame(dataMatt_hg38_gr), aes(x = "all", y=alpha_geomean)) +
  geom_violin() +
  geom_boxplot(width = .2) +
  geom_jitter() +
  theme_minimal(base_size = 14) 






#####################################################################################
## NEXT Find correlated regions = at least 5CpGs in 50bp which are hypervariable (>50%) ##

