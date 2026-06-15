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
if (!exists("table3layers")) load(here("gitignore/fullTable3layers.Rda"))
totalSiteswGeomMean <- table3layers[!is.na(table3layers$alpha_geomean), ]$chr_pos
top90SNPrm <- table3layers[!is.na(table3layers$alpha_geomean) & 
                             (table3layers$alpha_geomean >= .9), ]$chr_pos

############################################
if (!file.exists(here("gitignore/fullTable3layers.Rda"))){
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
  
  ## Endoderm &  Mesoderm
  if (!file.exists(here::here(
    "B_MultiTissues/dataOut/figures/correlations/correlation_meso-endoFullvsReduced6gp.pdf"))){
    
    ## Use data table to handle large data
    setDT(endo6gp)
    setDT(endo)
    
    x <- endo6gp[, .(name, alpha_6gp = alpha)]; setkey(x, name)
    y <- endo[, .(name, alpha_endo = alpha)]; setkey(y, name)
    
    m <- x[y, nomatch = 0]   # keeps matched names only
    mycor <- cor(m$alpha_6gp, m$alpha_endo, use = "complete.obs")
    set.seed(1234)
    p_endo <- ggplot(m[sample(nrow(m), 100000),], aes(x = alpha_6gp, y = alpha_endo))+
      geom_point(pch = 21, alpha = 0.1) +
      theme_minimal(base_size = 14) +
      ylim(c(0,1)) +
      annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
      labs(title = "Probability of being hypervariable in Loyfer WGBS endoderm cell types",
           subtitle = "(100k random CpG plotted)",
           x = "Pr(hv) calculated on a subset of cell types (N=6)",
           y = "Pr(hv) calculated on all cell types (N=21)")
    
    ## Meso
    setDT(meso6gp)
    setDT(meso)
    
    x <- meso6gp[, .(name, alpha_6gp = alpha)]; setkey(x, name)
    y <- meso[, .(name, alpha_meso = alpha)]; setkey(y, name)
    
    m <- x[y, nomatch = 0]   # keeps matched names only
    mycor <- cor(m$alpha_6gp, m$alpha_meso, use = "complete.obs")
    
    set.seed(1234)
    p_meso <- ggplot(m[sample(nrow(m), 100000),], aes(x = alpha_6gp, y = alpha_meso))+
      geom_point(pch = 21, alpha = 0.1) +
      theme_minimal(base_size = 14) +
      ylim(c(0,1)) +
      annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
      labs(title = "Probability of being hypervariable in Loyfer WGBS mesoderm cell types",
           subtitle = "(100k random CpG plotted)",
           x = "Pr(hv) calculated on a subset of cell types (N=6)", 
           y = "Pr(hv) calculated on all cell types (N=19)")
    
    ggplot2::ggsave(
      filename = here::here(paste0("B_MultiTissues/dataOut/figures/correlations/correlation_meso-endoFullvsReduced6gp.pdf")),
      plot = cowplot::plot_grid(p_endo, p_meso, labels = c("A", "B")), width = 17, height = 8)
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
  
  if (sum(grepl("chr", seqlevels(table3layers))) == 0){
    seqlevels(table3layers) <- paste0("chr", seqlevels(table3layers))
  }
  
  ### SAVED ###
  save(table3layers, file = here("gitignore/fullTable3layers.Rda"))
}

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

sets <- list(
  mQTLcontrols = makeGRfromMyCpGPos(vec = mQTLcontrols_hg38, setname = "mQTLcontrols"),
  HarrisSIV = HarrisSIV_hg38_GR,
  VanBaakSIV = VanBaakSIV_hg38_GR,
  VanBaakESS = VanBaakESS_hg38_GR,
  KesslerSIV = KesslerSIV_GRanges_hg38,
  CoRSIV = corSIV_GRanges_hg38,
  hvCpG = DerakhshanhvCpGs_hg38_GR
)

# ── ME overlap ────────────────────────────────────────────
MEsetdt            <- make_MEsetdt(sets, geomMeanGR)
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
# mQTLcontrols  0.211 0.00634 69972    0.199    0.224
# CoRSIV        0.315 0.00141 69972    0.312    0.318
# HarrisSIV     0.361 0.00974 69972    0.342    0.380
# hvCpG         0.526 0.00655 69972    0.513    0.539
# KesslerSIV    0.404 0.00684 69972    0.391    0.418
# VanBaakESS    0.474 0.01070 69972    0.453    0.495
# VanBaakSIV    0.633 0.02110 69972    0.591    0.674
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

############################################################ 
## Now, summarising for the regions, one value per region ##
############################################################
MEsetdt_regionMean <- make_MEsetdt_regionMean(sets, geomMeanGR)
MEsetdt_regionMean[, ME := relevel(factor(ME), ref = "mQTLcontrols")]

# Set controls as baseline
MEsetdt_regionMean[, ME := relevel(factor(ME), ref = "mQTLcontrols")]

p1 <- ggplot(MEsetdt_regionMean, aes(x = ME, y = alpha_geomean)) +
  geom_jitter(data = MEsetdt,
              aes(fill=ME), pch=21, size = 3, alpha = .05)+
  geom_violin(aes(col=ME))+
  geom_boxplot(aes(col=ME), width = .1) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylab("Pr(hv) (geometric mean)")

## Statistical comparisons of alpha between MEs

# Fit the model with controls as baseline
fit <- lm(alpha_geomean ~ ME, data = MEsetdt_regionMean)

# Get estimated marginal means and contrasts vs baseline
emm <- emmeans(fit, ~ ME)
contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = "mQTLcontrols", adjust = "sidak") %>%
  as.data.frame()

emm
# ME           emmean      SE    df lower.CL upper.CL
# mQTLcontrols  0.211 0.00405 18292    0.203    0.219
# CoRSIV        0.262 0.00224 18292    0.258    0.267
# HarrisSIV     0.361 0.00621 18292    0.349    0.373
# hvCpG         0.526 0.00418 18292    0.517    0.534
# KesslerSIV    0.402 0.00989 18292    0.382    0.421
# VanBaakESS    0.474 0.00680 18292    0.461    0.487
# VanBaakSIV    0.633 0.01350 18292    0.606    0.659
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

pdf(here("B_MultiTissues/dataOut/figures/alphaComparisonBetweenMEtypes_meanPerRegion.pdf"), 
    width = 14, height = 4)
cowplot::plot_grid(p1,p2, rel_widths = c(1, .8))
dev.off()

############################
## decay curve per ME set ##
thresholds <- seq(10, 90, by = 10) / 100

# one row per threshold per ME set
prop_table <- rbindlist(lapply(thresholds, function(thr) {
  MEsetdt[, .(
    proportion = mean(alpha_geomean > thr, na.rm = TRUE),
    n_above    = sum(alpha_geomean > thr, na.rm = TRUE),
    n_total    = .N
  ), by = ME][, threshold := thr]
}))

# plot
p1 <- plot_decay_curve(MEsetdt,            "Based on one Pr(HV) per CpG")

# one row per threshold per ME set
prop_table <- rbindlist(lapply(thresholds, function(thr) {
  MEsetdt_regionMean[, .(
    proportion = mean(alpha_geomean > thr, na.rm = TRUE),
    n_above    = sum(alpha_geomean > thr, na.rm = TRUE),
    n_total    = .N
  ), by = ME][, threshold := thr]
}))

# plot
p2 <- plot_decay_curve(MEsetdt_regionMean, "Based on one Pr(HV) per region")

pdf(here("B_MultiTissues/dataOut/figures/decayAlpha.pdf"), width = 14, height = 4)
p1 + p2 + plot_layout(guides = "collect")
dev.off()

# ── Save key objects for S07 ──────────────────────────────────────────────────
saveRDS(MEsetdt,            here("gitignore/MEsetdt.rds"))
saveRDS(MEsetdt_regionMean, here("gitignore/MEsetdt_regionMean.rds"))
saveRDS(geomMeanGR,         here("gitignore/geomMeanGR.rds"))

##############################
## Save the top alpha > 90% ##
##############################

totalSiteswGeomMean <- table3layers[!is.na(table3layers$alpha_geomean), ]$chr_pos

top90SNPrm <- table3layers[!is.na(table3layers$alpha_geomean) & 
                             (table3layers$alpha_geomean >= .9), ]$chr_pos

message(paste0("Total CpG sites with non NA geometric mean: ", length(totalSiteswGeomMean)))
message(paste0("Total top90 CpG sites: ", length(top90SNPrm), " (", 
               round(length(top90SNPrm)/length(totalSiteswGeomMean)*100,2), "% of total)"))
# Total CpG sites with non NA geometric mean: 21522541
# 196333 (0.91% of total)

####################################################################################
## Test enrichement of the most likely germ layer-universal hvCpG in previous MEs ##
####################################################################################
if (!exists("listGR")){
  listGR <- list(top90 = makeGRfromMyCpGPos(vec = top90SNPrm, setname = "top90CpGs"),
                 allButTop90 = makeGRfromMyCpGPos(
                   setdiff(totalSiteswGeomMean, top90SNPrm), "allButTop90"))}

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
  listGR$allButTop90$featureType <- ifelse(allButTop90Anno@members[, "prom"] == 1, "promoter",
                                           ifelse(allButTop90Anno@members[, "exon"] == 1, "exon",
                                                  ifelse(allButTop90Anno@members[, "intron"] == 1, "intron",
                                                         "intergenic")))
  
  ## 1. Build counts for subset vs background 
  
  # Feature levels in fixed order
  feat_levels <- c("promoter", "exon", "intron", "intergenic")
  
  # Count CpGs per feature
  bg_counts   <- table(factor(listGR$allButTop90$featureType,   levels = feat_levels))
  sub_counts  <- table(factor(listGR$top90$featureType, levels = feat_levels))
  
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
  enrich_df
  # feature subset_n     bg_n subset_prop    bg_prop odds_ratio       p_value         p_adj    log2_or
  # 1   promoter    17168  3171023  0.08744327 0.14869137  0.5461601  1.000000e+00  1.000000e+00 -0.8726043
  # 2       exon     7081  1194839  0.03606628 0.05602679  0.6281852  1.000000e+00  1.000000e+00 -0.6707382
  # 3     intron   116681 12025864  0.59430152 0.56390072  1.1341610 3.351803e-165 6.703607e-165  0.1816254 ***
  # 4 intergenic    55403  4934482  0.28218893 0.23138113  1.3093824  0.000000e+00  0.000000e+00  0.3888865 ***
  
  ## Enrichment in intron & intergenic regions of the top 90% hvCpGs ***
}

###############################
## Test GO of top candidates ##
###############################

# totalSiteswGeomMean <- table3layers[!is.na(table3layers$alpha_geomean), ]$chr_pos
# top90SNPrm <- table3layers[!is.na(table3layers$alpha_geomean) & 
#                              (table3layers$alpha_geomean >= .9), ]$chr_pos

# Method. ClusterProfiler

## 1. Keep CpGs in regions where at least 2 CpGs are in 50bp distance to each other
## 2. annotate with associated genes (in gene body or +/- 10kb from TSS)
## 3. run GO term enrichment with clusterProfiler::enrichGO

minimum_CpG_per_cluster = 2

## Create universe
universe <- annotateCpGs_txdb(
  clusterCpGs(totalSiteswGeomMean, max_gap = 50, min_size = minimum_CpG_per_cluster),
  tss_window = 10000)

print(paste0("Gene universe contains ", length(universe), " genes"))
## "Gene universe contains 32717 genes"

## Annotate 
# resAnnot_top90SNPrm <- CpG_GO_pipeline(
#   top90SNPrm, universe = universe, 
#   max_gap = 50, min_size = minimum_CpG_per_cluster, tss_window = 10000)
# # Reduced from 196333 to 7111 clustered CpGs
# # Found 1641 Entrez genes
# # Running GO enrichment...

# Original (no length control)
resAnnot_top90SNPrm <- CpG_GO_pipeline(
  top90SNPrm, universe = universe,
  max_gap = 50, min_size = minimum_CpG_per_cluster, tss_window = 10000)
# Reduced from 196333 to 7111 clustered CpGs
# 2169 genes were dropped because they have exons located on both strands of the same reference sequence or on more than one
# reference sequence, so cannot be represented by a single genomic range.
# Found 1641 Entrez genes
# Running GO enrichment...

# Length-controlled
resAnnot_top90SNPrm_lenCtrl <- CpG_GO_pipeline_lengthControlled(
  top90SNPrm, universe = universe,
  max_gap = 50, min_size = minimum_CpG_per_cluster, tss_window = 10000,
  control_length = TRUE)
# Found 1641 Entrez genes
# Controlling for gene length...
# Median gene length — foreground: 5,514 bp, universe: 3,100 bp, ratio: 1.78
# Length-matched universe: 16779 genes (was 32717)
# Running GO enrichment...

# Compare which terms survive
terms_original    <- resAnnot_top90SNPrm$BP@result %>% filter(p.adjust < 0.05) %>% pull(ID)
terms_lenCtrl     <- resAnnot_top90SNPrm_lenCtrl$BP@result %>% filter(p.adjust < 0.05) %>% pull(ID)

cat("BP terms before length control:", length(terms_original), "\n")
cat("BP terms after  length control:", length(terms_lenCtrl),  "\n")
cat("Terms lost:", length(setdiff(terms_original, terms_lenCtrl)), "\n")
cat("Terms gained:", length(setdiff(terms_lenCtrl, terms_original)), "\n")

df_all <- purrr::imap_dfr(resAnnot_top90SNPrm_lenCtrl, function(er, ont_name) {
  if (is.null(er) || nrow(er@result) == 0) return(tibble())
  as_tibble(er@result) |> 
    mutate(group_raw = "top90SNPrm", ontology = ont_name)
}) |> bind_rows()

df_all <- df_all |>
  mutate(
    group = "top90SNPrm",
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

write.csv(df_sig, file = here("B_MultiTissues/dataOut/df_sig_GOtop90SNPrm.csv"),
          quote = F, row.names = F)

# Plot
p <- ggplot(df_sig, aes(x = group, y = Description)) +
  geom_point(aes(size = FoldEnrichment, color = p.adjust), alpha = 0.9) +
  scale_size_continuous(name = "Fold Enrichment", 
                        range = c(1.5, 8), breaks = c(2, 2.5, 3, 3.5)) +  
  scale_color_viridis_c(name = "FDR", option = "plasma", direction = -1) +
  facet_wrap(ontology ~ ., scales = "free", space = "free_x") +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = NULL, 
       title = "GO Enrichment: top90SNPrm (FDR < 0.05)") +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold")
  )

print(p)

#######################
## Enrichement in TE ##
#######################

# UCSC RepeatMasker annotations (Oct2022) for Human (hg38) from AnnotationHub
ah <- AnnotationHub()
query(ah, c("UCSC", "RepeatMasker", "Homo sapiens"))

# Retrieve the desired resource, UCSC RepeatMasker annotations for hg38:
rmskhg38 <- ah[["AH111333"]]

# Filter for ERV1 & ERVK
te_regions <- rmskhg38[mcols(rmskhg38)$repFamily %in% c("ERV1", "ERVK")]

# View summary
table(mcols(te_regions)$repFamily)
length(te_regions)  # Total TE regions

top90SNPrm_GR <- makeGRfromMyCpGPos(top90SNPrm, "top90SNPrm")
totalSiteswGeomMean_GR <- makeGRfromMyCpGPos(totalSiteswGeomMean, "totalSiteswGeomMean")

# Strict background = non-hvCpG sites only
bg_only_GR <- makeGRfromMyCpGPos(totalSiteswGeomMean[!totalSiteswGeomMean %in% top90SNPrm], "bg_only")

# Test: of all CpGs, are hvCpGs more likely to fall inside an ERV1/ERVK region than background CpGs?
fisher_test_erv <- function(family, target, mybackground,
                            nameTarget = "foreground") {
  
  te_family_gr <- te_regions[mcols(te_regions)$repFamily == family]  # fine as-is
  
  fg_in  <- sum(overlapsAny(target,       te_family_gr, ignore.strand = TRUE))
  fg_out <- length(target) - fg_in
  bg_in  <- sum(overlapsAny(mybackground, te_family_gr, ignore.strand = TRUE))
  bg_out <- length(mybackground) - bg_in
  
  contingency <- matrix(c(fg_in, fg_out,
                          bg_in, bg_out),
                        nrow = 2, byrow = TRUE,
                        dimnames = list(c(nameTarget, "background"),
                                        c("in_ERV", "not_in_ERV")))
  test <- fisher.test(contingency, alternative = "greater")
  list(family      = family,
       contingency = contingency,
       pvalue      = test$p.value,
       odds_ratio  = test$estimate)
}

# Call with explicit nameTarget
erv1_test <- fisher_test_erv("ERV1", top90SNPrm_GR, bg_only_GR,
                             nameTarget = "top90SNPrm")
erv1_test
# $contingency
#             in_ERV not_in_ERV
# top90SNPrm   7265     189068
# background 632733   20693475
# 
# $pvalue
# [1] 1.30206e-75
# 
# $odds_ratio
# odds ratio 
# 1.256695
ervk_test <- fisher_test_erv("ERVK", top90SNPrm_GR, bg_only_GR,
                             nameTarget = "top90SNPrm")
ervk_test
# $contingency
#             in_ERV not_in_ERV
# top90SNPrm   1043     195290
# background  84328   21241880
# 
# $pvalue
# [1] 7.030883e-20
# 
# $odds_ratio
# odds ratio 
# 1.345317 
