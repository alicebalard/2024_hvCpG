#################################################
## Plot results of algorithm ran on atlas data ##
#################################################
## New run: SNPs MAF 1% excluded

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

## Load array results
if (!exists("resArray")) {
  resArray <- readRDS(here("B_MultiTissues/dataOut/resArraymean0.9p10.8p0.RDS"))
}

#######################
## Data in WGBS atlas:

## from the CS cluster: 
## /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/output_atlas_general/sample_metadata.tsv
sample_groups <- read.table(
  here("B_MultiTissues/resultsDir_gitIgnored/Atlas/atlas_general/sample_metadata.tsv"), 
  sep = "\t", header = T)

sample_groups %>% group_by(dataset) %>% summarise(n = n()) %>% 
  ggplot(aes(x = n)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of number of samples per dataset",
    x = "Number of samples",
    y = "Count of datasets"
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 1))

SupTab1_Loyfer2023 <- read.csv(here("B_MultiTissues/dataIn/SupTab1_Loyfer2023.csv"))
SupTab1_Loyfer2023$group <- paste(SupTab1_Loyfer2023$Source.Tissue, SupTab1_Loyfer2023$Cell.type, sep = " - ")
table(table(SupTab1_Loyfer2023$group)[table(SupTab1_Loyfer2023$group) >=3])
# 3  4  5  6 10 
# 33  9  2  1  1 

##################################
## Save all data in RDS objects ##
##################################
savePrepedAtlasFile <- function(
    file, p0, p1,
    res       = "fullres_",
    atlas_dir = here("B_MultiTissues/resultsDir_gitIgnored/Atlas"),
    out_dir   = here("gitignore/resultsAtlasPrepared")) {
  
  # `file` is the bare subdir name, e.g. "12_endo"
  search_dir <- file.path(atlas_dir, file)                      # where the batches live
  out_path   <- file.path(out_dir,
                          paste0(res, p0, "p0_", p1, "p1_", file, ".rds"))
  
  # ── Completeness check ────────────────────────────────────────────────────
  pattern   <- paste0("^results_.*", p0, "p0_", p1, "p1\\.rds$")
  rds_files <- base::dir(search_dir, pattern = pattern,
                         recursive = TRUE, full.names = TRUE)
  n_files   <- length(rds_files)
  
  if (n_files == 0) {
    message("SKIP ", file, " [", p0, "/", p1, "] - no matching files.")
    return(invisible(NULL))
  }
  
  # Check 1: file count matches highest Atlas_batch number
  batch_dirs <- unique(dirname(rds_files))
  batch_nums <- as.integer(regmatches(batch_dirs,
                                      regexpr("(?<=Atlas_batch)\\d+", batch_dirs, perl = TRUE)))
  expected_n <- max(batch_nums, na.rm = TRUE)
  
  if (n_files != expected_n) {
    message("SKIP ", file, " - expected ", expected_n,
            " files but found ", n_files, ".")
    return(invisible(NULL))
  }
  
  # Check 2: last batch must not be exactly 250000 CpGs (truncation sentinel)
  last_file   <- sort(rds_files) |> tail(1)                     # <- fixed sort
  cpg_in_name <- as.integer(regmatches(last_file,
                                       regexpr("(?<=_)\\d+(?=CpGs)", last_file, perl = TRUE)))
  
  if (!is.na(cpg_in_name) && cpg_in_name == 250000) {
    message("SKIP ", file, " - last batch has 250000 CpGs (truncated?): ",
            basename(last_file))
    return(invisible(NULL))
  }
  
  message("OK   ", file, " - ", n_files, "/", expected_n,
          " files, last batch has ", cpg_in_name, " CpGs.")
  
  if (file.exists(out_path)) {
    message("     already prepared - skipping."); return(invisible(NULL))
  }
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  system.time(Atlas_dt <- prepAtlasdt(file, p0, p1, atlas_dir = atlas_dir))
  saveRDS(Atlas_dt, file = out_path)
  message("     saved: ", out_path)
}

subdirs <- list.files(here("B_MultiTissues/resultsDir_gitIgnored/Atlas/"))
subdirs <- subdirs[!grepl("^PREVIOUS", subdirs)]      # bare names only

for (subdir in subdirs) {
  savePrepedAtlasFile(subdir, p0 = "0_8",  p1 = "0_65")
  savePrepedAtlasFile(subdir, p0 = "0_8",  p1 = "0_9")
  savePrepedAtlasFile(subdir, p0 = "0_55", p1 = "0_65")
}

prevDir <- here("B_MultiTissues/resultsDir_gitIgnored/Atlas/PREVIOUSalgosumloglikoverinds")
subdirs <- list.files(prevDir)
subdirs <- subdirs[!grepl("^Atlas10X_tissueAnalysis", subdirs)]

for (subdir in subdirs) {
  savePrepedAtlasFile(subdir, p0 = "0_8",  p1 = "0_65", res = "previousres_sumlog_", atlas_dir = prevDir)
  savePrepedAtlasFile(subdir, p0 = "0_8",  p1 = "0_9",  res = "previousres_sumlog_", atlas_dir = prevDir)
  savePrepedAtlasFile(subdir, p0 = "0_55", p1 = "0_65", res = "previousres_sumlog_", atlas_dir = prevDir)
}

## New July 2026 (post bug, post sum --> mean)
# Saved: /home/alice/Documents/Research/GIT/2024_hvCpG/gitignore/resultsAtlasPrepared/

################################################################################
## Compare algorithm with sum of with mean per dataset of the logliks         ##
################################################################################

Atlas_dt <- readRDS(
  here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_9p1_atlas_general.rds"))
nrow(Atlas_dt) # 21522541

Atlas_dt_sum <- readRDS(
  here("gitignore/resultsAtlasPrepared/previousres_sumlog_0_8p0_0_9p1_atlas_general.rds"))
nrow(Atlas_dt_sum) # 21522541

# Set key if not already set
setkey(Atlas_dt, name)
setkey(Atlas_dt_sum, name)

# Merge all three
merged_dt <- Atlas_dt[Atlas_dt_sum,
                      .(name, alpha_mean = alpha, alpha_sum = i.alpha),
                      on = "name", nomatch = NULL]

set.seed(1234)
merged_dt[sample(.N, 100000)] |>
  ggplot(aes(x = alpha_mean, y = alpha_sum)) +
  geom_point(alpha = 0.05, size = 0.3) +   # small/transparent for density
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  labs(x = "Pr(hv) atlas general - mean", y = "Pr(hv) atlas general - sum") +
  theme_minimal()

## Ranking is preserved. The tight linear band confirms both settings agree on which
# CpGs are most variable

###################
## rmMultSamples ## 
###################

# Some individuals have multiple cells sampled. Does that affect our results? NOPE
plotrmtest <- makeCompPlot(minplot = 1000000, # plot all
  X = readRDS(here::here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_9p1_atlas_general.rds")),
  Y = readRDS(here::here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_9p1_02_rmMultSamples.rds")),
  whichAlphaX = "alpha",
  whichAlphaY = "alpha",
  title = "Effect of keeping only one cell type per individual (1M CPGs plotted)",
  xlab = "Pr(hv) calculated on WGBS atlas datasets",
  ylab = "Pr(hv) calculated on WGBS atlas datasets keeping one sample/individual only")

ggplot2::ggsave(
  filename = here::here(paste0("B_MultiTissues/dataOut/figures/script03/rmMultipleSamples.pdf")),
  plot = plotrmtest, width = 10, height = 10)

################################################################################
## Load layer-specific analyses                                               ##
################################################################################

endo <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_9p1_12_endo.rds"))
meso <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_9p1_13_meso.rds"))
ecto <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_9p1_14_ecto.rds"))

################################
## Test: are 6 groups enough? ##
################################

endo6gp <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_9p1_12_2_endo6gp.rds"))
meso6gp <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_9p1_13_2_meso6gp.rds"))

if (!file.exists(here::here(
  "B_MultiTissues/dataOut/figures/script03/correlation_endomesoFullvsReduced6gp.pdf"))){
  ## Use data table to handle large data
  setDT(endo6gp)
  setDT(endo)
  
  x <- endo6gp[, .(name, alpha_6gp = alpha)]; setkey(x, name)
  y <- endo[, .(name, alpha_endo = alpha)]; setkey(y, name)
  
  m <- x[y, nomatch = 0]   # keeps matched names only
  mycor <- cor(m$alpha_6gp, m$alpha_endo, use = "complete.obs")
  set.seed(1234)
  p1 <- ggplot(m[sample(nrow(m), 100000),], aes(x = alpha_6gp, y = alpha_endo))+
    geom_point(pch = 21, alpha = 0.1) +
    geom_smooth(method = "lm") +
    theme_minimal(base_size = 14) +
    ylim(c(0,1)) +
    annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
    labs(title = "Probability of being hypervariable in WGBS atlas endoderm cell types",
         subtitle = "(100k random CpG plotted)",
         x = "Pr(hv) calculated on a subset of cell types (N=6)",
         y = "Pr(hv) calculated on all cell types (N=21)")
  
  setDT(meso6gp)
  setDT(meso)
  
  x <- meso6gp[, .(name, alpha_6gp = alpha)]; setkey(x, name)
  y <- meso[, .(name, alpha_meso = alpha)]; setkey(y, name)
  
  m <- x[y, nomatch = 0]   # keeps matched names only
  mycor <- cor(m$alpha_6gp, m$alpha_meso, use = "complete.obs")
  
  p2 <- ggplot(m[sample(nrow(m), 100000),], aes(x = alpha_6gp, y = alpha_meso))+
    geom_point(pch = 21, alpha = 0.1) +
    geom_smooth(method = "lm") +
    theme_minimal(base_size = 14) +
    ylim(c(0,1)) +
    annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
    labs(title = "Probability of being hypervariable in WGBS atlas mesoderm cell types",
         subtitle = "(100k random CpG plotted)",
         x = "Pr(hv) calculated on a subset of cell types (N=6)", 
         y = "Pr(hv) calculated on all cell types (N=19)")
  
  ggplot2::ggsave(
    filename = here::here(paste0("B_MultiTissues/dataOut/figures/script03/correlation_endomesoFullvsReduced6gp.pdf")),
    plot = cowplot::plot_grid(p1, p2, labels = c("A", "B")), width = 17, height = 8)
  
  rm(x, y, m, mycor, p1, p2)
}

################################################################################
## Compare algo ran on all datasets with geometric mean                       ##
################################################################################

endoGR <- GRanges(seqnames = endo$chr,
                  ranges = IRanges(start = endo$pos, end = endo$pos),
                  alpha_endo = endo$alpha)

ectoGR <- GRanges(seqnames = ecto$chr,
                  ranges = IRanges(start = ecto$pos, end = ecto$pos),
                  alpha_ecto = ecto$alpha)

mesoGR <- GRanges(seqnames = meso$chr,
                  ranges = IRanges(start = meso$pos, end = meso$pos),
                  alpha_meso = meso$alpha)

Atlas_dt <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_9p1_atlas_general.rds"))

allLayersGR <- GRanges(seqnames = Atlas_dt$chr,
                       ranges = IRanges(start = Atlas_dt$pos, end = Atlas_dt$pos),
                       alpha_allLayers = Atlas_dt$alpha)

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
# save(table3layers, file = 
#        here(paste0("gitignore/fullTable3layers_", format(Sys.Date(), "%d_%m_%y"), ".Rda")))

##############################
## Save the top alpha > 99% ##
##############################

totalSiteswGeomMean <- table3layers[!is.na(table3layers$alpha_geomean), ]$chr_pos

top99SNPrm <- table3layers[!is.na(table3layers$alpha_geomean) & 
                             (table3layers$alpha_geomean >= .99), ]$chr_pos

message(paste0("Total CpG sites with non NA geometric mean: ", length(totalSiteswGeomMean)))
message(paste0("Total top99 CpG sites: ", length(top99SNPrm), " (", 
               round(length(top99SNPrm)/length(totalSiteswGeomMean)*100,2), "% of total)"))
# Total CpG sites with non NA geometric mean: 21522541
# Total top90 CpG sites: 385370 (1.79% of total)

saveRDS(top99SNPrm, file = here("gitignore/top99SNPrm_july26.RDS"))
## To use for testFetalSIV_ingp5.R

################################################################################
################################## CHECKPOINT ##################################
################################################################################
load(here(paste0("gitignore/fullTable3layers_23_07_26.Rda")))

df <- as.data.frame(table3layers)

set.seed(1234)
p <- ggplot(df[sample(nrow(df), 100000),],
            aes(x = alpha_geomean, y = alpha_allLayers)) +
  geom_point(pch = 21, alpha = 0.1) +
  geom_abline(slope = 1) +
  theme_minimal(base_size = 14) +
  labs(title = "Pr(hv) using all layers vs geometric mean of each 3 layers",
       x = "Geometric mean Pr(hv) on the three layers",
       y = "Pr(hv) on all cell types",
       subtitle = "(100k random CpG plotted)")

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/script03/geomMean_vs_all.png"),
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

plotManhattangeomMean <- plotManhattanFromdt(table3layersdt, plotDerakhshan = FALSE, 
                                             centro = centro)
saveRDS(plotManhattangeomMean, file = here("gitignore/plotManhattangeomMean.RDS"))

############################################
## Compare with Maria's results           ##
############################################

###################################################################
## Calculate proba hvCpG minus matching control: is it always +? ##

data <- read.table(
  here("B_MultiTissues/01_dataPrep/prepDatasetsMaria_LSHTMserver/cistrans_GoDMC_hvCpG_matched_control.txt"), header = T)

x = dico$chrpos_hg38[match(data$hvCpG_name, dico$CpG)]
y = dico$chrpos_hg38[match(data$controlCpG_name, dico$CpG)]

# Build mapping from hvCpG -> control
pairs <- data.frame(
  hvCpG = x,
  control = y,
  stringsAsFactors = FALSE
)

# Merge hvCpG alphas
res <- table3layersdt[!is.na(group)]

hv_alpha <- res[, c("chr_pos", "alpha")]
colnames(hv_alpha) <- c("hvCpG", "alpha_hvCpG")

# Merge control alphas
ctrl_alpha <- res[, c("chr_pos", "alpha")]
colnames(ctrl_alpha) <- c("control", "alpha_control")

# Join everything
merged <- pairs %>%
  left_join(hv_alpha, by = "hvCpG") %>%
  left_join(ctrl_alpha, by = "control") %>%
  mutate(diffAlpha=alpha_hvCpG-alpha_control)

merged <- merged %>%
  mutate(chr = str_extract(hvCpG, "^chr[0-9XYM]+"))%>%
  filter(!is.na(diffAlpha))

# DifferenceOfProbabilityForhvCpG-matching_controlInAtlas
pdiffhv_controls <- ggplot(merged, aes(x="diff", y=diffAlpha))+
  geom_jitter(data=merged[merged$diffAlpha>=0,], col="black", alpha=.5)+
  geom_jitter(data=merged[merged$diffAlpha<0,], fill="yellow",col="black",pch=21, alpha=.5)+
  geom_violin(width=.5, fill = "grey", alpha=.8) +
  geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
  theme_minimal(base_size = 14)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), title = element_text(size=10))+
  ggtitle("Pr(hvCpG) minus Pr(matching control) in atlas")+
  ylab("Difference of probability")

plotManhattangeomMean2 <- plotManhattanFromdt(table3layersdt, plotDerakhshan = TRUE, centro = centro)

compPlotArrayAtlas <- makeCompPlot(
  X = resArray,
  Y = here::here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_9p1_atlas_general.rds"),
  whichAlphaX = "alpha",
  whichAlphaY = "alpha",
  title = "Array_vs_Atlas",
  xlab = "Pr(hv) calculated on array datasets",
  ylab = "Pr(hv) calculated on WGBS atlas datasets", minplot = 24000000)

compPlotArrayAtlasMESO <- makeCompPlot(
  X = resArray,
  Y = here::here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_9p1_13_meso.rds"),
  whichAlphaX = "alpha",
  whichAlphaY = "alpha",
  title = "Array_vs_Atlas",
  xlab = "Pr(hv) calculated on array datasets",
  ylab = "Pr(hv) calculated on WGBS atlas datasets MESODERM derived", minplot = 24000000)

ggplot2::ggsave(
  filename = here::here(
    "B_MultiTissues/dataOut/figures/script03/CompareWithResultsDerakhshan.png"),
  plot = cowplot::plot_grid(
    cowplot::plot_grid(plotManhattangeomMean2, pdiffhv_controls,
                       rel_widths = c(3,1), labels = c("A", "B")),
    cowplot::plot_grid(compPlotArrayAtlas, compPlotArrayAtlasMESO, labels = c("C", "D"), ncol = 2),
    nrow = 2, rel_heights = c(1,2)),
  width = 15, height = 10, dpi = 300, bg = "white")

##########################################
## What are the gaps in Manhattan plot? ##
# Compute the gap between consecutive CpGs on the same chromosome
table3layersdt[, gap := pos - data.table::shift(pos), by = chr]

# Identify large gaps (>= 500k bp)
gaps_dt <- table3layersdt[gap >= 500000, .(
  chr,
  gap_start = data.table::shift(pos),
  gap_end = pos,
  gap_size = gap
)]

# Drop first NA (since shift introduces one per chromosome)
gaps_dt[!is.na(gap_size)]

# chr gap_start   gap_end gap_size
# <fctr>     <int>     <int>    <int>
#   1:      1        NA 123979940  1102558
# 2:      1 123979940 124778153   798022
# 3:      1 124778153 143184605 18000029
# 4:      2 143184605  91406021  1003516
# 5:      5  91406021  49591807  2282986
# 6:      9  49591807  60518620 15013475
# 7:     14  60518620  18223731  2127270
# 8:     16  18223731  46380693  8100344
# 9:     19  46380693  26051961  1143713
# 10:     19  26051961  27240939  1188969
# 11:     21  27240939  12966132  2151768
# 12:     22  12966132  15154526  2250197
# 13:      X  15154526  60188774   701154
# 14:      Y  60188774   5585195   541347
# 15:      Y   5585195   7252770   503589
# 16:      Y   7252770   8945914   807694
# 17:      Y   8945914   9591421   645507
# 18:      Y   9591421  12916601   778842
# 19:      Y  12916601  13638567   721948
# 20:      Y  13638567  15801485  2076134
# 21:      Y  15801485  17064099  1262614
# 22:      Y  17064099  20534098  3131980
# 23:      Y  20534098  26432094  5897996
# 24:      Y  26432094  56677947 30006397

####################################
## Mitochondrial DNAm variability ##
####################################
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-023-09541-9?utm_source=chatgpt.com
## Near absence of 5mC, so expected

ggplot() +
  geom_point(data = table3layersdt[table3layersdt$chr == "M",], 
             aes(x = pos2, y = alpha),
             color = "black", size = 1, alpha = .5)+
  theme_classic() + theme(legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,.1)) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

###################################
## Y chromosome DNAm variability ##
###################################
ggplot() +
  geom_point(data = table3layersdt[table3layersdt$chr == "Y",], 
             aes(x = pos2, y = alpha),
             color = "black", size = 1, alpha = .5)+
  theme_classic() + theme(legend.position = "none") +
  geom_hline(yintercept = .7, linetype = 3)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

# High alpha in 3 regions

#####################################################
## Test enrichment of features for high alpha      ##
#####################################################

# Create GRanges of 23036026 CpGs
gr_cpg <- GRanges(
  seqnames = paste0("chr", table3layersdt$chr),
  ranges = IRanges(start = table3layersdt$pos, end = table3layersdt$pos),
  alpha = table3layersdt$alpha
)

# Import bed file
bed_features <- genomation::readTranscriptFeatures(here("gitignore/hg38_GENCODE_V47.bed"))

# restrict to autosomes and chr X
gr_cpg <- gr_cpg[gr_cpg@seqnames %in% c(paste0("chr", 1:22), "chrX"),]

# Annotate CpGs and see which regions have higher alpha (takes long)
anno_result <- genomation::annotateWithGeneParts(
  target = gr_cpg, feature = bed_features)

anno_result@perc.of.OlapFeat
# promoter     exon   intron 
# 97.66919 79.59063 92.61569 

## Add info from annotation to our GRange object
gr_cpg$featureType <- ifelse(anno_result@members[, "prom"] == 1, "promoter",
                             ifelse(anno_result@members[, "exon"] == 1, "exon",
                                    ifelse(anno_result@members[, "intron"] == 1, "intron", "intergenic")))

mcols(gr_cpg) %>% as.data.frame() %>%
  dplyr::group_by(featureType) %>%
  dplyr::summarise(meanAlpha = mean(alpha, na.rm=T),
                   medianAlpha = median(alpha, na.rm=T))
# # A tibble: 4 × 3
# featureType meanAlpha medianAlpha
# <chr>           <dbl>       <dbl>
# 1 exon           0.0979 0.00000134 
# 2 intergenic     0.152  0.00000231 ********* The higghest
# 3 intron         0.106  0.00000137 
# 4 promoter       0.105  0.000000902

dt <- as.data.table(mcols(gr_cpg))[!is.na(alpha)]
dt[, band := cut(alpha, c(-Inf, 0.1, 0.5, 0.7, 0.9, Inf),
                 labels = c("~0","0.1–0.5","0.5–0.7","0.7–0.9", "≥0.9"))]

frac <- dt[, .N, by = .(featureType, band)][, pct := 100*N/sum(N), by = featureType]

pfeatures <- ggplot(frac[band != "~0"], aes(featureType, pct, fill = band)) +
  geom_col(position = "dodge") +
  labs(y = "% of CpGs", x = NULL, fill = "pr(hv)") +
  scale_fill_brewer(palette = "OrRd") + theme_minimal(base_size = 14)

pfeatures

# ## ── Statistical testing ──────────────────────────────────────────────────────
# intergenic and intronic blocks have the highest hvCpG rates

## Work from `dt` (which carries featureType + alpha), NOT table3layersdt,
## and add the genomic coordinates so we can build independence blocks.
dt <- as.data.table(mcols(gr_cpg))
dt[, `:=`(chr = as.character(seqnames(gr_cpg)),
          pos = start(gr_cpg))]
dt <- dt[!is.na(alpha)]

## Define hypervariable by a threshold, then aggregate to ~independent units.
## 24M CpGs are NOT independent (adjacent CpGs co-methylate); testing per-CpG
## gives meaningless p-values. Collapse to 100 kb blocks so observations are
## roughly independent, and summarise the QUANTITY OF INTEREST — the fraction of
## hvCpGs per block — rather than mean(alpha), which buries the signal under the
## 60% near-zero pile.
dt[, hv := alpha >= 0.9]
dt[, block := paste0(chr, "_", pos %/% 1e5)]

agg <- dt[, .(frac_hv    = mean(hv),      # fraction of hvCpGs in this block
              mean_alpha = mean(alpha),   # kept for comparison
              n          = .N),
          by = .(featureType, block)][n >= 20]   # drop sparse blocks (unstable rates)

## Omnibus: does the hvCpG rate differ across features?
## Rank-based (no distributional assumption); n = number of blocks, not CpGs.
kruskal.test(frac_hv ~ featureType, data = agg)

## Effect size — THE ACTUAL ANSWER. At this n every p-value is ~0; epsilon^2
## says how much of the variance in block hv-rate the feature explains.
rstatix::kruskal_effsize(agg, frac_hv ~ featureType)

## Pairwise, Holm-corrected: which features differ (e.g. is intergenic > others?)
rstatix::dunn_test(agg, frac_hv ~ featureType, p.adjust.method = "holm")

## Report medians per feature alongside the test (the numbers to quote):
agg[, .(median_frac_hv = median(frac_hv),
        mean_frac_hv   = mean(frac_hv),
        n_blocks       = .N), by = featureType][order(-median_frac_hv)]

# Across 77,794 100-kb blocks, the fraction of hypervariable CpGs (p(hv) ≥ 0.9) differed by genomic feature (Kruskal–Wallis χ²(3) = 7454, p < 2.2×10⁻¹⁶), with a moderate effect size (η²_H = 0.096) — i.e. genomic feature explains ~10% of the between-block variance in hvCpG rate. Intergenic and intronic regions showed the highest hvCpG rates (median 1.7% and 1.4% of CpGs per block), exonic regions the lowest (median 0%), and promoters an intermediate but strongly right-skewed distribution (median 0.6%, mean 1.8%), consistent with a minority of promoter regions — likely CpG-island promoters — carrying disproportionately many hypervariable sites. All pairwise differences were significant after Holm correction (Dunn test, all p_adj < 10⁻⁵⁰), though given the large number of blocks, the effect size and the rate estimates are more informative than the p-values.

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
# $family
# [1] "ERV1"
# 
# $contingency
# in_ERV not_in_ERV
# top90SNPrm  14501     370869
# background 625497   20511674
# 
# $pvalue
# [1] 8.552909e-173
# 
# $odds_ratio
# odds ratio 
# 1.282191 

ervk_test <- fisher_test_erv("ERVK", top90SNPrm_GR, bg_only_GR,
                             nameTarget = "top90SNPrm")
ervk_test
# $family
# [1] "ERVK"
# 
# $contingency
# in_ERV not_in_ERV
# top90SNPrm   2213     383157
# background  83158   21054013
# 
# $pvalue
# [1] 3.783231e-62
# 
# $odds_ratio
# odds ratio 
# 1.462265 

# Calculate proportions
prop_data <- data.frame(
  TE_Family = c("ERV1", "ERVK"),
  Group = rep(c("top90SNPrm", "Background"), each = 2),
  In_TE_pct = c(14501/385370*100, 625497/21137171*100,
                2213/385370*100, 83158/21137171*100),
  OR = c(1.28, 1.28, 1.46, 1.46),
  P_val = c("8.55e-173", "8.55e-173", "3.78e-62", "3.78e-62")
)

# Main bar plot
plotTE <- ggplot(prop_data, aes(x = TE_Family, y = In_TE_pct, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.2f%%", In_TE_pct)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(
    name = "CpG group",  
    values = c("top90SNPrm" = "steelblue", "Background" = "gray70"),
    labels = c("Background (the remaining CpGs)", "Top CpGs with Pr(hv)>.9")  # Custom legend labels (was in wrong order)
  ) +
  labs(
    title = "Percentage of CpG Sites in Transposable Elements",
    y = "Percentage in TE",
    x = "",
    fill = ""  # ← Remove this or it will override scale_fill_manual(name=)
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 11),  # Optional: style legend title
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) + 
  annotate("text", x = 0.7, y = max(prop_data$In_TE_pct) * 0.9,
           label = "ERV1: OR = 1.28, p < 0.0001\nERVK: OR = 1.46, p < 0.0001",
           size = 3.5, fontface = "italic", color = "darkred")

#################################
## plot manhattan, feature, TE ##
#################################
ggplot2::ggsave(
  filename = here::here(
    "B_MultiTissues/dataOut/figures/script03/MappingVariability.png"),
  plot = plot_grid(plotManhattangeomMean,
                   plot_grid(
                     pfeatures, plotTE, 
                     ncol= 2,
                     labels = c("B", "C")
                   ), nrow = 2,
                   labels = c("A", "")),
  width = 15, height = 10,
  dpi = 300, bg = "white")

# **************************************************************************** #
###########################################
## Compare our results with previous MEs ##
###########################################

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
  VanBaakSIV = VanBaakSIV_hg38_GR,
  VanBaakESS = VanBaakESS_hg38_GR,
  KesslerSIV = KesslerSIV_GRanges_hg38,
  GunasekaraCorSIV = corSIV_GRanges_hg38,
  DerakhshanhvCpGs = DerakhshanhvCpGs_hg38_GR
)

group_cols <- c(
  "background"         = "#999999",
  "mQTLcontrols"       = "#000000",
  "HarrisSIV"          = RColorBrewer::brewer.pal(8, "Set2")[1],
  "KesslerSIV"         = RColorBrewer::brewer.pal(8, "Set2")[2],
  "DerakhshanhvCpGs"   = RColorBrewer::brewer.pal(8, "Set2")[3],
  "GunasekaraCorSIV"   = RColorBrewer::brewer.pal(8, "Set2")[4],
  "VanBaakESS"         = RColorBrewer::brewer.pal(8, "Set2")[5],
  "top99SNPrm"         = RColorBrewer::brewer.pal(8, "Set2")[6],
  "VanBaakSIV"         = RColorBrewer::brewer.pal(8, "Set2")[7]
)

# ── ME overlap ────────────────────────────────────────────
MEsetdt <- make_MEsetdt(sets, geomMeanGR)

MEsetdt <- na.omit(MEsetdt) ## 69979

# Set controls as baseline
MEsetdt[, ME := relevel(factor(ME), ref = "mQTLcontrols")]

## Shared colour palette, built the same way plot_decay_curve() already does internally
# me_levels    <- levels(factor(MEsetdt$ME))
# other_levels <- setdiff(me_levels, "mQTLcontrols")
# set2_cols    <- RColorBrewer::brewer.pal(max(length(other_levels), 3), "Set2")
# meColours    <- c(mQTLcontrols = "black",
#                   setNames(set2_cols[seq_along(other_levels)], other_levels))

## Statistical comparisons of alpha between MEs
fit <- lm(alpha_geomean ~ ME, data = MEsetdt)
emm <- emmeans(fit, ~ ME)
contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = "mQTLcontrols", adjust = "sidak") %>%
  as.data.frame()

contrasts <- contrasts %>%
  mutate(ME = contrast,
         ME_name = sub(" - mQTLcontrols$", "", contrast),   # match colour to the ME group being compared
         lower = estimate - 1.96 * SE,
         upper = estimate + 1.96 * SE)

pMEalpha <- ggplot(MEsetdt, aes(x = ME, y = alpha_geomean)) +
  geom_jitter(aes(colour = ME), size = 3, alpha = .05) +
  geom_violin(aes(colour = ME)) +
  geom_boxplot(aes(colour = ME), width = .1) +
  scale_colour_manual(values = group_cols, name = "CpG set") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylab("Pr(hv) (geometric mean)")

pcontrast <- ggplot(contrasts, aes(x = ME, y = estimate, colour = ME_name)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_colour_manual(values = group_cols, name = "CpG set") +
  coord_flip() +
  labs(y = "Difference in Pr(hv) (geometric mean) vs mQTLcontrols", x = "",
       title = "Comparison of previous CPG sets groups to mQTLcontrols",
       subtitle = "lm with multiple comparison correction (Sidak)") +
  theme_minimal() +
  theme(legend.position = "none")

## pdecay - legend inside the plot, bottom-left corner
pdecay <- plot_decay_curve(MEsetdt, "Decay curve of Pr(HV) per threshold") +
  scale_colour_manual(values = group_cols, name = "CpG set") +
  theme(legend.position = "inside",
        legend.position.inside = c(0, 0),
        legend.justification = c(0, 0),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3))

# ── Save key objects for S07 ──────────────────────────────────────────────────
saveRDS(MEsetdt,            here("gitignore/MEsetdt.rds"))
saveRDS(geomMeanGR,         here("gitignore/geomMeanGR.rds"))

####################################################################################
## Test enrichement of the most likely germ layer-universal hvCpG in previous MEs ##
####################################################################################
if (!exists("listGR")){
  listGR <- list(top90 = makeGRfromMyCpGPos(vec = top99SNPrm, setname = "top99SNPrm"),
                 allButTop90 = makeGRfromMyCpGPos(
                   setdiff(totalSiteswGeomMean, top99SNPrm), "allButTop99"))}

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

plot_top99CpGsEnrichME <- ggplot(res_plot2, aes(x = quadrant_ord, y = log2OR, fill = signif)) +
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

################################################################################
## Load SIV plots calculated in fetalSIV folder script (in ing-p5)            ##
################################################################################
plots <- readRDS(here("gitignore/intercorrelationSIVfetal_sepSIV.rds"))

pinterlayer_corr <- ggplot(plots$interlayer_corr, aes(x=group, y=interlayer_r, group = group, fill = group))+
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_manual(values = group_cols) +
  theme_minimal(base_size = 14) +
  labs(y = "Mean inter-germ layer correlation\n(Pearson's r)")+
  theme(axis.title.x = element_blank(), legend.position = "none") 

pinterindividual_var <- ggplot(plots$CpG_summary, aes(x = interindividual_var, color = group)) +
  geom_density(alpha = 0.5)+
  scale_colour_manual(values = group_cols) +
  theme_minimal(base_size = 14) +
  labs(x = "Interindividual variation")

pbinned <- ggplot(plots$binned_summary_boot,
                  aes(x = bin, y = median_r, color = group, fill = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(
    aes(ymin = low, ymax = high),
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Interindividual variation",
    y = "Inter-germ layer correlation \n(median ± bootstrap CI)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

upperRow <- plot_grid(
  plot_grid(pMEalpha, pcontrast, nrow = 2, labels = c("A", "B")),
  pdecay,
  ncol = 2,
  rel_widths = c(1, 1),
  labels = c("", "C")
)

SIV_plot <- plot_grid(
  pinterlayer_corr + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  plot_grid(pinterindividual_var, pbinned, ncol = 1, align = "v",labels = c("D", "E")),
  ncol = 2,
  rel_widths = c(1, 1), labels = c("C", "")
)

final_plot <- plot_grid(
  upperRow,
  SIV_plot,
  ncol = 1,
  rel_heights = c(1, 1.2)
)

pdf(here("B_MultiTissues/dataOut/figures/script03/CompareWithpreviousMEs.pdf"),
    width = 18, height = 14)
print(final_plot)
dev.off()

### TBC here

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

# Convert gt to a patchwork-compatible element
table_element <- patchwork::wrap_elements(
  full = grid::grid.draw(gt::gt_as_gtable(summary_table_gt))
)

### TBC


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
length(top90SNPrm) # 385370

## Map on arrays
matches <- match(x = top90SNPrm, table = dico$chrpos_hg38)

Pos <- dico[na.omit(matches), ]

table(Pos$array)
# 450k 450k and EPIC          EPIC 
# 275          3608          2538 

#################################
## Test enrichment of features ##
#################################

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

<<<<<<< HEAD
# Method. ClusterProfiler
=======
  <<<<<<< HEAD
# Method. ClusterProfiler
=======
  <<<<<<< HEAD
# Method 1. ClusterProfiler
=======
  # Method. ClusterProfiler
  >>>>>>> f4378a16865649083ff37e95e78135fd1036eeb7
>>>>>>> 6ff651eada8a7269cd0538e27365b7ce01a915a4
>>>>>>> 44923781579b28d6862056d849b5e5c1f3e87b32

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

<<<<<<< HEAD
write.csv(df_sig, file = here("B_MultiTissues/dataOut/df_sig_GOtop90SNPrm.csv"),
          quote = F, row.names = F)

=======
  <<<<<<< HEAD
write.csv(df_sig, file = here("B_MultiTissues/dataOut/df_sig_GOtop90SNPrm.csv"),
          quote = F, row.names = F)

=======
  <<<<<<< HEAD
=======
  write.csv(df_sig, file = here("B_MultiTissues/dataOut/df_sig_GOtop90SNPrm.csv"),
            quote = F, row.names = F)

>>>>>>> f4378a16865649083ff37e95e78135fd1036eeb7
>>>>>>> 6ff651eada8a7269cd0538e27365b7ce01a915a4
>>>>>>> 44923781579b28d6862056d849b5e5c1f3e87b32
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
<<<<<<< HEAD
=======
  <<<<<<< HEAD
=======
  <<<<<<< HEAD
write.csv(df_sig, file = here("B_MultiTissues/dataOut/df_sig_GOtop90SNPrm.csv"),
          quote = F, row.names = F)






# ###########################################
# ## Test different p0 and p1 in raw alpha ##
# ###########################################
# 
# ## Sub test: just one small batch
# # R0 <- readRDS(here("B_MultiTissues/resultsDir_gitIgnored/Atlas/atlas_general/Atlas_batch001/results_atlas_general_250000CpGs_0_8p0_0_65p1.rds"))
# # R1 <- readRDS(here("B_MultiTissues/resultsDir_gitIgnored/Atlas/atlas_general/Atlas_batch001/results_atlas_general_250000CpGs_0_55p0_0_65p1.rds"))
# # R2 <- readRDS(here("B_MultiTissues/resultsDir_gitIgnored/Atlas/atlas_general/Atlas_batch001/results_atlas_general_250000CpGs_0_9p0_0_65p1.rds"))
# R3 <- readRDS(here("B_MultiTissues/resultsDir_gitIgnored/Atlas/atlas_general/Atlas_batch001/results_atlas_general_250000CpGs_0_8p0_0_9p1.rds"))
# 
# # R0 <- data.table(alpha = as.vector(R0), param = "0.8p0_0.65p1", chr_pos = rownames(R0))
# # R1 <- data.table(alpha = as.vector(R1), param = "0.55p0_0.65p1", chr_pos = rownames(R1))
# # R2 <- data.table(alpha = as.vector(R2), param = "0.9p0_0.65p1", chr_pos = rownames(R2))
# R3 <- data.table(alpha = as.vector(R3), param = "0.8p0_0.9p1", chr_pos = rownames(R3))
# 
# # setkey(R0, chr_pos)
# # setkey(R1, chr_pos)
# # setkey(R2, chr_pos)
# setkey(R3, chr_pos)

# # Merge all into one wide table
# merged_dt <- R0[R1, .(chr_pos, alpha_R0 = alpha, alpha_R1 = i.alpha),
#                 on = "chr_pos", nomatch = NULL]
# merged_dt <- merged_dt[R2, .(chr_pos, alpha_R0, alpha_R1, alpha_R2 = i.alpha),
#                        on = "chr_pos", nomatch = NULL]
# merged_dt <- merged_dt[R3, .(chr_pos, alpha_R0, alpha_R1, alpha_R2, alpha_R3 = i.alpha),
#                        on = "chr_pos", nomatch = NULL]
# 
# # Reshape to long: one row per CpG per comparison vs R0
# plot_long <- merged_dt %>%
#   pivot_longer(cols = c(alpha_R1, alpha_R2, alpha_R3),
#                names_to = "param", values_to = "alpha_other") %>%
#   mutate(param = factor(param,
#                         levels = c("alpha_R1", "alpha_R2", "alpha_R3"),
#                         labels = c("p0=0.55, p1=0.65",
#                                    "p0=0.90, p1=0.65",
#                                    "p0=0.80, p1=0.90")))
# 
# # Plot: 3 panels, x = R0 baseline, y = alternative
# ggplot(plot_long, aes(x = alpha_R0, y = alpha_other)) +
#   geom_point(alpha = 0.15, size = 0.4) +
#   geom_abline(slope = 1, intercept = 0,
#               colour = "red", linetype = "dashed", linewidth = 0.5) +
#   facet_wrap(~ param, ncol = 3) +
#   labs(x     = "Pr(hv) baseline (p0=0.80, p1=0.65)",
#        y     = "Pr(hv) alternative parameters",
#        title = "Sensitivity to parameter choice (batch 001, 250 CpGs)") +
#   theme_minimal(base_size = 11) +
#   theme(strip.text = element_text(face = "bold"))

