#################################################
## Plot results of algorithm ran on atlas data ##
#################################################
## Previously:
# 9% CpGs sequenced in the 3 germ layers are SNPs (1559985 out of 18195276)
# 52% CpGs sequenced in the 3 germ layers with pr(hv) >=90% are SNPs (94627 out of 182710)
## --> Variability could simply come from SNPs!! We need to exclude them.
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
  resArray <- readRDS(here("B_MultiTissues/dataOut/resArray.RDS"))
}

##################################
## Save all data in RDS objects ##
##################################
savePrepedAtlasFile <- function(file, p0, p1) {
  out_path <- here(paste0("gitignore/resultsAtlasPrepared/fullres_", 
                          p0, "p0_", p1, "p1_", file, ".rds"))
  
  if (file.exists(out_path)) {
    message("File ", file, " already prepared - skipping.")
    return(invisible(NULL))
  }
  
  # Check files exist before attempting
  parent_dir <- here(paste0("B_MultiTissues/resultsDir_gitIgnored/Atlas/", file))
  n_files <- length(base::dir(parent_dir, 
                              pattern = paste0(p0, "p0_", p1, "p1.rds$"),
                              recursive = TRUE))
  if (n_files == 0) {
    message("Skipping ", file, " - no matching files found (incomplete run?).")
    return(invisible(NULL))
  }
  
  system.time(Atlas_dt <- prepAtlasdt(file, p0, p1))
  saveRDS(Atlas_dt, file = out_path)
  message("Saved: ", out_path)
}

for (subdir in list.files(here("B_MultiTissues/resultsDir_gitIgnored/Atlas/"))) {
  savePrepedAtlasFile(file = subdir, p0 = "0_8", p1 = "0_65")
}
## NB: the non atlas_general are INCOMPLETE --> rm and rerun when all finished!!

## Different p0 and p1 tested
for (subdir in list.files(here("B_MultiTissues/resultsDir_gitIgnored/Atlas/"))) {
  savePrepedAtlasFile(file = subdir, p0 = "0_8", p1 = "0_9")
}

###########################################
## Test different p0 and p1 in raw alpha ##
###########################################

Atlas_dt <- readRDS(
  "/home/alice/Documents/GIT/2024_hvCpG/gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_atlas_general.rds")
nrow(Atlas_dt) # 21.522.541

if (exists("doIprepAtlas") && isTRUE(doIprepAtlas)) {
  stop("stop here to only prepare atlas_dt")
}

Atlas_dt_80p090p1 <- readRDS(
  "/home/alice/Documents/GIT/2024_hvCpG/gitignore/resultsAtlasPrepared/fullres_0_8p0_0_9p1_atlas_general.rds")

# Set key if not already set
setkey(Atlas_dt, name)
setkey(Atlas_dt_80p090p1, name)

# Merge only columns needed for plotting
merged_dt <- Atlas_dt[Atlas_dt_80p090p1, 
                      .(name, 
                        alpha_general  = alpha,   
                        alpha_80p090p1 = i.alpha),
                      on = "name",
                      nomatch = NULL] # inner join

# Check
nrow(merged_dt)
head(merged_dt)

set.seed(1234)

merged_dt[sample(.N, 100000)] |>
  ggplot(aes(x = alpha_general, y = alpha_80p090p1)) +
  geom_point(alpha = 0.05, size = 0.3) +   # small/transparent for density
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  labs(x = "Pr(hv) p0=80%, p1=65%", y = "Pr(hv) p0=80%, p1=90%") +
  theme_minimal()

## Ranking is preserved. The tight linear band confirms both settings agree on which 
# CpGs are most variable --> hvCpG list is robust to this parameter choice

#######################
## Data in WGBS atlas:

## from the CS cluster: sample_groups <- h5read("/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/all_matrix_noscale.h5","sample_groups")
# sample_groups <- readRDS(here("05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/sample_groups.RDS"))
# 
# ggplot(data.frame(table(sample_groups)), aes(x = Freq)) +
#   geom_histogram(bins = 100, fill = "steelblue", color = "white") +
#   theme_minimal(base_size = 14) +
#   labs(
#     title = "Distribution of number of samples per dataset",
#     x = "Number of samples",
#     y = "Count of datasets"
#   ) +
#   scale_x_continuous(breaks = seq(0, 10, by = 1))
# 
# table(sample_groups)

SupTab1_Loyfer2023 <- read.csv(here("B_MultiTissues/dataIn/SupTab1_Loyfer2023.csv"))
SupTab1_Loyfer2023$group <- paste(SupTab1_Loyfer2023$Source.Tissue, SupTab1_Loyfer2023$Cell.type, sep = " - ")
table(table(SupTab1_Loyfer2023$group)[table(SupTab1_Loyfer2023$group) >=3])
# 3  4  5  6 10 
# 33  9  2  1  1 

###########################################
## Histogram of coverage across datasets ##
###########################################

# t5 <- read.table(here("04_prepAtlas/CpG_coverage_freqtable5X.tsv"), header = T)
# t10 <- read.table(here("04_prepAtlas/CpG_coverage_freqtable10X.tsv"), header = T)
# 
# # Add coverage type label
# t5$coverage <- "≥5"
# t10$coverage <- "≥10"
# 
# # Combine into one data.table
# t_combined <- rbind(t5, t10)
# 
# # Filter out CpGs with zero dataset coverage if needed
# t_combined <- t_combined[t_combined$datasets_covered_in > 0, ]
# 
# t_combined <- t_combined %>% group_by(coverage) %>%
#   dplyr::arrange(-dplyr::row_number(datasets_covered_in)) %>%
#   mutate(nCpGcum = cumsum(num_CpGs))
# 
# options(scipen=0)
# # Plot
# pdf(here("05_hvCpGalgorithm", "figures", "freqCpGperdataset.pdf"), width = 14, height = 4)
# ggplot(t_combined, aes(x = as.factor(datasets_covered_in), y = nCpGcum, fill = coverage)) +
#   geom_col(position = "dodge") +
#   scale_y_continuous(
#     breaks = seq(0, 100000000, by = 10000000),  # 10 million steps
#     labels = label_number(scale = 1e-6, suffix = "M")
#   ) +  scale_fill_manual(values = c("≥5" = "steelblue", "≥10" = "firebrick")) +
#   labs(title = "Nbr of CpG covered across X or more datasets",
#        x = ">= X datasets",
#        y = "Number of CpGs with ≥3 samples covered",
#        fill = "Coverage threshold") +
#   theme_minimal(base_size = 14) +
#   guides(fill = guide_legend(position = "inside")) +
#   theme(legend.position.inside = c(.2,.5),
#         legend.box = "horizontal",
#         legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
#         legend.key = element_rect(fill = "white", color = NA))
# dev.off()
# 
# t_combined[t_combined$coverage %in% "≥10" &  t_combined$datasets_covered_in %in% 46,"num_CpGs"] /
#   sum(t_combined[t_combined$coverage %in% "≥10","num_CpGs"])
# 
# ## 84% of all CpGs (23/27.5M) are covered in 46 cell types
# rm(t_combined, t5, t10)

####################
## Plot Manhattan ##
####################

plotManhattan1 <- plotManhattanFromdt(Atlas_dt, plotDerakhshan = FALSE)
ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/Manhattan/ManhattanAlphaPlot_atlas.png"),
  plot = plotManhattan1, width = 14, height = 4,
  dpi = 300, bg = "white")


## Only previous MEs: 

# 1. Convert GRanges to data.table for foverlaps
putativeME_dt <- as.data.table(putativeME_GR)[, .(
  chr   = sub("chr", "", seqnames),   # "chr1" -> "1" to match Atlas_dt
  start = start,
  end   = end,
  set = set
)]

# 2. Add a point range to Atlas_dt (foverlaps needs two position columns)
Atlas_dt[, pos_end := pos]   # point interval: start == end

# 3. Set keys for foverlaps
setkey(Atlas_dt,  chr, pos, pos_end)
setkey(putativeME_dt, chr, start, end)

# 4. Overlap join — returns only Atlas_dt rows falling inside a corSIV region
Atlas_putativeME <- foverlaps(Atlas_dt, putativeME_dt,
                              by.x = c("chr", "pos", "pos_end"),
                              by.y = c("chr", "start", "end"),
                              type = "within",
                              nomatch = NULL)   # NULL = only matched rows (like inner join)

Atlas_dt[, pos_end := NULL]

plotManhattan3 <- plotManhattanFromdt(Atlas_putativeME, colorBySet = TRUE,
                                      plotDerakhshan = FALSE, transp = .2) +
  theme(legend.position = "none")

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/Manhattan/ManhattanAlphaPlot_atlas_prevMEs.png"),
  plot = plotManhattan3, width = 14, height = 10,
  dpi = 300, bg = "white")

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
res <- Atlas_dt[!is.na(group)]

hv_alpha <- res[, c("name", "alpha")]
colnames(hv_alpha) <- c("hvCpG", "alpha_hvCpG")

# Merge control alphas
ctrl_alpha <- res[, c("name", "alpha")]
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
p <- ggplot(merged, aes(x="diff", y=diffAlpha))+
  geom_jitter(data=merged[merged$diffAlpha>=0,], col="black", alpha=.5)+
  geom_jitter(data=merged[merged$diffAlpha<0,], fill="yellow",col="black",pch=21, alpha=.5)+
  geom_violin(width=.5, fill = "grey", alpha=.8) +
  geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
  theme_minimal(base_size = 14)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), title = element_text(size=10))+
  ggtitle("P(hvCpG) minus P(matching control) in atlas")+
  ylab("Difference of probability")

##########################################
## What are the gaps in Manhattan plot? ##
# Compute the gap between consecutive CpGs on the same chromosome
Atlas_dt[, gap := pos - data.table::shift(pos), by = chr]

# Identify large gaps (>= 500k bp)
gaps_dt <- Atlas_dt[gap >= 500000, .(
  chr,
  gap_start = data.table::shift(pos),
  gap_end = pos,
  gap_size = gap
)]

# Drop first NA (since shift introduces one per chromosome)
gaps_dt[!is.na(gap_size)]

####################################
## Mitochondrial DNAm variability ##
####################################
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-023-09541-9?utm_source=chatgpt.com
## Near absence of 5mC, so expected

ggplot() +
  geom_point(data = Atlas_dt[Atlas_dt$chr == "M",], 
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
  geom_point(data = Atlas_dt[Atlas_dt$chr == "Y",], 
             aes(x = pos2, y = alpha),
             color = "black", size = 1, alpha = .5)+
  theme_classic() + theme(legend.position = "none") +
  geom_hline(yintercept = .7, linetype = 3)+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

# High alpha in 3 regions: chrY:5,043,848-6,534,238, chrY:10,107,290-11,747,410, chrY:56,822,399-56,841,336

#####################################################
## III. Test enrichment of features for high alpha ##
#####################################################

# Create GRanges of 23036026 CpGs
gr_cpg <- GRanges(
  seqnames = paste0("chr", Atlas_dt$chr),
  ranges = IRanges(start = Atlas_dt$pos, end = Atlas_dt$pos),
  alpha = Atlas_dt$alpha
)

# Import bed file
bed_features <- genomation::readTranscriptFeatures(here("gitignore/hg38_GENCODE_V47.bed"))

# Annotate CpGs and see which regions have higher alpha (takes long)
anno_result <- genomation::annotateWithGeneParts(
  target = gr_cpg, feature = bed_features)

anno_result@perc.of.OlapFeat
# promoter     exon   intron 
# 96.59648 73.64187 89.42734

## Add info from annotation to our GRange object
gr_cpg$featureType <- ifelse(anno_result@members[, "prom"] == 1, "promoter",
                             ifelse(anno_result@members[, "exon"] == 1, "exon",
                                    ifelse(anno_result@members[, "intron"] == 1, "intron", "intergenic")))

mcols(gr_cpg) %>% as.data.frame() %>%
  dplyr::group_by(featureType) %>%
  dplyr::summarise(meanAlpha = mean(alpha),
                   medianAlpha = median(alpha))
# featureType meanAlpha medianAlpha
# <chr>           <dbl>       <dbl>
#   1 exon            0.156      0.0805
# 2 intergenic      0.203      0.135 ********* The higghest
# 3 intron          0.162      0.0834
# 4 promoter        0.155      0.0661 ********* The lowest! Very conserved

## Kruskal–Wallis test: are the groups different in alpha values?
kruskal.test(alpha ~ featureType, data = mcols(gr_cpg))
# Kruskal-Wallis rank sum test
# 
# data:  alpha by featureType
# Kruskal-Wallis chi-squared = 248241, df = 3, p-value < 2.2e-16

## Post-hoc pairwise comparison
# pairwise_results <- pairwise.wilcox.test(
#   mcols(gr_cpg)$alpha,
#   mcols(gr_cpg)$featureType,
#   p.adjust.method = "fdr"
# )
# pairwise_results

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  mcols(gr_cpg)$alpha and mcols(gr_cpg)$featureType 
# 
# exon   intergenic intron
# intergenic <2e-16 -          -     
#   intron     <2e-16 <2e-16     -     
#   promoter   <2e-16 <2e-16     <2e-16
# P value adjustment method: fdr 

# visualize methylation levels by region
p <- ggplot(mcols(gr_cpg), aes(x = featureType, y = alpha, fill = featureType)) +
  geom_violin()+
  geom_boxplot(outlier.size = 0.5, alpha = 0.8, width = .3) +
  theme_minimal(base_size = 14) +
  labs(y = "p(hv)") + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none", axis.title.x = element_blank())

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/barplotFeaturesbyAlpha.png"),
  plot = p, width = 6, height = 4,
  dpi = 300, bg = "white")