#################################################
## Plot results of algorithm ran on atlas data ##
#################################################
library(here)

source(here("05_hvCpGalgorithm", "quiet_library.R"))

## Load array results
resArray <- readRDS(here("05_hvCpGalgorithm/dataOut/resArray.RDS"))

## Add previous MEs including Maria's results
source(here("05_hvCpGalgorithm/runAlgo_myDatasets/exploreResults/prepPreviousSIV.R"))

## This code does:
### I. Histogram of coverage across datasets
### II. Load data & Manhattan plot
### III. Test enrichment of features for high alpha
### IV. Test for enrichment in other putative MEs for hvCpGs with alpha > threshold

## Data in WGBS atlas:
## from the CS cluster: sample_groups <- h5read("/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/all_matrix_noscale.h5","sample_groups")
sample_groups <- readRDS(here("05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/sample_groups.RDS"))

ggplot(data.frame(table(sample_groups)), aes(x = Freq)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of number of samples per dataset",
    x = "Number of samples",
    y = "Count of datasets"
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 1))

table(sample_groups)

##############################################
## I. Histogram of coverage across datasets ##
##############################################

t5 <- read.table(here("04_prepAtlas/CpG_coverage_freqtable5X.tsv"), header = T)
t10 <- read.table(here("04_prepAtlas/CpG_coverage_freqtable10X.tsv"), header = T)

# Add coverage type label
t5$coverage <- "≥5"
t10$coverage <- "≥10"

# Combine into one data.table
t_combined <- rbind(t5, t10)

# Filter out CpGs with zero dataset coverage if needed
t_combined <- t_combined[t_combined$datasets_covered_in > 0, ]

t_combined <- t_combined %>% group_by(coverage) %>%
  dplyr::arrange(-dplyr::row_number(datasets_covered_in)) %>%
  mutate(nCpGcum = cumsum(num_CpGs))

options(scipen=0)
# Plot
pdf(here("05_hvCpGalgorithm", "figures", "freqCpGperdataset.pdf"), width = 14, height = 4)
ggplot(t_combined, aes(x = as.factor(datasets_covered_in), y = nCpGcum, fill = coverage)) +
  geom_col(position = "dodge") +
  scale_y_continuous(
    breaks = seq(0, 100000000, by = 10000000),  # 10 million steps
    labels = label_number(scale = 1e-6, suffix = "M")
  ) +  scale_fill_manual(values = c("≥5" = "steelblue", "≥10" = "firebrick")) +
  labs(title = "Nbr of CpG covered across X or more datasets",
       x = ">= X datasets",
       y = "Number of CpGs with ≥3 samples covered",
       fill = "Coverage threshold") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(.2,.5),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA))
dev.off()

t_combined[t_combined$coverage %in% "≥10" &  t_combined$datasets_covered_in %in% 46,"num_CpGs"] /
  sum(t_combined[t_combined$coverage %in% "≥10","num_CpGs"])

## 84% of all CpGs (23/27.5M) are covered in 46 cell types
rm(t_combined, t5, t10)

####################################
## II. Load data & Manhattan plot ##
####################################

prepAtlasdt <- function(){
  # Define parent folder containing all "Atlas_batchXXX" folders
  parent_dir <- here("05_hvCpGalgorithm/resultsDir/Atlas/Atlas10X/")
  
  # Get list of relevant RData files
  rdata_files <- dir(parent_dir, pattern = "results_Atlas10X_[0-9]+CpGs_0_8p0_0_65p1\\.RData$", 
                     recursive = TRUE, full.names = TRUE)
  
  ## Check if all batches have ran
  length(rdata_files) == 231
  
  all_cpg_values <- numeric()
  pb <- progress_bar$new(total = length(rdata_files), format = "📦 :current/:total [:bar] :percent")
  
  for (file in rdata_files) {
    e <- new.env()
    load(file, envir = e)
    obj <- e[[ls(e)[1]]]
    if (is.matrix(obj)) obj <- obj[, 1]
    all_cpg_values <- c(all_cpg_values, obj)
    pb$tick()
  }
  
  # Create data.table from named vector
  dt <- data.table(
    name =   sub("-[0-9]+$", "", names(all_cpg_values)), # just keep the C position instead of C + 1
    alpha = as.numeric(all_cpg_values)
  )
  
  rm(e, pb, all_cpg_values, obj, file, parent_dir, rdata_files)
  
  #######################################################################
  # Parse "chr_pos" in name into chr, start_pos, end_pos. NB: takes a couple of minutes
  dt[, c("chr", "pos") := tstrsplit(name, "_", fixed = TRUE)]
  
  # Convert to integer/numeric if not already
  dt[, pos := as.integer(pos)]
  
  ## Check chromosomes present:
  message("Chromosomes in the dataset:")
  table(unique(dt$chr))
  
  # Convert chr from "chrN" to  factor
  dt[, chr := sub("chr", "", chr)]
  dt[, chr := factor(chr, levels = as.character(c(1:22, "X", "Y", "M")))]
  
  ## Check chromosomes order:
  message("Chromosomes in the dataset:")
  table(unique(dt$chr))
  
  ## Mark group membership in dt
  dt[, group := NA_character_]
  dt[name %in% DerakhshanhvCpGs_hg38, group := "hvCpG_Derakhshan"]
  dt[name %in% mQTLcontrols_hg38, group := "mQTLcontrols"]
  
  # Compute cumulative position offsets for Manhattan plot
  setorder(dt, chr, pos)
  
  offsets <- dt[, .(max_pos = max(pos, na.rm = TRUE)), by = chr]
  offsets[, cum_offset := c(0, head(cumsum(as.numeric(max_pos)), -1))]
  
  dt <- merge(dt, offsets[, .(chr, cum_offset)], by = "chr", all.x = TRUE, sort = FALSE)
  
  # Convert to integer/numeric if not already
  dt[, cum_offset := as.numeric(cum_offset)]
  dt[, pos2 := pos + cum_offset]
  
  return(dt)
}

system.time(Atlas_dt <- prepAtlasdt())

if (exists("doIprepAtlas") && isTRUE(doIprepAtlas)) {
  stop("stop here to only prepare atlas_dt")
}

# Compute chromosome centers for x-axis labeling
df2 <- Atlas_dt[, .(center = mean(range(pos2, na.rm = TRUE))), by = chr]
df2 <- merge(data.frame(chr = factor(c(1:22, "X", "Y", "M"), levels=as.character(c(1:22, "X", "Y", "M")))),
             df2, by = "chr", all.x = TRUE, sort = TRUE)
df2 <- na.omit(df2)

# Compute chromosome boundaries
df_bounds <- Atlas_dt[, .(min_pos = min(pos2, na.rm = TRUE), 
                          max_pos = max(pos2, na.rm = TRUE)), by = chr]

# Midpoints between chromosomes = where to draw dotted lines
df_bounds[, next_start := data.table::shift(min_pos, n = 1, type = "lead")]
vlines <- df_bounds[!is.na(next_start), .(xintercept = (max_pos + next_start)/2)]

plot <- ggplot() +
  # background cloud
  geom_point_rast(data = Atlas_dt[is.na(group)], 
                  aes(x = pos2, y = alpha),
                  color = "black", size = 0.01, alpha = 0.01, raster.dpi = 72) +
  # hvCpG highlights
  geom_point(data = Atlas_dt[group == "hvCpG_Derakhshan"],
             aes(x = pos2, y = alpha),
             color = "#DC3220", size = 1, alpha = 0.7) +
  # mQTL controls highlights
  geom_point(data = Atlas_dt[group == "mQTLcontrols"],
             aes(x = pos2, y = alpha),
             color = "#005AB5", size = 1, alpha = 0.7) +
  # Add dotted separators
  geom_vline(data = vlines, aes(xintercept = xintercept),
             linetype = 3, color = "grey60") +
  theme_classic() + theme(legend.position = "none") +
  scale_x_continuous(breaks = df2$center, labels = as.character(df2$chr), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

# Save as PDF — rasterization improves performance and file size
CairoPDF(here("05_hvCpGalgorithm/figures/ManhattanAlphaPlot_previoushvCpGplotted_atlas.pdf"), width = 15, height = 3)
print(plot)
dev.off()

## Without the layer of previous hvCpG and controls plotted:
plot <- ggplot() +
  geom_point_rast(data = Atlas_dt, 
                  aes(x = pos2, y = alpha),
                  color = "black", size = 0.01, alpha = 0.01, raster.dpi = 72) +
  # Add dotted separators
  geom_vline(data = vlines, aes(xintercept = xintercept),
             linetype = 3, color = "grey60") +
  theme_classic() + theme(legend.position = "none") +
  scale_x_continuous(breaks = df2$center, labels = as.character(df2$chr), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

# Save as PDF — rasterization improves performance and file size
CairoPDF(here("05_hvCpGalgorithm/figures/ManhattanAlphaPlot_atlas.pdf"), width = 15, height = 3)
print(plot)
dev.off()

###################################################################
## Calculate proba hvCpG minus matching control: is it always +? ##

data <- read.table(here("03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt"), header = T)

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

pdf(here("05_hvCpGalgorithm/figures/DifferenceOfProbabilityForhvCpG-matching_controlInAtlas.pdf"),
    width = 4, height = 5)
ggplot(merged, aes(x="diff", y=diffAlpha))+
  geom_jitter(data=merged[merged$diffAlpha>=0,], col="black", alpha=.5)+
  geom_jitter(data=merged[merged$diffAlpha<0,], fill="yellow",col="black",pch=21, alpha=.5)+
  geom_violin(width=.5, fill = "grey", alpha=.8) +
  geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
  theme_minimal(base_size = 14)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), title = element_text(size=10))+
  ggtitle("P(hvCpG) minus P(matching control) in atlas")+
  ylab("Difference of probability")
dev.off()

##################################
## Save names high alpha points ##
##################################
table(Atlas_dt$alpha >= 0.7)
# FALSE     TRUE 
# 22335423   700603 

## Map on arrays
matches <- match(
  x = unlist(Atlas_dt[Atlas_dt$alpha >= 0.7, "name"]),
  table = dico$chrpos_hg38
)

highAlphaPos <- dico[na.omit(matches), ]

table(highAlphaPos$array)
# 450k    450k and EPIC      EPIC 
# 340          4483          3431 

# 4483+340 = 4823 on the 450k array

saveRDS(highAlphaPos, here("05_hvCpGalgorithm/runAlgo_myDatasets/exploreResults/fetalSIV/highAlphaPos_atlas0.7.RDS"))

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
## Load the results at https://genome.ucsc.edu/cgi-bin/hgTracks

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
  scale_x_continuous(breaks = df2[df2$chr == "M","center"],
                     labels = as.character(df2[df2$chr == "M","chr"]), expand = c(0, 0)) +
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
  scale_x_continuous(breaks = df2[df2$chr == "Y","center"],
                     labels = as.character(df2[df2$chr == "Y","chr"]), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

Atlas_dt[Atlas_dt$chr == "Y" & Atlas_dt$alpha > 0.7,]
# High alpha in 3 regions: chrY:5,043,848-6,534,238, chrY:10,107,290-11,747,410, chrY:56,822,399-56,841,336

#####################################################
## Find regions of high alpha using sliding window ##
#####################################################
# 
# high_alpha_thresh <- 0.7
# Atlas_dt[, high_alpha := alpha > high_alpha_thresh]
# 
# # Function to compute high-alpha fraction in each window
# find_high_alpha_windows <- function(dt, window_size, step_size) {
#   res <- list()
#   for (ch in unique(dt$chr)) {
#     chr_dt <- dt[chr == ch]
#     starts <- seq(min(chr_dt$start_pos), max(chr_dt$end_pos), by = step_size)
#     for (s in starts) {
#       e <- s + window_size - 1
#       w <- chr_dt[start_pos >= s & end_pos <= e]
#       if (nrow(w) > 0) {
#         frac <- sum(w$high_alpha) / nrow(w)
#         res[[length(res)+1]] <- data.table(chr = ch, window_start = s, window_end = e, frac_high_alpha = frac)
#       }
#     }
#   }
#   rbindlist(res)
# }
# 
# # system.time(high_alpha_windows <- find_high_alpha_windows(Atlas_dt, window_size = 100000, step_size = 50000))
## 100000, 50000: 10 minutes to run
# saveRDS(object = high_alpha_windows, 
#         file = here("05_hvCpGalgorithm/runAlgo_myDatasets/exploreResults/high_alpha_windows.RDS"))
# 
# # Keep windows where >50% of CpGs have high alpha
# hotspots <- high_alpha_windows[frac_high_alpha > 0.5]
# hotspots
# chr window_start window_end frac_high_alpha
# <char>        <num>      <num>           <num>
#   1:      1       517452     617451       1.0000000
# 2:      3     92911806   93011805       0.6666667
# 3:      6     32496889   32596888       0.5670498
# 4:      6     32546889   32646888       0.5853659
# 5:      6     32596889   32696888       0.5498721
# 6:      9     61210596   61310595       0.6666667
# 7:     15     21850602   21950601       0.6000000
# 8:     16     37260230   37360229       1.0000000
# 9:     16     37310230   37410229       1.0000000


# subAtlasdt <- Atlas_dt[Atlas_dt$chr ==3 & 
#                          Atlas_dt$start_pos >= 92911806 & 
#                          Atlas_dt$end_pos <= 93011805,]

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
bed_features <- genomation::readTranscriptFeatures(
  "~/Documents/Project_hvCpG/hg38_GENCODE_V47.bed")

# Annotate CpGs and see which regions have higher alpha
anno_result <- genomation::annotateWithGeneParts(
  target = gr_cpg, feature = bed_features)

anno_result@perc.of.OlapFeat
# promoter     exon   intron 
# 96.61240 74.37305 89.72931 

## Add info from annotation to our GRange object
gr_cpg$featureType <- ifelse(anno_result@members[, "prom"] == 1, "promoter",
                             ifelse(anno_result@members[, "exon"] == 1, "exon",
                                    ifelse(anno_result@members[, "intron"] == 1, "intron", "intergenic")))

############################
## find exon-intron junction
txdb <- txdbmaker::makeTxDbFromGFF("~/Documents/Project_hvCpG/gencode.v38.annotation.gtf", format = "gtf")

exons <- GenomicFeatures::exons(txdb)
introns <- GenomicFeatures::intronsByTranscript(txdb)

# Flatten intron list
introns_flat <- unlist(introns, use.names = FALSE)

# Define a small window (±25 bp) around exon–intron boundaries
junction_window <- 25

# Create exon-intron junction regions (exon ends + intron starts)
exon_ends <- resize(exons, width = 1, fix = "end")
intron_starts <- resize(introns_flat, width = 1, fix = "start")

# Combine and expand around the junction
junctions <- suppressWarnings(reduce(c(
  flank(exon_ends, width = junction_window, both = TRUE),
  flank(intron_starts, width = junction_window, both = TRUE)
)))

# Identify CpGs overlapping exon–intron junctions
is_junction <- countOverlaps(gr_cpg, junctions) > 0
region_type[is_junction] <- "ex-intr junction"

# add annotation to GRanges object
gr_cpg$preciseFeatureType <- region_type

## check that the only diff is ex-intr junction
table(gr_cpg[gr_cpg$featureType != gr_cpg$preciseFeatureType,]$preciseFeatureType) 

#########################################################################
## Find the first promoter and the first exon by absolute distance to TSS
gr_cpg$dist2TSS <- anno_result@dist.to.TSS$dist.to.feature
gr_cpg$TSSname <- anno_result@dist.to.TSS$feature.name

mcols(gr_cpg) <- mcols(gr_cpg) %>% as.data.frame() %>% 
  dplyr::group_by(TSSname, featureType) %>% 
  dplyr::mutate(pos = min_rank(abs(dist2TSS))) %>% data.frame()

mcols(gr_cpg) <- mcols(gr_cpg) %>% as.data.frame() %>% 
  mutate(
    preciseFeatureType = case_when(
      featureType == "promoter" & pos == 1 ~ "first promoter",
      featureType == "exon" & pos == 1 ~ "first exon",
      TRUE ~ as.character(preciseFeatureType)  # keep existing otherwise
    )
  )

mcols(gr_cpg) %>% as.data.frame() %>%
  dplyr::group_by(preciseFeatureType) %>%
  dplyr::summarise(meanAlpha = mean(alpha),
                   medianAlpha = median(alpha))
# preciseFeatureType meanAlpha medianAlpha
# <chr>                  <dbl>       <dbl>
# 1 ex-intr junction       0.128      0.0476 ********* The lowest! Very conserved
# 2 exon                   0.148      0.0703
# 3 first exon             0.169      0.0994
# 4 first promoter         0.163      0.0784
# 5 intergenic             0.196      0.121 ********* The higghest
# 6 intron                 0.158      0.0755
# 7 promoter               0.147      0.0591

mcols(gr_cpg) %>% as.data.frame() %>%
  dplyr::group_by(featureType) %>%
  dplyr::summarise(meanAlpha = mean(alpha),
                   medianAlpha = median(alpha))
# featureType meanAlpha medianAlpha
# <chr>           <dbl>       <dbl>
# 1 exon            0.148      0.0710
# 2 intergenic      0.196      0.121 ********* The higghest
# 3 intron          0.158      0.0755
# 4 promoter        0.147      0.0587 ********* The lowest

# visualize methylation levels by region
pdf(here("05_hvCpGalgorithm/figures/barplotFeaturesbyAlpha.pdf"), width = 6, height = 4)
ggplot(mcols(gr_cpg), aes(x = preciseFeatureType, y = alpha, fill = preciseFeatureType)) +
  geom_violin()+
  geom_boxplot(outlier.size = 0.5, alpha = 0.8, width = .3) +
  theme_minimal(base_size = 14) +
  labs(y = "p(hv)") + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none", axis.title.x = element_blank())
dev.off()

# saveRDS(Atlas_dt, "~/Documents/Project_hvCpG/Atlas_dt-3nov25.RDS")


# Promoters/enhancers: FANTOM5 database: https://fantom.gsc.riken.jp/5/
#   
#   I think it would also be interesting to look at chromatin states using ChromHMM  to get a sense of regions that are transcriptionally active etc (as we’ve done in several of our papers):
#   ChromHMM: https://compbio.mit.edu/ChromHMM/
#   R package to work with ChromHMM: https://www.bioconductor.org/packages/release/bioc/vignettes/segmenter/inst/doc/segmenter.html
# 
# Tissue-specificity is a major complicating factor of course. Even though we’re dealing with MEs, functional relevance could vary according to cell type, as we’ve found with PAX8, POMC and LY6S!



######################################################################################
## NB: to save time, we can start directly here during development of the analysis! ##
# Atlas_dt <- readRDS("~/Documents/Project_hvCpG/Atlas_dt-3nov25.RDS")


mcols(gr_cpg) %>% as.data.frame() %>%
  dplyr::group_by(preciseFeatureType) %>%
  dplyr::summarise(meanAlpha = mean(alpha))

## Kruskal–Wallis test: are the groups different in alpha values?
kruskal.test(alpha ~ preciseFeatureType, data = mcols(gr_cpg))
# Kruskal-Wallis chi-squared = 256663, df = 6, p-value < 2.2e-16

## Post-hoc pairwise comparison
pairwise_results <- pairwise.wilcox.test(
  mcols(gr_cpg)$alpha,
  mcols(gr_cpg)$preciseFeatureType,
  p.adjust.method = "fdr"
)
pairwise_results
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
#                 ex-intr junction  exon    first exon first promoter intergenic intron 
#   exon           < 2e-16          -       -          -              -          -      
#   first exon     < 2e-16          < 2e-16 -          -              -          -      
#   first promoter < 2e-16          3.5e-05 < 2e-16    -              -          -      
#   intergenic     < 2e-16          < 2e-16 < 2e-16    < 2e-16        -          -      
#   intron         < 2e-16          < 2e-16 < 2e-16    < 2e-16        < 2e-16    -      
#   promoter       < 2e-16          < 2e-16 < 2e-16    < 2e-16        < 2e-16    < 2e-16
# 
# P value adjustment method: fdr 

#################################################################################
## Test for enrichment in other putative MEs for hvCpGs with alpha > threshold ##
#################################################################################

# Parse with regex all the cpg tested
parsed <- str_match(Atlas_dt$name, "(chr[0-9XYM]+)_(\\d+)")

# Build GRanges
allcpg_GR <- GRanges(
  seqnames = parsed[,2],
  ranges   = IRanges(start = as.numeric(parsed[,3]),
                     end   = as.numeric(parsed[,3]))
); rm(parsed)

###########################################
## meQTL (vmeQTL) identified in MZ twins by Jordana Bell 
vmeQTL_hg19probes <- readxl::read_xlsx(here("05_hvCpGalgorithm/dataPrev/vmeQTL_vCpG_359pair_sig_Zhang2025.xlsx"))
vmeQTL_hg38 <- na.omit(dico$chrpos_hg38[match(vmeQTL_hg19probes$vCpG, dico$CpG)]) ; rm(vmeQTL_hg19probes)

# a vector of 1773 SIV from Harris 2012 (HarrisSIV_hg38)
HarrisSIV_hg38

# one of 1579 ESS from Van Baak 2018 (VanBaakESS_hg38)
VanBaakESS_hg38

# a GRange object for Kessler 2018 676 SIV regions (KesslerSIV_GRanges_hg38)
overlaps <- findOverlaps(query = KesslerSIV_GRanges_hg38, subject = allcpg_GR)
# Extract the overlapping ranges
CpG_overlapping     <- allcpg_GR[subjectHits(overlaps)]
KesslerSIV_hg38 <- na.omit(paste0(CpG_overlapping@seqnames, "_", CpG_overlapping@ranges))
length(KesslerSIV_hg38) # 2700

# a GRange object for Gunasekara 2019 9926 corSIV regions (corSIV_GRanges_hg38)
overlaps <- findOverlaps(query = corSIV_GRanges_hg38, subject = allcpg_GR)
# Extract the overlapping ranges
CpG_overlapping     <- allcpg_GR[subjectHits(overlaps)]
corSIV_hg38 <- na.omit(paste0(CpG_overlapping@seqnames, "_", CpG_overlapping@ranges))
length(corSIV_hg38) # 70222

# a vector of 3644 hvCpG from Derakhshan 2022 (DerakhshanhvCpGs_hg38)
DerakhshanhvCpGs_hg38

# a vector for matching mQTL controls (mQTLcontrols_hg38)
mQTLcontrols_hg38

# Silver2022_SoCCpGs_10WGBS 
SoCCpGs_hg38

#################################
## Check overlaps with upset plot
sets <- list(
  vmeQTL = vmeQTL_hg38,
  HarrisSIV = HarrisSIV_hg38,
  VanBaakESS = VanBaakESS_hg38,
  KesslerSIV = KesslerSIV_hg38,
  CoRSIV = corSIV_hg38,
  SoCCpGs = SoCCpGs_hg38,
  hvCpG = DerakhshanhvCpGs_hg38,
  mQTLcontrols = mQTLcontrols_hg38
)

# Create the plot in a base graphics device
pdf(NULL)  # draw to null device to avoid displaying
upset(fromList(sets), nsets = 7, order.by = "freq")
grid_plot <- grid.grab()  # Capture as a grid object
dev.off()

# Now save the captured grid object to a real PDF
pdf(here("05_hvCpGalgorithm/figures/upsetPreviousME.pdf"), width = 12, height = 5)
grid.draw(grid_plot)
dev.off()

#######################################################
## Check alpha for the different MEs: are they high? ##

MEsetdt <- rbindlist(lapply(names(sets), function(nm) {
  Atlas_dt[.(sets[[nm]]), on = .(name), .(name, alpha)][, ME := nm]
}))

MEsetdt <- na.omit(MEsetdt) ## 33478 so far (2 sept)

## Control as baseline
MEsetdt[, ME := relevel(factor(ME), ref = "mQTLcontrols")]

p1 <- ggplot(MEsetdt, aes(x = ME, y = alpha)) +
  geom_jitter(data = MEsetdt,
              aes(fill=ME), pch=21, size = 3, alpha = .1)+
  geom_violin(aes(col=ME))+
  geom_boxplot(aes(col=ME), width = .1) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylab("Probability of being a hvCpG")

p1

## Statistical comparisons of alpha between MEs

# 1️⃣ Fit the model with mQTLcontrols as baseline
fit <- lm(alpha ~ ME, data = MEsetdt)

# 2️⃣ Get estimated marginal means and contrasts vs baseline
emm <- emmeans(fit, ~ ME)
contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = "mQTLcontrols", adjust = "sidak") %>%
  as.data.frame()

# 3️⃣ Prepare for plotting
contrasts <- contrasts %>%
  mutate(ME = contrast,  # rename for clarity
         lower = estimate - 1.96*SE,
         upper = estimate + 1.96*SE)

# 4️⃣ Plot
p2 <- ggplot(contrasts, aes(x = ME, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    y = "Difference in alpha vs mQTLcontrols",
    x = "ME group",
    title = "Comparison of ME groups to mQTLcontrols",
    subtitle = "lm with multiple comparison correction (Sidak)"
  ) +
  theme_minimal()

pdf(here("05_hvCpGalgorithm/figures/alphaComparisonBetweenMEtypes.pdf"), width = 14, height = 4)
cowplot::plot_grid(p1,p2, rel_widths = c(1, .8))
dev.off()

###############
## FIND PAX8
# Chromosome 2, NC_000002.12 (113215997..113278921,

# PAX8 as a GRanges
gr_PAX8 <- c(bed_features$promoters[grep("ENST00000429538.8", bed_features$promoters$name)],
  bed_features$exons[grep("ENST00000429538.8", bed_features$exons$name)],
  bed_features$introns[grep("ENST00000429538.8", bed_features$introns$name)],
  bed_features$TSSes[grep("ENST00000429538.8", bed_features$TSSes$name)])

# Find overlaps between CpGs and the gene region
hits <- findOverlaps(gr_cpg, gr_PAX8)

# Extract matching rows from original df
df_hits <- gr_cpg[queryHits(hits), ]

# Add all annotations
df_hits <- Atlas_dt[match(paste0(df_hits@seqnames, "_", df_hits@ranges), Atlas_dt$name)]

# Determine limits and breaks
x_min <- floor(min(df_hits$pos) / 5000) * 5000
x_max <- ceiling(max(df_hits$pos) / 5000) * 5000
breaks_seq <- seq(x_min, x_max, by = 5000)

ggplot(df_hits, aes(x = pos, y = alpha)) +
  geom_smooth(col = "black") +
  geom_point(aes(fill = region_type), pch = 21) +
  theme_minimal(base_size = 14) +
  ylab("p(hv)") +
  scale_x_continuous(
    breaks = breaks_seq,
    labels = function(x) paste0(formatC(x / 1000, format = "f", digits = 0), "k")
  ) +
  xlab("Genomic position (chr2)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
