#############################################
## Overlap plot: Atlas (x) vs Array (y)    ##
#############################################
library(here)
source(here("05_hvCpGalgorithm/quiet_library.R"))

source(here("05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/prephvCpGandControls.R"))
hvCpGandControls <- prephvCpGandControls(codeDir = "~/Documents/GIT/2024_hvCpG/")

################################ 
## --- Prepare Array data --- ##
################################ 

source(here("05_hvCpGalgorithm/runAlgo_myDatasets/exploreResults/S02_analyseResultsArray_local.R"))

################################ 
## --- Prepare Atlas data --- ##
################################

mergeAtlasRunBatches <- function(parent_dir = "05_hvCpGalgorithm/resultsDir/Atlas10X/", # parent folder containing all "Atlas_batchXXX" folders
                                 analysis = "Atlas10X", alphaname = "alpha_atlas"){
  # Get list of relevant RData files
  rdata_files = dir(here(parent_dir), pattern = paste0("^results_", analysis, ".*CpGs_0_8p0_0_65p1\\.RData$"),
                    recursive = TRUE, full.names = TRUE
  )
  
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
    name = names(all_cpg_values),
    alpha = as.numeric(all_cpg_values)
  )
  
  atlas_dt <- dt[, .(chrpos = name, tmp = alpha)]
  setnames(atlas_dt, "tmp", alphaname)
  return(atlas_dt)
}

atlas_dt <- mergeAtlasRunBatches(parent_dir = "05_hvCpGalgorithm/resultsDir/Atlas10X/", 
                                 analysis = "Atlas10X", alphaname = "alpha_atlas")

###################################### 
## --- Merge Array & Atlas data --- ##
######################################

res_Alpha_Atlas <- dplyr::full_join(resCompArray, atlas_dt)

######################################### 
## --- Test 1: matching pos error? --- ##
#########################################

## Positions to test
## high alpha in array, low in atlas
## high alpha in atlas, low in array
pos2test <- c(res_Alpha_Atlas[res_Alpha_Atlas$alpha_array_all > 0.9 &
                                    res_Alpha_Atlas$alpha_atlas < 0.1 &
                                    !is.na(res_Alpha_Atlas$alpha_atlas ),][1,"chrpos"],
              res_Alpha_Atlas[res_Alpha_Atlas$alpha_array_all <0.1 &
                                    res_Alpha_Atlas$alpha_atlas > 0.9 &
                                    !is.na(res_Alpha_Atlas$alpha_atlas ),][1,"chrpos"])

# Extract only the desired columns from the matrix
# Load the names of CpGs
datadir <-"~/Documents/Project_hvCpG/10X/all_matrix_noscale.h5"
cpg_names <- h5read(datadir, "cpg_names")

col_indices <- match(pos2test, cpg_names)

subset_matrix <- h5read(datadir, "matrix", index = list(NULL, col_indices))

# Check result
dim(subset_matrix)
apply(subset_matrix, 2, hist)

## seems correct

## Chr1 seems to work more, are the positions matching? Test chr1 and chr2
head(res_Alpha_Atlas)

compChr1 <- res_Alpha_Atlas[res_Alpha_Atlas$chr %in% "chr1",]
summary(lm(alpha_array_all~alpha_atlas, data = compChr1))
ggplot(compChr1, aes(x=alpha_array_all, y=alpha_atlas))+
  geom_point()+
  geom_smooth(method = "lm")

compChr2 <- res_Alpha_Atlas[res_Alpha_Atlas$chr %in% "chr2",]
summary(lm(alpha_array_all~alpha_atlas, data = compChr2))
ggplot(compChr2, aes(x=alpha_array_all, y=alpha_atlas))+
  geom_point()+
  geom_smooth(method = "lm")

## THERE MUST BE A CODE ERROR all positions must be mixed!

################################## 
## --- Test 2: low N error? --- ##
##################################

p1 <- ggplot(res_Alpha_Atlas, aes(x=alpha_array_all, y=alpha_array_3ind, fill = group, col = group)) +
  geom_point(data = res_Alpha_Atlas[is.na(res_Alpha_Atlas$group),], pch = 21, alpha = 0.05) +
  geom_point(data = res_Alpha_Atlas[!is.na(res_Alpha_Atlas$group),], pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.18,0.85),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(title = "Probability of being hypervariable",
       x = "P(hv) considering all array data",
       y = "P(hv) considering 3 individuals per dataset")

p2 <- ggplot(res_Alpha_Atlas, aes(x=alpha_atlas, y=alpha_array_all, fill = group, col = group)) +
  geom_point(data = res_Alpha_Atlas[is.na(res_Alpha_Atlas$group),], pch = 21, alpha = 0.05) +
  geom_point(data = res_Alpha_Atlas[!is.na(res_Alpha_Atlas$group),], pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  labs(title = "Probability of being hypervariable",
       x = "P(hv) considering all atlas data",
       y = "P(hv) considering all array data")

p3 <- ggplot(res_Alpha_Atlas, aes(x=alpha_atlas, y=alpha_array_3ind, fill = group, col = group)) +
  geom_point(data = res_Alpha_Atlas[is.na(res_Alpha_Atlas$group),], pch = 21, alpha = 0.05) +
  geom_point(data = res_Alpha_Atlas[!is.na(res_Alpha_Atlas$group),], pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  labs(title = "Probability of being hypervariable",
       x = "P(hv) considering all atlas data",
       y = "P(hv) considering array with 3 individuals per dataset only")

# --- Turn off legends inside plots ---
p1_clean <- p1 + theme(legend.position = "none")
p2_clean <- p2 + theme(legend.position = "none")
p3_clean <- p3 + theme(legend.position = "none")

# --- Extract one legend (e.g. from p1) ---
legend <- cowplot::get_legend(
  p1 + theme(legend.position = "bottom",
             legend.box = "horizontal",
             legend.title = element_blank(),
             legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
             legend.key = element_rect(fill = "white", color = NA))
)

# --- Arrange plots with legend as 4th panel ---
pdf(here("05_hvCpGalgorithm/figures/test2_compAtlasvsArrayalland3ind.pdf"), width = 10, height = 10)
cowplot::plot_grid(p1_clean, p2_clean, p3_clean, legend,
                   ncol = 2)  # grid layout: 2 cols × 2 rows
dev.off()

## Not the issue

############################################# 
## --- Test 3: bulk vs purified cells? --- ##
#############################################

## Run only cd4+ cd8+ in both and see

## Array CD4+ CD8+
load("/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Arrays/results_Arrays_CD4_CD8___352914CpGs_0_8p0_0_65p1.RData")

resArrayCD4CD8 <- results_Arrays_CD4_CD8___352914CpGs_0_8p0_0_65p1 %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "cpgprobe") %>%
  dplyr::rename(alpha_array_CD4CD8 = alpha)

resArrayCD4CD8$chrpos = hvCpGandControls$dictionary$hg38[
  match(resArrayCD4CD8$cpgprobe,
        hvCpGandControls$dictionary$illu450k)]

## Atlas CD4+ CD8+ on array positions
atlas_dt_CD4CD8 <- mergeAtlasRunBatches(parent_dir = "05_hvCpGalgorithm/resultsDir/10X_CD4+CD8+/", 
                                 analysis = "Atlas10X_CD4_CD8", alphaname = "alpha_atlas_CD4CD8")

## Merge
res_Alpha_Atlas_CD4CD8 <- full_join(resArrayCD4CD8, data.frame(atlas_dt_CD4CD8))
res_Alpha_Atlas_CD4CD8 <- left_join(res_Alpha_Atlas_CD4CD8, res_Alpha_Atlas)

sum(!is.na(res_Alpha_Atlas_CD4CD8$alpha_array_CD4CD8) & 
      !is.na(res_Alpha_Atlas_CD4CD8$alpha_atlas_CD4CD8))
# 155635 covered in both in enough samples/coverage

## Indicate the hvCpG of Maria and controls
res_Alpha_Atlas_CD4CD8$group <- "background"
res_Alpha_Atlas_CD4CD8$group[res_Alpha_Atlas_CD4CD8$chrpos %in% 
                      hvCpGandControls$DerakhshanhvCpGs_names] <- "hvCpG_Derakhshan"
res_Alpha_Atlas_CD4CD8$group[res_Alpha_Atlas_CD4CD8$chrpos %in% 
                      hvCpGandControls$mQTLcontrols_names] <- "mQTLcontrols"
res_Alpha_Atlas_CD4CD8$group <- as.factor(res_Alpha_Atlas_CD4CD8$group)

resplot <- subset(res_Alpha_Atlas_CD4CD8, !is.na(group))
resplot$group <- factor(resplot$group, levels = c("background","hvCpG_Derakhshan", "mQTLcontrols"))

p <- ggplot(resplot, aes(x = alpha_array_CD4CD8,
                         y = alpha_atlas_CD4CD8,
                         fill = group,
                         colour = group)) +
  geom_point(shape = 21, alpha = 0.1) +
  geom_point(data = resplot[!resplot$group %in% "background",],
             shape = 21, alpha = 0.4) +
  scale_fill_manual(values = c("background"       = "grey",
                               "hvCpG_Derakhshan" = "#DC3220",
                               "mQTLcontrols"     = "#005AB5")) +
  scale_colour_manual(values = c("background"       = "grey",
                                 "hvCpG_Derakhshan" = "#DC3220",
                                 "mQTLcontrols"     = "#005AB5")) +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.7,0.7),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(title = "Probability of being hypervariable",
       x = "P(hv) of array data CD4+ CD8+",
       y = "P(hv) of atlas data CD4+ CD8+")

# Now add marginal densities by group
p_with_dens <- ggMarginal(
  p,
  type = "density",
  groupFill = TRUE,
  groupColour = TRUE,
  alpha = 0.4
)

p1 <- p_with_dens

## Compare array all +-vs array CD4CD8
x <- full_join(resArrayCD4CD8, resCompArray)

p2 <- ggplot(x, aes(x=alpha_array_all, y=alpha_array_CD4CD8, fill = group, col = group)) +
  geom_point(data = x[is.na(x$group),], pch = 21, alpha = 0.05) +
  geom_point(data = x[!is.na(x$group),], pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.3,0.8),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(title = "Probability of being hypervariable",
       x = "P(hv) considering all array data",
       y = "P(hv) considering array CD4+ CD8+ groups only")

pdf(here("05_hvCpGalgorithm/figures/test3_compArrayAtlasCD4CD8.pdf"), width = 14, height = 7)
cowplot::plot_grid(p2, p1, ncol = 2)
dev.off()

################################# 
## --- Test 4: Sex effect? --- ##
#################################

## Male only atlas
atlas_dt_male <- mergeAtlasRunBatches(parent_dir = "05_hvCpGalgorithm/resultsDir/10X_males/", 
                                 analysis = "Atlas10X_males", alphaname = "alpha_atlas_males")

## Female only atlas
atlas_dt_female <- mergeAtlasRunBatches(parent_dir = "05_hvCpGalgorithm/resultsDir/10X_females/", 
                                      analysis = "Atlas10X_females", alphaname = "alpha_atlas_females")


atlas_dt_male_female <- full_join(atlas_dt_male, atlas_dt_female)

## Compare array all +-vs array CD4CD8
x <- left_join(atlas_dt_male_female, res_Alpha_Atlas)

p1 <- ggplot(x, aes(x=alpha_atlas_males, y=alpha_atlas_females, fill = group, col = group)) +
  geom_point(data = x[is.na(x$group),], pch = 21, alpha = 0.05) +
  geom_point(data = x[!is.na(x$group),], pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.3,0.8),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(x = "P(hv) considering male atlas data (12ds)",
       y = "P(hv) considering female atlas data (16ds)")

p2 <- ggplot(x, aes(x=alpha_atlas_males, y=alpha_atlas, fill = group, col = group)) +
  geom_point(data = x[is.na(x$group),], pch = 21, alpha = 0.05) +
  geom_point(data = x[!is.na(x$group),], pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.3,0.8),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(x = "P(hv) considering male atlas data (12ds)",
       y = "P(hv) considering all atlas data (46ds)")

p3 <- ggplot(x, aes(x=alpha_atlas_females, y=alpha_atlas, fill = group, col = group)) +
  geom_point(data = x[is.na(x$group),], pch = 21, alpha = 0.05) +
  geom_point(data = x[!is.na(x$group),], pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.3,0.8),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(x = "P(hv) considering female atlas data (16ds)",
       y = "P(hv) considering all atlas data (46ds)")

# --- Turn off legends inside plots ---
p1_clean <- p1 + theme(legend.position = "none")
p2_clean <- p2 + theme(legend.position = "none")
p3_clean <- p3 + theme(legend.position = "none")

# --- Extract one legend (e.g. from p1) ---
legend <- cowplot::get_legend(
  p1 + theme(legend.position = "bottom",
             legend.box = "horizontal",
             legend.title = element_blank(),
             legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
             legend.key = element_rect(fill = "white", color = NA))
)

# Combine plots
plots_grid <- plot_grid(
  p1_clean, p2_clean,
  p3_clean, legend,
  ncol = 2
)

# Add an overall title on top
final_plot <- ggdraw() +
  draw_label(
    "Probability of being hypervariable P(hv) in male, female, and combined WGBS atlas data",
    fontface = "bold",
    x = 0.5, y = 0.98, hjust = 0.5, vjust = 1, size = 18
  ) +
  draw_plot(plots_grid, x = 0, y = 0, width = 1, height = 0.95)

pdf(here("05_hvCpGalgorithm/figures/test4_compAtlasBothSexes.pdf"), width = 12, height = 12)
print(final_plot)
dev.off()

##############################################
## --- Test 5: batch correction effect? --- ##
##############################################

## Compare results array with either raw data uncorrected or corrected
load("/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Arrays_noCorrectionInRaw/results_Arrays_noCorrectionInRaw__406334CpGs_0_8p0_0_65p1.RData")

resArrayNoCor <- results_Arrays_noCorrectionInRaw__406334CpGs_0_8p0_0_65p1
rm(results_Arrays_noCorrectionInRaw__406334CpGs_0_8p0_0_65p1)

resArrayNoCor <- resArrayNoCor %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "cpgprobe") %>%
  dplyr::rename(alpha_array_nocor = alpha)

resArrayNoCor$chrpos = hvCpGandControls$dictionary$hg38[
  match(resArrayNoCor$cpgprobe,
        hvCpGandControls$dictionary$illu450k)]

resCommon_Array_Atlas_rawArray <- full_join(res_Alpha_Atlas, resArrayNoCor)

p1 <- ggplot(resCommon_Array_Atlas_rawArray, 
             aes(x=alpha_array_all, y=alpha_array_nocor, fill = group, col = group)) +
  geom_point(data = resCommon_Array_Atlas_rawArray[is.na(resCommon_Array_Atlas_rawArray$group),],
             pch = 21, alpha = 0.05) +
  geom_point(data = resCommon_Array_Atlas_rawArray[!is.na(resCommon_Array_Atlas_rawArray$group),],
             pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.3,0.8),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(title = "Probability of being hypervariable",
       x = "P(hv) considering array data after full correction",
       y = "P(hv) considering array data without correction")

p2 <- ggplot(resCommon_Array_Atlas_rawArray, 
             aes(x=alpha_atlas, y=alpha_array_nocor, fill = group, col = group)) +
  geom_point(data = resCommon_Array_Atlas_rawArray[is.na(resCommon_Array_Atlas_rawArray$group),],
             pch = 21, alpha = 0.05) +
  geom_point(data = resCommon_Array_Atlas_rawArray[!is.na(resCommon_Array_Atlas_rawArray$group),],
             pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.3,0.8),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(title = "Probability of being hypervariable",
       x = "P(hv) considering atlas WGBS data",
       y = "P(hv) considering array data without correction")

# --- Turn off legends inside plots ---
p1_clean <- p1 + theme(legend.position = "none")
p2_clean <- p2 + theme(legend.position = "none")

# --- Extract one legend (e.g. from p1) ---
legend <- cowplot::get_legend(
  p1 + theme(legend.position = "bottom",
             legend.box = "horizontal",
             legend.title = element_blank(),
             legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
             legend.key = element_rect(fill = "white", color = NA))
)

# --- Arrange plots with legend as 4th panel ---
pdf(here("05_hvCpGalgorithm/figures/test5_compArrayVsnocor.pdf"), width = 10, height = 7)
cowplot::plot_grid(p1_clean, p2_clean, legend,
                   ncol = 2, rel_heights = c(.9, .3))  # grid layout: 2 cols × 2 rows
dev.off()

## tbc





###############
## FIND PAX8
# Chromosome 2, NC_000002.12 (113215997..113278921,

atlas_dt <- atlas_dt %>%
  mutate(
    chr   = str_extract(chrpos, "chr[0-9XY]+"),
    start = as.numeric(str_extract(chrpos, "(?<=_)[0-9]+")),
    end   = as.numeric(str_extract(chrpos, "(?<=-)[0-9]+"))
  )

# 1. Convert to GRanges object
gr_atlas <- GRanges(
  seqnames = atlas_dt$chr,
  ranges = IRanges(start = atlas_dt$pos, end = atlas_dt$pos),
  mcols = atlas_dt
)

# 2. PAX8 as a GRanges
gr_PAX8 <- GRanges(
  seqnames = "chr2",
  ranges = IRanges(start = 113215997, end = 113278921)
)

# 3. Find overlaps between CpGs and the gene region
hits <- findOverlaps(gr_resCommonAlphaAtlas, gr_PAX8)

# 4. Extract matching rows from original df
df_hits <- gr_resCommonAlphaAtlas[queryHits(hits), ]

as.data.frame(df_hits)
