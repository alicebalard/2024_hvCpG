#############################################
## Overlap plot: Atlas (x) vs Array (y)    ##
#############################################
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

# Define parent folder containing all "Atlas_batchXXX" folders
parent_dir <- "05_hvCpGalgorithm/resultsDir/Atlas10X/"

# Get list of relevant RData files
rdata_files <- dir(here(parent_dir), pattern = "results_Atlas10X_100000CpGs_0_8p0_0_65p1\\.RData$", 
                   recursive = TRUE, full.names = TRUE)

## Check if all batches have ran
length(rdata_files) ## 230

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

rm(e, pb, all_cpg_values, obj, file, parent_dir, rdata_files)

atlas_dt <- dt[, .(name, alpha_atlas = alpha)]
setnames(atlas_dt, "name", "chrpos")

###################################### 
## --- Merge Array & Atlas data --- ##
######################################

resCommonAlphaAtlas <- dplyr::left_join(resCompArray, atlas_dt)

######################################### 
## --- Test 1: matching pos error? --- ##
#########################################

## Positions to test
## high alpha in array, low in atlas
## high alpha in atlas, low in array
pos2test <- c(resCommonAlphaAtlas[resCommonAlphaAtlas$alpha_array_all > 0.9 &
                                    resCommonAlphaAtlas$alpha_atlas < 0.1 &
                                    !is.na(resCommonAlphaAtlas$alpha_atlas ),][1,"chrpos"],
              resCommonAlphaAtlas[resCommonAlphaAtlas$alpha_array_all <0.1 &
                                    resCommonAlphaAtlas$alpha_atlas > 0.9 &
                                    !is.na(resCommonAlphaAtlas$alpha_atlas ),][1,"chrpos"])

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

################################## 
## --- Test 2: low N error? --- ##
##################################

p1 <- ggplot(resCommonAlphaAtlas, aes(x=alpha_array_all, y=alpha_array_3ind, fill = group, col = group)) +
  geom_point(data = resCommonAlphaAtlas[is.na(resCommonAlphaAtlas$group),], pch = 21, alpha = 0.05) +
  geom_point(data = resCommonAlphaAtlas[!is.na(resCommonAlphaAtlas$group),], pch = 21, alpha = 0.4) +
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

p2 <- ggplot(resCommonAlphaAtlas, aes(x=alpha_atlas, y=alpha_array_all, fill = group, col = group)) +
  geom_point(data = resCommonAlphaAtlas[is.na(resCommonAlphaAtlas$group),], pch = 21, alpha = 0.05) +
  geom_point(data = resCommonAlphaAtlas[!is.na(resCommonAlphaAtlas$group),], pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  labs(title = "Probability of being hypervariable",
       x = "P(hv) considering all atlas data",
       y = "P(hv) considering all array data")

p3 <- ggplot(resCommonAlphaAtlas, aes(x=alpha_atlas, y=alpha_array_3ind, fill = group, col = group)) +
  geom_point(data = resCommonAlphaAtlas[is.na(resCommonAlphaAtlas$group),], pch = 21, alpha = 0.05) +
  geom_point(data = resCommonAlphaAtlas[!is.na(resCommonAlphaAtlas$group),], pch = 21, alpha = 0.4) +
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
pdf(here("05_hvCpGalgorithm/figures/compAtlasvsArrayalland3ind.pdf"), width = 10, height = 10)
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

## Atlas CD4+ CD8+ on array positions
# Define parent folder containing all "Atlas_batchXXX" folders
parent_dir <- here("05_hvCpGalgorithm/resultsDir/10X_CD4+CD8+/")

# Get list of relevant RData files
rdata_files <- dir(parent_dir, pattern = "^results_Atlas10X.*CpGs_0_8p0_0_65p1\\.RData$",
  recursive = TRUE, full.names = TRUE
)

## Check if all batches have ran
length(rdata_files) 

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

rm(e, pb, all_cpg_values, obj, file, parent_dir, rdata_files)

## Merge
resArrayCD4CD8 <- data.frame(
    alpha_array_CD4CD8 = as.vector(results_Arrays_CD4_CD8___352914CpGs_0_8p0_0_65p1),
  cpgProbe = rownames(results_Arrays_CD4_CD8___352914CpGs_0_8p0_0_65p1), row.names = NULL)

resArrayCD4CD8$chrpos = hvCpGandControls$dictionary$hg38[
  match(resArrayCD4CD8$cpgProbe,hvCpGandControls$dictionary$illu450k)]

nrow(resArrayCD4CD8) # 352914
nrow(dt) # 216791

resbothCD4CD8 <- na.omit(full_join(resArrayCD4CD8, 
          data.frame(alpha_atlas_CD4CD8 = dt$alpha, chrpos = dt$name)))
nrow(resbothCD4CD8) # 155635 covered in both in enough samples/coverage

## Indicate the hvCpG of Maria and controls
resbothCD4CD8$group <- "background"
resbothCD4CD8$group[resbothCD4CD8$chrpos %in% 
                      hvCpGandControls$DerakhshanhvCpGs_names] <- "hvCpG_Derakhshan"
resbothCD4CD8$group[resbothCD4CD8$chrpos %in% 
                      hvCpGandControls$mQTLcontrols_names] <- "mQTLcontrols"
resbothCD4CD8$group <- as.factor(resbothCD4CD8$group)

resplot <- subset(resbothCD4CD8, !is.na(group))
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

pdf(here("05_hvCpGalgorithm/figures/CD4CD8ArrayAtlas.pdf"), width = 7, height = 7)
p_with_dens
dev.off()

## Compare array all +-vs array CD4CD8
x <- full_join(resArrayCD4CD8, resCompArray)

pdf(here("05_hvCpGalgorithm/figures/compArrayallVsCD4CD8.pdf"), width = 6, height = 6)
ggplot(x, aes(x=alpha_array_all, y=alpha_array_CD4CD8, fill = group, col = group)) +
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
dev.off()






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
