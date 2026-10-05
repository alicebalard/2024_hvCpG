######################################################
## Explore differences between array and atlas data ##
######################################################

#####################################################################
## Prepare
library(here)
## Load libraries
source(here("B_MultiTissues", "quiet_library.R"))

## Load functions
if (!exists("functionsLoaded")) {
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))}

# if (!exists("resArray3ind")) {
#   resArray3ind <- readRDS(here("B_MultiTissues/dataOut/resArray3ind_0_8p0_0_65p1.RDS"))
# }
# 
# if (!exists("table3layers_coveredIn3")) {
#   load(here(paste0("gitignore/table3layers_coveredIn3_26_08_26.Rda")))
# }

if (!exists("dico")) source(here("B_MultiTissues/03_exploreResults/makeProbes2GenDictionary.R"))
setDT(dico)

## Atlas tissue
if (!exists("tissueAtlas")) {
  parent_dir_atlas <- here("B_MultiTissues/resultsDir_gitIgnored/Atlas/Atlas10X_tissueAnalysis/")
  rds_files_atlas <- list.files(parent_dir_atlas, pattern = "\\p1.rds$", recursive = TRUE, full.names = TRUE)
  tissueAtlas <- do.call(rbind, lapply(rds_files_atlas, readRDS))
}

## Array tissue
if (!exists("tissueArray")) {
  parent_dir_array <- here("B_MultiTissues/resultsDir_gitIgnored/Arrays/tissue/")
  rds_files_array <- list.files(parent_dir_array, pattern = "\\p1.rds$", recursive = TRUE, full.names = TRUE)
  tissueArray <- readRDS(rds_files_array)
}

######################################################
## Goal: find CpGs the two platforms score very differently (high in one, low
## in the other), then explain WHICH TISSUES drive the discordance, using the
## per-tissue score matrices (tissueArray / tissueAtlas) bridged by `dico`.
######################################################

# ══════════════════════════════════════════════════════════════════════════════
# Join the two platforms on hg38 chr_pos and rank discordance
# ══════════════════════════════════════════════════════════════════════════════
# The array is keyed on probe ids; map probe -> chrpos_hg38 via dico so it shares
# the atlas coordinate space. (resArray3ind already has chrpos == hg38 chr_pos.)

# atlas <- as.data.table(table3layers_coveredIn3)[, .(chr_pos, atlas_score = logBF_per_ds_allLayers)]
# arr   <- as.data.table(resArray3ind)[, .(chr_pos = chrpos, array_score = logBF_per_ds)]
# 
# rm(table3layers_coveredIn3, resArray3ind)
# 
# m <- merge(atlas, arr, by = "chr_pos")
# saveRDS(m, here("gitignore/m_S07.RDS"))

m <- readRDS(here("gitignore/m_S07.RDS"))

message(nrow(m), " CpGs scored by BOTH platforms")
# 280664 CpGs scored by BOTH platforms

## signed discordance + a strict "opposite extremes" flag
m[, diff := atlas_score - array_score]
qa <- quantile(m$atlas_score, c(.5,.9), na.rm=TRUE)
qr <- quantile(m$array_score, c(.5,.9), na.rm=TRUE)
m[, class := fifelse(atlas_score >= qa[2] & array_score <= qr[1], "atlas-high / array-low",
                     fifelse(array_score >= qr[2] & atlas_score <= qa[1], "array-high / atlas-low",
                             fifelse(array_score >= qa[2] & atlas_score >= qa[2], "array-high / atlas-high",
                             "concordant")))]
print(m[, .N, by = class])
#                      class      N
# 1:              concordant 257772
# 2:  array-high / atlas-low   1165
# 3: array-high / atlas-high  21312
# 4:  atlas-high / array-low    415

classes <- c(
  "array-high / atlas-low",
  "atlas-high / array-low",
  "array-high / atlas-high"
)

reference_class <- "array-high / atlas-low"

# Select CpGs belonging to the classes of interest
selected_m <- unique(
  m[
    class %chin% classes,
    .(chr_pos, class)
  ],
  by = c("chr_pos", "class")
)

# # Extract the corresponding CpGs from tissueAtlas
# X <- tissueAtlas[rownames(tissueAtlas) %in% selected_m$chr_pos, , drop = FALSE]
# X <- as.data.table(X, keep.rownames = "cpg_site")
# 
# saveRDS(X, "../../gitignore/X.RDS")

X <- readRDS("../../gitignore/X.RDS")

# Wide -> long
X_long <- melt(X,id.vars = "cpg_site",variable.name = "tissue",
  value.name = "value", variable.factor = FALSE)

# Attach class
X_long <- merge(X_long, selected_m, by.x = "cpg_site", by.y = "chr_pos", all = FALSE)

# Attach germ layer
Loyfer <- read.csv(here("B_MultiTissues/dataIn/SupTab1_Loyfer2023.csv"))
Loyfer$tissue <- paste0(Loyfer$Source.Tissue, " - ", Loyfer$Cell.type)
X_long$GermLayer <- Loyfer$Germ.layer[match(X_long$tissue, Loyfer$tissue)]

# Calculate tissue-level statistics separately by class
tissue_summary <- X_long[!is.na(value), .(n = .N, 
                                          GermLayer=GermLayer,
    mean = mean(value),
    sd = if (.N > 1L) sd(value) else NA_real_,
    se = if (.N > 1L) sd(value) / sqrt(.N) else NA_real_,
    ci_lower = if (.N > 1L) {
      mean(value) -
        qt(0.975, df = .N - 1L) *
        sd(value) / sqrt(.N)
    } else {
      NA_real_
    },
    ci_upper = if (.N > 1L) {
      mean(value) +
        qt(0.975, df = .N - 1L) *
        sd(value) / sqrt(.N)
    } else {
      NA_real_
    }
  ),
  by = .(tissue, class)
]

tissue_order <- tissue_summary[class == reference_class][order(mean), tissue]

# Add tissues not present in the reference class at the end
tissue_order <- c(tissue_order,
  setdiff(unique(tissue_summary$tissue), tissue_order))

class_means <- X_long[!is.na(value),.(overall_mean = mean(value),
    n = .N),by = class]

p <- ggplot(head(tissue_summary, 100), aes(x = tissue,y = mean, colour = GermLayer,
    group = class)) +  
  geom_point(aes(shape = class), 
             position = position_dodge(width = 0.7), size = 2.5) +
  geom_errorbar(aes(ymin = ci_lower,ymax = ci_upper),
    position = position_dodge(width = 0.7), width = 0.25) +
  geom_hline(data = class_means,aes(yintercept = overall_mean,
      colour = class),
    linetype = "dashed",linewidth = 0.6,inherit.aes = FALSE) +
  coord_flip() +
  scale_colour_manual(
    values = c("array-high / atlas-low" = "#D55E00",
      "atlas-high / array-low" = "#0072B2",
      "array-high / atlas-high" = "grey40")) +
  labs(x = NULL,y = "Mean value", colour = "Class",
       title = "Mean tissue values by discordance class",
    subtitle = paste0("Tissues ordered by ", reference_class,
      "; dashed lines show class-wide means")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom")

print(p)

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/script07/DiffArrayAtlasbyTissue.png"),
  plot = p, width = 8, height = 8, dpi = 300, bg = "white")
