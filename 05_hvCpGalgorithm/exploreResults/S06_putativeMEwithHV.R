###################################################################
## Plot secondary analyses results of algorithm ran on atlas data ##
####################################################################
library(here)

## Load libraries
if (!exists("libLoaded")) {
  source(here("05_hvCpGalgorithm", "quiet_library.R"))}

## Load functions
if (!exists("functionsLoaded")) {
  source(here("05_hvCpGalgorithm/exploreResults", "functions.R"))}

## Load array results for comparison
if (!exists("resCompArray")) {
  source(here("05_hvCpGalgorithm/exploreResults/S02_analyseResultsArray_local.R"))}

## Load the set of previously tested MEs & vmeQTL
if (!exists("previousSIVprepared")) {
  source(here("05_hvCpGalgorithm/exploreResults/prepPreviousSIV.R"))}

##################################
## Save all data in RDS objects ##
##################################

for (file in list.files(here("05_hvCpGalgorithm/resultsDir/Atlas/"))){
  if (!file.exists(here(paste0("gitignore/fullres_", file)))){
    ## Add previous MEs including Maria's results if not sourced yet
    if (!exists("KesslerSIV_GRanges_hg38")) {
      source(here("05_hvCpGalgorithm/exploreResults/prepPreviousSIV.R"))
    }
    
    system.time(Atlas_dt <- prepAtlasdt(file))
    
    print("Number of CpG tested:")
    print(nrow(Atlas_dt))
    
    print(paste0("Saving results for ", file, "in ", here(paste0("gitignore/fullres_", file))))
    saveRDS(Atlas_dt, file = here(paste0("gitignore/fullres_", file)))
    print("Saved")
  }
}

#####################################
## one run per developmental layer ##
#####################################

endo = readRDS(here::here("gitignore/fullres_10X_12_endo"))
meso = readRDS(here::here("gitignore/fullres_10X_13_meso"))
ecto = readRDS(here::here("gitignore/fullres_10X_14_ecto"))
allLayers = readRDS(here::here("gitignore/fullres_Atlas10X"))

WGBS_Array_datasets <- read.csv(here("05_hvCpGalgorithm/figures/WGBS_Array_datasets.csv"))

table(WGBS_Array_datasets[WGBS_Array_datasets$assay %in% "atlas", "Germ.layer"])
# ectoderm endoderm mesoderm 
# 6       21       19 

##############################
## new candidate locus Matt ##
##############################

## Matt's data are in hg19
dataMatt <- readxl::read_xlsx(here("gitignore/DEGCAGS_intersect_repeats_Alice.xlsx"))

head(dataMatt)

# --- Download and import hg19 → hg38 chain file ---
chain_dir <- here("05_hvCpGalgorithm/dataPrev")
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

endo_gr <- GRanges(seqnames = paste0("chr", endo$chr),
                   ranges = IRanges(start = endo$pos, end = endo$pos), 
                   alpha = endo$alpha)
ecto_gr <- GRanges(seqnames = paste0("chr", ecto$chr),
                   ranges = IRanges(start = ecto$pos, end = ecto$pos), 
                   alpha = ecto$alpha)
meso_gr <- GRanges(seqnames = paste0("chr", meso$chr),
                   ranges = IRanges(start = meso$pos, end = meso$pos), 
                   alpha = meso$alpha)
all_gr <- GRanges(seqnames = paste0("chr", allLayers$chr),
                   ranges = IRanges(start = allLayers$pos, end = allLayers$pos), 
                   alpha = allLayers$alpha)

# Find overlapping ranges
overlaps <- findOverlaps(dataMatt_hg38_gr, endo_gr)
dataMatt_hg38_gr[queryHits(overlaps),]$alpha_endo <- endo_gr[subjectHits(overlaps),]$alpha

# Find overlapping ranges
overlaps <- findOverlaps(dataMatt_hg38_gr, meso_gr)
dataMatt_hg38_gr[queryHits(overlaps),]$alpha_meso <- meso_gr[subjectHits(overlaps),]$alpha

# Find overlapping ranges
overlaps <- findOverlaps(dataMatt_hg38_gr, ecto_gr)
dataMatt_hg38_gr[queryHits(overlaps),]$alpha_ecto <- ecto_gr[subjectHits(overlaps),]$alpha

# Find overlapping ranges
overlaps <- findOverlaps(dataMatt_hg38_gr, all_gr)
dataMatt_hg38_gr[queryHits(overlaps),]$alpha_all <- all_gr[subjectHits(overlaps),]$alpha

hist(dataMatt_hg38_gr$alpha_endo, breaks = 100)
hist(dataMatt_hg38_gr$alpha_ecto, breaks = 100)
hist(dataMatt_hg38_gr$alpha_meso, breaks = 100)
hist(dataMatt_hg38_gr$alpha_all, breaks = 100)

dataMatt_hg38_gr$is50pc <- (dataMatt_hg38_gr$alpha_endo > .5 & dataMatt_hg38_gr$alpha_meso > .5 & 
                              dataMatt_hg38_gr$alpha_ecto > .5)
dataMatt_hg38_gr$is40pc <- (dataMatt_hg38_gr$alpha_endo > .4 & dataMatt_hg38_gr$alpha_meso > .4 & 
                              dataMatt_hg38_gr$alpha_ecto > .4)

write.csv(dataMatt_hg38_gr, here("gitignore/DEGCAGS_intersect_repeats_Alice_withalpha.csv"), 
          row.names = F, quote = F)

#######################################
## Save intersection for alpha > 90% ##
#######################################

# 1) Overlap among all names (unfiltered)
sets_unfilt <- lapply(list(endo, ecto, meso, allLayers), function(df) df$name)
overlap <- Reduce(intersect, sets_unfilt)

message(paste0("There are ", length(overlap), " overlapping CpGs between these groups (unfiltered)."))
# There are 17474840 overlapping CpGs between these groups (unfiltered).

# 2) Apply cutoff AND keep only overlapping names (so sets are comparable on the same universe)
sets <- lapply(list(endo, ecto, meso, allLayers),
               function(df) df$name[df$alpha >= 0.9 & df$name %in% overlap])
names(sets) <- c("endo", "ecto", "meso", "allLayers")

# 3) Select the intersection
topIntersect90 <- Reduce(intersect, sets)

# Exactly in 1 group
only_endo <- setdiff(sets$endo, union(sets$ecto, sets$meso))
only_ecto <- setdiff(sets$ecto, union(sets$endo, sets$meso))
only_meso <- setdiff(sets$meso, union(sets$endo, sets$ecto))
exactly_one <- unique(c(only_endo, only_ecto, only_meso))

# Exactly in 2 groups
endo_ecto <- intersect(sets$endo, sets$ecto)
endo_ecto <- setdiff(endo_ecto, sets$meso)
endo_meso <- intersect(sets$endo, sets$meso)
endo_meso <- setdiff(endo_meso, sets$ecto)
ecto_meso <- intersect(sets$ecto, sets$meso)
ecto_meso <- setdiff(ecto_meso, sets$endo)
exactly_two <- unique(c(endo_ecto, endo_meso, ecto_meso))

# Combined: sites in exactly 1 or 2 groups
sites_1or2 <- unique(c(exactly_one, exactly_two))

message(paste0("Sites with alpha>=0.9 in exactly 1 group: ", length(exactly_one)))
message(paste0("Sites with alpha>=0.9 in exactly 2 groups: ", length(exactly_two)))
message(paste0("Total sites in 1 or 2 groups: ", length(sites_1or2)))

## Save
saveRDS(overlap, file = here("gitignore/overlapLayers.RDS"))
saveRDS(sites_1or2, file = here("gitignore/overlap.1or2.RDS"))
saveRDS(topIntersect90, file = here("gitignore/topIntersect90.RDS"))

##############################################
## What is the overlap for different alpha? ##
##############################################

# Compute overlap across any number of groups and plot a Venn diagram
plotMyVenn <- function(cutoff, ...) {
  groups <- list(...)                         # list of data.frames (each has name, alpha)
  
  # 1) Overlap among all names (unfiltered)
  sets_unfilt <- lapply(groups, function(df) df$name)
  overlap <- Reduce(intersect, sets_unfilt)
  
  message(paste0("There are ", length(overlap), " overlapping CpGs between these groups (unfiltered)."))
  
  # 2) Apply cutoff AND keep only overlapping names (so sets are comparable on the same universe)
  sets <- lapply(groups, function(df) df$name[df$alpha >= cutoff & df$name %in% overlap])
  
  # 3) Optional: add names to sets if you passed named arguments
  if (!is.null(dots <- match.call(expand.dots = FALSE)$...) && length(names(dots))) {
    nm <- names(dots)
    if (any(nzchar(nm))) names(sets) <- ifelse(nzchar(nm), nm, paste0("group", seq_along(sets)))
  }
  
  # 4) Plot
  p <- ggVennDiagram(sets, label_alpha = 0) +
    scale_fill_gradient2(low = "white", mid = "yellow", high = "red") +
    ggtitle(paste0("Pr(hv) ≥ ", cutoff), 
            subtitle = paste0("Sequenced CpGs N = ", length(overlap)))
  
  return(p)
}

p1 <- plotMyVenn(0.5, endo = endo, meso = meso, ecto = ecto, all = allLayers)
p2 <- plotMyVenn(0.75, endo = endo, meso = meso, ecto = ecto, all = allLayers)
p3 <- plotMyVenn(0.90, endo = endo, meso = meso, ecto = ecto, all = allLayers)

pdf(here("05_hvCpGalgorithm/figures/vennGermLayers/Venn_top90.pdf"), width = 16, height = 8)
cowplot::plot_grid(p1 + theme(legend.position = "none"),
                   p2 + theme(legend.position = "none"),
                   p3 + theme(legend.position = "none"), 
                   nrow = 1, 
                   labels = "Ectoderm (N=6 cell types), mesoderm (N=19), endoderm (N=21), and all combined")
dev.off()

####################################################################################
## Test enrichement of the most likely germ layer-universal hvCpG in previous MEs ##
####################################################################################
total <- allLayers$name[allLayers$name %in% endo$name]
total <- total[total %in% meso$name]
total <- total[total %in% ecto$name]

top_cpgs <- intersect(
  intersect(allLayers$name[allLayers$alpha > 0.9], endo$name[endo$alpha > 0.9]),
  intersect(meso$name[meso$alpha > 0.9],ecto$name[ecto$alpha > 0.9]))
print(paste("Found", length(top_cpgs), "overlapping high-alpha CpGs")) # 174494

listGR <- list(top90 = makeGRfromMyCpGPos(vec = top_cpgs, setname = "topCpGs"),
  allButTop90 = makeGRfromMyCpGPos(total[!total %in% top_cpgs], "allButTop90"))

# ---- Run it (ME sets in putativeME_GR$set will be tested separately)
res_quadrants <- test_enrichment_quadrants(listGR, putativeME_GR, me_col = "set")

# Order quadrants within each facet by log2OR
res_plot2 <- res_quadrants %>%
  mutate(
    log2OR = log2(odds_ratio),
    signif  = p_adj_BH < 0.05
  ) %>%
  group_by(ME_set) %>%
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
    subtitle = "Bars ordered by effect within ME set; dashed line = OR = 1"
  ) +
  facet_wrap(~ ME_set, scales = "free_x", nrow = 1) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

pdf(here("05_hvCpGalgorithm/figures/topCpGsEnrichME.pdf"), width = 9, height = 6)
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
    set = "Metastable Epiallele Set",
    n_total = "Total MEs",
    n_in_top90 = "In Top90",
    pc_in_top90 = "%",
    n_in_allButTop90 = "In Rest",
    pc_in_allButTop90 = "%",
    n_in_both = "In Both",
    pc_in_both = "%",
    n_in_neither = "In Neither",
    pc_in_neither = "%"
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_column_labels()
  ) %>%
  tab_options(
    table.font.size = 11,
    data_row.padding = px(3)
  )
