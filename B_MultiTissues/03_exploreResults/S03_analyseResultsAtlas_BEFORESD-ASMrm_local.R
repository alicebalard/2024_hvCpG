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

## --> Jump top "checkpoint" if table3layers_coveredIn3 has been saved in gitignore
table3layers_coveredIn3_saved = TRUE

if (table3layers_coveredIn3_saved == FALSE){
  
  ## Load array results
  if (!exists("resArray")) {
    resArray <- readRDS(here("B_MultiTissues/dataOut/resArray_0_8p0_0_65p1.RDS"))
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
    out_dir   = here("gitignore/resultsAtlasPrepared"),
    a0ornot   = TRUE
  ) {
    search_dir <- file.path(atlas_dir, file)
    out_path   <- file.path(out_dir,
                            paste0(res, p0, "p0_", p1, "p1_", file, ".rds"))
    
    pattern <- if (a0ornot)
      paste0("^results_.*", p0, "p0_", p1, "p1_.*a0\\.rds$")
    else
      paste0("^results_.*", p0, "p0_", p1, "p1\\.rds$")
    
    rds_files <- base::dir(search_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
    if (length(rds_files) == 0) {
      message("SKIP ", file, " [", p0, "/", p1, "] - no matching files."); return(invisible(NULL))
    }
    
    # completeness: one file per batch, up to the highest batch number, no gaps
    batch_nums <- as.integer(regmatches(dirname(rds_files),
                                        regexpr("(?<=Atlas_batch)\\d+", dirname(rds_files), perl = TRUE)))
    expected_n <- max(batch_nums, na.rm = TRUE)
    if (length(rds_files) != expected_n || !all(sort(batch_nums) == seq_len(expected_n))) {
      message("SKIP ", file, " - found ", length(rds_files), " files, batches run 1..",
              expected_n, " (missing: ",
              paste(setdiff(seq_len(expected_n), batch_nums), collapse = ","), ")")
      return(invisible(NULL))
    }
    
    # last batch should be the short remainder, not a full 250000 (truncation check)
    cpgs <- as.integer(regmatches(rds_files,
                                  regexpr("(?<=_)\\d+(?=CpGs)", rds_files, perl = TRUE)))
    ord  <- order(batch_nums)
    if (!is.na(tail(cpgs[ord], 1)) && tail(cpgs[ord], 1) == 250000) {
      message("SKIP ", file, " - last batch has 250000 CpGs (truncated?).")
      return(invisible(NULL))
    }
    
    message("OK   ", file, " - ", length(rds_files), " batches.")
    if (file.exists(out_path)) { message("     already prepared - skipping."); return(invisible(NULL)) }
    
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    Atlas_dt <- prepAtlasdt(file, p0, p1, atlas_dir = atlas_dir, mypattern = pattern)
    saveRDS(Atlas_dt, out_path)
    message("     saved: ", out_path)
  }
  
  subdirs <- list.files(here("B_MultiTissues/resultsDir_gitIgnored/Atlas/"))
  subdirs <- subdirs[!grepl("^PREVIOUS", subdirs)]
  for (subdir in subdirs) {
    savePrepedAtlasFile(subdir, p0 = "0_8", p1 = "0_65")
  }
  
  ## New August 2026
  # Saved: /home/alice/Documents/Research/GIT/2024_hvCpG/gitignore/resultsAtlasPrepared/
  
  ################################################################################
  ## Load layer-specific analyses                                               ##
  ################################################################################
  
  endo <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_12_endo.rds"))
  meso <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_meso.rds"))
  ecto <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_14_ecto.rds"))
  
  ################################################################################
  ## Compare algo ran on all datasets with geometric mean                       ##
  ################################################################################
  
  endoGR <- GRanges(seqnames = endo$chr,
                    ranges = IRanges(start = endo$pos, end = endo$pos),
                    logBF_per_ds_endo = endo$logBF_per_ds)
  
  ectoGR <- GRanges(seqnames = ecto$chr,
                    ranges = IRanges(start = ecto$pos, end = ecto$pos),
                    logBF_per_ds_ecto = ecto$logBF_per_ds)
  
  mesoGR <- GRanges(seqnames = meso$chr,
                    ranges = IRanges(start = meso$pos, end = meso$pos),
                    logBF_per_ds_meso = meso$logBF_per_ds)
  
  Atlas_dt <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_atlas_general.rds"))
  
  allLayersGR <- GRanges(seqnames = Atlas_dt$chr,
                         ranges = IRanges(start = Atlas_dt$pos, end = Atlas_dt$pos),
                         logBF_per_ds_allLayers = Atlas_dt$logBF_per_ds)
  
  ####################################################################
  ## Create a table with all CpG sites & pr(hv) for each germ layer ##
  ####################################################################
  
  ## 1. Create union of all unique CpG positions
  table3layers_coveredIn3 <- union(union(allLayersGR, union(ectoGR, mesoGR)), endoGR)
  
  ## 2. Use findOverlaps to map logBF_per_ds values back
  # we want endoGR[i] -> table3layers_coveredIn3[endoHits[[i]]] for each i, etc.
  
  endoHits <- findOverlaps(endoGR, table3layers_coveredIn3, select = "first")
  ectoHits <- findOverlaps(ectoGR, table3layers_coveredIn3, select = "first")
  mesoHits <- findOverlaps(mesoGR, table3layers_coveredIn3, select = "first")
  allLayersHits <- findOverlaps(allLayersGR, table3layers_coveredIn3, select = "first")
  
  # initialize columns with NA
  mcols(table3layers_coveredIn3)$logBF_per_ds_endo <- NA_real_
  mcols(table3layers_coveredIn3)$logBF_per_ds_ecto <- NA_real_
  mcols(table3layers_coveredIn3)$logBF_per_ds_meso <- NA_real_
  mcols(table3layers_coveredIn3)$logBF_per_ds_allLayers <- NA_real_
  
  # copy only the hits
  mcols(table3layers_coveredIn3)$logBF_per_ds_endo[endoHits]   <- mcols(endoGR)$logBF_per_ds_endo
  mcols(table3layers_coveredIn3)$logBF_per_ds_ecto[ectoHits]   <- mcols(ectoGR)$logBF_per_ds_ecto
  mcols(table3layers_coveredIn3)$logBF_per_ds_meso[mesoHits]   <- mcols(mesoGR)$logBF_per_ds_meso
  mcols(table3layers_coveredIn3)$logBF_per_ds_allLayers[allLayersHits]   <- mcols(allLayersGR)$logBF_per_ds_allLayers
  
  ## Add chr_pos column to identify positions
  table3layers_coveredIn3$chr_pos <- paste0("chr", table3layers_coveredIn3@seqnames, "_", table3layers_coveredIn3@ranges@start)
  
  ################################################################################
  ### SAVED ###
  # save(table3layers_coveredIn3, file = 
  #        here(paste0("gitignore/fulltable3layers_coveredIn3PreGeomMean_", format(Sys.Date(), "%d_%m_%y"), ".Rda")))
  
  # load("../../gitignore/fulltable3layers_coveredIn3PreGeomMean_10_08_26.Rda")
  
  ################################################################################
  ## Arithmetic sum to stay in interpretable log-BF units
  m <- as.matrix(mcols(table3layers_coveredIn3)[, c("logBF_per_ds_endo",
                                                    "logBF_per_ds_ecto",
                                                    "logBF_per_ds_meso")])
  
  table3layers_coveredIn3$logBF_per_ds_mean3layers <- rowMeans(m, na.rm = FALSE)
  
  ### SAVED ###
  save(table3layers_coveredIn3, file =
         here(paste0("gitignore/fulltable3layers_coveredIn3_", format(Sys.Date(), "%d_%m_%y"), ".Rda")))
  
  load(here(paste0("gitignore/fulltable3layers_10_08_26.Rda")))
  
  ## We select only the sites covered in the 3 germ layer analyses and will use logBF_per_ds
  m <- as.matrix(mcols(table3layers)[, c("logBF_per_ds_endo",
                                         "logBF_per_ds_ecto",
                                         "logBF_per_ds_meso")])
  
  table3layers$n_layers <- rowSums(!is.na(m))
  table(table3layers$n_layers)
  #     1        2        3 
  # 916.558  2.111.121 21.522.541 
  
  ### SAVED ###
  table3layers_coveredIn3 <- table3layers[table3layers$n_layers %in% 3,]
  save(table3layers_coveredIn3, file =
         here(paste0("gitignore/table3layers_coveredIn3_", format(Sys.Date(), "%d_%m_%y"), ".Rda")))
  
}

################################################################################
################################## CHECKPOINT ##################################
################################################################################
load(here(paste0("gitignore/table3layers_coveredIn3_10_08_26.Rda")))

################################################
## Save the top 99% quantile for logBF_per_ds ##
################################################

top99q <- quantile(table3layers_coveredIn3$logBF_per_ds_allLayers, probs = 0.99, na.rm = FALSE)

top99q_CpGs <- table3layers_coveredIn3[table3layers_coveredIn3$logBF_per_ds_allLayers >= top99q, ]$chr_pos

message(paste0("Total CpG sites: ", length(table3layers_coveredIn3)))
message(paste0("Total top99q CpG sites: ", length(top99q_CpGs), " (",
               round(length(top99q_CpGs)/length(table3layers_coveredIn3)*100,2), "% of total)"))
# Total CpG sites: 21522541
# Total top99q CpG sites: 215226 (1% of total)

###################################
## Enrichement in SD-ASM regions ##
###################################
NOTtop99q <- table3layers_coveredIn3$chr_pos[!table3layers_coveredIn3$chr_pos %in% top99q_CpGs]
source(here("B_MultiTissues/01_dataPrep/SD-ASMprep.R"))

## Are hypervariable CpGs more likely to sit in SD-ASM regions than other covered CpGs?

# --- build GRanges for foreground / background CpGs ---
mk_gr <- function(cpg) {
  chr <- sub("_.*", "", cpg)
  pos <- as.integer(sub(".*_", "", cpg))
  GRanges(chr, IRanges(pos, pos))
}
fg_gr <- mk_gr(top99q_CpGs)   # top 1%
bg_gr <- mk_gr(NOTtop99q)     # the rest of covered-in-3 (disjoint from fg)

# --- generic overlap-enrichment test against a region set ---
sdasm_fisher <- function(region_gr, label) {
  fg_in  <- sum(overlapsAny(fg_gr, region_gr, ignore.strand = TRUE))
  bg_in  <- sum(overlapsAny(bg_gr, region_gr, ignore.strand = TRUE))
  fg_out <- length(fg_gr) - fg_in
  bg_out <- length(bg_gr) - bg_in
  
  m  <- matrix(c(fg_in, fg_out, bg_in, bg_out), nrow = 2, byrow = TRUE)
  ft <- fisher.test(m, alternative = "two.sided")   # two-sided: detect enrichment OR depletion
  
  data.table(
    set        = label,
    fg_in      = fg_in, fg_total = length(fg_gr),
    bg_in      = bg_in, bg_total = length(bg_gr),
    fg_frac    = fg_in / length(fg_gr),
    bg_frac    = bg_in / length(bg_gr),
    odds_ratio = unname(ft$estimate),
    conf_low   = ft$conf.int[1],
    conf_high  = ft$conf.int[2],
    pvalue     = ft$p.value
  )
}

# --- 1. overall SD-ASM (any classification) ---
res_all <- sdasm_fisher(reduce(SDASM_GR), "SD-ASM (all)")

# --- 2. per classification ---
cats <- unique(SDASM_GR$classification)
res_cat <- rbindlist(lapply(cats, function(cl)
  sdasm_fisher(reduce(SDASM_GR[SDASM_GR$classification == cl]), cl)))

res <- rbind(res_all, res_cat)
res[, p.adj := p.adjust(pvalue, "BH")]
res[order(-odds_ratio)]
#                        set fg_in fg_total   bg_in bg_total     fg_frac      bg_frac odds_ratio   conf_low conf_high        pvalue         p.adj
# 1:                   ubiq.   408   215226    2862 21307315 0.001895682 0.0001343201 14.1432937 12.7128788 15.692384 1.427789e-296 2.379648e-296
# 2:                   other 10468   215226  566545 21307315 0.048637246 0.0265892253  1.8716496  1.8347200  1.909029  0.000000e+00  0.000000e+00
# 3:            SD-ASM (all) 19464   215226 1256398 21307315 0.090435170 0.0589655712  1.5867418  1.5632898  1.610523  0.000000e+00  0.000000e+00
# 4: tissue-specific demeth.  8191   215226  633116 21307315 0.038057670 0.0297135514  1.2919267  1.2633532  1.320942 6.146815e-105 7.683519e-105
# 5:            denovo meth.   397   215226   53875 21307315 0.001844573 0.0025284744  0.7290201  0.6587763  0.804765  4.897757e-11  4.897757e-11

# a positive OR means "hv CpGs are enriched for cis-genetic control relative to typical covered CpGs"

################################################################################
## What about layer specific hvCpGs?                                          ##
################################################################################

# ── Load ───────────────────────────────────────────────────────────────────
if (!exists("endo")){
  endo     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_12_endo.rds"))
  meso     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_meso.rds"))
  ecto     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_14_ecto.rds"))
  analyses <- list(endo = endo, meso = meso, ecto = ecto)
}

# ── Make wide table ─────────────────────────────────────────────────────────────
# keep only name + the score, renamed per layer, from each table
slim <- Map(function(dt, nm) {
  dt <- as.data.table(dt)
  setNames(dt[, .(name, logBF_per_ds)], c("name", nm))
}, analyses, names(analyses))

# key once, then join — much cheaper than merging full tables
wide <- Reduce(function(a, b) a[b, on = "name"], slim)

## Keep only CpGs tested in all 3 layers
wide <- wide[rowSums(is.na(wide)) == 0]
nrow(wide) # 21.522.541

################################################################################
# In the top 1% per layer analysis, and in bottom 50% for the other layers

classify_wide <- function(wide, hv_q = 0.99, nothv_q = 0.50) {
  q <- lapply(c("endo", "meso", "ecto"), function(layer) {
    c(hv    = unname(quantile(wide[[layer]], hv_q,    na.rm = TRUE)),
      notHV = unname(quantile(wide[[layer]], nothv_q, na.rm = TRUE)))
  })
  names(q) <- c("endo", "meso", "ecto")
  
  wide[, `:=`(
    HV_endo = endo > q$endo["hv"], HV_meso = meso > q$meso["hv"], HV_ecto = ecto > q$ecto["hv"],
    notHV_endo = endo < q$endo["notHV"], notHV_meso = meso < q$meso["notHV"], notHV_ecto = ecto < q$ecto["notHV"]
  )]
  
  wide[, category := fcase(
    HV_meso & HV_endo & HV_ecto,           "ME",
    HV_meso & notHV_endo & notHV_ecto,     "Meso_specific",
    HV_endo & notHV_meso & notHV_ecto,     "Endo_specific",
    HV_ecto & notHV_meso & notHV_endo,     "Ecto_specific",
    notHV_meso & notHV_endo & notHV_ecto,  "constitutive",
    default =                              "ambiguous"
  )]
  attr(wide, "thresholds") <- q
  wide
}

wide   <- classify_wide(wide);   print(table(wide$category))
# ambiguous  constitutive Ecto_specific Endo_specific            ME Meso_specific 
# 14310278       7148471           390            17         63296            89

# ══════════════════════════════════════════════════════════════════════════════
# SD-ASM enrichment per germ-layer category (ME / *_specific / constitutive / ...)
# Same test as sdasm_fisher, but foreground/background passed as arguments so we
# can ask, for each category: is THIS category enriched in SD-ASM vs all other CpGs?
# ══════════════════════════════════════════════════════════════════════════════

# mk_gr() already defined earlier in the script; reuse it.

# --- generalised: foreground and background CpG-GRanges are arguments ---
sdasm_fisher2 <- function(region_gr, fg_gr, bg_gr, label) {
  fg_in  <- sum(overlapsAny(fg_gr, region_gr, ignore.strand = TRUE))
  bg_in  <- sum(overlapsAny(bg_gr, region_gr, ignore.strand = TRUE))
  fg_out <- length(fg_gr) - fg_in
  bg_out <- length(bg_gr) - bg_in
  
  m  <- matrix(c(fg_in, fg_out, bg_in, bg_out), nrow = 2, byrow = TRUE)
  ft <- fisher.test(m, alternative = "two.sided")
  
  data.table(
    fg_in      = fg_in, fg_total = length(fg_gr),
    bg_in      = bg_in, bg_total = length(bg_gr),
    fg_frac    = fg_in / length(fg_gr),
    bg_frac    = bg_in / length(bg_gr),
    odds_ratio = unname(ft$estimate),
    conf_low   = ft$conf.int[1],
    conf_high  = ft$conf.int[2],
    pvalue     = ft$p.value
  )
}

# --- pre-reduce the SD-ASM region sets ONCE (overall + per classification) ---
sdasm_sets <- c(
  list(`SD-ASM (all)` = reduce(SDASM_GR)),
  lapply(split(SDASM_GR, SDASM_GR$classification), reduce)
)

# --- germ-layer categories to test (skip 'ambiguous' as the residual bucket) ---
# `wide` carries chr_pos + category (from classify_wide)
setDT(wide)
cats_to_test <- setdiff(unique(wide$category), "ambiguous")

# --- loop categories × SD-ASM sets ---
res_list <- lapply(cats_to_test, function(cat) {
  fg_cpg <- wide[category == cat,  name]
  bg_cpg <- wide[category != cat & !is.na(category), name]
  
  fg_gr <- mk_gr(fg_cpg)
  bg_gr <- mk_gr(bg_cpg)
  
  rbindlist(lapply(names(sdasm_sets), function(sdclass) {
    out <- sdasm_fisher2(sdasm_sets[[sdclass]], fg_gr, bg_gr,
                         label = paste(cat, sdclass, sep = " | "))
    out[, `:=`(category = cat, sdasm_class = sdclass)]
    out
  }))
})

res_cat <- rbindlist(res_list)
res_cat[, p.adj := p.adjust(pvalue, "BH")]                        # correct across ALL tests
setcolorder(res_cat, c("category", "sdasm_class", "odds_ratio",
                       "conf_low", "conf_high", "fg_in", "fg_total",
                       "bg_in", "bg_total", "fg_frac", "bg_frac",
                       "pvalue", "p.adj"))

# view: strongest enrichment first, within each category
res_cat[order(category, -odds_ratio)]

# quick pivot of just the odds ratios (category × SD-ASM class)
dcast(res_cat, category ~ sdasm_class, value.var = "odds_ratio")

## Plot
pd <- copy(res_cat)

# recode category names for display
cat_rename <- c(
  constitutive  = "constitutive",
  ME            = "systemic hv",
  Endo_specific = "endoderm hv",
  Meso_specific = "mesoderm hv",
  Ecto_specific = "ectoderm hv"
)
pd[, category := factor(cat_rename[category], levels = cat_rename)]

pd[, sdasm_class := factor(sdasm_class,
                           levels = c("SD-ASM (all)", "ubiq.", "tissue-specific demeth.",
                                      "denovo meth.", "other"))]

# per-category N (fg_total is constant within a category) -> strip label "Category (N=…)"
pd[, cat_lab := sprintf("%s\n(N = %s)", category,
                        format(fg_total, big.mark = ",", trim = TRUE))]

# order facets by each category's overall SD-ASM OR
cat_order <- pd[sdasm_class == "SD-ASM (all)"][order(odds_ratio), cat_lab]
pd[, cat_lab := factor(cat_lab, levels = cat_order)]

pd[, sig := ifelse(p.adj < 0.05, "FDR < 0.05", "n.s.")]
pd[, or_plot := fifelse(odds_ratio == 0, NA_real_, odds_ratio)]

# lock the y order (reverse so the first level sits at the TOP of each panel)
pd[, sdasm_class := factor(sdasm_class,
                           levels = rev(c("SD-ASM (all)", "ubiq.", "tissue-specific demeth.",
                                          "denovo meth.", "other")))]

plot <- ggplot(pd, aes(x = or_plot, y = sdasm_class, colour = sig)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = pmax(conf_low, 1e-2), xmax = pmin(conf_high, 1e3)),
                height = 0.25, na.rm = TRUE) +
  geom_point(na.rm = TRUE) +
  scale_size_continuous(trans = "log10", range = c(1.5, 4), name = "n CpGs") +
  facet_grid(cat_lab ~ ., switch = "y") +      # drop scales="free_y", space="free_y"
  scale_y_discrete(drop = FALSE) +             # keep all 5 levels even if a panel lacks a point
  scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10),
                labels = c("0.1", "0.3", "1", "3", "10")) +
  scale_colour_manual(values = c("FDR < 0.05" = "#DC3220", "n.s." = "grey60"),
                      name = NULL) +
  labs(x = "Odds ratio (category vs rest, log scale)", y = NULL,
       title = "SD-ASM enrichment of hypervariability categories") +
  theme_bw(base_size = 12) +
  theme(strip.placement   = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        panel.grid.minor  = element_blank())

plot
ggplot2::ggsave(
  filename = here::here(paste0("B_MultiTissues/dataOut/figures/script03/SD-ASMenrichment.pdf")),
  plot = plot, width = 10, height = 6)
