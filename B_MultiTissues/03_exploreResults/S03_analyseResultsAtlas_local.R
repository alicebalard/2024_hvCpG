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
  
  ###################
  ## rmMultSamples ## 
  ###################
  
  if (!file.exists(here(paste0("B_MultiTissues/dataOut/figures/script03/rmMultipleSamples.pdf")))){
    # Some individuals have multiple cells sampled. Does that affect our results? NOPE
    plotrmtest <- makeCompPlot(
      minplot = 1000000, # plot all
      X = readRDS(here::here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_02_atlas_general.rds")),
      Y = readRDS(here::here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_rmMultSamples.rds")),
      whichAlphaX = "logBF_per_ds",
      whichAlphaY = "logBF_per_ds",
      title = "Effect of keeping only one cell type per individual (1M CPGs plotted)",
      xlab = "logBF_per_ds calculated on WGBS atlas datasets",
      ylab = "logBF_per_ds calculated on WGBS atlas datasets keeping one sample/individual only")
    
    ggplot2::ggsave(
      filename = here::here(paste0("B_MultiTissues/dataOut/figures/script03/rmMultipleSamples.pdf")),
      plot = plotrmtest, width = 10, height = 10)
  }
  
  ################################################################################
  ## Load layer-specific analyses                                               ##
  ################################################################################
  
  endo <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_12_endo.rds"))
  meso <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_meso.rds"))
  ecto <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_14_ecto.rds"))
  
  ################################
  ## Test: are 6 groups enough? ##
  ################################
  
  endo6gp <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_12_2_endo6gp.rds"))
  meso6gp <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_2_meso6gp.rds"))
  
  if (!file.exists(here::here(
    "B_MultiTissues/dataOut/figures/script03/correlation_endomesoFullvsReduced6gp.pdf"))){
    ## Use data table to handle large data
    setDT(endo6gp)
    setDT(endo)
    
    x <- endo6gp[, .(name, logBF_per_ds_6gp = logBF_per_ds)]; setkey(x, name)
    y <- endo[, .(name, logBF_per_ds_endo = logBF_per_ds)]; setkey(y, name)
    
    m <- x[y, nomatch = 0]   # keeps matched names only
    mycor <- cor(m$logBF_per_ds_6gp, m$logBF_per_ds_endo, use = "complete.obs")
    set.seed(1234)
    p1 <- ggplot(m[sample(nrow(m), 100000),], aes(x = logBF_per_ds_6gp, y = logBF_per_ds_endo))+
      geom_point(pch = 21, alpha = 0.1) +
      geom_smooth() +
      geom_smooth(method = "lm") +
      theme_minimal(base_size = 14) +
      ylim(c(0,1)) +
      annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
      labs(title = "Score of hypervariability in WGBS atlas endoderm cell types",
           subtitle = "(100k random CpG plotted)",
           x = "logBF per dataset calculated on a subset of cell types (N=6)",
           y = "logBF per dataset calculated on all cell types (N=21)")
    
    setDT(meso6gp)
    setDT(meso)
    
    x <- meso6gp[, .(name, logBF_per_ds_6gp = logBF_per_ds)]; setkey(x, name)
    y <- meso[, .(name, logBF_per_ds_meso = logBF_per_ds)]; setkey(y, name)
    
    m <- x[y, nomatch = 0]   # keeps matched names only
    mycor <- cor(m$logBF_per_ds_6gp, m$logBF_per_ds_meso, use = "complete.obs")
    
    p2 <- ggplot(m[sample(nrow(m), 100000),], aes(x = logBF_per_ds_6gp, y = logBF_per_ds_meso))+
      geom_point(pch = 21, alpha = 0.1) +
      geom_smooth() +
      geom_smooth(method = "lm") +
      theme_minimal(base_size = 14) +
      ylim(c(0,1)) +
      annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
      labs(title = "Score of hypervariability in WGBS atlas mesoderm cell types",
           subtitle = "(100k random CpG plotted)",
           x = "logBF per dataset calculated on a subset of cell types (N=6)",
           y = "logBF per dataset calculated on all cell types (N=19)")
    
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
  
  ################################################################################
  ## Plot the difference between the mean between layers and the logBF per ds 
  ## calculated on all datasets together
  ################################################################################
  
  if (!file.exists(here("B_MultiTissues/dataOut/figures/script03/mean3layers_vs_all.png"))){
    df <- as.data.frame(table3layers)
    
    set.seed(1234)
    p <- ggplot(df[sample(nrow(df), 100000),],
                aes(x = logBF_per_ds_mean3layers, y = logBF_per_ds_allLayers)) +
      geom_point(pch = 21, alpha = 0.1) +
      geom_abline(slope = 1) +
      theme_minimal(base_size = 14) +
      labs(title = "Hypervariability score using either all layers jointly or separately",
           x = "Average of the three layers LogBF per dataset",
           y = "LogBF per dataset calculated on all cell types",
           subtitle = "(100k random CpG plotted)")
    
    ggplot2::ggsave(
      filename = here::here("B_MultiTissues/dataOut/figures/script03/mean3layers_vs_all.png"),
      plot = p, width = 8, height = 8,
      dpi = 300, bg = "white")
    
    rm(df, p)
  }
  
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

##################################
## Make GR object as data.table ##
##################################

table3layers_coveredIn3dt <- as.data.table(table3layers_coveredIn3)
names(table3layers_coveredIn3dt)[names(table3layers_coveredIn3dt) %in% "seqnames"] <- "chr"
names(table3layers_coveredIn3dt)[names(table3layers_coveredIn3dt) %in% "start"] <- "pos"

# Compute cumulative position offsets for Manhattan plot
setorder(table3layers_coveredIn3dt, chr, pos)

offsets <- table3layers_coveredIn3dt[, .(max_pos = max(pos, na.rm = TRUE)), by = chr]
offsets[, cum_offset := c(0, head(cumsum(as.numeric(max_pos)), -1))]

table3layers_coveredIn3dt <- merge(table3layers_coveredIn3dt,
                                   offsets[, .(chr, cum_offset)], 
                                   by = "chr", all.x = TRUE, sort = FALSE)

## Mark group membership in dt
table3layers_coveredIn3dt[, group := NA_character_]
table3layers_coveredIn3dt[chr_pos %in% DerakhshanhvCpGs_hg38, group := "hvCpG_Derakhshan"]
table3layers_coveredIn3dt[chr_pos %in% mQTLcontrols_hg38, group := "mQTLcontrols"]

# Convert to integer/numeric if not already
table3layers_coveredIn3dt[, cum_offset := as.numeric(cum_offset)]
table3layers_coveredIn3dt[, pos2 := pos + cum_offset]

#################################################
## Figure "Compare with Derakhshan's results"  ##
#################################################

if (!file.exists(here("B_MultiTissues/dataOut/figures/script03/CompareWithResultsDerakhshan.png"))){
  
  ## Load array results
  if (!exists("resArray")) {
    resArray <- readRDS(here("B_MultiTissues/dataOut/resArray_0_8p0_0_65p1.RDS"))
  }
  
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
  
  # Merge hvCpG logBF_per_ds
  res <- table3layers_coveredIn3dt[!is.na(group)]
  
  message("We shift the score so it starts at zero")
  res[["logBF_per_ds"]] <- res[["logBF_per_ds_allLayers"]] + abs(min(res[["logBF_per_ds_allLayers"]]))
  
  hv_logBF_per_ds <- res[, c("chr_pos", "logBF_per_ds")]
  colnames(hv_logBF_per_ds) <- c("hvCpG", "logBF_per_ds_hvCpG")
  
  # Merge control logBF_per_ds
  ctrl_logBF_per_ds <- res[, c("chr_pos", "logBF_per_ds")]
  colnames(ctrl_logBF_per_ds) <- c("control", "logBF_per_ds_control")
  
  # Join everything
  merged <- pairs %>%
    left_join(hv_logBF_per_ds, by = "hvCpG") %>%
    left_join(ctrl_logBF_per_ds, by = "control") %>%
    mutate(difflogBF_per_ds=logBF_per_ds_hvCpG-logBF_per_ds_control)
  
  merged <- merged %>%
    mutate(chr = str_extract(hvCpG, "^chr[0-9XYM]+"))%>%
    filter(!is.na(difflogBF_per_ds))
  
  # DifferenceOfProbabilityForhvCpG-matching_controlInAtlas
  pdiffhv_controls <- ggplot(merged, aes(x="", y=difflogBF_per_ds))+
    geom_jitter(data=merged[merged$difflogBF_per_ds>=0,], col="black", alpha=.5)+
    geom_jitter(data=merged[merged$difflogBF_per_ds<0,], fill="yellow",col="black",pch=21, alpha=.5)+
    geom_violin(width=.5, fill = "grey", alpha=.8) +
    geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
    theme_minimal(base_size = 16)+
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
    ggtitle("Difference of score between Derakhshan\nhvCpGs and matching controls")+
    ylab("Difference of logBF per ds shifted to start at zero")
  
  plotManhattan_Derakh <- plotManhattanFromdt(
    table3layers_coveredIn3dt[!chr %in% c("M", "Y"),], plotDerakhshan = TRUE, centro = centro) +
    ggtitle("Hypervariability score along the genome")
  
  compPlotArrayAtlas <- makeCompPlot(
    X = resArray,
    Y = here::here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_atlas_general.rds"),
    whichX = "logBF_per_ds", whichY = "logBF_per_ds",
    title = "Array vs Atlas",
    xlab = "logBF per ds (array datasets)",
    ylab = "logBF per ds (WGBS atlas datasets)",
    minplot = 1e5)
  
  compPlotArrayAtlasMESO <- makeCompPlot(
    X = resArray,
    Y = here::here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_meso.rds"),
    whichX = "logBF_per_ds", whichY = "logBF_per_ds",
    title = "Array vs Atlas (mesoderm)",
    xlab = "logBF per ds (array datasets)",
    ylab = "logBF per ds (WGBS atlas datasets, mesoderm-derived)",
    minplot = 1e5)
  
  row2 <- cowplot::plot_grid(pdiffhv_controls, compPlotArrayAtlas, compPlotArrayAtlasMESO, 
                             labels = c("B", "C","D"), ncol = 3)
  
  ggsave(here("B_MultiTissues/dataOut/figures/script03/CompareWithResultsDerakhshan.png"),
         plot = cowplot::plot_grid(plotManhattan_Derakh, row2,
                                   nrow = 2, rel_heights = c(1,1.5), labels = c("A", "")),
         width = 18, height = 11, dpi = 300, bg = "white")
}

################################################################################
## Figure 3: MappingVariability                                               ##
################################################################################

####################################
## Plot Manhattan of logBF_per_ds ##

if (!file.exists(here(paste0("gitignore/plotManhattan_noDerakh.RDS")))){
  plotManhattan_noDerakh <- plotManhattanFromdt(
    table3layers_coveredIn3dt[!chr %in% c("M", "Y"),], plotDerakhshan = FALSE, centro = centro) +
    ggtitle("Hypervariability score along the genome")
  saveRDS(plotManhattan_noDerakh, here(paste0("gitignore/plotManhattan_noDerakh.RDS")))
}

##########################################
## What are the gaps in Manhattan plot? ##
# Compute the gap between consecutive CpGs on the same chromosome
table3layers_coveredIn3dt[, gap := pos - data.table::shift(pos), by = chr]

# Identify large gaps (>= 500k bp)
gaps_dt <- table3layers_coveredIn3dt[gap >= 500000, .(
  chr,
  gap_start = data.table::shift(pos),
  gap_end = pos,
  gap_size = gap
)]

# Drop first NA (since shift introduces one per chromosome)
gaps_dt[!is.na(gap_size)]

# chr gap_start   gap_end gap_size
# <fctr>     <int>     <int>    <int>
#   1:      1        NA 124793275  2292432
# 2:      1 124793275 143184605 18000029
# 3:      2 143184605  91406100  1003595
# 4:      5  91406100  49592147  2283407
# 5:      9  49592147  60518620 15041410
# 6:     14  60518620  18223731  2127270
# 7:     16  18223731  46380693  8100344
# 8:     19  46380693  27240939  2332752
# 9:     21  27240939   6070102   677715
# 10:     21   6070102  12966132  2151768
# 11:     22  12966132  15158090  2253761
# 12:      X  15158090  60274012   786392
# 13:      X  60274012  61918466   981691
# 14:      Y  61918466   5043719   914090
# 15:      Y   5043719   5765474   721626
# 16:      Y   5765474   6533721   768247
# 17:      Y   6533721   7395909   860254
# 18:      Y   7395909   8138220   742288
# 19:      Y   8138220  10087019  1948799
# 20:      Y  10087019  13725309  1975980
# 21:      Y  13725309  17064135  3338784
# 22:      Y  17064135  26436653  9337728
# 23:      Y  26436653  56677947 30022457

####################################
## Mitochondrial DNAm variability ##
####################################
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-023-09541-9?utm_source=chatgpt.com
## Near absence of 5mC, so expected

ggplot() +
  geom_point(data = table3layers_coveredIn3dt[table3layers_coveredIn3dt$chr == "M",], 
             aes(x = pos2, y = logBF_per_ds_allLayers),
             color = "black", size = 1, alpha = .5)+
  theme_classic() + theme(legend.position = "none") +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

###################################
## Y chromosome DNAm variability ##
###################################
ggplot() +
  geom_point(data = table3layers_coveredIn3dt[table3layers_coveredIn3dt$chr == "Y",], 
             aes(x = pos2, y = logBF_per_ds_allLayers),
             color = "black", size = 1, alpha = .5)+
  theme_classic() + theme(legend.position = "none") +
  labs(x = "Chromosome", y = "Probability of being a hvCpG")+
  theme_minimal(base_size = 14)

# High logBF_per_ds in 3 regions

#######################################################
## Test enrichment of features for high logBF_per_ds ##
#######################################################

if (!file.exists(here(paste0("gitignore/pfeatures.RDS")))){
  
  # Create GRanges
  gr_cpg <- GRanges(
    seqnames = paste0("chr", table3layers_coveredIn3dt$chr),
    ranges = IRanges(start = table3layers_coveredIn3dt$pos, end = table3layers_coveredIn3dt$pos),
    logBF_per_ds = table3layers_coveredIn3dt$logBF_per_ds_allLayers
  )
  length(gr_cpg) # 21.522.541 CpGs
  
  # restrict to autosomes and chr X
  gr_cpg <- gr_cpg[gr_cpg@seqnames %in% c(paste0("chr", 1:22), "chrX"),]
  length(gr_cpg) # 21.517.622 CpGs
  
  # Import bed file
  bed_features <- genomation::readTranscriptFeatures(here("gitignore/hg38_GENCODE_V47.bed"))
  
  # Annotate CpGs and see which regions have higher logBF_per_ds_allLayers (takes long)
  anno_result <- genomation::annotateWithGeneParts(
    target = gr_cpg, feature = bed_features)
  
  saveRDS(anno_result, "../../gitignore/anno_result.RDS")
  
  anno_result@perc.of.OlapFeat
  # promoter     exon   intron 
  # 96.57276 73.63779 89.41976 9 
  
  ## Add info from annotation to our GRange object
  gr_cpg$featureType <- ifelse(anno_result@members[, "prom"] == 1, "promoter",
                               ifelse(anno_result@members[, "exon"] == 1, "exon",
                                      ifelse(anno_result@members[, "intron"] == 1, "intron", "intergenic")))
  
  mcols(gr_cpg) %>% as.data.frame() %>%
    dplyr::group_by(featureType) %>%
    dplyr::summarise(meanlogBF_per_ds = mean(logBF_per_ds, na.rm=T),
                     medianlogBF_per_ds = median(logBF_per_ds, na.rm=T))
  # # A tibble: 4 × 3
  # featureType meanlogBF_per_ds medianlogBF_per_ds
  # <chr>                  <dbl>              <dbl>
  # 1 exon                  -0.307             -0.399
  # 2 intergenic            -0.190             -0.248
  # 3 intron                -0.291             -0.379
  # 4 promoter              -0.326             -0.432
  
  dt <- as.data.table(mcols(gr_cpg))[!is.na(logBF_per_ds)]
  
  # global quantile breaks:
  probs <- seq(0, 1, 0.1)
  qs    <- quantile(dt$logBF_per_ds, probs, na.rm = TRUE)
  
  dt[, band := cut(logBF_per_ds, breaks = qs, include.lowest = TRUE,
                   labels = c("Bottom 10%", "10–20th", "20–30th", "30–40th",
                              "40–50th", "50–60th", "60–70th", "70–80th", "80–90th", "Top 10%"))]
  
  frac <- dt[!is.na(band), .N, by = .(featureType, band)][
    , pct := 100 * N / sum(N), by = featureType]     # denominator = all CpGs in the feature
  
  pfeatures <- ggplot(frac, aes(featureType, pct, fill = band)) +
    geom_col(position = "dodge") +
    labs(y = "% of CpGs in feature", x = NULL, fill = "logBF per ds\n(percentile)") +
    scale_fill_viridis_d() +
    theme_minimal(base_size = 14) +
    geom_hline(yintercept = 10, linetype = 3, linewidth = 0.3) +
    ggtitle("Distribution of CpG hypervariability across feature types")
  pfeatures
  
  saveRDS(pfeatures, here(paste0("gitignore/pfeatures.RDS")))
}

# Distribution of CpG hypervariability across feature types. CpGs were binned into quantiles 
# of the per-dataset log Bayes factor (logBF per ds) computed genome-wide, and for each
# feature type the percentage of its CpGs falling in each upper quantile is shown 
# (the bottom quintile is omitted; percentages are of all CpGs in the feature). 
# Bars exceeding 10% (the expected share under no enrichment) indicate feature types enrichement.

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

saveRDS(top99q_CpGs, file = here("gitignore/top99q_CpGs_august26.RDS"))
## To use for testFetalSIV_ingp5.R

#######################
## Enrichement in TE ##
#######################

if (!file.exists(here("gitignore/TEplot.RDS"))){
  # UCSC RepeatMasker annotations (Oct2022) for Human (hg38) from AnnotationHub
  ah <- AnnotationHub()
  query(ah, c("UCSC", "RepeatMasker", "Homo sapiens"))
  
  # Retrieve the desired resource, UCSC RepeatMasker annotations for hg38:
  rmskhg38 <- ah[["AH111333"]]
  
  table(mcols(rmskhg38)$repFamily)
  
  te_classes <- c("LINE", "SINE", "LTR", "DNA", "RC", "Retroposon")
  te_regions <- rmskhg38[mcols(rmskhg38)$repClass %in% te_classes]
  
  # View summary
  table(mcols(te_regions)$repFamily)
  length(te_regions)  # Total TE regions
  
  top99q_CpGs_GR <- makeGRfromMyCpGPos(top99q_CpGs, "top99q_CpGs")
  totalSites_GR  <- makeGRfromMyCpGPos(table3layers_coveredIn3dt$chr_pos, "totalSites")
  
  # Strict background = non-hvCpG sites only (must be disjoint from top99q_CpGs)
  bg_only_GR <- makeGRfromMyCpGPos(table3layers_coveredIn3dt$chr_pos[
    !table3layers_coveredIn3dt$chr_pos %in% top99q_CpGs], "bg_only")
  
  # Enrichment/depletion of `target` CpGs vs `background` CpGs inside a set of repeat regions.
  fisher_test_te <- function(te_gr, target, background,
                             label = "TE", nameTarget = "foreground") {
    fg_in  <- sum(overlapsAny(target,     te_gr, ignore.strand = TRUE))
    bg_in  <- sum(overlapsAny(background, te_gr, ignore.strand = TRUE))
    fg_out <- length(target)     - fg_in
    bg_out <- length(background) - bg_in
    
    contingency <- matrix(c(fg_in, fg_out, bg_in, bg_out),
                          nrow = 2, byrow = TRUE,
                          dimnames = list(c(nameTarget, "background"),
                                          c(paste0("in_", label), paste0("not_in_", label))))
    test <- fisher.test(contingency, alternative = "two.sided")   # two-sided -> finite CI + detects depletion
    list(label       = label,
         contingency = contingency,
         pvalue      = test$p.value,
         odds_ratio  = unname(test$estimate),
         conf_low    = unname(test$conf.int[1]),
         conf_high   = unname(test$conf.int[2]))
  }
  
  # columns to pull from each result into a table
  res_cols <- c("label", "pvalue", "odds_ratio", "conf_low", "conf_high")
  
  ## ---- All TEs pooled ----
  allTE_test <- fisher_test_te(te_regions, top99q_CpGs_GR, bg_only_GR,
                               label = "TE", nameTarget = "top99q_CpGs")
  allTE_test
  # $label
  # [1] "TE"
  # 
  # $contingency
  # in_TE not_in_TE
  # top99q_CpGs   116812     98414
  # background  11104516  10202799
  # 
  # $pvalue
  # [1] 1.321875e-88
  # 
  # $odds_ratio
  # [1] 1.090531
  # 
  # $conf_low
  # [1] 1.081324
  # 
  # $conf_high
  # [1] 1.099918
  
  ## ---- Per repClass ----
  te_by_class <- split(te_regions, mcols(te_regions)$repClass)
  class_res <- lapply(names(te_by_class), function(cl)
    fisher_test_te(te_by_class[[cl]], top99q_CpGs_GR, bg_only_GR,
                   label = cl, nameTarget = "top99q_CpGs"))
  class_dt <- data.table::rbindlist(lapply(class_res, `[`, res_cols))
  class_dt[, p.adj := p.adjust(pvalue, "BH")]
  class_dt[order(-odds_ratio)]
  # label        pvalue odds_ratio  conf_low conf_high         p.adj
  # <char>         <num>      <num>     <num>     <num>         <num>
  # 1:       LINE  0.000000e+00  1.5007019 1.4839041 1.5176769  0.000000e+00
  # 2:         RC  6.083903e-02  1.4566919 0.9573108 2.1265146  7.300683e-02
  # 3:        LTR 4.911794e-294  1.3151295 1.2965659 1.3338986 9.823587e-294
  # 4:        DNA  1.773608e-20  1.1265118 1.0988141 1.1547253  2.660413e-20
  # 5: Retroposon  2.870427e-01  1.0732033 0.9382938 1.2222976  2.870427e-01
  # 6:       SINE  0.000000e+00  0.7523114 0.7448020 0.7598449  0.000000e+00
  
  ## ---- Per repFamily (drop tiny/uncertain families so ORs are stable) ----
  fam_tab  <- table(mcols(te_regions)$repFamily)
  fam_keep <- names(fam_tab)[fam_tab >= 1000 & !grepl("\\?$", names(fam_tab))]
  
  fam_by  <- split(te_regions, mcols(te_regions)$repFamily)
  fam_res <- lapply(fam_keep, function(f)
    fisher_test_te(fam_by[[f]], top99q_CpGs_GR, bg_only_GR,
                   label = f, nameTarget = "top99q_CpGs"))
  
  res_dt <- data.table::rbindlist(lapply(fam_res, `[`, res_cols))
  res_dt[, p.adj := p.adjust(pvalue, "BH")]           # correct across families
  res_dt[order(-odds_ratio)]
  
  # label        pvalue odds_ratio  conf_low conf_high         p.adj
  # <char>         <num>      <num>     <num>     <num>         <num>
  #   1:            L1  0.000000e+00  1.6187370 1.5991236 1.6385663  0.000000e+00
  # 2:      Helitron  6.083903e-02  1.4566919 0.9573108 2.1265146  1.303693e-01
  # 3:          ERVK  3.462892e-30  1.4196923 1.3399359 1.5030963  1.731446e-29
  # 4:     ERVL-MaLR 9.454857e-108  1.3025413 1.2732764 1.3323807 9.454857e-107
  # 5:     MULE-MuDR  8.998519e-02  1.2861086 0.9528308 1.6991305  1.799704e-01
  # 6:          ERV1  5.404470e-98  1.2817104 1.2532836 1.3106481  4.053352e-97
  # 7:          ERVL  2.510347e-49  1.2775112 1.2378570 1.3181502  1.506208e-48
  # 8: TcMar-Mariner  1.477872e-03  1.2458473 1.0874009 1.4210912  4.030559e-03
  # 9:           hAT  9.892430e-02  1.2328956 0.9520584 1.5712670  1.854831e-01
  # 10:    hAT-Tip100  1.637336e-05  1.2060434 1.1082202 1.3102336  4.912008e-05
  # 11: hAT-Blackjack  1.943938e-02  1.1808174 1.0256216 1.3530594  4.859846e-02
  # 12:         Gypsy  4.636022e-02  1.1520838 0.9992608 1.3218340  1.069851e-01
  # 13:           LTR  4.661397e-01  1.1205208 0.7745040 1.5693880  5.593677e-01
  # 14:  TcMar-Tigger  2.262577e-07  1.1199032 1.0732461 1.1681062  8.484663e-07
  # 15:     TcMar-Tc2  3.016151e-01  1.1198648 0.8946401 1.3848991  4.112933e-01
  # 16:   hAT-Charlie  5.688018e-08  1.1067295 1.0671752 1.1474161  2.437722e-07
  # 17:           SVA  2.870427e-01  1.0732033 0.9382938 1.2222976  4.100610e-01
  # 18:            L2  7.620033e-06  1.0609828 1.0338610 1.0885584  2.540011e-05
  # 19:           CR1  3.179310e-01  1.0415198 0.9599120 1.1282027  4.146926e-01
  # 20:          tRNA  9.031703e-01  1.0011871 0.5818425 1.6079355  1.000000e+00
  # 21:         RTE-X  1.000000e+00  0.9999340 0.8493747 1.1696736  1.000000e+00
  # 22:     5S-Deu-L2  1.000000e+00  0.9935032 0.5773893 1.5955355  1.000000e+00
  # 23:           MIR  1.825613e-01  0.9807094 0.9529342 1.0090565  2.956557e-01
  # 24:      tRNA-RTE  5.277820e-01  0.8850044 0.6189119 1.2273666  6.089792e-01
  # 25:      PiggyBac  4.116123e-01  0.8522759 0.5820681 1.2053285  5.145153e-01
  # 26:      RTE-BovB  1.872486e-01  0.8426688 0.6496571 1.0753592  2.956557e-01
  # 27:      Penelope  1.000000e+00  0.8389777 0.2276490 2.1623876  1.000000e+00
  # 28:           Alu  0.000000e+00  0.7406438 0.7330074 0.7483663  0.000000e+00
  # 29:        hAT-Ac  2.114208e-01  0.7276029 0.4230976 1.1675562  3.171312e-01
  # 30:           DNA  1.094098e-01  0.5917307 0.2832258 1.0909646  1.930761e-01
  
  res_dt[label %in% c("ERV1", "ERVK")]                # your two candidates
  
  # label       pvalue odds_ratio conf_low conf_high        p.adj
  # <char>        <num>      <num>    <num>     <num>        <num>
  # 1:   ERV1 5.404470e-98   1.281710 1.253284  1.310648 4.053352e-97
  # 2:   ERVK 3.462892e-30   1.419692 1.339936  1.503096 1.731446e-29
  
  ## ---- Plot: OR with 95% CI, ordered, coloured by significance ----
  res_plot <- res_dt[order(odds_ratio)]
  res_plot[, label := factor(label, levels = label)]  # lock the OR order
  res_plot[, sig := ifelse(p.adj < 0.05, "FDR < 0.05", "n.s.")]
  
  TEplot <- ggplot(res_plot, aes(odds_ratio, label, colour = sig)) +
    geom_vline(xintercept = 1, linetype = 3) +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.25) +
    geom_point(size = 3) +
    scale_x_log10(breaks = c(0,.5, 1, 1.5,2)) +                                 # OR is multiplicative -> log axis
    scale_colour_manual(values = c("FDR < 0.05" = "#DC3220", "n.s." = "grey60")) +
    labs(x = "Odds ratio (top99q hvCpG vs background, log scale)",
         y = NULL, colour = NULL,
         title = "TE family enrichment among 1% top hypervariable CpGs") +
    theme_minimal(base_size = 12)
  
  saveRDS(TEplot, file = here("gitignore/TEplot.RDS"))
}

#################################
## plot manhattan, feature, TE ##
#################################

if (!file.exists(here("B_MultiTissues/dataOut/figures/script03/MappingVariability.png"))){
  plotManhattan_noDerakh <- readRDS(here(paste0("gitignore/plotManhattan_noDerakh.RDS")))
  pfeatures <- readRDS(here(paste0("gitignore/pfeatures.RDS")))
  TEplot <- readRDS(here("gitignore/TEplot.RDS"))
  
  ggplot2::ggsave(
    filename = here::here(
      "B_MultiTissues/dataOut/figures/script03/MappingVariability.png"),
    plot = plot_grid(plotManhattan_noDerakh,
                     plot_grid(
                       pfeatures, TEplot, 
                       ncol= 2,
                       labels = c("B", "C")
                     ), nrow = 2,
                     labels = c("A", "")),
    width = 18, height = 10,
    dpi = 300, bg = "white")
}

###########################################
## Compare our results with previous MEs ##
###########################################

if (!file.exists(here("B_MultiTissues/dataOut/figures/script03/CompareWithpreviousMEs.png"))){
  
  ## Use the GR object with analyses in the 3 layers
  
  # Fix chromosome names in geomMeanGR (1 -> chr1)
  seqlevels(table3layers_coveredIn3) <- paste0("chr", seqlevels(table3layers_coveredIn3))
  
  sets <- list(
    mQTLcontrols = makeGRfromMyCpGPos(vec = mQTLcontrols_hg38, setname = "mQTLcontrols"),
    HarrisSIV = HarrisSIV_hg38_GR,
    VanBaakSIV = VanBaakSIV_hg38_GR,
    VanBaakESS = VanBaakESS_hg38_GR,
    KesslerSIV = KesslerSIV_GRanges_hg38,
    GunasekaraCorSIV = corSIV_GRanges_hg38,
    DerakhshanhvCpGs = DerakhshanhvCpGs_hg38_GR
  )
  
  ## Associate a colour to a group
  group_cols <- c(
    "background"         = "#999999",
    "mQTLcontrols"       = "#000000",
    "HarrisSIV"          = RColorBrewer::brewer.pal(8, "Set2")[1],
    "KesslerSIV"         = RColorBrewer::brewer.pal(8, "Set2")[2],
    "DerakhshanhvCpGs"   = RColorBrewer::brewer.pal(8, "Set2")[3],
    "GunasekaraCorSIV"   = RColorBrewer::brewer.pal(8, "Set2")[4],
    "VanBaakESS"         = RColorBrewer::brewer.pal(8, "Set2")[5],
    "top99q"         = RColorBrewer::brewer.pal(8, "Set2")[6],
    "VanBaakSIV"         = RColorBrewer::brewer.pal(8, "Set2")[7]
  )
  
  # ── ME overlap ────────────────────────────────────────────
  MEsetdt <- make_MEsetdt(sets, GR = table3layers_coveredIn3)
  
  MEsetdt <- na.omit(MEsetdt) ## 69979
  
  # Set controls as baseline
  MEsetdt[, ME := relevel(factor(ME), ref = "mQTLcontrols")]
  
  ## Statistical comparisons of alpha between MEs
  fit <- lm(logBF_per_ds ~ ME, data = MEsetdt)
  emm <- emmeans(fit, ~ ME)
  contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = "mQTLcontrols", adjust = "sidak") %>%
    as.data.frame()
  
  contrasts <- contrasts %>%
    mutate(ME = contrast,
           ME_name = sub(" - mQTLcontrols$", "", contrast),   # match colour to the ME group being compared
           lower = estimate - 1.96 * SE,
           upper = estimate + 1.96 * SE)
  
  pMElogBF_per_ds <- ggplot(MEsetdt, aes(x = ME, y = logBF_per_ds)) +
    geom_jitter(aes(colour = ME), size = 3, alpha = .2) +
    geom_violin(aes(colour = ME)) +
    geom_boxplot(aes(colour = ME), width = .1) +
    scale_colour_manual(values = group_cols, name = "CpG set") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    ylab("Hypervariability score (logBF per ds)")
  
  pcontrast <- ggplot(contrasts, aes(x = ME, y = estimate, colour = ME_name)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_colour_manual(values = group_cols, name = "CpG set") +
    coord_flip() +
    labs(y = "Difference in hypervariability score (logBF per ds) vs mQTLcontrols", x = "") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ## pdecay - legend inside the plot, bottom-left corner
  pdecay <- plot_decay_curve(MEsetdt) +
    scale_colour_manual(values = group_cols, name = "CpG set")
  
  # ── Save key objects for S07 ──────────────────────────────────────────────────
  saveRDS(MEsetdt,            here("gitignore/MEsetdt.rds"))
  
  ####################################################################################
  ## Test enrichement of the most likely germ layer-universal hvCpG in previous MEs ##
  ####################################################################################
  if (!exists("listGR")){
    listGR <- list(top99q = makeGRfromMyCpGPos(vec = top99q_CpGs, setname = "top99q"),
                   allButTop99q = makeGRfromMyCpGPos(
                     setdiff(table3layers_coveredIn3$chr_pos, top99q_CpGs), "allButTop99q"))
  }
  
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
  
  plot_top99qCpGsEnrichME <- ggplot(res_plot2, aes(x = quadrant_ord, y = log2OR, fill = signif)) +
    geom_col(width = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c("black", "grey")) +
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
  
  ## Saved earlier:
  # saveRDS(top99q_CpGs, file = here("gitignore/top99q_CpGs_august26.RDS"))
  ## Then run on ingp5: testFetalSIV_ingp5.R
  
  plots <- readRDS(here("gitignore/intercorrelationSIVfetal_sepSIV.rds"))
  
  pinterlayer_corr <- ggplot(plots$interlayer_corr, aes(x=group, y=interlayer_r, group = group, fill = group))+
    geom_violin(width=1.4) +
    geom_boxplot(width=0.3, fill="white") +
    scale_fill_manual(values = group_cols) +
    theme_minimal(base_size = 14) +
    labs(y = "Mean inter-germ layer correlation\n(Pearson's r)")
  
  pinterindividual_var <- ggplot(plots$CpG_summary, aes(x = interindividual_var, fill = group)) +
    geom_density(alpha = .8)+
    scale_fill_manual(values = group_cols) +
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
    )
  
  upperRow <- plot_grid(
    plot_grid(
      pMElogBF_per_ds, 
      pcontrast + theme_minimal(base_size = 17) + 
        theme(plot.title = element_text(size=16), legend.position = "none"),
      nrow = 2, 
      labels = c("A. Distribution of the hypervariability score for each CpG set", 
                 "B. Comparison of previous CPG sets groups to mQTLcontrols"),
      label_size = 18, label_x = 0, hjust = 0),
    pdecay + theme_minimal(base_size = 20) + 
      theme(legend.position = "inside",
            legend.position.inside = c(.7, .6),
            legend.justification = c(0, 0),
            legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3)),
    ncol = 2,
    rel_widths = c(1, 1),
    labels = c("", "C. Decay curve of hypervariability score (logBF per ds) per threshold"), 
    label_size = 18, label_x = 0, hjust = 0
  )
  
  SIV_plot <- plot_grid(
    pinterlayer_corr + theme_minimal(base_size = 20) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1),
            axis.title.x = element_blank(), legend.position = "none"),
    plot_grid(pinterindividual_var +
                labs(fill = "CpG set"), 
              pbinned + theme_minimal(base_size = 16) +
                theme(axis.text.x = element_text(angle = 20, hjust = 1))+
                labs(colour = "CpG set", fill = "CpG set"), 
              ncol = 1, align = "v",
              labels = c("E. Densities of interindividual variation per CpG within the fetal data, by set", 
                         "F. Inter-germ-layer correlation per interindividual variation, binned"),
              label_size = 18, label_x = 0, hjust = 0),
    ncol = 2,
    rel_widths = c(1, 1),
    labels = c("D. Mean inter-germ-layer correlation for each CpG set", ""),
    label_size = 18, label_x = 0, hjust = 0
  )
  
  final_plot <- plot_grid(
    upperRow,
    SIV_plot,
    ncol = 1,
    rel_heights = c(1, 1)
  )
  
  ggplot2::ggsave(
    filename = here::here(
      "B_MultiTissues/dataOut/figures/script03/CompareWithpreviousMEs.png"),
    plot = final_plot, width = 22, height = 20,  dpi = 300, bg = "white")
}

############################################################
## How many of each putative ME is actually in the top90? ##
############################################################

# Overlaps of ME ranges with CpG sites (top99q and allButTop99q)
hits_top <- findOverlaps(putativeME_GR, listGR$top99q, ignore.strand = TRUE)
hits_all <- findOverlaps(putativeME_GR, listGR$allButTop99q, ignore.strand = TRUE)

# Logical flags per ME range (query index)
overlap_top <- logical(length(putativeME_GR))
overlap_top[unique(queryHits(hits_top))] <- TRUE

overlap_all <- logical(length(putativeME_GR))
overlap_all[unique(queryHits(hits_all))] <- TRUE

# Build a small data.frame with one row per ME range
df_me <- as.data.frame(mcols(putativeME_GR)) |>
  mutate(
    in_top99q       = overlap_top,
    in_allButTop99q = overlap_all
  )

# Summaries per set
summary_df <- df_me |>
  group_by(set) |>
  summarise(
    n_total = n(),
    n_in_top99q = sum(in_top99q),
    pc_in_top99q = n_in_top99q/n_total*100,
    n_in_allButTop99q = sum(in_allButTop99q),
    pc_in_allButTop99q = n_in_allButTop99q/n_total*100,
    n_in_both = sum(in_top99q & in_allButTop99q),   # non 0 if region rather than CpG
    pc_in_both = n_in_both/n_total*100,
    n_in_neither = n_total - n_in_top99q - n_in_allButTop99q + n_in_both,
    pc_in_neither = n_in_neither/n_total*100)

## Format pretty
summary_df %>%
  mutate(across(starts_with("pc_"), ~ scales::percent(.x / 100, accuracy = 0.1))) %>%
  gt() %>%
  fmt_number(columns = starts_with("n_"), decimals = 0) %>%
  cols_label(
    set = "Putative ME set",
    n_total = "Total CpG sites or regions",
    n_in_top99q = "N in top 99q CpGs",
    pc_in_top99q = "%",
    n_in_allButTop99q = "N in all but top 99q CpGs",
    pc_in_allButTop99q = "%",
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

getEnrichCentroTelo <- function(){
  top <- top99q_CpGs
  hv_dt <- parse_cpg_names(top)
  totalSites <- table3layers_coveredIn3$chr_pos
  total_dt <- parse_cpg_names(totalSites)  # all background sites
  
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
  
  enrich_test(hv_in_centro,  bg_in_centro,  top, totalSites, "Centromere")
  enrich_test(hv_in_subtelo, bg_in_subtelo, top, totalSites, "Subtelomere (1Mb)")
}

getEnrichCentroTelo()
# ── Centromere ──
# hvCpGs in region:      3413 / 215226 (1.59%)
# Background in region:  190918 / 21522541 (0.89%)
# Fold enrichment:       1.79
# Fisher p (one-sided):  3.072442e-210
# Odds ratio:            1.8
# 
# ── Subtelomere (1Mb) ──
# hvCpGs in region:      4840 / 215226 (2.25%)
# Background in region:  482013 / 21522541 (2.24%)
# Fold enrichment:       1
# Fisher p (one-sided):  0.3887386
# Odds ratio:            1

#######################
## Test GO of top 1% ##
#######################

totalSites <- table3layers_coveredIn3$chr_pos
top99q_CpGs

# Method. ClusterProfiler
## 1. Keep CpGs in regions where at least 2 CpGs are in 50bp distance to each other
## 2. annotate with associated genes (in gene body or +/- 10kb from TSS)
## 3. run GO term enrichment with clusterProfiler::enrichGO

minimum_CpG_per_cluster = 2

## Create universe
universe <- annotateCpGs_txdb(
  clusterCpGs(totalSites, max_gap = 50, min_size = minimum_CpG_per_cluster),
  tss_window = 10000)

print(paste0("Gene universe contains ", length(universe), " genes"))
## "Gene universe contains 32717 genes"

## Annotate
# Length-controlled
resAnnot_top99q_lenCtrl <- CpG_GO_pipeline_lengthControlled(
  top99q_CpGs, universe = universe,
  max_gap = 50, min_size = minimum_CpG_per_cluster, tss_window = 10000,
  control_length = TRUE)
# Found 1641 Entrez genes
# Controlling for gene length...
# Median gene length — foreground: 5,514 bp, universe: 3,100 bp, ratio: 1.78
# Length-matched universe: 16779 genes (was 32717)
# Running GO enrichment...

# Which terms survive
terms_lenCtrl <- resAnnot_top99q_lenCtrl$BP@result %>% filter(p.adjust < 0.05) %>% pull(ID)

cat("BP terms after  length control:", length(terms_lenCtrl),  "\n")
# BP terms after length control: 15 

df_all <- purrr::imap_dfr(resAnnot_top99q_lenCtrl, function(er, ont_name) {
  if (is.null(er) || nrow(er@result) == 0) return(tibble())
  as_tibble(er@result) |>
    mutate(group_raw = "top90SNPrm", ontology = ont_name)
}) |> bind_rows()

df_all <- df_all |>
  mutate(
    group = "top99q",
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

write.csv(df_sig, file = here("B_MultiTissues/dataOut/df_sig_GOtop99q.csv"),
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
       title = "GO Enrichment: top99q (FDR < 0.05)") +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold")
  )

p

ggplot2::ggsave(
  filename = here::here(
    "B_MultiTissues/dataOut/figures/script03/GOplottop99q.png"),
  plot = p, width = 10, height = 5,  dpi = 300, bg = "white")
