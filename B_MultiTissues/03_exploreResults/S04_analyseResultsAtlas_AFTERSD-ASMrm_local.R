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

variant <- "SNP_SDASMrm"
prep_dir <- here("gitignore/resultsAtlasPrepared", variant)
fr <- function(f) readRDS(file.path(prep_dir, f)) # short reader

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
    here("B_MultiTissues/dataIn/sample_metadata.tsv"), 
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
  
  subdirs <- list.files(file.path(here("B_MultiTissues/resultsDir_gitIgnored/Atlas"), variant))
  subdirs <- subdirs[!grepl("^PREVIOUS", subdirs)]
  for (subdir in subdirs) {
    savePrepedAtlasFile(subdir, p0 = "0_8", p1 = "0_65", variant = variant)
  }
  
  # Saved: /home/alice/Documents/Research/GIT/2024_hvCpG/gitignore/resultsAtlasPrepared/SNP_SDASMrm
  
  ###################
  ## rmMultSamples ## 
  ###################
  
  if (!file.exists(here(paste0("B_MultiTissues/dataOut/figures/script04/rmMultipleSamples.pdf")))){
    # Some individuals have multiple cells sampled. Does that affect our results? NOPE
    plotrmtest <- makeCompPlot(
      minplot = 100000, # plot all
      X = file.path(prep_dir, "fullres_0_8p0_0_65p1_atlas_general.rds"),
      Y = file.path(prep_dir, "fullres_0_8p0_0_65p1_02_rmMultSamples.rds"),
      whichX = "logBF_per_ds",
      whichY = "logBF_per_ds",
      title = "Effect of keeping only one cell type per individual (100k CPGs plotted)",
      xlab = "Hypervariability score calculated on WGBS atlas datasets",
      ylab = "Hypervariability score calculated on WGBS atlas datasets keeping 1 cell type/individual only")
    
    ggplot2::ggsave(
      filename = here::here(paste0("B_MultiTissues/dataOut/figures/script04/rmMultipleSamples.pdf")),
      plot = plotrmtest, width = 9, height = 9)
  }
  
  ################################################################################
  ## Load layer-specific analyses                                               ##
  ################################################################################
  
  endo <- fr("fullres_0_8p0_0_65p1_12_endo.rds")
  meso <- fr("fullres_0_8p0_0_65p1_13_meso.rds")
  ecto <- fr("fullres_0_8p0_0_65p1_14_ecto.rds")
  
  ################################
  ## Test: are 6 groups enough? ##
  ################################
  
  endo6gp <- fr("fullres_0_8p0_0_65p1_12_2_endo6gp.rds")
  meso6gp <- fr("fullres_0_8p0_0_65p1_13_2_meso6gp.rds")
  # ecto6gp is ecto
  
  if (!file.exists(here::here(
    "B_MultiTissues/dataOut/figures/script04/correlation_endomesoFullvsReduced6gp.pdf"))){
    
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
      geom_smooth(method = "lm") +
      theme_minimal(base_size = 14) +
      annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
      labs(title = "Hypervariability score in WGBS atlas endoderm cell types",
           subtitle = "(100k random CpG plotted)",
           x = "Hypervariability score calculated on a subset of cell types (N=6)",
           y = "Hypervariability score calculated on all cell types (N=21)")
    
    setDT(meso6gp)
    setDT(meso)
    
    x <- meso6gp[, .(name, logBF_per_ds_6gp = logBF_per_ds)]; setkey(x, name)
    y <- meso[, .(name, logBF_per_ds_meso = logBF_per_ds)]; setkey(y, name)
    
    m <- x[y, nomatch = 0]   # keeps matched names only
    mycor <- cor(m$logBF_per_ds_6gp, m$logBF_per_ds_meso, use = "complete.obs")
    
    p2 <- ggplot(m[sample(nrow(m), 100000),], aes(x = logBF_per_ds_6gp, y = logBF_per_ds_meso))+
      geom_point(pch = 21, alpha = 0.1) +
      geom_smooth(method = "lm") +
      theme_minimal(base_size = 14) +
      annotate("text", x = .2, y = .9, label = sprintf("Pearson correlation: r = %.2f\n", mycor)) +
      labs(title = "Hypervariability score in WGBS atlas mesoderm cell types",
           subtitle = "(100k random CpG plotted)",
           x = "Hypervariability score calculated on a subset of cell types (N=6)",
           y = "Hypervariability score calculated on all cell types (N=19)")
    
    ggplot2::ggsave(
      filename = here::here(paste0("B_MultiTissues/dataOut/figures/script04/correlation_endomesoFullvsReduced6gp.pdf")),
      plot = cowplot::plot_grid(p1, p2, labels = c("A", "B")), width = 17, height = 8)
    
    rm(x, y, m, mycor, p1, p2)
  }
  
  ######################################################################
  ## Compare algo ran on all datasets with mean                       ##
  ######################################################################
  
  endoGR <- GRanges(seqnames = endo$chr,
                    ranges = IRanges(start = endo$pos, end = endo$pos),
                    logBF_per_ds_endo = endo$logBF_per_ds)
  
  ectoGR <- GRanges(seqnames = ecto$chr,
                    ranges = IRanges(start = ecto$pos, end = ecto$pos),
                    logBF_per_ds_ecto = ecto$logBF_per_ds)
  
  mesoGR <- GRanges(seqnames = meso$chr,
                    ranges = IRanges(start = meso$pos, end = meso$pos),
                    logBF_per_ds_meso = meso$logBF_per_ds)
  
  Atlas_dt <- fr("fullres_0_8p0_0_65p1_atlas_general.rds")
  allLayersGR <- GRanges(seqnames = Atlas_dt$chr,
                         ranges = IRanges(start = Atlas_dt$pos, end = Atlas_dt$pos),
                         logBF_per_ds_allLayers = Atlas_dt$logBF_per_ds)
  
  endo6gpGR <- GRanges(seqnames = endo6gp$chr,
                    ranges = IRanges(start = endo6gp$pos, end = endo6gp$pos),
                    logBF_per_ds_endo6gp = endo6gp$logBF_per_ds)
  
  meso6gpGR <- GRanges(seqnames = meso6gp$chr,
                    ranges = IRanges(start = meso6gp$pos, end = meso6gp$pos),
                    logBF_per_ds_meso6gp = meso6gp$logBF_per_ds)
  
  ####################################################################
  ## Create a table with all CpG sites & score for each germ layer ##
  ###################################################################
  
  ## 1. Create union of all unique CpG positions
  table3layers <- union(union(union(union(allLayersGR, union(ectoGR, mesoGR)), endoGR), endo6gpGR), meso6gpGR)
  
  ## 2. Use findOverlaps to map logBF_per_ds values back
  # we want endoGR[i] -> table3layers[endoHits[[i]]] for each i, etc.
  
  endoHits <- findOverlaps(endoGR, table3layers, select = "first")
  ectoHits <- findOverlaps(ectoGR, table3layers, select = "first")
  mesoHits <- findOverlaps(mesoGR, table3layers, select = "first")
  allLayersHits <- findOverlaps(allLayersGR, table3layers, select = "first")
  endo6gpHits <- findOverlaps(endo6gpGR, table3layers, select = "first")
  meso6gpHits <- findOverlaps(meso6gpGR, table3layers, select = "first")
  
  # initialize columns with NA
  mcols(table3layers)$logBF_per_ds_endo <- NA_real_
  mcols(table3layers)$logBF_per_ds_ecto <- NA_real_
  mcols(table3layers)$logBF_per_ds_meso <- NA_real_
  mcols(table3layers)$logBF_per_ds_allLayers <- NA_real_
  mcols(table3layers)$logBF_per_ds_endo6gp <- NA_real_
  mcols(table3layers)$logBF_per_ds_meso6gp <- NA_real_
  
  # copy only the hits
  mcols(table3layers)$logBF_per_ds_endo[endoHits]   <- mcols(endoGR)$logBF_per_ds_endo
  mcols(table3layers)$logBF_per_ds_ecto[ectoHits]   <- mcols(ectoGR)$logBF_per_ds_ecto
  mcols(table3layers)$logBF_per_ds_meso[mesoHits]   <- mcols(mesoGR)$logBF_per_ds_meso
  mcols(table3layers)$logBF_per_ds_allLayers[allLayersHits]   <- mcols(allLayersGR)$logBF_per_ds_allLayers
  mcols(table3layers)$logBF_per_ds_endo6gp[endo6gpHits]   <- mcols(endo6gpGR)$logBF_per_ds_endo6gp
  mcols(table3layers)$logBF_per_ds_meso6gp[meso6gpHits]   <- mcols(meso6gpGR)$logBF_per_ds_meso6gp
  
  ## Add chr_pos column to identify positions
  table3layers$chr_pos <- paste0("chr", table3layers@seqnames,
                                 "_", table3layers@ranges@start)
  
  ################################################################################
  ### SAVED ###
  save(table3layers, file =
         here(paste0("gitignore/fulltable3layersPreMean_", format(Sys.Date(), "%d_%m_%y"), ".Rda")))
  
  # load(here("gitignore/fulltable3layersPreMean_26_08_26.Rda"))
  
  ################################################################################
  ## Arithmetic sum to stay in interpretable log-BF units
  m <- as.matrix(mcols(table3layers)[, c("logBF_per_ds_endo",
                                         "logBF_per_ds_ecto",
                                         "logBF_per_ds_meso")])
  
  table3layers$logBF_per_ds_mean3layers <- rowMeans(m, na.rm = FALSE)
  
  ### SAVED ###
  save(table3layers, file =
         here(paste0("gitignore/fulltable3layers_", format(Sys.Date(), "%d_%m_%y"), ".Rda")))
  
  # load(here(paste0("gitignore/fulltable3layers_26_08_26.Rda")))
  
  ################################################################################
  ## Plot the difference between the mean between layers and the logBF per ds 
  ## calculated on all datasets together
  ################################################################################
  
  if (!file.exists(here("B_MultiTissues/dataOut/figures/script04/mean3layers_vs_all.png"))){
    df <- as.data.frame(table3layers)
    
    set.seed(1234)
    p <- ggplot(df[sample(nrow(df), 100000),],
                aes(x = logBF_per_ds_mean3layers, y = logBF_per_ds_allLayers)) +
      geom_point(pch = 21, alpha = 0.1) +
      geom_abline(slope = 1) +
      theme_minimal(base_size = 14) +
      labs(title = "Hypervariability score using either all layers jointly or separately",
           x = "Average of the three layers hypervariability score",
           y = "Hypervariability score calculated on all cell types",
           subtitle = "(100k random CpG plotted)")
    
    ggplot2::ggsave(
      filename = here::here("B_MultiTissues/dataOut/figures/script04/mean3layers_vs_all.png"),
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
  #     0        1        2        3 
  # 103677   857187  1960134 20246679
  
  ## Select only rows covered in the 3 layers!
  table3layers_coveredIn3 <- table3layers[table3layers$n_layers %in% 3,]
  table(table3layers_coveredIn3$n_layers)
  ## 20.246.679
  
  ### SAVED ###
  save(table3layers_coveredIn3, file =
         here(paste0("gitignore/table3layers_coveredIn3_", format(Sys.Date(), "%d_%m_%y"), ".Rda")))
}

################################################################################
################################## CHECKPOINT ##################################
################################################################################
load(here(paste0("gitignore/table3layers_coveredIn3_26_08_26.Rda")))

################################################
## Get the top 99% quantile for logBF_per_ds ##
################################################

top99q <- quantile(table3layers_coveredIn3$logBF_per_ds_allLayers, probs = 0.99, na.rm = FALSE)
top99q_CpGs <- table3layers_coveredIn3[table3layers_coveredIn3$logBF_per_ds_allLayers >= top99q, ]$chr_pos

message(paste0("Total CpG sites: ", length(table3layers_coveredIn3)))
message(paste0("Total top99q CpG sites: ", length(top99q_CpGs), " (",
               round(length(top99q_CpGs)/length(table3layers_coveredIn3)*100,2), "% of total)"))
# Total CpG sites: 20246679
# Total top99q CpG sites: 202467 (1% of total)

if (!file.exists(here(paste0("gitignore/top99q_CpGs_", variant, ".RDS")))){
  saveRDS(top99q_CpGs, here(paste0("gitignore/top99q_CpGs_", variant, ".RDS")))  
}

## To use for testFetalSIV_ingp5.R

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

if (!file.exists(here("B_MultiTissues/dataOut/figures/script04/CompareWithResultsDerakhshan.png"))){
  
  ## Load array results
  if (!exists("resArray")) {
    resArray <- readRDS(here("B_MultiTissues/dataOut/resArray_0_8p0_0_65p1.RDS"))
  }
  
  if (!exists("resArray3ind")) {
    resArray3ind <- readRDS(here("B_MultiTissues/dataOut/resArray3ind_0_8p0_0_65p1.RDS"))
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
    ylab("Difference of hypervariability score shifted to start at zero")
  
  plotManhattan_Derakh <- plotManhattanFromdt(
    table3layers_coveredIn3dt[!chr %in% c("M", "Y"),], plotDerakhshan = TRUE, centro = centro) +
    ggtitle("Hypervariability score along the genome") +
    ylab("Hypervariability score")
  
  compPlotArrayAtlas <- makeCompPlot(
    X = resArray,
    Y = file.path(prep_dir, "fullres_0_8p0_0_65p1_atlas_general.rds"),
    whichX = "logBF_per_ds", whichY = "logBF_per_ds",
    title = "Derakhshan datasets vs Loyfer datasets",
    xlab = "Hypervariability score (Derakhshan datasets)",
    ylab = "Hypervariability score (Loyfer datasets)",
    minplot = 1e5)
  
  compPlotArrayAtlasRED <- makeCompPlot(
    X = resArray3ind,
    Y = file.path(prep_dir, "fullres_0_8p0_0_65p1_atlas_general.rds"),
    whichX = "logBF_per_ds", whichY = "logBF_per_ds",
    title = "Reduced Derakhshan datasets vs Loyfer datasets",
    xlab = "Hypervariability score (Derakhshan datasets reduced to 3 ind/ds)",
    ylab = "Hypervariability score (Loyfer datasets)",
    minplot = 1e5)
  
  row2 <- cowplot::plot_grid(pdiffhv_controls, compPlotArrayAtlas, compPlotArrayAtlasRED, 
                             labels = c("B", "C","D"), ncol = 3)
  
  ggsave(here("B_MultiTissues/dataOut/figures/script04/CompareWithResultsDerakhshan.png"),
         plot = cowplot::plot_grid(plotManhattan_Derakh, row2,
                                   nrow = 2, rel_heights = c(1,1.5), labels = c("A", "")),
         width = 20, height = 12, dpi = 300, bg = "white")
}

################################################################################
## Figure 3: MappingVariability                                               ##
################################################################################

####################################
## Plot Manhattan of logBF_per_ds ##

if (!file.exists(here(paste0("gitignore/plotManhattan_noDerakh_", variant, ".RDS")))){
  plotManhattan_noDerakh <- plotManhattanFromdt(
    table3layers_coveredIn3dt[!chr %in% c("M", "Y"),], plotDerakhshan = FALSE, centro = centro) +
    ggtitle("Hypervariability score along the genome")
  saveRDS(plotManhattan_noDerakh, here(paste0("gitignore/plotManhattan_noDerakh_", variant, ".RDS")))
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
  labs(x = "Chromosome", y = "Hypervariability score")+
  theme_minimal(base_size = 14)

###################################
## Y chromosome DNAm variability ##
###################################
ggplot() +
  geom_point(data = table3layers_coveredIn3dt[table3layers_coveredIn3dt$chr == "Y",], 
             aes(x = pos2, y = logBF_per_ds_allLayers),
             color = "black", size = 1, alpha = .5)+
  theme_classic() + theme(legend.position = "none") +
  labs(x = "Chromosome", y = "Hypervariability score")+
  theme_minimal(base_size = 14)

# logBF_per_ds in 3 regions

#######################################################
## Test enrichment of features for high logBF_per_ds ##
#######################################################

if (!file.exists(here(paste0("gitignore/pfeatures_", variant, ".RDS")))){
  
  # Create GRanges
  gr_cpg <- GRanges(
    seqnames = paste0("chr", table3layers_coveredIn3dt$chr),
    ranges = IRanges(start = table3layers_coveredIn3dt$pos, end = table3layers_coveredIn3dt$pos),
    logBF_per_ds = table3layers_coveredIn3dt$logBF_per_ds_allLayers
  )
  length(gr_cpg) # 20.246.679 CpGs
  
  # restrict to autosomes and chr X
  gr_cpg <- gr_cpg[gr_cpg@seqnames %in% c(paste0("chr", 1:22), "chrX"),]
  length(gr_cpg) # 20.241.760 CpGs
  
  # Import bed file
  bed_features <- genomation::readTranscriptFeatures(here("gitignore/hg38_GENCODE_V47.bed"))
  
  # Annotate CpGs and see which regions have higher logBF_per_ds_allLayers (takes long)
  anno_result <- genomation::annotateWithGeneParts(
    target = gr_cpg, feature = bed_features)
  
  saveRDS(anno_result, here("gitignore/anno_result.RDS"))
  
  anno_result@perc.of.OlapFeat
  # promoter     exon   intron 
  # 93.03910 70.03343 87.28536 
  
  ## Add info from annotation to our GRange object
  gr_cpg$featureType <- ifelse(anno_result@members[, "prom"] == 1, "promoter",
                               ifelse(anno_result@members[, "exon"] == 1, "exon",
                                      ifelse(anno_result@members[, "intron"] == 1, "intron", "intergenic")))
  
  mcols(gr_cpg) %>% as.data.frame() %>%
    dplyr::group_by(featureType) %>%
    dplyr::summarise(meanlogBF_per_ds = mean(logBF_per_ds, na.rm=T),
                     medianlogBF_per_ds = median(logBF_per_ds, na.rm=T))
  # featureType    meanlogBF_per_ds medianlogBF_per_ds
  # 1 exon                  -0.313             -0.408
  # 2 intergenic            -0.184             -0.244
  # 3 intron                -0.291             -0.383
  # 4 promoter              -0.343             -0.459
  
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
    labs(y = "% of CpGs in feature", x = NULL, fill = "Hypervariability score\n(percentile)") +
    scale_fill_viridis_d() +
    theme_minimal(base_size = 14) +
    geom_hline(yintercept = 10, linetype = 3, linewidth = 0.3) +
    ggtitle("Distribution of CpG hypervariability across feature types")
  pfeatures
  
  saveRDS(pfeatures, here(paste0("gitignore/pfeatures_", variant, ".RDS")))
}

# Distribution of CpG hypervariability across feature types. CpGs were binned into quantiles 
# of the per-dataset log Bayes factor (logBF per ds) computed genome-wide, and for each
# feature type the percentage of its CpGs falling in each upper quantile is shown 
# (the bottom quintile is omitted; percentages are of all CpGs in the feature). 
# Bars exceeding 10% (the expected share under no enrichment) indicate feature types enrichement.

#######################
## Enrichement in TE ##
#######################

if (!file.exists(here(paste0("gitignore/TEplot_", variant, ".RDS")))){
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
  # top99q_CpGs   112737     89730
  # background  10669103   9375109
  # 
  # $pvalue
  # [1] 1.061446e-107
  # 
  # $odds_ratio
  # [1] 1.103994
  # 
  # $conf_low
  # [1] 1.094269
  # 
  # $conf_high
  # [1] 1.113781
  
  ## ---- Per repClass ----
  te_by_class <- split(te_regions, mcols(te_regions)$repClass)
  class_res <- lapply(names(te_by_class), function(cl)
    fisher_test_te(te_by_class[[cl]], top99q_CpGs_GR, bg_only_GR,
                   label = cl, nameTarget = "top99q_CpGs"))
  class_dt <- data.table::rbindlist(lapply(class_res, `[`, res_cols))
  class_dt[, p.adj := p.adjust(pvalue, "BH")]
  class_dt[order(-odds_ratio)]
  #         label        pvalue odds_ratio  conf_low conf_high         p.adj
  # 1:       LINE  0.000000e+00  1.5634939 1.5459085 1.5812874  0.000000e+00
  # 2:        LTR 7.486682e-282  1.3169069 1.2978621 1.3362254 1.497336e-281
  # 3:        DNA  6.510516e-25  1.1447071 1.1160691 1.1739134  9.765774e-25
  # 4: Retroposon  2.377514e-01  1.0827551 0.9457392 1.2342082  2.853017e-01
  # 5:         RC  9.039149e-01  1.0294752 0.6086499 1.6319471  9.039149e-01
  # 6:       SINE  0.000000e+00  0.7350373 0.7274819 0.7426105  0.000000e+00
  
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
  
  #            label        pvalue odds_ratio  conf_low conf_high         p.adj
  # 1:            L1  0.000000e+00  1.6895785 1.6690312 1.7103443  0.000000e+00
  # 2:          ERVK  1.992181e-20  1.3537627 1.2722759 1.4391576  9.960903e-20
  # 3:     ERVL-MaLR 1.488889e-128  1.3401095 1.3097615 1.3709956 1.488889e-127
  # 4:     MULE-MuDR  5.065536e-02  1.3340302 0.9914094 1.7577955  1.013107e-01
  # 5:          ERVL  4.344403e-56  1.3047879 1.2637651 1.3469158  2.606642e-55
  # 6:           hAT  3.944642e-02  1.2987802 1.0007350 1.6584729  8.452803e-02
  # 7: TcMar-Mariner  1.941548e-04  1.2943674 1.1300442 1.4760304  5.295130e-04
  # 8:          ERV1  2.168918e-70  1.2436195 1.2147891 1.2729879  1.626688e-69
  # 9:    hAT-Tip100  3.246848e-06  1.2304593 1.1286943 1.3390571  9.740544e-06
  # 10:           LTR  2.596615e-01  1.2171044 0.8460100 1.6968093  3.709449e-01
  # 11:         Gypsy  2.601682e-02  1.1731130 1.0156331 1.3482781  6.504205e-02
  # 12: hAT-Blackjack  3.219628e-02  1.1705158 1.0110940 1.3481605  7.429911e-02
  # 13:  TcMar-Tigger  5.071898e-10  1.1470843 1.0989709 1.1967884  2.173670e-09
  # 14:   hAT-Charlie  1.319619e-08  1.1152870 1.0743518 1.1573815  4.948571e-08
  # 15:         RTE-X  2.524237e-01  1.0920580 0.9309300 1.2732332  3.709449e-01
  # 16:     TcMar-Tc2  4.467066e-01  1.0854058 0.8596229 1.3527137  5.413684e-01
  # 17:           SVA  2.377514e-01  1.0827551 0.9457392 1.2342082  3.709449e-01
  # 18:            L2  5.969915e-07  1.0699628 1.0420265 1.0985182  1.989972e-06
  # 19:           CR1  1.986627e-01  1.0552915 0.9707259 1.1452053  3.317800e-01
  # 20:      Helitron  9.039149e-01  1.0294752 0.6086499 1.6319471  1.000000e+00
  # 21:           MIR  3.158246e-01  0.9848103 0.9559255 1.0143601  4.306700e-01
  # 22:     5S-Deu-L2  1.000000e+00  0.9611612 0.5366395 1.5903342  1.000000e+00
  # 23:      tRNA-RTE  8.716565e-01  0.9519133 0.6691719 1.3145428  1.000000e+00
  # 24:          tRNA  1.000000e+00  0.9476658 0.5290862 1.5679598  1.000000e+00
  # 25:      Penelope  1.000000e+00  0.9209277 0.2497598 2.3751213  1.000000e+00
  # 26:      PiggyBac  4.511404e-01  0.8565219 0.5810600 1.2180101  5.413684e-01
  # 27:      RTE-BovB  1.587532e-01  0.8293667 0.6337110 1.0666768  2.976622e-01
  # 28:        hAT-Ac  3.852743e-01  0.7919817 0.4604471 1.2711185  5.025317e-01
  # 29:           Alu  0.000000e+00  0.7227254 0.7150546 0.7304433  0.000000e+00
  # 30:           DNA  1.990680e-01  0.6453533 0.3088561 1.1901137  3.317800e-01
  
  ## ---- Plot: OR with 95% CI, ordered, coloured by significance ----
  res_plot <- res_dt[order(odds_ratio)]
  res_plot[, label := factor(label, levels = label)]  # lock the OR order
  res_plot[, sig := ifelse(p.adj < 0.05, "FDR < 0.05", "n.s.")]
  
  TEplot <- ggplot(res_plot, aes(odds_ratio, label, colour = sig)) +
    geom_vline(xintercept = 1, linetype = 3) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), height = 0.25, orientation = "y") +
    geom_point(size = 3) +
    scale_x_log10(breaks = c(0,.5, 1, 1.5,2)) +                                 # OR is multiplicative -> log axis
    scale_colour_manual(values = c("FDR < 0.05" = "#DC3220", "n.s." = "grey60")) +
    labs(x = "Odds ratio (top99q hvCpG vs background, log scale)",
         y = NULL, colour = NULL,
         title = "TE family enrichment among 1% top hypervariable CpGs") +
    theme_minimal(base_size = 12)
  
  saveRDS(TEplot, here(paste0("gitignore/TEplot_", variant, ".RDS")))
}

#################################
## plot manhattan, feature, TE ##
#################################

if (!file.exists(here("B_MultiTissues/dataOut/figures/script04/MappingVariability.png"))){
  plotManhattan_noDerakh <- readRDS(here(paste0("gitignore/plotManhattan_noDerakh_", variant, ".RDS")))
  pfeatures <- readRDS(here(paste0("gitignore/pfeatures_", variant, ".RDS")))
  TEplot <- readRDS(here(paste0("gitignore/TEplot_", variant, ".RDS")))
  
  ggplot2::ggsave(
    filename = here::here(
      "B_MultiTissues/dataOut/figures/script04/MappingVariability.png"),
    plot = plot_grid(plotManhattan_noDerakh + ylab("Hypervariability score"),
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

if (!file.exists(here("B_MultiTissues/dataOut/figures/script04/CompareWithpreviousMEs.png"))){
  
  ## Use the GR object with analyses in the 3 layers
  
  # Fix chromosome names in geomMeanGR (1 -> chr1)
  if (sum(grepl("chr", seqlevels(table3layers_coveredIn3))) == 0){
    seqlevels(table3layers_coveredIn3) <- paste0("chr", seqlevels(table3layers_coveredIn3))
  }
  
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
    "top99q"             = RColorBrewer::brewer.pal(8, "Set2")[6],
    "VanBaakSIV"         = RColorBrewer::brewer.pal(8, "Set2")[7]
  )
  
  # ── shared spacing refinements ───────────────────────────────────────────────
  # Added AFTER each theme_minimal()/theme_classic() so it is not overwritten.
  # Top margin reserves room for the cowplot panel labels; axis-title margins
  # push titles off the tick text.
  spacing <- theme(
    plot.margin  = margin(t = 24, r = 10, b = 12, l = 10),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  )
  
  # ── ME overlap ────────────────────────────────────────────
  MEsetdt <- make_MEsetdt(sets, GR = table3layers_coveredIn3)
  
  MEsetdt <- na.omit(MEsetdt)
  nrow(MEsetdt) ## 52.244
  
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
  
  # after fit/emmeans/contrasts are computed on the original MEsetdt,
  # reorder ME by mean score FOR THE PLOT
  MEsetdt_plot <- copy(MEsetdt)
  MEsetdt_plot[, ME := forcats::fct_reorder(ME, logBF_per_ds, .fun = median, na.rm = TRUE)]
  
  # recompute N labels on the reordered data (levels now in mean order)
  n_labels <- MEsetdt_plot[, .(n = .N), by = ME]
  y_top    <- max(MEsetdt_plot$logBF_per_ds, na.rm = TRUE)
  
  pMElogBF_per_ds <- ggplot(MEsetdt_plot, aes(x = ME, y = logBF_per_ds)) +
    geom_jitter(aes(colour = ME), size = 3, alpha = .2) +
    geom_violin(aes(colour = ME)) +
    geom_boxplot(aes(colour = ME), width = .1) +
    geom_text(data = n_labels,
              aes(x = ME, y = y_top, label = format(n, big.mark = ",")),
              vjust = -0.4, size = 5.5, fontface = "bold", inherit.aes = FALSE) +
    scale_colour_manual(values = group_cols, name = "CpG set") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    ylab("Hypervariability score")
  
  contrasts_plot <- contrasts %>%
    mutate(ME_name = forcats::fct_reorder(ME_name, estimate, .desc = TRUE))
  
  pcontrast <- ggplot(contrasts_plot, aes(x = ME_name, y = estimate, colour = ME_name)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_colour_manual(values = group_cols, name = "CpG set") +
    coord_flip() +
    labs(y = "Difference in hypervariability score vs mQTLcontrols", x = "") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ## pdecay - legend inside the plot, bottom-left corner
  pdecay <- plot_decay_curve(MEsetdt) +
    scale_colour_manual(values = group_cols, name = "CpG set")
  
  # ── Save key objects for S07 ──────────────────────────────────────────────────
  saveRDS(MEsetdt, here(paste0("gitignore/MEsetdt_", variant, ".rds")))
  
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
  # saveRDS(top99q_CpGs, file = here(paste0("gitignore/top99q_CpGs_", variant, ".RDS"))
  
  ## Then run on ingp5: testFetalSIV_ingp5.R
  
  plots <- readRDS(here("gitignore/intercorrelationSIVfetal_sepSIV.rds"))
  
  # per-group N for panel D
  nD <- as.data.table(plots$interlayer_corr)[, .(n = .N), by = group]
  yD_top <- max(plots$interlayer_corr$interlayer_r, na.rm = TRUE)
  
  pinterlayer_corr <- ggplot(plots$interlayer_corr,
                             aes(x = group, y = interlayer_r, group = group, fill = group)) +
    geom_violin(width = 1.4) +
    geom_boxplot(width = 0.3, fill = "white") +
    geom_text(data = nD,
              aes(x = group, y = yD_top, label = format(n, big.mark = ",")),
              vjust = -0.4, size = 5.5, fontface = "bold", inherit.aes = FALSE) +
    scale_fill_manual(values = group_cols) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
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
      pMElogBF_per_ds + spacing + theme(axis.title.x = element_blank()),
      pcontrast + theme_minimal(base_size = 18) + spacing +
        theme(plot.title = element_text(size=16), legend.position = "none"),
      nrow = 2,
      labels = c("A. Distribution of the hypervariability score for each CpG set",
                 "B. Comparison of previous CPG sets groups to mQTLcontrols"),
      label_size = 16, label_x = 0, hjust = 0, label_y = 0.98, vjust = 1),
    pdecay + theme_minimal(base_size = 20) + spacing + 
      theme(legend.position = "inside",
            legend.position.inside = c(.7, .6),
            legend.justification = c(0, 0),
            legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3)),
    ncol = 2,
    rel_widths = c(1, 1),
    labels = c("", "C. Decay curve of hypervariability score per percentile"),
    label_size = 16, label_x = 0, hjust = 0, label_y = 0.98, vjust = 1
  )
  
  SIV_plot <- plot_grid(
    pinterlayer_corr + theme_minimal(base_size = 20) + spacing +
      theme(axis.text.x = element_text(angle = 20, hjust = 1),
            axis.title.x = element_blank(), legend.position = "none"),
    plot_grid(pinterindividual_var + theme_minimal(base_size = 14) + spacing +
                labs(fill = "CpG set"),
              pbinned + theme_minimal(base_size = 16) + spacing +
                theme(axis.text.x = element_text(angle = 20, hjust = 1))+
                labs(colour = "CpG set", fill = "CpG set"),
              ncol = 1, align = "v",
              labels = c("E. Densities of interindividual variation per CpG within the fetal data, by set",
                         "F. Inter-germ-layer correlation per interindividual variation, binned"),
              label_size = 16, label_x = 0, hjust = 0, label_y = 0.98, vjust = 1),
    ncol = 2,
    rel_widths = c(1, 1),
    labels = c("D. Mean inter-germ-layer correlation for each CpG set", ""),
    label_size = 16, label_x = 0, hjust = 0, label_y = 0.98, vjust = 1
  )
  
  final_plot <- plot_grid(
    upperRow,
    SIV_plot,
    ncol = 1,
    rel_heights = c(1, 1)
  ) + theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
  
  ggplot2::ggsave(
    filename = here::here(
      "B_MultiTissues/dataOut/figures/script04/CompareWithpreviousMEs.png"),
    plot = final_plot, width = 24, height = 20,  dpi = 300, bg = "white")
}

############################################################
## How many of each putative ME is actually in the top99q? ##
############################################################

# Universe of covered CpGs, each already labelled top99q vs not
# (both are single-CpG GRanges built from your covered-in-3 sites)
top_gr  <- listGR$top99q          # top 1% CpGs
rest_gr <- listGR$allButTop99q    # the other 99%

# For each ME set, find which COVERED CpGs fall inside that set's regions,
# then classify those CpGs as top99q or not.
me_sets <- split(putativeME_GR, putativeME_GR$set)

summary_df <- rbindlist(lapply(names(me_sets), function(s) {
  regions <- reduce(me_sets[[s]], ignore.strand = TRUE)   # collapse overlapping regions
  
  # covered CpGs inside this set's regions
  top_in  <- sum(overlapsAny(top_gr,  regions, ignore.strand = TRUE))
  rest_in <- sum(overlapsAny(rest_gr, regions, ignore.strand = TRUE))
  n_cov   <- top_in + rest_in                              # total covered CpGs in the set
  
  data.table(
    set           = s,
    n_cpg_covered = n_cov,                                 # covered CpGs the set overlaps
    n_top99q      = top_in,
    pc_top99q     = if (n_cov) 100 * top_in  / n_cov else NA,
    n_rest        = rest_in,
    pc_rest       = if (n_cov) 100 * rest_in / n_cov else NA
  )
}))

summary_df$set <- as.character(summary_df$set)
summary_df$set[summary_df$set %in% "Gunasekara corSIV"] <- "Gunasekara CoRSIV"

## add fold-enrichment vs the genome-wide top99q rate (baseline ≈ 1%)
baseline <- length(top_gr) / (length(top_gr) + length(rest_gr))   # ~0.01
summary_df[, fold := (pc_top99q / 100) / baseline]

## Format pretty
summary_df %>%
  mutate(
    across(starts_with("pc_"), ~ scales::percent(.x / 100, accuracy = 0.1)),
    fold = scales::number(fold, accuracy = 0.1, suffix = "×")
  ) %>%
  gt() %>%
  fmt_number(columns = starts_with("n_"), decimals = 0) %>%
  cols_label(
    set           = "Previously published CpG set",
    n_cpg_covered = "Covered CpGs in set",
    n_top99q      = "N in top 1% hvCpGs",
    pc_top99q     = "%",
    n_rest        = "N in rest (non-top 1% hvCpGs)",
    pc_rest       = "%",
    fold          = "Fold enrichment"
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

getEnrichCentroTelo()
# ── Centromere ──
# hvCpGs in region:      3541 / 202467 (1.75%)
# Background in region:  187301 / 20246679 (0.93%)
# Fold enrichment:       1.88
# Fisher p (one-sided):  2.290384e-257
# Odds ratio:            1.91
# 
# ── Subtelomere (1Mb) ──
# hvCpGs in region:      3841 / 202467 (1.9%)
# Background in region:  417898 / 20246679 (2.06%)
# Fold enrichment:       0.92
# Fisher p (one-sided):  1
# Odds ratio:            0.92

#######################
## Test GO of top 1% ##
#######################

totalSites <- table3layers_coveredIn3$chr_pos
if (!exists("top99q_CpGs")){
  top99q <- quantile(table3layers_coveredIn3$logBF_per_ds_allLayers, probs = 0.99, na.rm = FALSE)
  top99q_CpGs <- table3layers_coveredIn3[table3layers_coveredIn3$logBF_per_ds_allLayers >= top99q, ]$chr_pos
}

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
## "Gene universe contains 32326 genes"

fg_bypass  <- unique(cpg_cluster_count_per_gene(top99q_CpGs)$entrez_id)
fg_wrapper <- annotateCpGs_txdb(clusterCpGs(top99q_CpGs, 50, 2), 10000)
length(intersect(fg_bypass, fg_wrapper)) / length(union(fg_bypass, fg_wrapper))  # want ≈ 1
# 1 
## The Jaccard is 1.0 --> the cpg_cluster_count_per_gene counts genes under exactly 
# the same body∪promoter rule as the annotation wrapper. The CpG-count covariate is valid.

# WGBS-correct control
res_cpg <- CpG_GO_pipeline_lengthControlled(
  top99q_CpGs, universe = universe,
  control_method = "cpg_count", all_sites = totalSites)

# Clustering CpGs...
# Reduced from 202467 to 9321 clustered CpGs
# Annotating genes...
# Found 1977 Entrez genes
# Controlling for CpG-cluster count per gene (WGBS)...
# Median clustered CpGs/gene — foreground: 1204.0, universe: 327.0, ratio: 3.68
# Matched universe: 10427 genes (was 32326)
# Foreground genes in matched universe: 1977 / 1977 (100.0%)
# Running GO enrichment...

# for the sensitivity table, also run the other two
res_len  <- CpG_GO_pipeline_lengthControlled(top99q_CpGs, universe = universe,
                                             control_method = "length")
# Clustering CpGs...
# Reduced from 202467 to 9321 clustered CpGs
# Annotating genes...
# Found 1977 Entrez genes
# Controlling for gene length (bp)...
# Median gene length — foreground: 5,434 bp, universe: 3,137 bp, ratio: 1.73
# Matched universe: 11862 genes (was 32326)
# Foreground genes in matched universe: 1977 / 1977 (100.0%)
# Running GO enrichment...

res_none <- CpG_GO_pipeline_lengthControlled(top99q_CpGs, universe = universe,
                                             control_method = "none")
# Clustering CpGs...
# Reduced from 202467 to 9321 clustered CpGs
# Annotating genes...
# Found 1977 Entrez genes
# Matched universe: 32326 genes (was 32326)
# Foreground genes in matched universe: 1977 / 1977 (100.0%)
# Running GO enrichment...

sapply(list(none = res_none, length = res_len, cpg = res_cpg),
       function(r) sum(r$BP@result$p.adjust < 0.05, na.rm = TRUE))
# none length    cpg 
# 17      3      3 

lapply(list(none = res_none, length = res_len, cpg = res_cpg),
       function(r) r$BP@result %>% filter(p.adjust < 0.05) %>% pull(Description) %>% head(20))

# $none
# [1] "homophilic cell-cell adhesion"                          
# [2] "synapse assembly"                                       
# [3] "cell junction assembly"                                 
# [4] "negative chemotaxis"                                    
# [5] "axon guidance"                                          
# [6] "neuron projection guidance"                             
# [7] "axonogenesis"                                           
# [8] "neuron recognition"                                     
# [9] "cell morphogenesis involved in neuron differentiation"  
# [10] "axon development"                                       
# [11] "regulation of small GTPase mediated signal transduction"
# [12] "regulation of postsynaptic membrane potential"          
# [13] "regulation of membrane potential"                       
# [14] "negative regulation of cAMP/PKA signal transduction"    
# [15] "regulation of axonogenesis"                             
# [16] "regulation of synapse organization"                     
# [17] "regulation of cytosolic calcium ion concentration"      
# 
# $length
# [1] "cell junction organization"    "homophilic cell-cell adhesion" "synapse organization"         
# 
# $cpg
# [1] "homophilic cell-cell adhesion" "cell junction organization"    "negative chemotaxis"   

## --> Likely, the broad neurodevelopmental signature is largely a CpG-density artefact, not an ME signal.

## cpg_count matching can under-power the protocadherin signal specifically,
# because the PCDH clusters are so CpG-dense that few comparable genes exist to match them
# so if adhesion terms weaken here, that's not proof they're artefactual.

## Direct test of prothocadherin

pcdh <- GRanges("chr5", IRanges(140710000, 141510000))     # hg38 PCDH clusters
top_gr <- makeGRfromMyCpGPos(top99q_CpGs, "top99q")
bg_gr  <- makeGRfromMyCpGPos(totalSites,  "bg")
obs <- sum(overlapsAny(top_gr, pcdh))
exp <- length(top_gr) * mean(overlapsAny(bg_gr, pcdh))
c(observed = obs, expected = exp, fold = obs / exp)
# observed   expected       fold 
# 161.000000  99.380103   1.620043 

in_top <- sum(overlapsAny(top_gr, pcdh))
in_bg  <- sum(overlapsAny(bg_gr,  pcdh))
mat <- matrix(c(in_top, length(top_gr) - in_top,
                in_bg,  length(bg_gr)  - in_bg), nrow = 2, byrow = TRUE)
fisher.test(mat, alternative = "greater")   # p-value + CI for the 1.34x
# Fisher's Exact Test for Count Data
# data:  mat
# p-value = 1.029e-08
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#   1.41476     Inf
# sample estimates:
#   odds ratio 
# 1.620518 

# --> top99q CpGs are 1.62× enriched at the PCDH clusters (161 observed vs 99 expected)
# a naive GO enrichment shows a broad neurodevelopmental signature (17 BP terms),
# but this is largely attributable to the higher CpG density of the foreground genes
# (3.7× the universe); after matching the background on per-gene CpG count,
# only homophilic cell-cell adhesion, cell junction organization, and
# negative chemotaxis remain, likely driven by enrichment at the clustered
# protocadherin locus on chr5 (Fisher's exact test OR = 1.62, p < 0.0001), a
# a known systemically-variable ME region.

## ── Leave-one-out: does removing the protocadherin cluster drop the terms? ──

# PCDH-cluster entrez IDs (hg38 chr5 ~140.7–141.5 Mb)
pcdh_region <- GRanges("chr5", IRanges(140710000, 141510000))
genes_gr <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
pcdh_entrez <- genes_gr$gene_id[overlapsAny(genes_gr, pcdh_region, ignore.strand = TRUE)]

res_cpg_noPCDH <- CpG_GO_pipeline_lengthControlled(
  top99q_CpGs, universe = universe,
  control_method = "cpg_count", all_sites = totalSites,
  exclude_genes = pcdh_entrez)

# Clustering CpGs...
# Reduced from 202467 to 9321 clustered CpGs
# Annotating genes...
# Found 1977 Entrez genes
# Excluded 12 foreground genes (e.g. PCDH cluster); 1977 -> 1965
# Controlling for CpG-cluster count per gene (WGBS)...

# terms lost by removing PCDH = attributable to the cluster
bp_with    <- res_cpg$BP@result        %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::pull(Description)
bp_without <- res_cpg_noPCDH$BP@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::pull(Description)
setdiff(bp_with, bp_without)
# [1] "homophilic cell-cell adhesion" "negative chemotaxis"     

bp_without
# [1] "cell junction organization"

res_cpg$BP@result %>% filter(Description == "negative chemotaxis") %>% pull(geneID)
# ROBO2/UNC5C/SEMA3E/ROBO1/SEMA5B/NRG3/SEMA5A/EFNA5/SEMA6A/SLIT2/SEMA6D/ITGB3/NTN1/SLIT3

res_cpg_noPCDH$BP@result  %>% dplyr::filter(p.adjust < 0.05)
# ID                Description GeneRatio  BgRatio RichFactor FoldEnrichment   zScore
# GO:0034330 GO:0034330 cell junction organization  104/1044 411/6133  0.2530414       1.486497 4.624471
# pvalue   p.adjust     qvalue
# GO:0034330 6.874005e-06 0.02492514 0.02492514
# geneID
# GO:0034330 NOS1AP/BMP6/NEURL1/ROBO2/ITSN1/CNTN2/LZTS1/EPHB2/NRXN3/SLITRK2/CDH22/PECAM1/LRP4/IQSEC1/KIRREL1/ADGRL2/EGFLAM/STRN/PPFIA3/TMEM108/CDH13/SHANK1/SEMA3E/CHRNA1/DRD2/SORBS1/FN1/ITGA6/CNTN4/DOCK4/PKP2/LAMC1/ERBB4/CNTN5/ARHGAP4/GABRA5/NRG3/ABCC8/CSMD2/DUSP22/GABRA1/HAPLN4/GPR158/IL1RAPL2/EFNA5/ASIC2/SDK1/UBE3B/GABRB3/PLXNC1/THBS2/GRID1/PKP1/KCNK13/ILDR1/IL1RAPL1/KLK8/NLGN1/EPHB1/MTMR2/ERC2/SRPX2/LIMS1/CDH5/MTSS1/GAP43/ARHGAP6/RDX/TENM3/DOCK1/PRMT8/BCR/CLDN11/TENM2/ITGB3/GRID2/DNER/GABRG1/SHISA6/CACNB2/CDKL5/SYN3/MAPRE2/LAMA1/PRICKLE1/PEAK1/OPCML/PIP5K1A/TJP2/PAK2/LRRC4C/EPHB3/HMCN2/NTN1/EGLN1/SLC8A2/LARGE1/PDLIM5/SHANK2/PKHD1/APBB2/CDH8/TENM4/CACNA2D3
# Count
# GO:0034330   104

res_cpg_noPCDH$CC@result  %>% dplyr::filter(p.adjust < 0.05)
# ID           Description GeneRatio  BgRatio RichFactor FoldEnrichment   zScore       pvalue
# GO:0097060 GO:0097060     synaptic membrane   65/1089 237/6426  0.2742616       1.618370 4.381414 2.583329e-05
# GO:0009986 GO:0009986          cell surface   78/1089 306/6426  0.2549020       1.504132 4.081609 6.848117e-05
# GO:0045211 GO:0045211 postsynaptic membrane   47/1089 165/6426  0.2848485       1.680841 4.001917 1.269869e-04
# p.adjust     qvalue
# GO:0097060 0.01304581 0.01304581
# GO:0009986 0.01729150 0.01729150
# GO:0045211 0.02137612 0.02137612
# geneID
# GO:0097060                                                                                  CNTN2/CNIH3/EPHB2/PCDH9/NRXN3/SLC1A2/SLITRK2/LRP4/IQSEC1/ADGRL2/STRN/GRM1/TMEM108/SHANK1/GRM7/SYNE1/CHRNA1/DRD2/ERBB4/PICALM/GRIP1/CNTN5/GABRA5/ABCC8/CSMD2/GABRA1/SORCS2/GPR158/SCN10A/ASIC2/PRKCG/UNC13C/GABRB3/GRID1/NCAM2/CLTA/IL1RAPL1/CHRM5/NLGN1/SLC6A6/MTMR2/ERC2/SRPX2/DMD/HCN1/NETO1/KCNA1/LRRC7/CLMP/TENM2/ITGB3/GRID2/GABRG1/SHISA6/KCNA2/RGS7/LRRC4C/SLC1A7/DGKB/SHANK2/GSG1L/CDH8/CACNG5/CACNA2D3/GRIN2A
# GO:0009986 GPC6/GPC5/ROBO2/GFRA1/PTPRT/UNC5C/CNTN2/EPHB2/TSPAN8/AJAP1/SLC1A2/ADAMTS7/ANTXRL/PECAM1/SLAMF8/LRP4/TRPC4/CLEC4C/NRROS/SIRPA/CDH13/ENTPD1/CD200R1L/CHRNA1/PCSK6/ITGA6/IGSF5/HAVCR2/MYO18A/CR1/ROBO1/PICALM/ATP5PO/FCGR3A/PKD1L3/DSCAML1/TPO/COL23A1/MUC17/TFRC/IL1RAPL2/EFNA5/GABRB3/P2RX7/TMC1/IL1RAPL1/NLGN1/LMO7/SLIT2/KCNN2/SRPX2/DMD/ITGB8/CDH5/HCN1/SEMA6D/KCNA1/TMX3/TF/SLC46A2/CLMP/ITGB3/ENOX2/MELTF/PDGFRA/STAB2/ANXA4/GFRA2/TSPEAR/SCNN1A/ZPLD1/PKHD1/TMEM8B/FCGR2A/ABCG1/BTNL8/ASTN1/GRIN2A
# GO:0045211                                                                                                                                                                                                 CNTN2/CNIH3/EPHB2/SLITRK2/IQSEC1/ADGRL2/STRN/GRM1/TMEM108/SHANK1/GRM7/SYNE1/CHRNA1/DRD2/ERBB4/PICALM/GRIP1/GABRA5/CSMD2/GABRA1/SORCS2/GPR158/ASIC2/GABRB3/GRID1/IL1RAPL1/CHRM5/NLGN1/SLC6A6/DMD/HCN1/NETO1/LRRC7/CLMP/TENM2/ITGB3/GRID2/GABRG1/SHISA6/KCNA2/RGS7/LRRC4C/DGKB/SHANK2/GSG1L/CACNG5/GRIN2A
# Count
# GO:0097060    65
# GO:0009986    78
# GO:0045211    47

res_cpg_noPCDH$MF@result  %>% dplyr::filter(p.adjust < 0.05)
# ID                                Description GeneRatio  BgRatio RichFactor FoldEnrichment
# GO:0005216 GO:0005216            monoatomic ion channel activity   57/1081 201/6345  0.2835821       1.664504
# GO:0022836 GO:0022836                     gated channel activity   43/1081 148/6345  0.2905405       1.705347
# GO:0015267 GO:0015267                           channel activity   57/1081 213/6345  0.2676056       1.570729
# GO:0022803 GO:0022803 passive transmembrane transporter activity   57/1081 213/6345  0.2676056       1.570729
# zScore       pvalue   p.adjust     qvalue
# GO:0005216 4.338167 3.357596e-05 0.02350317 0.02350317
# GO:0022836 3.934406 1.690876e-04 0.03369974 0.03369974
# GO:0015267 3.839326 1.925700e-04 0.03369974 0.03369974
# GO:0022803 3.839326 1.925700e-04 0.03369974 0.03369974
# geneID
# GO:0005216 TPTE/PIEZO2/CLIC5/TRPM3/CACNA1H/OTOP1/TRPC4/RYR3/CLCN1/TPTE2/RYR2/KCNIP1/ANO2/CHRNA1/KCNJ5/KCNG2/ATP5PO/CACNA1F/GABRA5/PKD1L3/ABCC8/GABRA1/TRPC4AP/SCN10A/ASIC2/GABRB3/KCNK3/P2RX7/GRID1/ANO5/TMC1/KCNK13/TMEM175/FXYD3/ANO3/KCNN2/HCN1/CLCN5/ANO7/KCNA1/ITPR1/CACNA1G/RYR1/GRID2/GABRG1/KCNMB2/CLIC2/CACNB2/KCNA2/TMEM63A/TRPC1/SLC1A7/SCNN1A/CATSPER3/CACNG5/CACNA2D3/GRIN2A
# GO:0022836                                                                                         PIEZO2/CACNA1H/RYR3/CLCN1/RYR2/ANO2/CHRNA1/KCNJ5/KCNG2/CACNA1F/GABRA5/PKD1L3/ABCC8/GABRA1/SCN10A/ASIC2/GABRB3/KCNK3/P2RX7/GRID1/ANO5/TMC1/KCNK13/KCNN2/HCN1/CLCN5/ANO7/KCNA1/ITPR1/CACNA1G/RYR1/GRID2/GABRG1/KCNMB2/CACNB2/KCNA2/TMEM63A/SLC1A7/SCNN1A/CATSPER3/CACNG5/CACNA2D3/GRIN2A
# GO:0015267 TPTE/PIEZO2/CLIC5/TRPM3/CACNA1H/OTOP1/TRPC4/RYR3/CLCN1/TPTE2/RYR2/KCNIP1/ANO2/CHRNA1/KCNJ5/KCNG2/ATP5PO/CACNA1F/GABRA5/PKD1L3/ABCC8/GABRA1/TRPC4AP/SCN10A/ASIC2/GABRB3/KCNK3/P2RX7/GRID1/ANO5/TMC1/KCNK13/TMEM175/FXYD3/ANO3/KCNN2/HCN1/CLCN5/ANO7/KCNA1/ITPR1/CACNA1G/RYR1/GRID2/GABRG1/KCNMB2/CLIC2/CACNB2/KCNA2/TMEM63A/TRPC1/SLC1A7/SCNN1A/CATSPER3/CACNG5/CACNA2D3/GRIN2A
# GO:0022803 TPTE/PIEZO2/CLIC5/TRPM3/CACNA1H/OTOP1/TRPC4/RYR3/CLCN1/TPTE2/RYR2/KCNIP1/ANO2/CHRNA1/KCNJ5/KCNG2/ATP5PO/CACNA1F/GABRA5/PKD1L3/ABCC8/GABRA1/TRPC4AP/SCN10A/ASIC2/GABRB3/KCNK3/P2RX7/GRID1/ANO5/TMC1/KCNK13/TMEM175/FXYD3/ANO3/KCNN2/HCN1/CLCN5/ANO7/KCNA1/ITPR1/CACNA1G/RYR1/GRID2/GABRG1/KCNMB2/CLIC2/CACNB2/KCNA2/TMEM63A/TRPC1/SLC1A7/SCNN1A/CATSPER3/CACNG5/CACNA2D3/GRIN2A
# Count
# GO:0005216    57
# GO:0022836    43
# GO:0015267    57
# GO:0022803    57