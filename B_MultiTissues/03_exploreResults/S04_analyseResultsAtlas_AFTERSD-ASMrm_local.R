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
  
  # Fix chromosome names in geomMeanGR (1 -> chr1)
  if (sum(grepl("chr", seqlevels(table3layers_coveredIn3))) == 0){
    seqlevels(table3layers_coveredIn3) <- paste0("chr", seqlevels(table3layers_coveredIn3))
  }
  
  ### SAVED ###
  save(table3layers_coveredIn3, file =
         here(paste0("gitignore/table3layers_coveredIn3_", format(Sys.Date(), "%d_%m_%y"), ".Rda")))
}

################################################################################
################################## CHECKPOINT ##################################
################################################################################
load(here(paste0("gitignore/table3layers_coveredIn3_26_08_26.Rda")))

##################################
## Distribution of logBF_per_ds ##
##################################

if (!file.exists(here("B_MultiTissues/dataOut/figures/script04/DistributionProba.pdf"))){
  makeplotdist <- function(x = "logBF_per_ds_allLayers"){
    qa <- quantile(mcols(table3layers_coveredIn3)[[x]],
                   probs = c(0.5, 0.9, 0.95, 0.98, 0.99), na.rm = TRUE)
    
    qa_df <- data.frame(quantile = c("50th percentile", "90th percentile", 
                                     "95th percentile", "98th percentile", "99th percentile"),
                        value = as.numeric(qa))
    
    pdist <- ggplot(data.frame(table3layers_coveredIn3), aes(x = .data[[x]])) +
      geom_histogram(bins = 60, fill = "grey80", colour = "white") +
      geom_vline(data = qa_df, aes(xintercept = value, colour = quantile),
                 linetype = "dashed",linewidth = 0.8) +
      geom_text(data = qa_df, aes(
        x = value, y = Inf, label = paste0(quantile, " = ",format(round(value, 3), nsmall = 3)),
        colour = quantile), angle = 90, vjust = 1.2, hjust = 1.1, show.legend = FALSE) +
      labs(x = "Atlas score", y = "Number of CpGs", title = "Distribution of atlas scores",
           subtitle = x) +
      theme_bw() + theme(legend.position = "none")
    
    return(pdist)
  }
  
  ggsave(here("B_MultiTissues/dataOut/figures/script04/DistributionProba.pdf"),
         plot = cowplot::plot_grid(
           makeplotdist(x = "logBF_per_ds_allLayers"), 
           makeplotdist(x = "logBF_per_ds_meso"),
           makeplotdist(x = "logBF_per_ds_endo"),
           makeplotdist(x = "logBF_per_ds_ecto"),
           nrow = 2, labels = c("A", "B", "C", "D")),
         width = 12, height = 6, dpi = 300, bg = "white")
}

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

## What is the top 1% in the 3 layers (overlap)?
top99q_endo <- quantile(table3layers_coveredIn3$logBF_per_ds_endo, probs = 0.99, na.rm = FALSE)
top99q_ecto <- quantile(table3layers_coveredIn3$logBF_per_ds_ecto, probs = 0.99, na.rm = FALSE)
top99q_meso <- quantile(table3layers_coveredIn3$logBF_per_ds_meso, probs = 0.99, na.rm = FALSE)

top99q_in3layersoverlap_CpGs <- table3layers_coveredIn3[
  table3layers_coveredIn3$logBF_per_ds_endo >= top99q_endo &
    table3layers_coveredIn3$logBF_per_ds_ecto >= top99q_ecto &
    table3layers_coveredIn3$logBF_per_ds_meso >= top99q_meso, ]$chr_pos

message(paste0("Total top99q in3layersoverlap CpG sites: ", length(top99q_in3layersoverlap_CpGs), " (",
               round(length(top99q_in3layersoverlap_CpGs)/length(table3layers_coveredIn3)*100,2), "% of total)"))
# Total top99q in3layersoverlap CpG sites: 60424 (0.3% of total)
length(intersect(top99q_in3layersoverlap_CpGs, top99q_CpGs)) # 60424

if (!file.exists(here(paste0("gitignore/top99q_CpGs_", variant, ".RDS")))){
  saveRDS(top99q_CpGs, here(paste0("gitignore/top99q_CpGs_", variant, ".RDS")))  
}

if (!file.exists(here(paste0("gitignore/top99q_in3layersoverlap_CpGs_", variant, ".RDS")))){
  saveRDS(top99q_in3layersoverlap_CpGs, here(paste0("gitignore/top99q_in3layersoverlap_CpGs_", variant, ".RDS")))  
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
    table3layers_coveredIn3dt, plotDerakhshan = FALSE, centro = centro) 
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
# 1:      1        NA 124793275  2292432
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

if (!file.exists(here(paste0("gitignore/TEplot_plotclass_", variant, ".RDS")))){
  
  # UCSC RepeatMasker annotations (Oct2022) for Human (hg38) from AnnotationHub
  ah <- AnnotationHub()
  query(ah, c("UCSC", "RepeatMasker", "Homo sapiens"))
  
  # Retrieve the desired resource, UCSC RepeatMasker annotations for hg38:
  rmskhg38 <- ah[["AH111333"]]
  
  te_classes <- c("LINE", "SINE", "LTR", "DNA", "RC", "Retroposon")
  te_regions <- rmskhg38[mcols(rmskhg38)$repClass %in% te_classes]
  
  # View summary
  table(mcols(te_regions)$repFamily, mcols(te_regions)$repClass)
  #                   DNA    LINE     LTR      RC Retroposon    SINE
  # 5S-Deu-L2           0       0       0       0          0    2588
  # Alu                 0       0       0       0          0 1282821
  # CR1                 0   69603       0       0          0       0
  # Crypton             3       0       0       0          0       0
  # Crypton-A           1       0       0       0          0       0
  # DNA              2430       0       0       0          0       0
  # Dong-R4             0     548       0       0          0       0
  # ERV1                0       0  187201       0          0       0
  # ERV1?               0       0    1288       0          0       0
  # ERVK                0       0   11903       0          0       0
  # ERVL                0       0  171950       0          0       0
  # ERVL-MaLR           0       0  367470       0          0       0
  # ERVL?               0       0    2326       0          0       0
  # Gypsy               0       0   17427       0          0       0
  # Gypsy?              0       0    7548       0          0       0
  # hAT              8888       0       0       0          0       0
  # hAT-Ac           4596       0       0       0          0       0
  # hAT-Blackjack   20062       0       0       0          0       0
  # hAT-Charlie    270942       0       0       0          0       0
  # hAT-Tag1          244       0       0       0          0       0
  # hAT-Tip100      47824       0       0       0          0       0
  # hAT-Tip100?      2041       0       0       0          0       0
  # hAT?             1543       0       0       0          0       0
  # Helitron            0       0       0    1820          0       0
  # I-Jockey            0       5       0       0          0       0
  # Kolobok             4       0       0       0          0       0
  # L1                  0 1031524       0       0          0       0
  # L1-Tx1              0       1       0       0          0       0
  # L2                  0  486431       0       0          0       0
  # LTR                 0       0    3438       0          0       0
  # Merlin             63       0       0       0          0       0
  # MIR                 0       0       0       0          0  616589
  # MULE-MuDR        2078       0       0       0          0       0
  # MULE-MuDR?          7       0       0       0          0       0
  # Penelope            0    1128       0       0          0       0
  # PIF-Harbinger      33       0       0       0          0       0
  # PiggyBac         2253       0       0       0          0       0
  # PiggyBac?         225       0       0       0          0       0
  # RTE-BovB            0    9297       0       0          0       0
  # RTE-X               0   15944       0       0          0       0
  # SVA                 0       0       0       0       5974       0
  # TcMar             172       0       0       0          0       0
  # TcMar-Mariner   16871       0       0       0          0       0
  # TcMar-Pogo         35       0       0       0          0       0
  # TcMar-Tc1          15       0       0       0          0       0
  # TcMar-Tc2        8479       0       0       0          0       0
  # TcMar-Tigger   123237       0       0       0          0       0
  # TcMar?            358       0       0       0          0       0
  # tRNA                0       0       0       0          0    2295
  # tRNA-Deu            0       0       0       0          0     643
  # tRNA-RTE            0       0       0       0          0    5695
  length(te_regions)  # Total TE regions
  
  top99q_CpGs_GR <- makeGRfromMyCpGPos(top99q_CpGs, "top99q_CpGs")
  top99q_in3layersoverlap_CpGs_GR <- makeGRfromMyCpGPos(top99q_in3layersoverlap_CpGs, "top99q_in3layersoverlap_CpGs")
  totalSites_GR  <- makeGRfromMyCpGPos(table3layers_coveredIn3dt$chr_pos, "totalSites")
  
  makeTEplot <- function(top = top99q_CpGs_GR, whichtop = "top99q_CpGs",
                         mode = c("shell","nested"), famorclass = c("family","class")){
    mode       <- match.arg(mode)
    famorclass <- match.arg(famorclass)
    
    bg_only_GR <- makeGRfromMyCpGPos(table3layers_coveredIn3dt$chr_pos[
      !table3layers_coveredIn3dt$chr_pos %in% top], "bg_only")
    res_cols <- c("label","pvalue","odds_ratio","conf_low","conf_high")
    
    ## ── choose the grouping column, the split, and the keys ONCE ────────────
    grp_col <- if (famorclass == "family") mcols(te_regions)$repFamily
    else                        mcols(te_regions)$repClass
    gr_by   <- split(te_regions, grp_col)
    
    grp_tab <- table(grp_col)
    if (famorclass == "family") {
      keep <- names(grp_tab)[grp_tab >= 1000 & !grepl("\\?$", names(grp_tab))]
    } else {
      keep <- names(grp_tab)                       # keep all 6 TE classes
    }
    
    ## windows
    if (mode == "nested") {
      wins  <- list("TE body"=0L, "TE+1kb"=1000L, "TE+3kb"=3000L, "TE+10kb"=10000L)
      build <- function(gr, w) trim(gr + wins[[w]])
    } else {
      bounds <- list("TE body"=c(0,0), "0-1kb"=c(0,1000),
                     "1-3kb"=c(1000,3000), "3-10kb"=c(3000,10000),
                     "10-15kb"=c(10000,15000))
      build <- function(gr, w) {
        b <- bounds[[w]]
        if (b[1]==b[2]) reduce(trim(gr))
        else GenomicRanges::setdiff(reduce(trim(gr + b[2])), reduce(trim(gr + b[1])))
      }
      wins <- bounds
    }
    
    run_one_window <- function(w) {
      res <- lapply(keep, function(f)
        fisher_test_te(build(gr_by[[f]], w), top, bg_only_GR,
                       label = f, nameTarget = whichtop))
      dt <- data.table::rbindlist(lapply(res, `[`, res_cols))
      dt[, window := w]
      dt
    }
    
    res_dt <- data.table::rbindlist(lapply(names(wins), run_one_window))
    res_dt[, window := factor(window, levels = names(wins))]
    res_dt[, p.adj  := p.adjust(pvalue, "BH")]
    res_dt[, sig    := p.adj < 0.05]
    res_dt[, log2OR := log2(odds_ratio)]            # <- needed for ordering & heatmap
    
    ref <- res_dt[window == "TE body"][order(log2OR), unique(as.character(label))]
    res_dt[, label := factor(as.character(label), levels = ref)]
    res_dt[]
  }
  
  res_dt_class <- makeTEplot(top99q_CpGs_GR, "top99q_CpGs", mode = "shell", famorclass = "class")
  res_dt_fam <- makeTEplot(top99q_CpGs_GR, "top99q_CpGs", mode = "shell", famorclass = "family")
  
  plotclass <- ggplot(res_dt_class, aes(window, label, fill = log2OR)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(sig, "*", "")), vjust = 0.75, size = 6) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#DC3220",
                         midpoint = 0, name = expression(log[2]~OR)) +
    labs(x = NULL, y = NULL,
         title = "TE class enrichment of top 1% hypervariable CpGs\nby distance from element") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  plotfam <- ggplot(res_dt_fam, aes(window, label, fill = log2OR)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(sig, "*", "")), vjust = 0.75, size = 6) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#DC3220",
                         midpoint = 0, name = expression(log[2]~OR)) +
    labs(x = NULL, y = NULL,
         title = "TE family enrichment of top 1% hypervariable CpGs\nby distance from element") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  saveRDS(res_dt_class, here("gitignore/TEplot_res_dt_class.RDS"))
  saveRDS(res_dt_fam, here("gitignore/TEplot_res_dt_fam.RDS"))
  saveRDS(plotclass, here(paste0("gitignore/TEplot_plotclass_", variant, ".RDS")))
  saveRDS(plotfam, here(paste0("gitignore/TEplot_plotfam_", variant, ".RDS")))
}

########################################################################
## KZFP BINDING-SITE enrichment of top 1% hypervariable CpGs          ##
##                                                                    ##
## Data : KRABopedia hg38 peaks (Imbeault 2017, GSE78099), 293T,      ##
##        one BED per KZFP experiment (_peaks.bed = binding intervals).##
## Model: same Fisher target-vs-background overlap as the TE script.  ##
##                                                                    ##
## Three tests:                                                       ##
##   (1) pooled "bound by any KZFP" + distance shells (decay)         ##
##   (2) DECISIVE: hv enriched at KZFP-bound vs unbound TEs           ##
##       (conditional on being in a TE -> removes the TE confound)    ##
##   (3) per-factor: which KZFPs' targets are most hypervariable      ##
########################################################################

if (!file.exists(here(paste0("gitignore/KZFP_bindingSite_", variant, ".RDS")))) {
  
  peak_dir <- here("gitignore/hg38_peaks_un_filt_Krabopedia6oct26")
  
  ## ============================================================
  ## 0. LOAD KZFP BINDING SITES  (keep ONLY KZFPs)
  ## ============================================================
  files  <- list.files(peak_dir, pattern = "_peaks\\.bed$", full.names = TRUE)
  
  # factor = token before the FIRST underscore. Gene symbols contain no "_",
  # so this keeps ZNF324B / ZNF75D / ZNF585A etc. intact while pooling
  # ZNF100_rep1_pubM, ZNF100_rep2_pubM, ... -> ZNF100.
  factor <- sub("_.*$", "", basename(files))
  
  # --- keep ONLY KZFP binding sites: drop the 3 non-KZFP entries ---
  non_kzfp <- c("293TexoT",  # 293T control track (not a factor)
                "PRDM9",      # PR/SET meiotic ZNF, no KRAB domain
                "ZSCAN18")    # SCAN-only zinc finger, no KRAB domain
  keep   <- !(factor %in% non_kzfp)
  files  <- files[keep];  factor <- factor[keep]
  
  # drop low-confidence ChIPs (several are huge & noisy, e.g. ZNF677)
  lc     <- grepl("lowconf", basename(files))
  files  <- files[!lc];   factor <- factor[!lc]
  
  stopifnot(!any(grepl("^293T|^PRDM9$|^ZSCAN18$", factor)))   # sanity
  message(length(unique(factor)), " KZFP factors from ", length(files), " files")
  # 277 KZFP factors from 298 files
  
  # read, tag by factor, pool each factor's runs/replicates into ONE set
  read_bed <- function(f) { g <- import(f, format = "BED"); seqlevelsStyle(g) <- "UCSC"; g }
  peaks_by_factor <- tapply(files, factor, function(fs)
    reduce(trim(do.call(c, lapply(fs, read_bed)))))
  peaks_by_factor <- GRangesList(lapply(
    split(files, factor),
    function(fs) reduce(trim(do.call(c, lapply(fs, read_bed))))
  ))  
  
  # drop failed ChIPs (too few peaks). Inspect sort(lengths(.)) to tune.
  MIN_PEAKS <- 100
  peaks_by_factor <- peaks_by_factor[lengths(peaks_by_factor) >= MIN_PEAKS]
  message(length(peaks_by_factor), " KZFP factors kept after >= ", MIN_PEAKS, " peaks")
  # 228 KZFP factors kept after >= 100 peaks
  
  # union across all KZFPs = "bound by any KZFP"
  kzfp_bound <- reduce(unlist(peaks_by_factor))
  message("union KZFP-bound regions: ", length(kzfp_bound))
  # union KZFP-bound regions: 1289764
  
  ## ============================================================
  ## 1. COMMON OBJECTS
  ## ============================================================
  top        <- top99q_CpGs_GR
  bg_only_GR <- makeGRfromMyCpGPos(table3layers_coveredIn3dt$chr_pos[
    !table3layers_coveredIn3dt$chr_pos %in% top99q_CpGs], "bg_only")
  
  res_cols <- c("label", "pvalue", "odds_ratio", "conf_low", "conf_high")
  
  # safe wrapper: reduce/trim the region, tolerate an empty set, return a 1-row dt
  te_or_na <- function(region, target, background, label) {
    region <- reduce(trim(region))
    if (length(region) == 0L)
      return(data.table(label = label, pvalue = NA_real_, odds_ratio = NA_real_,
                        conf_low = NA_real_, conf_high = NA_real_))
    r <- fisher_test_te(region, target, background,
                        label = label, nameTarget = "top99q_CpGs")
    as.data.table(r[res_cols])
  }
  
  ## ============================================================
  ## 2. TEST 1 — pooled "bound by any KZFP" with distance shells
  ## ============================================================
  shells <- list("peak"   = c(0, 0),    "0-1kb"   = c(0, 1000),
                 "1-3kb"  = c(1000, 3000), "3-10kb" = c(3000, 10000),
                 "10-15kb"= c(10000, 15000))
  build  <- function(gr, b){
    if (b[1] == b[2]) reduce(trim(gr))
    else GenomicRanges::setdiff(reduce(trim(gr + b[2])), reduce(trim(gr + b[1])))}
  
  res_shell <- rbindlist(lapply(names(shells), function(w) {
    dt <- te_or_na(build(kzfp_bound, shells[[w]]), top, bg_only_GR,
                   label = "KZFP-bound")
    dt[, window := w]; dt
  }))
  res_shell[, window := factor(window, levels = names(shells))]
  res_shell[, p.adj  := p.adjust(pvalue, "BH")]
  res_shell[, sig    := p.adj < 0.05]
  print(res_shell[])
  #         label        pvalue odds_ratio  conf_low conf_high  window         p.adj    sig
  # 1: KZFP-bound  0.000000e+00  0.7791534 0.7692930 0.7891024    peak  0.000000e+00   TRUE
  # 2: KZFP-bound  1.523958e-01  0.9936093 0.9848967 1.0023874   0-1kb  1.523958e-01  FALSE
  # 3: KZFP-bound 1.463525e-113  1.1239430 1.1127375 1.1352422   1-3kb 3.658813e-113   TRUE
  # 4: KZFP-bound  1.103181e-95  1.2170571 1.1952168 1.2392084  3-10kb  1.838635e-95   TRUE
  # 5: KZFP-bound  2.819590e-04  1.2161231 1.0945366 1.3475939 10-15kb  3.524487e-04   TRUE
  
  ## ============================================================
  ## 3. TEST 2 — DECISIVE: hv at KZFP-bound vs unbound TEs
  ##    (conditional on being in a TE -> TE confound removed)
  ## ============================================================
  
  if(!exists(x = "te_regions")){
    # UCSC RepeatMasker annotations (Oct2022) for Human (hg38) from AnnotationHub
    ah <- AnnotationHub()
    query(ah, c("UCSC", "RepeatMasker", "Homo sapiens"))
    
    # Retrieve the desired resource, UCSC RepeatMasker annotations for hg38:
    rmskhg38 <- ah[["AH111333"]]
    
    te_classes <- c("LINE", "SINE", "LTR", "DNA", "RC", "Retroposon")
    te_regions <- rmskhg38[mcols(rmskhg38)$repClass %in% te_classes]
  }
  
  te_all     <- reduce(trim(te_regions))
  in_te      <- overlapsAny(totalSites_GR, te_all, ignore.strand = TRUE)  # same order as chr_pos
  te_chr_pos <- table3layers_coveredIn3dt$chr_pos[in_te]
  
  te_top_GR  <- makeGRfromMyCpGPos(intersect(te_chr_pos, top99q_CpGs), "te_top")  # hv,     in TE
  te_bg_GR   <- makeGRfromMyCpGPos(setdiff(te_chr_pos,   top99q_CpGs), "te_bg")   # non-hv, in TE
  
  te_bound   <- subsetByOverlaps(te_regions, kzfp_bound)
  te_unbound <- te_regions[!overlapsAny(te_regions, kzfp_bound)]
  
  # (a) the single conditional OR — the headline statistic
  res_cond <- te_or_na(te_bound, te_top_GR, te_bg_GR,
                       "hv at KZFP-bound vs unbound TE (TE-conditional)")
  res_cond[, p.adj := p.adjust(pvalue, "BH")]
  print(res_cond[])
  #                                                label       pvalue odds_ratio conf_low conf_high        p.adj
  #   1: hv at KZFP-bound vs unbound TE (TE-conditional) 2.834361e-60   1.106766  1.09347   1.12025 2.834361e-60
  
  # (b) bound vs unbound TEs, each vs the genome background (for the plot)
  res_split <- rbindlist(list(
    te_or_na(te_bound,   top, bg_only_GR, "TE: KZFP-bound"),
    te_or_na(te_unbound, top, bg_only_GR, "TE: unbound")))
  res_split[, p.adj := p.adjust(pvalue, "BH")][, sig := p.adj < 0.05]
  print(res_split[])
  #             label        pvalue odds_ratio conf_low conf_high         p.adj    sig
  # 1: TE: KZFP-bound 5.426718e-134   1.146341 1.134055  1.158678 1.085344e-133   TRUE
  # 2:    TE: unbound  1.180616e-02   1.011907 1.002622  1.021269  1.180616e-02   TRUE
  
  peak_core_te <- te_or_na(subsetByOverlaps(kzfp_bound, te_all),
                           te_top_GR, te_bg_GR,
                           "hv at KZFP footprint (TE-conditional)")
  print(peak_core_te)
  #                                   label       pvalue odds_ratio  conf_low conf_high
  # 1: hv at KZFP footprint (TE-conditional) 4.077268e-35  0.9005643 0.8854686 0.9158105
  
  ## ============================================================
  ## 4. TEST 3 — per-factor enrichment (peak-level)
  ## ============================================================
  per_kzfp <- rbindlist(lapply(names(peaks_by_factor), function(f)
    te_or_na(peaks_by_factor[[f]], top, bg_only_GR, f)))
  per_kzfp[, p.adj  := p.adjust(pvalue, "BH")]
  per_kzfp[, sig    := p.adj < 0.05]
  per_kzfp[, log2OR := log2(odds_ratio)]
  setorder(per_kzfp, -log2OR)
  print(head(per_kzfp, 25))
  
  ## test ZNF808 for Matt
  per_kzfp[per_kzfp$label %in% "ZNF808",] # not significant
  
  ## ============================================================
  ## 5. PLOTS  (house style: #DC3220 significant, grey60 n.s.)
  ## ============================================================
  sig_cols <- c(`TRUE` = "#DC3220", `FALSE` = "grey60")
  sig_labs <- c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s.")
  
  ## (1) pooled distance-decay
  p_shell <- ggplot(res_shell, aes(window, odds_ratio)) +
    geom_hline(yintercept = 1, linetype = 3) +
    geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.2) +
    geom_point(aes(colour = sig), size = 3) +
    scale_colour_manual(values = sig_cols, labels = sig_labs, name = NULL) +
    scale_y_log10() +
    labs(x = NULL, y = "Odds ratio (top1% hvCpG vs background)",
         title = "hvCpG enrichment at KZFP binding sites, by distance") +
    theme_minimal(base_size = 12)
  
  ## (2) bound vs unbound TE
  res_split[, label := factor(label, levels = c("TE: unbound", "TE: KZFP-bound"))]
  p_split <- ggplot(res_split, aes(odds_ratio, label)) +
    geom_vline(xintercept = 1, linetype = 3) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), height = 0.2, orientation = "y") +
    geom_point(aes(colour = sig), size = 3) +
    scale_colour_manual(values = sig_cols, labels = sig_labs, name = NULL) +
    scale_x_log10() +
    labs(x = "Odds ratio (vs background)", y = NULL,
         title = "hvCpG enrichment: KZFP-bound vs unbound TEs",
         subtitle = sprintf("TE-conditional OR = %.2f (FDR %.1e)",
                            res_cond$odds_ratio, res_cond$p.adj)) +
    theme_minimal(base_size = 12)
  
  ## (3) per-factor forest — most extreme factors for legibility
  N_SHOW   <- 30
  show_fac <- rbind(head(per_kzfp, N_SHOW), tail(per_kzfp, 5))   # top enriched + a few depleted
  show_fac <- unique(show_fac)
  show_fac[, label := factor(label, levels = label[order(log2OR)])]
  p_fac <- ggplot(show_fac, aes(odds_ratio, label)) +
    geom_vline(xintercept = 1, linetype = 3) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), height = 0.3, orientation = "y") +
    geom_point(aes(colour = sig), size = 2.5) +
    scale_colour_manual(values = sig_cols, labels = sig_labs, name = NULL) +
    scale_x_log10() +
    labs(x = "Odds ratio (top1% hvCpG vs background)", y = NULL) +
    theme_minimal(base_size = 10)
  
  ## assemble if cowplot is available
  if (requireNamespace("cowplot", quietly = TRUE)) {
    fig <- cowplot::plot_grid(
      cowplot::plot_grid(p_shell, p_split, ncol = 1, rel_heights = c(1, 0.8)),
      p_fac, ncol = 2, rel_widths = c(1, 1))
    ggplot2::ggsave(
      here("B_MultiTissues/dataOut/figures/script04/KZFP_bindingSite_enrichment.png"),
      fig, width = 12, height = 8, dpi = 300, bg = "white")
  }
  
  ## ============================================================
  ## 6. SAVE
  ## ============================================================
  saveRDS(list(res_shell = res_shell, res_cond = res_cond, res_split = res_split,
               per_kzfp = per_kzfp, n_factors = length(peaks_by_factor)),
          here(paste0("gitignore/KZFP_bindingSite_", variant, ".RDS")))
  saveRDS(list(p_shell = p_shell, p_split = p_split, p_fac = p_fac),
          here(paste0("gitignore/KZFP_bindingSite_plots_", variant, ".RDS")))
}

#################################
## plot manhattan, feature, TE ##
#################################

if (!file.exists(here("B_MultiTissues/dataOut/figures/script04/MappingVariability.png"))){
  plotManhattan_noDerakh <- readRDS(here(paste0("gitignore/plotManhattan_noDerakh_", variant, ".RDS")))
  
  pfeatures <- readRDS(here(paste0("gitignore/pfeatures_", variant, ".RDS")))
  row1 <- plot_grid(plotManhattan_noDerakh+ ylab("Hypervariability score"),
                    pfeatures + theme_minimal(base_size = 12) +
                      ggtitle(""),
                    ncol = 2, rel_widths = c(2,1),
                    labels = c("A", "B"))
  
  TEplotfam <- readRDS(here("gitignore/TEplot_plotfam_SNP_SDASMrm.RDS"))
  TEplotclass <- readRDS(here("gitignore/TEplot_plotclass_SNP_SDASMrm.RDS"))
  row2 <- plot_grid(TEplotclass, TEplotfam, ncol = 2,
                    labels = c("C", "D"))
  
  KZFPPlots <- readRDS(here(paste0("gitignore/KZFP_bindingSite_plots_", variant, ".RDS")))
  row3 <- plot_grid(KZFPPlots$p_shell+ labs(x = NULL, y = "Odds ratio\n(top1% hvCpG vs background)"),
                    KZFPPlots$p_split, ncol = 2,
                    labels = c("E", "F"))
  
  ggplot2::ggsave(
    filename = here::here(
      "B_MultiTissues/dataOut/figures/script04/MappingVariability.pdf"),
    plot = plot_grid(
      row1 , row2, row3, nrow = 3, rel_heights = c(1,2,1)),
    width = 15, height = 12, dpi = 300, bg = "white")
  
  ggplot2::ggsave(
    filename = here::here(
      "B_MultiTissues/dataOut/figures/script04/MappingVariability.png"),
    plot = plot_grid(
      row1 , row2, row3, nrow = 3, rel_heights = c(1,2,1)),
    width = 15, height = 12, dpi = 300, bg = "white")
}

#######################
## Test GO of top 1% ##
#######################

totalSites <- table3layers_coveredIn3$chr_pos
if (!exists("top99q_CpGs")){
  top99q <- readRDS(here(paste0("gitignore/top99q_CpGs_", variant, ".RDS")))
  top99q_in3layersoverlap_CpGs <- readRDS(here(paste0("gitignore/top99q_in3layersoverlap_CpGs_", variant, ".RDS")))
}

# Method. ClusterProfiler
## 1. Keep CpGs in regions where at least 2 CpGs are in 50bp distance to each other
## 2. annotate with associated genes (in gene body or +/- 10kb from TSS)
## 3. run GO term enrichment with clusterProfiler::enrichGO

minimum_CpG_per_cluster = 2

getGOtop <- function(top){
  
  ## Create universe
  universe <- suppressWarnings(suppressMessages(
    annotateCpGs_txdb(
      clusterCpGs(totalSites, max_gap = 50, min_size = minimum_CpG_per_cluster),
      tss_window = 10000)
  ))
  
  print(paste0("Gene universe contains ", length(universe), " genes"))
  ## "Gene universe contains 32326 genes"
  
  fg_bypass  <- unique(cpg_cluster_count_per_gene(top)$entrez_id)
  fg_wrapper <- annotateCpGs_txdb(clusterCpGs(top, 50, 2), 10000)
  message("Jaccard:")
  print(length(intersect(fg_bypass, fg_wrapper)) / length(union(fg_bypass, fg_wrapper)))  # want ≈ 1
  ## The Jaccard is 1.0 --> the cpg_cluster_count_per_gene counts genes under exactly 
  # the same body∪promoter rule as the annotation wrapper. The CpG-count covariate is valid.
  
  # WGBS-correct control
  res_cpg <- CpG_GO_pipeline_lengthControlled(
    top, universe = universe, min_size = minimum_CpG_per_cluster,
    control_method = "cpg_count", all_sites = totalSites)
  
  # for the sensitivity table, also run the other two
  res_len  <- CpG_GO_pipeline_lengthControlled(top, universe = universe,
                                               control_method = "length")
  
  res_none <- CpG_GO_pipeline_lengthControlled(top, universe = universe,
                                               control_method = "none")
  
  ## cpg_count matching can under-power the protocadherin signal specifically,
  # because the PCDH clusters are so CpG-dense that few comparable genes exist to match them
  # so if adhesion terms weaken here, that's not proof they're artefactual.
  
  print("Direct test of prothocadherin")
  
  pcdh <- GRanges("chr5", IRanges(140710000, 141510000))     # hg38 PCDH clusters
  top_gr <- makeGRfromMyCpGPos(top, "top")
  bg_gr  <- makeGRfromMyCpGPos(totalSites,  "bg")
  obs <- sum(overlapsAny(top_gr, pcdh))
  exp <- length(top_gr) * mean(overlapsAny(bg_gr, pcdh))
  print(c(observed = obs, expected = exp, fold = obs / exp))
  in_top <- sum(overlapsAny(top_gr, pcdh))
  in_bg  <- sum(overlapsAny(bg_gr,  pcdh))
  mat <- matrix(c(in_top, length(top_gr) - in_top,
                  in_bg,  length(bg_gr)  - in_bg), nrow = 2, byrow = TRUE)
  
  message("Fisher test:")
  print(fisher.test(mat, alternative = "greater"))
  
  message("Leave-one-out: does removing the protocadherin cluster drop the terms?")
  # PCDH-cluster entrez IDs (hg38 chr5 ~140.7–141.5 Mb)
  pcdh_region <- GRanges("chr5", IRanges(140710000, 141510000))
  genes_gr <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
  pcdh_entrez <- genes_gr$gene_id[overlapsAny(genes_gr, pcdh_region, ignore.strand = TRUE)]
  
  res_cpg_noPCDH <- CpG_GO_pipeline_lengthControlled(
    top, universe = universe,
    control_method = "cpg_count", all_sites = totalSites,
    exclude_genes = pcdh_entrez)
  
  return(list(res_cpg = res_cpg, res_len = res_len, res_none = res_none,
              res_cpg_noPCDH = res_cpg_noPCDH))
}

if(!file.exists(here::here(paste0("B_MultiTissues/dataOut/figures/script04/GOplottop1pc.pdf")))){
  GO_top99q_CpGs <- getGOtop(top = top99q_CpGs)
  
  go_dt <- rbindlist(lapply(names(GO_top99q_CpGs)[1:3], function(method) {
    res <- GO_top99q_CpGs[[method]]
    rbindlist(lapply(names(res), function(ontology) {
      if (is.null(res[[ontology]]) || nrow(as.data.frame(res[[ontology]])) == 0)
        return(NULL)
      x <- as.data.table(as.data.frame(res[[ontology]]))
      x[, `:=`(method = method, ontology = ontology)]
      x
    }), fill = TRUE)
  }), fill = TRUE)
  
  ## GO plot top 10 terms by ontology
  go_top <- go_dt[order(p.adjust),
                  head(.SD, 10), by = .(method, ontology)]
  
  go_top[, Description := factor(Description,
                                 levels = rev(unique(Description)))]
  
  go_top$method[go_top$method == "res_none"] <- "no correction"
  go_top$method[go_top$method == "res_len"] <- "gene-length correction"
  go_top$method[go_top$method == "res_cpg"] <- "CpG density correction"
  
  p <- ggplot(go_top, aes(x = method, y = Description,size = Count,
                          colour = -log10(p.adjust))) +
    geom_point() +
    facet_wrap(~ontology, scales = "free_y", ncol = 2) +
    scale_colour_viridis_c(option = "plasma") +
    labs(x = NULL, y = NULL,
         colour = "-log10 adjusted P",
         size = "Gene count", title = "GO enrichment comparison") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8))
  
  ggplot2::ggsave(
    filename = here::here(paste0("B_MultiTissues/dataOut/figures/script04/GOplottop1pc.pdf")),
    plot = p, width = 10, height = 8)
  ggplot2::ggsave(
    filename = here::here(paste0("B_MultiTissues/dataOut/figures/script04/GOplottop1pc.png")),
    plot = p, width = 10, height = 8)
  ## --> Likely, the broad neurodevelopmental signature is largely a CpG-density artefact, not an ME signal.
  
  # [1] "Direct test of prothocadherin"
  # observed   expected       fold 
  # 161.000000  99.380103   1.620043 
  # Fisher test:
  #   
  #   Fisher's Exact Test for Count Data
  # 
  # data:  mat
  # p-value = 1.029e-08
  # alternative hypothesis: true odds ratio is greater than 1
  # 95 percent confidence interval:
  #  1.41476     Inf
  # sample estimates:
  # odds ratio 
  #   1.620518 
  
  # --> top99q CpGs are 1.62× enriched at the PCDH clusters (161 observed vs 99 expected)
  # a naive GO enrichment shows a broad neurodevelopmental signature (17 BP terms),
  # but this is largely attributable to the higher CpG density of the foreground genes
  # (3.7× the universe); after matching the background on per-gene CpG count,
  # only homophilic cell-cell adhesion, cell junction organization, and
  # negative chemotaxis remain, likely driven by enrichment at the clustered
  # protocadherin locus on chr5 (Fisher's exact test OR = 1.62, p < 0.0001), a
  # a known systemically-variable ME region.
  
  # terms lost by removing PCDH = attributable to the cluster
  bp_with    <- GO_top99q_CpGs$res_cpg$BP@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::pull(Description)
  bp_without <- GO_top99q_CpGs$res_cpg_noPCDH$BP@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::pull(Description)
  setdiff(bp_with, bp_without)
  # [1] "homophilic cell-cell adhesion" "negative chemotaxis"     
  
  cc_with    <- GO_top99q_CpGs$res_cpg$CC@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::pull(Description)
  cc_without <- GO_top99q_CpGs$res_cpg_noPCDH$CC@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::pull(Description)
  setdiff(cc_with, cc_without)
  # 0
  
  mf_with    <- GO_top99q_CpGs$res_cpg$MF@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::pull(Description)
  mf_without <- GO_top99q_CpGs$res_cpg_noPCDH$MF@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::pull(Description)
  setdiff(mf_with, mf_without)
  # 0
}

## =============================================================================
## KEGG ENRICHMENT — same CpG-density-matched design as the GO block
## =============================================================================

## also ran with minimum_CpG_per_cluster <- 2 (in gitignore)
redoKegg = FALSE
if (redoKegg == TRUE) {
  
  minimum_CpG_per_cluster <- 2
  
  # universe = all genes annotated from covered CpGs (entrez), as in GO
  universe <- annotateCpGs_txdb(
    clusterCpGs(totalSites, max_gap = 50, min_size = minimum_CpG_per_cluster),
    tss_window = 10000)
  
  # foreground genes (same clustering + annotation as GO)
  ensg <- annotateCpGs_txdb(clusterCpGs(top99q_CpGs, 50, minimum_CpG_per_cluster), 10000)
  
  # CpG-density-matched background (reuse the helper function), then enrichKEGG
  cpg_count_dt <- cpg_cluster_count_per_gene(totalSites, 50, minimum_CpG_per_cluster, 10000)
  uni_matched  <- match_universe_by_covariate(ensg, universe, cpg_count_dt, "n_cpg")
  
  kk <- enrichKEGG(gene          = as.character(ensg),
                   universe      = as.character(uni_matched),
                   organism      = "hsa",
                   keyType       = "ncbi-geneid",   # entrez ids
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)
  # translate entrez -> symbol in the geneID column for readable plots
  setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  # 0 enriched terms found
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
    "top 1% hv"             = RColorBrewer::brewer.pal(8, "Set2")[6],## name dep on the plot
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
  
  res_quadrants$quadrant[res_quadrants$quadrant %in% "top99q"] <- "top 1% hv"
  res_quadrants$quadrant[res_quadrants$quadrant %in% "allButTop99q"] <- "all but top 1% hv"
  
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
      title = "CpG set enrichment in top 1% hypervariable CpGs by group (vs other group)",
      subtitle = "2x2 Fisher's exact test"
    ) +
    facet_wrap(~ CpG_set, scales = "free_x", nrow = 1) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none" , ## if all significant
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold")
    )
  
  ggplot2::ggsave(
    filename = here::here(
      "B_MultiTissues/dataOut/figures/script04/plot_top99qCpGsEnrichME.pdf"),
    plot = plot_top99qCpGsEnrichME, width = 8, height = 6,  dpi = 300, bg = "white")
  ggplot2::ggsave(
    filename = here::here(
      "B_MultiTissues/dataOut/figures/script04/plot_top99qCpGsEnrichME.png"),
    plot = plot_top99qCpGsEnrichME, width = 8, height = 6,  dpi = 300, bg = "white")
  
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
  
  ## rename
  renameTop <- function(x){
    levels(x)[levels(x)=="top99q"] <- "top 1% hv"
    x <- droplevels(x)
  }
  
  plots$interlayer_corr$group <- renameTop(plots$interlayer_corr$group)
  nD$group <- renameTop(nD$group)
  plots$CpG_summary$group <- renameTop(plots$CpG_summary$group)
  plots$binned_summary_boot$group <- renameTop(plots$binned_summary_boot$group)
  
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
  
  BASE <- 18
  fig_theme <- theme_minimal(base_size = BASE) +
    theme(
      plot.title   = element_text(size = BASE, hjust = 0, margin = margin(b = 6)),
      axis.title   = element_text(size = BASE),
      axis.text    = element_text(size = BASE - 2),
      legend.title = element_text(size = BASE),
      legend.text  = element_text(size = BASE - 2),
      plot.margin  = margin(t = 20, r = 10, b = 10, l = 10)   # replaces `spacing`
    )
  
  row1 <- plot_grid(
    pMElogBF_per_ds + fig_theme + 
      theme(axis.title.x = element_blank(), 
            legend.position = "none",
            axis.text.x = element_text(angle = 20, hjust = 1)),
    pcontrast + fig_theme + theme(legend.position = "none"),
    ncol = 2, labels = c("A", "B"),
    label_size = 16, label_x = 0, hjust = 0, label_y = 0.98, vjust = 1)
  
  row2 <- plot_grid(pdecay + fig_theme +
              theme(legend.position = "inside",
                    legend.position.inside = c(.7, .6),
                    legend.justification = c(0, 0),
                    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3)),
            
            pinterlayer_corr + fig_theme +
              theme(axis.text.x = element_text(angle = 30, hjust = 1),
                    axis.title.x = element_blank(), legend.position = "none"),
            ncol = 2, labels = c("C", "D"),
            label_size = 16, label_x = 0, hjust = 0, label_y = 0.98, vjust = 1)

 row3 <- plot_grid(
    pinterindividual_var + fig_theme + labs(fill = "CpG set") ,
    pbinned + fig_theme +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(colour = "CpG set", fill = "CpG set") ,
    ncol = 2, align = "v", labels = c("E", "F"),
    label_size = 16, label_x = 0, hjust = 0, label_y = 0.98, vjust = 1)

final_plot <- plot_grid(row1, row2, row3, ncol = 1) +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

## Titles:
# A. Distribution of the hypervariability score for each CpG set
# B. Comparison of previous CpG set groups to mQTLcontrols
# C. Decay curve of hypervariability score per percentile
# D. Mean inter-germ-layer correlation for each CpG set
# E. Densities of interindividual variation per CpG (fetal data), by set
# F. Inter-germ-layer correlation per interindividual variation, binned

ggplot2::ggsave(
  filename = here::here(
    "B_MultiTissues/dataOut/figures/script04/CompareWithpreviousMEs.pdf"),
  plot = final_plot, width = 20, height = 16,  dpi = 300, bg = "white")

ggplot2::ggsave(
  filename = here::here(
    "B_MultiTissues/dataOut/figures/script04/CompareWithpreviousMEs.png"),
  plot = final_plot, width = 20, height = 16,  dpi = 300, bg = "white")
}

############################################################
## How many of each putative ME is actually in the top99q? ##
############################################################

if (!exists("listGR")){
  listGR <- list(top99q = makeGRfromMyCpGPos(vec = top99q_CpGs, setname = "top99q"),
                 allButTop99q = makeGRfromMyCpGPos(
                   setdiff(table3layers_coveredIn3$chr_pos, top99q_CpGs), "allButTop99q"))
}

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
tab <- summary_df %>%
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
library(dplyr)
library(gt)

tab <- summary_df %>%
  mutate(
    across(
      starts_with("pc_"),
      ~ scales::percent(.x / 100, accuracy = 0.1)
    ),
    fold = scales::number(
      fold,
      accuracy = 0.1,
      suffix = "×"
    )
  ) %>%
  gt() %>%
  fmt_number(
    columns = starts_with("n_"),
    decimals = 0
  ) %>%
  cols_label(
    set = "Previously published CpG set",
    n_cpg_covered = "Covered CpGs in set",
    n_top99q = "N in top 1% hvCpGs",
    pc_top99q = "%",
    n_rest = "N in rest (non-top 1% hvCpGs)",
    pc_rest = "%",
    fold = "Fold enrichment"
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_column_labels()
  ) %>%
  tab_options(
    table.font.size = 13,
    data_row.padding = px(3)
  )

gtsave(tab,
  filename = here("B_MultiTissues/dataOut/figures/script04/topCpGsEnrichME_table.pdf"))

######################################################################
## Check enrichement of telomeres and centromeres for high geomMean ##
######################################################################

# Centromeres (from cytoBand - acen bands)
# wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz \
# | zcat | awk '$5=="acen"' > centromeres_hg38.bed

# getEnrichCentroTelo()
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
