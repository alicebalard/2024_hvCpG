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
  top99q_in3layersoverlap_CpGs_GR <- makeGRfromMyCpGPos(top99q_in3layersoverlap_CpGs, "top99q_in3layersoverlap_CpGs")
  totalSites_GR  <- makeGRfromMyCpGPos(table3layers_coveredIn3dt$chr_pos, "totalSites")
  
  makeTEplot <- function(top = top99q_CpGs_GR, whichtop = "top99q_CpGs"){
    
    # Strict background = non-hvCpG sites only (must be disjoint from top)
    bg_only_GR <- makeGRfromMyCpGPos(table3layers_coveredIn3dt$chr_pos[
      !table3layers_coveredIn3dt$chr_pos %in% top], "bg_only")
    
    # Enrichment/depletion of `target` CpGs vs `background` CpGs inside a set of repeat regions.
    
    # columns to pull from each result into a table
    res_cols <- c("label", "pvalue", "odds_ratio", "conf_low", "conf_high")
    
    ## ---- All TEs pooled ----
    allTE_test <- fisher_test_te(te_regions, top, bg_only_GR,
                                 label = "TE", nameTarget = "top99q_CpGs")
    print(allTE_test)
    
    ## ---- Per repClass ----
    te_by_class <- split(te_regions, mcols(te_regions)$repClass)
    class_res <- lapply(names(te_by_class), function(cl)
      fisher_test_te(te_by_class[[cl]], top, bg_only_GR,
                     label = cl, nameTarget = "top99q_CpGs"))
    class_dt <- data.table::rbindlist(lapply(class_res, `[`, res_cols))
    class_dt[, p.adj := p.adjust(pvalue, "BH")]
    print(class_dt[order(-odds_ratio)])
    
    ## ---- Per repFamily (drop tiny/uncertain families so ORs are stable) ----
    fam_tab  <- table(mcols(te_regions)$repFamily)
    fam_keep <- names(fam_tab)[fam_tab >= 1000 & !grepl("\\?$", names(fam_tab))]
    
    fam_by  <- split(te_regions, mcols(te_regions)$repFamily)
    fam_res <- lapply(fam_keep, function(f)
      fisher_test_te(fam_by[[f]], top, bg_only_GR,
                     label = f, nameTarget = "top99q_CpGs"))
    
    res_dt <- data.table::rbindlist(lapply(fam_res, `[`, res_cols))
    res_dt[, p.adj := p.adjust(pvalue, "BH")]           # correct across families
    print(res_dt[order(-odds_ratio)])
    
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
    
    saveRDS(TEplot, here(paste0("gitignore/TEplot_", whichtop, "_", variant, ".RDS")))
  }
  makeTEplot(top99q_CpGs_GR, "top99q_CpGs")
  # $label
  # [1] "TE"
  # 
  # $contingency
  # in_TE not_in_TE
  # top99q_CpGs   112737     89730
  # background  10781840   9464839
  # 
  # $pvalue
  # [1] 1.332613e-105
  # 
  # $odds_ratio
  # [1] 1.102905
  # 
  # $conf_low
  # [1] 1.09325
  # 
  # $conf_high
  # [1] 1.112685
  # 
  # label        pvalue odds_ratio  conf_low conf_high         p.adj
  # <char>         <num>      <num>     <num>     <num>         <num>
  # 1:       LINE  0.000000e+00  1.5553095 1.5377996  1.573005  0.000000e+00
  # 2:        LTR 7.649540e-276  1.3128470 1.2938973  1.332009 1.529908e-275
  # 3:        DNA  2.030086e-24  1.1430595 1.1144968  1.172222  3.045129e-24
  # 4: Retroposon  2.380273e-01  1.0818589 0.9449651  1.233174  2.856327e-01
  # 5:         RC  9.039334e-01  1.0291717 0.6084856  1.631416  9.039334e-01
  # 6:       SINE  0.000000e+00  0.7371553 0.7295801  0.744744  0.000000e+00
  # label        pvalue odds_ratio  conf_low conf_high         p.adj
  # <char>         <num>      <num>     <num>     <num>         <num>
  # 1:            L1  0.000000e+00  1.6787164 1.6583288 1.6993290  0.000000e+00
  # 2:          ERVK  5.391939e-20  1.3489968 1.2677907 1.4341111  2.695969e-19
  # 3:     ERVL-MaLR 8.698016e-126  1.3356109 1.3053668 1.3664271 8.698016e-125
  # 4:     MULE-MuDR  5.108937e-02  1.3295892 0.9881317 1.7518940  1.021787e-01
  # 5:          ERVL  6.407374e-55  1.3008409 1.2599006 1.3427905  3.844425e-54
  # 6:           hAT  4.674471e-02  1.2949115 0.9977691 1.6535247  1.001672e-01
  # 7: TcMar-Mariner  2.302027e-04  1.2905693 1.1267510 1.4717433  6.278255e-04
  # 8:          ERV1  6.445483e-69  1.2406186 1.2118375 1.2699249  4.834112e-68
  # 9:    hAT-Tip100  3.763599e-06  1.2276315 1.1260981 1.3359825  1.129080e-05
  # 10:           LTR  2.602610e-01  1.2144678 0.8441980 1.6930860  3.852471e-01
  # 11:         Gypsy  2.876436e-02  1.1710860 1.0138878 1.3459350  7.191091e-02
  # 12: hAT-Blackjack  3.243543e-02  1.1685235 1.0093529 1.3458349  7.485099e-02
  # 13:  TcMar-Tigger  8.088091e-10  1.1454019 1.0973642 1.1949964  3.466325e-09
  # 14:   hAT-Charlie  1.903610e-08  1.1140045 1.0731271 1.1560324  7.138539e-08
  # 15:         RTE-X  2.696730e-01  1.0910537 0.9300923 1.2720495  3.852471e-01
  # 16:     TcMar-Tc2  4.470745e-01  1.0844786 0.8588974 1.3516349  5.412316e-01
  # 17:           SVA  2.380273e-01  1.0818589 0.9449651 1.2331745  3.758325e-01
  # 18:            L2  7.738175e-07  1.0692152 1.0413072 1.0977500  2.579392e-06
  # 19:           CR1  1.988895e-01  1.0547082 0.9701952 1.1445814  3.314825e-01
  # 20:      Helitron  9.039334e-01  1.0291717 0.6084856 1.6314159  1.000000e+00
  # 21:           MIR  3.194412e-01  0.9849600 0.9560782 1.0145126  4.356017e-01
  # 22:     5S-Deu-L2  1.000000e+00  0.9615346 0.5368616 1.5909042  1.000000e+00
  # 23:      tRNA-RTE  8.716416e-01  0.9523713 0.6695017 1.3151591  1.000000e+00
  # 24:          tRNA  1.000000e+00  0.9481621 0.5293797 1.5687321  1.000000e+00
  # 25:      Penelope  1.000000e+00  0.9216565 0.2499683 2.3768415  1.000000e+00
  # 26:      PiggyBac  4.510264e-01  0.8577526 0.5819027 1.2197413  5.412316e-01
  # 27:      RTE-BovB  1.586413e-01  0.8307844 0.6347992 1.0684900  2.974525e-01
  # 28:        hAT-Ac  3.852737e-01  0.7936326 0.4614145 1.2737448  5.025309e-01
  # 29:           Alu  0.000000e+00  0.7248988 0.7171715 0.7326592  0.000000e+00
  # 30:           DNA  1.987925e-01  0.6476503 0.3099628 1.1943424  3.314825e-01
  # label        pvalue odds_ratio  conf_low conf_high         p.adj
  # <char>         <num>      <num>     <num>     <num>         <num>
  
  makeTEplot(top99q_in3layersoverlap_CpGs_GR, "top99q_in3layersoverlap_CpGs")
  # $label
  # [1] "TE"
  # 
  # $contingency
  # in_TE not_in_TE
  # top99q_CpGs    35261     25163
  # background  10781840   9464839
  # 
  # $pvalue
  # [1] 6.478401e-140
  # 
  # $odds_ratio
  # [1] 1.230133
  # 
  # $conf_low
  # [1] 1.210341
  # 
  # $conf_high
  # [1] 1.250265
  # 
  # label       pvalue odds_ratio   conf_low conf_high        p.adj
  # <char>        <num>      <num>      <num>     <num>        <num>
  # 1:       LINE 0.000000e+00  1.5392355 1.50767987 1.5713301 0.000000e+00
  # 2:        LTR 9.661159e-99  1.3407969 1.30589380 1.3765104 2.898348e-98
  # 3: Retroposon 3.129908e-02  1.2778213 1.01266477 1.5913013 3.755890e-02
  # 4:        DNA 9.961784e-05  1.0987018 1.04792619 1.1513696 1.494268e-04
  # 5:       SINE 3.745399e-68  0.8529081 0.83757039 0.8684888 7.490798e-68
  # 6:         RC 2.635788e-01  0.3831443 0.04636648 1.3861755 2.635788e-01
  # label       pvalue odds_ratio   conf_low conf_high        p.adj
  # <char>        <num>      <num>      <num>     <num>        <num>
  #   1:            L1 0.000000e+00  1.6496875 1.61306606 1.6869590 0.000000e+00
  # 2:           LTR 2.299307e-01  1.3952303 0.72005481 2.4412298 5.306093e-01
  # 3:    hAT-Tip100 2.387713e-05  1.3866956 1.19295680 1.6031119 1.023306e-04
  # 4: hAT-Blackjack 1.038024e-02  1.3856883 1.07748324 1.7549754 3.114072e-02
  # 5:     ERVL-MaLR 7.345924e-44  1.3575755 1.30223758 1.4146648 7.345924e-43
  # 6:          ERVK 3.806708e-07  1.3513168 1.20513500 1.5104462 1.903354e-06
  # 7:          ERVL 2.619105e-20  1.3254995 1.25063454 1.4037191 1.571463e-19
  # 8:     MULE-MuDR 2.966080e-01  1.3103329 0.73260252 2.1641110 5.932160e-01
  # 9:     5S-Deu-L2 4.804499e-01  1.2887827 0.47221508 2.8116678 6.910699e-01
  # 10:           SVA 3.129908e-02  1.2778213 1.01266477 1.5913013 8.536113e-02
  # 11:          ERV1 2.069362e-27  1.2725703 1.21967519 1.3272144 1.552021e-26
  # 12: TcMar-Mariner 8.290316e-02  1.2436934 0.95929918 1.5862117 2.072579e-01
  # 13:         RTE-X 3.704673e-01  1.1300288 0.84089451 1.4868134 6.537659e-01
  # 14:  TcMar-Tigger 4.439191e-03  1.1221660 1.03627527 1.2133016 1.479730e-02
  # 15:            L2 2.048230e-04  1.0953816 1.04406535 1.1486259 7.680864e-04
  # 16:         Gypsy 4.837489e-01  1.0931941 0.82533034 1.4204854 6.910699e-01
  # 17:           CR1 4.119120e-01  1.0626189 0.91059748 1.2328826 6.865199e-01
  # 18:      tRNA-RTE 8.822943e-01  1.0350088 0.53431629 1.8102144 9.803269e-01
  # 19:   hAT-Charlie 7.853974e-01  1.0094805 0.93923189 1.0836184 9.062278e-01
  # 20:           MIR 3.188433e-01  0.9722543 0.92022980 1.0264913 5.978311e-01
  # 21:           hAT 1.000000e+00  0.9344619 0.51047466 1.5695110 1.000000e+00
  # 22:           DNA 1.000000e+00  0.8680658 0.23621208 2.2273309 1.000000e+00
  # 23:           Alu 4.993578e-69  0.8482141 0.83250811 0.8641875 7.490367e-68
  # 24:          tRNA 1.000000e+00  0.8472131 0.23054482 2.1736911 1.000000e+00
  # 25:      PiggyBac 7.589880e-01  0.8344222 0.38122165 1.5858478 9.062278e-01
  # 26:      RTE-BovB 5.191092e-01  0.8214371 0.48655567 1.2992922 7.078761e-01
  # 27:     TcMar-Tc2 4.546537e-01  0.8175384 0.48421783 1.2930953 6.910699e-01
  # 28:        hAT-Ac 5.467680e-01  0.6257022 0.17034465 1.6044730 7.131757e-01
  # 29:      Helitron 2.635788e-01  0.3831443 0.04636648 1.3861755 5.648118e-01
  # 30:      Penelope 6.446330e-01  0.0000000 0.00000000 2.8601722 8.057913e-01
}

#################################
## plot manhattan, feature, TE ##
#################################

if (!file.exists(here("B_MultiTissues/dataOut/figures/script04/MappingVariability.png"))){
  plotManhattan_noDerakh <- readRDS(here(paste0("gitignore/plotManhattan_noDerakh_", variant, ".RDS")))
  pfeatures <- readRDS(here(paste0("gitignore/pfeatures_", variant, ".RDS")))
  TEplot <- readRDS(here("gitignore/TEplot_top99q_CpGs_SNP_SDASMrm.RDS"))
  TEplot2 <- readRDS(here("gitignore/TEplot_top99q_in3layersoverlap_CpGs_SNP_SDASMrm.RDS"))
  
  bottomrow <- plot_grid(pfeatures, TEplot, 
                         TEplot2 + labs(title = "1% top overlap 3 layers") , ncol = 3,
                         labels = c("B", "C", "D"))
  ggplot2::ggsave(
    filename = here::here(
      "B_MultiTissues/dataOut/figures/script04/MappingVariability.pdf"),
    plot = plot_grid(
      plotManhattan_noDerakh + ylab("Hypervariability score"), bottomrow,
      labels = c("A", ""), nrow = 2),
    width = 20, height = 10, dpi = 300, bg = "white")
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
  
  p <- ggplot(go_top, aes(x = method, y = Description,size = Count,
                          colour = -log10(p.adjust))) +
    geom_point() +
    facet_wrap(~ontology, scales = "free_y") +
    scale_colour_viridis_c(option = "plasma") +
    labs(x = NULL, y = NULL,
         colour = "-log10 adjusted P",
         size = "Gene count", title = "GO enrichment comparison") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8))
  
  ggplot2::ggsave(
    filename = here::here(paste0("B_MultiTissues/dataOut/figures/script04/GOplottop1pc.pdf")),
    plot = p, width = 14, height = 4)
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
  bp_without
  # [1] "cell junction organization"
}

######################################
## And for top 1% overlap 3 layers? ##
######################################

rerunGO2 = FALSE
if (rerunGO2){
  GO_top99q_in3layersoverlap_CpGs <- getGOtop(top = top99q_in3layersoverlap_CpGs)
  ## No GO terms significant (too few CpGs)  
}
