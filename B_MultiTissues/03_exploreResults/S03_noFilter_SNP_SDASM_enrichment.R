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

variant <- "noFilter"    # was "SNPrm"  -> reads .../Atlas/noFilter results
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
  
  # Saved:
  # /home/alice/Documents/Research/GIT/2024_hvCpG/gitignore/resultsAtlasPrepared/noFilter/fullres_0_8p0_0_65p1_12_endo.rds
  # /home/alice/Documents/Research/GIT/2024_hvCpG/gitignore/resultsAtlasPrepared/noFilter/fullres_0_8p0_0_65p1_13_meso.rds
  # /home/alice/Documents/Research/GIT/2024_hvCpG/gitignore/resultsAtlasPrepared/noFilter/fullres_0_8p0_0_65p1_14_ecto.rds
  # /home/alice/Documents/Research/GIT/2024_hvCpG/gitignore/resultsAtlasPrepared/noFilter/fullres_0_8p0_0_65p1_atlas_general.rds
  
  ################################################################################
  ## Load layer-specific analyses                                               ##
  ################################################################################
  
  endo <- fr("fullres_0_8p0_0_65p1_12_endo.rds")
  meso <- fr("fullres_0_8p0_0_65p1_13_meso.rds")
  ecto <- fr("fullres_0_8p0_0_65p1_14_ecto.rds")
  
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
  
  Atlas_dt <- fr("fullres_0_8p0_0_65p1_atlas_general.rds")
  
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
  
  ## We select only the sites covered in the 3 germ layer analyses and will use logBF_per_ds
  m <- as.matrix(mcols(table3layers_coveredIn3)[, c("logBF_per_ds_endo",
                                                    "logBF_per_ds_ecto",
                                                    "logBF_per_ds_meso")])
  
  table3layers_coveredIn3$logBF_per_ds_mean3layers <- rowMeans(m, na.rm = FALSE)
  
  ## We select only the sites covered in the 3 germ layer analyses and will use logBF_per_ds
  table3layers_coveredIn3$n_layers <- rowSums(!is.na(m))
  
  table(table3layers_coveredIn3$n_layers)
  # 1        2        3 
  # 964581  2250919 23036026 

  ### SAVED ###
  table3layers_coveredIn3 <- table3layers_coveredIn3[table3layers_coveredIn3$n_layers %in% 3,]
  save(table3layers_coveredIn3, file =
         here(paste0("gitignore/table3layers_coveredIn3_noFilter_", format(Sys.Date(), "%d_%m_%y"), ".Rda")))
}

################################################################################
################################## CHECKPOINT ##################################
################################################################################

load(here("gitignore/table3layers_coveredIn3_noFilter_11_11_26.Rda"))

################################################
## Get the top 99% quantile for logBF_per_ds ##
################################################

top99q <- quantile(table3layers_coveredIn3$logBF_per_ds_allLayers, probs = 0.99, na.rm = FALSE)

top99q_CpGs <- table3layers_coveredIn3[table3layers_coveredIn3$logBF_per_ds_allLayers >= top99q, ]$chr_pos

message(paste0("Total CpG sites: ", length(table3layers_coveredIn3)))
message(paste0("Total top99q CpG sites: ", length(top99q_CpGs), " (",
               round(length(top99q_CpGs)/length(table3layers_coveredIn3)*100,2), "% of total)"))
# Total CpG sites: 23036026
# Total top99q CpG sites: 230361 (1% of total)

###################################
## Enrichement in SD-ASM regions ##
###################################
NOTtop99q <- table3layers_coveredIn3$chr_pos[!table3layers_coveredIn3$chr_pos %in% top99q_CpGs]
source(here("gitignore/sdasm_all_chr_pos.txt")) 

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
#               set fg_in fg_total   bg_in bg_total    fg_frac    bg_frac odds_ratio conf_low conf_high        pvalue
#   1: SD-ASM (all) 16207   230361 1347478 22805665 0.07035479 0.05908523   1.205167 1.185887  1.224695 4.136448e-109

# --- 2. per classification ---
cats <- unique(SDASM_GR$classification)
res_cat <- rbindlist(lapply(cats, function(cl)
  sdasm_fisher(reduce(SDASM_GR[SDASM_GR$classification == cl]), cl)))

res <- rbind(res_all, res_cat)
res[, p.adj := p.adjust(pvalue, "BH")]
res[order(-odds_ratio)]
#                        set fg_in fg_total   bg_in bg_total     fg_frac      bg_frac odds_ratio  conf_low conf_high        pvalue         p.adj
# 1:                   ubiq.   234   230361    3269 22805665 0.001015797 0.0001433416  7.0913056 6.1846892 8.1012938 2.217986e-111 1.108993e-110
# 2: tissue-specific demeth.  8536   230361  679451 22805665 0.037054883 0.0297930799  1.2531186 1.2260053 1.2807032  4.435065e-86  7.391776e-86
# 3:            SD-ASM (all) 16207   230361 1347478 22805665 0.070354791 0.0590852317  1.2051672 1.1858871 1.2246949 4.136448e-109 1.034112e-108
# 4:                   other  7206   230361  608505 22805665 0.031281337 0.0266821862  1.1779323 1.1503377 1.2060783  3.256213e-40  3.256213e-40
# 5:            denovo meth.   231   230361   56253 22805665 0.001002774 0.0024666240  0.4059567 0.3551735 0.4619666  1.533142e-57  1.916427e-57

# a positive OR means "hv CpGs are enriched for cis-genetic control relative to typical covered CpGs"

###########################################################
## Enrichment of top 1% hv CpGs in SNPs (1000G, MAF>0.01) ##
##   Only meaningful on the "noFilter" run, where SNP CpGs ##
##   are still present. Membership test on chr_pos mirrors ##
##   exactly how the pipeline EXCLUDES SNP CpGs.           ##
###########################################################
snp_file  <- here("gitignore/snps_maf_0.01_chr_pos.txt")      # produced by prepSNP1000GP.sh
snp_sites <- data.table::fread(snp_file, header = FALSE)$V1   # "chrX_pos" strings

## foreground = top 1% hv ; background = the rest of covered-in-3 (already defined above)
fg_tot <- length(top99q_CpGs); bg_tot <- length(NOTtop99q)
fg_in_snp <- sum(top99q_CpGs %in% snp_sites)                  # top-1% hv CpGs that ARE a SNP position
bg_in_snp <- sum(NOTtop99q  %in% snp_sites)                   # same for the background

## 2x2 Fisher: is the top 1% enriched for SNP-overlapping CpGs vs the rest?
m_snp  <- matrix(c(fg_in_snp, fg_tot - fg_in_snp,
                   bg_in_snp, bg_tot - bg_in_snp), nrow = 2, byrow = TRUE)
ft_snp <- fisher.test(m_snp, alternative = "two.sided")

res_snp <- data.table(
  set        = "SNP (1000G MAF>0.01)",
  fg_in = fg_in_snp, fg_total = fg_tot,
  bg_in = bg_in_snp, bg_total = bg_tot,
  fg_frac = fg_in_snp / fg_tot, bg_frac = bg_in_snp / bg_tot,
  odds_ratio = unname(ft_snp$estimate),
  conf_low = ft_snp$conf.int[1], conf_high = ft_snp$conf.int[2],
  pvalue = ft_snp$p.value)
print(res_snp)
#                       set fg_in fg_total bg_in bg_total    fg_frac     bg_frac odds_ratio conf_low conf_high pvalue
#   1: SNP (1000G MAF>0.01)  2538   230361 37525 22805665 0.01101749 0.001645425    6.75919 6.489171  7.039012      0

###########################################################
## Combined summary: how much of the UNFILTERED top 1% hv ##
## sits in SNPs and in SD-ASM, with enrichment ORs        ##
###########################################################
## res_all = SD-ASM(all) Fisher result computed just above (fg_in = top-1% in SD-ASM)
summary_tab <- rbind(
  res_snp[, .(set, n_top1pc_in = fg_in, pct_of_top1pc = 100*fg_in/fg_total,
              odds_ratio, pvalue)],
  res_all[, .(set, n_top1pc_in = fg_in, pct_of_top1pc = 100*fg_in/fg_total,
              odds_ratio, pvalue)])
print(summary_tab)
#                     set n_top1pc_in pct_of_top1pc odds_ratio        pvalue
# 1: SNP (1000G MAF>0.01)        2538      1.101749   6.759190  0.000000e+00
# 2:         SD-ASM (all)       16207      7.035479   1.205167 4.136448e-109

cat(sprintf(
  "\nUNFILTERED top 1%% hv CpGs: %s\n  in SNPs   : %s (%.1f%%)  OR = %.2f\n  in SD-ASM : %s (%.1f%%)  OR = %.2f\n",
  format(fg_tot, big.mark=","),
  format(fg_in_snp,  big.mark=","), 100*fg_in_snp/fg_tot,  res_snp$odds_ratio,
  format(res_all$fg_in, big.mark=","), 100*res_all$fg_in/fg_tot, res_all$odds_ratio))

# UNFILTERED top 1% hv CpGs: 230,361
# in SNPs   : 2,538 (1.1%)  OR = 6.76
# in SD-ASM : 16,207 (7.0%)  OR = 1.21

