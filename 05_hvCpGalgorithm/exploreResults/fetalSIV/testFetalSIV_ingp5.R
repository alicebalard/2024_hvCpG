setwd("/home/alice/2024_hvCpG")
library(here)

source(here("05_hvCpGalgorithm/quiet_library.R"))
source(here("05_hvCpGalgorithm/exploreResults/prepPreviousSIV.R"))

#################################################################################################
## Load positions to test for SIV (calculated in 05_hvCpGalgorithm/exploreResults/S03 and S07) ##
#################################################################################################

### 1/ based on high alpha
x <- readRDS(here("05_hvCpGalgorithm/exploreResults/fetalSIV/x.RDS"))

### 2/ based on high alpha when test is run also in each germ layer, independently
topIntersect90_pos <- readRDS(here("05_hvCpGalgorithm/exploreResults/fetalSIV/topIntersect90_pos.RDS"))

## Prepare fetal data

fetalData <- readRDS("/mnt/ing-s1/ICH_fetal/fetal_EPIC_processed/Beta_matrices/Foetal_ICH_beta_mat_adjusted_EPIC.RDS")
nrow(fetalData) # 746,492 sites after preprocessing by Maria

fetalData_long <- data.frame(fetalData) %>%
  tibble::rownames_to_column("CpG") %>% # move CpG IDs into a column
  pivot_longer(
    cols = -CpG,
    names_to = "sample_full",
    values_to = "beta"
  ) %>%
  dplyr::mutate(
    sample = str_extract(sample_full, "^S[0-9]+"),      # extract sample ID (e.g. S48)
    tissue = str_extract(sample_full, "(?<=_)[A-Za-z]+(?=_)"),  # between underscores
    layer  = str_extract(sample_full, "[A-Za-z]+$")     # after last underscore
  ) %>%
  dplyr::select(CpG, sample, tissue, layer, beta)

fetalData_long$chrpos_hg38 <- dico$chrpos_hg38[match(fetalData_long$CpG, dico$CpG)]

## Prepare array background

# 372 571 array background CpGs (excluding previously identified SIV CpGs and hvCpGs) (dark grey)
# van Baak et al. 2018
# Harris et al. 2014
# Kessler et al. 2018
# Gunasekara et al. 2019

# Create GRanges for probes (1 bp each)
dico_GRanges_hg38 <- GRanges(
  seqnames = dico$chr_hg38,
  ranges = IRanges(start = dico$pos_hg38, end = dico$pos_hg38),
  names = dico$CpG
)

## Find overlap between Kessler & Gunasekara SIV range and the array probes
KesslerSIV_hg38 <- dico[queryHits(findOverlaps(dico_GRanges_hg38, KesslerSIV_GRanges_hg38)), "chrpos_hg38"]
length(KesslerSIV_hg38)

corSIV_hg38 <- dico[queryHits(findOverlaps(dico_GRanges_hg38, corSIV_GRanges_hg38)), "chrpos_hg38"]

prevDetSIV <- c(VanBaakESS_hg38, HarrisSIV_hg38, KesslerSIV_hg38, corSIV_hg38)

arrayBkg <- unique(setdiff(unique(fetalData_long$chrpos_hg38), 
                           union(prevDetSIV, x$chrpos_hg38)))
length(arrayBkg) # 741,559

arrayBkg2 <- unique(setdiff(unique(fetalData_long$chrpos_hg38), 
                            union(prevDetSIV, topIntersect90_pos$chrpos_hg38)))
length(arrayBkg2) # 742,416

## Maria arrives at 372 571 array background CpGs because she uses only overlaps with array450k 
table(dico[match(arrayBkg, dico$chrpos_hg38),"array"])
# 450k and EPIC          EPIC 
#   380362        361196

table(dico[match(arrayBkg2, dico$chrpos_hg38),"array"])
# 450k and EPIC          EPIC 
# 380864        361551 

makePlotSIV <- function(x=highAlphaPos_atlas0.7){
  fetalData_subset_highVar <- fetalData_long[fetalData_long$chrpos_hg38 %in% x$chrpos_hg38,]
  fetalData_subset_backgrd <- fetalData_long[fetalData_long$chrpos_hg38 %in% arrayBkg,]
  fetalData_subset_SIV <- fetalData_long[fetalData_long$chrpos_hg38 %in% prevDetSIV,]
  
  ## Sanity check that it's different CpGs
  intersect(fetalData_subset_highVar$CpG, fetalData_subset_backgrd$CpG)
  intersect(fetalData_subset_highVar$chrpos_hg38, fetalData_subset_backgrd$chrpos_hg38)
  
  # For the nine individuals with available endoderm–mesoderm samples, we calculated 
  # the Pearson r between germ layer methylation values for each hvCpG and repeated this
  # for individuals with endoderm–ectoderm and mesoderm–ectoderm samples. The inter-germ 
  # layer correlation was then defined as the average Pearson r across these three comparisons.
  
  ## NB: absolute values (just the strength of correlation)
  interlayer_corr_highVar <- fetalData_subset_highVar %>%
    dplyr::select(sample, layer, CpG, beta) %>%
    pivot_wider(names_from = layer, values_from = beta) %>%
    group_by(CpG) %>%
    summarise(
      r_Endo_Meso = cor(Endo, Meso, use = "pairwise.complete.obs"),
      r_Endo_Ecto = cor(Endo, Ecto, use = "pairwise.complete.obs"),
      r_Meso_Ecto = cor(Meso, Ecto, use = "pairwise.complete.obs"),
      interlayer_r = mean(c(r_Endo_Meso, r_Endo_Ecto, r_Meso_Ecto), na.rm = TRUE)
    ) %>% mutate(group = "highVar")
  
  mean(interlayer_corr_highVar$interlayer_r) # 0.57
  
  interlayer_corr_backgrd <- fetalData_subset_backgrd %>%
    dplyr::select(sample, layer, CpG, beta) %>%
    pivot_wider(names_from = layer, values_from = beta) %>%
    group_by(CpG) %>%
    summarise(
      r_Endo_Meso = cor(Endo, Meso, use = "pairwise.complete.obs"),
      r_Endo_Ecto = cor(Endo, Ecto, use = "pairwise.complete.obs"),
      r_Meso_Ecto = cor(Meso, Ecto, use = "pairwise.complete.obs"),
      interlayer_r = mean(c(r_Endo_Meso, r_Endo_Ecto, r_Meso_Ecto), na.rm = TRUE)
    ) %>% mutate(group = "background")
  
  mean(interlayer_corr_backgrd$interlayer_r)
  
  interlayer_corr_SIV <- fetalData_subset_SIV %>%
    dplyr::select(sample, layer, CpG, beta) %>%
    pivot_wider(names_from = layer, values_from = beta) %>%
    group_by(CpG) %>%
    summarise(
      r_Endo_Meso = cor(Endo, Meso, use = "pairwise.complete.obs"),
      r_Endo_Ecto = cor(Endo, Ecto, use = "pairwise.complete.obs"),
      r_Meso_Ecto = cor(Meso, Ecto, use = "pairwise.complete.obs"),
      interlayer_r = mean(c(r_Endo_Meso, r_Endo_Ecto, r_Meso_Ecto), na.rm = TRUE)
    ) %>% mutate(group = "SIV")
  
  mean(interlayer_corr_SIV$interlayer_r)
  
  interlayer_corr <- dplyr::full_join(dplyr::full_join(interlayer_corr_highVar, interlayer_corr_backgrd),
                                      interlayer_corr_SIV)
  
  p1 <- ggplot(interlayer_corr, aes(x=group, y=interlayer_r, group = group, fill = group))+
    geom_violin(width=1.4) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    scale_fill_viridis(discrete = TRUE) +
    theme_minimal(base_size = 14) +
    labs(y = "Mean inter-germ layer correlation\n(Person r)")+
    theme(axis.title.x = element_blank(), legend.position = "none") 
  
  p1
  
  # We calculated interindividual variation using the same metric as van Baak et al. (27):
  # for each CpG, we took the mean methylation value across the two germ-layer derived 
  # tissues for every individual (giving 27 values for each CpG) and defined interindividual 
  # variation of the CpG as the range of these means.
  
  interindividual_var <- fetalData_long %>%
    group_by(sample, CpG) %>%
    summarise(mean_beta = mean(beta, na.rm = TRUE), .groups = "drop") %>%
    group_by(CpG) %>%
    summarise(interindividual_var = max(mean_beta, na.rm = TRUE) - min(mean_beta, na.rm = TRUE))
  
  CpG_summary <- interlayer_corr %>%
    left_join(interindividual_var, by = "CpG")
  
  table(CpG_summary$group)
  # background    highVar        SIV 
  # 741621       1295       3994 
  
  p2 <- ggplot(CpG_summary, aes(x = interindividual_var, color = group)) +
    geom_density(alpha = 0.5)+
    scale_colour_viridis(discrete = TRUE) +
    theme_minimal(base_size = 14) +
    labs(x = "Interindividual variation")
  
  p2
  
  # Bin interindividual_var into 0.1 intervals and bootstrap for CI
  library(boot)
  
  # Function to compute bootstrap CI for median
  boot_median_ci <- function(x, nboot = 1000, conf = 0.95) {
    x <- x[!is.na(x)]
    if (length(x) < 5) return(c(median = NA, low = NA, high = NA))  # skip small bins
    
    bootfun <- function(data, idx) median(data[idx], na.rm = TRUE)
    b <- boot(x, statistic = bootfun, R = nboot)
    ci <- boot.ci(b, type = "perc", conf = conf)
    
    if (!is.null(ci) && "percent" %in% names(ci)) {
      c(median = median(x, na.rm = TRUE), low = ci$percent[4], high = ci$percent[5])
    } else {
      c(median = median(x, na.rm = TRUE), low = NA, high = NA)
    }
  }
  
  # Apply bootstrap per (group, bin)
  binned_summary_boot <- CpG_summary %>%
    mutate(bin = cut(
      interindividual_var,
      breaks = seq(0, max(interindividual_var, na.rm = TRUE) + 0.1, by = 0.1),
      include.lowest = TRUE
    )) %>%
    group_by(group, bin) %>%
    summarise(
      boot_res = list(boot_median_ci(interlayer_r)),
      .groups = "drop"
    ) %>%
    mutate(
      median_r = sapply(boot_res, `[[`, "median"),
      low = sapply(boot_res, `[[`, "low"),
      high = sapply(boot_res, `[[`, "high")
    ) %>%
    dplyr::select(-boot_res)
  
  # Plot
  p3 <- ggplot(binned_summary_boot,
               aes(x = bin, y = median_r, color = group, fill = group)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(
      aes(ymin = low, ymax = high),
      width = 0.2,
      position = position_dodge(width = 0.5)
    ) +
    scale_color_viridis(discrete = TRUE) +
    scale_fill_viridis(discrete = TRUE) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Interindividual variation",
      y = "Inter-germ layer correlation \n(median ± bootstrap CI)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  p3
  
  final_plot <- plot_grid(p1, plot_grid(p2, p3, ncol = 1, align = "v"),
                          ncol = 2, rel_widths = c(.7, 1))
  print(final_plot)  # CRITICAL: print the result to save as pdf
  return(final_plot)
}

pdf(here("05_hvCpGalgorithm/figures/intercorrelationSIVfetal_highAlphaAtlaspos70pc.pdf"),
    width = 10, height = 7)
makePlotSIV()
dev.off()

pdf(here("05_hvCpGalgorithm/figures/intercorrelationSIVfetal_topIntersect90_pos.pdf"),
    width = 10, height = 7)
makePlotSIV(x = topIntersect90_pos)
dev.off()
