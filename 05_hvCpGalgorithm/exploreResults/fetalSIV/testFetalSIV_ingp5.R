setwd("/home/alice/2024_hvCpG")
library(here)
library(ggVennDiagram)
library(boot)

source(here("05_hvCpGalgorithm/quiet_library.R"))
source(here("05_hvCpGalgorithm/exploreResults/prepPreviousSIV.R"))

########################
## Prepare fetal data ##
########################

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

#################################################################################################
## Load positions to test for SIV (calculated in 05_hvCpGalgorithm/exploreResults/S03 and S07) ##
#################################################################################################

### 1/ based on high alpha
highAlphaPos_atlas0.7 <- readRDS(here("05_hvCpGalgorithm/exploreResults/fetalSIV/highAlphaPos_atlas0.7.RDS"))
fetalData_subset_highVar0.7 <- fetalData_long[
  fetalData_long$chrpos_hg38 %in% highAlphaPos_atlas0.7$chrpos_hg38,]

### 2/ based on high alpha when test is run also in each germ layer, independently
topIntersect90_pos <- readRDS(here("05_hvCpGalgorithm/exploreResults/fetalSIV/topIntersect90_pos.RDS"))
fetalData_subset_topIntersect90 <- fetalData_long[
  fetalData_long$chrpos_hg38 %in% topIntersect90_pos$chrpos_hg38,]

#######################################
## Prepare previously identified SIV ##
#######################################
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
fetalData_subset_prevSIV <- fetalData_long[fetalData_long$chrpos_hg38 %in% prevDetSIV,]

##############################
## Prepare array background ##
##############################
fetalData_subset_backgrd <- fetalData_long[
  !fetalData_long$chrpos_hg38 %in% 
    c(fetalData_subset_highVar0.7$chrpos_hg38, 
      fetalData_subset_topIntersect90$chrpos_hg38,
      fetalData_subset_prevSIV$chrpos_hg38),]

## Check overlap on a Venn diagram

cpgs <- list(backgrd = fetalData_subset_backgrd$CpG, highVar0.7 = fetalData_subset_highVar0.7$CpG,
             topIntersect90 = fetalData_subset_topIntersect90$CpG, prevSIV = fetalData_subset_prevSIV$CpG)

ggVennDiagram(cpgs, label_alpha = 0, label = "count") +
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red")+
  theme(legend.position = "none")

#########################
## Shape data for plot ##
#########################

## Derakhshan 2022: For the nine individuals with available endoderm–mesoderm samples, we calculated 
# the Pearson r between germ layer methylation values for each hvCpG and repeated this
# for individuals with endoderm–ectoderm and mesoderm–ectoderm samples. The inter-germ 
# layer correlation was then defined as the average Pearson r across these three comparisons.

## NB: absolute values (just the strength of correlation)
getinterlayer_corr <- function(fetalData_subset, name){
  interlayer_corr = fetalData_subset %>%
    dplyr::select(sample, layer, CpG, beta) %>%
    pivot_wider(names_from = layer, values_from = beta) %>%
    group_by(CpG) %>%
    summarise(
      r_Endo_Meso = cor(Endo, Meso, use = "pairwise.complete.obs"),
      r_Endo_Ecto = cor(Endo, Ecto, use = "pairwise.complete.obs"),
      r_Meso_Ecto = cor(Meso, Ecto, use = "pairwise.complete.obs"),
      interlayer_r = mean(c(r_Endo_Meso, r_Endo_Ecto, r_Meso_Ecto), na.rm = TRUE)
    ) %>% mutate(group = name)
  
  print(mean(interlayer_corr$interlayer_r))
  return(interlayer_corr)
}

interlayer_corr_backgrd <- getinterlayer_corr(fetalData_subset_backgrd, "background")
# mean: 0.057
interlayer_corr_prevSIV <- getinterlayer_corr(fetalData_subset_prevSIV, "prevSIV")
# mean: 0.34
interlayer_corr_highVar0.7 <- getinterlayer_corr(fetalData_subset_highVar0.7, "highVar0.7")
# mean: 0.50
interlayer_corr_topIntersect90 <- getinterlayer_corr(fetalData_subset_topIntersect90, "topIntersect90")
# mean: 0.74

interlayer_corr <- interlayer_corr_backgrd |>
  dplyr::full_join(interlayer_corr_prevSIV) |>
  dplyr::full_join(interlayer_corr_highVar0.7) |>
  dplyr::full_join(interlayer_corr_topIntersect90) |>
  mutate(group = forcats::fct_relevel(group,
    "background", "prevSIV", "highVar0.7", "topIntersect90"))

p1 <- ggplot(interlayer_corr, aes(x=group, y=interlayer_r, group = group, fill = group))+
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal(base_size = 14) +
  labs(y = "Mean inter-germ layer correlation\n(Pearson's r)")+
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
# background        prevSIV     highVar0.7 topIntersect90 
# 741621           3994           1295             50 

p2 <- ggplot(CpG_summary, aes(x = interindividual_var, color = group)) +
  geom_density(alpha = 0.5)+
  scale_colour_viridis(discrete = TRUE) +
  theme_minimal(base_size = 14) +
  labs(x = "Interindividual variation")

p2

###################################################################
# Bin interindividual_var into 0.1 intervals and bootstrap for CI #
###################################################################
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

pdf(here("05_hvCpGalgorithm/figures/intercorrelationSIVfetal.pdf"),
    width = 10, height = 7)
final_plot
dev.off()