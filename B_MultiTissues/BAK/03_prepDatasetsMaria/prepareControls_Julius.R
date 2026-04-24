## Maria's paper:
# "mQTL-matched controls:	3722	CpGs reported by the GoDMC meta-GWAS (51) with 
# the same number of mQTL associations and similar mean % variance explained by an mQTL, 
# requiring each control CpG to be present in at least as many datasets as the hvCpG."

## Julius: Sure, I’ve simplified my script to make it a bit easier for you, I directly load all the data into the R session. 
# The script still includes some internal checks to e.g. check the number of controls at each filtering stage, 
# because my original script focussed on all CpGs but they had to include ≥1 trans mQTL (the current script just matches all hvCpGs)
# on cis and trans effects, irrespective of the requirement for ≥1 trans mQTL. I’ve also attached a link to the website where you 
# can download the data (I recommend reading through the readme file). 
## Note, I only match based on clumped SNPs, when matching on all SNPs (irrespective of clumped or not)
# matching on SNP frequency becomes a lot harder because some CpGs have lots of associations.
# There are also two files mentioned in the script which lists previously discovered MEs from Noah, SIV and ESS CpGs (from Van Baak) 
# and SoC CpGs from Matt’s hotspots paper, not sure if you need this? Either way, I’ve attached that as well.
# The link from which to download the data: http://mqtldb.godmc.org.uk/downloads

## Short script to match hvCpGs to controls 
library(dplyr)
library(tidyverse)
library(meffil)

# Julius homemade function:
'%!in%' <- function(x,y)!('%in%'(x,y))

setwd("/home/alice/2024_hvCpG/03_prepDatasetsMaria/") # set path

## 🔍 1. Load hvCpG (Maria’s hyper-variable CpGs)
hvCpGs_names.vec <- read.csv("~/2024_hvCpG/03_prepDatasetsMaria/Derakhshan2022_ST5_hvCpG.txt") %>%
  dplyr::select(cpg_name = CpG) %>%
  pull(cpg_name)
length(hvCpGs_names.vec) # 4143

# 🧬 2. Load GoDMC summary‐level mQTL data
# (http://mqtldb.godmc.org.uk/downloads) --> dl then rm after script run because large file
## NB: this is based on hg19

# ## Done in bash to avoid loading massive file in R:
# zcat assoc_meta_all.csv.gz | awk -F, '
# NR==1 {
#   # Build header index map
#   for (i=1; i<=NF; i++) {
#     h = $i
#     gsub(/^"|"$/, "", h)
#     hdr[h] = i
#   }
#   print
#   next
# }
# {
#   clumped = $hdr["clumped"]
#   cistrans = $hdr["cistrans"]
#   pval     = $hdr["pval"]
#   # Remove quotes
#   gsub(/^"|"$/, "", clumped)
#   gsub(/^"|"$/, "", cistrans)
#   gsub(/^"|"$/, "", pval)
#   p = pval + 0  # force numeric
#   
#   if (clumped == "TRUE" &&
#      ((cistrans == "TRUE"  && p < 1e-8) ||
#       (cistrans == "FALSE" && p < 1e-14)))
#     print
# }' > filtered_goDMC.csv
#   
filtered_goDMC <- read.csv("filtered_goDMC.csv")

table(hvCpGs_names.vec %in% filtered_goDMC$cpg)

GoDMC.df.FULL <- filtered_goDMC
nrow(GoDMC.df.FULL) # 271724

# 📊 3. Compute variance explained per mQTL
# Formula for variance explained: # 2 * MAF * (1 - MAF) * (beta^2) 
# Labels each as “cis” or “trans"
names(GoDMC.df.FULL)

GoDMC.df <- GoDMC.df.FULL %>%  
  dplyr::select(cpg_name=cpg, snp, beta_a1, pval, allele1, allele2, freq_a1, cistrans) %>% 
  mutate(MAF = ifelse(freq_a1 < 0.5, freq_a1, 1-freq_a1)) %>%
  mutate(var_explained = 2 * MAF * (1 - MAF) * (beta_a1 ^ 2)) %>% 
  mutate(cistrans = ifelse(cistrans == F, "trans", "cis")) %>% 
  dplyr::select(cpg_name, snp, pval, var_explained, cistrans) 

# 🧮 4. Count QTL frequency per CpG
# Here calculate the frequency of mQTL associations for each CpG site
snp_freq_per_cpg.df <- GoDMC.df %>%
  dplyr::group_by(cpg_name, cistrans) %>%
  dplyr::summarise(snp_freq = n_distinct(snp)) %>%
  ungroup() %>% 
  pivot_wider(names_from = cistrans, values_from = snp_freq, values_fill = list(snp_freq = 0)) %>% # if no snps for cis or trans make frequency 0
  dplyr::rename(cis_snp_freq = cis, trans_snp_freq = trans)

# 🎯 5. Subset data for hvCpGs
hvCpG_sumstats.df <- GoDMC.df %>% 
  filter(cpg_name %in% unique(hvCpGs_names.vec))

length(unique(hvCpG_sumstats.df$cpg_name)) # 3722, like in Maria's paper = all good

# 📏 6. Identify “potential controls”
# here filter the GoDMC.df for controls that are not in the same cluster as hvCpGs, at least 10 kb away 

#meffil::meffil.list.featuresets()
annot_450k.df <- meffil::meffil.get.features("450k") %>% 
  dplyr::rename(cpg_name = name) %>% 
  dplyr::select(cpg_name, cpg_chr = chromosome, cpg_pos = position) %>% 
  mutate(cpg_chr = as.character(cpg_chr), 
         cpg_pos = as.numeric(cpg_pos))

GoDMC.df <- GoDMC.df %>% 
  left_join(., annot_450k.df, by = "cpg_name")

# check that all positions annotated - nrow should be 0, all correct
GoDMC.df %>% filter(is.na(cpg_pos)) %>% nrow(.)
GoDMC.df %>% filter(is.na(cpg_chr))  %>% nrow(.)

GoDMC.df$cpg_clust <- bumphunter::clusterMaker(chr =GoDMC.df$cpg_chr, pos = GoDMC.df$cpg_pos,
                                               maxGap = 1e4, assumeSorted = F)

#filter out hvCpG clusters 
hvCpG_clust.vec <- GoDMC.df %>% 
  filter(cpg_name %in% hvCpGs_names.vec) %>%  
  distinct(cpg_clust) %>% 
  pull(cpg_clust)

## Here create dataframe of potential controls. 
potential_controls.df <- GoDMC.df %>% 
  filter(cpg_clust %!in% hvCpG_clust.vec) 

length(unique(potential_controls.df$cpg_name)) # check number of potential controls: 154452
length(intersect(potential_controls.df$cpg_name, hvCpGs_names.vec)) # no hvCpGs intersect potential controls 

# 🔄 7. Join SNP frequency & variance info
# Both hvCpGs and potential controls get their cis/trans SNP frequency and mean variance explained.
hvCpG_sumstats.df <- hvCpG_sumstats.df %>% 
  inner_join(., snp_freq_per_cpg.df, by = "cpg_name") 

potential_controls.df <- potential_controls.df %>% 
  inner_join(., snp_freq_per_cpg.df, by = "cpg_name") 

## !!NEW!! verify that the pote

# 🤝 8. Match each hvCpG to controls
# For each hvCpG find eligible controls with the same cis and trans snp frequency - useful to check number of eligible controls based on 
# snp frequency for each hvCpG. 
# For each hvCpG:
# 1. Filter controls that match exactly on cis/trans SNP counts,
# 2. Further narrow to controls with mean variance explained (cis & trans) within ±0.005 of the hvCpG,
# 3. If multiple are eligible, choose the one minimizing the total variance difference,
# 4. Ensure each control is used only once (no repeats).

# this will take a couple of seconds, can speed this up if want e.g. mclapply

hvCpG_sumstats.list <- hvCpG_sumstats.df %>% 
  split(., .$cpg_name) 

hvCpG_w_match_control_snp_freq.list <- lapply(hvCpG_sumstats.list, function(me_input.df){
  hvCpG_select_snp_freq.df <- me_input.df %>% 
    distinct(cpg_name, cis_snp_freq, trans_snp_freq)
  
  # select controls with the same snp frequencies 
  control_match_snp_freq.vec <- potential_controls.df %>% 
    distinct(cpg_name, cis_snp_freq, trans_snp_freq) %>% 
    filter(cis_snp_freq == hvCpG_select_snp_freq.df$cis_snp_freq & 
             trans_snp_freq == hvCpG_select_snp_freq.df$trans_snp_freq) %>% 
    distinct(cpg_name) %>% 
    pull(cpg_name)
}) 

## Check which hvcpgs dont have matching control based on mqtl frequency
hvCpG_nomatch.list <- sapply(unique(names(hvCpG_w_match_control_snp_freq.list)), function(hvCpG_name.value){
  length(unique(hvCpG_w_match_control_snp_freq.list[[hvCpG_name.value]]))
})

hvCpG_nomatch_names.vec <- names( which(hvCpG_nomatch.list == 0))
length(hvCpG_nomatch_names.vec) ## 1

# filter out hvCpG without control 
hvCpG_sumstats.list <- hvCpG_sumstats.list[which(names(hvCpG_sumstats.list) %!in% hvCpG_nomatch_names.vec)]

## now calculate the difference in cis and trans variacne explained 
# function to calculate the mean variance explained in cis and in trans

hvCpG_var_explained.mat <- hvCpG_sumstats.df %>% 
  group_by(cpg_name, cistrans) %>% 
  dplyr::summarize(meanvar = mean(var_explained, na.rm = T)) %>% 
  ungroup()  %>%
  complete(cpg_name, cistrans = c("cis", "trans"), 
           fill = list(meanvar = 0)) %>% # if no snps for cis or trans make mean variance 0
  pivot_wider(names_from = cpg_name, values_from = meanvar) %>% 
  column_to_rownames(., var = "cistrans") %>% 
  as.matrix() 

controls_var_explained.mat <- potential_controls.df %>% 
  group_by(cpg_name, cistrans) %>% 
  dplyr::summarize(meanvar = mean(var_explained, na.rm = T)) %>% 
  ungroup()  %>%
  complete(cpg_name, cistrans = c("cis", "trans"), 
           fill = list(meanvar = 0)) %>% # if no snps for cis or trans make mean variance 0
  pivot_wider(names_from = cpg_name, values_from = meanvar)  %>% 
  column_to_rownames(., var = "cistrans") %>% 
  as.matrix() 

# make vector of controls already selected
already_selected_controls.vec <- rep(NA, length(hvCpG_sumstats.list ))

# calculate number of eligible controls for each hvCpG
hvCpG_eligibility.df <- data.frame(
  hvCpG_name = names(hvCpG_sumstats.list),
  num_eligible_controls = sapply(names(hvCpG_sumstats.list), function(hvCpG_name.value) {
    length(unique(hvCpG_w_match_control_snp_freq.list[[hvCpG_name.value]]))
  })
)

# sort hvCpGs by increasing number of eligible controls
sorted_mes.vec <- hvCpG_eligibility.df %>%
  arrange(num_eligible_controls) %>%
  pull(hvCpG_name)

# iterate to match hvCpGs - can take a minute 
hvCpG_w_matched_control.df <- sapply(sorted_mes.vec, function(hvCpG_name.value) {
  me_input.mat <- hvCpG_var_explained.mat[, hvCpG_name.value, drop = F]
  eligible_control_names.vec <- hvCpG_w_match_control_snp_freq.list[[hvCpG_name.value]]
  eligible_control_names.vec <- eligible_control_names.vec[which(eligible_control_names.vec %!in% already_selected_controls.vec)] # dont select already selected controls 
  
  eligible_controls.mat <- controls_var_explained.mat[, eligible_control_names.vec, drop = F]
  
  var_explained_threshold.value <- 5e-3 # threshold for mean variance explained
  lower_bounds <- me_input.mat - var_explained_threshold.value
  upper_bounds <- me_input.mat + var_explained_threshold.value
  
  # Repeat bounds across columns
  lower_bounds.mat <- matrix(rep(lower_bounds, ncol(eligible_controls.mat)), 
                             nrow = 2, byrow = FALSE)
  upper_bounds.mat <- matrix(rep(upper_bounds, ncol(eligible_controls.mat)), 
                             nrow = 2, byrow = FALSE)
  
  
  # matrix for TRUE/FALSE if within bounds 
  within_bounds.mat <- eligible_controls.mat >= lower_bounds.mat & 
    eligible_controls.mat <= upper_bounds.mat
  
  # Filter: both cis and trans must be within bounds
  cols_to_keep <- colSums(within_bounds.mat) == 2
  
  filtered_controls.mat <- eligible_controls.mat[, cols_to_keep, drop = FALSE]
  
  # get the controls
  if (ncol(filtered_controls.mat) == 0) {
    selected_control.value <- NA  # No eligible control
  } else if (ncol(filtered_controls.mat) == 1) {
    selected_control.value <- colnames(filtered_controls.mat)[1]  # get name of only eligible control
  } else {
    # if multiple controls within threshold and calculate minimum difference 
    var_diff.mat <- apply(filtered_controls.mat, MARGIN = 2, function(input.col) {
      abs(input.col - me_input.mat)
    })
    
    # sum differences across cis and trans
    var_diff.vec <- colSums(var_diff.mat)
    
    # identify control with the minimum variance difference
    min_var_diff.cols <- names(var_diff.vec)[var_diff.vec == min(var_diff.vec)]
    
    # unlikely ot happen but if multiple controls with same difference in variance randomly choose 1
    selected_control.value <- if (length(min_var_diff.cols) > 1) {
      sample(min_var_diff.cols, 1)
    } else {
      min_var_diff.cols
    }
  }
  # add selected control to already_selected_controls.vec if it is not NA
  if (!is.na(selected_control.value)) {
    already_selected_controls.vec <<- c(already_selected_controls.vec, selected_control.value)
  }
  
  # output dataframe
  matched_control_output.df <- data.frame(
    hvCpG_name = hvCpG_name.value,
    controlCpG_name = selected_control.value
  )
  
  return(matched_control_output.df)
  
}, simplify = F, USE.NAMES = T)   %>%
  data.table::rbindlist(., idcol = F) %>%
  as.data.frame()

nr_hvCpGs_prematch.value <- nrow(hvCpG_w_matched_control.df)

# below: filter out hvCpGs that could not be matched 
hvCpG_w_matched_control.df <- hvCpG_w_matched_control.df %>% 
  filter(!is.na(controlCpG_name))

nr_hvCpGs_postmatch.value <- nrow(hvCpG_w_matched_control.df)

# get number of hvCpGs that cannot be matched 
hvCpG_loss.value <- nr_hvCpGs_postmatch.value - nr_hvCpGs_prematch.value
hvCpG_loss.value # 77

# check how many hvCpGs could be matched and number of controls - should be equal : 3644
hvCpG_w_matched_control.df %>% 
  summarise(hvCpG_freq = n_distinct(hvCpG_name), 
            controlCpG_freq = n_distinct(controlCpG_name))

# check the mqtl frequencies match 
snp_per_cpg_controlmatch.df <- GoDMC.df %>% 
  dplyr::filter(cpg_name %in% c(hvCpG_w_matched_control.df$hvCpG_name,
                         hvCpG_w_matched_control.df$controlCpG_name)) %>%
  group_by(cpg_name, cistrans) %>%
  dplyr::summarise(snp_freq = n_distinct(snp)) %>%
  ungroup() %>% 
  pivot_wider(names_from = cistrans, values_from = snp_freq, values_fill = list(snp_freq = 0)) %>% 
  dplyr::rename(cis_snp_freq = cis, trans_snp_freq = trans) %>% 
  mutate(hvCpG_or_control = ifelse(cpg_name %in% hvCpG_w_matched_control.df$hvCpG_name, "hvCpG", 
                                   "controlCpG"), 
         total_snp_freq = cis_snp_freq + trans_snp_freq) %>% 
  pivot_longer(., cols = c("cis_snp_freq", "trans_snp_freq", "total_snp_freq"),
               names_to = "snp_freq_cat", values_to = "freq")  %>% 
  mutate(facet_label = case_when(snp_freq_cat == "cis_snp_freq" ~ "Cis", 
                                 snp_freq_cat == "trans_snp_freq" ~ "Trans", 
                                 snp_freq_cat == "total_snp_freq" ~  "Total", 
                                 TRUE ~ NA_character_), 
         facet_label = factor(facet_label, 
                              levels = c("Total", "Cis", "Trans")))

# 📊 9. Validate matching

snp_per_cpg_controlmatch.df %>% 
  ggplot() + 
  geom_histogram(aes(x = freq, fill = hvCpG_or_control), position = "dodge") + 
  facet_wrap(~ facet_label) +
  labs(
    title = "SNP frequency of hvCpGs and matched controls",
    x = "SNP Frequency",
    y = "CpG frequency",
    fill = "CpG Category"
  )

# check that the cis and trans variances match --> it does!

cistrans_meanvar.df <- GoDMC.df %>% 
  dplyr::filter(cpg_name %in% c(hvCpG_w_matched_control.df$hvCpG_name,
                         hvCpG_w_matched_control.df$controlCpG_name)) %>%
  group_by(cpg_name, cistrans) %>%
  dplyr::summarise(mean_var_explained = mean(var_explained, na.rm = T)) %>%
  ungroup() %>% 
  pivot_wider(names_from = cistrans, values_from = mean_var_explained, values_fill = list(mean_var_explained = 0)) %>% 
  dplyr::rename(cis_meanvar = cis, trans_meanvar = trans) 

all_meanvar.df <- GoDMC.df %>% 
  dplyr::filter(cpg_name %in% c(hvCpG_w_matched_control.df$hvCpG_name,
                         hvCpG_w_matched_control.df$controlCpG_name)) %>%
  group_by(cpg_name) %>%
  dplyr::summarise(all_meanvar_explained = mean(var_explained, na.rm = T)) %>%
  ungroup() %>% 
  left_join(., cistrans_meanvar.df, by = "cpg_name") %>% 
  mutate(hvCpG_or_control = ifelse(cpg_name %in% hvCpG_w_matched_control.df$hvCpG_name, "hvCpG", 
                                   "controlCpG")) %>% 
  pivot_longer(., cols = c("cis_meanvar", "trans_meanvar", "all_meanvar_explained"),
               names_to = "cis_trans_var_explained", 
               values_to = "var_explained") %>% 
  mutate(facet_label = case_when(cis_trans_var_explained == "cis_meanvar" ~ "Cis", 
                                 cis_trans_var_explained == "trans_meanvar" ~ "Trans", 
                                 cis_trans_var_explained == "all_meanvar_explained" ~  "Total", 
                                 TRUE ~ NA_character_), 
         facet_label = factor(facet_label, 
                              levels = c("Total", "Cis", "Trans")))


# density plot
all_meanvar.df %>% 
  ggplot() + 
  geom_density(aes(x = var_explained, fill = hvCpG_or_control),  
               alpha = 0.3) + 
  facet_wrap(~ facet_label, scales = "free") +
  labs(
    title = "Mean variance explained of hvCpGs and matched controls",
    x = "Mean variance explained",
    y = "CpG frequency",
    fill = "CpG Category"
  )

# histogram 
all_meanvar.df %>% 
  ggplot() + 
  geom_histogram(aes(x = var_explained, fill = hvCpG_or_control),  
                 alpha = 1, position = "dodge") + 
  facet_wrap(~ facet_label, scales = "free") +
  labs(
    title = "Mean variance explained of hvCpGs and matched controls",
    x = "Mean variance explained",
    y = "CpG frequency",
    fill = "CpG Category"
  )

nrow(hvCpG_w_matched_control.df) ## 3644 match/control hvCpG from Maria

#save the dataframe of hvCpGs and matched controls 
write_delim(hvCpG_w_matched_control.df, file = "/home/alice/2024_hvCpG/03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt",col_names = T,quote = "needed")
