#######################################################
## Plot results of algorithm ran at the tissue level ##
#######################################################
library(here)
source(here("05_hvCpGalgorithm", "quiet_library.R"))

## Atlas
parent_dir_atlas <- here("05_hvCpGalgorithm/resultsDir/Atlas/Atlas10X_tissueAnalysis/")
rds_files_atlas <- list.files(parent_dir_atlas, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
all_medsd_lambda_atlas <- read.table(here("04_prepAtlas/all_medsd_lambda.tsv"), sep = "\t", header = T)

## Array
parent_dir_array <- here("05_hvCpGalgorithm/resultsDir/Arrays/tissue/")
rds_files_array <- list.files(parent_dir_array, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
all_medsd_lambda_array <- read.table(here("03_prepDatasetsMaria/all_medsd_lambda.tsv"), sep = "\t", header = T)

## Both
WGBS_Array_datasets <- read.csv(here("05_hvCpGalgorithm/figures/WGBS_Array_datasets.csv"))

#################
## Violin plot ##
#################

makeViolin <- function(size = 10, rds_files, fill = "nothing", all_medsd_lambda, mytitle){
  set.seed(1234)
  mat1 = readRDS(rds_files[1])
  mat1 = mat1[sample(nrow(mat1), size = size), ,drop = F]
  pb <- progress_bar$new(total = length(rds_files), format = "📦 :current/:total [:bar] :percent")
  for (file in rds_files) {
    mat = readRDS(file)  # rows ~100000, cols ~46
    ## Select a sub sample because it's too big otherwise
    mat1 = rbind(mat1, mat[sample(nrow(mat), size = size), ,drop = F])
    pb$tick()
  }
  mat1 = reshape2::melt(mat1)
  names(mat1) = c("pos", "dataset", "prob")
  
  ## Exponentialise because it is on a log scale
  mat1$prob <- exp(mat1$prob)
  
  ## Add germ layer
  if (fill == "Germ layer"){
    SupTab1_Loyfer2023 = read.csv(here("05_hvCpGalgorithm/dataPrev/SupTab1_Loyfer2023.csv"))
    mat1$Germ.layer = SupTab1_Loyfer2023$Germ.layer[
      match(mat1$dataset, paste0(SupTab1_Loyfer2023$Source.Tissue, " - ", SupTab1_Loyfer2023$Cell.type))]
  }
  
  ## add median SD and lambda
  mat1 = merge(mat1, all_medsd_lambda)
  
  ## Add mitosis
  mat1$mitosis_rate = WGBS_Array_datasets$mitosis_rate[match(mat1$dataset, WGBS_Array_datasets$dataset)]
  
  ## Plot on 2 axes p(hv) and lambda
  range_prob = range(mat1$prob, na.rm = TRUE)
  range_lambda = range(mat1$lambda, na.rm = TRUE)
  
  # slope and intercept for mapping
  a = diff(range_lambda) / diff(range_prob)
  b = range_lambda[1] - a * range_prob[1]
  
  # Map lambda to the primary y-scale
  mat1$lambda_y <- (mat1$lambda - b) / a
  
  # Reorder factor levels ascending by median prob 
  mat1 <- mat1 %>%
    mutate(dataset = forcats::fct_reorder(dataset, lambda, .desc = FALSE, .na_rm = TRUE))
  
  p = ggplot(mat1, aes(x = dataset, y = prob))
  
  if (fill == "Germ layer"){
    p = p + geom_violin(aes(fill = Germ.layer), width = 1.2)
    summaryMat1 = mat1 %>% 
      dplyr::group_by(dataset, Germ.layer) %>%
      dplyr::summarise(medianPrhv = median(prob, na.rm = T),
                       lambda = median(lambda, na.rm = T)) # all the same
  } else if (fill == "mitosis_rate"){
    p = p + geom_violin(aes(fill = mitosis_rate), width = 1.2)
    summaryMat1 = mat1 %>% 
      dplyr::group_by(dataset, mitosis_rate) %>%
      dplyr::summarise(medianPrhv = median(prob, na.rm = T),
                       lambda = median(lambda, na.rm = T)) # all the same
  } else {
    p = p + geom_violin(width = 1.2)
    summaryMat1 = mat1 %>% 
      dplyr::group_by(dataset) %>%
      dplyr::summarise(medianPrhv = median(prob, na.rm = T),
                       lambda = median(lambda, na.rm = T)) # all the same
  }
  
  p = p +
    scale_fill_viridis_d() +
    geom_boxplot(width = .2, fill = "white", outliers = F) +
    geom_hline(yintercept = median(mat1$prob, na.rm = T), linetype = 3) +
    theme_minimal(base_size = 14) +
    labs(title = mytitle,
         x = NULL, y = "p(hv)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(
      name = "Probability of being hv",
      sec.axis = sec_axis(~ a * . + b, name = "lambda")
    ) +
    geom_point(aes(y = lambda_y), color = "red", size = 3, show.legend = FALSE) +
    geom_hline(aes(yintercept = median(lambda_y, na.rm = T)), linetype = 3, col = "red") # median lambda
  
  return(list(p=p, summaryMat1=summaryMat1))
}

#########################
## mitosis_rate effect ##
#########################

## Atlas
plot_violin_tissues_atlas <- makeViolin(
  rds_files = rds_files_atlas, size = 10, fill = "mitosis_rate",
  all_medsd_lambda = all_medsd_lambda_atlas, mytitle = "atlas")

## Array
plot_violin_tissues_array <- makeViolin(
  rds_files = rds_files_array, size = 10, fill = "mitosis_rate",
  all_medsd_lambda = all_medsd_lambda_array, mytitle = "arrays")

pdf(here("05_hvCpGalgorithm/figures/tissuePlot_arrayAtlas_mitosis.pdf"), width = 15, height = 12)
cowplot::plot_grid(plot_violin_tissues_atlas$p, plot_violin_tissues_array$p, nrow = 2)
dev.off()



# library(tidyverse)
# library(rstatix)
# library(ggpubr)

plot_violin_tissues_atlas$summaryMat1 %>%
  group_by(mitosis_rate) %>%
  rstatix::get_summary_stats(lambda, type = "median_iqr")

ggboxplot(
  plot_violin_tissues_atlas$summaryMat1, x = "mitosis_rate", y = "lambda", 
  ylab = "lambda", xlab = "mitosis_rate", add = "jitter"
)

wilcox.test(plot_violin_tissues_atlas$summaryMat1$mitosis_rate,
            plot_violin_tissues_atlas$summaryMat1$lamda)











# A high lambda value indicates that a dataset contains a pronounced subset of CpGs with exceptionally high inter‑sample methylation variability relative to the typical CpG, pointing to locus‑specific biological heterogeneity (or, less desirably, technical outliers) rather than uniform variability across the methylome.

## Atlas
p1 <- makeViolin(
  rds_files = rds_files_atlas, size = 10, fill = "Germ layer",
  all_medsd_lambda = all_medsd_lambda_atlas, mytitle = "atlas")

## Array
p2 <- makeViolin(
  rds_files = rds_files_array, size = 10,
  all_medsd_lambda = all_medsd_lambda_array, mytitle = "arrays")

pdf(here("05_hvCpGalgorithm/figures/tissuePlot_arrayAtlas_germLayer.pdf"), width = 15, height = 12)
cowplot::plot_grid(p1$p, p2$p, nrow = 2)
dev.off()

plot_violin_tissues_atlas$summaryMat1$method = "WGBS atlas of purified cells"
plot_violin_tissues_array$summaryMat1$method = "methylation array of tissues"

sumMatlambdaprhv <- rbind(plot_violin_tissues_array$summaryMat1, plot_violin_tissues_atlas$summaryMat1)
print(sumMatlambdaprhv, n = 100)

# Compute residuals from linear model
model <- lm(lambda ~ medianPrhv, data = sumMatlambdaprhv)
sumMatlambdaprhv$residual <- residuals(model)

# Select top 10 positive and negative residuals
clean_df <- as.data.frame(lapply(sumMatlambdaprhv, as.vector))

# Arrange by residual and take top 10
top_residuals <- clean_df %>%
  arrange(desc(abs(residual))) %>%
  dplyr::slice(1:10)

# Plot
ggplot(sumMatlambdaprhv, aes(x = medianPrhv, y = lambda)) +
  geom_point(aes(col = mitosis), size = 3) +
  scale_color_viridis_d() +
  geom_abline(intercept = model$coefficients[1],
              slope = model$coefficients[2],
              color = "blue", linetype = "dashed")+
  geom_label_repel(data = top_residuals,
                   aes(label = dataset),
                   box.padding = 0.3,
                   max.overlaps = Inf) +
  theme_bw() +
  labs(title = "Residuals: medianPrhv vs lambda",
       x = "Median Pr(hv)",
       y = "lambda") +
  theme(legend.position = "top")

wilcox.test(sumMatlambdaprhv$medianPrhv[sumMatlambdaprhv$method %in% "WGBS atlas of purified cells"],
            sumMatlambdaprhv$medianPrhv[sumMatlambdaprhv$method %in% "methylation array of tissues"])

wilcox.test(sumMatlambdaprhv$lambda[sumMatlambdaprhv$method %in% "WGBS atlas of purified cells"],
            sumMatlambdaprhv$lambda[sumMatlambdaprhv$method %in% "methylation array of tissues"])
## non significant

p1 <- ggplot(sumMatlambdaprhv, aes(x = method, y = lambda)) +
  geom_violin() +
  geom_boxplot(width = .2) +
  geom_jitter(size = 3) +
  theme_bw() + labs(x = NULL)

p2 <- ggplot(sumMatlambdaprhv, aes(x = method, y = medianPrhv)) +
  geom_violin() +
  geom_boxplot(width = .2) +
  geom_jitter(size = 3) +
  theme_bw() + labs(x = NULL)

cowplot::plot_grid(p1,p2)

##################
## Summary plot ##
##################
makeSummaryDT <- function(rds_files, selectSomeCpGs = FALSE, CpGs = NA){
  
  n_files <- length(rds_files)
  
  # Initialize accumulators
  first_mat <- readRDS(rds_files[1])
  datasets <- colnames(first_mat)
  n_cols <- length(datasets)
  
  sum_vals <- rep(0, n_cols)
  sum_sq <- rep(0, n_cols)
  count_vals <- rep(0, n_cols)
  min_vals <- rep(Inf, n_cols)
  max_vals <- rep(-Inf, n_cols)
  
  pb <- progress_bar$new(total = n_files, format = "📦 :current/:total [:bar] :percent")
  
  for (file in rds_files) {
    mat <- readRDS(file)  # rows ~100000, cols ~46
    
    if (selectSomeCpGs){
      ## Select only some
      mat <- mat[sub("-.*", "", rownames(mat)) %in% CpGs, ,drop = F]
    }
    
    # Update accumulators
    sum_vals <- sum_vals + colSums(mat, na.rm = TRUE)
    sum_sq <- sum_sq + colSums(mat^2, na.rm = TRUE)
    count_vals <- count_vals + colCounts(!is.na(mat))
    
    min_vals <- pmin(min_vals, colMins(mat, na.rm = TRUE))
    max_vals <- pmax(max_vals, colMaxs(mat, na.rm = TRUE))
    
    pb$tick()
  }
  
  # Compute final stats
  means <- sum_vals / count_vals
  sds <- sqrt((sum_sq / count_vals) - (means^2))
  
  # Compute 95% CI
  alpha <- 0.05
  t_val <- qt(1 - alpha/2, df = count_vals - 1)  # t critical value
  ci_lower <- means - t_val * (sds / sqrt(count_vals))
  ci_upper <- means + t_val * (sds / sqrt(count_vals))
  
  summary_dt <- data.table(
    dataset = datasets,
    min = min_vals,
    max = max_vals,
    mean = means,
    sd = sds,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    meanexp = exp(means),
    sdexp = exp(sds),
    ci_lowerexp = exp(ci_lower),
    ci_upperexp = exp(ci_upper)
  )
  
  ## Add information
  SupTab1_Loyfer2023 = read.csv(here("05_hvCpGalgorithm/dataPrev/SupTab1_Loyfer2023.csv"))
  summary_dt$Germ.layer = SupTab1_Loyfer2023$Germ.layer[
    match(summary_dt$dataset, paste0(SupTab1_Loyfer2023$Source.Tissue, " - ", SupTab1_Loyfer2023$Cell.type))]
  
  summary_dt$isBlood <- ifelse(grepl("Blood", summary_dt$dataset), "Blood", "Other")
  
  return(summary_dt)
}

summary_dt <- makeSummaryDT(rds_files_atlas)

p0 <- ggplot(summary_dt, aes(x = dataset, y = meanexp)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lowerexp, ymax = ci_upperexp)) +
  theme_minimal(base_size = 14) +
  labs(title = "Mean probability of being hypervariable, with 95% CI", x = "Dataset", y = "Mean p(hv)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Statistical differences between samples?
kruskal.test(mean ~ dataset, data = summary_dt)

## Statistical differences between germ layers in terms of variability?

## Differences between germ layers
kruskal.test(mean ~ Germ.layer, data = summary_dt)

## Blood vs rest
kruskal.test(mean ~ isBlood, data = summary_dt)

####################################################################################
## That's for all CpGs. What about looking at the top 100k and bottom 100k alpha? ##
####################################################################################
top100k <- readRDS(here("05_hvCpGalgorithm/dataOut/top100k_4Hamdan.rds"))
bottom100k <- readRDS(here("05_hvCpGalgorithm/dataOut/bottom100k_4Hamdan.rds"))

summary_dt_top100 <- makeSummaryDT(rds_files_atlas, selectSomeCpGs = TRUE, CpGs = top100k$name)
summary_dt_bottom100 <- makeSummaryDT(rds_files_atlas, selectSomeCpGs = TRUE, CpGs = bottom100k$name)

summary_dt_all <- rbind(summary_dt %>% mutate(group="allCpGs"), 
                        summary_dt_top100 %>% mutate(group="top100kCpGs"), 
                        summary_dt_bottom100 %>% mutate(group="bottom100kCpGs"))

summary_dt_all <- summary_dt_all %>%
  group_by(group) %>%
  arrange(desc(mean), .by_group = TRUE) %>%
  mutate(top5 = row_number() <= 5) %>%
  arrange(mean, .by_group = TRUE) %>%
  mutate(bottom5 = row_number() <= 5)

summary_dt_all$grouptop <- ifelse(summary_dt_all$top5, "top5", ifelse(summary_dt_all$bottom5, "bottom5", NA))

ggplot(summary_dt_all, aes(x = dataset, y = meanexp, col = group)) +
  geom_point(aes(group=grouptop, fill = grouptop), col = "white", size = 5, pch = 21) +
  scale_fill_manual(values = c("blue", "red", "lightgrey")) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lowerexp, ymax = ci_upperexp)) +
  theme_minimal(base_size = 14) +
  labs(title = "Mean probability of being hypervariable, with 95% CI", x = "Dataset", y = "Mean p(hv)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_viridis_d()

## Statistical differences between samples?
kruskal.test(mean ~ dataset, data = summary_dt_top100)
kruskal.test(mean ~ dataset, data = summary_dt_bottom100)

## Statistical differences between germ layers in terms of variability?

## Differneces between germ layers
kruskal.test(mean ~ Germ.layer, data = summary_dt_top100)
kruskal.test(mean ~ Germ.layer, data = summary_dt_bottom100)

## Blood vs rest
kruskal.test(mean ~ isBlood, data = summary_dt_top100)
kruskal.test(mean ~ isBlood, data = summary_dt_bottom100)

ggplot(summary_dt_top100, aes(x = isBlood, y = mean, group = isBlood)) +
  geom_boxplot() + 
  geom_label_repel(data = summary_dt_top100[summary_dt_top100$mean > -0.15 | summary_dt_top100$mean < -0.30,], aes(label = dataset)) +
  geom_jitter(height = 0, width = .2) +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank())

######################################################################
## What about looking at the top alpha for array and low for atlas? ##
######################################################################

############### 95+50-
CpGArray95moreAtlas50less <- readRDS(here("05_hvCpGalgorithm/dataOut/CpGArray95moreAtlas50less.RDS"))$chrpos

summary_dt_CpGArray95moreAtlas50less <- makeSummaryDT(rds_files_atlas, selectSomeCpGs = TRUE, CpGs = CpGArray95moreAtlas50less)

ggplot(summary_dt_CpGArray95moreAtlas50less, aes(x = dataset, y = meanexp, col = Germ.layer)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lowerexp, ymax = ci_upperexp)) +
  theme_minimal(base_size = 14) +
  labs(title = "Mean probability of being hypervariable, with 95% CI", 
       subtitle = "Array95+Atlas50-", x = "Dataset", y = "Mean p(hv)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  scale_color_viridis_d()

## Differences between germ layers
kruskal.test(mean ~ Germ.layer, data = summary_dt_CpGArray95moreAtlas50less)

## Blood vs rest
kruskal.test(mean ~ isBlood, data = summary_dt_CpGArray95moreAtlas50less)

ggplot(summary_dt_CpGArray95moreAtlas50less, aes(x = isBlood, y = mean, group = isBlood)) +
  geom_boxplot() + 
  geom_jitter(height = 0, width = .2) +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank())

############### 95+50+ 
CpGArray95moreAtlas50more <- readRDS(here("05_hvCpGalgorithm/dataOut/CpGArray95moreAtlas50more.RDS"))$chrpos

summary_dt_CpGArray95moreAtlas50more <- makeSummaryDT(rds_files_atlas, selectSomeCpGs = TRUE, CpGs = CpGArray95moreAtlas50more)

ggplot(summary_dt_CpGArray95moreAtlas50more, aes(x = dataset, y = meanexp, col = Germ.layer)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lowerexp, ymax = ci_upperexp)) +
  theme_minimal(base_size = 14) +
  labs(title = "Mean probability of being hypervariable, with 95% CI", 
       subtitle = "Array95+Atlas50+", x = "Dataset", y = "Mean p(hv)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  scale_color_viridis_d()

## Statistical differences between germ layers in terms of variability?

## Differneces between germ layers
kruskal.test(mean ~ Germ.layer, data = summary_dt_CpGArray95moreAtlas50more)

## Blood vs rest
kruskal.test(mean ~ isBlood, data = summary_dt_CpGArray95moreAtlas50more)

ggplot(summary_dt_CpGArray95moreAtlas50more, aes(x = isBlood, y = mean, group = isBlood)) +
  geom_boxplot() + 
  geom_jitter(height = 0, width = .2) +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank())

############### 50+50+ 
CpGArray50moreAtlas50more <- readRDS(here("05_hvCpGalgorithm/dataOut/CpGArray50moreAtlas50more.RDS"))$chrpos

summary_dt_CpGArray50moreAtlas50more <- makeSummaryDT(rds_files_atlas, selectSomeCpGs = TRUE, CpGs = CpGArray50moreAtlas50more)

ggplot(summary_dt_CpGArray50moreAtlas50more, aes(x = dataset, y = meanexp, col = Germ.layer)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lowerexp, ymax = ci_upperexp)) +
  theme_minimal(base_size = 14) +
  labs(title = "Mean probability of being hypervariable, with 95% CI", 
       subtitle = "Array50+Atlas50+", x = "Dataset", y = "Mean p(hv)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  scale_color_viridis_d()
