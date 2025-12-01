#######################################################
## Plot results of algorithm ran at the tissue level ##
#######################################################
library(here)
source(here("05_hvCpGalgorithm", "quiet_library.R"))

parent_dir <- here("gitignore/Atlas10X_tissueAnalysis_NONGITED/")

rds_files <- list.files(parent_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

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
      mat <- mat[sub("-.*", "", rownames(mat)) %in% CpGs,]
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

summary_dt <- makeSummaryDT(rds_files)

p0 <- ggplot(summary_dt, aes(x = dataset, y = meanexp)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lowerexp, ymax = ci_upperexp)) +
  theme_minimal(base_size = 14) +
  labs(title = "Mean probability of being hypervariable, with 95% CI", x = "Dataset", y = "Mean p(hv)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Statistical differences between samples?
kruskal.test(mean ~ dataset, data = summary_dt)

## Statistical differences between germ layers in terms of variability?

## Differneces between germ layers
kruskal.test(mean ~ Germ.layer, data = summary_dt)

## Blood vs rest
kruskal.test(mean ~ isBlood, data = summary_dt)

####################################################################################
## That's for all CpGs. What about looking at the top 100k and bottom 100k alpha? ##
####################################################################################
top100k <- readRDS(here("05_hvCpGalgorithm/dataOut/top100k_4Hamdan.rds"))
bottom100k <- readRDS(here("05_hvCpGalgorithm/dataOut/bottom100k_4Hamdan.rds"))

summary_dt_top100 <- makeSummaryDT(rds_files, selectSomeCpGs = TRUE, CpGs = top100k$name)
summary_dt_bottom100 <- makeSummaryDT(rds_files, selectSomeCpGs = TRUE, CpGs = bottom100k$name)

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

