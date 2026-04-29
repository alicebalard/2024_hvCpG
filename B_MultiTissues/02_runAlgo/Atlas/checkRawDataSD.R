checkRawDataSD <- function(cpgNamesVec){
  samples <- h5read("/home/alice/Documents/Project_hvCpG/10X/all_matrix_noscale.h5",
                    "samples")
  
  sample_groups <- h5read("/home/alice/Documents/Project_hvCpG/10X/all_matrix_noscale.h5",
                          "sample_groups")
  
  indexVec = which(cpg_names_all %in% cpgNamesVec)
  subRawData <- h5read("/home/alice/Documents/Project_hvCpG/10X/all_matrix_noscale.h5",
                       "matrix", index = list(NULL, indexVec))
  
  colnames(subRawData) = cpg_names_all[indexVec]
  subRawData <- as.data.frame(subRawData)
  subRawData$samples = samples
  subRawData$sample_groups = sample_groups
  
  # Exclude metadata columns
  num_cols <- setdiff(names(subRawData), c("samples", "sample_groups"))
  
  # Step 1: compute SD per group per CpG site
  df_long <- subRawData %>%
    select(all_of(num_cols), sample_groups) %>%
    pivot_longer(cols = -sample_groups, names_to = "CpG", values_to = "value") %>%
    group_by(sample_groups, CpG) %>%
    summarise(sd = sd(value, na.rm = TRUE), .groups = "drop")
  
  # Step 2: median of SDs per group
  med <- df_long %>%
    group_by(CpG) %>%
    summarise(median_sd = median(sd, na.rm = TRUE))
  
  # Step 3: histogram of SD distributions per group
  ggplot(df_long, aes(x = sd)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    facet_wrap(~ CpG, scales = "free_y") +
    theme_minimal() +
    labs(title = "Distribution of SD per CpG",
         x = "Standard deviation across samples",
         y = "Count")+
    geom_vline(data = med, aes(xintercept = median_sd))
}