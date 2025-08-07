source("../05_hvCpGalgorithm/hvCpG_algorithm_detection_v4scan.R")

library(patchwork)

exploreAlgo <- function(x,title){
  alpha_trace <<- list()
  
  runAndSave(
    analysis = "Maria",
    cpgPos_vec = x,
    resultDir = "~/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Mariads/RES/",
    NCORES = 1,
    p0 = 0.80,
    p1 = 0.65,
    overwrite = TRUE
  ) 
  
  load("~/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Mariads/RES/results_Maria_1CpGs_0_8p0_0_65p1.RData")
  message("alpha:")
  resalgo = results_Maria_1CpGs_0_8p0_0_65p1
  print(resalgo)

  unlog <- function(x) {
    odds <- 2^x
    beta <- odds / (1 + odds)
    beta
  }
  
  ## Prepare data in the environment:
  prep = prepData("Maria")
  metadata = prep$metadata
  medsd_lambdas = prep$medsd_lambdas
  cpg_names_all = prep$cpg_names_all
  source_M_1CpG = prep$source_M_1CpG

  cpgRaw = unlog(source_M_1CpG(x))
  cpgRaw <- as.data.frame(cpgRaw)
  cpgRaw$sample <- rownames(cpgRaw)

  # Reshape to long format: sample, CpG, value
  long_df <- reshape2::melt(cpgRaw, id.vars = "sample", variable.name = "CpG", value.name = "value")
  # Merge with metadata to get dataset info
  long_df <- left_join(long_df, metadata, by = "sample")

  # Now plot: dataset on x, value on y
  p1 <- ggplot(long_df, aes(x = dataset, y = value)) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.7, size = 2) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Dataset") +
    theme(axis.text.x = element_text(size =6)) +
    facet_wrap(~ CpG, scales = "free_y", ncol = 3)

  sd_long_df <- long_df %>% group_by(dataset, CpG) %>% summarise(sdMethyl = sd(value, na.rm = T))
  sd_long_df <- merge(sd_long_df, medsd_lambdas)
  sd_long_df$thr <- sd_long_df$median_sd / sd_long_df$lambda

  sd_long_df <- sd_long_df[sd_long_df$dataset %in% long_df$dataset,]

  p2 <- ggplot(sd_long_df, aes(x = dataset)) +
    geom_segment(aes(
      x = dataset,
      xend = dataset,
      y = thr,
      yend = sdMethyl,
      color = sdMethyl > thr
    ),
    arrow = arrow(length = unit(0.15, "cm")),
    position = position_jitter(width = 0.2)) +

    scale_color_manual(values = c("FALSE" = "red", "TRUE" = "green")) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    xlab("Dataset") +
    coord_cartesian(ylim = c(0, 1)) +
    facet_wrap(~ CpG, scales = "free_y", ncol = 3) +
    guides(color = "none") + # Remove legend if desired
    ggtitle(paste0(title, "\n alpha = ", round(as.numeric(resalgo), 3)))

  alpha_df <- do.call(rbind, lapply(alpha_trace, as.data.frame))
  p3 <- ggplot(alpha_df, aes(x = alpha, y = loglik)) +
    geom_point() +
    geom_line() +
    theme_minimal() +
    labs(title = "Alpha Optimization Trace",
         x = expression(alpha),
         y = "Log-likelihood")

  layout <- (p1 / p3) | p2
  pdf("figures/testMaria1hvCpGexplo.pdf", width = 10, height = 10)
  layout
  dev.off()
  print(table(sd_long_df$sdMethyl > sd_long_df$thr))
}

# a hvCpG:
exploreAlgo(x = match("cg03396347",  cpg_names_all), title = "a hvCpG from Maria")
