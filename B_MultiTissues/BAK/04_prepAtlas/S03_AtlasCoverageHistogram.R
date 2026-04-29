library(ggplot2)

hist_data <- read.table("../04_prepAtlas/coverage_histogram.tsv", header = TRUE, sep = "\t")

## Compute the median
weighted_median <- function(x, w) {
  # force weights to numeric to avoid integer overflow
  w <- as.numeric(w)
  
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cum_w <- cumsum(w) / sum(w)
  x[which(cum_w >= 0.5)[1]]
}

median_val <- weighted_median(hist_data$coverage, hist_data$frequency)
median_val # 39

ggplot(hist_data, aes(x = coverage, y = frequency)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Coverage histogram across beta files in Atlas data",
    x = "Coverage",
    y = "Frequency"
  ) + 
  scale_x_continuous(breaks=seq(0, 300, 10))+
  geom_vline(xintercept = median_val, col = "red")+
  theme_minimal()
