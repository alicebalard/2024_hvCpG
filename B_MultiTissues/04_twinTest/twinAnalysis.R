library(data.table)
library(ggplot2)

## ---------------- parameters ---------------------
tabRDS  <- "/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/table3layers_coveredIn3dt.RDS"
covdir  <- "/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01Methcall/Blood"  # blood twin .cov(.gz)
covpat  <- "\\.cov(\\.gz)?$"
MINCOV  <- 10L     # "sufficient" per-sample read depth at a CpG
MINSAMP <- NA      # min #samples that must reach MINCOV; NA = require ALL samples

## ---------------- (1) top 1% CpGs ----------------
dt  <- readRDS(tabRDS)

top <- data.table(chr_pos = dt$chr_pos[dt$percentile_logBF_per_ds_allLayers >= 99])
setkey(top, chr_pos)
cat("top-1% CpGs:", nrow(top), "\n")          # expect 202467

## ---------------- (2) coverage at those CpGs in each blood twin ----------------
files <- list.files(covdir, pattern = covpat, full.names = TRUE)
stopifnot(length(files) > 0)
nsamp <- length(files)
cat("blood samples:", nsamp, "\n")

read_cov_at_top <- function(f) {
    x <- fread(f, header = FALSE,
               col.names = c("chr","start","end","meth_pct","M","U"))
    x[, chr_pos := paste0(chr, "_", start)]
    x[, cov := M + U]
    x[top, on = "chr_pos", nomatch = 0L][        # keep only the target CpGs
      , .(chr_pos, cov, sample = sub(covpat, "", basename(f)))]
}
covlong <- rbindlist(lapply(files, read_cov_at_top))

## ---------------- per-CpG summary (absent sample = 0x) ----------------
perCpG <- covlong[, .(n_samples_present = .N,
                      n_cov_ge_min      = sum(cov >= MINCOV),
                      total_cov         = sum(cov),
                      median_cov_present= as.numeric(median(cov))),
                  by = chr_pos]
perCpG <- perCpG[top, on = "chr_pos"]                      # add CpGs seen in no file
perCpG[is.na(n_samples_present),
       `:=`(n_samples_present = 0L, n_cov_ge_min = 0L,
            total_cov = 0, median_cov_present = 0)]
perCpG[, mean_cov := total_cov / nsamp]                    # mean depth across ALL samples

## "sufficient": reaches MINCOV in enough samples
need <- if (is.na(MINSAMP)) nsamp else MINSAMP
perCpG[, sufficient := n_cov_ge_min >= need]

## ---------------- answers ----------------
cat("\nCpGs with sufficient coverage (>=", MINCOV, "x in >=", need, "samples): ",
    perCpG[sufficient == TRUE, .N], " / ", nrow(perCpG), "\n", sep = "")
cat("match rate (CpGs seen in >=1 blood file): ",
    round(100 * perCpG[n_samples_present > 0, .N] / nrow(perCpG), 1), "%\n", sep = "")
cat("\nper-CpG mean coverage across samples:\n")
cat("  average:", round(mean(perCpG$mean_cov), 2),
    " median:", round(median(perCpG$mean_cov), 2), "\n")

## ---------------- histogram ----------------
p <- ggplot(perCpG, aes(mean_cov)) +
    geom_histogram(bins = 60, fill = "#2166AC", colour = "white") +
    geom_vline(xintercept = MINCOV, linetype = 2, colour = "#DC3220") +
    scale_x_continuous(limits = c(0, quantile(perCpG$mean_cov, 0.99))) +  # trim tail
    labs(x = "mean read depth across blood twin samples",
         y = "number of top-1% CpGs",
         title = "Coverage of top-1% hypervariable CpGs in blood twins") +
    theme_minimal(base_size = 13)
ggsave("blood_twin_coverage_top1pct.png", p, width = 8, height = 5, dpi = 300, bg = "white")

