#####################################################################
## S07_testTargetMatt.R
## Check Pr(HV) for Matt's target regions (LY6S-VMR/MER11C, ACTL8/LTR41)
## Produces: ACTL8_LTR41.pdf
#####################################################################

library(here)
source(here("B_MultiTissues", "quiet_library.R"))
if (!exists("functionsLoaded"))     source(here("B_MultiTissues/03_exploreResults", "functions.R"))
if (!exists("previousSIVprepared")) source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))
if (!exists("putativeME_GR"))       source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))
if (!exists("table3layers"))        load(here("gitignore/fullTable3layers_23_07_26.Rda"))
if (!exists("MEsetdt"))             MEsetdt            <- readRDS(here("gitignore/MEsetdt.rds"))
if (!exists("MEsetdt_regionMean"))  MEsetdt_regionMean <- readRDS(here("gitignore/MEsetdt_regionMean.rds"))
if (!exists("geomMeanGR"))          geomMeanGR         <- readRDS(here("gitignore/geomMeanGR.rds"))

if (is.null(table3layers$percentile)) {
  table3layers$percentile <- ecdf(table3layers$alpha_geomean)(table3layers$alpha_geomean) * 100
}

# ── Load inter-layer correlation from fetal EPIC data ────────────────────────
interlayer_corr <- readRDS(here("B_MultiTissues/dataOut/interlayer_corr_all.RDS"))
interlayer_corr$chr_pos <- dico$chrpos_hg38[match(interlayer_corr$CpG, dico$CpG)]
interlayer_corr$interlayer_r <- rowMeans(
  abs(interlayer_corr[, c("r_Endo_Meso", "r_Endo_Ecto", "r_Meso_Ecto")]),
  na.rm = TRUE)
interlayer_corr$percentile_r <- ecdf(interlayer_corr$interlayer_r)(
  interlayer_corr$interlayer_r) * 100

interlayer_corr_clean <- interlayer_corr[!is.na(interlayer_corr$chr_pos), ]
pos_check <- as.integer(sub(".*_", "", interlayer_corr_clean$chr_pos))
interlayer_corr_clean <- interlayer_corr_clean[!is.na(pos_check), ]

interlayer_gr <- GRanges(
  seqnames = sub("_.*", "", interlayer_corr_clean$chr_pos),
  ranges   = IRanges(
    start = as.integer(sub(".*_", "", interlayer_corr_clean$chr_pos)),
    width = 1))
mcols(interlayer_gr) <- interlayer_corr_clean

# ══════════════════════════════════════════════════════════════════════════════
# Genomic coordinates for Matt's target regions
# ══════════════════════════════════════════════════════════════════════════════

LY6SVMR_hg19  <- GRanges("chr8", IRanges(144120106, 144120706), name = "LY6S-VMR")
LY6SVMR_hg38  <- unlist(liftOver(LY6SVMR_hg19, chain))
LY6SVMR_hg38$name <- "LY6S-VMR"
LY6S_AS1_hg38 <- GRanges("chr8", IRanges(143039209, 143054303), name = "LY6S_AS1")
MER11C_hg38   <- GRanges("chr8", IRanges(143040739, 143041811), name = "MER11C")

LTR41_hg19 <- GRanges("chr1",
                       IRanges(start = c(18081648, 18085651),
                               end   = c(18082190, 18086109)),
                       name = c("LTR41_1", "LTR41_2"))
LTR41_hg38 <- unlist(liftOver(LTR41_hg19, chain))
LTR41_hg38$name <- c("LTR41_1", "LTR41_2")
ACTL8_hg38 <- GRanges("chr1", IRanges(17755333, 17827063), name = "ACTL8")

SUSPECT_ALICE <- c(GRanges("chr1", IRanges(17757500, 17757900), name = "SUSPECT_1"),
                   GRanges("chr1", IRanges(17760700, 17761100), name = "SUSPECT_2"))

# ── Region GRanges (plot windows) ────────────────────────────────────────────
region1_gr <- c(
  LY6SVMR_hg38, MER11C_hg38,
  GRanges("chr8",
          IRanges(start = min(start(LY6SVMR_hg38), start(MER11C_hg38)) - 3000,
                  end   = max(end(LY6SVMR_hg38),   end(MER11C_hg38))   + 3000),
          name = NA))

region2_gr <- c(
  LTR41_hg38, ACTL8_hg38, SUSPECT_ALICE,
  GRanges("chr1",
          IRanges(start = min(start(LTR41_hg38)) - 3000,
                  end   = max(start(LTR41_hg38)) + 3000),
          name = NA))

annot_gr_region1 <- GRanges(
  seqnames   = "chr8",
  ranges     = IRanges(
    start = c(start(LY6S_AS1_hg38), start(LY6SVMR_hg38), start(MER11C_hg38)),
    end   = c(end(LY6S_AS1_hg38),   end(LY6SVMR_hg38),   end(MER11C_hg38))),
  name       = c("LY6S_AS1", "LY6S-VMR", "MER11C"),
  annot_type = c("gene",     "geneVMR",  "TE"))

annot_gr_region2 <- GRanges(
  seqnames   = "chr1",
  ranges     = IRanges(
    start = c(start(LTR41_hg38), start(ACTL8_hg38), start(SUSPECT_ALICE)),
    end   = c(end(LTR41_hg38),   end(ACTL8_hg38),   end(SUSPECT_ALICE))),
  name       = c("LTR41_1", "LTR41_2", "ACTL8", "SUSPECT_1", "SUSPECT_2"),
  annot_type = c("TE",      "TE",      "gene",  "geneVMR",   "geneVMR"))

# ══════════════════════════════════════════════════════════════════════════════
# CpG extraction for pchuckle (Matt's regions only)
# ══════════════════════════════════════════════════════════════════════════════

matt_regions_to_extract <- c(
  region1_gr[is.na(region1_gr$name)],
  region2_gr[is.na(region2_gr$name)])

cpg_extract_list <- suppressWarnings(rbindlist(lapply(
  seq_along(matt_regions_to_extract), function(i) {
    gr   <- matt_regions_to_extract[i]
    hits <- findOverlaps(table3layers, gr)
    data.table(
      chr_pos = table3layers$chr_pos[queryHits(hits)],
      region  = if (!is.null(gr$name)) gr$name else paste0("region_", i))
  })))

unique_cpgs <- unique(cpg_extract_list$chr_pos)
unique_cpgs <- unique_cpgs[!is.na(unique_cpgs)]
message(sprintf("Total unique CpGs to extract: %d", length(unique_cpgs)))

out_path <- here("B_MultiTissues/dataOut",
                 sprintf("targetCpGs_Matt_%s.txt", format(Sys.time(), "%Y%m%d")))
writeLines(unique_cpgs, out_path)
message("Written to: ", out_path)

## git push, then in pchuckle:
# source /share/apps/source_files/python/python-3.13.0a6.source
# cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# python3 S00_extractRawMethylationForTargetCpG.py \
# --cpg_list  .../targetCpGs_Matt_YYYYMMDD.txt \
# --cpg_bed   .../CpG.bed.gz \
# --beta_files ".../GSM*.hg38.beta" \
# --meta      .../SupTab1_Loyfer2023_amended.csv \
# --output    .../gitignore/methylation_targetRegions.tsv \
# --minCov    10

# ══════════════════════════════════════════════════════════════════════════════
# Load extracted methylation
# ══════════════════════════════════════════════════════════════════════════════

meth <- fread(here("gitignore/methylation_targetRegions.tsv"))
setDT(meth)

loyfer_meta <- fread(here("B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv"))
loyfer_meta[, source_tissue_celltype := paste0(`Source Tissue`, " - ", `Cell type`)]
tissue_to_layer <- unique(loyfer_meta[, .(source_tissue_celltype, germ_layer = `Germ layer`)])
meth <- merge(meth, tissue_to_layer, by = "source_tissue_celltype", all.x = TRUE)
message(sprintf("Rows with missing germ_layer: %d", sum(is.na(meth$germ_layer))))

if (!"pos" %in% names(meth)) meth[, pos := as.integer(sub(".*_", "", cpg_site))]
if (!"chr" %in% names(meth)) meth[, chr := sub("_.*", "", cpg_site)]

germ_colours <- c(Endo = "#1D9E75", Meso = "#D85A30", Ecto = "#185FA5")

# ══════════════════════════════════════════════════════════════════════════════
# Inter-layer correlation (Loyfer, Endo×Meso only, n=4 patients, low power)
# ══════════════════════════════════════════════════════════════════════════════

compute_interlayer_r_loyfer <- function(meth_sub, min_patients = 3) {
  agg  <- meth_sub[, .(methylation = mean(methylation, na.rm = TRUE)),
                   by = .(cpg_site, patient_id, germ_layer)]
  wide <- dcast(agg, cpg_site + patient_id ~ germ_layer, value.var = "methylation")
  for (col in c("Endo","Meso","Ecto")) {
    if (!col %in% names(wide)) wide[, (col) := NA_real_]
  }
  wide[, {
    idx_EM <- !is.na(Endo) & !is.na(Meso)
    idx_EE <- !is.na(Endo) & !is.na(Ecto)
    idx_ME <- !is.na(Meso) & !is.na(Ecto)
    n_EM <- sum(idx_EM); n_EE <- sum(idx_EE); n_ME <- sum(idx_ME)
    r_EM <- if (n_EM >= min_patients) suppressWarnings(cor(Endo[idx_EM], Meso[idx_EM])) else NA_real_
    r_EE <- if (n_EE >= min_patients) suppressWarnings(cor(Endo[idx_EE], Ecto[idx_EE])) else NA_real_
    r_ME <- if (n_ME >= min_patients) suppressWarnings(cor(Meso[idx_ME], Ecto[idx_ME])) else NA_real_
    list(n_Endo_Meso = n_EM, n_Endo_Ecto = n_EE, n_Meso_Ecto = n_ME,
         r_Endo_Meso = r_EM, r_Endo_Ecto = r_EE, r_Meso_Ecto = r_ME,
         interlayer_r = mean(abs(c(r_EM, r_EE, r_ME)), na.rm = TRUE),
         low_power    = n_EM < 8)
  }, by = cpg_site]
}

interlayer_loyfer <- compute_interlayer_r_loyfer(meth, min_patients = 3)
message(sprintf("CpGs with Endo x Meso r: %d / %d",
                interlayer_loyfer[!is.na(r_Endo_Meso), .N], nrow(interlayer_loyfer)))

# ══════════════════════════════════════════════════════════════════════════════
# Plot functions (plot_raw_meth, plot_percpg_interlayer_corr, plot_region)
# are defined in functions.R — source it before running
# ══════════════════════════════════════════════════════════════════════════════

# ── Regional plots ────────────────────────────────────────────────────────────
p1 <- plot_region(region_gr = region1_gr,
                  annot_gr  = annot_gr_region1,
                  meth      = meth,
                  title     = "Region 1: LY6S-VMR / MER11C (chr8)")

p2 <- plot_region(region_gr = region2_gr,
                  annot_gr  = annot_gr_region2,
                  meth      = meth,
                  title     = "Region 2: ACTL8 / LTR41 (chr1)")

# ── Decay curves ──────────────────────────────────────────────────────────────
makePlotDecayTarget <- function(window_gr) {
  hits_target <- findOverlaps(window_gr, geomMeanGR)
  target_dt   <- data.table(
    alpha_geomean = geomMeanGR$alpha_geomean[subjectHits(hits_target)],
    ME            = paste0(window_gr$name, collapse = "_"))
  MEsetdt_with_target <- rbind(MEsetdt, target_dt)
  MEsetdt_with_target[, ME := relevel(factor(ME), ref = "mQTLcontrols")]
  plot_decay_curve(MEsetdt_with_target, title = "Decay curve including target region")
}

p4 <- makePlotDecayTarget(region1_gr[!is.na(region1_gr$name)])
p5 <- makePlotDecayTarget(region2_gr[grep("LTR41_2", region2_gr$name), ])

# ── LTR41 genome-wide comparison ─────────────────────────────────────────────
library(AnnotationHub)
ah       <- AnnotationHub()
rmskhg38 <- ah[["AH111333"]]

LTR41_all_hg38 <- rmskhg38[mcols(rmskhg38)$repName == "LTR41"]
message("Total LTR41 elements: ", length(LTR41_all_hg38))

target_idx <- which(
  as.character(seqnames(LTR41_all_hg38)) == "chr1" &
    start(LTR41_all_hg38) %in% start(LTR41_hg38))

target_cols      <- c("LTR41_1" = "#E69F00", "LTR41_2" = "#CC79A7")
mean_ltr41_width <- mean(width(LTR41_all_hg38))
ltr41_widths     <- width(LTR41_hg38)

fit_loess_ribbon <- function(dt, x_col, y_col, span = 0.5, n_out = 300,
                              max_n = 20000) {
  df <- dt[!is.na(get(y_col)), .(x = get(x_col), y = get(y_col))]
  if (nrow(df) < 10) return(NULL)
  if (nrow(df) > max_n) { set.seed(42); df <- df[sample(.N, max_n)] }
  fit   <- loess(y ~ x, data = df, span = span)
  x_seq <- seq(min(df$x), max(df$x), length.out = n_out)
  pred  <- predict(fit, newdata = data.frame(x = x_seq), se = TRUE)
  data.table(x = x_seq, y = pred$fit,
             ymin = pred$fit - 1.96 * pred$se.fit,
             ymax = pred$fit + 1.96 * pred$se.fit)
}

make_ltr41_plot <- function(smooth_mean, smooth_target, smooth_indiv,
                             y_lab, y_limits, y_breaks,
                             x_limits = NULL, extra_layers = NULL,
                             title = NULL, subtitle = NULL) {
  p <- ggplot() +
    geom_line(data = smooth_indiv,
              aes(x = x, y = y, group = ltr_idx),
              colour = "grey70", alpha = 0.2, linewidth = 0.2) +
    geom_ribbon(data = smooth_mean,
                aes(x = x, ymin = ymin, ymax = ymax),
                fill = "grey70", alpha = 0.3) +
    geom_line(data = smooth_mean, aes(x = x, y = y),
              colour = "grey30", linewidth = 1) +
    geom_ribbon(data = smooth_target,
                aes(x = x, ymin = ymin, ymax = ymax, fill = ltr_name),
                alpha = 0.25) +
    geom_line(data = smooth_target,
              aes(x = x, y = y, colour = ltr_name), linewidth = 1.2) +
    scale_colour_manual(values = target_cols, name = "Target LTR41") +
    scale_fill_manual(values   = target_cols, name = "Target LTR41") +
    geom_vline(xintercept = 0,                linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = mean_ltr41_width, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = ltr41_widths[1],
               colour = target_cols["LTR41_1"], linetype = "dotted", linewidth = 0.8) +
    geom_vline(xintercept = ltr41_widths[2],
               colour = target_cols["LTR41_2"], linetype = "dotted", linewidth = 0.8) +
    annotate("text", x = 0,                y = max(y_limits),
             label = "TE start", hjust = 0, vjust = 1.5, size = 2.5, colour = "grey40") +
    annotate("text", x = mean_ltr41_width, y = max(y_limits),
             label = sprintf("TE end\n(mean=%d bp)", round(mean_ltr41_width)),
             hjust = 1, vjust = 1.5, size = 2.5, colour = "grey40") +
    scale_x_continuous("Position relative to LTR41 start (bp)",
                       labels = function(x) paste0(round(x / 1e3, 1), " kb")) +
    scale_y_continuous(y_lab, limits = y_limits, breaks = y_breaks) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  if (!is.null(x_limits))     p <- p + coord_cartesian(xlim = x_limits)
  if (!is.null(extra_layers)) p <- p + extra_layers
  if (!is.null(title))        p <- p + ggtitle(title, subtitle = subtitle)
  p
}

LTR41_windows <- GRanges(
  seqnames = seqnames(LTR41_all_hg38),
  ranges   = IRanges(start = pmax(start(LTR41_all_hg38) - 3000, 1),
                     end   = end(LTR41_all_hg38) + 3000))

hits_all <- findOverlaps(table3layers, LTR41_windows)
cpg_dt   <- as.data.table(table3layers[queryHits(hits_all)])[
  , .(start, alpha_geomean, percentile, chr_pos)]
cpg_dt[, ltr_idx   := subjectHits(hits_all)]
cpg_dt[, ltr_start := start(LTR41_all_hg38)[ltr_idx]]
cpg_dt[, rel_pos   := start - ltr_start]
cpg_dt[, is_target := ltr_idx %in% target_idx]
cpg_dt[, ltr_name  := fcase(
  ltr_idx == target_idx[1], "LTR41_1",
  ltr_idx == target_idx[2], "LTR41_2",
  default = "other")]

ltr41_cpg_dt <- cpg_dt

set.seed(42)
non_target_idx <- setdiff(unique(ltr41_cpg_dt$ltr_idx), target_idx)
sample_idx     <- sample(non_target_idx, min(200, length(non_target_idx)))

message("Computing LTR41 smoothers...")
smooth_mean_alpha <- fit_loess_ribbon(ltr41_cpg_dt[is_target == FALSE], "rel_pos", "alpha_geomean")
smooth_mean_pct   <- fit_loess_ribbon(ltr41_cpg_dt[is_target == FALSE], "rel_pos", "percentile")

smooth_target_alpha <- rbindlist(lapply(unique(target_idx), function(ti) {
  dt <- ltr41_cpg_dt[ltr_idx == ti]
  s  <- fit_loess_ribbon(dt, "rel_pos", "alpha_geomean")
  if (is.null(s)) return(NULL)
  s[, ltr_name := dt$ltr_name[1]]; s}))

smooth_target_pct <- rbindlist(lapply(unique(target_idx), function(ti) {
  dt <- ltr41_cpg_dt[ltr_idx == ti]
  s  <- fit_loess_ribbon(dt, "rel_pos", "percentile")
  if (is.null(s)) return(NULL)
  s[, ltr_name := dt$ltr_name[1]]; s}))

smooth_indiv_alpha <- rbindlist(lapply(sample_idx, function(ti) {
  dt <- ltr41_cpg_dt[ltr_idx == ti & !is.na(alpha_geomean)]
  if (nrow(dt) < 5) return(NULL)
  fit   <- loess(alpha_geomean ~ rel_pos, data = dt, span = 0.5)
  x_seq <- seq(min(dt$rel_pos), max(dt$rel_pos), length.out = 100)
  data.table(x = x_seq, y = predict(fit, x_seq), ltr_idx = ti)}))

smooth_indiv_pct <- rbindlist(lapply(sample_idx, function(ti) {
  dt <- ltr41_cpg_dt[ltr_idx == ti & !is.na(percentile)]
  if (nrow(dt) < 5) return(NULL)
  fit   <- loess(percentile ~ rel_pos, data = dt, span = 0.5)
  x_seq <- seq(min(dt$rel_pos), max(dt$rel_pos), length.out = 100)
  data.table(x = x_seq, y = predict(fit, x_seq), ltr_idx = ti)}))

message("LTR41 smoothers done.")

zoom_levels <- list("1kb" = 1000, "2kb" = 2000, "3kb" = 3000)

ltr41_plots <- lapply(names(zoom_levels), function(zoom_name) {
  pad      <- zoom_levels[[zoom_name]]
  x_limits <- c(-pad, mean_ltr41_width + pad)
  p_alpha <- make_ltr41_plot(
    smooth_mean_alpha, smooth_target_alpha, smooth_indiv_alpha,
    y_lab = "Pr(HV) geomean", y_limits = c(0,1), y_breaks = c(0,.25,.5,.75,1),
    x_limits = x_limits,
    title    = sprintf("LTR41 +/- %s", zoom_name),
    subtitle = sprintf("grey = all LTR41 (n=%d, n=%d shown) | coloured = our 2 LTR41s",
                       length(LTR41_all_hg38), length(sample_idx)))
  p_pct <- make_ltr41_plot(
    smooth_mean_pct, smooth_target_pct, smooth_indiv_pct,
    y_lab = "Percentile (geomean)", y_limits = c(0,100),
    y_breaks = c(0,25,50,75,95,100), x_limits = x_limits,
    extra_layers = geom_hline(yintercept = 95,
                              linetype = "dashed", colour = "firebrick"))
  p_alpha / p_pct + plot_layout(guides = "collect") & theme(legend.position = "right")
})

legend_plot <- make_ltr41_plot(
  smooth_mean_alpha, smooth_target_alpha, smooth_indiv_alpha,
  y_lab = "Pr(HV) geomean", y_limits = c(0,1), y_breaks = c(0,.25,.5,.75,1),
  x_limits = c(-3000, mean_ltr41_width + 3000))
legend <- cowplot::get_legend(legend_plot)
ltr41_plots_nolegend <- lapply(ltr41_plots, function(x) x & theme(legend.position = "none"))

# ── LD heatmaps ───────────────────────────────────────────────────────────────
ld_region1_gr <- GRanges("chr8",
  IRanges(start = min(start(LY6SVMR_hg38), start(MER11C_hg38)) - 3000,
          end   = max(end(LY6SVMR_hg38),   end(MER11C_hg38))   + 3000))
ld_region2_gr <- GRanges("chr1",
  IRanges(start = min(start(LTR41_hg38)) - 3000,
          end   = max(end(LTR41_hg38))   + 3000))

feature_ticks_region1 <- data.table(
  pos   = c(start(LY6SVMR_hg38), end(LY6SVMR_hg38),
            start(MER11C_hg38),  end(MER11C_hg38)),
  label = c("VMR start","VMR end","MER11C start","MER11C end"))

feature_ticks_region2 <- data.table(
  pos   = c(start(SUSPECT_ALICE[1]), end(SUSPECT_ALICE[1]),
            start(SUSPECT_ALICE[2]), end(SUSPECT_ALICE[2]),
            start(LTR41_hg38[1]),   end(LTR41_hg38[1]),
            start(LTR41_hg38[2]),   end(LTR41_hg38[2])),
  label = c("SUSPECT_1 start","SUSPECT_1 end","SUSPECT_2 start","SUSPECT_2 end",
            "LTR41_1 start","LTR41_1 end","LTR41_2 start","LTR41_2 end"))

ld_region1 <- make_region_ld_plots(ld_region1_gr, meth, "LY6S-VMR / MER11C",
                                    feature_ticks_region1)
ld_region2 <- make_region_ld_plots(ld_region2_gr, meth, "ACTL8 / LTR41",
                                    feature_ticks_region2)

# ── Final figure assembly ─────────────────────────────────────────────────────
p6 <- cowplot::plot_grid(
  cowplot::plot_grid(plotlist = ltr41_plots_nolegend, ncol = 3,
                     labels = c("e.i","e.ii","e.iii")),
  legend, ncol = 2, rel_widths = c(1, 0.08))

p_ld1 <- ld_region1$meth
p_ld2 <- ld_region2$meth

row1 <- cowplot::plot_grid(p1, p2,       ncol = 2, labels = c("a","b"))
row2 <- cowplot::plot_grid(p_ld1, p_ld2, ncol = 2, labels = c("c","d"))
row3_labelled <- cowplot::plot_grid(
  cowplot::plot_grid(NULL, labels = "e"),
  p6, ncol = 1, rel_heights = c(0.05, 1))
row4 <- cowplot::plot_grid(p4, p5, ncol = 2, labels = c("f","g"))

final_plot <- cowplot::plot_grid(
  row1, row2, row3_labelled, row4,
  ncol        = 1,
  rel_heights = c(3, 1.8, 2, 1))

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/ACTL8_LTR41.pdf"),
  plot     = final_plot, width = 14, height = 22)

message("Saved: ACTL8_LTR41.pdf")
