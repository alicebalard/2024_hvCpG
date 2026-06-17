#####################################################################
## Check pr(hv) geom means for target regions vs random background ##
#####################################################################

#####################################################################
## Prepare
library(here)
## Load libraries
source(here("B_MultiTissues", "quiet_library.R"))

## Load functions
if (!exists("functionsLoaded")) {
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))}

## Add previous MEs including Maria's results
## Load the set of previously tested MEs & vmeQTL
if (!exists("previousSIVprepared")) {
  source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))}
#####################################################################

# objects needed from S06 — load only what's required
if (!exists("putativeME_GR")) source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))
if (!exists("table3layers")) load(here("gitignore/fullTable3layers.Rda"))
if (!exists("MEsetdt")) MEsetdt <- readRDS(here("gitignore/MEsetdt.rds"))
if (!exists("MEsetdt_regionMean")) MEsetdt_regionMean <- readRDS(here("gitignore/MEsetdt_regionMean.rds"))
if (!exists("geomMeanGR")) geomMeanGR <- readRDS(here("gitignore/geomMeanGR.rds"))
################################################################################

# Compute the percentile rank of each alpha_geomean value
# (what % of all values are <= this value)
table3layers$percentile <- ecdf(table3layers$alpha_geomean)(table3layers$alpha_geomean) * 100
# Percentile = 95 means the site is in the top 5%
# top X% = percentile >= (100 - X)

## Interlayer correlation obtained from fetal script
interlayer_corr <- readRDS(here("B_MultiTissues/dataOut/interlayer_corr_all.RDS"))
interlayer_corr$chr_pos <- dico$chrpos_hg38[match(interlayer_corr$CpG, dico$CpG)]

## to check!!! commit
interlayer_corr$interlayer_r <- rowMeans(
  abs(interlayer_corr[, c("r_Endo_Meso", "r_Endo_Ecto", "r_Meso_Ecto")]),
  na.rm = TRUE
)
interlayer_corr$percentile_r <- ecdf(interlayer_corr$interlayer_r)(interlayer_corr$interlayer_r) * 100

# ══════════════════════════════════════════════════════════════════════════════
# Generic plotting function for any genomic region
# ══════════════════════════════════════════════════════════════════════════════

plot_region <- function(region_gr,    # GRanges with named elements + NA for the window
                        annot_colours, # named vector of colours for named elements
                        title = "") {
  
  # ── Extract the query window (the NA-named range) ──────────────────────────
  window_gr <- region_gr[is.na(region_gr$name)]
  annot_dt  <- as.data.table(region_gr[!is.na(region_gr$name)])
  annot_dt[, mid := (start + end) / 2]  # midpoint for label
  
  
  # ── Overlap table3layers with window ──────────────────────────────────────
  hits_t3 <- findOverlaps(table3layers, window_gr)
  t3_sub  <- as.data.table(table3layers[queryHits(hits_t3)])[
    , .(seqnames, start, chr_pos,
        alpha_endo, alpha_ecto, alpha_meso,
        alpha_allLayers, alpha_geomean, percentile)]
  
  # ── Overlap interlayer_corr with window ───────────────────────────────────
  hits_ic <- findOverlaps(interlayer_gr, window_gr)
  ic_sub  <- as.data.table(mcols(interlayer_gr[queryHits(hits_ic)]))
  
  # ── Join ──────────────────────────────────────────────────────────────────
  result <- merge(t3_sub, ic_sub, by = "chr_pos", all = TRUE)
  result <- result[!is.na(result$start)]
  
  # abs r values
  r_cols <- c("r_Endo_Meso", "r_Endo_Ecto", "r_Meso_Ecto", "interlayer_r")
  result[, (r_cols) := lapply(.SD, abs), .SDcols = r_cols]
  
  # ── Melt ──────────────────────────────────────────────────────────────────
  alpha_cols <- c("alpha_endo", "alpha_ecto", "alpha_meso", "alpha_geomean")
  long_alpha <- melt(result, id.vars = "start",
                     measure.vars = alpha_cols,
                     variable.name = "track", value.name = "value")
  long_r     <- melt(result, id.vars = "start",
                     measure.vars = r_cols,
                     variable.name = "track", value.name = "value")
  
  x_min <- start(window_gr)
  x_max <- end(window_gr)
  
  # ── Annotation helper ─────────────────────────────────────────────────────
  add_annot <- function(p, show_labels = FALSE) {
    p <- p +
      geom_rect(data = annot_dt,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = name),
                inherit.aes = FALSE, alpha = 0.12) +
      scale_fill_manual(values = annot_colours)
    
    if (show_labels) {
      p <- p +
        geom_text(data = annot_dt,
                  aes(x = mid, y = Inf, label = name),
                  inherit.aes = FALSE,
                  vjust = 1.4, size = 3, colour = "black", fontface = "bold")
    }
    p
  }
  
  # ── Themes ────────────────────────────────────────────────────────────────
  th <- theme_bw(base_size = 10) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position  = "right")
  
  th_bottom <- theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "none")
  
  # ── P1: Alpha tracks — labels here only ───────────────────────────────────
  p_alpha <- add_annot(
    ggplot(long_alpha, aes(x = start, y = value, colour = track)) +
      geom_smooth(alpha = 0.1, na.rm = TRUE, span = 0.3) +
      geom_point(size = 1.2, alpha = 0.4, na.rm = TRUE) +
      scale_colour_manual(values = c(
        alpha_endo    = "#1D9E75",
        alpha_ecto    = "#185FA5",
        alpha_meso    = "#D85A30",
        alpha_geomean = "black")) +
      scale_y_continuous("Pr(HV)", limits = c(0, 1)) +
      coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
      ggtitle(title) +
      th,
    show_labels = TRUE   # <-- labels only on top panel
  )
  
  # ── P2: Percentile geomean ────────────────────────────────────────────────
  p_pct <- add_annot(
    ggplot(result, aes(x = start, y = percentile)) +
      geom_smooth(alpha = 0.1, na.rm = TRUE, span = 0.3) +
      geom_point(size = 1.5, alpha = 0.4, na.rm = TRUE, colour = "grey50") +
      geom_hline(yintercept = 95, linetype = "dashed", colour = "firebrick") +
      scale_y_continuous("Percentile\n(geomean)") +
      coord_cartesian(xlim = c(x_min, x_max)) +
      th_bottom
  )
  
  # ── P3: Interlayer r ──────────────────────────────────────────────────────
  p_r <- add_annot(
    ggplot(long_r, aes(x = start, y = value, colour = track)) +
      geom_point(size = 1.2, alpha = 0.7, na.rm = TRUE) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      scale_colour_manual(values = c(
        r_Endo_Meso  = "#1D9E75",
        r_Endo_Ecto  = "#185FA5",
        r_Meso_Ecto  = "#D85A30",
        interlayer_r = "black")) +
      scale_y_continuous("|Interlayer r|", limits = c(0, 1)) +
      coord_cartesian(xlim = c(x_min, x_max)) +
      th
  )
  
  # ── P4: Percentile r ──────────────────────────────────────────────────────
  p_pct_r <- add_annot(
    ggplot(result, aes(x = start, y = percentile_r)) +
      geom_point(size = 1.5, alpha = 0.4, na.rm = TRUE, colour = "grey50") +
      geom_hline(yintercept = 95, linetype = "dashed", colour = "firebrick") +
      scale_y_continuous("Percentile\n(interlayer r)") +
      scale_x_continuous("Position (hg38)") +
      coord_cartesian(xlim = c(x_min, x_max)) +
      th_bottom
  )
  
  # ── P2-P4: no labels ──────────────────────────────────────────────────────
  p_pct   <- add_annot(p_pct,   show_labels = FALSE)
  p_r     <- add_annot(p_r,     show_labels = FALSE)
  p_pct_r <- add_annot(p_pct_r, show_labels = FALSE)
  
  # ── Stack vertically ──────────────────────────────────────────────────────
  p_alpha / p_pct / p_r / p_pct_r +
    plot_layout(heights = c(1, 1, 1, 1))
}

# ══════════════════════════════════════════════════════════════════════════════
# Prepare interlayer_gr once (shared across regions)
# ══════════════════════════════════════════════════════════════════════════════
interlayer_corr_clean <- interlayer_corr[!is.na(interlayer_corr$chr_pos), ]
interlayer_gr <- GRanges(
  seqnames = sub("_.*",  "", interlayer_corr_clean$chr_pos),
  ranges   = IRanges(
    start = as.integer(sub(".*_", "", interlayer_corr_clean$chr_pos)),
    width = 1)
)
mcols(interlayer_gr) <- interlayer_corr_clean

# ══════════════════════════════════════════════════════════════════════════════
# Region 1: LY6S-VMR / MER11C  (hg38 coordinates)
# LY6S-VMR hg19 chr8:144120106-144120706 → liftover to hg38
# MER11C   hg38 chr8:143040739-143041811
# Window: -3000 before LY6S-VMR start, +3000 after MER11C end
# ══════════════════════════════════════════════════════════════════════════════
LY6SVMR_hg19 <- GRanges("chr8", IRanges(144120106, 144120706), name = "LY6S-VMR")
LY6SVMR_hg38 <- unlist(liftOver(LY6SVMR_hg19, chain))
LY6SVMR_hg38$name <- "LY6S-VMR"

MER11C_hg38 <- GRanges("chr8", IRanges(143040739, 143041811), name = "MER11C")

region1_gr <- c(
  LY6SVMR_hg38,
  MER11C_hg38,
  GRanges("chr8",
          IRanges(
            start = min(start(LY6SVMR_hg38), start(MER11C_hg38)) - 3000,
            end   = max(end(LY6SVMR_hg38),   end(MER11C_hg38))   + 3000),
          name = NA)
)

p1 <- plot_region(
  region_gr     = region1_gr,
  annot_colours = c("LY6S-VMR" = "#CC79A7", "MER11C" = "#56B4E9"),
  title         = "Region 1: LY6S-VMR / MER11C (chr8)"
)

# ══════════════════════════════════════════════════════════════════════════════
# Region 2: ACTL8 / LTR41 
# Window: -3000 before first LTR41 start, +3000 after ACTL8 start
# ══════════════════════════════════════════════════════════════════════════════
LTR41_hg19 <- GRanges("chr1",
                      IRanges(start = c(18081648, 18085651),
                              end   = c(18082190, 18086109)),
                      name = c("LTR41_1", "LTR41_2"))
LTR41_hg38 <- unlist(liftOver(LTR41_hg19, chain))
LTR41_hg38$name <- c("LTR41_1", "LTR41_2")

ACTL8_hg38 <- GRanges("chr1", IRanges(17755333, 17827063), name = "ACTL8")

region2_gr <- c(
  LTR41_hg38,
  ACTL8_hg38,
  GRanges("chr1",
          IRanges(
            start = min(start(LTR41_hg38)) - 3000,
            end   = max(start(LTR41_hg38)) + 3000),
          name = NA)
)

p2 <- plot_region(
  region_gr     = region2_gr,
  annot_colours = c("LTR41_1" = "#56B4E9",
                    "LTR41_2" = "#56B4E9",
                    "ACTL8"   = "#CC79A7"),
  title         = "Region 2: ACTL8 / LTR41 (chr1)"
)

# ══════════════════════════════════════════════════════════════════════════════
# Region 3: ACTL8 / LTR41 
# Window: -10000 before first LTR41 start, +10000 after ACTL8 start
# ══════════════════════════════════════════════════════════════════════════════
region3_gr <- c(
  LTR41_hg38,
  ACTL8_hg38,
  GRanges("chr1",
          IRanges(
            start = min(start(LTR41_hg38)) - 10000,
            end   = max(start(LTR41_hg38)) + 10000),
          name = NA)
)

p3 <- plot_region(
  region_gr     = region3_gr,
  annot_colours = c("LTR41_1" = "#56B4E9",
                    "LTR41_2" = "#56B4E9",
                    "ACTL8"   = "#CC79A7"),
  title         = "Region 2: ACTL8 / LTR41 (chr1) extended"
)

####################################
## Compare with other sets of MEs ##
####################################

makePlotDecayTarget <- function(window_gr){
  hits_target <- findOverlaps(
    window_gr,
    geomMeanGR
  )
  
  target_dt <- data.table(
    alpha_geomean = geomMeanGR$alpha_geomean[subjectHits(hits_target)],
    ME = paste0(window_gr$name, collapse='_')
  )
  
  # ── Append to MEsetdt and replot ──────────────────────────────────────────────
  MEsetdt_with_target <- rbind(MEsetdt, target_dt)
  MEsetdt_with_target[, ME := relevel(factor(ME), ref = "mQTLcontrols")]
  
  plot_decay_curve(MEsetdt_with_target, 
                   title = "Decay curve including target region")
}

p4 <- makePlotDecayTarget(region1_gr[!is.na(region1_gr$name)])
p5 <- makePlotDecayTarget(region2_gr[grep("LTR41_2", region2_gr$name),])

ggplot2::ggsave(
  filename = here::here(paste0("B_MultiTissues/dataOut/figures/ACTL8_LTR41.pdf")),
  plot = cowplot::plot_grid(
    cowplot::plot_grid(p4, p5, ncol = 2, labels = c("a", "b")),
    cowplot::plot_grid(p1,p2, p3, ncol = 3, labels = c("c", "d", "e")),
    nrow = 2., rel_heights = c(1,3)),
  width = 15, height = 12
)

