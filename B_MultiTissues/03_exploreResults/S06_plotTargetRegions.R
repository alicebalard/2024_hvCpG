#####################################################################
## Check hypervariability for target regions
#####################################################################

#####################################################################
## Prepare
library(here)
## Load libraries
source(here("B_MultiTissues", "quiet_library.R"))
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## Load functions
if (!exists("functionsLoaded")) {
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))}

## Add previous MEs including Maria's results
if (!exists("previousSIVprepared")) {
  source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))}
#####################################################################

## Object created in S04 (GRanges of CpGs covered in all 3 layers, with per-layer scores)
load(here("gitignore/table3layers_coveredIn3_26_08_26.Rda"))

## Table 3 layers, before subset for those covered in 3 layers analyses
load(here(paste0("gitignore/fulltable3layers_26_08_26.Rda")))
if (sum(grepl("chr", seqlevels(table3layers))) == 0){
  seqlevels(table3layers) <- paste0("chr", seqlevels(table3layers))
}

dt <- as.data.table(table3layers_coveredIn3)
dt_full <- as.data.table(table3layers)

dt$percentile_logBF_per_ds_allLayers <- ecdf(
  dt$logBF_per_ds_allLayers)(dt$logBF_per_ds_allLayers) * 100
dt_full$percentile_logBF_per_ds_allLayers <- ecdf(
  dt_full$logBF_per_ds_allLayers)(dt_full$logBF_per_ds_allLayers) * 100

## ── colour per feature TYPE (region backdrops) ───────────────────────────────
type_colours <- c(
  gene = "grey10",
  TE   = "green",
  VMR  = "#FFB90F",
  other = "lightgrey"
)

## ── colour per GENOMIC FEATURE (the CpG point fills) ─────────────────────────
feature_colours <- c(
  promoter   = "#D73027",   # red
  exon       = "#4575B4",   # blue
  intron     = "#91BFDB",   # light blue
  intergenic = "grey85"
)

## ── CpG GRanges from dt (once) ───────────────────────────────────────────────
if (!exists("dt_gr")) {
  dt_gr <- GRanges(sub("_.*", "", dt$chr_pos),
                   IRanges(as.integer(sub(".*_", "", dt$chr_pos)), width = 1),
                   percentile = dt$percentile_logBF_per_ds_allLayers,
                   chr_pos    = dt$chr_pos)
}
if (!exists("dt_gr_full")) {
  dt_gr_full <- GRanges(sub("_.*", "", dt_full$chr_pos),
                        IRanges(as.integer(sub(".*_", "", dt_full$chr_pos)), width = 1),
                        percentile = dt_full$percentile_logBF_per_ds_allLayers,
                        chr_pos    = dt_full$chr_pos)
}

## ── build the txdb feature sets ONCE (reused for point annotation) ────────────
if (!exists("txdb")) txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
if (!exists("exon_gr"))   exon_gr   <- reduce(exons(txdb))
if (!exists("intron_gr")) intron_gr <- reduce(unlist(intronsByTranscript(txdb)))
if (!exists("prom_gr"))   prom_gr   <- reduce(promoters(txdb, upstream = 2000, downstream = 200))

## classify each CpG (a 1bp GRanges) into promoter/exon/intron/intergenic.
## precedence: promoter > exon > intron > intergenic.
classify_feature <- function(gr) {
  f <- rep("intergenic", length(gr))
  f[overlapsAny(gr, intron_gr, ignore.strand = TRUE)] <- "intron"
  f[overlapsAny(gr, exon_gr,   ignore.strand = TRUE)] <- "exon"
  f[overlapsAny(gr, prom_gr,   ignore.strand = TRUE)] <- "promoter"
  factor(f, levels = c("promoter","exon","intron","intergenic"))
}

## ── optional: TEs in a window (unchanged) ────────────────────────────────────
get_TEs_in_window <- function(chr, x_min, x_max,
                              te_classes = c("LINE","SINE","LTR","DNA","RC","Retroposon"),
                              te_names = NULL, te_family = NULL) {
  if (!exists("rmskhg38", envir = .GlobalEnv)) {
    message("Loading hg38 RepeatMasker (AH111333) — once...")
    ah <- AnnotationHub::AnnotationHub()
    assign("rmskhg38", ah[["AH111333"]], envir = .GlobalEnv)
  }
  rmsk <- get("rmskhg38", envir = .GlobalEnv)
  win  <- GRanges(chr, IRanges(x_min, x_max))
  te <- rmsk[mcols(rmsk)$repClass %in% te_classes]
  if (!is.null(te_family)) te <- te[mcols(te)$repFamily %in% te_family]
  if (!is.null(te_names))  te <- te[mcols(te)$repName   %in% te_names]
  te <- subsetByOverlaps(te, win, ignore.strand = TRUE)
  if (!length(te)) return(NULL)
  GRanges(seqnames(te), ranges(te), name = mcols(te)$repName, type = "TE")
}

## ═════════════════════════════════════════════════════════════════════════════
## Plot: percentile along a region.
##   - feature backdrops (gene/TE/VMR) shaded by `type`
##   - CpG POINTS: pch 21, size 3, black outline, FILL = genomic feature
##     (promoter / exon / intron / intergenic)
##   - NO gene-model track below
## ═════════════════════════════════════════════════════════════════════════════
`%||%` <- function(a, b) if (is.null(a)) b else a

plot_region_percentile <- function(region_gr, dt_gr, flank = 5000,
                                   add_TEs = FALSE, te_names = NULL, te_family = NULL,
                                   xlim = NULL, title = NULL,
                                   type_colours    = get("type_colours",    envir = .GlobalEnv),
                                   feature_colours = get("feature_colours", envir = .GlobalEnv)) {
  chr   <- as.character(seqnames(region_gr))[1]
  if (!is.null(xlim)) { x_min <- xlim[1]; x_max <- xlim[2] }
  else { x_min <- min(start(region_gr)) - flank; x_max <- max(end(region_gr)) + flank }
  
  ## optionally add TEs to the shaded backdrops
  if (add_TEs) {
    tes <- get_TEs_in_window(chr, x_min, x_max, te_names = te_names, te_family = te_family)
    if (!is.null(tes)) {
      keep <- function(g) GRanges(seqnames(g), ranges(g), name = g$name, type = g$type)
      region_gr <- c(keep(region_gr), keep(tes))
    }
  }
  
  ## covered CpGs in the window, classified by genomic feature
  win <- GRanges(chr, IRanges(x_min, x_max))
  hit <- dt_gr[overlapsAny(dt_gr, win, ignore.strand = TRUE)]
  if (!length(hit)) { message("No covered CpGs in window for ", title); return(NULL) }
  cpg_dt <- data.table(pos = start(hit), percentile = hit$percentile,
                       feature = classify_feature(hit))
  
  ## feature backdrops
  feat <- as.data.table(region_gr)[, .(start, end, name = region_gr$name, type = region_gr$type)]
  feat <- feat[order(start)]
  feat[, `:=`(start = pmax(start, x_min), end = pmin(end, x_max))]
  feat <- feat[end > start]
  
  ggplot(cpg_dt, aes(pos, percentile)) +
    geom_rect(data = feat, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = 0, ymax = 100, fill = type),
              alpha = 0.15) +
    geom_rect(data = feat, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = 101, ymax = 106, fill = type),
              alpha = 0.4) +
    geom_text(data = feat, inherit.aes = FALSE,
              aes(x = (start + end) / 2, y = 103.5, label = name),
              size = 2.4, angle = 45, colour = "firebrick", fontface = "bold") +
    geom_line(col = "grey", alpha = .5) +
    ## CpG points: pch 21, size 3, black outline, fill = genomic feature
    geom_point(aes(fill = feature), shape = 21, colour = "black", size = 3, stroke = 0.4) +
    geom_hline(yintercept = c(95, 99), linetype = c("dashed", "dotted"), colour = "firebrick") +
    scale_x_continuous("Position (hg38)",
                       labels = function(x) paste0(round(x / 1e3, 1), " kb"),
                       limits = c(x_min, x_max), expand = c(0.01, 0)) +
    scale_y_continuous("Hypervariability score percentile\n(all layers)",
                       limits = c(0, 106), breaks = c(0, 25, 50, 75, 95, 100)) +
    ## ONE fill scale must serve both the backdrops (type) and points (feature):
    ## combine the two palettes so every level has a colour.
    scale_fill_manual(values = c(type_colours, feature_colours), name = NULL) +
    ggtitle(title %||% paste(unique(feat$name), collapse = " / ")) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

## ═════════════════════════════════════════════════════════════════════════════
## Region definitions
## ═════════════════════════════════════════════════════════════════════════════
LY6SVMR_hg19  <- GRanges("chr8", IRanges(144120106, 144120706), name = "LY6S-VMR")
LY6SVMR_hg38  <- unlist(liftOver(LY6SVMR_hg19, chain))
LY6SVMR_hg38$name <- "LY6S-VMR"; LY6SVMR_hg38$type <- "VMR"

LY6S_AS1_hg38 <- GRanges("chr8", IRanges(143039209, 143054303), name = "LY6S_AS1", type = "gene")
MER11C_hg38   <- GRanges("chr8", IRanges(143040739, 143041811), name = "MER11C",  type = "TE")

LTR41_hg19 <- GRanges("chr1",
                      IRanges(start = c(18081648, 18085651),
                              end   = c(18082190, 18086109)),
                      name = c("LTR41_1", "LTR41_2"))
LTR41_hg38 <- unlist(liftOver(LTR41_hg19, chain))
LTR41_hg38$name <- c("LTR41_1", "LTR41_2"); LTR41_hg38$type <- "TE"

ACTL8_hg38 <- GRanges("chr1", IRanges(17755333, 17827063), name = "ACTL8", type = "gene")

## ═════════════════════════════════════════════════════════════════════════════
## Plot  (now single-panel: plot_region_percentile directly, no gene model)
## ═════════════════════════════════════════════════════════════════════════════
region1 <- c(LY6SVMR_hg38, LY6S_AS1_hg38, MER11C_hg38)

p_region1 <- plot_region_percentile(
  region1, dt_gr, flank = 3000, add_TEs = FALSE,
  title = "chr8: LY6S-VMR / LY6S-AS1 / MER11C")

p_region1_zoom <- plot_region_percentile(
  region1, dt_gr, flank = 1000, add_TEs = FALSE,
  title = "chr8: LY6S-VMR / LY6S-AS1 / MER11C  (zoomed in)",
  xlim = c(143038689, 143045000))

region2 <- c(LTR41_hg38, ACTL8_hg38)

p_region2 <- plot_region_percentile(
  region2, dt_gr, flank = 3000, add_TEs = FALSE,
  title = "chr1: LTR41 / ACTL8")

p_region2_zoom <- plot_region_percentile(
  region2, dt_gr, flank = 1000, add_TEs = FALSE,
  title = "chr1: LTR41 / ACTL8 (zoomed in)",
  xlim = c(17755153, 17762000))

ggsave(here("B_MultiTissues/dataOut/figures/script06/hypervariabilityTargetsMatt.png"),
       plot_grid(p_region1, p_region1_zoom, p_region2, p_region2_zoom, nrow = 2),
       width = 18, height = 8, dpi = 300, bg = "white")

############ Candidate MEs
plot_POMC <- plot_region_percentile(
  GRanges("chr2", IRanges(25160860, 25168580), name = "POMC", type = "gene"),
  dt_gr, flank = 1000, add_TEs = TRUE, title = "chr2: POMC")

plot_CYP2E1 <- plot_region_percentile(
  GRanges("chr10", IRanges(133527363, 133539123), name = "CYP2E1", type = "gene"),
  dt_gr, flank = 1000, add_TEs = TRUE, title = "chr10: CYP2E1")

plot_PAX8 <- plot_region_percentile(
  GRanges("chr2", IRanges(113215997, 113278921), name = "PAX8", type = "gene"),
  dt_gr, flank = 1000, add_TEs = TRUE, title = "chr2: PAX8")

ggsave(here("B_MultiTissues/dataOut/figures/script06/hypervariabilityMEs.png"),
       plot_grid(plot_POMC, plot_CYP2E1, plot_PAX8, nrow = 3),
       width = 12, height = 9, dpi = 300, bg = "white")
