#####################################################################
## Check hypervariability for target regions
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

## Object created in S04 (GRanges of CpGs covered in all 3 layers, with per-layer scores)
load(here("gitignore/table3layers_coveredIn3_26_08_26.Rda"))

## Table 3 layers, before subset for those covered in 3 layers analyses
load(here(paste0("gitignore/fulltable3layers_26_08_26.Rda")))
# Fix chromosome names in geomMeanGR (1 -> chr1)
if (sum(grepl("chr", seqlevels(table3layers))) == 0){
  seqlevels(table3layers) <- paste0("chr", seqlevels(table3layers))
}

dt <- as.data.table(table3layers_coveredIn3)           # one row per covered-in-3 CpG
dt_full <- as.data.table(table3layers)           # one row per covered-in-3 CpG

## Use percentile for logBF_per_ds, more interpretable to look at candidates
dt$percentile_logBF_per_ds_allLayers <- ecdf(
  dt$logBF_per_ds_allLayers)(dt$logBF_per_ds_allLayers) * 100

## Use percentile for logBF_per_ds, more interpretable to look at candidates
dt_full$percentile_logBF_per_ds_allLayers <- ecdf(
  dt_full$logBF_per_ds_allLayers)(dt_full$logBF_per_ds_allLayers) * 100

## ── colour per TYPE ──────────────────────────────────────────────────────────
type_colours <- c(
  gene = "grey10",   
  TE   = "#C1FFC1", 
  VMR  = "#FFB90F",
  other = "lightgrey"
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

## ── (3) optional: all TEs in a window, from hg38 RepeatMasker (same as S04) ────
## Returns a GRanges of TEs overlapping the window, with name = repName, type="TE".
## Loads rmskhg38 once and caches it in the global env.
get_TEs_in_window <- function(chr, x_min, x_max,
                              te_classes = c("LINE","SINE","LTR","DNA","RC","Retroposon"),
                              te_names   = NULL,      # e.g. "LTR41", or c("LTR41","MER11C")
                              te_family  = NULL) {    # e.g. "ERVL", or a repFamily
  if (!exists("rmskhg38", envir = .GlobalEnv)) {
    message("Loading hg38 RepeatMasker (AH111333) — once...")
    ah <- AnnotationHub::AnnotationHub()
    assign("rmskhg38", ah[["AH111333"]], envir = .GlobalEnv)
  }
  rmsk <- get("rmskhg38", envir = .GlobalEnv)
  win  <- GRanges(chr, IRanges(x_min, x_max))
  
  te <- rmsk[mcols(rmsk)$repClass %in% te_classes]
  if (!is.null(te_family)) te <- te[mcols(te)$repFamily %in% te_family]   # filter by family
  if (!is.null(te_names))  te <- te[mcols(te)$repName   %in% te_names]     # filter by exact name
  te <- subsetByOverlaps(te, win, ignore.strand = TRUE)
  if (!length(te)) return(NULL)
  GRanges(seqnames(te), ranges(te),
          name = mcols(te)$repName, type = "TE")
}

## ═════════════════════════════════════════════════════════════════════════════
## Plot: percentile along a region, features coloured by TYPE
##   region_gr : GRanges with mcols `name` and `type`
##   add_TEs   : if TRUE, auto-add all RepeatMasker TEs in the window
## ═════════════════════════════════════════════════════════════════════════════
`%||%` <- function(a, b) if (is.null(a)) b else a

plot_region_percentile <- function(region_gr, dt_gr, flank = 5000,
                                   add_TEs = FALSE, te_names = NULL, te_family = NULL,
                                   xlim = NULL,                       # <- optional c(min, max) to zoom
                                   title = NULL,
                                   type_colours = get("type_colours", envir = .GlobalEnv)) {
  chr   <- as.character(seqnames(region_gr))[1]
  # window: explicit xlim if given, else region +/- flank
  if (!is.null(xlim)) {
    x_min <- xlim[1]; x_max <- xlim[2]
  } else {
    x_min <- min(start(region_gr)) - flank
    x_max <- max(end(region_gr))   + flank
  }
  
  ## optionally pull in all TEs in the window and append to the features
  if (add_TEs) {
    tes <- get_TEs_in_window(chr, x_min, x_max,
                             te_names = te_names, te_family = te_family)   # <- pass filter
    if (!is.null(tes)) {
      keep <- function(g) GRanges(seqnames(g), ranges(g), name = g$name, type = g$type)
      region_gr <- c(keep(region_gr), keep(tes))
    }
  }
  
  ## covered CpGs in the window
  win <- GRanges(chr, IRanges(x_min, x_max))
  hit <- dt_gr[overlapsAny(dt_gr, win, ignore.strand = TRUE)]
  if (!length(hit)) { message("No covered CpGs in window for ", title); return(NULL) }
  cpg_dt <- data.table(pos = start(hit), percentile = hit$percentile)
  
  ## feature rectangles, with type (and a stacked y-row per feature so overlapping
  ## TEs/genes don't hide each other in the annotation strip)
  feat <- as.data.table(region_gr)[, .(start, end,
                                       name = region_gr$name,
                                       type = region_gr$type)]
  feat <- feat[order(start)]
  feat[, `:=`(start = pmax(start, x_min), end = pmin(end, x_max))]   # <- clip to window
  feat <- feat[end > start]                                          # drop features fully outside
  feat[, row := seq_len(.N)]                    # each feature its own shaded band height
  
  ggplot(cpg_dt, aes(pos, percentile)) +
    ## shaded features, coloured by type, spanning full y so they read as backdrops
    geom_rect(data = feat, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = 0, ymax = 100, fill = type),
              alpha = 0.8) +
    ## a thin labelled bar per feature at the top, to name each without clutter
    geom_rect(data = feat, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = 101, ymax = 106, fill = type),
              alpha = 0.4) +
    geom_text(data = feat, inherit.aes = FALSE,
              aes(x = (start + end) / 2, y = 103.5, label = name),
              size = 2.4, colour = "firebrick", fontface = "bold") +
    geom_line(col = "grey", alpha = .5) +
    geom_point(aes(col = percentile), size = 2) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 50)+
    guides(col="none") + # no legend
    geom_hline(yintercept = c(95, 99), linetype = c("dashed", "dotted"),
               colour = "firebrick") +
    scale_x_continuous("Position (hg38)",
                       labels = function(x) paste0(round(x / 1e3, 1), " kb"),
                       limits = c(x_min, x_max), expand = c(0.01, 0)) +
    scale_y_continuous("Hypervariability score percentile\n(all layers)",
                       limits = c(0, 106), breaks = c(0, 25, 50, 75, 95, 100)) +
    scale_fill_manual(values = type_colours, name = "Feature type") +
    ggtitle(title %||% paste(unique(feat$name), collapse = " / ")) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

## ── gene-model track for a window: promoter / exon / intron as a black bar ─────
gene_model_track <- function(chr, x_min, x_max,
                             txdb = TxDb.Hsapiens.UCSC.hg38.knownGene) {
  win <- GRanges(chr, IRanges(x_min, x_max))
  
  # exons and introns overlapping the window
  ex  <- subsetByOverlaps(reduce(exons(txdb)), win, ignore.strand = TRUE)
  intr <- subsetByOverlaps(reduce(unlist(intronsByTranscript(txdb))), win, ignore.strand = TRUE)
  # promoters (TSS -2000/+200) overlapping the window
  prom <- subsetByOverlaps(reduce(promoters(txdb, upstream = 2000, downstream = 200)),
                           win, ignore.strand = TRUE)
  
  # assemble into one data.table of features with a 'part' label
  parts <- rbindlist(list(
    if (length(intr)) data.table(start = start(intr), end = end(intr), part = "intron"),
    if (length(ex))   data.table(start = start(ex),   end = end(ex),   part = "exon"),
    if (length(prom)) data.table(start = start(prom), end = end(prom), part = "promoter")
  ), fill = TRUE)
  if (!nrow(parts)) return(NULL)
  # clip to window
  parts[, `:=`(start = pmax(start, x_min), end = pmin(end, x_max))]
  
  ggplot(parts) +
    # introns as a thin blocks, exons as thick blocks, promoters as a distinct block
    geom_rect(data = parts[part == "intron"],
              aes(xmin = start, xmax = end, ymin = 0.9, ymax = 1.1),
              fill = "grey") +
    geom_rect(data = parts[part == "exon"],
              aes(xmin = start, xmax = end, ymin = 0.7, ymax = 1.3),
              fill = "black") +
    geom_rect(data = parts[part == "promoter"],
              aes(xmin = start, xmax = end, ymin = 0.55, ymax = 1.45),
              fill = "firebrick") +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0.01, 0)) +
    scale_y_continuous(NULL, limits = c(0.4, 1.6), breaks = NULL) +
    labs(x = NULL) +
    scale_x_continuous("Position (hg38)",
                       labels = function(x) paste0(round(x / 1e3, 1), " kb"),
                       limits = c(x_min, x_max), expand = c(0.01, 0)) +
    scale_y_continuous(NULL, limits = c(0.4, 1.6), breaks = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin  = margin(0, 5, 0, 5))
}

## build both panels for a region and stack with shared x
plot_region_with_genemodel <- function(region_gr, dt_gr, flank = 3000,
                                       add_TEs = FALSE, te_names = NULL, te_family = NULL,
                                       xlim = NULL, title = NULL) {
  chr   <- as.character(seqnames(region_gr))[1]
  if (!is.null(xlim)) { x_min <- xlim[1]; x_max <- xlim[2] }
  else { x_min <- min(start(region_gr)) - flank; x_max <- max(end(region_gr)) + flank }
  
  p_main <- plot_region_percentile(region_gr, dt_gr, flank = flank,
                                   add_TEs = add_TEs, te_names = te_names, te_family = te_family,
                                   xlim = xlim, title = title) +           # <- pass xlim
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  p_track <- gene_model_track(chr, x_min, x_max)                            # uses the same bounds
  if (is.null(p_track)) return(p_main)
  p_main / p_track + plot_layout(heights = c(6, 1))
}

## ═════════════════════════════════════════════════════════════════════════════
## Region definitions — each subregion carries a `type` attribute
## (used to colour the fill and to group features in the legend).
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
## Plot
## ═════════════════════════════════════════════════════════════════════════════

region1 <- c(LY6SVMR_hg38, LY6S_AS1_hg38, MER11C_hg38)

p_region1 <- plot_region_with_genemodel(
  region1, dt_gr, flank = 3000,
  add_TEs = FALSE, te_names = "LTR41",
  title = "chr1: LY6S-VMR / LY6S-AS1 / MER11C")

p_region1_zoom <- plot_region_with_genemodel(
  region1, dt_gr, flank = 1000,
  add_TEs = FALSE, te_names = "LTR41",
  title = "chr1: LY6S-VMR / LY6S-AS1 / MER11C  (zoomed in)",
  xlim = c(143038689, 143045000))

region2 <- c(LTR41_hg38, ACTL8_hg38)

p_region2 <- plot_region_with_genemodel(
  region2, dt_gr, flank = 3000,
  add_TEs = FALSE, te_names = "LTR41",
  title = "chr1: LTR41 / ACTL8")

p_region2_zoom <- plot_region_with_genemodel(
  region2, dt_gr, flank = 1000,
  add_TEs = FALSE, te_names = "LTR41",
  title = "chr1: LTR41 / ACTL8 (zoomed in)",
  xlim = c(17755153, 17762000))

ggsave(here("B_MultiTissues/dataOut/figures/script06/hypervariabilityTargetsMatt.png"),
       plot_grid(p_region1, p_region1_zoom, p_region2, p_region2_zoom, nrow = 2),
       width = 18, height = 10, dpi = 300, bg = "white")

############ Candidate MEs 

# POMC
# NC_000002.12 (25160860..25168580, complement)
plot_POMC <- plot_region_with_genemodel(
  GRanges("chr2",
          IRanges(start = c(25160860),
                  end   = c(25168580)),
          name = "POMC", type = "gene"), dt_gr, flank = 1000,
  add_TEs = TRUE, title = "chr2: POMC")

plot_CYP2E1 <- plot_region_with_genemodel(
  GRanges("chr10",
          IRanges(start = 133527363,
                  end   = 133539123),
          name = "CYP2E1", type = "gene"), dt_gr, flank = 1000,
  add_TEs = TRUE, title = "chr10: CYP2E1")

ggsave(here("B_MultiTissues/dataOut/figures/script06/hypervariabilityMEs.png"),
       plot_grid(plot_POMC, plot_CYP2E1, nrow = 2),
       width = 10, height = 9, dpi = 300, bg = "white")
