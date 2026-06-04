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

load(here("gitignore/fullTable3layers.Rda"))
seqlevels(table3layers) <- paste0("chr", seqlevels(table3layers))

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

############################
## Define candidate sites ##
############################

## Matt's data are in hg19
LTR41table_gr_hg19 <- GRanges(seqnames = c("chr1", "chr1"), 
                              ranges = IRanges(start = c(18081648, 18085651),
                                               end = c(18082190, 18086109)),
                              name = c("LTR41_1", "LTR41_2"))

## ACTL8 position in hg38: chr1:17755333-17827063 (+)

# --- Liftover (hg19 → hg38) ---
LTR41table_gr_hg38 <- unlist(liftOver(LTR41table_gr_hg19, chain))

LTR41table_gr_hg38 <- c(
  LTR41table_gr_hg38,
  GRanges(seqnames = "chr1", 
          ranges = IRanges(start = 17755333, end = 17827063),
          name = "ACTL8"))

## Focus 10k around ACTL8
LTR41table_gr_hg38 <- c(
  LTR41table_gr_hg38,
  GRanges(seqnames = "chr1", 
          ranges = IRanges(start = min(LTR41table_gr_hg38@ranges@start) - 3000,
                           end = max(LTR41table_gr_hg38@ranges@start+
                                       LTR41table_gr_hg38@ranges@width) + 3000),
          name = NA))

# ── 1. Overlap table3layers with LTR41 region ─────────────────────────────────
hits_t3 <- findOverlaps(table3layers, LTR41table_gr_hg38[is.na(LTR41table_gr_hg38$name)])
t3_sub  <- as.data.table(table3layers[queryHits(hits_t3)])
# keep only metadata columns (drop GRanges coordinate columns if not needed)
t3_sub  <- t3_sub[, .(seqnames, start, chr_pos, 
                      alpha_endo, alpha_ecto, alpha_meso, 
                      alpha_allLayers, alpha_geomean, percentile)]

# ── 2. Overlap interlayer_corr with LTR41 region ──────────────────────────────
# interlayer_corr is a tibble with chr_pos → convert to GRanges first
interlayer_corr_clean <- interlayer_corr[!is.na(interlayer_corr$chr_pos), ]

interlayer_gr <- GRanges(
  seqnames = sub("_.*", "", interlayer_corr_clean$chr_pos),
  ranges   = IRanges(
    start = as.integer(sub(".*_", "", interlayer_corr_clean$chr_pos)),
    width = 1
  )
)
mcols(interlayer_gr) <- interlayer_corr_clean

hits_ic  <- findOverlaps(interlayer_gr, LTR41table_gr_hg38[is.na(LTR41table_gr_hg38$name)])
ic_sub   <- as.data.table(mcols(interlayer_gr[queryHits(hits_ic)]))

# ── Join on chr_pos ────────────────────────────────────────────────────────
result <- merge(t3_sub, ic_sub, by = "chr_pos", all = TRUE)
result <- result[!is.na(result$start),]

# ── result → GRanges ──────────────────────────────────────────────────────────
result_gr <- GRanges(
  seqnames = result$seqnames,
  ranges   = IRanges(start = result$start, width = 1)
)
mcols(result_gr)$chr_pos <- result$chr_pos

# ── overlap with putativeME_GR ────────────────────────────────────────────────
hits <- findOverlaps(result_gr, putativeME_GR)
length(hits) # no hits

# ── melt to long for the alpha / r tracks ─────────────────────────────────────
alpha_cols <- c("alpha_endo", "alpha_ecto", "alpha_meso", "alpha_geomean")
r_cols     <- c("r_Endo_Meso", "r_Endo_Ecto", "r_Meso_Ecto", "interlayer_r")

long_alpha <- melt(result, id.vars = c("start"),
                   measure.vars = alpha_cols,
                   variable.name = "track", value.name = "value")

long_r     <- melt(result, id.vars = c("start"),
                   measure.vars = r_cols,
                   variable.name = "track", value.name = "value")

# ── colour for set annotation ─────────────────────────────────────────────────
set_colours <- c("Derakhshan" = "#E69F00",
                 "Gunasekara" = "#56B4E9",
                 "other"      = "#009E73")

# ── shared theme ──────────────────────────────────────────────────────────────
th <- theme_bw(base_size = 10) +
  theme(axis.title.x  = element_blank(),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position  = "none")

th_bottom <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "none")

# ── LTR41 / ACTL8 annotation rectangles ──────────────────────────────────────
annot <- as.data.table(LTR41table_gr_hg38)[!is.na(name)]

add_annot <- function(p) {
  p + geom_rect(data = annot,
                aes(xmin = start, xmax = end,
                    ymin = -Inf, ymax = Inf,
                    fill = name),
                inherit.aes = FALSE,
                alpha = 0.12) +
    scale_fill_manual(values = c("LTR41_1" = "#CC79A7",
                                 "LTR41_2" = "#CC79A7",
                                 "ACTL8"   = "grey60"))
}

# ── 1. Alpha tracks (0–1) ─────────────────────────────────────────────────────
p_alpha <- ggplot(long_alpha, aes(x = start, y = value, colour = track)) +
  geom_smooth(alpha = 0.2, na.rm = TRUE, span = 0.1) +
  geom_point(size = 1.2, alpha = 0.7, na.rm = TRUE) +
  scale_colour_manual(values = c(
    alpha_endo      = "#1D9E75",
    alpha_ecto      = "#185FA5",
    alpha_meso      = "#D85A30",
    alpha_geomean   = "black")) +
  scale_y_continuous("Pr(HV)", limits = c(0, 1)) +
  th +
  theme(legend.position = "right")
p_alpha <- add_annot(p_alpha)

# ── 2. Percentile (alpha_geomean) ─────────────────────────────────────────────
p_pct <- ggplot(result, aes(x = start, y = percentile)) +
  geom_smooth(alpha = 0.2, na.rm = TRUE, span = 0.1) +
  geom_point(aes(colour = !is.na(set)), size = 1.5, na.rm = TRUE) +
  geom_hline(yintercept = 95, linetype = "dashed", colour = "firebrick") +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#E69F00")) +
  scale_y_continuous("Percentile\n(geomean)") +
  th
p_pct <- add_annot(p_pct)

# ── 3. Interlayer r tracks ────────────────────────────────────────────────────
p_r <- ggplot(long_r[long_r$track %in% "interlayer_r",],
              aes(x = start, y = value, colour = track)) +
  geom_point(size = 1.2, alpha = 0.7, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(
    r_Endo_Meso  = "#1D9E75",
    r_Endo_Ecto  = "#185FA5",
    r_Meso_Ecto  = "#D85A30",
    interlayer_r = "black")) +
  scale_y_continuous("Interlayer r", limits = c(-1, 1)) +
  th +
  theme(legend.position = "right")
p_r <- add_annot(p_r)

# ── 4. Percentile_r ───────────────────────────────────────────────────────────
p_pct_r <- ggplot(result, aes(x = start, y = percentile_r)) +
  geom_point(aes(colour = !is.na(set)), size = 1.5, na.rm = TRUE) +
  geom_hline(yintercept = 95, linetype = "dashed", colour = "firebrick") +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#E69F00")) +
  scale_y_continuous("Percentile\n(interlayer r)") +
  scale_x_continuous("Position (hg38)") +
  th_bottom
p_pct_r <- add_annot(p_pct_r)

# ── 5. Stack ──────────────────────────────────────────────────────────────────
(p_alpha / p_pct / p_r / p_pct_r) +
  plot_layout(heights = c(2, 1, 2, 1))

# ggplot2::ggsave(
#   filename = here::here(paste0("B_MultiTissues/dataOut/figures/SIV/allEPICfetus_corprhv_rinter.pdf")),
#   plot = p3, width = 5, height = 5
# )
