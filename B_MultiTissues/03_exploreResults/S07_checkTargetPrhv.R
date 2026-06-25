#####################################################################
## Check pr(hv) geom means for target regions vs random background ##
#####################################################################

library(here)
source(here("B_MultiTissues", "quiet_library.R"))
if (!exists("functionsLoaded"))    source(here("B_MultiTissues/03_exploreResults", "functions.R"))
if (!exists("previousSIVprepared")) source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))
if (!exists("putativeME_GR"))       source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))
if (!exists("table3layers"))        load(here("gitignore/fullTable3layers.Rda"))
if (!exists("MEsetdt"))             MEsetdt            <- readRDS(here("gitignore/MEsetdt.rds"))
if (!exists("MEsetdt_regionMean"))  MEsetdt_regionMean <- readRDS(here("gitignore/MEsetdt_regionMean.rds"))
if (!exists("geomMeanGR"))          geomMeanGR         <- readRDS(here("gitignore/geomMeanGR.rds"))

if (is.null(table3layers$percentile)) {
  table3layers$percentile <- ecdf(table3layers$alpha_geomean)(table3layers$alpha_geomean) * 100
}

interlayer_corr <- readRDS(here("B_MultiTissues/dataOut/interlayer_corr_all.RDS"))
interlayer_corr$chr_pos <- dico$chrpos_hg38[match(interlayer_corr$CpG, dico$CpG)]
interlayer_corr$interlayer_r <- rowMeans(
  abs(interlayer_corr[, c("r_Endo_Meso", "r_Endo_Ecto", "r_Meso_Ecto")]),
  na.rm = TRUE)
interlayer_corr$percentile_r <- ecdf(interlayer_corr$interlayer_r)(
  interlayer_corr$interlayer_r) * 100

# ── Build interlayer_gr once (used inside plot_region) ───────────────────────
interlayer_corr_clean <- interlayer_corr[!is.na(interlayer_corr$chr_pos), ]

# check for non-standard chr_pos format
pos_check <- as.integer(sub(".*_", "", interlayer_corr_clean$chr_pos))
interlayer_corr_clean <- interlayer_corr_clean[!is.na(pos_check), ]

interlayer_gr <- GRanges(
  seqnames = sub("_.*", "", interlayer_corr_clean$chr_pos),
  ranges   = IRanges(
    start = as.integer(sub(".*_", "", interlayer_corr_clean$chr_pos)),
    width = 1
  )
)
mcols(interlayer_gr) <- interlayer_corr_clean

# ══════════════════════════════════════════════════════════════════════════════
# plot_region
# peak_anchor_gr : GRanges used ONLY to define the detection zone
#                  (should cover the TE + VMR but NOT the full gene body
#                   so the zone is tight around the features of interest)
# annot_gr       : GRanges for display (coloured bands + labels in plot)
# ══════════════════════════════════════════════════════════════════════════════

plot_region <- function(region_gr, annot_gr, peak_anchor_gr, meth,
                        title                     = "",
                        peak_span                 = 0.2){
  window_gr  <- region_gr[is.na(region_gr$name)]
  
  # ── x limits: set ONCE as plain integers ─────────────────────────────────
  x_min      <- as.integer(start(window_gr))
  x_max      <- as.integer(end(window_gr))
  chr_window <- as.character(seqnames(window_gr))
  
  ext_annot_dt <- as.data.table(annot_gr)
  ext_annot_dt[, mid := (start + end) / 2]
  
  # ── Order: gene first (background), then geneVMR, then TE (foreground) ────
  ext_annot_dt[, annot_order := fcase(
    annot_type == "gene",    1L,
    annot_type == "geneVMR", 2L,
    annot_type == "TE",      3L,
    default = 4L
  )]
  setorder(ext_annot_dt, annot_order)
  
  # ── Clip rectangles to plot window so gene body doesn't extend beyond ──────
  ext_annot_dt[, start := pmax(start, x_min)]
  ext_annot_dt[, end   := pmin(end,   x_max)]
  
  # recompute mid after clipping
  ext_annot_dt[, mid := (start + end) / 2]
  
  # ── Plot window data ───────────────────────────────────────────────────────
  hits_t3 <- findOverlaps(table3layers, window_gr)
  t3_sub  <- as.data.table(table3layers[queryHits(hits_t3)])[
    , .(seqnames, start, chr_pos, alpha_geomean, percentile)]
  hits_ic <- findOverlaps(interlayer_gr, window_gr)
  ic_sub  <- as.data.table(mcols(interlayer_gr[queryHits(hits_ic)]))
  result  <- merge(t3_sub, ic_sub, by = "chr_pos", all.x = TRUE)
  result  <- result[!is.na(result$start)]
  r_cols  <- c("r_Endo_Meso", "r_Endo_Ecto", "r_Meso_Ecto", "interlayer_r")
  result[, (r_cols) := lapply(.SD, abs), .SDcols = r_cols]
  
  # ── meth: parse positions ONCE, filter to window ──────────────────────────
  if (!"pos" %in% names(meth)) meth[, pos := as.integer(sub(".*_", "", cpg_site))]
  if (!"chr" %in% names(meth)) meth[, chr := sub("_.*", "", cpg_site)]
  meth_window <- meth[chr == chr_window & pos >= x_min & pos <= x_max]
  
  # replace peaks with empty table so add_annot doesn't break:
  peaks <- data.table(peak_start = integer(), peak_end = integer(),
                      n_cpgs = integer(), width_bp = integer(),
                      mean_smoothed = numeric())
  
  # ── Colours ───────────────────────────────────────────────────────────────
  annot_type_colours <- c(TE = "#b7d7e0", gene = "lightgrey", geneVMR = "#ebc30e")
  
  # ── X axis ────────────────────────────────────────────────────────────────
  x_breaks <- pretty(c(x_min, x_max), n = 6)
  x_breaks <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
  x_labels <- if ((x_max - x_min) < 500000) {
    function(x) paste0(round(x / 1e3, 1), " kb")
  } else {
    function(x) paste0(round(x / 1e6, 2), " Mb")
  }
  shared_x_top    <- scale_x_continuous(breaks = x_breaks, labels = x_labels,
                                        expand = c(0.01, 0))
  shared_x_bottom <- scale_x_continuous(name = "Position (hg38)", breaks = x_breaks,
                                        labels = x_labels, expand = c(0.01, 0))
  shared_coord      <- coord_cartesian(xlim = c(x_min, x_max))
  shared_coord_clip <- coord_cartesian(xlim = c(x_min, x_max), clip = "off")
  
  # ── Themes ────────────────────────────────────────────────────────────────
  th <- theme_bw(base_size = 10) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_line(), panel.grid.minor = element_blank(),
          legend.position = "none")
  th_bottom <- theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # ── Annotation helpers ────────────────────────────────────────────────────
  add_labels <- function(p, show_labels = FALSE) {
    if (nrow(peaks) > 0) {
      p <- p + geom_rect(
        data = peaks,
        aes(xmin = peak_start, xmax = peak_end, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = NA, colour = "#D55E00", linewidth = 0.8, linetype = "solid")
    }
    if (show_labels) {
      if (nrow(ext_annot_dt) > 0) {
        p <- p + geom_text(
          data = ext_annot_dt,
          aes(x = mid, y = Inf, label = name),
          inherit.aes = FALSE,
          vjust = 1.4, size = 2.8, colour = "black", fontface = "bold")
      }
    }
    p
  }
  
  # ── Panels ────────────────────────────────────────────────────────────────
  p_alpha <- add_labels(
    ggplot(result, aes(x = start, y = alpha_geomean)) +
      # ── background FIRST ───────────────────────────────────────────────
      geom_rect(data = ext_annot_dt,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
                    fill = annot_type),
                inherit.aes = FALSE, alpha = 0.8) +
      scale_fill_manual(values = annot_type_colours, guide = "none") +
      # ── data on top ────────────────────────────────────────────────────
      geom_smooth(colour = "black", fill = "grey80",
                  alpha = 0.2, na.rm = TRUE, span = peak_span) +
      geom_point(size = 1.2, alpha = 0.6, colour = "black", na.rm = TRUE) +
      scale_y_continuous("Pr(HV)\n(geomean)", limits = c(0, 1)) +
      shared_x_top + shared_coord_clip + ggtitle(title) + th,
    show_labels = TRUE)
  
  # ── P2: Percentile geomean — fix y to 0-100 ──────────────────────────────
  p_pct <- add_labels(
    ggplot(result, aes(x = start, y = percentile)) +
      # ── background FIRST ───────────────────────────────────────────────
      geom_rect(data = ext_annot_dt,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
                    fill = annot_type),
                inherit.aes = FALSE, alpha = 0.8) +
      scale_fill_manual(values = annot_type_colours, guide = "none") +
      # ── data on top ────────────────────────────────────────────────────
      geom_smooth(colour = "black", fill = "grey80",
                  alpha = 0.1, na.rm = TRUE, span = peak_span) +
      geom_point(size = 1.5, alpha = 0.6, colour = "black", na.rm = TRUE) +
      geom_hline(yintercept = 95, linetype = "dashed", colour = "firebrick") +
      scale_y_continuous("Percentile\n(geomean)",
                         limits = c(0, 100),
                         breaks = c(0, 25, 50, 75, 95, 100)) +
      shared_x_top + shared_coord + th)
  
  # ── Parse position from cpg_site ─────────────────────────────────────────
  if (!"pos" %in% names(meth)) {
    meth[, pos := as.integer(sub(".*_", "", cpg_site))]
  }
  if (!"chr" %in% names(meth)) {
    meth[, chr := sub("_.*", "", cpg_site)]
  }
  
  # ── Filter meth to plot window ────────────────────────────────────────────
  chr_window  <- as.character(seqnames(window_gr))
  meth_window <- meth[chr == chr_window & pos >= x_min & pos <= x_max]
  
  # filter meth to the plot window
  meth_window <- meth[
    pos >= x_min & pos <= x_max
  ]
  
  # need pos column
  if (!"pos" %in% names(meth)) {
    meth[, pos := as.integer(sub(".*_", "", cpg_site))]
  }
  
  chr_window <- sub("_.*", "", result$chr_pos[1])   # e.g. "chr8"
  meth_window <- meth[
    grepl(paste0("^", chr_window, "_"), cpg_site) &
      as.integer(sub(".*_", "", cpg_site)) >= x_min   &
      as.integer(sub(".*_", "", cpg_site)) <= x_max
  ]
  
  # ── P3: raw methylation — add_annot works fine (genomic x axis) ───────────
  p_raw <- add_labels(
    plot_raw_meth(copy(meth_window), x_min, x_max,
                  annot_dt     = ext_annot_dt,
                  annot_colours = annot_type_colours) +
      shared_x_top +
      geom_smooth(aes(group = germ_layer),
                  colour = "black", fill = "grey60",
                  alpha = 0.2, na.rm = TRUE, span = peak_span,
                  method = "loess", se = TRUE) +
      th + theme(legend.position = "right"))
  
  # ── Diagnose inter-layer coverage ─────────────────────────────────────────
  layer_per_patient <- meth_window[, .(n_layers = uniqueN(germ_layer)), by = patient_id]
  message(sprintf("  Patients with >=2 germ layers in window: %d / %d",
                  sum(layer_per_patient$n_layers >= 2),
                  nrow(layer_per_patient)))
  message("  Layer combinations present:")
  print(meth_window[, .N, by = germ_layer])
  
  # ── P4: inter-layer correlation — NO add_annot (categorical x axis) ───────
  p_corr <- add_labels(
    plot_percpg_interlayer_corr(copy(meth_window), x_min, x_max,
                                annot_dt      = ext_annot_dt,
                                annot_colours = annot_type_colours) +
      geom_smooth(colour = "black", fill = "grey80",
                  alpha = 0.2, na.rm = TRUE, span = peak_span) +
      shared_x_top +
      coord_cartesian(xlim = c(x_min, x_max)) +
      scale_y_continuous("|Pearson r|\n(inter-layer)",
                         limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1)))
  
  p_alpha / p_pct / p_raw / p_corr +
    plot_layout(heights = c(1, 1, 3, 1))
}

# ══════════════════════════════════════════════════════════════════════════════
# Genomic coordinates
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

SUSPECT_ALICE <- c(GRanges("chr1", IRanges(17757500, 17757900),name = "SUSPECT_1"),
                   GRanges("chr1", IRanges(17760700, 17761100),name = "SUSPECT_2"))

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

# ── Display annotation GRanges ────────────────────────────────────────────────
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
    end   = c(end(LTR41_hg38),   end(ACTL8_hg38), end(SUSPECT_ALICE))),
  name       = c("LTR41_1", "LTR41_2", "ACTL8", "SUSPECT_1", "SUSPECT_2"),
  annot_type = c("TE",      "TE",      "gene", "geneVMR", "geneVMR"))

# ══════════════════════════════════════════════════════════════════════════════
# Extract all CpGs in target regions → chr_pos format for S00 script
# ══════════════════════════════════════════════════════════════════════════════

# All regions to extract CpGs from, with a label
target_regions <- list(
  LY6SVMR      = LY6SVMR_hg38,
  MER11C       = MER11C_hg38,
  LY6S_AS1     = LY6S_AS1_hg38,
  LTR41_1      = LTR41_hg38[1],
  LTR41_2      = LTR41_hg38[2],
  ACTL8        = ACTL8_hg38,
  region1_window = region1_gr[is.na(region1_gr$name)],
  region2_window = region2_gr[is.na(region2_gr$name)]
)

# For each region, find overlapping CpGs in table3layers and extract chr_pos
cpg_list <- rbindlist(lapply(names(target_regions), function(nm) {
  gr   <- target_regions[[nm]]
  hits <- findOverlaps(table3layers, gr)
  cpgs <- table3layers$chr_pos[queryHits(hits)]
  data.table(chr_pos = cpgs, region = nm)
}))

# Unique CpGs across all regions (S00 needs one per line, no duplicates)
unique_cpgs <- unique(cpg_list$chr_pos)
unique_cpgs <- unique_cpgs[!is.na(unique_cpgs)]

message(sprintf("Total unique CpGs across all regions: %d", length(unique_cpgs)))
print(cpg_list[, .N, by = region])  # how many per region

# ── Write the cpg_list file for S00 ──────────────────────────────────────────
out_path <- here("B_MultiTissues/dataOut/targetCpGs_regions.txt")
writeLines(unique_cpgs, out_path)
message("Written to: ", out_path)

## In pchuckle:
# source /share/apps/source_files/python/python-3.13.0a6.source
# cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# python3 S00_extractRawMethylationForTargetCpG.py \
# --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/targetCpGs_regions.txt \
# --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/methylation_targetRegions.tsv \
# --minCov    10

meth <- fread(here("B_MultiTissues/dataOut/methylation_targetRegions.tsv"))
head(meth)

# ══════════════════════════════════════════════════════════════════════════════
# Prepare meth data for regional plots
# ══════════════════════════════════════════════════════════════════════════════

setDT(meth)

# Germ layer colours consistent with rest of script
germ_colours <- c(Endo = "#1D9E75", Meso = "#D85A30", Ecto = "#185FA5")

# ── (1) Raw methylation per person, coloured by germ layer ───────────────────

plot_raw_meth <- function(meth_sub, x_min, x_max, annot_dt,
                          annot_colours, title = "") {
  meth_sub[, pos := as.integer(sub(".*_", "", cpg_site))]
  
  ggplot(meth_sub, aes(x = pos, y = methylation, group = patient_id)) +
    # ── background rectangles FIRST ──────────────────────────────────────
    geom_rect(data = annot_dt,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
                  fill = annot_type),
              inherit.aes = FALSE, alpha = 0.8) +
    scale_fill_manual(values = annot_colours, guide = "none") +
    # ── data on top ───────────────────────────────────────────────────────
    geom_point(size = 0.8, alpha = 0.05) +
    scale_y_continuous("Methylation", limits = c(0, 1)) +
    scale_x_continuous(breaks = pretty(c(x_min, x_max), n = 6),
                       labels = function(x) paste0(round(x / 1e3, 1), " kb"),
                       expand = c(0.01, 0)) +
    coord_cartesian(xlim = c(x_min, x_max)) +
    facet_grid(germ_layer ~ .) +
    theme_bw(base_size = 10) +
    theme(legend.position = "right", panel.grid.minor = element_blank())
}

# ── (2) Intra-individual inter-layer Pearson correlation ─────────────────────
# For each patient who has samples from >= 2 different germ layers,
# correlate methylation across CpG positions between layers

# ── Compute per-CpG interlayer correlation from raw meth ─────────────────────
compute_percpg_interlayer_corr <- function(meth_sub) {
  
  meth_sub[, pos := as.integer(sub(".*_", "", cpg_site))]
  
  # mean per (patient, germ_layer, pos) — collapses multiple tissues same layer
  meth_agg <- meth_sub[, .(methylation = mean(methylation, na.rm = TRUE)),
                       by = .(patient_id, germ_layer, pos, cpg_site)]
  
  # wide: one column per germ layer
  wide <- dcast(meth_agg, patient_id + pos + cpg_site ~ germ_layer,
                value.var = "methylation")
  
  layers_present <- intersect(c("Endo", "Meso", "Ecto"), names(wide))
  if (length(layers_present) < 2) return(NULL)
  
  layer_pairs <- combn(layers_present, 2, simplify = FALSE)
  
  rbindlist(lapply(layer_pairs, function(pair) {
    l1 <- pair[1]; l2 <- pair[2]
    
    # per-CpG correlation across patients
    wide[, {
      idx <- !is.na(get(l1)) & !is.na(get(l2))
      if (sum(idx) >= 3) {
        list(r    = cor(get(l1)[idx], get(l2)[idx], method = "pearson"),
             n    = sum(idx),
             pair = paste(l1, l2, sep = "-"))
      } else {
        list(r = NA_real_, n = sum(idx), pair = paste(l1, l2, sep = "-"))
      }
    }, by = .(pos, cpg_site)]
  }))
}

# ── Plot per-CpG interlayer r along genomic position ─────────────────────────
plot_percpg_interlayer_corr <- function(meth_sub, x_min, x_max,
                                        annot_dt, annot_colours) {
  corr_dt <- compute_percpg_interlayer_corr(meth_sub)
  
  if (is.null(corr_dt) || nrow(corr_dt) == 0 || all(is.na(corr_dt$r))) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "Insufficient data\nfor per-CpG correlation",
                      hjust = 0.5, size = 2.8, colour = "grey50") +
             xlim(0, 1) + ylim(0, 1) + theme_void())
  }
  
  pair_colours <- c("Endo-Meso" = "#1D9E75",
                    "Endo-Ecto" = "#185FA5",
                    "Meso-Ecto" = "#D85A30")
  
  ggplot(corr_dt[!is.na(r)], aes(x = pos, y = abs(r), colour = pair)) +
    # ── background rectangles FIRST ──────────────────────────────────────
    geom_rect(data = annot_dt,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
                  fill = annot_type),
              inherit.aes = FALSE, alpha = 0.8) +
    scale_fill_manual(values = annot_colours, guide = "none") +
    # ── data on top ───────────────────────────────────────────────────────
    geom_point(size = 1.2, alpha = 0.7) +
    scale_colour_manual(values = pair_colours,
                        limits = names(pair_colours), name = "Layer pair") +
    scale_y_continuous("|Pearson r|\n(inter-layer)",
                       limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_continuous(breaks = pretty(c(x_min, x_max), n = 6),
                       labels = function(x) paste0(round(x / 1e3, 1), " kb"),
                       expand = c(0.01, 0)) +
    coord_cartesian(xlim = c(x_min, x_max)) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
          legend.position = "none")
}

# ══════════════════════════════════════════════════════════════════════════════
# Plots
# ══════════════════════════════════════════════════════════════════════════════
p1 <- do.call(plot_region, c(list(
  region_gr      = region1_gr,
  annot_gr       = annot_gr_region1,
  meth           = meth,
  title          = "Region 1: LY6S-VMR / MER11C (chr8)")))

p2 <- do.call(plot_region, c(list(
  region_gr      = region2_gr,
  meth           = meth,
  annot_gr       = annot_gr_region2,
  title          = "Region 2: ACTL8 / LTR41 (chr1)")))

# ══════════════════════════════════════════════════════════════════════════════
# Decay curves
# ══════════════════════════════════════════════════════════════════════════════
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

# ==========================================
## Compare our LTR41s with average LTR41s
# ==========================================
library(AnnotationHub)
ah        <- AnnotationHub()
rmskhg38  <- ah[["AH111333"]]

LTR41_all_hg38 <- rmskhg38[mcols(rmskhg38)$repName == "LTR41"]
message("Total LTR41 elements: ", length(LTR41_all_hg38))

# ── which LTR41s are our targets ─────────────────────────────────────────────
target_idx <- which(
  as.character(seqnames(LTR41_all_hg38)) == "chr1" &
    start(LTR41_all_hg38) %in% start(LTR41_hg38)
)

target_cols      <- c("LTR41_1" = "#E69F00", "LTR41_2" = "#CC79A7")
mean_ltr41_width <- mean(width(LTR41_all_hg38))
ltr41_widths     <- width(LTR41_hg38)

# ── Loess ribbon fitter ───────────────────────────────────────────────────────
fit_loess_ribbon <- function(dt, x_col, y_col, span = 0.5, n_out = 300,
                             max_n = 20000) {
  df <- dt[!is.na(get(y_col)), .(x = get(x_col), y = get(y_col))]
  if (nrow(df) < 10) return(NULL)
  if (nrow(df) > max_n) { set.seed(42); df <- df[sample(.N, max_n)] }
  fit   <- loess(y ~ x, data = df, span = span)
  x_seq <- seq(min(df$x), max(df$x), length.out = n_out)
  pred  <- predict(fit, newdata = data.frame(x = x_seq), se = TRUE)
  data.table(x    = x_seq,
             y    = pred$fit,
             ymin = pred$fit - 1.96 * pred$se.fit,
             ymax = pred$fit + 1.96 * pred$se.fit)
}

# ── Plot function ─────────────────────────────────────────────────────────────
make_ltr41_plot <- function(smooth_mean, smooth_target, smooth_indiv,
                            y_lab, y_limits, y_breaks,
                            x_limits = NULL,
                            extra_layers = NULL,
                            title = NULL, subtitle = NULL) {
  p <- ggplot() +
    geom_line(data = smooth_indiv,
              aes(x = x, y = y, group = ltr_idx),
              colour = "grey70", alpha = 0.2, linewidth = 0.2) +
    geom_ribbon(data = smooth_mean,
                aes(x = x, ymin = ymin, ymax = ymax),
                fill = "grey70", alpha = 0.3) +
    geom_line(data = smooth_mean,
              aes(x = x, y = y),
              colour = "grey30", linewidth = 1) +
    geom_ribbon(data = smooth_target,
                aes(x = x, ymin = ymin, ymax = ymax, fill = ltr_name),
                alpha = 0.25) +
    geom_line(data = smooth_target,
              aes(x = x, y = y, colour = ltr_name),
              linewidth = 1.2) +
    scale_colour_manual(values = target_cols, name = "Target LTR41") +
    scale_fill_manual(values   = target_cols, name = "Target LTR41") +
    geom_vline(xintercept = 0,
               linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = mean_ltr41_width,
               linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = ltr41_widths[1],
               colour = target_cols["LTR41_1"],
               linetype = "dotted", linewidth = 0.8) +
    geom_vline(xintercept = ltr41_widths[2],
               colour = target_cols["LTR41_2"],
               linetype = "dotted", linewidth = 0.8) +
    annotate("text", x = 0, y = max(y_limits),
             label = "TE start", hjust = 0, vjust = 1.5,
             size = 2.5, colour = "grey40") +
    annotate("text", x = mean_ltr41_width, y = max(y_limits),
             label = sprintf("TE end\n(mean=%d bp)", round(mean_ltr41_width)),
             hjust = 1, vjust = 1.5, size = 2.5, colour = "grey40") +
    scale_x_continuous("Position relative to LTR41 start (bp)",
                       labels = function(x) paste0(round(x / 1e3, 1), " kb")) +
    scale_y_continuous(y_lab, limits = y_limits, breaks = y_breaks) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  
  if (!is.null(x_limits))
    p <- p + coord_cartesian(xlim = x_limits)
  if (!is.null(extra_layers)) p <- p + extra_layers
  if (!is.null(title))        p <- p + ggtitle(title, subtitle = subtitle)
  p
}

# ── Build dataset + smoothers for each zoom level ────────────────────────────
zoom_levels <- list(
  "1kb"  = 1000,
  "2kb"  = 2000,
  "3kb"  = 3000
)

# Build once with max padding (3kb), reuse for all zoom levels
LTR41_windows <- GRanges(
  seqnames = seqnames(LTR41_all_hg38),
  ranges   = IRanges(
    start = pmax(start(LTR41_all_hg38) - 3000, 1),
    end   = end(LTR41_all_hg38) + 3000
  )
)

message("Extracting CpGs for ", length(LTR41_all_hg38), " LTR41 elements...")
hits_all <- findOverlaps(table3layers, LTR41_windows)

cpg_dt <- as.data.table(table3layers[queryHits(hits_all)])[
  , .(start, alpha_geomean, percentile, chr_pos)]
cpg_dt[, ltr_idx   := subjectHits(hits_all)]
cpg_dt[, ltr_start := start(LTR41_all_hg38)[ltr_idx]]
cpg_dt[, rel_pos   := start - ltr_start]
cpg_dt[, is_target := ltr_idx %in% target_idx]
cpg_dt[, ltr_name  := fcase(
  ltr_idx == target_idx[1], "LTR41_1",
  ltr_idx == target_idx[2], "LTR41_2",
  default = "other"
)]

ltr41_cpg_dt <- cpg_dt
message(sprintf("  %d CpG observations across %d elements (%d unique CpGs)",
                nrow(ltr41_cpg_dt),
                uniqueN(ltr41_cpg_dt$ltr_idx),
                uniqueN(ltr41_cpg_dt$chr_pos)))

# subsample for individual grey lines
set.seed(42)
non_target_idx <- setdiff(unique(ltr41_cpg_dt$ltr_idx), target_idx)
sample_idx     <- sample(non_target_idx, min(200, length(non_target_idx)))

# ── Pre-compute smoothers once on full ±3kb range ────────────────────────────
message("Computing smoothers...")

smooth_mean_alpha <- fit_loess_ribbon(
  ltr41_cpg_dt[is_target == FALSE], "rel_pos", "alpha_geomean")
smooth_mean_pct   <- fit_loess_ribbon(
  ltr41_cpg_dt[is_target == FALSE], "rel_pos", "percentile")

smooth_target_alpha <- rbindlist(lapply(unique(target_idx), function(ti) {
  dt <- ltr41_cpg_dt[ltr_idx == ti]
  s  <- fit_loess_ribbon(dt, "rel_pos", "alpha_geomean")
  if (is.null(s)) return(NULL)
  s[, ltr_name := dt$ltr_name[1]]; s
}))

smooth_target_pct <- rbindlist(lapply(unique(target_idx), function(ti) {
  dt <- ltr41_cpg_dt[ltr_idx == ti]
  s  <- fit_loess_ribbon(dt, "rel_pos", "percentile")
  if (is.null(s)) return(NULL)
  s[, ltr_name := dt$ltr_name[1]]; s
}))

smooth_indiv_alpha <- rbindlist(lapply(sample_idx, function(ti) {
  dt <- ltr41_cpg_dt[ltr_idx == ti & !is.na(alpha_geomean)]
  if (nrow(dt) < 5) return(NULL)
  fit   <- loess(alpha_geomean ~ rel_pos, data = dt, span = 0.5)
  x_seq <- seq(min(dt$rel_pos), max(dt$rel_pos), length.out = 100)
  data.table(x = x_seq, y = predict(fit, x_seq), ltr_idx = ti)
}))

smooth_indiv_pct <- rbindlist(lapply(sample_idx, function(ti) {
  dt <- ltr41_cpg_dt[ltr_idx == ti & !is.na(percentile)]
  if (nrow(dt) < 5) return(NULL)
  fit   <- loess(percentile ~ rel_pos, data = dt, span = 0.5)
  x_seq <- seq(min(dt$rel_pos), max(dt$rel_pos), length.out = 100)
  data.table(x = x_seq, y = predict(fit, x_seq), ltr_idx = ti)
}))

message("Smoothers done.")

# ── Build 3 zoom-level pairs ──────────────────────────────────────────────────
ltr41_plots <- lapply(names(zoom_levels), function(zoom_name) {
  pad      <- zoom_levels[[zoom_name]]
  x_limits <- c(-pad, mean_ltr41_width + pad)
  
  p_alpha <- make_ltr41_plot(
    smooth_mean   = smooth_mean_alpha,
    smooth_target = smooth_target_alpha,
    smooth_indiv  = smooth_indiv_alpha,
    y_lab    = "Pr(HV) geomean",
    y_limits = c(0, 1),
    y_breaks = c(0, 0.25, 0.5, 0.75, 1),
    x_limits = x_limits,
    title    = sprintf("LTR41 ± %s", zoom_name),
    subtitle = sprintf("grey = all LTR41 (n=%d, n=%d shown) | coloured = our 2 LTR41s",
                       length(LTR41_all_hg38), length(sample_idx))
  )
  
  p_pct <- make_ltr41_plot(
    smooth_mean   = smooth_mean_pct,
    smooth_target = smooth_target_pct,
    smooth_indiv  = smooth_indiv_pct,
    y_lab         = "Percentile (geomean)",
    y_limits      = c(0, 100),
    y_breaks      = c(0, 25, 50, 75, 95, 100),
    x_limits      = x_limits,
    extra_layers  = geom_hline(yintercept = 95,
                               linetype = "dashed", colour = "firebrick")
  )
  
  p_alpha / p_pct + plot_layout(guides = "collect") &
    theme(legend.position = "right")
})

# ── Combine 3 zoom levels side by side + legend ─────────────────────────

# ── Extract legend from a single ggplot panel (not the patchwork) ─────────────
legend_plot <- make_ltr41_plot(
  smooth_mean   = smooth_mean_alpha,
  smooth_target = smooth_target_alpha,
  smooth_indiv  = smooth_indiv_alpha,
  y_lab    = "Pr(HV) geomean",
  y_limits = c(0, 1),
  y_breaks = c(0, 0.25, 0.5, 0.75, 1),
  x_limits = c(-3000, mean_ltr41_width + 3000)
)

legend <- cowplot::get_legend(legend_plot)

# ── Remove legend from all patchwork plots ────────────────────────────────────
ltr41_plots_nolegend <- lapply(ltr41_plots, function(x) {
  x & theme(legend.position = "none")
})

# ══════════════════════════════════════════════════════════════════════════════
# LD-style pairwise methylation correlation heatmaps
#
# CONCEPT: In genetics, LD (r²) measures how often two SNP alleles co-occur
# across individuals. Here we do the exact analogue for methylation:
# for each pair of CpGs, we compute Pearson r of their methylation values
# across 135 patients. High r = the two CpGs tend to be methylated (or
# unmethylated) together across individuals → epigenetic co-regulation.
# Negative r = when one is methylated, the other tends not to be → competing
# states. This is displayed as an upper-triangle heatmap ordered by genomic
# position, exactly like a standard LD block plot.
# ══════════════════════════════════════════════════════════════════════════════

# ── Region windows ────────────────────────────────────────────────────────────
ld_region1_gr <- GRanges("chr8",
                         IRanges(start = min(start(LY6SVMR_hg38), start(MER11C_hg38)) - 3000,
                                 end   = max(end(LY6SVMR_hg38),   end(MER11C_hg38))   + 3000))

ld_region2_gr <- GRanges("chr1",
                         IRanges(start = min(start(LTR41_hg38)) - 3000,
                                 end   = max(end(LTR41_hg38))   + 3000))

feature_ticks_region1 <- data.table(
  pos   = c(start(LY6SVMR_hg38), end(LY6SVMR_hg38),
            start(MER11C_hg38),  end(MER11C_hg38)),
  label = c("VMR start", "VMR end", "MER11C start", "MER11C end")
)

feature_ticks_region2 <- data.table(
  pos   = c(start(SUSPECT_ALICE[1]), end(SUSPECT_ALICE[1]),
            start(SUSPECT_ALICE[2]), end(SUSPECT_ALICE[2]),
            start(LTR41_hg38[1]), end(LTR41_hg38[1]),
            start(LTR41_hg38[2]), end(LTR41_hg38[2])),
  label = c("SUSPECT_1 start", "SUSPECT_1 end",
            "SUSPECT_2 start", "SUSPECT_2 end",
            "LTR41_1 start", "LTR41_1 end",
            "LTR41_2 start", "LTR41_2 end")
)

# ── build_meth_matrix ─────────────────────────────────────────────────────────
# Extracts raw methylation from `meth` for all CpGs in the window,
# averages across tissues per patient (so each patient contributes one value
# per CpG regardless of how many tissues they have), then returns a
# CpG × patient matrix sorted by genomic position.
# Rows with fewer than min_patients non-NA values are dropped.
build_meth_matrix <- function(window_gr, meth, min_patients = 5) {
  
  x_min      <- as.integer(start(window_gr))
  x_max      <- as.integer(end(window_gr))
  chr_window <- as.character(seqnames(window_gr))
  
  if (!"pos" %in% names(meth)) meth[, pos := as.integer(sub(".*_", "", cpg_site))]
  if (!"chr" %in% names(meth)) meth[, chr := sub("_.*", "", cpg_site)]
  
  mw <- meth[chr == chr_window & pos >= x_min & pos <= x_max]
  message(sprintf("  Window %s:%d-%d — %d rows, %d CpGs, %d patients",
                  chr_window, x_min, x_max,
                  nrow(mw), uniqueN(mw$cpg_site), uniqueN(mw$patient_id)))
  if (nrow(mw) == 0) return(NULL)
  
  agg  <- mw[, .(methylation = mean(methylation, na.rm = TRUE)),
             by = .(cpg_site, pos, patient_id)]
  wide <- dcast(agg, pos + cpg_site ~ patient_id, value.var = "methylation")
  setorder(wide, pos)
  mat  <- as.matrix(wide[, -(1:2)])
  keep <- rowSums(!is.na(mat)) >= min_patients
  
  list(mat       = mat[keep, , drop = FALSE],
       positions = wide$pos[keep],
       cpg_sites = wide$cpg_site[keep])
}

# ── compute_pairwise_cor ──────────────────────────────────────────────────────
# THE CORE LD-EQUIVALENT CALCULATION.
# Takes a CpG × sample matrix and computes all pairwise Pearson correlations
# between rows (CpGs) across columns (samples = patients or layers).
# This is exactly r from LD analysis, applied to methylation instead of alleles:
#   r(CpG_i, CpG_j) = cor(meth[i, ], meth[j, ], use="pairwise.complete.obs")
# Rows with near-zero variance (CpGs uniformly methylated across all samples)
# are removed first as they produce undefined correlations.
compute_pairwise_cor <- function(mat) {
  rv  <- apply(mat, 1, var, na.rm = TRUE)
  mat <- mat[!is.na(rv) & rv > 1e-10, , drop = FALSE]
  if (nrow(mat) < 2) return(NULL)
  cor_mat <- suppressWarnings(
    cor(t(mat), use = "pairwise.complete.obs", method = "pearson"))
  keep    <- rowSums(!is.na(cor_mat)) > 0
  cor_mat <- cor_mat[keep, keep, drop = FALSE]
  if (nrow(cor_mat) < 2) return(NULL)
  cor_mat
}

# ── trim_positions ────────────────────────────────────────────────────────────
# Removes genomic positions corresponding to zero-variance CpG rows,
# keeping positions in sync with the rows that survive into the cor matrix.
trim_positions <- function(mat, positions) {
  rv   <- apply(mat, 1, var, na.rm = TRUE)
  keep <- !is.na(rv) & rv > 1e-10
  positions[keep]
}

# ── make_cpg_density_strip ────────────────────────────────────────────────────
# Builds a barplot showing how many CpGs fall in each genomic bin across the
# window. This is shown below the heatmap so the reader can interpret
# dense vs sparse regions of the rank-ordered x/y axis.
make_cpg_density_strip <- function(positions, feature_ticks = NULL, n_bins = 50) {
  if (length(positions) == 0 || !is.finite(min(positions))) return(ggplot() + theme_void())
  
  pos_min <- min(positions); pos_max <- max(positions)
  breaks  <- seq(pos_min, pos_max, length.out = n_bins + 1)
  counts  <- hist(positions, breaks = breaks, plot = FALSE)$counts
  dt      <- data.table(bin   = seq_len(n_bins), count = counts,
                        pos   = (breaks[-length(breaks)] + breaks[-1]) / 2)
  
  label_pos  <- pretty(positions, n = 5)
  label_pos  <- label_pos[label_pos >= pos_min & label_pos <= pos_max]
  label_bins <- sapply(label_pos, function(p) which.min(abs(dt$pos - p)))
  label_text <- paste0(round(label_pos / 1e3, 1), " kb")
  
  p <- ggplot(dt, aes(x = bin, y = count, fill = count)) +
    geom_col(width = 1) +
    scale_fill_gradient(low = "grey85", high = "grey20", guide = "none") +
    scale_x_continuous(breaks = label_bins, labels = label_text, expand = c(0,0)) +
    scale_y_continuous(expand = c(0, 0), name = "CpG\ndensity") +
    theme_bw(base_size = 7) +
    theme(axis.text.x  = element_text(size = 6, angle = 45, hjust = 1),
          axis.text.y  = element_text(size = 5),
          axis.title.x = element_blank(),
          panel.grid   = element_blank(),
          plot.margin  = margin(0,0,0,0))
  
  if (!is.null(feature_ticks) && nrow(feature_ticks) > 0) {
    ft_bins <- sapply(feature_ticks$pos, function(fp) which.min(abs(dt$pos - fp)))
    p <- p + geom_vline(xintercept = ft_bins, colour = "black",
                        linetype = "dashed", linewidth = 0.4)
  }
  p
}

# ── make_ld_heatmap ───────────────────────────────────────────────────────────
# Renders the pairwise correlation matrix as an upper-triangle heatmap.
# X and Y axes use rank indices (1 to n CpGs) rather than genomic coordinates
# so every CpG gets equal visual space regardless of genomic spacing —
# exactly as in standard LD plots. Axis tick labels are converted back to
# genomic kb positions for readability. Feature boundaries are shown as
# dashed lines at the nearest CpG rank. Colour scale is diverging:
# blue = negative r (anti-correlated), white = r≈0, orange/red = positive r.
make_ld_heatmap <- function(cor_mat, positions, title, feature_ticks = NULL) {
  
  if (is.null(cor_mat) || nrow(cor_mat) < 2) {
    return(ggplot() +
             annotate("text", x=0.5, y=0.5,
                      label = sprintf("Insufficient data\n(n=%d CpGs)",
                                      if (is.null(cor_mat)) 0L else nrow(cor_mat)),
                      hjust=0.5, size=3, colour="grey50") +
             theme_void() + ggtitle(title))
  }
  
  n   <- nrow(cor_mat)
  idx <- which(upper.tri(cor_mat, diag = TRUE), arr.ind = TRUE)
  dt  <- data.table(i = idx[,1], j = idx[,2], r = cor_mat[idx])
  dt  <- dt[!is.na(r)]
  
  label_idx  <- round(seq(1, n, length.out = min(6, n)))
  label_text <- paste0(round(positions[label_idx] / 1e3, 1), " kb")
  
  p <- ggplot(dt, aes(x = i, y = j, fill = r)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_gradientn(
      colours  = c("blue", "dodgerblue", "white", "orange", "red"),
      rescaler = ~ scales::rescale_mid(.x, mid = 0),
      limits   = c(-1, 1),
      breaks   = c(-1, -0.5, 0, 0.5, 1),
      name     = "Pearson r",
      na.value = "grey90"
    ) +
    scale_x_continuous(breaks = label_idx, labels = label_text, expand = c(0,0)) +
    scale_y_continuous(breaks = label_idx, labels = label_text, expand = c(0,0)) +
    coord_fixed() +
    ggtitle(title) +
    theme_bw(base_size = 9) +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y     = element_text(size = 7),
          axis.title      = element_blank(),
          panel.grid      = element_blank(),
          legend.position = "right")
  
  if (!is.null(feature_ticks) && nrow(feature_ticks) > 0) {
    ft_ranks <- sapply(feature_ticks$pos, function(p) which.min(abs(positions - p)))
    ft <- data.table(rank = ft_ranks, label = feature_ticks$label)
    p  <- p +
      geom_vline(data = ft, aes(xintercept = rank), colour = "grey",
                 linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE) +
      geom_hline(data = ft, aes(yintercept = rank), colour = "grey",
                 linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE) +
      annotate("text", x = ft_ranks, y = n, label = ft$label,
               angle = 90, hjust = 1, vjust = -0.3, size = 2, colour = "black")
  }
  p
}

# ── make_ld_panel ─────────────────────────────────────────────────────────────
# Stacks the LD heatmap on top of the CpG density strip (15:1 height ratio),
# giving a compact panel that shows both the correlation structure and
# where CpGs are actually located in the genome.
make_ld_panel <- function(cor_mat, positions, title, feature_ticks = NULL) {
  p_heat <- make_ld_heatmap(cor_mat, positions, title, feature_ticks)
  if (length(positions) >= 2 && is.finite(min(positions))) {
    p_strip <- make_cpg_density_strip(positions, feature_ticks)
    return(p_heat / p_strip + plot_layout(heights = c(15, 1)))
  }
  p_heat
}

# ── make_region_ld_plots ──────────────────────────────────────────────────────
# raw methylation LD — 135 patients × CpGs
# Returns a named list of three patchwork panels.
make_region_ld_plots <- function(window_gr, meth, region_name,
                                 feature_ticks = NULL, min_samples = 5) {
  
  x_min      <- as.integer(start(window_gr))
  x_max      <- as.integer(end(window_gr))
  chr_window <- as.character(seqnames(window_gr))
  
  if (!"pos" %in% names(meth)) meth[, pos := as.integer(sub(".*_", "", cpg_site))]
  if (!"chr" %in% names(meth)) meth[, chr := sub("_.*", "", cpg_site)]
  mw <- meth[chr == chr_window & pos >= x_min & pos <= x_max]
  
  message(sprintf("  %s: %d rows, %d CpGs, %d patients",
                  region_name, nrow(mw), uniqueN(mw$cpg_site),
                  uniqueN(mw$patient_id)))
  if (nrow(mw) == 0) { message("  No data"); return(list()) }
  
  # raw methylation
  mw_agg    <- mw[, .(methylation = mean(methylation, na.rm=TRUE)),
                  by = .(cpg_site, pos, patient_id)][!is.na(methylation)]
  meth_wide <- dcast(mw_agg, pos + cpg_site ~ patient_id, value.var = "methylation")
  setorder(meth_wide, pos)
  meth_mat       <- as.matrix(meth_wide[, -(1:2)])
  keep_m         <- rowSums(!is.na(meth_mat)) >= min_samples
  meth_mat       <- meth_mat[keep_m, , drop = FALSE]
  positions_meth <- meth_wide$pos[keep_m]
  cor_meth       <- compute_pairwise_cor(meth_mat)
  positions_meth <- trim_positions(meth_mat, positions_meth)
  if (!is.null(cor_meth))
    positions_meth <- positions_meth[seq_len(nrow(cor_meth))]
  message(sprintf("  Raw meth: %d CpGs in cor matrix", length(positions_meth)))
  
  list(
    meth = make_ld_panel(cor_meth, positions_meth,
                         title         = sprintf("%s\nRaw methylation", region_name),
                         feature_ticks = feature_ticks)
  )
}

# ── Run ───────────────────────────────────────────────────────────────────────
ld_region1 <- make_region_ld_plots(ld_region1_gr, meth, "LY6S-VMR / MER11C",
                                   feature_ticks_region1)
ld_region2 <- make_region_ld_plots(ld_region2_gr, meth, "ACTL8 / LTR41",
                                   feature_ticks_region2)

##################
## Final figure ##
##################

# Regional overview (p1, p2) — what does Pr(HV) look like at each locus
# LD heatmaps (ld_region1, ld_region2) — methylation co-variation structure
# LTR41 genome-wide comparison (ltr41_plots) — how our LTR41s compare to all others
# Decay curves (p4, p5) — ME signature strength

# ── Combine: 3 zoom plots + shared legend on the right ───────────────────────
p6 <- cowplot::plot_grid(
  cowplot::plot_grid(plotlist = ltr41_plots_nolegend, ncol = 3,
                     labels = c("e.i", "e.ii", "e.iii")),
  legend,
  ncol       = 2,
  rel_widths = c(1, 0.08)
)

# ── Extract single LD plots from lists ───────────────────────────────────────
p_ld1 <- ld_region1$meth
p_ld2 <- ld_region2$meth

# ── Row 1: regional Pr(HV) overview ──────────────────────────────────────────
row1 <- cowplot::plot_grid(
  p1, p2,
  ncol   = 2,
  labels = c("a", "b")
)

# ── Row 2: LD heatmaps ───────────────────────────────────────────────────────
row2 <- cowplot::plot_grid(
  p_ld1, p_ld2,
  ncol   = 2,
  labels = c("c", "d")
)

# ── Row 3: decay curves ───────────────────────────────────────────────────────
row3 <- cowplot::plot_grid(
  p4, p5,
  ncol   = 2,
  labels = c("e", "f") # wait, used below — reorder
)

# ── Final assembly ────────────────────────────────────────────────────────────
# Relabel to keep narrative order:
# a/b = region plots, c/d = LD heatmaps, e = LTR41 genome-wide, f/g = decay

row1 <- cowplot::plot_grid(p1, p2,    ncol=2, labels=c("a","b"))
row2 <- cowplot::plot_grid(p_ld1, p_ld2, ncol=2, labels=c("c","d"))
row3 <- p6   # already has labels e.i/ii/iii internally; add outer label
row3_labelled <- cowplot::plot_grid(
  cowplot::plot_grid(NULL, labels="e"),  # blank label holder
  row3,
  ncol       = 1,
  rel_heights = c(0.05, 1)
)
row4 <- cowplot::plot_grid(p4, p5, ncol=2, labels=c("f","g"))

final_plot <- cowplot::plot_grid(
  row1,
  row2,
  row3_labelled,
  row4,
  ncol        = 1,
  rel_heights = c(
    3,    # row1: region plots are tall (4 panels stacked each)
    1.8,  # row2: LD heatmaps — square
    2,  # row3: LTR41 comparison — 2 panels stacked × 3 zooms
    1   # row4: decay curves — compact
  )
)

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/ACTL8_LTR41.pdf"),
  plot     = final_plot,
  width    = 14,
  height   = 22
)
