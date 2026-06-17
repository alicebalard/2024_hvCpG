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

# ── Region GRanges (plot windows) ────────────────────────────────────────────
region1_gr <- c(
  LY6SVMR_hg38, MER11C_hg38,
  GRanges("chr8",
          IRanges(start = min(start(LY6SVMR_hg38), start(MER11C_hg38)) - 3000,
                  end   = max(end(LY6SVMR_hg38),   end(MER11C_hg38))   + 3000),
          name = NA))

region2_gr <- c(
  LTR41_hg38, ACTL8_hg38,
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
    start = c(start(LTR41_hg38), start(ACTL8_hg38)),
    end   = c(end(LTR41_hg38),   end(ACTL8_hg38))),
  name       = c("LTR41_1", "LTR41_2", "ACTL8"),
  annot_type = c("TE",      "TE",      "gene"))

# ── Peak anchor GRanges (detection zone only — no full gene body) ─────────────
# Region 1: VMR + MER11C only; gene body excluded to avoid biasing baseline
peak_anchor_region1 <- GRanges(
  seqnames   = "chr8",
  ranges     = IRanges(
    start = c(start(LY6SVMR_hg38), start(MER11C_hg38)),
    end   = c(end(LY6SVMR_hg38),   end(MER11C_hg38))),
  name       = c("LY6S-VMR", "MER11C"),
  annot_type = c("geneVMR",  "TE"))

# Region 2: LTR41 elements only; ACTL8 gene body excluded
peak_anchor_region2 <- GRanges(
  seqnames   = "chr1",
  ranges     = IRanges(start = start(LTR41_hg38), end = end(LTR41_hg38)),
  name       = c("LTR41_1", "LTR41_2"),
  annot_type = c("TE", "TE"))

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
  region2_window = region2_gr[is.na(region2_gr$name)],
  region3_window = region3_gr[is.na(region3_gr$name)]
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
          legend.position = "none")    # theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
    #       axis.ticks.x = element_line(), panel.grid.minor = element_blank(),
    #       legend.position = "right")
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

ggplot2::ggsave(
  filename = here::here("B_MultiTissues/dataOut/figures/ACTL8_LTR41.pdf"),
  plot = cowplot::plot_grid(
    cowplot::plot_grid(p1, p2, ncol = 2, labels = c("a", "b")),
    cowplot::plot_grid(p4, p5, ncol = 2, labels = c("c", "d")),
    nrow = 2, rel_heights = c(3, 1)),
  width = 12, height = 12
)
