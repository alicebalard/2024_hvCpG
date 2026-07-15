#####################################################################
## S13_plotMattTargets.R
## Per-locus plots (same panels as S11's plot_candidate_locus) for
## Matt's target REGIONS, defined by explicit coordinates (not symbols).
##
## Self-contained: needs only functions.R, table3layers, wideFull,
## the liftOver `chain` (from prepPreviousSIV.R), and an extracted
## methylation TSV covering these regions' CpGs.
##
## Produces: figures/lociAlice/locus_Matt_<REGION>.png
#####################################################################

library(here)
source(here("B_MultiTissues", "quiet_library.R"))
if (!exists("functionsLoaded"))     source(here("B_MultiTissues/03_exploreResults", "functions.R"))
if (!exists("previousSIVprepared")) source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))  # for `chain`
if (!exists("table3layers"))        load(here("gitignore/fullTable3layers.Rda"))
if (!exists("wideFull"))            wideFull <- readRDS(here("gitignore/wide_script10_3layers_full.RDS"))

if (is.null(table3layers$percentile)) {
  table3layers$percentile <- ecdf(table3layers$alpha_geomean)(table3layers$alpha_geomean) * 100
}

library(GenomicRanges)
library(GenomeInfoDb)

# ══════════════════════════════════════════════════════════════════════════════
# 0. Target coordinates  (from S07; hg19 ones are lifted to hg38)
# ══════════════════════════════════════════════════════════════════════════════
LY6SVMR_hg19 <- GRanges("chr8", IRanges(144120106, 144120706), name = "LY6S-VMR")
LY6SVMR_hg38 <- unlist(liftOver(LY6SVMR_hg19, chain)); LY6SVMR_hg38$name <- "LY6S-VMR"
LY6S_AS1_hg38 <- GRanges("chr8", IRanges(143039209, 143054303), name = "LY6S_AS1")
MER11C_hg38   <- GRanges("chr8", IRanges(143040739, 143041811), name = "MER11C")

LTR41_hg19 <- GRanges("chr1",
                      IRanges(start = c(18081648, 18085651),
                              end   = c(18082190, 18086109)),
                      name = c("LTR41_1", "LTR41_2"))
LTR41_hg38 <- unlist(liftOver(LTR41_hg19, chain)); LTR41_hg38$name <- c("LTR41_1", "LTR41_2")
ACTL8_hg38 <- GRanges("chr1", IRanges(17755333, 17827063), name = "ACTL8")

library(readxl)
dmr_dt <- as.data.table(read_excel(
  here("B_MultiTissues/dataIn/DEGCAGS_intersect_repeats_Alice_withalpha_GITIGNORE.xlsx")))
# coordinates in seqnames/start/end are already hg38 (hg19_* kept separately)
dmr_gr <- GRanges(dmr_dt$seqnames,
                  IRanges(dmr_dt$start, dmr_dt$end),
                  pct_change = dmr_dt$`X.change`,
                  padj       = dmr_dt$padj,
                  TE_family  = dmr_dt$TE_family)
seqlevelsStyle(dmr_gr) <- "UCSC"

# collect every target as one flat GRanges (one row per "locus" to plot)
targets_gr <- c(LY6SVMR_hg38, MER11C_hg38, LY6S_AS1_hg38,
                LTR41_hg38, ACTL8_hg38)
seqlevelsStyle(targets_gr) <- "UCSC"
targets_gr <- keepStandardChromosomes(targets_gr, pruning.mode = "coarse")

# ── shape targets into the table getFullPlotsTarget expects ──────────────────
# plot_candidate_locus treats `gene`/gene_start/gene_end as the plot window,
# so here "gene" = region name, and gene_start/end = region bounds.
all_candidates_Matt <- data.table(
  gene       = targets_gr$name,
  chr        = as.character(seqnames(targets_gr)),
  gene_start = start(targets_gr),
  gene_end   = end(targets_gr))
# small features (VMR, TEs) are tiny — pad so the window is meaningful
all_candidates_Matt[, span := gene_end - gene_start]
all_candidates_Matt[span < 2000, `:=`(                 # widen features < 2 kb
  gene_start = gene_start - 1500, gene_end = gene_end + 1500)]
all_candidates_Matt[, span := NULL]

# ── combined REGION windows: gene kept only `gene_keep` bp from its start ────
pad       <- 5000   # flank on each side
gene_keep <- 5000   # how much of the long gene to keep, from its start

## Region 1: LY6S (first 5 kb) + LY6S-VMR + MER11C, ±5 kb
ly6s_tss  <- start(LY6S_AS1_hg38)                     # gene assumed + strand; see note
r1_left   <- min(start(LY6SVMR_hg38), start(MER11C_hg38), ly6s_tss) - pad
r1_right  <- max(end(LY6SVMR_hg38),  end(MER11C_hg38),  ly6s_tss + gene_keep) + pad

region1 <- data.table(
  gene       = "Region1_LY6SVMR_MER11C",
  chr        = "chr8",
  gene_start = r1_left,
  gene_end   = r1_right)

## Region 2: ACTL8 (first 5 kb) + LTR41_1 + LTR41_2, ±5 kb
actl8_tss <- start(ACTL8_hg38)
r2_left   <- min(start(LTR41_hg38), actl8_tss) - pad    # 5 kb upstream of the most upstream LTR
r2_right  <- max(end(LTR41_hg38),   actl8_tss + gene_keep) + pad

region2 <- data.table(
  gene       = "Region2_ACTL8_LTR41",
  chr        = "chr1",
  gene_start = r2_left,
  gene_end   = r2_right)

message(sprintf("Region1: chr8:%s-%s  (%.1f kb)",
                format(r1_left, big.mark=","), format(r1_right, big.mark=","),
                (r1_right - r1_left)/1e3))
message(sprintf("Region2: chr1:%s-%s  (%.1f kb)",
                format(r2_left, big.mark=","), format(r2_right, big.mark=","),
                (r2_right - r2_left)/1e3))

# keep the individual targets AND add the two combined regions
all_candidates_Matt <- rbind(all_candidates_Matt, region1, region2, fill = TRUE)

candidates_gr_Matt <- unique(GRanges(
  seqnames = all_candidates_Matt$chr,
  ranges   = IRanges(all_candidates_Matt$gene_start - 5000,
                     all_candidates_Matt$gene_end   + 5000),
  name     = all_candidates_Matt$gene))

candidates_gr_Matt <- unique(GRanges(
  seqnames = all_candidates_Matt$chr,
  ranges   = IRanges(all_candidates_Matt$gene_start - 5000,
                     all_candidates_Matt$gene_end   + 5000),
  name     = all_candidates_Matt$gene))

# ══════════════════════════════════════════════════════════════════════════════
# 1. Write CpG list to extract (then run python extractor on the cluster)
# ══════════════════════════════════════════════════════════════════════════════
getCpG4pchuckle <- function(candidates_gr, layer){
  cpg_extract_list <- suppressWarnings(rbindlist(lapply(
    seq_along(candidates_gr), function(i) {
      gr <- candidates_gr[i]; hits <- findOverlaps(table3layers, gr)
      if (length(hits) == 0) return(NULL)
      data.table(chr_pos = table3layers$chr_pos[queryHits(hits)], gene = gr$name)
    })))
  unique_cpgs <- unique(na.omit(cpg_extract_list$chr_pos))
  message(sprintf("Total unique CpGs to extract: %d", length(unique_cpgs)))
  print(cpg_extract_list[, .N, by = gene][order(-N)])
  out_path <- here("B_MultiTissues/dataOut",
                   sprintf("targetCpGs_Alice_%s%s.txt", layer, format(Sys.time(), "%Y%m%d")))
  writeLines(unique_cpgs, out_path); message("Written to: ", out_path); out_path
}
matt_cpg_file <- getCpG4pchuckle(candidates_gr_Matt, "Matt")
# Written to: /home/alice/Documents/Research/GIT/2024_hvCpG/B_MultiTissues/dataOut/targetCpGs_Alice_Matt20260709.txt

## Then on pchuckle (adapt the dated filename printed above):
# source /share/apps/source_files/python/python-3.13.0a6.source
# cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# python3 S00_extractRawMethylationForTargetCpG.py \
# --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/targetCpGs_Alice_Matt20260709.txt \
# --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/methylation_targetCpGs_Alice_Matt20260709.tsv \
# --minCov    10

# ══════════════════════════════════════════════════════════════════════════════
# 2. Plotting driver  (identical logic to S11's getFullPlotsTarget)
# ══════════════════════════════════════════════════════════════════════════════
getFullPlotsTarget <- function(meth, candidates_gr, all_candidates,
                               layer, gene_coords_dt = NULL, flank = 5000,
                               features_gr = NULL, style = "full",
                               trunc_kb = NULL, fetal_dt = NULL,
                               show_raw_meth = TRUE, dmr_gr = NULL){
  
  setDT(meth); setDT(all_candidates)
  
  if (is.null(gene_coords_dt)) {
    stopifnot(all(c("gene","chr","gene_start","gene_end") %in% names(all_candidates)))
    gene_coords_dt <- unique(all_candidates[, .(
      gene_name = gene, seqnames = chr, start = gene_start, end = gene_end,
      strand = "*", width = gene_end - gene_start + 1)])
    setorder(gene_coords_dt, gene_name, -width)
    gene_coords_dt <- unique(gene_coords_dt, by = "gene_name")
    message(sprintf("Built gene_coords_dt for %d regions", nrow(gene_coords_dt)))
  }
  
  loyfer_meta <- fread(here("B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv"))
  loyfer_meta[, source_tissue_celltype := paste0(`Source Tissue`, " - ", `Cell type`)]
  tissue_to_layer <- unique(loyfer_meta[, .(source_tissue_celltype, germ_layer = `Germ layer`)])
  meth <- merge(meth, tissue_to_layer, by = "source_tissue_celltype", all.x = TRUE)
  message(sprintf("Rows with missing germ_layer: %d", sum(is.na(meth$germ_layer))))
  if (!"pos" %in% names(meth)) meth[, pos := as.integer(sub(".*_", "", cpg_site))]
  if (!"chr" %in% names(meth)) meth[, chr := sub("_.*", "", cpg_site)]
  
  compute_interlayer_r_loyfer <- function(meth, min_patients = 3) {
    agg  <- meth[, .(methylation = mean(methylation, na.rm = TRUE)),
                 by = .(cpg_site, patient_id, germ_layer)]
    wide <- dcast(agg, cpg_site + patient_id ~ germ_layer, value.var = "methylation")
    for (col in c("Endo","Meso","Ecto")) if (!col %in% names(wide)) wide[, (col) := NA_real_]
    wide[, {
      idx_EM <- !is.na(Endo) & !is.na(Meso); idx_EE <- !is.na(Endo) & !is.na(Ecto); idx_ME <- !is.na(Meso) & !is.na(Ecto)
      n_EM <- sum(idx_EM); n_EE <- sum(idx_EE); n_ME <- sum(idx_ME)
      r_EM <- if (n_EM >= min_patients) suppressWarnings(cor(Endo[idx_EM], Meso[idx_EM])) else NA_real_
      r_EE <- if (n_EE >= min_patients) suppressWarnings(cor(Endo[idx_EE], Ecto[idx_EE])) else NA_real_
      r_ME <- if (n_ME >= min_patients) suppressWarnings(cor(Meso[idx_ME], Ecto[idx_ME])) else NA_real_
      
      list(n_Endo_Meso = n_EM, n_Endo_Ecto = n_EE, n_Meso_Ecto = n_ME,
           r_Endo_Meso = r_EM, r_Endo_Ecto = r_EE, r_Meso_Ecto = r_ME,
           interlayer_r = abs(mean(c(r_EM, r_EE, r_ME), na.rm = TRUE)), ## match fetal data but not all have all 3 layers
           low_power    = n_EM < 8)
      
    }, by = cpg_site]
  }
  interlayer_loyfer <- compute_interlayer_r_loyfer(meth, min_patients = 3)
  message(sprintf("Inter-layer r computed for %d CpGs", nrow(interlayer_loyfer)))
  
  compute_intralayer_r_loyfer <- function(meth, layer_name, min_obs = 3) {
    layer_meth <- meth[germ_layer == layer_name]
    eligible   <- layer_meth[, .(n = uniqueN(source_tissue_celltype)), by = patient_id][n > 1, patient_id]
    if (length(eligible) == 0) { message(sprintf("  No patients with >1 tissue in %s", layer_name)); return(NULL) }
    message(sprintf("  %s: %d eligible patients with >1 tissue", layer_name, length(eligible)))
    all_pairs <- rbindlist(lapply(eligible, function(pid) {
      pat <- layer_meth[patient_id == pid]; tiss <- unique(pat$source_tissue_celltype)
      if (length(tiss) < 2) return(NULL)
      rbindlist(lapply(combn(tiss, 2, simplify = FALSE), function(pr) {
        d1 <- pat[source_tissue_celltype == pr[1], .(cpg_site, m1 = methylation)]
        d2 <- pat[source_tissue_celltype == pr[2], .(cpg_site, m2 = methylation)]
        merge(d1, d2, by = "cpg_site")[, `:=`(patient_id = pid, tissue1 = pr[1], tissue2 = pr[2])]
      }), fill = TRUE)
    }), fill = TRUE)
    if (is.null(all_pairs) || nrow(all_pairs) == 0) return(NULL)
    all_pairs[!is.na(m1) & !is.na(m2), {
      if (.N >= min_obs) list(r_intralayer = suppressWarnings(cor(m1, m2, method = "pearson")), n_obs = .N, layer = layer_name)
      else               list(r_intralayer = NA_real_, n_obs = .N, layer = layer_name)
    }, by = cpg_site]
  }
  intralayer_loyfer <- rbindlist(lapply(c("Endo","Meso","Ecto"), function(lay) {
    message(sprintf("Computing intra-layer r for: %s", lay))
    compute_intralayer_r_loyfer(meth, lay, min_obs = 3)
  }), fill = TRUE)
  intralayer_wide <- dcast(intralayer_loyfer[!is.na(r_intralayer)], cpg_site ~ layer, value.var = "r_intralayer")
  setnames(intralayer_wide,
           intersect(c("Endo","Meso","Ecto"), names(intralayer_wide)),
           paste0("r_intra_", intersect(c("Endo","Meso","Ecto"), names(intralayer_wide))),
           skip_absent = TRUE)
  
  seqlevelsStyle(candidates_gr) <- "UCSC"; candidates_gr <- keepStandardChromosomes(candidates_gr, pruning.mode = "coarse")
  seqlevelsStyle(table3layers)  <- "UCSC"; table3layers  <- keepStandardChromosomes(table3layers,  pruning.mode = "coarse")
  
  candidates_cpgs <- rbindlist(lapply(seq_along(candidates_gr), function(i) {
    gr <- candidates_gr[i]; hits <- findOverlaps(table3layers, gr)
    if (length(hits) == 0) return(NULL)
    dt <- as.data.table(table3layers[queryHits(hits)])[, .(chr_pos, alpha_endo, alpha_meso, alpha_ecto, alpha_geomean)]
    dt[, gene := candidates_gr$name[i]]; dt
  }), fill = TRUE)
  candidates_cpgs_r <- merge(candidates_cpgs,
                             interlayer_loyfer[, .(cpg_site, r_Endo_Meso, interlayer_r, n_Endo_Meso, low_power)],
                             by.x = "chr_pos", by.y = "cpg_site", all.x = TRUE)
  candidates_cpgs_full <- merge(candidates_cpgs_r, intralayer_wide,
                                by.x = "chr_pos", by.y = "cpg_site", all.x = TRUE)
  for (col in c("r_intra_Endo","r_intra_Meso","r_intra_Ecto"))
    if (!col %in% names(candidates_cpgs_full)) candidates_cpgs_full[, (col) := NA_real_]
  message(sprintf("candidates_cpgs_full: %d rows, %d CpGs, %d regions",
                  nrow(candidates_cpgs_full), uniqueN(candidates_cpgs_full$chr_pos),
                  uniqueN(candidates_cpgs_full$gene)))
  
  dir.create(here("B_MultiTissues/dataOut/figures/lociAlice"), showWarnings = FALSE, recursive = TRUE)
  all_windows <- GRanges(seqnames = gene_coords_dt$seqnames,
                         ranges = IRanges(gene_coords_dt$start - flank, gene_coords_dt$end + flank))
  table3layers_sub <- suppressWarnings(subsetByOverlaps(table3layers, all_windows))
  message(sprintf("table3layers subset: %d CpGs", length(table3layers_sub)))
  
  meth[, pos := as.integer(sub(".*_", "", cpg_site))]; meth[, chr := sub("_.*", "", cpg_site)]
  message("Pre-splitting meth by region...")
  meth_by_gene <- lapply(unique(candidates_cpgs_full$gene), function(gn) {
    if (!gn %in% gene_coords_dt$gene_name) return(NULL)
    gc <- gene_coords_dt[gene_name == gn]
    meth[chr == as.character(gc$seqnames) & pos >= gc$start - flank & pos <= gc$end + flank]
  })
  names(meth_by_gene) <- unique(candidates_cpgs_full$gene)
  
  for (gn in unique(candidates_cpgs_full$gene)) {
    if (!gn %in% gene_coords_dt$gene_name) { message(sprintf("Skipping %s: no coords", gn)); next }
    cat_label <- layer[1]; if (is.na(cat_label)) cat_label <- "unknown"
    out_file  <- here(sprintf("B_MultiTissues/dataOut/figures/script13/locus_%s_%s.png", cat_label, gn))
    if (file.exists(out_file)) { message(sprintf("Skipping %s: exists", gn)); next }
    meth_gn <- meth_by_gene[[gn]]
    if (is.null(meth_gn) || nrow(meth_gn) == 0) { message(sprintf("Skipping %s: no meth", gn)); next }
    tryCatch({
      if (style == "acd") {
        p <- plot_locus_acd_annotated(
          gene_name_arg = gn, candidates_cpgs_full = candidates_cpgs_full,
          interlayer_loyfer = interlayer_loyfer, meth = meth_gn,
          gene_coords_dt = gene_coords_dt, table3layers = table3layers_sub,
          features_gr = features_gr, flank = flank,
          trunc_kb = trunc_kb, fetal_dt = fetal_dt,
          show_raw_meth = show_raw_meth, , dmr_gr = dmr_gr)
        out_file <- here(sprintf("B_MultiTissues/dataOut/figures/script13/acd_%s_%s.png", cat_label, gn))
        ggplot2::ggsave(out_file, p, width = 15, height = 13, dpi = 300, bg = "white")
      } else {
        p <- plot_candidate_locus(
          gene_name_arg = gn, candidates_cpgs_full = candidates_cpgs_full,
          interlayer_loyfer = interlayer_loyfer, meth = meth_gn,
          gene_coords_dt = gene_coords_dt, table3layers = table3layers_sub, flank = flank)
        ggplot2::ggsave(out_file, p, width = 15, height = 13, dpi = 300, bg = "white")
      }
      message(sprintf("Saved: %s", basename(out_file))); rm(p); gc()
    }, error = function(e) message(sprintf("Failed for %s: %s", gn, e$message)))
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 2b. Trimmed locus figure: panels a (Pr(HV)) + c (correlation) + d (LD),
#     with a gene-model (promoter/exon/intron) + feature/LTR annotation track.
# ══════════════════════════════════════════════════════════════════════════════

## ── Gene model with strand-aware truncation n kb after the TSS ───────────────
## trunc_kb: NULL (no truncation), a single numeric (applies to all genes),
##           or a named numeric vector, e.g. c(ACTL8 = 20, LY6S = 5),
##           optionally with a "default" entry: c(default = 10, ACTL8 = 30).
get_gene_model <- function(chr, win_start, win_end,
                           promoter_up = 2000,
                           trunc_kb    = NULL) {
  edb <- EnsDb.Hsapiens.v86
  w   <- GRanges(sub("^chr", "", chr), IRanges(win_start, win_end))
  ex  <- tryCatch(exonsBy(edb, by = "gene", filter = GRangesFilter(w, type = "any")),
                  error = function(e) NULL)
  if (is.null(ex) || length(ex) == 0) return(NULL)
  gsym <- mapIds(edb, names(ex), "SYMBOL", "GENEID")
  
  trunc_for <- function(sym) {
    if (is.null(trunc_kb)) return(NA_real_)
    if (is.null(names(trunc_kb))) return(as.numeric(trunc_kb[1]))
    if (sym %in% names(trunc_kb)) return(as.numeric(trunc_kb[[sym]]))
    if ("default" %in% names(trunc_kb)) return(as.numeric(trunc_kb[["default"]]))
    NA_real_
  }
  
  segs <- rbindlist(lapply(seq_along(ex), function(i) {
    e   <- reduce(ex[[i]]); e <- e[order(start(e))]
    gid <- names(ex)[i]; sym <- gsym[gid]; if (is.na(sym)) sym <- gid
    g_start <- min(start(e)); g_end <- max(end(e)); str <- as.character(strand(e))[1]
    tss <- if (str == "-") g_end else g_start
    
    n_kb <- trunc_for(sym)
    if (is.finite(n_kb)) {
      cut_start <- if (str == "-") tss - n_kb * 1000 else tss
      cut_end   <- if (str == "-") tss               else tss + n_kb * 1000
    } else { cut_start <- -Inf; cut_end <- Inf }
    
    exdt <- data.table(gene = sym, type = "exon",
                       start = start(e), end = end(e), strand = str)
    indt <- if (length(e) > 1)
      data.table(gene = sym, type = "intron",
                 start = end(e)[-length(e)] + 1, end = start(e)[-1] - 1, strand = str) else NULL
    prdt <- data.table(gene = sym, type = "promoter",
                       start = if (str == "-") tss + 1 else tss - promoter_up,
                       end   = if (str == "-") tss + promoter_up else tss - 1,
                       strand = str)
    
    body <- rbindlist(list(exdt, indt), fill = TRUE)
    if (!is.null(body) && nrow(body)) {                     # clip body to truncation
      body[start < cut_start, start := cut_start]
      body[end   > cut_end,   end   := cut_end]
      body <- body[end >= start]
    }
    out <- rbindlist(list(body, prdt), fill = TRUE)          # promoter never truncated
    out[, `:=`(tss = tss,
               cut_at = if (is.finite(n_kb)) (if (str == "-") cut_start else cut_end) else NA_real_,
               truncated = is.finite(n_kb))]
    out
  }), fill = TRUE)
  
  segs[, chr := chr]
  segs[start < win_start, start := win_start]
  segs[end   > win_end,   end   := win_end]
  segs <- segs[end >= start]
  segs[, type := factor(type, levels = c("promoter", "exon", "intron"))]
  segs[]
}

## ── Gene-model track with a fill legend ─────────────────────────────────────
gene_model_track <- function(gm, shared_x, shared_coord, base_th, hide_x,
                             show_cut = TRUE, show_x = FALSE) {
  seg_cols <- c(promoter = "#4575B4", exon = "grey20", intron = "grey55")
  if (is.null(gm) || !nrow(gm)) {
    return(ggplot() + shared_x + shared_coord + base_th + (if (show_x) NULL else hide_x) +
             scale_y_continuous("genes", breaks = NULL, limits = c(0, 1)))
  }
  gm[, y := as.integer(factor(gene))]
  hh <- c(promoter = 0.30, exon = 0.25, intron = 0.04)      # half-heights
  gm[, h := hh[as.character(type)]]
  
  p <- ggplot(gm) +
    geom_rect(aes(xmin = start, xmax = end, ymin = y - h, ymax = y + h, fill = type),
              alpha = 0.95) +
    geom_text(data = gm[, .(x = mean(c(min(start), max(end)))), by = .(gene, y)],
              aes(x = x, y = y + 0.45, label = gene), size = 2.6) +
    scale_fill_manual("Gene model", values = seg_cols, drop = FALSE) +
    shared_x + shared_coord +
    scale_y_continuous("genes", breaks = NULL, limits = c(0.4, max(gm$y) + 0.7)) +
    base_th + (if (show_x) NULL else hide_x) +
    theme(axis.ticks.y = element_blank(), legend.position = "right",
          legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 7),
          legend.title = element_text(size = 8))
  
  cuts <- unique(gm[truncated == TRUE, .(gene, y, cut_at)])
  if (show_cut && nrow(cuts))
    p <- p + geom_segment(data = cuts,
                          aes(x = cut_at, xend = cut_at, y = y - 0.35, yend = y + 0.35),
                          linetype = 2, linewidth = 0.35, colour = "grey20")
  p
}

plot_locus_acd_annotated <- function(gene_name_arg, candidates_cpgs_full,
                                     interlayer_loyfer, meth, gene_coords_dt,
                                     table3layers, features_gr = NULL,
                                     flank = 5000, min_samples = 5,
                                     peak_span = 0.3, max_ld_cpg = 600L,
                                     trunc_kb = NULL, fetal_dt = NULL, 
                                     show_raw_meth = TRUE, dmr_gr = NULL) {
  
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("plot_locus_acd_annotated() needs 'patchwork'.")
  
  layer_pal <- c(Endo = "#009E73", Meso = "#56B4E9", Ecto = "#CC79A7")
  corr_pal  <- c("Inter-layer |r|" = "grey55", "Intra Endo" = "#009E73",
                 "Intra Meso" = "#56B4E9", "Intra Ecto" = "#CC79A7")
  me_yellow <- "#E6AB02"
  feat_col  <- "#D55E00"
  
  gc <- gene_coords_dt[gene_name == gene_name_arg]
  if (nrow(gc) == 0) stop(sprintf("Gene '%s' not found", gene_name_arg))
  gene_start <- gc$start; gene_end <- gc$end; chr <- as.character(gc$seqnames)
  x_min <- gene_start - flank; x_max <- gene_end + flank
  
  kb <- function(x) paste0(round(x / 1e3, 1), " kb")
  shared_x     <- scale_x_continuous("Position (hg38)", labels = kb, expand = c(0.01, 0))
  shared_coord <- coord_cartesian(xlim = c(x_min, x_max))
  base_th <- theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 7.5))
  hide_x  <- theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  gene_band <- annotate("rect", xmin = gene_start, xmax = gene_end,
                        ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.6)
  
  # ── feature overlay (vertical bands) for a & c ──────────────────────────────
  feat_dt <- NULL
  if (!is.null(features_gr) && length(features_gr)) {
    fd <- as.data.table(features_gr)[as.character(seqnames) == chr]
    fd <- fd[end >= x_min & start <= x_max]
    if (nrow(fd)) feat_dt <- fd
  }
  feature_bands <- if (!is.null(feat_dt))
    geom_rect(data = feat_dt,
              aes(xmin = pmax(start, x_min), xmax = pmin(end, x_max),
                  ymin = -Inf, ymax = Inf),
              fill = feat_col, alpha = 0.10, inherit.aes = FALSE) else NULL
  
  # DMR CpGs in this window
  dmr_dt_w <- NULL
  if (!is.null(dmr_gr) && length(dmr_gr)) {
    dd <- as.data.table(dmr_gr)[as.character(seqnames) == chr]
    dd <- dd[end >= x_min & start <= x_max]
    if (nrow(dd)) dmr_dt_w <- dd
  }
  dmr_lines <- if (!is.null(dmr_dt_w))
    geom_vline(data = dmr_dt_w, aes(xintercept = start),
               colour = "#B2182B", linetype = 3, linewidth = 0.4,
               inherit.aes = FALSE) else NULL
  
  # feat_dt used by the LD annotation band: sub-gene features only
  # (genes are already shown by the gene-model track, and would mask exon/intron)
  gene_level <- c("ACTL8", "LY6S_AS1")
  feat_dt_ld <- if (!is.null(feat_dt)) feat_dt[!name %in% gene_level] else NULL
  if (!is.null(dmr_dt_w)) {
    dmr_feat   <- dmr_dt_w[, .(seqnames, start, end, name = "DMR")]
    feat_dt_ld <- rbindlist(list(feat_dt_ld, dmr_feat), fill = TRUE)
  }
  if (!is.null(feat_dt_ld) && !nrow(feat_dt_ld)) feat_dt_ld <- NULL
  
  # ── panel a ─────────────────────────────────────────────────────────────────
  cpg_sub <- copy(as.data.table(candidates_cpgs_full)[gene == gene_name_arg])
  cpg_sub[, pos := as.integer(sub(".*_", "", chr_pos))]
  if ("alpha_geomean" %in% names(cpg_sub)) cpg_sub[is.na(alpha_geomean), alpha_geomean := 0]
  cpg_sub[, is_me := alpha_geomean >= 0.5]
  
  pa <- ggplot(cpg_sub, aes(pos, alpha_geomean)) +
    gene_band + feature_bands + dmr_lines +
    geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3) +
    geom_smooth(method = "loess", span = peak_span, se = FALSE, linewidth = 0.7, alpha = 0.15) +
    geom_point(aes(colour = is_me), size = 1, alpha = 0.85) +
    scale_colour_manual(values = c(`TRUE` = me_yellow, `FALSE` = "grey60"), guide = "none") +
    shared_x + shared_coord +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
    base_th + hide_x +
    labs(y = "Pr(HV)\ngeomean",
         title = sprintf("%s  (%s:%s-%s, flank \u00b1%gkb)", gene_name_arg, chr,
                         format(gene_start, big.mark = ","),
                         format(gene_end, big.mark = ","), flank / 1e3),
         subtitle = "a. ME test: yellow = alpha_geomean \u2265 0.5") +
    theme(plot.title = element_text(face = "bold", size = 11))
  
  # ── panel c ─────────────────────────────────────────────────────────────────
  il <- copy(as.data.table(interlayer_loyfer))
  if ("cpg_site" %in% names(il)) il[, pos := as.integer(sub(".*_", "", cpg_site))]
  
  inter_dt <- if (all(c("pos","interlayer_r") %in% names(il)))
    il[pos >= x_min & pos <= x_max, .(pos, r = interlayer_r, grp = "Inter-layer |r|")] else NULL
  
  intra_parts <- list()
  
  if ("r_intra_Endo" %in% names(cpg_sub)) intra_parts[["Intra Endo"]] <- cpg_sub[, .(pos, r = abs(r_intra_Endo), grp = "Intra Endo")]
  if ("r_intra_Meso" %in% names(cpg_sub)) intra_parts[["Intra Meso"]] <- cpg_sub[, .(pos, r = abs(r_intra_Meso), grp = "Intra Meso")]
  if ("r_intra_Ecto" %in% names(cpg_sub)) intra_parts[["Intra Ecto"]] <- cpg_sub[, .(pos, r = abs(r_intra_Ecto), grp = "Intra Ecto")]
  
  corr_dt <- rbindlist(c(list(inter_dt), intra_parts), use.names = TRUE, fill = TRUE)
  corr_dt <- corr_dt[!is.na(r)]
  pc <- NULL
  if (nrow(corr_dt)) {
    corr_dt[, grp := factor(grp, levels = names(corr_pal))]
    pc <- ggplot(corr_dt, aes(pos, r, colour = grp, fill = grp)) +
      gene_band + feature_bands + dmr_lines +
      geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.25, colour = "grey40") +
      geom_point(size = 0.5, alpha = 0.30) +
      geom_smooth(method = "loess", span = peak_span, se = TRUE, linewidth = 0.7, alpha = 0.15) +
      scale_colour_manual("Correlation", values = corr_pal, drop = FALSE) +
      scale_fill_manual("Correlation", values = corr_pal, drop = FALSE) +
      shared_x + shared_coord + scale_y_continuous(limits = c(0, 1)) +
      base_th + hide_x +
      labs(y = "|Pearson r|",
           subtitle = "c. Grey = inter-layer |r| (only endo-meso, 4 people) | Coloured = intra-layer |r| (within individuals)")
  }
  
  # ── gene-model track (promoter / exon / intron) ─────────────────────────────
  gm     <- tryCatch(get_gene_model(chr, x_min, x_max, trunc_kb = trunc_kb),
                     error = function(e) NULL)
  
  p_gene <- gene_model_track(gm, shared_x, shared_coord, base_th, hide_x,
                             show_x = TRUE)
  
  # ── feature / LTR track (bottom, shows x-axis) ──────────────────────────────
  p_feat <- ggplot() + shared_x + shared_coord +
    scale_y_continuous(NULL, breaks = NULL, limits = c(0, 1)) +
    base_th + hide_x                                  # now hides its x-axis
  
  if (!is.null(feat_dt)) {
    p_feat <- p_feat +
      geom_rect(data = feat_dt, aes(xmin = pmax(start, x_min), xmax = pmin(end, x_max),
                                    ymin = 0.30, ymax = 0.70, fill = name), alpha = 0.75) +
      geom_text(data = feat_dt, aes(x = (pmax(start, x_min) + pmin(end, x_max)) / 2,
                                    y = 0.88, label = name), size = 2.4) +
      scale_fill_brewer(palette = "Set2", guide = "none")
  }
  
  if (!is.null(dmr_dt_w)) {
    p_feat <- p_feat +
      geom_point(data = dmr_dt_w, aes(x = start, y = 0.12),
                 shape = 25, size = 1.8, fill = "#B2182B", colour = "#B2182B",
                 inherit.aes = FALSE) +
      annotate("text", x = x_min, y = 0.12, label = "DMR CpGs",
               hjust = -0.05, size = 2.2, colour = "#B2182B")
  }
  
  # ── panel d: index-based LD triangle + density inset + feature index lines ──
  make_ld <- function(md_sub, label, ms) {
    if (is.null(md_sub) || !nrow(md_sub)) return(NULL)
    keep <- md_sub[, .(n = uniqueN(sample_id)), by = .(cpg_site, pos)][n >= ms]
    md_sub <- md_sub[cpg_site %in% keep$cpg_site]
    lev <- sort(unique(md_sub$pos)); if (length(lev) < 3) return(NULL)
    if (length(lev) > max_ld_cpg) lev <- lev[round(seq(1, length(lev), length.out = max_ld_cpg))]
    md_sub <- md_sub[pos %in% lev]
    w <- dcast(md_sub, sample_id ~ pos, value.var = "methylation")
    mat <- as.matrix(w[, -1, with = FALSE]); mat <- mat[, order(as.integer(colnames(mat))), drop = FALSE]
    if (ncol(mat) < 3) return(NULL)
    cpos <- as.integer(colnames(mat))
    cormat <- suppressWarnings(stats::cor(mat, use = "pairwise.complete.obs")); n <- ncol(cormat)
    
    # section labels along the CpG-index axis  (moved here: needs `n`)
    ann  <- annotate_positions(cpos, feat_dt_ld, gm)
    runs <- label_runs(ann)
    
    generic <- c("intergenic", "promoter", "exon", "intron")
    runs <- runs[n >= 2 | !(label %in% generic)]
    
    band_h <- max(4, round(0.03 * n))
    
    long <- as.data.table(as.table(cormat)); setnames(long, c("c1","c2","r"))
    idx <- setNames(seq_len(n), colnames(cormat))
    long[, i := idx[as.character(c1)]]; long[, j := idx[as.character(c2)]]
    long <- long[j > i & !is.na(r)]; if (!nrow(long)) return(NULL)
    br <- unique(round(seq(1, n, length.out = 6)))
    fl <- NULL
    if (!is.null(feat_dt_ld)) {  
      bnd <- sort(unique(c(feat_dt_ld$start, feat_dt_ld$end)))
      bnd <- bnd[bnd >= min(cpos) & bnd <= max(cpos)]
      if (length(bnd)) fl <- sapply(bnd, function(b) which.min(abs(cpos - b)))
    }
    
    g <- ggplot(long, aes(i, j, fill = r)) +
      geom_raster() +
      # coloured annotation band below the triangle
      geom_rect(data = runs, inherit.aes = FALSE,
                aes(xmin = i_start - 0.5, xmax = i_end + 0.5,
                    ymin = -band_h, ymax = -1),
                fill = "grey85", colour = "white", linewidth = 0.2) +
      scale_fill_gradient2("Pearson r", low = "#2166AC", mid = "white",
                           high = "#D73027", midpoint = 0, limits = c(-1, 1)) +
      scale_x_continuous(breaks = runs$i_mid, labels = runs$label, expand = c(0, 0)) +
      scale_y_continuous(breaks = br, labels = kb(cpos[br]),
                         limits = c(-band_h, n + 1), expand = c(0, 0)) +
      coord_fixed() + base_th +
      theme(axis.title = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
      labs(subtitle = label)
    if (!is.null(fl))
      g <- g + geom_vline(xintercept = fl, colour = feat_col, linewidth = 0.3, alpha = 0.6) +
      geom_hline(yintercept = fl, colour = feat_col, linewidth = 0.3, alpha = 0.6)
    g
    
  }
  
  meth_dt <- copy(as.data.table(meth)); mdw <- data.table()
  if (nrow(meth_dt) > 0 && all(c("methylation","cpg_site","patient_id") %in% names(meth_dt))) {
    if (!"pos" %in% names(meth_dt)) meth_dt[, pos := as.integer(sub(".*_", "", cpg_site))]
    mdw <- meth_dt[pos >= x_min & pos <= x_max & !is.na(methylation)]
    if (nrow(mdw)) {
      mdw[, sample_id := if ("source_tissue_celltype" %in% names(mdw))
        paste(patient_id, source_tissue_celltype, sep = "|") else as.character(patient_id)]
      mdw <- mdw[, .(methylation = mean(methylation, na.rm = TRUE)), by = .(sample_id, cpg_site, pos)]
    }
  }
  ld_all <- if (nrow(mdw)) make_ld(mdw, sprintf("d. Co-methylation pairwise r \u2014 %s \u00b1%gkb (all samples)",
                                                gene_name_arg, flank / 1e3), min_samples) else NULL
  pd_dens <- NULL
  if (nrow(mdw)) {
    pd_dens <- ggplot(data.table(pos = sort(unique(mdw$pos))), aes(pos)) +
      geom_histogram(bins = 60, fill = "grey40") +
      scale_x_continuous(labels = kb, limits = c(x_min, x_max), expand = c(0, 0)) +
      base_th + theme(axis.title = element_text(size = 6), axis.text = element_text(size = 5),
                      plot.background = element_rect(fill = "white", colour = NA)) +
      labs(x = NULL, y = "CpG\ndensity")
  }
  pd <- if (!is.null(ld_all) && !is.null(pd_dens))
    ld_all + patchwork::inset_element(pd_dens, left = 0.50, bottom = 0.02, right = 1.0, top = 0.32, align_to = "panel")
  else ld_all
  
  # fetal correlation track (panel b)
  pf <- if (!is.null(fetal_dt))
    fetal_corr_track(fetal_dt, chr, x_min, x_max, shared_x, shared_coord, base_th, hide_x,
                     peak_span) else NULL
  
  # raw methylation track (panel e)
  pr <- if (isTRUE(show_raw_meth))
    raw_meth_track(meth, chr, x_min, x_max, shared_x, shared_coord, base_th, hide_x,
                   peak_span) else NULL
  
  # ── assemble: a | b(fetal) | c | e(raw meth) | genes | features   ||   d ────
  left_panels <- Filter(Negate(is.null), list(pa, pf, pc, pr, p_feat, p_gene))
  hts <- c(3.2, 3.0, 3.2, 3.0, 1.1, 1.4)[seq_along(left_panels)]
  left <- patchwork::wrap_plots(left_panels, ncol = 1, heights = hts)
  if (!is.null(pd)) patchwork::wrap_plots(left, pd, ncol = 2, widths = c(1.15, 1))
  else left
  
}

# ══════════════════════════════════════════════════════════════════════════════
# 2c. Extra tracks: fetal inter-layer r, raw methylation, LD annotation band
# ══════════════════════════════════════════════════════════════════════════════

## Fetal inter-germ-layer correlation (from testFetalSIV_ingp5.R)
##   interlayer_corr_all: CpG (EPIC probe), interlayer_r, percentile_r
##   `dico` (from prepPreviousSIV.R) maps CpG -> chrpos_hg38
get_fetal_corr_dt <- function(path = here("B_MultiTissues/dataOut/interlayer_corr_all.RDS")) {
  fc <- as.data.table(readRDS(path))
  fc[, chrpos_hg38 := dico$chrpos_hg38[match(CpG, dico$CpG)]]
  fc <- fc[!is.na(chrpos_hg38) & !is.na(interlayer_r)]
  fc[, `:=`(chr = sub("_.*", "", chrpos_hg38),
            pos = as.integer(sub(".*_", "", chrpos_hg38)))]
  fc[]
}

## Track 1: fetal intra-individual mean inter-germ-layer correlation
fetal_corr_track <- function(fetal_dt, chr, x_min, x_max,
                             shared_x, shared_coord, base_th, hide_x,
                             peak_span = 0.3) {
  chr_w <- chr                                   # avoid masking by the column
  d <- fetal_dt[chr == chr_w & pos >= x_min & pos <= x_max]
  p <- ggplot(d, aes(pos, interlayer_r)) +
    shared_x + shared_coord +
    scale_y_continuous(limits = c(0, 1)) +
    base_th + hide_x +
    labs(y = "Fetal inter-layer r\n(mean |Pearson|)",
         subtitle = "b. Fetal EPIC: intra-individual mean inter-germ-layer correlation")
  if (nrow(d) >= 1) {
    p <- p + geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "grey40") +
      geom_point(size = 1.1, colour = "#6A3D9A", alpha = 0.8)
  } else {
    p <- p + annotate("text", x = (x_min + x_max)/2, y = 0.5,
                      label = "no fetal EPIC probes in window", size = 2.8, colour = "grey40")
  }
  p
}

## Track 2: raw methylation along the chromosome, coloured by germ layer
raw_meth_track <- function(meth, chr, x_min, x_max,
                           shared_x, shared_coord, base_th, hide_x,
                           peak_span = 0.3) {
  layer_pal <- c(Endo = "#009E73", Meso = "#56B4E9", Ecto = "#CC79A7")
  md <- as.data.table(meth)
  if (!"pos" %in% names(md)) md[, pos := as.integer(sub(".*_", "", cpg_site))]
  md <- md[pos >= x_min & pos <= x_max & !is.na(methylation) &
             germ_layer %in% names(layer_pal)]
  p <- ggplot(md, aes(pos, methylation, colour = germ_layer, fill = germ_layer)) +
    shared_x + shared_coord +
    scale_colour_manual("Germ layer", values = layer_pal, drop = FALSE) +
    scale_fill_manual("Germ layer",   values = layer_pal, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1)) +
    base_th + hide_x +
    labs(y = "Raw methylation\n(beta)",
         subtitle = "e. Raw methylation per sample along the locus")
  if (nrow(md)) {
    p <- p + geom_point(size = 0.35, alpha = 0.20)
    if (nrow(md) >= 20)
      p <- p + geom_smooth(method = "loess", span = peak_span, se = FALSE, linewidth = 0.7)
  }
  p
}

## Per-CpG annotation label: named feature > gene-model segment > intergenic
annotate_positions <- function(cpos, feat_dt, gm) {
  
  lab <- rep("intergenic", length(cpos))
  if (!is.null(gm) && nrow(gm)) {
    for (ty in c("intron", "exon", "promoter")) {          # promoter wins over exon/intron
      sub <- gm[type == ty]
      if (nrow(sub))
        for (k in seq_len(nrow(sub)))
          lab[cpos >= sub$start[k] & cpos <= sub$end[k]] <- ty
    }
  }
  if (!is.null(feat_dt) && nrow(feat_dt))                   # named features win over everything
    for (k in seq_len(nrow(feat_dt)))
      lab[cpos >= feat_dt$start[k] & cpos <= feat_dt$end[k]] <- as.character(feat_dt$name[k])
  lab
}

## Contiguous runs of identical labels -> rects + midpoints on the index axis
label_runs <- function(lab) {
  r <- rle(lab)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  data.table(label = r$values, i_start = starts, i_end = ends,
             i_mid = (starts + ends) / 2, n = r$lengths)
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Run (after the python extraction produced the TSV)
# ══════════════════════════════════════════════════════════════════════════════
matt_meth_tsv <- here("gitignore", "methylation_targetCpGs_Alice_Matt20260709.tsv")

matt_features_gr <- c(LY6SVMR_hg38, LY6S_AS1_hg38, MER11C_hg38,
                      LTR41_hg38, ACTL8_hg38)
seqlevelsStyle(matt_features_gr) <- "UCSC"

fetal_dt <- get_fetal_corr_dt()      # once, before the calls

getFullPlotsTarget(
  meth = fread(matt_meth_tsv),
  candidates_gr  = candidates_gr_Matt[candidates_gr_Matt$name %in%
                                        c("Region1_LY6SVMR_MER11C","Region2_ACTL8_LTR41")],
  all_candidates = all_candidates_Matt[gene %in%
                                         c("Region1_LY6SVMR_MER11C","Region2_ACTL8_LTR41")],
  layer = "Matt", features_gr = matt_features_gr, style = "acd",
  trunc_kb = NULL, flank = 0, fetal_dt = fetal_dt, dmr_gr = dmr_gr)

message("S13 complete.")

