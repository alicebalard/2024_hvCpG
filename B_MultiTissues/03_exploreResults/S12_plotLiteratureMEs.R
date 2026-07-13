#####################################################################
## S12_plotLiteratureMEs.R
## Produce per-locus plots (same panels as S11) for a hand-picked
## list of literature metastable-epiallele (ME) genes.
##
## Self-contained: does NOT depend on S10/S11 objects. Needs only
##   - functions.R (plot_candidate_locus)
##   - table3layers, wideFull  (as in S11)
##   - a methylation TSV extracted for these genes' CpGs
##
## Produces: B_MultiTissues/dataOut/figures/lociAlice/locus_MElit_<GENE>.png
#####################################################################

library(here)
source(here("B_MultiTissues", "quiet_library.R"))
if (!exists("functionsLoaded")) source(here("B_MultiTissues/03_exploreResults", "functions.R"))
if (!exists("table3layers"))    load(here("gitignore/fullTable3layers.Rda"))
if (!exists("wideFull"))        wideFull <- readRDS(here("gitignore/wide_script10_3layers_full.RDS"))

if (is.null(table3layers$percentile)) {
  table3layers$percentile <- ecdf(table3layers$alpha_geomean)(table3layers$alpha_geomean) * 100
}

library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(GenomeInfoDb)

# ══════════════════════════════════════════════════════════════════════════════
# 0. The literature ME gene list
# ══════════════════════════════════════════════════════════════════════════════
lit_MEs <- c("VTRNA2-1","PM20D1","POMC","CYP2E1","PAX8","SPATC1L",
             "HLA-DQA1","HLA-A","HLA-DPA1","ZXDA","GPRASP1","MMGT1",
             "IDS","CNKSR2","IRS4","ZMAT1")

# some canonical MEs carry alias symbols in EnsDb — map the awkward ones
alias_map <- c("VTRNA2-1" = "VTRNA2-1")   # add fallbacks here if lookup misses, e.g. "nc886","MIR886"

# ══════════════════════════════════════════════════════════════════════════════
# 1. Gene coordinates (hg38) from EnsDb  → gene_coords_dt + all_candidates + GRanges
# ══════════════════════════════════════════════════════════════════════════════
edb <- EnsDb.Hsapiens.v86

get_gene_coords <- function(symbols) {
  g <- genes(edb, filter = GeneNameFilter(symbols))
  seqlevelsStyle(g) <- "UCSC"                                   # -> chr1, chr2 ...
  g <- keepStandardChromosomes(g, pruning.mode = "coarse")
  dt <- as.data.table(g)[, .(gene = gene_name,
                             chr = as.character(seqnames),
                             gene_start = start, gene_end = end)]
  # if a symbol maps to several loci, keep the widest span
  dt[, w := gene_end - gene_start]
  setorder(dt, gene, -w)
  unique(dt, by = "gene")[, w := NULL][]
}

coords <- get_gene_coords(unique(c(lit_MEs, alias_map)))
found   <- coords$gene
missing <- setdiff(lit_MEs, found)
if (length(missing))
  message("No EnsDb coords for: ", paste(missing, collapse = ", "),
          "  (check aliases, e.g. VTRNA2-1 = nc886/MIR886)")

# `all_candidates` in the exact shape getFullPlotsTarget expects
# (needs columns: gene, chr, gene_start, gene_end)
all_candidates_MElit <- copy(coords)

# GRanges (gene ± 5 kb) for CpG extraction
candidates_gr_MElit <- unique(GRanges(
  seqnames = all_candidates_MElit$chr,
  ranges   = IRanges(start = all_candidates_MElit$gene_start - 5000,
                     end   = all_candidates_MElit$gene_end   + 5000),
  name     = all_candidates_MElit$gene))

# ══════════════════════════════════════════════════════════════════════════════
# 2. Write CpG list to extract (run the python extractor on the cluster, as S11)
# ══════════════════════════════════════════════════════════════════════════════
getCpG4pchuckle <- function(candidates_gr, layer){
  cpg_extract_list <- suppressWarnings(rbindlist(lapply(
    seq_along(candidates_gr), function(i) {
      gr   <- candidates_gr[i]
      hits <- findOverlaps(table3layers, gr)
      if (length(hits) == 0) return(NULL)
      data.table(chr_pos = table3layers$chr_pos[queryHits(hits)], gene = gr$name)
    })))
  unique_cpgs <- unique(na.omit(cpg_extract_list$chr_pos))
  message(sprintf("Total unique CpGs to extract: %d", length(unique_cpgs)))
  print(cpg_extract_list[, .N, by = gene][order(-N)])
  out_path <- here("B_MultiTissues/dataOut",
                   sprintf("targetCpGs_Alice_%s%s.txt", layer, format(Sys.time(), "%Y%m%d")))
  writeLines(unique_cpgs, out_path)
  message("Written to: ", out_path)
  out_path
}

me_cpg_file <- getCpG4pchuckle(candidates_gr_MElit, "MElit")

## Then on pchuckle (adapt the dated filename printed above):
# source /share/apps/source_files/python/python-3.13.0a6.source
# cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# python3 S00_extractRawMethylationForTargetCpG.py \
# --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/targetCpGs_Alice_MElit20260708.txt \
# --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/methylation_targetCpGs_Alice_MElit20260708.tsv \
# --minCov    10

# ══════════════════════════════════════════════════════════════════════════════
# 3. Plotting driver  (identical logic to S11's getFullPlotsTarget)
# ══════════════════════════════════════════════════════════════════════════════
getFullPlotsTarget <- function(meth, candidates_gr, all_candidates,
                               layer, gene_coords_dt = NULL, flank = 5000){
  
  setDT(meth); setDT(all_candidates)
  
  # gene_coords_dt straight from the candidate table (has gene_start/gene_end)
  if (is.null(gene_coords_dt)) {
    stopifnot(all(c("gene","chr","gene_start","gene_end") %in% names(all_candidates)))
    gene_coords_dt <- unique(all_candidates[, .(
      gene_name = gene, seqnames = chr,
      start = gene_start, end = gene_end,
      strand = "*", width = gene_end - gene_start + 1)])
    setorder(gene_coords_dt, gene_name, -width)
    gene_coords_dt <- unique(gene_coords_dt, by = "gene_name")
    message(sprintf("Built gene_coords_dt for %d genes", nrow(gene_coords_dt)))
  }
  
  loyfer_meta <- fread(here("B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv"))
  loyfer_meta[, source_tissue_celltype := paste0(`Source Tissue`, " - ", `Cell type`)]
  tissue_to_layer <- unique(loyfer_meta[, .(source_tissue_celltype, germ_layer = `Germ layer`)])
  meth <- merge(meth, tissue_to_layer, by = "source_tissue_celltype", all.x = TRUE)
  message(sprintf("Rows with missing germ_layer: %d", sum(is.na(meth$germ_layer))))
  
  if (!"pos" %in% names(meth)) meth[, pos := as.integer(sub(".*_", "", cpg_site))]
  if (!"chr" %in% names(meth)) meth[, chr := sub("_.*", "", cpg_site)]
  
  # ── inter-layer r ───────────────────────────────────────────────────────────
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
           interlayer_r = mean(abs(c(r_EM, r_EE, r_ME)), na.rm = TRUE),
           low_power    = n_EM < 8)
    }, by = cpg_site]
  }
  interlayer_loyfer <- compute_interlayer_r_loyfer(meth, min_patients = 3)
  message(sprintf("Inter-layer r computed for %d CpGs", nrow(interlayer_loyfer)))
  
  # ── intra-layer r ───────────────────────────────────────────────────────────
  compute_intralayer_r_loyfer <- function(meth, layer_name, min_obs = 3) {
    layer_meth <- meth[germ_layer == layer_name]
    eligible   <- layer_meth[, .(n = uniqueN(source_tissue_celltype)), by = patient_id][n > 1, patient_id]
    if (length(eligible) == 0) { message(sprintf("  No patients with >1 tissue in %s", layer_name)); return(NULL) }
    message(sprintf("  %s: %d eligible patients with >1 tissue", layer_name, length(eligible)))
    all_pairs <- rbindlist(lapply(eligible, function(pid) {
      pat  <- layer_meth[patient_id == pid]; tiss <- unique(pat$source_tissue_celltype)
      if (length(tiss) < 2) return(NULL)
      tpairs <- combn(tiss, 2, simplify = FALSE)
      rbindlist(lapply(tpairs, function(pr) {
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
  
  # ── build candidates_cpgs_full ──────────────────────────────────────────────
  seqlevelsStyle(candidates_gr) <- "UCSC"; candidates_gr <- keepStandardChromosomes(candidates_gr, pruning.mode = "coarse")
  seqlevelsStyle(table3layers)  <- "UCSC"; table3layers  <- keepStandardChromosomes(table3layers,  pruning.mode = "coarse")
  
  candidates_cpgs <- rbindlist(lapply(seq_along(candidates_gr), function(i) {
    gr <- candidates_gr[i]; hits <- findOverlaps(table3layers, gr)
    if (length(hits) == 0) return(NULL)
    dt <- as.data.table(table3layers[queryHits(hits)])[, .(chr_pos, alpha_endo, alpha_meso, alpha_ecto, alpha_geomean)]
    dt[, gene := candidates_gr$name[i]]; dt
  }), fill = TRUE)
  
  candidates_cpgs_r <- merge(candidates_cpgs,
                             interlayer_loyfer[, .(cpg_site, r_Endo_Meso, n_Endo_Meso, low_power)],
                             by.x = "chr_pos", by.y = "cpg_site", all.x = TRUE)
  candidates_cpgs_full <- merge(candidates_cpgs_r, intralayer_wide,
                                by.x = "chr_pos", by.y = "cpg_site", all.x = TRUE)
  for (col in c("r_intra_Endo","r_intra_Meso","r_intra_Ecto"))
    if (!col %in% names(candidates_cpgs_full)) candidates_cpgs_full[, (col) := NA_real_]
  
  message(sprintf("candidates_cpgs_full: %d rows, %d CpGs, %d genes",
                  nrow(candidates_cpgs_full), uniqueN(candidates_cpgs_full$chr_pos),
                  uniqueN(candidates_cpgs_full$gene)))
  
  # ── per-locus plots ─────────────────────────────────────────────────────────
  dir.create(here("B_MultiTissues/dataOut/figures/lociAlice"), showWarnings = FALSE, recursive = TRUE)
  standard_chrs <- paste0("chr", c(1:22,"X","Y","M"))
  all_windows <- GRanges(seqnames = gene_coords_dt$seqnames,
                         ranges = IRanges(start = gene_coords_dt$start - flank,
                                          end   = gene_coords_dt$end   + flank))
  table3layers_sub <- suppressWarnings(subsetByOverlaps(table3layers, all_windows))
  message(sprintf("table3layers subset: %d CpGs", length(table3layers_sub)))
  
  meth[, pos := as.integer(sub(".*_", "", cpg_site))]; meth[, chr := sub("_.*", "", cpg_site)]
  
  message("Pre-splitting meth by gene...")
  meth_by_gene <- lapply(unique(candidates_cpgs_full$gene), function(gn) {
    if (!gn %in% gene_coords_dt$gene_name) return(NULL)
    gc <- gene_coords_dt[gene_name == gn]
    meth[chr == as.character(gc$seqnames) & pos >= gc$start - flank & pos <= gc$end + flank]
  })
  names(meth_by_gene) <- unique(candidates_cpgs_full$gene)
  
  for (gn in unique(candidates_cpgs_full$gene)) {
    if (!gn %in% gene_coords_dt$gene_name) { message(sprintf("Skipping %s: no coords", gn)); next }
    cat_label <- layer[1]; if (is.na(cat_label)) cat_label <- "unknown"
    out_file  <- here(sprintf("B_MultiTissues/dataOut/figures/script12/locus_%s_%s.png", cat_label, gn))
    if (file.exists(out_file)) { message(sprintf("Skipping %s: exists", gn)); next }
    meth_gn <- meth_by_gene[[gn]]
    if (is.null(meth_gn) || nrow(meth_gn) == 0) { message(sprintf("Skipping %s: no meth", gn)); next }
    
    tryCatch({
      p <- plot_candidate_locus(
        gene_name_arg        = gn,
        candidates_cpgs_full = candidates_cpgs_full,
        interlayer_loyfer    = interlayer_loyfer,
        meth                 = meth_gn,
        gene_coords_dt       = gene_coords_dt,
        table3layers         = table3layers_sub,
        flank                = flank)
      ggplot2::ggsave(out_file, p, width = 13, height = 10, dpi = 300, bg = "white")
      message(sprintf("Saved: locus_%s_%s.png", cat_label, gn))
      rm(p); gc()
    }, error = function(e) message(sprintf("Failed for %s: %s", gn, e$message)))
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. Run  (after the python extraction has produced the TSV)
# ══════════════════════════════════════════════════════════════════════════════
me_meth_tsv <- here("gitignore","methylation_targetCpGs_Alice_MElit20260708.tsv")

getFullPlotsTarget(
  meth           = fread(me_meth_tsv),
  candidates_gr  = candidates_gr_MElit,
  all_candidates = all_candidates_MElit,
  layer          = "MElit")

message("S12 complete.")