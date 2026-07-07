#####################################################################
## S11_testCandidatesS10.R
## Test Alice's candidates: ME positive controls + Endo/Meso/Ecto
## candidates identified in S10 (data-driven + literature)
## Depends on: all_candidates.RDS (from S10_findTargetGenes.R)
## Produces: candidates_alice.pdf + per-locus PDFs
#####################################################################

library(here)
source(here("B_MultiTissues", "quiet_library.R"))
if (!exists("functionsLoaded"))  source(here("B_MultiTissues/03_exploreResults", "functions.R"))
if (!exists("table3layers"))     load(here("gitignore/fullTable3layers.Rda"))
if (!exists("wideFull"))         wideFull <- readRDS(here("gitignore/wide_script10_3layers_full.RDS"))

if (is.null(table3layers$percentile)) {
  table3layers$percentile <- ecdf(table3layers$alpha_geomean)(table3layers$alpha_geomean) * 100
}

library(EnsDb.Hsapiens.v86)

# ══════════════════════════════════════════════════════════════════════════════
# Load candidate list from S10
# ══════════════════════════════════════════════════════════════════════════════

all_candidates_meso <- read.csv(
  here("B_MultiTissues/dataOut/meso_hits_genesClusters.csv"))

all_candidates_endo <- read.csv(
  here("B_MultiTissues/dataOut/endo_hits_genesClusters.csv"))

# ── Build GRanges for extraction ─────────────────────────────────────────────
candidates_gr_meso <- unique(GRanges(
  seqnames = all_candidates_meso$chr,
  ranges   = IRanges(start = all_candidates_meso$gene_start - 5000,
                     end   = all_candidates_meso$gene_end   + 5000),
  name     = all_candidates_meso$gene))

candidates_gr_endo <- unique(GRanges(
  seqnames = all_candidates_endo$chr,
  ranges   = IRanges(start = all_candidates_endo$gene_start - 5000,
                     end   = all_candidates_endo$gene_end   + 5000),
  name     = all_candidates_endo$gene))

# ══════════════════════════════════════════════════════════════════════════════
# CpG extraction for pchuckle
# ══════════════════════════════════════════════════════════════════════════════

getCpG4pchuckle <- function(candidates_gr, layer){
  cpg_extract_list <- suppressWarnings(rbindlist(lapply(
    seq_along(candidates_gr), function(i) {
      gr   <- candidates_gr[i]
      hits <- findOverlaps(table3layers, gr)
      if (length(hits) == 0) return(NULL)
      data.table(chr_pos = table3layers$chr_pos[queryHits(hits)],
                 gene    = gr$name)
    })))
  
  unique_cpgs <- unique(cpg_extract_list$chr_pos)
  unique_cpgs <- unique_cpgs[!is.na(unique_cpgs)]
  message(sprintf("Total unique CpGs to extract: %d", length(unique_cpgs)))
  print(cpg_extract_list[, .N, by = gene][order(-N)])
  
  out_path <- here("B_MultiTissues/dataOut",
                   sprintf(paste0("targetCpGs_Alice_", layer, "%s.txt"), 
                           format(Sys.time(), "%Y%m%d")))
  writeLines(unique_cpgs, out_path)
  message("Written to: ", out_path)
}

getCpG4pchuckle(candidates_gr_meso, "meso")
# Written to: /home/alice/Documents/Research/GIT/2024_hvCpG/B_MultiTissues/dataOut/targetCpGs_Alice_meso20260707.txt
getCpG4pchuckle(candidates_gr_endo, "endo")
# Written to: /home/alice/Documents/Research/GIT/2024_hvCpG/B_MultiTissues/dataOut/targetCpGs_Alice_endo20260707.txt

## git push, then in pchuckle interactive:
# source /share/apps/source_files/python/python-3.13.0a6.source
# cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# python3 S00_extractRawMethylationForTargetCpG.py \
# --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/targetCpGs_Alice_meso20260707.txt \
# --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/methylation_targetCpGs_Alice_meso20260707.txt.tsv \
# --minCov    10

# python3 S00_extractRawMethylationForTargetCpG.py \
# --cpg_list  /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/targetCpGs_Alice_endo20260707.txt \
# --cpg_bed   /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
# --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
# --meta      /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
# --output    /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/methylation_targetCpGs_Alice_endo20260707.txt.tsv \
# --minCov    10










# ══════════════════════════════════════════════════════════════════════════════
# Load extracted methylation
# ══════════════════════════════════════════════════════════════════════════════

meth <- fread(here("gitignore/methylation_targetRegions_20260706.tsv"))
setDT(meth)

loyfer_meta <- fread(here("B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv"))
loyfer_meta[, source_tissue_celltype := paste0(`Source Tissue`, " - ", `Cell type`)]
tissue_to_layer <- unique(loyfer_meta[, .(source_tissue_celltype, germ_layer = `Germ layer`)])
meth <- merge(meth, tissue_to_layer, by = "source_tissue_celltype", all.x = TRUE)
message(sprintf("Rows with missing germ_layer: %d", sum(is.na(meth$germ_layer))))

if (!"pos" %in% names(meth)) meth[, pos := as.integer(sub(".*_", "", cpg_site))]
if (!"chr" %in% names(meth)) meth[, chr := sub("_.*", "", cpg_site)]

# ══════════════════════════════════════════════════════════════════════════════
# Inter-layer correlation (Endo x Meso, n=4 patients, low power)
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
message(sprintf("Inter-layer r computed for %d CpGs", nrow(interlayer_loyfer)))
# Inter-layer r computed for 44584 CpGs

# ══════════════════════════════════════════════════════════════════════════════
# Intra-layer correlation (within individuals, across tissues of same layer)
# ══════════════════════════════════════════════════════════════════════════════

compute_intralayer_r_loyfer <- function(meth_sub, layer_name, min_obs = 3) {
  layer_meth <- meth_sub[germ_layer == layer_name]
  eligible   <- layer_meth[, .(n = uniqueN(source_tissue_celltype)),
                            by = patient_id][n > 1, patient_id]
  if (length(eligible) == 0) {
    message(sprintf("  No patients with >1 tissue in %s", layer_name))
    return(NULL)
  }
  message(sprintf("  %s: %d eligible patients with >1 tissue", layer_name, length(eligible)))
  all_pairs <- rbindlist(lapply(eligible, function(pid) {
    pat  <- layer_meth[patient_id == pid]
    tiss <- unique(pat$source_tissue_celltype)
    if (length(tiss) < 2) return(NULL)
    tpairs <- combn(tiss, 2, simplify = FALSE)
    rbindlist(lapply(tpairs, function(pr) {
      d1 <- pat[source_tissue_celltype == pr[1], .(cpg_site, m1 = methylation)]
      d2 <- pat[source_tissue_celltype == pr[2], .(cpg_site, m2 = methylation)]
      dm <- merge(d1, d2, by = "cpg_site")
      dm[, patient_id := pid]
      dm[, tissue1    := pr[1]]
      dm[, tissue2    := pr[2]]
      dm
    }), fill = TRUE)
  }), fill = TRUE)
  if (is.null(all_pairs) || nrow(all_pairs) == 0) return(NULL)
  all_pairs[!is.na(m1) & !is.na(m2), {
    if (.N >= min_obs)
      list(r_intralayer = suppressWarnings(cor(m1, m2, method = "pearson")),
           n_obs = .N, layer = layer_name)
    else
      list(r_intralayer = NA_real_, n_obs = .N, layer = layer_name)
  }, by = cpg_site]
}

intralayer_loyfer <- rbindlist(lapply(c("Endo","Meso","Ecto"), function(lay) {
  message(sprintf("Computing intra-layer r for: %s", lay))
  compute_intralayer_r_loyfer(meth, lay, min_obs = 3)
}), fill = TRUE)

intralayer_wide <- dcast(intralayer_loyfer[!is.na(r_intralayer)],
                          cpg_site ~ layer, value.var = "r_intralayer")
setnames(intralayer_wide,
         intersect(c("Endo","Meso","Ecto"), names(intralayer_wide)),
         paste0("r_intra_", intersect(c("Endo","Meso","Ecto"), names(intralayer_wide))),
         skip_absent = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# Build candidates_cpgs_full: CpG-level data with all metrics
# ══════════════════════════════════════════════════════════════════════════════

library(GenomeInfoDb)

# put both on the same naming style and drop alt/haplotype scaffolds
seqlevelsStyle(candidates_gr) <- "UCSC"
candidates_gr <- keepStandardChromosomes(candidates_gr, pruning.mode = "coarse")

# table3layers is presumably already UCSC-standard, but make it explicit
seqlevelsStyle(table3layers) <- "UCSC"
table3layers  <- keepStandardChromosomes(table3layers,  pruning.mode = "coarse")

# ── extract CpGs per gene from table3layers ───────────────────────────────────
candidates_cpgs <- rbindlist(lapply(seq_along(candidates_gr), function(i) {
  gr   <- candidates_gr[i]
  hits <- findOverlaps(table3layers, gr)
  if (length(hits) == 0) return(NULL)
  dt <- as.data.table(table3layers[queryHits(hits)])[
    , .(chr_pos,
        alpha_endo,    alpha_meso,    alpha_ecto,
        alpha_geomean)]               # ← add this
  dt[, gene := candidates_gr$name[i]]
  dt
}), fill = TRUE)

# join inter-layer r
candidates_cpgs_r <- merge(
  candidates_cpgs,
  interlayer_loyfer[, .(cpg_site, r_Endo_Meso, n_Endo_Meso, low_power)],
  by.x = "chr_pos", by.y = "cpg_site", all.x = TRUE)

# join intra-layer r
candidates_cpgs_full <- merge(
  candidates_cpgs_r,
  intralayer_wide,
  by.x = "chr_pos", by.y = "cpg_site", all.x = TRUE)

# ensure intra-r columns exist
for (col in c("r_intra_Endo","r_intra_Meso","r_intra_Ecto")) {
  if (!col %in% names(candidates_cpgs_full))
    candidates_cpgs_full[, (col) := NA_real_]
}

message(sprintf("candidates_cpgs_full: %d rows, %d CpGs, %d genes",
                nrow(candidates_cpgs_full),
                uniqueN(candidates_cpgs_full$chr_pos),
                uniqueN(candidates_cpgs_full$gene)))
# candidates_cpgs_full: 44584 rows, 44584 CpGs, 49 genes

# ══════════════════════════════════════════════════════════════════════════════
# Per-locus plots using plot_candidate_locus() from functions.R
# ══════════════════════════════════════════════════════════════════════════════

dir.create(here("B_MultiTissues/dataOut/figures/lociAlice"), showWarnings = FALSE)

# ── Before the loop: pre-restrict to candidate regions ────────────────────────
# avoids searching all 21M CpGs for each gene inside plot_candidate_locus
standard_chrs <- paste0("chr", c(1:22,"X","Y","M"))

# build one merged window covering all candidate genes ± flank
all_windows <- GRanges(
  seqnames = gene_coords_dt$seqnames,
  ranges   = IRanges(start = gene_coords_dt$start - 5000,
                     end   = gene_coords_dt$end   + 5000))

# subset table3layers once
table3layers_sub <- suppressWarnings(
  subsetByOverlaps(table3layers, all_windows))
message(sprintf("table3layers subset: %d CpGs (vs %d total)",
                length(table3layers_sub), length(table3layers)))
# table3layers subset: 44584 CpGs (vs 24550220 total)

# subset meth once
meth[, pos := as.integer(sub(".*_", "", cpg_site))]
meth[, chr := sub("_.*",  "", cpg_site)]

candidate_regions_dt <- as.data.table(as.data.frame(all_windows))
meth_sub <- meth[chr %in% standard_chrs]  # already small for candidate genes

# ── Pre-split meth by gene to avoid passing full meth into every iteration ────
message("Pre-splitting meth by gene...")
meth_by_gene <- lapply(unique(candidates_cpgs_full$gene), function(gn) {
  if (!gn %in% gene_coords_dt$gene_name) return(NULL)
  gc    <- gene_coords_dt[gene_name == gn]
  x_min <- gc$start - 5000
  x_max <- gc$end   + 5000
  chr_w <- as.character(gc$seqnames)
  meth[chr == chr_w & pos >= x_min & pos <= x_max]
})
names(meth_by_gene) <- unique(candidates_cpgs_full$gene)
message("Done pre-splitting.")

for (gn in unique(candidates_cpgs_full$gene)) {
  if (!gn %in% gene_coords_dt$gene_name) {
    message(sprintf("Skipping %s: no gene coordinates", gn)); next
  }
  
  cat_label <- all_candidates[SYMBOL == gn, candidate_category]
  if (length(cat_label) == 0 || is.na(cat_label[1])) cat_label <- "unknown"
  cat_label <- cat_label[1]
  
  out_file <- here(sprintf(
    "B_MultiTissues/dataOut/figures/lociAlice/locus_%s_%s.png",
    cat_label, gn))
  if (file.exists(out_file)) {
    message(sprintf("Skipping %s: already exists", gn)); next
  }
  
  meth_gn <- meth_by_gene[[gn]]
  if (is.null(meth_gn) || nrow(meth_gn) == 0) {
    message(sprintf("Skipping %s: no meth data", gn)); next
  }
  
  tryCatch({
    p <- plot_candidate_locus(
      gene_name_arg        = gn,
      candidates_cpgs_full = candidates_cpgs_full,
      interlayer_loyfer    = interlayer_loyfer,
      meth                 = meth_gn,          # ← pre-subset meth
      gene_coords_dt       = gene_coords_dt,
      table3layers         = table3layers_sub,
      flank                = 5000)
    
    # ── save as PNG not PDF — dramatically smaller files ─────────────────────
    ggplot2::ggsave(
      filename = out_file,
      plot     = p,
      width    = 13, height = 10,
      dpi      = 300,
      bg       = "white")
    message(sprintf("Saved: locus_%s_%s.png", cat_label, gn))
    
    # free memory explicitly after each gene
    rm(p); gc()
    
  }, error = function(e) {
    message(sprintf("Failed for %s: %s", gn, e$message))
  })
}

message("S11 complete.")
