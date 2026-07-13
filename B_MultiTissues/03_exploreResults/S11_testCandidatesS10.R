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

getFullPlotsTarget <- function(meth, candidates_gr, all_candidates,
                               layer, gene_coords_dt = NULL){
  
  setDT(meth)
  setDT(all_candidates)                      # in case it came from read.csv
  
  # ── Build gene_coords_dt from the candidate CSV itself ──────────────────────
  if (is.null(gene_coords_dt)) {
    stopifnot(all(c("gene","chr","gene_start","gene_end") %in% names(all_candidates)))
    gene_coords_dt <- unique(all_candidates[, .(
      gene_name = gene,
      seqnames  = chr,
      start     = gene_start,
      end       = gene_end,
      strand    = "*",
      width     = gene_end - gene_start + 1)])
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
  
  # ══════════════════════════════════════════════════════════════════════════════
  # Inter-layer correlation (Endo x Meso, n=4 patients, low power)
  # ══════════════════════════════════════════════════════════════════════════════
  
  compute_interlayer_r_loyfer <- function(meth, min_patients = 3) {
    agg  <- meth[, .(methylation = mean(methylation, na.rm = TRUE)),
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
  
  compute_intralayer_r_loyfer <- function(meth, layer_name, min_obs = 3) {
    layer_meth <- meth[germ_layer == layer_name]
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
    
    cat_label <- layer 
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
}

getFullPlotsTarget(
  meth          = fread(here("gitignore/methylation_targetCpGs_Alice_endo20260707.txt.tsv")),
  candidates_gr = candidates_gr_endo,
  all_candidates = all_candidates_endo,
  layer         = "endo")

getFullPlotsTarget(
  meth          = fread(here("gitignore/methylation_targetCpGs_Alice_meso20260707.txt.tsv")),
  candidates_gr = candidates_gr_meso,
  all_candidates = all_candidates_meso,
  layer         = "meso")


### ME litt (to improve!)
library(data.table)

build_candidate_enrichment <- function(gene_cpgs,           # per-gene CpG table (long)
                                       hit_col,              # e.g. "is_ME_hit" / "is_meso_hit"
                                       uni,                  # full CpG universe (cpg_site, category)
                                       hits_full = NULL,     # optional: hit table w/ r_own,r_cross
                                       literature_genes = NULL, # char vec of named MEs to append
                                       genes_gr = NULL,      # for gene_length / cpgs_per_kb
                                       cluster_gap = 100) {
  
  setDT(gene_cpgs)
  stopifnot(hit_col %in% names(gene_cpgs))
  
  # ── genome-wide background rate: P(a CpG is a hit) ──────────────────────────
  # numerator = total hits genome-wide; denominator = all tested CpGs
  n_hits_total <- uniqueN(gene_cpgs[get(hit_col) == TRUE, cpg_site])
  n_uni        <- nrow(uni)
  p0           <- n_hits_total / n_uni      # background hit probability
  
  # ── per-gene counts ─────────────────────────────────────────────────────────
  per_gene <- gene_cpgs[, .(
    n_cpgs     = sum(get(hit_col)),                 # hits in this gene
    total_cpgs = uniqueN(cpg_site),                 # all CpGs in gene+flank
    chr        = chr[1],
    span_bp    = max(pos) - min(pos) + 1
  ), by = .(SYMBOL = gene)][total_cpgs > 0]
  
  # ── binomial enrichment vs background ───────────────────────────────────────
  per_gene[, expected   := total_cpgs * p0]
  per_gene[, enrichment := n_cpgs / expected]
  per_gene[, p_binom := mapply(function(k, n)
    binom.test(k, n, p = p0, alternative = "greater")$p.value,
    n_cpgs, total_cpgs)]
  per_gene[, p_adj := p.adjust(p_binom, method = "BH")]
  
  # ── gene length / density ───────────────────────────────────────────────────
  if (!is.null(genes_gr)) {
    gl <- data.table(SYMBOL = genes_gr$symbol,
                     gene_length = width(genes_gr))
    gl <- gl[!is.na(SYMBOL)][, .(gene_length = max(gene_length)), by = SYMBOL]
    per_gene <- merge(per_gene, gl, by = "SYMBOL", all.x = TRUE)
  } else per_gene[, gene_length := span_bp]
  per_gene[, cpgs_per_kb := n_cpgs / (gene_length / 1000)]
  
  # ── cluster peak (max hits within cluster_gap bp) ───────────────────────────
  clpk <- gene_cpgs[get(hit_col) == TRUE][
    order(chr, pos),
    {
      if (.N == 0) .(cluster_peak = 0L)
      else {
        cl <- cumsum(c(TRUE, diff(pos) > cluster_gap))   # new cluster when gap exceeded
        .(cluster_peak = max(tabulate(cl)))
      }
    }, by = .(SYMBOL = gene)]
  per_gene <- merge(per_gene, clpk, by = "SYMBOL", all.x = TRUE)
  per_gene[is.na(cluster_peak), cluster_peak := 0L]
  
  # ── optional correlation summaries from the hit table ───────────────────────
  if (!is.null(hits_full) && all(c("r_own","r_cross") %in% names(hits_full))) {
    hf <- as.data.table(hits_full)
    # need a gene mapping on hits_full; expect a `gene` column, else skip
    if ("gene" %in% names(hf)) {
      rs <- hf[, .(mean_r_own   = mean(r_own,   na.rm = TRUE),
                   mean_r_cross = mean(r_cross, na.rm = TRUE)),
               by = .(SYMBOL = gene)]
      per_gene <- merge(per_gene, rs, by = "SYMBOL", all.x = TRUE)
    }
  }
  for (col in c("mean_r_own","mean_r_cross"))
    if (!col %in% names(per_gene)) per_gene[, (col) := NA_real_]
  per_gene[, spec_score := NA_real_]   # fill with your own definition if wanted
  
  per_gene[, `:=`(source = "data_driven", candidate_category = "ME")]
  
  # ── append literature genes not recovered (all-NA stats, as in your table) ──
  if (!is.null(literature_genes)) {
    missing_lit <- setdiff(literature_genes, per_gene$SYMBOL)
    if (length(missing_lit)) {
      lit_dt <- data.table(SYMBOL = missing_lit)
      lit_dt[, source := "literature"][, candidate_category := "ME"]
      per_gene <- rbind(per_gene, lit_dt, fill = TRUE)
    }
    # mark the recovered literature ones too
    per_gene[SYMBOL %in% literature_genes, source := "literature"]
  }
  
  setorder(per_gene, -enrichment, p_adj, na.last = TRUE)
  per_gene[]
}

lit_MEs <- c("VTRNA2-1","PM20D1","POMC","CYP2E1","PAX8","SPATC1L",
             "HLA-DQA1","HLA-A","HLA-DPA1","ZXDA","GPRASP1","MMGT1",
             "IDS","CNKSR2","IRS4","ZMAT1")

all_candidates_ME <- build_candidate_enrichment(
  gene_cpgs        = ME_res$gene_cpgs,      # from annotate_layer_hits(ME_hits,...)
  hit_col          = "is_ME_hit",
  uni              = uni,
  hits_full        = ME_hits,               # if it carries gene + r_own/r_cross
  literature_genes = lit_MEs,
  genes_gr         = genes_gr,
  cluster_gap      = 100)

all_candidates_ME <- read.csv(
  here("B_MultiTissues/dataOut/all_candidates.csv"))
all_candidates_ME <- all_candidates_ME[all_candidates_ME$source %in% "literature" &
                                         all_candidates_ME$candidate_category %in% "ME",]

candidates_gr_ME <- unique(GRanges(
  seqnames = all_candidates_MEslitt$chr,
  ranges   = IRanges(start = all_candidates_meso$gene_start - 5000,
                     end   = all_candidates_meso$gene_end   + 5000),
  name     = all_candidates_meso$gene))

getFullPlotsTarget(
  meth          = fread(here("gitignore/methylation_layerspecific_and_ME.tsv")),
  candidates_gr = candidates_gr_ME,
  all_candidates = all_candidates_ME,
  layer         = "MElitterature") ## tbc

meth          = fread(here("gitignore/methylation_layerspecific_and_ME.tsv"))

meth




message("S11 complete.")
