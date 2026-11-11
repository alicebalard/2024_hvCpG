#############################################################################################
## hvCpG twin validation: BLOOD + ADIPOSE in one script
##   Set TISSUE below to "Blood" or "Adipose"; everything else adapts automatically.
##
##   PART I   : matched top-10% vs mQTL control CpG sets (done in a previous script)
##   PART II  : per-sample beta matrix for the matched CpGs (cached)
##   PART III : twin structure (sample sheet -> pair_id + zygosity)
##   PART IV  : cell-composition covariates  <-- THIS DIFFERS BY TISSUE
##                Blood   : REFERENCE-BASED  (Houseman proportions, sorted-blood reference)
##                Adipose : REFERENCE-FREE   (PCA surrogate; no sorted-adipose reference exists)
##   PART V   : sample / pair accounting
##   Unified test runner + driver (top-1%, VanBaak SIV, VanBaak ESS; raw & composition-adjusted)
##   + per-set hierarchical clustering of samples (do twins cluster by pair?)
#############################################################################################

## ---- external dependencies: libraries, helper functions, previous ME sets ----
library(here)
if (!exists("libLoaded"))          source(here("B_MultiTissues", "quiet_library.R"))
if (!exists("functionsLoaded"))    source(here("B_MultiTissues/03_exploreResults", "functions.R"))
if (!exists("previousSIVprepared"))source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))

## ===========================================================================================
# TISSUE SWITCH — the ONE line to change to run the other tissue
## ===========================================================================================
TISSUE <- "Blood"
# TISSUE <- "Adipose"

## ---------------- shared parameters (tissue-independent) ----------------
tabRDS  <- here("gitignore/table3layers_coveredIn3dt.RDS")
covpat  <- "_5X\\.cov(\\.gz)?$"
figdir  <- here("B_MultiTissues/dataOut/figures/twin"); dir.create(figdir, FALSE, TRUE)
MINCOV  <- 5L; MINSAMP <- 10L; BG_MULT <- 15L; SEED <- 1234
N_SV    <- 5L            # (Adipose only) number of reference-free surrogate PCs to remove (tune: 2–10)

## RefSeq accession -> UCSC chromosome name (the .cov files use NC_000001.11 etc.)
refseq2ucsc <- c(
  NC_000001.11="chr1",  NC_000002.12="chr2",  NC_000003.12="chr3",  NC_000004.12="chr4",
  NC_000005.10="chr5",  NC_000006.12="chr6",  NC_000007.14="chr7",  NC_000008.11="chr8",
  NC_000009.12="chr9",  NC_000010.11="chr10", NC_000011.10="chr11", NC_000012.12="chr12",
  NC_000013.11="chr13", NC_000014.9="chr14",  NC_000015.10="chr15", NC_000016.10="chr16",
  NC_000017.11="chr17", NC_000018.10="chr18", NC_000019.10="chr19", NC_000020.11="chr20",
  NC_000021.9="chr21",  NC_000022.11="chr22", NC_000023.11="chrX",  NC_000024.10="chrY",
  NC_012920.1="chrM")

## a .cov filename -> sample id (MUST equal twins$sample and props$sample, i.e. the EGAN accession)
sample_from_file <- function(f) sub("_R1.*$", "", basename(f))

covdir       <- here(paste0("gitignore/", TISSUE))
matchfile <- here(paste0("gitignore/matched_CpGs_top_vs_control_", TISSUE, ".rds"))
betalongfile <- here(paste0("gitignore/betalongfile_", TISSUE, ".rds"))
samplecsv    <- here(paste0("B_MultiTissues/dataIn/sample_file_", tolower(TISSUE), "WGBS.csv"))
resultsfile  <- here(paste0("gitignore/twin_results_", TISSUE, ".rds"))
cachefile    <- here(paste0("gitignore/match_acc_", TISSUE, ".rds"))

## for blood
propsfile    <- here("B_MultiTissues/dataOut/blood_cellprops_cell_proportions.tsv")

if (TISSUE == "Blood") {
  alias_prefix <- "WB_"                                   # "WB_MZ1:1" -> "MZ1"
} else if (TISSUE == "Adipose") {
  alias_prefix <- "AT_"                                   # "AT_MZ1:1" -> "MZ1"
} else stop("TISSUE must be 'Blood' or 'Adipose'")

#############################################################################################
## PART I — GoDMC mQTL-matched hv / control pairs, filtered to the twin coverage
##   We take the PRE-matched pairs from cistrans_GoDMC_mytop10pchv_matched_control
## (each row = one top-10% hvCpG + its mQTL-matched control), map BOTH CpGs to hg38 
## (via `dico` from prepPreviousSIV), and keep only pairs where BOTH members are covered >=MINCOV (5x)
## in >=MINSAMP twin samples.
##   Output `out`: one row per CpG, labelled top1pct (hv) / control, with the CpG-pair kept.
##   NB: these controls are matched on mQTL properties, NOT on coverage/methylation, so the
##       two sets will NOT have identical mean_cov / mean_meth (that's expected here).
#############################################################################################
if (file.exists(matchfile)) {
  out <- readRDS(matchfile)  # already built — reuse
} else {
  ## 1. read the pre-matched pairs; map both CpGs (450k probe IDs) to hg38 chr_pos
  mc <- as.data.table(read.table(
    here("B_MultiTissues/01_dataPrep/cistrans_GoDMC_mytop10pchv_matched_control.txt"), header = TRUE))
  mc[, `:=`(hv_pos   = dico$chrpos_hg38[match(hvCpG_name,      dico$CpG)],
            ctrl_pos = dico$chrpos_hg38[match(controlCpG_name, dico$CpG)])]
  mc <- mc[!is.na(hv_pos) & !is.na(ctrl_pos)]                   # drop pairs not mappable to hg38
  mc[, cpg_pair := .I]                                          # keep the hv<->control pairing
  cat("GoDMC pairs with both CpGs mapped to hg38:", nrow(mc), "\n")
  # GoDMC pairs with both CpGs mapped to hg38: 16289 
  
  ## 2. coverage of every involved CpG in the twin .cov (per CpG: in how many samples >=5x)
  allpos <- data.table(chr_pos = unique(c(mc$hv_pos, mc$ctrl_pos))); setkey(allpos, chr_pos)
  files  <- list.files(covdir, pattern = covpat, full.names = TRUE)
  acc <- copy(allpos)[, `:=`(n_samp = 0L, total_cov = 0, total_M = 0)]; setkey(acc, chr_pos)
  for (f in files) {
    x <- fread(f, header = FALSE, select = c(1,2,5,6), col.names = c("chr","start","M","U"))
    x[, chr := refseq2ucsc[as.character(chr)]]; x <- x[!is.na(chr)]
    x[, chr_pos := paste0(chr, "_", start)]
    x <- x[allpos, on = "chr_pos", nomatch = 0L]
    x[, cov := M + U][cov >= MINCOV]                            # the _5X files are already >=5x
    s <- x[, .(np = .N, tc = sum(cov), tM = sum(M)), by = chr_pos]
    acc[s, on = "chr_pos", `:=`(n_samp = n_samp + np, total_cov = total_cov + tc, total_M = total_M + tM)]
    rm(x, s)
  }
  nsamp <- length(files)
  acc[, `:=`(mean_cov  = total_cov / nsamp,
             mean_meth = fifelse(total_cov > 0, total_M / total_cov, NA_real_),
             usable    = n_samp >= MINSAMP)]                    # >=5x in >=MINSAMP samples
  
  ## 3. keep only pairs where BOTH the hvCpG and its matched control are usable
  usable_pos <- acc[usable == TRUE, chr_pos]
  mc <- mc[hv_pos %in% usable_pos & ctrl_pos %in% usable_pos]
  cat("pairs with BOTH covered >=", MINCOV, "x in >=", MINSAMP, "samples: ", nrow(mc), "\n", sep = "")
  ## Blood: pairs with BOTH covered >=5x in >=10samples: 1668
  
  ## 4. build `out` (long: one labelled row per CpG; cpg_pair keeps the hv<->control link)
  out <- rbind(
    acc[match(mc$hv_pos,   chr_pos)][, `:=`(set = "top1pct", cpg_pair = mc$cpg_pair)],
    acc[match(mc$ctrl_pos, chr_pos)][, `:=`(set = "control", cpg_pair = mc$cpg_pair)]
  )[, .(chr_pos, mean_cov, mean_meth, n_samp, set, cpg_pair)]
  saveRDS(out, matchfile)
}
## sanity: counts per set (mean_cov / mean_meth need NOT match — controls are mQTL-matched)
print(out[, .(n = .N, mean_cov = mean(mean_cov), mean_meth = mean(mean_meth)), by = set])

## Blood
#        set     n mean_cov mean_meth
# 1: top1pct  1668 4.525580 0.6100107
# 2: control  1668 4.630118 0.6288850

## Adipose
#        set     n mean_cov mean_meth
# 1: top1pct   169 2.960842  0.514025
# 2: control   169 2.927950  0.650956

# Check that the mean-methylation distributions overlap well
ggplot(out, aes(mean_meth, colour=set)) + geom_density()
## Roughly, but the control is bimodal 0-1 while the hv have more intermediate values

#############################################################################################
## PART II — per-sample beta matrix for the matched CpGs (cached)
##   Re-read the .cov files keeping per-sample beta (M/coverage) at the GoDMC hv+control CpGs.
#############################################################################################
targ  <- unique(out[, .(chr_pos, set)]); setkey(targ, chr_pos)
files <- list.files(covdir, pattern = covpat, full.names = TRUE)
## PART VII needs per-sample coverage (to remove binomial sequencing noise): an old cache
## without a 'cov' column is rebuilt once.
betalong <- if (file.exists(betalongfile)) readRDS(betalongfile) else NULL
if (is.null(betalong) || !"cov" %in% names(betalong)) {
  betalong <- rbindlist(lapply(files, function(f) {
    x <- fread(f, header = FALSE, select = c(1,2,5,6), col.names = c("chr","start","M","U"))
    x[, chr := refseq2ucsc[as.character(chr)]]; x <- x[!is.na(chr)]
    x[, chr_pos := paste0(chr, "_", start)]
    x <- x[targ, on = "chr_pos", nomatch = 0L]
    x <- x[, cov := M + U][cov >= MINCOV]                        # (assign: filter now applied)
    x[, .(chr_pos, sample = sample_from_file(f), beta = M / cov, cov)]  # per-sample methylation + depth
  }))
  saveRDS(betalong, betalongfile)
}

betalong
## Blood
#            chr_pos          sample      beta   cov
# 1: chr10_100065428 EGAN00001342016 0.5000000     6
# 2:  chr10_10093262 EGAN00001342016 1.0000000     7
# 3: chr10_107165464 EGAN00001342016 0.6000000     5
# 4: chr10_110006146 EGAN00001342016 1.0000000     5
# 5: chr10_112838756 EGAN00001342016 1.0000000     6
# ---                                                
#   43587:   chrX_71053404 EGAN00001342042 0.1666667     6
# 43588:   chrX_75785267 EGAN00001342042 0.5833333    12
# 43589:   chrX_75785296 EGAN00001342042 0.5000000    10
# 43590:    chrX_9352941 EGAN00001342042 0.0000000     9
# 43591:   chrX_96337963 EGAN00001342042 0.6363636    11

## Adipose
#            chr_pos          sample      beta
# 1: chr10_121063891 EGAN00001341982 0.5714286
# 2:  chr10_18005634 EGAN00001341982 0.0000000
# 3:  chr10_58457785 EGAN00001341982 0.0000000
# 4:  chr10_60354605 EGAN00001341982 0.7500000
# 5:  chr10_64011988 EGAN00001341982 0.0000000
# ---                                          
#   4839:   chr8_32341667 EGAN00001342015 0.8571429
# 4840:   chr8_56682350 EGAN00001342015 0.2857143
# 4841:   chr8_80174059 EGAN00001342015 0.5714286
# 4842:  chr9_113057084 EGAN00001342015 0.0000000
# 4843:  chrX_100593990 EGAN00001342015 0.7142857

#############################################################################################
## PART III — twin structure
##   sample sheet: sample (EGAN) -> pair_id (e.g. MZ1) + zygosity (MZ/DZ), one row per sample.
##   Only the alias prefix differs between tissues (WB_ for blood, AT_ for adipose).
#############################################################################################
s <- fread(samplecsv)
twins <- unique(s[, .(sample = sample_accession_id, alias = sample_alias)])
twins[, pair_id  := sub(":.*$", "", sub(paste0("^", alias_prefix), "", alias))]  # "<pref>MZ1:1" -> "MZ1"
twins[, zygosity := sub("^([A-Z]{2}).*$", "\\1", pair_id)]                        # "MZ1" -> "MZ"
twins <- twins[!zygosity %in% "ST", .(sample, pair_id, zygosity)]                 # drop non-twin ST sample
print(twins[, .N, by = zygosity])                                                # MZ / DZ sample counts

## Blood
#    zygosity     N
# 1:       MZ    14
# 2:       DZ    12

## Adipose
# zygosity     N
# 1:       MZ    14
# 2:       DZ    12

#############################################################################################
## PART IV — cell-composition covariates  (THE PART THAT DIFFERS BY TISSUE)
##   Goal in both cases: a 'props' table (sample + covariate columns) and a 'celltypes' vector
##   naming those covariates, so the correction can regress them out of each CpG.
##
##   Blood   -> REFERENCE-BASED: sorted-blood reference exists (Jaffe/Reinius); read Houseman
##              proportions estimated per sample (script S03). The 6 proportions sum to ~1, so
##              the correction drops one (Gran) to avoid collinearity.
##   Adipose -> REFERENCE-FREE : no sorted-adipose reference; estimate latent between-sample
##              axes directly from the data (PCA on the CONTROL CpGs). The top PCs capture the
##              dominant genome-wide structure — largely cell composition — and are regressed
##              out the same way. These PCs also absorb genetic/ancestry/age/batch structure, so
##              this is a conservative / possibly over-correction (watch var_ratio; if it
##              collapses toward 1, lower N_SV).
#############################################################################################
if (TISSUE == "Blood") {
  ## ---- REFERENCE-BASED composition: Houseman deconvolution proportions (script S03) ----
  ## Keep every sample that HAS a usable estimate: drop only too_few_cpgs / collinear_reference
  ## (NA estimates); keep high_sampling_error & low_cpg_proportion (usable, just noisier).
  props_all <- fread(propsfile)
  props     <- props_all[!(flag_reason %like% "too_few_cpgs|collinear_reference")]
  celltypes <- c("CD8T","CD4T","NK","Bcell","Mono","Gran")        # the 6 sorted blood cell types
  cat("proportions kept:", nrow(props), "of", nrow(props_all), "samples\n")
  # proportions kept: 23 of 27 samples
} else {
  ## ---- REFERENCE-FREE composition surrogate: PCA on the CONTROL CpGs ----
  ## Built from the control CpGs only (independent of the hv foreground) so we don't circularly
  ## remove the hypervariability we test. Packaged to look exactly like blood's 'props'.
  ctrl_cpgs <- out[set == "control", chr_pos]
  bmat <- dcast(betalong[chr_pos %in% ctrl_cpgs], chr_pos ~ sample, value.var = "beta")
  rn <- bmat$chr_pos; bmat <- as.matrix(bmat[, -1]); rownames(bmat) <- rn
  
  ## adipose coverage is sparse: a fixed ">=80% of samples" filter can empty the matrix.
  ## Instead drop CpGs covered in NO sample, then keep the best-covered ones (never empty).
  present <- rowMeans(!is.na(bmat))
  bmat <- bmat[present > 0, , drop = FALSE]; present <- present[present > 0]
  if (nrow(bmat) < 10) stop("too few covered control CpGs for a composition surrogate (", nrow(bmat), ")")
  keep <- order(-present)[seq_len(min(5000L, nrow(bmat)))]       # most-covered control CpGs
  bmat <- bmat[keep, , drop = FALSE]
  cat("surrogate basis:", nrow(bmat), "control CpGs (best-covered)\n")
  
  ## impute the remaining NAs with the CpG (row) mean
  rmu <- rowMeans(bmat, na.rm = TRUE); na_idx <- which(is.na(bmat), arr.ind = TRUE)
  if (nrow(na_idx)) bmat[na_idx] <- rmu[na_idx[, 1]]
  
  ## drop zero-variance CpGs, then take the most-variable ones for speed/signal
  rv   <- rowSums((bmat - rowMeans(bmat))^2) / (ncol(bmat) - 1)
  bmat <- bmat[rv > 0, , drop = FALSE]; rv <- rv[rv > 0]
  bmat <- bmat[order(-rv)[seq_len(min(20000L, nrow(bmat)))], , drop = FALSE]
  
  ## PCA; never request more components than samples allow
  pcmod <- prcomp(t(bmat), center = TRUE, scale. = FALSE)         # samples x features -> PCs
  k     <- min(N_SV, ncol(pcmod$x))                              # cap at available PCs (<= nsamples-1)
  sv    <- pcmod$x[, seq_len(k), drop = FALSE]                   # sample scores on top PCs
  cat("surrogate: ", k, " PCs explain ",
      round(sum(pcmod$sdev[1:k]^2) / sum(pcmod$sdev^2) * 100, 1),
      "% of control-CpG variance\n", sep = "")
  props <- data.table(sample = rownames(sv))
  props[, (paste0("PC", seq_len(k))) := as.data.table(sv)]
  celltypes <- paste0("PC", seq_len(k))                          # the "covariates" to regress out
}

# Blood: proportions kept: 23 of 27 samples

# Adipose: 
# surrogate basis: 169 control CpGs (best-covered)
# surrogate: 5 PCs explain 34.5% of control-CpG variance

#############################################################################################
## PART V — sample / pair accounting: where does N drop?
##   (For adipose the PCA surrogate covers every sample, so "dropped by props" is usually empty.)
#############################################################################################
samp_B <- unique(betalong$sample); samp_twins <- twins$sample; samp_props <- props$sample
cat("beta:", length(samp_B), "| twins:", length(samp_twins), "| props:", length(samp_props), "\n")

# Adipose: beta: 33 | twins: 26 | props: 33 
cat("twins dropped by props:", paste(setdiff(samp_twins, samp_props), collapse = ", "), "\n")

# Adipose: NA

keep <- Reduce(intersect, list(samp_B, samp_twins, samp_props))  # samples surviving all three
print(twins[sample %in% keep][, .N, by = pair_id][N == 2][, .(complete_pairs = .N),
                                                          by = sub("[0-9]+$","",pair_id)])
## Blood
#       sub complete_pairs
# 1:     MZ              6
# 2:     DZ              4

## Adipose
#       sub complete_pairs
# 1:     MZ              7
# 2:     DZ              6

#############################################################################################
## PART VI — within-pair discordance test, run DIRECTLY on the GoDMC hv-vs-control set
##   Test (PRIMARY and only): within-pair methylation discordance, hv vs mQTL control,
##   split by zygosity (MZ vs DZ). The mQTL-matched control absorbs the cis-genetic component,
##   so residual excess discordance of hv over control — especially any present in MZ co-twins —
##   is variation genotype cannot explain (stochastic). An additive-genetic effect would instead
##   appear as larger excess in DZ than in MZ pairs.
##   restrict_same = TRUE uses the samples with usable composition covariates, so the raw and
##   adjusted runs sit on the SAME samples (composition correction is then the only difference).
#############################################################################################
run_twin_tests <- function(name = "top10%hv (mQTL-ctrl)", correct = TRUE, restrict_same = TRUE) {
  
  ## label each CpG hv(=ME) / control; drop any CpG that is ambiguously both across pairs
  so  <- unique(out[, .(chr_pos, set)])
  dup <- so$chr_pos[duplicated(so$chr_pos)]
  so  <- so[!chr_pos %in% dup]
  so[, lab := ifelse(set == "top1pct", "ME", "control")]
  set_of <- setNames(so$lab, so$chr_pos)
  
  dt <- copy(betalong)[chr_pos %in% names(set_of)]  # per-sample beta at matched CpGs
  dt[, set := factor(set_of[chr_pos], levels = c("control", "ME"))]
  dt <- merge(dt, twins, by = "sample")
  ## same samples as the corrected run (samples with usable composition covariates)
  if (restrict_same && exists("props")) {
    props_ok <- props$sample[complete.cases(props[, ..celltypes])]
    dt <- dt[sample %in% props_ok]
  }
  cat(sprintf("== %s | %s | %s ==  samples: %d (%d pairs)\n",
              name, basename(covdir), if (correct) "adjusted" else "raw",
              uniqueN(dt$sample), uniqueN(dt$pair_id)))
  
  ## composition correction: regress beta on the covariates, keep residual + CpG mean
  if (correct) {
    if (!exists("props")) stop("correct=TRUE but no 'props' in the environment")
    dt <- merge(dt, props[, c("sample", ..celltypes)], by = "sample")
    covn <- setdiff(celltypes, "Gran")   # blood: drop Gran (sum-to-1). adipose PCs: nothing to drop.
    dt[, beta := {
      P <- as.matrix(.SD); ok <- complete.cases(P) & !is.na(beta); r <- rep(NA_real_, .N)
      if (sum(ok) >= length(covn) + 2L) {
        fit <- .lm.fit(cbind(1, P[ok, , drop = FALSE]), beta[ok]); r[ok] <- fit$residuals + mean(beta[ok])
      }; r
    }, by = chr_pos, .SDcols = covn]
  }
  
  ## ---- within-pair discordance, hv vs control, split by zygosity ----
  ## per CpG x complete twin pair: within-pair |Δbeta|; averaged per pair per set; the per-pair
  ## EXCESS = hv - control. Pair is the unit of inference.
  wp <- dt[!is.na(beta), .(d = if (.N == 2L) abs(diff(beta)) else NA_real_),
           by = .(chr_pos, pair_id, zygosity, set)][!is.na(d)]
  pp <- dcast(wp[, .(mean_d = mean(d)), by = .(pair_id, zygosity, set)],
              pair_id + zygosity ~ set, value.var = "mean_d")
  pp[, excess := ME - control]
  t2a <- wilcox.test(pp$ME, pp$control, paired = TRUE, alternative = "greater")  # hv > control
  t2b <- wilcox.test(excess ~ zygosity, data = pp)                               # MZ vs DZ excess
  zyg <- pp[, .(mean_excess = mean(excess, na.rm = TRUE)), by = zygosity]
  ncpg <- wp[, .(n = uniqueN(chr_pos)), by = set]                                # CpGs contributing
  
  summary_row <- data.table(
    set = name, tissue = basename(covdir), corrected = correct,
    n_hv = ncpg[set == "ME", n], n_ctrl = ncpg[set == "control", n],
    n_pairs = nrow(pp), n_MZ = pp[zygosity == "MZ", .N], n_DZ = pp[zygosity == "DZ", .N],
    p_ME_gt_control = t2a$p.value,
    excess_MZ = zyg[zygosity == "MZ", mean_excess], excess_DZ = zyg[zygosity == "DZ", mean_excess],
    p_MZ_vs_DZ = t2b$p.value)
  print(summary_row)
  list(name = name, correct = correct, pp = pp, summary = summary_row,
       dt = dt[, .(chr_pos, sample, beta, set, pair_id, zygosity)])   # kept for clustering
}

## ======================================================================
## plotting (A: per-pair hv vs control | B: excess discordance by zygosity)
## ======================================================================
plot_set_tests <- function(res, outfile = NULL, bare = FALSE, panels = FALSE) {
  pp <- copy(res$pp); s <- res$summary; nm <- res$name
  cols_zyg <- c(MZ = "#0072B2", DZ = "#D55E00")
  lim <- range(c(pp$control, pp$ME))
  st <- if (s$p_ME_gt_control > 0.001) sprintf("p = %.3g", s$p_ME_gt_control) else "p < 0.001"
  
  pA <- ggplot(pp, aes(control, ME, colour = zygosity)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 3) + coord_equal(xlim = lim, ylim = lim) +
    scale_colour_manual(values = cols_zyg, name = NULL) +
    labs(x = "within-pair discordance, matched mQTL-control",
         y = "within-pair discordance, top 10% hv CpGs",
         title = "Per pair: hv vs control", subtitle = st) +
    theme_minimal(base_size = 12) + theme(legend.position = c(.9, .1))
  pB <- ggplot(pp, aes(zygosity, excess, colour = zygosity)) +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, fatten = 1.5) +
    geom_jitter(width = 0.08, size = 3) + scale_colour_manual(values = cols_zyg, guide = "none") +
    annotate("text", x = 1.5, y = max(pp$excess), label = sprintf("p = %.3g", s$p_MZ_vs_DZ),
             size = 3.2, vjust = 1) +
    labs(x = NULL, y = "excess discordance (hv − control)", title = "Excess discordance by zygosity") +
    theme_minimal(base_size = 12)
  
  ## panels = TRUE: return the two ggplots separately, for the multi-letter combined figure
  if (panels) return(list(scatter = pA, excess = pB))
  
  ## otherwise behave as before (standalone 2-panel block)
  fig <- cowplot::plot_grid(pA, pB, nrow = 1)
  if (!bare) {
    fig <- fig + patchwork::plot_annotation(
      title = nm,
      subtitle = sprintf("%s · %d twin pairs (%d MZ, %d DZ)",
                         s$tissue, nrow(pp), pp[zygosity == "MZ", .N], pp[zygosity == "DZ", .N]))
  }
  if (!is.null(outfile)) ggsave(outfile, fig, width = 9, height = 5, dpi = 300, bg = "white")
  fig
}

## ======================================================================
## Violin of ALL within-pair discordances (not the per-pair average)
##   Reuses res$dt (per-sample beta already labelled set / pair_id / zygosity)
##   and recomputes the per-CpG x per-pair |Δbeta| exactly as the test does
##   (same line as `wp` inside run_twin_tests). One point = one CpG in one pair.
##   x = set (control vs hv); points coloured by twin pair; MZ / DZ in facets.
## ======================================================================
plot_discordance_violin <- function(res, outfile = NULL) {
  
  ## 1. rebuild the per-CpG, per-pair discordances from the returned long table.
  ##    For each CpG x complete pair (exactly 2 co-twins), |Δbeta| = |beta1 - beta2|.
  wp <- res$dt[!is.na(beta),
               .(d = if (.N == 2L) abs(diff(beta)) else NA_real_),
               by = .(chr_pos, pair_id, zygosity, set)][!is.na(d)]
  
  ## 2. readable x-axis labels; keep control on the left, hv on the right
  wp[, grp := factor(ifelse(set == "ME", "top 10% hv", "control"),
                     levels = c("control", "top 10% hv"))]
  
  ## 3. violin = full distribution of all discordances; jittered points coloured by pair
  g <- ggplot(wp, aes(grp, d)) +
    geom_violin(aes(fill = grp), colour = NA, alpha = 0.25, scale = "width") +      # all CpGs' spread
    geom_jitter(aes(colour = pair_id), width = 0.18, size = 0.4, alpha = 0.3) +     # every CpG, by pair
    stat_summary(fun = median, geom = "point", shape = 95, size = 8, colour = "black") +  # median tick
    scale_fill_manual(values = c(control = "grey70", "top 10% hv" = "#B2182B"), guide = "none") +
    facet_wrap(~ zygosity) +                                                        # MZ vs DZ panels
    labs(x = NULL, y = "within-pair discordance  (per CpG, per pair)",
         colour = "twin pair",
         title = sprintf("%s — all within-pair discordances", res$name),
         subtitle = sprintf("%s · %s top 10pc hv CpGs with matched mQTL controls CpGs",
                            res$summary$tissue, length(unique(wp$chr_pos))/2)) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +       # legible legend
    theme_minimal(base_size = 12)
  
  if (!is.null(outfile)) ggsave(outfile, g, width = 10, height = 5, dpi = 300, bg = "white")
  g
}

## ======================================================================
## driver — the GoDMC hv-vs-control comparison, raw & composition-adjusted
## ======================================================================
res_adj <- run_twin_tests(correct = TRUE)    # composition-adjusted
res_raw <- run_twin_tests(correct = FALSE)   # raw
res_all <- rbindlist(list(res_adj$summary, res_raw$summary), fill = TRUE)
print(res_all); saveRDS(res_all, resultsfile)
#                 set tissue corrected  n_hv n_ctrl n_pairs  n_MZ  n_DZ p_ME_gt_control  excess_MZ  excess_DZ p_MZ_vs_DZ
# top10%hv (mQTL-ctrl)  Blood      TRUE  1668   1668      10     6     4    0.0009765625 0.04576979 0.03373396  0.1142857
# top10%hv (mQTL-ctrl)  Blood     FALSE  1668   1668      10     6     4    0.0009765625 0.05645989 0.06200971  0.4761905

#                 set  tissue corrected  n_hv n_ctrl n_pairs  n_MZ  n_DZ p_ME_gt_control  excess_MZ  excess_DZ p_MZ_vs_DZ
# top10%hv (mQTL-ctrl) Adipose      TRUE   168    166      13     7     6    0.0001220703 0.06641051 0.06941351  0.8356643
# top10%hv (mQTL-ctrl) Adipose     FALSE   168    166      13     7     6    0.0001220703 0.07135265 0.07997177  0.7307692

## Plot discordance averaged per twin pair
plot_set_tests(res_adj, here(paste0("B_MultiTissues/dataOut/figures/twin/", TISSUE, "_CompositionCorrected.pdf")))
plot_set_tests(res_raw, here(paste0("B_MultiTissues/dataOut/figures/twin/", TISSUE, "_raw.pdf")))

## Plot discordance not averaged (violin plot)
plot_discordance_violin(res_adj, here(paste0("B_MultiTissues/dataOut/figures/twin/", TISSUE, "_discordance_violin.pdf")))

## Save object (run for both tissues) for combined plot:
saveRDS(res_adj, here(paste0("gitignore/res_adj_", TISSUE, ".rds")))

#############################################################################################
## PART VII — twin variance decomposition (simplified from Hannon et al. 2018, PLoS Genet
##            14:e1007544, doi:10.1371/journal.pgen.1007544)
##   Hannon et al. split each site's methylation variance into additive genetic (A), shared
##   environment (C) and unique environment / stochastic + measurement error (E) by comparing
##   MZ and DZ co-twin similarity (ACE models, 426 MZ + 306 DZ pairs). They found A is HIGHEST at
##   highly variable sites and at intermediate methylation (0.2-0.8). Our hvCpGs are, by
##   construction, more variable and more intermediate than their controls, so if hypervariability
##   were mostly genetic we would expect A(hv) > A(control). A(hv) <= A(control), with the extra
##   variance of hv sites falling mainly into E, supports a non-genetic origin.
##
##   With 6-7 MZ and 4-6 DZ pairs a per-CpG ACE fit is impossible, so we use the moment
##   (Falconer / within-pair mean-square) version, pooled over each CpG set:
##     V    = variance of beta across all twins                    (per CpG)
##     W_MZ = mean over MZ pairs of (b1 - b2)^2 / 2                (= E: MZ share all genes)
##     W_DZ = mean over DZ pairs of (b1 - b2)^2 / 2                (= E + A/2)
##   averaged over the CpGs of a set:
##     absolute : sA = 2 (W_DZ - W_MZ) ;  sE = W_MZ ;  sC = V - sA - sE
##     proportion: A = sA / V ;  E = sE / V ;  C = 1 - A - E
##   Binomial sequencing noise is removed first: each beta carries variance b(1-b)/(cov-1),
##   subtracted from V and from each pair's (b1-b2)^2/2 (otherwise noise inflates E).
##   Extra variance of hv over control, and the share of it that is non-genetic / genetic:
##     extra_E = (sE_hv - sE_ctrl) / (V_hv - V_ctrl) ; extra_A = (sA_hv - sA_ctrl) / (V_hv - V_ctrl)
##   Uncertainty: 95% CI from a bootstrap over twin pairs (resampled within zygosity) —
##   pairs, not CpGs, are the independent units. Estimates are unconstrained (can fall
##   outside 0-1). meth_range = c(0.2, 0.8) repeats it on intermediate CpGs only (Hannon's
##   definition), removing the methylation-level difference between hv and controls.
#############################################################################################
variance_decomp <- function(res, B = 1000L, meth_range = c(0, 1), seed = SEED) {
  ## per-sample beta (raw or composition-corrected, from run_twin_tests) + depth for the noise
  d <- merge(res$dt, betalong[, .(chr_pos, sample, braw = beta, cov)], by = c("chr_pos", "sample"))
  d <- d[!is.na(beta)]
  d[, noise := braw * (1 - braw) / (cov - 1)]                     # binomial variance of beta
  keep_cpg <- d[, .(m = mean(braw)), by = chr_pos][m >= meth_range[1] & m <= meth_range[2], chr_pos]
  d <- d[chr_pos %in% keep_cpg]
  
  ## one row per CpG x complete twin pair
  pr <- d[, if (.N == 2L) .(b1 = beta[1], b2 = beta[2], n1 = noise[1], n2 = noise[2]),
          by = .(chr_pos, set, pair_id, zygosity)]
  pr[, g := .GRP, by = .(chr_pos, set)]
  gset <- pr[, .(set = as.character(set[1])), keyby = g]$set
  pr[, `:=`(s1 = b1 + b2, s2 = b1^2 + b2^2, sn = n1 + n2,
            wd = ((b1 - b2)^2 - n1 - n2) / 2,                      # noise-corrected within-pair var
            mz = as.numeric(zygosity == "MZ"))]
  pairs <- unique(pr[, .(pair_id, zygosity)])
  
  ## estimator for a vector of pair weights (all 1 = observed data; bootstrap counts otherwise)
  est <- function(wpair) {
    wt <- wpair[pr$pair_id]
    M  <- rowsum(cbind(2 * wt, wt * pr$s1, wt * pr$s2, wt * pr$sn,
                       wt * pr$mz, wt * pr$mz * pr$wd, wt * (1 - pr$mz), wt * (1 - pr$mz) * pr$wd),
                 pr$g, reorder = TRUE)
    V    <- (M[, 3] - M[, 2]^2 / M[, 1]) / (M[, 1] - 1) - M[, 4] / M[, 1]   # total var - noise
    W_MZ <- M[, 6] / M[, 5]; W_DZ <- M[, 8] / M[, 7]
    st   <- gset[as.integer(rownames(M))]
    ok   <- is.finite(V) & is.finite(W_MZ) & is.finite(W_DZ)
    agg  <- function(s) { i <- ok & st == s
    v <- mean(V[i]); sE <- mean(W_MZ[i]); sA <- 2 * (mean(W_DZ[i]) - sE)
    c(V = v, sA = sA, sE = sE, sC = v - sA - sE, A = sA / v, E = sE / v, C = 1 - (sA + sE) / v, n = sum(i)) }
    h <- agg("ME"); k <- agg("control"); dV <- h[["V"]] - k[["V"]]
    c(setNames(h, paste0(names(h), "_hv")), setNames(k, paste0(names(k), "_ctrl")),
      A_hv_minus_ctrl = h[["A"]] - k[["A"]], E_hv_minus_ctrl = h[["E"]] - k[["E"]],
      extra_A = (h[["sA"]] - k[["sA"]]) / dV, extra_E = (h[["sE"]] - k[["sE"]]) / dV)
  }
  
  pt <- est(setNames(rep(1, nrow(pairs)), pairs$pair_id))
  set.seed(seed)
  pid  <- split(pairs$pair_id, pairs$zygosity)
  boot <- t(replicate(B, {
    smp <- unlist(lapply(pid, function(p) sample(p, length(p), replace = TRUE)), use.names = FALSE)
    est(c(table(factor(smp, levels = pairs$pair_id))))
  }))
  ci  <- apply(boot, 2, quantile, c(.025, .975), na.rm = TRUE)
  tab <- data.table(tissue = basename(covdir), corrected = res$correct,
                    meth_range = paste(meth_range, collapse = "-"),
                    quantity = names(pt), estimate = pt, lo95 = ci[1, ], hi95 = ci[2, ])
  tab[, p_boot_le0 := colMeans(boot <= 0, na.rm = TRUE)[quantity]]   # one-sided bootstrap p (H1: > 0)
  list(table = tab, boot = boot)
}

vd_adj <- variance_decomp(res_adj)                              # composition-adjusted (main)
vd_raw <- variance_decomp(res_raw)                              # raw
vd_int <- variance_decomp(res_adj, meth_range = c(0.2, 0.8))    # intermediate CpGs only
vd_all <- rbindlist(list(vd_adj$table, vd_raw$table, vd_int$table))
print(vd_all[!grepl("^n_", quantity)], digits = 3)
saveRDS(vd_all, here(paste0("gitignore/twin_varDecomp_", TISSUE, ".rds")))
## Key rows: A_hv vs A_ctrl (expectation from Hannon: A_hv > A_ctrl if hv were genetic),
##           A_hv_minus_ctrl (CI), extra_E (share of hv's extra variance that is non-genetic).

## ============================================================
## FINAL figure: Blood (A, B, C) over Adipose (D, E, F)
##   A/D = violin of all discordances ; B/E = per-pair scatter ; C/F = excess by zygosity
##   One bold line over each row names the tissue + hv CpG count.
## ============================================================
runFinal <- T   # set TRUE once res_adj_Blood.rds AND res_adj_Adipose.rds both exist

if (runFinal) {
  res_blood   <- readRDS(here("gitignore/res_adj_Blood.rds"))      # blood, composition-corrected
  res_adipose <- readRDS(here("gitignore/res_adj_Adipose.rds"))    # adipose, surrogate-corrected
  
  ## build one tissue row: violin + scatter + excess, each lettered, under a single header line
  build_row <- function(res, lab3) {
    v  <- plot_discordance_violin(res) +
      labs(title = "All within-pair discordances", subtitle = NULL)   # tissue text goes in the header
    ps <- plot_set_tests(res, panels = TRUE)                                 # list(scatter, excess)
    row <- cowplot::plot_grid(v, ps$scatter, ps$excess,
                              labels = lab3, label_size = 16,
                              nrow = 1, rel_widths = c(1.5, 1, 1))           # violin is 2 facets wide
    header <- cowplot::ggdraw() +
      cowplot::draw_label(
        sprintf("%s · %s top 10%% hv CpGs with matched mQTL control CpGs",
                res$summary$tissue, format(res$summary$n_hv, big.mark = ",")),
        fontface = "bold", x = 0.01, hjust = 0, size = 13)
    cowplot::plot_grid(header, row, ncol = 1, rel_heights = c(0.08, 1))      # header line over the row
  }
  
  rowBlood   <- build_row(res_blood,   c("A", "B", "C"))
  rowAdipose <- build_row(res_adipose, c("D", "E", "F"))
  final <- cowplot::plot_grid(rowBlood, rowAdipose, ncol = 1)               # two stacked rows
  
  ggsave(here("B_MultiTissues/dataOut/figures/twin/Twin_final_blood_adipose.pdf"),
         final, width = 18, height = 11, dpi = 300, bg = "white")
  ggsave(here("B_MultiTissues/dataOut/figures/twin/Twin_final_blood_adipose.png"),
         final, width = 18, height = 11, dpi = 300, bg = "white")
}

d <- merge(res_raw$dt, betalong[, .(chr_pos, sample, cov)], by = c("chr_pos", "sample"))
d[, noise := beta * (1 - beta) / (cov - 1)]
wp2 <- d[, if (.N == 2L) .(wd = ((beta[1] - beta[2])^2 - noise[1] - noise[2]) / 2),
         by = .(chr_pos, set, pair_id, zygosity)]
pp2 <- dcast(wp2[, .(w = mean(wd)), by = .(pair_id, zygosity, set)],
             pair_id + zygosity ~ set, value.var = "w")
wilcox.test(pp2$ME, pp2$control, paired = TRUE, alternative = "greater")
pp2[, .(ME = mean(ME), control = mean(control)), by = zygosity]
