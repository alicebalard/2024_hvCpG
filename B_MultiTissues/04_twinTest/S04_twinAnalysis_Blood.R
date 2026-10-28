#############################################################################################
## hvCpG blood-twin validation
##   PART I   : build coverage+methylation-matched top-1% vs control CpG sets (cached)
##   PART II  : per-sample methylation (beta) matrix for the matched CpGs (cached)
##   PART III : twin sheet + blood cell-composition proportions
##   PART IV  : sample / pair accounting (where does N drop?)
##   PART V   : composition-adjusted beta
##   PART VI  : Test 1 — inter-individual variance (out-of-sample hypervariability)
##   PART VII : Test 2 — within-pair discordance, MZ vs DZ (stochastic vs genetic)
#############################################################################################
library(here); library(data.table); library(ggplot2); library(patchwork)

## ---------------- shared parameters ----------------
tabRDS     <- here("gitignore/table3layers_coveredIn3dt.RDS")         # 20M-CpG hvCpG scores
covdir     <- here("gitignore/Blood")                                  # blood-twin .cov files
covpat     <- "_5X\\.cov(\\.gz)?$"                                     # the 5x-prefiltered files
targetfile <- here("gitignore/match_target.rds")                       # cached top + candidate pool
cachefile  <- here("gitignore/match_acc.rds")                          # cached streamed accumulator
matchfile  <- here("gitignore/matched_CpGs_top_vs_control.rds")        # PART I output
betalongfile  <- here("gitignore/betalongfile.rds")        # PART II output
propsfile  <- here("B_MultiTissues/dataOut/blood_cellprops_cell_proportions.tsv")
samplecsv  <- here("B_MultiTissues/dataIn/sample_file_bloodWGBS.csv")

MINCOV  <- 5L      # a CpG is "covered" in a sample at >= this depth
MINSAMP <- 10L     # a CpG is "usable" if covered in >= this many samples
BG_MULT <- 15L     # candidate controls drawn = BG_MULT x (#top CpGs); large enough to fill every bin
SEED    <- 1234
celltypes <- c("CD8T","CD4T","NK","Bcell","Mono","Gran")

## chromosomes in the .cov files are RefSeq accessions; our CpG ids use chr1 etc.
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

#############################################################################################
## PART I — matched top-1% vs control CpG sets
##   Stream the 5x .cov once over (top-1% hvCpGs) + (BG_MULT x random non-top CpGs),
##   accumulate coverage & methylation, keep the usable ones (>=5x in >=10 samples),
##   then draw controls that match the top set on (mean coverage, mean methylation).
##   Cached: reruns skip the heavy streaming unless the signature changes.
#############################################################################################
if (file.exists(matchfile)) {
  out <- readRDS(matchfile)   # already built — reuse
} else {
  
  ## (1) target set = top-1% + candidate controls (load the 20M-row table once)
  sig <- list(mincov = MINCOV, bg_mult = BG_MULT, seed = SEED)
  if (file.exists(targetfile) && identical(readRDS(targetfile)$sig, sig)) {
    tb <- readRDS(targetfile); target <- tb$target; n_top <- tb$n_top
  } else {
    dtab    <- readRDS(tabRDS)[, .(chr_pos, pct = percentile_logBF_per_ds_allLayers)]
    top_pos <- dtab[pct >= 99, chr_pos] # top 1%
    n_top   <- length(top_pos)
    set.seed(SEED)
    bg_pos  <- sample(dtab[pct < 99, chr_pos], BG_MULT * n_top)   # candidate controls
    target  <- data.table(chr_pos = c(top_pos, bg_pos),
                          is_top  = c(rep(TRUE, n_top), rep(FALSE, length(bg_pos))))
    setkey(target, chr_pos)
    saveRDS(list(target = target, n_top = n_top, sig = sig), targetfile)
    rm(dtab, top_pos, bg_pos); gc()
  }
  cat("target CpGs:", nrow(target), "(top:", n_top, ")\n")
  
  ## (2) one streaming pass: per-CpG summed coverage + methylation across samples
  sig_acc <- c(sig, list(n_target = nrow(target),
                         tsig = paste(target$chr_pos[1], target$chr_pos[nrow(target)])))
  if (file.exists(cachefile)) {
    cache <- readRDS(cachefile)
    if (!identical(cache$sig, sig_acc)) stop("acc cache stale — delete ", cachefile)
    acc <- cache$acc; done <- cache$done
  } else {
    acc  <- copy(target)[, `:=`(n_samp = 0L, total_cov = 0, total_M = 0)]
    setkey(acc, chr_pos); done <- character(0)
  }
  files    <- list.files(covdir, pattern = covpat, full.names = TRUE)
  newfiles <- files[!(basename(files) %in% done)]
  for (f in newfiles) {
    x <- fread(f, header = FALSE, select = c(1,2,5,6), col.names = c("chr","start","M","U"))
    x[, chr := refseq2ucsc[as.character(chr)]]
    x <- x[!is.na(chr)]
    x[, chr_pos := paste0(chr, "_", start)]
    x <- x[target, on = "chr_pos", nomatch = 0L]              # keep target CpGs only
    x[, cov := M + U][cov >= MINCOV]
    s <- x[, .(np = .N, tc = sum(cov), tM = sum(M)), by = chr_pos]
    acc[s, on = "chr_pos",
        `:=`(n_samp = n_samp + np, total_cov = total_cov + tc, total_M = total_M + tM)]
    done <- c(done, basename(f))
    saveRDS(list(acc = acc, done = done, sig = sig_acc), cachefile)   # resume-safe checkpoint
    rm(x, s); gc(FALSE)
  }
  nsamp <- length(done)
  acc[, `:=`(mean_cov  = total_cov / nsamp,
             mean_meth = fifelse(total_cov > 0, total_M / total_cov, NA_real_))]
  acc[, usable := n_samp >= MINSAMP]                          # >=5x in >=10 samples
  
  ## (3) usable sets
  fg <- acc[is_top == TRUE  & usable]
  bg <- acc[is_top == FALSE & usable & !is.na(mean_meth)]
  cat("usable top:", nrow(fg), "| usable candidate controls:", nrow(bg), "\n")
  
  ## (4) match on (mean_cov, mean_meth): bin both, keep min(#top,#control) per bin => equal N
  cov_breaks  <- unique(quantile(fg$mean_cov, seq(0, 1, 0.1)))
  meth_breaks <- seq(0, 1, by = 0.05)
  cellf <- function(cv, mt) paste(findInterval(cv, cov_breaks,  all.inside = TRUE),
                                  findInterval(mt, meth_breaks, all.inside = TRUE), sep = "_")
  fg[, cell := cellf(mean_cov, mean_meth)]
  bg[, cell := cellf(mean_cov, mean_meth)]
  k <- merge(fg[, .(n_fg = .N), by = cell], bg[, .(n_bg = .N), by = cell], by = "cell", all.x = TRUE)
  k[is.na(n_bg), n_bg := 0L][, keep := pmin(n_fg, n_bg)]
  set.seed(SEED)
  fg <- fg[sample(.N)]; fg[, rk := rowid(cell)]
  bg <- bg[sample(.N)]; bg[, rk := rowid(cell)]
  top_matched  <- fg[k, on = "cell"][rk <= keep]
  ctrl_matched <- bg[k, on = "cell"][rk <= keep]
  
  ## (5) save two equal-size, usable, matched lists
  out <- rbind(top_matched [, .(chr_pos, mean_cov, mean_meth, n_samp, set = "top1pct")],
               ctrl_matched[, .(chr_pos, mean_cov, mean_meth, n_samp, set = "control")])
  saveRDS(out, matchfile)
}
## sanity: the two sets should have near-identical mean coverage & mean methylation
print(out[, .(n = .N, mean_cov = mean(mean_cov), mean_meth = mean(mean_meth)), by = set])
#        set     n mean_cov mean_meth
# 1: top1pct 99644 4.840099 0.5530534
# 2: control 99644 5.000096 0.5540202

#############################################################################################
## PART II — per-sample methylation (beta) matrix for the matched CpGs
##   Re-read the .cov files, this time keeping per-sample beta (M/coverage) at the matched CpGs.
##   Cached because a bit heavy
#############################################################################################
targ  <- unique(out[, .(chr_pos, set)]); setkey(targ, chr_pos)
files <- list.files(covdir, pattern = covpat, full.names = TRUE)

if (file.exists(betalongfile)) {
  betalong <- readRDS(betalongfile)   # already built — reuse
} else {
  betalong <- rbindlist(lapply(files, function(f) {
    x <- fread(f, header = FALSE, select = c(1,2,5,6), col.names = c("chr","start","M","U"))
    x[, chr := refseq2ucsc[as.character(chr)]]
    x <- x[!is.na(chr)]
    x[, chr_pos := paste0(chr, "_", start)]
    x <- x[targ, on = "chr_pos", nomatch = 0L]
    x[, cov := M + U][cov >= MINCOV]
    x[, .(chr_pos, sample = sample_from_file(f), beta = M / cov)]   # per-sample methylation
  }))
  saveRDS(betalong, betalongfile)
}

set_of <- setNames(targ$set, targ$chr_pos) # top1pct / control per CpG

#############################################################################################
## PART III — twin structure + blood cell composition
#############################################################################################
## twin sheet: sample (EGAN) -> pair_id (e.g. MZ1) + zygosity (MZ/DZ), one row per sample
s <- fread(samplecsv)
twins <- unique(s[, .(sample = sample_accession_id, alias = sample_alias)])
twins[, pair_id  := sub(":.*$", "", sub("^WB_", "", alias))]      # "WB_MZ1:1" -> "MZ1"
twins[, zygosity := sub("^([A-Z]{2}).*$", "\\1", pair_id)]        # "MZ1" -> "MZ"
twins <- twins[!zygosity %in% "ST", .(sample, pair_id, zygosity)] # drop the non-twin ST sample
print(twins[, .N, by = pair_id][N != 2])                         # should be empty (each pair = 2)
print(twins[, .N, by = zygosity])                                # MZ / DZ sample counts
# zygosity     N
# 1:       MZ    14
# 2:       DZ    12

## cell proportions (Houseman deconvolution, prepared in script S03). Keep every sample that HAS an estimate:
## drop only too_few_cpgs / collinear_reference (NA estimates); keep high_sampling_error &
## low_cpg_proportion (usable, just noisier) so we don't needlessly lose whole twin pairs.
props_all <- fread(propsfile)
props     <- props_all[!(flag_reason %like% "too_few_cpgs|collinear_reference")]
cat("proportions kept:", nrow(props), "of", nrow(props_all), "samples\n")
# proportions kept: 23 of 27 samples

#############################################################################################
## PART IV — sample / pair accounting: where does N drop?
##   Every inner join below is a potential drop; this shows exactly who is lost and where.
#############################################################################################
samp_B     <- unique(betalong$sample)
samp_twins <- twins$sample
samp_props <- props$sample
cat("beta matrix:", length(samp_B), "| twins:", length(samp_twins),
    "| props(kept):", length(samp_props), "\n")
# beta matrix: 27 | twins: 26 | props(kept): 23 

cat("in cov but not twin sheet :", paste(setdiff(samp_B, samp_twins), collapse = ", "), "\n")
# in cov but not twin sheet : EGAN00001342042 (the ST file, not a twin)

cat("twins with no cov file     :", paste(setdiff(samp_twins, samp_B), collapse = ", "), "\n")
# none

cat("twins dropped by props     :", paste(setdiff(samp_twins, samp_props), collapse = ", "), "\n")
# twins dropped by props     : EGAN00001342026, EGAN00001342038, EGAN00001342039, EGAN00001342040 

keep <- Reduce(intersect, list(samp_B, samp_twins, samp_props))  # samples surviving all three
surv <- twins[sample %in% keep]
print(surv[, .N, by = pair_id][N == 2][, .(complete_pairs = .N), by = sub("[0-9]+$","",pair_id)])
#       sub complete_pairs
# 1:     MZ              6
# 2:     DZ              4

#############################################################################################
## PART V — composition-adjusted beta
##   Assemble beta + set + pair/zygosity + proportions, then regress each CpG's beta across
##   samples on the cell proportions (drop Gran to avoid the sum-to-1 collinearity) and keep
##   the residual + the CpG mean. This removes blood-composition effects from each sample.
#############################################################################################
dt <- as.data.table(betalong)             # one row per (chr_pos, sample): the measured beta
dt[, set := set_of[chr_pos]]              # label each CpG top1pct / control
dt <- merge(dt, twins, by = "sample")     # add pair_id, zygosity for that sample
dt <- merge(dt, props[, c("sample", ..celltypes)], by = "sample")  # add the 6 cell proportions

covn <- setdiff(celltypes, "Gran")        # use 5 of the 6 proportions as predictors
# Why drop one cell type: the six proportions are fractions of a whole, so they (nearly) sum to 1.
# If we put all six plus an intercept into a regression, the intercept equals the sum of the six columns
# = perfect collinearity, and the fit is singular. Dropping one type (Gran here, as the reference)
# removes the redundancy; the intercept absorbs the baseline. 

dt[, beta_adj := {
  P  <- as.matrix(.SD)                       # samples x 5 proportions, for this CpG
  ok <- complete.cases(P) & !is.na(beta)     # samples with full proportions AND a measured beta
  r  <- rep(NA_real_, .N)                    # output vector, one slot per sample-row, prefilled NA
  if (sum(ok) >= length(covn) + 2L) {
    fit <- .lm.fit(cbind(1, P[ok, , drop = FALSE]), beta[ok])
    r[ok] <- fit$residuals + mean(beta[ok])  # add CpG mean back (keep beta scale)
  }
  r
}, by = chr_pos, .SDcols = covn] # by = chr_pos runs the { … } block independently for every CpG

#############################################################################################
## PART VI — Test 1: are top-1% more variable BETWEEN individuals? (hypervariability, out-of-sample)
##   Mean is matched between sets, so the signal is variance. Compare per-CpG variance of the
##   composition-adjusted beta, top vs control.
#############################################################################################
pc <- dt[!is.na(beta_adj),
         .(set = set[1], v = var(beta_adj), mu = mean(beta_adj), n = .N), by = chr_pos][n >= MINSAMP]
print(pc[, .(median_var = median(v), median_mu = median(mu)), by = set])   # mu matched; v should not be
#        set median_var median_mu
# 1: control 0.01462524 0.5692676
# 2: top1pct 0.03806299 0.5666787

wilcox.test(v ~ set, data = pc, alternative = "less") # control < top1pct  (NB direction!)
# Wilcoxon rank sum test with continuity correction
# 
# data:  v by set
# W = 1447233304, p-value < 2.2e-16
# alternative hypothesis: true location shift is less than 0

cat("variance ratio top/control:",
    pc[, median(v), by = set][, V1[set=="top1pct"]/V1[set=="control"]], "\n")
# variance ratio top/control: 2.602555 

#############################################################################################
## PART VII — Test 2: within-pair discordance, and is it the same in MZ and DZ?
##   |Δbeta_adj| within each complete pair, averaged per pair over top and over control CpGs.
##   pair = unit of inference (honest given few pairs). Coverage is matched, so sampling noise
##   cancels in the top-vs-control and MZ-vs-DZ contrasts.
##     ME (stochastic)  : top > control AND MZ excess ~ DZ excess (both > 0)
##     genetic (mQTL)   : DZ excess >> MZ excess (MZ co-twins concordant)
#############################################################################################
wp <- dt[!is.na(beta_adj),
         .(d = if (.N == 2L) abs(diff(beta_adj)) else NA_real_),   # need BOTH twins at this CpG
         by = .(chr_pos, pair_id, zygosity, set)][!is.na(d)]

pp <- dcast(wp[, .(mean_d = mean(d)), by = .(pair_id, zygosity, set)],
            pair_id + zygosity ~ set, value.var = "mean_d")
pp[, excess := top1pct - control]  # extra discordance at top vs control
print(pp)
# Key: <pair_id, zygosity>
#    pair_id zygosity   control   top1pct       excess
# 1:     DZ1       DZ 0.1093775 0.1708473  0.061469806
# 2:     DZ2       DZ 0.1208263 0.1993430  0.078516675
# 3:     DZ3       DZ 0.1550658 0.2348242  0.079758415
# 4:     DZ4       DZ 0.1255058 0.2003904  0.074884580
# 5:     MZ1       MZ 0.1990551 0.1654870 -0.033568078
# 6:     MZ2       MZ 0.1374220 0.1198821 -0.017539835
# 7:     MZ3       MZ 0.1396272 0.1643299  0.024702640
# 8:     MZ4       MZ 0.1415203 0.1763687  0.034848451
# 9:     MZ5       MZ 0.1654368 0.1862900  0.020853201
# 10:    MZ9       MZ 0.1799872 0.1782695 -0.001717655

## (a) are top CpGs more discordant within pairs than matched controls? (paired across pairs)
print(wilcox.test(pp$top1pct, pp$control, paired = TRUE, alternative = "greater"))
# Wilcoxon signed rank exact test
# 
# data:  pp$top1pct and pp$control
# V = 47, p-value = 0.02441
# alternative hypothesis: true location shift is greater than 0

## (b) is that excess the same in MZ and DZ?
print(wilcox.test(excess ~ zygosity, data = pp))
# Wilcoxon rank sum exact test
# 
# data:  excess by zygosity
# W = 24, p-value = 0.009524
# alternative hypothesis: true location shift is not equal to 0

print(pp[, .(mean_top = mean(top1pct), mean_ctrl = mean(control), mean_excess = mean(excess)),
         by = zygosity])
#    zygosity  mean_top mean_ctrl mean_excess
# 1:       DZ 0.2013512 0.1276939 0.073657369
# 2:       MZ 0.1651045 0.1605081 0.004596454

####################
## Visualisation: ##
####################

## ---- Panel A: inter-individual variance, top vs control (Test 1) ----

cols_set <- setNames(c("grey60", "#B2182B"), c("control", "top1pct"))

## counts entering the tests
ncpg  <- pc[, .N, by = set] # CpGs per group (Test 1)
xlab  <- setNames(sprintf("%s\n(n = %s)", ncpg$set, format(ncpg$N, big.mark = ",")), ncpg$set)
npair <- nrow(pp)               # twin pairs (Tests 2a/2b)

## A — inter-individual variance (Test 1)
pA <- ggplot(pc, aes(set, v, fill = set)) +
  geom_violin(colour = NA, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = cols_set, guide = "none") +
  scale_x_discrete(labels = xlab) +
  annotate("text", x = 1.5, y = max(pc$v),
           label = sprintf("%.1f× more variable\n(p < 2.2e-16)",
                           with(pc, median(v[set=="top1pct"])/median(v[set=="control"]))),
           vjust = 1, size = 3.3) +
  labs(x = NULL, y = "per-CpG variance of beta_adj",
       title = "A  Inter-individual variance") +
  theme_minimal(base_size = 12)

## ---- Panel B: per-pair top vs control discordance (Test 2a + the MZ/DZ split) ----
lim <- range(c(pp$control, pp$top1pct))
pB <- ggplot(pp, aes(control, top1pct, colour = zygosity)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 3) +
  coord_equal(xlim = lim, ylim = lim) +
  scale_colour_manual(values = c(MZ = "#0072B2", DZ = "#D55E00"), name = NULL) +
  labs(x = "within-pair discordance, control CpGs",
       y = "within-pair discordance, top-1% CpGs",
       title = "B  Per pair: top vs control") +
  theme_minimal(base_size = 12) + theme(legend.position = c(.15,.9))

## ---- Panel C: excess discordance by zygosity (Test 2b) ----
pC <- ggplot(pp, aes(zygosity, excess, colour = zygosity)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, fatten = 1.5) +
  geom_jitter(width = 0.08, size = 3) +
  scale_colour_manual(values = c(MZ = "#0072B2", DZ = "#D55E00"), guide = "none") +
  annotate("text", x = 1.5, y = max(pp$excess), label = "p = 0.0095", size = 3.3, vjust = 1) +
  labs(x = NULL, y = "excess discordance\n(top-1% − control)",
       title = "C  Excess discordance by zygosity") +
  theme_minimal(base_size = 12)

fig <- pA + pB + pC + plot_layout(widths = c(1, 1.1, 0.9))
ggsave(here("B_MultiTissues/dataOut/figures/twin/twin_hvCpG_validation.pdf"),
       fig, width = 14, height = 4.8)

## Interpretation

# Test 1: hypervariability replicates out-of-sample. At matched mean methylation 
# (0.567 vs 0.569), the top-1% CpGs are 2.6× more variable between individuals 
# than their matched controls (median variance 0.038 vs 0.015, p < 2.2×10⁻¹⁶). 
# Our hypervariability calls, made in multi-tissue discovery data, hold up independently in blood WGBS.

# Test 2: but the variability looks genetic, not metastable. 

# (a) top > control within pairs (paired across 10 pairs, p = 0.024): 
# top-1% CpGs are more discordant within twin pairs than matched controls,
# consistent with them being more variable.

# (b) that excess is driven almost entirely by DZ pairs (MZ vs DZ excess, p = 0.0095): 
#   DZ pairs show a large excess discordance at top-1% (+0.074), MZ pairs essentially none 
# (+0.005). Look at the per-pair numbers — all four DZ pairs are strongly positive (0.061–0.080); 
# the six MZ pairs scatter around zero (−0.034 to +0.035).

#############################################################################################
## Run the twin validation tests (variance + MZ/DZ discordance) on EXTERNAL putative-ME sets
## Prerequisites already in the environment:
##   - from the main script: files, refseq2ucsc, twins, props, celltypes,
##                            MINCOV, MINSAMP, SEED, sample_from_file(), targetfile
##   - from prepPreviousSIV.R (sourced): VanBaakSIV_hg38, VanBaakESS_hg38 (chr_pos vectors, hg38)
#############################################################################################

## Load the set of previously tested MEs & vmeQTL
if (!exists("previousSIVprepared")) {
  source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))}

## Load libraries
if (!exists("libLoaded")) {
  source(here("B_MultiTissues", "quiet_library.R"))}

## Load functions
if (!exists("functionsLoaded")) {
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))}

## background universe to draw matched controls from = the random non-top CpGs already sampled
## in the main target cache (3M genomic CpGs, excludes the top-1%). Cheap, no 20M-table reload.
universe <- readRDS(targetfile)$target[is_top == FALSE, chr_pos]

## ---- one reusable function: everything for a single foreground CpG set ----
run_set_tests <- function(fg_pos, name, bg_mult = 15L, min_fg = 50L) {
  
  ## normalise / clean the foreground ids (expects "chr1_12345"; fix ":"/"-" just in case)
  fg_pos <- unique(na.omit(as.character(fg_pos)))
  fg_pos <- gsub("[:\\-]", "_", fg_pos)
  
  ## target = ME CpGs + a random matched-control candidate pool (not in the ME set)
  set.seed(SEED)
  bgpool <- sample(setdiff(universe, fg_pos), bg_mult * length(fg_pos))
  targ <- data.table(chr_pos = c(fg_pos, bgpool),
                     is_fg   = c(rep(TRUE, length(fg_pos)), rep(FALSE, length(bgpool))))
  setkey(targ, chr_pos)
  
  ## ONE streaming pass over the blood .cov: per-sample coverage + methylation at target CpGs
  bl <- rbindlist(lapply(files, function(f) {
    x <- fread(f, header = FALSE, select = c(1,2,5,6), col.names = c("chr","start","M","U"))
    x[, chr := refseq2ucsc[as.character(chr)]]; x <- x[!is.na(chr)]
    x[, chr_pos := paste0(chr, "_", start)]
    x <- x[targ, on = "chr_pos", nomatch = 0L]
    x[, cov := M + U][cov >= MINCOV]
    x[, .(chr_pos, sample = sample_from_file(f), M, cov, beta = M / cov)]
  }))
  
  ## per-CpG summaries -> usable (>=5x in >=MINSAMP samples)
  summ <- bl[, .(n_samp = .N, mean_cov = mean(cov), mean_meth = sum(M)/sum(cov)), by = chr_pos]
  summ <- merge(summ, targ, by = "chr_pos")
  summ[, usable := n_samp >= MINSAMP]
  fg <- summ[is_fg == TRUE  & usable]
  bg <- summ[is_fg == FALSE & usable & !is.na(mean_meth)]
  cat(sprintf("\n== %s ==\n  foreground: %d given, %d with blood coverage, %d usable\n",
              name, length(fg_pos), summ[is_fg == TRUE, .N], nrow(fg)))
  if (nrow(fg) < min_fg) { message("  too few usable CpGs — skipped"); return(NULL) }
  
  ## match controls on (mean_cov, mean_meth); equal N per cell
  cov_breaks  <- unique(quantile(fg$mean_cov, seq(0, 1, 0.1)))
  meth_breaks <- seq(0, 1, 0.05)
  cf <- function(cv, mt) paste(findInterval(cv, cov_breaks,  all.inside = TRUE),
                               findInterval(mt, meth_breaks, all.inside = TRUE), sep = "_")
  fg[, cell := cf(mean_cov, mean_meth)]; bg[, cell := cf(mean_cov, mean_meth)]
  k <- merge(fg[, .(n_fg = .N), by = cell], bg[, .(n_bg = .N), by = cell], by = "cell", all.x = TRUE)
  k[is.na(n_bg), n_bg := 0L][, keep := pmin(n_fg, n_bg)]
  set.seed(SEED)
  fg <- fg[sample(.N)]; fg[, rk := rowid(cell)]
  bg <- bg[sample(.N)]; bg[, rk := rowid(cell)]
  topm  <- fg[k, on = "cell"][rk <= keep]
  ctrlm <- bg[k, on = "cell"][rk <= keep]
  set_of <- setNames(c(rep("ME", nrow(topm)), rep("control", nrow(ctrlm))),
                     c(topm$chr_pos, ctrlm$chr_pos))
  cat(sprintf("  matched set: %d ME + %d control (cov %.2f/%.2f, meth %.3f/%.3f)\n",
              nrow(topm), nrow(ctrlm), mean(topm$mean_cov), mean(ctrlm$mean_cov),
              mean(topm$mean_meth), mean(ctrlm$mean_meth)))
  
  ## per-sample beta for the matched CpGs, + twin/props, + composition adjustment
  dt <- bl[chr_pos %in% names(set_of)]
  dt[, set := factor(set_of[chr_pos], levels = c("control", "ME"))]
  dt <- merge(dt, twins, by = "sample")
  dt <- merge(dt, props[, c("sample", ..celltypes)], by = "sample")
  covn <- setdiff(celltypes, "Gran")
  dt[, beta_adj := {
    P <- as.matrix(.SD); ok <- complete.cases(P) & !is.na(beta); r <- rep(NA_real_, .N)
    if (sum(ok) >= length(covn) + 2L) {
      fit <- .lm.fit(cbind(1, P[ok, , drop = FALSE]), beta[ok]); r[ok] <- fit$residuals + mean(beta[ok])
    }; r
  }, by = chr_pos, .SDcols = covn]
  
  ## ---- Test 1: inter-individual variance, ME vs matched control ----
  pc <- dt[!is.na(beta_adj),
           .(set = set[1], v = var(beta_adj), mu = mean(beta_adj), n = .N), by = chr_pos][n >= MINSAMP]
  w1  <- wilcox.test(v ~ set, data = pc, alternative = "less")   # control < ME
  vr  <- pc[, median(v), by = set][, V1[set=="ME"]/V1[set=="control"]]
  
  ## ---- Test 2: within-pair discordance, ME vs control, split by zygosity ----
  wp <- dt[!is.na(beta_adj),
           .(d = if (.N == 2L) abs(diff(beta_adj)) else NA_real_),
           by = .(chr_pos, pair_id, zygosity, set)][!is.na(d)]
  pp <- dcast(wp[, .(mean_d = mean(d)), by = .(pair_id, zygosity, set)],
              pair_id + zygosity ~ set, value.var = "mean_d")
  pp[, excess := ME - control]
  t2a <- wilcox.test(pp$ME, pp$control, paired = TRUE, alternative = "greater")  # ME > control
  t2b <- wilcox.test(excess ~ zygosity, data = pp)                               # MZ vs DZ excess
  
  ## tidy summary row
  zyg <- pp[, .(mean_excess = mean(excess)), by = zygosity]
  summary_row <- data.table(
    set = name, n_matched = nrow(topm),
    var_ratio = vr, var_p = w1$p.value,
    p_ME_gt_control = t2a$p.value,
    excess_MZ = zyg[zygosity=="MZ", mean_excess],
    excess_DZ = zyg[zygosity=="DZ", mean_excess],
    p_MZ_vs_DZ = t2b$p.value)
  print(summary_row)
  list(name = name, pc = pc, pp = pp, summary = summary_row)
}

## ---- run for the Van Baak sets ----
res_VanBaakSIV <- run_set_tests(VanBaakSIV_hg38, "VanBaak SIV")
res_VanBaakESS <- run_set_tests(VanBaakESS_hg38, "VanBaak ESS")

## combined summary table across sets (add your own top-1% row for comparison)
res_all <- rbindlist(lapply(list(res_VanBaakSIV, res_VanBaakESS), function(r) if (!is.null(r)) r$summary),
                     fill = TRUE)
print(res_all)
#            set n_matched var_ratio      var_p p_ME_gt_control   excess_MZ     excess_DZ p_MZ_vs_DZ
# 1: VanBaak SIV        54  1.653965 0.02505658      0.46093750 -0.01900826  0.0284081555 0.25714286
# 2: VanBaak ESS       283  1.136995 0.01595351      0.01855469  0.01893762 -0.0000491144 0.01904762
saveRDS(list(SIV = res_VanBaakSIV, ESS = res_VanBaakESS), here("gitignore/prevME_twin_tests.rds"))

###################
## Visualisation ##
###################

plot_set_tests <- function(res, outfile = NULL) {
  pc <- copy(res$pc); pp <- copy(res$pp); s <- res$summary; nm <- res$name
  
  ## nice labels: "control" vs the set name instead of "ME"
  pc[, grp := factor(ifelse(set == "ME", nm, "control"), levels = c("control", nm))]
  cols_grp <- setNames(c("grey60", "#B2182B"), c("control", nm))
  cols_zyg <- c(MZ = "#0072B2", DZ = "#D55E00")
  
  ## counts entering the tests
  ncpg  <- pc[, .N, by = grp]                                   # CpGs per group (Test 1)
  xlab  <- setNames(sprintf("%s\n(n = %s)", ncpg$grp, format(ncpg$N, big.mark = ",")), ncpg$grp)
  npair <- nrow(pp)                                             # twin pairs (Tests 2a/2b)
  
  ## A — inter-individual variance (Test 1)
  pA <- ggplot(pc, aes(grp, v, fill = grp)) +
    geom_violin(colour = NA, alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    scale_x_discrete(labels = xlab) +                           # <- n under each group
    scale_fill_manual(values = cols_grp, guide = "none") +
    annotate("text", x = 1.5, y = max(pc$v),
             label = sprintf("%.1f× more variable\n(p = %.2g)", s$var_ratio, s$var_p),
             vjust = 1, size = 3.2) +
    labs(x = NULL, y = "per-CpG variance of beta_adj", title = "A  Inter-individual variance") +
    theme_minimal(base_size = 12)
  
  ## B — per-pair control vs ME discordance (Test 2a + zygosity split)
  lim <- range(c(pp$control, pp$ME))
  pB <- ggplot(pp, aes(control, ME, colour = zygosity)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 3) + coord_equal(xlim = lim, ylim = lim) +
    scale_colour_manual(values = cols_zyg, name = NULL) +
    labs(x = "within-pair discordance, control",
         y = paste0("within-pair discordance, ", nm),
         title = "B  Per pair: ME vs control") +
    theme_minimal(base_size = 12) + theme(legend.position = c(.15, .9))
  
  ## C — excess discordance by zygosity (Test 2b)
  pC <- ggplot(pp, aes(zygosity, excess, colour = zygosity)) +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, fatten = 1.5) +
    geom_jitter(width = 0.08, size = 3) +
    scale_colour_manual(values = cols_zyg, guide = "none") +
    annotate("text", x = 1.5, y = max(pp$excess), label = sprintf("p = %.3g", s$p_MZ_vs_DZ),
             size = 3.2, vjust = 1) +
    labs(x = NULL, y = "excess discordance (ME − control)",
         title = "C  Excess discordance by zygosity") +
    theme_minimal(base_size = 12)
  
  fig <- pA + pB + pC + plot_layout(widths = c(1, 1.1, 0.9)) +
    patchwork::plot_annotation(
      title = nm,
      subtitle = sprintf("%s matched CpGs per group · %d twin pairs (%d MZ, %d DZ)",
                         format(ncpg$N[1], big.mark = ","), npair,
                         pp[zygosity == "MZ", .N], pp[zygosity == "DZ", .N]))
  
  if (!is.null(outfile)) ggsave(outfile, fig, width = 14, height = 4.8)
  fig
}

## one figure per set
plot_set_tests(res_VanBaakSIV, here("B_MultiTissues/dataOut/figures/twin/twin_VanBaakSIV.pdf"))
plot_set_tests(res_VanBaakESS, here("B_MultiTissues/dataOut/figures/twin/twin_VanBaakESS.pdf"))
