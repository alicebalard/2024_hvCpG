# add to the S01
# 
# cd $METHCALL_OUTDIR
# instead of 
# 
# cd 01Methcall/Blood
# 
# for f in *.bismark.cov.gz; do
# base=${f%.bismark.cov.gz}
# zcat "$f" | awk -F'\t' '($5+$6) >= 5' > "${base}_5X.cov"
# done
# 
# 
# 
# DATADIR="/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/Busche2015"
# EGAN00001342042 still running
# 
# 
# I have data for extracting twins (mono and dizygotic, MZ and DZ) in script
# S01_alignCallMeth.sh. It is done.
# 
# I also have a table saved as a RDS object in:
#   
#   /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/table3layers_coveredIn3dt.RDS
# that contains when opened in R:
#   
#   Index: <chr_pos>
#   chr   pos   end width strand logBF_per_ds_endo logBF_per_ds_ecto logBF_per_ds_meso logBF_per_ds_allLayers
# <fctr> <int> <int> <int> <fctr>             <num>             <num>             <num>                  <num>
#   1:   chr1 17452 17452     1      *       -0.04351713        -0.5127993       -0.09070382             -0.1152128
# 2:   chr1 17478 17478     1      *       -0.66693206        -0.7742006       -0.75360608             -0.7150446
# 3:   chr1 17483 17483     1      *       -0.77952104        -0.7235260       -0.65861379             -0.7208407
# 4:   chr1 17492 17492     1      *       -0.78485244        -0.7959979       -0.77317453             -0.7807970
# 5:   chr1 17562 17562     1      *       -0.67376480        -0.1162556       -0.73765576             -0.6249438
# ---                                                                                                             
#   20246675:   chrM 16427 16427     1      *       -0.76611492        -0.8092041       -0.78422530             -0.7777684
# 20246676:   chrM 16449 16449     1      *       -0.77896241        -0.8102262       -0.78454732             -0.7842658
# 20246677:   chrM 16454 16454     1      *       -0.76012998        -0.8100064       -0.78424614             -0.7750473
# 20246678:   chrM 16495 16495     1      *       -0.79843223        -0.8107169       -0.78475610             -0.7940611
# 20246679:   chrM 16542 16542     1      *       -0.79789180        -0.8098529       -0.78450459             -0.7935754
# logBF_per_ds_endo6gp logBF_per_ds_meso6gp    chr_pos logBF_per_ds_mean3layers n_layers percentile_logBF_per_ds_allLayers
# <num>                <num>     <char>                    <num>    <num>                             <num>
#   1:         -0.090381521            0.7125974 chr1_17452               -0.2156734        3                      68.707598910
# 2:         -0.209438572           -0.2523771 chr1_17478               -0.7315796        3                       5.926571958
# 3:          0.211521887           -0.3124793 chr1_17483               -0.7205536        3                       5.087940595
# 4:         -0.794229508           -0.3687644 chr1_17492               -0.7846750        3                       0.031269326
# 5:         -0.391897218           -0.7440644 chr1_17562               -0.5092254        3                      20.162447382
# ---                                                                                                                         
#   20246675:          0.022283014            0.4930459 chrM_16427               -0.7865148        3                       0.069700320
# 20246676:          0.558905432           -0.3465653 chrM_16449               -0.7912453        3                       0.011231472
# 20246677:         -0.413283469           -0.5471890 chrM_16454               -0.7847942        3                       0.125087181
# 20246678:          0.545002981           -0.4953212 chrM_16495               -0.7979684        3                       0.001052024
# 20246679:         -0.005184287           -0.6119597 chrM_16542               -0.7974164        3                       0.002276917
# cum_offset  group       pos2
# <num> <char>      <num>
#   1:          0   <NA>      17452
# 2:          0   <NA>      17478
# 3:          0   <NA>      17483
# 4:          0   <NA>      17492
# 5:          0   <NA>      17562
# ---                             
#   20246675: 3086384434   <NA> 3086400861
# 20246676: 3086384434   <NA> 3086400883
# 20246677: 3086384434   <NA> 3086400888
# 20246678: 3086384434   <NA> 3086400929
# 20246679: 3086384434   <NA> 3086400976
# 
# From there, I would like to:
#   
#   (1) Select the top 1% chr_pos on percentile_logBF_per_ds_allLayers with this:
#   table3layers_coveredIn3dt$chr_pos[
#     table3layers_coveredIn3dt$percentile_logBF_per_ds_allLayers >= 99] (total is 202467)
# 
# (2) find how many of them have a sufficient coverage in twins for blood cells (find the average and median coverage, do a histogram)
# 
# (2) extract methylation for all twins of the top 1% using S02_extractTargetMeth.sh
# 
# (3) extract methylation for all twins of a random set of chr_pos contained in table3layers_coveredIn3dt$chr_pos, maybe 500 000 for a start. Then I want to select 
# 
# 
# EGAN00001342041_R1_val_1_bismark_bt2_pe







library(here); library(data.table)

## ---------------- parameters ----------------
tabRDS     <- here("gitignore/table3layers_coveredIn3dt.RDS")
covdir     <- here("gitignore/Blood")
covpat     <- "_5X\\.cov(\\.gz)?$"                 # the 5x-prefiltered blood files
targetfile <- here("gitignore/match_target.rds")   # top + candidate pool (reproducible)
cachefile  <- here("gitignore/match_acc.rds")       # streamed cov+meth accumulator

refseq2ucsc <- c(
  NC_000001.11="chr1",  NC_000002.12="chr2",  NC_000003.12="chr3",
  NC_000004.12="chr4",  NC_000005.10="chr5",  NC_000006.12="chr6",
  NC_000007.14="chr7",  NC_000008.11="chr8",  NC_000009.12="chr9",
  NC_000010.11="chr10", NC_000011.10="chr11", NC_000012.12="chr12",
  NC_000013.11="chr13", NC_000014.9="chr14",  NC_000015.10="chr15",
  NC_000016.10="chr16", NC_000017.11="chr17", NC_000018.10="chr18",
  NC_000019.10="chr19", NC_000020.11="chr20", NC_000021.9="chr21",
  NC_000022.11="chr22", NC_000023.11="chrX",  NC_000024.10="chrY",
  NC_012920.1="chrM")

MINCOV  <- 5L     # a CpG is "covered" in a sample at >= this depth
MINSAMP <- 10L    # "usable" = covered in >= this many samples
BG_MULT <- 3L     # candidate controls drawn = BG_MULT x (#top CpGs)
SEED    <- 1234

## ================================================================
## (1) TARGET = top-1% hvCpGs + 3x random non-top CpGs (load big RDS once)
## ================================================================
sig <- list(mincov = MINCOV, bg_mult = BG_MULT, seed = SEED)
if (file.exists(targetfile) && identical(readRDS(targetfile)$sig, sig)) {
  tb <- readRDS(targetfile); target <- tb$target; n_top <- tb$n_top
} else {
  dt      <- readRDS(tabRDS)[, .(chr_pos, pct = percentile_logBF_per_ds_allLayers)]
  top_pos <- dt[pct >= 99, chr_pos]                      # top 1%
  n_top   <- length(top_pos)
  set.seed(SEED)
  bg_pos  <- sample(dt[pct < 99, chr_pos], BG_MULT * n_top)   # 3x candidate controls
  target  <- data.table(chr_pos = c(top_pos, bg_pos),
                        is_top  = c(rep(TRUE, n_top), rep(FALSE, length(bg_pos))))
  setkey(target, chr_pos)
  saveRDS(list(target = target, n_top = n_top, sig = sig), targetfile)
  rm(dt, top_pos, bg_pos); gc()
}
cat("target CpGs:", nrow(target), "(top:", n_top, ")\n")
# target CpGs: 809868 (top: 202467 )

## ================================================================
## (2) ONE streaming pass: coverage + methylation for every target CpG
## ================================================================
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
cat(length(newfiles), "new /", length(files), "files\n")

for (f in newfiles) {
  print(f)
  x <- fread(f, header = FALSE, select = c(1, 2, 5, 6),
             col.names = c("chr", "start", "M", "U"))
  x[, chr := refseq2ucsc[as.character(chr)]]
  x <- x[!is.na(chr)]
  x[, chr_pos := paste0(chr, "_", start)]
  x <- x[target, on = "chr_pos", nomatch = 0L]        # keep top + candidate controls only
  x[, cov := M + U]
  x <- x[cov >= MINCOV]                               # >=5x (redundant on 5X files, safe on full)
  s <- x[, .(np = .N, tc = sum(cov), tM = sum(M)), by = chr_pos]
  acc[s, on = "chr_pos",
      `:=`(n_samp = n_samp + np, total_cov = total_cov + tc, total_M = total_M + tM)]
  done <- c(done, basename(f))
  saveRDS(list(acc = acc, done = done, sig = sig_acc), cachefile)   # resume-safe checkpoint
  rm(x, s); gc(FALSE)
}

nsamp <- length(done)
acc[, `:=`(mean_cov  = total_cov / nsamp,                            # mean depth (absent = 0)
           mean_meth = fifelse(total_cov > 0, total_M / total_cov, NA_real_))]  # pooled beta
acc[, usable := n_samp >= MINSAMP]                                   # >=5x in >=10 samples

## ================================================================
## (3) usable sets for each category
## ================================================================
fg <- acc[is_top == TRUE  & usable]                      # usable top-1% CpGs
bg <- acc[is_top == FALSE & usable & !is.na(mean_meth)]  # usable candidate controls
cat("usable top:", nrow(fg), " | usable candidate controls:", nrow(bg), "\n")

## ================================================================
## (4) match on (mean_cov, mean_meth); keep EQUAL numbers in both
## ================================================================
cov_breaks  <- unique(quantile(fg$mean_cov, seq(0, 1, 0.1)))  # 10 coverage bands
meth_breaks <- seq(0, 1, by = 0.05)                           # 20 methylation bands
cellf <- function(cv, mt) paste(findInterval(cv, cov_breaks,  all.inside = TRUE),
                                findInterval(mt, meth_breaks, all.inside = TRUE), sep = "_")
fg[, cell := cellf(mean_cov, mean_meth)]
bg[, cell := cellf(mean_cov, mean_meth)]

## per cell, we can form min(#top, #control) matched pairs
k <- merge(fg[, .(n_fg = .N), by = cell],
           bg[, .(n_bg = .N), by = cell], by = "cell", all.x = TRUE)
k[is.na(n_bg), n_bg := 0L]
k[, keep := pmin(n_fg, n_bg)]

set.seed(SEED)
fg <- fg[sample(.N)]; fg[, rk := rowid(cell)]            # random order within each cell
bg <- bg[sample(.N)]; bg[, rk := rowid(cell)]
top_matched  <- fg[k, on = "cell"][rk <= keep]           # keep `keep` top per cell
ctrl_matched <- bg[k, on = "cell"][rk <= keep]           # keep `keep` controls per cell

cat("final matched set size (each):", nrow(top_matched),
    " | top dropped for lack of in-cell control:", nrow(fg) - nrow(top_matched), "\n")

## ================================================================
## (5) outputs: two equal-size, usable, matched CpG lists
## ================================================================
usable_top     <- top_matched [, .(chr_pos, mean_cov, mean_meth, n_samp, set = "top1pct")]
usable_control <- ctrl_matched[, .(chr_pos, mean_cov, mean_meth, n_samp, set = "control")]
stopifnot(nrow(usable_top) == nrow(usable_control))

out <- rbind(usable_top, usable_control)
saveRDS(out, here("gitignore/matched_CpGs_top_vs_control.rds"))

## checks: the two groups should have near-identical cov + meth
print(out[, .(n = .N, mean_cov = mean(mean_cov), mean_meth = mean(mean_meth)), by = set])
## cells where the control pool ran short (reduce the final n)
print(k[n_bg < n_fg][order(-(n_fg - n_bg))])