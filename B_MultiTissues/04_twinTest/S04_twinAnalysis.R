library(here); library(data.table)

#############################################################################################
## PART I. Select top 1% hvCpG and same number of random control matching methylation mean ##
## accross samples. Min cov = 5X, in at least 10 twins.
## Output saved in here("gitignore/matched_CpGs_top_vs_control.rds")
#############################################################################################
# Twin methylation extracted in cluster in: S01_alignCallMeth.sh
## After S01 methyl call has ran, in cluster run the following to select only 5X
## CpG positions:
# /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/Busche2015/01Methcall/Blood
# for f in *.bismark.cov.gz; do
# base=${f%.bismark.cov.gz}
# zcat "$f" | awk -F'\t' '($5+$6) >= 5' > "${base}_5X.cov"
# done

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
BG_MULT <- 15L     # candidate controls drawn = BG_MULT x (#top CpGs) (7 was too small)
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
cat("target CpGs:", nrow(target), "(top:", n_top, ")")
# target CpGs: 3239472 (top: 202467 )

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
cat(length(newfiles), "new /", length(files), "files")
# 27 new / 27 files

for (f in newfiles) {
  print(f)
  x <- fread(f, header = FALSE, select = c(1, 2, 5, 6),
             col.names = c("chr", "start", "M", "U"))
  x[, chr := refseq2ucsc[as.character(chr)]] # correct chromosome name
  x <- x[!is.na(chr)] # remove extra chromosomes
  x[, chr_pos := paste0(chr, "_", start)] # same format in both
  x <- x[target, on = "chr_pos", nomatch = 0L] # keep top + candidate controls only
  x[, cov := M + U] # calculate coverage
  x <- x[cov >= MINCOV] # safety step (redundant on 5X files, safe on full)
  ## calc total accross samples for later mean calculation by sample
  s <- x[, .(np = .N, tc = sum(cov), tM = sum(M)), by = chr_pos]
  acc[s, on = "chr_pos",
      `:=`(n_samp = n_samp + np, total_cov = total_cov + tc, total_M = total_M + tM)]
  done <- c(done, basename(f))
  saveRDS(list(acc = acc, done = done, sig = sig_acc), cachefile)   # resume-safe checkpoint
  rm(x, s); gc(FALSE)
}

nsamp <- length(done) # 27
acc[, `:=`(mean_cov  = total_cov / nsamp,  # mean depth (absent = 0)
           mean_meth = fifelse(total_cov > 0, total_M / total_cov, NA_real_))]  # pooled beta
acc[, usable := n_samp >= MINSAMP] # >=5x in >=10 samples

## ================================================================
## (3) usable sets for each category
## ================================================================
fg <- acc[is_top == TRUE  & usable]                      # usable top-1% CpGs
bg <- acc[is_top == FALSE & usable & !is.na(mean_meth)]  # usable candidate controls
cat("usable top:", nrow(fg), " | usable candidate controls:", nrow(bg), "\n")
# usable top: 99644  | usable candidate controls: 1495851 

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
# final matched set size (each): 99644  | top dropped for lack of in-cell control: 0 

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
# set     n mean_cov mean_meth
# <char> <int>    <num>     <num>
# 1: top1pct 99644 4.840099 0.5530534
# 2: control 99644 5.000096 0.5540202

## cells where the control pool ran short (reduce the final n)
print(k[n_bg < n_fg][order(-(n_fg - n_bg))])
# Key: <cell>
#   Empty data.table (0 rows and 4 cols): cell,n_fg,n_bg,keep

#######################################################################################
## PART II. Comparison of top1pct and controls for (1) methylation and (2) hv scores ## 
#######################################################################################

######################################
## A. Correct for blood composition ##
######################################
## Done in bash script S03_bloodComposition.sh

######################################
## B. Correct for blood composition ##
######################################

chr_pos mean_cov  mean_meth n_samp     set
<char>    <num>      <num>  <int>  <char>
  1:  chr2_137865182 8.148148 0.03636364     22 top1pct
2:   chr5_55720156 7.111111 0.04687500     18 top1pct
3:  chr4_174840068 7.259259 0.03061224     20 top1pct
4:   chr1_85991691 8.185185 0.03619910     20 top1pct
