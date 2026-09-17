library(data.table)

##############################################################################
## Step 1 — extract target CpGs from coverage files (done in S08.1 in bash) ##
##############################################################################

meth <- fread("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/twins_endo_target_meth_all.tsv")
meth[, beta := meth_pct / 100]

#########################################################################
## Step 2 — R: attach twin structure and build within-pair differences ##
#########################################################################

# --- map sample_accession_id -> zygosity + pair + twin from the sample file ---
## copied from /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/Busche2015/sample_file_adiposeWGBS.csv
samp <- fread("/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataIn/sample_file_adiposeWGBS.csv")
# alias like "AT_MZ1:1" -> zygosity=MZ, pair=1, twin=1
samp[, alias := sub("^AT_", "", sample_alias)]                       # "MZ1:1"
samp[, `:=`(zygosity = sub("([A-Z]+).*", "\\1", alias),              # MZ / DZ
            pair     = sub("[A-Z]+([0-9]+):.*", "\\1", alias),       # 1
            twin     = sub(".*:", "", alias))]                       # 1 / 2
samp[, pair_id := paste0(zygosity, pair)]                            # MZ1, DZ5
meta <- unique(samp[, .(sample = sample_accession_id, zygosity, pair_id, twin)])

meth <- merge(meth, meta, by = "sample")

# --- one beta per (CpG, sample); pivot twins side by side per pair ---
# keep only CpGs with BOTH twins of a pair covered (else no within-pair diff)
wide <- dcast(meth, chr_pos + pair_id + zygosity ~ twin, value.var = "beta")
setnames(wide, c("1","2"), c("twinA","twinB"), skip_absent = TRUE)
wide <- wide[!is.na(twinA) & !is.na(twinB)]

# within-pair absolute methylation difference
wide[, abs_diff := abs(twinA - twinB)]

#############################################################################
## Step 3 — the test: MZ vs DZ within-pair difference, per CpG and overall ##
#############################################################################

## overall (pooled across CpGs) — more power than per-CpG at your N
mz <- wide[zygosity == "MZ", abs_diff]
dz <- wide[zygosity == "DZ", abs_diff]

wilcox.test(mz, dz, alternative = "less")   # H1: MZ diff < DZ diff (genetic effect)
data.table(group = c("MZ","DZ"),
           n = c(length(mz), length(dz)),
           median_absdiff = c(median(mz), median(dz)),
           mean_absdiff   = c(mean(mz),  mean(dz)))

## per-CpG (if you want each target scored)
per_cpg <- wide[, {
  m <- abs_diff[zygosity=="MZ"]; d <- abs_diff[zygosity=="DZ"]
  if (length(m) >= 3 && length(d) >= 3)
    .(n_mz=length(m), n_dz=length(d),
      med_mz=median(m), med_dz=median(d),
      p = wilcox.test(m, d, alternative="less")$p.value)
  else .(n_mz=length(m), n_dz=length(d), med_mz=median(m), med_dz=median(d), p=NA_real_)
}, by = chr_pos]
print(per_cpg[order(p)])
