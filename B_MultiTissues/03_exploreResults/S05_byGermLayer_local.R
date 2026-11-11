#####################################################################
# S05 — Per-germ-layer hypervariability categories
#####################################################################

#####################################################################
## Prepare
library(here)
## Load libraries
source(here("B_MultiTissues", "quiet_library.R"))

## Load functions
if (!exists("functionsLoaded")) {
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))}

## Add previous MEs including Maria's results
## Load the set of previously tested MEs & vmeQTL
# if (!exists("previousSIVprepared")) {
# source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))}
#####################################################################

variant  <- "SNP_SDASMrm"
prep_dir <- here("gitignore/resultsAtlasPrepared", variant)

## Object created in S04 (GRanges of CpGs covered in all 3 layers, with per-layer scores)
load(here("gitignore/table3layers_coveredIn3_26_08_26.Rda"))
dt <- as.data.table(table3layers_coveredIn3)           # one row per covered-in-3 CpG
setnames(dt, "chr_pos", "chr_pos", skip_absent = TRUE) # ensure chr_pos exists

# category_colours <- c(
#   "top1pcEctoOnly" = paletteer::paletteer_d("nationalparkcolors::Everglades")[1],
#   "top1pcEndoOnly" = paletteer::paletteer_d("nationalparkcolors::Everglades")[2],
#   "top1pcMesoOnly" = paletteer::paletteer_d("nationalparkcolors::Everglades")[3],
#   "top1pc3layers" = "orange",
#   "constitutive" = "grey80",
#   "background" = "black"
# )

layers_colours <- c(
  "ecto" = paletteer::paletteer_d("nationalparkcolors::Everglades")[1],
  "endo" = paletteer::paletteer_d("nationalparkcolors::Everglades")[2],
  "meso" = paletteer::paletteer_d("nationalparkcolors::Everglades")[3]
)

## ==============================================================================
## Definition of layer-specific HV CpGs
## Top 1% in one layer and bottom 50% in the other two layers
## ==============================================================================

q_top <- function(x) {
  quantile(x, probs = 0.99, na.rm = TRUE, names = FALSE)
}

q_bot <- function(x) {
  quantile(x, probs = 0.50, na.rm = TRUE, names = FALSE)
}

top_endo  <- q_top(dt$logBF_per_ds_endo)
top_meso  <- q_top(dt$logBF_per_ds_meso)
top_ecto  <- q_top(dt$logBF_per_ds_ecto)

bot_endo  <- q_bot(dt$logBF_per_ds_endo)
bot_meso  <- q_bot(dt$logBF_per_ds_meso)
bot_ecto  <- q_bot(dt$logBF_per_ds_ecto)

## Make category for every CpG
dt[, `:=`(
  hv_none =
    logBF_per_ds_endo <= bot_endo &
    logBF_per_ds_meso <= bot_meso &
    logBF_per_ds_ecto <= bot_ecto,
  hv_endo_only =
    logBF_per_ds_endo >= top_endo &
    logBF_per_ds_meso <= bot_meso &
    logBF_per_ds_ecto <= bot_ecto,
  hv_meso_only =
    logBF_per_ds_meso >= top_meso &
    logBF_per_ds_endo <= bot_endo &
    logBF_per_ds_ecto <= bot_ecto,
  hv_ecto_only =
    logBF_per_ds_ecto >= top_ecto &
    logBF_per_ds_endo <= bot_endo &
    logBF_per_ds_meso <= bot_meso,
  hv_endomeso =
    logBF_per_ds_endo >= top_endo &
    logBF_per_ds_meso >= top_meso &
    logBF_per_ds_ecto <= bot_ecto,
  hv_endoecto = 
    logBF_per_ds_endo >= top_endo &
    logBF_per_ds_meso <= bot_meso &
    logBF_per_ds_ecto >= top_ecto,
  hv_mesoecto = 
    logBF_per_ds_endo <= bot_endo &
    logBF_per_ds_meso >= top_meso &
    logBF_per_ds_ecto >= top_ecto,
  hv_mesoectoendo = ## // MEs
    logBF_per_ds_endo >= top_endo &
    logBF_per_ds_meso >= top_meso &
    logBF_per_ds_ecto >= top_ecto
)]

hv_counts <- data.table(
  category = c(
    "In no layer",
    "Mesoderm only",
    "Endoderm only",
    "Ectoderm only",
    "Endoderm + mesoderm",
    "Endoderm + ectoderm",
    "Mesoderm + ectoderm",
    "Endoderm + mesoderm + ectoderm"
  ),
  variable = c(
    "hv_none",
    "hv_meso_only",
    "hv_endo_only",
    "hv_ecto_only",
    "hv_endomeso",
    "hv_endoecto",
    "hv_mesoecto",
    "hv_mesoectoendo"
  ),
  n = c(
    sum(dt$hv_none, na.rm = TRUE),
    sum(dt$hv_meso_only, na.rm = TRUE),
    sum(dt$hv_endo_only, na.rm = TRUE),
    sum(dt$hv_ecto_only, na.rm = TRUE),
    sum(dt$hv_endomeso, na.rm = TRUE),
    sum(dt$hv_endoecto, na.rm = TRUE),
    sum(dt$hv_mesoecto, na.rm = TRUE),
    sum(dt$hv_mesoectoendo, na.rm = TRUE)
  )
)

hv_counts
#                          category        variable       n
# 1:                    In no layer         hv_none 6723418
# 2:                  Mesoderm only    hv_meso_only      63
# 3:                  Endoderm only    hv_endo_only      18
# 4:                  Ectoderm only    hv_ecto_only     297
# 5:            Endoderm + mesoderm     hv_endomeso      39
# 6:            Endoderm + ectoderm     hv_endoecto      18
# 7:            Mesoderm + ectoderm     hv_mesoecto       3
# 8: Endoderm + mesoderm + ectoderm hv_mesoectoendo   60424 

## =============================================================================
## Test overlap between the top in each layer
## =============================================================================

if (!file.exists(here("B_MultiTissues/dataOut/figures/script05/tophvvsbottombylayer.pdf"))){
  
  ## 1. Pairwise scatter plots — shared fill scale + one common legend
  pairs <- list(c("logBF_per_ds_endo","logBF_per_ds_meso","Hypervariability score endoderm","Hypervariability score mesoderm"),
                c("logBF_per_ds_endo","logBF_per_ds_ecto","Hypervariability score endoderm","Hypervariability score ectoderm"),
                c("logBF_per_ds_meso","logBF_per_ds_ecto","Hypervariability score mesoderm","Hypervariability score ectoderm"))
  
  make_scatter <- function(p) {
    xcol <- p[1]; ycol <- p[2]; xlab <- p[3]; ylab <- p[4]
    x <- dt[[xcol]]; y <- dt[[ycol]]
    ok <- is.finite(x) & is.finite(y)
    ct <- cor.test(x[ok], y[ok], method = "spearman")
    lab <- sprintf("Spearman rho = %.2f\n%s\nn = %s",
                   unname(ct$estimate),
                   if (ct$p.value < 0.001) "p < 0.001" else sprintf("p = %.2g", ct$p.value),
                   format(sum(ok), big.mark = ","))
    ggplot(dt[ok, ], aes(x = .data[[xcol]], y = .data[[ycol]])) +
      geom_hex(bins = 120) +                                  # layer 1 -> counts read from here
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
      geom_smooth(method = "lm", se = FALSE,
                  colour = "#DC3220", linewidth = 0.8) +      # correlation line
      annotate("label", x = -Inf, y = Inf, hjust = 0, vjust = 1,
               label = lab, size = 3.2, label.size = NA,
               fill = "white", alpha = 0.7) +                 # stats box (readable on viridis)
      labs(x = xlab, y = ylab) +
      theme_minimal(base_size = 12)
  }
  
  sc <- lapply(pairs, make_scatter)              # build the 3 panels ONCE
  
  ## shared colour limits, read off the SAME objects (no second round of binning)
  counts <- unlist(lapply(sc, function(g) ggplot_build(g)$data[[1]]$count))
  lims   <- c(1, max(counts, na.rm = TRUE))
  
  sc <- lapply(sc, function(g)
    g + scale_fill_viridis_c(trans = "log10", name = "CpGs", limits = lims))
  
  legend <- cowplot::get_legend(sc[[1]] + theme(legend.position = "right"))
  sc     <- lapply(sc, function(g) g + theme(legend.position = "none"))
  scatter_row  <- cowplot::plot_grid(plotlist = sc, labels = c("A", "B","C"), nrow = 1)
  scatterplots <- cowplot::plot_grid(scatter_row, legend, nrow = 1, rel_widths = c(1, 0.12))
  
  ## 2. for the CpGs that are top-1% in layer X, where do they sit in the other layers' rank distributions?
  
  # within-layer percentile rank (0–1), comparable across layers
  pr <- function(x) data.table::frank(x, na.last = "keep", ties.method = "average") /
    sum(!is.na(x))
  dt[, `:=`(pr_endo = pr(logBF_per_ds_endo),
            pr_meso = pr(logBF_per_ds_meso),
            pr_ecto = pr(logBF_per_ds_ecto))]
  
  layers <- c("endo","meso","ecto")
  
  makeridgeplotGermLayers <- function(thr){
    long <- rbindlist(lapply(layers, function(f) {
      foc <- dt[get(paste0("pr_", f)) >= 1-thr/100] # top thr% in focal layer
      rbindlist(lapply(setdiff(layers, f), function(o)
        data.table(focal = f, other = o,
                   pr_other = foc[[paste0("pr_", o)]])))
    }))
    
    long[, panel := paste0("top", thr, "% ", toupper(focal), " → rank in ", toupper(other))]
    
    ggplot(long, aes(pr_other, panel, fill = focal)) +
      geom_density_ridges(scale = 1.3, alpha = .7, rel_min_height = 0.01,
                          quantile_lines = F) +
      geom_vline(xintercept = 0.99, linetype = 2, colour = "grey40") +
      scale_x_continuous("Percentile rank in the OTHER layer", labels = scales::percent) +
      scale_fill_manual(values = layers_colours)+
      theme_ridges(font_size = 14)+
      theme(legend.position = "none", axis.title.y = element_blank())
    ## title: Where each layer's top", thr, "% hvCpGs fall in the other layers
  }
  
  p1 <- makeridgeplotGermLayers(1)
  p2 <- makeridgeplotGermLayers(5)
  p3 <- makeridgeplotGermLayers(10)
  
  ridgeplots <- cowplot::plot_grid(p1, p2, p3,
                                   labels = c("D", "E", "F"), nrow = 1)
  
  ggsave(here("B_MultiTissues/dataOut/figures/script05/tophvvsbottombylayer.pdf"),
         cowplot::plot_grid(scatterplots, ridgeplots, nrow = 2),
         width = 16, height = 8, dpi = 300, bg = "white")
}

################################################################################
## PART A — CpG lists per category (+ controls) for raw-methylation extraction ##
################################################################################
set.seed(1234)

## the six layer-restricted categories are small -> take ALL their CpGs
cat_cpgs <- list(
  meso_only = dt[hv_meso_only == TRUE, chr_pos],   # hv in Meso only
  endo_only = dt[hv_endo_only == TRUE, chr_pos],   # hv in Endo only
  ecto_only = dt[hv_ecto_only == TRUE, chr_pos],   # hv in Ecto only
  endomeso  = dt[hv_endomeso  == TRUE, chr_pos],   # hv in Endo + Meso (ingression branch)
  endoecto  = dt[hv_endoecto  == TRUE, chr_pos],   # hv in Endo + Ecto
  mesoecto  = dt[hv_mesoecto  == TRUE, chr_pos]    # hv in Meso + Ecto
)
## controls: 1000 random CpGs from each non-layer-specific class
cat_cpgs$control_none    <- sample(dt[hv_none == TRUE,        chr_pos], 1000)  # neg ctrl (low hv everywhere)
cat_cpgs$control_3layers <- sample(dt[hv_mesoectoendo == TRUE, chr_pos], 1000) # pos ctrl (systemic hv = MEs)

stopifnot(!any(duplicated(unlist(cat_cpgs))))                 # each CpG in one list only

## lookup table + the union to extract
cpg_category <- rbindlist(lapply(names(cat_cpgs),
                                 function(c) data.table(cpg_site = cat_cpgs[[c]], category = c)))
writeLines(unique(cpg_category$cpg_site), here("B_MultiTissues/dataOut/CpG2extractS05.txt"))
message(sprintf("Wrote %d CpGs to extract", uniqueN(cpg_category$cpg_site)))
# Wrote 2438 CpGs to extract

## --- Extraction (run on pchuckle, then copy the .tsv back) --------------------
# In pchuckle (after git pull):
# source /share/apps/source_files/python/python-3.13.0a6.source
# cd /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/03_exploreResults
# python3 S00_extractRawMethylationForTargetCpG.py \
#   --cpg_list /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/dataOut/CpG2extractS05.txt \
#   --cpg_bed /SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/wgbs_tools/references/hg38/CpG.bed.gz \
#   --beta_files "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/GSM*.hg38.beta" \
#   --meta /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv \
#   --output /SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/gitignore/CpG2extractS05.tsv --minCov 10
## git push and pull in local

################################################################################
## PART B — load extracted methylation, tidy, keep multi-tissue patients       ##
################################################################################
meth <- fread(here("gitignore/CpG2extractS05.tsv"))
meth <- unique(meth, by = c("cpg_site","patient_id","source_tissue_celltype"))  # drop replicate runs

loyfer <- fread(here("B_MultiTissues/01_dataPrep/SupTab1_Loyfer2023_amended.csv"))
loyfer[, source_tissue_celltype := paste0(`Source Tissue`, " - ", `Cell type`)]
meth <- merge(meth, unique(loyfer[, .(source_tissue_celltype, germ_layer = `Germ layer`)]),
              by = "source_tissue_celltype", all.x = TRUE)
meth <- merge(meth, cpg_category, by = "cpg_site", all.x = TRUE)

multi      <- meth[, uniqueN(source_tissue_celltype), by = patient_id][V1 > 1, patient_id]
meth_multi <- meth[patient_id %in% multi]                    # patients sampled in >1 tissue

## Multi-tissue patient accounting (how much concordance power?)

## 1. patients sampled in >1 tissue, and which germ layers they span
per_patient <- meth_multi[, .(
  n_tissues   = uniqueN(source_tissue_celltype),
  germ_layers = paste(sort(unique(germ_layer)), collapse = "+")
), by = patient_id]
message(sprintf("%d patients with >1 tissue", nrow(per_patient)))
# 32 patients with >1 tissue

print(per_patient[, .N, by = germ_layers][order(-N)])          # patients per layer combination
# germ_layers     N
# 1:        Endo    13
# 2:        Meso     9
# 3:        Ecto     6
# 4:   Endo+Meso     4

## 2. patients with >=2 tissues WITHIN one germ layer  -> powers same-layer r (r_endo/meso/ecto)
same_layer <- meth_multi[, .(n = uniqueN(source_tissue_celltype)),
                         by = .(patient_id, germ_layer)][n >= 2]
print(same_layer[, .(n_patients = uniqueN(patient_id)), by = germ_layer][order(germ_layer)])
# germ_layer n_patients
# 1:       Ecto          6
# 2:       Endo         15
# 3:       Meso          9

## 3. patients spanning TWO layers-> powers any cross-layer comparison, per pair
pair_n <- function(a, b)
  meth_multi[germ_layer %in% c(a, b), uniqueN(germ_layer), by = patient_id][V1 == 2, .N]
message(sprintf("Donors with both: Endo+Meso %d | Endo+Ecto %d | Meso+Ecto %d",
                pair_n("Endo","Meso"), pair_n("Endo","Ecto"), pair_n("Meso","Ecto")))
# Donors with both: Endo+Meso 4 | Endo+Ecto 0 | Meso+Ecto 0

################################################################################
## PART C — within-layer systemic concordance per CpG                          ##
##   For a germ layer: per patient with >=2 tissues of that layer, pick 2       ##
##   tissues at random; for each CpG correlate those 2 tissues ACROSS patients. ##
##   High |r| = the between-individual variation is shared across that layer's  ##
##   tissues (systemic within the layer) -> early, clonally-propagated signal.  ##
################################################################################
layer_concordance <- function(meth, layer, seed = 1234) {
  d <- meth[germ_layer == layer, .(cpg_site, patient_id, source_tissue_celltype, methylation)]
  pts <- d[, uniqueN(source_tissue_celltype), by = patient_id][V1 >= 2, patient_id]
  d <- d[patient_id %in% pts]
  if (!nrow(d)) return(data.table(cpg_site = character()))
  set.seed(seed)
  picks <- unique(d[, .(patient_id, source_tissue_celltype)])[
    , .(source_tissue_celltype = sample(unique(source_tissue_celltype), 2),
        slot = c("A","B")), by = patient_id]                 # 2 random tissues per patient
  d <- merge(d, picks, by = c("patient_id","source_tissue_celltype"))
  w <- dcast(d, cpg_site + patient_id ~ slot, value.var = "methylation", fun.aggregate = mean)
  w[, {
    ok <- is.finite(A) & is.finite(B); n <- sum(ok)
    r  <- if (n >= 3 && sd(A[ok]) > 0 && sd(B[ok]) > 0) abs(cor(A[ok], B[ok])) else NA_real_
    .(abs_r = r, n_patients = n)
  }, by = cpg_site]
}

conc <- Reduce(function(a,b) merge(a,b, by = "cpg_site", all = TRUE), list(
  layer_concordance(meth_multi, "Endo")[, .(cpg_site, r_endo = abs_r, n_endo = n_patients)],
  layer_concordance(meth_multi, "Meso")[, .(cpg_site, r_meso = abs_r, n_meso = n_patients)],
  layer_concordance(meth_multi, "Ecto")[, .(cpg_site, r_ecto = abs_r, n_ecto = n_patients)]))
conc <- merge(conc, cpg_category, by = "cpg_site", all.x = TRUE)

## own-layer and other-layer concordance per CpG (mean |r| over the relevant layers)
## NB. only test categories whose defining layer(s) are actually sampled:
##   single-layer: fully powered;  endomeso: weak (4 donors);
##   endoecto / mesoecto: NOT testable (0 donors span those layer pairs) -> excluded
own_layers <- list(meso_only = "meso", endo_only = "endo", ecto_only = "ecto",
                   endomeso = c("endo","meso"))

conc[, own := mapply(function(cat,re,rm,rc){
  L <- own_layers[[cat]]; if (is.null(L)) NA_real_ else mean(c(endo=re,meso=rm,ecto=rc)[L], na.rm=TRUE)
}, category, r_endo, r_meso, r_ecto)]

conc[, other := mapply(function(cat,re,rm,rc){
  L <- own_layers[[cat]]; if (is.null(L)) NA_real_
  else mean(c(endo=re,meso=rm,ecto=rc)[setdiff(c("endo","meso","ecto"),L)], na.rm=TRUE)
}, category, r_endo, r_meso, r_ecto)]

## controls have no single "own" layer -> use the mean over all three layers
conc[category %in% c("control_none","control_3layers"),
     own := rowMeans(.SD, na.rm=TRUE), .SDcols=c("r_endo","r_meso","r_ecto")]

################################################################################
## PART D — tests                                                              ##
##  (i)  own-layer concordance of a category  vs  background, in the same layers ##
##  (ii) within category: own-layer  vs  other-layer concordance (specificity)  ##
################################################################################
test_tab <- rbindlist(lapply(names(own_layers), function(cat){
  cols <- paste0("r_", own_layers[[cat]])
  o   <- conc[category == cat, own];            o   <- o[is.finite(o)]
  bg  <- conc[category == "control_none",    rowMeans(.SD, na.rm=TRUE), .SDcols=cols]; bg  <- bg[is.finite(bg)]
  pos <- conc[category == "control_3layers", rowMeans(.SD, na.rm=TRUE), .SDcols=cols]; pos <- pos[is.finite(pos)]
  pr  <- conc[category == cat, .(own, other)][is.finite(own) & is.finite(other)]
  data.table(
    category = cat, n = length(o), n_bg = length(bg), n_pos = length(pos),
    median_own = if (length(o))   median(o)   else NA_real_,
    median_bg  = if (length(bg))  median(bg)  else NA_real_,
    median_pos = if (length(pos)) median(pos) else NA_real_,
    ## (i) own-layer concordance higher than background? (often NA: bg undefined, see note)
    p_vs_bg  = if (length(o)>=3 && length(bg)>=3)  wilcox.test(o, bg,  alternative="greater")$p.value else NA_real_,
    ## (ii) is it LESS concordant than systemic MEs? (both sets vary -> always computable)
    p_vs_pos = if (length(o)>=3 && length(pos)>=3) wilcox.test(o, pos, alternative="less")$p.value   else NA_real_,
    ## (iii) SPECIFICITY: own-layer vs other-layer, paired within the same CpGs (no background needed)
    p_own_gt_other = if (nrow(pr)>=3) wilcox.test(pr$own, pr$other, paired=TRUE, alternative="greater")$p.value else NA_real_)
}))
test_tab[, `:=`(p_vs_bg = p.adjust(p_vs_bg,"BH"),
                p_vs_pos = p.adjust(p_vs_pos,"BH"),
                p_own_gt_other = p.adjust(p_own_gt_other,"BH"))]
print(test_tab)
#     category     n  n_bg n_pos median_own median_bg median_pos p_vs_bg      p_vs_pos p_own_gt_other
# 1: meso_only    63     0   980  0.6794972        NA  0.9487253      NA  5.502890e-26   9.909766e-08
# 2: endo_only    18     0   980  0.6697061        NA  0.9410914      NA  8.959902e-11   6.582642e-02
# 3: ecto_only   297     0   980  0.3951251        NA  0.9411471      NA 3.725522e-102   1.165062e-05
# 4:  endomeso     5     0   980  0.9189015        NA  0.9364590      NA  4.624837e-02   9.375000e-02

################################################################################
## PART E — list germ-layer-specific candidates + annotate                     ##
##   Candidates = CpGs in the testable restricted categories. Concordance is    ##
##   carried over where it could be computed; annotation added where available. ##
################################################################################
suppressPackageStartupMessages({
  library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db); library(GenomicFeatures)
})

cand_cats  <- c("meso_only","endo_only","ecto_only","endomeso")
candidates <- rbindlist(lapply(cand_cats, function(c)
  data.table(cpg_site = cat_cpgs[[c]], category = c)))
## carry over concordance (own / other / per-layer r) where it was computed
candidates <- merge(candidates,
                    conc[, .(cpg_site, r_endo, r_meso, r_ecto, own, other)],
                    by = "cpg_site", all.x = TRUE)

## ---- feature sets (built once) ----------------------------------------------
txdb      <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr  <- genes(txdb); genes_gr$symbol <- mapIds(org.Hs.eg.db, names(genes_gr), "SYMBOL", "ENTREZID")
tss_gr    <- promoters(genes_gr, upstream = 0, downstream = 1)
prom_gr   <- promoters(txdb, upstream = 2000, downstream = 200)
exon_gr   <- reduce(exons(txdb))
intron_gr <- reduce(unlist(intronsByTranscript(txdb)))
if (!exists("te_regions")) {                       # reuse the TE block's object if present
  library(AnnotationHub); ah <- AnnotationHub(); rmsk <- ah[["AH111333"]]
  te_regions <- rmsk[mcols(rmsk)$repClass %in% c("LINE","SINE","LTR","DNA","RC","Retroposon")]
}

## ---- annotate each candidate CpG --------------------------------------------
cpg_to_gr <- function(cpg) GRanges(sub("_.*", "", cpg),
                                   IRanges(as.integer(sub(".*_", "", cpg)), width = 1), cpg_site = cpg)
g <- cpg_to_gr(candidates$cpg_site)
candidates[, `:=`(chr = as.character(seqnames(g)), pos = start(g))]

## nearest gene + distance
nr <- distanceToNearest(g, genes_gr, ignore.strand = TRUE)
candidates[, `:=`(gene = NA_character_, dist_to_gene = NA_integer_)]
candidates[queryHits(nr), `:=`(gene = genes_gr$symbol[subjectHits(nr)],
                               dist_to_gene = mcols(nr)$distance)]
## distance to nearest TSS
nt <- distanceToNearest(g, tss_gr, ignore.strand = TRUE)
candidates[, dist_to_tss := NA_integer_]
candidates[queryHits(nt), dist_to_tss := mcols(nt)$distance]
## genomic feature: promoter > exon > intron > intergenic
candidates[, feature := "intergenic"]
candidates[overlapsAny(g, intron_gr, ignore.strand = TRUE), feature := "intron"]
candidates[overlapsAny(g, exon_gr,   ignore.strand = TRUE), feature := "exon"]
candidates[overlapsAny(g, prom_gr,   ignore.strand = TRUE), feature := "promoter"]
## transposable element overlap (first TE per CpG, deterministic)
te <- findOverlaps(g, te_regions, ignore.strand = TRUE); te <- te[!duplicated(queryHits(te))]
candidates[, `:=`(in_TE = FALSE, TE_class = NA_character_, TE_family = NA_character_)]
candidates[queryHits(te), `:=`(in_TE = TRUE,
                               TE_class  = as.character(mcols(te_regions)$repClass [subjectHits(te)]),
                               TE_family = as.character(mcols(te_regions)$repFamily[subjectHits(te)]))]

## optional concordance flag (own-layer |r| >= 0.5; adjust threshold as you like)
candidates[, concordant := !is.na(own) & own >= 0.5]

## ---- tidy + save ------------------------------------------------------------
setcolorder(candidates, c("category","cpg_site","chr","pos","gene","dist_to_gene",
                          "feature","dist_to_tss","in_TE","TE_class","TE_family",
                          "own","other","r_endo","r_meso","r_ecto","concordant"))
setorder(candidates, category, -own, na.last = TRUE)   # most concordant first within category

fwrite(candidates, here("B_MultiTissues/dataOut/S05_layerSpecific_candidates_annotated.csv"))
message(sprintf("Wrote %d candidate CpGs", nrow(candidates)))
print(candidates[, .(n = .N, n_concordant = sum(concordant),
                     n_in_TE = sum(in_TE)), by = category])
## quick peek at the strongest per category
print(candidates[concordant == TRUE, head(.SD, 5),
                 by = category, .SDcols = c("cpg_site","gene","feature","in_TE","own")])

################################################################################
## PART F — best candidates: concordant WITHIN own layer, discordant in others ##
##   Per CpG: own_r must be HIGH in ALL its own layer(s); oth_r must be LOW (or  ##
##   undefined = doesn't even vary) in the other layer(s). A Fisher-z test       ##
##   compares own vs the strongest competing other layer.                        ##
################################################################################
sel <- conc[category %in% names(own_layers)]

## per CpG: own_r = MIN concordance over own layers (all must be high);
##          oth_r = MAX concordance over other layers (even the best must be low)
sel <- sel[, {
  L  <- own_layers[[category]]; O <- setdiff(c("endo","meso","ecto"), L)
  rv <- c(endo = r_endo, meso = r_meso, ecto = r_ecto)
  nv <- c(endo = n_endo, meso = n_meso, ecto = n_ecto)
  own_r <- if (all(is.na(rv[L]))) NA_real_ else min(rv[L], na.rm = TRUE)
  own_n <- if (all(is.na(nv[L]))) NA_real_ else min(nv[L], na.rm = TRUE)
  i     <- if (all(is.na(rv[O]))) NA_integer_ else which.max(rv[O])
  oth_r <- if (is.na(i)) NA_real_ else rv[O][i]       # NA = CpG doesn't vary in any other layer
  oth_n <- if (is.na(i)) NA_real_ else nv[O][i]
  .(own_r, own_n, oth_r, oth_n)
}, by = .(cpg_site, category)]

## Fisher-z test: is own-layer correlation GREATER than the best other-layer one?
## (only computable when the other layer actually varies, i.e. oth_r is defined)
sel[, z := (atanh(pmin(own_r, .999)) - atanh(pmin(oth_r, .999))) /
      sqrt(1/(own_n - 3) + 1/(oth_n - 3))]
sel[own_n > 3 & oth_n > 3, p_spec := pnorm(z, lower.tail = FALSE)]
sel[, p_spec_adj := p.adjust(p_spec, "BH")]

## ---- selection --------------------------------------------------------------
HI <- 0.5   # min own-layer |r| to count as "concordant within the layer"   (tune)
LO <- 0.3   # max other-layer |r| to count as "discordant in the others"    (tune)
sel[, best := !is.na(own_r) & own_r >= HI &
      (is.na(oth_r) | oth_r <= LO)]          # high own AND (no variation | low) in others

print(sel[, .(n = .N, n_best = sum(best),
              n_sig = sum(!is.na(p_spec_adj) & p_spec_adj < 0.05)), by = category])
best_candidates <- sel[best == TRUE][order(category, -own_r)]
print(best_candidates[, .(cpg_site, category, own_r, oth_r, p_spec_adj)])

## merge the specificity result onto the annotated candidates and write the best ones
candidates <- merge(candidates,
                    sel[, .(cpg_site, own_r, oth_r, p_spec_adj, best)],
                    by = "cpg_site", all.x = TRUE)

final2send <- candidates[best == TRUE][order(category, -own_r)]

fwrite(final2send,
       here("B_MultiTissues/dataOut/S05_bestLayerSpecific_candidates_annotated.csv"))
