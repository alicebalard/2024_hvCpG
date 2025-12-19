###################################################################
## Plot secondary analyses results of algorithm ran on atlas data ##
####################################################################
library(here)

if (!exists("libLoaded")) {
  source(here("05_hvCpGalgorithm", "quiet_library.R"))
}

if (!exists("functionsLoaded")) {
  source(here("05_hvCpGalgorithm/exploreResults", "functions.R"))
}

#if (!exists("resCompArray")) {
#  source(here("05_hvCpGalgorithm/exploreResults/S02_analyseResultsArray_local.R"))  
#}

##################################
## Save all data in RDS objects ##
##################################

for (file in list.files(here("05_hvCpGalgorithm/resultsDir/Atlas/"))){
  if (!file.exists(here(paste0("gitignore/fullres_", file)))){
    ## Add previous MEs including Maria's results if not sourced yet
    if (!exists("KesslerSIV_GRanges_hg38")) {
      source(here("05_hvCpGalgorithm/exploreResults/prepPreviousSIV.R"))
    }
    
    system.time(Atlas_dt <- prepAtlasdt(file))
    
    print("Number of CpG tested:")
    print(nrow(Atlas_dt))
    
    print(paste0("Saving results for ", file, "in ", here(paste0("gitignore/fullres_", file))))
    saveRDS(Atlas_dt, file = here(paste0("gitignore/fullres_", file)))
    print("Saved")
  }
}

##################
## 1_byDevLayer ## --> not very useful
##################

if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_0_vs_1_byDevLayer.pdf")))){
  makeCompPlot(
    X = readRDS(here::here("gitignore/fullres_Atlas10X")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_1_byDevLayer")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_0_vs_1_byDevLayer",
    xlab = "Pr(hv) calculated on WGBS atlas datasets",
    ylab = "Pr(hv) calculated on WGBS atlas datasets cut in 3 groups (germ layers)")
}

#####################
## 2_rmMultSamples ## 
#####################

## Some individuals have multiple cells sampled. Does that affect our results? NOPE
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_0_vs_2_rmMultSamples.pdf")))){
  makeCompPlot(
    X = readRDS(here::here("gitignore/fullres_Atlas10X")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_2_rmMultSamples")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_0_vs_2_rmMultSamples",
    xlab = "Pr(hv) calculated on WGBS atlas datasets",
    ylab = "Pr(hv) calculated on WGBS atlas datasets keeping one sample/individual only")
}

###########################
## 3_correspMariaTissues ## --> very similar than full atlas
###########################

if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_0_vs_3_correspMariaTissues.pdf")))){
  makeCompPlot(
    X = readRDS(here("gitignore/fullres_Atlas10X")),
    Y = readRDS(here("gitignore/fullres_Atlas10X_3_correspMariaTissues")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_0_vs_3_correspMariaTissues",
    xlab = "Pr(hv) calculated on WGBS atlas datasets",
    ylab = "Pr(hv) calculated on WGBS atlas datasets with only cells found in array")
}

################
## Sex effect ##
################

## 1/ male
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_4_vs_5_maleEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here("gitignore/fullres_Atlas10X_4_maleOnly")),
    Y = readRDS(here("gitignore/fullres_Atlas10X_5_allButMaleOnly")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_4_vs_5_maleEffect",
    xlab = "Pr(hv) calculated on WGBS atlas males only datasets",
    ylab = "Pr(hv) calculated on WGBS atlas all but males only datasets")
}

## 2/ female
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_6_vs_7_femaleEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here("gitignore/fullres_Atlas10X_6_femaleOnly")),
    Y = readRDS(here("gitignore/fullres_Atlas10X_7_allButfemaleOnly")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_6_vs_7_femaleEffect",
    xlab = "Pr(hv) calculated on WGBS atlas females only datasets",
    ylab = "Pr(hv) calculated on WGBS atlas all but females only datasets")
}

################
## 8_byTissue ##
################

if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_0_vs_8_byTissue.pdf")))){
  ## Cut by tissue rather than by cell type. Is is closer to array data?
  makeCompPlot(
    X = readRDS(here::here("gitignore/fullres_Atlas10X")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_8_byTissue")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_0_vs_8_byTissue",
    xlab = "Pr(hv) calculated on WGBS atlas datasets (cut by cell types)",
    ylab = "Pr(hv) calculated on WGBS atlas datasets cut by tissues")
}

if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Array_vs_8_byTissue.pdf")))){
  makeCompPlot(
    X = resCompArray,
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_8_byTissue")),
    whichAlphaX = "alpha_array_all",
    whichAlphaY = "alpha",          
    title = "Array_vs_8_byTissue",
    xlab = "Pr(hv) calculated on array datasets",
    ylab = "Pr(hv) calculated on WGBS atlas datasets cut by tissues")
}

#########################
## 9_immune cells only ##
#########################

if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_0_vs_9_immuneOnly.pdf")))){
  makeCompPlot(
    
    X = readRDS(here::here("gitignore/fullres_Atlas10X")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_0_vs_9_immuneOnly",
    xlab = "Pr(hv) calculated on WGBS atlas datasets",
    ylab = "Pr(hv) calculated on WGBS atlas datasets, immune cells only")
}

#############################
## 10_immune cells removed ##
#############################

## Hypothesis testing:

## Prediction 1:
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_10_vs_9_immuneEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here::here("gitignore/fullres_Atlas10X_10_noImmune")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_10_vs_9_immuneEffect",
    xlab = "Pr(hv) calculated on WGBS atlas datasets, immune cells removed",
    ylab = "Pr(hv) calculated on WGBS atlas datasets, only immune cells")
}

system.time(Z_inner_immvsnoimm <- makeZ_inner(
  X = readRDS(here::here("gitignore/fullres_Atlas10X_10_noImmune")),
  Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
  whichAlphaX = "alpha",
  whichAlphaY = "alpha")) ## takes 2 minutes

print("Z_inner_immvsnoimm created")

stable <- Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X < 0.5 & Z_inner_immvsnoimm$alpha_Y < 0.5]
cellUniversal <- Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X > 0.5 & Z_inner_immvsnoimm$alpha_Y > 0.5]
immune <- Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X < 0.5 & Z_inner_immvsnoimm$alpha_Y > 0.5]
rest <-  Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X > 0.5 & Z_inner_immvsnoimm$alpha_Y < 0.5]

message(paste0("We found ", 
               length(stable), " (", round(length(stable)/length(Z_inner_immvsnoimm$name)*100), "%) stable CpGs, ",
               length(cellUniversal), " (", round(length(cellUniversal)/length(Z_inner_immvsnoimm$name)*100), "%) cell universal hvCpGs, ", 
               length(immune), " (", round(length(immune)/length(Z_inner_immvsnoimm$name)*100), "%) immune cell hvCpGs and ",
               length(rest), " (", round(length(rest)/length(Z_inner_immvsnoimm$name)*100), "%) hvCpGs undetected in immune cells, out of a total of ",
               length(Z_inner_immvsnoimm$name), " CpGs investigated"))

## We found 19672473 (86%) stable CpGs, 1047081 (5%) cell universal hvCpGs and 1464480 (6%) immune cell hvCpGs, out of a total of 22805115 CpGs investigated

print("Test enrichment in the 2 important categories:")

## Info:
# https://www.youtube.com/watch?v=_ycOp3P4AG0

# gr: GRanges of CpG positions (width=1)
# max_gap: maximum distance allowed between adjacent CpGs within a region
# min_cpg: minimum number of CpGs per merged region
# min_width: minimum merged width (bp)

mergeCpGsByGap <- function(ids, max_gap = 200, min_cpg = 3, min_width = 50) {
  ## Make GRanges object
  chr = sub("_.+$", "", ids, perl = TRUE)
  pos = as.integer(sub("^.*_", "", ids, perl = TRUE))
  gr = GRanges(
    seqnames = Rle(chr),
    ranges   = IRanges(start = pos, width = 1),
    genome   = "hg38")
  
  ## Sanity check: keep standard chromosomes
  gr = keepStandardChromosomes(gr, pruning.mode = "coarse")
  seqlevelsStyle(gr) <- "UCSC"  
  
  # Ensure sorted input
  gr <- sort(gr)
  
  # Reduce with a minimum gap width rule:
  # IRanges::reduce merges consecutive ranges if the gap BETWEEN them is <= max_gap
  merged <- reduce(gr, min.gapwidth = max_gap + 1)
  
  # Count CpGs per merged region
  ov <- findOverlaps(gr, merged, ignore.strand = TRUE)
  cpg_counts <- tabulate(subjectHits(ov), nbins = length(merged))
  merged$CpG_count <- cpg_counts
  merged$width <- width(merged)
  
  # Filter by counts and width
  keep <- (merged$CpG_count >= min_cpg) & (merged$width >= min_width)
  
  ## rm metadata which pose problem later on
  gr = merged[keep]
  mcols(gr) <- NULL
  
  return(gr)
}

for (x in c("cellUniversal", "immune")){
  gr <- mergeCpGsByGap(ids = get(x))
  ## rm metadata which pose problem later on
  mcols(gr) <- NULL
  
  if (!file.exists(here(paste0("05_hvCpGalgorithm/exploreResults/annotations_rGREAT/", x, ".RDS")))){
    system.time(res <- great(gr, "GO:BP", "hg38", cores = 10))
    saveRDS(res, file = here(paste0("05_hvCpGalgorithm/exploreResults/annotations_rGREAT/", x, ".RDS")))
    print(paste0("Annotation for ", x, " done!"))
  }
}

## Plot simplify GO for the regions
for (x in c("cellUniversal", "immune")){
  if (!file.exists(here(paste0("05_hvCpGalgorithm/exploreResults/annotations_rGREAT/plot_", x, "_heatmap.pdf")))){
    res = readRDS(here(paste0("05_hvCpGalgorithm/exploreResults/annotations_rGREAT/", x, ".RDS")))
    tab = res@table
    print("Pruning poorly informative and redundant enriched terms...")
    system.time(simplifiedGOterms <- compEpiTools::simplifyGOterms(
      goterms=tab$id, maxOverlap= 0.1, ontology='BP', go2allEGs = org.Hs.egGO2ALLEGS))
    print(paste0("Length simplifiedGOterms: ", length(simplifiedGOterms)))
    saveRDS(simplifiedGOterms, file = here(paste0("05_hvCpGalgorithm/exploreResults/annotations_rGREAT/simplifiedGOterms_", x, ".RDS")))
    print("Simplify and plot...")
    pdf(here(paste0("05_hvCpGalgorithm/exploreResults/annotations_rGREAT/plot_", x, "_heatmap.pdf")), width = 8, height = 6)
    simplifyEnrichment::simplifyGO(mat = simplifiedGOterms)
    dev.off()
  }
}

print("Enrichment done!")

simplifiedGOterms_immune <- readRDS(here("05_hvCpGalgorithm/exploreResults/annotations_rGREAT/simplifiedGOterms_immune.RDS"))
simplifiedGOterms_cellUniversal <- readRDS(here("05_hvCpGalgorithm/exploreResults/annotations_rGREAT/simplifiedGOterms_cellUniversal.RDS"))

simplifiedGOterms_immune@table


###################
## Prediction 2: ##
###################
Z_inner_all <- makeZ_inner(
  X = resCompArray,
  Y = readRDS(here::here("gitignore/fullres_Atlas10X")),
  whichAlphaX = "alpha_array_all",
  whichAlphaY = "alpha")

Z_inner_noimm <- makeZ_inner(
  X = resCompArray,
  Y = readRDS(here::here("gitignore/fullres_Atlas10X_10_noImmune")),
  whichAlphaX = "alpha_array_all",
  whichAlphaY = "alpha")

Z_inner_onlyimm <- makeZ_inner(
  X = resCompArray,
  Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
  whichAlphaX = "alpha_array_all",
  whichAlphaY = "alpha")

c_noimm <- cor.test(Z_inner_noimm$alpha_X, Z_inner_noimm$alpha_Y)
c_all <- cor.test(Z_inner_all$alpha_X, Z_inner_all$alpha_Y)
c_onlyimm <- cor.test(Z_inner_onlyimm$alpha_X, Z_inner_onlyimm$alpha_Y)

slope_noimm <- lm(data = Z_inner_noimm, alpha_Y ~ alpha_X)$coefficients[["alpha_X"]]
slope_all <- lm(data = Z_inner_all,  alpha_Y ~ alpha_X)$coefficients[["alpha_X"]]
slope_onlyimm <- lm(data = Z_inner_onlyimm,  alpha_Y ~ alpha_X)$coefficients[["alpha_X"]]

p <- data.frame(x=seq(0,1,0.01)) %>%
  dplyr::mutate(noImmun = slope_noimm * x,
                all = slope_all * x,
                onlyimm = slope_onlyimm * x) %>%
  pivot_longer(cols = c(noImmun, all, onlyimm)) %>%
  ggplot(aes(x = x, y = value, group = name, col = name)) +
  geom_abline(slope = 1, linetype = 3) +
  geom_line(linewidth = 2) +
  scale_color_viridis_d() +
  theme_minimal(base_size = 14) +
  labs(x = "Pr(hv) array datasets", y = "Pr(hv) WGBS datasets", col = "WGBS datasets") +
  ylim(c(0,1))

ggplot2::ggsave(
  filename = here::here("05_hvCpGalgorithm/figures/correlations/correlation_prediction2.pdf"),
  plot = p, width = 9, height = 7
)
