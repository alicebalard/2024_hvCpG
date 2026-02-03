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

if (!exists("resCompArray")) {
  source(here("05_hvCpGalgorithm/exploreResults/S02_analyseResultsArray_local.R"))
}

###########################################
## Prepare putative MEs GRanges objects
vmeQTL_hg19probes <- readxl::read_xlsx(here("05_hvCpGalgorithm/dataPrev/vmeQTL_vCpG_359pair_sig_Zhang2025.xlsx"))
vmeQTL_hg38 <- na.omit(dico$chrpos_hg38[match(vmeQTL_hg19probes$vCpG, dico$CpG)]) ; rm(vmeQTL_hg19probes)
vmeQTL_hg38_GR <- makeGRfromMyCpGPos(vmeQTL_hg38, "vmeQTL (Zhang2025)")

HarrisSIV_hg38_GR <- makeGRfromMyCpGPos(HarrisSIV_hg38, "HarrisSIV")

VanBaakESS_hg38_GR <- makeGRfromMyCpGPos(VanBaakESS_hg38, "VanBaakESS")

KesslerSIV_GRanges_hg38$set <- "KesslerSIV"

corSIV_GRanges_hg38$set <- "Gunasekara 2019 corSIV"

DerakhshanhvCpGs_hg38_GR <- makeGRfromMyCpGPos(DerakhshanhvCpGs_hg38, "hvCpG Derakhshan 2022")

SoCCpGs_hg38_GR <- makeGRfromMyCpGPos(SoCCpGs_hg38, "Silver2022_SoCCpGs")

putativeME_GR <- c(vmeQTL_hg38_GR, HarrisSIV_hg38_GR, VanBaakESS_hg38_GR, KesslerSIV_GRanges_hg38,
                   corSIV_GRanges_hg38, DerakhshanhvCpGs_hg38_GR, SoCCpGs_hg38_GR)

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

#######################################
## 1 one run per developmental layer ##
#######################################

endo = readRDS(here::here("gitignore/fullres_10X_12_endo"))
meso = readRDS(here::here("gitignore/fullres_10X_13_meso"))
ecto = readRDS(here::here("gitignore/fullres_10X_14_ecto"))
allLayers = readRDS(here::here("gitignore/fullres_Atlas10X"))

WGBS_Array_datasets <- read.csv(here("05_hvCpGalgorithm/figures/WGBS_Array_datasets.csv"))

table(WGBS_Array_datasets[WGBS_Array_datasets$assay %in% "atlas", "Germ.layer"])
# ectoderm endoderm mesoderm 
# 6       21       19 

## What is the overlap for different alpha?

library(ggVennDiagram)
library(ggplot2)

# Compute overlap across any number of groups and plot a Venn diagram
plotMyVenn <- function(cutoff, ...) {
  groups <- list(...)                         # list of data.frames (each has name, alpha)
  
  # 1) Overlap among all names (unfiltered)
  sets_unfilt <- lapply(groups, function(df) df$name)
  overlap <- Reduce(intersect, sets_unfilt)
  
  message(paste0("There are ", length(overlap), " overlapping CpGs between these groups (unfiltered)."))
  
  # 2) Apply cutoff AND keep only overlapping names (so sets are comparable on the same universe)
  sets <- lapply(groups, function(df) df$name[df$alpha >= cutoff & df$name %in% overlap])
  
  # 3) Optional: add names to sets if you passed named arguments
  if (!is.null(dots <- match.call(expand.dots = FALSE)$...) && length(names(dots))) {
    nm <- names(dots)
    if (any(nzchar(nm))) names(sets) <- ifelse(nzchar(nm), nm, paste0("group", seq_along(sets)))
  }
  
  # 4) Plot
  p <- ggVennDiagram(sets, label_alpha = 0) +
    scale_fill_gradient2(low = "white", mid = "yellow", high = "red") +
    ggtitle(paste0("Pr(hv) ≥ ", cutoff), 
            subtitle = paste0("Sequenced CpGs N = ", length(overlap)))
  
  return(p)
}

p1 <- plotMyVenn(0.5, endo = endo, meso = meso, ecto = ecto, all = allLayers)
p2 <- plotMyVenn(0.75, endo = endo, meso = meso, ecto = ecto, all = allLayers)
p3 <- plotMyVenn(0.90, endo = endo, meso = meso, ecto = ecto, all = allLayers)

cowplot::plot_grid(p1 + theme(legend.position = "none"),
                   p2 + theme(legend.position = "none"),
                   p3 + theme(legend.position = "none"), 
                   nrow = 1, 
                   labels = "Ectoderm (N=6 cell types), mesoderm (N=19), endoderm (N=21), and all combined")

## Test enrichement of the most likely germ layer-universal hvCpG in previous MEs
total <- allLayers$name[allLayers$name %in% endo$name]
total <- total[total %in% meso$name]
total <- total[total %in% ecto$name]

top_cpgs <- intersect(
  intersect(allLayers$name[allLayers$alpha > 0.9], endo$name[endo$alpha > 0.9]),
  intersect(meso$name[meso$alpha > 0.9],ecto$name[ecto$alpha > 0.9]))
print(paste("Found", length(top_cpgs), "overlapping high-alpha CpGs"))

listGR <- list(top90 = makeGRfromMyCpGPos(vec = top_cpgs, setname = "topCpGs"),
  allButTop90 = makeGRfromMyCpGPos(total[!total %in% top_cpgs], "allButTop90"))

# --- Helper: safe Fisher with edge cases (all zeros, etc.)
.safe_fisher <- function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                dimnames = list(c("this_quadrant","others"), c("inME","notME")))
  # If any row or column totals are zero, Fisher’s test is not defined
  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    return(list(or = NA_real_, p = NA_real_))
  } else {
    ft <- fisher.test(mat)
    return(list(or = unname(ft$estimate), p = ft$p.value))
  }
}

# --- Main: test enrichment of ME for each quadrant vs the other three combined
test_enrichment_quadrants <- function(quad_list, putativeME_GR, me_col = "set") {
  # If no 'set' column, treat all ME as one group "ALL"
  if (!(me_col %in% names(mcols(putativeME_GR)))) {
    me_sets <- "ALL"
    putativeME_GR$..tmp_set.. <- "ALL"
    me_col <- "..tmp_set.."
  } else {
    me_sets <- unique(as.character(mcols(putativeME_GR)[[me_col]]))
  }
  
  # Total elements per quadrant (each element counted once)
  totals <- vapply(quad_list, length, integer(1))
  
  # Iterate over ME sets
  out <- lapply(me_sets, function(me_name) {
    me_subset <- putativeME_GR[mcols(putativeME_GR)[[me_col]] == me_name]
    
    # Count how many elements in each quadrant overlap the ME subset
    inME <- vapply(quad_list, function(gr) {
      ov <- findOverlaps(gr, me_subset, ignore.strand = TRUE)
      length(unique(queryHits(ov)))  # number of quadrant elements hitting at least one ME
    }, integer(1))
    
    # Build quadrant-vs-others 2x2 tests
    bind_rows(lapply(names(quad_list), function(q) {
      a <- inME[[q]]
      b <- totals[[q]] - a
      c <- sum(inME[names(inME) != q])
      d <- sum(totals[names(totals) != q]) - c
      
      fs <- .safe_fisher(a, b, c, d)
      
      tibble(
        ME_set            = me_name,
        quadrant          = q,
        inME              = a,
        total             = totals[[q]],
        others_inME       = c,
        others_total      = sum(totals) - totals[[q]],
        pct_inME          = 100 * a / totals[[q]],
        pct_inME_others   = 100 * c / (sum(totals) - totals[[q]]),
        odds_ratio        = fs$or,
        p_value           = fs$p
      )
    }))
  }) %>% bind_rows() %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
  
  out
}

# ---- Run it (ME sets in putativeME_GR$set will be tested separately)
res_quadrants <- test_enrichment_quadrants(listGR, putativeME_GR, me_col = "set")

# Order quadrants within each facet by log2OR
res_plot2 <- res_quadrants %>%
  mutate(
    log2OR = log2(odds_ratio),
    signif  = p_adj_BH < 0.05
  ) %>%
  group_by(ME_set) %>%
  mutate(quadrant_ord = reorder(quadrant, log2OR)) %>%
  ungroup()

plot <- ggplot(res_plot2, aes(x = quadrant_ord, y = log2OR, fill = signif)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("grey", "black")) +
  labs(
    x = NULL,
    y = expression(log[2]~"(odds ratio)"),
    title = "ME enrichment by group (vs other group)",
    subtitle = "Bars ordered by effect within ME set; dashed line = OR = 1"
  ) +
  facet_wrap(~ ME_set, scales = "free_x", nrow = 1) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

pdf(here("05_hvCpGalgorithm/figures/topCpGsEnrichME.pdf"), width = 12, height = 6)
plot
dev.off()

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

fullres_Atlas10X_5_femaleOnly6gp
fullres_Atlas10X_6_bothsexes6gp

## 1/ male
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_4_vs_6_maleEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here("gitignore/fullres_Atlas10X_4_maleOnly")),
    Y = readRDS(here("gitignore/fullres_Atlas10X_6_bothsexes6gp")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_4_vs_6_maleEffect",
    xlab = "Pr(hv) on WGBS atlas with 6 datasets of only males",
    ylab = "Pr(hv) on WGBS atlas with 6 datasets of mix males/females")
}

## 2/ female
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_5_vs_6_femaleEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here("gitignore/fullres_Atlas10X_5_femaleOnly6gp")),
    Y = readRDS(here("gitignore/fullres_Atlas10X_6_bothsexes6gp")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_5_vs_6_femaleEffect",
    xlab = "Pr(hv) on WGBS atlas with 6 datasets of only females",
    ylab = "Pr(hv) on WGBS atlas with 6 datasets of mix males/females")
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
if (!file.exists(file.path(here::here("05_hvCpGalgorithm/figures/correlations/correlation_Atlas_11_vs_9_immuneEffect.pdf")))){
  makeCompPlot(
    X = readRDS(here::here("gitignore/fullres_Atlas10X_11_noImmune_sample11groups")),
    Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
    whichAlphaX = "alpha",
    whichAlphaY = "alpha",          
    title = "Atlas_11_vs_9_immuneEffect",
    xlab = "Pr(hv) calculated on WGBS atlas datasets, immune cells removed (11 gp)",
    ylab = "Pr(hv) calculated on WGBS atlas datasets, only immune cells (11 gp)")
}

system.time(Z_inner_immvsnoimm <- makeZ_inner(
  X = readRDS(here::here("gitignore/fullres_Atlas10X_11_noImmune_sample11groups")),
  Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly")),
  whichAlphaX = "alpha",
  whichAlphaY = "alpha")) ## takes 2 minutes

print("Z_inner_immvsnoimm created")

stable <- Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X < 0.5 & Z_inner_immvsnoimm$alpha_Y < 0.5]
cellUniversal <- Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X > 0.5 & Z_inner_immvsnoimm$alpha_Y > 0.5]
immune <- Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X < 0.5 & Z_inner_immvsnoimm$alpha_Y > 0.5]
notimmune <-  Z_inner_immvsnoimm$name[Z_inner_immvsnoimm$alpha_X > 0.5 & Z_inner_immvsnoimm$alpha_Y < 0.5]

message(paste0("We found ", 
               length(stable), " (", round(length(stable)/length(Z_inner_immvsnoimm$name)*100), "%) stable CpGs, ",
               length(cellUniversal), " (", round(length(cellUniversal)/length(Z_inner_immvsnoimm$name)*100), "%) cell universal hvCpGs, ", 
               length(immune), " (", round(length(immune)/length(Z_inner_immvsnoimm$name)*100), "%) immune cell hvCpGs and ",
               length(notimmune), " (", round(length(notimmune)/length(Z_inner_immvsnoimm$name)*100), "%) hvCpGs undetected in immune cells, out of a total of ",
               length(Z_inner_immvsnoimm$name), " CpGs investigated"))

# We found 20268230 (84%) stable CpGs, 1190281 (5%) cell universal hvCpGs, 1507272 (6%) 
# immune cell hvCpGs and 1123906 (5%) hvCpGs undetected in immune cells, out of a total of 24089689 CpGs investigated

## Prediction 1: Using only the immune cells in WGBS Atlas should give a result closer to what is observed using the array

Z_inner_all <- makeZ_inner(
  X = resCompArray,
  Y = readRDS(here::here("gitignore/fullres_Atlas10X")),
  whichAlphaX = "alpha_array_all",
  whichAlphaY = "alpha")

Z_inner_noimm <- makeZ_inner(
  X = resCompArray,
  Y = readRDS(here::here("gitignore/fullres_Atlas10X_11_noImmune_sample11groups")),
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

# Combine datasets
Z_inner_all$group <- "all"
Z_inner_noimm$group <- "noImm"
Z_inner_onlyimm$group <- "onlyImm"

combined <- rbind(Z_inner_all, Z_inner_noimm, Z_inner_onlyimm)

# Fit model with interaction
model <- lm(alpha_Y ~ alpha_X * group, data = combined)
summary(model)

library(emmeans)
emm <- emtrends(model, ~ group, var = "alpha_X")
pairs(emm)

# contrast        estimate      SE     df t.ratio p.value
# all - noImm     -0.00505 0.00166 988312  -3.042  0.0066
# all - onlyImm   -0.10615 0.00169 988312 -62.974  <.0001
# noImm - onlyImm -0.10110 0.00163 988312 -62.005  <.0001

## Difference with controls in both cases?
data <- read.table(here("03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt"), header = T)

x = dico$chrpos_hg38[match(data$hvCpG_name, dico$CpG)]
y = dico$chrpos_hg38[match(data$controlCpG_name, dico$CpG)]

# Build mapping from hvCpG -> control
pairs <- data.frame(
  hvCpG = x,
  control = y,
  stringsAsFactors = FALSE
)

# Merge hvCpG alphas
res <- X[!is.na(group)]

hv_alpha <- res[, c("name", "alpha")]
colnames(hv_alpha) <- c("hvCpG", "alpha_hvCpG")

# Merge control alphas
ctrl_alpha <- res[, c("name", "alpha")]
colnames(ctrl_alpha) <- c("control", "alpha_control")

# Join everything
merged <- pairs %>%
  left_join(hv_alpha, by = "hvCpG") %>%
  left_join(ctrl_alpha, by = "control") %>%
  mutate(diffAlpha=alpha_hvCpG-alpha_control)

merged1 <- merged %>%
  mutate(chr = str_extract(hvCpG, "^chr[0-9XYM]+"))%>%
  filter(!is.na(diffAlpha))

p1 <- ggplot(merged1, aes(x="diff", y=diffAlpha))+
  geom_jitter(data=merged[merged$diffAlpha>=0,], col="black", alpha=.5)+
  geom_jitter(data=merged[merged$diffAlpha<0,], fill="yellow",col="black",pch=21, alpha=.5)+
  geom_violin(width=.5, fill = "grey", alpha=.8) +
  geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
  theme_minimal(base_size = 14)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), title = element_text(size=10))+
  ggtitle("P(hvCpG) minus P(matching control)", subtitle = "no immune cells")+
  ylab("Difference of probability") +
  coord_cartesian(ylim = c(-1,1))

# Merge hvCpG alphas
res <- Y[!is.na(group)]

hv_alpha <- res[, c("name", "alpha")]
colnames(hv_alpha) <- c("hvCpG", "alpha_hvCpG")

# Merge control alphas
ctrl_alpha <- res[, c("name", "alpha")]
colnames(ctrl_alpha) <- c("control", "alpha_control")

# Join everything
merged <- pairs %>%
  left_join(hv_alpha, by = "hvCpG") %>%
  left_join(ctrl_alpha, by = "control") %>%
  mutate(diffAlpha=alpha_hvCpG-alpha_control)

merged2 <- merged %>%
  mutate(chr = str_extract(hvCpG, "^chr[0-9XYM]+"))%>%
  filter(!is.na(diffAlpha))

p2 <- ggplot(merged2, aes(x="diff", y=diffAlpha))+
  geom_jitter(data=merged[merged$diffAlpha>=0,], col="black", alpha=.5)+
  geom_jitter(data=merged[merged$diffAlpha<0,], fill="yellow",col="black",pch=21, alpha=.5)+
  geom_violin(width=.5, fill = "grey", alpha=.8) +
  geom_boxplot(width=0.1, color="black", fill = "grey", alpha=0.8) +
  theme_minimal(base_size = 14)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), title = element_text(size=10))+
  ggtitle("P(hvCpG) minus P(matching control)", subtitle =  "immune cells only")+
  ylab("Difference of probability") +
  coord_cartesian(ylim = c(-1,1))

cowplot::plot_grid(p1, p2, ncol = 2)

wilcox.test(merged1$diffAlpha, merged2$diffAlpha)

mean(merged1$diffAlpha); median(merged1$diffAlpha)
mean(merged2$diffAlpha); median(merged$diffAlpha, na.rm = T)

##############################################################
## Are these hvCpG detected also in individual germ layers? ##
##############################################################

cellUniversal_GR <- makeGRfromMyCpGPos(cellUniversal, "cellUniversal")
immune_GR <- makeGRfromMyCpGPos(immune, "immune")
notimmune_GR <- makeGRfromMyCpGPos(notimmune, "notimmune")
stable_GR <- makeGRfromMyCpGPos(stable, "stable")

gr_list <- list(cellUniversal = cellUniversal_GR, 
                immune = immune_GR, 
                notimmune = notimmune_GR, 
                stable = stable_GR, 
                DerakhshanhvCpGs = DerakhshanhvCpGs_hg38_GR) ## test previous hvCpGs

lapply(names(gr_list), function(name){
  x <- gr_list[[name]]
  gr_names <- paste0(as.character(seqnames(x)), "_", start(x))
  p1 <- plotMyVenn(0.5, endo = endo[endo$name %in% gr_names,],
                   meso = meso[meso$name %in% gr_names,],
                   ecto = ecto[ecto$name %in% gr_names,],
                   all = allLayers[allLayers$name %in% gr_names,])
  p2 <- plotMyVenn(0.75, endo = endo[endo$name %in% gr_names,],
                   meso = meso[meso$name %in% gr_names,],
                   ecto = ecto[ecto$name %in% gr_names,],
                   all = allLayers[allLayers$name %in% gr_names,])
  p3 <- plotMyVenn(0.9, endo = endo[endo$name %in% gr_names,],
                   meso = meso[meso$name %in% gr_names,],
                   ecto = ecto[ecto$name %in% gr_names,],
                   all = allLayers[allLayers$name %in% gr_names,])
  pdf(file = paste0("../../05_hvCpGalgorithm/figures/vennGermLayers/Venn_", 
                    name, "_GR.pdf"), width = 13, height = 6)
  print(cowplot::plot_grid(p1 + theme(legend.position = "none"),
                           p2 + theme(legend.position = "none"),
                           p3 + theme(legend.position = "none"), 
                           nrow = 1, labels = name)) 
  dev.off()
})

# Prediction 2. The putative ME are those CpG sites which variability is high no matter which collection of cells one uses (“cell-universal” hvCpGs)

###########################################
## Find them in each quadrants
hits <- GenomicRanges::findOverlaps(cellUniversal_GR, putativeME_GR)
me_cellUniversal <- putativeME_GR[unique(subjectHits(hits))]

hits <- GenomicRanges::findOverlaps(immune_GR, putativeME_GR)
me_immune <- putativeME_GR[unique(subjectHits(hits))]

hits <- GenomicRanges::findOverlaps(notimmune_GR, putativeME_GR)
me_notimmune <- putativeME_GR[unique(subjectHits(hits))]

hits <- GenomicRanges::findOverlaps(stable_GR, putativeME_GR)
me_stable <- putativeME_GR[unique(subjectHits(hits))]

###########################################
# Compute enrichment for one ME set

# # Dependencies
# library(GenomicRanges)
# library(dplyr)
# library(tibble)

quads <- list(
  stable        = stable_GR,
  immune        = immune_GR,
  notimmune     = notimmune_GR,
  cellUniversal = cellUniversal_GR
)

# --- Helper: safe Fisher with edge cases (all zeros, etc.)
.safe_fisher <- function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                dimnames = list(c("this_quadrant","others"), c("inME","notME")))
  # If any row or column totals are zero, Fisher’s test is not defined
  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    return(list(or = NA_real_, p = NA_real_))
  } else {
    ft <- fisher.test(mat)
    return(list(or = unname(ft$estimate), p = ft$p.value))
  }
}

# --- Main: test enrichment of ME for each quadrant vs the other three combined
test_enrichment_quadrants <- function(quad_list, putativeME_GR, me_col = "set") {
  # If no 'set' column, treat all ME as one group "ALL"
  if (!(me_col %in% names(mcols(putativeME_GR)))) {
    me_sets <- "ALL"
    putativeME_GR$..tmp_set.. <- "ALL"
    me_col <- "..tmp_set.."
  } else {
    me_sets <- unique(as.character(mcols(putativeME_GR)[[me_col]]))
  }
  
  # Total elements per quadrant (each element counted once)
  totals <- vapply(quad_list, length, integer(1))
  
  # Iterate over ME sets
  out <- lapply(me_sets, function(me_name) {
    me_subset <- putativeME_GR[mcols(putativeME_GR)[[me_col]] == me_name]
    
    # Count how many elements in each quadrant overlap the ME subset
    inME <- vapply(quad_list, function(gr) {
      ov <- findOverlaps(gr, me_subset, ignore.strand = TRUE)
      length(unique(queryHits(ov)))  # number of quadrant elements hitting at least one ME
    }, integer(1))
    
    # Build quadrant-vs-others 2x2 tests
    bind_rows(lapply(names(quad_list), function(q) {
      a <- inME[[q]]
      b <- totals[[q]] - a
      c <- sum(inME[names(inME) != q])
      d <- sum(totals[names(totals) != q]) - c
      
      fs <- .safe_fisher(a, b, c, d)
      
      tibble(
        ME_set            = me_name,
        quadrant          = q,
        inME              = a,
        total             = totals[[q]],
        others_inME       = c,
        others_total      = sum(totals) - totals[[q]],
        pct_inME          = 100 * a / totals[[q]],
        pct_inME_others   = 100 * c / (sum(totals) - totals[[q]]),
        odds_ratio        = fs$or,
        p_value           = fs$p
      )
    }))
  }) %>% bind_rows() %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
  
  out
}

# ---- Run it (ME sets in putativeME_GR$set will be tested separately)
res_quadrants <- test_enrichment_quadrants(quads, putativeME_GR, me_col = "set")

# Order quadrants within each facet by log2OR
res_plot2 <- res_quadrants %>%
  mutate(
    log2OR       = log2(odds_ratio),
    signif  = p_adj_BH < 0.05
  ) %>%
  group_by(ME_set) %>%
  mutate(quadrant_ord = reorder(quadrant, log2OR)) %>%
  ungroup()

plot <- ggplot(res_plot2, aes(x = quadrant_ord, y = log2OR, fill = signif)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("grey", "black")) +
  labs(
    x = "Quadrant",
    y = expression(log[2]~"(odds ratio)"),
    title = "ME enrichment by quadrant (vs other three quadrants)",
    subtitle = "Bars ordered by effect within ME set; dashed line = OR = 1"
  ) +
  facet_wrap(~ ME_set, scales = "free_x", nrow = 1) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

pdf(here("05_hvCpGalgorithm/figures/quadrantsEnrichME.pdf"), width = 12, height = 6)
plot
dev.off()

## Prediction 3. CpGs that are only variable using immune cells only (“immune variable hvCpGs”) or using non-immune cells capture variability acquired later in life

###################
## Where is MHC? ##
###################
# https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13
chr6pos = Z_inner_immvsnoimm$name[grep("chr6", Z_inner_immvsnoimm$name)]
chr6pos = as.integer(sub(".*_", "", chr6pos))

MHCpos <- paste0("chr6_", chr6pos[chr6pos >= 28510120 & chr6pos <= 33480577])
length(MHCpos) # 56491

print(paste0(length(stable[stable %in% MHCpos]), " MHC CpGs are in the stable quadrant, ",
             length(immune[immune %in% MHCpos]), " in the immune one, ",
             length(notimmune[notimmune %in% MHCpos]), " in the not immune one, ",
             length(cellUniversal[cellUniversal %in% MHCpos]), " in the cell universal one"))

conting <- data.frame(category = c("stable", "immune", "notimmune", "cellUniversal"),
                      all = c(length(stable), length(immune),
                              length(notimmune), length(cellUniversal)),
                      MHC = c(length(stable[stable %in% MHCpos]), length(immune[immune %in% MHCpos]),
                              length(notimmune[notimmune %in% MHCpos]), length(cellUniversal[cellUniversal %in% MHCpos])))

conting$nonMHC <- conting$all - conting$MHC

tbl <- as.matrix(conting[, c("MHC", "nonMHC")])
rownames(tbl) <- conting$category

tbl

chisq.test(tbl)

fisher <- sapply(conting$category, function(cat) {
  MHC_cat <- conting$MHC[conting$category == cat]
  nonMHC_cat <- conting$nonMHC[conting$category == cat]
  
  MHC_other <- sum(conting$MHC) - MHC_cat
  nonMHC_other <- sum(conting$nonMHC) - nonMHC_cat
  
  data.frame(estimate = fisher.test(matrix(c(MHC_cat, nonMHC_cat,
                                             MHC_other, nonMHC_other),
                                           nrow = 2))$estimate,
             p.adj = p.adjust(fisher.test(matrix(c(MHC_cat, nonMHC_cat,
                                                   MHC_other, nonMHC_other),
                                                 nrow = 2))$p.value, method = "BH"))
  
})

fisher

## Plot
X = readRDS(here::here("gitignore/fullres_Atlas10X_11_noImmune_sample11groups"))
Y = readRDS(here::here("gitignore/fullres_Atlas10X_9_immuneOnly"))

MHC_noImmune = X[X$name %in% MHCpos,]
MHC_immune = Y[Y$name %in% MHCpos,]

Z_inner_MHC <- makeZ_inner(X = MHC_noImmune, Y = MHC_immune, whichAlphaX = "alpha", whichAlphaY = "alpha")

ggplot(Z_inner_MHC, aes(alpha_X, alpha_Y)) +
  geom_point(pch = 21, alpha = 0.05) +
  geom_abline(slope = 1, linetype = 3) +
  geom_smooth(linetype = 3) +
  geom_smooth(method = "lm", fill = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "MHC CpGs",
       x = "Pr(hv) calculated on WGBS atlas datasets, immune cells removed (11 gp)",
       y = "Pr(hv) calculated on WGBS atlas datasets, only immune cells (11 gp)")

## Manhattan plot of both cases
cowplot::plot_grid(plotManhattanFromdt(MHC_noImmune, transp = .3),
                   plotManhattanFromdt(MHC_immune, transp = .3),
                   nrow = 2, labels = c("immune excluded", "only immune"))

CairoPDF(here("05_hvCpGalgorithm/figures/Manhattan/ManhattanAlphaPlot_previoushvCpGplotted_atlas_immunevsnon.pdf"), width = 15, height = 7)
cowplot::plot_grid(plotManhattanFromdt(X),
                   plotManhattanFromdt(Y),
                   nrow = 2, labels = c("immune excluded", "only immune"))
dev.off()

####################
## Annotate genes ##
####################
print("Test enrichment by categories:")

## 1. Keep CpGs in regions where at least 5 CpGs are in 50bp distance to each other
## 2. annotate with associated genes (in gene body or +/- 10kb from TSS)
## 3. run GO term enrichment with clusterProfiler::enrichGO

## Create universe
universe <- annotateCpGs_txdb(
  clusterCpGs(Z_inner_immvsnoimm$name, max_gap = 50, min_size = 5),
  tss_window = 10000
)

print(paste0("Gene universe contains ", length(universe), " genes"))
## "Gene universe contains 30190 genes"

## Annotate the different CpGs categories 
resAnnot_cellUniversal <- CpG_GO_pipeline(cellUniversal, universe = universe)
resAnnot_immune <- CpG_GO_pipeline(immune, universe = universe)
resAnnot_notimmune <- CpG_GO_pipeline(notimmune, universe = universe)

resAnnot <- list(resAnnot_cellUniversal=resAnnot_cellUniversal,
                 resAnnot_immune=resAnnot_immune,
                 resAnnot_notimmune=resAnnot_notimmune)

# 1) Flatten the nested list: groups -> ontologies -> result rows
# resAnnot names: resAnnot_cellUniversal, resAnnot_immune, resAnnot_notimmune
grp_label_map <- c(
  resAnnot_cellUniversal = "cellUniversal",
  resAnnot_immune        = "immune",
  resAnnot_notimmune     = "not immune"
)

df_all <- purrr::imap(resAnnot, function(ont_list, grp_name) {
  purrr::imap_dfr(ont_list, function(er, ont_name) {
    if (is.null(er) || nrow(er@result) == 0) return(tibble())
    as_tibble(er@result) %>%
      mutate(group_raw = grp_name,  ontology = ont_name)
  })
}) %>% bind_rows()

# 2) Harmonize labels / factors
df_all <- df_all %>%
  mutate(
    group    = grp_label_map[group_raw],
    group    = factor(group, levels = c("cellUniversal", "immune", "not immune")),
    ontology = factor(ontology, levels = c("BP", "MF", "CC"))
  )

# 3) Keep only terms that are significant AND have counts AND more than 5 genes
df_sig <- df_all %>%
  filter(!is.na(p.adjust) & p.adjust < 0.01) %>%
  filter(!is.na(Description) & !is.na(FoldEnrichment) & !is.na(Count) & Count > 10 & FoldEnrichment > 2)

# 4) Optional: reorder Description within each ontology by (high FE, low FDR)
df_sig <- df_sig %>%
  group_by(ontology) %>%
  mutate(Description = fct_reorder2(Description, FoldEnrichment, -p.adjust, .fun = max, .desc = TRUE)) %>%
  ungroup()

# 5) Plot: X = group, Y = Description, facet by ontology, size = FE, color = p.adjust
p <- ggplot(df_sig, aes(x = group, y = Description)) +
  geom_point(aes(size = FoldEnrichment, color = p.adjust), alpha = 0.9) +
  scale_size_continuous(name = "FoldEnrichment", range = c(1.5, 8)) +
  # smaller FDR should look darker; direction = -1 handles that
  scale_color_viridis_c(name = "FDR (p.adjust)", option = "plasma", direction = -1) +
  facet_grid(group ~ ontology, scales = "free", space = "free_y") +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = NULL,
       title = "GO enrichment across groups and ontologies (FDR < 0.05)") +
  theme(
    legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
    legend.background     = element_rect(fill = "#ebebeb", color = "#ebebeb"),
    legend.key            = element_rect(fill = "#ebebeb", color = "#ebebeb"),
    legend.position       = "top",
    axis.text.y           = element_text(size = 8),
    axis.text.x           = element_text(size = 8, angle = 45, hjust = 1),
    strip.text            = element_text(face = "bold")
  )

print(p)

p <- ggplot(df_sig[df_sig$ontology %in% "BP",], aes(x = group, y = Description)) +
  geom_point(aes(size = FoldEnrichment, color = p.adjust), alpha = 0.9) +
  scale_size_continuous(name = "FoldEnrichment", range = c(1.5, 8)) +
  # smaller FDR should look darker; direction = -1 handles that
  scale_color_viridis_c(name = "FDR (p.adjust)", option = "plasma", direction = -1) +
  facet_grid(group ~ ontology, scales = "free", space = "free_y") +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = NULL,
       title = "GO enrichment across groups and ontologies (FDR < 0.05)") +
  theme(
    legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
    legend.background     = element_rect(fill = "#ebebeb", color = "#ebebeb"),
    legend.key            = element_rect(fill = "#ebebeb", color = "#ebebeb"),
    legend.position       = "top",
    axis.text.y           = element_text(size = 8),
    axis.text.x           = element_text(size = 8, angle = 45, hjust = 1),
    strip.text            = element_text(face = "bold")
  )

pdf(here("05_hvCpGalgorithm/figures/GOenrichment.pdf"), width = 10, height = 7)
p
dev.off()

## Nothing very clear...
