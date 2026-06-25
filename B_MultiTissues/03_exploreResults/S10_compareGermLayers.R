#####################################################################
# Step 1 - Find inter-individual hypervariable CpGs per germ layer
# Step 2 — Test systemic intra-individual correlations
#####################################################################

## Prepare
library(here)
## Load libraries
source(here("B_MultiTissues", "quiet_library.R"))

## Load functions
if (!exists("functionsLoaded")) {
  source(here("B_MultiTissues/03_exploreResults", "functions.R"))}
####################################################################

# High Pr(hv) CpGs within a germ layer group are candidates set after that layer committed:
# they escaped the systemic lock and drifted tissue-specifically.
# Pr(hv)(Endo) >> Pr(hv)(Meso) ≈ Pr(hv)(Ecto) → set post-endoderm commitment
# Pr(hv)(Meso) ≈ Pr(hv)(Endo) >> Pr(hv)(Ecto) → set post-mesoderm commitment
# Pr(hv)(Ecto) >> the others → set late, post-ectoderm commitment
# All three Pr(hv) high and correlated → pre-gastrulation ME
#
# Step 2 — Test systemic intra-individual correlations
# For the CpGs flagged in Step 1, ask: does the methylation level at this CpG correlate across tissues within the same individual?
# Take donors who have samples from multiple germ layers (this table has several multi-tissue donors — check the PatientID column). For each CpG, compute the Pearson or Spearman correlation of methylation between, say, a liver sample (Endo) and a blood sample (Meso) across individuals.
#
# High cross-germ-layer correlation = systemic = set pre-gastrulation → ME candidate
# High within-germ-layer correlation only = set post-commitment in that lineage → germ-layer DMR (your Class ③ from the first figure)
# Low correlation everywhere = tissue-specific or experience-dependent
#
# The key contrast that tests the timing hypothesis
# PatternInterpretationHigh Pr(hv), high cross-layer rPre-gastrulation ME — earliestLow Pr(hv),
# high within-Endo r, low cross-layer rSet at endoderm commitmentLow Pr(hv), high within-Meso r,
# low cross-layer rSet at mesoderm commitmentLow Pr(hv), high within-Ecto r, low cross-layer 
# rSet at ectoderm commitment — latestLow Pr(hv) everywhere, no correlationConstitutive (always on/off)

###################################################################
#── 1. Load ───────────────────────────────────────────────────────────────────
# endo     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_12_endo.rds"))
# meso     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_meso.rds"))
# ecto     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_14_ecto.rds"))
# 
# analyses <- list(endo = endo, meso = meso, ecto = ecto)
# 
# # ── 2. Wide table ─────────────────────────────────────────────────────────────
# wide <- Reduce(
#   function(a, b) merge(a, b, by = "name"),
#   Map(function(dt, nm) setnames(copy(dt), "alpha", nm), analyses, names(analyses))
# )
# saveRDS(wide, here("gitignore/wide_script10_3layers.RDS"))
# rm(endo, meso, ecto, analyses)

wide <- readRDS(here("gitignore/wide_script10_3layers6gp.RDS"))

## Define a cutoff
HV <- 0.90 
notHV <- 0.5 

# ── 1. Categories based on 3 "one vs two others" comparisons ─────────────────
wide[, `:=`(
  HV_meso     = meso > HV,
  HV_endo     = endo > HV,
  HV_ecto     = ecto > HV,
  notHV_meso  = meso < notHV,
  notHV_endo  = endo < notHV,
  notHV_ecto  = ecto < notHV
)]

wide[, category := fcase(
  # ME: HV in all 3
  HV_meso & HV_endo & HV_ecto,
  "ME",
  
  # layer-specific: HV in one, clearly NOT HV in the other two
  HV_meso & notHV_endo & notHV_ecto, "Meso_specific",
  HV_endo & notHV_meso & notHV_ecto, "Endo_specific",
  HV_ecto & notHV_meso & notHV_endo, "Ecto_specific",
  
  # clearly not HV anywhere → constitutive
  notHV_meso & notHV_endo & notHV_ecto, "constitutive",
  
  # anything else: one layer HV but others not clearly low, or intermediate
  default = "ambiguous"
)]

print(table(wide$category))
# ambiguous  constitutive Ecto_specific Endo_specific            ME Meso_specific 
# 4026474      17287114         69748          1674        131844          5687 

# ── 2. Three plots: each layer vs its "two others" pool ───────────────────────
category_colours <- c(
  ME             = "#E69F00",
  Meso_specific  = "#56B4E9",
  Endo_specific  = "#009E73",
  Ecto_specific  = "#CC79A7",
  ambiguous      = "grey90", 
  constitutive   = "black"
)

set.seed(1234)

make_plot <- function(data, x_col, y_col, x_lab, y_lab, title) {
  data <- data[sample(nrow(data), 100000), ]
  
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], 
                   colour = category, shape = category)) +
    geom_point(data = data[category == "constitutive"],
               alpha = 0.3, size = 0.3) +
    geom_point(data = data[category == "ambiguous"],
               alpha = 0.4, size = 0.5) +
    geom_point(data = data[!category %in% c("constitutive", "ambiguous")],
               alpha = 0.8, size = 1.5) +
    geom_hline(yintercept = HV,    linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = HV,    linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = notHV, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = notHV, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
    scale_x_continuous(limits = c(0, 1), name = x_lab) +
    scale_y_continuous(limits = c(0, 1), name = y_lab) +
    scale_colour_manual(values = category_colours, drop = FALSE) +
    scale_shape_manual(values = c(
      ME             = 16,   # filled circle
      Meso_specific  = 16,
      Endo_specific  = 16,
      Ecto_specific  = 16,
      ambiguous      = 1,    # empty circle
      constitutive   = 4     # cross
    ), drop = FALSE) +
    ggtitle(title) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
}

p_mesovsendo <- make_plot(wide,
                          x_col = "meso",     y_col = "endo",
                          x_lab = "Pr(HV) meso",
                          y_lab = "Pr(HV) endo",
                          title = "Meso vs endo")

p_mesovsecto <- make_plot(wide,
                          x_col = "meso",     y_col = "ecto",
                          x_lab = "Pr(HV) meso",
                          y_lab = "Pr(HV) ecto",
                          title = "Meso vs ecto")

p_endovsecto <- make_plot(wide,
                          x_col = "endo",     y_col = "ecto",
                          x_lab = "Pr(HV) endo",
                          y_lab = "Pr(HV) ecto",
                          title = "Endo vs ecto")

# shared legend
legend_p <- ggplot(wide[sample(.N, 10000)],
                   aes(x = meso, y = endo, colour = category)) +
  geom_point(size = 3) +
  scale_colour_manual(values = category_colours, drop = FALSE, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_void() +
  theme(legend.position = "right")

(p_mesovsendo | p_mesovsecto | p_endovsecto | cowplot::get_legend(legend_p)) +
  plot_layout(widths = c(1, 1, 1, 0.35))

############################
## Test of SIV of targets ##
############################

## Test of raw data between people
Loyfer <- read.csv("../dataIn/SupTab1_Loyfer2023.csv")

dupPeople <- Loyfer[Loyfer$PatientID %in% Loyfer$PatientID[duplicated(Loyfer$PatientID)],]

# from the table above, build a summary of patients by their germ layer combination
layer_summary <- as.data.table(dupPeople)[
  , .(
    n_tissues   = .N,
    germ_layers = paste(sort(unique(Germ.layer)), collapse = "+")
  ),
  by = PatientID
]

# cross-tab: how many patients have each combination
result <- layer_summary[, .N, by = .(germ_layers, n_tissues)]
setorder(result, germ_layers, n_tissues)
result
# germ_layers n_tissues     N
# <char>     <int> <int>
#   1:        Ecto         2     5
# 2:        Ecto         3     1
# 3:        Ecto         4     1
# 4:        Endo         2     9
# 5:        Endo         3     5
# 6:   Endo+Meso         2     1
# 7:   Endo+Meso         3     1
# 8:   Endo+Meso         4     1
# 9:   Endo+Meso         5     1
# 10:        Meso         2     3
# 11:        Meso         3     1
# 12:        Meso         5     1
# 13:        Meso         6     2
# 14:        Meso         7     2

# Test if the meso-specific candidate are correlated within mesodermal t

####################
## GO enrichement ##
####################

