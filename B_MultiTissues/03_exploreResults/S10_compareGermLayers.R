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
### Same with 6 groups only per germ layer to exclude a power issue 
## that would explain the difference in ectoderm vs the rest

#── 1. Load ───────────────────────────────────────────────────────────────────
# endo6gp     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_12_2_endo6gp.rds"))
# meso6gp     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_13_2_meso6gp.rds"))
# ecto6gp     <- readRDS(here("gitignore/resultsAtlasPrepared/fullres_0_8p0_0_65p1_14_ecto.rds"))
#  
# analyses <- list(endo = endo6gp, meso = meso6gp, ecto = ecto6gp)
#  
# # ── 2. Wide table ─────────────────────────────────────────────────────────────
# wide <- Reduce(
#   function(a, b) merge(a, b, by = "name"),
#   Map(function(dt, nm) setnames(copy(dt), "alpha", nm), analyses, names(analyses))
# )
# 
# saveRDS(wide, here("gitignore/wide_script10_3layers6gp.RDS"))
# rm(endo6gp, meso6gp, ecto6gp, analyses)
wide <- readRDS(here("gitignore/wide_script10_3layers6gp.RDS"))

## Define a cutoff
HV <- 0.90 # top 10%
notHV <- 0.5 # bottom 50% are clearly not hv

# ── 1. Categories based on 3 "one vs two others" comparisons ─────────────────
wide[, `:=`(
  HV_meso     = meso     > HV,   # HV within mesoderm
  HV_endo     = endo     > HV,   # HV within endoderm
  HV_ecto     = ecto     > HV   # HV within ectoderm
)]

wide[, category := fcase(
  # HV all 3 layers → pre-gastrulation ME
  HV_meso & HV_endo & HV_ecto,
  "ME",
  
  # meso-specific stochastic
  HV_meso & !HV_endo & !HV_ecto,
  "Meso_specific",
  
  # endo-specific stochastic
  HV_endo & !HV_meso & !HV_ecto,
  "Endo_specific",
  
  # ecto-specific stochastic
  HV_ecto & !HV_meso & !HV_endo,
  "Ecto_specific",
  
  default = "constitutive"
)]

print(table(wide$category))
# constitutive Ecto_specific Endo_specific     ME     Meso_specific 
# 21729450        187031        155443        122180        133357 

# ── 2. Three plots: each layer vs its "two others" pool ───────────────────────
category_colours <- c(
  ME               = "#E69F00",   
  Meso_specific    = "#56B4E9",   
  Endo_specific    = "#009E73",   
  Ecto_specific    = "#CC79A7",   
  constitutive     = "grey80"
)

set.seed(1234)

make_plot <- function(data, x_col, y_col, x_lab, y_lab, title) {
  data = data[sample(nrow(data), 100000), ]
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], colour = category)) +
    geom_point(data = data[category == "constitutive"],
               alpha = 0.15, size = 0.3) +
    geom_point(data = data[category != "constitutive"],
               alpha = 0.6,  size = 1.5) +
    geom_hline(yintercept = HV, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = HV, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    scale_x_continuous(limits = c(0, 1), name = x_lab) +
    scale_y_continuous(limits = c(0, 1), name = y_lab) +
    scale_colour_manual(values = category_colours, drop = FALSE) +
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
legend_p <- ggplot(wide[sample(.N, 1000)],
                   aes(x = meso, y = endo, colour = category)) +
  geom_point(size = 3) +
  scale_colour_manual(values = category_colours, drop = FALSE, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_void() +
  theme(legend.position = "right")

(p_mesovsendo | p_mesovsecto | p_endovsecto | cowplot::get_legend(legend_p)) +
  plot_layout(widths = c(1, 1, 1, 0.35))

### With all data (power issue)
# # ── 1. Load ───────────────────────────────────────────────────────────────────
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
# 
# saveRDS(wide, here("gitignore/wide_script10_3layersfull.RDS"))
# rm(ecto, endo, meso, analyses)
# wide <- readRDS(here("gitignore/wide_script10_3layersfull.RDS"))
# 
# ## Define a cutoff
# HV <- 0.9 # top 10%
# notHV <- 0.5 # bottom 50% are clearly not hv
# 
# # ── 1. Categories based on 3 "one vs two others" comparisons ─────────────────
# wide[, `:=`(
#   HV_meso     = meso     > HV,   # HV within mesoderm
#   HV_endo     = endo     > HV,   # HV within endoderm
#   HV_ecto     = ecto     > HV   # HV within ectoderm
# )]
# 
# wide[, category := fcase(
#   # HV all 3 layers → pre-gastrulation ME
#   HV_meso & HV_endo & HV_ecto,
#   "ME",
#   
#   # meso-specific stochastic
#   HV_meso & !HV_endo & !HV_ecto,
#   "Meso_specific",
#   
#   # endo-specific stochastic
#   HV_endo & !HV_meso & !HV_ecto,
#   "Endo_specific",
#   
#   # ecto-specific stochastic
#   HV_ecto & !HV_meso & !HV_endo,
#   "Ecto_specific",
#   
#    default = "constitutive"
# )]
# 
# print(table(wide$category))
# # constitutive Ecto_specific Endo_specific     ME      Meso_specific 
# # 21121234        179727         28385        131844         61351 
# 
# # ── 2. Three plots: each layer vs its "two others" pool ───────────────────────
# category_colours <- c(
#   ME               = "#E69F00",   
#   Meso_specific    = "#56B4E9",   
#   Endo_specific    = "#009E73",   
#   Ecto_specific    = "#CC79A7",   
#   constitutive     = "grey80"
# )
# 
# set.seed(1234)
# 
# make_plot <- function(data, x_col, y_col, x_lab, y_lab, title) {
#   data = data[sample(nrow(data), 100000), ]
#   ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], colour = category)) +
#     geom_point(data = data[category == "constitutive"],
#                alpha = 0.15, size = 0.3) +
#     geom_point(data = data[category != "constitutive"],
#                alpha = 0.6,  size = 1.5) +
#     geom_hline(yintercept = HV, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
#     geom_vline(xintercept = HV, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
#     scale_x_continuous(limits = c(0, 1), name = x_lab) +
#     scale_y_continuous(limits = c(0, 1), name = y_lab) +
#     scale_colour_manual(values = category_colours, drop = FALSE) +
#     ggtitle(title) +
#     theme_bw(base_size = 11) +
#     theme(legend.position = "none")
# }
# 
# p_mesovsendo <- make_plot(wide,
#                     x_col = "meso",     y_col = "endo",
#                     x_lab = "Pr(HV) meso",
#                     y_lab = "Pr(HV) endo",
#                     title = "Meso vs endo")
# 
# p_mesovsecto <- make_plot(wide,
#                           x_col = "meso",     y_col = "ecto",
#                           x_lab = "Pr(HV) meso",
#                           y_lab = "Pr(HV) ecto",
#                           title = "Meso vs ecto")
# 
# p_endovsecto <- make_plot(wide,
#                           x_col = "endo",     y_col = "ecto",
#                           x_lab = "Pr(HV) endo",
#                           y_lab = "Pr(HV) ecto",
#                           title = "Endo vs ecto")
# 
# # shared legend
# legend_p <- ggplot(wide[sample(.N, 1000)],
#                    aes(x = meso, y = endo, colour = category)) +
#   geom_point(size = 3) +
#   scale_colour_manual(values = category_colours, drop = FALSE, name = NULL) +
#   guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
#   theme_void() +
#   theme(legend.position = "right")
# 
# (p_mesovsendo | p_mesovsecto | p_endovsecto | cowplot::get_legend(legend_p)) +
#   plot_layout(widths = c(1, 1, 1, 0.35))
