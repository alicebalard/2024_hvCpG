#####################################################################
## Check pr(hv) geom means for target regions vs random background ##
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
if (!exists("previousSIVprepared")) {
  source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))}
#####################################################################

load(here("gitignore/fullTable3layers.Rda"))
seqlevels(table3layers) <- paste0("chr", seqlevels(table3layers))

############################
## Define candidate sites ##
############################

## Matt's data are in hg19
dataMatt <- readxl::read_xlsx(here("gitignore/DEGCAGS_intersect_repeats_Alice.xlsx"))

LTR41table <- dataMatt[grepl("LTR41", dataMatt$TE_family),]

# --- Liftover (hg19 → hg38) ---
LTR41table_gr <- GRanges(
  seqnames = LTR41table$chromosome,
  ranges = IRanges(start = LTR41table$start, end = LTR41table$end),
  hg19_chr = LTR41table$chromosome, hg19_start = LTR41table$start, hg19_end = LTR41table$end,
  `%change` = LTR41table$`%change`, padj = LTR41table$padj,
  TE_chromosome = LTR41table$TE_chromosome, TE_start = LTR41table$TE_start, 
  TE_end = LTR41table$TE_end, TE_family = LTR41table$TE_family,TE_type = LTR41table$TE_type
)

mapped <- unlist(liftOver(LTR41table_gr, chain))

# Find overlaps
hits <- findOverlaps(mapped, table3layers)

# Extract both sides and combine metadata
mapped_hits      <- mapped[queryHits(hits)]
table3layers_hits <- table3layers[subjectHits(hits)]

# Combine 
result_df <- as.data.frame(mapped_hits) %>%
  cbind(as.data.frame(mcols(table3layers_hits)))

chr_order <- paste0("chr", c(1:22, "X", "Y"))

result_df <- result_df %>%
  mutate(seqnames = factor(seqnames,
                           levels = chr_order[chr_order %in% unique(seqnames)]))
## Save for SIV test in fetal script
saveRDS(result_df, "fetalSIV/LTR41table.RDS")

# Get genome-wide means for reference lines
genome_means <- table3layers %>%
  as.data.frame() %>%
  summarise(across(c(alpha_endo, alpha_meso, alpha_ecto, alpha_geomean), 
                   ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "layer", values_to = "genome_mean") %>%
  mutate(layer = factor(layer,
                        levels = c("alpha_endo", "alpha_meso", "alpha_ecto", "alpha_geomean"),
                        labels = c("Endo", "Meso", "Ecto", "Geom.")))

plotLTR41 <- result_df %>%
  arrange(seqnames, start) %>%
  mutate(chr_pos = factor(chr_pos, levels = unique(chr_pos)))%>%
  pivot_longer(cols = c(alpha_endo, alpha_ecto, alpha_meso, alpha_geomean),
               names_to = "layer", values_to = "alpha") %>%
  mutate(layer = factor(layer,
                        levels = c("alpha_endo", "alpha_meso", "alpha_ecto", "alpha_geomean"),
                        labels = c("Endo", "Meso", "Ecto", "Geom."))) %>%
  ggplot(aes(x = layer, y = alpha, fill = layer)) +
  geom_col(width = 0.7) +
  geom_hline(data = genome_means,
             aes(yintercept = genome_mean, colour = layer),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = c(Endo = "#5C8AE0", Meso = "#4CAF82",
                               Ecto = "#E05C5C", "Geom." = "black"),
                    guide = "none") +
  scale_colour_manual(values = c(Endo = "#5C8AE0", Meso = "#4CAF82",
                                 Ecto = "#E05C5C", "Geom." = "black")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  facet_wrap(~ chr_pos, ncol = 6) +
  labs(x = NULL, y = "Pr(hvCpG)",
       caption = "Dashed lines = genome-wide mean per layer") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        strip.text     = element_text(size = 7))

ggplot2::ggsave(
  filename = here::here(paste0("B_MultiTissues/dataOut/figures/plotLTR41_prhv.pdf")),
  plot = plotLTR41, width = 8, height = 8
)

####################################
## Formal test (if enough points) ##
####################################

# Build GRanges from geometric mean
geomMeanGR <- GRanges(seqnames = table3layers@seqnames,
                      ranges = IRanges(start = table3layers@ranges@start, 
                                       end = table3layers@ranges@start),
                      alpha_geomean = table3layers$alpha_geomean)

geomMeanGR <- geomMeanGR[!is.na(geomMeanGR$alpha_geomean)]

sets <- list(
  mQTLcontrols = makeGRfromMyCpGPos(vec = mQTLcontrols_hg38, setname = "mQTLcontrols"),
  HarrisSIV = HarrisSIV_hg38_GR,
  VanBaakESS = VanBaakESS_hg38_GR,
  KesslerSIV = KesslerSIV_GRanges_hg38,
  CoRSIV = corSIV_GRanges_hg38,
  hvCpG = DerakhshanhvCpGs_hg38_GR,
  LTR41 = makeGRfromMyCpGPos(result_df$chr_pos, "LTR41")
)

# Now do overlap join for each set
MEsetdt <- rbindlist(lapply(names(sets), function(nm) {
  hits <- findOverlaps(sets[[nm]], geomMeanGR)
  data.table(
    alpha_geomean = geomMeanGR$alpha_geomean[subjectHits(hits)],
    ME    = nm
  )
}))

MEsetdt <- na.omit(MEsetdt) ## 69732

# Set controls as baseline
MEsetdt[, ME := relevel(factor(ME), ref = "mQTLcontrols")]

p1 <- ggplot(MEsetdt, aes(x = ME, y = alpha_geomean)) +
  geom_jitter(data = MEsetdt,
              aes(fill=ME), pch=21, size = 3, alpha = .05)+
  geom_violin(aes(col=ME))+
  geom_boxplot(aes(col=ME), width = .1) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylab("Pr(hv) (geometric mean)")

## Statistical comparisons of alpha between MEs

# Fit the model with controls as baseline
fit <- lm(alpha_geomean ~ ME, data = MEsetdt)

# Get estimated marginal means and contrasts vs baseline
emm <- emmeans(fit, ~ ME)
contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = "mQTLcontrols", adjust = "sidak") %>%
  as.data.frame()

emm
# ME           emmean      SE    df lower.CL upper.CL
# mQTLcontrols  0.211 0.00635 69709    0.199    0.224
# CoRSIV        0.315 0.00141 69709    0.312    0.318
# HarrisSIV     0.361 0.00974 69709    0.342    0.380
# hvCpG         0.526 0.00656 69709    0.513    0.539
# KesslerSIV    0.404 0.00684 69709    0.391    0.418
# VanBaakESS    0.474 0.01070 69709    0.453    0.495
# 
# Confidence level used: 0.95 

contrasts <- contrasts %>%
  mutate(ME = contrast,  # rename for clarity
         lower = estimate - 1.96*SE,
         upper = estimate + 1.96*SE)

p2 <- ggplot(contrasts, aes(x = ME, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    y = "Difference in Pr(hv) (geometric mean) vs mQTLcontrols",
    x = "",
    title = "Comparison of ME groups to mQTLcontrols",
    subtitle = "lm with multiple comparison correction (Sidak)"
  ) +
  theme_minimal()

cowplot::plot_grid(p1,p2, rel_widths = c(1, .8))


