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

# Compute the percentile rank of each alpha_geomean value
# (what % of all values are <= this value)
table3layers$percentile <- ecdf(table3layers$alpha_geomean)(table3layers$alpha_geomean) * 100
# Percentile = 95 means the site is in the top 5%
# top X% = percentile >= (100 - X)

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

## Focus 10k around ACTL8 primary transcript (chr1:17755333-17827063)
target_window <- GRanges("chr1", IRanges(17755333-1000, 17827063+1000))  # region ±10kb
hits <- subsetByOverlaps(mapped, target_window, ignore.strand = TRUE)

# Find overlaps with full results
hits <- findOverlaps(hits, table3layers)

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

## Interlayer correlation obtained from fetal script
interlayer_corr <-
  readRDS(here("B_MultiTissues/dataOut/interlayer_corr_all.RDS"))

interlayer_corr$chr_pos <- dico$chrpos_hg38[match(interlayer_corr$CpG, dico$CpG)]

# Join interlayer_r to result_df by chr_pos
result_df_annot <- result_df %>%
  left_join(interlayer_corr %>%
              dplyr::select(chr_pos, interlayer_r, percentile_r) %>%
              distinct(chr_pos, .keep_all = TRUE),
            by = "chr_pos") %>%
  arrange(seqnames, start) %>%
  mutate(chr_pos_r = factor(
    paste0(chr_pos, "\nr = ", round(interlayer_r, 2)),
    levels = unique(paste0(chr_pos, "\nr = ", round(interlayer_r, 2)))
  ))

result_df_annot %>% dplyr::select(chr_pos, alpha_geomean, percentile, interlayer_r, percentile_r)

plot_df <- result_df_annot %>%
  arrange(seqnames, start) %>%
  mutate(chr_num = as.integer(factor(seqnames, levels = paste0("chr", c(1:22,"X","Y"))))) %>%
  dplyr::select(seqnames, start, chr_pos, alpha_geomean, percentile, interlayer_r, percentile_r, chr_num) %>%
  pivot_longer(cols = c(alpha_geomean, percentile, interlayer_r, percentile_r),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("alpha_geomean", "percentile", "interlayer_r", "percentile_r"),
                         labels = c("Pr(hv) (geom. mean)", "Percentile rank Pr(hv)", "Inter-layer r", "Percentile rank inter-layer r")))

# Add percentile reference line (e.g. 90th percentile threshold)
ref_lines <- data.frame(
  metric = c("Pr(hv) (geom. mean)", "Percentile rank Pr(hv)", "Inter-layer r", "Percentile rank inter-layer r"),
  ref    = c(genome_means$genome_mean[genome_means$layer == "Geom."],
             90,
             mean(result_df_annot$interlayer_r),
             90)
)

# Build offsets (unchanged)
chr_offsets <- result_df_annot %>%
  arrange(seqnames, start) %>%
  distinct(seqnames) %>%
  mutate(
    chr_num = as.integer(factor(seqnames, levels = paste0("chr", c(1:22,"X","Y")))),
    offset  = (chr_num - 1) * 5e8
  )

plot_df <- plot_df %>%
  left_join(chr_offsets, by = "seqnames") %>%
  mutate(x_pos = start + offset)

chr_labels <- plot_df %>%
  group_by(seqnames, offset) %>%
  summarise(x_mid = mean(start + offset), .groups = "drop")

# Replace metric labels with numbered versions to avoid any string matching issues
plot_df <- plot_df %>%
  mutate(metric_num = as.integer(metric))  # 1,2,3,4 in level order

# Custom labeller
metric_labels <- c(
  "1" = "Pr(hv) (geom. mean)",
  "2" = "Percentile rank Pr(hv)",
  "3" = "Inter-layer r",
  "4" = "Percentile rank inter-layer r"
)

# Match ref_lines to numeric too
ref_lines$metric_num <- as.integer(factor(ref_lines$metric,
                                          levels = levels(plot_df$metric)))

plot_df <- plot_df %>%
  mutate(metric_num = factor(as.integer(metric), 
                             levels = c(4, 3, 2, 1))) %>%
  ggplot(aes(x = x_pos, y = value, colour = seqnames)) +
  geom_point(size = 3) +
  geom_hline(data = ref_lines, aes(yintercept = ref),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = chr_offsets$offset[-1] - 2.5e8,
             colour = "grey85", linewidth = 0.4) +
  facet_wrap(~ metric_num, ncol = 1, scales = "free_y",
             as.table = FALSE,
             labeller = labeller(metric_num = metric_labels)) +
  scale_x_continuous(
    breaks = plot_df %>% distinct(chr_pos, x_pos) %>% pull(x_pos),
    labels = plot_df %>%
      distinct(chr_pos, x_pos) %>%
      mutate(label = format(as.integer(sub(".*_", "", chr_pos)),
                            big.mark = ",", scientific = FALSE)) %>%
      pull(label),
    expand = c(0.05, 0)
  ) +
  scale_colour_viridis_d(guide = "none") +
  labs(x = NULL, y = NULL,
       title = "LTR41 CpGs: Pr(hv), percentile rank and inter-germ layer correlation",
       caption = "Dashed lines: 90th percentile (percentile ranks) / genome wide mean (Pr(hv) and r)") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    panel.grid.major.x = element_line(colour = "grey92"),
    panel.grid.minor.x = element_blank(),
    strip.text         = element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = here::here(paste0("B_MultiTissues/dataOut/figures/SIV/plotLTR41_alongChr.pdf")),
  plot = plot_df, width = 8, height = 8
)

plotLTR41 <- result_df_annot %>%
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
  facet_wrap(~ chr_pos_r, ncol = 6) +
  labs(x = NULL, y = "Pr(hvCpG)",
       caption = "Dashed lines = genome-wide mean per layer. r = mean inter-layer correlation in fetus data") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        strip.text         = element_text(size = 7))

ggplot2::ggsave(
  filename = here::here(paste0("B_MultiTissues/dataOut/figures/SIV/plotLTR41_prhv.pdf")),
  plot = plotLTR41, width = 8, height = 8
)

result_df_annot$CpG <- dico$CpG[match(result_df_annot$chr_pos, dico$chrpos_hg38)]

p2 <- ggplot(result_df_annot, aes(x = alpha_geomean, y = interlayer_r, 
                                  label = paste(chr_pos, CpG, sep = "\n"))) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_label(aes(col = seqnames), size = 2) +
  theme_bw() +
  theme(legend.position = "none")

ggplot2::ggsave(
  filename = here::here(paste0("B_MultiTissues/dataOut/figures/SIV/plotLTR41_corprhv_r.pdf")),
  plot = p2, width = 5, height = 5
)

####################################################################################################
## Correlations pr(hv)(WGBS Loyfer) and inter-germ layer correlation (fetus) for all data on EPIC ##
####################################################################################################
interlayer_corr_all <- readRDS(here("B_MultiTissues/dataOut/interlayer_corr_all.RDS"))

## Add hg38 coordinates
interlayer_corr_all$chr_pos <- dico$chrpos_hg38[match(interlayer_corr_all$CpG, dico$CpG)]

interlayer_corr_all$alpha_geomean <- table3layers$alpha_geomean[
  match(interlayer_corr_all$chr_pos, table3layers$chr_pos)]

p3 <- ggplot(interlayer_corr_all, aes(x = alpha_geomean, y = abs(interlayer_r))) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = .05) +
  geom_smooth(col="red")+
  geom_smooth(method = "lm") +
  labs(x="Pr(hv) (geometric mean on the 3 germ layers analyses)",
       y="Inter-germ layer correlation (absolute value)")+
  theme_bw() +
  theme(legend.position = "none")

ggplot2::ggsave(
  filename = here::here(paste0("B_MultiTissues/dataOut/figures/SIV/allEPICfetus_corprhv_rinter.pdf")),
  plot = p3, width = 5, height = 5
)
