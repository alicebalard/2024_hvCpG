library(here)
if (!exists("libLoaded"))          source(here("B_MultiTissues", "quiet_library.R"))
if (!exists("functionsLoaded"))    source(here("B_MultiTissues/03_exploreResults", "functions.R"))
if (!exists("previousSIVprepared"))source(here("B_MultiTissues/03_exploreResults/prepPreviousSIV.R"))

tabRDS  <- readRDS(here("gitignore/table3layers_coveredIn3dt.RDS"))

Blood_evCpGs <- readxl::read_xlsx(here("B_MultiTissues/dataIn/Jimenez_13059_2020_2223_MOESM2_ESM.xlsx"), sheet = 2)
AT_Blood_evCpGs <- readxl::read_xlsx(here("B_MultiTissues/dataIn/Jimenez_13059_2020_2223_MOESM2_ESM.xlsx"), sheet = 3)

nrow(Blood_evCpGs) # 333
nrow(AT_Blood_evCpGs) # 154

## Get the probes from Planterose Jiménez et al. in hg38 format

x1 <- dico$chrpos_hg38[match(Blood_evCpGs$CpG, dico$CpG)] # full
x2 <- dico$chrpos_hg38[match(AT_Blood_evCpGs$CpG, dico$CpG)] # replicated in 2 tissues

## Get the hv scores for these probes
y1 <- tabRDS[tabRDS$chr_pos %in% x1,]
y2 <- tabRDS[tabRDS$chr_pos %in% x2,]
nrow(y1) # 190
nrow(y2) # 91

## Plot
med1   <- median(y1$percentile_logBF_per_ds_allLayers, na.rm = TRUE)
med2   <- median(y2$percentile_logBF_per_ds_allLayers, na.rm = TRUE)
med1 # 90.39304
med2 # 95.37517

## Prepare plot with both full and replicated
cols <- c("Blood" = "grey60", "Replicated in blood & adipose tissue" = "#B2182B")

## one long table with a grouping column
dat <- rbind(
  data.frame(score = y1$percentile_logBF_per_ds_allLayers, set = "Blood"),
  data.frame(score = y2$percentile_logBF_per_ds_allLayers, set = "Replicated in blood & adipose tissue")
)
dat$set <- factor(dat$set, levels = names(cols))

## medians, one per set, for the dashed lines
meds <- data.frame(set = factor(names(cols), levels = names(cols)),
                   med = c(med1, med2))

ggplot(dat, aes(x = score, fill = set)) +
  geom_histogram(breaks = seq(0, 100, 1), position = "identity",
                 alpha = .5, colour = "white", linewidth = .1) +
  geom_vline(data = meds, aes(xintercept = med, colour = set),
             linetype = "dashed", linewidth = .6, show.legend = FALSE) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_colour_manual(values = cols, guide = "none") +   # same palette for the median lines
  labs(x = "hvCpG score (percentile)", y = "Number of CpGs",
       title = "Hypervariability scores of evCpGs (Planterose Jiménez et al.)") +
  scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0.01, 0.01))) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11))

ggsave(here("B_MultiTissues/dataOut/figures/twin/evCpG_hvscore_hist.pdf"),
       width = 6, height = 4.5, dpi = 300, bg = "white")
