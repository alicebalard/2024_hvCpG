#############################################
## Overlap plot: Atlas (x) vs Array (y)    ##
#############################################
library(here)
source(here("05_hvCpGalgorithm/quiet_library.R"))

source(here("05_hvCpGalgorithm/runAlgo_myDatasets/Atlas/prephvCpGandControls.R"))
hvCpGandControls <- prephvCpGandControls(codeDir = "~/Documents/GIT/2024_hvCpG/")

################################ 
## --- Prepare Array data --- ##
################################ 

source(here("05_hvCpGalgorithm/runAlgo_myDatasets/exploreResults/S02_analyseResultsArray_local.R"))

################################
## --- Prepare Atlas data --- ##
################################
doIprepAtlas = TRUE
source(here("05_hvCpGalgorithm/runAlgo_myDatasets/exploreResults/S03_analyseResultsAtlas_local.R"))

###################################### 
## --- Merge Array & Atlas data --- ##
######################################

resCompArray$name <- resCompArray$chrpos ## key to merge
res_Alpha_Atlas <- dplyr::left_join(resCompArray, Atlas_dt, by = "name") %>%
  dplyr::select("alpha_array_all", "alpha_array_3ind", "alpha", "chrpos", "group.x") %>% 
  dplyr::rename(alpha_atlas=alpha, group=group.x)

######################### 
## --- Correlation --- ##
#########################

mod <- lm(alpha_atlas ~ alpha_array_all, data = res_Alpha_Atlas)
correlCoef <- mod$coefficients["alpha_array_all"]
summary(mod)

p1 <- ggplot(res_Alpha_Atlas, aes(x=alpha_array_all, y=alpha_atlas)) +
  geom_point(pch = 21, alpha = 0.05) +
  geom_abline(slope = 1, linetype = 3) +
  geom_smooth(linetype = 3)+
  geom_smooth(method = "lm", fill = "black") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.18,0.85),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(title = "Proba(hypervariable) at common CpGs",
       x = "P(hv) using array datasets",
       y = "P(hv) using WGBS atlas datasets") + 
  scale_x_continuous(breaks = seq(0, 1, by = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = .1))

pdf(here("05_hvCpGalgorithm/figures/correlation_array_atlas.pdf"), width = 5, height = 5)
p1
dev.off()

## What are alpha array for different cutoffs of alpha atlas?

df_plot <- bind_rows(
  res_Alpha_Atlas %>% filter(alpha_atlas > 0.8) %>% mutate(threshold = "> 0.8"),
  res_Alpha_Atlas %>% filter(alpha_atlas > 0.7) %>% mutate(threshold = "> 0.7"),
  res_Alpha_Atlas %>% filter(alpha_atlas > 0.6) %>% mutate(threshold = "> 0.6")
)

pdf(here("05_hvCpGalgorithm/figures/histCutoff-arrayAtlas.pdf"), width = 8, height = 5)
ggplot(df_plot, aes(x = alpha_array_all, fill=threshold)) +
  geom_histogram(bins = 30,color = "white") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("yellow", "orange", "red"))+
  xlab("p(hv) in array") +
  ylab("Count")
dev.off()

table(df_plot$threshold)

######################### 
## --- Sex effect? --- ##
#########################

## Male only atlas
atlas_dt_male <- mergeAtlasRunBatches(parent_dir = "05_hvCpGalgorithm/resultsDir/10X_males/", 
                                 analysis = "Atlas10X_males", alphaname = "alpha_atlas_males")

## Female only atlas
atlas_dt_female <- mergeAtlasRunBatches(parent_dir = "05_hvCpGalgorithm/resultsDir/10X_females/", 
                                      analysis = "Atlas10X_females", alphaname = "alpha_atlas_females")


atlas_dt_male_female <- full_join(atlas_dt_male, atlas_dt_female)

## Compare array all +-vs array CD4CD8
x <- left_join(atlas_dt_male_female, res_Alpha_Atlas)

p1 <- ggplot(x, aes(x=alpha_atlas_males, y=alpha_atlas_females, fill = group, col = group)) +
  geom_point(data = x[is.na(x$group),], pch = 21, alpha = 0.05) +
  geom_point(data = x[!is.na(x$group),], pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.3,0.8),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(x = "P(hv) considering male atlas data (12ds)",
       y = "P(hv) considering female atlas data (16ds)")

p2 <- ggplot(x, aes(x=alpha_atlas_males, y=alpha_atlas, fill = group, col = group)) +
  geom_point(data = x[is.na(x$group),], pch = 21, alpha = 0.05) +
  geom_point(data = x[!is.na(x$group),], pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.3,0.8),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(x = "P(hv) considering male atlas data (12ds)",
       y = "P(hv) considering all atlas data (46ds)")

p3 <- ggplot(x, aes(x=alpha_atlas_females, y=alpha_atlas, fill = group, col = group)) +
  geom_point(data = x[is.na(x$group),], pch = 21, alpha = 0.05) +
  geom_point(data = x[!is.na(x$group),], pch = 21, alpha = 0.4) +
  geom_smooth(method = "lm", fill = "black") +
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"), 
                    labels = c("hvCpG (Derakhshan)", "mQTL controls", "background")) +
  scale_colour_manual(values = c("#DC3220", "#005AB5", "grey"),guide = "none") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.3,0.8),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(x = "P(hv) considering female atlas data (16ds)",
       y = "P(hv) considering all atlas data (46ds)")

# --- Turn off legends inside plots ---
p1_clean <- p1 + theme(legend.position = "none")
p2_clean <- p2 + theme(legend.position = "none")
p3_clean <- p3 + theme(legend.position = "none")

# --- Extract one legend (e.g. from p1) ---
legend <- cowplot::get_legend(
  p1 + theme(legend.position = "bottom",
             legend.box = "horizontal",
             legend.title = element_blank(),
             legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
             legend.key = element_rect(fill = "white", color = NA))
)

# Combine plots
plots_grid <- plot_grid(
  p1_clean, p2_clean,
  p3_clean, legend,
  ncol = 2
)

# Add an overall title on top
final_plot <- ggdraw() +
  draw_label(
    "Probability of being hypervariable P(hv) in male, female, and combined WGBS atlas data",
    fontface = "bold",
    x = 0.5, y = 0.98, hjust = 0.5, vjust = 1, size = 18
  ) +
  draw_plot(plots_grid, x = 0, y = 0, width = 1, height = 0.95)

pdf(here("05_hvCpGalgorithm/figures/test4_compAtlasBothSexes.pdf"), width = 12, height = 12)
print(final_plot)
dev.off()


## tbc




