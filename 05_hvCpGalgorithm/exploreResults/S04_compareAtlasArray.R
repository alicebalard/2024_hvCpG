#############################################
## Overlap plot: Atlas (x) vs Array (y)    ##
#############################################
library(here)
source(here("05_hvCpGalgorithm", "quiet_library.R"))
source(here("05_hvCpGalgorithm/exploreResults", "functions.R"))

## Prepare Array data 
source(here("05_hvCpGalgorithm/exploreResults/S02_analyseResultsArray_local.R"))

X <- resCompArray
Y <- readRDS(here("gitignore/fullres_Atlas10X"))

# Ensure data.table types
setDT(X); setDT(Y)

# Keep only what you need
X <- X[, .(name = chrpos, alpha_array = alpha_array_all)]
Y <- Y[, .(name, alpha_atlas = alpha)]

# Set keys for efficient joins (optional but helpful on very large data)
setkey(X, name)
setkey(Y, name)

# Inner join (only names present in both)
Z_inner <- X[Y, nomatch = 0]

c <- cor.test(Z_inner$alpha_array, Z_inner$alpha_atlas)

p1 <- ggplot(Z_inner, aes(x=alpha_array, y=alpha_atlas)) +
  geom_point(pch = 21, alpha = 0.05) +
  geom_abline(slope = 1, linetype = 3) +
  geom_smooth(linetype = 3)+
  geom_smooth(method = "lm", fill = "black") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  annotate("text", x=.2, y=.8, label= paste0("R2 : ", round(c$estimate,2))) + 
  theme(legend.position.inside = c(0.18,0.85),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(title = "Pr(hv) at overlapping CpGs",
       x = "Pr(hv) calculated on array datasets",
       y = "Pr(hv) calculated on WGBS atlas datasets") + 
  scale_x_continuous(breaks = seq(0, 1, by = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = .1))

pdf(here("05_hvCpGalgorithm/figures/correlation_Arrray_Atlas10X.pdf"), width = 8, height = 8)
p1
dev.off()

#################################################################
## Is the slope different if I take array reduced to 3 ind/ds? ##
#################################################################

# Keep only what you need
X2 <- resCompArray
setDT(X2)
X2 <- X2[, .(name = chrpos, alpha_array = alpha_array_3ind)]

# Set keys for efficient joins (optional but helpful on very large data)
setkey(X2, name)

# Inner join (only names present in both)
Z_inner2 <- X2[Y, nomatch = 0]

c <- cor.test(Z_inner2$alpha_array, Z_inner2$alpha_atlas)

p1 <- ggplot(Z_inner2, aes(x=alpha_array, y=alpha_atlas)) +
  geom_point(pch = 21, alpha = 0.05) +
  geom_abline(slope = 1, linetype = 3) +
  geom_smooth(linetype = 3)+
  geom_smooth(method = "lm", fill = "black") +
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  annotate("text", x=.2, y=.8, label= paste0("R2 : ", round(c$estimate,2))) + 
  theme(legend.position.inside = c(0.18,0.85),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(title = "Pr(hv) at overlapping CpGs",
       x = "Pr(hv) calculated on array datasets with ONLY 3 individuals per dataset",
       y = "Pr(hv) calculated on WGBS atlas datasets") + 
  scale_x_continuous(breaks = seq(0, 1, by = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = .1))

pdf(here("05_hvCpGalgorithm/figures/correlation_Arrray3inds_Atlas10X.pdf"), width = 8, height = 8)
p1
dev.off()

################################################################
## What are alpha array for different cutoffs of alpha atlas? ##
################################################################

df_plot <- bind_rows(
  Z_inner %>% filter(alpha_atlas > 0.8) %>% mutate(threshold = "> 0.8"),
  Z_inner %>% filter(alpha_atlas > 0.7) %>% mutate(threshold = "> 0.7"),
  Z_inner %>% filter(alpha_atlas > 0.6) %>% mutate(threshold = "> 0.6")
)

pdf(here("05_hvCpGalgorithm/figures/histCutoff-arrayAtlas.pdf"), width = 8, height = 5)
ggplot(df_plot, aes(x = alpha_array, fill=threshold)) +
  geom_histogram(bins = 30,color = "white") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("yellow", "orange", "red"))+
  xlab("p(hv) in array") +
  ylab("Count")
dev.off()

table(df_plot$threshold)

########################################################################
## What are the positions with a high alpha in array but low in WGBS? ##
########################################################################

saveRDS(object = Z_inner[!is.na(Z_inner$alpha_array) & !is.na(Z_inner$alpha_atlas) &
                           Z_inner$alpha_array > 0.95 & Z_inner$alpha_atlas < 0.5,] %>% 
          dplyr::select(name) %>% dplyr::rename(chrpos = name) %>% as.vector(), here("05_hvCpGalgorithm/dataOut/CpGArray95moreAtlas50less.RDS"))

saveRDS(object = Z_inner[!is.na(Z_inner$alpha_array) & !is.na(Z_inner$alpha_atlas) &
                           Z_inner$alpha_array > 0.95 & Z_inner$alpha_atlas > 0.5,] %>% 
          dplyr::select(name) %>% dplyr::rename(chrpos = name) %>% as.vector(), here("05_hvCpGalgorithm/dataOut/CpGArray95moreAtlas50more.RDS"))

saveRDS(object = Z_inner[!is.na(Z_inner$alpha_array) & !is.na(Z_inner$alpha_atlas) &
                           Z_inner$alpha_array > 0.5 & Z_inner$alpha_atlas > 0.5,] %>% 
          dplyr::select(name) %>% dplyr::rename(chrpos = name) %>% as.vector(), here("05_hvCpGalgorithm/dataOut/CpGArray50moreAtlas50more.RDS"))
