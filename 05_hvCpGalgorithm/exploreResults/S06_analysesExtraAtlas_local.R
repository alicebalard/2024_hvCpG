####################################################################
## Plot secondary analyses results of algorithm ran on atlas data ##
####################################################################
library(here)

source(here("05_hvCpGalgorithm", "quiet_library.R"))
source(here("05_hvCpGalgorithm/exploreResults", "functions.R"))

##################################
## Save all data in RDS objects ##
##################################
saveAgain = FALSE
if (saveAgain){
  ## Add previous MEs including Maria's results
  source(here("05_hvCpGalgorithm/exploreResults/prepPreviousSIV.R"))
  
  for (file in list.files(here("05_hvCpGalgorithm/resultsDir/Atlas/"))){
    system.time(Atlas_dt <- prepAtlasdt(file))
    
    print("Number of CpG tested:")
    print(nrow(Atlas_dt)) # 23036026
    
    print(paste0("Saving results for ", file, "in ", here(paste0("gitignore/fullres_", file))))
    saveRDS(Atlas_dt, file = here(paste0("gitignore/fullres_", file)))
    print("Saved")
  }
}

##################
## 1_byDevLayer ## --> not very useful
##################

#####################
## 2_rmMultSamples ## 
#####################

## Some individuals have multiple cells sampled. Does that affect our results? NOPE

list.files(here("gitignore"))

X <- readRDS(here("gitignore/fullres_Atlas10X"))
Y <- readRDS(here("gitignore/fullres_Atlas10X_2_rmMultSamples"))

# Ensure data.table types
setDT(X); setDT(Y)

# Keep only what you need
X <- X[, .(name, alpha_X = alpha)]
Y <- Y[, .(name, alpha_Y = alpha)]

# Set keys for efficient joins (optional but helpful on very large data)
setkey(X, name)
setkey(Y, name)

# Inner join (only names present in both)
Z_inner <- X[Y, nomatch = 0]

c <- cor.test(Z_inner$alpha_X, Z_inner$alpha_Y)

set.seed(1234)
Z_inner_1M <- Z_inner[sample(nrow(Z_inner), min(nrow(Z_inner), 1000000)),]

p1 <- ggplot(Z_inner_1M, aes(x=alpha_X, y=alpha_Y)) +
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
  labs(title = "Effect of multiple samples per individual",
       x = "Pr(hv) calculated on WGBS atlas datasets with all samples",
       y = "Pr(hv) calculated on WGBS atlas datasets removing multiple samples in individuals") + 
  scale_x_continuous(breaks = seq(0, 1, by = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = .1))

pdf(here("05_hvCpGalgorithm/figures/correlation_Atlas10X_0_vs_2_rmMultSamples.pdf"), width = 8, height = 8)
p1
dev.off()

###########################
## 3_correspMariaTissues ## --> very similar than full atlas
###########################

list.files(here("gitignore"))

## Prepare Array data 
source(here("05_hvCpGalgorithm/exploreResults/S02_analyseResultsArray_local.R"))
X <- resCompArray
Y <- readRDS(here("gitignore/fullres_Atlas10X_3_correspMariaTissues"))

# Ensure data.table types
setDT(X); setDT(Y)

# Keep only what you need
X <- X[, .(name = chrpos, alpha_X = alpha_array_all)]
Y <- Y[, .(name, alpha_Y = alpha)]

# Set keys for efficient joins (optional but helpful on very large data)
setkey(X, name)
setkey(Y, name)

# Inner join (only names present in both)
Z_inner <- X[Y, nomatch = 0]

c <- cor.test(Z_inner$alpha_X, Z_inner$alpha_Y)

set.seed(1234)
Z_inner_1M <- Z_inner[sample(nrow(Z_inner), min(nrow(Z_inner), 1000000)),]

p1 <- ggplot(Z_inner_1M, aes(x=alpha_X, y=alpha_Y)) +
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
  labs(x = "Pr(hv) calculated on array datasets",
       y = "Pr(hv) calculated on WGBS atlas datasets with only cell types found in array's datasets") + 
  scale_x_continuous(breaks = seq(0, 1, by = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = .1))

pdf(here("05_hvCpGalgorithm/figures/correlation_Atlas10X_0_vs_3_correspMariaTissues.pdf"), width = 8, height = 8)
p1
dev.off()

################
## Sex effect ##
################

list.files(here("gitignore"))

## 1/ male

X <- readRDS(here("gitignore/fullres_Atlas10X_4_maleOnly"))
Y <- readRDS(here("gitignore/fullres_Atlas10X_5_allButMaleOnly"))

# Ensure data.table types
setDT(X); setDT(Y)

# Keep only what you need
X <- X[, .(name, alpha_X = alpha)]
Y <- Y[, .(name, alpha_Y = alpha)]

# Set keys for efficient joins (optional but helpful on very large data)
setkey(X, name)
setkey(Y, name)

# Inner join (only names present in both)
Z_inner <- X[Y, nomatch = 0]

c <- cor.test(Z_inner$alpha_X, Z_inner$alpha_Y)

# split chrX vs autosomes
Z_inner <- Z_inner[!grepl("chrM", Z_inner$name)& !grepl("chrY", Z_inner$name),]
Z_inner$chrtype <- "autosome"
Z_inner[grepl("chrX", Z_inner$name),"chrtype"] <- "chrX"

c1 <- cor.test(Z_inner$alpha_X[Z_inner$chrtype  %in% "autosome"], Z_inner$alpha_Y[Z_inner$chrtype  %in% "autosome"])
c2 <- cor.test(Z_inner$alpha_X[Z_inner$chrtype  %in% "chrX"], Z_inner$alpha_Y[Z_inner$chrtype  %in% "chrX"])

c1
c2

set.seed(1234)
Z_inner_1M <- Z_inner[sample(nrow(Z_inner), min(nrow(Z_inner), 1000000)),]

p1 <- ggplot(Z_inner_1M, aes(x=alpha_X, y=alpha_Y)) +
  geom_point(pch = 21, alpha = 0.05) +
  geom_abline(slope = 1, linetype = 3) +
  geom_smooth(linetype = 3)+
  geom_smooth(method = "lm", fill = "black") +
  facet_grid(.~chrtype)+
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.18,0.85),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(title = "Effect of sex (males only vs rest)",
       x = "Pr(hv) calculated on WGBS atlas datasets with only males",
       y = "Pr(hv) calculated on all other WGBS atlas datasets") + 
  scale_x_continuous(breaks = seq(0, 1, by = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = .1))

pdf(here("05_hvCpGalgorithm/figures/correlation_Atlas10X_0_vs_4-5_maleEffect.pdf"), width = 16, height = 8)
p1
dev.off()

## 2/ female

X <- readRDS(here("gitignore/fullres_Atlas10X_6_femaleOnly"))
Y <- readRDS(here("gitignore/fullres_Atlas10X_7_allButfemaleOnly"))

# Ensure data.table types
setDT(X); setDT(Y)

# Keep only what you need
X <- X[, .(name, alpha_X = alpha)]
Y <- Y[, .(name, alpha_Y = alpha)]

# Set keys for efficient joins (optional but helpful on very large data)
setkey(X, name)
setkey(Y, name)

# Inner join (only names present in both)
Z_inner <- X[Y, nomatch = 0]

c <- cor.test(Z_inner$alpha_X, Z_inner$alpha_Y)

# split chrX vs autosomes
Z_inner <- Z_inner[!grepl("chrM", Z_inner$name)& !grepl("chrY", Z_inner$name),]
Z_inner$chrtype <- "autosome"
Z_inner[grepl("chrX", Z_inner$name),"chrtype"] <- "chrX"

c1 <- cor.test(Z_inner$alpha_X[Z_inner$chrtype  %in% "autosome"], Z_inner$alpha_Y[Z_inner$chrtype  %in% "autosome"])
c2 <- cor.test(Z_inner$alpha_X[Z_inner$chrtype  %in% "chrX"], Z_inner$alpha_Y[Z_inner$chrtype  %in% "chrX"])

c1
c2

set.seed(1234)
Z_inner_1M <- Z_inner[sample(nrow(Z_inner), min(nrow(Z_inner), 1000000)),]

p1 <- ggplot(Z_inner_1M, aes(x=alpha_X, y=alpha_Y)) +
  geom_point(pch = 21, alpha = 0.05) +
  geom_abline(slope = 1, linetype = 3) +
  geom_smooth(linetype = 3)+
  geom_smooth(method = "lm", fill = "black") +
  facet_grid(.~chrtype)+
  theme_minimal(base_size = 14) +
  guides(fill = guide_legend(position = "inside"))+
  theme(legend.position.inside = c(0.18,0.85),
        legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(title = "Effect of sex (females only vs rest)",
       x = "Pr(hv) calculated on WGBS atlas datasets with only females",
       y = "Pr(hv) calculated on all other WGBS atlas datasets") + 
  scale_x_continuous(breaks = seq(0, 1, by = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = .1))

pdf(here("05_hvCpGalgorithm/figures/correlation_Atlas10X_0_vs_6-7_femaleEffect.pdf"), width = 16, height = 8)
p1
dev.off()

################
## 8_byTissue ##
################

## Cut by tissue rather than by cell type. Is is closer to array data?

## Prepare Array data 
source(here("05_hvCpGalgorithm/exploreResults/S02_analyseResultsArray_local.R"))
X <- resCompArray
Y <- readRDS(here("gitignore/fullres_Atlas10X_8_byTissue"))

# Ensure data.table types
setDT(X); setDT(Y)

# Keep only what you need
X <- X[, .(name = chrpos, alpha_X = alpha_array_all)]
Y <- Y[, .(name, alpha_Y = alpha)]

# Set keys for efficient joins (optional but helpful on very large data)
setkey(X, name)
setkey(Y, name)

# Inner join (only names present in both)
Z_inner <- X[Y, nomatch = 0]

c <- cor.test(Z_inner$alpha_X, Z_inner$alpha_Y)

set.seed(1234)
Z_inner_1M <- Z_inner[sample(nrow(Z_inner), min(nrow(Z_inner), 1000000)),]

p1 <- ggplot(Z_inner_1M, aes(x=alpha_X, y=alpha_Y)) +
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
  labs(title = "Effect of tissue vs cell type grouping",
       x = "Pr(hv) calculated on array datasets",
       y = "Pr(hv) calculated on WGBS atlas datasets separated by tissue (and not cell type)") + 
  scale_x_continuous(breaks = seq(0, 1, by = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = .1))

pdf(here("05_hvCpGalgorithm/figures/correlation_Array_vs_Atlas_8_byTissue.pdf"), width = 8, height = 8)
p1
dev.off()

## Cut by tissue rather than by cell type. Is is closer to array data?

## Prepare Array data 
source(here("05_hvCpGalgorithm/exploreResults/S02_analyseResultsArray_local.R"))
X <- readRDS(here("gitignore/fullres_Atlas10X"))
Y <- readRDS(here("gitignore/fullres_Atlas10X_8_byTissue"))

# Ensure data.table types
setDT(X); setDT(Y)

# Keep only what you need
X <- X[, .(name, alpha_X = alpha)]
Y <- Y[, .(name, alpha_Y = alpha)]

# Set keys for efficient joins (optional but helpful on very large data)
setkey(X, name)
setkey(Y, name)

# Inner join (only names present in both)
Z_inner <- X[Y, nomatch = 0]

c <- cor.test(Z_inner$alpha_X, Z_inner$alpha_Y)

set.seed(1234)
Z_inner_1M <- Z_inner[sample(nrow(Z_inner), min(nrow(Z_inner), 1000000)),]

p1 <- ggplot(Z_inner_1M, aes(x=alpha_X, y=alpha_Y)) +
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
  labs(title = "Effect of tissue vs cell type grouping",
       x = "Pr(hv) calculated on WGBS atlas datasets separated by cell type",
       y = "Pr(hv) calculated on WGBS atlas datasets separated by tissue") + 
  scale_x_continuous(breaks = seq(0, 1, by = .1)) +
  scale_y_continuous(breaks = seq(0, 1, by = .1))

pdf(here("05_hvCpGalgorithm/figures/correlation_Atlas10X_0_vs_8_byTissue.pdf"), width = 8, height = 8)
p1
dev.off()
