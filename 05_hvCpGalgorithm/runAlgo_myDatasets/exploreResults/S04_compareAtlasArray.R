#############################################
## Overlap plot: Atlas (x) vs Array (y)    ##
#############################################
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

# Define parent folder containing all "Atlas_batchXXX" folders
parent_dir <- "05_hvCpGalgorithm/resultsDir/Atlas10X/"

# Get list of relevant RData files
rdata_files <- dir(parent_dir, pattern = "results_Atlas10X_100000CpGs_0_8p0_0_65p1\\.RData$", 
                   recursive = TRUE, full.names = TRUE)

## Check if all batches have ran
length(rdata_files) ## 4 so far (27 aug)

all_cpg_values <- numeric()
pb <- progress_bar$new(total = length(rdata_files), format = "📦 :current/:total [:bar] :percent")

for (file in rdata_files) {
  e <- new.env()
  load(file, envir = e)
  obj <- e[[ls(e)[1]]]
  if (is.matrix(obj)) obj <- obj[, 1]
  all_cpg_values <- c(all_cpg_values, obj)
  pb$tick()
}

# Create data.table from named vector
dt <- data.table(
  name = names(all_cpg_values),
  alpha = as.numeric(all_cpg_values)
)

rm(e, pb, all_cpg_values, obj, file, parent_dir, rdata_files)

atlas_dt <- dt[, .(name, alpha_atlas = alpha)]
setnames(atlas_dt, "name", "chrpos")

## --- Merge by chrpos ---
merged <- merge(array_dt, atlas_dt, by = "chrpos", all = FALSE)

## Do the alpha correlate?
mod <- lm(alpha_array~alpha_atlas, data = merged)

summary(mod)
# alpha_atlas 0.710523   0.016648   42.68   <2e-16 ***
## --- Plot ---
pdf("05_hvCpGalgorithm/figures/Atlas_vs_Array_overlap.pdf", width = 8, height = 6)
ggplot(merged, aes(x = alpha_atlas, y = alpha_array)) +
  # background cloud
  geom_point(data = merged[is.na(merged$group), ],
             color = "grey60", alpha = 0.2, size = 0.5) +
  # hvCpGs (red)
  geom_point(data = merged[merged$chrpos %in% hvCpGandControls$DerakhshanhvCpGs_names, ],
             color = "#DC3220", alpha = 0.8, size = 0.8) +
  # mQTL controls (blue)
  geom_point(data = merged[merged$chrpos %in% hvCpGandControls$mQTLcontrols_names, ],
             color = "#005AB5", alpha = 0.8, size = 0.8) +
  # any relationship between x and y
  geom_smooth(method = "lm", aes(col=group, fill = group))+
  scale_color_manual(values = c("#DC3220", "#005AB5", "grey60"))+
  scale_fill_manual(values = c("#DC3220", "#005AB5","grey60"))+
  # linear relationship between both alpha values
  geom_abline(intercept = mod$coefficients["(Intercept)"], 
              slope = mod$coefficients["alpha_atlas"], linetype = 2)+
  theme_minimal(base_size = 14) +
  labs(
    x = "Atlas alpha (probability hvCpG)",
    y = "Array alpha (probability hvCpG)",
    title = "Overlap of hvCpG results between Atlas and Array"
  ) +
  # theme(legend.position = "none")+
  coord_cartesian(xlim = 0:1, ylim= 0:1)
dev.off()

## Plot the difference between both
merged <- merged %>% mutate(diffArrayMinusAlpha=alpha_array- alpha_atlas)

pdf("05_hvCpGalgorithm/figures/Atlas_vs_Array_overlap_fig2.pdf", width = 5, height = 5)
ggplot(merged, aes(x=group, y=diffArrayMinusAlpha))+
  geom_jitter(data=merged[!merged$group %in% c("hvCpG_Derakhshan", "mQTLcontrols"),], col="grey", alpha=.1)+
  geom_jitter(data=merged[merged$group %in% "hvCpG_Derakhshan",], col="#DC3220", alpha=.1)+
  geom_jitter(data=merged[merged$group %in% "mQTLcontrols",], col="#005AB5", alpha=.1)+
  scale_fill_manual(values = c("#DC3220", "#005AB5", "grey"))+
  geom_violin(aes(fill = group), width=.5, alpha=.8) +
  geom_boxplot(aes(fill = group), width=0.1, color="black", alpha=0.8) +
  theme_minimal(base_size = 14)+
  theme(legend.position =  "none", axis.title.x = element_blank(), title = element_text(size=10))+
  ggtitle("P(array) minus P(atlas)")+
  ylab("Difference of probability of being hypervariable")
dev.off()

##############################################
## What is the most variable point on chr1? ##
##############################################
merged[merged$alpha_array>0.99 & merged$alpha_atlas>0.99,"chrpos"]
## chr1_18784484-18784485

## Explore
cpg_names_all <- rhdf5::h5read("/home/alice/Documents/Project_hvCpG/10X/all_matrix_noscale.h5", "cpg_names")

x <- match(merged[merged$alpha_array>0.99 & merged$alpha_atlas>0.99,"chrpos"], cpg_names_all)

subRawData <- h5read("/home/alice/Documents/Project_hvCpG/10X/all_matrix_noscale.h5",
                     "matrix", index = list(NULL, x))
metadata <- read.table("/home/alice/Documents/Project_hvCpG/10X/sample_metadata.tsv", sep ="\t", header = TRUE)

colnames(subRawData) = cpg_names_all[x]
subRawData <- as.data.frame(subRawData)

subRawData$samples = metadata$sample
subRawData$sample_groups = metadata$dataset

subRawData <- melt(subRawData)

# Create alternating colors for sample_groups
group_levels <- unique(subRawData$sample_groups)
color_values <- rep(c("black", "grey60"), length.out = length(group_levels))

ggplot(subRawData, aes(x = sample_groups, y = value, color = sample_groups)) +
  geom_point() +
  scale_color_manual(values = color_values) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"  # remove legend if not needed
  ) +
  labs(x = "Sample Groups", y = "Value")

# Exclude metadata columns
num_cols <- setdiff(names(subRawData), c("samples", "sample_groups"))

# Compute SD per group per CpG site
df_long <- subRawData %>%
  dplyr::group_by(sample_groups, variable) %>%
  dplyr::summarise(sd = sd(value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::rename(CpG = variable)

# Histogram of SD distributions per group
ggplot(df_long, aes(x = sd)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ CpG, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution of SD per CpG",
       x = "Standard deviation across samples",
       y = "Count")



