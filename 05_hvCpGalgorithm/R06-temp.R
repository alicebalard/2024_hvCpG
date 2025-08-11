## 7th August
## Run on local laptop for hvCpGs and controls before poster

##########################
## Source the algorithm ##
source("hvCpG_algorithm_detection_v4scan.R")

###################
## Load the cpgs ##
cpg_names_all <- h5read("/home/alice/Documents/Project_hvCpG/10X/all_scaled_matrix.h5", "cpg_names")

###################################
## Which are covered in 26 cells ##
cpg_46 <- read.table("/home/alice/Documents/Project_hvCpG/10X/selected_cpgs_min3_in46_datasets.txt")$V1

#############################################
## Needed to perform liftover hg19 to hg38 ##
library(rtracklayer)
library(GenomicRanges)

# Download the chain file
chain_url <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
chain_gz <- "dataPrev/hg19ToHg38.over.chain.gz"
chain_file <- "dataPrev/hg19ToHg38.over.chain"
if (!file.exists(chain_file)) {
  download.file(chain_url, chain_gz)
  R.utils::gunzip(chain_gz, destname = chain_file, remove = FALSE)
}
chain <- import.chain(chain_file)

## Manifest illumina450k to check arrays
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

###########################
## Derakhshan2022_hvCpGs ##
DerakhshanhvCpGs <- readxl::read_excel("dataPrev/Derakhshan2022_4143hvCpGs_450k.xlsx", sheet = 6, skip = 3)

DerakhshanhvCpGs_GRanges <- GRanges(
  seqnames = anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"chr"],
  ranges = IRanges(start = ifelse(anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"strand"] %in% "+",
                                  anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"pos"],
                                  anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"pos"] - 1),
                   end = ifelse(anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"strand"] %in% "+",
                                anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"pos"] + 1,
                                anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"pos"])),
  strand = anno450k[match(DerakhshanhvCpGs$CpG, anno450k$Name),"strand"])

DerakhshanhvCpGs_GRanges_hg38 <- unlist(liftOver(DerakhshanhvCpGs_GRanges, chain))
length(DerakhshanhvCpGs_GRanges_hg38) # 3430

DerakhshanhvCpGs_names <- paste0(DerakhshanhvCpGs_GRanges_hg38@seqnames, "_", DerakhshanhvCpGs_GRanges_hg38@ranges)

# Restrict to CpGs in cpg_46
DerakhshanhvCpGs_names_filtered <- DerakhshanhvCpGs_names[DerakhshanhvCpGs_names %in% cpg_46]

# Find positions in cpg_names_all
DerakhshanhvCpGs_positions <- match(DerakhshanhvCpGs_names_filtered, cpg_names_all)

##########################################
## Matching genetic controls to hvCpGs  ##
mQTLcontrols <- read.table("../03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt", header = T)

mQTLcontrols_GRanges <- GRanges(
  seqnames = anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"chr"],
  ranges = IRanges(start = ifelse(anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"strand"] %in% "+",
                                  anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"pos"],
                                  anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"pos"] - 1),
                   end = ifelse(anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"strand"] %in% "+",
                                anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"pos"] + 1,
                                anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"pos"])),
  strand = anno450k[match(mQTLcontrols$controlCpG_name, anno450k$Name),"strand"])

mQTLcontrols_GRanges_hg38 <- unlist(liftOver(mQTLcontrols_GRanges, chain))

length(mQTLcontrols_GRanges_hg38) # 3644

mQTLcontrols_names <- paste0(mQTLcontrols_GRanges_hg38@seqnames, "_", mQTLcontrols_GRanges_hg38@ranges)

# Restrict to CpGs in cpg_46
mQTLcontrols_names_filtered <- mQTLcontrols_names[mQTLcontrols_names %in% cpg_46]

# Find positions in cpg_names_all
mQTLcontrols_positions <- match(mQTLcontrols_names_filtered, cpg_names_all)

########################################
## Run algorithm on hvCpG and control ##

# system.time(runAndSave(
#   analysis = "testLocalPC",
#   cpgPos_vec = DerakhshanhvCpGs_positions,
#   resultDir = "/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X/",
#   NCORES = 10,
#   p0 = 0.80,
#   p1 = 0.65, 
#   batch_size = 200)
# )

# system.time(runAndSave(
#   analysis = "testLocalPC",
#   cpgPos_vec = mQTLcontrols_positions,
#   resultDir = "/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X/",
#   NCORES = 10,
#   p0 = 0.80,
#   p1 = 0.65, 
#   batch_size = 200
# ))

#####################
## Explore results ##

load("/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X/results_testLocalPC_1663CpGs_0_8p0_0_65p1.RData")
reshvCpG = results_testLocalPC_1663CpGs_0_8p0_0_65p1
rm(results_testLocalPC_1663CpGs_0_8p0_0_65p1)

load("/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X/results_testLocalPC_1474CpGs_0_8p0_0_65p1.RData")
resControls = results_testLocalPC_1474CpGs_0_8p0_0_65p1
rm(results_testLocalPC_1474CpGs_0_8p0_0_65p1)

df = rbind(data.frame(alpha = reshvCpG, type = "hvCpG Derakhshan"),
           data.frame(alpha = resControls, type = "mQTL controls"))

ggplot(df, aes(x = type, y = alpha)) +
  geom_jitter(aes(fill=type), pch=21, size = 3, alpha = .1)+ 
  geom_violin(aes(col=type))+
  geom_boxplot(aes(col=type), width = .1) + 
  theme_minimal(base_size = 14) +                                                    
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylab("Probability of being a hvCpG") 

summary(lm(alpha~type, df))

# Create matching names for GRanges
gr_names <- paste0(
  as.character(seqnames(DerakhshanhvCpGs_GRanges_hg38)), "_",
  start(DerakhshanhvCpGs_GRanges_hg38), "-", end(DerakhshanhvCpGs_GRanges_hg38)
)

# Match alpha values to GRanges
DerakhshanhvCpGs_GRanges_hg38$alpha <- reshvCpG[gr_names]

###################################################################
## Compare values for Derakshan hvCpGs from arrays vs from Atlas ##

load("resultsDir/Mariads/results_MariasarraysREDUCED_3samples_15datasets_6906CpGs_0_8p0_0_65p1.RData")
Rred3array <- results_MariasarraysREDUCED_3samples_15datasets_6906CpGs_0_8p0_0_65p1

Rred3array_GRanges <- GRanges(
  seqnames = anno450k[match(rownames(Rred3array), anno450k$Name),"chr"],
  ranges = IRanges(start = ifelse(anno450k[match(rownames(Rred3array), anno450k$Name),"strand"] %in% "+",
                                  anno450k[match(rownames(Rred3array), anno450k$Name),"pos"],
                                  anno450k[match(rownames(Rred3array), anno450k$Name),"pos"] - 1),
                   end = ifelse(anno450k[match(rownames(Rred3array), anno450k$Name),"strand"] %in% "+",
                                anno450k[match(rownames(Rred3array), anno450k$Name),"pos"] + 1,
                                anno450k[match(rownames(Rred3array), anno450k$Name),"pos"])),
  strand = anno450k[match(rownames(Rred3array), anno450k$Name),"strand"],
  alpha = Rred3array)

Rred3array_GRanges$alpha <- Rred3array_GRanges$alpha.alpha
  
Rred3array_GRanges_hg38 <- unlist(liftOver(Rred3array_GRanges, chain))

# Find overlapping CpGs
hits <- findOverlaps(Rred3array_GRanges_hg38, DerakhshanhvCpGs_GRanges_hg38)

# Extract alpha values
alpha_dt <- data.table(
  alpha_x = Rred3array_GRanges_hg38$alpha[queryHits(hits)],
  alpha_y = DerakhshanhvCpGs_GRanges_hg38$alpha[subjectHits(hits)]
)

# Remove rows with NA values if needed
alpha_dt <- na.omit(alpha_dt)

cor_val <- cor(alpha_dt$alpha_x, alpha_dt$alpha_y)
ggplot(alpha_dt, aes(x = alpha_x, y = alpha_y)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 0.05, y = 0.95, hjust = 0, label = paste("r =", round(cor_val, 2))) +
  labs(
    x = "Alpha (Rred3array)",
    y = "Alpha (Derakhshan)"
  ) +
  theme_minimal(base_size = 14)

