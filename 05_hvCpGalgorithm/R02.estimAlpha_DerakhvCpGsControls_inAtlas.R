#!/usr/bin/env Rscript
# 12th August
# Run algo v5 for hvCpGs and controls, without transformation

quiet_library <- function(pkg) {
  # Check if installed
  installed <- requireNamespace(pkg, quietly = TRUE)
  
  # Install if missing
  if (!installed) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      suppressMessages(suppressWarnings(
        install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
      ))
    }
    
    cran_pkgs <- suppressMessages(available.packages(repos = "https://cloud.r-project.org"))
    if (pkg %in% rownames(cran_pkgs)) {
      suppressMessages(suppressWarnings(
        install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
      ))
    } else {
      suppressMessages(suppressWarnings(
        BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
      ))
    }
  }
  
  # Load silently
  suppressPackageStartupMessages(
    suppressMessages(
      suppressWarnings(
        library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      )
    )
  )
  
  # Get version after loading
  version <- as.character(utils::packageVersion(pkg))
  
  # Print only our message
  cat(sprintf("Load package %s v%s\n", pkg, version))
}

quiet_library_all <- function(pkgs) {
  invisible(lapply(pkgs, quiet_library))
}

quiet_library_all(c(
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "GenomicRanges", "rtracklayer" ## Needed to perform liftover hg19 to hg38 ##
))

###############################
message("Define parameters...")
args <- commandArgs(trailingOnly = TRUE)
dataDir <- args[1] ## scratch datadir
nslots <- as.integer(args[2]) ## number of cores

message(paste0("We work in", dataDir, " with ", nslots, " cores."))

codeDir = "/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/"
resDir = file.path(codeDir, "05_hvCpGalgorithm/resultsDir/Atlas10X/")

##################################
message("Source the algorithm...")
source(file.path(codeDir, "05_hvCpGalgorithm/hvCpG_algorithm_detection_v5batches.R"))

###########################
message("Load the cpgs...")
cpg_names_all <- h5read(file.path(dataDir, "all_matrix_noscale.h5"), "cpg_names")

###################################
## Which are covered in 26 cells ##
cpg_46 <- read.table(file.path(dataDir, "selected_cpgs_min3_in46_datasets.txt"))$V1

#############################################
message("Download the chain file for liftover...")
chain_url <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
chain_gz <- file.path(codeDir, "05_hvCpGalgorithm/dataPrev/hg19ToHg38.over.chain.gz")
chain_file <- file.path(codeDir, "05_hvCpGalgorithm/dataPrev/hg19ToHg38.over.chain")
if (!file.exists(chain_file)) {
  download.file(chain_url, chain_gz)
  R.utils::gunzip(chain_gz, destname = chain_file, remove = FALSE)
}
chain <- import.chain(chain_file)

## Manifest illumina450k to check arrays
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

###########################################
message("Prepare Derakhshan2022_hvCpGs...")
DerakhshanhvCpGs <- readxl::read_excel(file.path(codeDir, "05_hvCpGalgorithm/dataPrev/Derakhshan2022_4143hvCpGs_450k.xlsx"), sheet = 6, skip = 3)

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
length(DerakhshanhvCpGs_GRanges_hg38) == 4143

DerakhshanhvCpGs_names <- paste0(DerakhshanhvCpGs_GRanges_hg38@seqnames, "_", DerakhshanhvCpGs_GRanges_hg38@ranges)

# Restrict to CpGs in cpg_46
DerakhshanhvCpGs_names_filtered <- DerakhshanhvCpGs_names[DerakhshanhvCpGs_names %in% cpg_46]

# Find positions in cpg_names_all
DerakhshanhvCpGs_positions <- match(DerakhshanhvCpGs_names_filtered, cpg_names_all)

###############################################
message("Matching genetic controls to hvCpGs..")
mQTLcontrols <- read.table(file.path(codeDir, "03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt"), header = T)

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

length(mQTLcontrols_GRanges_hg38) == 3644

mQTLcontrols_names <- paste0(mQTLcontrols_GRanges_hg38@seqnames, "_", mQTLcontrols_GRanges_hg38@ranges)

# Restrict to CpGs in cpg_46
mQTLcontrols_names_filtered <- mQTLcontrols_names[mQTLcontrols_names %in% cpg_46]

# Find positions in cpg_names_all
mQTLcontrols_positions <- match(mQTLcontrols_names_filtered, cpg_names_all)

########################################################
message(paste0("Run algorithm on the ", length(DerakhshanhvCpGs_positions), " hvCpG covered in 46 cells..."))

system.time(runAndSave(
   analysis = "hvCpGsMariav5",
   cpgPos_vec = DerakhshanhvCpGs_positions,
   resultDir = resDir,
   NCORES = nslots,
   p0 = 0.80,
   p1 = 0.65, 
   batch_size = 1000,
   dataDir = dataDir)
 )

###########################################################
message(paste0("Run algorithm on the ", length(mQTLcontrols_positions), " controls covered in 46 cells..."))
        
 system.time(runAndSave(
   analysis = "mQTLcontrolsv5",
   cpgPos_vec = mQTLcontrols_positions,
   resultDir = resDir,
   NCORES = nslots,
   p0 = 0.80,
   p1 = 0.65, 
   batch_size = 1,
   dataDir = dataDir, overwrite = TRUE)
   )

#######################
#### Explore results ##
##
##load("/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X/results_testLocalPC_1663CpGs_0_8p0_0_65p1.RData")
##reshvCpG = results_testLocalPC_1663CpGs_0_8p0_0_65p1
##rm(results_testLocalPC_1663CpGs_0_8p0_0_65p1)
##
##load("/home/alice/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas10X/results_testLocalPC_1474CpGs_0_8p0_0_65p1.RData")
##resControls = results_testLocalPC_1474CpGs_0_8p0_0_65p1
##rm(results_testLocalPC_1474CpGs_0_8p0_0_65p1)
##
##df = rbind(data.frame(alpha = reshvCpG, type = "hvCpG Derakhshan"),
##           data.frame(alpha = resControls, type = "mQTL controls"))
##
##ggplot(df, aes(x = type, y = alpha)) +
##  geom_jitter(aes(fill=type), pch=21, size = 3, alpha = .1)+ 
##  geom_violin(aes(col=type))+
##  geom_boxplot(aes(col=type), width = .1) + 
##  theme_minimal(base_size = 14) +                                                    
##  theme(legend.position = "none", axis.title.x = element_blank()) +
##  ylab("Probability of being a hvCpG") 
##
##summary(lm(alpha~type, df))
##
### Create matching names for GRanges
##gr_names <- paste0(
##  as.character(seqnames(DerakhshanhvCpGs_GRanges_hg38)), "_",
##  start(DerakhshanhvCpGs_GRanges_hg38), "-", end(DerakhshanhvCpGs_GRanges_hg38)
##)
##
### Match alpha values to GRanges
##DerakhshanhvCpGs_GRanges_hg38$alpha <- reshvCpG[gr_names]
##
#####################################################################
#### Compare values for Derakshan hvCpGs from arrays vs from Atlas ##
##
##load("resultsDir/Mariads/results_MariasarraysREDUCED_3samples_15datasets_6906CpGs_0_8p0_0_65p1.RData")
##Rred3array <- results_MariasarraysREDUCED_3samples_15datasets_6906CpGs_0_8p0_0_65p1
##
##Rred3array_GRanges <- GRanges(
##  seqnames = anno450k[match(rownames(Rred3array), anno450k$Name),"chr"],
##  ranges = IRanges(start = ifelse(anno450k[match(rownames(Rred3array), anno450k$Name),"strand"] %in% "+",
##                                  anno450k[match(rownames(Rred3array), anno450k$Name),"pos"],
##                                  anno450k[match(rownames(Rred3array), anno450k$Name),"pos"] - 1),
##                   end = ifelse(anno450k[match(rownames(Rred3array), anno450k$Name),"strand"] %in% "+",
##                                anno450k[match(rownames(Rred3array), anno450k$Name),"pos"] + 1,
##                                anno450k[match(rownames(Rred3array), anno450k$Name),"pos"])),
##  strand = anno450k[match(rownames(Rred3array), anno450k$Name),"strand"],
##  alpha = Rred3array)
##
##Rred3array_GRanges$alpha <- Rred3array_GRanges$alpha.alpha
##  
##Rred3array_GRanges_hg38 <- unlist(liftOver(Rred3array_GRanges, chain))
##
### Find overlapping CpGs
##hits <- findOverlaps(Rred3array_GRanges_hg38, DerakhshanhvCpGs_GRanges_hg38)
##
### Extract alpha values
##alpha_dt <- data.table(
##  alpha_x = Rred3array_GRanges_hg38$alpha[queryHits(hits)],
##  alpha_y = DerakhshanhvCpGs_GRanges_hg38$alpha[subjectHits(hits)]
##)
##
### Remove rows with NA values if needed
##alpha_dt <- na.omit(alpha_dt)
##
##cor_val <- cor(alpha_dt$alpha_x, alpha_dt$alpha_y)
##ggplot(alpha_dt, aes(x = alpha_x, y = alpha_y)) +
##  geom_point(alpha = 0.6) +
##  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
##  annotate("text", x = 0.05, y = 0.95, hjust = 0, label = paste("r =", round(cor_val, 2))) +
##  labs(
##    x = "Alpha (Rred3array)",
##    y = "Alpha (Derakhshan)"
##  ) +
##  theme_minimal(base_size = 14)
##
##
