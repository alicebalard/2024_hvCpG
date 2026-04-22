library(here)

## Load libraries
if (!exists("libLoaded")) {
  source(here("05_hvCpGalgorithm", "quiet_library.R"))}

## Load functions
if (!exists("functionsLoaded")) {
  source(here("05_hvCpGalgorithm/exploreResults", "functions.R"))}


## TBC

## Bayes factor!

# Statistical test




##########################
## Test candidate genes ##
##########################
# +/-10kb of TSS
# CSGALNACT1
# chr8:19669045-20022919

# 
# ## Matt's data are in hg19
# dataMatt <- readxl::read_xlsx(here("gitignore/DEGCAGS_intersect_repeats_Alice.xlsx"))
# 
# head(dataMatt)

endo_gr <- GRanges(seqnames = paste0("chr", endo$chr),
                   ranges = IRanges(start = endo$pos, end = endo$pos), 
                   alpha = endo$alpha)
ecto_gr <- GRanges(seqnames = paste0("chr", ecto$chr),
                   ranges = IRanges(start = ecto$pos, end = ecto$pos), 
                   alpha = ecto$alpha)
meso_gr <- GRanges(seqnames = paste0("chr", meso$chr),
                   ranges = IRanges(start = meso$pos, end = meso$pos), 
                   alpha = meso$alpha)
all_gr <- GRanges(seqnames = paste0("chr", allLayers$chr),
                  ranges = IRanges(start = allLayers$pos, end = allLayers$pos), 
                  alpha = allLayers$alpha)

######################################################
## A function to output all the alphas per position ##
######################################################

getAlphasCandidates <- function(
    df, ENDO = endo_gr, ECTO = ecto_gr, MESO = meso_gr, ALL = all_gr){
  ## e.g.  df = data.frame(chromosome = "chr8", start = 19669045, end = 20022919, strand = "-", genome = "hg38")
  
  # --- Create GR object --- #
  candidate_gr <- GRanges(
    seqnames = df$chromosome,
    ranges = IRanges(start = df$start, end = df$end),
    strand = df$strand,
    genome = df$genome,
    gene = df$gene,
    alpha_endo = NA,  alpha_meso = NA,  alpha_ecto = NA,  alpha_all = NA)
  
  # --- Liftover (hg19 → hg38) IF NEEDED --- #
  if (candidate_gr$genome != "hg38"){
    ## Download and import hg19 → hg38 chain file
    chain_dir <- here("05_hvCpGalgorithm/dataPrev")
    chain_gz <- file.path(chain_dir, "hg19ToHg38.over.chain.gz")
    chain_file <- file.path(chain_dir, "hg19ToHg38.over.chain")
    
    if (!file.exists(chain_file)) {
      message("⬇️  Downloading UCSC liftOver chain file...")
      dir.create(chain_dir, showWarnings = FALSE, recursive = TRUE)
      download.file(
        url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        destfile = chain_gz,
        quiet = TRUE
      )
      R.utils::gunzip(chain_gz, destname = chain_file, remove = FALSE)
    }
    chain <- import.chain(chain_file)
    
    mapped <- liftOver(candidate_gr, chain)
    
    # Keep one-to-one mappings only
    keep <- lengths(mapped) == 1
    candidate_gr <- unlist(mapped[keep])
  }
  
  # One data frame per layer
  endo_df <- data.frame(
    candidate_seqnames = seqnames(candidate_gr)[queryHits(findOverlaps(candidate_gr, ENDO))],
    candidate_start = start(candidate_gr)[queryHits(findOverlaps(candidate_gr, ENDO))],
    candidate_end = end(candidate_gr)[queryHits(findOverlaps(candidate_gr, ENDO))],
    candidate_strand = as.character(strand(candidate_gr)[queryHits(findOverlaps(candidate_gr, ENDO))]),
    match_seqnames = seqnames(ENDO)[subjectHits(findOverlaps(candidate_gr, ENDO))],
    match_position = start(ENDO)[subjectHits(findOverlaps(candidate_gr, ENDO))],
    layer = "endo",
    alpha = ENDO$alpha[subjectHits(findOverlaps(candidate_gr, ENDO))]
  )
  
  meso_df <- data.frame(
    candidate_seqnames = seqnames(candidate_gr)[queryHits(findOverlaps(candidate_gr, MESO))],
    candidate_start = start(candidate_gr)[queryHits(findOverlaps(candidate_gr, MESO))],
    candidate_end = end(candidate_gr)[queryHits(findOverlaps(candidate_gr, MESO))],
    candidate_strand = as.character(strand(candidate_gr)[queryHits(findOverlaps(candidate_gr, MESO))]),
    match_seqnames = seqnames(MESO)[subjectHits(findOverlaps(candidate_gr, MESO))],
    match_position = start(MESO)[subjectHits(findOverlaps(candidate_gr, MESO))],
    layer = "meso",
    alpha = MESO$alpha[subjectHits(findOverlaps(candidate_gr, MESO))]
  )
  
  ecto_df <- data.frame(
    candidate_seqnames = seqnames(candidate_gr)[queryHits(findOverlaps(candidate_gr, ECTO))],
    candidate_start = start(candidate_gr)[queryHits(findOverlaps(candidate_gr, ECTO))],
    candidate_end = end(candidate_gr)[queryHits(findOverlaps(candidate_gr, ECTO))],
    candidate_strand = as.character(strand(candidate_gr)[queryHits(findOverlaps(candidate_gr, ECTO))]),
    match_seqnames = seqnames(ECTO)[subjectHits(findOverlaps(candidate_gr, ECTO))],
    match_position = start(ECTO)[subjectHits(findOverlaps(candidate_gr, ECTO))],
    layer = "ecto",
    alpha = ECTO$alpha[subjectHits(findOverlaps(candidate_gr, ECTO))]
  )
  
  all_df <- data.frame(
    candidate_seqnames = seqnames(candidate_gr)[queryHits(findOverlaps(candidate_gr, ALL))],
    candidate_start = start(candidate_gr)[queryHits(findOverlaps(candidate_gr, ALL))],
    candidate_end = end(candidate_gr)[queryHits(findOverlaps(candidate_gr, ALL))],
    candidate_strand = as.character(strand(candidate_gr)[queryHits(findOverlaps(candidate_gr, ALL))]),
    match_seqnames = seqnames(ALL)[subjectHits(findOverlaps(candidate_gr, ALL))],
    match_position = start(ALL)[subjectHits(findOverlaps(candidate_gr, ALL))],
    layer = "all",
    alpha = ALL$alpha[subjectHits(findOverlaps(candidate_gr, ALL))]
  )
  
  # Combine ALL layers
  result_df <- rbind(endo_df, meso_df, ecto_df, all_df)
  
  return(result_df)
}

result_df_CSGALNACT1 <- getAlphasCandidates(
  data.frame(chromosome = "chr8", start = 19669045, end = 20022919, strand = "-", 
             genome = "hg38", gene = "CSGALNACT1"))

ggplot(result_df, aes(x=match_position, y = alpha, col = layer))+
  geom_point()+
  theme_bw()

ggplot(result_df, aes(x=layer, y = alpha, col = layer))+
  geom_violin()+
  geom_boxplot(width = .2)+
  theme_bw()

## BS & 95% CI like in Maria's paper TO DO


## TBC


##############################
## new candidate locus Matt ##
##############################

## Matt's data are in hg19
dataMatt <- readxl::read_xlsx(here("gitignore/DEGCAGS_intersect_repeats_Alice.xlsx"))

head(dataMatt)

# --- Download and import hg19 → hg38 chain file ---
chain_dir <- here("05_hvCpGalgorithm/dataPrev")
chain_gz <- file.path(chain_dir, "hg19ToHg38.over.chain.gz")
chain_file <- file.path(chain_dir, "hg19ToHg38.over.chain")

if (!file.exists(chain_file)) {
  message("⬇️  Downloading UCSC liftOver chain file...")
  dir.create(chain_dir, showWarnings = FALSE, recursive = TRUE)
  download.file(
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    destfile = chain_gz,
    quiet = TRUE
  )
  R.utils::gunzip(chain_gz, destname = chain_file, remove = FALSE)
}
chain <- import.chain(chain_file)

# --- Liftover (hg19 → hg38) ---
dataMatt_gr <- GRanges(
  seqnames = dataMatt$chromosome,
  ranges = IRanges(start = dataMatt$start, end = dataMatt$end),
  alpha_endo = NA,  alpha_meso = NA,  alpha_ecto = NA,  alpha_all = NA,
  hg19_chr = dataMatt$chromosome, hg19_start = dataMatt$start, hg19_end = dataMatt$end,
  `%change` = dataMatt$`%change`, padj = dataMatt$padj,
  TE_chromosome = dataMatt$TE_chromosome, TE_start = dataMatt$TE_start, 
  TE_end = dataMatt$TE_end, TE_family = dataMatt$TE_family,TE_type = dataMatt$TE_type
)

mapped <- liftOver(dataMatt_gr, chain)

# Keep one-to-one mappings only
keep <- lengths(mapped) == 1
dataMatt_hg38_gr <- unlist(mapped[keep])

endo_gr <- GRanges(seqnames = paste0("chr", endo$chr),
                   ranges = IRanges(start = endo$pos, end = endo$pos), 
                   alpha = endo$alpha)
ecto_gr <- GRanges(seqnames = paste0("chr", ecto$chr),
                   ranges = IRanges(start = ecto$pos, end = ecto$pos), 
                   alpha = ecto$alpha)
meso_gr <- GRanges(seqnames = paste0("chr", meso$chr),
                   ranges = IRanges(start = meso$pos, end = meso$pos), 
                   alpha = meso$alpha)
all_gr <- GRanges(seqnames = paste0("chr", allLayers$chr),
                  ranges = IRanges(start = allLayers$pos, end = allLayers$pos), 
                  alpha = allLayers$alpha)

# Find overlapping ranges
overlaps <- findOverlaps(dataMatt_hg38_gr, endo_gr)
dataMatt_hg38_gr[queryHits(overlaps),]$alpha_endo <- endo_gr[subjectHits(overlaps),]$alpha

# Find overlapping ranges
overlaps <- findOverlaps(dataMatt_hg38_gr, meso_gr)
dataMatt_hg38_gr[queryHits(overlaps),]$alpha_meso <- meso_gr[subjectHits(overlaps),]$alpha

# Find overlapping ranges
overlaps <- findOverlaps(dataMatt_hg38_gr, ecto_gr)
dataMatt_hg38_gr[queryHits(overlaps),]$alpha_ecto <- ecto_gr[subjectHits(overlaps),]$alpha

# Find overlapping ranges
overlaps <- findOverlaps(dataMatt_hg38_gr, all_gr)
dataMatt_hg38_gr[queryHits(overlaps),]$alpha_all <- all_gr[subjectHits(overlaps),]$alpha

hist(dataMatt_hg38_gr$alpha_endo, breaks = 100)
hist(dataMatt_hg38_gr$alpha_ecto, breaks = 100)
hist(dataMatt_hg38_gr$alpha_meso, breaks = 100)
hist(dataMatt_hg38_gr$alpha_all, breaks = 100)

dataMatt_hg38_gr$is50pc <- (dataMatt_hg38_gr$alpha_endo > .5 & dataMatt_hg38_gr$alpha_meso > .5 & 
                              dataMatt_hg38_gr$alpha_ecto > .5)
dataMatt_hg38_gr$is40pc <- (dataMatt_hg38_gr$alpha_endo > .4 & dataMatt_hg38_gr$alpha_meso > .4 & 
                              dataMatt_hg38_gr$alpha_ecto > .4)

write.csv(dataMatt_hg38_gr, here("gitignore/DEGCAGS_intersect_repeats_Alice_withalpha.csv"), 
          row.names = F, quote = F)




