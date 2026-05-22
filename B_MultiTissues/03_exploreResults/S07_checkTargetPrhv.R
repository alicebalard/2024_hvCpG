
##############################
## new candidate locus Matt ##
##############################

## Matt's data are in hg19
dataMatt <- readxl::read_xlsx(here("gitignore/DEGCAGS_intersect_repeats_Alice.xlsx"))

# --- Download and import hg19 → hg38 chain file ---
chain_dir <- here("B_MultiTissues/dataIn")
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

# Initialise the column with NA first
dataMatt_hg38_gr$alpha_geomean <- NA_real_
# Find overlapping ranges

# Check and fix chromosome names
if (!any(grepl("chr", as.character(seqlevels(table3layers))))) {
  seqlevels(table3layers) <- paste0("chr", seqlevels(table3layers))
}

overlaps <- findOverlaps(dataMatt_hg38_gr, table3layers)
dataMatt_hg38_gr[queryHits(overlaps),]$alpha_geomean <- 
  table3layers[subjectHits(overlaps),]$alpha_geomean

ggplot(as.data.frame(dataMatt_hg38_gr), aes(x = seqnames, y=alpha_geomean)) +
  geom_violin() +
  geom_boxplot(width = .2) +
  geom_jitter() +
  theme_minimal(base_size = 14) 
