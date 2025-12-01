# Load quiet library and required packages
source(here("05_hvCpGalgorithm/quiet_library.R"))

message("🧬 Building dictionary linking 450k and EPIC probes to hg19 and hg38 coordinates...")

# --- 1️⃣ Load 450k annotation ---
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno450k <- data.frame(
  CpG = rownames(anno450k),
  chr = anno450k$chr,
  pos = anno450k$pos,
  gene = anno450k$UCSC_RefGene_Name,
  relation_to_island = anno450k$Relation_to_Island,
  strand = anno450k$strand,
  stringsAsFactors = FALSE
)

# --- 2️⃣ Load EPIC annotation ---
annoEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annoEPIC <- data.frame(
  CpG = rownames(annoEPIC),
  chr = annoEPIC$chr,
  pos = annoEPIC$pos,
  gene = annoEPIC$UCSC_RefGene_Name,
  relation_to_island = annoEPIC$Relation_to_Island,
  strand = annoEPIC$strand,
  stringsAsFactors = FALSE
)

# --- 3️⃣ Combine both arrays ---
anno_combined <- bind_rows(
  anno450k %>% mutate(array = "450k"),
  annoEPIC %>% mutate(array = "EPIC")
)

anno_combined$chrpos_hg19 <- paste0(anno_combined$chr, "_", anno_combined$pos)

# --- 4️⃣ Download and import hg19 → hg38 chain file ---
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

# --- 5️⃣ Liftover (hg19 → hg38) ---
message("🔄 Performing liftover to hg38...")

hg19_gr <- GRanges(
  seqnames = anno_combined$chr,
  ranges = IRanges(start = anno_combined$pos, width = 1),
  names = anno_combined$CpG
)

mapped <- liftOver(hg19_gr, chain)

# Keep one-to-one mappings only
keep <- lengths(mapped) == 1
hg38_gr <- unlist(mapped[keep])

# --- 6️⃣ Creat9 and hg38 coordinate format ---

dico <- anno_combined[keep, ]

# Keep both hg19 and add hg38 coordinates
dico$chr_hg38 <- as.character(seqnames(hg38_gr))
dico$pos_hg38 <- start(hg38_gr)

# chr_pos string for hg38
dico$chrpos_hg38 <- paste0(dico$chr_hg38, "_", dico$pos_hg38)

# Find the ones present in both array
dup <- dico$chrpos_hg38[duplicated(dico$chrpos_hg38)]

dico[dico$chrpos_hg38 %in% dup,"array"] <- "450k and EPIC"
dico <- unique(dico)

table(dico$array)
# 450k     450k and EPIC      EPIC 
# 33041        452305        413313 

# --- 7️⃣R m junk ---
rm("anno_combined", "anno450k", "annoEPIC", "chain_dir", "chain_file", "chain_gz",
   "dup","hg19_gr", "hg38_gr", "keep", "mapped")
