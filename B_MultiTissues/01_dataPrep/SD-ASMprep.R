#####################################################################
## Prepare
library(here)
## Load libraries
source(here("B_MultiTissues", "quiet_library.R"))
#####################################################################

createSDASM_GR <- function(){
  sup5Rosenski2026 <- readxl::read_excel(
    here("B_MultiTissues/dataIn/41467_2026_71693_MOESM3_ESM.xlsx"), 
    sheet = 4, skip = 2) |> as.data.table()
  
  length(unique(sup5Rosenski2026$region))
  # 33574 bimodal regions with 1 or more SNPs
  
  table(sup5Rosenski2026$`SD-ASM classification`)
  # denovo meth.    other  tissue-specific demeth.           ubiq. 
  # 1435            22847                25714               361 
  
  # 1. unique (region, classification) combinations
  reduced <- unique(sup5Rosenski2026[, .(region, `SD-ASM classification`)])
  setnames(reduced, "SD-ASM classification", "classification")
  
  # 2. check counts against the paper
  table(reduced$classification)
  # denovo meth.      other    tissue-specific demeth.       ubiq. 
  # 1173              15546                   19548          166
  
  ## Make GRanges object from regions
  gr <- GRanges(sup5Rosenski2026$region)
  mcols(gr)$classification <- sup5Rosenski2026$`SD-ASM classification`
  return(unique(gr))
}

SDASM_GR <- createSDASM_GR()

## ─────────────────────────────────────────────────────────────────────────────
## Append to SD-ASMprep.R
## Expand SD-ASM REGIONS -> individual CpG positions ("chr_pos"), so the list can
## be fed to the beta-matrix pipeline via --exclude_sites (which excludes by
## chr_pos, not by interval). Uses the SAME CpG.bed the pipeline uses, so the
## coordinates line up exactly.
## ─────────────────────────────────────────────────────────────────────────────

writeSDASM_blacklist <- function(sdasm_gr,
                                 cpg_bed,            # same CpG.bed.gz the pipeline uses
                                 out_file,
                                 classifications = NULL) {  # NULL = ALL SD-ASM regions
  stopifnot(file.exists(cpg_bed))
  
  # optionally restrict to certain SD-ASM classes (default: all)
  gr <- sdasm_gr
  if (!is.null(classifications)) {
    gr <- gr[mcols(gr)$classification %in% classifications]
  }
  gr <- GenomicRanges::reduce(gr)          # collapse overlapping/abutting regions
  
  # load reference CpG positions as a GRanges of width-1 sites
  bed <- data.table::fread(cpg_bed, header = FALSE)      # cols: chr, pos (0- or 1-based? see note)
  data.table::setnames(bed, 1:2, c("chr", "pos"))
  cpg_gr <- GenomicRanges::GRanges(bed$chr,
                                   IRanges::IRanges(bed$pos, bed$pos))
  
  # CpGs falling inside any SD-ASM region
  hit <- GenomicRanges::overlapsAny(cpg_gr, gr, ignore.strand = TRUE)
  black <- paste0(bed$chr[hit], "_", bed$pos[hit])
  
  writeLines(black, out_file)
  message(sprintf("Wrote %s SD-ASM CpGs (of %s regions) -> %s",
                  format(length(black), big.mark = ","),
                  format(length(gr),    big.mark = ","),
                  out_file))
  invisible(black)
}

## Example: write the FULL SD-ASM blacklist (all classifications)
## Point cpg_bed at the SAME file used in run_pipeline_atlas.sh (--cpg_bed).
writeSDASM_blacklist(
  SDASM_GR,
  cpg_bed  = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/hg38/CpG.bed.gz",
  out_file = here::here("gitignore/sdasm_all_chr_pos.txt")
)

## Optional: exclude ONLY specific classes (e.g. keep the confound test symmetric)
# writeSDASM_blacklist(SDASM_GR, cpg_bed = "/.../CpG.bed.gz",
#   out_file = here::here("gitignore/sdasm_ubiq_chr_pos.txt"),
#   classifications = "ubiq.")
