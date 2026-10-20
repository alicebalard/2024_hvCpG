#!/usr/bin/env Rscript
###############################################################################
# build_reference_jaffe.R
#
# One-off preparation of the reference file used by ont_cell_proportions.R.
#
# It takes the 500 CpG coefficients of the Jaffe model for adult venous blood
# (FlowSorted.Blood.450k.JaffeModelPars; CD8T, CD4T, NK, Bcell, Mono, Gran),
# looks up the genomic position of each CpG from the Illumina 450k annotation
# (hg19), converts the positions to hg38 with a UCSC chain file if required,
# and writes a table in the layout that ont_cell_proportions.R reads.
#
# The 450k annotation gives one position per probe, and it is not always clear
# whether that position is the C of the CpG or the G on the opposite strand.
# If a genome package is supplied with --genome_pkg, the script checks the
# sequence and records the position of the C on the forward strand for every
# CpG, and drops any site that is not a CpG in the target genome. This check is
# strongly recommended.
#
# USAGE
#   Rscript build_reference_jaffe.R --build hg38 \
#       --chain hg19ToHg38.over.chain \
#       --genome_pkg BSgenome.Hsapiens.UCSC.hg38 \
#       --output jaffe_hg38_reference.tsv
#
#   Options (defaults in brackets)
#     --build       hg38 or hg19                                          [hg38]
#     --chain       UCSC chain file, hg19ToHg38.over.chain; needed for hg38
#     --genome_pkg  BSgenome package matching --build, for the CpG check
#     --output      file to write                  [jaffe_reference.tsv]
#
# REQUIREMENTS (Bioconductor)
#   FlowSorted.Blood.450k, IlluminaHumanMethylation450kanno.ilmn12.hg19, minfi,
#   GenomicRanges, rtracklayer, and the BSgenome package if --genome_pkg is used.
#   The chain file can be downloaded from
#   https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
#
# OUTPUT
#   A tab-separated table with the columns cpg, chrom, cpg_start, followed by
#   one column per cell type. A short account of how many CpGs were kept at
#   each step is printed to the console.
#
# NOTE ON USE
#   The Jaffe coefficients come from sorted blood of adults measured on the
#   450k array. They are a reasonable starting reference for adult peripheral
#   blood. Their suitability for children, for populations of African ancestry
#   and for Nanopore-derived methylation should be checked against an
#   independent estimate for a subset of samples.
###############################################################################

suppressPackageStartupMessages({
  library(FlowSorted.Blood.450k)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(minfi)
  library(GenomicRanges)
  library(rtracklayer)
})

# --- command line ----------------------------------------------------------

defaults <- list(build = "hg38",
                 chain = NA_character_,
                 genome_pkg = NA_character_,
                 output = "jaffe_reference.tsv")

a <- commandArgs(trailingOnly = TRUE)
opts <- defaults
i <- 1
while (i <= length(a)) {
  key <- sub("^--", "", a[i])
  if (!grepl("^--", a[i]) || !(key %in% names(defaults))) stop("Unknown option: ", a[i])
  if (i == length(a)) stop("Option --", key, " needs a value")
  opts[[key]] <- a[i + 1]
  i <- i + 2
}
if (!(opts$build %in% c("hg38", "hg19"))) stop("--build must be hg38 or hg19")
if (opts$build == "hg38" && is.na(opts$chain)) stop("--chain is required for --build hg38")

# --- coefficients and hg19 positions ----------------------------------------

coef <- as.matrix(FlowSorted.Blood.450k.JaffeModelPars)
message("Jaffe coefficients: ", nrow(coef), " CpGs, cell types: ",
        paste(colnames(coef), collapse = ", "))

ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
absent <- setdiff(rownames(coef), rownames(ann))
if (length(absent) > 0) stop(length(absent), " CpGs are missing from the 450k annotation")
ann <- ann[rownames(coef), ]

sites <- data.frame(cpg = rownames(coef),
                    chrom = as.character(ann$chr),
                    pos = as.integer(ann$pos),
                    stringsAsFactors = FALSE)

# --- convert to hg38 if required ---------------------------------------------

if (opts$build == "hg38") {
  gr <- GRanges(seqnames = sites$chrom,
                ranges = IRanges(start = sites$pos, width = 1),
                strand = "+",
                cpg = sites$cpg)
  lifted <- liftOver(gr, import.chain(opts$chain))

  # Keep CpGs that map to exactly one place, on the same strand orientation.
  one_hit <- lengths(lifted) == 1L
  message("Mapped uniquely to hg38: ", sum(one_hit), " of ", length(gr))
  lifted <- unlist(lifted[one_hit])
  same_orientation <- as.character(strand(lifted)) == "+"
  message("Same orientation after conversion: ", sum(same_orientation),
          " of ", length(lifted))
  lifted <- lifted[same_orientation]

  sites <- data.frame(cpg = as.character(mcols(lifted)$cpg),
                      chrom = as.character(seqnames(lifted)),
                      pos = start(lifted),
                      stringsAsFactors = FALSE)
}

# --- check against the genome sequence ---------------------------------------

if (!is.na(opts$genome_pkg)) {
  if (!requireNamespace(opts$genome_pkg, quietly = TRUE)) {
    stop("Package not installed: ", opts$genome_pkg)
  }
  genome <- getExportedValue(opts$genome_pkg, "Hsapiens")

  sites <- sites[sites$chrom %in% seqlevels(genome), ]
  dinucleotide <- function(chrom, start) {
    toupper(as.character(Biostrings::getSeq(
      genome, GRanges(chrom, IRanges(start = start, width = 2)))))
  }
  at_pos    <- dinucleotide(sites$chrom, sites$pos)
  one_before <- dinucleotide(sites$chrom, sites$pos - 1L)

  sites$cpg_start <- ifelse(at_pos == "CG", sites$pos,
                     ifelse(one_before == "CG", sites$pos - 1L, NA_integer_))
  message("Annotated position is the C of a CpG: ", sum(at_pos == "CG"),
          "; the G of a CpG: ", sum(at_pos != "CG" & one_before == "CG"),
          "; not a CpG in this genome (dropped): ", sum(is.na(sites$cpg_start)))
  sites <- sites[!is.na(sites$cpg_start), ]
} else {
  warning("No genome supplied, so the annotated position is taken to be the C of ",
          "the CpG without checking.")
  sites$cpg_start <- sites$pos
}

# --- write -------------------------------------------------------------------

out <- data.frame(cpg = sites$cpg,
                  chrom = sites$chrom,
                  cpg_start = sites$cpg_start,
                  coef[sites$cpg, , drop = FALSE],
                  check.names = FALSE,
                  stringsAsFactors = FALSE)
write.table(out, opts$output, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", nrow(out), " CpGs to ", opts$output)
