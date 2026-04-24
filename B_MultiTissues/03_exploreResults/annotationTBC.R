## For a vector of CpGs in the format chromosome_position 
## CpGvec <- c("chr1_17452", "chr1_17478", "chr1_17483")
## 1. Keep CpGs in regions where at least 5 CpGs are in 50bp distance to each other
## 2. annotate with genes
## 3. run GO term enrichment with clusterProfiler::enrichGO


library(data.table)
# library(GenomicRanges)
# library(ensembldb)
# library(EnsDb.Hsapiens.v105)  
library(AnnotationHub)
ah <- AnnotationHub()
# snapshotDate(): 2025-04-08
query(ah, c("EnsDb", "v105"))

## Query for available H.Sapiens EnsDb databases
ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb"))

genes(ahDb)


library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

## STEP 1 — ULTRA‑FAST DENSE CpG CLUSTERING (5 CpGs min within 50 bp)
clusterCpGs <- function(CpGvec, max_gap = 50, min_size = 5) {
  dt <- data.table(
    raw = CpGvec,
    chr = sub("_.*", "", CpGvec),
    pos = as.integer(sub(".*_", "", CpGvec))
  )
  setkey(dt, chr, pos)
  
  # gap to previous CpG
  dt[, gap := pos - data.table::shift(pos), by = chr]
  
  # run ID increments whenever gap > max_gap OR different chromosome
  dt[, run_id := cumsum(is.na(gap) | gap > max_gap), by = chr]
  
  # Count CpGs in each run
  dt[, run_size := .N, by = .(chr, run_id)]
  
  # Keep only large runs
  dt[run_size >= min_size, raw]
}

test <- clusterCpGs(head(stable, 1000))

## STEP 2 — FAST GENE ANNOTATION (OFFLINE)


library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

annotateCpGs_txdb <- function(CpGs, tss_window = 10000) {
  
  if (length(CpGs) == 0) return(character(0))
  
  chr <- sub("_.*", "", CpGs)
  pos <- as.integer(sub(".*_", "", CpGs))
  gr  <- GRanges(chr, IRanges(pos, pos))
  
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  txdb
  # genes_txdb <- genes(txdb)
  # promoters_txdb <- promoters(txdb, upstream = tss_window, downstream = tss_window)
  # 
  # # overlaps with gene bodies
  # o1 <- findOverlaps(gr, genes_txdb)
  # g1 <- genes_txdb$gene_id[subjectHits(o1)]
  # 
  # # overlaps with promoters
  # o2 <- findOverlaps(gr, promoters_txdb)
  # g2 <- promoters_txdb$gene_id[subjectHits(o2)]
  # 
  # unique(c(g1, g2))
}

test2 <- annotateCpGs_txdb(test)

## STEP 3 — GO ENRICHMENT (clusterProfiler)

runGO <- function(entrez_ids, universe = NULL) {
  enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    keyType       = "ENTREZID",
    universe      = universe,
    pAdjustMethod = "BH",
    readable      = TRUE
  )
}

## 🚀 FULL PIPELINE FUNCTION 

CpG_GO_pipeline <- function(CpGvec,
                            max_gap = 50, min_size = 5,
                            tss_window = 10000,
                            universe = NULL) {
  
  message("Clustering CpGs...")
  CpGclustered <- clusterCpGs(CpGvec, max_gap, min_size)
  message(sprintf("Reduced from %d to %d clustered CpGs",
                  length(CpGvec), length(CpGclustered)))
  
  if (length(CpGclustered) == 0) {
    warning("No CpG clusters found.")
    return(NULL)
  }
  
  message("Annotating genes...")
  ensg <- annotateCpGs(CpGclustered, tss_window)
  
  message(sprintf("Found %d ENSEMBL genes", length(ensg)))
  
  message("Converting to Entrez IDs...")
  entrez <- convertToEntrez(ensg)
  message(sprintf("Converted to %d Entrez IDs", length(entrez)))
  
  message("Running GO enrichment...")
  runGO(entrez, universe)
}


  
  
  
  