prephvCpGandControls <- function(codeDir, cpg_names_all, cpg_46){

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
  
  data <- read.table(file.path(codeDir, "03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt"), header = T)
  
  message("Make a dictionary to link 450k probe names to hg19 and hg38 positions...")
  vecCpGs <- unlist(data); names(vecCpGs) = NULL
  idx <- match(vecCpGs, anno450k$Name)
  hg19_GRanges <- GRanges(
    seqnames = anno450k$chr[idx],
    ranges = IRanges(start = anno450k$pos[idx], width = 1),
    strand = anno450k$strand[idx]
  )
  hg38_GRanges <- unlist(liftOver(hg19_GRanges, chain))
  
  dictionary <- data.frame(
    illu450k = vecCpGs,
    hg19 = paste0(
      hg19_GRanges@seqnames, "_",
      ifelse(hg19_GRanges@strand %in% "+", hg19_GRanges@ranges@start, hg19_GRanges@ranges@start -1), "-",
      ifelse(hg19_GRanges@strand %in% "+", hg19_GRanges@ranges@start +1, hg19_GRanges@ranges@start)),
    hg38 = paste0(
      hg38_GRanges@seqnames, "_",
      ifelse(hg38_GRanges@strand %in% "+", hg38_GRanges@ranges@start, hg38_GRanges@ranges@start -1), "-",
      ifelse(hg38_GRanges@strand %in% "+", hg38_GRanges@ranges@start +1, hg38_GRanges@ranges@start)))
  
  message("Prepare Derakhshan2022_hvCpGs...")
  DerakhshanhvCpGs_names <- dictionary[match(data$hvCpG_name, dictionary$illu450k), "hg38"]
  # Restrict to CpGs in cpg_46
  DerakhshanhvCpGs_names_filtered <- DerakhshanhvCpGs_names[DerakhshanhvCpGs_names %in% cpg_46]
  
  message("Matching genetic controls to hvCpGs..")
  mQTLcontrols_names <- dictionary[match(data$controlCpG_name, dictionary$illu450k), "hg38"]
  # Restrict to CpGs in cpg_46
  mQTLcontrols_names_filtered <- mQTLcontrols_names[mQTLcontrols_names %in% cpg_46]
  
  return(list(dictionary = dictionary,
              DerakhshanhvCpGs_names_filtered = DerakhshanhvCpGs_names_filtered, 
              mQTLcontrols_names_filtered = mQTLcontrols_names_filtered))
}
