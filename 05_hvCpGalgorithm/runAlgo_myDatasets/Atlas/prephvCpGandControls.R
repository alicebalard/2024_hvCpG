prephvCpGandControls <- function(codeDir, cpg_46=FALSE){
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
  
  message("Make a dictionary to link 450k probe names to hg19 and hg38 positions...")
  options(scipen = 999)   # large penalty against scientific notation
  
  # Preserve probe IDs when building hg19_GRanges
  hg19_GRanges <- GRanges(
    seqnames = anno450k$chr,
    ranges = IRanges(start = anno450k$pos, width = 1),
    strand = anno450k$strand,
    names = anno450k$Name   # <-- keep probe names
  )
  
  mapped <- liftOver(hg19_GRanges, chain)
  keep <- lengths(mapped) == 1
  
  hg19_unique <- hg19_GRanges[keep]
  hg38_unique <- unlist(mapped[keep])
  
  # Subset anno450k by probe IDs directly
  anno450k_unique <- anno450k[anno450k$Name %in% hg19_unique$names, ]
  
  dictionary <- data.frame(
    illu450k = anno450k_unique$Name,
    hg19 = paste0(
      hg19_unique@seqnames, "_",
      ifelse(hg19_unique@strand %in% "+", hg19_unique@ranges@start, hg19_unique@ranges@start -1), "-",
      ifelse(hg19_unique@strand %in% "+", hg19_unique@ranges@start +1, hg19_unique@ranges@start)),
    hg38 = paste0(
      hg38_unique@seqnames, "_",
      ifelse(hg38_unique@strand %in% "+", hg38_unique@ranges@start, hg38_unique@ranges@start -1), "-",
      ifelse(hg38_unique@strand %in% "+", hg38_unique@ranges@start +1, hg38_unique@ranges@start)))
  
  data <- read.table(file.path(codeDir, "03_prepDatasetsMaria/cistrans_GoDMC_hvCpG_matched_control.txt"), header = T)
  message("Prepare Derakhshan2022_hvCpGs...")
  DerakhshanhvCpGs_names <- dictionary[match(data$hvCpG_name, dictionary$illu450k), "hg38"]
  
  message("Matching genetic controls to hvCpGs..")
  mQTLcontrols_names <- dictionary[match(data$controlCpG_name, dictionary$illu450k), "hg38"]
  
  if (cpg_46 == FALSE){
    return(list(chain = chain, dictionary = dictionary,
                DerakhshanhvCpGs_names = DerakhshanhvCpGs_names, 
                mQTLcontrols_names = mQTLcontrols_names))
  } else {
    # Restrict to CpGs in cpg_46
    DerakhshanhvCpGs_names_filtered <- DerakhshanhvCpGs_names[DerakhshanhvCpGs_names %in% cpg_46]
    mQTLcontrols_names_filtered <- mQTLcontrols_names[mQTLcontrols_names %in% cpg_46]
    
    return(list(chain = chain, dictionary = dictionary, 
                DerakhshanhvCpGs_names_filtered = DerakhshanhvCpGs_names_filtered, 
                mQTLcontrols_names_filtered = mQTLcontrols_names_filtered))
  }
}
