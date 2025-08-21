prephvCpGandControls <- function(codeDir, cpg_names_all, cpg_46){
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
  
  return(list(DerakhshanhvCpGs_names_filtered = DerakhshanhvCpGs_names_filtered, 
              mQTLcontrols_names_filtered = mQTLcontrols_names_filtered))
}
