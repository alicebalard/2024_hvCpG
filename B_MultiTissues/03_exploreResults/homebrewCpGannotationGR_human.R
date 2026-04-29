#####################################################################
## Annotate a vector of CpGs using a bed12 and a gff3 annotation files
#####################################################################

# Usage:
# myHomebrewCpGannotation(CpGvec = c("chr1_26116565", "chr2_2635914"),
#                         myannotBed12 = annotBed12, myannotGff3 = annotGff3)

## Issue in the original "annotateWithGeneParts" from genomation function:
## involved associating all CpGs by proximity to a TSS first, then assigning
## gene ID/feature type. However, if a CpG is at the end of a large gene, it could
## be closer to the TSS of the next gene along, so would get associated incorrectly.

## Solution: 2 cases, one for the genic part, overlapping directly with the annotation
## via intersection; then the TSS proximity method is only applied for intergenic sites.

## Charley Yen's (https://github.com/eugeniecyen) improvements on my original function:
## - removing a filter by distance <10kb from TSS before splitting into the 2 cases,
## as this might also remove some genic sites on large genes.
## - having to set unique.prom=FALSE for promoters to come up in readTranscriptFeatures
## to load bed12 annotation.
## - for some reason, the gene ID association to DMS was getting jumbled up in the
## final output table. This seemed to be occurring with the old looping code, but
## no longer with the adapted homebrew function

# set rerun = T when change list to reload
myHomebrewCpGannotation <- function(CpGvec, myannotBed12, myannotGff3, rerun=F){
  myGRanges = getGRange_corrected(CpGvec, myannotBed12)
  message(paste0("We have ", length(myGRanges), " positions extracted from the bed12 file."))
  return(myGRanges)
}

getGRange_corrected <- function(CpGvec, myannotBed12){
  # Change the vector into a GRange:
  GRangeOBJ = makeGRangesFromDataFrame(data.frame(chr=sapply(strsplit(CpGvec, "_"), `[`, 1),
                                                  start=sapply(strsplit(CpGvec, "_"), `[`, 2),
                                                  end=sapply(strsplit(CpGvec, "_"), `[`, 2),
                                                  DMS=CpGvec), keep.extra.columns = T)
  
  # Annotate with features
  GRangeOBJ_annot = genomation::annotateWithGeneParts(target = as(GRangeOBJ,"GRanges"),
                                          feature = myannotBed12)
  
  ## We assign the feature type to GRangeOBJ
  GRangeOBJ$featureType = ifelse(GRangeOBJ_annot@members[,1]==1, "promoters",
                                 ifelse(GRangeOBJ_annot@members[,2]==1, "exons",
                                        ifelse(GRangeOBJ_annot@members[,3]==1, "introns", "intergenic")))
  
  ## CASE 1: genic --> intersection with BED12 file
  GRangeOBJ1=GRangeOBJ[!GRangeOBJ$featureType %in% "intergenic"]
  
  # Add empty column called geneInfo to fill with gene IDs
  GRangeOBJ1$feature.name <- NA
  
  add_geneInfo_genic <- function(x, GRangeOBJ, annotBed12){
    ov = GenomicRanges::findOverlaps(
      annotBed12[[x]],
      GRangeOBJ[GRangeOBJ$featureType %in% x,])
    ## Add gene annotation to subject GRanges (i.e. left join)
    mcols(GRangeOBJ[GRangeOBJ$featureType %in% x,])[subjectHits(ov), "feature.name"] =
      mcols(annotBed12[[x]])[queryHits(ov), "name"]
    return(GRangeOBJ)
  }
  
  GRangeOBJ_ex = add_geneInfo_genic("exons", GRangeOBJ1, myannotBed12) # Add gene names for sites on exons
  GRangeOBJ_ex_in = add_geneInfo_genic("introns", GRangeOBJ_ex, myannotBed12) # Add gene names for sites on introns
  GRangeOBJ_genic = add_geneInfo_genic("promoters", GRangeOBJ_ex_in, myannotBed12) # Add gene names for sites on promoters
  
  ## CASE 2: intergenic --> original method
  GRangeOBJ2=GRangeOBJ[GRangeOBJ$featureType %in% "intergenic"]
  
  ### Extract annotation from the overall annotation object
  GRangeOBJ_annot_2 = GRangeOBJ_annot@dist.to.TSS[GRangeOBJ$featureType %in% "intergenic",]
  
  GRangeOBJ_annot_2$feature.name
  
  ## Filter out the DMS non associated with a gene:
  # Heckwolf 2020: To be associated to a gene, the DMS had to be either inside the gene or,
  # if intergenic, not further than 10 kb away from the TSS.
  rows2rm = which((GRangeOBJ2$dist.to.TSS>10000 |
                     GRangeOBJ2$dist.to.TSS < -10000) &
                    GRangeOBJ2$featureType %in% "intergenic")
  if (!sjmisc::is_empty(rows2rm)){
    GRangeOBJ2 = GRangeOBJ2[-rows2rm,]
    GRangeOBJ_annot_2 = GRangeOBJ_annot_2[-rows2rm,]
  }
  
  GRangeOBJ_intergenic = GRangeOBJ2
  GRangeOBJ_intergenic$feature.name = GRangeOBJ_annot_2$feature.name
  
  ## Merge both cases
  GRangeOBJ_both = c(GRangeOBJ_genic, GRangeOBJ_intergenic)
  
  # Change recursively the gene names to keep only ID
  getName <- function(x) {sub(";.*", "", sub(".*ID=", "", x))}
  
  for (i in 1:length(GRangeOBJ_both)){
    GRangeOBJ_both$feature.name[i] <- getName(GRangeOBJ_both$feature.name[i])
  }
  
  return(GRangeOBJ_both)
}
