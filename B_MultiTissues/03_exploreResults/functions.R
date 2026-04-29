############################################
# loads functions useful in multiple scripts

## makeVennArrayReduced
## prepAtlasdt
## plotManhattanFromdt
## makeZ_inner
## makeCompPlot
## makeGRfromMyCpGPos
## GO old functions (commented out)
## .safe_fisher & test_enrichment_quadrants --> test enrichment of target CpGs 
## for each quadrant vs the other three combined
### plotMyVenn: Compute overlap across any number of groups and plot a Venn diagram

makeVennArrayReduced <- function(df_circles){
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = df_circles,
                          aes(x, y, group = group),
                          fill = "#4981BF", alpha = 0.3,
                          colour = "grey40", linewidth = 0.5) +
    # Exclusive regions
    ggplot2::annotate("text", x = -0.9,  y =  0.55, label = fmt(n_only_full,      total), size = 3, lineheight = 0.9) +
    ggplot2::annotate("text", x =  0.9,  y =  0.55, label = fmt(n_only_2ind,      total), size = 3, lineheight = 0.9) +
    ggplot2::annotate("text", x =  0.0,  y = -1.05, label = fmt(n_only_3ind,      total), size = 3, lineheight = 0.9) +
    # Pairwise only
    ggplot2::annotate("text", x =  0.0,  y =  0.75, label = fmt(n_full_2ind_only, total), size = 3, lineheight = 0.9) +
    ggplot2::annotate("text", x = -0.55, y = -0.25, label = fmt(n_full_3ind_only, total), size = 3, lineheight = 0.9) +
    ggplot2::annotate("text", x =  0.55, y = -0.25, label = fmt(n_2ind_3ind_only, total), size = 3, lineheight = 0.9) +
    # Triple
    ggplot2::annotate("text", x =  0.0,  y =  0.15, label = fmt(v$n_all,          total), size = 3, fontface = "bold", lineheight = 0.9) +
    # Set labels with total size
    ggplot2::annotate("text", x = -1.1,  y =  1.5,  label = paste0("Full array\nn=",    v$n_full), size = 3.5, fontface = "bold") +
    ggplot2::annotate("text", x =  1.1,  y =  1.5,  label = paste0("Array 2ind/ds\nn=", v$n_2ind), size = 3.5, fontface = "bold") +
    ggplot2::annotate("text", x =  0.0,  y = -1.7,  label = paste0("Array 3ind/ds\nn=", v$n_3ind), size = 3.5, fontface = "bold") +
    ggplot2::coord_fixed() +
    ggplot2::theme_void()
}

prepAtlasdt <- function(dir = "Atlas10X"){
  # Define parent folder containing all "Atlas_batchXXX" folders
  parent_dir = here(paste0("05_hvCpGalgorithm/resultsDir/Atlas/", dir))
  
  # Get list of relevant RDS files
  rds_files <- dir(parent_dir, pattern = ".rds$", 
                   recursive = TRUE, full.names = TRUE)
  rds_files
  all_cpg_values <- numeric()
  pb <- progress_bar$new(total = length(rds_files), format = "📦 :current/:total [:bar] :percent")
  
  for (file in rds_files) {
    obj <- readRDS(file)
    if (is.matrix(obj)) obj <- obj[, 1]
    all_cpg_values <- c(all_cpg_values, obj)
    pb$tick()
  }
  
  # Create data.table from named vector
  dt <- data.table(
    name =   sub("-[0-9]+$", "", names(all_cpg_values)), # just keep the C position instead of C + 1
    alpha = as.numeric(all_cpg_values)
  )

  #######################################################################
  # Parse "chr_pos" in name into chr, start_pos, end_pos. NB: takes a couple of minutes
  dt[, c("chr", "pos") := tstrsplit(name, "_", fixed = TRUE)]
  
  # Convert to integer/numeric if not already
  dt[, pos := as.integer(pos)]
  
  ## Check chromosomes present:
  message("Chromosomes in the dataset:")
  table(unique(dt$chr))
  
  # Convert chr from "chrN" to  factor
  dt[, chr := sub("chr", "", chr)]
  dt[, chr := factor(chr, levels = as.character(c(1:22, "X", "Y", "M")))]
  
  ## Check chromosomes order:
  message("Chromosomes in the dataset:")
  table(unique(dt$chr))
  
  ## Mark group membership in dt
  dt[, group := NA_character_]
  dt[name %in% DerakhshanhvCpGs_hg38, group := "hvCpG_Derakhshan"]
  dt[name %in% mQTLcontrols_hg38, group := "mQTLcontrols"]
  
  # Compute cumulative position offsets for Manhattan plot
  setorder(dt, chr, pos)
  
  offsets <- dt[, .(max_pos = max(pos, na.rm = TRUE)), by = chr]
  offsets[, cum_offset := c(0, head(cumsum(as.numeric(max_pos)), -1))]
  
  dt <- merge(dt, offsets[, .(chr, cum_offset)], by = "chr", all.x = TRUE, sort = FALSE)
  
  # Convert to integer/numeric if not already
  dt[, cum_offset := as.numeric(cum_offset)]
  dt[, pos2 := pos + cum_offset]
  
  return(dt)
}

plotManhattanFromdt <- function(dt, transp = 0.01){
  # Compute chromosome centers for x-axis labeling
  df2 <- dt[, .(center = mean(range(pos2, na.rm = TRUE))), by = chr]
  df2 <- merge(data.frame(chr = factor(c(1:22, "X", "Y", "M"), levels=as.character(c(1:22, "X", "Y", "M")))),
               df2, by = "chr", all.x = TRUE, sort = TRUE)
  df2 <- na.omit(df2)
  
  # Compute chromosome boundaries
  df_bounds <- dt[, .(min_pos = min(pos2, na.rm = TRUE), 
                      max_pos = max(pos2, na.rm = TRUE)), by = chr]
  
  # Midpoints between chromosomes = where to draw dotted lines
  df_bounds[, next_start := data.table::shift(min_pos, n = 1, type = "lead")]
  vlines <- df_bounds[!is.na(next_start), .(xintercept = (max_pos + next_start)/2)]
  
  ggplot() +
    # background cloud
    geom_point_rast(data = dt[is.na(group)], 
                    aes(x = pos2, y = alpha),
                    color = "black", size = 0.01, alpha = transp, raster.dpi = 72) +
    # hvCpG highlights
    geom_point(data = dt[group == "hvCpG_Derakhshan"],
               aes(x = pos2, y = alpha),
               color = "#DC3220", size = 1, alpha = 0.7) +
    # mQTL controls highlights
    geom_point(data = dt[group == "mQTLcontrols"],
               aes(x = pos2, y = alpha),
               color = "#005AB5", size = 1, alpha = 0.7) +
    # Add dotted separators
    geom_vline(data = vlines, aes(xintercept = xintercept),
               linetype = 3, color = "grey60") +
    theme_classic() + theme(legend.position = "none") +
    scale_x_continuous(breaks = df2$center, labels = as.character(df2$chr), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Chromosome", y = "Probability of being a hvCpG")+
    theme_minimal(base_size = 14)
}

# ---- Inner helper for makeCompPlot: build Z_inner ----
makeZ_inner <- function(X, Y, whichAlphaX = NULL, whichAlphaY = NULL) {
  setDT(X); setDT(Y)
  
  # Determine X side
  is_array_X <- any(grepl("array", names(X)))
  if (is_array_X) {
    if (is.null(whichAlphaX)) {
      stop("X looks like an *array* table (columns contain 'array'). ",
           "Please provide whichAlphaX, e.g. 'alpha_array_all'.")
    }
    colX <- if (is.character(whichAlphaX)) whichAlphaX else deparse(substitute(whichAlphaX))
    stopifnot("chrpos" %in% names(X), colX %in% names(X))
    X <- X[, .(name = chrpos, alpha_X = get(colX))]
  } else {
    stopifnot(all(c("name", "alpha") %in% names(X)))
    X <- X[, .(name, alpha_X = alpha)]
  }
  
  # Determine Y side
  is_array_Y <- any(grepl("array", names(Y)))
  if (is_array_Y) {
    if (is.null(whichAlphaY)) {
      stop("Y looks like an *array* table (columns contain 'array'). ",
           "Please provide whichAlphaY.")
    }
    colY <- if (is.character(whichAlphaY)) whichAlphaY else deparse(substitute(whichAlphaY))
    stopifnot("chrpos" %in% names(Y), colY %in% names(Y))
    Y <- Y[, .(name = chrpos, alpha_Y = get(colY))]
  } else {
    stopifnot(all(c("name", "alpha") %in% names(Y)))
    Y <- Y[, .(name, alpha_Y = alpha)]
  }
  
  # Ensure same type for join column
  X[, name := as.character(name)]
  Y[, name := as.character(name)]
  
  # Explicit inner join on 'name'
  Z_inner <- X[Y, on = "name", nomatch = 0]
  return(Z_inner)
}

makeCompPlot <- function(X, Y, title, xlab, ylab,
                         whichAlphaX = NULL, whichAlphaY = NULL,
                         minplot = 100000) {
  
  # ---- Build Z_inner *and assign it* ----
  Z_inner <- makeZ_inner(X, Y, whichAlphaX = whichAlphaX, whichAlphaY = whichAlphaY)
  
  # ---- Plot & save ----
  c <- cor.test(Z_inner$alpha_X, Z_inner$alpha_Y)
  
  set.seed(1234)
  if(nrow(Z_inner) > minplot){
    Z_inner_plot = Z_inner[sample(nrow(Z_inner), minplot),]  
  } else {
    Z_inner_plot = Z_inner  
  }
  p1 <- ggplot(Z_inner_plot, aes(alpha_X, alpha_Y)) +
    geom_point(pch = 21, alpha = 0.05) +
    geom_abline(slope = 1, linetype = 3) +
    geom_smooth(linetype = 3) +
    geom_smooth(method = "lm", fill = "black") +
    theme_minimal(base_size = 14) +
    annotate("text", x = .2, y = .8, label = paste0("R : ", round(c$estimate, 2))) +
    labs(title = title, x = xlab, y = ylab)
  
  # Make sure the folder exists
  dir.create(here::here("B_MultiTissues/dataOut/figures/correlations"),
             recursive = TRUE, showWarnings = FALSE)
  
  ggplot2::ggsave(
    filename = here::here(paste0("B_MultiTissues/dataOut/figures/correlations/correlation_", title, ".pdf")),
    plot = p1, width = 8, height = 8
  )
  
  invisible(Z_inner)
}

makeGRfromMyCpGPos <- function(vec, setname){# Parse with regex all the cpg tested
  parsed = str_match(vec, "(chr[0-9XYM]+)_(\\d+)")
  
  # Build GRanges
  GR <- GRanges(
    seqnames = parsed[,2],
    ranges   = IRanges(start = as.numeric(parsed[,3]),
                       end   = as.numeric(parsed[,3]))
  )
  GR$set = setname
  return(GR)
}

# ## For a vector of CpGs in the format chromosome_position 
# ## CpGvec <- c("chr1_17452", "chr1_17478", "chr1_17483")
# ## 1. Keep CpGs in regions where at least 5 CpGs are in 50bp distance to each other
# ## 2. annotate with associated genes
# ## 3. run GO term enrichment with clusterProfiler::enrichGO
# 
# ## STEP 1 — ULTRA‑FAST DENSE CpG CLUSTERING (5 CpGs min within 50 bp)
# clusterCpGs <- function(CpGvec, max_gap = 50, min_size = 5) {
#   dt <- data.table(
#     raw = CpGvec,
#     chr = sub("_.*", "", CpGvec),
#     pos = as.integer(sub(".*_", "", CpGvec))
#   )
#   setkey(dt, chr, pos)
#   
#   # gap to previous CpG
#   dt[, gap := pos - data.table::shift(pos), by = chr]
#   
#   # run ID increments whenever gap > max_gap OR different chromosome
#   dt[, run_id := cumsum(is.na(gap) | gap > max_gap), by = chr]
#   
#   # Count CpGs in each run
#   dt[, run_size := .N, by = .(chr, run_id)]
#   
#   # Keep only large runs
#   dt[run_size >= min_size, raw]
# }
# 
# ## STEP 2 — FAST GENE ANNOTATION (OFFLINE)
# annotateCpGs_txdb <- function(CpGs, tss_window = 10000) {
#   
#   if (length(CpGs) == 0) return(character(0))
#   
#   chr <- sub("_.*", "", CpGs)
#   pos <- as.integer(sub(".*_", "", CpGs))
#   gr  <- GRanges(chr, IRanges(pos, pos))
#   
#   # Trim to seqinfo bounds to avoid out-of-bound warnings
#   gr <- GenomicRanges::trim(gr)
#   
#   txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#   genes_txdb <- GenomicFeatures::genes(txdb)
#   promoters_txdb <- GenomicFeatures::promoters(txdb, upstream = tss_window, downstream = tss_window)
#   
#   # overlaps with gene bodies
#   o1 <- findOverlaps(gr, genes_txdb)
#   g1 <- genes_txdb$gene_id[subjectHits(o1)]
#   
#   # overlaps with promoters
#   o2 <- findOverlaps(gr, promoters_txdb)
#   g2 <- promoters_txdb$gene_id[subjectHits(o2)]
#   
#   return(unique(c(g1, g2)))
# }
# 
# ## STEP 3 — GO ENRICHMENT (clusterProfiler)
# runGO <- function(entrez_ids, universe = NULL, myont) {
#   clusterProfiler::enrichGO(
#     gene          = entrez_ids,
#     OrgDb         = org.Hs.eg.db,
#     ont           = myont,
#     keyType       = "ENTREZID",
#     universe      = universe,
#     pAdjustMethod = "BH",
#     readable      = TRUE
#   )
# }
# 
# ## 🚀 FULL PIPELINE FUNCTION 
# CpG_GO_pipeline <- function(CpGvec,
#                             max_gap = 50, min_size = 5,
#                             tss_window = 10000,
#                             universe = NULL) {
#   
#   message("Clustering CpGs...")
#   CpGclustered <- clusterCpGs(CpGvec, max_gap, min_size)
#   message(sprintf("Reduced from %d to %d clustered CpGs",
#                   length(CpGvec), length(CpGclustered)))
#   
#   if (length(CpGclustered) == 0) {
#     warning("No CpG clusters found.")
#     return(NULL)
#   }
#   
#   message("Annotating genes...")
#   ensg <- annotateCpGs_txdb(CpGclustered, tss_window)
#   
#   message(sprintf("Found %d Entrez genes", length(ensg)))
#   
#   message("Running GO enrichment...")
#   list = lapply(c("BP", "MF", "CC"), function(x) {runGO(ensg, universe, x)})
#   names(list) = c("BP", "MF", "CC")
#   return(list)
# }
# 
# ###############
# ## Check GO slim terms for easy interpretation
# ## GO subsets (also known as GO slims) are condensed versions of the GO containing a subset of the terms.
# # dl the GO slim Developed by GO Consortium for the Alliance of Genomes Resources
# # download.file(url = "https://current.geneontology.org/ontology/subsets/goslim_agr.obo",
# # destfile = here("gitignore/goslim_agr.obo"))
# slim <- GSEABase::getOBOCollection(here("gitignore/goslim_agr.obo"))
# 
# getGOslim <- function(x){
#   res = x@result
#   onto = x@ontology
#   
#   # Create the GOCollection 
#   go_collection = GSEABase::GOCollection(ids = res[res$p.adjust < 0.05 & res$Count >= 10, "ID"])
#   
#   # Perform the GO slim mapping
#   slimdf = GSEABase::goSlim(idSrc = go_collection, 
#                             slimCollection = slim,
#                             ontology = onto)
#   
#   # Map the original GO terms to the slim terms
#   mappedIds <- function(df, collection, OFFSPRING) {
#     map <- as.list(OFFSPRING[rownames(df)])
#     mapped <- lapply(map, intersect, ids(collection))
#     df[["mapped_go_terms"]] <- vapply(unname(mapped), paste, collapse = ";", character(1L))
#     
#     # Get GO term names
#     go_names <- AnnotationDbi::select(GO.db, keys = unlist(mapped), columns = "TERM", keytype = "GOID")
#     go_names_list <- split(go_names$TERM, go_names$GOID)
#     df[["mapped_go_names"]] <- vapply(mapped, function(x) {
#       paste(go_names_list[x], collapse = ";")
#     }, character(1L))
#     
#     # Get GO slim term full names
#     slim_names <- AnnotationDbi::select(GO.db, keys = rownames(df), columns = "TERM", keytype = "GOID")
#     df[["go_slim_full_name"]] <- slim_names$TERM
#     
#     df
#   }
#   
#   # Use the appropriate OFFSPRING database based on the ontology
#   offspring_db <- switch(onto,
#                          "BP" = GO.db::GOBPOFFSPRING,
#                          "CC" = GO.db::GOCCOFFSPRING,
#                          "MF" = GO.db::GOMFOFFSPRING)
#   
#   slimdf_with_terms <- mappedIds(slimdf, go_collection, offspring_db)
#   
#   slimdf_with_terms %>%
#     filter(Count != 0) %>%
#     dplyr::mutate(GO.category = onto)
# }
# 
# makeGOslimPlot_all <- function(dfGOslim_all, posleg = "top") {
#   if (nrow(dfGOslim_all) == 0) {
#     stop("dfGOslim_all is empty — no slim categories to plot.")
#   }
#   
#   ggplot(dfGOslim_all, aes(x = group, y = Term)) +
#     geom_point(aes(size = Percent)) +
#     scale_size_continuous(
#       name = "% of GO terms in this GO slim category",
#       range = c(1, 8)
#     ) +
#     theme_bw() +
#     ylab("") + xlab("") +
#     theme(
#       legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
#       legend.background     = element_rect(fill = "#ebebeb", color = "#ebebeb"),
#       legend.key            = element_rect(fill = "#ebebeb", color = "#ebebeb"),
#       legend.position       = posleg,
#       axis.text.y           = element_text(size = 12),
#       axis.text.x           = element_text(size = 12, angle = 45, hjust = 1)
#     ) +
#     facet_grid(group ~ fct_inorder(GO.category), scales = "free", space = "free") +
#     coord_flip()
# }

# test enrichment of ME for each quadrant vs the other three combined

# --- Helper: safe Fisher with edge cases (all zeros, etc.)
.safe_fisher <- function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                dimnames = list(c("this_quadrant","others"), c("inME","notME")))
  # If any row or column totals are zero, Fisher’s test is not defined
  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
    return(list(or = NA_real_, p = NA_real_))
  } else {
    ft <- fisher.test(mat)
    return(list(or = unname(ft$estimate), p = ft$p.value))
  }
}

test_enrichment_quadrants <- function(quad_list, putativeME_GR, me_col = "set") {
  # If no 'set' column, treat all ME as one group "ALL"
  if (!(me_col %in% names(mcols(putativeME_GR)))) {
    me_sets <- "ALL"
    putativeME_GR$..tmp_set.. <- "ALL"
    me_col <- "..tmp_set.."
  } else {
    me_sets <- unique(as.character(mcols(putativeME_GR)[[me_col]]))
  }
  
  # Total elements per quadrant (each element counted once)
  totals <- vapply(quad_list, length, integer(1))
  
  # Iterate over ME sets
  out <- lapply(me_sets, function(me_name) {
    me_subset <- putativeME_GR[mcols(putativeME_GR)[[me_col]] == me_name]
    
    # Count how many elements in each quadrant overlap the ME subset
    inME <- vapply(quad_list, function(gr) {
      ov <- findOverlaps(gr, me_subset, ignore.strand = TRUE)
      length(unique(queryHits(ov)))  # number of quadrant elements hitting at least one ME
    }, integer(1))
    
    # Build quadrant-vs-others 2x2 tests
    bind_rows(lapply(names(quad_list), function(q) {
      a <- inME[[q]]
      b <- totals[[q]] - a
      c <- sum(inME[names(inME) != q])
      d <- sum(totals[names(totals) != q]) - c
      
      fs <- .safe_fisher(a, b, c, d)
      
      tibble(
        CpG_set            = me_name,
        quadrant          = q,
        inME              = a,
        total             = totals[[q]],
        others_inME       = c,
        others_total      = sum(totals) - totals[[q]],
        pct_inME          = 100 * a / totals[[q]],
        pct_inME_others   = 100 * c / (sum(totals) - totals[[q]]),
        odds_ratio        = fs$or,
        p_value           = fs$p
      )
    }))
  }) %>% bind_rows() %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
  
  out
}

# Compute overlap across any number of groups and plot a Venn diagram
plotMyVenn <- function(cutoff, ...) {
  groups <- list(...)                         # list of data.frames (each has name, alpha)
  
  # 1) Overlap among all names (unfiltered)
  sets_unfilt <- lapply(groups, function(df) df$name)
  overlap <- Reduce(intersect, sets_unfilt)
  
  message(paste0("There are ", length(overlap), " overlapping CpGs between these groups (unfiltered)."))
  
  # 2) Apply cutoff AND keep only overlapping names (so sets are comparable on the same universe)
  sets <- lapply(groups, function(df) df$name[df$alpha >= cutoff & df$name %in% overlap])
  
  # 3) Optional: add names to sets if you passed named arguments
  if (!is.null(dots <- match.call(expand.dots = FALSE)$...) && length(names(dots))) {
    nm <- names(dots)
    if (any(nzchar(nm))) names(sets) <- ifelse(nzchar(nm), nm, paste0("group", seq_along(sets)))
  }
  
  # 4) Plot
  p <- ggVennDiagram(sets, label_alpha = 0) +
    scale_fill_gradient2(low = "white", mid = "yellow", high = "red") +
    ggtitle(paste0("Pr(hv) ≥ ", cutoff), 
            subtitle = paste0("Sequenced CpGs N = ", length(overlap)))
  
  return(p)
}

functionsLoaded = TRUE
