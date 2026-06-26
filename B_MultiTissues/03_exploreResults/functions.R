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
## Functions extracted from S06 for reuse in downstream scripts:
### make_MEsetdt
### make_MEsetdt_regionMean
### plot_decay_curve
### plot_decay_curve_layered (by germ layers)
### compute_percpg_interlayer_corr: Compute per-CpG interlayer correlation from raw meth

makeVennArrayReduced <- function(df_circles, v, counts, fmt_fn){
  size = 4
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = df_circles,
                          aes(x, y, group = group),
                          fill = "#4981BF", alpha = 0.3,
                          colour = "grey40", linewidth = 0.5) +
    # Exclusive regions
    ggplot2::annotate("text", x = -0.9,  y =  0.55, label = fmt_fn(counts$n_only_full,      v$n_union), size = size, lineheight = 0.9) +
    ggplot2::annotate("text", x =  0.9,  y =  0.55, label = fmt_fn(counts$n_only_2ind,      v$n_union), size = size, lineheight = 0.9) +
    ggplot2::annotate("text", x =  0.0,  y = -1.05, label = fmt_fn(counts$n_only_3ind,      v$n_union), size = size, lineheight = 0.9) +
    # Pairwise only
    ggplot2::annotate("text", x =  0.0,  y =  0.75, label = fmt_fn(counts$n_full_2ind_only, v$n_union), size = size, lineheight = 0.9) +
    ggplot2::annotate("text", x = -0.55, y = -0.25, label = fmt_fn(counts$n_full_3ind_only, v$n_union), size = size, lineheight = 0.9) +
    ggplot2::annotate("text", x =  0.55, y = -0.25, label = fmt_fn(counts$n_2ind_3ind_only, v$n_union), size = size, lineheight = 0.9) +
    # Triple
    ggplot2::annotate("text", x =  0.0,  y =  0.15, label = fmt_fn(v$n_all,          v$n_union), size = size, lineheight = 0.9) +
    # Set labels with total size
    ggplot2::annotate("text", x = -1,  y =  1.5,  label = "Full array", size = size, fontface = "bold") +
    ggplot2::annotate("text", x =  .9,  y =  1.5,  label = "Array 2 ind/ds", size = size, fontface = "bold") +
    ggplot2::annotate("text", x =  0.0,  y = -1.7,  label = "Array 3 ind/ds", size = size, fontface = "bold") +
    ggplot2::coord_fixed() +
    ggplot2::theme_void()
}

prepAtlasdt <- function(subdir, p0, p1) {
  parent_dir <- here(paste0("B_MultiTissues/resultsDir_gitIgnored/Atlas/", subdir))
  rds_files  <- base::dir(parent_dir, pattern = paste0(p0, "p0_", p1, "p1.rds$"),
                          recursive = TRUE, full.names = TRUE)
  
  if (length(rds_files) == 0)
    stop("No files matched in: ", parent_dir)
  
  all_cpg_values <- numeric()
  pb <- progress_bar$new(total = length(rds_files), 
                         format = "📦 :current/:total [:bar] :percent")
  
  for (file in rds_files) {
    obj <- readRDS(file)
    if (is.matrix(obj)) obj <- setNames(obj[, 1], rownames(obj))  # fix order
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
  
  # Convert chr from "chrN" to  factor
  dt[, chr := sub("chr", "", chr)]
  dt[, chr := factor(chr, levels = as.character(c(1:22, "X", "Y", "M")))]
  
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

plotManhattanFromdt <- function(dt, transp = 0.01, plotDerakhshan = TRUE,
                                colorBySet = FALSE){
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
  
  p <- ggplot() +
    # Add dotted separators
    geom_vline(data = vlines, aes(xintercept = xintercept),
               linetype = 3, color = "black", linewidth = 1) +
    theme_classic() + theme(legend.position = "none") +
    scale_x_continuous(breaks = df2$center, labels = as.character(df2$chr), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Chromosome", y = "Pr(hv)")+
    theme_minimal(base_size = 14)
  
  if (plotDerakhshan == TRUE){
    p <- p +  
      # background cloud
      geom_point_rast(data = dt[is.na(group)], 
                      aes(x = pos2, y = alpha),
                      color = "black", size = 0.01, alpha = transp, raster.dpi = 72) +
      geom_point(data = dt[group == "hvCpG_Derakhshan"],
                 aes(x = pos2, y = alpha),
                 color = "#DC3220", size = 1, alpha = 0.7) +
      # mQTL controls highlights
      geom_point(data = dt[group == "mQTLcontrols"],
                 aes(x = pos2, y = alpha),
                 color = "#005AB5", size = 1, alpha = 0.7)
  }
  if (colorBySet == TRUE){
    p <- p +
      geom_point(data = dt,
                 aes(x = pos2, y = alpha, color = set),
                 alpha = transp, size = 1) +
      facet_wrap(.~set, nrow = 5)
  } else {
    p <- p +
      # background cloud
      geom_point_rast(data = dt, 
                      aes(x = pos2, y = alpha),
                      color = "black", size = 0.01, alpha = transp, raster.dpi = 72) 
  }
  return(p)
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
                         minplot = 100000, drawline = TRUE) {
  
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
    theme_minimal(base_size = 14) +
    annotate("text", x = .2, y = .8, label = paste0("R : ", round(c$estimate, 2))) +
    labs(title = title, x = xlab, y = ylab)
  
  if (drawline == TRUE){
    p1 <- p1 +     
      geom_smooth(linetype = 3) +
      geom_smooth(method = "lm", fill = "black")
  }
  
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

########
## GO ##
########

## For a vector of CpGs in the format chromosome_position
## CpGvec <- c("chr1_17452", "chr1_17478", "chr1_17483")
## 1. Keep CpGs in regions where at least 5 CpGs are in 50bp distance to each other
## 2. annotate with associated genes CONTROLLING FOR GENE LENGTH
## 3. run GO term enrichment with clusterProfiler::enrichGO

## ULTRA‑FAST DENSE CpG CLUSTERING (5 CpGs min within 50 bp)
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

## FAST GENE ANNOTATION (OFFLINE)
annotateCpGs_txdb <- function(CpGs, tss_window = 10000) {

  if (length(CpGs) == 0) return(character(0))

  chr <- sub("_.*", "", CpGs)
  pos <- as.integer(sub(".*_", "", CpGs))
  gr  <- GRanges(chr, IRanges(pos, pos))

  # Trim to seqinfo bounds to avoid out-of-bound warnings
  gr <- GenomicRanges::trim(gr)

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genes_txdb <- GenomicFeatures::genes(txdb)
  promoters_txdb <- GenomicFeatures::promoters(txdb, upstream = tss_window, downstream = tss_window)

  # overlaps with gene bodies
  o1 <- findOverlaps(gr, genes_txdb)
  g1 <- genes_txdb$gene_id[subjectHits(o1)]

  # overlaps with promoters
  o2 <- findOverlaps(gr, promoters_txdb)
  g2 <- promoters_txdb$gene_id[subjectHits(o2)]

  return(unique(c(g1, g2)))
}

# ── Helper: get gene lengths from TxDb ──────────────────────────────────────
get_gene_lengths <- function() {
  txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
  exons <- GenomicFeatures::exonsBy(txdb, by = "gene")
  # Sum of non-overlapping exon widths per gene
  lengths <- sum(width(GenomicRanges::reduce(exons)))
  data.table(entrez_id = names(lengths), gene_length = as.integer(lengths))
}

# ── Length-matched universe subsampling ─────────────────────────────────────
length_match_universe <- function(foreground_genes, universe_genes,
                                  gene_lengths_dt, n_bins = 20, seed = 42) {
  set.seed(seed)
  
  fg_dt  <- data.table(entrez_id = foreground_genes)
  uni_dt <- data.table(entrez_id = universe_genes)
  
  # Merge lengths
  fg_dt  <- merge(fg_dt,  gene_lengths_dt, by = "entrez_id", all.x = TRUE)
  uni_dt <- merge(uni_dt, gene_lengths_dt, by = "entrez_id", all.x = TRUE)
  
  # Bin gene lengths
  breaks <- quantile(fg_dt$gene_length, probs = seq(0, 1, length.out = n_bins + 1),
                     na.rm = TRUE)
  breaks[1]         <- 0
  breaks[n_bins + 1] <- Inf
  
  fg_dt[,  len_bin := cut(gene_length, breaks, include.lowest = TRUE)]
  uni_dt[, len_bin := cut(gene_length, breaks, include.lowest = TRUE)]
  
  # For each bin, sample from universe to match foreground count
  fg_counts  <- fg_dt[, .N, by = len_bin]
  matched_bg <- uni_dt[!entrez_id %in% foreground_genes][  # exclude foreground
    , .SD[sample(.N, min(.N, fg_counts[len_bin == .BY$len_bin, N] * 10))],
    # 10x more background than foreground per bin = good power
    by = len_bin
  ]
  
  unique(c(foreground_genes, matched_bg$entrez_id))
}

## GO ENRICHMENT core function (clusterProfiler)
runGO <- function(entrez_ids, universe = NULL, myont) {
  clusterProfiler::enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = myont,
    keyType       = "ENTREZID",
    universe      = universe,
    pAdjustMethod = "BH",
    readable      = TRUE
  )
}

## 🚀 FULL PIPELINE FUNCTIONS
# ── Updated CpG_GO_pipeline with length control ──────────────────────────────
CpG_GO_pipeline_lengthControlled <- function(CpGvec,
                                             max_gap    = 50,
                                             min_size   = minimum_CpG_per_cluster,
                                             tss_window = 10000,
                                             universe   = NULL,
                                             control_length = TRUE) {
  message("Clustering CpGs...")
  CpGclustered <- clusterCpGs(CpGvec, max_gap, min_size)
  message(sprintf("Reduced from %d to %d clustered CpGs",
                  length(CpGvec), length(CpGclustered)))
  
  if (length(CpGclustered) == 0) { warning("No CpG clusters found."); return(NULL) }
  
  message("Annotating genes...")
  ensg <- annotateCpGs_txdb(CpGclustered, tss_window)
  message(sprintf("Found %d Entrez genes", length(ensg)))
  
  if (control_length) {
    message("Controlling for gene length...")
    gene_lengths_dt <- get_gene_lengths()
    
    # Check: are foreground genes longer than universe?
    fg_len  <- gene_lengths_dt[entrez_id %in% ensg, median(gene_length, na.rm = TRUE)]
    uni_len <- gene_lengths_dt[entrez_id %in% universe, median(gene_length, na.rm = TRUE)]
    message(sprintf("  Median gene length — foreground: %s bp, universe: %s bp, ratio: %.2f",
                    format(fg_len, big.mark = ","),
                    format(uni_len, big.mark = ","),
                    fg_len / uni_len))
    
    universe_matched <- length_match_universe(ensg, universe, gene_lengths_dt)
    message(sprintf("  Length-matched universe: %d genes (was %d)",
                    length(universe_matched), length(universe)))
  } else {
    universe_matched <- universe
  }
  
  message("Running GO enrichment...")
  result <- lapply(c("BP", "MF", "CC"), function(x) {
    runGO(ensg, universe_matched, x)
  })
  names(result) <- c("BP", "MF", "CC")
  return(result)
}

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

## functions_S06.R
## Functions extracted from S06 for reuse in downstream scripts

make_MEsetdt <- function(sets, geomMeanGR) {
  MEsetdt <- rbindlist(lapply(names(sets), function(nm) {
    hits <- findOverlaps(sets[[nm]], geomMeanGR)
    data.table(
      alpha_geomean = geomMeanGR$alpha_geomean[subjectHits(hits)],
      ME = nm
    )
  }))
  na.omit(MEsetdt)
}

make_MEsetdt_regionMean <- function(sets, geomMeanGR) {
  MEsetdt <- rbindlist(lapply(names(sets), function(nm) {
    hits <- findOverlaps(sets[[nm]], geomMeanGR)
    dt <- data.table(
      region_idx    = queryHits(hits),
      alpha_geomean = geomMeanGR$alpha_geomean[subjectHits(hits)],
      ME = nm
    )
    dt[, .(alpha_geomean = mean(alpha_geomean, na.rm = TRUE)),
       by = .(region_idx, ME)]
  }))
  na.omit(MEsetdt)
}

plot_decay_curve <- function(MEsetdt, title = "Decay curve") {
  thresholds <- seq(10, 90, by = 10) / 100
  prop_table <- rbindlist(lapply(thresholds, function(thr) {
    MEsetdt[, .(
      proportion = mean(alpha_geomean > thr, na.rm = TRUE),
      n_above    = sum(alpha_geomean > thr, na.rm = TRUE),
      n_total    = .N
    ), by = ME][, threshold := thr]
  }))
  
  # Build colour vector: black for reference, Set2 for the rest
  me_levels <- unique(prop_table$ME)
  other_levels <- setdiff(me_levels, "mQTLcontrols")
  set2_cols <- RColorBrewer::brewer.pal(max(length(other_levels), 3), "Set2")
  my_cols <- c(mQTLcontrols = "black",
               setNames(set2_cols[seq_along(other_levels)], other_levels))
  
  ggplot(prop_table, aes(x = threshold, y = proportion, colour = ME)) +
    geom_line() +
    geom_point() +
    scale_x_continuous("Pr(HV) threshold", breaks = thresholds) +
    scale_y_continuous("Proportion above threshold", labels = scales::percent) +
    scale_colour_manual(values = my_cols) +
    theme_bw() +
    ggtitle(title)
}

plot_decay_curve_layered <- function(MEsetdt, title = "Decay curve by layer") {
  thresholds <- seq(0, 100, by = 10) / 100
  proportion <- seq(0, 100, by = 10) / 100
  
  # melt the 3 alpha columns to long format
  long <- melt(MEsetdt,
               id.vars = "ME",
               measure.vars = c("alpha_endo", "alpha_meso", "alpha_ecto"),
               variable.name = "layer", value.name = "alpha")
  
  prop_table <- rbindlist(lapply(thresholds, function(thr) {
    long[, .(
      proportion = mean(alpha > thr, na.rm = TRUE),
      n_above    = sum(alpha > thr, na.rm = TRUE),
      n_total    = .N
    ), by = .(ME, layer)][, threshold := thr]
  }))
  
  # colours
  me_levels    <- unique(prop_table$ME)
  other_levels <- setdiff(me_levels, "mQTLcontrols")
  set2_cols    <- RColorBrewer::brewer.pal(max(length(other_levels), 3), "Set2")
  my_cols      <- c(mQTLcontrols = "black",
                    setNames(set2_cols[seq_along(other_levels)], other_levels))
  
  ggplot(prop_table, aes(x = threshold, y = proportion, colour = ME)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ layer, nrow = 1,
               labeller = labeller(layer = c(
                 alpha_endo = "Endoderm",
                 alpha_meso = "Mesoderm",
                 alpha_ecto = "Ectoderm"))) +
    scale_x_continuous("Pr(HV) threshold", breaks = thresholds) +
    scale_y_continuous("Proportion above threshold", breaks = proportion, labels = scales::percent) +
    scale_colour_manual(values = my_cols) +
    theme_bw() +
    ggtitle(title)
}

plot_residuals_layered <- function(MEsetdt, title = "Excess Pr(HV) vs mQTLcontrols by layer") {
  thresholds <- seq(0, 100, by = 10) / 100
  
  long <- melt(MEsetdt,
               id.vars = "ME",
               measure.vars = c("alpha_endo", "alpha_meso", "alpha_ecto"),
               variable.name = "layer", value.name = "alpha")
  
  prop_table <- rbindlist(lapply(thresholds, function(thr) {
    long[, .(
      proportion = mean(alpha > thr, na.rm = TRUE)
    ), by = .(ME, layer)][, threshold := thr]
  }))
  
  # subtract control proportion per layer × threshold
  ctrl <- prop_table[ME == "mQTLcontrols", .(layer, threshold, ctrl_proportion = proportion)]
  prop_table <- merge(prop_table, ctrl, by = c("layer", "threshold"))
  prop_table[, residual := proportion - ctrl_proportion]
  prop_table <- prop_table[ME != "mQTLcontrols"]
  
  # layer factor order
  layer_labels <- c(alpha_endo = "Endoderm",
                    alpha_meso = "Mesoderm",
                    alpha_ecto = "Ectoderm")
  prop_table[, layer_label := factor(layer_labels[as.character(layer)],
                                     levels = c("Endoderm", "Mesoderm", "Ectoderm"))]
  
  # colours
  me_levels    <- unique(prop_table$ME)
  other_levels <- setdiff(me_levels, "mQTLcontrols")
  set2_cols    <- RColorBrewer::brewer.pal(max(length(other_levels), 3), "Set2")
  my_cols      <- setNames(set2_cols[seq_along(other_levels)], other_levels)
  
  ggplot(prop_table, aes(x = threshold, y = residual, colour = ME)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.4) +
    facet_wrap(~ layer_label, nrow = 1) +
    scale_x_continuous("Pr(HV) threshold", breaks = thresholds) +
    scale_y_continuous("Excess proportion vs mQTLcontrols",
                       labels = scales::percent) +
    scale_colour_manual(values = my_cols) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right") +
    ggtitle(title)
}

# ── Compute per-CpG interlayer correlation from raw meth
compute_percpg_interlayer_corr <- function(meth_sub) {
  
  meth_sub[, pos := as.integer(sub(".*_", "", cpg_site))]
  
  # mean per (patient, germ_layer, pos) — collapses multiple tissues same layer
  meth_agg <- meth_sub[, .(methylation = mean(methylation, na.rm = TRUE)),
                       by = .(patient_id, germ_layer, pos, cpg_site)]
  
  # wide: one column per germ layer
  wide <- dcast(meth_agg, patient_id + pos + cpg_site ~ germ_layer,
                value.var = "methylation")
  
  layers_present <- intersect(c("Endo", "Meso", "Ecto"), names(wide))
  if (length(layers_present) < 2) return(NULL)
  
  layer_pairs <- combn(layers_present, 2, simplify = FALSE)
  
  rbindlist(lapply(layer_pairs, function(pair) {
    l1 <- pair[1]; l2 <- pair[2]
    
    # per-CpG correlation across patients
    wide[, {
      idx <- !is.na(get(l1)) & !is.na(get(l2))
      if (sum(idx) >= 3) {
        list(r    = cor(get(l1)[idx], get(l2)[idx], method = "pearson"),
             n    = sum(idx),
             pair = paste(l1, l2, sep = "-"))
      } else {
        list(r = NA_real_, n = sum(idx), pair = paste(l1, l2, sep = "-"))
      }
    }, by = .(pos, cpg_site)]
  }))
}

functionsLoaded = TRUE
