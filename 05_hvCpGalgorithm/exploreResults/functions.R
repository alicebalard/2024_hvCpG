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
                         minplot = 1000000) {
  
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
  dir.create(here::here("05_hvCpGalgorithm/figures/correlations"),
             recursive = TRUE, showWarnings = FALSE)
  
  ggplot2::ggsave(
    filename = here::here(paste0("05_hvCpGalgorithm/figures/correlations/correlation_", title, ".pdf")),
    plot = p1, width = 8, height = 8
  )
  
  invisible(Z_inner)
}

functionsLoaded = TRUE
