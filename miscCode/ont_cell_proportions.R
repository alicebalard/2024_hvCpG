#!/usr/bin/env Rscript
###############################################################################
# ont_cell_proportions.R
#
# Estimates white blood cell proportions from DNA methylation called on
# Oxford Nanopore (ONT) data, and flags samples whose coverage is too low for
# the estimate to be trusted.
#
# METHOD
#   Constrained projection of Houseman et al. (2012, BMC Bioinformatics 13:86).
#   For each sample, the observed methylation fractions at a set of reference
#   CpGs are regressed on the cell-type-specific methylation of those same
#   CpGs (the reference matrix), under two constraints: no proportion may be
#   negative, and the proportions may not sum to more than 1. Only the CpGs
#   that have adequate coverage in a given sample are used for that sample.
#
#   The default reference is the Jaffe et al. (2014) set of 500 CpGs derived
#   from sorted adult venous blood (Reinius et al. 2012, 450k array), with six
#   cell types: CD8T, CD4T, NK, Bcell, Mono, Gran. It is prepared once, in
#   Nanopore-compatible coordinates, by build_reference_jaffe.R. Any other
#   reference with the same layout can be supplied instead (see REFERENCE).
#
# INPUT (choose one with --input_type)
#   bedmethyl (default)
#       --input is a tab-separated sample sheet with a header and two columns,
#       'sample' and 'file' (paths are relative to the directory the script is
#       run from). Each file is a modkit bedMethyl table (optionally gzipped),
#       as written by 'modkit pileup'. The columns used are:
#       1 chromosome, 2 start (0-based), 4 modification code, 6 strand,
#       10 valid coverage, 12 number of modified calls.
#   long
#       --input is a single tab-separated table with a header and the columns
#         sample   sample identifier
#         chrom    chromosome (with or without the 'chr' prefix)
#         pos      1-based position of the CpG cytosine (or of its guanine on
#                  the opposite strand)
#         n_valid  number of reads with a usable call at that position
#         n_mod    number of reads called as methylated   (preferred), or
#         meth     methylation level, in the units set by --scale
#       An optional 'strand' column is accepted. Supply one row per position
#       and strand. Use this format for output from callers other than modkit.
#
# REFERENCE (--reference)
#   Tab-separated, with a header. The first three columns are
#       cpg        CpG identifier (for example the array probe name)
#       chrom      chromosome, in the same genome build as the ONT alignments
#       cpg_start  1-based position of the C of the CpG on the forward strand
#   Every remaining column is a cell type and holds the methylation level
#   (0 to 1) of that CpG in that cell type. Reads from both strands of a CpG
#   (positions cpg_start and cpg_start + 1) are pooled before use.
#
# OUTPUT (files are named from --output)
#   <output>_cell_proportions.tsv  One row per sample: cell proportions, an
#                                  overall flag, and the coverage measures
#                                  described below.
#   <output>_coverage_qc.pdf       CpG coverage per sample and sampling error.
#   <output>_run_settings.txt      Settings and software versions for the
#                                  methods section.
#
# HOW LOW COVERAGE IS FLAGGED
#   A CpG counts as 'used' in a sample if it has at least --min_depth valid
#   reads. A sample is flagged (flag = FLAG, with the reasons listed in
#   flag_reason) if any of the following holds.
#     too_few_cpgs          Fewer than --min_cpgs CpGs are used (and always if
#                           fewer than one more than the number of cell types).
#                           No estimate is produced for these samples.
#     collinear_reference   The CpGs used cannot separate the cell types (the
#                           reference matrix for that subset is rank deficient).
#                           No estimate is produced.
#     low_cpg_proportion    The proportion of reference CpGs used is below
#                           --min_prop.
#     high_sampling_error   The largest standard error among the cell types
#                           exceeds --max_se. The standard error comes from
#                           re-drawing the methylated read counts at every used
#                           CpG from a binomial distribution (--nboot times)
#                           and re-fitting. It reflects read-depth noise and
#                           how well the available CpGs separate the cell
#                           types. It does not capture differences between the
#                           reference and the sample population, or between
#                           array and Nanopore chemistry.
#   Samples with reasons other than the first two still receive an estimate,
#   unless --mask_flagged is given. The default thresholds are pragmatic
#   starting values and have not been calibrated on Nanopore data. Check them
#   by thinning reads or CpGs from your best-covered samples and seeing where
#   the estimates start to drift.
#
# OTHER COLUMNS IN THE OUTPUT
#   sum_proportions    Sum of the estimated proportions. The constraint allows
#                      this to be below 1. Values well below 1 suggest the
#                      reference describes the sample poorly.
#   fit_rmse           Root mean square difference between observed and fitted
#                      methylation at the used CpGs.
#   n_cpgs_covered     Reference CpGs with at least one valid read.
#   n_cpgs_used        Reference CpGs with at least --min_depth valid reads.
#   median_depth_used  Median valid read depth over the used CpGs.
#   condition_number   Condition number of the reference matrix restricted to
#                      the used CpGs. Large values mean the cell types are hard
#                      to tell apart with that subset of CpGs.
#   se_<cell type>     Standard error described above.
#
# USAGE
#   Rscript ont_cell_proportions.R --input samples.tsv \
#       --reference jaffe_hg38_reference.tsv --output results/cohort1
#
#   Options (defaults in brackets)
#     --input_type    bedmethyl or long                               [bedmethyl]
#     --scale         percent or fraction, for 'meth' in long input     [percent]
#     --mod_codes     modification codes to count; use m,h to count
#                     5mC plus 5hmC, which is what bisulfite arrays
#                     measure                                                 [m]
#     --min_depth     minimum valid reads for a CpG to be used               [10]
#     --min_cpgs      minimum number of CpGs used per sample                 [50]
#     --min_prop      minimum proportion of reference CpGs used            [0.25]
#     --max_se        largest acceptable standard error                    [0.05]
#     --nboot         resamples for the standard error (0 switches off)     [100]
#     --seed          random seed for the resampling                          [1]
#     --unconstrained fit without the non-negativity and sum constraints
#     --mask_flagged  set the estimates to NA for every flagged sample
#     --help          print this text
#
#   The functions can also be used from an R session:
#     source("ont_cell_proportions.R")
#     opts <- default_options()
#     opts$input <- "samples.tsv"; opts$reference <- "ref.tsv"
#     res <- run_cell_proportions(opts)
#
# REQUIREMENTS
#   R with the packages data.table and quadprog (both on CRAN).
#
# REFERENCES
#   Houseman EA et al. (2012) BMC Bioinformatics 13:86.
#   Jaffe AE, Irizarry RA (2014) Genome Biology 15:R31.
#   Reinius LE et al. (2012) PLoS ONE 7:e41361.
###############################################################################

suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("quadprog", quietly = TRUE)) {
  stop("The 'quadprog' package is required: install.packages('quadprog')")
}


# ---------------------------------------------------------------------------
# Settings and command line
# ---------------------------------------------------------------------------

default_options <- function() {
  list(
    input         = NA_character_,
    input_type    = "bedmethyl",
    reference     = NA_character_,
    output        = "cell_proportions",
    scale         = "percent",
    mod_codes     = "m",
    min_depth     = 10,
    min_cpgs      = 50,
    min_prop      = 0.25,
    max_se        = 0.05,
    nboot         = 100,
    seed          = 1,
    unconstrained = FALSE,
    mask_flagged  = FALSE,
    help          = FALSE
  )
}

# Reads --name value pairs from the command line. Options whose default is
# TRUE/FALSE are switches and take no value.
parse_args <- function(defaults) {
  a <- commandArgs(trailingOnly = TRUE)
  opts <- defaults
  i <- 1
  while (i <= length(a)) {
    if (!grepl("^--", a[i])) stop("Unexpected argument: ", a[i])
    key <- sub("^--", "", a[i])
    if (!(key %in% names(defaults))) stop("Unknown option: --", key)
    if (is.logical(defaults[[key]])) {
      opts[[key]] <- TRUE
      i <- i + 1
    } else {
      if (i == length(a)) stop("Option --", key, " needs a value")
      val <- a[i + 1]
      if (is.numeric(defaults[[key]])) {
        val <- suppressWarnings(as.numeric(val))
        if (is.na(val)) stop("Option --", key, " must be a number")
      }
      opts[[key]] <- val
      i <- i + 2
    }
  }
  opts
}

check_options <- function(opts) {
  if (is.na(opts$input) || is.na(opts$reference)) {
    stop("Both --input and --reference are required (see --help)")
  }
  if (!file.exists(opts$input)) stop("Input not found: ", opts$input)
  if (!file.exists(opts$reference)) stop("Reference not found: ", opts$reference)
  if (!(opts$input_type %in% c("bedmethyl", "long"))) {
    stop("--input_type must be 'bedmethyl' or 'long'")
  }
  if (!(opts$scale %in% c("percent", "fraction"))) {
    stop("--scale must be 'percent' or 'fraction'")
  }
  if (opts$min_depth < 1) stop("--min_depth must be at least 1")
  if (opts$min_prop < 0 || opts$min_prop > 1) stop("--min_prop must be between 0 and 1")
  if (opts$nboot < 0) stop("--nboot cannot be negative")
  invisible(opts)
}


# ---------------------------------------------------------------------------
# Reference
# ---------------------------------------------------------------------------

# Adds the 'chr' prefix if it is missing, so that 'chr1' and '1' both work.
norm_chrom <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

read_reference <- function(path) {
  ref <- fread(path, data.table = TRUE)
  needed <- c("cpg", "chrom", "cpg_start")
  if (!all(needed %in% names(ref))) {
    stop("Reference must contain the columns: ", paste(needed, collapse = ", "))
  }
  cell_types <- setdiff(names(ref), needed)
  if (length(cell_types) < 2) {
    stop("Reference needs at least two cell-type columns after cpg, chrom, cpg_start")
  }
  if (anyDuplicated(ref$cpg)) stop("Duplicated CpG identifiers in the reference")

  ref[, chrom := norm_chrom(chrom)]
  ref[, cpg_start := as.integer(cpg_start)]
  ref <- ref[complete.cases(ref[, c("cpg", "chrom", "cpg_start", cell_types), with = FALSE])]

  X <- as.matrix(ref[, cell_types, with = FALSE])
  rownames(X) <- ref$cpg
  # Regression-based references can stray slightly outside 0 to 1, so only
  # values that look like percentages are rejected.
  if (max(X) > 1.5) {
    stop("Reference values must be methylation fractions between 0 and 1, not percentages")
  }
  list(X = X,
       sites = ref[, .(cpg, chrom, cpg_start)],
       cell_types = cell_types)
}


# ---------------------------------------------------------------------------
# Reading methylation calls
# ---------------------------------------------------------------------------

# When one position has several rows (for example separate 5mC and 5hmC rows
# from modkit), the valid coverage is the same on each, so it is taken once
# and the modified counts are added.
collapse_codes <- function(dt) {
  keys <- intersect(c("sample", "chrom", "pos", "strand"), names(dt))
  dt[, .(n_valid = max(n_valid), n_mod = sum(n_mod)), by = keys]
}

read_bedmethyl <- function(file, mod_codes) {
  if (!file.exists(file)) stop("bedMethyl file not found: ", file)
  keep_cols <- c(1, 2, 4, 6, 10, 12)
  if (grepl("\\.gz$", file)) {
    dt <- fread(cmd = paste("gzip -dc", shQuote(file)), header = FALSE,
                select = keep_cols, showProgress = FALSE)
  } else {
    dt <- fread(file, header = FALSE, select = keep_cols, showProgress = FALSE)
  }
  setnames(dt, c("V1", "V2", "V4", "V6", "V10", "V12"),
           c("chrom", "start", "code", "strand", "n_valid", "n_mod"))

  codes <- strsplit(mod_codes, ",")[[1]]
  dt <- dt[code %in% codes]
  # bedMethyl starts are 0-based; the reference positions are 1-based.
  dt[, pos := as.integer(start) + 1L]
  dt[, chrom := norm_chrom(chrom)]
  dt <- dt[n_valid > 0, .(chrom, pos, strand, n_valid, n_mod)]
  collapse_codes(dt)
}

read_long_table <- function(path, scale) {
  dt <- fread(path, data.table = TRUE)
  needed <- c("sample", "chrom", "pos", "n_valid")
  if (!all(needed %in% names(dt))) {
    stop("Long-format input must contain the columns: ", paste(needed, collapse = ", "),
         ", and either n_mod or meth")
  }
  if ("n_mod" %in% names(dt)) {
    # use the counts as supplied
  } else if ("meth" %in% names(dt)) {
    frac <- if (scale == "percent") dt$meth / 100 else dt$meth
    dt[, n_mod := round(frac * n_valid)]
  } else {
    stop("Long-format input needs a column called n_mod or meth")
  }
  if (!("strand" %in% names(dt))) dt[, strand := "."]

  dt[, sample := as.character(sample)]
  dt[, chrom := norm_chrom(chrom)]
  dt[, pos := as.integer(pos)]
  dt <- dt[!is.na(n_valid) & n_valid > 0 & !is.na(n_mod),
           .(sample, chrom, pos, strand, n_valid, n_mod)]
  dt[, n_mod := pmin(n_mod, n_valid)]
  collapse_codes(dt)
}

# Assigns each call to a reference CpG and pools the two strands. A CpG
# occupies positions cpg_start (the C) and cpg_start + 1 (the G on the other
# strand). Two CpGs cannot share a position, so the assignment is unambiguous.
aggregate_to_cpg <- function(dt, sites) {
  lookup <- rbind(
    sites[, .(cpg, chrom, pos = cpg_start)],
    sites[, .(cpg, chrom, pos = cpg_start + 1L)]
  )
  m <- merge(dt, lookup, by = c("chrom", "pos"))
  m[, .(n_valid = sum(n_valid), n_mod = sum(n_mod)), by = .(sample, cpg)]
}

read_all_input <- function(opts, sites) {
  if (opts$input_type == "bedmethyl") {
    sheet <- read.delim(opts$input, stringsAsFactors = FALSE, colClasses = "character")
    if (!all(c("sample", "file") %in% names(sheet))) {
      stop("The sample sheet needs the columns 'sample' and 'file'")
    }
    if (anyDuplicated(sheet$sample)) stop("Duplicated sample names in the sample sheet")
    pieces <- vector("list", nrow(sheet))
    for (i in seq_len(nrow(sheet))) {
      message("Reading ", sheet$sample[i], " (", i, " of ", nrow(sheet), ")")
      raw <- read_bedmethyl(sheet$file[i], opts$mod_codes)
      raw[, sample := sheet$sample[i]]
      pieces[[i]] <- aggregate_to_cpg(raw, sites)
    }
    return(list(data = rbindlist(pieces), samples = as.character(sheet$sample)))
  }
  raw <- read_long_table(opts$input, opts$scale)
  list(data = aggregate_to_cpg(raw, sites), samples = unique(raw$sample))
}


# ---------------------------------------------------------------------------
# Estimation
# ---------------------------------------------------------------------------

# Projects the observed methylation y (one value per CpG) onto the reference
# matrix X (CpGs x cell types). With constrained = TRUE the proportions are
# non-negative and sum to at most 1. Returns NA if the fit fails.
project_one <- function(y, X, constrained = TRUE) {
  K <- ncol(X)
  XtX <- crossprod(X)
  Xty <- as.vector(crossprod(X, y))
  tryCatch({
    if (constrained) {
      Amat <- cbind(rep(-1, K), diag(K))   # -sum(b) >= -1, and b >= 0
      b0   <- c(-1, rep(0, K))
      quadprog::solve.QP(XtX, Xty, Amat, b0)$solution
    } else {
      as.vector(solve(XtX, Xty))
    }
  }, error = function(e) rep(NA_real_, K))
}

# Standard error of each proportion from binomial re-drawing of the read
# counts. A half-read is added to each side when forming the success
# probability so that CpGs observed at exactly 0 or 1 are not treated as
# noise-free.
sampling_se <- function(n_valid, n_mod, Xo, nboot, constrained) {
  p <- (n_mod + 0.5) / (n_valid + 1)
  n <- round(n_valid)
  draws <- matrix(NA_real_, nrow = nboot, ncol = ncol(Xo))
  for (b in seq_len(nboot)) {
    y_b <- rbinom(length(p), size = n, prob = p) / n
    draws[b, ] <- project_one(y_b, Xo, constrained)
  }
  apply(draws, 2, sd, na.rm = TRUE)
}

# Fits one sample. 'used' holds the CpGs that passed the depth filter.
estimate_sample <- function(used, X, opts) {
  K <- ncol(X)
  out <- list(status = "ok",
              est = rep(NA_real_, K), se = rep(NA_real_, K),
              rmse = NA_real_, cond = NA_real_)

  if (nrow(used) < max(opts$min_cpgs, K + 1)) {
    out$status <- "too_few_cpgs"
    return(out)
  }
  Xo <- X[used$cpg, , drop = FALSE]
  if (qr(Xo)$rank < K) {
    out$status <- "collinear_reference"
    return(out)
  }

  y <- used$n_mod / used$n_valid
  b <- project_one(y, Xo, constrained = !opts$unconstrained)
  if (anyNA(b)) {
    out$status <- "collinear_reference"
    return(out)
  }

  out$est  <- b
  out$rmse <- sqrt(mean((y - drop(Xo %*% b))^2))
  out$cond <- kappa(Xo, exact = TRUE)
  if (opts$nboot > 0) {
    out$se <- sampling_se(used$n_valid, used$n_mod, Xo, opts$nboot,
                          constrained = !opts$unconstrained)
  }
  out
}


# ---------------------------------------------------------------------------
# Quality control plot and run record
# ---------------------------------------------------------------------------

plot_qc <- function(res, opts, min_needed, file) {
  grDevices::pdf(file, width = 11, height = 5)
  tryCatch({
    par(mfrow = c(1, 2))
    flagged <- res$flag != "PASS"
    cols <- ifelse(flagged, "firebrick", "grey30")

    ord <- order(res$n_cpgs_used)
    plot(seq_along(ord), res$n_cpgs_used[ord], pch = 19, col = cols[ord],
         xlab = "Samples, ordered by CpGs used",
         ylab = "Reference CpGs used",
         main = "CpG coverage per sample",
         ylim = c(0, max(res$n_ref_cpgs[1], res$n_cpgs_used, na.rm = TRUE)))
    abline(h = min_needed, lty = 2)
    abline(h = opts$min_prop * res$n_ref_cpgs[1], lty = 3)
    legend("topleft", bty = "n", pch = c(19, 19, NA, NA),
           lty = c(NA, NA, 2, 3), col = c("grey30", "firebrick", "black", "black"),
           legend = c("PASS", "FLAG", "minimum CpGs", "minimum proportion"))

    if (any(!is.na(res$max_se))) {
      plot(res$n_cpgs_used, res$max_se, pch = 19, col = cols,
           xlab = "Reference CpGs used",
           ylab = "Largest standard error across cell types",
           main = "Sampling error")
      abline(h = opts$max_se, lty = 2)
    } else {
      plot.new()
      text(0.5, 0.5, "No standard errors available")
    }
  }, error = function(e) message("QC plot could not be drawn: ", conditionMessage(e)),
  finally = grDevices::dev.off())
}

write_run_settings <- function(opts, ref, file) {
  lines <- c(
    paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("R version:", R.version.string),
    paste("data.table version:", as.character(packageVersion("data.table"))),
    paste("quadprog version:", as.character(packageVersion("quadprog"))),
    paste("reference CpGs:", nrow(ref$X)),
    paste("cell types:", paste(ref$cell_types, collapse = ", ")),
    "",
    "settings:",
    paste0("  ", names(opts), " = ", unlist(lapply(opts, as.character)))
  )
  writeLines(lines, file)
}


# ---------------------------------------------------------------------------
# Main routine
# ---------------------------------------------------------------------------

run_cell_proportions <- function(opts) {
  check_options(opts)
  set.seed(opts$seed)

  ref <- read_reference(opts$reference)
  X <- ref$X
  K <- ncol(X)
  n_ref <- nrow(X)
  message("Reference: ", n_ref, " CpGs, ", K, " cell types (",
          paste(ref$cell_types, collapse = ", "), ")")

  inp <- read_all_input(opts, ref$sites)
  samples <- inp$samples
  by_sample <- if (nrow(inp$data) > 0) split(inp$data, by = "sample", keep.by = FALSE) else list()
  message("Reference CpGs with data in at least one sample: ",
          uniqueN(inp$data$cpg), " of ", n_ref)

  n <- length(samples)
  est_mat <- se_mat <- matrix(NA_real_, n, K,
                              dimnames = list(NULL, ref$cell_types))
  qc <- data.frame(n_cpgs_covered = integer(n), n_cpgs_used = integer(n),
                   median_depth_used = NA_real_, fit_rmse = NA_real_,
                   condition_number = NA_real_)
  flag_reason <- character(n)

  min_needed <- max(opts$min_cpgs, K + 1)

  for (i in seq_len(n)) {
    d <- by_sample[[samples[i]]]
    if (is.null(d)) d <- data.table(cpg = character(), n_valid = numeric(), n_mod = numeric())
    used <- d[n_valid >= opts$min_depth]

    fit <- estimate_sample(used, X, opts)

    qc$n_cpgs_covered[i] <- nrow(d)
    qc$n_cpgs_used[i] <- nrow(used)
    if (nrow(used) > 0) qc$median_depth_used[i] <- median(used$n_valid)
    qc$fit_rmse[i] <- fit$rmse
    qc$condition_number[i] <- fit$cond
    est_mat[i, ] <- fit$est
    se_mat[i, ] <- fit$se

    reasons <- character()
    if (fit$status != "ok") reasons <- c(reasons, fit$status)
    if (nrow(used) / n_ref < opts$min_prop) reasons <- c(reasons, "low_cpg_proportion")
    if (fit$status == "ok" && opts$nboot > 0 && any(!is.na(fit$se)) &&
        max(fit$se, na.rm = TRUE) > opts$max_se) {
      reasons <- c(reasons, "high_sampling_error")
    }
    flag_reason[i] <- if (length(reasons) == 0) "" else paste(reasons, collapse = ";")
  }

  flag <- ifelse(flag_reason == "", "PASS", "FLAG")
  if (opts$mask_flagged) est_mat[flag == "FLAG", ] <- NA_real_

  max_se <- apply(se_mat, 1, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
  colnames(se_mat) <- paste0("se_", ref$cell_types)

  res <- data.frame(
    sample = samples,
    flag = flag,
    flag_reason = flag_reason,
    est_mat,
    sum_proportions = rowSums(est_mat),
    fit_rmse = qc$fit_rmse,
    n_ref_cpgs = n_ref,
    n_cpgs_covered = qc$n_cpgs_covered,
    n_cpgs_used = qc$n_cpgs_used,
    prop_cpgs_used = qc$n_cpgs_used / n_ref,
    median_depth_used = qc$median_depth_used,
    condition_number = qc$condition_number,
    se_mat,
    max_se = max_se,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  out_dir <- dirname(opts$output)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  write.table(res, paste0(opts$output, "_cell_proportions.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  plot_qc(res, opts, min_needed, paste0(opts$output, "_coverage_qc.pdf"))
  write_run_settings(opts, ref, paste0(opts$output, "_run_settings.txt"))

  # Summary for the console
  n_flag <- sum(flag == "FLAG")
  message("\n", n, " samples: ", n - n_flag, " PASS, ", n_flag, " FLAG")
  if (n_flag > 0) {
    why <- table(unlist(strsplit(flag_reason[flag == "FLAG"], ";")))
    for (r in names(why)) message("  ", r, ": ", why[[r]])
  }
  message("CpGs used per sample: median ", median(qc$n_cpgs_used),
          ", range ", min(qc$n_cpgs_used), " to ", max(qc$n_cpgs_used))
  if (n_flag > n / 2) {
    message("More than half the samples are flagged. Check the thresholds against ",
            "the coverage plot before using them as they stand.")
  }
  message("Results written with prefix ", opts$output)

  invisible(res)
}


print_usage <- function() {
  cat("Usage: Rscript ont_cell_proportions.R --input FILE --reference FILE [--output PREFIX] [options]\n",
      "See the comments at the top of this script for the formats and all options.\n",
      "Options: --input_type --scale --mod_codes --min_depth --min_cpgs --min_prop\n",
      "         --max_se --nboot --seed --unconstrained --mask_flagged --help\n", sep = "")
}

main <- function() {
  opts <- parse_args(default_options())
  if (opts$help) {
    print_usage()
    quit(status = 0)
  }
  run_cell_proportions(opts)
}

# Runs only when the file is executed with Rscript, not when it is source()d.
if (sys.nframe() == 0L) main()
