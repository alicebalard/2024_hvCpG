#!/usr/bin/env python3
"""
cpg_utils.py — Shared core library for CpG methylation data preparation
========================================================================
Contains all reusable logic shared across atlas and array-based pipelines:
  - CpG name loading (BED, HDF5)
  - Beta file loading (.beta WGBS binary format)
  - RDS matrix loading (Illumina array data via pyreadr)
  - Coverage-aware beta matrix construction (chunked, memory-efficient)
  - Per-dataset median SD and lambda — computed AFTER final CpG selection
  - HDF5 / TSV output writing
  - Metadata filtering and sampling helpers

Design principle — statistics ordering
---------------------------------------
Median SD and lambda are ALWAYS computed on the final filtered CpG set
(i.e. after the min_datasets site-selection mask has been applied).
compute_stats_on_filtered_matrix() is the single point of truth for this.

Author: Alice Balard
"""

import os
import re
import gzip
import numpy as np
import pandas as pd
import h5py
import bottleneck as bn


# ══════════════════════════════════════════════════════════════════════════════
#  CpG coordinate / name loading
# ══════════════════════════════════════════════════════════════════════════════

def load_cpg_names_from_bed(cpg_bed):
    """
    Load CpG site names from a (optionally gzipped) BED file.
    Each row must have at least two tab-separated columns: chrom and pos.
    Returns strings formatted as 'chr1_12345'.
    """
    cpg_names = []
    opener = gzip.open(cpg_bed, "rt") if cpg_bed.endswith(".gz") else open(cpg_bed, "rt")
    with opener as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom, pos = parts[0], int(parts[1])
            cpg_names.append(f"{chrom}_{pos}")
    print(f"  Loaded {len(cpg_names):,} CpG sites from {os.path.basename(cpg_bed)}")
    return cpg_names


def load_exclude_sites(path):
    """
    Load a list of CpG sites to EXCLUDE from any analysis.

    The file must contain one site per line in 'chr_pos' format
    (e.g. 'chr1_123456'), matching the identifiers used throughout
    the pipeline.  Blank lines and lines starting with '#' are ignored.

    Sites in this list are masked out at the earliest possible point:
      - In build_beta_matrix_chunked: the exclude_mask is applied before
        Pass 1, so excluded sites never accumulate dataset-coverage counts
        and cannot appear in the final matrix.
      - In prepare_arrays.py: the mask is applied before the NA filter.

    Parameters
    ----------
    path : str  -- path to the exclusion list (plain text, one site per line)

    Returns
    -------
    frozenset of str  -- site identifiers to exclude
    """
    sites = set()
    opener = gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")
    with opener as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            sites.add(line)
    print(f"  Loaded {len(sites):,} sites to exclude from: {os.path.basename(path)}")
    return frozenset(sites)


def build_exclude_mask(cpg_names, exclude_sites):
    """
    Build a boolean array marking sites to EXCLUDE.

    Parameters
    ----------
    cpg_names     : list of str  -- full reference CpG list (length = nr_sites)
    exclude_sites : frozenset    -- output of load_exclude_sites()

    Returns
    -------
    np.ndarray bool, shape (nr_sites,)
        True  = exclude this site
        False = keep this site
    """
    if not exclude_sites:
        return np.zeros(len(cpg_names), dtype=bool)
    mask = np.array([c in exclude_sites for c in cpg_names], dtype=bool)
    n_excl = int(mask.sum())
    print(f"  build_exclude_mask: {n_excl:,} / {len(cpg_names):,} sites will be excluded.")
    return mask


def load_cpg_names_from_h5(h5_path, dataset_key="cpg_names"):
    """Load CpG site names stored inside an HDF5 file."""
    with h5py.File(h5_path, "r") as f:
        raw   = f[dataset_key][:]
        names = [s.decode("utf-8") if isinstance(s, bytes) else s for s in raw]
    print(f"  Loaded {len(names):,} CpG names from {os.path.basename(h5_path)}")
    return names


# ══════════════════════════════════════════════════════════════════════════════
#  Beta file I/O  (.beta WGBS binary format)
# ══════════════════════════════════════════════════════════════════════════════

def load_beta_file(path):
    """
    Load a WGBS .beta file -> (beta_values, coverage) arrays.
    Format: flat uint8 binary, interleaved [meth, cov] pairs.
    Positions with cov == 0 returned as NaN in beta.
    """
    arr  = np.fromfile(path, dtype=np.uint8).reshape((-1, 2))
    meth = arr[:, 0].astype(np.float32)
    cov  = arr[:, 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = np.where(cov == 0, np.nan, meth / cov).astype(np.float32)
    return beta, cov


def build_sample_to_path_map(beta_files, id_pattern=r"-([A-Za-z0-9_]+)\.hg38\.beta$"):
    """
    Build a mapping  short_sample_id -> beta_file_path.
    id_pattern: regex with ONE capture group extracting the sample ID.
    Default matches Loyfer atlas naming: '-SAMPLEID.hg38.beta'
    """
    sample_to_path = {}
    dup_count = 0
    for fn in beta_files:
        m = re.search(id_pattern, os.path.basename(fn))
        if m:
            key = m.group(1)
            if key in sample_to_path:
                dup_count += 1
            sample_to_path[key] = fn
    if dup_count:
        print(f"  Warning: {dup_count} duplicate sample IDs; keeping last occurrence.")
    print(f"  Mapped {len(sample_to_path):,} unique sample IDs -> beta files.")
    return sample_to_path


# ══════════════════════════════════════════════════════════════════════════════
#  Metadata helpers — Analysis group construction
# ══════════════════════════════════════════════════════════════════════════════

def add_analysis_group(df, tissue_col="Source Tissue", cell_col="Cell type", sep=" - "):
    """Add 'Analysis group' = tissue + sep + cell_type."""
    for col in [tissue_col, cell_col]:
        if col not in df.columns:
            raise ValueError(f"Missing required column: '{col}'")
    df = df.copy()
    df["Analysis group"] = df[tissue_col].astype(str) + sep + df[cell_col].astype(str)
    return df


# ══════════════════════════════════════════════════════════════════════════════
#  Metadata helpers — row-level filters
# ══════════════════════════════════════════════════════════════════════════════

def filter_by_column_value(df, column, value, keep=True):
    """Keep (keep=True) or drop (keep=False) rows where column == value."""
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in metadata.")
    mask   = df[column] == value
    result = df[mask] if keep else df[~mask]
    action = "kept" if keep else "dropped"
    print(f"  {action} rows where {column}=={value!r} -> {len(result):,} remain.")
    return result


def deduplicate_by_column(df, column, keep="first"):
    """
    Drop duplicate rows so each value of `column` appears only once.
    Used for removing multiple samples per individual (script 2).
      e.g. deduplicate_by_column(df, 'PatientID')
    keep: 'first' (default) or 'last'
    """
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found for deduplication.")
    before = len(df)
    df = df.drop_duplicates(subset=column, keep=keep)
    print(f"  deduplicate_by_column('{column}'): "
          f"removed {before - len(df):,} rows -> {len(df):,} remain.")
    return df


# ══════════════════════════════════════════════════════════════════════════════
#  Metadata helpers — group-level filters
# ══════════════════════════════════════════════════════════════════════════════

def filter_groups_by_column(df, group_col="Analysis group",
                             filter_col="sex", required_value="M"):
    """Keep only groups where ALL rows in filter_col equal required_value."""
    filtered = df.groupby(group_col).filter(
        lambda g: (g[filter_col].astype(str).str.strip() == str(required_value)).all()
    )
    print(f"  filter_groups_by_column: {filtered[group_col].nunique()} groups remain "
          f"(all {filter_col}=={required_value!r}).")
    return filtered


def filter_groups_mixed_sex(df, group_col="Analysis group", sex_col="sex"):
    """
    Keep only groups containing BOTH at least one 'M' and one 'F'.
    Used for bothsexes analysis (script 6).
    """
    filtered = df.groupby(group_col).filter(
        lambda g: (
            ("M" in g[sex_col].astype(str).str.strip().values) and
            ("F" in g[sex_col].astype(str).str.strip().values)
        )
    )
    print(f"  filter_groups_mixed_sex: {filtered[group_col].nunique()} mixed-sex groups remain.")
    return filtered


def filter_valid_groups(meta, group_col="Analysis group",
                         sample_col="Sample name", min_samples=3):
    """Return {group_name: [sample_names]} for groups with >= min_samples."""
    counts = meta.groupby(group_col).size()
    valid  = counts[counts >= min_samples].index.tolist()
    print(f"  {len(valid)} groups pass min_samples={min_samples} (of {len(counts)} total).")
    return {g: meta.loc[meta[group_col] == g, sample_col].tolist() for g in valid}


# ══════════════════════════════════════════════════════════════════════════════
#  Metadata helpers — sampling
# ══════════════════════════════════════════════════════════════════════════════

def normalize_sex(s):
    """Normalise sex values to 'M' or 'F'; return pd.NA if unrecognised."""
    if pd.isna(s):
        return pd.NA
    s = str(s).strip().upper()
    if s in {"M", "MALE"}:
        return "M"
    if s in {"F", "FEMALE"}:
        return "F"
    return pd.NA


def sample_n_groups(df, n, group_col="Analysis group", min_per_group=3, seed=42):
    """
    Randomly sample up to n groups that each have >= min_per_group rows.
    Used by scripts 5, 6, 11, 12.2, 13.2.

    Parameters
    ----------
    n             : maximum number of groups to keep
    min_per_group : minimum rows required for eligibility
    seed          : random seed for reproducibility
    """
    counts      = df.groupby(group_col).size()
    eligible    = counts[counts >= min_per_group].index.tolist()
    n_to_sample = min(n, len(eligible))
    rng         = np.random.default_rng(seed)
    sampled     = list(rng.choice(eligible, size=n_to_sample, replace=False))
    result      = df[df[group_col].isin(sampled)]
    print(f"  sample_n_groups: sampled {len(sampled)}/{len(eligible)} eligible groups "
          f"(>={min_per_group} samples each) -> {len(result):,} rows.")
    return result


def sample_pairs(df, pair_type, group_col="Analysis group",
                  sex_col="sex", n_groups=None, seed=42):
    """
    For each eligible Analysis group, sample exactly 2 individuals
    according to pair_type:
      'MM' -- 2 males   (group needs >= 2 males)
      'FF' -- 2 females (group needs >= 2 females)
      'MF' -- 1 male + 1 female (group needs >= 1 of each)

    Implements scripts 15 (MM), 16 (FF), 17 (MF).

    Parameters
    ----------
    pair_type : 'MM', 'FF', or 'MF'
    n_groups  : cap on number of groups (None = all eligible)
    seed      : random seed
    """
    if pair_type not in {"MM", "FF", "MF"}:
        raise ValueError(f"pair_type must be 'MM', 'FF', or 'MF'; got '{pair_type}'")

    rng = np.random.default_rng(seed)
    df  = df.copy()
    df[sex_col] = df[sex_col].map(normalize_sex)

    eligible_groups = []
    for grp, g in df.groupby(group_col, dropna=False):
        males   = g[g[sex_col] == "M"]
        females = g[g[sex_col] == "F"]
        if pair_type == "MM" and len(males) >= 2:
            eligible_groups.append(grp)
        elif pair_type == "FF" and len(females) >= 2:
            eligible_groups.append(grp)
        elif pair_type == "MF" and len(males) >= 1 and len(females) >= 1:
            eligible_groups.append(grp)

    if n_groups is not None and n_groups < len(eligible_groups):
        eligible_groups = list(rng.choice(eligible_groups, size=n_groups, replace=False))

    selected = []
    for grp in eligible_groups:
        g       = df[df[group_col] == grp]
        males   = g[g[sex_col] == "M"]
        females = g[g[sex_col] == "F"]
        if pair_type == "MM":
            idx = list(rng.choice(males.index, size=2, replace=False))
        elif pair_type == "FF":
            idx = list(rng.choice(females.index, size=2, replace=False))
        else:  # MF
            idx = list(rng.choice(males.index, size=1, replace=False)) + \
                  list(rng.choice(females.index, size=1, replace=False))
        selected.append(df.loc[idx])

    if selected:
        result = pd.concat(selected).sort_values(group_col)
    else:
        result = pd.DataFrame(columns=df.columns)

    print(f"  sample_pairs({pair_type}): {result[group_col].nunique()} groups, "
          f"{len(result)} samples.")
    return result


def shorten_sample_ids(samples_per_group, sep="-", part=-1):
    """Strip a prefix from sample names: 'GSE123-SAMPLEID' -> 'SAMPLEID'."""
    return {g: [s.split(sep)[part] for s in samples]
            for g, samples in samples_per_group.items()}


# ══════════════════════════════════════════════════════════════════════════════
#  Statistics — always on the FINAL filtered CpG set
# ══════════════════════════════════════════════════════════════════════════════

def compute_median_and_lambda(row_sds, lambda_percentile=95.0):
    """
    Compute median SD and lambda = perc(lambda_percentile) / median.
    Always call on the final filtered CpG set, never on raw unfiltered data.
    """
    if row_sds is None or np.all(np.isnan(row_sds)):
        return np.nan, np.nan
    median_sd    = float(np.nanmedian(row_sds))
    perc_upper   = float(np.nanpercentile(row_sds, lambda_percentile))
    lambda_value = (
        perc_upper / median_sd
        if median_sd and np.isfinite(median_sd) and median_sd > 0
        else np.nan
    )
    return median_sd, lambda_value


def compute_stats_on_filtered_matrix(matrix, sample_groups, lambda_percentile=95.0):
    """
    Compute per-dataset median SD and lambda on the ALREADY-FILTERED matrix.

    Must be called AFTER the final CpG site-selection mask has been applied.
    This is enforced by every entry-point script.

    Parameters
    ----------
    matrix        : np.ndarray (n_selected_cpgs x n_samples), float32
    sample_groups : list of str -- one group label per column
    """
    groups        = list(dict.fromkeys(sample_groups))
    group_medians = {}
    group_lambdas = {}
    for group in groups:
        col_idx  = [i for i, g in enumerate(sample_groups) if g == group]
        submat   = matrix[:, col_idx]
        row_sds  = bn.nanstd(submat, axis=1)
        med, lam = compute_median_and_lambda(row_sds, lambda_percentile)
        group_medians[group] = med
        group_lambdas[group] = lam
        print(f"    {group}: median_sd={med:.4f}, lambda={lam:.4f}")
    return group_medians, group_lambdas


# ══════════════════════════════════════════════════════════════════════════════
#  Core WGBS matrix builder  (two-pass, chunked, coverage-aware)
# ══════════════════════════════════════════════════════════════════════════════

def build_beta_matrix_chunked(
    samples_per_group_short, sample_to_path, nr_sites,
    min_cov=10, min_samples_per_group=3, min_datasets=1,
    chunk_size=100_000, lambda_percentile=95.0,
    exclude_mask=None,
):
    """
    Two-pass WGBS beta matrix builder.

    Pass 1  -- per-site coverage counting only (no statistics computed).
    Selection -- keep CpGs covered in >= min_datasets groups.
    Pass 2  -- build final float32 matrix (selected sites only).
    Stats   -- compute_stats_on_filtered_matrix() on the Pass 2 output.

    Parameters
    ----------
    exclude_mask : np.ndarray bool, shape (nr_sites,), optional
        Sites marked True are excluded from the earliest possible point.
        Their dataset_pass_count slots are permanently zeroed so they
        can never satisfy >= min_datasets and will not appear in the
        output matrix or statistics.
        Build with: build_exclude_mask(cpg_names, load_exclude_sites(path))

    Returns
    -------
    final_matrix, group_medians, group_lambdas,
    all_sample_names, all_sample_groups, final_mask
    """
    dataset_pass_count = np.zeros(nr_sites, dtype=np.uint16)
    all_sample_names  = []
    all_sample_groups = []

    # Zero-out excluded sites before Pass 1 so they can never
    # accumulate coverage counts and will be absent from final_mask.
    if exclude_mask is not None and np.any(exclude_mask):
        dataset_pass_count[exclude_mask] = 0  # already 0; made explicit

    # Pass 1
    print("\n  -- Pass 1: counting per-site dataset coverage --")
    for group, samples in samples_per_group_short.items():
        group_paths = []
        for s in samples:
            path = sample_to_path.get(s)
            if path is None:
                print(f"    Warning: no beta file for sample '{s}' -- skipping.")
                continue
            group_paths.append(path)
            all_sample_names.append(s)
            all_sample_groups.append(group)

        if not group_paths:
            print(f"    Warning: {group}: no valid samples, skipped.")
            continue

        for start in range(0, nr_sites, chunk_size):
            end      = min(start + chunk_size, nr_sites)
            mat_cols = []
            for path in group_paths:
                mm      = np.memmap(path, dtype=np.uint8, mode="r").reshape(-1, 2)
                meth_sl = mm[start:end, 0].astype(np.float32, copy=False)
                cov_sl  = mm[start:end, 1]
                col     = np.empty(end - start, dtype=np.float32)
                np.divide(meth_sl, cov_sl, out=col, where=(cov_sl != 0))
                col[cov_sl < min_cov] = np.nan
                mat_cols.append(col)

            mat_chunk    = np.column_stack(mat_cols)
            valid_counts = np.sum(~np.isnan(mat_chunk), axis=1)
            covered      = valid_counts >= min_samples_per_group
            # Apply increment only to non-excluded sites
            if exclude_mask is not None:
                covered &= ~exclude_mask[start:end]
            dataset_pass_count[start:end][covered] += 1

    # Site selection — excluded sites have count==0 so they fail this test
    final_mask = dataset_pass_count >= min_datasets
    if exclude_mask is not None:
        final_mask &= ~exclude_mask   # belt-and-braces: ensure excluded never selected
    n_sel  = int(np.sum(final_mask))
    n_excl = int(np.sum(exclude_mask)) if exclude_mask is not None else 0
    n_cols = len(all_sample_names)
    print(f"\n  Selected {n_sel:,} CpGs "
          f"(>={min_samples_per_group} samples in >={min_datasets} datasets"
          + (f", {n_excl:,} sites excluded by exclusion list" if n_excl else "") + ").")

    # Pass 2
    print("  -- Pass 2: building filtered beta matrix --")
    final_matrix = np.full((n_sel, n_cols), np.nan, dtype=np.float32)
    col_idx = 0
    for group, samples in samples_per_group_short.items():
        for s in samples:
            path = sample_to_path.get(s)
            if path is None:
                col_idx += 1
                continue
            beta, cov = load_beta_file(path)
            if len(beta) != nr_sites:
                raise ValueError(
                    f"Size mismatch '{s}': got {len(beta):,} CpGs, expected {nr_sites:,}"
                )
            beta[cov < min_cov] = np.nan
            final_matrix[:, col_idx] = beta[final_mask]
            col_idx += 1

    # Statistics on final filtered matrix
    print("\n  -- Computing median SD and lambda on final filtered CpGs --")
    group_medians, group_lambdas = compute_stats_on_filtered_matrix(
        final_matrix, all_sample_groups, lambda_percentile
    )

    return (final_matrix, group_medians, group_lambdas,
            all_sample_names, all_sample_groups, final_mask)


# ══════════════════════════════════════════════════════════════════════════════
#  HDF5 / TSV output
# ══════════════════════════════════════════════════════════════════════════════

def write_h5_matrix(output_path, matrix, cpg_names, sample_names,
                     sample_groups, compression="gzip"):
    """Write beta matrix + metadata to HDF5."""
    with h5py.File(output_path, "w") as h5f:
        h5f.create_dataset("matrix",        data=matrix,
                           dtype="float32",  compression=compression)
        h5f.create_dataset("cpg_names",     data=np.array(cpg_names, dtype="S"))
        h5f.create_dataset("samples",       data=sample_names,
                           dtype=h5py.string_dtype(encoding="utf-8"))
        h5f.create_dataset("sample_groups", data=sample_groups,
                           dtype=h5py.string_dtype(encoding="utf-8"))
    print(f"  Saved HDF5 ({matrix.shape[0]:,} CpGs x {matrix.shape[1]:,} samples) -> {output_path}")


def write_medsd_lambda_tsv(output_path, group_medians, group_lambdas):
    """Save per-dataset median SD and lambda to TSV."""
    pd.DataFrame({
        "dataset":   list(group_medians.keys()),
        "median_sd": list(group_medians.values()),
        "lambda":    [group_lambdas[k] for k in group_medians],
    }).to_csv(output_path, sep="\t", index=False)
    print(f"  Saved median_sd / lambda -> {output_path}")


def write_metadata_tsv(output_path, sample_names, sample_groups):
    """Save sample-level metadata to TSV."""
    pd.DataFrame({"sample": sample_names, "dataset": sample_groups}).to_csv(
        output_path, sep="\t", index=False
    )
    print(f"  Saved sample metadata -> {output_path}")
