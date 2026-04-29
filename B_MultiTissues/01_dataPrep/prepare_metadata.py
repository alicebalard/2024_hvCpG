#!/usr/bin/env python3
"""
prepare_metadata.py — Flexible metadata preparation for CpG methylation pipelines
==================================================================================
Single script replacing all 18 numbered prepare_metadata_*.py scripts.

Every combination of options maps to one of the original scripts:

  Script 1  (byDevLayer)            --group_col "Germ layer"
  Script 2  (rmMultSamples)         --dedup_col "PatientID"
  Script 3  (correspMariaTissues)   --keep_col "Found in Maria DS?" --keep_val True
  Script 4  (maleOnly)              --sex_filter M
  Script 5  (femaleOnly6gp)         --sex_filter F --sample_n_groups 6
  Script 6  (bothsexes6gp)          --mixed_sex --sample_n_groups 6
  Script 8  (byTissue)              --group_col "Group simplified"
  Script 9  (immuneOnly)            --keep_col "Immune?" --keep_val True
  Script 10 (noImmune)              --exclude_immune
  Script 11 (noImmune_sample11gp)   --exclude_immune --sample_n_groups 11
  Script 12 (endo)                  --germ_layer Endo
  Script 12.2 (endo6gp)             --germ_layer Endo --sample_n_groups 6
  Script 13 (meso)                  --germ_layer Meso
  Script 13.2 (meso6gp)             --germ_layer Meso --sample_n_groups 6
  Script 14 (ecto)                  --germ_layer Ecto
  Script 15 (pairs_MM)              --pairs MM
  Script 16 (pairs_FF)              --pairs FF
  Script 17 (pairs_MF)              --pairs MF

All options can be freely combined and applied in sequence.

Usage examples
--------------
# Script 1 equivalent — group by developmental layer:
python prepare_metadata.py --meta SupTab1.csv --group_col "Germ layer"

# Script 2 — remove multiple samples per individual:
python prepare_metadata.py --meta SupTab1.csv --dedup_col "PatientID"

# Script 4 — male-only groups:
python prepare_metadata.py --meta SupTab1.csv --sex_filter M

# Script 5 — female-only, sample 6 groups:
python prepare_metadata.py --meta SupTab1.csv --sex_filter F --sample_n_groups 6

# Script 6 — both sexes, sample 6 mixed-sex groups:
python prepare_metadata.py --meta SupTab1.csv --mixed_sex --sample_n_groups 6

# Script 11 — no immune, sample 11 groups:
python prepare_metadata.py --meta SupTab1.csv --exclude_immune --sample_n_groups 11

# Script 12 — endoderm only:
python prepare_metadata.py --meta SupTab1.csv --germ_layer Endo

# Script 12.2 — endoderm, sample 6 groups:
python prepare_metadata.py --meta SupTab1.csv --germ_layer Endo --sample_n_groups 6

# Script 15 — paired MM:
python prepare_metadata.py --meta SupTab1.csv --pairs MM

# Script 17 — paired MF, limit to 10 groups:
python prepare_metadata.py --meta SupTab1.csv --pairs MF --n_groups 10

Author: Alice Balard
"""

import argparse
import pandas as pd
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cpg_utils import (
    add_analysis_group,
    set_analysis_group_from_column,
    filter_by_column_value,
    filter_groups_by_column,
    filter_groups_mixed_sex,
    deduplicate_by_column,
    sample_n_groups,
    sample_pairs,
)

# ──────────────────────────────────────────────
#  Arguments
# ──────────────────────────────────────────────

parser = argparse.ArgumentParser(
    description="Prepare and filter sample metadata for CpG methylation pipelines.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)

parser.add_argument("--meta",   required=True, help="Input metadata CSV file.")
parser.add_argument("--output", required=False, help="Output CSV path (auto-named if omitted).")

# ── Analysis group construction ───────────────────────────────────────────────
group_def = parser.add_mutually_exclusive_group()
group_def.add_argument(
    "--group_col", default=None, metavar="COL",
    help="Use a single existing column as 'Analysis group' directly "
         "(e.g. 'Germ layer' for script 1, 'Group simplified' for script 8). "
         "Mutually exclusive with --tissue_col / --cell_col.")
group_def.add_argument(
    "--tissue_col", default=None, metavar="COL",
    help="Tissue column for 'Analysis group' = tissue + ' - ' + cell_type "
         "(default: 'Source Tissue'). Mutually exclusive with --group_col.")

parser.add_argument("--cell_col", default="Cell type",
                    help="Cell-type column (used with --tissue_col; default: 'Cell type').")
parser.add_argument("--sample_col", default="Sample name",
                    help="Column for sample identifiers (default: 'Sample name').")

# ── Row-level filters ─────────────────────────────────────────────────────────
parser.add_argument("--exclude_immune", action="store_true",
                    help="Drop rows where 'Immune?' == True. (Scripts 10, 11)")
parser.add_argument("--keep_col", action="append", default=[], metavar="COL",
                    help="Keep rows where this column == --keep_val. Repeatable.")
parser.add_argument("--keep_val", action="append", default=[], metavar="VAL",
                    help="Value for --keep_col filter. Repeatable.")
parser.add_argument("--drop_col", action="append", default=[], metavar="COL",
                    help="Drop rows where this column == --drop_val. Repeatable.")
parser.add_argument("--drop_val", action="append", default=[], metavar="VAL",
                    help="Value for --drop_col filter. Repeatable.")
parser.add_argument("--germ_layer", default=None, metavar="LAYER",
                    help="Keep rows where 'Germ layer' == LAYER. "
                         "One of: Endo, Meso, Ecto. (Scripts 12–14)")
parser.add_argument("--dedup_col", default=None, metavar="COL",
                    help="Deduplicate keeping first row per unique value of COL. "
                         "E.g. 'PatientID' to keep one sample per individual. (Script 2)")

# ── Group-level filters ───────────────────────────────────────────────────────
parser.add_argument("--sex_filter", default=None, metavar="SEX",
                    help="Keep groups where ALL samples have this sex: M or F. (Scripts 4, 5)")
parser.add_argument("--mixed_sex", action="store_true",
                    help="Keep only groups with BOTH at least one M and one F. (Script 6)")
parser.add_argument("--keep_groups_col", default=None, metavar="COL",
                    help="Keep groups where ALL rows have this column == --keep_groups_val.")
parser.add_argument("--keep_groups_val", default=None, metavar="VAL",
                    help="Value for --keep_groups_col filter.")

# ── Sampling ──────────────────────────────────────────────────────────────────
parser.add_argument("--sample_n_groups", type=int, default=None, metavar="N",
                    help="Randomly sample up to N groups (each with >= --min_per_group rows). "
                         "Used by scripts 5, 6, 11, 12.2, 13.2.")
parser.add_argument("--min_per_group", type=int, default=3,
                    help="Minimum rows per group to be eligible for sampling (default: 3).")
parser.add_argument("--pairs", default=None, metavar="TYPE",
                    choices=["MM", "FF", "MF"],
                    help="Sample exactly 2 individuals per group according to pair type: "
                         "MM (2 males), FF (2 females), MF (1 male + 1 female). (Scripts 15-17)")
parser.add_argument("--n_groups", type=int, default=None, metavar="N",
                    help="Cap on number of groups used with --pairs (default: all eligible).")
parser.add_argument("--seed", type=int, default=42,
                    help="Random seed for all sampling operations (default: 42).")

args = parser.parse_args()

# Validate paired args
if len(args.keep_col) != len(args.keep_val):
    parser.error("--keep_col and --keep_val must be used in equal numbers.")
if len(args.drop_col) != len(args.drop_val):
    parser.error("--drop_col and --drop_val must be used in equal numbers.")
if bool(args.keep_groups_col) != bool(args.keep_groups_val):
    parser.error("--keep_groups_col and --keep_groups_val must be used together.")
if args.mixed_sex and args.sex_filter:
    parser.error("--mixed_sex and --sex_filter are mutually exclusive.")
if args.pairs and args.sample_n_groups:
    parser.error("--pairs and --sample_n_groups are mutually exclusive "
                 "(use --n_groups to cap groups in --pairs mode).")

# ──────────────────────────────────────────────
#  Load
# ──────────────────────────────────────────────

print(f"\nLoading: {args.meta}")
df = pd.read_csv(args.meta)
print(f"  {len(df):,} rows, {df.shape[1]} columns.")

# ──────────────────────────────────────────────
#  Step 1 — Build Analysis group
# ──────────────────────────────────────────────

if args.group_col:
    # Use a single column directly (scripts 1, 8)
    df = set_analysis_group_from_column(df, args.group_col)
    print(f"  Analysis group <- '{args.group_col}' "
          f"({df['Analysis group'].nunique()} unique groups).")
else:
    # Combine tissue + cell type (all other scripts)
    tissue_col = args.tissue_col if args.tissue_col else "Source Tissue"
    df = add_analysis_group(df, tissue_col=tissue_col, cell_col=args.cell_col)
    print(f"  Analysis group <- '{tissue_col}' + '{args.cell_col}' "
          f"({df['Analysis group'].nunique()} unique groups).")

# ──────────────────────────────────────────────
#  Step 2 — Deduplication (script 2)
# ──────────────────────────────────────────────

if args.dedup_col:
    df = deduplicate_by_column(df, args.dedup_col)

# ──────────────────────────────────────────────
#  Step 3 — Row-level filters
# ──────────────────────────────────────────────

# --exclude_immune (scripts 10, 11)
if args.exclude_immune:
    if "Immune?" not in df.columns:
        print("  Warning: --exclude_immune requested but 'Immune?' column not found -- skipping.")
    else:
        before = len(df)
        df = df[df["Immune?"] != True]
        print(f"  --exclude_immune: removed {before - len(df):,} rows -> {len(df):,} remain.")

# --germ_layer (scripts 12, 12.2, 13, 13.2, 14)
if args.germ_layer:
    if "Germ layer" not in df.columns:
        parser.error("--germ_layer requires a 'Germ layer' column in the metadata.")
    df = filter_by_column_value(df, "Germ layer", args.germ_layer, keep=True)

# Generic --keep_col / --keep_val pairs (scripts 3, 9)
for col, val in zip(args.keep_col, args.keep_val):
    if col in df.columns:
        try:
            val = df[col].dtype.type(val)
        except (ValueError, AttributeError):
            pass
    df = filter_by_column_value(df, col, val, keep=True)

# Generic --drop_col / --drop_val pairs
for col, val in zip(args.drop_col, args.drop_val):
    if col in df.columns:
        try:
            val = df[col].dtype.type(val)
        except (ValueError, AttributeError):
            pass
    df = filter_by_column_value(df, col, val, keep=False)

# ──────────────────────────────────────────────
#  Step 4 — Group-level filters
# ──────────────────────────────────────────────

# --sex_filter M|F (scripts 4, 5)
if args.sex_filter:
    df["sex"] = df["sex"].astype(str).str.strip()
    df = filter_groups_by_column(df, filter_col="sex", required_value=args.sex_filter)

# --mixed_sex (script 6)
if args.mixed_sex:
    df["sex"] = df["sex"].astype(str).str.strip()
    df = filter_groups_mixed_sex(df)

# --keep_groups_col / --keep_groups_val
if args.keep_groups_col:
    df = filter_groups_by_column(df,
                                  filter_col=args.keep_groups_col,
                                  required_value=args.keep_groups_val)

# ──────────────────────────────────────────────
#  Step 5 — Sampling
# ──────────────────────────────────────────────

# --sample_n_groups N (scripts 5, 6, 11, 12.2, 13.2)
if args.sample_n_groups:
    df = sample_n_groups(df, n=args.sample_n_groups,
                          min_per_group=args.min_per_group, seed=args.seed)

# --pairs MM|FF|MF (scripts 15, 16, 17)
if args.pairs:
    df = sample_pairs(df, pair_type=args.pairs,
                       n_groups=args.n_groups, seed=args.seed)

# ──────────────────────────────────────────────
#  Auto-name output
# ──────────────────────────────────────────────

if args.output:
    meta_out = args.output
else:
    parts = []
    if args.group_col:
        parts.append(args.group_col.replace(" ", "_"))
    if args.dedup_col:
        parts.append(f"dedup_{args.dedup_col}".replace(" ", "_"))
    if args.exclude_immune:
        parts.append("noImmune")
    if args.germ_layer:
        parts.append(args.germ_layer.lower())
    for c, v in zip(args.keep_col, args.keep_val):
        parts.append(f"keep_{c}_{v}".replace(" ", "_"))
    if args.sex_filter:
        parts.append(f"{args.sex_filter}Only")
    if args.mixed_sex:
        parts.append("mixedSex")
    if args.keep_groups_col:
        parts.append(f"{args.keep_groups_col}_{args.keep_groups_val}".replace(" ", "_"))
    if args.sample_n_groups:
        parts.append(f"sample{args.sample_n_groups}gp")
    if args.pairs:
        parts.append(f"pairs{args.pairs}")
    suffix   = "_".join(parts) if parts else "filtered"
    meta_out = os.path.splitext(args.meta)[0] + f".{suffix}.csv"

df.to_csv(meta_out, index=False)
print(f"\nSaved {len(df):,} rows -> {meta_out}")
print(f"  Final groups: {df['Analysis group'].nunique()}")
