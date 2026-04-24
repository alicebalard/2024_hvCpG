#!/usr/bin/env python3
"""
prepare_metadata.py — Flexible metadata preparation for CpG methylation pipelines
==================================================================================
Adds an 'Analysis group' column and applies optional filtering rules.

Supported filters (all combinable):
  --exclude_immune          Drop rows where 'Immune?' == True
  --sex_filter M|F          Keep only groups where ALL samples have given sex
  --keep_col / --keep_val   Generic keep: rows where column == value (repeatable)
  --drop_col / --drop_val   Generic drop: rows where column == value (repeatable)
  --keep_groups_col / --keep_groups_val
                            Keep groups where ALL rows have column == value

Usage examples
--------------
# All groups:
python prepare_metadata.py --meta SupTab1.csv --output meta_all.csv

# No immune cells:
python prepare_metadata.py --meta SupTab1.csv --output meta_noImmune.csv \\
    --exclude_immune

# Male-only groups:
python prepare_metadata.py --meta SupTab1.csv --output meta_maleOnly.csv \\
    --sex_filter M

# CD4+/CD8+ T cells only:
python prepare_metadata.py --meta SupTab1.csv --output meta_TcellSubsets.csv \\
    --keep_col "Cell type" --keep_val "CD4+ T cells" \\
    --keep_col "Cell type" --keep_val "CD8+ T cells"

# Drop blood tissue:
python prepare_metadata.py --meta SupTab1.csv \\
    --drop_col "Source Tissue" --drop_val "Blood"

Author: Alice Balard
"""

import argparse
import pandas as pd
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cpg_utils import add_analysis_group, filter_by_column_value, filter_groups_by_column

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

# Column name overrides
parser.add_argument("--tissue_col", default="Source Tissue",
                    help="Column for tissue (default: 'Source Tissue').")
parser.add_argument("--cell_col",   default="Cell type",
                    help="Column for cell type (default: 'Cell type').")
parser.add_argument("--sample_col", default="Sample name",
                    help="Column for sample name (default: 'Sample name').")

# Filtering
parser.add_argument("--exclude_immune",   action="store_true",
                    help="Drop rows where 'Immune?' == True.")
parser.add_argument("--sex_filter",       default=None, metavar="SEX",
                    help="Keep groups where ALL samples have this sex (M or F).")
parser.add_argument("--keep_col", action="append", default=[], metavar="COL",
                    help="Column to filter (use with --keep_val). Repeatable.")
parser.add_argument("--keep_val", action="append", default=[], metavar="VAL",
                    help="Value to keep in --keep_col column. Repeatable.")
parser.add_argument("--drop_col", action="append", default=[], metavar="COL",
                    help="Column to filter (use with --drop_val). Repeatable.")
parser.add_argument("--drop_val", action="append", default=[], metavar="VAL",
                    help="Value to drop in --drop_col column. Repeatable.")
parser.add_argument("--keep_groups_col", default=None, metavar="COL",
                    help="Keep groups where ALL rows have this column…")
parser.add_argument("--keep_groups_val", default=None, metavar="VAL",
                    help="…equal to this value (pair with --keep_groups_col).")

args = parser.parse_args()

# Validate paired args
if len(args.keep_col) != len(args.keep_val):
    parser.error("--keep_col and --keep_val must be used in pairs.")
if len(args.drop_col) != len(args.drop_val):
    parser.error("--drop_col and --drop_val must be used in pairs.")
if bool(args.keep_groups_col) != bool(args.keep_groups_val):
    parser.error("--keep_groups_col and --keep_groups_val must be used together.")

# ──────────────────────────────────────────────
#  Load and process
# ──────────────────────────────────────────────

print(f"\n📂 Loading: {args.meta}")
df = pd.read_csv(args.meta)
print(f"   {len(df):,} rows, {df.shape[1]} columns.")

# 1. Add Analysis group
df = add_analysis_group(df, tissue_col=args.tissue_col, cell_col=args.cell_col)
print(f"   Added 'Analysis group' ({df['Analysis group'].nunique()} unique groups).")

# 2. --exclude_immune
if args.exclude_immune:
    if "Immune?" not in df.columns:
        print("  ⚠️  --exclude_immune: 'Immune?' column not found — skipping.")
    else:
        before = len(df)
        df = df[df["Immune?"] != True]
        print(f"  --exclude_immune: removed {before - len(df):,} rows → {len(df):,} remain.")

# 3. --keep_col / --keep_val  (OR within same column, AND across different columns)
for col, val in zip(args.keep_col, args.keep_val):
    # Try coercing value to match column dtype
    if col in df.columns:
        try:
            val = df[col].dtype.type(val)
        except (ValueError, AttributeError):
            pass
    df = filter_by_column_value(df, col, val, keep=True)

# 4. --drop_col / --drop_val
for col, val in zip(args.drop_col, args.drop_val):
    if col in df.columns:
        try:
            val = df[col].dtype.type(val)
        except (ValueError, AttributeError):
            pass
    df = filter_by_column_value(df, col, val, keep=False)

# 5. --sex_filter
if args.sex_filter:
    df = filter_groups_by_column(df, filter_col="sex", required_value=args.sex_filter)

# 6. --keep_groups_col / --keep_groups_val
if args.keep_groups_col:
    df = filter_groups_by_column(df,
                                 filter_col=args.keep_groups_col,
                                 required_value=args.keep_groups_val)

# ──────────────────────────────────────────────
#  Auto-name output if not given
# ──────────────────────────────────────────────

if args.output:
    meta_out = args.output
else:
    parts = []
    if args.exclude_immune:        parts.append("noImmune")
    if args.sex_filter:            parts.append(f"{args.sex_filter}only")
    if args.keep_groups_col:       parts.append(f"{args.keep_groups_col}_{args.keep_groups_val}")
    for c, v in zip(args.keep_col, args.keep_val):
        parts.append(f"keep_{c}_{v}".replace(" ", "_"))
    suffix   = "_".join(parts) if parts else "filtered"
    meta_out = os.path.splitext(args.meta)[0] + f".{suffix}.csv"

df.to_csv(meta_out, index=False)
print(f"\n✅ Saved {len(df):,} rows → {meta_out}")
print(f"   Final groups: {df['Analysis group'].nunique()}")
