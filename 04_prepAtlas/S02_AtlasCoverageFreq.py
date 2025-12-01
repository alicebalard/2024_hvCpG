import numpy as np
import glob
import os
from collections import Counter

# Directory with your beta files
beta_dir = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/betaFiles/"
beta_files = glob.glob(os.path.join(beta_dir, "*.beta"))

# Global counter for coverage frequencies
coverage_counter = Counter()

for fn in beta_files:
    content = np.fromfile(fn, dtype=np.uint8).reshape((-1, 2))
    coverage = content[:, 1]  # second column = coverage
    # Update counts efficiently
    unique, counts = np.unique(coverage, return_counts=True)
    coverage_counter.update(dict(zip(unique, counts)))

# Convert to numpy arrays for saving/plotting
coverages = np.array(sorted(coverage_counter.keys()), dtype=np.int32)
frequencies = np.array([coverage_counter[c] for c in coverages], dtype=np.int64)

# Save as tab-separated file for R
out_file = "/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/04_prepAtlas/coverage_histogram.tsv"
np.savetxt(out_file, np.column_stack([coverages, frequencies]),
           fmt="%d\t%d", header="coverage\tfrequency", comments="")

print(f"✅ Coverage histogram written to {out_file}")
