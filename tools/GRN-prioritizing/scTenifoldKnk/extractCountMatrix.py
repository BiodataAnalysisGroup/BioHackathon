#!/usr/bin/env python

import argparse
from pathlib import Path
import scanpy as sc
import pandas as pd

#
# Parse the arguments
#
parser = argparse.ArgumentParser(description = "Load AnnData object, extract the count matrix and save it as a CSV.")

# Single-cell data input argument
parser.add_argument(
    "-sc", 
    "--sc_data", 
    help = "File path to the single-cell data as AnnData object.", 
    required = True
)

# Single-cell count matrix CSV output file path
parser.add_argument(
    "-o", 
    "--out_path", 
    help = "Single-cell count matrix CSV output file path.", 
    required = True
)

# Parse all arguments
args = parser.parse_args()

#
# Read and write the data
# 
print("Reading single-cell data...")
ann_data = sc.read_h5ad(Path(args.sc_data))

print("Writing single-cell count matrix into a CSV...")
count_matrix = ann_data.raw.X.toarray()
pd.DataFrame(data = count_matrix, 
             index = count_matrix.obs_names, 
             columns = count_matrix.raw.var_names).to_csv(Path(args.out_path))