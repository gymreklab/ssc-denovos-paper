#!/usr/bin/env python3
"""
Usage: ./merge_allele_counts.py [files]
"""

import sys
import pandas as pd
import numpy as np

files = sys.argv[1:]

# Read in each dataframe for each chrom
dflist = []
for f in files:
    df = pd.read_csv(f, sep="\t", names=["chrom","pos","allele","count"])
    dflist.append(df)

# Concatenate all the dfs
data = pd.concat(dflist, ignore_index=True)

# counts
counts = data.groupby(["chrom","pos","allele"], as_index=False).agg({"count": np.sum})

# output
counts.to_csv(sys.stdout, sep="\t", index=False)
