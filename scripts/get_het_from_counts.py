#!/usr/bin/env python3
"""
Usage: ./get_het_from_counts.py <countsfile>

Counts file is sorted by chrom, pos and has:
chrom pos allele count [freq]
"""
import numpy as np
import sys

try:
    countsfile = sys.argv[1]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

def OutputHet(chrom, pos, allele_counts):
    counts = np.array([item for item in allele_counts.values()])
    freqs = counts/np.sum(counts)
    het = 1-sum([item**2 for item in freqs])
    sys.stdout.write("\t".join([chrom, str(pos), str(het), str(np.sum(counts))])+"\n")

current_chrom = None
current_pos = None
allele_counts = {}
with open(countsfile, "r") as f:
    for line in f:
        if "pos" in line: continue # header
        chrom, pos = line.strip().split()[0:2]
        pos = int(pos)
        allele, count = line.strip().split()[2:4]
        allele = int(allele)
        count = int(count)
        if chrom == current_chrom and pos == current_pos:
            allele_counts[allele] = allele_counts.get(allele, 0) + count
        else:
            if current_chrom is not None: OutputHet(current_chrom, current_pos, allele_counts)
            current_chrom = chrom
            current_pos = pos
            allele_counts = {allele: count}
