#!/usr/bin/env python

"""
Usage: ./compute_allele_counts.py <vcf> <parentids>
"""

import vcf
import sys

try:
    vcffile = sys.argv[1]
    parentfile = sys.argv[2]
except:
    sys.stderr(__doc__)
    sys.exit(1)

parents = [line.strip() for line in open(parentfile, "r").readlines()]

def GetEnclReads(sample):
    enclreads = {}
    encl = sample["ENCLREADS"]
    if encl == "NULL": return {}
    for item in encl.split("|"):
        al, count = item.split(",")
        enclreads[int(al)] = int(count)
    return enclreads

def FilterCall(sample):
    enclreads = GetEnclReads(sample)
    # Check has at least 10 enclosing reads
    if sum(enclreads.values()) < 10: return True
    # Check each allele supported by >=3 reads
    a1, a2 = [int(item) for item in sample["REPCN"]]
    if enclreads.get(a1, 0) < 3: return True
    if enclreads.get(a2, 0) < 3: return True
    # Check >= 90% of enclosing reads match gt
    badreads = 0
    for a in enclreads.keys():
        if a != a1 and a != a2: badreads += enclreads[a]
    if badreads*1.0/sum(enclreads.values()) > 0.1: return True
    return False

reader = vcf.Reader(open(vcffile, "rb"))
for record in reader:
    allele_counts = {} # allele->count
    for sample in record:
        if not sample.called: continue
        if not sample.sample in parents: continue
        if FilterCall(sample):
            continue
        gt = sample["REPCN"]
        for a in gt:
            a = int(a)
            allele_counts[a] = allele_counts.get(a, 0) + 1
    for a in allele_counts:
        items = [record.CHROM, record.POS, a, allele_counts[a]]
        sys.stdout.write("\t".join([str(item) for item in items])+"\n")
