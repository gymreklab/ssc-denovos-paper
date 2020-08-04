#!/usr/bin/env python3

import gzip
import numpy as np
import pandas as pd
import scipy.stats
import sys

ALLMUTFILE = "/storage/ileena/ssc-gangstr-denovos/denovos_GW_priors_Jan20/combined/SSC_allphases_011720_denovos_GW_priors_Jan20.sorted.all_mutations.tab.gz"
FILTMUTFILE = "/storage/ileena/ssc-gangstr-denovos/denovos_GW_priors_Jan20/combined/SSC_allphases_011720_denovos_GW_priors_Jan20.final_qc_mutations.tab"

# Force all chroms to be ints
def GetChrom(x):
    if "chr" in str(x): return int(x[3:])
    else: return int(x)

# Keep track of column indices
COL_PT = 6
COL_POO = 11
COL_POO_MAT = 23
COL_POO_PAT = 24
COL_CHILD_GT = 16
COL_MAT_GT = 17
COL_PAT_GT = 18

# Read in filtered mutations
filtmut = pd.read_csv(FILTMUTFILE, sep="\t")
filtmut["chrom"] = filtmut["chrom"].apply(GetChrom)

# Filter problematic families
rmfams = [14151, 12434, 12281, 13673, 13351, 13355, 13143]
filtmut = filtmut[~filtmut["family"].isin(rmfams)]
usefams = set(filtmut["family"])

# Write header
sys.stdout.write("\t".join(["chrom","pos", \
                            "all_unaff","all_aff","mut_unaff","mut_aff", \
                            "newalleles_unaff","newalleles_aff", \
                            "fisher_p","fisher_odds", \
                            "assoc_p", \
                            "tdt_counts_unaff_mother","tdt_counts_unaff_father","tdt_counts_unaff_combined",
                            "tdt_counts_aff_mother","tdt_counts_aff_father","tdt_counts_aff_combined",
                            "tdt_p_mother","tdt_odds_mother", \
                            "tdt_p_father","tdt_odds_father",
                            "tdt_p_combined","tdt_odds_combined"])+"\n")

def TestTDT(tdt_unaff, tdt_aff):
    tdt_table = [[tdt_unaff[1], tdt_aff[1]], \
                 [tdt_unaff[2], tdt_aff[2]]]
    return scipy.stats.fisher_exact(tdt_table)

def ComputeLocusStats(proc_lines, filtmut):
    if len(proc_lines) == 0: return

    # Figure out chrom and pos. First check them
    chroms = set([GetChrom(l.split()[0]) for l in proc_lines])
    poss = set([l.split()[1] for l in proc_lines])
    if len(chroms) > 1:
        sys.stderr.write("Error. found more than one chrom in proc_lines %s"%chroms)
        sys.exit(1)
    if len(poss) > 1:
        sys.stderr.write("Error. found more than one pos in proc_lines %s"%chroms)
        sys.exit(1)
    chrom = GetChrom(proc_lines[0].split()[0])
    pos = int(proc_lines[0].split()[1])

    # Subset mutations file based on this
    muts = filtmut[(filtmut["chrom"]==chrom) & (filtmut["pos"]==pos)]

    # initialize values
    num_calls_aff = 0
    num_calls_unaff = 0
    num_mut_aff = muts[muts["phenotype"]==2].shape[0]
    num_mut_unaff = muts[muts["phenotype"]==1].shape[0]
    if num_mut_aff > 0:
        new_alleles_aff = ",".join(list(muts[muts["phenotype"]==2].sort_values("newallele")["newallele"].apply(str)))
    else:
        new_alleles_aff = "."
    if num_mut_unaff > 0:
        new_alleles_unaff = ",".join(list(muts[muts["phenotype"]==1].sort_values("newallele")["newallele"].apply(str)))
    else:
        new_alleles_unaff = "."
    tdt_unaff_mother = [0,0,0]
    tdt_aff_mother = [0,0,0]
    tdt_unaff_father = [0,0,0]
    tdt_aff_father = [0,0,0]
    alleles_unaff = []
    alleles_aff = []

    # Go through lines and add to values
    for line in proc_lines:
        items = line.strip().split()
        # Update numcall counts
        if int(items[COL_PT]) == 1:
            num_calls_unaff += 1
        else:
            num_calls_aff += 1
        if int(items[COL_POO]) == 1:
            # Update TDT and allele counts
            tdt_mother = int(items[COL_POO_MAT])
            tdt_father = int(items[COL_POO_PAT])
            child_gt = sum([int(item) for item in items[COL_CHILD_GT].split(",")])
            # Fix TDT if all have same genotype (TODO bug in cookiemonstr)
            if (items[COL_CHILD_GT] == items[COL_MAT_GT]) and (items[COL_CHILD_GT] == items[COL_PAT_GT]):
                tdt_mother = 0
                tdt_father = 0
            if int(items[COL_PT]) == 1:
                tdt_unaff_mother[tdt_mother] += 1
                tdt_unaff_father[tdt_father] += 1
                alleles_unaff.append(child_gt)
            else:
                tdt_aff_mother[tdt_mother] += 1
                tdt_aff_father[tdt_mother] += 1
                alleles_aff.append(child_gt)

    # Compute locus-level mutation burden
    fisher_table = [[num_mut_unaff, num_mut_aff], \
                    [num_calls_unaff-num_mut_unaff, num_calls_aff-num_mut_aff]]
    fisher_odds, fisher_p = scipy.stats.fisher_exact(fisher_table, alternative="greater")

    # Compute association p-values
    try:
        assoc_p = scipy.stats.mannwhitneyu(alleles_unaff, alleles_aff, alternative="two-sided")[1]
    except ValueError:
        assoc_p = np.nan

    # Compute TDT stats
    tdt_unaff_combined = [tdt_unaff_mother[i]+tdt_unaff_father[i] for i in range(3)]
    tdt_aff_combined = [tdt_aff_mother[i]+tdt_aff_father[i] for i in range(3)]
    tdt_mother_odds, tdt_mother_p = TestTDT(tdt_unaff_mother, tdt_aff_mother)
    tdt_father_odds, tdt_father_p = TestTDT(tdt_unaff_father, tdt_aff_father)
    tdt_combined_odds, tdt_combined_p = TestTDT(tdt_unaff_combined, tdt_aff_combined)

    # Output results
    output_vals = [chrom, pos, \
                   num_calls_unaff, num_calls_aff, num_mut_unaff, num_mut_aff, \
                   new_alleles_unaff, new_alleles_aff, \
                   fisher_p, fisher_odds, \
                   assoc_p, \
                   ":".join([str(item) for item in tdt_unaff_mother]), \
                   ":".join([str(item) for item in tdt_unaff_father]), \
                   ":".join([str(item) for item in tdt_unaff_combined]), \
                   ":".join([str(item) for item in tdt_aff_mother]), \
                   ":".join([str(item) for item in tdt_aff_father]), \
                   ":".join([str(item) for item in tdt_aff_combined]), \
                   tdt_mother_p, tdt_mother_odds, \
                   tdt_father_p, tdt_father_odds, \
                   tdt_combined_p, tdt_combined_odds]

    sys.stdout.write("\t".join([str(item) for item in output_vals])+"\n")

proc_lines = [] # keep track of list of lines to process for the locus
current_pos = -1
with gzip.open(ALLMUTFILE, "rt") as f:
    for line in f:
        if line.startswith("chrom"): continue # header
        if int(line.split()[4]) not in usefams: continue # filtered fam
        pos = int(line.split()[1])
        if pos == current_pos:
            proc_lines.append(line)
        else:
            if current_pos != -1:
                ComputeLocusStats(proc_lines, filtmut)
            proc_lines = [line]
            current_pos = pos
