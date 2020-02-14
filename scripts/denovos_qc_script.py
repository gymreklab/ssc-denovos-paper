#!/usr/bin/env python3

"""
QC filtering script from raw de novo output -> final list of mutations

"""

import argparse
import numpy as np
import pandas as pd
import scipy.stats
import sys

def load_data(fn):
    df = pd.read_csv(fn, sep="\t", header=0)
    df["chrom"]  = df["chrom"].map(lambda x: x.lstrip("chr"))
    df.family = df.family.astype(int)
    df.chrom = df.chrom.astype(int)
    df.pos = df.pos.astype(int)
    df.sort_values(by=["family", "chrom", "pos"], inplace=True)
    return df

# Further het dropout filtering
def FilterHetDropout(x, require_encl):
    if abs(x["mutsize"])<5: return False # don't worry about short ones
    # For longer ones, require more enclosing reads
    # If both het, we're good
    if not x["pat_hom"] and not x["mat_hom"]: return False
    # If both hom, or POO is hom, require new threshold encl
    if (x["pat_hom"] and x["mat_hom"]) or \
        (x["poocase"]==2 and x["pat_hom"]) or \
        (x["poocase"]==3 and x["mat_hom"]):
        return (x["encl_child"] < require_encl)
    # If neither of those cases is true but one of parents is homozygous
    if (x["pat_hom"]) or (x["mat_hom"]):
        return (x["encl_child"] < require_encl)
    return False # most of these are unknown POO. look ok

# Fix but in getting mutation size. 
# Returns tuple of updated mutsize,new allele
def FixMutationSize(x):
    if x["poocase"] not in [2,3]: return (x["mutsize"], x["newallele"])
    # Mutation from father. If both child alleles in mother, reassess
    child_gt = x["child_gt"].split(",")
    parent2_gt = None # novel allele
    if x["poocase"] == 2:
        if x["child_gt"] == x["mat_gt"]:
            parent2_gt = x["pat_gt"].split(",")
        else: return (x["mutsize"], x["newallele"])
    # Mutation from mother. If both child alleles in father, reassess
    if x["poocase"] == 3:
        if x["child_gt"] == x["pat_gt"]:
            parent2_gt = x["mat_gt"].split(",")
        else: return (x["mutsize"], x["newallele"])
    # Check all possible mut sizes
    mindiff = 10000
    newallele = 0
    for possible_newallele in child_gt:
        for parent_allele in parent2_gt:
            if abs(int(possible_newallele)-int(parent_allele)) < abs(mindiff):
                mindiff = int(possible_newallele)-int(parent_allele)
                newallele = int(possible_newallele)
    return (mindiff, newallele)

# Test if biased toward insertions or deletions
def GetDirectionPval(x):
    n = len(x)
    k = len([item for item in x if item<0])
    return scipy.stats.binom_test(k, n=n, p=0.5, alternative="two-sided")

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--mutfile", help="File with mutations passing posterior threshold e.g. ....filtered_mutations.tab.gz", type=str, required=True)
    parser.add_argument("--locussummary", help="File with locus summary e.g. ...locus_summary.tab.gz", type=str, required=False) # TO ADD
    parser.add_argument("--removefams", help="Text file with families to remove e.g. due to sample contamination, high mutation rates, etc.", type=str, required=False, default="/storage/ileena/denovos5/metadata/ssc_remove_families.txt")
    parser.add_argument("--outmutfile", help="Output file with mutations passing QC filters.", type=str, required=True)
    parser.add_argument("--filter-outlier-child", type=int, required=False, dest="CHILD_THRESH", default=5, help="Use to exclude children with # mutations above standard deviations from mean (Default 5 SD).")
    parser.add_argument("--filter-outlier-loc-sd", type=int, required=False, dest="LOCI_SD_THRESH",  help="Use to exclude TR loci with # mutations above standard deviations from mean.")
    parser.add_argument("--filter-outlier-loc-num", type=int, required=False, dest="LOCI_NUM_THRESH",  help="Use to exclude TR loci with # mutations above number threshold.")
    parser.add_argument("--filter-homozygous-calls", action='store_true', required=False, dest="HCALL",  help="Use to de novo calls where the child is homozygous for the de novo allele.")
    parser.add_argument("--filter-direction-bias", type=float, required=False, dest="LOCI_DIR_P", help="Filter calls with a ins/del bias with p< this threshold.")
    parser.add_argument("--fix-mutsize", action="store_true", help="Check that mutation size and new allele are set correctly.")
    parser.add_argument("--filter-hetdropout", type=int, required=False, dest="LOCI_HETDROPOUT", help="For long new alleles de novo from homozygous parents, require this many enclosing reads in the child")

    args = parser.parse_args()

    ### Load mutations list
    print("############### Loading mutations list #########")
    denovos = load_data(args.mutfile)
    outcols = denovos.columns
    print("Number of denovos in list", len(denovos))

    ### Remove families
    print("############### (FAMILY FILTER) Removing contaminated families #########")
    rmfams = pd.read_csv(args.removefams, sep="\t", header=None, names=["family"])
    rmfams.family  = rmfams.family.astype(int)
    denovos = denovos[~denovos.family.isin(rmfams.family.tolist())]
    print(len(rmfams.family), "families removed.")
    print("Number of denovos remaining in list", len(denovos))

    ### Remove trios
    print("############### (FAMILY FILTER) Removing trios #########")
    rmtrios = denovos[["family", "child"]].drop_duplicates().groupby(["family"]).size().reset_index()
    rmtrios = rmtrios[rmtrios[0]<2].family.tolist()
    denovos = denovos[~denovos.family.isin(rmtrios)]
    print(len(rmtrios), "trio families removed.")
    print("Number of denovos remaining in list", len(denovos))

    ### Remove outlier families
    print("############### (FAMILY FILTER) Removing outlier families #########")
    count = denovos.groupby(["family", "child"]).size().reset_index()
    if args.CHILD_THRESH > 0: child_threshold = count[0].mean() + args.CHILD_THRESH*count[0].std()
    else: child_threshold=None
    if child_threshold != None:
        rmoutlierfams = count[count[0]>child_threshold].family.tolist()
        denovos = denovos[~denovos.family.isin(rmoutlierfams)]
        print("{} outlier families removed with # muts >= {} SD + mean ({})".format(len(rmoutlierfams), args.CHILD_THRESH, child_threshold))
        print(rmoutlierfams)
        print("Number of denovos remaining in list", len(denovos))

    ### Remove call if both sibling have mutation
    print("############### (CALL FILTER) Remove call if both sibilings have a mutation #########")
    count = denovos.groupby(["chrom", "pos", "family"]).child.count().reset_index()
    remove2 = count[count.child>1].copy()
    remove2["both_sibs"] = True
    denovos = pd.merge(denovos, remove2[["chrom", "pos", "family", "both_sibs"]], on=["chrom", "pos", "family"], how="left")
    denovos = denovos[denovos.both_sibs!=True]
    print(remove2.shape[0], "calls with both sibling mutations removed.")
    print("Number of denovos remaining in list", len(denovos))

    ### Remove call if homozygous in child
    print("############### (CALL FILTER) Remove calls homozygous for new allele in child #########")
    if args.HCALL:
        denovos["child_hcall"] = denovos.child_gt.apply(lambda x: len(set(x.split(",")))==1)
        denovos = denovos[~denovos.child_hcall]
        print(len(denovos.child_hcall), "homozygous denovo calls in child removed.")
        print("Number of denovos remaining in list", len(denovos))

    #### Fix mutation size and new allele
    print("############### (CALL BUG FIX) Fixing mutation size and new allele #########")
    if args.fix_mutsize:
        mutdata = denovos.apply(FixMutationSize, 1)
        denovos["mutsize"] = [item[0] for item in mutdata]
        denovos["newallele"] = [item[1] for item in mutdata]

    #### Apply further filtering on enclosing reads to avoid het dropout
    print("############### (CALL FILTER) Stricter filtering for possible het dropouts #########")
    denovos["mat_hom"] = denovos["mat_gt"].apply(lambda x: len(set(x.split(",")))==1)
    denovos["pat_hom"] = denovos["pat_gt"].apply(lambda x: len(set(x.split(",")))==1)
    if args.LOCI_HETDROPOUT is not None:
        denovos["possible.het.do"] = denovos.apply(lambda x: FilterHetDropout(x, args.LOCI_HETDROPOUT), 1)
        print("{} denovos removed due to possible het dropout and child encl < {}.".format(np.sum(denovos["possible.het.do"]), args.LOCI_HETDROPOUT))
        denovos = denovos[~denovos["possible.het.do"]]
        print("Number of denovos remaining in list", len(denovos))

    #### Remove loci with a strong bias toward insertions vs. deletions
    print("############### (LOCUS FILTER) Remove loci with ins/del bias #########")
    if args.LOCI_DIR_P is not None:
        dirp = denovos.groupby(["chrom","pos"], as_index=False).agg({"mutsize": GetDirectionPval})
        dirp.columns = ["chrom","pos","dir.P"]
        denovos = pd.merge(denovos, dirp, on=["chrom","pos"])
        denovos = denovos[denovos["dir.P"]>args.LOCI_DIR_P]
        print("{} loci with direction bias p < {} removed.".format(dirp[dirp["dir.P"]<=args.LOCI_DIR_P].shape[0], args.LOCI_DIR_P))
        print("Number of denovos remaining in list", len(denovos))

    #### Remove outlier loci
    print("############### (LOCUS FILTER) Remove outlier loci #########")
    count = denovos.groupby(["chrom", "pos"]).child.count().reset_index()
    if args.LOCI_SD_THRESH != None: tr_threshold = count.child.mean() + args.LOCI_SD_THRESH*count.child.std()
    elif args.LOCI_NUM_THRESH != None: tr_threshold = args.LOCI_NUM_THRESH
    else: tr_threshold = None
    if tr_threshold != None:
        remove1 = count[count.child>tr_threshold].copy()
        remove1["extreme_loci"] = True
        denovos = pd.merge(denovos, remove1[["chrom", "pos", "extreme_loci"]], on=["chrom", "pos"], how="left")
        denovos = denovos[denovos.extreme_loci!=True]
        print("{} loci with > {} children with mutations removed.".format(remove1.shape[0], tr_threshold))
        print("Number of denovos remaining in list", len(denovos))

    ## Write output final list of de novos
    print("############### Output final de novos list #########")
    denovos[outcols].to_csv(args.outmutfile, sep="\t", index=False, header=True)

if __name__ == "__main__":
    main()
