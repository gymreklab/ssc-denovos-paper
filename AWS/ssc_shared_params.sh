#!/bin/bash

# Extra params
THREADS=4

# GangSTR params:
#Reference set of regions to genotype - Yes "chr" in hg38
REFBED=s3://gangstr/hg38/hg38_ver16.bed.gz
REFFASTA=s3://sscwgs-hg38/documentation/GRCh38_full_analysis_set_plus_decoy_hla.fa

# DumpSTR params:
FILTERREGIONS=s3://gangstr-refs/dumpstr/GRCh38GenomicSuperDup.sorted.bed.gz
MINDP=20
MAXDP=1000
HWEP=0.00001
MINCALLRATE=0.8
READLEN=150
NREADS=2
MINSAMPLECALLRATE=0.85

# CookieMonSTR - shared params:
AWSBUCKET=ssc-gangstr-denovos
FAMFILE=s3://ssc-gangstr-denovos/other/ssc_4phases_ids.ped
MAXALLELES=100
MINCOV=10
MINSCORE=0.9
MINSPANCOV=10
MINSUPPREADS=2
PTHRESH=0.8
MUTMODEL=s3://ssc-gangstr-denovos/denovos_mean_priors/mean_str_mutrates_hg38_ver16.bed

# CookieMonSTR - NAIVE params AND MODEL params:
MODELPTHRESH=0.5
DEFAULTPRIOR=-3
MINENCCHILD=3
MINENCPERPARENT=0.05
MINENCMATCH=0.9
MINTOTALENC=10

# Other params:
OTHERDIR=s3://ssc-gangstr-denovos/other
