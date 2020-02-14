#!/bin/bash

# GangSTR params:
#The output path for the STR genotypes
GANGSTRDIR=s3://ssc-gangstr-denovos/vcf/phase4
PHASE=phase4

# DumpSTR params:
DUMPSTRDIR=s3://ssc-gangstr-denovos/vcf/filtered2/phase4

# MergeSTR params:
MERGEDIR=s3://ssc-gangstr-denovos/vcf/merged/${PHASE}

# CookieMonSTR params:
DENOVOVCF=s3://ssc-gangstr-denovos/vcf/filtered_PASS_only_Jan20/filtered_PASS_only_Jan20/${PHASE}
DENOVODIR=s3://ssc-gangstr-denovos/denovos_mean_priors_Nov19/${PHASE}
DENOVOMODELFILTERDIR=s3://ssc-gangstr-denovos/denovos_GW_priors_Jan20/${PHASE}
DENOVONAIVEDIR=s3://ssc-gangstr-denovos/denovos_naive_Jan20/${PHASE}
