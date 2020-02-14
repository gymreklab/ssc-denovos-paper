#!/bin/bash

# GangSTR params:
#The output path for the STR genotypes
PHASE=phase3_2
GANGSTRDIR=s3://ssc-gangstr-denovos/vcf/${PHASE}

# DumpSTR params:
DUMPSTRDIR=s3://ssc-gangstr-denovos/vcf/filtered2/${PHASE}

# MergeSTR params:
MERGEDIR=s3://ssc-gangstr-denovos/vcf/merged/${PHASE}

# CookieMonSTR params:
DENOVOVCF=s3://ssc-gangstr-denovos/vcf/filtered_PASS_only_Jan20/filtered_PASS_only_Jan20/${PHASE}
DENOVODIR=s3://ssc-gangstr-denovos/denovos_mean_priors_Nov19/${PHASE}
DENOVOMODELFILTERDIR=s3://ssc-gangstr-denovos/denovos_GW_priors_Jan20/${PHASE}
DENOVONAIVEDIR=s3://ssc-gangstr-denovos/denovos_naive_Jan20/${PHASE}
