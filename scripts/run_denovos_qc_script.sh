#!/bin/bash


# python3 -W ignore denovos_qc_script.py \
# --mutfile /storage/ileena/ssc-gangstr-denovos/denovos_GW_priors_Jan20/SSC_allphases_011720_denovos_GW_priors_Jan20.filtered_mutations.tab.gz \
# --outmutfile /storage/ileena/ssc-gangstr-denovos/denovos_GW_priors_Jan20/SSC_allphases_011720_denovos_GW_priors_Jan20.final_qc_mutations.tab \
# --filter-outlier-child 4 \
# --filter-outlier-loc-num 25 \
# --filter-homozygous-calls \
# --fix-mutsize \
# --filter-direction-bias 0.0001 \
# --filter-hetdropout 6

python3 -W ignore denovos_qc_script.py \
 --mutfile /storage/ileena/ssc-gangstr-denovos/denovos_naive_Jan20/SSC_allphases_011720_denovos_naive_Jan20.filtered_mutations.tab.gz \
 --outmutfile /storage/ileena/ssc-gangstr-denovos/denovos_naive_Jan20/SSC_allphases_011720_denovos_naive_Jan20.final_qc_mutations.tab \
 --filter-outlier-child 4 \
 --filter-outlier-loc-num 25 \
 --filter-homozygous-calls \
 --fix-mutsize \
 --filter-direction-bias 0.0001 \
 --filter-hetdropout 6
