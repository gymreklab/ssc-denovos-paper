#!/bin/bash

set -e
OUTDIR=/storage/mgymrek/ssc-denovos/allele-freqs/
mkdir -p ${OUTDIR}/intermediate

# Compute allele counts separately per phase
PARENTS=/storage/mgymrek/workspace/ssc-denovos/metadata/all_parent_ids.txt
for phase in 1 2 3_1 3_2
do
    #PARENTS=/storage/mgymrek/workspace/ssc-denovos/metadata/all_parent_sample_IDs_phase${phase}.txt
    for chrom in $(seq 1 22)
    do
	# Output chrom, pos, allele, count for each phase
	VCF=/storage/ileena/ssc-gangstr-denovos/vcf/merged/phase${phase}/phase${phase}_${chrom}.filtered.vcf.gz
	cmd="./compute_allele_counts.py $VCF $PARENTS > ${OUTDIR}/intermediate/phase${phase}_chr${chrom}_allele_counts.tab"
	echo $cmd
    done
done #| xargs -I% -n1 -P10 sh -c "%"

# TODO - this part. redo with phase 3_1 added
# And chrom, pos, heterozygosity, count
for chrom in $(seq 1 22)
do
    echo $chrom
    FILE1=${OUTDIR}/intermediate/phase1_chr${chrom}_allele_counts.tab
    FILE2=${OUTDIR}/intermediate/phase2_chr${chrom}_allele_counts.tab
#    FILE3=${OUTDIR}/intermediate/phase3_1_chr${chrom}_allele_counts.tab
    FILE4=${OUTDIR}/intermediate/phase3_2_chr${chrom}_allele_counts.tab
    # Convert to chrom, pos, allele, count, freq by merging across phases
#    ./merge_allele_counts.py $FILE1 $FILE2 $FILE3 $FILE4 > ${OUTDIR}/intermediate/merged_chr${chrom}_allele_counts.tab
    # Get het for each locus
    ./get_het_from_counts.py ${OUTDIR}/intermediate/merged_chr${chrom}_allele_counts.tab > ${OUTDIR}/intermediate/merged_chr${chrom}_het.tab
done

# Combine by chrom
#echo "chrom,pos,allele,count" | sed 's/,/\t/g' > ${OUTDIR}/SSC_merged_allele_counts.tab
#cat ${OUTDIR}/intermediate/merged_chr*_allele_counts.tab | grep -v chrom >> ${OUTDIR}/SSC_merged_allele_counts.tab
echo "chrom,pos,het,total" | sed 's/,/\t/g' > ${OUTDIR}/SSC_merged_het.tab
cat ${OUTDIR}/intermediate/merged_chr*_het.tab >> ${OUTDIR}/SSC_merged_het.tab
