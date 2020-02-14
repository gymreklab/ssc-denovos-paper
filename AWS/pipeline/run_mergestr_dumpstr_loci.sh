#!/bin/bash

#Input values
CHROM=$1
SSC_PARAMS=${2}
echo "[run_mergestr.sh]: ${CHROM}"

### Set up
DATADIR=/scratch # This is where we have all the EBS storage space mounted
mkdir -p ${DATADIR}/${CHROM}/datafiles
mkdir -p ${DATADIR}/${CHROM}/results
mkdir -p ${DATADIR}/${CHROM}/tmp
rm ${DATADIR}/${CHROM}/results/*
rm ${DATADIR}/${CHROM}/datafiles/*
rm ${DATADIR}/${CHROM}/tmp/*

aws s3 cp s3://ssc-gangstr-denovos/scripts/ssc_shared_params.sh ${DATADIR}/${FAMID}/ssc_shared_params.sh || die "Error copying ssc_shared_params.sh"
source ${DATADIR}/${FAMID}/ssc_shared_params.sh
aws s3 cp s3://ssc-gangstr-denovos/scripts/${SSC_PARAMS} ${DATADIR}/${FAMID}/${SSC_PARAMS} || die "Error copying ${SSC_PARAMS}"
source ${DATADIR}/${FAMID}/${SSC_PARAMS}

### Download VCF files needed
# SSC GangSTR VCF files
aws s3 cp ${DUMPSTRDIR}/  ${DATADIR}/${CHROM}/datafiles/ --recursive  --exclude "*" --include "*_${CHROM}.filtered.vcf*"

### Download data files needed for DumpSTR
# SEGDUPS Filter file
aws s3 cp ${FILTERREGIONS} ${DATADIR}/${CHROM}/datafiles/filter.bed.gz || die "Error copying ${FILTERREGIONS}"
aws s3 cp ${FILTERREGIONS}.tbi ${DATADIR}/${CHROM}/datafiles/filter.bed.gz.tbi || die "Error copying ${FILTERREGIONS}"

### Remove VCFs for families not passing call rate
aws s3 cp ${DUMPSTRDIR}/${PHASE}_remove_families_callrate.txt  ${DATADIR}/${FAMID}/datafiles/${PHASE}_remove_families_callrate.txt
while read FAMID;
do
  rm ${DATADIR}/${CHROM}/datafiles/${FAMID}_*.filtered.vcf.gz
done <${DATADIR}/${FAMID}/datafiles/${PHASE}_remove_families_callrate.txt

### Run MergeSTR per chrom
MERGEFILE=${PHASE}_${CHROM}
rm ${DATADIR}/${CHROM}/datafiles/${MERGEFILE}.filtered.vcf.gz  #avoid duplicate
rm ${DATADIR}/${CHROM}/datafiles/${MERGEFILE}.sorted.vcf.gz #avoid duplicate
VCFLIST=$(ls ${DATADIR}/${CHROM}/datafiles/*_${CHROM}.filtered.vcf.gz | awk '{print $1};' | tr '\n' ',' | sed 's/,$//')

echo "[run_mergestr.sh]: Starting mergeSTR for ${VCFLIST}"
mergeSTR --vcfs ${VCFLIST} --out ${DATADIR}/${CHROM}/results/${MERGEFILE}
df -h ${DATADIR}/${CHROM}/ ## DEBUG
df -h ${DATADIR}/${CHROM}/tmp/  ## DEBUG
cat ${DATADIR}/${CHROM}/results/${MERGEFILE}.vcf | vcf-sort -t ${DATADIR}/${CHROM}/tmp/  -p 4 | bgzip -c > ${DATADIR}/${CHROM}/results/${MERGEFILE}.sorted.vcf.gz
tabix -p vcf ${DATADIR}/${CHROM}/results/${MERGEFILE}.sorted.vcf.gz
aws s3 cp ${DATADIR}/${CHROM}/results/${MERGEFILE}.sorted.vcf.gz ${MERGEDIR}/
aws s3 cp ${DATADIR}/${CHROM}/results/${MERGEFILE}.sorted.vcf.gz.tbi ${MERGEDIR}/

### Run DumpSTR - Locus filters
cmd="echo [run_dumpstr.sh]: Running DumpSTR ${MERGEFILE}.sorted.vcf.gz"
cmd="${cmd}; dumpSTR --vcf ${DATADIR}/${CHROM}/results/${MERGEFILE}.sorted.vcf.gz \
--filter-regions ${DATADIR}/${CHROM}/datafiles/filter.bed.gz \
--filter-regions-names SEGDUP \
--min-locus-hwep ${HWEP} \
--min-locus-callrate ${MINCALLRATE} \
--out ${DATADIR}/${CHROM}/results/${MERGEFILE}.filtered"
cmd="${cmd} && cat ${DATADIR}/${CHROM}/results/${MERGEFILE}.filtered.vcf | vcf-sort -t ${DATADIR}/${CHROM}/tmp/ | bgzip -c > ${DATADIR}/${CHROM}/results/${MERGEFILE}.filtered.vcf.gz"
cmd="${cmd} && tabix -p vcf ${DATADIR}/${CHROM}/results/${MERGEFILE}.filtered.vcf.gz"

### Upload merged VCF file to S3
cmd="${cmd} && aws s3 cp ${DATADIR}/${CHROM}/results/${MERGEFILE}.filtered.vcf.gz ${MERGEDIR}/"
cmd="${cmd} && aws s3 cp ${DATADIR}/${CHROM}/results/${MERGEFILE}.filtered.vcf.gz.tbi ${MERGEDIR}/"
cmd="${cmd} && aws s3 cp ${DATADIR}/${CHROM}/results/${MERGEFILE}.filtered.loclog.tab  ${MERGEDIR}/"
cmd="${cmd} && aws s3 cp ${DATADIR}/${CHROM}/results/${MERGEFILE}.filtered.samplog.tab ${MERGEDIR}/"

echo $cmd | xargs -n1 -I% -P${THREADS} sh -c "%"

# Get BED file of loci to remove
vcf-query ${DATADIR}/${CHROM}/results/${MERGEFILE}.filtered.vcf.gz -f '%CHROM\t%POS\t%INFO/END\t%FILTER\n' |  awk '{if ($4 != "PASS") print $0 ;}' >  ${DATADIR}/${CHROM}/results/${MERGEFILE}_merged_loci_remove.bed
aws s3 cp  ${DATADIR}/${CHROM}/results/${MERGEFILE}_merged_loci_remove.bed ${MERGEDIR}/

### Cleanup before moving on to next job
rm ${DATADIR}/${CHROM}/results/*
rm ${DATADIR}/${CHROM}/datafiles/*
rm ${DATADIR}/${CHROM}/tmp/*
exit 0
