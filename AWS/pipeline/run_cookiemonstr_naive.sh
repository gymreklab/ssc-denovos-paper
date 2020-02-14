#!/bin/bash

#CookieMonSTR run on filtered GangSTR VCF files on s3. Run after QC for all phases. Concats all CHRs and Runs per family.

#Input values
FAMID=$1
SSC_PARAMS=${2}
echo "[run_cookiemonstr.sh]: ${FAMID}"

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

# #############
DATADIR=/scratch # This is where we have all the EBS storage space mounted
# Set up
mkdir -p ${DATADIR}/${FAMID}/datafiles
mkdir -p ${DATADIR}/${FAMID}/results
mkdir -p ${DATADIR}/${FAMID}/tmp
rm ${DATADIR}/${FAMID}/results/*
rm ${DATADIR}/${FAMID}/datafiles/*
rm ${DATADIR}/${FAMID}/tmp/*

aws s3 cp s3://ssc-gangstr-denovos/scripts/ssc_shared_params.sh ${DATADIR}/${FAMID}/ssc_shared_params.sh || die "Error copying ssc_shared_params.sh"
source ${DATADIR}/${FAMID}/ssc_shared_params.sh
aws s3 cp s3://ssc-gangstr-denovos/scripts/${SSC_PARAMS} ${DATADIR}/${FAMID}/${SSC_PARAMS} || die "Error copying ${SSC_PARAMS}"
source ${DATADIR}/${FAMID}/${SSC_PARAMS}

### First, download data files needed for CookieMonSTR - .filtered.vcf.gz
# SSC filtered VCF files
aws s3 cp ${DUMPSTRDIR}/  ${DATADIR}/${FAMID}/datafiles/  --recursive  --exclude "*" --include "${FAMID}_*.filtered.vcf.gz*" || die "Error copying ${FAMID}_*.filtered.vcf.gz"
# Other filters
aws s3 cp  ${FAMFILE} ${DATADIR}/${FAMID}/datafiles/fam.ped



for chrom in {1..22}
do
  #Get filtered VCF file is created
  aws s3 cp ${DENOVOVCF}/${FAMID}_${chrom}.filtered.PASS_only.recode.vcf.gz ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.PASS_only.recode.vcf.gz
  aws s3 cp ${DENOVOVCF}/${FAMID}_${chrom}.filtered.PASS_only.recode.vcf.gz.tbi ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.PASS_only.recode.vcf.gz.tbi

  #If does not exist, then create filtered VCF
  if test ! -f "${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.PASS_only.recode.vcf.gz"; then

    # Get list of loci to remove
    aws s3 cp s3://ssc-gangstr-denovos/vcf/merged/${chrom}_merged_loci_remove.bed  ${DATADIR}/${FAMID}/results/${chrom}_merged_loci_remove.bed
    ##If does not exist, then create file
    if test ! -f "${DATADIR}/${FAMID}/results/${chrom}_merged_loci_remove.bed"; then
        ## Remove dumpSTR loci from all phases for all families
        aws s3 cp s3://ssc-gangstr-denovos/vcf/merged/phase1/phase1_${chrom}_merged_loci_remove.bed ${DATADIR}/${FAMID}/datafiles/
        aws s3 cp s3://ssc-gangstr-denovos/vcf/merged/phase2/phase2_${chrom}_merged_loci_remove.bed ${DATADIR}/${FAMID}/datafiles/
        aws s3 cp s3://ssc-gangstr-denovos/vcf/merged/phase3_1/phase3_1_${chrom}_merged_loci_remove.bed ${DATADIR}/${FAMID}/datafiles/
        aws s3 cp s3://ssc-gangstr-denovos/vcf/merged/phase3_2/phase3_2_${chrom}_merged_loci_remove.bed ${DATADIR}/${FAMID}/datafiles/
        aws s3 cp s3://ssc-gangstr-denovos/vcf/merged/phase4/phase4_${chrom}_merged_loci_remove.bed ${DATADIR}/${FAMID}/datafiles/
        cat ${DATADIR}/${FAMID}/datafiles/phase*_${chrom}_merged_loci_remove.bed | sort -k1,1 -k2,2n -T ${DATADIR}/${FAMID}/tmp/ | datamash groupby 1,2,3 collapse 4 >  ${DATADIR}/${FAMID}/results/${chrom}_merged_loci_remove.bed
        cat ${DATADIR}/${FAMID}/datafiles/phase*_${chrom}_merged_loci_remove.bed | sort -k1,1 -k2,2n -T ${DATADIR}/${FAMID}/tmp/ | datamash groupby 1,2,3 collapse 4 >  ${DATADIR}/${FAMID}/results/${chrom}_merged_loci_remove.bed
        aws s3 cp ${DATADIR}/${FAMID}/results/${chrom}_merged_loci_remove.bed s3://ssc-gangstr-denovos/vcf/merged/${chrom}_merged_loci_remove.bed
    fi
    ### Remove problematic loci
    vcftools --exclude-bed ${DATADIR}/${FAMID}/results/${chrom}_merged_loci_remove.bed --gzvcf ${DATADIR}/${FAMID}/datafiles/${FAMID}_${chrom}.filtered.vcf.gz  --temp ${DATADIR}/${FAMID}/tmp/ --recode --recode-INFO-all --out ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.PASS_only
    cat ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.PASS_only.recode.vcf | vcf-sort -t ${DATADIR}/${FAMID}/tmp/ | bgzip -c > ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.PASS_only.recode.vcf.gz
    tabix -p vcf ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.PASS_only.recode.vcf.gz

    aws s3 cp ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.PASS_only.recode.vcf.gz ${DENOVOVCF}/
    aws s3 cp ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.PASS_only.recode.vcf.gz.tbi ${DENOVOVCF}/
  fi

  # Run CookieMonSTR
  cmd="echo [run_cookiemonstr.sh]: Running CookieMonSTR chr${chrom}"
  cmd="${cmd}; CookieMonSTR \
  --strvcf  ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.PASS_only.recode.vcf.gz \
  --fam ${DATADIR}/${FAMID}/datafiles/fam.ped \
  --max-num-alleles ${MAXALLELES} \
  --round-alleles --include-invariant
  --output-all-loci \
  --gangstr \
  --require-all-children \
  --naive \
  --naive-expansions-frr 3,8 \
  --min-num-encl-child ${MINENCCHILD} --max-perc-encl-parent ${MINENCPERPARENT} --min-encl-match ${MINENCMATCH} --min-total-encl ${MINTOTALENC} \
  --out ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.denovos"
  echo $cmd | xargs -n1 -I% -P${THREADS} sh -c "%"
done

# Get header
cat ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.denovos.all_mutations.tab | grep chrom > ${DATADIR}/${FAMID}/results/${FAMID}_allchrs.denovos.all_mutations.tab
cat ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.denovos.locus_summary.tab | grep chrom > ${DATADIR}/${FAMID}/results/${FAMID}_allchrs.denovos.locus_summary.tab
# Join all CHRs
for chrom in {1..22}
do
    cat ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.denovos.all_mutations.tab | grep -v chrom >> ${DATADIR}/${FAMID}/results/${FAMID}_allchrs.denovos.all_mutations.tab
    cat ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.denovos.locus_summary.tab | grep -v chrom >> ${DATADIR}/${FAMID}/results/${FAMID}_allchrs.denovos.locus_summary.tab
done

## Upload results to S3
aws s3 cp  ${DATADIR}/${FAMID}/results/${FAMID}_allchrs.denovos.all_mutations.tab ${DENOVONAIVEDIR}/${FAMID}_allchrs.denovos.all_mutations.tab || die "Error uploading ${FAMID}_allchrs.denovos.all_mutations.tab "
aws s3 cp  ${DATADIR}/${FAMID}/results/${FAMID}_allchrs.denovos.locus_summary.tab ${DENOVONAIVEDIR}/${FAMID}_allchrs.denovos.locus_summary.tab || die "Error uploading ${FAMID}_allchrs.denovos.locus_summary.tab "

# QC all CHRs denovos
echo "Running CookieMonSTR filtering ${FAMID}"
grep chrom ${DATADIR}/${FAMID}/results/${FAMID}_allchrs.denovos.all_mutations.tab > ${DATADIR}/${FAMID}/results/${FAMID}_allchrs.denovos.filtered_mutations.tab
awk '{if ($8==1 || $8==-1)  print $0};' ${DATADIR}/${FAMID}/results/${FAMID}_allchrs.denovos.all_mutations.tab >> ${DATADIR}/${FAMID}/results/${FAMID}_allchrs.denovos.filtered_mutations.tab
aws s3 cp  ${DATADIR}/${FAMID}/results/${FAMID}_allchrs.denovos.filtered_mutations.tab ${DENOVONAIVEDIR}/${FAMID}_allchrs.denovos.filtered_mutations.tab || die "Error uploading ${FAMID}_allchrs.denovos.filtered_mutations.tab "


### Cleanup before moving on to next job
rm ${DATADIR}/${FAMID}/results/*
rm ${DATADIR}/${FAMID}/datafiles/*
rm ${DATADIR}/${FAMID}/tmp/*
exit 0
