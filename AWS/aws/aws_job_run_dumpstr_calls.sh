#!/bin/bash

# Run this script on snorlax inside script dir
# Run command ./aws_job_run_dumpstr_calls_calls.sh
# Creates job for 1 family VCF running dumpSTR for removing bad calls per sample


SSC_FILE='../files/SSC_phase2_crams_input.tab'
SSC_PARAMS='sscp2_params.sh'

source ../ssc_shared_params.sh
source ../${SSC_PARAMS}
aws s3 cp ../ssc_shared_params.sh s3://ssc-gangstr-denovos/scripts/ || echo "ssc_shared_params.sh upload to s3 failed"
aws s3 cp ../${SSC_PARAMS} s3://ssc-gangstr-denovos/scripts/ || echo "${SSC_PARAMS} upload to s3 failed"
aws s3 cp ../pipeline/run_dumpstr_calls.sh s3://ssc-gangstr-denovos/scripts/run_dumpstr_calls.sh || echo "run_dumpstr_calls.sh upload to s3 failed"


while read -r FAMID CRAMSLIST; do
   # sorted_vcf_not_missing=true
   # # Check .sorted.vcf.gz VCF files exists - if missing, then do NOT run dumpSTR
   # #for chrom in $(seq 1 22) X Y
   # for chrom in $(seq 1 22)
   # do
   #   VCF=vcf/${PHASE}/${FAMID}_${chrom}.sorted.vcf.gz
   #   aws s3api head-object --bucket ssc-gangstr-denovos --key ${VCF} || sorted_vcf_not_missing=false
   #  if (! ${sorted_vcf_not_missing}); then
   #    echo "File missing: s3://ssc-gangstr-denovos/${VCF}"
   #    break
   #  fi
   # done
   # if (${sorted_vcf_not_missing}); then
       ### if .sorted.vcf.gz VCF files exsits
       ### Then, check call rate file b/c we don't want to run dumpster if it has already been run
   STARTCHR=0
   FILTEREDVCF=vcf/filtered2/${PHASE}/${FAMID}.filtered.callrate.tab
   #if file does not exist, set start chrom number to chrom
   aws s3api head-object --bucket ssc-gangstr-denovos --key ${FILTEREDVCF} || STARTCHR=1
   # call submit job if any of the chroms were missing
   if [ ${STARTCHR} -ne 0 ]; then
       echo "Run dumpstr for ${FAMID} chr${STARTCHR}-22"
       aws batch submit-job \
           --job-name qc-ssc-${FAMID} \
           --job-queue ssc-denovo \
           --job-definition test-ileena:2 \
           --container-overrides 'command=["run_dumpstr_calls.sh",'"${FAMID}"', '"${STARTCHR}"', '"${SSC_PARAMS}"'],environment=[{name="BATCH_FILE_TYPE",value="script"},{name="BATCH_FILE_S3_URL",value="s3://ssc-gangstr-denovos/scripts/run_dumpstr_calls.sh"}]'
    fi
    # fi
done <${SSC_FILE}
