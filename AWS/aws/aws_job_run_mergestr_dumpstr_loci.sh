#!/bin/bash

# Run this script on snorlax inside script dir
# Run this script after dumpSTR VCF filtering has been done for all families in phase

# Select phase to run
SSC_PARAMS='sscp1_params.sh'

source ../ssc_shared_params.sh
source ../${SSC_PARAMS}
aws s3 cp ../ssc_shared_params.sh s3://ssc-gangstr-denovos/scripts/ || echo "ssc_shared_params.sh upload to s3 failed"
aws s3 cp ../${SSC_PARAMS} s3://ssc-gangstr-denovos/scripts/ || echo "${SSC_PARAMS} upload to s3 failed"
aws s3 cp ../pipeline/run_mergestr_dumpstr_loci.sh s3://ssc-gangstr-denovos/scripts/run_mergestr_dumpstr_loci.sh || echo "run_mergestr_dumpstr_loci.sh upload to s3 failed"

for chrom in {1..22}
do
 RUN=false
  ###if results file does not exist, then run job for chrom
  VCF=/vcf/merged/${PHASE}/${PHASE}_${chrom}_merged_loci_remove.bed
  aws s3api head-object --bucket ssc-gangstr-denovos-denovos --key ${VCF} || RUN=true
  if (${RUN}); then
   echo "Run mergeSTR + dumpSTR loci for chr${chrom}"
   aws batch submit-job \
       --job-name ssc-merge-${chrom} \
       --job-queue ssc-500GB \
       --job-definition ssc-denovos:4 \
       --container-overrides 'command=["run_mergestr_dumpstr_loci.sh",'"${chrom}"', '"${SSC_PARAMS}"'],environment=[{name="BATCH_FILE_TYPE",value="script"},{name="BATCH_FILE_S3_URL",value="s3://ssc-gangstr-denovos/scripts/run_mergestr_dumpstr_loci.sh"}]'
  fi
done
