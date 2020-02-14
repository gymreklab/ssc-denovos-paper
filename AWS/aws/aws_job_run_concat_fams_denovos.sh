#!/bin/bash

# Run this script on snorlax inside script dir
# Run command ./aws_job_concat_fams_denovos.sh


SSC_PARAMS='sscp3_2_params.sh'  ### CHANGE
DENOVORES="denovos_naive_Jan20"

source ../ssc_shared_params.sh
source ../${SSC_PARAMS}
aws s3 cp ../ssc_shared_params.sh s3://ssc-gangstr-denovos/scripts/ || echo "ssc_shared_params.sh upload to s3 failed"
aws s3 cp ../${SSC_PARAMS} s3://ssc-gangstr-denovos/scripts/ || echo "${SSC_PARAMS} upload to s3 failed"
aws s3 cp ../pipeline/concat_fams_denovos.sh s3://ssc-gangstr-denovos/scripts/concat_fams_denovos.sh || echo "concat_fams_denovos.sh upload to s3 failed"


aws batch submit-job \
   --job-name ssc-concact \
   --job-queue ssc-denovo \
   --job-definition ssc-denovos:5  \
   --container-overrides 'command=["concat_fams_denovos.sh", '"${SSC_PARAMS}"', '"${DENOVORES}"'],environment=[{name="BATCH_FILE_TYPE",value="script"},{name="BATCH_FILE_S3_URL",value="s3://ssc-gangstr-denovos/scripts/concat_fams_denovos.sh"}]'
