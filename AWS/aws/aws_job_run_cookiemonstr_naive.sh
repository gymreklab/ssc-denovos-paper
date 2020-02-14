#!/bin/bash

# Run this script on snorlax inside script dir
# Run command ./aws_job_run_cookiemonstr_naive.sh


SSC_FILE='../files/SSC_phase4_crams_input.tab' ### CHANGE
SSC_PARAMS='sscp4_params.sh'  ### CHANGE

source ../ssc_shared_params.sh
source ../${SSC_PARAMS}
aws s3 cp ../ssc_shared_params.sh s3://ssc-gangstr-denovos/scripts/ || echo "ssc_shared_params.sh upload to s3 failed"
aws s3 cp ../${SSC_PARAMS} s3://ssc-gangstr-denovos/scripts/ || echo "${SSC_PARAMS} upload to s3 failed"
aws s3 cp ../pipeline/run_cookiemonstr_naive.sh s3://ssc-gangstr-denovos/scripts/run_cookiemonstr_naive.sh|| echo "run_cookiemonstr_naive.sh upload to s3 failed"

while read -r FAMID CRAMSLIST; do
     RUN=false
     ### Check .denovos.all_mutations.tab files exsits - if exist, then do not re-run CookieMonSTR
     FILE=denovos_naive_Jan20/${PHASE}/${FAMID}_allchrs.denovos.filtered_mutations.tab
     aws s3api head-object --bucket ssc-gangstr-denovos --key ${FILE} || RUN=true
     if (${RUN}); then
       echo "Run cookiemonstr for ${FAMID}"
       aws batch submit-job \
           --job-name denovo-ssc-${FAMID} \
           --job-queue ssc-denovo-4 \
           --job-definition ssc-denovos:5  \
           --container-overrides 'command=["run_cookiemonstr_naive.sh",'"${FAMID}"', '"${SSC_PARAMS}"'],environment=[{name="BATCH_FILE_TYPE",value="script"},{name="BATCH_FILE_S3_URL",value="s3://ssc-gangstr-denovos/scripts/run_cookiemonstr_naive.sh"}]'
     fi
done <${SSC_FILE}
