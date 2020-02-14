 #!/bin/bash python3

# Run this script on snorlax inside script dir
# Run command ./aws_job_run_gangstr.sh

SSC_FILE='../files/SSC_phase31_crams_input.tab'
SSC_PARAMS='sscp3_1_params.sh'

# Upload scripts
source ../${SSC_PARAMS}
aws s3 cp ../${SSC_PARAMS} s3://ssc-denovos/scripts/ || echo "${SSC_PARAMS} upload to s3 failed"
aws s3 cp ../pipeline/run_gangstr.sh s3://ssc-denovos/scripts/ || echo "run_gangstr.sh upload to s3 failed"
aws s3 cp ../src/decrypt.py s3://ssc-denovos/scripts/ || echo "decrypt.py upload to s3 failed"
aws s3 cp ../files/encrypt_code.txt s3://ssc-denovos/scripts/encrypt_code.txt ||  echo "encrypt_code.txt upload to s3 failed"


SSC_ACCESS_KEY=$(cat ~/.aws/credentials | grep -A 2 ssc2 | grep id | cut -f 2 -d '=' | cut -f 2 -d' ' )
SSC_SECRET_ACCESS_KEY=$(cat ~/.aws/credentials | grep -A 2 ssc2 | grep secret | cut -f 2 -d '=' | cut -f 2 -d' ' )
MESSAGE=$(cat ../files/encrypt_code.txt)
ENC_SSC_ACCESS_KEY=$(python3 ../src/encrypt.py ${MESSAGE} ${SSC_ACCESS_KEY})
ENC_SSC_SECRET_ACCESS_KEY=$(python3 ../src/encrypt.py ${MESSAGE} ${SSC_SECRET_ACCESS_KEY})

# Read SSC family list
while read -r FAMID CRAMSLIST; do
  # Check VCF files exsits - if missing, then run gangSTR
   STARTCHR=0
   for chrom in $(seq 1 22)
   do
     VCF=vcf/${PHASE}/${FAMID}_${chrom}.sorted.vcf.gz
     #if file does not exist, set start chrom number to chrom
     aws s3api head-object --bucket ssc-denovos --key ${VCF} || STARTCHR=${chrom}
     if [ ${STARTCHR} -ne 0 ]; then
       break
     fi
   done

   # call submit job if any of the chroms were missing
    if [ ${STARTCHR} -ne 0 ]; then
    echo "Run GangSTR for ${FAMID} chr${STARTCHR}"
    aws batch submit-job \
        --job-name ssc-${FAMID} \
        --job-queue ssc-denovos \
        --job-definition ssc-denovos:13 \
        --container-overrides 'command=["run_gangstr.sh",'"${CRAMSLIST}"','"${FAMID}"','"${ENC_SSC_ACCESS_KEY}"','"${ENC_SSC_SECRET_ACCESS_KEY}"','"${STARTCHR}"', '"${SSC_PARAMS}"'],environment=[{name="BATCH_FILE_TYPE",value="script"},{name="BATCH_FILE_S3_URL",value="s3://ssc-denovos/scripts/run_gangstr.sh"}]'
   fi
done <${SSC_FILE}
