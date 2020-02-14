#!/bin/bash

#DumpSTR run on GangSTR VCF files on s3. Uses Level 1 filters. Runs per family.

#Input values
FAMID=$1
STARTCHR=${2:-1}
SSC_PARAMS=${3}

echo "[run_dumpstr_calls.sh]: ${FAMID}"

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
mkdir -p ${DATADIR}/${FAMID}/tmp/
rm ${DATADIR}/${FAMID}/results/*
rm ${DATADIR}/${FAMID}/datafiles/*
rm ${DATADIR}/${FAMID}/tmp/*

aws s3 cp s3://ssc-gangstr-denovos/scripts/ssc_shared_params.sh ${DATADIR}/${FAMID}/ssc_shared_params.sh || die "Error copying ssc_shared_params.sh"
source ${DATADIR}/${FAMID}/ssc_shared_params.sh
aws s3 cp s3://ssc-gangstr-denovos/scripts/${SSC_PARAMS} ${DATADIR}/${FAMID}/${SSC_PARAMS} || die "Error copying ${SSC_PARAMS}"
source ${DATADIR}/${FAMID}/${SSC_PARAMS}

### First, download data files needed for DumpSTR - SSC GangSTR VCF files
aws s3 cp ${GANGSTRDIR}/  ${DATADIR}/${FAMID}/datafiles/  --recursive  --exclude "*" --include "${FAMID}_*.sorted.vcf*" || die "Error copying ${FAMID}_*.sorted.vcf*"
aws s3 cp ${DUMPSTRDIR}/${PHASE}_remove_families_callrate.txt  ${DATADIR}/${FAMID}/datafiles/${PHASE}_remove_families_callrate.txt || echo "" >  ${DATADIR}/${FAMID}/datafiles/${PHASE}_remove_families_callrate.txt

if [ ${STARTCHR} == "X" ]; then
  ARRAY=(X Y)
elif [ ${STARTCHR} == "Y" ]; then
  ARRAY=(Y)
else
  #ARRAY=($(seq ${STARTCHR} 22) X Y)
  ARRAY=($(seq ${STARTCHR} 22))
fi

for chrom in "${ARRAY[@]}"
do
    cmd="echo [run_dumpstr_calls.sh]: Running DumpSTR chr${chrom}"

    ### Second, run DumpSTR - Level 1 filters
    cmd="${cmd}; dumpSTR \
    --vcf ${DATADIR}/${FAMID}/datafiles/${FAMID}_${chrom}.sorted.vcf.gz \
    --min-call-DP ${MINDP} \
    --max-call-DP  ${MAXDP} \
    --filter-spanbound-only \
    --filter-badCI \
    --require-support ${NREADS} \
    --readlen ${READLEN} \
    --out ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}"
    cmd="${cmd} && cat ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.vcf | vcf-sort -t ${DATADIR}/${FAMID}/tmp/ | bgzip -c > ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.vcf.gz"
    cmd="${cmd} && tabix -p vcf ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.vcf.gz"

    ### Third, upload results to S3
    cmd="${cmd} && aws s3 cp ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.vcf.gz ${DUMPSTRDIR}/"
    cmd="${cmd} && aws s3 cp ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.filtered.vcf.gz.tbi ${DUMPSTRDIR}/"
    cmd="${cmd} && aws s3 cp ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.loclog.tab  ${DUMPSTRDIR}/${FAMID}_${chrom}.filtered.loclog.tab"
    cmd="${cmd} && aws s3 cp ${DATADIR}/${FAMID}/results/${FAMID}_${chrom}.samplog.tab ${DUMPSTRDIR}/${FAMID}_${chrom}.filtered.samplog.tab"

    echo $cmd
done | xargs -n1 -I% -P${THREADS} sh -c "%"

### Third, create call rate file
NCHRSCOMPLETED=`ls  ${DATADIR}/${FAMID}/results/ | grep .samplog.tab | wc -l`
if [ "${NCHRSCOMPLETED}" -eq 22 ]
then
  TOTALCALLS=$(grep FILTER:PASS ${DATADIR}/${FAMID}/results/${FAMID}_*.loclog.tab | cut -f 2 | datamash sum 1)
  cat  ${DATADIR}/${FAMID}/results/${FAMID}_*.samplog.tab |  grep -v ^sample | datamash --sort  groupby 1 count 1 sum 2 > ${DATADIR}/${FAMID}/results/callrate.tab
  awk -v TOTALCALLS=${TOTALCALLS} '{print $1, $2, $3, TOTALCALLS, $3/TOTALCALLS} ' ${DATADIR}/${FAMID}/results/callrate.tab |sed  "1isample\tnumchroms\tnumcalls\ttotalcalls\tcallrate" > ${DATADIR}/${FAMID}/results/${FAMID}.filtered.callrate.tab
  aws s3 cp ${DATADIR}/${FAMID}/results/${FAMID}.filtered.callrate.tab ${DUMPSTRDIR}/${FAMID}.filtered.callrate.tab

  ### Fourth, if any sample fails call rate then remove family
  NSAMPLESFAIL=`cat ${DATADIR}/${FAMID}/results/${FAMID}.filtered.callrate.tab | awk -v MINSAMPLECALLRATE=${MINSAMPLECALLRATE} '{if($4 < MINSAMPLECALLRATE) print $1}' | wc -l`
  if [ "${NSAMPLESFAIL}" -gt 0 ]
  then
      cat ${DATADIR}/${FAMID}/datafiles/${PHASE}_remove_families_callrate.txt | sed  "1i${FAMID}"| sort -u >  ${DATADIR}/${FAMID}/results/new_remove_families.txt
      aws s3 cp ${DATADIR}/${FAMID}/results/new_remove_families.txt ${DUMPSTRDIR}/${PHASE}_remove_families_callrate.txt
  fi
fi

### Cleanup before moving on to next job
rm ${DATADIR}/${FAMID}/results/*
rm ${DATADIR}/${FAMID}/datafiles/*
rm ${DATADIR}/${FAMID}/tmp/*
exit 0
