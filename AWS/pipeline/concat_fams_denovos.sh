#!/bin/bash

#Concats all de novo mutation files for all families in 1 phase.

#Input values
SSC_PARAMS=${1}
DENOVORES=${2}

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

# #############
DATADIR=/scratch # This is where we have all the EBS storage space mounted
# Set up
rm -r ${DATADIR}/*
mkdir -p ${DATADIR}/datafiles
mkdir -p ${DATADIR}/results
mkdir -p ${DATADIR}/tmp

# Copy Files
aws s3 cp s3://ssc-gangstr-denovos/scripts/ssc_shared_params.sh ${DATADIR}/ssc_shared_params.sh || die "Error copying ssc_shared_params.sh"
source ${DATADIR}/ssc_shared_params.sh
aws s3 cp s3://ssc-gangstr-denovos/scripts/${SSC_PARAMS} ${DATADIR}/${SSC_PARAMS} || die "Error copying ${SSC_PARAMS}"
source ${DATADIR}/${SSC_PARAMS}
mkdir -p ${DATADIR}/datafiles/${PHASE}
mkdir -p ${DATADIR}/results/${PHASE}


cd ${DATADIR}

for FILE in "filtered_mutations"
do
  FAMFILE=_allchrs.denovos.${FILE}.tab
  OUTFILE=SSC_${PHASE}_011720_${DENOVORES}.${FILE}.tab

  # Download SSC de novo mutation files
  #aws s3 cp s3://ssc-gangstr-denovos/denovos_GW_priors_Jan20/phase1/ . --recursive --exclude "*"  --include "*.filtered_mutations.tab*" --dryrun
  RESDIR=s3://ssc-gangstr-denovos/${DENOVORES}/${PHASE}
  aws s3 cp ${RESDIR}/ ${DATADIR}/datafiles/${PHASE}/ --recursive  --exclude "*"  --include "*.${FILE}.tab*"

  #Add header
  echo 'chrom	pos	period	prior	family	child	phenotype	posterior	newallele	mutsize	inparents	poocase	isnew	case_count	ctrl_count	unk_count	child_gt	mat_gt	pat_gt	encl_child	encl_mother	encl_father	encl_parent	long_mother	long_father' > ${DATADIR}/results/${PHASE}/${OUTFILE}

  for RESULT in ${DATADIR}/datafiles/${PHASE}/*${FAMFILE}
  do
    cat  ${RESULT} | grep -v ^chrom >> ${DATADIR}/results/${PHASE}/${OUTFILE}
  done

  gzip -f ${DATADIR}/results/${PHASE}/${OUTFILE}
  echo "Created file: ${DATADIR}/results/${PHASE}/${OUTFILE}.gz"
  aws s3 cp ${DATADIR}/results/${PHASE}/${OUTFILE}.gz ${RESDIR}/${OUTFILE}.gz

  rm ${DATADIR}/datafiles/${PHASE}/*
done


rm -r ${DATADIR}/*
