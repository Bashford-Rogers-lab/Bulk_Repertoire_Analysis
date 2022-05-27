#!/bin/bash
#$ -cwd

# If there's an error, fail the whole script
set -e -o pipefail

# Set permissions so any user in the group can 
# read/write what it's created by the script
umask 002

echo "********************************************************"

# Job Arguments
# Extracting Outputdir
DEPENDANCIES=$1
SAMPLES_FILE_POST=$2
LAYOUTS_FILE=$3
CODE_DIRECTORY=$4
IMGT_MUTATION=$5

## Job arguments 
OUTPUTDIR=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$8}" $SAMPLES_FILE_POST)
GENE=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$4}" $SAMPLES_FILE_POST)    # Output folder path

## Location of IDs! 
SAMPLE_FILE=COMMANDLOGS/${SAMPLES_FILE_POST}_SAMPLES.txt
IDS_FILE=COMMANDLOGS/${SAMPLES_FILE_POST}_IDS.txt

echo "* Job Details"
echo "********************************************************"
echo "SGE job ID       : "$JOB_ID
echo "SGE task ID      : "$SGE_TASK_ID
echo "Run on host      : "`hostname`
echo "Operating system : "`uname -s`
echo "Username         : "`whoami`
echo "Started at       : "`date`
echo
echo "********************************************************"
echo "* Job Parameters"
echo "********************************************************"
echo "OUTPUTDIR       : ${OUTPUTDIR}"
echo "Sample File     : ${SAMPLE_FILE}"
echo "IDs File        : ${IDS_FILE}"
echo "********************************************************"

# "catch exit status 1" grep wrapper
# exit status 1 when no lines matched  - this is causing the script to fail with set -e hence the catch error
# exit status 0 when lines matched
# exit status 2 when error
c1grep() { grep "$@" || test $? = 1; }

## Load Recquired Modules
module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
echo "Loaded R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0 Module"

echo "RUNNING ISOTYPE ANALYSIS PLOTS"
cd ${CODE_DIRECTORY}
CMD="Rscript ${CODE_DIRECTORY}ISOTYPER_ANALYSIS.R -o ${OUTPUTDIR} -s ${SAMPLES_FILE_POST} -g ${GENE} -l ${LAYOUTS_FILE} -r ${IMGT_MUTATION}"
echo ${CMD}
eval "${CMD}" 
echo "DONE"

## IF JOB RUN SUCESSFULLY SAVE TO SAMPLE COUNTER FILE 
NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_SUMMARYISOTYPER.txt"
eval "${NWCMD}"

# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0