#!/bin/bash

#$ -cwd
#$ -N BCR_TCR_p2
#$ -q short.qc
#$ -pe shmem 1
#$ -e COMMANDLOGS/
#$ -o COMMANDLOGS/

# If there's an error, fail the whole script
set -e -o pipefail

# Set permissions so any user in the group can 
# read/write what it's created by the script
umask 002

SAMPLES_FILE=$1
STAGE=$3

TASKS=$(cat ${SAMPLES_FILE} | wc -l)
TASKS=$((TASKS+1))

FILE="job_$SAMPLES_FILE_$3.txt"
LENGTH=$(cat ${FILE} | wc -l) 

if [[ "$LENGTH" == "$TASKS" ]]; then
	echo "YES"
else 
	echo "NO"
fi 

exit 0