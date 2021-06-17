#!/bin/bash
## Pipeline Manager for BCR_TCR_ANALYSIS Pipeline 
## Author: Lauren Overend: lauren.overend@oriel.ox.ac.uk
## Lab: Bashford Rogers
## May 2021

# Job Arguments
# File one containing the libraries
# File two containing the PCR multiplexed samples 
SAMPLES_FILE_PRE=$1
SAMPLES_FILE_POST=$2
STAGE=$3
RUNNAME=$4
BATCH_FILE=$5

## Make a directory specific for this sample RUN -> this will be the location of all the LOG Files
mkdir COMMANDLOGS/${SAMPLES_FILE_POST}

## RUN DEMULTIPLEXING STAGE 1 USING SAMPLES_FILE_PRE

if [[ "$STAGE" == "1" ]]; then 
	echo "RUNNING STAGES 1-5 and RS"
	## Calculate number of samples to run array job on
	TASKS=$(cat ${SAMPLES_FILE_PRE} | wc -l)
	TASKS=$((TASKS+1))

	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	JOB_ID=$(qsub -t 1-${TASKS} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh $SAMPLES_FILE_PRE 1 $SAMPLES_FILE_POST)
	echo "${SAMPLES_FILE_PRE},TASK1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}

	##Part two uses all jobs from post file. 
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))

	JOB_ID=$(qsub -t 1-${TASKS} -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK2,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	JOB_ID=$(qsub -t 1-${TASKS} -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${SAMPLES_FILE_POST} 3)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK3,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	JOB_ID=$(qsub -t 1-${TASKS} -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${SAMPLES_FILE_POST} 4)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK4,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	#Part five uses only 'one'
	JOB_ID=$(qsub -t 1 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${SAMPLES_FILE_POST} 5)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK5,ARRAYSIZE:1,JOBID:${JOB_ID}"

	#Part five uses only 'one' and runs the R script 
	JOB_ID=$(qsub -t 1 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${SAMPLES_FILE_POST} RS ${RUNNAME})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK.RS,ARRAYSIZE:1,JOBID:${JOB_ID},RUNNAME:${RUNNAME}"

	## Submitted all jobs stages 1-5+RS
	exit 0
fi

## Running the isotyper/subsampling functions 
if [[ "$STAGE" == "2" ]]; then 
	echo "RUNNING STAGES 6 and ISO1"
	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	JOB_ID=$(qsub -t 1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh $SAMPLES_FILE_PRE 6)
	echo "${SAMPLES_FILE_PRE},TASK6,ARRAYSIZE:1},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}
	
	##Part two uses all jobs from post file.
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))

	#Run Jobs and capture job submission IDS
	JOB_ID=$(qsub -t 1-${TASKS} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse -hold_jid ${JOB_ID} BCR_TCR_Wrapper_Cluster.sh $SAMPLES_FILE_PRE ISO1)
	echo "${SAMPLES_FILE_PRE},TASK.ISO1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}
	
	exit 0

fi 


## running the jaccard functions 
if [[ "$STAGE" == "3" ]]; then 
	echo "RUNNING STAGES CONSENSUS and JACCARD"
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	JOB_ID=$(qsub -t 1-${TASKS} -terse BCR_TCR_Wrapper_Cluster.sh $SAMPLES_FILE_PRE CONSENSUS)
	echo "${SAMPLES_FILE_PRE},TASK.CONSENSUS,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}
	
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
	JOB_ID=$(qsub -t 1 -terse -hold_jid ${JOB_ID} BCR_TCR_Wrapper_Cluster_HighMemory.sh $SAMPLES_FILE_PRE JACCARD ${RUNNAME} ${BATCH_FILE})
	echo "${SAMPLES_FILE_PRE},TASK.JACCARD,ARRAYSIZE:1,JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}
	
	exit 0
fi 



