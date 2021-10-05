#!/bin/bash
## Pipeline Manager for BCR_TCR_ANALYSIS Pipeline 
## Author: Lauren Overend: lauren.overend@oriel.ox.ac.uk
## Lab: Bashford Rogers
## May 2021
# Job Arguments
# File one containing the libraries
# File two containing the PCR multiplexed samples 
DEPENDANCIES=$1
SAMPLES_FILE_PRE=$2
SAMPLES_FILE_POST=$3
STAGE=$4
RUNNAME=$5
BATCH_FILE=$6

## Make a directory specific for this sample RUN -> this will be the location of all the LOG Files
mkdir COMMANDLOGS/${SAMPLES_FILE_POST}

## Make an IDs File
IDS=$(awk -F '\t' "{ print \$1 }" $SAMPLES_FILE_POST)
for item in "$IDS"
	do  
		printf "${item}\n" > COMMANDLOGS/${SAMPLES_FILE_POST}_IDS.txt
	done 
# Make a Samples File
SAMPLE==$(awk -F '\t' "{ print \$2 }" $SAMPLES_FILE_POST)
for item in "$SAMPLE"
	do  
		printf "${item}\n" > COMMANDLOGS/${SAMPLES_FILE_POST}_SAMPLES.txt
	done




## RUN DEMULTIPLEXING STAGE 1 USING SAMPLES_FILE_PRE

if [[ "$STAGE" -eq 1 ]]; then 
	echo "RUNNING STAGES 1-5 and RS"
	## Calculate number of samples to run array job on
	TASKS=$(cat ${SAMPLES_FILE_PRE} | wc -l)
	TASKS=$((TASKS+1))

	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_1.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} $SAMPLES_FILE_PRE 1 $SAMPLES_FILE_POST)
	echo "${SAMPLES_FILE_PRE},TASK1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}

	##Part two uses all jobs from post file. 
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_2.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_2 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK2,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_3.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_3 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 3)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK3,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_4.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_4 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 4)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK4,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_5.txt
	#Part five uses only 'one'
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1 -N BCR_TCR_PIPELINE_5 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 5)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK5,ARRAYSIZE:1,JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_RS.txt
	#Part five uses only 'one' and runs the R script 
	JOB_ID=$(qsub -pe shmem 8 -q short.qc -t 1 -N BCR_TCR_PIPELINE_RS -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK.RS,ARRAYSIZE:1,JOBID:${JOB_ID},RUNNAME:${RUNNAME},BATCH_FILE:${BATCH_FILE}"

	## Submitted all jobs stages 1-5+RS
	exit 0
fi


## running all and the jaccard functions (without consensus)
if [[ "$STAGE" -eq 2 ]]; then 
	echo "RUNNING STAGES 1-5 and RS"
	## Calculate number of samples to run array job on
	TASKS=$(cat ${SAMPLES_FILE_PRE} | wc -l)
	TASKS=$((TASKS+1))

	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_1.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} $SAMPLES_FILE_PRE 1 $SAMPLES_FILE_POST)
	echo "${SAMPLES_FILE_PRE},TASK1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}

	##Part two uses all jobs from post file. 
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_2.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_2 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK2,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_3.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_3 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 3)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK3,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_4.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_4 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 4)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK4,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	#Part five uses only 'one'
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_5.txt
	JOB_ID5=$(qsub -pe shmem 1 -q short.qc -t 1 -hold_jid ${JOB_ID} -N BCR_TCR_PIPELINE_5 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 5)
	JOB_ID5=${JOB_ID5%.*}
	echo "${SAMPLES_FILE_POST},TASK5,ARRAYSIZE:1,JOBID:${JOB_ID5}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_RS.txt
	#Part five uses only 'one' and runs the R script 
	JOB_ID=$(qsub -pe shmem 8 -q short.qc -t 1 -hold_jid ${JOB_ID5} -N BCR_TCR_PIPELINE_RS -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK.RS,ARRAYSIZE:1,JOBID:${JOB_ID},RUNNAME:${RUNNAME},BATCH_FILE:${BATCH_FILE}"
	
	echo "RUNNING STAGES JACCARD"
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_JACCARD.txt
	JOB_ID1=$(qsub -pe shmem 10 -q long.qc -t 1 -N BCR_TCR_PIPELINE_JACCARD1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse -hold_jid ${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 1)
	JOB_ID2=$(qsub -pe shmem 10 -q long.qc -t 1 -N BCR_TCR_PIPELINE_JACCARD2 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse -hold_jid ${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 2)
	 
	exit 0
fi 


## running all and the jaccard functions + consensus 
if [[ "$STAGE" -eq 3 ]]; then 
	echo "RUNNING STAGES 1-5, RS and CONSENSUS JACCARD"
	## Calculate number of samples to run array job on
	TASKS=$(cat ${SAMPLES_FILE_PRE} | wc -l)
	TASKS=$((TASKS+1))

	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_1.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} $SAMPLES_FILE_PRE 1 $SAMPLES_FILE_POST)
	echo "${SAMPLES_FILE_PRE},TASK1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}

	##Part two uses all jobs from post file. 
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_2.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_2 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK2,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_3.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_3 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 3)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK3,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_4.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_4 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 4)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK4,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	#Part five uses only 'one'
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_5.txt
	JOB_ID5=$(qsub -pe shmem 1 -q short.qc -t 1 -hold_jid ${JOB_ID} -N BCR_TCR_PIPELINE_5 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 5)
	JOB_ID5=${JOB_ID5%.*}
	echo "${SAMPLES_FILE_POST},TASK5,ARRAYSIZE:1,JOBID:${JOB_ID5}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_RS.txt
	#Part five uses only 'one' and runs the R script 
	JOB_ID=$(qsub -pe shmem 8 -q short.qc -t 1 -hold_jid ${JOB_ID5} -N BCR_TCR_PIPELINE_RS -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK.RS,ARRAYSIZE:1,JOBID:${JOB_ID},RUNNAME:${RUNNAME},BATCH_FILE:${BATCH_FILE}"
	
	echo "RUNNING STAGES CONSENSUS and JACCARD"
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_JACCARD.txt
	JOB_ID1=$(qsub -pe shmem 10 -q long.qc -t 1 -N BCR_TCR_PIPELINE_JACCARD1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse -hold_jid ${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 1)
	JOB_ID2=$(qsub -pe shmem 10 -q long.qc -t 1 -N BCR_TCR_PIPELINE_JACCARD2 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse -hold_jid ${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 2)
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_CONSENSUS.txt
	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	JOB_ID=$(qsub -pe shmem 1 -q long.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_CONSENSUS -hold_jid ${JOB_ID5} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} CONSENSUS)
	echo "${SAMPLES_FILE_PRE},TASK.CONSENSUS,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}
	
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
    # Running the JACCARD with UMI correction 
	JOB_ID3=$(qsub -pe shmem 10 -q long.qc -t 1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -N BCR_TCR_PIPELINE_JACCARD3 -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse -hold_jid ${JOB_ID} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 3)

	echo "${SAMPLES_FILE_PRE},TASK.JACCARD1-3,ARRAYSIZE:1,JOBID:${JOB_ID}"
	exit 0
fi 


## Running the isotyper/subsampling functions 
if [[ "$STAGE" -eq 4 ]]; then 
	echo "RUNNING STAGES 6 and ISO1"
	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_6.txt
	JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1 -N BCR_TCR_PIPELINE_6 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 6)
	echo "${SAMPLES_FILE_PRE},TASK6,ARRAYSIZE:1},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}
	
	##Part two uses all jobs from post file.
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1.txt
	#Run Jobs and capture job submission IDS
	JOB_ID1=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_ISO1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse -hold_jid ${JOB_ID} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ISO1)
	echo "${SAMPLES_FILE_PRE},TASK.ISO1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID1}"
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_PRODUCTIVE.txt
	JOB_ID2=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_ISO1_PRODUCTIVE -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse -hold_jid ${JOB_ID} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ISO1_PRODUCTIVE)
	echo "${SAMPLES_FILE_PRE},TASK.ISO1_PRODUCTIVE,ARRAYSIZE:${TASKS},JOBID:${JOB_ID2}"
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_NON_PRODUCTIVE.txt
	JOB_ID3=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -N BCR_TCR_PIPELINE_ISO1_UNPRODUCTIVE -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse -hold_jid ${JOB_ID} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ISO1_NON_PRODUCTIVE)
	echo "${SAMPLES_FILE_PRE},TASK.ISO1_NONPRODUCTIVE,ARRAYSIZE:${TASKS},JOBID:${JOB_ID3}"
	
	exit 0

fi 


if [[ "$STAGE" -eq 5 ]]; then 
	echo "RUNNING STAGES ISO1"
	#Run Jobs and capture job submission IDS
	# Part one uses pre-file

	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1.txt
	#Run Jobs and capture job submission IDS
	JOB_ID1=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_ISO1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse  BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ISO1)
	echo "${SAMPLES_FILE_PRE},TASK.ISO1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID1}"
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_PRODUCTIVE.txt
	JOB_ID2=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N BCR_TCR_PIPELINE_ISO1_PRODUCTIVE -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ISO1_PRODUCTIVE)
	echo "${SAMPLES_FILE_PRE},TASK.ISO1_PRODUCTIVE,ARRAYSIZE:${TASKS},JOBID:${JOB_ID2}"
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_NON_PRODUCTIVE.txt
	JOB_ID3=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -e COMMANDLOGS/${SAMPLES_FILE_POST}/ -N BCR_TCR_PIPELINE_ISO1_UNPRODUCTIVE -o COMMANDLOGS/${SAMPLES_FILE_POST}/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ISO1_NON_PRODUCTIVE)
	echo "${SAMPLES_FILE_PRE},TASK.ISO1_NONPRODUCTIVE,ARRAYSIZE:${TASKS},JOBID:${JOB_ID3}"
	
	
	exit 0

fi 

