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
IMGT_MUTATION=$7


CODE_DIRECTORY='/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/'

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
	## Make a directory specific for this sample RUN -> this will be the location of all the LOG Files
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}
	mkdir COMMANDLOGS/${SAMPLES_FILE_POST}
	CMD="rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_*"
	eval "${CMD}"
	CMD="rm COMMANDLOGS/job_${SAMPLES_FILE_PRE}_*"
	eval "${CMD}"
	
	## Calculate number of samples to run array job on
	TASKS=$(cat ${SAMPLES_FILE_PRE} | wc -l)
	TASKS=$((TASKS+1))

	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_1.txt
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/1/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/1/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} $SAMPLES_FILE_PRE 1 $SAMPLES_FILE_POST)
	echo "${SAMPLES_FILE_PRE},TASK1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}

	##Part two uses all jobs from post file. 
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_2.txt
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_2 --dependency=afterok:${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/2/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/2/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK2,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_3.txt
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_3 --dependency=afterok:${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/3/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/3/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 3)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK3,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_4.txt
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_4 --dependency=afterok:${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/4/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/4/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 4)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK4,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_5.txt
	#Part five uses only 'one'
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_5 --dependency=afterok:${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/5/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/5/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 5)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK5,ARRAYSIZE:1,JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_RS.txt
	#Part five uses only 'one' and runs the R script 
	JOB_ID=$(sbatch --cpus-per-task=8 -p short --array=1 --job-name=BCR_TCR_PIPELINE_RS --dependency=afterok:${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK.RS,ARRAYSIZE:1,JOBID:${JOB_ID},RUNNAME:${RUNNAME},BATCH_FILE:${BATCH_FILE}"

	## Submitted all jobs stages 1-5+RS
	exit 0
fi


## running all and the jaccard functions (without consensus)
if [[ "$STAGE" -eq 2 ]]; then 
	echo "RUNNING STAGES 1-5 and RS and Jaccard"
	## Make a directory specific for this sample RUN -> this will be the location of all the LOG Files
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}
	mkdir COMMANDLOGS/${SAMPLES_FILE_POST}
	CMD="rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_*"
	eval "${CMD}"
	CMD="rm COMMANDLOGS/job_${SAMPLES_FILE_PRE}_*"
	eval "${CMD}"
	## Calculate number of samples to run array job on
	TASKS=$(cat ${SAMPLES_FILE_PRE} | wc -l)
	TASKS=$((TASKS+1))

	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_1.txt
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/1/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/1/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} $SAMPLES_FILE_PRE 1 $SAMPLES_FILE_POST)
	echo "${SAMPLES_FILE_PRE},TASK1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}

	##Part two uses all jobs from post file. 
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_2.txt
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_2 --dependency=afterok:${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/2/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/2/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK2,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_3.txt
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_3 --dependency=afterok:${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/3/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/3/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 3)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK3,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_4.txt
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_4 --dependency=afterok:${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/4/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/4/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 4)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK4,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	##Part five uses only 'one'
	JOB_ID5=$(sbatch --cpus-per-task=1 -p short --array=1 --job-name=BCR_TCR_PIPELINE_5 --dependency=afterok:${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/5/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/5/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 5)
	JOB_ID5=${JOB_ID5%.*}
	echo "${SAMPLES_FILE_POST},TASK5,ARRAYSIZE:1,JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_RS.txt
	#Part five uses only 'one' and runs the R script 
	JOB_ID=$(sbatch --cpus-per-task=8 -p short --array=1 --job-name=BCR_TCR_PIPELINE_RS --dependency=afterok:${JOB_ID5} -e COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK.RS,ARRAYSIZE:1,JOBID:${JOB_ID},RUNNAME:${RUNNAME},BATCH_FILE:${BATCH_FILE}"

	echo "RUNNING STAGES JACCARD"
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_JACCARD.txt
	JOB_ID1=$(sbatch --cpus-per-task=10 -p long --array=1 --job-name=BCR_TCR_PIPELINE_JACCARD1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ --parsable --dependency=afterok:${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 1)
	JOB_ID2=$(sbatch --cpus-per-task=10 -p long --array=1 --job-name=BCR_TCR_PIPELINE_JACCARD2 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ --parsable --dependency=afterok:${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 2)
	exit 0
fi 

## Run stage JACCARD and RS - useful if you change batch file descriptor or want to run with different pdf dimensions etc!
if [[ "$STAGE" == "R" ]]; then 
	echo "RUNNING STAGES RS and Jaccard"
	## Make a directory specific for this sample RUN -> this will be the location of all the LOG Files
	rm -r COMMANDLOGS/RS
	rm -r COMMANDLOGS/JI
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_RS.txt
	#Part five uses only 'one' and runs the R script 
	JOB_ID=$(sbatch --cpus-per-task=8 -p short -array=1 --job-name=BCR_TCR_PIPELINE_RS -e COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK.RS,ARRAYSIZE:1,JOBID:${JOB_ID},RUNNAME:${RUNNAME},BATCH_FILE:${BATCH_FILE}"

	echo "RUNNING STAGES JACCARD"	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_JACCARD.txt
	JOB_ID1=$(sbatch --cpus-per-task=10 -p long --array=1 --job-name=BCR_TCR_PIPELINE_JACCARD1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 1)
	JOB_ID2=$(sbatch --cpus-per-task=10 -p long --array=1 --job-name=BCR_TCR_PIPELINE_JACCARD2 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 2)
	exit 0
fi 




## running all and the jaccard functions + consensus 
if [[ "$STAGE" -eq 3 ]]; then 
	echo "RUNNING STAGES 1-5, RS and CONSENSUS JACCARD"
	## Make a directory specific for this sample RUN -> this will be the location of all the LOG Files
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}
	mkdir COMMANDLOGS/${SAMPLES_FILE_POST}
	CMD="rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_*"
	eval "${CMD}"
	CMD="rm COMMANDLOGS/job_${SAMPLES_FILE_PRE}_*"
	eval "${CMD}"

	## Calculate number of samples to run array job on
	TASKS=$(cat ${SAMPLES_FILE_PRE} | wc -l)
	TASKS=$((TASKS+1))

	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_1.txt
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/1/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/1/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} $SAMPLES_FILE_PRE 1 $SAMPLES_FILE_POST)
	echo "${SAMPLES_FILE_PRE},TASK1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}

	##Part two uses all jobs from post file. 
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_2.txt
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_2 --dependency=afterok:${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/2/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/2/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK2,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_3.txt
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_3 --dependency=afterok:${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/3/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/3/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 3)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK3,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_4.txt
	JOB_ID=$(sbatch --cpus-per-task=1 -p short --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_4 --dependency=afterok:${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/4/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/4/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 4)
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK4,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	#Part five uses only 'one'
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_5.txt
	JOB_ID5=$(sbatch --cpus-per-task=1 -p short --array=1 --job-name=BCR_TCR_PIPELINE_5 -e COMMANDLOGS/${SAMPLES_FILE_POST}/5/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/5/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 5)
	JOB_ID5=${JOB_ID5%.*}
	echo "${SAMPLES_FILE_POST},TASK5,ARRAYSIZE:1,JOBID:${JOB_ID5}"

	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_RS.txt
	#Part five uses only 'one' and runs the R script 
	JOB_ID=$(sbatch --cpus-per-task=8 -p short --array=1 --dependency=afterok:${JOB_ID5} --job-name=BCR_TCR_PIPELINE_RS -e COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE})
	JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK.RS,ARRAYSIZE:1,JOBID:${JOB_ID},RUNNAME:${RUNNAME},BATCH_FILE:${BATCH_FILE}"
	
	echo "RUNNING STAGES CONSENSUS and JACCARD"
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_JACCARD.txt
	JOB_ID1=$(sbatch --cpus-per-task=10 -p long --array=1 --job-name=BCR_TCR_PIPELINE_JACCARD1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ --parsable --dependency=afterok:${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 1)
	JOB_ID2=$(sbatch --cpus-per-task=10 -p long --array=1 --job-name=BCR_TCR_PIPELINE_JACCARD2 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ --parsable --dependency=afterok:${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 2)
	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_CONSENSUS.txt
	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	JOB_ID=$(sbatch --cpus-per-task=1 -p long --array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_CONSENSUS --dependency=afterok:${JOB_ID5} -e COMMANDLOGS/${SAMPLES_FILE_POST}/CONSENSES/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/CONSENSES/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} CONSENSUS)
	echo "${SAMPLES_FILE_POST},TASK.CONSENSUS,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}
	
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
    # Running the JACCARD with UMI correction 
	JOB_ID3=$(sbatch --cpus-per-task=10 -p long --array=1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ --job-name=BCR_TCR_PIPELINE_JACCARD3 -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ --parsable --dependency=afterok:${JOB_ID} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 3)

	echo "${SAMPLES_FILE_PRE},TASK.JACCARD1-3,ARRAYSIZE:1,JOBID:${JOB_ID}"
	exit 0
fi 


## Running the isotyper/subsampling functions 
if [[ "$STAGE" -eq 4 ]]; then 
	echo "RUNNING STAGES 6 and ISO1"
	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_6*.txt
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1*.txt
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_PRODUCTIVE*.txt
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_UNPRODUCTIVE*.txt
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO_COMPLETE*.txt
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_CAT*.txt
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}/6
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}/ISO1
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}/ISO1_PRODUCTIVE
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}/ISO1_UNPRODUCTIVE
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}/ISO_COMPLETE
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}/CAT
	
	JOB_ID=$(sbatch --cpus-per-task=6 -p long -array=1 --job-name=BCR_TCR_PIPELINE_6 -e COMMANDLOGS/${SAMPLES_FILE_POST}/6/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/6/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 6)
	echo "${SAMPLES_FILE_POST},TASK6,ARRAYSIZE:1},JOBID:${JOB_ID}"
	JOB_ID=${JOB_ID%.*}
	
	##Part two uses all jobs from post file.
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))

	# Run ISO1_PRODUCTIVE_NONPRODUCTIVE sequentially (requires fewer slots to become availible)
	JOB_ID1=$(sbatch --cpus-per-task=1 -p long -array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_ISO1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ISO_COMPLETE/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ISO_COMPLETE/ --parsable --dependency=afterok:${JOB_ID} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ISO1_COMPLETE)
	echo "${SAMPLES_FILE_POST},TASK.ISO1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID1}"
	JOB_ID1=${JOB_ID1%.*}
	
	# RUN CAT (run IMGT mutation = YES)
	JOB_ID4=$(sbatch --cpus-per-task=16 -p long -array=1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/ --job-name=BCR_TCR_PIPELINE_CAT -o COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/ --parsable --dependency=afterok:${JOB_ID1} cat_IMGT_files.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ${BATCH_FILE} ${CODE_DIRECTORY} YES)
	echo "${SAMPLES_FILE_POST},TASKCAT,ARRAYSIZE:1,JOBID:${JOB_ID4}"
	exit 0
fi 

if [[ "$STAGE" -eq 5 ]]; then 
	echo "RUNNING STAGES ISO1"
	#Run Jobs and capture job submission IDS
	# Part one uses pre-file
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}/ISO1
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}/ISO1_PRODUCTIVE
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}/ISO1_UNPRODUCTIVE
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}/ISO_COMPLETE
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1*.txt
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_PRODUCTIVE*.txt
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_UNPRODUCTIVE*.txt
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO_COMPLETE*.txt
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_CAT*.txt
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}/CAT
	
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))

	# Run ISO1_PRODUCTIVE_NONPRODUCTIVE sequentially (requires fewer slots to become availible)
	JOB_ID1=$(sbatch --cpus-per-task=1 -p long -array=1-${TASKS} --job-name=BCR_TCR_PIPELINE_ISO1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ISO_COMPLETE/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ISO_COMPLETE/ --parsable BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ISO1_COMPLETE)
	JOB_ID1=${JOB_ID1%.*}

	# RUN CAT (run IMGT mutation = YES)
	JOB_ID4=$(sbatch --cpus-per-task=16 -p long -array=1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/ --job-name=BCR_TCR_PIPELINE_CAT -o COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/ --parsable --dependency=afterok:${JOB_ID1} cat_IMGT_files.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ${BATCH_FILE} ${CODE_DIRECTORY} YES)
	echo "${SAMPLES_FILE_POST},TASKCAT,ARRAYSIZE:1,JOBID:${JOB_ID4}"
	exit 0
fi 

if [[ "$STAGE" -eq 6 ]]; then 
	echo "RUNNING ISOTYPER SUMMARY"
	# RUN CAT
	JOB_ID=$(sbatch --cpus-per-task=10 -p long -array=1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/ --job-name=BCR_TCR_PIPELINE_CAT -o COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/ --parsable BCR_TCR_Wrapper_2.sh NO ${SAMPLES_FILE_POST} ${BATCH_FILE} ${CODE_DIRECTORY} ${IMGT_MUTATION})
	echo "${SAMPLES_FILE_POST},TASKCAT,ARRAYSIZE:1,JOBID:${JOB_ID}"
	exit 0

fi 

#########################################################################################################################
######################################################################################################################
## ALL PREPROCESSING STAGES DONE!!!!!!
## Running Tensor/ Module Reduction 