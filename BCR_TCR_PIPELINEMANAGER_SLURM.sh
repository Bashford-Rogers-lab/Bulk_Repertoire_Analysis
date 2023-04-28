#!/bin/bash

###----------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------
## Pipeline Manager for BCR_TCR_ANALYSIS Pipeline using SLURM!!!!!
## Author: Lauren Overend: lauren.overend@oriel.ox.ac.uk
## Lab: Bashford Rogers
## May 2021 - updated May 2023 for SLURM 


###----------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------
## Job Arguments - Please provide these files: 
## DEPENDENCIES: Run dependencies - This will check that previous jobs ran successfully to completion - RECOMEND YES or Y
## SAMPLES_FILE_PRE: File one (pre) containing the libraries aka test_pre.txt 
## SAMPLES_FILE_POST:File two containing the PCR multiplexed samples  aka test_post.txt
## STAGE: Which part of this pipeline do you want to run. Minimum is 1 then upload to IMGT then 4
## RUNNAME: e.g. "SEPSIS_BCR" 
## BATCH_FILE: Layouts file detailing sample organisation e.g. lanes/libraries etc 
## TECHNICAL_SAMPLES: subset of BATCH_FILE detailing which are technical samples 
## SAMPLES_EXCLUDE: Only include if there are any samples you want to exclude all together - e.g. you discover they have significant comorbidities 
DEPENDANCIES=$1
SAMPLES_FILE_PRE=$2
SAMPLES_FILE_POST=$3
STAGE=$4
RUNNAME=$5
BATCH_FILE=$6
TECHNICAL_SAMPLES=$7
SAMPLES_EXCLUDE=$8

## HARD CODED YOU MAY NEED TO EDIT 
## MOST PEOPLE WILL WANT NORMALISE NORMALISE_VGENEUSAGE==
CODE_DIRECTORY='/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/'
NORMALISE_VGENEUSAGE=TRUE
## If this is BCR set to either YES or NO - if TCR it doesnt matter as it wont run 
IMGT_MUTATION=YES 
GENE=$(awk -F '\t' "{if (NR==1) print \$4}" $SAMPLES_FILE_POST)  

## Set memory for Isotyper summary script 
if [[ "$GENE" == "IGH" ]]; then
	echo "BCR Repertoire ANALYSIS"
	NODES=10
fi

if [[ "$GENE" == "TCR" ]]; then
	echo "TCR Repertoire ANALYSIS"
	NODES=1
fi

echo "RUNNING SLURM PIPELINE"

###----------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------
### IF BATCH, TECHNICAL FILE OR SAMPLES EXCLUDE ARE NOT SUPPLIED SET TO DEFAULT
if [ -z "$6" ]
  then
    echo "No Batch file supplied defaults to FALSE"
	BATCH_FILE=FALSE
fi

if [ -z "$7" ]
  then
    echo "No technical file supplied defaults to FALSE"
	TECHNICAL_SAMPLES=FALSE
fi

if [ -z "$5" ]
  then
    echo "No runname supplied. Defaults to: TRY"
	RUNNAME=REPSEQ
fi

if [ -z "$8" ]
  then
    echo "No Samples to Exclude supplied. Defaults to: NA"
	SAMPLES_EXCLUDE=NA
fi

###----------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------
### This is important for dependency checking
### CREATE SAMPLES ID FILE 

IDS=$(awk -F '\t' "{ print \$1 }" $SAMPLES_FILE_POST)
for item in "$IDS"
	do  
		printf "${item}\n" > COMMANDLOGS/${SAMPLES_FILE_POST}_IDS.txt
	done 
	
### CREATE LIBRARIES ID FILE 
SAMPLE==$(awk -F '\t' "{ print \$2 }" $SAMPLES_FILE_POST)
for item in "$SAMPLE"
	do  
		printf "${item}\n" > COMMANDLOGS/${SAMPLES_FILE_POST}_SAMPLES.txt
	done


###----------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------
## OPTION 1:     RUN BASIC PREPROCESSING AND R QC ANALYSIS 
##               AT THE END OF THIS YOU WILL NEED TO UPLOAD THE FULLY REDUCED FILE TO IMGT FOR ANNOTATION 
##               YOU WILL ALSO NEED TO SET SUBSAMPLE DEPTH FOR ISOTYPER 

if [[ "$STAGE" -eq 1 ]]; then 
	echo "RUNNING STAGES 1-5 and RS"
	
	# SETUP: Make a directory specific for this sample RUN -> this will be the location of all the LOG Files
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}
	mkdir COMMANDLOGS/${SAMPLES_FILE_POST}
	CMD="rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_*"
	eval "${CMD}"
	CMD="rm COMMANDLOGS/job_${SAMPLES_FILE_PRE}_*"
	eval "${CMD}"
	
	# TASK NUMBER: CALCULATE NUMBER OF LIBRARIES 
	TASKS=$(cat ${SAMPLES_FILE_PRE} | wc -l)
	TASKS=$((TASKS+1))

	# RUN STAGE 1
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_1.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/1/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short --cpus-per-task=1 -J STAGE_1 -e ${OUT}/STAGE_1_%a.e -o ${OUT}/STAGE_1_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_PRE} 1 ${SAMPLES_FILE_POST})
	echo "${SAMPLES_FILE_PRE},TASK1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	# TASK NUMBER: CALCULATE NUMBER OF DEMULTIPLEXED SAMPLES
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
	# RUN STAGE 2
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_2.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/2/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_2 -e ${OUT}/STAGE_2_%a.e -o ${OUT}/STAGE_2_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	#JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N STAGE_2 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/2/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/2/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	echo "${SAMPLES_FILE_POST},TASK2,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	# RUN STAGE 3
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_3.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/3/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_3 -e ${OUT}/STAGE_3_%a.e -o ${OUT}/STAGE_3_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 3)
	#JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N STAGE_3 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/3/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/3/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 3)
	echo "${SAMPLES_FILE_POST},TASK3,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	
	# RUN STAGE 4 
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_4.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/4/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_4 -e ${OUT}/STAGE_4_%a.e -o ${OUT}/STAGE_4_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 4)
	#JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N STAGE_4 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/4/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/4/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 4)
	#JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK4,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	
	# RUN STAGE 5: uses only 'one' task 
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_5.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/5/
	mkdir ${OUT}
	#Part five uses only 'one task e.g. not an array'
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_5 -e ${OUT}/STAGE_5_%a.e -o ${OUT}/STAGE_5_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 5)
	#JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1 -N STAGE_5 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/5/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/5/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 5)
	echo "${SAMPLES_FILE_POST},TASK5,ARRAYSIZE:1,JOBID:${JOB_ID}"

	# RUN STAGE RS: uses only 'one' task 
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_RS.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/RS/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=8 -J STAGE_RS -e ${OUT}/STAGE_RS_%a.e -o ${OUT}/STAGE_RS_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh  ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE} ${TECHNICAL_SAMPLES})
	#JOB_ID=$(qsub -pe shmem 8 -q short.qc -t 1 -N STAGE_RS -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE} ${TECHNICAL_SAMPLES})
	echo "${SAMPLES_FILE_POST},TASK.RS,ARRAYSIZE:1,JOBID:${JOB_ID},RUNNAME:${RUNNAME},BATCH_FILE:${BATCH_FILE}"

	## Submitted all jobs stages 1-5+RS
	exit 0
fi

###----------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------
## OPTION 2:     RUN BASIC PREPROCESSING AND R QC ANALYSIS + JACCARD (sample overlap) 
##               AT THE END OF THIS YOU WILL NEED TO UPLOAD THE FULLY REDUCED FILE TO IMGT FOR ANNOTATION 
##               YOU WILL ALSO NEED TO SET SUBSAMPLE DEPTH FOR ISOTYPER 

if [[ "$STAGE" -eq 2 ]]; then 
	echo "RUNNING STAGES 1-5 and RS and Jaccard"
	
	# SETUP: Make a directory specific for this sample RUN -> this will be the location of all the LOG Files
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}
	mkdir COMMANDLOGS/${SAMPLES_FILE_POST}
	CMD="rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_*"
	eval "${CMD}"
	CMD="rm COMMANDLOGS/job_${SAMPLES_FILE_PRE}_*"
	eval "${CMD}"
	
	# TASK NUMBER: CALCULATE NUMBER OF LIBRARIES 
	TASKS=$(cat ${SAMPLES_FILE_PRE} | wc -l)
	TASKS=$((TASKS+1))

	# RUN STAGE 1
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_1.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/1/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short --cpus-per-task=1 -J STAGE_1 -e ${OUT}/STAGE_1_%a.e -o ${OUT}/STAGE_1_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_PRE} 1 ${SAMPLES_FILE_POST})
	echo "${SAMPLES_FILE_PRE},TASK1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	# TASK NUMBER: CALCULATE NUMBER OF DEMULTIPLEXED SAMPLES
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
	# RUN STAGE 2
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_2.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/2/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_2 -e ${OUT}/STAGE_2_%a.e -o ${OUT}/STAGE_2_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	#JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N STAGE_2 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/2/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/2/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	echo "${SAMPLES_FILE_POST},TASK2,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	# RUN STAGE 3
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_3.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/3/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_3 -e ${OUT}/STAGE_3_%a.e -o ${OUT}/STAGE_3_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 3)
	#JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N STAGE_3 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/3/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/3/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 3)
	echo "${SAMPLES_FILE_POST},TASK3,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"
	
	# RUN STAGE 4
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_4.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/4/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_4 -e ${OUT}/STAGE_4_%a.e -o ${OUT}/STAGE_4_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 4)
	#JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N STAGE_4 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/4/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/4/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 4)
	#JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK4,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	# RUN STAGE 5: uses only 'one' task 
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_5.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/5/
	mkdir ${OUT}
	#Part five uses only 'one task e.g. not an array'
	JOB5_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_5 -e ${OUT}/STAGE_5_%a.e -o ${OUT}/STAGE_5_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 5)
	#JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1 -N STAGE_5 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/5/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/5/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 5)
	echo "${SAMPLES_FILE_POST},TASK5,ARRAYSIZE:1,JOBID:${JOB_ID}"

	# RUN STAGE RS: uses only 'one' task 
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_RS.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/RS/
	mkdir ${OUT}
	#Part five uses only 'one' and runs the R script 
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOB5_ID} --cpus-per-task=8 -J STAGE_RS -e ${OUT}/STAGE_RS_%a.e -o ${OUT}/STAGE_RS_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE} ${TECHNICAL_SAMPLES})
	#JOB_ID=$(qsub -pe shmem 8 -q short.qc -t 1 -N STAGE_RS -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE} ${TECHNICAL_SAMPLES})
	echo "${SAMPLES_FILE_POST},TASK.RS,ARRAYSIZE:1,JOBID:${JOB_ID},RUNNAME:${RUNNAME},BATCH_FILE:${BATCH_FILE}"

	# RUN STAGE JACCARD: uses only 'one' task 
	# Note this will be submitted in parallel to RS analysis 
	echo "RUNNING STAGES JACCARD "
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_JACCARD.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/JI/
	mkdir ${OUT}
	JOBA_ID1=$(sbatch --parsable -p short -d afterok:${JOB5_ID} --cpus-per-task=10 -J STAGE_J1 -e ${OUT}/STAGE_JACCARD1_%a.e -o ${OUT}/STAGE_JACCARD1_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 1)
	#JOB_ID1=$(qsub -pe shmem 10 -q long.qc -t 1 -N STAGE_JACCARD1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -terse -hold_jid ${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 1)
	JOBA_ID2=$(sbatch --parsable -p short -d afterok:${JOB5_ID} --cpus-per-task=10 -J STAGE_J2 -e ${OUT}/STAGE_JACCARD2_%a.e -o ${OUT}/STAGE_JACCARD2_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 2)
	#JOB_ID2=$(qsub -pe shmem 10 -q long.qc -t 1 -N STAGE_JACCARD2 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -terse -hold_jid ${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 2)
	exit 0
fi 

###----------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------
## OPTION R:     RERUN QC AND JACCARD - FOR EXAMPLE IF YOU MADE A MISTAKE IN YOUR LAYOUTS FILE 

if [[ "$STAGE" == "R" ]]; then 
	echo "RUNNING STAGES RS and Jaccard"
	
	# SETUP: Make a directory specific for this sample RUN -> this will be the location of all the LOG Files
	rm -r COMMANDLOGS/RS
	rm -r COMMANDLOGS/JI
	
	# RUN STAGE RS: uses only 'one' task 
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_RS.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/RS/
	mkdir ${OUT}
	#Part five uses only 'one' and runs the R script 
	JOBA_ID=$(sbatch --parsable -p short --cpus-per-task=8 -J STAGE_RS -e ${OUT}/STAGE_RS_%a.e -o ${OUT}/STAGE_RS_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE} ${TECHNICAL_SAMPLES})
	#JOB_ID=$(qsub -pe shmem 8 -q short.qc -t 1 -N STAGE_RS -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE} ${TECHNICAL_SAMPLES})
	echo "${SAMPLES_FILE_POST},TASK.RS,ARRAYSIZE:1,JOBID:${JOB_ID},RUNNAME:${RUNNAME},BATCH_FILE:${BATCH_FILE}"

	# RUN STAGE JACCARD: uses only 'one' task 
	echo "RUNNING STAGES JACCARD"	
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_JACCARD.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/JI/
	mkdir ${OUT}
	JOBA_ID1=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=10 -J STAGE_J1 -e ${OUT}/STAGE_JACCARD1_%a.e -o ${OUT}/STAGE_JACCARD1_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 1)
	#JOB_ID1=$(qsub -pe shmem 10 -q long.qc -t 1 -N STAGE_JACCARD1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -terse -hold_jid ${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 1)
	JOBA_ID2=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=10 -J STAGE_J2 -e ${OUT}/STAGE_JACCARD2_%a.e -o ${OUT}/STAGE_JACCARD2_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 2)
	#JOB_ID2=$(qsub -pe shmem 10 -q long.qc -t 1 -N STAGE_JACCARD2 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -terse -hold_jid ${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 2)
	exit 0
fi 



###----------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------
## OPTION 3:     RUN BASIC PREPROCESSING AND R QC ANALYSIS + JACCARD (sample overlap) + CONSENSUS JACCARD (UMI OVERLAP)
##               AT THE END OF THIS YOU WILL NEED TO UPLOAD THE FULLY REDUCED FILE TO IMGT FOR ANNOTATION 
##               YOU WILL ALSO NEED TO SET SUBSAMPLE DEPTH FOR ISOTYPER 

if [[ "$STAGE" -eq 3 ]]; then 
	echo "RUNNING STAGES 1-5, RS and CONSENSUS JACCARD"
	
	# SETUP: Make a directory specific for this sample RUN -> this will be the location of all the LOG Files
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}
	mkdir COMMANDLOGS/${SAMPLES_FILE_POST}
	CMD="rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_*"
	eval "${CMD}"
	CMD="rm COMMANDLOGS/job_${SAMPLES_FILE_PRE}_*"
	eval "${CMD}"
	
	# TASK NUMBER: CALCULATE NUMBER OF LIBRARIES 
	TASKS=$(cat ${SAMPLES_FILE_PRE} | wc -l)
	TASKS=$((TASKS+1))

	# RUN STAGE 1
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_1.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/1/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short --cpus-per-task=1 -J STAGE_1 -e ${OUT}/STAGE_1_%a.e -o ${OUT}/STAGE_1_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_PRE} 1 ${SAMPLES_FILE_POST})
	echo "${SAMPLES_FILE_PRE},TASK1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	# TASK NUMBER: CALCULATE NUMBER OF DEMULTIPLEXED SAMPLES
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))
	
	# RUN STAGE 2
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_2.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/2/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_2 -e ${OUT}/STAGE_2_%a.e -o ${OUT}/STAGE_2_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	#JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N STAGE_2 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/2/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/2/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 2 ${SAMPLES_FILE_PRE})
	echo "${SAMPLES_FILE_POST},TASK2,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	# RUN STAGE 3
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_3.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/3/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_3 -e ${OUT}/STAGE_3_%a.e -o ${OUT}/STAGE_3_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 3)
	#JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N STAGE_3 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/3/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/3/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 3)
	echo "${SAMPLES_FILE_POST},TASK3,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	# RUN STAGE 4
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_4.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/4/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_4 -e ${OUT}/STAGE_4_%a.e -o ${OUT}/STAGE_4_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 4)
	#JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N STAGE_4 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/4/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/4/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 4)
	#JOB_ID=${JOB_ID%.*}
	echo "${SAMPLES_FILE_POST},TASK4,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	# RUN STAGE 5
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_5.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/5/
	mkdir ${OUT}
	#Part five uses only 'one task e.g. not an array'
	JOB5_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_5 -e ${OUT}/STAGE_5_%a.e -o ${OUT}/STAGE_5_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 5)
	#JOB_ID=$(qsub -pe shmem 1 -q short.qc -t 1 -N STAGE_5 -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/5/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/5/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 5)
	echo "${SAMPLES_FILE_POST},TASK5,ARRAYSIZE:1,JOBID:${JOB_ID}"

	# RUN STAGE RS
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_RS.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/RS/
	mkdir ${OUT}
	#Part five uses only 'one' and runs the R script 
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOB5_ID} --cpus-per-task=8 -J STAGE_RS -e ${OUT}/STAGE_RS_%a.e -o ${OUT}/STAGE_RS_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE} ${TECHNICAL_SAMPLES})
	#JOB_ID=$(qsub -pe shmem 8 -q short.qc -t 1 -N STAGE_RS -hold_jid ${JOB_ID} -e COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/RS/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} RS ${RUNNAME} ${BATCH_FILE} ${TECHNICAL_SAMPLES})
	echo "${SAMPLES_FILE_POST},TASK.RS,ARRAYSIZE:1,JOBID:${JOB_ID},RUNNAME:${RUNNAME},BATCH_FILE:${BATCH_FILE}"

	# RUN STAGE  JACCARD
	# Note this will be submitted in parallel to RS analysis 
	echo "RUNNING STAGES JACCARD "
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_JACCARD.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/JI/
	mkdir ${OUT}
	JOBA_ID1=$(sbatch --parsable -p short -d afterok:${JOB5_ID} --cpus-per-task=10 -J STAGE_J1 -e ${OUT}/STAGE_JACCARD1_%a.e -o ${OUT}/STAGE_JACCARD1_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 1)
	#JOB_ID1=$(qsub -pe shmem 10 -q long.qc -t 1 -N STAGE_JACCARD1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -terse -hold_jid ${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 1)
	JOBA_ID2=$(sbatch --parsable -p short -d afterok:${JOB5_ID} --cpus-per-task=10 -J STAGE_J2 -e ${OUT}/STAGE_JACCARD2_%a.e -o ${OUT}/STAGE_JACCARD2_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 2)
	#JOB_ID2=$(qsub -pe shmem 10 -q long.qc -t 1 -N STAGE_JACCARD2 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -terse -hold_jid ${JOB_ID5} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 2)
	
	# RUN STAGE  CONSENSUS
	echo "RUNNING STAGES UMI CONSENSUS"
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_CONSENSUS.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/CONSENSUS/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOB5_ID} --cpus-per-task=1 -J STAGE_CO -e ${OUT}/STAGE_CONSENSUS_%a.e -o ${OUT}/STAGE_CONSENSUS_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} CONSENSUS)
	#JOB_ID=$(qsub -pe shmem 1 -q long.qc -t 1-${TASKS} -N STAGE_CONSENSUS -hold_jid ${JOB_ID5} -e COMMANDLOGS/${SAMPLES_FILE_POST}/CONSENSUS/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/CONSENSES/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} CONSENSUS)
	echo "${SAMPLES_FILE_POST},TASK.CONSENSUS,ARRAYSIZE:${TASKS},JOBID:${JOB_ID}"

	# RUN STAGE  JACCARD CONSENSUS 
	# Note this needs to wait until RUNNING STAGE UMI CONSENUS IS COMPLETE 
	echo "RUNNING STAGES UMI CONSENSUS"
	rm COMMANDLOGS/job_${SAMPLES_FILE_POST}_CONSENSUS.txt
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/JI/
    # Running the JACCARD with UMI correction 
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=10 -J STAGE_CJ -e ${OUT}/STAGE_CONSENSUS_JI_%a.e -o ${OUT}/STAGE_CONSENSUS_JI_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 3)
	#JOB_ID3=$(qsub -pe shmem 10 -q long.qc -t 1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -N STAGE_JACCARD3 -o COMMANDLOGS/${SAMPLES_FILE_POST}/JI/ -terse -hold_jid ${JOB_ID} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} JACCARD ${RUNNAME} ${BATCH_FILE} 3)
	echo "${SAMPLES_FILE_PRE},TASK.JACCARD1-3,ARRAYSIZE:1,JOBID:${JOB_ID}"
	exit 0
fi 

###----------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------
## OPTION 4:     (POST IMGT) UNZIP IMGT OUTPUTFILES, RUN ISOTYPER SCRIPT AND COMPILE MATRIX 


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
	
	# RUN STAGE 6
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/6/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short --cpus-per-task=6 -J STAGE_6 -e ${OUT}/STAGE_6_%a.e -o ${OUT}/STAGE_6_%a.o --array 1 BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 6)
	#JOB_ID=$(qsub -pe shmem 6 -q short.qc -t 1 -N STAGE_6 -e COMMANDLOGS/${SAMPLES_FILE_POST}/6/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/6/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} 6)
	echo "${SAMPLES_FILE_POST},TASK6,ARRAYSIZE:1,JOBID:${JOB_ID}"
		
	# TASK NUMBER: CALCULATE NUMBER OF DEMULTIPLEXED SAMPLES
	TASKS=$(cat ${SAMPLES_FILE_POST} | wc -l)
	TASKS=$((TASKS+1))

	# RUN STAGE: ISO1_PRODUCTIVE_NONPRODUCTIVE sequentially (requires fewer slots to become availible)
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/ISO_COMPLETE/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=1 -J STAGE_ISO1 -e ${OUT}/STAGES_ISO1_%a.e -o ${OUT}/STAGE_ISO1_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ISO1_COMPLETE)
	#JOB_ID1=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N STAGE_ISO1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ISO_COMPLETE/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ISO_COMPLETE/ -terse -hold_jid ${JOB_ID} BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ISO1_COMPLETE)
	echo "${SAMPLES_FILE_POST},TASK.ISO1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID1}"
	JOB_ID1=${JOB_ID1%.*}
	
	# RUN STAGE CAT 
	# Needs more memory if running BCR IMGT 
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=5 -J STAGE_CAT -e ${OUT}/STAGE_CAT_%a.e -o ${OUT}/STAGE_CAT_%a.o --array 1 cat_IMGT_files_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ${BATCH_FILE} ${CODE_DIRECTORY} ${IMGT_MUTATION} ${SAMPLES_EXCLUDE} ${NORMALISE_VGENEUSAGE})
	#JOB_ID4=$(qsub -pe shmem 10 -q short.qc -t 1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/ -N STAGE_CAT -o COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/ -terse -hold_jid ${JOB_ID1} cat_IMGT_files.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ${BATCH_FILE} ${CODE_DIRECTORY} YES ${SAMPLES_EXCLUDE})
	echo "${SAMPLES_FILE_POST},TASKCAT,ARRAYSIZE:1,JOBID:${JOB_ID4}"
	exit 0
fi 

###----------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------
## OPTION 5:     (POST IMGT, POST FILES BEING UNZIPPED) RUN ISOTYPER SCRIPT AND COMPILE MATRIX 
##                USEFUL IF YOU WANT TO RERUN THE ISOTYPER SCRIPT WITH A DIFFERENT READ DEPTH 
##                READ DEPTH WILL BE INCLUDED IN SUMMARY FILE NAME, SO MAKE SURE NOT TO GET CONFUSED BETWEEN MULTIPLE RUNS / DELETE THE EXTRA FILES BEFOREHAND

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
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/ISO_COMPLETE/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short --cpus-per-task=1 -J STAGE_ISO1 -e ${OUT}/STAGE_ISO1_%a.e -o ${OUT}/STAGE_ISO1_%a.o --array 1-${TASKS} BCR_TCR_Wrapper_Cluster_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ISO1_COMPLETE)
	#JOB_ID1=$(qsub -pe shmem 1 -q short.qc -t 1-${TASKS} -N STAGE_ISO1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/ISO_COMPLETE/ -o COMMANDLOGS/${SAMPLES_FILE_POST}/ISO_COMPLETE/ -terse BCR_TCR_Wrapper_Cluster.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ISO1_COMPLETE)
	echo "${SAMPLES_FILE_POST},TASK.ISO1,ARRAYSIZE:${TASKS},JOBID:${JOB_ID1}"

	# RUN CAT 
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short -d afterok:${JOBA_ID} --cpus-per-task=5 -J STAGE_CAT -e ${OUT}/STAGE_CAT_%a.e -o ${OUT}/STAGE_CONSENSUS_CAT_%a.o --array 1 cat_IMGT_files_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ${BATCH_FILE} ${CODE_DIRECTORY} ${IMGT_MUTATION} ${SAMPLES_EXCLUDE} ${NORMALISE_VGENEUSAGE})
	#JOB_ID4=$(qsub -pe shmem 10 -q short.qc -t 1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/ -N STAGE_CAT -o COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/ -terse -hold_jid ${JOB_ID1} cat_IMGT_files.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ${BATCH_FILE} ${CODE_DIRECTORY} YES ${SAMPLES_EXCLUDE})
	echo "${SAMPLES_FILE_POST},TASKCAT,ARRAYSIZE:1,JOBID:${JOB_ID4}"
	exit 0
fi 

###----------------------------------------------------------------------------------------
###----------------------------------------------------------------------------------------
## OPTION 5:     (POST IMGT, POST FILES BEING UNZIPPED, POST ISOTYPER) COMPILE MATRIX 
##                USEFUL IF YOU WANT TO RERUN WITH DIFFERENT FILTERING THRESHOLDS AKAK MORE THAN x MISSINGNESS 
##                RECOMMENDED TO JUST USE MY DEFAULT PARAMATERS
##                WARNING USING A SMALL NUMBER OF SAMPLES WE MAY NEED TO CUSTOMISE SCRIPT SO IT DOESNT REMOVE FEATURES WITH LESS THAN N UNIQUE READS (DEFAULT IS 8)
##                ^ IF YOU HAVE 7 SAMPLES THIS WOULD CAUSE A PROBLEM  

if [[ "$STAGE" -eq 6 ]]; then 
	echo "RUNNING ISOTYPER SUMMARY"
	rm -r COMMANDLOGS/${SAMPLES_FILE_POST}/CAT
	
	# RUN CAT
	OUT=COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/
	mkdir ${OUT}
	JOBA_ID=$(sbatch --parsable -p short --cpus-per-task=5 -J STAGE_CAT -e ${OUT}/STAGE_CAT_%a.e -o ${OUT}/STAGE_CAT_%a.o --array 1 cat_IMGT_files_SLURM.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ${BATCH_FILE} ${CODE_DIRECTORY} ${IMGT_MUTATION} ${SAMPLES_EXCLUDE} ${NORMALISE_VGENEUSAGE})
	#JOB_ID=$(qsub -pe shmem 10 -q long.qc -t 1 -e COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/ -N STAGE_CAT -o COMMANDLOGS/${SAMPLES_FILE_POST}/CAT/ -terse BCR_TCR_Wrapper_2.sh ${DEPENDANCIES} ${SAMPLES_FILE_POST} ${BATCH_FILE} ${CODE_DIRECTORY} YES ${SAMPLES_EXCLUDE})
	echo "${SAMPLES_FILE_POST},TASKCAT,ARRAYSIZE:1,JOBID:${JOB_ID}"
	exit 0

fi 

#########################################################################################################################
######################################################################################################################
## ALL PREPROCESSING STAGES DONE!!!!!!
## Running Tensor/ Module Reduction
