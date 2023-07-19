#!/bin/bash
#$ -cwd

# If there's an error, fail the whole script
#set -e -o pipefail

# Set permissions so any user in the group can 
# read/write what it's created by the script
umask 002

echo `date`: Executing task ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID} on `hostname` as user ${USER} 
echo "********************************************************"
# Load software modules
module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load networkx/2.2-foss-2019b-Python-2.7.16
echo "Loaded networkx/2.2-foss-2019b-Python-2.7.16 Module"

# Job Arguments
# File one containing the libraries
# File two containing the PCR multiplexed samples 
DEPENDANCIES=$1
SAMPLES_FILE_POST=$2
TASK=$3
## Runname or in the case of part 2 it will be the samples file pre 
RUNNAME=$4
BATCH_FILE=$5
JACCARD_TASK=$6
TECHNICAL_SAMPLES=$6
CPU_TYPE=$MODULE_CPU_TYPE

# "catch exit status 1" grep wrapper
# exit status 1 when no lines matched  - this is causing the script to fail with set -e hence the catch error
# exit status 0 when lines matched
# exit status 2 when error
c1grep() { grep "$@" || test $? = 1; }

echo "Running SLURM WRAPPER with ERROR CATCH"

## Set up module environement for R scripts
if [[ "$TASK" == "RS" || "$TASK" == "JACCARD" ]]; then
module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0

if [[ "$CPU_TYPE" == "skylake" ]]; then
echo "CPU TYPE: ${CPU_TYPE}"
echo "exporting R packages from /well/immune-rep/shared/CODE/Communal_R_packages/4.0/skylake"
export R_LIBS=/well/immune-rep/shared/CODE/Communal_R_packages/4.0/skylake
fi 

if [[ "$CPU_TYPE" == "ivybridge" ]]; then
echo "CPU TYPE: ${CPU_TYPE}"
echo "exporting R packages from /well/immune-rep/shared/CODE/Communal_R_packages/4.0/ivybridge"
export R_LIBS=/well/immune-rep/shared/CODE/Communal_R_packages/4.0/ivybridge
fi 

echo "Loaded R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0 Module"
fi

echo "********************************************************"
echo "********************************************************"
echo "* Job Dependancies"
echo "********************************************************"

#------------------------------------------------------------
## Determining Pior Task Dependancies: 
## For numerical tasks: 
PRIORTASK=$((TASK-1))

## For named tasks:  
if [[ "$TASK" == "RS" ]]; then 
	PRIORTASK=5
fi
if [[ "$TASK" == "ISO1"  || "$TASK" == "TCRISO1" || "$TASK" == "ISO1_NON_PRODUCTIVE" || "$TASK" == "ISO1_PRODUCTIVE" || "$TASK" == "ISO1_COMPLETE" ]]; then 
	PRIORTASK=6
fi
if [[ "$TASK" == "CONSENSUS" ]]; then 
	PRIORTASK=5
fi

if [[ "$TASK" == "JACCARD" ]]; then 
	PRIORTASK="CONSENSUS"
fi

if [[ "$JACCARD_TASK" == "1" || "$JACCARD_TASK" == "2" ]]; then 
	PRIORTASK=5
fi

if [ -z "$4" ]
  then
    echo "No Batch file Supplied. Defaulting to FALSE"
	BATCH_FILE="FALSE"	
fi

echo ${PRIORTASK}

#------------------------------------------------------------
## Determining which sample file (pre/post was used) and file length. 
## Stages using 'post-file' (2 onwards and Consensus)
if [[ "$DEPENDANCIES" == "YES" || "$DEPENDANCIES" == "Y" || "$DEPENDANCIES" == "y" || "$DEPENDANCIES" == "yes" ]]; then
if [[ "$PRIORTASK" -ge 2 || "$PRIORTASK" == "-1" || "$PRIORTASK" == "CONSENSUS" ]]; then 
FILE="COMMANDLOGS/job_${SAMPLES_FILE_POST}_${PRIORTASK}.txt"
echo ${FILE}
	if [ ! -e "$FILE" ]; then
		echo "DEPENDANCIES: Final File from previous job does NOT exist - something has gone wrong and no samples have run!"
		exit 888
	else 
		echo "DEPENDANCIES: Final File from previous job exists: ${FILE}"
	fi 
LENGTHJOBS=$(cat ${FILE} | wc -l)
SAMPLES=$(cat ${SAMPLES_FILE_POST} | wc -l)
SAMPLES=$((SAMPLES+1))
fi 

## Stages using 'pre-file' (stage 1)
if [[ "$PRIORTASK" == "1" ]]; then
FILE="COMMANDLOGS/job_${RUNNAME}_${PRIORTASK}.txt"
echo ${FILE}
	if [ ! -e "$FILE" ]; then
		echo "DEPENDANCIES: Final File from previous job does NOT exist - something has gone wrong and no samples have run!"
		exit 888
	else 
		echo "DEPENDANCIES: Final File from previous job exists: ${FILE}"
	fi 
LENGTHJOBS=$(cat ${FILE} | wc -l)
SAMPLES=$(cat ${RUNNAME} | wc -l)
SAMPLES=$((SAMPLES+1))
fi 
fi 

if [[ "$DEPENDANCIES" == "NO" || "$DEPENDANCIES" == "N" || "$DEPENDANCIES" == "n" || "$DEPENDANCIES" == "no" ]]; then 
echo "CHECK FOR DEPENDANCIES HAD BEEN TURNED OFF BY USER"
fi  
#------------------------------------------------------------
## Checking whether prior jobs ran sucessfully - if not terminate script. 
## Print out Failed Samples to a failure file to easily identify.

# Task Arguments for checks
ID=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$1}" $SAMPLES_FILE_POST) 
# All Samples for array job
SAMPLE=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$2}" $SAMPLES_FILE_POST) 
# All IDs for array job
IDS=$(awk -F '\t' "{ print \$1 }" $SAMPLES_FILE_POST)

## Samples and ID files
SAMPLE_FILE=COMMANDLOGS/${SAMPLES_FILE_POST}_SAMPLES.txt
IDS_FILE=COMMANDLOGS/${SAMPLES_FILE_POST}_IDS.txt

if [[ "$DEPENDANCIES" == "YES" || "$DEPENDANCIES" == "Y" ]]; then 
if [[ "$PRIORTASK" -lt 5 &&  "$PRIORTASK" -ge 1 || "$PRIORTASK" == "-1" ||"$PRIORTASK" == "CONSENSUS" ]]; then
	if [[ "$PRIORTASK" -ge 2 && "$PRIORTASK" -lt 4 ]]; then 
		if fgrep -qx "$ID" $FILE
		then
			echo "ID ${ID} PRESENT in Outs file of PRIOR TASK ${PRIORTASK}"
		else 
			echo "Warning ${ID} is NOT Present in OUTS file of PRIOR TASK ${PRIORTASK} - check sample naming"
			#FAILED="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_${PRIORTASK}_FAILED_SAMPLES.txt"
			#eval "${FAILED}" 
		fi 
	fi 
	
	if [[ "$PRIORTASK" -eq 4 ]]; then
		echo "Looking for failed samples.." 
		#c1grep -v -f ${FILE} ${IDS_FILE} > COMMANDLOGS/job_${SAMPLES_FILE_POST}_${PRIORTASK}_FAILED_SAMPLES.txt
		echo "DONE"
	fi 
		
	if [[ "$PRIORTASK" -eq 1 ]]; then 
		if fgrep -qx "$SAMPLE" $FILE
		then
			echo "ID ${SAMPLE} PRESENT in Outs file of PRIOR TASK ${PRIORTASK}"
		#else 
			#FAILED="echo ${SAMPLE} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_${PRIORTASK}_FAILED_SAMPLES.txt"
			#eval "${FAILED}" 
		fi 
	fi 
		
	if [[ $LENGTHJOBS -ne $SAMPLES ]]; then 
		echo "WARNING: DEPENDANCIES: Not all dependancies ran sucessfully to completion"
		echo "ERROR: NO FURTHER ANALYSIS WILL BE PERFORMED"
		exit 999
	else 
		echo "DEPENDANCIES: All dependancies ran sucessfully to completion"
		echo "SUCCESS: ANALYSIS STAGE ${TASK} WILL BE PERFORMED"
	fi
else   
    if [[ "$PRIORTASK" -eq 0 ]]; then
		echo "DEPENDANCIES: NO DEPENDANCIES REQUIRED"
		echo "SUCCESS: ANALYSIS STAGE ${TASK} WILL BE PERFORMED"
    fi 
	 if [[ "$PRIORTASK" -eq 6 && $LENGTHJOBS -eq "1" ]]; then
		echo "DEPENDANCIES: STAGE RAN SUCESSFULLY TO COMPLETION (Post-IMGT)"
		echo "SUCCESS: ANALYSIS STAGE ${TASK} WILL BE PERFORMED"
    fi 
    if [[ "$PRIORTASK" -eq 5 && $LENGTHJOBS -ne "1" ]]; then
		echo "DEPENDANCIES: STAGE 5 DID not run sucessfully to completion"
		echo "ERROR: NO FURTHER ANALYSIS WILL BE PERFORMED"
		exit 999
    fi 
	   if [[ "$PRIORTASK" -eq 6 && $LENGTHJOBS -ne "1" ]]; then
		echo "DEPENDANCIES: STAGE 6 DID not run sucessfully to completion (Post-IMGT)"
		echo "ERROR: NO FURTHER ANALYSIS WILL BE PERFORMED"
		exit 999
    fi 
	if [[ "$PRIORTASK" -eq 5 && $LENGTHJOBS -eq "1" ]]; then
		echo "DEPENDANCIES: STAGE 5 ran sucessfully to completion"
		echo "SUCCESS: ANALYSIS STAGE ${TASK} WILL BE PERFORMED"
	fi
fi 
fi


#------------------------------------
# Task Arguments for TASK 1
ID=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$1}" $SAMPLES_FILE_POST) 
SAMPLE=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$2}" $SAMPLES_FILE_POST)      # Sample ID
INFO=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$3}" $SAMPLES_FILE_POST)     # 10X FastQ sample files path
GENE=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$4}" $SAMPLES_FILE_POST)    # Output folder path
FASTQDIRECTORY=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$5}" $SAMPLES_FILE_POST)
INITIALPAIRING=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$6}" $SAMPLES_FILE_POST)
FINALPAIRING=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$7}" $SAMPLES_FILE_POST)
OUTPUTDIR=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$8}" $SAMPLES_FILE_POST)
PLATFORM=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$9}" $SAMPLES_FILE_POST)
SPECIES=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$10}" $SAMPLES_FILE_POST)
CONSTANTPRIMER=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$11}" $SAMPLES_FILE_POST)
TYPE=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$12}" $SAMPLES_FILE_POST)
BARCODE=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$13}" $SAMPLES_FILE_POST)
RECEPTOR=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$14}" $SAMPLES_FILE_POST)
ORF_FILTER=$(awk -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$15}" $SAMPLES_FILE_POST)
FILE_LENGTH=$(awk '{print NF}' $SAMPLES_FILE_POST | sort -nu | tail -n 1)

OTHER="${TYPE},${BARCODE},${RECEPTOR}"

## If ORF Filter not provided (when full 14 other columns provided - defaults to TRUE 
if [[ -z "$ORF_FILTER" && "$FILE_LENGTH" -eq 14 ]] ;then 
	ORF_FILTER="True"
fi 

if [[ "$FILE_LENGTH" -lt 14 && "$FILE_LENGTH" -gt 11 ]]; then 
	ORF_FILTER=$(awk  -F '\t' "{if (NR==$SLURM_ARRAY_TASK_ID) print \$NF}" $SAMPLES_FILE_POST)
	RECEPTOR="STANDARD"
fi 

if [[ ! "$ORF_FILTER" =~ (True|ORF) ]]; then 
	ORF_FILTER="True"
fi 

if [[ "$ORF_FILTER" == "TRUE" ]]; then 
	ORF_FILTER="True"
fi 

if [[ "$FILE_LENGTH" -le 12 ]]; then 
	OTHER=""
fi 
	

echo "********************************************************"
echo "* Job Details"
echo "********************************************************"
echo "SLURM job ID     : "$SLURM_ARRAY_JOB_ID
echo "SLURM task ID    : "$SLURM_ARRAY_TASK_ID
echo "Run on host      : "`hostname`
echo "Operating system : "`uname -s`
echo "Username         : "`whoami`
echo "Started at       : "`date`
echo "RUNNING BCR_TCR_PIPELINE STAGE: "$TASK 
echo "Authors: Rachael Bashford-Rogers and Lauren Overend"
echo "********************************************************"
echo "* Job Parameters"
echo "********************************************************"
echo "No. Inputs      : ${FILE_LENGTH}"
echo "ID     		  : ${ID}"
echo "SAMPLE          : ${SAMPLE}"
echo "INFO     		  : ${INFO}"
echo "GENE            : ${GENE}"
echo "FASTQDIRECTORY  : ${FASTQDIRECTORY}"
echo "INITIALPAIRING  : ${INITIALPAIRING}"
echo "FINALPAIRING    : ${FINALPAIRING}"
echo "OUTPUTDIR       : ${OUTPUTDIR}"
echo "PLATFORM        : ${PLATFORM}"
echo "SPECIES         : ${SPECIES}"
echo "CONSTANTPRIMER  : ${CONSTANTPRIMER}"
echo "TYPE  		  : ${TYPE}"
echo "BARCODE         : ${BARCODE}"
echo "RECEPTOR        : ${RECEPTOR}"
echo "OTHER           : ${OTHER}"
echo "ORF Filtering   : ${ORF_FILTER}"
echo "Run Name        : ${RUNNAME}"
echo "Batch File      : ${BATCH_FILE}"
echo "Sample File     : ${SAMPLE_FILE}"
echo "IDs File        : ${IDS_FILE}"
echo "********************************************************"


# PRINT JOB 1
if [[ "$TASK" == 1 ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Read_processing_and_quality.py ${OUTPUTDIR} ${ID} ${SAMPLE} ${GENE} ${INITIALPAIRING} ${SPECIES} ${FASTQDIRECTORY} 200 ${CONSTANTPRIMER} ${PLATFORM} 1 $OTHER ${RECEPTOR}"
elif [[ "$TASK" == 2 ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Read_processing_and_quality.py ${OUTPUTDIR} ${ID} ${SAMPLE} ${GENE} ${INITIALPAIRING} ${SPECIES} ${FASTQDIRECTORY} 200 ${CONSTANTPRIMER} ${PLATFORM} 2 $OTHER ${RECEPTOR} ${ORF_FILTER}"
elif [[ "$TASK" == "2UJ" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Read_processing_and_quality_unjoined.py ${OUTPUTDIR} ${ID} ${SAMPLE} ${GENE} ${INITIALPAIRING} ${SPECIES} ${FASTQDIRECTORY} 200 ${CONSTANTPRIMER} ${PLATFORM} 2 $OTHER ${RECEPTOR}"
elif [[ "$TASK" == 3 ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Read_processing_and_quality.py ${OUTPUTDIR} ${ID} ${SAMPLE} ${GENE} ${INITIALPAIRING} ${SPECIES} ${FASTQDIRECTORY} 200 ${CONSTANTPRIMER} ${PLATFORM} 3 $OTHER ${RECEPTOR}"
elif [[ "$TASK" == 4 ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Generate_repertoire_statistics.py ${OUTPUTDIR}ORIENTATED_SEQUENCES/ANNOTATIONS/ ${ID} ${OUTPUTDIR}ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_${ID}.fasta ${OUTPUTDIR}ORIENTATED_SEQUENCES/Filtered_ORFs_sequences_all_${ID}.fasta ${GENE} ${SPECIES} ${OUTPUTDIR}ORIENTATED_SEQUENCES/NETWORKS/Cluster_identities_${ID}.txt ANNOTATE,STATISTICS ${RECEPTOR}"
elif [[ "$TASK" == 5 ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Get_batch_information.py $2"
elif [[ "$TASK" == 6 ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Combine_extract_IMGT_information.py $2 ${OUTPUTDIR}"
elif [[ "$TASK" == "RS" ]]; then
CMD="Rscript AnalysisStages1to4.R -o ${OUTPUTDIR} -r ${RUNNAME} -g ${GENE} -b ${BATCH_FILE} -t ${TECHNICAL_SAMPLES}"
elif [[ "$TASK" == "JACCARD" ]]; then
CMD="Rscript AnalysisJaccard.R -o ${OUTPUTDIR} -r ${RUNNAME} -g ${GENE} -b ${BATCH_FILE} -t ${JACCARD_TASK}"
elif [[ "$TASK" == "ISO1" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/IsoTyper_2.0.py ${ID} ${ID} ${OUTPUTDIR} ${SPECIES} ${RECEPTOR} $2"
CMD2="cp ${OUTPUTDIR}ORIENTATED_SEQUENCES/ANNOTATIONS/Sampling_depth_per_isotype_${SAMPLES_FILE_POST} ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/"
elif [[ "$TASK" == "ISO1_PRODUCTIVE" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/IsoTyper_2.0.py ${ID}_productive ${ID}_productive ${OUTPUTDIR} ${SPECIES} ${RECEPTOR} $2"
CMD2="cp ${OUTPUTDIR}ORIENTATED_SEQUENCES/ANNOTATIONS/Sampling_depth_per_isotype_${SAMPLES_FILE_POST} ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/"
elif [[ "$TASK" == "ISO1_NON_PRODUCTIVE" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/IsoTyper_2.0.py ${ID}_unproductive ${ID}_unproductive ${OUTPUTDIR} ${SPECIES} ${RECEPTOR} $2"
CMD2="cp ${OUTPUTDIR}ORIENTATED_SEQUENCES/ANNOTATIONS/Sampling_depth_per_isotype_${SAMPLES_FILE_POST} ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/"
elif [[ "$TASK" == "ISO1_COMPLETE" ]]; then
CMDA="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/IsoTyper_2.0.py ${ID} ${ID} ${OUTPUTDIR} ${SPECIES} ${RECEPTOR} $2"
CMDB="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/IsoTyper_2.0.py ${ID}_productive ${ID}_productive ${OUTPUTDIR} ${SPECIES} ${RECEPTOR} $2"
CMDC="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/IsoTyper_2.0.py ${ID}_unproductive ${ID}_unproductive ${OUTPUTDIR} ${SPECIES} ${RECEPTOR} $2"
CMD2="cp ${OUTPUTDIR}ORIENTATED_SEQUENCES/ANNOTATIONS/Sampling_depth_per_isotype_${SAMPLES_FILE_POST} ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/"
elif [[ "$TASK" == "NONISO1" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Non_isotyper_1.0.py ${ID} ${ID} ${OUTPUTDIR} ${SPECIES}"
elif [[ "$TASK" == "TCRISO1" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Non_isotyper_1.0.py ${ID} ${ID} ${OUTPUTDIR} ${SPECIES} ${RECEPTOR}"
elif [[ "$TASK" == "CSR" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Class_switch_recombination_analysis.py ${OUTPUTDIR}ORIENTATED_SEQUENCES/ANNOTATIONS/ ${ID} ${OUTPUTDIR}ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_${ID}.fasta ${OUTPUTDIR}ORIENTATED_SEQUENCES/Filtered_ORFs_sequences_all_${ID}.fasta ${GENE} ${SPECIES} ${OUTPUTDIR}ORIENTATED_SEQUENCES/NETWORKS/Cluster_identities_${ID}.txt 1"
elif [[ "$TASK" == "SUBSAMPLE" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Subsampling_networks.py ${OUTPUTDIR}ORIENTATED_SEQUENCES/ANNOTATIONS/ ${ID} ${OUTPUTDIR}ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_${ID}.fasta ${OUTPUTDIR}ORIENTATED_SEQUENCES/NETWORKS/Edges_${ID}.txt" 
elif [[ "$TASK" == "CONSENSUS" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Get_consensus_counts.py ${OUTPUTDIR} ${ID}" 
else
echo "PARAMETER NOT RECOGNISED."
fi 

 #########################################################
## RUN THE JOB 1
if [[ "$TASK" -ne "ISO1_COMPLETE" ]]; then
# PRINT JOB TO LOGFILE
echo "********************************************************"
echo "["`date`"] Running TCR/BCR"
echo "********************************************************"
echo "Command : ${CMD}"
echo "********************************************************"
echo
eval "${CMD}"
status=$?
### Has the sample/file run successfully 
if [ $status -eq 0 ]; then
	echo "********************************************************"
    echo "COMMAND EXECUTED SUCCESSFULLY"
	echo "********************************************************"
	NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_${TASK}.txt"
	eval "${NWCMD}"
else
	echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    echo "ERROR: COMMAND FAILED WITH ERROR STATUS: $status"
	echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
	NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_${TASK}_FAILED_SAMPLESX.txt"
	eval "${NWCMD}"
fi
fi

##### If its the R script
if [[ "$TASK" == "RS" || "$TASK" == "JACCARD" ]]; then
# PRINT JOB TO LOGFILE
echo "********************************************************"
echo "["`date`"] Running TCR/BCR"
echo "********************************************************"
echo "Command : ${CMD}"
echo "********************************************************"
echo
eval "${CMD}"
status=$?
if [ $status -eq 0 ]; then
    echo "********************************************************"
    echo "COMMAND EXECUTED SUCCESSFULLY"
	echo "********************************************************"
	NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_${TASK}.txt"
	eval "${NWCMD}"
else
	echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    echo "ERROR: COMMAND FAILED WITH ERROR STATUS: $status"
	echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
	NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_${TASK}_FAILED_SAMPLESX.txt"
	eval "${NWCMD}"
fi
fi
#########################################################


#########################################################
if [[ "$TASK" == "ISO1_COMPLETE" ]]; then
##ALL
echo "********************************************************"
echo "["`date`"] Running TCR/BCR"
echo "********************************************************"
echo "Command : ${CMDA}"
echo "********************************************************"
eval "${CMDA}" 
status=$?
if [ $status -eq 0 ]; then
	echo "********************************************************"
    echo "COMMAND EXECUTED SUCCESSFULLY"
	echo "********************************************************"
	NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1.txt"
	eval "${NWCMD}"
else
	echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    echo "ERROR: COMMAND FAILED WITH ERROR STATUS: $status"
	echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
	NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_FAILED_SAMPLESX.txt"
	eval "${NWCMD}"
fi
eval "${CMD2}"
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
echo "********************************************************"

## PRODUCTIVE
echo "Command : ${CMDB}"
echo "********************************************************"
eval "${CMDB}" 
status=$?
if [ $status -eq 0 ]; then
	echo "********************************************************"
    echo "COMMAND EXECUTED SUCCESSFULLY"
	echo "********************************************************"
	NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_PRODUCTIVE.txt"
	eval "${NWCMD}"
else
	echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    echo "ERROR: COMMAND FAILED WITH ERROR STATUS: $status"
	echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
	NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_PRODUCTIVE_FAILED_SAMPLESX.txt"
	eval "${NWCMD}"
fi
eval "${CMD2}" 
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"


## NON_PRODUCTIVE
echo "********************************************************"
echo "Command : ${CMDC}"
echo "********************************************************"
eval "${CMDC}" 
status=$?
if [ $status -eq 0 ]; then
	echo "********************************************************"
    echo "COMMAND EXECUTED SUCCESSFULLY"
	echo "********************************************************"
	NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_NON_PRODUCTIVE.txt"
	eval "${NWCMD}"
else
	echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    echo "ERROR: COMMAND FAILED WITH ERROR STATUS: $status"
	echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
	NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISO1_NON_PRODUCTIVE_FAILED_SAMPLESX.txt"
	eval "${NWCMD}"
fi
eval "${CMD2}" 
echo "********************************************************"
echo "["`date`"] Done"
echo "******************************************************"
fi
#########################################################

	
## Make read writeable by everyone 
## Sometimes you get errors if its accessing files that another job has already edited to lets not save the output 
echo "Updating Permissions"
chmod -R g+wrx ${OUTPUTDIR} 2>/dev/null


# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0