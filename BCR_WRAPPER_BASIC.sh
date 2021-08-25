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

# Load software modules
module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load networkx/2.2-foss-2019b-Python-2.7.16


# Job Arguments
# File one containing the libraries
# File two containing the PCR multiplexed samples 
SAMPLES_FILE_POST=$1
TASK=$2
## Runname or in the case of part 2 and 1  it will be the samples file pre/post 
RUNNAME=$3
BATCH_FILE=$4
JACCARD_TASK=$5

#------------------------------------

## Set up module environement for R scripts
if [[ "$TASK" == "RS" || "$TASK" == "JACCARD" ]]; then
module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
fi


# Task Arguments for TASK 1
ID=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$1}" $SAMPLES_FILE_POST) 
SAMPLE=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$2}" $SAMPLES_FILE_POST)      # Sample ID
INFO=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$3}" $SAMPLES_FILE_POST)     # 10X FastQ sample files path
GENE=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$4}" $SAMPLES_FILE_POST)    # Output folder path
FASTQDIRECTORY=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$5}" $SAMPLES_FILE_POST)
INITIALPAIRING=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$6}" $SAMPLES_FILE_POST)
FINALPAIRING=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$7}" $SAMPLES_FILE_POST)
OUTPUTDIR=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$8}" $SAMPLES_FILE_POST)
PLATFORM=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$9}" $SAMPLES_FILE_POST)
SPECIES=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$10}" $SAMPLES_FILE_POST)
CONSTANTPRIMER=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$11}" $SAMPLES_FILE_POST)
TYPE=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$12}" $SAMPLES_FILE_POST)
BARCODE=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$13}" $SAMPLES_FILE_POST)
RECEPTOR=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$14}" $SAMPLES_FILE_POST)
ORF_FILTER=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$15}" $SAMPLES_FILE_POST)
FILE_LENGTH=$(awk '{print NF}' $SAMPLES_FILE_POST | sort -nu | tail -n 1)

OTHER="${TYPE},${BARCODE},${RECEPTOR}"

## If ORF Filter not provided (when full 14 other columns provided - defaults to TRUE 
if [[ -z "$ORF_FILTER" && "$FILE_LENGTH" -eq 14 ]] ;then 
	ORF_FILTER="True"
fi 

if [[ "$FILE_LENGTH" -lt 14 && "$FILE_LENGTH" -gt 11 ]]; then 
	ORF_FILTER=$(awk  -F '\t' "{if (NR==$SGE_TASK_ID) print \$NF}" $SAMPLES_FILE_POST)
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
echo "OTHER           : ${OTHER}"
echo "RECEPTOR        : ${RECEPTOR}"
echo "ORF Filtering   : ${ORF_FILTER}"
echo "Run Name        : ${RUNNAME}"
echo "********************************************************"



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
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Get_batch_information.py $1"
elif [[ "$TASK" == 6 ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Combine_extract_IMGT_information.py $1 ${OUTPUTDIR}"
elif [[ "$TASK" == "RS" ]]; then
CMD="Rscript AnalysisStages1to4.R -o ${OUTPUTDIR} -r ${RUNNAME} -g ${GENE} -b ${BATCH_FILE}"
elif [[ "$TASK" == "JACCARD" ]]; then
CMD="Rscript AnalysisJaccard.R -o ${OUTPUTDIR} -r ${RUNNAME} -g ${GENE} -b ${BATCH_FILE} -t ${JACCARD_TASK}"
elif [[ "$TASK" == "ISO1" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/IsoTyper_2.0.py ${ID} ${ID} ${OUTPUTDIR} ${SPECIES} ${RECEPTOR} $1"
CMD2="cp ${OUTPUTDIR}ORIENTATED_SEQUENCES/ANNOTATIONS/Sampling_depth_per_isotype_${SAMPLES_FILE_POST} ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/"
elif [[ "$TASK" == "ISO1" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/IsoTyper_2.0.py ${ID} ${ID} ${OUTPUTDIR} ${SPECIES} ${RECEPTOR} $1"
CMD2="cp ${OUTPUTDIR}ORIENTATED_SEQUENCES/ANNOTATIONS/Sampling_depth_per_isotype_${SAMPLES_FILE_POST} ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/"
elif [[ "$TASK" == "ISO1_PRODUCTIVE" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/IsoTyper_2.0.py ${ID}_productive ${ID}_productive ${OUTPUTDIR} ${SPECIES} ${RECEPTOR} $1"
CMD2="cp ${OUTPUTDIR}ORIENTATED_SEQUENCES/ANNOTATIONS/Sampling_depth_per_isotype_${SAMPLES_FILE_POST} ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/"
elif [[ "$TASK" == "ISO1_NON_PRODUCTIVE" ]]; then
CMD="python /well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/IsoTyper_2.0.py ${ID}_unproductive ${ID}_unproductive ${OUTPUTDIR} ${SPECIES} ${RECEPTOR} $1"
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

 
# PRINT JOB TO LOGFILE
echo "********************************************************"
echo "["`date`"] Running TCR/BCR"
echo "********************************************************"
echo "Command : ${CMD}"
echo "********************************************************"
echo

## RUN THE JOB 1
eval "${CMD}"


# Move the file if running isotyper 
if [[ "$TASK" == "ISO1" || "$TASK" == "ISO1_PRODUCTIVE" || "$TASK" == "ISO1_NON_PRODUCTIVE" ]]; then 
	eval "${CMD2}"
fi

# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0