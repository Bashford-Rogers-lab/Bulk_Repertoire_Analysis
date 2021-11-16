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

echo "********************************************************"
echo "********************************************************"
echo "* Job Dependancies"
echo "********************************************************"

## Specifying dependancies
PRIORTASK1="ISO1"
PRIORTASK2="ISO1_PRODUCTIVE"
PRIORTASK3="ISO1_NON_PRODUCTIVE"

SAMPLES=$(cat ${SAMPLES_FILE_POST} | wc -l)
SAMPLES=$((SAMPLES+1))

# Testing for sample success!
if [[ "$DEPENDANCIES" == "YES" || "$DEPENDANCIES" == "Y" || "$DEPENDANCIES" == "y" || "$DEPENDANCIES" == "yes" ]]; then
######################
FILE1="COMMANDLOGS/job_${SAMPLES_FILE_POST}_${PRIORTASK1}.txt"
echo ${FILE1}
	if [ ! -e "$FILE1" ]; then
		echo "DEPENDANCIES: Final File ${PRIORTASK1} does NOT exist - something has gone wrong and no samples have run!"
		exit 888
	else 
		echo "DEPENDANCIES: Final File ${PRIORTASK1} from previous job exists: ${FILE1}"
	fi 
LENGTHJOBS1=$(cat ${FILE1} | wc -l)
###########################
FILE2="COMMANDLOGS/job_${SAMPLES_FILE_POST}_${PRIORTASK2}.txt"
echo ${FILE2}
	if [ ! -e "$FILE2" ]; then
		echo "DEPENDANCIES: Final File ${PRIORTASK2} NOT exist - something has gone wrong and no samples have run!"
		exit 888
	else 
		echo "DEPENDANCIES: Final File ${PRIORTASK2} from previous job exists: ${FILE2}"
	fi 
LENGTHJOBS2=$(cat ${FILE2} | wc -l)
###########################
FILE3="COMMANDLOGS/job_${SAMPLES_FILE_POST}_${PRIORTASK3}.txt"
echo ${FILE3}
	if [ ! -e "$FILE3" ]; then
		echo "DEPENDANCIES: Final File ${PRIORTASK3} does NOT exist - something has gone wrong and no samples have run!"
		exit 888
	else 
		echo "DEPENDANCIES: Final File ${PRIORTASK3} from previous job exists: ${FILE3}"
	fi 
LENGTHJOBS3=$(cat ${FILE3} | wc -l)
##
echo "Looking for failed samples.." 
c1grep -v -f ${FILE1} ${IDS_FILE} > COMMANDLOGS/job_${SAMPLES_FILE_POST}_${PRIORTASK1}_FAILED_SAMPLES.txt
c1grep -v -f ${FILE2} ${IDS_FILE} > COMMANDLOGS/job_${SAMPLES_FILE_POST}_${PRIORTASK2}_FAILED_SAMPLES.txt
c1grep -v -f ${FILE3} ${IDS_FILE} > COMMANDLOGS/job_${SAMPLES_FILE_POST}_${PRIORTASK3}_FAILED_SAMPLES.txt
echo "DONE"
##
## Note we will only error the isotyper script if ISO1 (productive+ nonproductive is not successful
## If ISO1 is successful, then any failures in productive/non productive are most likely due to a sample having 0 reads in either category hence we can still continue with analysis! 
if [[ $LENGTHJOBS1 -ne $SAMPLES ]]; then 
	echo "DEPENDANCIES: Final File ${PRIORTASK1}: Not all dependancies ran sucessfully to completion"
	echo "ERROR: NO FURTHER ANALYSIS WILL BE PERFORMED"
	echo "NOTE: This is **most likely due to a sample containing 0 total reads** for isotyoer analysis - investigate further!!"
	exit 999
else 
	echo "DEPENDANCIES: All dependancies ran sucessfully to completion"
	echo "SUCCESS: ANALYSIS STAGE ${TASK} WILL BE PERFORMED"
fi
##
if [[ $LENGTHJOBS2 -ne $SAMPLES ]]; then 
	echo "DEPENDANCIES: Final File ${PRIORTASK2}: Not all dependancies ran sucessfully to completion"
	echo "WARNING: ANALYSIS STAGE ${TASK} WILL BE PERFORMED ON INCOMPLETE SAMPLE SET"
	echo "Note: This is **most likely** due to a sample containing 0 PRODUCTIVE READs for isotyper analysis - investigate further!!"
else 
	echo "DEPENDANCIES: All dependancies ran sucessfully to completion"
	echo "SUCCESS: ANALYSIS STAGE ${TASK} WILL BE PERFORMED"
fi
##
if [[ $LENGTHJOBS3 -ne $SAMPLES ]]; then 
	echo "DEPENDANCIES: Final File ${PRIORTASK2}: Not all dependancies ran sucessfully to completion"
	echo "WARNING: ANALYSIS STAGE ${TASK} WILL BE PERFORMED ON INCOMPLETE SAMPLE SET"
	echo "Note: This is **most likely** due to a sample containing 0 NON-PRODUCTIVE READS for isotyper analysis - investigate further!!"
else 
	echo "DEPENDANCIES: All dependancies ran sucessfully to completion"
	echo "SUCCESS: ANALYSIS STAGE ${TASK} WILL BE PERFORMED"
fi	
fi 

## Creating Relevant Subdirectories 
mkdir ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE || true
mkdir ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE || true
mkdir ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL || true

for file in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*_productive*
do
  mv "$file" ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE || true
done

for file in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*_unproductive*
do
  mv "$file" ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE || true
done

for file in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*.txt
do
  mv "$file" ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL || true
done

# Moving Sampling depth file up one level 
mv ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Sampling_depth_per_isotype_* ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ || true

#Sorting the output files into named directories to make them easier to cat 
##Move them one at a time 
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}CDR3_VJ_genes || true
	mv ${d}CDR3_VJ_genes_*  ${d}CDR3_VJ_genes || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}CDR3_charge || true
	mv ${d}CDR3_charge_*  ${d}CDR3_charge || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}CDR3_grouped_length || true
	mv ${d}CDR3_grouped_length_*  ${d}CDR3_grouped_length || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}CDR3_length_array || true
	mv ${d}CDR3_length_array_*  ${d}CDR3_length_array || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}CDR3_length_distribution || true
	mv ${d}CDR3_length_distribution_*  ${d}CDR3_length_distribution || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}CDR3_lengths_overall || true
	mv ${d}CDR3_lengths_overall_*  ${d}CDR3_lengths_overall || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}CDR3_n_unique_regions_per_isotype_group || true
	mv ${d}CDR3_n_unique_regions_per_isotype_group_*  ${d}CDR3_n_unique_regions_per_isotype_group || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}CDR_charge || true
	mv ${d}CDR_charge_*  ${d}CDR_charge || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Cluster_cassification_per_cluster_size || true
	mv ${d}Cluster_cassification_per_cluster_size_*  ${d}Cluster_cassification_per_cluster_size || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Cluster_expansion_isotype || true
	mv ${d}Cluster_expansion_isotype_*  ${d}Cluster_expansion_isotype || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED || true
	mv ${d}Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_*  ${d}Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Cluster_per_sequence_network_parameters || true
	mv ${d}Cluster_per_sequence_network_parameters_*  ${d}Cluster_per_sequence_network_parameters || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Isotype_counts_shared || true
	mv ${d}Isotype_counts_shared_*  ${d}Isotype_counts_shared || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Isotype_overlapping_frequencies || true
	mv ${d}Isotype_overlapping_frequencies_*  ${d}Isotype_overlapping_frequencies || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Isotype_usages_SUBSAMPLED || true
	mv ${d}Isotype_usages_SUBSAMPLED_*  ${d}Isotype_usages_SUBSAMPLED || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}J_gene_grouped_isotype_frequency || true
	mv ${d}J_gene_grouped_isotype_frequency_*  ${d}J_gene_grouped_isotype_frequency || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Mutations_VJ_genes || true
	mv ${d}Mutations_VJ_genes_*  ${d}Mutations_VJ_genes || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}SHM_Indel_summary || true
	mv ${d}SHM_Indel_summary_*  ${d}SHM_Indel_summary || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}SHM_Mutation_summmary_selection || true
	mv ${d}SHM_Mutation_summmary_selection_*  ${d}SHM_Mutation_summmary_selection || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}SHM_Mutations_per_expanded_cluster || true
	mv ${d}SHM_Mutations_per_expanded_cluster_*  ${d}SHM_Mutations_per_expanded_cluster || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}SHM_Unmutated_sequences || true
	mv ${d}SHM_Unmutated_sequences_*  ${d}SHM_Unmutated_sequences || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}SHM_per_cluster_Mutational_development_classifiation || true
	mv ${d}SHM_per_cluster_Mutational_development_classifiation_*  ${d}SHM_per_cluster_Mutational_development_classifiation || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Secondary_rearrangements_file_SAMPLED || true
	mv ${d}Secondary_rearrangements_file_SAMPLED_*  ${d}Secondary_rearrangements_file_SAMPLED || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Secondary_rearrangements_RAW || true
	mv ${d}Secondary_rearrangements_RAW_*  ${d}Secondary_rearrangements_RAW || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Secondary_rearrangements_clone_sizes || true
	mv ${d}Secondary_rearrangements_clone_sizes_*  ${d}Secondary_rearrangements_clone_sizes || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Secondary_rearrangements_V_genes || true
	mv ${d}Secondary_rearrangements_V_genes_*  ${d}Secondary_rearrangements_V_genes || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Secondary_rearrangements || true
	mv ${d}Secondary_rearrangements_*.txt  ${d}Secondary_rearrangements || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}V_gene_IGHV4_34_quantification || true
	mv ${d}V_gene_IGHV4_34_quantification_*.txt  ${d}V_gene_IGHV4_34_quantification || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}V_gene_grouped_isotype_frequency || true
	mv ${d}V_gene_grouped_isotype_frequency_*.txt  ${d}V_gene_grouped_isotype_frequency || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}V_gene_isotype_frequency || true
	mv ${d}V_gene_isotype_frequency_*.txt  ${d}V_gene_isotype_frequency || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}V_gene_per_cluster_VJ_gene_usage_by_cluster_classification || true
	mv ${d}V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_*.txt  ${d}V_gene_per_cluster_VJ_gene_usage_by_cluster_classification || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}V_gene_summary_cluster_file || true
	mv ${d}V_gene_summary_cluster_file_*.txt  ${d}V_gene_summary_cluster_file || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Isotype_normalised_overlap_frequencies_uniq || true
	mv ${d}Isotye_normalised_overlap_frequencies_uniq_*.txt  ${d}Isotype_normalised_overlap_frequencies_uniq || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Classification_per_sequence || true
	mv ${d}Classification_per_sequence_*.txt  ${d}Classification_per_sequence || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}Classification_per_cluster || true
	mv ${d}Classification_per_cluster_*.txt  ${d}Classification_per_cluster || true
done
##
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	mkdir ${d}IMGT_trimmed_sequences || true
	mv ${d}IMGT_trimmed_sequences_*.txt  ${d}IMGT_trimmed_sequences || true
done

## Cat the files to form 'ALL' summary files 
for d in ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/*/ ; do
	part2=$(basename "${d}")
	echo $part2
	echo $d
	cat ${d}CDR3_VJ_genes/CDR3_VJ_genes_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_CDR3_VJ_genes_${part2}.txt || true
	cat ${d}CDR3_charge/CDR3_charge_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_CDR3_charge_${part2}.txt || true
	cat ${d}CDR3_grouped_length/CDR3_grouped_length_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_CDR3_grouped_length_${part2}.txt || true
	cat ${d}CDR3_length_array/CDR3_length_array_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_CDR3_length_array_${part2}.txt || true
	cat ${d}CDR3_length_distribution/CDR3_length_distribution_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_CDR3_length_distribution_${part2}.txt || true
	cat ${d}CDR3_lengths_overall/CDR3_lengths_overall_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_CDR3_lengths_overall_${part2}.txt || true
	cat ${d}CDR3_n_unique_regions_per_isotype_group/CDR3_n_unique_regions_per_isotype_group_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_CDR3_n_unique_regions_per_isotype_group_${part2}.txt || true
	cat ${d}CDR_charge/CDR_charge_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_CDR_charge_${part2}.txt || true 
	cat ${d}Cluster_cassification_per_cluster_size/Cluster_cassification_per_cluster_size_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_cassification_per_cluster_size_${part2}.txt || true
	cat ${d}Cluster_expansion_isotype/Cluster_expansion_isotype_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_expansion_isotype_${part2}.txt || true
	cat ${d}Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED/Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_${part2}.txt || true
	cat ${d}Cluster_per_sequence_network_parameters/Cluster_per_sequence_network_parameters_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_per_sequence_network_parameters_${part2}.txt || true
	cat ${d}Isotype_counts_shared/Isotype_counts_shared_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Isotype_counts_shared_${part2}.txt || true
	cat ${d}Isotype_usages_SUBSAMPLED/Isotype_usages_SUBSAMPLED_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Isotype_usages_SUBSAMPLED_${part2}.txt || true
	cat ${d}J_gene_grouped_isotype_frequency/J_gene_grouped_isotype_frequency_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_J_gene_grouped_isotype_frequency_${part2}.txt || true
	cat ${d}Mutations_VJ_genes/Mutations_VJ_genes_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Mutations_VJ_genes_${part2}.txt || true
	cat ${d}SHM_Indel_summary/SHM_Indel_summary_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_SHM_Indel_summary_${part2}.txt || true
	cat ${d}SHM_Mutation_summmary_selection/SHM_Mutation_summmary_selection_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_SHM_Mutation_summmary_selection_${part2}.txt || true
	cat ${d}SHM_Mutations_per_expanded_cluster/SHM_Mutations_per_expanded_cluster_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_SHM_Mutations_per_expanded_cluster_${part2}.txt || true
	cat ${d}SHM_Unmutated_sequences/SHM_Unmutated_sequences_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_SHM_Unmutated_sequences_${part2}.txt || true
	cat ${d}SHM_per_cluster_Mutational_development_classifiation/SHM_per_cluster_Mutational_development_classifiation_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_SHM_per_cluster_Mutational_development_classifiation_${part2}.txt || true
	cat ${d}Secondary_rearrangements/Secondary_rearrangements_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Secondary_rearrangements_${part2}.txt || true
	cat ${d}Secondary_rearrangements_RAW/Secondary_rearrangements_RAW_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Secondary_rearrangements_RAW_${part2}.txt || true
	cat ${d}Secondary_rearrangements_file_SAMPLED/Secondary_rearrangements_file_SAMPLED_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Secondary_rearrangements_file_SAMPLED_${part2}.txt || true
	cat ${d}V_gene_IGHV4_34_quantification/V_gene_IGHV4_34_quantification_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_IGHV4_34_quantification_${part2}.txt || true
	cat ${d}V_gene_grouped_isotype_frequency/V_gene_grouped_isotype_frequency_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_grouped_isotype_frequency_${part2}.txt || true
	cat ${d}V_gene_isotype_frequency/V_gene_isotype_frequency_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_isotype_frequency_${part2}.txt || true
	cat ${d}V_gene_per_cluster_VJ_gene_usage_by_cluster_classification/V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_${part2}.txt || true
	cat ${d}V_gene_summary_cluster_file/V_gene_summary_cluster_file_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_summary_cluster_file_${part2}.txt || true
	cat ${d}Secondary_rearrangements_clone_sizes/Secondary_rearrangements_clone_sizes_* > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Secondary_rearrangements_clone_sizes_${part2}.txt || true 
	cat ${d}Isotype_normalised_overlap_frequencies_uniq/Isotye_normalised_overlap_frequencies_uniq_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Isotype_normalised_overlap_frequencies_uniq_${part2}.txt || true
	cat ${d}Isotype_overlapping_frequencies/Isotype_overlapping_frequencies_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Isotype_overlapping_frequencies_${part2}.txt || true
	cat ${d}Classification_per_sequence/Classification_per_sequence_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Classification_per_sequence_${part2}.txt || true
	cat ${d}Classification_per_cluster/Classification_per_cluster_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/All_Classification_per_clusterR_${part2}.txt || true	
	cat ${d}IMGT_trimmed_sequences/IMGT_trimmed_sequences_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/IMGT_trimmed_sequences_${part2}.txt || true	
	cat ${d}Secondary_rearrangements_V_genes/Secondary_rearrangements_V_genes_*.txt > ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/Secondary_rearrangements_V_genes_${part2}.txt || true	
done

## IF JOB RUN SUCESSFULLY SAVE TO SAMPLE COUNTER FILE 
NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_CAT.txt"
eval "${NWCMD}"

echo "DONE PART 1"
module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
echo "Loaded R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0 Module"

echo "RUNNING ISOTYPE ANALYSIS PLOTS"
CMD="Rscript ISOTYPER_ANALYSIS.R -o ${OUTPUTDIR} -s ${SAMPLES_FILE_POST} -g ${GENE} -l ${LAYOUTS_FILE}"
echo ${CMD}
eval "${CMD}" 
echo "DONE"

## IF JOB RUN SUCESSFULLY SAVE TO SAMPLE COUNTER FILE 
NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_ISOPLOTTING.txt"
eval "${NWCMD}"

# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0