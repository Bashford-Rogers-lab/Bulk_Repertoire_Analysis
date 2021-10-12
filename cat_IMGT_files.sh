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

## Job arguments 
OUTPUTDIR=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$8}" $SAMPLES_FILE_POST)

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
SAMPLES=$((SAMPLES1+1))

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
FILE2="COMMANDLOGS/job_${SAMPLES_FILE_POST}_${PRIORTASK3}.txt"
echo ${FILE2}
	if [ ! -e "$FILE2" ]; then
		echo "DEPENDANCIES: Final File ${PRIORTASK2} NOT exist - something has gone wrong and no samples have run!"
		exit 888
	else 
		echo "DEPENDANCIES: Final File ${PRIORTASK2} from previous job exists: ${FILE2}"
	fi 
LENGTHJOBS2=$(cat ${FILE2} | wc -l)
###########################
FILE3="COMMANDLOGS/job_${SAMPLES_FILE_POST}_${PRIORTASK2}.txt"
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
if [[ $LENGTHJOBS1 -ne $SAMPLES ]]; then 
	echo "DEPENDANCIES: Final File ${PRIORTASK1}: Not all dependancies ran sucessfully to completion"
	echo "ERROR: NO FURTHER ANALYSIS WILL BE PERFORMED"
	exit 999
else 
	echo "DEPENDANCIES: All dependancies ran sucessfully to completion"
	echo "SUCCESS: ANALYSIS STAGE ${TASK} WILL BE PERFORMED"
fi
##
if [[ $LENGTHJOBS2 -ne $SAMPLES ]]; then 
	echo "DEPENDANCIES: Final File ${PRIORTASK2}: Not all dependancies ran sucessfully to completion"
	echo "ERROR: NO FURTHER ANALYSIS WILL BE PERFORMED"
	exit 999
else 
	echo "DEPENDANCIES: All dependancies ran sucessfully to completion"
	echo "SUCCESS: ANALYSIS STAGE ${TASK} WILL BE PERFORMED"
fi
##
if [[ $LENGTHJOBS3 -ne $SAMPLES ]]; then 
	echo "DEPENDANCIES: Final File ${PRIORTASK2}: Not all dependancies ran sucessfully to completion"
	echo "ERROR: NO FURTHER ANALYSIS WILL BE PERFORMED"
	exit 999
else 
	echo "DEPENDANCIES: All dependancies ran sucessfully to completion"
	echo "SUCCESS: ANALYSIS STAGE ${TASK} WILL BE PERFORMED"
fi	
fi 

## Creating Relevant Subdirectories 
mkdir ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE
mv *_productive* ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE

mkdir ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE
mv *_unproductive* ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE

mkdir ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL
mv *.txt ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL

# Moving Sampling depth file up one level 
mv ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Sampling_depth_per_isotype_* ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/

#Making into Combined Files
#ALL
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/CDR3_VJ_genes_*.txt > All_CDR3_VJ_genes_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/CDR3_charge_*.txt > All_CDR3_charge_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/CDR3_grouped_length_*.txt > All_CDR3_grouped_length_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/CDR3_length_array_*.txt > All_CDR3_length_array_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/CDR3_length_distribution_*.txt > All_CDR3_length_distribution_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/CDR3_lengths_overall_*.txt > All_CDR3_lengths_overall_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/CDR3_n_unique_regions_per_isotype_group_*.txt > All_CDR3_n_unique_regions_per_isotype_group_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/CDR_charge_*.txt > All_CDR_charge_ALL.txt || true 
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Cluster_cassification_per_cluster_size_*.txt > All_Cluster_cassification_per_cluster_size_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Cluster_expansion_isotype_*.txt > All_Cluster_expansion_isotype_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_*.txt > All_Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Cluster_per_sequence_network_parameters_*.txt > All_Cluster_per_sequence_network_parameters_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Isotype_counts_shared_*.txt > All_Isotype_counts_shared_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Isotype_overlapping_frequencies_*.txt > All_Isotype_overlapping_frequencies_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Isotype_usages_SUBSAMPLED_*.txt > All_Isotype_usages_SUBSAMPLED_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/J_gene_grouped_isotype_frequency_*.txt > All_J_gene_grouped_isotype_frequency_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Mutations_VJ_genes_*.txt > All_Mutations_VJ_genes_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/SHM_Indel_summary_*.txt > All_SHM_Indel_summary_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/SHM_Mutation_summmary_selection_*.txt > All_SHM_Mutation_summmary_selection_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/SHM_Mutations_per_expanded_cluster_*.txt > All_SHM_Mutations_per_expanded_cluster_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/SHM_Unmutated_sequences_*.txt > All_SHM_Unmutated_sequences_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/SHM_per_cluster_Mutational_development_classifiation_*.txt > All_SHM_per_cluster_Mutational_development_classifiation_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Secondary_rearrangements_*.txt > All_Secondary_rearrangements_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Secondary_rearrangements_RAW_*.txt > All_Secondary_rearrangements_RAW_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Secondary_rearrangements_file_SAMPLED_*.txt > All_Secondary_rearrangements_file_SAMPLED_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/V_gene_IGHV4_34_quantification_*.txt > All_V_gene_IGHV4_34_quantification_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/V_gene_grouped_isotype_frequency_*.txt > All_V_gene_grouped_isotype_frequency_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/V_gene_isotype_frequency_*.txt > All_V_gene_isotype_frequency_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_*.txt > All_V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/V_gene_summary_cluster_file_*.txt > All_V_gene_summary_cluster_file_ALL.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Secondary_rearrangements_clone_sizes_* > All_Secondary_rearrangements_clone_sizes_ALL.txt || true 
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Isotye_normalised_overlap_frequencies_uniq_BCR*_ALL.txt >All_Isotype_normalised_overlap_frequencies_uniq_BCR.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL/Isotye_normalised_overlap_frequencies_uniq_TCR*_ALL.txt >All_Isotype_normalised_overlap_frequencies_uniq_TCR.txt || true

#PRODUCTIVE 
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/CDR3_VJ_genes_*_productive.txt > All_CDR3_VJ_genes_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/CDR3_charge_*_productive.txt > All_CDR3_charge_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/CDR3_grouped_length_*_productive.txt > All_CDR3_grouped_length_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/CDR3_length_array_*_productive.txt > All_CDR3_length_array_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/CDR3_length_distribution_*_productive.txt > All_CDR3_length_distribution_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/CDR3_lengths_overall_*_productive.txt > All_CDR3_lengths_overall_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/CDR3_n_unique_regions_per_isotype_group_*_productive.txt > All_CDR3_n_unique_regions_per_isotype_group_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/CDR_charge_*_productive.txt > All_CDR_charge_PRODUCTIVE.txt || true 
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Cluster_cassification_per_cluster_size_*_productive.txt > All_Cluster_cassification_per_cluster_size_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Cluster_expansion_isotype_*_productive.txt > All_Cluster_expansion_isotype_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_*_productive.txt > All_Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Cluster_per_sequence_network_parameters_*_productive.txt > All_Cluster_per_sequence_network_parameters_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Isotype_counts_shared_*_productive.txt > All_Isotype_counts_shared_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Isotype_overlapping_frequencies_*_productive.txt > All_Isotype_overlapping_frequencies_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Isotype_usages_SUBSAMPLED_*_productive.txt > All_Isotype_usages_SUBSAMPLED_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/J_gene_grouped_isotype_frequency_*_productive.txt > All_J_gene_grouped_isotype_frequency_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Mutations_VJ_genes_*_productive.txt > All_Mutations_VJ_genes_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/SHM_Indel_summary_*_productive.txt > All_SHM_Indel_summary_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/SHM_Mutation_summmary_selection_*_productive.txt > All_SHM_Mutation_summmary_selection_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/SHM_Mutations_per_expanded_cluster_*_productive.txt > All_SHM_Mutations_per_expanded_cluster_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/SHM_Unmutated_sequences_*_productive.txt > All_SHM_Unmutated_sequences_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/SHM_per_cluster_Mutational_development_classifiation_*_productive.txt > All_SHM_per_cluster_Mutational_development_classifiation_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Secondary_rearrangements_*_productive.txt > All_Secondary_rearrangements_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Secondary_rearrangements_RAW_*_productive.txt > All_Secondary_rearrangements_RAW_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Secondary_rearrangements_file_SAMPLED_*_productive.txt > All_Secondary_rearrangements_file_SAMPLED_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/V_gene_IGHV4_34_quantification_*_productive.txt > All_V_gene_IGHV4_34_quantification_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/V_gene_grouped_isotype_frequency_*_productive.txt > All_V_gene_grouped_isotype_frequency_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/V_gene_isotype_frequency_*_productive.txt > All_V_gene_isotype_frequency_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_*_productive.txt > All_V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/V_gene_summary_cluster_file_*_productive.txt > All_V_gene_summary_cluster_file_PRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Secondary_rearrangements_clone_sizes_*_productive.txt > All_Secondary_rearrangements_clone_sizes_PRODUCTIVE.txt || true 
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Isotye_normalised_overlap_frequencies_uniq_BCR*_productive.txt >All_Isotype_normalised_overlap_frequencies_uniq_BCR.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE/Isotye_normalised_overlap_frequencies_uniq_TCR*_productive.txt >All_Isotype_normalised_overlap_frequencies_uniq_TCR.txt || true

#UNPRODUCTIVE
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/CDR3_VJ_genes_*_unproductive.txt > All_CDR3_VJ_genes_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/CDR3_charge_*_unproductive.txt > All_CDR3_charge_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/CDR3_grouped_length_*_unproductive.txt > All_CDR3_grouped_length_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/CDR3_length_array_*_unproductive.txt > All_CDR3_length_array_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/CDR3_length_distribution_*_unproductive.txt > All_CDR3_length_distribution_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/CDR3_lengths_overall_*_unproductive.txt > All_CDR3_lengths_overall_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/CDR3_n_unique_regions_per_isotype_group_*_unproductive.txt > All_CDR3_n_unique_regions_per_isotype_group_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/CDR_charge_*_unproductive.txt > All_CDR_charge_UNPRODUCTIVE.txt || true 
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Cluster_cassification_per_cluster_size_*_unproductive.txt > All_Cluster_cassification_per_cluster_size_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Cluster_expansion_isotype_*_unproductive.txt > All_Cluster_expansion_isotype_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_*_unproductive.txt > All_Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Cluster_per_sequence_network_parameters_*_unproductive.txt > All_Cluster_per_sequence_network_parameters_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Isotype_counts_shared_*_unproductive.txt > All_Isotype_counts_shared_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Isotype_overlapping_frequencies_*_unproductive.txt > All_Isotype_overlapping_frequencies_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Isotype_usages_SUBSAMPLED_*_unproductive.txt > All_Isotype_usages_SUBSAMPLED_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/J_gene_grouped_isotype_frequency_*_unproductive.txt > All_J_gene_grouped_isotype_frequency_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Mutations_VJ_genes_*_unproductive.txt > All_Mutations_VJ_genes_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/SHM_Indel_summary_*_unproductive.txt > All_SHM_Indel_summary_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/SHM_Mutation_summmary_selection_*_unproductive.txt > All_SHM_Mutation_summmary_selection_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/SHM_Mutations_per_expanded_cluster_*_unproductive.txt > All_SHM_Mutations_per_expanded_cluster_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/SHM_Unmutated_sequences_*_unproductive.txt > All_SHM_Unmutated_sequences_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/SHM_per_cluster_Mutational_development_classifiation_*_unproductive.txt > All_SHM_per_cluster_Mutational_development_classifiation_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Secondary_rearrangements_*_unproductive.txt > All_Secondary_rearrangements_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Secondary_rearrangements_RAW_*_unproductive.txt > All_Secondary_rearrangements_RAW_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Secondary_rearrangements_file_SAMPLED_*_unproductive.txt > All_Secondary_rearrangements_file_SAMPLED_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/V_gene_IGHV4_34_quantification_*_unproductive.txt > All_V_gene_IGHV4_34_quantification_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/V_gene_grouped_isotype_frequency_*_unproductive.txt > All_V_gene_grouped_isotype_frequency_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/V_gene_isotype_frequency_*_unproductive.txt > All_V_gene_isotype_frequency_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_*_unproductive.txt > All_V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/V_gene_summary_cluster_file_*_unproductive.txt > All_V_gene_summary_cluster_file_UNPRODUCTIVE.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Secondary_rearrangements_clone_sizes_*_unproductive.txt > All_Secondary_rearrangements_clone_sizes_UNPRODUCTIVE.txt || true 
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Isotye_normalised_overlap_frequencies_uniq_BCR*_unproductive.txt >All_Isotype_normalised_overlap_frequencies_uniq_BCR.txt || true
cat ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE/Isotye_normalised_overlap_frequencies_uniq_TCR*_unproductive.txt >All_Isotype_normalised_overlap_frequencies_uniq_TCR.txt || true

## IF JOB RUN SUCESSFULLY SAVE TO SAMPLE COUNTER FILE 
NWCMD="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_CAT.txt"
eval "${NWCMD}"

# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0