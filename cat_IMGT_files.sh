#!/bin/bash
#$ -cwd

# If there's an error, fail the whole script
set -e -o pipefail

# Set permissions so any user in the group can 
# read/write what it's created by the script
umask 002

# Extracting Outputdir
OUTPUTDIR=$1

## Creating Relevant Subdirectories 
mkdir ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE
mv *_productive* ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/PRODUCTIVE

mkdir ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE
mv *_unproductive* ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/UNPRODUCTIVE

mkdir ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL
mv *.txt ${OUTPUTDIR}ORIENTATED_SEQUENCES/ISOTYPER/ALL

# Moving Sampling deapth file
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

# Loading R module to do summary Plots 
module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
echo "Loaded R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0 Module"


## Call R script

#ISOTYPER_ANALYSIS
### TRIAL



cat PRODUCTIVE/CDR3_VJ_genes_*.txt > All_CDR3_VJ_genes_PRODUCTIVE.txt || true
cat PRODUCTIVE/CDR3_charge_*.txt > All_CDR3_charge_PRODUCTIVE.txt || true
cat PRODUCTIVE/CDR3_grouped_length_*.txt > All_CDR3_grouped_length_PRODUCTIVE.txt || true
cat PRODUCTIVE/CDR3_length_array_*.txt > All_CDR3_length_array_PRODUCTIVE.txt || true
cat PRODUCTIVE/CDR3_length_distribution_*.txt > All_CDR3_length_distribution_PRODUCTIVE.txt || true
cat PRODUCTIVE/CDR3_lengths_overall_*.txt > All_CDR3_lengths_overall_PRODUCTIVE.txt || true
cat PRODUCTIVE/CDR3_n_unique_regions_per_isotype_group_*.txt > All_CDR3_n_unique_regions_per_isotype_group_PRODUCTIVE.txt || true
cat PRODUCTIVE/CDR_charge_*.txt > All_CDR_charge_PRODUCTIVE.txt || true 
cat PRODUCTIVE/Cluster_cassification_per_cluster_size_*.txt > All_Cluster_cassification_per_cluster_size_PRODUCTIVE.txt || true
cat PRODUCTIVE/Cluster_expansion_isotype_*.txt > All_Cluster_expansion_isotype_PRODUCTIVE.txt || true
cat PRODUCTIVE/Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_*.txt > All_Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_PRODUCTIVE.txt || true
cat PRODUCTIVE/Cluster_per_sequence_network_parameters_*.txt > All_Cluster_per_sequence_network_parameters_PRODUCTIVE.txt || true
cat PRODUCTIVE/Isotype_counts_shared_*.txt > All_Isotype_counts_shared_PRODUCTIVE.txt || true
cat PRODUCTIVE/Isotype_overlapping_frequencies_*.txt > All_Isotype_overlapping_frequencies_PRODUCTIVE.txt || true
cat PRODUCTIVE/Isotype_unexpanded_cluster_*.txt > All_Isotype_unexpanded_cluster_PRODUCTIVE.txt || true
cat PRODUCTIVE/Isotype_usages_SUBSAMPLED_*.txt > All_Isotype_usages_SUBSAMPLED_PRODUCTIVE.txt || true
cat PRODUCTIVE/J_gene_grouped_isotype_frequency_*.txt > All_J_gene_grouped_isotype_frequency_PRODUCTIVE.txt || true
cat PRODUCTIVE/Mutations_VJ_genes_*.txt > All_Mutations_VJ_genes_PRODUCTIVE.txt || true
cat PRODUCTIVE/SHM_Cluster_mutation_sharing_probability_*.txt > All_SHM_Cluster_mutation_sharing_probability_PRODUCTIVE.txt || true
cat PRODUCTIVE/SHM_Indel_summary_*.txt > All_SHM_Indel_summary_PRODUCTIVE.txt || true
cat PRODUCTIVE/SHM_Mutation_summmary_selection_*.txt > All_SHM_Mutation_summmary_selection_PRODUCTIVE.txt || true
cat PRODUCTIVE/SHM_Mutational_positions_summmary_per_gene_*.txt > All_SHM_Mutational_positions_summmary_per_gene_PRODUCTIVE.txt || true
cat PRODUCTIVE/SHM_Mutations_per_expanded_cluster_*.txt > All_SHM_Mutations_per_expanded_cluster_PRODUCTIVE.txt || true
cat PRODUCTIVE/SHM_Unmutated_sequences_*.txt > All_SHM_Unmutated_sequences_PRODUCTIVE.txt || true
cat PRODUCTIVE/SHM_per_cluster_Mutational_development_classifiation_*.txt > All_SHM_per_cluster_Mutational_development_classifiation_PRODUCTIVE.txt || true
cat PRODUCTIVE/Secondary_rearrangements_*.txt > All_Secondary_rearrangements_PRODUCTIVE.txt || true
cat PRODUCTIVE/Secondary_rearrangements_RAW_*.txt > All_Secondary_rearrangements_RAW_PRODUCTIVE.txt || true
cat PRODUCTIVE/Secondary_rearrangements_V_genes_*.txt > All_Secondary_rearrangements_V_genes_PRODUCTIVE.txt || true
cat PRODUCTIVE/Secondary_rearrangements_file_SAMPLED_*.txt > All_Secondary_rearrangements_file_SAMPLED_PRODUCTIVE.txt || true
cat PRODUCTIVE/V_gene_IGHV4_34_quantification_*.txt > All_V_gene_IGHV4_34_quantification_PRODUCTIVE.txt || true
cat PRODUCTIVE/V_gene_grouped_isotype_frequency_*.txt > All_V_gene_grouped_isotype_frequency_PRODUCTIVE.txt || true
cat PRODUCTIVE/V_gene_isotype_frequency_*.txt > All_V_gene_isotype_frequency_PRODUCTIVE.txt || true
cat PRODUCTIVE/V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_*.txt > All_V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_PRODUCTIVE.txt || true
cat PRODUCTIVE/V_gene_summary_cluster_file_*.txt > All_V_gene_summary_cluster_file_PRODUCTIVE.txt || true
cat PRODUCTIVE/Isotye_normalised_overlap_frequencies_uniq_BCR_*.txt >All_Isotype_normalised_overlap_frequencies_uniq_BCR_PRODUCTIVE.txt || true
cat PRODUCTIVE/Isotye_normalised_overlap_frequencies_uniq_TCR_*.txt >All_Isotype_normalised_overlap_frequencies_uniq_TCR_PRODUCTIVE.txt || true
cat PRODUCTIVE/Secondary_rearrangements_clone_sizes_* > All_Secondary_rearrangements_clone_sizes_PRODUCTIVE.txt || true 
