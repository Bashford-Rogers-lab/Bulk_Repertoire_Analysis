module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load R/4.1.0-foss-2021a 
########################
library(ggplot2) 
library(data.table) 
library(stringr)
library(ggpubr)
library(rmcorr) 
library(lme4)
library(lmerTest)
library(ggplot2) 
library(M3C) 
library(lmerTest)
###################################################

####-----------------------------------------------------------------------------------------------------
cohort1_file <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/logcpm_181_20600_cohort2.txt'
cohort2_file <- "/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/Logcpm_864_20416.txt"
cohort_assignment <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/SampleOrganisation/RepertoireRNAseqCohorts.txt'

## BCR
eigenvectors_file <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt'
outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/GeneExpPlots/'
gene_expression_file <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/REPERTOIRE_SAMPLES_GEX_EnsembeID.txt'

# TCRAB
eigenvectors_file <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Eigenvectors_No_Technical_TCRAB_PRODUCTIVE.txt'
outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/GeneExpPlots/'
gene_expression_file <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/REPERTOIRE_SAMPLES_GEX_EnsembeID.txt'

#TCRGD
eigenvectors_file <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Eigenvectors_No_Technical_TCRGD_PRODUCTIVE.txt'
gene_expression_file <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/REPERTOIRE_SAMPLES_GEX_EnsembeID.txt'
outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/GeneExpPlots/'

#TCRGDNEW
eigenvectors_file <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW/Eigenvectors_No_Technical_TCRGD_PRODUCTIVE.txt'
gene_expression_file <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/REPERTOIRE_SAMPLES_GEX_EnsembeID.txt'
outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW/GeneExpPlots/'

####-----------------------------------------------------------------------------------------------------
## First look for batch effect in Gene Expression DATA!
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/GeneExpression/batcheffect_assessment.R')
assess_batcheffect(cohort1_file, cohort2_file, outputdir, cohort_assignment)
###
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/GeneExpression/CorrelateGenes.R')
correlate_genes(cohort1_file, cohort2_file, eigenvectors_file, gene_df, outputdir, cohort_assignment, module_list, list_name)


####-----------------------------------------------------------------------------------------------------
## Quantitative Trait Analysis
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/GeneExpression/Function_QuantitativeTrait.R')
### Lets do this for all of the modules !!!!
for(i in c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25",  "Module_26",  "Module_27",  "Module_28",  "Module_29",  "Module_30",  "Module_31")){
	print(i)
	get_best_genes(gene_expression_file, eigenvectors_file, outputdir, i)
	}
	
#### Now lets do the GO analysis this needs WIFI
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/GeneExpression/Function_Gene_Ontology.R')
files_LME <- list.files(paste0(outputdir, 'ALL_GENES/GeneTables/'), full.name=TRUE)
files_LME_a <- grep("significant", files_LME, value=TRUE)

#############################
for(i in files_LME_a){
	gene_ontology_analysis(i, outputdir)
}



####-----------------------------------------------------------------------------------------------------
### Gene Lists 
## Housekeepers
gene_df <- read.delim('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/GeneExpression/housekeepers.txt', sep="\t", header=FALSE)
module_list <- "Module_3"
list_name <- "housekeepers"

## Mod 1 
gene_df <- read.delim('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/GeneExpression/Module_1_genes.txt', sep="\t", header=FALSE)
module_list <- "Module_1"
list_name <- "Module1_genelist"

## Mod 3 
gene_df <- read.delim('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/GeneExpression/Module_3_genes.txt', sep="\t", header=FALSE)
module_list <- "Module_3"
list_name <- "Module3_genelist"

