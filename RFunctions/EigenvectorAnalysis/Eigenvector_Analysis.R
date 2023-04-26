## Functions to correlate eigenvectors to clinical variables!!
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk
## September 2022
## Need to use new version of R !! for rmcorr package
module load R/4.1.0-foss-2021a 
#module load R/4.2.1-foss-2022a
############################################
library(broom)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(multcomp)
library(car)
library(AICcmodavg)
#library(MuMIn)

#############################################
##BCR
eigenvectors <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt'
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR'
metadata <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt'
metahealth <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Healthies_ClinData.txt'
type_receptor <- "BCR"
feature_assignment <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/Clustered_Features_assignment_BCR_PRODUCTIVE_NON_IMPUTED.txt'
imputed_data <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Imputed_DATA_FINAL_SCALED_BCR_PRODUCTIVE.txt'
loadings_use <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/Module_FeaturePCALoadings_BCR_1500_PRODUCTIVE.rds'
srs_assignment <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/SampleOrganisation/Repertoire_SRS_RNA_seq_assignment.txt'
cyber_sort  <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/SampleOrganisation/Repertoire_CYBERSORT_assignment.txt'
gene_expression_file <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/REPERTOIRE_SAMPLES_GEX_EnsembeID.txt'
protein  <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/Yuxin_PROTEOMIC_PROCESSED/REPERTOIRE_SAMPLES_IMMUNOGLOBULIN.txt'
final_data <- "/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Imputed_DATA_FINAL_BCR_PRODUCTIVE.txt"
iso_type <- "PRODUCTIVE"
meta <- '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/sepsis_meta_health.txt'
minClusterSizex <- 20
thresh_sd <- 1.4


#TCRAB
eigenvectors <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Eigenvectors_No_Technical_TCRAB_PRODUCTIVE.txt'
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB'
metadata <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt'
type_receptor <- "TCRAB"
metahealth <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Healthies_ClinData.txt'
feature_assignment <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Summary/Clustered_Features_assignment_TCRAB_PRODUCTIVE_NON_IMPUTED.txt'
imputed_data <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Imputed_DATA_FINAL_SCALED_TCRAB_PRODUCTIVE.txt'
loadings_use <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Summary/Module_FeaturePCALoadings_TCRAB_375_PRODUCTIVE.rds'
srs_assignment <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/SampleOrganisation/Repertoire_SRS_RNA_seq_assignment.txt'
cyber_sort  <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/SampleOrganisation/Repertoire_CYBERSORT_assignment.txt'
final_data <- "/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Imputed_DATA_FINAL_TCRAB_PRODUCTIVE.txt"
iso_type <- "PRODUCTIVE"
meta <- '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/sepsis_meta_health.txt'
minClusterSizex <- 20
thresh_sd <- 1.2

##TCRGD
eigenvectors <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Eigenvectors_No_Technical_TCRGD_PRODUCTIVE.txt'
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD'
metadata <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt'
type_receptor <- "TCRGD"
metahealth <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Healthies_ClinData.txt'
feature_assignment <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Summary/Clustered_Features_assignment_TCRGD_PRODUCTIVE_NON_IMPUTED.txt'
imputed_data <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Imputed_DATA_FINAL_SCALED_TCRGD_PRODUCTIVE.txt'
loadings_use <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Summary/Module_FeaturePCALoadings_TCRGD_200_PRODUCTIVE.rds'
srs_assignment <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/SampleOrganisation/Repertoire_SRS_RNA_seq_assignment.txt'

##TCRGD NEW 
eigenvectors <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW/Eigenvectors_No_Technical_TCRGD_PRODUCTIVE.txt'
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW'
metadata <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt'
type_receptor <- "TCRGD"
metahealth <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Healthies_ClinData.txt'
feature_assignment <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW/Summary/Clustered_Features_assignment_TCRGD_PRODUCTIVE_NON_IMPUTED.txt'
imputed_data <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW/Imputed_DATA_FINAL_SCALED_TCRGD_PRODUCTIVE.txt'
loadings_use <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW/Summary/Module_FeaturePCALoadings_TCRGD_80_PRODUCTIVE.rds'
srs_assignment <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/SampleOrganisation/Repertoire_SRS_RNA_seq_assignment.txt'
final_data <- "/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Imputed_DATA_FINAL_TCRGD_PRODUCTIVE.txt"
iso_type <- "PRODUCTIVE"
meta <- '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/sepsis_meta_health.txt'
minClusterSizex <- 20
thresh_sd <- 0.75

##-----------------------------------------------------------------------------
## This first function looks at all day 1 measures and correlates them to the clinical data provided 
#########################################################################################################
## Function to plot difference in mean values as part of some of the other functions!
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotTukey.R')
## SURVIVAL ANALYSIS
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/SurvivalCurve.R')
plot_survival(eigenvectors, metadata, outputdir)
## STEP 1: IDENTIFY COVARIATE:which would then be added as terms in model....
#source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/Covariate.R')
#p_covariates <- correlate_covariates(eigenvectors, metadata, outputdir,type_receptor, metahealth)
## STEP 2: Model selection
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/Model_Eigenvectors_Selection.R')
model_selection(eigenvectors, metadata, outputdir,type_receptor, metahealth)
## STEP 3: HvsD 
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/Model_Eigenvectors_HvsD.R')
p_hd_extended <- correlate_eigenvectors_hvd(eigenvectors, metadata, outputdir, type_receptor, metahealth, "NO") 
## STEP 4: MORTALITY
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/Model_Eigenvectors_Mortality.R')
p_hd_extended_mort <- correlate_eigenvectors_outcome(eigenvectors, metadata, outputdir, type_receptor, "NO")
## STEP 5: MODULE VIGNETTE PLOTS 
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotModuleVignette.R')
plot_module_vignette(eigenvectors, feature_assignment, imputed_data, outputdir)
## STEP 6: SRS
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/SRSmodel.R')
correlate_eigenvectors_srs(eigenvectors, metadata, outputdir, type_receptor, srs_assignment)
## STEP 7: WITHIN MODULE CORRELATION NETWORK 
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotCorrelationWithin.R')
plot_correlation_within(eigenvectors, outputdir, type_receptor)
#### STEP 8: CYBERSORT ANALYSIS 
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/CybersortAnalysis.R')
correlate_cybersort(eigenvectors, cyber_sort, outputdir, metadata, metahealth)
#### STEP 9: CLINICAL TRAIT ANALYSIS
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/Model_Discrete_Metadata.R')
correlate_eigenvectors_t1(eigenvectors, metadata, outputdir, type_receptor, NA) 
## STEP 10: CORRELATE TO Longitudinal Clinical variables 
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/Model_Continuous_Metadata.R')
correlate_eigenvectors_continuous(eigenvectors, metadata, outputdir, type_receptor)
## Cluster Samples 
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/SampleClustering.R')
cluster_samples(final_data, paste0(outputdir, "/"), type_receptor, iso_type, meta, minClusterSizex, thresh_sd)

##------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------
### Cross-Receptor Analysis!!!!!!!
### Lets now look at the correlation across a pair of features...
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotCorrelationAcross.R')
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/'
eigenvector_list <- list('/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Eigenvectors_No_Technical_TCRAB_PRODUCTIVE.txt', '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt')
type_receptor <- "BCR_TCRAB"
plot_correlation_across(eigenvector_list, outputdir, type_receptor)
####
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/'
eigenvector_list <- list('/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Eigenvectors_No_Technical_TCRAB_PRODUCTIVE.txt', '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW/Eigenvectors_No_Technical_TCRGD_PRODUCTIVE.txt')
type_receptor <- "TCRAB_TCRGD"
plot_correlation_across(eigenvector_list, outputdir, type_receptor)
###
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/'
eigenvector_list <- list('/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt', '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW/Eigenvectors_No_Technical_TCRGD_PRODUCTIVE.txt')
type_receptor <- "BCR_TCRGD"
plot_correlation_across(eigenvector_list, outputdir, type_receptor)## Functions to correlate eigenvectors to clinical variables!!
