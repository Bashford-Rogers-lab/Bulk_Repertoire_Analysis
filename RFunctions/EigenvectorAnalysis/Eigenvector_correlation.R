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

#TCRAB
eigenvectors <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Eigenvectors_No_Technical_TCRAB_PRODUCTIVE.txt'
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB'
metadata <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt'
type_receptor <- "TCRAB"
metahealth <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Healthies_ClinData.txt'
feature_assignment <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Summary/Clustered_Features_assignment_TCRAB_PRODUCTIVE_NON_IMPUTED.txt'
imputed_data <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Imputed_DATA_FINAL_SCALED_TCRAB_PRODUCTIVE.txt'
loadings_use <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Summary/Module_FeaturePCALoadings_TCRAB_375_PRODUCTIVE.rds'

##TCRGD
eigenvectors <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Eigenvectors_No_Technical_TCRGD_PRODUCTIVE.txt'
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD'
metadata <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt'
type_receptor <- "TCRGD"
metahealth <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Healthies_ClinData.txt'
feature_assignment <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Summary/Clustered_Features_assignment_TCRGD_PRODUCTIVE_NON_IMPUTED.txt'
imputed_data <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Imputed_DATA_FINAL_SCALED_TCRGD_PRODUCTIVE.txt'
loadings_use <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Summary/Module_FeaturePCALoadings_TCRGD_200_PRODUCTIVE.rds'

##-----------------------------------------------------------------------------
## This first function looks at all day 1 measures and correlates them to the clinical data provided 
#########################################################################################################
## Function to plot difference in mean values as part of some of the other functions!
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotTukey.R')
## Lets look at the survival curve 
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/SurvivalCurve.R')
plot_survival(eigenvectors, metadata, outputdir)
## STEP 1: IDENTIFY COVARIATE:which would then be added as terms in model....
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/Covariate.R')
p_covariates <- correlate_covariates(eigenvectors, metadata, outputdir,type_receptor, metahealth)
## Perform Model selection
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/ModelSelection.R')
model_selection(eigenvectors, metadata, outputdir,type_receptor, metahealth)
## STEP 2: DIFFERENCE BETWEEN HEALTH AND DISEASE USING COVARIATES SEX AND AGE 
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotHvsD.R')
p_hd_extended <- correlate_eigenvectors_hvd(eigenvectors, metadata, outputdir, type_receptor, metahealth, "NO") 

## STEP 3: DYNAMICS WITHIN SEPSIS AND DIFFERENCE BETWEEN MORTALITY GROUPS 
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/OutcomeNoRandomEffect.R')
p_hd_extended_mort <- correlate_eigenvectors_outcome(eigenvectors, metadata, outputdir, type_receptor, "NO")
#### STEP 4: Plot the module vignette
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotModuleVignette.R')
plot_module_vignette(eigenvectors, feature_assignment, imputed_data, outputdir)
#### Do for srs
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/SRSmodel.R')
correlate_eigenvectors_srs(eigenvectors, metadata, outputdir, type_receptor, srs_assignment)


##-----------------------------------------------------------------------------
## STEP 4: PULL OUT SOME INTERESTING MODULES THAT DIFFER BETWEEN HEALTH AND SEPSIS or between mortality groups 
if(colnames(p_hd_extended)[2]!="DAY"){
	q <- p_hd_extended[p_hd_extended$DISEASE_STATUS <0.05,]
	q <- unique(q$Module)
} else {
	q <- unique(p_hd_extended$MODULE)
}	
m  <- p_hd_extended_mort[p_hd_extended_mort$p_value <0.05 & p_hd_extended_mort$Effect %in% c("LATE.DEATH", "EARLY.DEATH"),]
modules_keep <- as.character(c(q, unique(m$MODULE)))

## STEP 5: LETS CORRELATE TO SOME CLINICAL VARIABLES INCLUDING COVARIATES!
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotDvsDay1.R')
correlate_eigenvectors_t1(eigenvectors, metadata, outputdir, type_receptor, NA) 
correlate_eigenvectors_t1(eigenvectors, metadata, outputdir, type_receptor, modules_keep) 
## STEP 6: CORRELATE TO MEASURES WHERE WE HAVE MULTIPLE LONGITUDINAL MEASURES
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotContinuous.R')
correlate_eigenvectors_continuous(eigenvectors, metadata, outputdir, type_receptor)

#########################################################################################################
#### Lets nicely plot the correlation within Type 
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotCorrelationWithin.R')
plot_correlation_within(eigenvectors, outputdir, type_receptor)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotModuleVignette.R')
plot_module_vignette(eigenvectors, feature_assignment, imputed_data, outputdir)

##------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------
###########################################
### Lets now look at the correlation across a pair of features...
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotCorrelationAcross.R')
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/'
eigenvector_list <- list('/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Eigenvectors_No_Technical_TCRAB_PRODUCTIVE.txt', '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt')
type_receptor <- "BCR_TCRAB"
plot_correlation_across(eigenvector_list, outputdir, type_receptor)
####
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/'
eigenvector_list <- list('/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Eigenvectors_No_Technical_TCRAB_PRODUCTIVE.txt', '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Eigenvectors_No_Technical_TCRGD_PRODUCTIVE.txt')
type_receptor <- "TCRAB_TCRGD"
plot_correlation_across(eigenvector_list, outputdir, type_receptor)
###
outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/'
eigenvector_list <- list('/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt', '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/Eigenvectors_No_Technical_TCRGD_PRODUCTIVE.txt')
type_receptor <- "BCR_TCRGD"
plot_correlation_across(eigenvector_list, outputdir, type_receptor)## Functions to correlate eigenvectors to clinical variables!!
