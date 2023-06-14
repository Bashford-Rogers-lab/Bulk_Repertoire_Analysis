## Module Reduction QC analysis of BCR and TCR results from RBR analysis Pipeline 
## Lauren Overend 
## lauren.overend@oriel.ox.ac.uk
## ISOTYPER R ANALYSIS 
#module purge
#module use -a /apps/eb/dev/ivybridge/modules/all
#module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0

library(optparse)
option_list <- list( 
  make_option(c("-o", "--outputdir"), action="store", type="character", default="NA", help="Path to BCR/TCR Outputdir"),
  make_option(c("-s", "--samplesfilepost"), action="store", type="character", help="Sample Files Post"),
  make_option(c("-g", "--gene"), action="store", type="character", help="GENE either TCR or BCR"),
  make_option(c("-l", "--layouts"), action="store", type="character", help="Path to batch file/layouts file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE)) 
options(bitmapType='cairo-png')


## Location of OUTS directory from BCR/TCR Run
outputdir <- opt$o
samplesfilepost <- opt$s
gene <- opt$g
path_to_layout <- opt$l

# Source Auxillary Functions
setwd("/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE")
my_aux_functions <- c("RFunctions/ModuleSelection")           
source_files <- list.files(my_aux_functions, "*.R$", full.names=TRUE)  # locate all .R files
for (f in source_files) {
    source(f)
}

#if(type_use %like% "TCRG" | type_use %like% "TCRD"){ 
	##FOR TCRGD
	#outputdir1<- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG/'
	#outputdir2<- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRD/'
	#outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/'
	#type_use <- "TCRGD"
	#productivity <- "PRODUCTIVE"
	#subsampled_deptha <- 200
	#subsampled_depthb  <- 80
	#plot_outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD/'
	#meta_data <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Meta_data_for_Cohort1and2.txt'
	#try_imp <- "YES"
	#try_imp <- "NO"
	#try_cor <- "NO"
	#run_imp <- "NO"
	#cluster_immunerep(type_use, subsampled_deptha,subsampled_depthb, productivity, outputdir1, outputdir2, plot_outputdir, meta_data, try_imp, try_cor, run_imp)
	#print("Done TCRGD")
#}

#########################################
## REPEAT with filtered TCRG
###################################
#if(type_use %like% "TCRG" | type_use %like% "TCRD"){ 
	##FOR TCRGD
	print("Running for TCRGD")
	outputdir1<- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG_NEW/'
	outputdir2<- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRD/'
	outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW/'
	type_use <- "TCRGD"
	productivity <- "PRODUCTIVE"
	subsampled_deptha <- 50
	subsampled_depthb  <- 50
	plot_outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW/'
	meta_data <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Meta_data_for_Cohort1and2.txt'
	try_imp <- "NO"
	#try_imp <- "NO"
	try_cor <- "NO"
	run_imp <- "NO"
	cluster_immunerep(type_use, subsampled_deptha,subsampled_depthb, productivity, outputdir1, outputdir2, plot_outputdir, meta_data, try_imp, try_cor, run_imp)
	print("Done TCRGD")
#}


#if(type_use %like% "TCRA" | type_use %like% "TCRB"){ 
	## Trial
	print("Running for TCRAB")
	type_use <- "TCRAB"
	#subsampled_depth_all<- 400.375
	outputdir1<- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRA/'
	outputdir2<- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRB/'
	productivity <- "PRODUCTIVE"
	subsampled_deptha <- 300
	subsampled_depthb  <- 300
	plot_outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/'
	meta_data <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Meta_data_for_Cohort1and2.txt'
	#try_imp <- "NO"
	try_imp <- "NO"
	try_cor <- "NO"
	run_imp <- "NO"
	cluster_immunerep(type_use, subsampled_deptha,subsampled_depthb, productivity, outputdir1, outputdir2, plot_outputdir, meta_data,try_imp, try_cor, run_imp)
	print("Done TCRAB")
#}


#if(type_use %like% "BCR"){ 
	##FOR BCR
	#print("Running for BCR")
	#subsampled_depth_all<- 1500
	#outputdir1<- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/'
	#outputdir2<- NA
	#outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/'
	#type_use <- "BCR"
	#productivity <- "PRODUCTIVE"
	#subsampled_deptha <- 1500
	#subsampled_depthb  <- NA
	#plot_outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/'
	#meta_data <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Meta_data_for_Cohort1and2.txt'
	#try_imp <- "YES"
	#try_imp <- "NO"
	#try_cor <- "NO"
	#run_imp <- "NO"
	#cluster_immunerep(type_use, subsampled_deptha,subsampled_depthb, productivity, outputdir1, outputdir2, plot_outputdir, meta_data, try_imp, try_cor, run_imp)
	#print("Done BCR")
#}



