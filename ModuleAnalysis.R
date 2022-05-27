## Module Reduction QC analysis of BCR and TCR results from RBR analysis Pipeline 
## Lauren Overend 
## lauren.overend@oriel.ox.ac.uk
## ISOTYPER R ANALYSIS 


option_list <- list( 
  make_option(c("-o", "--outputdir"), action="store", type="character", default="NA", help="Path to BCR/TCR Outputdir"),
  make_option(c("-s", "--samplesfilepost"), action="store", type="character", help="Sample Files Post"),
  make_option(c("-g", "--gene"), action="store", type="character", help="GENE either TCR or BCR"),
  make_option(c("-l", "--layouts"), action="store", type="character", help="Path to batch file/layouts file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE)) 
options(bitmapType='cairo-png')
setwd("/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE")

## Location of OUTS directory from BCR/TCR Run
outputdir <- opt$o
samplesfilepost <- opt$s
gene <- opt$g
path_to_layout <- opt$l

# Source Auxillary Functions
my_aux_functions <- c("RFunctions/ModuleSelection")           
source_files <- list.files(my_aux_functions, "*.R$", full.names=TRUE)  # locate all .R files
for (f in source_files) {
    source(f)
}

if(type_use %like% "TCRA" | type_use %like% "TCRB"){ 
	## Trial
	type_use <- "TCRAB"
	subsampled_depth_all<- 400
	outputdir1<- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/TCRA_CH2/'
	outputdir2<- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/TCRB_CH2/'
	productivity <- "PRODUCTIVE"
	subsampled_deptha <- 300
	subsampled_depthb  <- 400
	plot_outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/'
	meta_data <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Meta_data_for_Cohort1and2.txt'
	#new_eigenvectors <- read.delim(paste0(outputdir, "Eigenvectors_", type_use, "_", iso_type, ".txt"), sep='\t')
	cluster_immunerep(type_use, subsampled_deptha,subsampled_depthb, productivity, outputdir1, outputdir2, plot_outputdir, meta_data)
}


if(type_use %like% "TCRG" | type_use %like% "TCRD"){ 
	##FOR TCRGD
	subsampled_depth_all<- 150
	outputdir1<- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/TCRG_CH2/'
	outputdir2<- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/TCRD_CH2/'
	outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/'
	type_use <- "TCRGD"
	productivity <- "PRODUCTIVE"
	subsampled_deptha <- 300
	subsampled_depthb  <- 150
	plot_outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/'
	meta_data <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Meta_data_for_Cohort1and2.txt'
	cluster_immunerep(type_use, subsampled_deptha,subsampled_depthb, productivity, outputdir1, outputdir2, plot_outputdir, meta_data)
}


if(type_use %like% "BCR"){ 
	##FOR BCR
	subsampled_depth_all<- 1500
	outputdir1<- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/BCR_CH2/'
	outputdir2<- NA
	outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/BCR_CH2/'
	type_use <- "BCR"
	productivity <- "PRODUCTIVE"
	subsampled_deptha <- 1000
	subsampled_depthb  <- NA
	plot_outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/BCR_CH2/'
	meta_data <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Meta_data_for_Cohort1and2.txt'
	cluster_immunerep(type_use, subsampled_deptha,subsampled_depthb, productivity, outputdir1, outputdir2, plot_outputdir, meta_data)
}



