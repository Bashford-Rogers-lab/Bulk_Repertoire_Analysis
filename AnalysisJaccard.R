## Basic QC analysis of BCR and TCR results from RBR analysis Pipeline 
## Lauren Overend 
## lauren.overend@oriel.ox.ac.uk
# Rescomp Modules
#
library("optparse")
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggforce))
suppressMessages(library(Gviz))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(gtools))

parser <- OptionParser()
option_list <- list( 
  make_option(c("-o", "--outputdir"), action="store", type="character", default="NA", help="Path to BCR/TCR Outputdir"),
  make_option(c("-r", "--runname"), action="store", type="character", default="BCR_TCR_ANALYSIS", help="Runname for analysis [default]"),
  make_option(c("-g", "--gene"), action="store", type="character", help="GENE either TCR or BCR"),
  make_option(c("-b", "--batchfile"), action="store", type="character", default="FALSE", help="Location of Barcode Batch file"), 
  make_option(c("-t", "--task"), action="store", type="character", default=1, help="JACCARD_INDEX_CALCULATE") 
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE)) 

## Location of OUTS directory from BCR/TCR Run
results_outputdir <- opt$o
runname <- opt$r
gene <- opt$g
batch <- opt$b

#results_outputdir <- "/well/immune-rep/shared/MISEQ/TEST_PIPE/"
#runname <- "TEST" 
#gene <-"IGH" 
#batch <-"LEO_SEPSIS_BCR_ALL_layouts.txt"


# Source Auxillary Functions
my_aux_functions <- c("RFunctions/Stages1_4")           
source_files <- list.files(my_aux_functions, "*.R$", full.names=TRUE)  # locate all .R files
for (f in source_files) {
    source(f)
}

plot_dir <- paste0(results_outputdir, "/Plots")
stat_dir <- paste0(results_outputdir, "/Summary")

## CREATE PLOT DIR
if (!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

if (!dir.exists(stat_dir)){
  dir.create(stat_dir)
}

## Part Number 1: Invesitgate sample read detection and filtering
if(opt$b != "FALSE" & opt$b != "False" & opt$b != "false" ){
    
	if(opt$t == "1"){
	# Basic Jaccard
		print("Calculating JACCARD Statistics for BASIC VDJ \n")
		calculate_jaccard_matrix(results_outputdir, runname)
		print("\n Done \n")
		gc()
		calculate_jaccard_matrix_libhopcorrection(results_outputdir, runname, batch)
		print("\n Done \n")
		gc()
		calculate_jaccard_matrix_libcontam_correction(results_outputdir, runname, batch)
		print("\n Done \n")
		gc()
	}
	
	if(opt$t == "2" ){
		print("Calculating JACCARD Statistics for VDJ with 2 or more reads \n")
		#Jaccard 2 or more
		calculate_jaccard_matrix_filter(results_outputdir, runname)
		print("\n Done \n")
		gc()
		calculate_jaccard_matrix_libhopcorrection_filter(results_outputdir, runname, batch)
		print("\n Done \n")
		gc()
		calculate_jaccard_matrix_libcontam_correction_filter(results_outputdir, runname, batch)
		print("\n Done \n")
		gc()
	}
	
	#Jaccard_UMI_RAW
	if(opt$t == "3"){
		print("\n Calculating JACCARD Statistics for UMIs (RAW) \n")
		calculate_jaccard_matrix_UMI_RAW(results_outputdir, runname)
		print("\n Done \n")
		gc()
		calculate_jaccard_matrix_libhopcorrection_UMI_RAW(results_outputdir, runname, batch)
		print("\n Done \n")
		gc()
		calculate_jaccard_matrix_libcontam_correction_UMI_RAW(results_outputdir, runname, batch)
		print("\n Done \n")
		gc()
	}
	
	} else {
	if(opt$t==1){
		print("\n Calculating BASIC JACCARD Statistics on VDJ \n")
		calculate_jaccard_matrix(results_outputdir, runname)
		print("\n Done \n")
		gc()
	} 
	
	if(opt$t==2){
		print("\n Calculating JACCARD Statistics for VDJ with 2 or more reads \n")
			calculate_jaccard_matrix_filter(results_outputdir, runname)
		print("\n Done \n")
		gc()
	} 
	
	if(opt$t==3){
		print("\n Calculating JACCARD Statistics for UMIs (RAW) \n")
		calculate_jaccard_matrix_UMI_RAW(results_outputdir, runname)
		print("\n Done \n")
		gc()
	}
	
} 


## Create Summary Plots and Statistics on Jaccard
## Will only run if a batch file is provided  
locations_of_jaccard <- paste0(results_outputdir, "Summary")
path_jaccard <- paste0(results_outputdir, "Summary")
locations_of_jaccard <- list.files(locations_of_jaccard, full.names=TRUE)
locations_of_jaccard <- grep("JACCARD", locations_of_jaccard, value=TRUE)
if(batch != "FALSE" & batch != "False" & batch != "false" ){
	for(i in locations_of_jaccard){
			create_jaccard_plots(results_outputdir, i, batch)		
		}
}
	
print("Done")

## PIPELINE COMPLETE

