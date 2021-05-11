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
  make_option(c("-b", "--batchfile"), action="store", type="character", default="FALSE", help="Location of Barcode Batch file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE)) 

## Location of OUTS directory from BCR/TCR Run
results_outputdir <- opt$o
runname <- opt$r
gene <- opt$g
# Source Auxillary Functions
my_aux_functions <- c("RFunctions/")           
source_files <- list.files(my_aux_functions, "*.R$", full.names=TRUE)  # locate all .R files
for (f in source_files) {
    source(f)
}

## Part Number 1: Invesitgate sample read detection and filtering
if(opt$b != "FALSE"){
	calculate_jaccard_matrix(results_outputdir, runname)
	print("\n Done \n")
	gc()
	calculate_jaccard_matrix_libhopcorrection(results_outputdir, runname, opt$b)
	print("\n Done \n")
	gc()
	calculate_jaccard_matrix_libcontam_correction(results_outputdir, runname, opt$b)
	print("\n Done \n")
	gc()
	calculate_jaccard_matrix_filter(results_outputdir, runname)
	print("\n Done \n")
	gc()
	calculate_jaccard_matrix_libhopcorrection_filter(results_outputdir, runname, opt$b)
	print("\n Done \n")
	gc()
	calculate_jaccard_matrix_libcontam_correction_filter(results_outputdir, runname, opt$b)
	print("\n Done \n")
	gc()
	
	## RUN on UMIs
	
	calculate_jaccard_matrix_UMI(results_outputdir, runname)
	print("\n Done \n")
	gc()
	calculate_jaccard_matrix_libhopcorrection_UMI(results_outputdir, runname, opt$b)
	print("\n Done \n")
	gc()
	calculate_jaccard_matrix_libcontam_correction_UMI(results_outputdir, runname, opt$b)
	print("\n Done \n")
	gc()
	calculate_jaccard_matrix_filter_UMI(results_outputdir, runname)
	print("\n Done \n")
	gc()
	calculate_jaccard_matrix_libhopcorrection_filter_UMI(results_outputdir, runname, opt$b)
	print("\n Done \n")
	gc()
	calculate_jaccard_matrix_libcontam_correction_filter_UMI(results_outputdir, runname, opt$b)
	print("\n Done \n")
	gc()
	
	} else {
	calculate_jaccard_matrix(results_outputdir, runname)
	print("\n Done \n")
	gc()
	calculate_jaccard_matrix_filter(results_outputdir, runname)
	print("\n Done \n")
	gc()
	
	## RUN on UMIs
	
	calculate_jaccard_matrix_UMI(results_outputdir, runname)
	print("\n Done \n")
	gc()
	calculate_jaccard_matrix_filter_UMI(results_outputdir, runname)
	print("\n Done \n")
	gc()
	
	}
	
} 


## Create Summary Plots and Statistics on Jaccard 
locations_of_jaccard <- paste0(results_outputdir, "Summary")
path_jaccard <- paste0(results_outputdir, "Summary")
locations_of_jaccard <- list.files(locations_of_jaccard, full.names=TRUE)
locations_of_jaccard <- grep("JACCARD", locations_of_jaccard, value=TRUE)
for(i in locations_of_jaccard){
		create_jaccard_plots(results_outputdir, i, opt$b)
	}
	
print("Done")

## PIPELINE COMPLETE
