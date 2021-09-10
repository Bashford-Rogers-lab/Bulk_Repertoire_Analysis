## Basic QC analysis of BCR and TCR results from RBR analysis Pipeline 
## Lauren Overend 
## lauren.overend@oriel.ox.ac.uk
## ISOTYPER R ANALYSIS 

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

option_list <- list( 
  make_option(c("-o", "--outputdir"), action="store", type="character", default="NA", help="Path to BCR/TCR Outputdir"),
  make_option(c("-s", "--samplesfilepost"), action="store", type="character", help="Sample Files Post")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE)) 

setwd("/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE")
outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_JUNE21_UPDATE/BCR/'
samplesfilepost <- '/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/LEO_SEPSIS_BCR_FINAL_post.txt'


## Location of OUTS directory from BCR/TCR Run
outputdir <- opt$o
samplesfilepost <- opt$s

# Source Auxillary Functions
my_aux_functions <- c("RFunctions")           
source_files <- list.files(my_aux_functions, "*.R$", full.names=TRUE)  # locate all .R files
for (f in source_files) {
    source(f)
}

# Create Directory for the Isotyper Plots 
dir.create(paste0(outputdir, "Plots/ISOTYPER"))

# Create Summary Matrices 
summary_isotyper(outputdir, samplesfilepost, "ALL")
summary_isotyper(outputdir, samplesfilepost, "UNPRODUCTIVE")
summary_isotyper(outputdir, samplesfilepost, "PRODUCTIVE")


# Find all analysis matrices
# Plot the analysis matrices 

files <- list.files(paste0(outputdir, "Summary"), full.names=TRUE)
files <- grep("All_raw_values", files, value=TRUE)

for (i in 1:length(files)){
	subsampled_depth <- unlist(str_split(files[i], "/"))
	subsampled_depth <- subsampled_depth[length(subsampled_depth)]
	subsampled_depth <- unlist(str_split(subsampled_depth, "_"))
	depth <- subsampled_depth[(length(subsampled_depth)-1)]
	
	productivity <- subsampled_depth[(length(subsampled_depth))]
	productivity <- gsub(".txt", "", productivity)
	
	file_use <- files[i]
	info_file <- paste0(outputdir, "Summary/isotyper_metrics_", depth, "_", productivity, ".txt")
	print(paste0("Plotting Isotyper Results for: ", info_file))
	plot_isotyper(file_use, outputdir, info_file, productivity)
	#plot_isotyper_sepsis(file_use, outputdir, info_file, productivity)
}


# Plot comparison plots of all vs functional vs non_functional 
analysis_matrices_list <- files

# get subsample depth file
files <- list.files(paste0(outputdir, "Summary"), full.names=TRUE)
files <- grep("isotyper_metrics_", files, value=TRUE)
sig <- grep("significant", files, value=TRUE)
files <- files[!files %in% sig]
info_file <- files[1]

# Plot comparison
plot_isotyper_comparison(analysis_matrices_list, outputdir, info_file)

