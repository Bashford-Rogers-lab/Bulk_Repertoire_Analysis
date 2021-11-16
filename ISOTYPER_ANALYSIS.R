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
suppressMessages(library(purrr))
suppressMessages(library(reshape2))
suppressMessages(library(Hmisc))
suppressMessages(library(corrplot))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(matrixStats))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrastr))
suppressMessages(library(ggpubr))

option_list <- list( 
  make_option(c("-o", "--outputdir"), action="store", type="character", default="NA", help="Path to BCR/TCR Outputdir"),
  make_option(c("-s", "--samplesfilepost"), action="store", type="character", help="Sample Files Post"),
  make_option(c("-g", "--gene"), action="store", type="character", help="GENE either TCR or BCR"),
  make_option(c("-l", "--layouts"), action="store", type="character", help="Path to batch file/layouts file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE)) 

setwd("/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE")
samplesfilepost <- "/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/LEO_Samples_IMMUNOAGEING_TCR_INT_post.txt"
outputdir <- "/well/immune-rep/shared/MISEQ/IMMUNOAGEING/TCR_INT/"
cluster_nodes <- 10
path_to_layout <- '/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/BatchingLayouts_TCR_INT.txt'
path_to_outputdir <- outputdir
productivity <- "ALL"
iso_type <- "ALL"

## Location of OUTS directory from BCR/TCR Run
outputdir <- opt$o
samplesfilepost <- opt$s
gene <- opt$g
path_to_layout <- opt$l

# Source Auxillary Functions
my_aux_functions <- c("RFunctions/Isotyper")           
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
# Plot the analysis matrices from isotyper

files <- list.files(paste0(outputdir, "Summary"), full.names=TRUE)
files <- grep("isotyper_metrics_filtered_FINAL_METRICS_", files, value=TRUE)

for (i in 1:length(files)){
	subsampled_depth <- unlist(str_split(files[i], "/"))
	subsampled_depth <- subsampled_depth[length(subsampled_depth)]
	subsampled_depth <- unlist(str_split(subsampled_depth, "_"))
	depth <- subsampled_depth[(length(subsampled_depth)-1)]
	productivity <- subsampled_depth[(length(subsampled_depth))]
	productivity <- gsub(".txt", "", productivity)
	file_use <- files[i]
	info_file <- paste0(outputdir, "Summary/isotyper_metrics_summary_stats", depth, "_", productivity, ".txt")
	depth_file <- paste0(outputdir, "Summary/Read_Depths_", productivity, ".txt")
	if(productivity=="ALL"){
		info_use <- info_file
	}
	print(paste0("Plotting Isotyper Results for: ", info_file))
	plot_isotyper(file_use, outputdir, info_file, productivity, depth_file)
	#plot_isotyper_sepsis(file_use, outputdir, info_file, productivity)
}

# Plot comparison plots of all vs functional vs non_functional 
analysis_matrices_list <- files

# get subsample depth file
info <- info_use

# Plot comparison
plot_isotyper_comparison(analysis_matrices_list, outputdir, info)

# Plot SHM summary from IMGT
if(gene=="IGH"){
	imgt_mutation_statistics(outputdir, cluster_nodes = 11, "ALL", path_to_layout)
	imgt_mutation_statistics(outputdir,  cluster_nodes = 11, "UNPRODUCTIVE", path_to_layout)
	imgt_mutation_statistics(outputdir, cluster_nodes = 11,  "PRODUCTIVE", path_to_layout)
}

####
