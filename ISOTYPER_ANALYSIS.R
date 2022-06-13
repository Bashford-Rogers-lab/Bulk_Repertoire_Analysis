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
  make_option(c("-l", "--layouts"), action="store", type="character", help="Path to batch file/layouts file"),
  make_option(c("-r", "--runimgtmutation"), action="store", type="character", help="Run IMGT Mutation")

)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE)) 

######################### 
## For practising code
setwd('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/')
samplesfilepost <- "LEO_SEPSIS_BCR_ALL_post.txt"
outputdir <- "/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/"
cluster_nodes <- 10
path_to_layout <- "LEO_SEPSIS_BCR_ALL_layouts.txt"
path_to_outputdir <- outputdir
#productivity <- "ALL"
#iso_type <- "ALL"
run_mutation <- "NO"
#######################


## Location of OUTS directory from BCR/TCR Run
outputdir <- opt$o
samplesfilepost <- opt$s
gene <- opt$g
path_to_layout <- opt$l
run_mutation<- opt$r
# Source Auxillary Functions
my_aux_functions <- c("RFunctions/Isotyper")           
source_files <- list.files(my_aux_functions, "*.R$", full.names=TRUE)  # locate all .R files
for (f in source_files) {
    source(f)
}

## DO IMGT ANALYSIS

## This can take quite some time so if it has been run previously we don't want  to rerun it just to update the summary isotyper!
if(run_mutation %like% "Y" | run_mutation %like% "y"){
	if(gene %like% "IGH" | gene %like% "BCR"){
		imgt_mutation_statistics(outputdir, cluster_nodes = 11, "ALL", path_to_layout)
		print("IMGT Mutation Summary statistics for ALL complete")
		gc()
		imgt_mutation_statistics(outputdir,  cluster_nodes = 11, "UNPRODUCTIVE", path_to_layout)
		print("IMGT Mutation Summary statistics for UNPRODUCTIVE complete")
		gc()
		imgt_mutation_statistics(outputdir, cluster_nodes = 11,  "PRODUCTIVE", path_to_layout)
		print("IMGT Mutation Summary statistics for Sepsis PRODUCTIVE complete")
		gc()
	}
}

# Create Directory for the Isotyper Plots 
dir.create(paste0(outputdir, "Plots/ISOTYPER"))
# Create Summary Matrices 
summary_isotyper(outputdir, samplesfilepost, "ALL")
print("ISOTYPER Summary statistics for  ALL complete")
gc()
summary_isotyper(outputdir, samplesfilepost, "UNPRODUCTIVE")
print("ISOTYPER Summary statistics for  UNPRODUCTIVE complete")
gc()
summary_isotyper(outputdir, samplesfilepost, "PRODUCTIVE")
print("ISOTYPER Summary statistics for  PRODUCTIVE complete")
gc()

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
	iso_type <- productivity 
	if(productivity=="ALL"){
		info_use <- info_file
	}
	print(paste0("Plotting Isotyper Results for: ", info_file))
	plot_isotyper(file_use, outputdir, info_file, productivity, depth_file)
	if(outputdir %like% "SEPSIS"){
		print("Making Sepsis Specific Plots")
		plot_isotyper_sepsis(file_use, outputdir, info_file, productivity, depth_file)
	}
}

# Plot comparison plots of all vs functional vs non_functional 
analysis_matrices_list <- files

# get subsample depth file
info <- info_use

# Plot comparison
plot_isotyper_comparison(analysis_matrices_list, outputdir, info)

####
print("Done")