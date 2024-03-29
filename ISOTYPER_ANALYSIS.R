## Basic QC analysis of BCR and TCR results from RBR analysis Pipeline 
## Lauren Overend 
## lauren.overend@oriel.ox.ac.uk
## ISOTYPER R ANALYSIS 


### Set Location of librarys for SKYLAKE!!!!!!
#### THIS WILL NEED TO BE EDITED
print("Location of Packages:")
print(.libPaths())####
######

library("optparse")
suppressMessages(library(corrplot))
suppressMessages(library(data.table))
suppressMessages(library(doParallel))
suppressMessages(library(dplyr))
suppressMessages(library(foreach))
suppressMessages(library(ggforce))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrastr))
suppressMessages(library(gridExtra))
suppressMessages(library(gtools))
#suppressMessages(library(Gviz))
suppressMessages(library(Hmisc))
suppressMessages(library(matrixStats))
suppressMessages(library(moments))
suppressMessages(library(mousetrap))
suppressMessages(library(Peptides))
suppressMessages(library(plot3D))
suppressMessages(library(plyr))
suppressMessages(library(purrr))
suppressMessages(library(reshape2))
#suppressMessages(library(ShortRead))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))


option_list <- list( 
  make_option(c("-o", "--outputdir"), action="store", type="character", default="NA", help="Path to BCR/TCR Outputdir"),
  make_option(c("-s", "--samplesfilepost"), action="store", type="character", help="Sample Files Post"),
  make_option(c("-g", "--gene"), action="store", type="character", help="GENE either TCR or BCR"),
  make_option(c("-l", "--layouts"), action="store", type="character", help="Path to batch file/layouts file"),
  make_option(c("-r", "--runimgtmutation"), action="store", type="character", help="Run IMGT Mutation"),
  make_option(c("-n", "--normalise"), action="store", type="character", help="Do I normalise V gene Usage"),
  make_option(c("-e", "--excludesamples"), action="store", default=NA, type="character", help="Samples to Exclude")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE)) 

setwd('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/')
######################### 
## For practising code

## Location of OUTS directory from BCR/TCR Run
cluster_nodes <- 10
outputdir <- opt$o
samplesfilepost <- opt$s
gene <- opt$g
path_to_layout <- opt$l
run_mutation<- opt$r
exclude_samples <- opt$e
normalise <- opt$n

#samplesfilepost <- "LEO_SEPSIS_TCRG_ALL_post.txt"
#outputdir <- "/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG_NEW/"
#cluster_nodes <- 10
#path_to_layout <- "LEO_SEPSIS_TCRG_ALL_layouts.txt"
#path_to_outputdir <- outputdir
#productivity <- "PRODUCTIVE"
#iso_type <- "PRODUCTIVE"
#run_mutation <- "Y"
#gene <- "TCR"
#exclude_samples <- 'LEO_Samples_Exclude.txt'

#samplesfilepost <- "R4RA_Anna_BCR_post_V0_1.txt"
#outputdir <- "/well/immune-rep/shared/MISEQ/R4RA_FINAL/BCR/"
#cluster_nodes <-5
#path_to_layout <- "R4RA_Anna_Batch_file_BCR.txt"
#path_to_outputdir <- outputdir
#productivity <- "PRODUCTIVE"
#iso_type <- "PRODUCTIVE"
#run_mutation <- "Y"
#gene <- "IGH"
#normalise <- "FALSE"
#exclude_samples <- NA
#class_vdj <- "IGH"
#all_class <- "TCRA"
#exclude_samples <- 'LEO_Samples_Exclude.txt'

# Source Auxillary Functions
my_aux_functions <- c("RFunctions/Isotyper")           
source_files <- list.files(my_aux_functions, "*.R$", full.names=TRUE)  # locate all .R files
for (f in source_files) {
    source(f)
}


############### Make sure directories are there 
plot_dir <- paste0(outputdir, "/Plots")
stat_dir <- paste0(outputdir, "/Stats")

## CREATE PLOT DIR
if (!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

if (!dir.exists(stat_dir)){
  dir.create(stat_dir)
}
############### 

if(normalise == "TRUE"){
	print("WARNING Normalised IGHV gene usage will be calculated!!")
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
if(dir.exists(paste0(outputdir, "Plots/ISOTYPER"))==FALSE){
		dir.create(paste0(outputdir, "Plots/ISOTYPER"))
}	

# Create Summary Matrices 
summary_isotyper(outputdir, samplesfilepost, "ALL", exclude_samples, path_to_layout, normalise, run_mutation)
print("ISOTYPER Summary statistics for  ALL complete")
gc()
summary_isotyper(outputdir, samplesfilepost, "UNPRODUCTIVE", exclude_samples, path_to_layout, normalise, run_mutation)
print("ISOTYPER Summary statistics for  UNPRODUCTIVE complete")
gc()
summary_isotyper(outputdir, samplesfilepost, "PRODUCTIVE", exclude_samples, path_to_layout, normalise, run_mutation)
print("ISOTYPER Summary statistics for  PRODUCTIVE complete")
gc()

# Find all analysis matrices
# Plot the analysis matrices from isotyper

files <- list.files(paste0(outputdir, "Summary"), full.names=TRUE)
files <- grep("isotyper_metrics_filtered_FINAL_ALL", files, value=TRUE)
files <- grep("TENSOR_FORMAT", files, invert=TRUE, value=TRUE)

## in case the unproductive file is emtpy 
file_skip <- c()

for (i in 1:length(files)){
	subsampled_depth <- unlist(str_split(files[i], "/"))
	subsampled_depth <- subsampled_depth[length(subsampled_depth)]
	subsampled_depth <- unlist(str_split(subsampled_depth, "_"))
	depth <- subsampled_depth[(length(subsampled_depth)-1)]
	productivity <- subsampled_depth[(length(subsampled_depth))]
	productivity <- gsub(".txt", "", productivity)
	file_use <- files[i]
	info <- paste0(outputdir, "Summary/isotyper_metrics_summary_stats", depth, "_", productivity, ".txt")
	depth_file <- paste0(outputdir, "Summary/Read_Depths_", productivity, ".txt")
	iso_type <- productivity 
	
	skip_to_next <- FALSE
  
	# Note that print(b) fails since b doesn't exist
	tryCatch(read.delim(file_use, sep="\t"), error = function(e) { 
	skip_to_next <<- TRUE
	})
	
    if(skip_to_next) { 
		file_skip <- c(file_skip, file_use)
		next 
	} 
  
	## Test if file is empty
	if(productivity=="ALL"){
		info_use <- info
	}
	print(paste0("Plotting Isotyper Results for: ", info))
	plot_isotyper(file_use, outputdir, info, productivity, depth_file)
	if(outputdir %like% "SEPSIS"){
		print("Making Sepsis Specific Plots")
		suppressMessages(plot_isotyper_sepsis(file_use, outputdir, info, productivity, depth_file))
	}
	layouts <- read.delim(path_to_layout, header=TRUE)
	if(any(layouts$SampleID!=layouts$Barcode)){
		plot_isotyper_groups(file_use, outputdir, info, productivity, depth_file)
	}
}

print("Plot Comparison plots")
# Plot comparison plots of all vs functional vs non_functional 
analysis_matrices_list <- files[!files %in% file_skip]

# get subsample depth file
info <- info_use

# Plot comparison
if(length(file_skip)==0){
	plot_isotyper_comparison(analysis_matrices_list, outputdir, info)
}
####
print("Done")