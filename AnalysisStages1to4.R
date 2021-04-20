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
  make_option(c("-g", "--gene"), action="store", type="character", help="GENE either TCR or BCR")
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
if(gene=="IGH"){
	visualise_filtering_bcr(results_outputdir, runname)
	visualise_constant_region_bcr(results_outputdir, runname)
	visualise_vj_usage_bcr(results_outputdir, runname)
	visualise_isoptype_cluster_bcr(results_outputdir, runname)
	calculate_jaccard_matrix(results_outputdir)
} 

if(gene=="TCR"){
	visualise_filtering_tcr(results_outputdir, runname)
	visualise_vj_usage_tcr(results_outputdir, runname)
	visualise_isoptype_cluster_tcr(results_outputdir, runname)
	calculate_jaccard_matrix(results_outputdir)
} 

## PIPELINE COMPLETE
