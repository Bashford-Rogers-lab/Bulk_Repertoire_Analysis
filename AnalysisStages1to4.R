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
  make_option(c("-b", "--batchfile"), action="store", type="character", default="FALSE", help="Path to Layouts file"), 
  make_option(c("-t", "--technicalreplicates"), action="store", type="character", default="FALSE", help="Path to technical replicates file") 
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE)) 

## Location of OUTS directory from BCR/TCR Run
results_outputdir <- opt$o
runname <- opt$r
gene <- opt$g
layoutsx <- opt$b
technicals <- opt$t
# Source Auxillary Functions
my_aux_functions <- c("RFunctions/Stages1_4")           
source_files <- list.files(my_aux_functions, "*.R$", full.names=TRUE)  # locate all .R files
for (f in source_files) {
    source(f)
}

## Part Number 1: Invesitgate sample read detection and filtering
if(gene=="IGH"){
	if(opt$b != "FALSE" & opt$b != "False" & opt$b != "false"){
			visualise_filtering_bcr_layouts(path_to_outputdir=results_outputdir, run_name=runname, path_to_layout=layoutsx)
			visualise_constant_region_bcr_layouts(path_to_outputdir=results_outputdir, run_name=runname, path_to_layout=layoutsx)
	} else {
			visualise_filtering_bcr(path_to_outputdir=results_outputdir, run_name=runname)
			visualise_constant_region_bcr(results_outputdir, runname)
	}
	
	visualise_vj_usage_bcr(results_outputdir, runname)
	visualise_isoptype_cluster_bcr(results_outputdir, runname)
	calculate_rarefaction(results_outputdir)
	visualise_vj_QC(results_outputdir)
	if(!is.na(opt$t)){
		compare_technicals(results_outputdir, runname, technicals)
	} 
	
} 

if(gene=="TCR" || gene=="TRB" || gene=="TRA"|| gene=="TRG"|| gene=="TRD"){
	if(opt$b != "FALSE" & opt$b != "False" & opt$b != "false"){
			visualise_filtering_tcr_layouts(path_to_outputdir=results_outputdir, run_name=runname, path_to_layout=layoutsx)
			visualise_constant_region_tcr_layouts(path_to_outputdir=results_outputdir, run_name=runname, path_to_layout=layoutsx)
	} else {
			visualise_filtering_tcr(path_to_outputdir=results_outputdir, run_name=runname)
			visualise_constant_region_tcr(results_outputdir, runname)
	}
	visualise_vj_usage_tcr(results_outputdir, runname)
	visualise_isoptype_cluster_tcr(results_outputdir, runname)
	calculate_rarefaction(results_outputdir)
	visualise_vj_QC(results_outputdir)
	if(!is.na(opt$t)){
		compare_technicals(results_outputdir, runname, technicals)
	} 
} 

print("Done")

## PIPELINE COMPLETE

