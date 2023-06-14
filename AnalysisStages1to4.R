## Basic QC analysis of BCR and TCR results from RBR analysis Pipeline 
## Lauren Overend 
## lauren.overend@oriel.ox.ac.uk


#### Set Location of librarys for SKYLAKE!!!!!!
#### THIS WILL NEED TO BE EDITED
print("Location of Packages:")
print(.libPaths())

####

####
library("optparse")
#suppressMessages(library(Gviz))
suppressMessages(library(BBmisc))
suppressMessages(library(coin))
suppressMessages(library(cowplot))
suppressMessages(library(data.table))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))
suppressMessages(library(ggforce))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))
suppressMessages(library(gridExtra))
suppressMessages(library(gtools) )
suppressMessages(library(gtools))
suppressMessages(library(optparse))
suppressMessages(library(plyr))
suppressMessages(library(purrr))
suppressMessages(library(rstatix))
suppressMessages(library(seqinr))
suppressMessages(library(stringr))
suppressMessages(library(tidyverse))


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


######################################
## For debugging 
#results_outputdir <- '/well/immune-rep/shared/MISEQ/TEST_PIPE'
#runname <- 'TEST'
#gene <- 'IGH'
#layoutsx <- 'LEO_SEPSIS_BCR_ALL_layouts.txt'
#technicals <- 'LEO_SEPSIS_BCR_ALL_technicals.txt'

#results_outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG/'
#runname <- 'TEST'
#gene <- 'TCRG'
#layoutsx <- 'LEO_SEPSIS_TCRG_ALL_layouts.txt'
#technicals <- 'LEO_SEPSIS_TCRG_ALL_technicals.txt'

#results_outputdir <- '/well/immune-rep/shared/MISEQ/TEST_LEO_ANNA_BCR/'
#runname <- 'TEST'
#gene <- 'IGH'
#layoutsx <- 'LEO_ANNA_BATCH.txt'
#technicals <- 'LEO_ANNA_tech.txt'
#########################################


### ARGUMENTS 
## Location of OUTS directory from BCR/TCR Run
results_outputdir <- opt$o
runname <- opt$r
gene <- opt$g
layoutsx <- opt$b
technicals <- opt$t

## OUTPUTDIRS
plot_dir <- paste0(results_outputdir, "/Plots")
stat_dir <- paste0(results_outputdir, "/Stats")

## CREATE PLOT DIR
if (!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

if (!dir.exists(stat_dir)){
  dir.create(stat_dir)
}
## SOURCE FUNCTIONS 
my_aux_functions <- c("RFunctions/Stages1_4")           
source_files <- list.files(my_aux_functions, "*.R$", full.names=TRUE)  # locate all .R files
for (f in source_files) {
    #print(f)
	source(f)
}

########################################
### BCR QC 
########################################
if(gene=="IGH"){
	print("Running BCR QC")
	if(opt$b != "FALSE" & opt$b != "False" & opt$b != "false"){
	
	## STAGE 1: VISUALISE READ DEPTH AND PERCENTAGE PASSED FILTERING BY LIBRARY 
	print("stage 1")
	visualise_filtering_bcr_layouts_neat(results_outputdir, runname, layoutsx, plot_dir, stat_dir, cluster_nodes = 5, gene)
	print("stage 2")
	## STAGE 2: VISUALISE READ DEPTH AND PERCENTAGE PASSED FILTERING BY SAMPLE
	visualise_filtering_bcr_layouts(results_outputdir, runname, layoutsx)
	print("stage 3")
	## STAGE 3: VISUALISE CONSTANT REGION USAGE 
	visualise_constant_region_bcr_layouts(results_outputdir,runname, layoutsx, plot_dir)
	print("stage 4")
	## STAGE 4: VISUALISE CONSTANT REGION + V GENE USAGE FOR TECHNICALS 
	compare_technicals_neat(results_outputdir, runname, technicals, plot_dir, layoutsx, stat_dir, cluster_nodes = 5)
	print("stage 5")
	## STAGE 5: VJ USAGE FOR ALL SAMPLES 
	visualise_vj_usage_bcr(results_outputdir, runname)
	print("stage 6")
	## STAGE 6: CLUSTER SIZE
	visualise_isoptype_cluster_bcr(results_outputdir, runname)
	print("stage 7")
	## RAREFACTION
	calculate_rarefaction_neat(results_outputdir, gene, plot_dir)
	print("stage 8")
	## DETECTED_VGENES
	visualise_vj_QC_neat(results_outputdir, cluster_nodes = 5, gene, plot_dir)	
}
} 

########################################
### TCR QC 
########################################
if(gene=="TCR" || gene=="TRB" || gene=="TRA"|| gene=="TRG"|| gene=="TRD"){
	if(opt$b != "FALSE" & opt$b != "False" & opt$b != "false"){
	
	print("Running TCR QC")
	## STAGE 1: VISUALISE READ DEPTH AND PERCENTAGE PASSED FILTERING BY LIBRARY 
	visualise_filtering_bcr_layouts_neat(results_outputdir, runname, layoutsx, plot_dir, stat_dir, cluster_nodes = 5, gene)
	## STAGE 2: VISUALISE READ DEPTH AND PERCENTAGE PASSED FILTERING BY SAMPLE
	visualise_filtering_tcr_layouts(path_to_outputdir=results_outputdir, run_name=runname, path_to_layout=layoutsx)
	## STAGE 3: VISUALISE CONSTANT REGION USAGE 
	visualise_constant_region_tcr_layouts(path_to_outputdir=results_outputdir, run_name=runname, path_to_layout=layoutsx)
	## STAGE 4: VISUALISE CONSTANT REGION USAGE FOR TECHNICALS 
	compare_technicals_neat_tcr(results_outputdir, runname, technicals, plot_dir, layoutsx, stat_dir, gene, cluster_nodes = 5)
	## STAGE 5: VJ USAGE FOR ALL SAMPLES
	visualise_vj_usage_tcr(results_outputdir, runname)
	## STAGE 6: CLUSTER SIZE
	visualise_isoptype_cluster_tcr(results_outputdir, runname)
	## RAREFACTION
	calculate_rarefaction_neat(results_outputdir, gene, plot_dir)
	## DETECTED_VGENES
	visualise_vj_QC_neat(results_outputdir, cluster_nodes = 5, gene, plot_dir)
} 
}
print("Done")

## PIPELINE COMPLETE

