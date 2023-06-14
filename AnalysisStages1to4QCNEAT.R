## Code to generate neat QC plots for BCR  for thesis 
# Lauren Overend
## November 2022
## lauren.overend@oriel.ox.ac.uk

module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0

#### Packages 
library(tidyverse)
library(ggplot2)
library(foreach)
library(doParallel)
library(gridExtra)
library(stringr)
library(tidyverse)
library(ggforce)
library(Gviz)
library(cowplot)
library(gtools)
library(data.table)
library(ggpubr)
library(ggrepel)
	
	
### Variables
setwd('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE')

### BCR
path_to_outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR'
run_name <- 'SEPSIS_NEAT'
gene <- 'IGH'
path_to_layout <- 'LEO_SEPSIS_BCR_ALL_layouts.txt'
technical_replicates_file <- 'LEO_SEPSIS_BCR_ALL_technicals.txt'
plot_dir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/THESIS_PLOTS'
stat_dir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/THESIS_STATS'
cluster_nodes = 5
type_use <- "BCR"
iso_type <- "PRODUCTIVE"
if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
if (!dir.exists(stat_dir)) {dir.create(stat_dir)}
jaccard_matrix <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/JACCARDMATRIX_BASIC_AllSAMPLES_LEO_SEPSIS.txt'
chain <- "BCR"

source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/QC_readdepth.R')
visualise_filtering_bcr_layouts_neat(path_to_outputdir, run_name, path_to_layout, plot_dir, stat_dir, cluster_nodes = 5, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/GeneUsage.R')
compare_technicals_neat(path_to_outputdir, run_name, technical_replicates_file, plot_dir, path_to_layout, stat_dir, cluster_nodes = 5)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/PlotJaccardNeat.R')
create_jaccard_plots(paste0(path_to_outputdir, "/"), jaccard_matrix, path_to_layout, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/DetectedGenes.R')
visualise_vj_QC_neat(paste0(path_to_outputdir, "/"), cluster_nodes = 5, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/RarefectionNeat.R')
calculate_rarefaction_neat(paste0(path_to_outputdir, "/"), chain)



#################################################
## TCRA
path_to_outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRA'
run_name <- 'SEPSIS_NEAT'
gene <- 'TCR'
path_to_layout <- 'LEO_SEPSIS_TCRA_ALL_layouts.txt'
technical_replicates_file <- 'LEO_SEPSIS_TCRA_ALL_technicals.txt'
plot_dir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRA/THESIS_PLOTS'
stat_dir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRA/THESIS_STATS'
cluster_nodes = 5
type_use <- "TCRA"
iso_type <- "PRODUCTIVE"
if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
if (!dir.exists(stat_dir)) {dir.create(stat_dir)}
jaccard_matrix <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRA/Summary/JACCARDMATRIX_BASIC_AllSAMPLES_LEO_SEPSIS.txt'
chain="TCRA"

source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/QC_readdepth.R')
a_list <- visualise_filtering_bcr_layouts_neat(path_to_outputdir, run_name, path_to_layout, plot_dir, stat_dir, cluster_nodes = 5, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/GeneUsage.R')
a2_list <- compare_technicals_neat_tcr(path_to_outputdir, run_name, technical_replicates_file, plot_dir, path_to_layout, stat_dir, chain, cluster_nodes = 5)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/PlotJaccardNeat.R')
Ja <- create_jaccard_plots(paste0(path_to_outputdir, "/"), jaccard_matrix, path_to_layout, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/DetectedGenes.R')
a3_list <- visualise_vj_QC_neat(paste0(path_to_outputdir, "/"), cluster_nodes = 5, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/RarefectionNeat.R')
rar1 <- calculate_rarefaction_neat(paste0(path_to_outputdir, "/"), chain)

#################################################
## TCRB
path_to_outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRB'
run_name <- 'SEPSIS_NEAT'
gene <- 'TCR'
path_to_layout <- 'LEO_SEPSIS_TCRB_ALL_layouts.txt'
technical_replicates_file <- 'LEO_SEPSIS_TCRB_ALL_technicals.txt'
plot_dir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRB/THESIS_PLOTS'
stat_dir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRB/THESIS_STATS'
cluster_nodes = 5
type_use <- "TCRB"
iso_type <- "PRODUCTIVE"
if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
if (!dir.exists(stat_dir)) {dir.create(stat_dir)}
jaccard_matrix <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRB/Summary/JACCARDMATRIX_BASIC_AllSAMPLES_LEO_SEPSIS.txt'
chain="TCRB"

source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/QC_readdepth.R')
b_list <- visualise_filtering_bcr_layouts_neat(path_to_outputdir, run_name, path_to_layout, plot_dir, stat_dir, cluster_nodes = 5, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/GeneUsage.R')
b2_list <- compare_technicals_neat_tcr(path_to_outputdir, run_name, technical_replicates_file, plot_dir, path_to_layout, stat_dir, chain, cluster_nodes = 5)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/PlotJaccardNeat.R')
Jb <- create_jaccard_plots(paste0(path_to_outputdir, "/"), jaccard_matrix, path_to_layout, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/DetectedGenes.R')
b3_list <- visualise_vj_QC_neat(paste0(path_to_outputdir, "/"), cluster_nodes = 5, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/RarefectionNeat.R')
rar2 <- calculate_rarefaction_neat(paste0(path_to_outputdir, "/"), chain)


#################################################
## TCRG
path_to_outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG'
run_name <- 'SEPSIS_NEAT'
gene <- 'TCR'
path_to_layout <- 'LEO_SEPSIS_TCRG_ALL_layouts.txt'
technical_replicates_file <- 'LEO_SEPSIS_TCRG_ALL_technicals.txt'
plot_dir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG/THESIS_PLOTS'
stat_dir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG/THESIS_STATS'
cluster_nodes = 5
type_use <- "TCRG"
iso_type <- "PRODUCTIVE"
if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
if (!dir.exists(stat_dir)) {dir.create(stat_dir)}
jaccard_matrix <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG/Summary/JACCARDMATRIX_BASIC_AllSAMPLES_LEO_SEPSIS.txt'
chain="TCRG"

source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/QC_readdepth.R')
g_list <- visualise_filtering_bcr_layouts_neat(path_to_outputdir, run_name, path_to_layout, plot_dir, stat_dir, cluster_nodes = 5, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/GeneUsage.R')
g2_list <- compare_technicals_neat_tcr(path_to_outputdir, run_name, technical_replicates_file, plot_dir, path_to_layout, stat_dir, chain, cluster_nodes = 5)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/PlotJaccardNeat.R')
Jc <- create_jaccard_plots(paste0(path_to_outputdir, "/"), jaccard_matrix, path_to_layout, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/DetectedGenes.R')
g3_list <- visualise_vj_QC_neat(paste0(path_to_outputdir, "/"), cluster_nodes = 5, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/RarefectionNeat.R')
rar3 <- calculate_rarefaction_neat(paste0(path_to_outputdir, "/"), chain)

#################################################
## TCRD
path_to_outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRD'
run_name <- 'SEPSIS_NEAT'
gene <- 'TCR'
path_to_layout <- 'LEO_SEPSIS_TCRD_ALL_layouts.txt'
technical_replicates_file <- 'LEO_SEPSIS_TCRD_ALL_technicals.txt'
plot_dir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRD/THESIS_PLOTS'
stat_dir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRD/THESIS_STATS'
cluster_nodes = 5
type_use <- "TCRD"
iso_type <- "PRODUCTIVE"
if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
if (!dir.exists(stat_dir)) {dir.create(stat_dir)}
jaccard_matrix <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRD/Summary/JACCARDMATRIX_BASIC_AllSAMPLES_LEO_SEPSIS.txt'
chain="TCRD" 

source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/QC_readdepth.R')
d_list <- visualise_filtering_bcr_layouts_neat(path_to_outputdir, run_name, path_to_layout, plot_dir, stat_dir, cluster_nodes = 5, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/GeneUsage.R')
d2_list <- compare_technicals_neat_tcr(path_to_outputdir, run_name, technical_replicates_file, plot_dir, path_to_layout, stat_dir, chain, cluster_nodes = 5)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/PlotJaccardNeat.R')
Jd <- create_jaccard_plots(paste0(path_to_outputdir, "/"), jaccard_matrix, path_to_layout, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/DetectedGenes.R')
d3_list <- visualise_vj_QC_neat(paste0(path_to_outputdir, "/"), cluster_nodes = 5, chain)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/RarefectionNeat.R')
rar4 <- calculate_rarefaction_neat(paste0(path_to_outputdir, "/"), chain)

############ Merge all the QC Plots 
pdf(paste0('/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCR_thesis_plots/Filtering_QC.pdf'), width=15, height=12)
plot(plot_grid(a_list, b_list, g_list, d_list, ncol=1, align = "v", axis = "btlr"))
dev.off()

pdf(paste0('/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCR_thesis_plots/Vgene_usage_technical_QC.pdf'), width=13, height=13)
plot(plot_grid(a2_list[[1]], b2_list[[1]], g2_list[[1]], d2_list[[1]], ncol=1, align = "v", rel_heights=c(1,1,0.5,0.5), labels="AUTO"))
dev.off()

pdf(paste0('/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCR_thesis_plots/Constant_usage_technical_QC.pdf'), width=10, height=10)
plot(plot_grid(a2_list[[2]], b2_list[[2]], g2_list[[2]], d2_list[[2]], ncol=1, align = "v", rel_heights=c(1,1,1,1), labels="AUTO", axis = "btlr"))
dev.off()

pdf(paste0('/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCR_thesis_plots/DetectedVgenes.pdf'), width=10, height=10)
plot(plot_grid(a3_list, b3_list, g3_list, d3_list, ncol=1, align = "v", rel_heights=c(1,1,1,1), labels="AUTO", axis = "btlr"))
dev.off()

pdf(paste0('/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCR_thesis_plots/Rarefaction.pdf'), width=14.25, height=10)
plot(plot_grid(rar1, rar2, rar3, rar4, ncol=1, align = "v", rel_heights=c(1,1,1,1), labels="AUTO", axis = "btlr"))
dev.off()

pdf(paste0('/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCR_thesis_plots/Jaccard.pdf'), width=10, height=13)
plot(plot_grid(Ja, Jb, Jc, Jd, ncol=1, align = "v", rel_heights=c(1,1,1,1), labels="AUTO", axis = "btlr"))
dev.off()

########################################################################################################################
## Combined TCRAB
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/summary_imputation.R')
impt <- summary_imputation('/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB', '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Plots', '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Summary', "TCRAB", "PRODUCTIVE")

## COMBINED TCRGD 
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ThesisSummaryPlots/summary_imputation.R')
impt2 <- summary_imputation('/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW', '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW/Plots', '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRGD_NEW/Summary', "TCRGD", "PRODUCTIVE")

if(length(impt) > 3){
pdf(paste0('/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCR_thesis_plots/Imputation.pdf'), width=14, height=7)
plot(plot_grid(impt[[4]], impt[[5]], impt[[6]], impt[[7]], impt2[[4]], impt2[[5]], impt2[[6]], impt2[[7]], ncol=4, align = "v", axis = "btlr", labels="AUTO"))
dev.off()
}

pdf(paste0('/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCR_thesis_plots/Imputation_Method.pdf'), width=14, height=8)
plot(plot_grid(impt[[1]], impt2[[1]], ncol=2, align = "v", axis = "btlr", labels="AUTO"))
dev.off()