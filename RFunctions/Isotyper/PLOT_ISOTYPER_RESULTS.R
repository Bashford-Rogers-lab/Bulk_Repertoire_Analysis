suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggforce))
#suppressMessages(library(Gviz))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(gtools))
suppressMessages(library(ggrastr))
suppressMessages(library(purrr))
suppressMessages(library(ggpubr))
suppressMessages(library(doParallel))
#suppressMessages(library(ShortRead))

plot_isotyper_sepsis <- function(file_use=file_use, outputdir, info, iso_type, depth_file){
	
	analysis_matrices <- file_use
	analysis_matrices <- read.delim(analysis_matrices, sep="\t")
	info <- read.delim(info, sep="\t")
	depth_file <- read.delim(depth_file, sep="\t")
	
	depth_file$SampleIDforDepths <- gsub("_unproductive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("_productive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("BCR_", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("TCR_", "", depth_file$SampleIDforDepths)
	info$Metric <- gsub("__TCR", "", info$Metric)
	info$Metric <- gsub("__TCRA", "", info$Metric)
	info$Metric <- gsub("__TCRB", "", info$Metric)
	info$Metric <- gsub("__TCRG", "", info$Metric)
	info$Metric <- gsub("__TCRD", "", info$Metric)
	info$Metric <- gsub("__TRAC", "", info$Metric)
	info$Metric <- gsub("__TRBC", "", info$Metric)
	info$Metric <- gsub("__TRGC", "", info$Metric)
	info$Metric <- gsub("__TRDC", "", info$Metric)
	
	meta_data <- read.delim('/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Meta_data_for_Cohort1and2.txt', sep='\t', header=TRUE)
	
	cores <- 10
	cl <- parallel::makeCluster(cores)
	doParallel::registerDoParallel(cl)
	foreach(i = 1:length(colnames(analysis_matrices)), .export='sepsis_plots',   .packages=c("ggplot2", "data.table", "stringr", "ggrastr", "ggpubr", "cowplot", "gridExtra")) %dopar% {
				sepsis_plots(analysis_matrices, iso_type, outputdir, i, depth_file, info, meta_data)
	}
	stopCluster(cl)
	print("Done Sepsis Plots")
}

plot_isotyper <- function(file_use=file_use, outputdir, info, iso_type, depth_file){
	analysis_matrices <- file_use
	analysis_matrices <- read.delim(analysis_matrices, sep="\t")
	info <- read.delim(info, sep="\t")
	depth_file <- read.delim(depth_file, sep="\t")
	
	depth_file$SampleIDforDepths <- gsub("_unproductive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("_productive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("BCR_", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("TCR_", "", depth_file$SampleIDforDepths)
	info$Metric <- gsub("__TCR", "", info$Metric)
	info$Metric <- gsub("__TCRA", "", info$Metric)
	info$Metric <- gsub("__TCRB", "", info$Metric)
	info$Metric <- gsub("__TCRG", "", info$Metric)
	info$Metric <- gsub("__TCRD", "", info$Metric)
	info$Metric <- gsub("__TRAC", "", info$Metric)
	info$Metric <- gsub("__TRBC", "", info$Metric)
	info$Metric <- gsub("__TRGC", "", info$Metric)
	info$Metric <- gsub("__TRDC", "", info$Metric)
	## Parellise the plotting function!
	
	cores <- 10
	cl <- parallel::makeCluster(cores)
	doParallel::registerDoParallel(cl)
	foreach(i = 1:length(colnames(analysis_matrices)), .export='basic_plot',   .packages=c("ggplot2", "data.table", "ggrastr", "ggpubr", "cowplot", "gridExtra")) %dopar% {
				basic_plot(analysis_matrices, iso_type, outputdir, i, depth_file, info)
	}
	stopCluster(cl)
	print("Done Basic Plots")
}


plot_isotyper_groups <- function(file_use=file_use, outputdir, info, iso_type, depth_file){
	
	analysis_matrices <- file_use
	analysis_matrices <- read.delim(analysis_matrices, sep="\t")
	info <- read.delim(info, sep="\t")
	depth_file <- read.delim(depth_file, sep="\t")
	
	depth_file$SampleIDforDepths <- gsub("_unproductive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("_productive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("BCR_", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("TCR_", "", depth_file$SampleIDforDepths)
	info$Metric <- gsub("__TCR", "", info$Metric)
	info$Metric <- gsub("__TCRA", "", info$Metric)
	info$Metric <- gsub("__TCRB", "", info$Metric)
	info$Metric <- gsub("__TCRG", "", info$Metric)
	info$Metric <- gsub("__TCRD", "", info$Metric)
	info$Metric <- gsub("__TRAC", "", info$Metric)
	info$Metric <- gsub("__TRBC", "", info$Metric)
	info$Metric <- gsub("__TRGC", "", info$Metric)
	info$Metric <- gsub("__TRDC", "", info$Metric)
	
	cores <- 10
	cl <- parallel::makeCluster(cores)
	doParallel::registerDoParallel(cl)
	foreach(i = 1:length(colnames(analysis_matrices)), .export='group_plot', .packages=c("ggplot2", "data.table", "stringr", "ggrastr", "ggpubr", "cowplot", "gridExtra")) %dopar% {
		group_plot(analysis_matrices, iso_type, outputdir, i, depth_file, info)
	}
	stopCluster(cl)
	print("Done Grouped Plots")	
}


plot_isotyper_comparison <- function(analysis_matrices_list, outputdir, info){
	
	dfiles <- list.files(paste0(outputdir, "Summary"), full.names=TRUE)
	dfiles <- grep("Read_Depths_", dfiles, value=TRUE)
	dprod <- grep("_PRODUCTIVE", dfiles, value=TRUE)
	dprod <- read.delim(dprod, sep="\t")
	duprod <- grep("_UNPRODUCTIVE", dfiles, value=TRUE)
	duprod <- read.delim(duprod, sep="\t")
	dall <- grep("_ALL", dfiles, value=TRUE)
	dall <- read.delim(dall, sep="\t")
	dprod$SampleIDforDepths <- gsub("_productive", "", dprod$SampleIDforDepths)
	duprod$SampleIDforDepths <- gsub("_unproductive", "", duprod$SampleIDforDepths)
	
	analysis_matrices1 <- read.delim(analysis_matrices_list[1], sep="\t")
	analysis_matrices1$productivity <- analysis_matrices_list[1]
	analysis_matrices1$sample <- row.names(analysis_matrices1)

	analysis_matrices2 <- read.delim(analysis_matrices_list[2], sep="\t")
	analysis_matrices2$productivity <- analysis_matrices_list[2]
	analysis_matrices2$sample <- row.names(analysis_matrices2)
	
	analysis_matrices3 <- read.delim(analysis_matrices_list[3], sep="\t")
	analysis_matrices3$productivity <- analysis_matrices_list[3]
	analysis_matrices3$sample <- row.names(analysis_matrices3)
	
	colnames(analysis_matrices1) <- gsub("TCR_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("BCR_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("TCRA_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("TCRB_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("TCRG_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("TCRD_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("UNPRODUCTIVE_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("PRODUCTIVE_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("ALL_", "", colnames(analysis_matrices1))
	
	colnames(analysis_matrices2) <- gsub("TCR_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("BCR_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("TCRA_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("TCRB_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("TCRG_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("TCRD_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("UNPRODUCTIVE_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("PRODUCTIVE_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("ALL_", "", colnames(analysis_matrices2))
	
	colnames(analysis_matrices3) <- gsub("TCR_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("BCR_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("TCRA_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("TCRB_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("TCRG_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("TCRD_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("UNPRODUCTIVE_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("PRODUCTIVE_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("ALL_", "", colnames(analysis_matrices3))
	
	big <- dplyr::bind_rows(analysis_matrices1, analysis_matrices2, analysis_matrices3)
	depth_file <- read.delim(depth_file, sep="\t")
	depth_file$SampleIDforDepths <- gsub("_unproductive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("_productive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("BCR_", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("TCR_", "", depth_file$SampleIDforDepths)
	
	
	duprod$SampleIDforDepths <- gsub("_unproductive", "", duprod$SampleIDforDepths)
	duprod$SampleIDforDepths <- gsub("_productive", "", duprod$SampleIDforDepths)
	duprod$SampleIDforDepths <- gsub("BCR_", "", duprod$SampleIDforDepths)
	duprod$SampleIDforDepths <- gsub("TCR_", "", duprod$SampleIDforDepths)
	
	dall$SampleIDforDepths <- gsub("_unproductive", "", dall$SampleIDforDepths)
	dall$SampleIDforDepths <- gsub("_productive", "", dall$SampleIDforDepths)
	dall$SampleIDforDepths <- gsub("BCR_", "", dall$SampleIDforDepths)
	dall$SampleIDforDepths <- gsub("TCR_", "", dall$SampleIDforDepths)
	
	dprod$SampleIDforDepths <- gsub("_unproductive", "", dprod$SampleIDforDepths)
	dprod$SampleIDforDepths <- gsub("_productive", "", dprod$SampleIDforDepths)
	dprod$SampleIDforDepths <- gsub("BCR_", "", dprod$SampleIDforDepths)
	dprod$SampleIDforDepths <- gsub("TCR_", "", dprod$SampleIDforDepths)
	
	for(i in 1:length(big$productivity)){ 
		if(grepl("_UNPRODUCTIVE", big$productivity[i])){
			big$productivity[i] <- "UNPRODUCTIVE"
		} 
		if(grepl("_PRODUCTIVE", big$productivity[i])){
			big$productivity[i] <- "PRODUCTIVE"
		} 
		if(grepl("ALL", big$productivity[i])){
			big$productivity[i] <- "ALL"
		}
	} 

	cores <- 10
	cl <- parallel::makeCluster(cores)
	doParallel::registerDoParallel(cl)
	foreach(i = 1:(length(colnames(big))-2), .export='comparator_plot',   .packages=c("ggplot2", "data.table", "ggrastr", "ggpubr", "cowplot", "gridExtra")) %dopar% {
				comparator_plot(big, outputdir, i, dall, duprod, dprod)
	}
	stopCluster(cl)
	print("Done Comparator Plots")
}
			
			
			