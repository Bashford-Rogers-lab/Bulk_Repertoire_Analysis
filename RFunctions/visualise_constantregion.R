
# Function to visualise the Filtering Reports for all files processed through the RBR Bulk BCR/TCR pipline
# Author: Lauren Overend 

visualise_constant_region_bcr <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 5){
	library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(doParallel)
	library(gridExtra)
	path <- path_to_outputdir
	Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES/ANNOTATIONS")
	files <- list.files(path, full.name=TRUE)
	files <- grep('Constant_region_counts', files, value=TRUE)
	cl <- cluster_nodes
	registerDoParallel(cl)
	Constant_Results <- foreach(i = 1:length(files), .combine=rbind, .packages='tidyverse') %dopar% {
				filtered_path <- files[i]
				output <- read.delim(filtered_path, sep="\t", header=TRUE)
				output$percentage <- output$frequency / output$frequency[output$gene=="ALL"]  * 100
				output$totalreads <- output$frequency[output$gene=="ALL"]
				output$percentage[output$gene == "IGHD/M_mutated" | output$gene == "ALL" | output$gene == "class_switched" | output$gene == "IGHD/M_unmutated"] <- NA
				IgM_D <- sum(output$frequency[output$gene == "IGHD/M_mutated" | output$gene == "IGHD/M_unmutated"])
				output$percent_IGM <- output$frequency / IgM_D  * 100
				output$percent_IGM[output$gene != "IGHD/M_mutated" & output$gene != "IGHD/M_unmutated"] <- NA
				colnames(output)[1] <- "Sample"
				return(output)
	}
	if(dir.exists(paste0(path_to_outputdir, "/Plots"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots"))
	}
	Constant_Results_subset <- Constant_Results[Constant_Results$gene != "IGHD/M_mutated" & Constant_Results$gene != "ALL" & Constant_Results$gene != "class_switched" & Constant_Results$gene != "IGHD/M_unmutated",]
	Constant_Results_mutation <- Constant_Results[Constant_Results$gene == "IGHD/M_mutated"  | Constant_Results$gene == "IGHD/M_unmutated",]
	pdf(paste0(path_to_outputdir,'/Plots/Constant_Region_Counts_QC_', Run_name, '.pdf'), width=23, height=14)
	p1 <- ggplot(Constant_Results_subset, aes(x=Sample, y=percentage, fill=gene)) + geom_bar(position="stack", stat="identity") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="Isotype")
	p2 <- ggplot(Constant_Results_mutation, aes(x=Sample, y=percent_IGM, fill=gene)) + geom_bar(position="stack", stat="identity") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% IgM Reads") +labs(fill="Mutation Status")
	plot(plot_grid(p1, p2, ncol=1))
	dev.off()				
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Constant_Results, paste0(path_to_outputdir, "/Summary/Constant_Results_", Run_name, ".txt"), sep="\t")		
}

visualise_constant_region_tcr <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 5){
	library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(doParallel)
	library(gridExtra)
	path <- path_to_outputdir
	Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES/ANNOTATIONS")
	files <- list.files(path, full.name=TRUE)
	files <- grep('Constant_region_counts', files, value=TRUE)
	cl <- cluster_nodes
	registerDoParallel(cl)
	Constant_Results <- foreach(i = 1:length(files), .combine=rbind, .packages='tidyverse') %dopar% {
				filtered_path <- files[i]
				output <- read.delim(filtered_path, sep="\t", header=TRUE)
				output <- output[output$gene != "IGHD/M_mutated" & output$gene != "class_switched" & output$gene != "IGHD/M_unmutated",]
				output$percentage <- output$frequency / output$frequency[output$gene=="ALL"]  * 100
				output$totalreads <- output$frequency[output$gene=="ALL"]
				colnames(output)[1] <- "Sample"
				return(output)
	}
	if(dir.exists(paste0(path_to_outputdir, "/Plots"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots"))
	}
	Constant_Results_subset <- Constant_Results[Constant_Results$gene != "ALL",]
	pdf(paste0(path_to_outputdir,'/Plots/Constant_Region_Counts_QC_', Run_name, '.pdf'), width=23, height=14)
	p1 <- ggplot(Constant_Results_subset, aes(x=Sample, y=percentage, fill=gene)) + geom_bar(position="stack", stat="identity") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="Constant Region")
	plot(plot_grid(p1, ncol=1))
	dev.off()				
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Constant_Results, paste0(path_to_outputdir, "/Summary/Constant_Results_", Run_name, ".txt"), sep="\t")		
}
