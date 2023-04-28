
# Function to visualise the Filtering Reports for all files processed through the RBR Bulk BCR/TCR pipline
# Author: Lauren Overend 
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(gridExtra))
	
	
visualise_vj_usage_bcr <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 5){

	path <- path_to_outputdir
	Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES/ANNOTATIONS")
	files <- list.files(path, full.name=TRUE)
	files <- grep('Gene_frequencies', files, value=TRUE)
	cl <- cluster_nodes
	registerDoParallel(cl)
	VJ_Results <- foreach(i = 1:length(files), .combine=rbind, .packages='tidyverse') %dopar% {
				filtered_path <- files[i]
				output <- read.delim(filtered_path, sep="\t", header=FALSE)
				sample <- output[1,1]
				genes <- strsplit(output$V2, "|", fixed=TRUE)
				df  <- matrix(unlist(genes), ncol=2, byrow=TRUE)
				df <- data.frame(df)
				output$V5 <- df$X1
				output$V6 <- df$X2
			    colnames(output) <- c("Sample", "VJ_Segement", "Count", "V_family", "V_Gene", "J_gene")
				return(output)
	}
	frequency_variable_genes_id <- aggregate(VJ_Results$Count, by=list(V_Gene=VJ_Results$V_Gene, Sample=VJ_Results$Sample), FUN=sum)
	colnames(frequency_variable_genes_id)[3] <- "count"
	frequency_variable_genes_family <- aggregate(VJ_Results$Count, by=list(V_family=VJ_Results$V_family, Sample=VJ_Results$Sample), FUN=sum)
	colnames(frequency_variable_genes_family)[3] <- "count"
	frequency_variable_genes_J <- aggregate(VJ_Results$Count, by=list(J_gene=VJ_Results$J_gene, Sample=VJ_Results$Sample), FUN=sum)
	colnames(frequency_variable_genes_J)[3] <- "count"
	if(dir.exists(paste0(path_to_outputdir, "/Plots"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots"))
	}
	
	widthx <- 0.2112676 * length(VJ_Results$Sample)
	if (widthx > 120){
		widthx <- 120
	}
	if(widthx < 5){
		widthx <- 5
	}
	
	pdf(paste0(path_to_outputdir, '/Plots/V_GENE_USAGE_QC.pdf'), width=widthx, height=20)
	p1 <- ggplot(frequency_variable_genes_id, aes(x=Sample, y=count, fill=V_Gene)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="V gene")
	p2 <- ggplot(frequency_variable_genes_family, aes(x=Sample, y=count, fill=V_family)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="V family")
	p3 <- ggplot(frequency_variable_genes_J, aes(x=Sample, y=count, fill=J_gene)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="J gene")
	plot(plot_grid(p1, p2, ncol=1))
	plot(plot_grid(p3, ncol=1))
	dev.off()				
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(VJ_Results, paste0(path_to_outputdir, "/Summary/VJ_Results.txt"), sep="\t")
	write.table(frequency_variable_genes_id, paste0(path_to_outputdir, "/Summary/V_Gene_Usage_", Run_name, ".txt"), sep="\t")
	write.table(frequency_variable_genes_family, paste0(path_to_outputdir, "/Summary/V_family_Usage_", Run_name, ".txt"), sep="\t")
	write.table(frequency_variable_genes_J, paste0(path_to_outputdir, "/Summary/J_gene_Usage_", Run_name, ".txt"), sep="\t")
}


visualise_vj_usage_tcr <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 5){
	library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(doParallel)
	library(gridExtra)
	path <- path_to_outputdir
	Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES/ANNOTATIONS")
	files <- list.files(path, full.name=TRUE)
	files <- grep('Gene_frequencies', files, value=TRUE)
	cl <- cluster_nodes
	registerDoParallel(cl)
	VJ_Results <- foreach(i = 1:length(files), .combine=rbind, .packages='tidyverse') %dopar% {
				filtered_path <- files[i]
				output <- read.delim(filtered_path, sep="\t", header=FALSE)
				sample <- output[1,1]
				genes <- strsplit(output$V2, "|", fixed=TRUE)
				df  <- matrix(unlist(genes), ncol=2, byrow=TRUE)
				df <- data.frame(df)
				output$V5 <- df$X1
				output$V6 <- df$X2
			    colnames(output) <- c("Sample", "VJ_Segement", "Count", "V_family", "V_Gene", "J_gene")
				return(output)
	}
	frequency_variable_genes_id <- aggregate(VJ_Results$Count, by=list(V_Gene=VJ_Results$V_Gene, Sample=VJ_Results$Sample), FUN=sum)
	colnames(frequency_variable_genes_id)[3] <- "count"
	frequency_variable_genes_family <- aggregate(VJ_Results$Count, by=list(V_family=VJ_Results$V_family, Sample=VJ_Results$Sample), FUN=sum)
	colnames(frequency_variable_genes_family)[3] <- "count"
	frequency_variable_genes_J <- aggregate(VJ_Results$Count, by=list(J_gene=VJ_Results$J_gene, Sample=VJ_Results$Sample), FUN=sum)
	colnames(frequency_variable_genes_J)[3] <- "count"
	if(dir.exists(paste0(path_to_outputdir, "/Plots"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots"))
	}
	
	widthx <- 0.2112676 * length(VJ_Results$Sample)
	if (widthx > 120){
		widthx <- 120
	} 
	if(widthx < 5){
		widthx <- 5
	}
	pdf(paste0(path_to_outputdir, '/Plots/V_GENE_USAGE_QC.pdf'), width=widthx, height=20)
	p1 <- ggplot(frequency_variable_genes_id, aes(x=Sample, y=count, fill=V_Gene)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="V gene")
	p2 <- ggplot(frequency_variable_genes_family, aes(x=Sample, y=count, fill=V_family)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="V family")
	p3 <- ggplot(frequency_variable_genes_J, aes(x=Sample, y=count, fill=J_gene)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="J gene")
	plot(plot_grid(p1, p2, ncol=1))
	plot(plot_grid(p3, ncol=1))
	dev.off()				
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(VJ_Results, paste0(path_to_outputdir, "/Summary/VJ_Results.txt"), sep="\t")
	write.table(frequency_variable_genes_id, paste0(path_to_outputdir, "/Summary/V_Gene_Usage_", Run_name, ".txt"), sep="\t")
	write.table(frequency_variable_genes_family, paste0(path_to_outputdir, "/Summary/V_family_Usage_", Run_name, ".txt"), sep="\t")
	write.table(frequency_variable_genes_J, paste0(path_to_outputdir, "/Summary/J_gene_Usage_", Run_name, ".txt"), sep="\t")
}