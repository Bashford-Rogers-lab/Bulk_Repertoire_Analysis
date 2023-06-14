# Function to visualise the Filtering Reports for all files processed through the RBR Bulk BCR/TCR pipline
# Author: Lauren Overend 
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(gridExtra))
	
visualise_isoptype_cluster_bcr <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 5){

	path <- path_to_outputdir
	Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES/ANNOTATIONS")
	files <- list.files(path, full.name=TRUE)
	files <- grep('IsoTyper_chain_repertoire_statistics_file_', files, value=TRUE)
	cl <- cluster_nodes
	registerDoParallel(cl)
	Iso_Clustering <- foreach(i = 1:length(files), .combine=rbind, .packages='tidyverse') %dopar% {
				filtered_path <- files[i]
				output <- read.delim(filtered_path, sep="\t", header=TRUE)
			    colnames(output) <- c("Sample", "Isotype", "N_Reads", "N_Vertices", "Vertex_Gini_Index", "Cluster_Gini_Index", "Largest_cluster_percent", "Second_cluster_percent")
				return(output)
	}
	if(dir.exists(paste0(path_to_outputdir, "/Plots"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots"))
	}
	widthx <- 0.2112676 * length(unique(Iso_Clustering$Sample))
	if (widthx > 120){
		widthx <- 120
	}
	if(widthx < 5){
		widthx <- 5
	}
	
	pdf(paste0(path_to_outputdir,'/Plots/CLUSTERING_QC.pdf'),width=widthx, height=7.5)
	p1 <- ggplot(Iso_Clustering, aes(x=Sample, y=Largest_cluster_percent, fill=Isotype)) +geom_bar(stat="identity", position=position_dodge()) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Largest Cluster (%)") +labs(fill="Isotype")
	p2 <- ggplot(Iso_Clustering, aes(x=Sample, y=Second_cluster_percent, fill=Isotype)) +geom_bar(stat="identity" , position=position_dodge()) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Second Large Cluster (%)") +labs(fill="Isotype")
	p3 <- ggplot(Iso_Clustering, aes(x=Sample, y= Cluster_Gini_Index, fill=Isotype)) +geom_bar(stat="identity" , position=position_dodge()) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Cluster Gini Index") +labs(fill="Isotype")
	p4 <- ggplot(Iso_Clustering, aes(x=Sample, y= Vertex_Gini_Index, fill=Isotype)) +geom_bar(stat="identity" , position=position_dodge()) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Vertex Gini Index") +labs(fill="Isotype")
	plot(plot_grid(p1, p2, ncol=1))
	plot(plot_grid(p3, p4, ncol=1))
	dev.off()				
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Iso_Clustering, paste0(path_to_outputdir, "/Summary/CLUSTERING_QC.txt"), sep="\t")
}


visualise_isoptype_cluster_tcr <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 5){
	library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(doParallel)
	library(gridExtra)
	path <- path_to_outputdir
	Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES/ANNOTATIONS")
	files <- list.files(path, full.name=TRUE)
	files <- grep('IsoTyper_chain_repertoire_statistics_file_', files, value=TRUE)
	cl <- cluster_nodes
	registerDoParallel(cl)
	Iso_Clustering <- foreach(i = 1:length(files), .combine=rbind, .packages='tidyverse') %dopar% {
				filtered_path <- files[i]
				output <- read.delim(filtered_path, sep="\t", header=TRUE)
			    colnames(output) <- c("Sample", "Isotype", "N_Reads", "N_Vertices", "Vertex_Gini_Index", "Cluster_Gini_Index", "Largest_cluster_percent", "Second_cluster_percent")
				output$read_count <- sum(output$N_Reads)
				output$constant_percentage <- output$N_Reads / output$read_count * 100
				return(output)
	}
	if(dir.exists(paste0(path_to_outputdir, "/Plots"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots"))
	}

	max_identified <- aggregate(Iso_Clustering$N_Reads, by=list(Iso_Clustering$Isotype), FUN=sum)
	max_identified <- as.character(max_identified$Group.1[which.max(max_identified$x)])
	
	if(max_identified =="TRBC2"){
		max_identified <- c(max_identified, "TRBC1")
	}
	
	if(max_identified =="TRBC1"){
		max_identified <- c(max_identified, "TRBC2")
	}
	
	if(max_identified =="TRGC2"){
		max_identified <- c(max_identified, "TRGC1")
	}
	
	if(max_identified =="TRGC1"){
		max_identified <- c(max_identified, "TRGC2")
	}
	
	Iso_Clustering_subset <- Iso_Clustering[Iso_Clustering$Isotype==max_identified[1]| Iso_Clustering$Isotype==max_identified[length(max_identified)],] 
		
	if(length(max_identified)==2){
		plot_width <- 60
	} else {
		plot_width <- 30
	} 
	
	widthx <- 0.2112676 * length(unique(Iso_Clustering$Sample))
	if (widthx > 120){
		widthx <- 120
	} 
	if(widthx < 5){
		widthx <- 5
	}
	
	pdf(paste0(path_to_outputdir,'/Plots/CLUSTERING_QC.pdf'),width=widthx, height=7.5)
	p1 <- ggplot(Iso_Clustering_subset, aes(x=Sample, y=Largest_cluster_percent, fill=Isotype)) +geom_bar(stat="identity", position=position_dodge(), colour="black") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Largest Cluster (%)") +labs(fill="Constant Region") +facet_wrap(~Isotype, scales="free_y")
	p2 <- ggplot(Iso_Clustering_subset, aes(x=Sample, y=Second_cluster_percent, fill=Isotype)) +geom_bar(stat="identity" , position=position_dodge(), colour="black") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Second Large Cluster (%)") +labs(fill="Constant Region")+facet_wrap(~Isotype, scales="free_y")
	p3 <- ggplot(Iso_Clustering_subset, aes(x=Sample, y= Cluster_Gini_Index, fill=Isotype)) +geom_bar(stat="identity" , position=position_dodge(), colour="black") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Cluster Gini Index") +labs(fill="Constant Region")+facet_wrap(~Isotype, scales="free_y")
	p4 <- ggplot(Iso_Clustering_subset, aes(x=Sample, y= Vertex_Gini_Index, fill=Isotype)) +geom_bar(stat="identity" , position=position_dodge(), colour="black") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Vertex Gini Index") +labs(fill="Constant Region")+facet_wrap(~Isotype, scales="free_y")
	plot(plot_grid(p1, p2, ncol=1))
	plot(plot_grid(p3, p4, ncol=1))
	dev.off()		
	
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Iso_Clustering, paste0(path_to_outputdir, "/Summary/CLUSTERING_QC.txt"), sep="\t")
}