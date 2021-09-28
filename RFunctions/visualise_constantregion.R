
# Function to visualise the Filtering Reports for all files processed through the RBR Bulk BCR/TCR pipline
# Author: Lauren Overend 
library(stringr)

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
	
	widthx <- 0.1619718 * length(Constant_Results$totalreads)
	if (widthx > 120){
		widthx <- 120
	}
	
	
	pdf(paste0(path_to_outputdir,'/Plots/Constant_Region_Counts_QC_', Run_name, '.pdf'), width=widthx, height=14)
	p1 <- ggplot(Constant_Results_subset, aes(x=Sample, y=percentage, fill=gene)) + geom_bar(position="stack", stat="identity", color="black") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="Isotype")
	p2 <- ggplot(Constant_Results_mutation, aes(x=Sample, y=percent_IGM, fill=gene)) + geom_bar(position="stack", stat="identity", color="black") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% IgM_D Reads Reads") +labs(fill="Mutation Status") + scale_fill_discrete(labels = c("Mutated", "Unmutated"))
	plot(plot_grid(p1, p2, ncol=1))
	dev.off()
	
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Constant_Results, paste0(path_to_outputdir, "/Summary/Constant_Results_", Run_name, ".txt"), sep="\t")		
}


visualise_constant_region_bcr_layouts <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 5, path_to_layout=path_to_layout){
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
	
	
	Constant_Results$SampleBase <- Constant_Results$Sample
	Constant_Results$SampleBase <- gsub("TCRA_", "", Constant_Results$SampleBase)
	Constant_Results$SampleBase <- gsub("TCRB_", "", Constant_Results$SampleBase)
	Constant_Results$SampleBase <- gsub("TCRG_", "", Constant_Results$SampleBase)
	Constant_Results$SampleBase <- gsub("TCRD_", "", Constant_Results$SampleBase)
	Constant_Results$SampleBase <- gsub("TCR_", "", Constant_Results$SampleBase)
	Constant_Results$SampleBase <- gsub("TR_", "", Constant_Results$SampleBase)
	Constant_Results$SampleBase <- gsub("BCR_", "", Constant_Results$SampleBase)
	
	#Read in layouts file for batch information
	layouts <- read.delim(path_to_layout, sep="\t", header=TRUE)

	layouts$SampleID <- gsub("TCRA_", "", layouts$SampleID)
	layouts$SampleID <- gsub("TCRB_", "", layouts$SampleID)
	layouts$SampleID <- gsub("TCRG_", "", layouts$SampleID)
	layouts$SampleID <- gsub("TCRD_", "", layouts$SampleID)
	layouts$SampleID <- gsub("TCR_", "", layouts$SampleID)
	layouts$SampleID <- gsub("TR_", "", layouts$SampleID)
	layouts$SampleID <- gsub("BCR_", "", layouts$SampleID)
	
	layouts$Barcode <- gsub("TCRA_", "", layouts$Barcode)
	layouts$Barcode <- gsub("TCRB_", "", layouts$Barcode)
	layouts$Barcode <- gsub("TCRG_", "", layouts$Barcode)
	layouts$Barcode <- gsub("TCRD_", "", layouts$Barcode)
	layouts$Barcode <- gsub("TCR_", "", layouts$Barcode)
	layouts$Barcode <- gsub("TR_", "", layouts$Barcode)
	layouts$Barcode <- gsub("BCR_", "", layouts$Barcode)
	
	
	#Merge batch and data 
	Constant_Results <- merge(Constant_Results, layouts, by.x="SampleBase", by.y="SampleID")

	Constant_Results$Position <- as.character(Constant_Results$Position)
	Constant_Results$PCRBarcode <- as.character(Constant_Results$PCRBarcode)
	Constant_Results$Library <- as.character(Constant_Results$Library)
	Constant_Results$Plate <- as.character(Constant_Results$Plate)
	Constant_Results$Lane <- as.character(Constant_Results$Lane)
	
	if(any(!layouts$Barcode==layouts$SampleID)){	
		days <- data.frame(str_split_fixed(Constant_Results$SampleBase, "_", 2))
		days <- days$X2
		Constant_Results <- cbind(Constant_Results, days)
		Constant_Results$days <- gsub("T", "", Constant_Results$days)
		Constant_Results$days[Constant_Results$days != 1 & Constant_Results$days != 3 & Constant_Results$days != 5] <- "Control"
	} else {
		Constant_Results$days <- "NA"
	} 
	
	Constant_Results_subset <- Constant_Results[Constant_Results$gene != "IGHD/M_mutated" & Constant_Results$gene != "ALL" & Constant_Results$gene != "class_switched" & Constant_Results$gene != "IGHD/M_unmutated",]
	Constant_Results_mutation <- Constant_Results[Constant_Results$gene == "IGHD/M_mutated"  | Constant_Results$gene == "IGHD/M_unmutated",]
	
	widthx <- 0.1619718 * length(Constant_Results$Lane)
	if (widthx > 120){
		widthx <- 120
	}
	
	pdf(paste0(path_to_outputdir,'/Plots/Constant_Region_Counts_QC_ISOTYPEUSAGE', Run_name, '.pdf'), width=widthx, height=14)
	p1 <- ggplot(Constant_Results_subset, aes(x=Sample, y=percentage, fill=gene)) + geom_bar(position="stack", stat="identity", color="black") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="Isotype")
	p2 <- ggplot(Constant_Results_mutation, aes(x=Sample, y=percent_IGM, fill=gene)) + geom_bar(position="stack", stat="identity", color="black") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% IgHM Reads") +labs(fill="Mutation Status") + scale_fill_discrete(labels = c("Mutated", "Unmutated"))
	plot(plot_grid(p1, p2, ncol=1))
	dev.off()
	
	pdf(paste0(path_to_outputdir,'/Plots/Constant_Region_Counts_QC_ISOTYPEUSAGE_Layouts', Run_name, '.pdf'), width=30, height=14)
	p1 <- ggplot(Constant_Results_subset, aes(x=Library, y=percentage, fill=Lane)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Reads") +facet_wrap(~gene, scales = "free_x")
	p2 <- ggplot(Constant_Results_mutation, aes(x=Library, y=percent_IGM, fill=Lane)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% of IgM_D Reads") +facet_wrap(~gene, scales = "free_x")
	plot(p1)
	plot(p2)
	p1 <- ggplot(Constant_Results_subset, aes(x=Lane, y=percentage, fill=days)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Reads") +facet_wrap(~gene, scales = "free_x")
	p2 <- ggplot(Constant_Results_mutation, aes(x=Lane, y=percent_IGM, fill=days)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% of IgM_D Reads") +facet_wrap(~gene, scales = "free_x")
	plot(p1)
	plot(p2)
	dev.off()

	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Constant_Results, paste0(path_to_outputdir, "/Summary/Constant_Results_", Run_name, ".txt"), sep="\t")		
}

##-------------------------------------------------------------
##-------------------------------------------------------------
##-------------------------------------------------------------


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
	
	widthx <- 0.1619718 * length(Constant_Results$totalreads)
	if (widthx > 120){
		widthx <- 120
	}

	Constant_Results_subset <- Constant_Results[Constant_Results$gene != "ALL",]
	pdf(paste0(path_to_outputdir,'/Plots/Constant_Region_Counts_QC_', Run_name, '.pdf'), width=widthx, height=14)
	p1 <- ggplot(Constant_Results_subset, aes(x=Sample, y=percentage, fill=gene)) + geom_bar(position="stack", stat="identity", color="black") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="Constant Region")
	plot(plot_grid(p1, ncol=1))
	dev.off()				
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Constant_Results, paste0(path_to_outputdir, "/Summary/Constant_Results_", Run_name, ".txt"), sep="\t")		
}

visualise_constant_region_tcr_layouts <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 5, path_to_layout=path_to_layout){
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
	
	Constant_Results$SampleBase <- Constant_Results$Sample
	Constant_Results$SampleBase <- gsub("TCRA_", "", Constant_Results$SampleBase)
	Constant_Results$SampleBase <- gsub("TCRB_", "", Constant_Results$SampleBase)
	Constant_Results$SampleBase <- gsub("TCRG_", "", Constant_Results$SampleBase)
	Constant_Results$SampleBase <- gsub("TCRD_", "", Constant_Results$SampleBase)
	Constant_Results$SampleBase <- gsub("TCR_", "", Constant_Results$SampleBase)
	Constant_Results$SampleBase <- gsub("TR_", "", Constant_Results$SampleBase)
	Constant_Results$SampleBase <- gsub("BCR_", "", Constant_Results$SampleBase)
	
	#Read in layouts file for batch information
	layouts <- read.delim(path_to_layout, sep="\t", header=TRUE)

	layouts$Barcode <- gsub("TCRA_", "", layouts$Barcode)
	layouts$Barcode <- gsub("TCRB_", "", layouts$Barcode)
	layouts$Barcode <- gsub("TCRG_", "", layouts$Barcode)
	layouts$Barcode <- gsub("TCRD_", "", layouts$Barcode)
	layouts$Barcode <- gsub("TCR_", "", layouts$Barcode)
	layouts$Barcode <- gsub("TR_", "", layouts$Barcode)
	layouts$Barcode <- gsub("BCR_", "", layouts$Barcode)
	

	#Merge batch and data 
	Constant_Results <- merge(Constant_Results, layouts, by.x="SampleBase", by.y="SampleID")
	
	Constant_Results$Position <- as.character(Constant_Results$Position)
	Constant_Results$PCRBarcode <- as.character(Constant_Results$PCRBarcode)
	Constant_Results$Library <- as.character(Constant_Results$Library)
	Constant_Results$Plate <- as.character(Constant_Results$Plate)
	Constant_Results$Lane <- as.character(Constant_Results$Lane)
	
	if(any(!layouts$Barcode==layouts$SampleID)){	
		days <- data.frame(str_split_fixed(Constant_Results$SampleBase, "_", 2))
		days <- days$X2
		Constant_Results <- cbind(Constant_Results, days)
		Constant_Results$days <- gsub("T", "", Constant_Results$days)
		Constant_Results$days[Constant_Results$days != 1 & Constant_Results$days != 3 & Constant_Results$days != 5] <- "Control"
	} else {
		Constant_Results$days <- "NA"
	} 
	
	
	if(dir.exists(paste0(path_to_outputdir, "/Plots"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots"))
	}
	
	widthx <- 0.1619718 * length(Constant_Results$totalreads)
	if (widthx > 120){
		widthx <- 120
	}


	Constant_Results_subset <- Constant_Results[Constant_Results$gene != "ALL",]
	pdf(paste0(path_to_outputdir,'/Plots/Constant_Region_Counts_QC_CONSTANTUSAGE', Run_name, '.pdf'), width=widthx, height=14)
	p1 <- ggplot(Constant_Results_subset, aes(x=Sample, y=percentage, fill=gene)) + geom_bar(position="stack", stat="identity", color="black") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="Constant Region")
	plot(plot_grid(p1, ncol=1))
	dev.off()
	
	pdf(paste0(path_to_outputdir,'/Plots/Constant_Region_Counts_QC_CONSTANTUSAGE_Layouts', Run_name, '.pdf'), width=30, height=14)
	p1 <- ggplot(Constant_Results_subset, aes(x=Library, y=percentage, fill=Lane)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Reads") +facet_wrap(~gene, scales = "free_x")
	plot(p1)
		
	p1 <- ggplot(Constant_Results_subset, aes(x=Lane, y=percentage, fill=days)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Reads") +facet_wrap(~gene, scales = "free_x")
	plot(p1)
	dev.off()		

	
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Constant_Results, paste0(path_to_outputdir, "/Summary/Constant_Results_", Run_name, ".txt"), sep="\t")		
}

