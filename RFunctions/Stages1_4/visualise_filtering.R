# Function to visualise the Filtering Reports for all files processed through the RBR Bulk BCR/TCR pipline
# Author: Lauren Overend 
visualise_filtering_bcr <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 5){
	library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(doParallel)
	library(gridExtra)
	library(cowplot)
	path <- path_to_outputdir
	Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES")
	files <- list.files(path, full.name=TRUE)
	files <- grep('Filtering_report', files, value=TRUE)
	cl <- cluster_nodes
	registerDoParallel(cl)
	Filtering_Results <- foreach(i = 1:length(files), .combine=rbind, .packages='tidyverse') %dopar% {
				filtered_path <- files[i]
				output <- read.delim(filtered_path, sep="\t", header=TRUE)
				return(output)
	}
	colnames(Filtering_Results) <- c("Directory", "Sample", "Species", "Gene", "PercentReadsBCRRetained", "JoinedReads_perBatch", "JoinedReads_perSample", "ReadswithUMI", "ReadswithUniqueUMIs", "Reads_post_Vmatch", "Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueBCRS")
	filtering <- Filtering_Results[, c("Sample", "ReadswithUMI", "Reads_post_Vmatch", "Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueBCRS")]
	wide_format <- gather(filtering , "FilteringStage", "NoReads", -Sample)
	wide_format$FilteringStage <- factor(wide_format$FilteringStage, levels = c("ReadswithUMI", "Reads_post_Vmatch","Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueBCRS"))
	levels(wide_format$FilteringStage) <- c("ReadswithUMI", "Reads_post_Vmatch","Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueBCRS")   
	
	Filtering_Results$PercentageFunctional <- Filtering_Results$Reads_post_ORFfilt / Filtering_Results$Reads_post_Jmatch *100
	Filtering_Results$Duplication <- Filtering_Results$N_UniqueBCRS / Filtering_Results$Reads_post_ORFfilt *100
	Filtering_Results$PercentagePassedFiltering <- Filtering_Results$N_UniqueBCRS / Filtering_Results$JoinedReads_perSample *100
	if(dir.exists(paste0(path_to_outputdir, "/Plots"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots"))
	}
	
	widthx <- 0.2464789 * length(Filtering_Results$PercentagePassedFiltering)
	if (widthx > 120){
		widthx <- 120
		heightx <- 50
		heighty <- 20
	} else {
		heightx <- 25
		heighty <- 10
	}
		
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_QC_', Run_name, '.pdf'), width=widthx, height=heightx)
	p4 <- ggplot(Filtering_Results, aes(x=Sample, y=PercentagePassedFiltering)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % BCRs Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100)+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')
	p2 <- ggplot(Filtering_Results, aes(x=Sample, y=N_UniqueBCRS)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("BCRs: No. Unique Sequences Post-Filtering") +xlab("Sample") + geom_hline(yintercept=mean(Filtering_Results$N_UniqueBCRS), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) + sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) - sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')
	p1 <- ggplot(Filtering_Results, aes(x=Sample, y=PercentageFunctional)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ylim(0,100)
	p3 <- ggplot(Filtering_Results, aes(x=Sample, y=Duplication)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Unique VDJ nt sequences from Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted') +ylim(0,100)
	plot(plot_grid(p1,p3, p2, p4, ncol=1))
	p4 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentagePassedFiltering), y=PercentagePassedFiltering)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % BCRs Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100) +xlab("Sample")  +geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')
	p2 <- ggplot(Filtering_Results, aes(x=reorder(Sample, N_UniqueBCRS), y=N_UniqueBCRS)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("BCRs: No. Unique Sequences Post-Filtering") +xlab("Sample")  +geom_hline(yintercept=mean(Filtering_Results$N_UniqueBCRS), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) + sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) - sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')
	p1 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentageFunctional), y=PercentageFunctional)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +ylim(0,100)+xlab("Sample") +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')
	p3 <- ggplot(Filtering_Results, aes(x=reorder(Sample, Duplication), y=Duplication)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Unique VDJ nt sequences from Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+ylim(0,100)+xlab("Sample")
	plot(plot_grid(p1,p3, p2, p4, ncol=1))
	dev.off()	
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_Stages_', Run_name, '.pdf'), width=widthx, height=10)
	p1 <- ggplot(wide_format, aes(x=Sample, y=NoReads, fill=FilteringStage)) + geom_bar(stat="identity", position=position_dodge())+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("Number of Reads")
	plot(p1)
	dev.off()
	
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Filtering_Results, paste0(path_to_outputdir, "/Summary/Filtering_Results_", Run_name, ".txt"), sep="\t")		
}



visualise_filtering_bcr_layouts <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, path_to_layout = path_to_layout, cluster_nodes = 5){
	library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(doParallel)
	library(gridExtra)
	library(cowplot)
	path <- path_to_outputdir
	Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES")
	files <- list.files(path, full.name=TRUE)
	files <- grep('Filtering_report', files, value=TRUE)
	cl <- cluster_nodes
	registerDoParallel(cl)
	Filtering_Results <- foreach(i = 1:length(files), .combine=rbind, .packages='tidyverse') %dopar% {
				filtered_path <- files[i]
				output <- read.delim(filtered_path, sep="\t", header=TRUE)
				return(output)
	}
	colnames(Filtering_Results) <- c("Directory", "Sample", "Species", "Gene", "PercentReadsBCRRetained", "JoinedReads_perBatch", "JoinedReads_perSample", "ReadswithUMI", "ReadswithUniqueUMIs", "Reads_post_Vmatch", "Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueBCRS")
	filtering <- Filtering_Results[, c("Sample", "ReadswithUMI", "Reads_post_Vmatch", "Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueBCRS")]
	wide_format <- gather(filtering , "FilteringStage", "NoReads", -Sample)
	wide_format$FilteringStage <- factor(wide_format$FilteringStage, levels = c("ReadswithUMI", "Reads_post_Vmatch","Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueBCRS"))
	levels(wide_format$FilteringStage) <- c("ReadswithUMI", "Reads_post_Vmatch","Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueBCRS")   
	
	Filtering_Results$PercentageFunctional <- Filtering_Results$Reads_post_ORFfilt / Filtering_Results$Reads_post_Jmatch *100
	Filtering_Results$Duplication <- Filtering_Results$N_UniqueBCRS / Filtering_Results$Reads_post_ORFfilt *100
	Filtering_Results$PercentagePassedFiltering <- Filtering_Results$N_UniqueBCRS / Filtering_Results$JoinedReads_perSample *100
	
	#Create the base name
	Filtering_Results$SampleBase <- Filtering_Results$Sample
	Filtering_Results$SampleBase <- gsub("BCR_", "", Filtering_Results$SampleBase)
	
	#Read in layouts file for batch information
	layouts <- read.delim(path_to_layout, sep="\t", header=TRUE)
	layouts$SampleID <- gsub("BCR_", "", layouts$SampleID)

	
	#Merge batch and data 
	Filtering_Results <- merge(Filtering_Results, layouts, by.x="SampleBase", by.y="SampleID")

	
	if(dir.exists(paste0(path_to_outputdir, "/Plots"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots"))
	}
	
	Filtering_Results$Position <- as.character(Filtering_Results$Position)
	Filtering_Results$PCRBarcode <- as.character(Filtering_Results$PCRBarcode)
	Filtering_Results$Library <- as.character(Filtering_Results$Library)
	Filtering_Results$Plate <- as.character(Filtering_Results$Plate)
	Filtering_Results$Lane <- as.character(Filtering_Results$Lane)
	
	widthx <- 0.2464789 * length(Filtering_Results$Lane)
	if (widthx > 120){
		widthx <- 120
		heightx <- 50
		heighty <- 20
	} else {
		heightx <- 25
		heighty <- 10
	}
		
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_QC_', Run_name, '.pdf'), width=widthx, height=heightx)
	p4 <- ggplot(Filtering_Results, aes(x=Sample, y=PercentagePassedFiltering, shape=Lane, colour=Library)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % BCRs Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100)+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted') + guides(colour=guide_legend(ncol=4), shape=guide_legend(ncol=4))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	p2 <- ggplot(Filtering_Results, aes(x=Sample, y=N_UniqueBCRS, shape=Lane, colour=Library)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("BCRs: No. Unique Sequences Post-Filtering") +xlab("Sample") + geom_hline(yintercept=mean(Filtering_Results$N_UniqueBCRS), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) + sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) - sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted') + guides(colour=guide_legend(ncol=4), shape=guide_legend(ncol=4))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	p1 <- ggplot(Filtering_Results, aes(x=Sample, y=PercentageFunctional, shape=Lane, colour=Library)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ylim(0,100) + guides(colour=guide_legend(ncol=4), shape=guide_legend(ncol=4))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	p3 <- ggplot(Filtering_Results, aes(x=Sample, y=Duplication, shape=Lane, colour=Library)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Unique VDJ nt sequences from Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted') +ylim(0,100)+ guides(colour=guide_legend(ncol=4), shape=guide_legend(ncol=4))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	plot(plot_grid(p1,p3, p2, p4, ncol=1))
	p4 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentagePassedFiltering), shape=Lane, colour=Library, y=PercentagePassedFiltering)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % BCRs Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100) +xlab("Sample")  +geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted') + guides(colour=guide_legend(ncol=4), shape=guide_legend(ncol=4)) + scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane))))) 
	p2 <- ggplot(Filtering_Results, aes(x=reorder(Sample, N_UniqueBCRS), y=N_UniqueBCRS, shape=Lane, colour=Library)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("BCRs: No. Unique Sequences Post-Filtering") +xlab("Sample")  +geom_hline(yintercept=mean(Filtering_Results$N_UniqueBCRS), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) + sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) - sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted') + guides(colour=guide_legend(ncol=4), shape=guide_legend(ncol=4))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	p1 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentageFunctional), shape=Lane, colour=Library, y=PercentageFunctional)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +ylim(0,100)+xlab("Sample") +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ guides(colour=guide_legend(ncol=4), shape=guide_legend(ncol=4))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	p3 <- ggplot(Filtering_Results, aes(x=reorder(Sample, Duplication), y=Duplication, shape=Lane, colour=Library)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Unique VDJ nt sequences from Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+ylim(0,100)+xlab("Sample")+ guides(colour=guide_legend(ncol=4), shape=guide_legend(ncol=4))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	plot(plot_grid(p1,p3, p2, p4, ncol=1))
	dev.off()
	
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_QC_BATCHCOMPARISON_', Run_name, '.pdf'), width=15, height=heighty)
	# Plot the per batch statistics
	p4 <- ggplot(Filtering_Results, aes(x=Library, y=PercentagePassedFiltering, fill=Library)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % BCRs Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100)+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted') + guides(fill=guide_legend(ncol=4)) +facet_wrap(~Lane, scales = "free_x")
	p2 <- ggplot(Filtering_Results, aes(x=Library, y=N_UniqueBCRS, fill=Library)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("BCRs: No. Unique Sequences Post-Filtering")  + geom_hline(yintercept=mean(Filtering_Results$N_UniqueBCRS), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) + sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) - sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted') + guides(fill=guide_legend(ncol=4))+facet_wrap(~Lane, scales = "free_x")
	p1 <- ggplot(Filtering_Results, aes(x=Library, y=PercentageFunctional, fill=Library)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ylim(0,100) + guides(fill=guide_legend(ncol=4), shape=guide_legend(ncol=4))+facet_wrap(~Lane, scales = "free_x")
	p3 <- ggplot(Filtering_Results, aes(x=Library, y=Duplication, fill=Library)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Unique VDJ nt sequences from Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted') +ylim(0,100)+ guides(colour=guide_legend(ncol=4), fill=guide_legend(ncol=4))+facet_wrap(~Lane, scales = "free_x")
	plot(plot_grid(p1,p3, ncol=1))
	plot(plot_grid(p2, p4, ncol=1))
	
	p4 <- ggplot(Filtering_Results, aes(x=Lane, y=PercentagePassedFiltering, fill=Lane)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100)+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted') + guides(fill=guide_legend(ncol=4)) 
	p2 <- ggplot(Filtering_Results, aes(x=Lane, y=N_UniqueBCRS, fill=Lane)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("BCRs: No. Unique Sequences Post-Filtering") + geom_hline(yintercept=mean(Filtering_Results$N_UniqueBCRS), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) + sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) - sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted') + guides(fill=guide_legend(ncol=4))
	p1 <- ggplot(Filtering_Results, aes(x=Lane, y=PercentageFunctional, fill=Lane)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ylim(0,100) + guides(fill=guide_legend(ncol=4), shape=guide_legend(ncol=4))
	p3 <- ggplot(Filtering_Results, aes(x=Lane, y=Duplication, fill=Lane)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("BCRs: % Unique VDJ nt sequences from Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted') +ylim(0,100)+ guides(fill=guide_legend(ncol=4), shape=guide_legend(ncol=4))
	plot(plot_grid(p1,p3, p2, p4, ncol=2))
	dev.off()
	
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_Stages_', Run_name, '.pdf'), width=widthx, height=10)
	p1 <- ggplot(wide_format, aes(x=Sample, y=NoReads, fill=FilteringStage)) + geom_bar(stat="identity", position=position_dodge())+theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("Number of Reads") + guides(fill=guide_legend("Filtering Stage")) + scale_fill_discrete(labels = c("Reads with UMI", "Post V gene match", "Post J gene match", "Post ORF Filter", "Unique VDJ nt sequence"))
	plot(p1)
	dev.off()
	
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Filtering_Results, paste0(path_to_outputdir, "/Summary/Filtering_Results_", Run_name, ".txt"), sep="\t")		
}


visualise_filtering_tcr <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 5){
	library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(doParallel)
	library(gridExtra)
	library(cowplot)
	path <- path_to_outputdir
	Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES")
	files <- list.files(path, full.name=TRUE)
	files <- grep('Filtering_report', files, value=TRUE)
	cl <- cluster_nodes
	registerDoParallel(cl)
	Filtering_Results <- foreach(i = 1:length(files), .combine=rbind, .packages='tidyverse') %dopar% {
				filtered_path <- files[i]
				output <- read.delim(filtered_path, sep="\t", header=TRUE)
				return(output)
	}
	colnames(Filtering_Results) <- c("Directory", "Sample", "Species", "Gene", "PercentReadsTCRRetained", "JoinedReads_perBatch", "JoinedReads_perSample", "ReadswithUMI", "ReadswithUniqueUMIs", "Reads_post_Vmatch", "Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueTCRS")
	filtering <- Filtering_Results[, c("Sample", "ReadswithUMI", "Reads_post_Vmatch", "Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueTCRS")]
	wide_format <- gather(filtering , "FilteringStage", "NoReads", -Sample)
	wide_format$FilteringStage <- factor(wide_format$FilteringStage, levels = c("ReadswithUMI", "Reads_post_Vmatch","Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueTCRS"))
	levels(wide_format$FilteringStage) <- c("ReadswithUMI", "Reads_post_Vmatch","Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueTCRS")   
	
	Filtering_Results$PercentageFunctional <- Filtering_Results$Reads_post_ORFfilt / Filtering_Results$Reads_post_Jmatch *100
	Filtering_Results$Duplication <- Filtering_Results$N_UniqueTCRS / Filtering_Results$Reads_post_ORFfilt *100
	Filtering_Results$PercentagePassedFiltering <- Filtering_Results$N_UniqueTCRS / Filtering_Results$JoinedReads_perSample *100
	if(dir.exists(paste0(path_to_outputdir, "/Plots"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots"))
	}
	
	widthx <- 0.2464789 * length(Filtering_Results$PercentagePassedFiltering)
	if (widthx > 120){
		widthx <- 120
		heightx <- 50
		heighty <- 20
	} else {
		heightx <- 25
		heighty <- 10
	}
	
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_QC_', Run_name, '.pdf'), width=widthx, height=heightx)
	p4 <- ggplot(Filtering_Results, aes(x=Sample, y=PercentagePassedFiltering)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % TCRs Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100)+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')
	p2 <- ggplot(Filtering_Results, aes(x=Sample, y=N_UniqueTCRS)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("TCRs: No. Unique Sequences Post-Filtering") +xlab("Sample") + geom_hline(yintercept=mean(Filtering_Results$N_UniqueTCRS), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) + sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) - sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted')
	p1 <- ggplot(Filtering_Results, aes(x=Sample, y=PercentageFunctional)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ylim(0,100)
	p3 <- ggplot(Filtering_Results, aes(x=Sample, y=Duplication)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % Unique VDJ nt sequences from Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted') +ylim(0,100)
	plot(plot_grid(p1,p3, p2, p4, ncol=1))
	p4 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentagePassedFiltering), y=PercentagePassedFiltering)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % TCRs Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100) +xlab("Sample")  +geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')
	p2 <- ggplot(Filtering_Results, aes(x=reorder(Sample, N_UniqueTCRS), y=N_UniqueTCRS)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("TCRs: No. Unique Sequences Post-Filtering") +xlab("Sample")  +geom_hline(yintercept=mean(Filtering_Results$N_UniqueTCRS), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) + sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) - sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted')
	p1 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentageFunctional), y=PercentageFunctional)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +ylim(0,100)+xlab("Sample") +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')
	p3 <- ggplot(Filtering_Results, aes(x=reorder(Sample, Duplication), y=Duplication)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % Unique VDJ nt sequences from Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+ylim(0,100)+xlab("Sample")
	plot(plot_grid(p1,p3, p2, p4, ncol=1))
	dev.off()	
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_Stages_', Run_name, '.pdf'), width=widthx, height=10)
	p1 <- ggplot(wide_format, aes(x=Sample, y=NoReads, fill=FilteringStage)) + geom_bar(stat="identity", position=position_dodge())+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("Number of Reads")
	plot(p1)
	dev.off()
	
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Filtering_Results, paste0(path_to_outputdir, "/Summary/Filtering_Results_", Run_name, ".txt"), sep="\t")		
}



visualise_filtering_tcr_layouts <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 5, path_to_layout = path_to_layout){
	library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(doParallel)
	library(gridExtra)
	library(cowplot)
	path <- path_to_outputdir
	Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES")
	files <- list.files(path, full.name=TRUE)
	files <- grep('Filtering_report', files, value=TRUE)
	cl <- cluster_nodes
	registerDoParallel(cl)
	Filtering_Results <- foreach(i = 1:length(files), .combine=rbind, .packages='tidyverse') %dopar% {
				filtered_path <- files[i]
				output <- read.delim(filtered_path, sep="\t", header=TRUE)
				return(output)
	}
	colnames(Filtering_Results) <- c("Directory", "Sample", "Species", "Gene", "PercentReadsTCRRetained", "JoinedReads_perBatch", "JoinedReads_perSample", "ReadswithUMI", "ReadswithUniqueUMIs", "Reads_post_Vmatch", "Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueTCRS")
	filtering <- Filtering_Results[, c("Sample", "ReadswithUMI", "Reads_post_Vmatch", "Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueTCRS")]
	wide_format <- gather(filtering , "FilteringStage", "NoReads", -Sample)
	wide_format$FilteringStage <- factor(wide_format$FilteringStage, levels = c("ReadswithUMI", "Reads_post_Vmatch","Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueTCRS"))
	levels(wide_format$FilteringStage) <- c("ReadswithUMI", "Reads_post_Vmatch","Reads_post_Jmatch", "Reads_post_ORFfilt", "N_UniqueTCRS")   
	
	Filtering_Results$PercentageFunctional <- Filtering_Results$Reads_post_ORFfilt / Filtering_Results$Reads_post_Jmatch *100
	Filtering_Results$Duplication <- Filtering_Results$N_UniqueTCRS / Filtering_Results$Reads_post_ORFfilt *100
	Filtering_Results$PercentagePassedFiltering <- Filtering_Results$N_UniqueTCRS / Filtering_Results$JoinedReads_perSample *100
	if(dir.exists(paste0(path_to_outputdir, "/Plots"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots"))
	}
	
	Filtering_Results$PercentageFunctional <- Filtering_Results$Reads_post_ORFfilt / Filtering_Results$Reads_post_Jmatch *100
	Filtering_Results$Duplication <- Filtering_Results$N_UniqueTCRS / Filtering_Results$Reads_post_ORFfilt *100
	Filtering_Results$PercentagePassedFiltering <- Filtering_Results$N_UniqueTCRS / Filtering_Results$JoinedReads_perSample *100
	
	#Create the base name
	Filtering_Results$SampleBase <- Filtering_Results$Sample
	Filtering_Results$SampleBase <- gsub("TCRA_", "", Filtering_Results$SampleBase)
	Filtering_Results$SampleBase <- gsub("TCRB_", "", Filtering_Results$SampleBase)
	Filtering_Results$SampleBase <- gsub("TCRG_", "", Filtering_Results$SampleBase)
	Filtering_Results$SampleBase <- gsub("TCRD_", "", Filtering_Results$SampleBase)
	Filtering_Results$SampleBase <- gsub("TCR_", "", Filtering_Results$SampleBase)
	Filtering_Results$SampleBase <- gsub("TR_", "", Filtering_Results$SampleBase)

	#Read in layouts file for batch information
	layouts <- read.delim(path_to_layout, sep="\t", header=TRUE)

	
	layouts$SampleID <- gsub("TCRA_", "", layouts$SampleID)
	layouts$SampleID <- gsub("TCRB_", "", layouts$SampleID)
	layouts$SampleID <- gsub("TCRG_", "", layouts$SampleID)
	layouts$SampleID <- gsub("TCRD_", "", layouts$SampleID)
	layouts$SampleID <- gsub("TCR_", "", layouts$SampleID)
	layouts$SampleID <- gsub("TR_", "", layouts$SampleID)
	
	
	#Merge batch and data 
	Filtering_Results <- merge(Filtering_Results, layouts, by.x="SampleBase", by.y="SampleID")

	Filtering_Results$Position <- as.character(Filtering_Results$Position)
	Filtering_Results$PCRBarcode <- as.character(Filtering_Results$PCRBarcode)
	Filtering_Results$Library <- as.character(Filtering_Results$Library)
	Filtering_Results$Plate <- as.character(Filtering_Results$Plate)
	Filtering_Results$Lane <- as.character(Filtering_Results$Lane)
	
	widthx <- 0.2464789 * length(Filtering_Results$Lane)
	if (widthx > 120){
		widthx <- 120
		heightx <- 50
		heighty <- 30
	} else {
		heightx <- 25
		heighty <- 10
	}
	
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_QC_', Run_name, '.pdf'), width=widthx, height=heightx)
	p4 <- ggplot(Filtering_Results, aes(x=Sample, y=PercentagePassedFiltering, shape=Lane, colour=Library)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % TCRs Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100)+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted') + guides(colour=guide_legend(ncol=5), shape=guide_legend(ncol=5))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	p2 <- ggplot(Filtering_Results, aes(x=Sample, y=N_UniqueTCRS, shape=Lane, colour=Library)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("TCRs: No. Unique Sequences Post-Filtering") +xlab("Sample") + geom_hline(yintercept=mean(Filtering_Results$N_UniqueTCRS), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) + sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) - sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted') + guides(colour=guide_legend(ncol=5), shape=guide_legend(ncol=5))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	p1 <- ggplot(Filtering_Results, aes(x=Sample, y=PercentageFunctional, shape=Lane, colour=Library)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ylim(0,100) + guides(colour=guide_legend(ncol=5), shape=guide_legend(ncol=5))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	p3 <- ggplot(Filtering_Results, aes(x=Sample, y=Duplication, shape=Lane, colour=Library)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % Unique VDJ nt sequences from Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted') +ylim(0,100)+ guides(colour=guide_legend(ncol=5), shape=guide_legend(ncol=5))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	plot(plot_grid(p1,p3, p2, p4, ncol=1))
	p4 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentagePassedFiltering), shape=Lane, colour=Library, y=PercentagePassedFiltering)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % TCRs Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100) +xlab("Sample")  +geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted') + guides(colour=guide_legend(ncol=5), shape=guide_legend(ncol=5))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	p2 <- ggplot(Filtering_Results, aes(x=reorder(Sample, N_UniqueTCRS), y=N_UniqueTCRS, shape=Lane, colour=Library)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("TCRs: No. Unique Sequences Post-Filtering") +xlab("Sample")  +geom_hline(yintercept=mean(Filtering_Results$N_UniqueTCRS), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) + sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) - sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted') + guides(colour=guide_legend(ncol=5), shape=guide_legend(ncol=5))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	p1 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentageFunctional), shape=Lane, colour=Library, y=PercentageFunctional)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +ylim(0,100)+xlab("Sample") +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ guides(colour=guide_legend(ncol=5), shape=guide_legend(ncol=5))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	p3 <- ggplot(Filtering_Results, aes(x=reorder(Sample, Duplication), y=Duplication, shape=Lane, colour=Library)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % Unique VDJ nt sequences from Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+ylim(0,100)+xlab("Sample")+ guides(colour=guide_legend(ncol=5), shape=guide_legend(ncol=5))+ scale_shape_manual(values = c(c(0:length(unique(Filtering_Results$Lane)))))
	plot(plot_grid(p1,p3, p2, p4, ncol=1))
	dev.off()
	
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_QC_BATCHCOMPARISON_', Run_name, '.pdf'), width=30, height=heighty)
	# Plot the per batch statistics
	p4 <- ggplot(Filtering_Results, aes(x=Library, y=PercentagePassedFiltering, fill=Library)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % TCRs Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100)+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted') + guides(fill=guide_legend(ncol=5)) +facet_wrap(~Lane, scales = "free_x")
	p2 <- ggplot(Filtering_Results, aes(x=Library, y=N_UniqueTCRS, fill=Library)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("TCRs: No. Unique Sequences Post-Filtering")  + geom_hline(yintercept=mean(Filtering_Results$N_UniqueTCRS), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) + sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) - sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted') + guides(fill=guide_legend(ncol=5))+facet_wrap(~Lane, scales = "free_x")
	p1 <- ggplot(Filtering_Results, aes(x=Library, y=PercentageFunctional, fill=Library)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ylim(0,100) + guides(fill=guide_legend(ncol=5), shape=guide_legend(ncol=5))+facet_wrap(~Lane, scales = "free_x")
	p3 <- ggplot(Filtering_Results, aes(x=Library, y=Duplication, fill=Library)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % Unique VDJ nt sequences from Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted') +ylim(0,100)+ guides(colour=guide_legend(ncol=5), fill=guide_legend(ncol=5))+facet_wrap(~Lane, scales = "free_x")
	plot(plot_grid(p1,p3, ncol=1))
	plot(plot_grid(p2, p4, ncol=1))
	
	p4 <- ggplot(Filtering_Results, aes(x=Lane, y=PercentagePassedFiltering, fill=Lane)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % TCRs Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100)+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted') + guides(fill=guide_legend(ncol=5)) 
	p2 <- ggplot(Filtering_Results, aes(x=Lane, y=N_UniqueTCRS, fill=Lane)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("TCRs: No. Unique Sequences Post-Filtering") + geom_hline(yintercept=mean(Filtering_Results$N_UniqueTCRS), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) + sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) - sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted') + guides(fill=guide_legend(ncol=5))
	p1 <- ggplot(Filtering_Results, aes(x=Lane, y=PercentageFunctional, fill=Lane)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ylim(0,100) + guides(fill=guide_legend(ncol=5), shape=guide_legend(ncol=5))
	p3 <- ggplot(Filtering_Results, aes(x=Lane, y=Duplication, fill=Lane)) + geom_boxplot()  +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("TCRs: % Unique VDJ nt sequences from Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted') +ylim(0,100)+ guides(fill=guide_legend(ncol=5), shape=guide_legend(ncol=5))
	plot(plot_grid(p1,p3, p2, p4, ncol=2))
	dev.off()	
	
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_Stages_', Run_name, '.pdf'), width=widthx, height=10)
	p1 <- ggplot(wide_format, aes(x=Sample, y=NoReads, fill=FilteringStage)) + geom_bar(stat="identity", position=position_dodge())+theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("Number of Reads") + guides(fill=guide_legend("Filtering Stage")) + scale_fill_discrete(labels = c("Reads with UMI", "Post V gene match", "Post J gene match", "Post ORF Filter", "Unique VDJ nt sequence"))
	plot(p1)
	dev.off()
	
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Filtering_Results, paste0(path_to_outputdir, "/Summary/Filtering_Results_", Run_name, ".txt"), sep="\t")		
}


