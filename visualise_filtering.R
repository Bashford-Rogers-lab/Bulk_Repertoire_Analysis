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
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_QC_', Run_name, '.pdf'), width=30, height=15)
	p4 <- ggplot(Filtering_Results, aes(x=Sample, y=PercentagePassedFiltering)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% BCRs passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100)+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')
	p2 <- ggplot(Filtering_Results, aes(x=Sample, y=N_UniqueBCRS)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("Unique BCR Sequences Post-Filtering") +xlab("Sample") + geom_hline(yintercept=mean(Filtering_Results$N_UniqueBCRS), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) + sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) - sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')
	p1 <- ggplot(Filtering_Results, aes(x=Sample, y=PercentageFunctional)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% of BCRs which are Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ylim(0,100)
	p3 <- ggplot(Filtering_Results, aes(x=Sample, y=Duplication)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Unique BCRs in Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted') +ylim(0,100)
	plot(plot_grid(p1,p3, p2, p4, ncol=1))
	p4 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentagePassedFiltering), y=PercentagePassedFiltering)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("%BCRs passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100) +xlab("Sample")  +geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')
	p2 <- ggplot(Filtering_Results, aes(x=reorder(Sample, N_UniqueBCRS), y=N_UniqueBCRS)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("Unique BCR Sequences Post-Filtering") +xlab("Sample")  +geom_hline(yintercept=mean(Filtering_Results$N_UniqueBCRS), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) + sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) - sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')
	p1 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentageFunctional), y=PercentageFunctional)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% of BCRs which are Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +ylim(0,100)+xlab("Sample") +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')
	p3 <- ggplot(Filtering_Results, aes(x=reorder(Sample, Duplication), y=Duplication)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Unique BCRs in Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+ylim(0,100)+xlab("Sample")
	plot(plot_grid(p1,p3, p2, p4, ncol=1))
	dev.off()	
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_Stages_', Run_name, '.pdf'), width=30, height=10)
	p1 <- ggplot(wide_format, aes(x=Sample, y=NoReads, fill=FilteringStage)) + geom_bar(stat="identity", position=position_dodge())+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("Number of Reads")
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
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_QC_', Run_name, '.pdf'), width=30, height=15)
	p4 <- ggplot(Filtering_Results, aes(x=Sample, y=PercentagePassedFiltering)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% TCRs passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100)+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')
	p2 <- ggplot(Filtering_Results, aes(x=Sample, y=N_UniqueTCRS)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("Unique TCR Sequences Post-Filtering") +xlab("Sample") + geom_hline(yintercept=mean(Filtering_Results$N_UniqueTCRS), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) + sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) - sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted')
	p1 <- ggplot(Filtering_Results, aes(x=Sample, y=PercentageFunctional)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% of TCRs which are Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ylim(0,100)
	p3 <- ggplot(Filtering_Results, aes(x=Sample, y=Duplication)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Unique TCRs in Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted') +ylim(0,100)
	plot(plot_grid(p1,p3, p2, p4, ncol=1))
	p4 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentagePassedFiltering), y=PercentagePassedFiltering)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("%TCRs passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100) +xlab("Sample")  +geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')
	p2 <- ggplot(Filtering_Results, aes(x=reorder(Sample, N_UniqueTCRS), y=N_UniqueTCRS)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("Unique TCR Sequences Post-Filtering") +xlab("Sample")  +geom_hline(yintercept=mean(Filtering_Results$N_UniqueTCRS), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) + sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueTCRS) - sd(Filtering_Results$N_UniqueTCRS)), col="purple", linetype='dotted')
	p1 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentageFunctional), y=PercentageFunctional)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% of TCRs which are Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +ylim(0,100)+xlab("Sample") +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')
	p3 <- ggplot(Filtering_Results, aes(x=reorder(Sample, Duplication), y=Duplication)) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Unique TCRs in Filtered UMIs") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+ylim(0,100)+xlab("Sample")
	plot(plot_grid(p1,p3, p2, p4, ncol=1))
	dev.off()	
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_Stages_', Run_name, '.pdf'), width=30, height=10)
	p1 <- ggplot(wide_format, aes(x=Sample, y=NoReads, fill=FilteringStage)) + geom_bar(stat="identity", position=position_dodge())+theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("Number of Reads")
	plot(p1)
	dev.off()
	
	if(dir.exists(paste0(path_to_outputdir, "/Summary"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary"))
	}				
	write.table(Filtering_Results, paste0(path_to_outputdir, "/Summary/Filtering_Results_", Run_name, ".txt"), sep="\t")		
}


