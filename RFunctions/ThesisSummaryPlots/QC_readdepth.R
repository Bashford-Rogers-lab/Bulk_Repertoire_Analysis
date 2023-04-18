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
	library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(doParallel)
	library(gridExtra)
	library(cowplot)
	library(data.table)
#--------------------------------------------------------------------------------------------------------------------------------------
visualise_filtering_bcr_layouts_neat <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, path_to_layout = path_to_layout, plot_dir=plot_dir, stat_dir=stat_dir, cluster_nodes = 5, chain=chain){
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
	if(widthx < 5){
		widthx <- 5
	}
	Filtering_Results$Library <- as.numeric(Filtering_Results$Library )
	Filtering_Results$Library <- factor(Filtering_Results$Library, levels=c(1:max(Filtering_Results$Library)))   
	Filtering_Results$Lane <- as.numeric(Filtering_Results$Lane )
	Filtering_Results$Lane <- factor(Filtering_Results$Lane, levels=c(1:max(Filtering_Results$Lane)))
	Filtering_Results$DISEASE <-  "S"
	Filtering_Results$DISEASE[Filtering_Results$Sample  %like% "HV"] <-  "H"
	Filtering_Results$DISEASE[Filtering_Results$Sample  %like% "POSITIVE"] <-  "T"
	pdf(paste0(plot_dir,'/Filtering_QC_', run_name, '.pdf'), width=13, height=7)
	p1 <- ggplot(Filtering_Results, aes(x = Library, y=N_UniqueBCRS, fill=Lane)) + geom_boxplot(alpha=0.5) +facet_grid(~Lane,  scales="free_x", space="free_x") + theme_bw() + xlab("Sequencing Library") + ylab("Number of  Unique VDJs") +guides(fill="none")+ geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) + sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) - sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted') +geom_hline(yintercept=mean(Filtering_Results$N_UniqueBCRS), col="red", linetype='dotted')+geom_hline(yintercept=1500, col="green")+ggtitle(chain)
	p2 <- ggplot(Filtering_Results, aes(x = Library, y=PercentagePassedFiltering, fill=Lane)) + geom_boxplot(alpha=0.5)  +facet_grid(~Lane,  scales="free_x", space="free_x") + theme_bw() + xlab("Sequencing Library") + ylab("% of Reads Passed Filtering") +guides(fill="none")+ geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted') +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted')+ggtitle(chain)
	p3 <- ggplot(Filtering_Results, aes(x = DISEASE, y=N_UniqueBCRS, fill=DISEASE)) + geom_boxplot(alpha=0.5)+ theme_bw() + xlab("DISEASE") + ylab("Number of  Unique VDJs") +guides(fill="none")+ geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) + sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) - sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted') +geom_hline(yintercept=mean(Filtering_Results$N_UniqueBCRS), col="red", linetype='dotted')+geom_hline(yintercept=1500, col="green")+ggtitle(chain)
	p4 <- ggplot(Filtering_Results, aes(x = DISEASE, y=PercentagePassedFiltering, fill=DISEASE)) + geom_boxplot(alpha=0.5)  + theme_bw() + xlab("DISEASE") + ylab("% of Reads Passed Filtering") +guides(fill="none")+ geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted') +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted')+ggtitle(chain)
	plot(plot_grid(p1, p3, p2, p4, labels = c('A', 'B', 'C', 'D'), ncol=2, rel_widths = c(4,1.3), align = "v", axis = "btlr"))
	plot_list <- plot_grid(p1, p3, p2, p4, ncol=4, labels = c('A', 'B', 'C', 'D'), rel_widths = c(4,1.3,4,1.3), align = "h", axis = "btlr")
	dev.off()
	## Calculate Mean depth and stdeviation 
	stats <- Filtering_Results %>% arrange(Library) %>% group_by(Lane)  %>% summarize(MeanUniqueBCRs = round(mean(N_UniqueBCRS, na.rm=TRUE),2), SDUniqueBCRs = round(sd(N_UniqueBCRS, na.rm=TRUE),2), NSDUniqueBCRs = round(sd(N_UniqueBCRS, na.rm=TRUE)/mean(N_UniqueBCRS, na.rm=TRUE),2), MeanPassed = round(mean(PercentagePassedFiltering, na.rm=TRUE),2), SDPassed = round(sd(PercentagePassedFiltering, na.rm=TRUE),2), NSDPassed = round(sd(PercentagePassedFiltering, na.rm=TRUE)/mean(PercentagePassedFiltering, na.rm=TRUE),2))
	colnames(stats) <- c("Lane", "Mean no. VDJs", "SD no. VDJs", "NSD no. VDJs", "Mean % Passed Filtering","SD % Passed Filtering", "NSD % Passed Filtering")
	pdf(paste0(plot_dir,'/Filtering_QC_stats', run_name, '.pdf'), width=15, height=6)
	grid.table(stats, rows = NULL)
	dev.off()	
	write.table(stats, paste0(stat_dir, "/Filtering_QC_stats.txt"), sep=",", row.names=FALSE)
	stats <- Filtering_Results %>% arrange(Library)   %>% summarize(MeanUniqueBCRs = round(mean(N_UniqueBCRS, na.rm=TRUE),2), SDUniqueBCRs = round(sd(N_UniqueBCRS, na.rm=TRUE),2), NSDUniqueBCRs = round(sd(N_UniqueBCRS, na.rm=TRUE)/mean(N_UniqueBCRS, na.rm=TRUE),2), MeanPassed = round(mean(PercentagePassedFiltering, na.rm=TRUE),2), SDPassed = round(sd(PercentagePassedFiltering, na.rm=TRUE),2), NSDPassed = round(sd(PercentagePassedFiltering, na.rm=TRUE)/mean(PercentagePassedFiltering, na.rm=TRUE),2))
	write.table(stats, paste0(stat_dir, "/Filtering_QC_stats_summary.txt"), sep=",", row.names=FALSE)
	return(plot_list)
}