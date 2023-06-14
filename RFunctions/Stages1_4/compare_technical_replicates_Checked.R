# Function to visualise the metrics for technical replicate for bulk BCR/TCR pipline to look for batch effect 
# Author: Lauren Overend 
# lauren.overend@oriel.ox.ac.uk 

#technical_replicates_file<-'/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/LEO_SEPSIS_BCR_ALL_technicals.txt'
#path_to_outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG'

suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(gridExtra))
suppressMessages(library(stringr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggforce))
#suppressMessages(library(Gviz))
suppressMessages(library(cowplot))
suppressMessages(library(gtools))
	
	
compare_technicals <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, technical_replicates_file=technical_replicates_file, cluster_nodes = 5){

	
	path <- path_to_outputdir
	#Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES/ANNOTATIONS")
	files <- list.files(path, full.name=TRUE)
	all_files <- files
	
	# Read in location of technical replicate 
	technical_samples <- read.delim(technical_replicates_file, sep="\t") 
	
	## Extract only files for technical replicates 
	files <- data.frame(grep('IsoTyper_chain_repertoire_statistics_file_', all_files, value=TRUE))
	files$sample <- str_split_fixed(files[,1], "IsoTyper_chain_repertoire_statistics_file_", 2)[,2]
	files$sample <- gsub(".txt", "", files$sample)
	files <- files[files$sample %in% technical_samples$SampleID,]
	files <- as.vector(files[,1])
	
	## Get Some Summary Stats! (clustering stats)
	############## CLUSTERING SUMMARY STATS
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
	Iso_Clustering <- merge(technical_samples, Iso_Clustering, by.x="SampleID", by.y="Sample")
	
	## Work out if TCR or BCR
	if(any(Iso_Clustering$Isotype %like% "TR")){
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
		Iso_Clustering <- Iso_Clustering_subset		
	}	
	
	## plot how metrics vary across measures for technicals 
	pdf(paste0(path_to_outputdir,'/Plots/Clustering_Results_QC_TECHNICAL_REPLICATES.pdf'),width=20, height=10)	
	coeff <- max(Iso_Clustering$N_Reads)/max(Iso_Clustering$Largest_cluster_percent)
	if(max(Iso_Clustering$Largest_cluster_percent)<1){
			coeff <- max(Iso_Clustering$N_Reads)*max(Iso_Clustering$Largest_cluster_percent)
	}
	p1 <- ggplot(Iso_Clustering, aes(x=SampleID)) +geom_bar(aes(y=Largest_cluster_percent), stat="identity", position=position_dodge(), fill="lightblue", colour="black") + geom_point(aes(y=N_Reads/coeff))+theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Largest Cluster (%)") +facet_wrap(~Isotype, scales="free_y") + scale_y_continuous(name="Largest Cluster %",  sec.axis = sec_axis(trans=~.*coeff, name="Read Depth")) + theme(axis.title.y = element_text(color = "blue", size=13),axis.title.y.right = element_text(color = "black", size=13))+ggtitle("Technical Replicates: No Library Size Correction")
	p1.1 <- ggplot(Iso_Clustering, aes(x=as.character(Lane), y=Largest_cluster_percent, fill=as.character(Lane))) +geom_boxplot(alpha=0.5)+ geom_point(aes(colour=log(N_Reads)))+theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Largest Cluster (%)")+ scale_colour_gradient2(low = "blue", mid = "yellow",high = "red" ) +facet_wrap(~Isotype) +xlab("Lane") + labs(fill="Lane", colour="Log(Number of Reads)") +ggtitle("Technical Replicates: No Library Size Correction")
	plot(plot_grid(p1, p1.1, ncol=2))
	coeff <- max(Iso_Clustering$N_Reads)/max(Iso_Clustering$Second_cluster_percent)
	if(max(Iso_Clustering$Second_cluster_percent)<1){
			coeff <- max(Iso_Clustering$N_Reads)*max(Iso_Clustering$Second_cluster_percent)
	}
	p2 <- ggplot(Iso_Clustering, aes(x=SampleID)) +geom_bar(aes(y=Second_cluster_percent), stat="identity" , position=position_dodge(), colour="black", fill="lightblue")+ geom_point(aes(y=N_Reads/coeff)) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Second Large Cluster (%)")  +facet_wrap(~Isotype, scales="free_y") + scale_y_continuous(name="2nd Largest Cluster %",  sec.axis = sec_axis(trans=~.*coeff, name="Read Depth")) + theme(axis.title.y = element_text(color = "blue", size=13),axis.title.y.right = element_text(color = "black", size=13))+ggtitle("Technical Replicates: No Library Size Correction")
	p2.1 <- ggplot(Iso_Clustering, aes(x=as.character(Lane), y=Second_cluster_percent, fill=as.character(Lane))) +geom_boxplot(alpha=0.5)+ geom_point(aes(colour=log(N_Reads)))+theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("2nd Largest Cluster (%)")+ scale_colour_gradient2(low = "blue", mid = "yellow",high = "red" ) +facet_wrap(~Isotype) +xlab("Lane") + labs(fill="Lane", colour="Log(Number of Reads)") +ggtitle("Technical Replicates: No Library Size Correction")
	plot(plot_grid(p2, p2.1, ncol=2))
	coeff <- max(Iso_Clustering$N_Reads)/max(Iso_Clustering$Cluster_Gini_Index)
	if(max(Iso_Clustering$Cluster_Gini_Index)<1){
			coeff <- max(Iso_Clustering$N_Reads)*max(Iso_Clustering$Cluster_Gini_Index)
	}
	p3 <- ggplot(Iso_Clustering, aes(x=SampleID)) +geom_bar(aes(y= Cluster_Gini_Index), stat="identity" , position=position_dodge(), colour="black", fill="lightblue") + geom_point(aes(y=N_Reads/coeff)) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Cluster Gini Index") +facet_wrap(~Isotype, scales="free_y") + scale_y_continuous(name="Cluster Gini Index",  sec.axis = sec_axis(trans=~.*coeff, name="Read Depth")) + theme(axis.title.y = element_text(color = "blue", size=13),axis.title.y.right = element_text(color = "black", size=13))+ggtitle("Technical Replicates: No Library Size Correction")
	p3.1 <- ggplot(Iso_Clustering, aes(x=as.character(Lane), y=Cluster_Gini_Index, fill=as.character(Lane))) +geom_boxplot(alpha=0.5)+ geom_point(aes(colour=log(N_Reads)))+theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Cluster Gini Index")+ scale_colour_gradient2(low = "blue", mid = "yellow",high = "red" ) +facet_wrap(~Isotype) +xlab("Lane") + labs(fill="Lane", colour="Log(Number of Reads)") +ggtitle("Technical Replicates: No Library Size Correction")
	plot(plot_grid(p3, p3.1, ncol=2))
	coeff <- max(Iso_Clustering$N_Reads)/max(Iso_Clustering$Vertex_Gini_Index)
	if(max(Iso_Clustering$Vertex_Gini_Index)<1){
			coeff <- max(Iso_Clustering$N_Reads)*max(Iso_Clustering$Vertex_Gini_Index)
	}
	p4 <- ggplot(Iso_Clustering, aes(x=SampleID)) +geom_bar(aes(, y= Vertex_Gini_Index), stat="identity" , position=position_dodge(), colour="black", fill="lightblue") + geom_point(aes(y=N_Reads/coeff))+theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Vertex Gini Index") +facet_wrap(~Isotype, scales="free_y") + scale_y_continuous(name="Vertex Gini Index",  sec.axis = sec_axis(trans=~.*coeff, name="Read Depth")) + theme(axis.title.y = element_text(color = "blue", size=13),axis.title.y.right = element_text(color = "black", size=13))+ggtitle("Technical Replicates: No Library Size Correction")
	p4.1 <- ggplot(Iso_Clustering, aes(x=as.character(Lane), y=Vertex_Gini_Index, fill=as.character(Lane))) +geom_boxplot(alpha=0.5)+ geom_point(aes(colour=log(N_Reads)))+theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Vertex Gini Index")+ scale_colour_gradient2(low = "blue", mid = "yellow",high = "red" ) +facet_wrap(~Isotype) +xlab("Lane") + labs(fill="Lane", colour="Log(Number of Reads)") +ggtitle("Technical Replicates: No Library Size Correction")
	plot(plot_grid(p4, p4.1, ncol=2))
	dev.off()	
	write.table(Iso_Clustering, paste0(path_to_outputdir, "/Summary/Clustering_Results_QC_TECHNICAL_REPLICATES.txt"), sep="\t")

	#### V Gene Usage
	files <- data.frame(grep('Gene_frequencies', all_files, value=TRUE))
	files$sample <- str_split_fixed(files[,1], "Gene_frequencies_", 2)[,2]
	files$sample <- gsub(".txt", "", files$sample)
	files <- files[files$sample %in% technical_samples$SampleID,]
	files <- as.vector(files[,1])
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
	frequency_variable_genes_id <- merge(frequency_variable_genes_id, technical_samples, by.x="Sample", by.y="SampleID")
	frequency_variable_genes_family <- aggregate(VJ_Results$Count, by=list(V_family=VJ_Results$V_family, Sample=VJ_Results$Sample), FUN=sum)
	colnames(frequency_variable_genes_family)[3] <- "count"
	frequency_variable_genes_family <- merge(frequency_variable_genes_family, technical_samples, by.x="Sample", by.y="SampleID")
	frequency_variable_genes_J <- aggregate(VJ_Results$Count, by=list(J_gene=VJ_Results$J_gene, Sample=VJ_Results$Sample), FUN=sum)
	colnames(frequency_variable_genes_J)[3] <- "count"
	frequency_variable_genes_J <- merge(frequency_variable_genes_J, technical_samples, by.x="Sample", by.y="SampleID")
	pdf(paste0(path_to_outputdir, '/Plots/V_Gene_Usage_QC_TECHNICAL.pdf'), width=20, height=15)
	p1 <- ggplot(frequency_variable_genes_id, aes(x=Sample, y=count, fill=V_Gene)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="V gene") +facet_wrap(~Lane, scales="free_x")+ggtitle("Technical Replicates: No Library Size Correction")
	p2 <- ggplot(frequency_variable_genes_family, aes(x=Sample, y=count, fill=V_family)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="V family")+facet_wrap(~Lane, scales="free_x")+ggtitle("Technical Replicates: No Library Size Correction")
	p3 <- ggplot(frequency_variable_genes_J, aes(x=Sample, y=count, fill=J_gene)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="J gene")+facet_wrap(~Lane, scales="free_x")+ggtitle("Technical Replicates: No Library Size Correction")
	plot(plot_grid(p1, ncol=1))
	plot(plot_grid(p2, ncol=1))
	plot(plot_grid(p3, ncol=1))
	dev.off()

	#### Constant Gene Usage
	files <- data.frame(grep('Constant_region_counts', all_files, value=TRUE))
	files$sample <- str_split_fixed(files[,1], "Constant_region_counts_", 2)[,2]
	files$sample <- gsub(".txt", "", files$sample)
	files <- files[files$sample %in% technical_samples$SampleID,]
	files <- as.vector(files[,1])
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
	Constant_Results_subset <- Constant_Results[Constant_Results$gene != "IGHD/M_mutated" & Constant_Results$gene != "ALL" & Constant_Results$gene != "class_switched" & Constant_Results$gene != "IGHD/M_unmutated",]
	Constant_Results_subset <- merge(Constant_Results_subset, technical_samples, by.x="Sample", by.y="SampleID")
	Constant_Results_mutation <- Constant_Results[Constant_Results$gene == "IGHD/M_mutated"  | Constant_Results$gene == "IGHD/M_unmutated",]
	Constant_Results_mutation <- merge(Constant_Results_mutation, technical_samples, by.x="Sample", by.y="SampleID")	
	pdf(paste0(path_to_outputdir,'/Plots/Constant_Region_Counts_QC_ISOTYPEUSAGE_TECHNICAL.pdf'), width=10, height=10)
	p1 <- ggplot(Constant_Results_subset, aes(x=as.character(Lane), y=percentage, fill=as.character(Lane))) + geom_boxplot(alpha=0.5) + geom_point(aes(colour=log(frequency))) +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Reads") +facet_wrap(~gene, scales = "free_x")+ scale_colour_gradient2(low = "blue", mid = "yellow",high = "red" )+ labs(fill="Lane", colour="Log(Number of Reads)")+ggtitle("Technical Replicates: No Library Size Correction") +xlab("Lane")
	plot(p1)
	if(dim(Constant_Results_mutation)[1] != 0){
		p2 <- ggplot(Constant_Results_mutation, aes(x=as.character(Lane), y=percent_IGM, fill=as.character(Lane))) + geom_boxplot(alpha=0.5) + geom_point(aes(colour=log(frequency))) +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% of IgM_D Reads") +facet_wrap(~gene, scales = "free_x")+ scale_colour_gradient2(low = "blue", mid = "yellow",high = "red" )+ labs(fill="Lane", colour="Log(Number of Reads)")+ggtitle("Technical Replicates: No Library Size Correction")+xlab("Lane")
		plot(p2)
	}
	dev.off()
	
	## Filtering Statsfiles <- list.files(path, full.name=TRUE)
	path <- paste0(path_to_outputdir, "/ORIENTATED_SEQUENCES")
	files <- list.files(path, full.name=TRUE)
	files <- data.frame(grep('Filtering_report', files, value=TRUE))
	files$sample <- str_split_fixed(files[,1], "Filtering_report_", 2)[,2]
	files$sample <- gsub(".txt", "", files$sample)
	files <- files[files$sample %in% technical_samples$SampleID,]
	files <- as.vector(files[,1])
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
	Filtering_Results <- merge(Filtering_Results, technical_samples, by.x="Sample", by.y="SampleID")
	
	pdf(paste0(path_to_outputdir,'/Plots/Filtering_QC_TECHNICAL.pdf'), width=10, height=10)
	p4 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentagePassedFiltering), y=PercentagePassedFiltering, colour=as.character(Lane))) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Passed Filtering") +geom_hline(yintercept=mean(Filtering_Results$PercentagePassedFiltering), col="red", linetype='dotted') +ylim(0,100) +xlab("Sample")  +geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) + sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentagePassedFiltering) - sd(Filtering_Results$PercentagePassedFiltering)), col="purple", linetype='dotted')+ labs(colour="Lane")
	p2 <- ggplot(Filtering_Results, aes(x=reorder(Sample, N_UniqueBCRS), y=N_UniqueBCRS, colour=as.character(Lane))) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +geom_hline(yintercept=1000, col="blue") +ylab("No. Unique Sequences Post-Filtering") +xlab("Sample")  +geom_hline(yintercept=mean(Filtering_Results$N_UniqueBCRS), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) + sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$N_UniqueBCRS) - sd(Filtering_Results$N_UniqueBCRS)), col="purple", linetype='dotted')+ labs(colour="Lane")
	p1 <- ggplot(Filtering_Results, aes(x=reorder(Sample, PercentageFunctional), y=PercentageFunctional, colour=as.character(Lane))) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Functional") +geom_hline(yintercept=mean(Filtering_Results$PercentageFunctional), col="red", linetype='dotted') +ylim(0,100)+xlab("Sample") +geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) + sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$PercentageFunctional) - sd(Filtering_Results$PercentageFunctional)), col="purple", linetype='dotted')+ labs(colour="Lane")
	p3 <- ggplot(Filtering_Results, aes(x=reorder(Sample, Duplication), y=Duplication, colour=as.character(Lane))) + geom_point()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab(" % Unique VDJ sequences") +geom_hline(yintercept=mean(Filtering_Results$Duplication), col="red", linetype='dotted') +geom_hline(yintercept=(mean(Filtering_Results$Duplication) + sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(Filtering_Results$Duplication) - sd(Filtering_Results$Duplication)), col="purple", linetype='dotted')+ylim(0,100)+xlab("Sample")+ labs(colour="Lane")
	plot(plot_grid(p1,p3, p2, p4, ncol=2))
	dev.off()	
	## DONE
	
}


		

