# Function to visualise the metrics for technical replicate for bulk BCR/TCR pipline to look for batch effect 
# Author: Lauren Overend 
# lauren.overend@oriel.ox.ac.uk 
# neat format for thesis 

#Rscript AnalysisStages1to4.R -o /well/immune-rep/shared/MISEQ/TEST_LEO_ANNA_BCR/ -r TESTING -g IGH -b R4RA_Anna_Batch_file_BCR.txt -t R4RA_Anna_BCR_technical.txt
#path_to_outputdir='/well/immune-rep/shared/MISEQ/TEST_LEO_ANNA_BCR/'
#run_name='TESTING' 
#path_to_layout='R4RA_Anna_Batch_file_BCR.txt'
#technical_replicates_file="R4RA_Anna_BCR_technical.txt"

compare_technicals_neat <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, technical_replicates_file=technical_replicates_file, plot_dir=plot_dir, path_to_layout = path_to_layout, stat_dir=stat_dir, cluster_nodes = 5){
	
	path <- path_to_outputdir
	#Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES/ANNOTATIONS")
	files <- list.files(path, full.name=TRUE)
	all_files <- files
	
	# Read in location of technical replicate 
	technical_samples <- read.delim(technical_replicates_file, sep="\t") 
	
	##--------------------------------------------------------------
	## Extract only files for technical replicates 
	files <- data.frame(grep('IsoTyper_chain_repertoire_statistics_file_', all_files, value=TRUE))
	files$sample <- str_split_fixed(files[,1], "IsoTyper_chain_repertoire_statistics_file_", 2)[,2]
	files$sample <- gsub(".txt", "", files$sample)
	files <- files[files$sample %in% technical_samples$SampleID,]
	files <- as.vector(files[,1])
	
	if(length(files)>=1){
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
		VJ_Results$V_family <- gsub("/OR15", "",VJ_Results$V_family)
		VJ_Results$V_family <- gsub("/OR16", "",VJ_Results$V_family) 
		VJ_Results$V_family <- gsub("/OR21", "",VJ_Results$V_family) 	
		
		frequency_variable_genes_id <- aggregate(VJ_Results$Count, by=list(V_Gene=VJ_Results$V_Gene, Sample=VJ_Results$Sample), FUN=sum)
		colnames(frequency_variable_genes_id)[3] <- "count"
		frequency_variable_genes_id <- merge(frequency_variable_genes_id, technical_samples, by.x="Sample", by.y="SampleID")
		frequency_variable_genes_family <- aggregate(VJ_Results$Count, by=list(V_family=VJ_Results$V_family, Sample=VJ_Results$Sample), FUN=sum)
		
		colnames(frequency_variable_genes_family)[3] <- "count"
		frequency_variable_genes_family <- merge(frequency_variable_genes_family, technical_samples, by.x="Sample", by.y="SampleID")
		frequency_variable_genes_J <- aggregate(VJ_Results$Count, by=list(J_gene=VJ_Results$J_gene, Sample=VJ_Results$Sample), FUN=sum)
		colnames(frequency_variable_genes_J)[3] <- "count"
		frequency_variable_genes_J <- merge(frequency_variable_genes_J, technical_samples, by.x="Sample", by.y="SampleID")
		
		percentage <- frequency_variable_genes_J %>% dplyr::group_by(Sample, J_gene) %>% dplyr::summarise(cnt = count) %>%  dplyr::mutate(percent = cnt / sum(cnt)*100, totalcnt =sum(cnt))
		percentage2 <- percentage  %>% dplyr::group_by(J_gene)  %>% dplyr::summarize(Meanpercentage = round(mean(percent, na.rm=TRUE),2), SDpercentage = round(sd(percent, na.rm=TRUE),2), NSDpercentage=round(sd(percent, na.rm=TRUE)/mean(percent, na.rm=TRUE),2))
		colnames(percentage2)[1] <- "Gene"
		percentagex <- frequency_variable_genes_family %>% dplyr::group_by(Sample, V_family) %>% dplyr::summarise(cnt = count) %>%  dplyr::mutate(percent = cnt / sum(cnt)*100, totalcnt =sum(cnt))
		percentage2x <- percentagex  %>% dplyr::group_by(V_family)  %>% dplyr::summarize(Meanpercentage = round(mean(percent, na.rm=TRUE),2), SDpercentage = round(sd(percent, na.rm=TRUE),2), NSDpercentage=round(sd(percent, na.rm=TRUE)/mean(percent, na.rm=TRUE),2))
		colnames(percentage2x)[1] <- "Gene"
		results <- rbind(percentage2, percentage2x)
		colnames(results) <- c("Gene", "Mean % of Repertoire", "SD % of Repertoire", "NSD % of Repertoire")
		pdf(paste0(plot_dir,'/Filtering_Vgene_stats.pdf'), width=10, height=10)
		grid.table(results, rows = NULL)
		dev.off()
		write.table(results, paste0(stat_dir, "/V_gene_stats.txt"), sep=",", row.names=FALSE)

		###############################
		layouts <- read.delim(path_to_layout, sep="\t", header=TRUE)
		layouts$SampleID <- gsub("BCR_", "", layouts$SampleID)
		percentage <- merge(percentage, layouts, by.x="Sample", by.y="SampleID")
		percentagex <- merge(percentagex, layouts, by.x="Sample", by.y="SampleID")
		
		percentage$Library <- as.numeric(percentage$Library )
		percentage$Library <- factor(percentage$Library, levels=c(1:max(percentage$Library)))   
		percentage$Lane <- as.numeric(percentage$Lane )
		percentage$Lane <- factor(percentage$Lane, levels=c(1:max(percentage$Lane)))
		percentage$totalcnt <- as.numeric(percentage$totalcnt)
		percentagex$Library <- as.numeric(percentagex$Library )
		percentagex$Library <- factor(percentagex$Library, levels=c(1:max(percentagex$Library)))   
		percentagex$Lane <- as.numeric(percentagex$Lane )
		percentagex$Lane <- factor(percentagex$Lane, levels=c(1:max(percentagex$Lane)))
		percentagex$totalcnt <- as.numeric(percentagex$totalcnt)

		 
		pdf(paste0(plot_dir, "/V_Gene_Usage_QC_TECHNICAL.pdf"), width=10, height=5)
		p1 <- ggplot(percentage, aes(x=J_gene, y=percent)) +geom_boxplot() +geom_point(aes(colour=totalcnt))+ facet_grid(~Lane)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("J gene") + ylab("% of Repertoire")+guides(fill="none")+ scale_colour_viridis_c()+labs(colour="Total Count")
		p2 <- ggplot(percentagex, aes(x=V_family, y=percent)) +geom_boxplot() +geom_point(aes(colour=totalcnt)) + facet_grid(~Lane)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("V gene family") + ylab("% of Repertoire")+guides(fill="none")+ scale_colour_viridis_c()+labs(colour="Total Count")
		plot(plot_grid(p1, p2, labels = c('A', 'B'), ncol=1))
		dev.off()
		
		
		##------------------------------------------
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
		Constant_Results_subset$Library <- as.numeric(Constant_Results_subset$Library )
		Constant_Results_subset$Library <- factor(Constant_Results_subset$Library, levels=c(1:max(Constant_Results_subset$Library)))   
		Constant_Results_subset$Lane <- as.numeric(Constant_Results_subset$Lane )
		Constant_Results_subset$Lane <- factor(Constant_Results_subset$Lane, levels=c(1:max(Constant_Results_subset$Lane)))
		pdf(paste0(plot_dir, "/Constant_Region_Counts_QC_ISOTYPEUSAGE_TECHNICAL.pdf"), width=8, height=8)
		p1 <- ggplot(Constant_Results_subset, aes(x=Lane, y=percentage, fill=as.character(Lane))) + geom_boxplot(alpha=0.5) + geom_point(aes(colour=totalreads)) +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Reads") +facet_wrap(~gene, scales = "free_x")+ scale_colour_viridis_c()+ labs(fill="Lane", colour="Total VDJs")+ggtitle("Isotype Usage across Technical Replicates") +xlab("Lane") + guides(fill="none")
		plot(p1)
		dev.off()
		stats <- Constant_Results_subset %>% arrange(Lane) %>% dplyr::group_by(Lane, gene)  %>% dplyr::summarize(Meanpercentage = round(mean(percentage, na.rm=TRUE),2), SDpercentage = round(sd(percentage, na.rm=TRUE),2), NSDpercentage=round(sd(percentage, na.rm=TRUE)/mean(percentage, na.rm=TRUE),2))
		colnames(stats) <- c("Lane", "Isotype", "Mean % of Repertoire", "SD % of Repertoire", "NSD % of Repertoire")
		stats2 <- Constant_Results_subset %>% arrange(Lane) %>% dplyr::group_by(gene)  %>% dplyr::summarize(Meanpercentage = round(mean(percentage, na.rm=TRUE),2), SDpercentage = round(sd(percentage, na.rm=TRUE),2),  NSDpercentage=round(sd(percentage, na.rm=TRUE)/mean(percentage, na.rm=TRUE),2))
		colnames(stats2) <- c("Isotype", "Mean % of Repertoire", "SD % of Repertoire", "NSD % of Repertoire")
		pdf(paste0(plot_dir,'/Filtering_QC_stats_iso', run_name, '.pdf'), width=10, height=18)
		grid.table(stats, rows = NULL)
		dev.off()
		pdf(paste0(plot_dir,'/Filtering_QC_stats_iso2', run_name, '.pdf'), width=10, height=5)
		grid.table(stats2, rows = NULL)
		dev.off()
		
		Constant_Results_subsetx <- unique(Constant_Results_subset[, c("Sample", "totalreads")])
		Constant_Results_subsetx$Type <- "TECHNICAL"
		stasx <- Constant_Results_subsetx %>% dplyr::summarize(Meanpercentage = round(mean(totalreads, na.rm=TRUE),2), SDpercentage = round(sd(totalreads, na.rm=TRUE),2), NSDpercentage = round(sd(totalreads, na.rm=TRUE)/mean(totalreads, na.rm=TRUE),2))
		pdf(paste0(plot_dir, "/TECHNICAL.pdf"), width=4, height=4)
		p1 <- ggplot(Constant_Results_subsetx, aes(x=Type, y=totalreads)) + geom_boxplot(alpha=0.5) + geom_point() +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("Total Redas")+ggtitle("Technical Replicate Read Depth") +xlab("") + guides(fill="none")
		plot(p1)
		dev.off()
		####################################################
		
		write.table(stats2, paste0(stat_dir, "/TECHNICALC_stats.txt"), sep=",", row.names=FALSE)
		print("Mean number of reads for technicals")
		print(mean(unique(Constant_Results$totalreads)))
		} else {
		 print("Only one technical to compare - not running")
		}
}


compare_technicals_neat_tcr <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, technical_replicates_file=technical_replicates_file, plot_dir=plot_dir, path_to_layout = path_to_layout, stat_dir=stat_dir, chain=chain, cluster_nodes = 5){
	
	path <- path_to_outputdir
	#Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES/ANNOTATIONS")
	files <- list.files(path, full.name=TRUE)
	all_files <- files
	
	# Read in location of technical replicate 
	technical_samples <- read.delim(technical_replicates_file, sep="\t") 
	
	##--------------------------------------------------------------
	## Extract only files for technical replicates 
	files <- data.frame(grep('IsoTyper_chain_repertoire_statistics_file_', all_files, value=TRUE))
	files$sample <- str_split_fixed(files[,1], "IsoTyper_chain_repertoire_statistics_file_", 2)[,2]
	files$sample <- gsub(".txt", "", files$sample)
	files <- files[files$sample %in% technical_samples$SampleID,]
	files <- as.vector(files[,1])
	if(length(files)>=1){
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
		VJ_Results$V_family <- gsub("/OR15", "",VJ_Results$V_family)
		VJ_Results$V_family <- gsub("/OR16", "",VJ_Results$V_family) 
		VJ_Results$V_family <- gsub("/OR21", "",VJ_Results$V_family) 	
		
		if(chain == "TCRA"){
			VJ_Results=VJ_Results[VJ_Results$V_Gene %like% "TRA",]
			VJ_Results=VJ_Results[VJ_Results$J_gene %like% "TRA",]
		}
		if(chain == "TCRB"){
			VJ_Results=VJ_Results[VJ_Results$V_Gene %like% "TRB",]
			VJ_Results=VJ_Results[VJ_Results$J_gene %like% "TRB",]
		}
		if(chain == "TCRG"){
			VJ_Results=VJ_Results[VJ_Results$V_Gene %like% "TRG",]
			VJ_Results=VJ_Results[VJ_Results$J_gene %like% "TRG",]
		}
		if(chain == "TCRD"){
			VJ_Results=VJ_Results[VJ_Results$V_Gene %like% "DV",]
			VJ_Results=VJ_Results[VJ_Results$J_gene %like% "TRD",]
		}
		
		
		frequency_variable_genes_id <- aggregate(VJ_Results$Count, by=list(V_Gene=VJ_Results$V_Gene, Sample=VJ_Results$Sample), FUN=sum)
		colnames(frequency_variable_genes_id)[3] <- "count"
		frequency_variable_genes_id <- merge(frequency_variable_genes_id, technical_samples, by.x="Sample", by.y="SampleID")
		frequency_variable_genes_family <- aggregate(VJ_Results$Count, by=list(V_Gene=VJ_Results$V_Gene, Sample=VJ_Results$Sample), FUN=sum)
		
		colnames(frequency_variable_genes_family)[3] <- "count"
		frequency_variable_genes_family <- merge(frequency_variable_genes_family, technical_samples, by.x="Sample", by.y="SampleID")
		frequency_variable_genes_J <- aggregate(VJ_Results$Count, by=list(J_gene=VJ_Results$J_gene, Sample=VJ_Results$Sample), FUN=sum)
		colnames(frequency_variable_genes_J)[3] <- "count"
		frequency_variable_genes_J <- merge(frequency_variable_genes_J, technical_samples, by.x="Sample", by.y="SampleID")
		
		percentage <- frequency_variable_genes_J %>% dplyr::group_by(Sample, J_gene) %>% dplyr::summarise(cnt = count) %>%  dplyr::mutate(percent = cnt / sum(cnt)*100, totalcnt =sum(cnt))
		percentage2 <- percentage  %>% dplyr::group_by(J_gene)  %>% dplyr::summarize(Meanpercentage = round(mean(percent, na.rm=TRUE),2), SDpercentage = round(sd(percent, na.rm=TRUE),2), NSDpercentage=round(sd(percent, na.rm=TRUE)/mean(percent, na.rm=TRUE),2))
		colnames(percentage2)[1] <- "Gene"
		percentagex <- frequency_variable_genes_family %>% dplyr::group_by(Sample, V_Gene) %>% dplyr::summarise(cnt = count) %>%  dplyr::mutate(percent = cnt / sum(cnt)*100, totalcnt =sum(cnt))
		percentage2x <- percentagex  %>% dplyr::group_by(V_Gene)  %>% dplyr::summarize(Meanpercentage = round(mean(percent, na.rm=TRUE),2), SDpercentage = round(sd(percent, na.rm=TRUE),2), NSDpercentage=round(sd(percent, na.rm=TRUE)/mean(percent, na.rm=TRUE),2))
		colnames(percentage2x)[1] <- "Gene"
		results <- rbind(percentage2, percentage2x)
		colnames(results) <- c("Gene", "Mean % of Repertoire", "SD % of Repertoire", "NSD % of Repertoire")
		pdf(paste0(plot_dir,'/Filtering_Vgene_stats.pdf'), width=10, height=10)
		grid.table(results, rows = NULL)
		dev.off()
		write.table(results, paste0(stat_dir, "/V_gene_stats.txt"), sep=",", row.names=FALSE)

		###############################
		layouts <- read.delim(path_to_layout, sep="\t", header=TRUE)
		layouts$SampleID <- gsub("BCR_", "", layouts$SampleID)
		percentage <- merge(percentage, layouts, by.x="Sample", by.y="SampleID")
		percentagex <- merge(percentagex, layouts, by.x="Sample", by.y="SampleID")
		
		percentage$Library <- as.numeric(percentage$Library )
		percentage$Library <- factor(percentage$Library, levels=c(1:max(percentage$Library)))   
		percentage$Lane <- as.numeric(percentage$Lane )
		percentage$Lane <- factor(percentage$Lane, levels=c(1:max(percentage$Lane)))
		percentage$totalcnt <- as.numeric(percentage$totalcnt)
		percentagex$Library <- as.numeric(percentagex$Library )
		percentagex$Library <- factor(percentagex$Library, levels=c(1:max(percentagex$Library)))   
		percentagex$Lane <- as.numeric(percentagex$Lane )
		percentagex$Lane <- factor(percentagex$Lane, levels=c(1:max(percentagex$Lane)))
		percentagex$totalcnt <- as.numeric(percentagex$totalcnt)

		 
		pdf(paste0(plot_dir, "/V_Gene_Usage_QC_TECHNICAL.pdf"), width=14, height=6)
		p1 <- ggplot(percentage, aes(x=J_gene, y=percent, fill=as.factor(Lane), colour=as.factor(Lane))) +geom_boxplot(alpha=0.5) +theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("J gene") + ylab("% of Repertoire")+labs(colour="Lane", fill="Lane") +ggtitle(paste0(chain, " J gene"))
		p2 <- ggplot(percentagex, aes(x=V_Gene, y=percent, fill=as.factor(Lane), colour=as.factor(Lane))) +geom_boxplot(alpha=0.5)  +theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("V gene family") + ylab("% of Repertoire")+labs(colour="Lane", fill="Lane")+ggtitle(paste0(chain, " V gene"))
		plot(plot_grid(p1, p2, labels = c('A', 'B'), ncol=1))
		dev.off()
		
		px <- plot_grid(p1, p2, ncol=1)
		
		if(chain == "TCRG" | chain == "TCRD"){
			px <- plot_grid(p1, p2, ncol=2)
		}
		##------------------------------------------
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
					if(chain == "TCRA"){
					output=output[output$gene %like% "TRA" | output$gene %like% "ALL",]
					}
					if(chain == "TCRB"){
					output=output[output$gene %like% "TRB" | output$gene %like% "ALL",]
					}
					if(chain == "TCRG"){
					output=output[output$gene %like% "TRG" | output$gene %like% "ALL",]
					}
					if(chain == "TCRD"){
					output=output[output$gene %like% "TRD" | output$gene %like% "ALL",]
					}
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
		Constant_Results_subset$Library <- as.numeric(Constant_Results_subset$Library )
		Constant_Results_subset$Library <- factor(Constant_Results_subset$Library, levels=c(1:max(Constant_Results_subset$Library)))   
		Constant_Results_subset$Lane <- as.numeric(Constant_Results_subset$Lane )
		Constant_Results_subset$Lane <- factor(Constant_Results_subset$Lane, levels=c(1:max(Constant_Results_subset$Lane)))
		pdf(paste0(plot_dir, "/Constant_Region_Counts_QC_ISOTYPEUSAGE_TECHNICAL.pdf"), width=8, height=8)
		p1x <- ggplot(Constant_Results_subset, aes(x=Lane, y=percentage, fill=as.character(Lane))) + geom_boxplot(alpha=0.5) + geom_point(aes(colour=totalreads)) +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("% Reads") +facet_wrap(~gene, scales = "free_x")+ scale_colour_viridis_c()+ labs(fill="Lane", colour="Unique VDJs")+ggtitle(paste0(chain)) +xlab("Lane") + guides(fill="none")
		plot(p1x)
		dev.off()
		stats <- Constant_Results_subset %>% arrange(Lane) %>% dplyr::group_by(Lane, gene)  %>% dplyr::summarize(Meanpercentage = round(mean(percentage, na.rm=TRUE),2), SDpercentage = round(sd(percentage, na.rm=TRUE),2), NSDpercentage=round(sd(percentage, na.rm=TRUE)/mean(percentage, na.rm=TRUE),2))
		colnames(stats) <- c("Lane", "Isotype", "Mean % of Repertoire", "SD % of Repertoire", "NSD % of Repertoire")
		stats2 <- Constant_Results_subset %>% arrange(Lane) %>% dplyr::group_by(gene)  %>% dplyr::summarize(Meanpercentage = round(mean(percentage, na.rm=TRUE),2), SDpercentage = round(sd(percentage, na.rm=TRUE),2),  NSDpercentage=round(sd(percentage, na.rm=TRUE)/mean(percentage, na.rm=TRUE),2))
		colnames(stats2) <- c("Isotype", "Mean % of Repertoire", "SD % of Repertoire", "NSD % of Repertoire")
		pdf(paste0(plot_dir,'/Filtering_QC_stats_iso', run_name, '.pdf'), width=10, height=18)
		grid.table(stats, rows = NULL)
		dev.off()
		pdf(paste0(plot_dir,'/Filtering_QC_stats_iso2', run_name, '.pdf'), width=10, height=5)
		grid.table(stats2, rows = NULL)
		dev.off()
		
		Constant_Results_subsetx <- unique(Constant_Results_subset[, c("Sample", "totalreads")])
		Constant_Results_subsetx$Type <- "TECHNICAL"
		stasx <- Constant_Results_subsetx %>% dplyr::summarize(Meanpercentage = round(mean(totalreads, na.rm=TRUE),2), SDpercentage = round(sd(totalreads, na.rm=TRUE),2), NSDpercentage = round(sd(totalreads, na.rm=TRUE)/mean(totalreads, na.rm=TRUE),2))
		pdf(paste0(plot_dir, "/TECHNICAL.pdf"), width=4, height=4)
		p1 <- ggplot(Constant_Results_subsetx, aes(x=Type, y=totalreads)) + geom_boxplot(alpha=0.5) + geom_point() +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +ylab("Total Redas")+ggtitle("Technical Replicate Read Depth") +xlab("") + guides(fill="none")
		plot(p1)
		dev.off()
		####################################################
		
		write.table(stats2, paste0(stat_dir, "/TECHNICALC_stats.txt"), sep=",", row.names=FALSE)
		print("Mean number of reads for technicals")
		print(mean(unique(Constant_Results$totalreads)))
	} else {
		 print("Only one technical to compare - not running")
	}
}