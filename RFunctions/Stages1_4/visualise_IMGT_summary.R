# Function to visualise the Filtering Reports for all files processed through the RBR Bulk BCR/TCR pipeline
# Author: Lauren Overend 

imgt_summary_bcr <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 3, productivity=NA){
	library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(imgt_summary_bcrdoParallel)
	library(gridExtra)
	library(purrr)
	library(seqinr)
	path <- path_to_outputdir
	Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT")
	files <- list.files(path, full.name=TRUE)
	if(is.na(productivity)){
		files <- grep('Summary', files, value=TRUE)
		files_not <- grep('_productive_', files, value=TRUE)
		files_not2 <- grep('_unproductive', files, value=TRUE)
		files <- files[!files %in% c(files_not, files_not2)]
	} 
	
	if(productivity=="productive"){	
		files <- grep('_productive_', files, value=TRUE)
	} 
	
	if(productivity=="unproductive"){	
		files <- grep('_unproductive', files, value=TRUE)
	} 
	
	
	cl <- cluster_nodes
	registerDoParallel(cl)
	#for(i in 1:length(files)){
	IMGT_Summary <- foreach(i = 1:length(files), .combine=rbind, .packages=c('tidyverse', 'seqinr', 'purrr')) %dopar% {
				filtered_path <- files[i]
				print(i)
				`%notin%` <- Negate(`%in%`)
				sample_id <- unlist(str_split(filtered_path, 'IMGT_SPLIT/'))[2]
				sample_id <- unlist(str_split(sample_id, '_1_Summary.txt'))[1]
				sample_id <- unlist(str_split(sample_id, 'IMGT_'))[2]
				output <- read.delim(filtered_path, sep="\t", header=FALSE)
				output$Sample <- sample_id
				no_unique_bcrs <- length(output[,1])
				output$unique_bcrs <- no_unique_bcrs
				read_header <- output$V2
				ids <- strsplit(read_header, "__", fixed=TRUE)
				z <- map(ids, 1)
				ids  <- matrix(unlist(z), ncol=1, byrow=TRUE)
				output$read_id <- ids 
				path_to_reduced <- paste0(path_to_outputdir, "/ORIENTATED_SEQUENCES/NETWORKS")
				fastas <- list.files(path_to_reduced, full.name=TRUE)
				fastas <- grep(sample_id, fastas, value=TRUE)
				fastas <- grep("Att", fastas, value=TRUE)
				ids_2 <- read.delim(fastas, sep="\t", header=FALSE)
				ids_2 <- ids_2$V1
				long_ids <- ids_2
				long_ids_split <- strsplit(long_ids, "__", fixed=TRUE)
				a <- map(long_ids_split, 1)
				b <- map(long_ids_split, 2)
				df_a <- matrix(unlist(a), ncol=1, byrow=TRUE)
				df_b <- matrix(unlist(b), ncol=1, byrow=TRUE)
				igh_split <- strsplit(df_b, "|", fixed=TRUE)
				c <- map(igh_split, 1)
				x <- map(igh_split, 2)
				x <- x[1]
				col_ids <- unlist(strsplit(unlist(x), "_", fixed=TRUE))
				all_cols <- c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHEP2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM")
				
				df_c <- matrix(unlist(c), ncol=1, byrow=TRUE)
				igh_split_2 <- strsplit(df_c, "_", fixed=TRUE)
				df_d <- matrix(unlist(igh_split_2), ncol=length(col_ids), byrow=TRUE)
				df_d <- data.frame(df_d)
				colnames(df_d) <- col_ids
				
				if(any(all_cols %notin% col_ids)=="TRUE"){
					new_cols <- all_cols[all_cols %notin% col_ids]
					for(i in 1:length(new_cols)){
						name <- new_cols[i]
						df_d$V3 <- "0"
						colnames(df_d)[colnames(df_d) == 'V3'] <- name
					}
				}
				df_d <- data.frame(df_d)
				df_d <- as.data.frame(sapply(df_d, as.numeric))
				df_d <- df_d %>% select(order(colnames(df_d))) 
				df_d$count <- rowSums(df_d)
				colnames(df_a) <- "read_id"
				usage <- cbind(df_a, df_d)
				no_nonunique_bcrs <- sum(usage$count)
				output$no_nonunique_bcr <- no_nonunique_bcrs
				output_final <- merge(output, usage, by="read_id")
				return(output_final)	
	}
	stopImplicitCluster()
	IMGT_Summary$Sample <- as.factor(IMGT_Summary$Sample)
	names_1 <- c("read_id", "Sequence number","Sequence_ID","V_DOMAIN_Functionality", "V_GENE_allele", "V_REGION score", "V_REGION_identity_percent", "V_REGION_identity_nt", "V_REGION_identity_percent_with ins_del_events", "V_REGION_identity_nt(with ins/del events)", "J_GENE_and_allele", "J_REGION_score", "J_REGION_identity_percent", "J_REGION_identity_nt", "D_GENE_and_allele", "D_REGION_reading_frame", "CDR1_IMGT_length", "CDR2_IMGT_length","CDR3_IMGT_length", "CDR_IMGT_lengths", "FR_IMGT_lengths", "AA_JUNCTION", "JUNCTION_frame", "Orientation", "V_DOMAIN_Functionality_comment", "V_REGION_potential_ins_del", "J_GENE_and_allele_comment", "V_REGION_insertions", "V_REGION_deletions", "Sequence", "5prime_trimmed_n_nb", "3prime_trimmed_n_nb", "Analysed_sequence_length", "Sequence_analysis_category", "Sample", "No_unique_bcrs", "No_non_unique_bcr", "IGHA1", "IGHA2", "IGHD", "IGHE", "IGHEP2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM", "Sequence_count")
	colnames(IMGT_Summary) <- names_1
	
	IMGT_Summary$Sample <- gsub("BCR1.", "", IMGT_Summary$Sample)
	IMGT_Summary$Sample <- gsub("BCR_._", "", IMGT_Summary$Sample)
	IMGT_Summary$Sample <- gsub("BCR_.._", "", IMGT_Summary$Sample)


	if(dir.exists(paste0(path_to_outputdir, "/Plots/IMGT"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots/IMGT"))
	}
	
	
	
	
	# Extract a V gene
	#-----------------------------------------------------------------------------
	IMGT_Summary$V_GENE_allele[IMGT_Summary$V_GENE_allele==""] <- "NO_ASSIGNMENT"
	genes <- strsplit(IMGT_Summary$V_GENE_allele, ",", fixed=TRUE)
	w <- map(genes, 1) 
	df  <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df <- gsub("Homsap ", "", df)
	genes <- strsplit(df, " F", fixed=TRUE)
	w <- map(genes, 1) 
	df  <- matrix(unlist(w), ncol=1, byrow=TRUE)	
	# Get just V gene Family
    genes <- strsplit(df, "*", fixed=TRUE)   
	w <- map(genes, 1)
    df_1 <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df_3 <- cbind(df, df_1)
	genes <- strsplit(df_1, "-", fixed=TRUE)   
	w <- map(genes, 1)
    df_4 <- matrix(unlist(w), ncol=1, byrow=TRUE)
	
	df_5 <- cbind(df_3, df_4)
	
	colnames(df_5) <- c("V_GENE_ALLELE_selected", "V_Gene_selected", "V_Gene_Family")
	IMGT_Summary <- cbind(IMGT_Summary, df_5)
	
    # Extract a J gene 	
	#------------------------------------------------------------------------------			
	IMGT_Summary$J_GENE_and_allele[IMGT_Summary$J_GENE_and_allele==""] <- "NO_ASSIGNMENT"
	genes <- strsplit(IMGT_Summary$J_GENE_and_allele, ",", fixed=TRUE)
	w <- map(genes, 1) 
	df  <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df <- gsub("Homsap ", "", df)
	genes <- strsplit(df, " F", fixed=TRUE)
	w <- map(genes, 1) 
	df  <- matrix(unlist(w), ncol=1, byrow=TRUE)	
	# Get just V gene Family
    genes <- strsplit(df, "*", fixed=TRUE)   
	w <- map(genes, 1)
    df_1 <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df_3 <- cbind(df, df_1)
	colnames(df_3) <- c("J_GENE_ALLELE_selected", "J_Gene_selected")
	IMGT_Summary <- cbind(IMGT_Summary, df_3)
	
	# Extract a D gene 	
	#------------------------------------------------------------------------------			
	IMGT_Summary$D_GENE_and_allele[IMGT_Summary$D_GENE_and_allele==""] <- "NO_ASSIGNMENT"
	genes <- strsplit(IMGT_Summary$D_GENE_and_allele, ",", fixed=TRUE)
	w <- map(genes, 1) 
	df  <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df <- gsub("Homsap ", "", df)
	genes <- strsplit(df, " F", fixed=TRUE)
	w <- map(genes, 1) 
	df  <- matrix(unlist(w), ncol=1, byrow=TRUE)	
	# Get just V gene Family
    genes <- strsplit(df, "*", fixed=TRUE)   
	w <- map(genes, 1)
    df_1 <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df_3 <- cbind(df, df_1)
	colnames(df_3) <- c("D_GENE_ALLELE_selected", "D_Gene_selected")
	IMGT_Summary <- cbind(IMGT_Summary, df_3)
	
	##-------------------------------------------------------------------------			
	IMGT_Summary$CDR3_IMGT_length <- as.numeric(IMGT_Summary$CDR3_IMGT_length)
	IMGT_Summary$CDR1_IMGT_length <- as.numeric(IMGT_Summary$CDR1_IMGT_length)
	IMGT_Summary$CDR2_IMGT_length <- as.numeric(IMGT_Summary$CDR2_IMGT_length)
	IMGT_Summary$Analysed_sequence_length <- as.numeric(IMGT_Summary$Analysed_sequence_length)
	pdf(paste0(path_to_outputdir,'/Plots/IMGT/IMGT_Summary_', Run_name, "_", productivity, '.pdf'), width=23, height=14)
	p1 <- ggplot(IMGT_Summary, aes(x=Sample, y=as.numeric(CDR3_IMGT_length))) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") +geom_hline(yintercept=mean(IMGT_Summary$CDR3_IMGT_length[!is.na(IMGT_Summary$CDR3_IMGT_length)]), color="red") +ylab("CDR3 Length")
	p1.1 <- ggplot(IMGT_Summary, aes(x=as.numeric(CDR3_IMGT_length))) + geom_histogram(binwidth=1, color="black", fill="white") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("CDR3 IMGT length")+geom_vline(xintercept=mean(IMGT_Summary$CDR3_IMGT_length[!is.na(IMGT_Summary$CDR3_IMGT_length)]), color="red") 
	p2 <- ggplot(IMGT_Summary, aes(x=Sample, y=as.numeric(CDR1_IMGT_length))) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") +geom_hline(yintercept=mean(IMGT_Summary$CDR1_IMGT_length[!is.na(IMGT_Summary$CDR1_IMGT_length)]), color="red") +ylab("CDR1 Length")
	p2.1 <- ggplot(IMGT_Summary, aes(x=as.numeric(CDR1_IMGT_length))) + geom_histogram(binwidth=1, color="black", fill="white") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("CDR1 IMGT length")+geom_vline(xintercept=mean(IMGT_Summary$CDR1_IMGT_length[!is.na(IMGT_Summary$CDR1_IMGT_length)]), color="red") 
  	p3 <- ggplot(IMGT_Summary, aes(x=Sample, y=as.numeric(CDR2_IMGT_length))) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") +geom_hline(yintercept=mean(IMGT_Summary$CDR2_IMGT_length[!is.na(IMGT_Summary$CDR2_IMGT_length)]), color="red") +ylab("CDR2 Length")
	p3.1 <- ggplot(IMGT_Summary, aes(x=as.numeric(CDR2_IMGT_length))) + geom_histogram(binwidth=1, color="black", fill="white") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("CDR2 IMGT length")+geom_vline(xintercept=mean(IMGT_Summary$CDR2_IMGT_length[!is.na(IMGT_Summary$CDR2_IMGT_length)]), color="red") 
	plot(plot_grid(p2, p2.1, ncol=1))
	plot(plot_grid(p3, p3.1, ncol=1))   
	plot(plot_grid(p1, p1.1, ncol=1))
	p4 <- ggplot(IMGT_Summary, aes(x=Sample, y=as.numeric(Analysed_sequence_length))) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") +geom_hline(yintercept=mean(IMGT_Summary$Analysed_sequence_length[!is.na(IMGT_Summary$Analysed_sequence_length)]), color="red") +ylab("Analysed Sequence Length")
	p4.1 <- ggplot(IMGT_Summary, aes(x=as.numeric(Analysed_sequence_length))) + geom_histogram(binwidth=2, color="black", fill="white") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("Analysed Sequence Length")+geom_vline(xintercept=mean(IMGT_Summary$Analysed_sequence_length[!is.na(IMGT_Summary$Analysed_sequence_length)]), color="red") 
	plot(plot_grid(p4, p4.1, ncol=1))
	p4 <- ggplot(IMGT_Summary, aes(x=Sample, y=as.numeric(V_REGION_identity_percent))) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") +geom_hline(yintercept=mean(IMGT_Summary$V_REGION_identity_percent[!is.na(IMGT_Summary$V_REGION_identity_percent)]), color="red") +ylab("V REGION identity %")
	p4.1 <- ggplot(IMGT_Summary, aes(x=as.numeric(V_REGION_identity_percent))) + geom_histogram(binwidth=1, color="black", fill="white") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("V REGION identity %")+geom_vline(xintercept=mean(IMGT_Summary$V_REGION_identity_percent[!is.na(IMGT_Summary$V_REGION_identity_percent)]), color="red") 
	plot(plot_grid(p4, p4.1, ncol=1))
	p4 <- ggplot(IMGT_Summary, aes(x=Sample, y=as.numeric(J_REGION_identity_percent))) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") +geom_hline(yintercept=mean(IMGT_Summary$J_REGION_identity_percent[!is.na(IMGT_Summary$J_REGION_identity_percent)]), color="red") +ylab("J REGION identity %")
	p4.1 <- ggplot(IMGT_Summary, aes(x=as.numeric(J_REGION_identity_percent))) + geom_histogram(binwidth=1, color="black", fill="white") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("J REGION identity %")+geom_vline(xintercept=mean(IMGT_Summary$J_REGION_identity_percent[!is.na(IMGT_Summary$J_REGION_identity_percent)]), color="red") 
	plot(plot_grid(p4, p4.1, ncol=1))
	dev.off()
	
	V_genes <- data.frame(table(IMGT_Summary$Sample, IMGT_Summary$V_Gene_selected))
	colnames(V_genes) <- c("Sample", "IGV_V", "Freq_V")
	V_genes_family <- data.frame(table(IMGT_Summary$Sample, IMGT_Summary$V_Gene_Family)) 
	colnames(V_genes_family) <- c("Sample", "IGV_V_fam", "Freq_V_fam")
	J_genes <- data.frame(table(IMGT_Summary$Sample, IMGT_Summary$J_Gene_selected)) 
	colnames(J_genes) <- c("Sample", "IGV_J", "Freq_J")
	D_genes <- data.frame(table(IMGT_Summary$Sample, IMGT_Summary$D_Gene_selected)) 
	colnames(D_genes) <- c("Sample", "IGV_D", "Freq_D")
	
	summary_usage <- merge(V_genes, V_genes_family, by="Sample")
	summary_usage <- merge(summary_usage, J_genes, by="Sample")
	summary_usage <- merge(summary_usage, D_genes, by="Sample")
	
	pdf(paste0(path_to_outputdir, '/Plots/IMGT/Gene_Usage_QC_', Run_name,"_", productivity, '.pdf'), width=23, height=14)
	p1 <- ggplot(V_genes, aes(x=Sample, y=Freq_V, fill=IGV_V)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="V gene")
	p2 <- ggplot(V_genes_family, aes(x=Sample, y=Freq_V_fam, fill=IGV_V_fam)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="V family")
	p3 <- ggplot(J_genes, aes(x=Sample, y=Freq_J, fill=IGV_J)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="J gene")
	p4 <- ggplot(D_genes, aes(x=Sample, y=Freq_D, fill=IGV_D)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="D gene")
	plot(plot_grid(p1, ncol=1))
	plot(plot_grid(p2, ncol=1))
	plot(plot_grid(p3, ncol=1))
	plot(plot_grid(p4, ncol=1))
	dev.off()

	IMGT_counts <- IMGT_Summary[, c("No_unique_bcrs", "No_non_unique_bcr", "Sample")]
	IMGT_counts <- unique(IMGT_counts)
	pdf(paste0(path_to_outputdir, '/Plots/IMGT/Counts_QC_', Run_name,"_", productivity, '.pdf'), width=23, height=14)
	p1 <- ggplot(IMGT_counts, aes(x=Sample, y=as.numeric(No_unique_bcrs))) +geom_bar(stat="identity", fill="white", color="black")  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Number of Unique BCRs (IMGT))") +geom_hline(yintercept=1000, color="blue") +geom_hline(yintercept=(mean(IMGT_counts$No_unique_bcrs) + sd(IMGT_counts$No_unique_bcrs)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(IMGT_counts$No_unique_bcrs) - sd(IMGT_counts$No_unique_bcrs)), col="purple", linetype='dotted') +geom_hline(yintercept=mean(IMGT_counts$No_unique_bcrs), col="red", linetype='dotted')
	p2 <- ggplot(IMGT_counts, aes(x=Sample, y=as.numeric(No_non_unique_bcr))) +geom_bar(stat="identity", fill="white", color="black")  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Number of Unique UMIs (IMGT))") +geom_hline(yintercept=1000, color="blue") +geom_hline(yintercept=(mean(IMGT_counts$No_non_unique_bcr) + sd(IMGT_counts$No_non_unique_bcr)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(IMGT_counts$No_non_unique_bcr) - sd(IMGT_counts$No_non_unique_bcr)), col="purple", linetype='dotted') +geom_hline(yintercept=mean(IMGT_counts$No_non_unique_bcr), col="red", linetype='dotted')
	plot(plot_grid(p2, p1, ncol=1))
	dev.off()	
	
	IMGT_Summary$V_PRODUCTIVITY_EDITED  <- NA
	IMGT_Summary$V_PRODUCTIVITY_EDITED[IMGT_Summary$V_DOMAIN_Functionality == "productive (see comment)" | IMGT_Summary$V_DOMAIN_Functionality == "productive" ] <- "Productive" 
	IMGT_Summary$V_PRODUCTIVITY_EDITED[IMGT_Summary$V_DOMAIN_Functionality == "unproductive (see comment)" | IMGT_Summary$V_DOMAIN_Functionality == "unproductive" ] <- "Unproductive" 
	IMGT_Summary$V_PRODUCTIVITY_EDITED[IMGT_Summary$V_DOMAIN_Functionality == "No rearrangement found" | IMGT_Summary$V_DOMAIN_Functionality == "No results" ] <- "No Rearrangement/Result" 
	IMGT_Summary$V_PRODUCTIVITY_EDITED[IMGT_Summary$V_DOMAIN_Functionality == "rearranged sequence (but no junction found) (see comment)" ] <- "Rearrangement No Junction" 
	
	IMGT_functionality <- IMGT_Summary[, c("Sample", "V_PRODUCTIVITY_EDITED")]
	table_1 <- data.frame(table(IMGT_functionality))
	
	pdf(paste0(path_to_outputdir, '/Plots/IMGT/Functionality_QC_', Run_name,"_", productivity, '.pdf'), width=30, height=14)
    p1 <- ggplot(IMGT_Summary, aes(x=Sample, y=Analysed_sequence_length, color=V_PRODUCTIVITY_EDITED)) +geom_boxplot()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ ylab("Analysed Sequence Length")
	p2 <- ggplot(IMGT_Summary, aes(y=Analysed_sequence_length, x=V_PRODUCTIVITY_EDITED, color=V_Gene_Family)) +geom_boxplot()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("V Gene Functionality") + ylab("Analysed Sequence Length")
	p3 <- ggplot(table_1, aes(x=Sample, y=Freq, fill=V_PRODUCTIVITY_EDITED)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="Functionality")
	plot(plot_grid(p1, p2, p3,  ncol=1))
	dev.off()
	
	if(dir.exists(paste0(path_to_outputdir, "/Summary/IMGT"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary/IMGT"))
	}				
	write.table(IMGT_Summary, paste0(path_to_outputdir, "/Summary/IMGT/IMGT_SUMMARY_", Run_name,"_", productivity, ".txt"), sep="\t")
	write.table(summary_usage, paste0(path_to_outputdir, "/Summary/IMGT/IMGT_SUMMARY_GENE_USAGE_", Run_name,"_", productivity, ".txt"), sep="\t")
	return(IMGT_Summary)
}




# Function to visualise the Filtering Reports for all files processed through the RBR Bulk BCR/TCR pipeline
# Author: Lauren Overend 

imgt_summary_tcr <- function(path_to_outputdir = path_to_outputdir, run_name = run_name, cluster_nodes = 10){
	library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(doParallel)
	library(gridExtra)
	library(purrr)
	library(seqinr)
	path <- path_to_outputdir
	Run_name <- run_name
	path <- paste0(path, "/ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT")
	files <- list.files(path, full.name=TRUE)
	files <- grep('Summary', files, value=TRUE)
	cl <- cluster_nodes
	registerDoParallel(cl)
	IMGT_Summary <- foreach(i = 1:length(files), .combine=rbind, .packages=c('tidyverse', 'seqinr', 'purrr')) %dopar% {
				filtered_path <- files[i]
				print(i)
				`%notin%` <- Negate(`%in%`)
				sample_id <- unlist(str_split(filtered_path, 'IMGT_SPLIT/'))[2]
				sample_id <- unlist(str_split(sample_id, '_1_Summary.txt'))[1]
				sample_id <- unlist(str_split(sample_id, 'IMGT_'))[2]
				output <- read.delim(filtered_path, sep="\t", header=FALSE)
				output$Sample <- sample_id
				no_unique_bcrs <- length(output[,1])
				output$unique_bcrs <- no_unique_bcrs
				read_header <- output$V2
				ids <- strsplit(read_header, "__", fixed=TRUE)
				z <- map(ids, 1)
				ids  <- matrix(unlist(z), ncol=1, byrow=TRUE)
				output$read_id <- ids 
				path_to_reduced <- paste0(path_to_outputdir, "/ORIENTATED_SEQUENCES/NETWORKS")
				fastas <- list.files(path_to_reduced, full.name=TRUE)
				fastas <- grep(sample_id, fastas, value=TRUE)
				fastas <- grep("Att", fastas, value=TRUE)
				ids_2 <- read.delim(fastas, sep="\t", header=FALSE)
				long_ids <- ids_2$V1
				long_ids_split <- strsplit(long_ids, "__", fixed=TRUE)
				a <- map(long_ids_split, 1)
				b <- map(long_ids_split, 2)
				df_a <- matrix(unlist(a), ncol=1, byrow=TRUE)
				df_b <- matrix(unlist(b), ncol=1, byrow=TRUE)
				igh_split <- strsplit(df_b, "|", fixed=TRUE)
				c <- map(igh_split, 1)
				x <- map(igh_split, 2)
				x <- x[1]
				col_ids <- unlist(strsplit(unlist(x), "_", fixed=TRUE))
				all_cols <- c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2")
				
				df_c <- matrix(unlist(c), ncol=1, byrow=TRUE)
				igh_split_2 <- strsplit(df_c, "_", fixed=TRUE)
				df_d <- matrix(unlist(igh_split_2), ncol=length(col_ids), byrow=TRUE)
				df_d <- data.frame(df_d)
				colnames(df_d) <- col_ids
				
				if(any(all_cols %notin% col_ids)=="TRUE"){
					new_cols <- all_cols[all_cols %notin% col_ids]
					for(i in 1:length(new_cols)){
						name <- new_cols[i]
						df_d$V3 <- "0"
						colnames(df_d)[colnames(df_d) == 'V3'] <- name
					}
				}
				
				df_d <- data.frame(df_d)
				df_d <- as.data.frame(sapply(df_d, as.numeric))
				df_d <- df_d %>% select(order(colnames(df_d))) 
				df_d$count <- rowSums(df_d)
				colnames(df_a) <- "read_id"
				usage <- cbind(df_a, df_d)
				no_nonunique_bcrs <- sum(usage$count)
				output$no_nonunique_bcr <- no_nonunique_bcrs
				output_final <- merge(output, usage, by="read_id")
				return(output_final)	
	}
	stopImplicitCluster()
	IMGT_Summary$Sample <- as.factor(IMGT_Summary$Sample)
	names_1 <- c("read_id", "Sequence number","Sequence_ID","V_DOMAIN_Functionality", "V_GENE_allele", "V_REGION score", "V_REGION_identity_percent", "V_REGION_identity_nt", "V_REGION_identity_percent_with ins_del_events", "V_REGION_identity_nt(with ins/del events)", "J_GENE_and_allele", "J_REGION_score", "J_REGION_identity_percent", "J_REGION_identity_nt", "D_GENE_and_allele", "D_REGION_reading_frame", "CDR1_IMGT_length", "CDR2_IMGT_length","CDR3_IMGT_length", "CDR_IMGT_lengths", "FR_IMGT_lengths", "AA_JUNCTION", "JUNCTION_frame", "Orientation", "V_DOMAIN_Functionality_comment", "V_REGION_potential_ins_del", "J_GENE_and_allele_comment", "V_REGION_insertions", "V_REGION_deletions", "Sequence", "5prime_trimmed_n_nb", "3prime_trimmed_n_nb", "Analysed_sequence_length", "Sequence_analysis_category", "Sample", "No_unique_tcrs", "No_non_unique_tcr", "TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "Sequence_count")
	colnames(IMGT_Summary) <- names_1
	
	if(dir.exists(paste0(path_to_outputdir, "/Plots/IMGT"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Plots/IMGT"))
	}
	
	# Extract a V gene
	#-----------------------------------------------------------------------------
	IMGT_Summary$V_GENE_allele[IMGT_Summary$V_GENE_allele==""] <- "NO_ASSIGNMENT"
	genes <- strsplit(IMGT_Summary$V_GENE_allele, ",", fixed=TRUE)
	w <- map(genes, 1) 
	df  <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df <- gsub("Homsap ", "", df)
	df <- gsub(" [[:punct:]]F[[:punct:]]", "", df)
	df <- gsub(" F", "", df)
	# Get just V gene Family
    genes <- strsplit(df, "*", fixed=TRUE)   
	w <- map(genes, 1)
    df_1 <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df_3 <- cbind(df, df_1)
	genes <- strsplit(df_1, "-", fixed=TRUE)   
	w <- map(genes, 1)
    df_4 <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df_5 <- cbind(df_3, df_4)
	
	colnames(df_5) <- c("V_GENE_ALLELE_selected", "V_Gene_selected", "V_Gene_Family")
	IMGT_Summary <- cbind(IMGT_Summary, df_5)
	
    # Extract a J gene 	
	#------------------------------------------------------------------------------			
	IMGT_Summary$J_GENE_and_allele[IMGT_Summary$J_GENE_and_allele==""] <- "NO_ASSIGNMENT"
	genes <- strsplit(IMGT_Summary$J_GENE_and_allele, ",", fixed=TRUE)
	w <- map(genes, 1) 
	df  <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df <- gsub("Homsap ", "", df)
	df <- gsub(" [[:punct:]]F[[:punct:]]", "", df)
	df <- gsub(" F", "", df)
	# Get just V gene Family
    genes <- strsplit(df, "*", fixed=TRUE)   
	w <- map(genes, 1)
    df_1 <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df_3 <- cbind(df, df_1)
	colnames(df_3) <- c("J_GENE_ALLELE_selected", "J_Gene_selected")
	IMGT_Summary <- cbind(IMGT_Summary, df_3)
	
	# Extract a D gene 	
	#------------------------------------------------------------------------------			
	IMGT_Summary$D_GENE_and_allele[IMGT_Summary$D_GENE_and_allele==""] <- "NO_ASSIGNMENT"
	genes <- strsplit(IMGT_Summary$D_GENE_and_allele, ",", fixed=TRUE)
	w <- map(genes, 1) 
	df  <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df <- gsub("Homsap ", "", df)
	df <- gsub(" [[:punct:]]F[[:punct:]]", "", df)
	df <- gsub(" F", "", df)	
	# Get just V gene Family
    genes <- strsplit(df, "*", fixed=TRUE)   
	w <- map(genes, 1)
    df_1 <- matrix(unlist(w), ncol=1, byrow=TRUE)
	df_3 <- cbind(df, df_1)
	colnames(df_3) <- c("D_GENE_ALLELE_selected", "D_Gene_selected")
	IMGT_Summary <- cbind(IMGT_Summary, df_3)
	
	##-------------------------------------------------------------------------			
	IMGT_Summary$CDR3_IMGT_length <- as.numeric(IMGT_Summary$CDR3_IMGT_length)
	IMGT_Summary$CDR1_IMGT_length <- as.numeric(IMGT_Summary$CDR1_IMGT_length)
	IMGT_Summary$CDR2_IMGT_length <- as.numeric(IMGT_Summary$CDR2_IMGT_length)
	IMGT_Summary$Analysed_sequence_length <- as.numeric(IMGT_Summary$Analysed_sequence_length)
	pdf(paste0(path_to_outputdir,'/Plots/IMGT/IMGT_Summary_', Run_name, '.pdf'), width=23, height=14)
	p1 <- ggplot(IMGT_Summary, aes(x=Sample, y=as.numeric(CDR3_IMGT_length))) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") +geom_hline(yintercept=mean(IMGT_Summary$CDR3_IMGT_length[!is.na(IMGT_Summary$CDR3_IMGT_length)]), color="red") +ylab("CDR3 Length")
	p1.1 <- ggplot(IMGT_Summary, aes(x=as.numeric(CDR3_IMGT_length))) + geom_histogram(binwidth=1, color="black", fill="white") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("CDR3 IMGT length")+geom_vline(xintercept=mean(IMGT_Summary$CDR3_IMGT_length[!is.na(IMGT_Summary$CDR3_IMGT_length)]), color="red") 
	p2 <- ggplot(IMGT_Summary, aes(x=Sample, y=as.numeric(CDR1_IMGT_length))) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") +geom_hline(yintercept=mean(IMGT_Summary$CDR1_IMGT_length[!is.na(IMGT_Summary$CDR1_IMGT_length)]), color="red") +ylab("CDR1 Length")
	p2.1 <- ggplot(IMGT_Summary, aes(x=as.numeric(CDR1_IMGT_length))) + geom_histogram(binwidth=1, color="black", fill="white") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("CDR1 IMGT length")+geom_vline(xintercept=mean(IMGT_Summary$CDR1_IMGT_length[!is.na(IMGT_Summary$CDR1_IMGT_length)]), color="red") 
  	p3 <- ggplot(IMGT_Summary, aes(x=Sample, y=as.numeric(CDR2_IMGT_length))) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") +geom_hline(yintercept=mean(IMGT_Summary$CDR2_IMGT_length[!is.na(IMGT_Summary$CDR2_IMGT_length)]), color="red") +ylab("CDR2 Length")
	p3.1 <- ggplot(IMGT_Summary, aes(x=as.numeric(CDR2_IMGT_length))) + geom_histogram(binwidth=1, color="black", fill="white") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("CDR2 IMGT length")+geom_vline(xintercept=mean(IMGT_Summary$CDR2_IMGT_length[!is.na(IMGT_Summary$CDR2_IMGT_length)]), color="red") 
	plot(plot_grid(p2, p2.1, ncol=1))
	plot(plot_grid(p3, p3.1, ncol=1))   
	plot(plot_grid(p1, p1.1, ncol=1))
	p4 <- ggplot(IMGT_Summary, aes(x=Sample, y=as.numeric(Analysed_sequence_length))) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") +geom_hline(yintercept=mean(IMGT_Summary$Analysed_sequence_length[!is.na(IMGT_Summary$Analysed_sequence_length)]), color="red") +ylab("Analysed Sequence Length")
	p4.1 <- ggplot(IMGT_Summary, aes(x=as.numeric(Analysed_sequence_length))) + geom_histogram(binwidth=2, color="black", fill="white") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("Analysed Sequence Length")+geom_vline(xintercept=mean(IMGT_Summary$Analysed_sequence_length[!is.na(IMGT_Summary$Analysed_sequence_length)]), color="red") 
	plot(plot_grid(p4, p4.1, ncol=1))
	p4 <- ggplot(IMGT_Summary, aes(x=Sample, y=as.numeric(V_REGION_identity_percent))) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") +geom_hline(yintercept=mean(IMGT_Summary$V_REGION_identity_percent[!is.na(IMGT_Summary$V_REGION_identity_percent)]), color="red") +ylab("V REGION identity %")
	p4.1 <- ggplot(IMGT_Summary, aes(x=as.numeric(V_REGION_identity_percent))) + geom_histogram(binwidth=1, color="black", fill="white") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("V REGION identity %")+geom_vline(xintercept=mean(IMGT_Summary$V_REGION_identity_percent[!is.na(IMGT_Summary$V_REGION_identity_percent)]), color="red") 
	plot(plot_grid(p4, p4.1, ncol=1))
	p4 <- ggplot(IMGT_Summary, aes(x=Sample, y=as.numeric(J_REGION_identity_percent))) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") +geom_hline(yintercept=mean(IMGT_Summary$J_REGION_identity_percent[!is.na(IMGT_Summary$J_REGION_identity_percent)]), color="red") +ylab("J REGION identity %")
	p4.1 <- ggplot(IMGT_Summary, aes(x=as.numeric(J_REGION_identity_percent))) + geom_histogram(binwidth=1, color="black", fill="white") +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("J REGION identity %")+geom_vline(xintercept=mean(IMGT_Summary$J_REGION_identity_percent[!is.na(IMGT_Summary$J_REGION_identity_percent)]), color="red") 
	plot(plot_grid(p4, p4.1, ncol=1))
	dev.off()
	
	V_genes <- data.frame(table(IMGT_Summary$Sample, IMGT_Summary$V_Gene_selected))
	colnames(V_genes) <- c("Sample", "IGV_V", "Freq_V")
	V_genes_family <- data.frame(table(IMGT_Summary$Sample, IMGT_Summary$V_Gene_Family)) 
	colnames(V_genes_family) <- c("Sample", "IGV_V_fam", "Freq_V_fam")
	J_genes <- data.frame(table(IMGT_Summary$Sample, IMGT_Summary$J_Gene_selected)) 
	colnames(J_genes) <- c("Sample", "IGV_J", "Freq_J")
	D_genes <- data.frame(table(IMGT_Summary$Sample, IMGT_Summary$D_Gene_selected)) 
	colnames(D_genes) <- c("Sample", "IGV_D", "Freq_D")
	
	summary_usage <- merge(V_genes, V_genes_family, by="Sample")
	summary_usage <- merge(summary_usage, J_genes, by="Sample")
	summary_usage <- merge(summary_usage, D_genes, by="Sample")
	
	pdf(paste0(path_to_outputdir, '/Plots/IMGT/Gene_Usage_QC_', Run_name, '.pdf'), width=23, height=14)
	p1 <- ggplot(V_genes, aes(x=Sample, y=Freq_V, fill=IGV_V)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="V gene")
	p2 <- ggplot(V_genes_family, aes(x=Sample, y=Freq_V_fam, fill=IGV_V_fam)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="V family")
	p3 <- ggplot(J_genes, aes(x=Sample, y=Freq_J, fill=IGV_J)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="J gene")
	p4 <- ggplot(D_genes, aes(x=Sample, y=Freq_D, fill=IGV_D)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="D gene")
	plot(plot_grid(p1, ncol=1))
	plot(plot_grid(p2, ncol=1))
	plot(plot_grid(p3, ncol=1))
	plot(plot_grid(p4, ncol=1))
	dev.off()

	IMGT_counts <- IMGT_Summary[, c("No_unique_tcrs", "No_non_unique_tcr", "Sample")]
	IMGT_counts <- unique(IMGT_counts)
	pdf(paste0(path_to_outputdir, '/Plots/IMGT/Counts_QC_', Run_name, '.pdf'), width=23, height=14)
	p1 <- ggplot(IMGT_counts, aes(x=Sample, y=as.numeric(No_unique_tcrs))) +geom_bar(stat="identity", fill="white", color="black")  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Number of Unique TCRs (IMGT))") +geom_hline(yintercept=1000, color="blue") +geom_hline(yintercept=(mean(IMGT_counts$No_unique_tcrs) + sd(IMGT_counts$No_unique_tcrs)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(IMGT_counts$No_unique_tcrs) - sd(IMGT_counts$No_unique_tcrs)), col="purple", linetype='dotted') +geom_hline(yintercept=mean(IMGT_counts$No_unique_tcrs), col="red", linetype='dotted')
	p2 <- ggplot(IMGT_counts, aes(x=Sample, y=as.numeric(No_non_unique_tcr))) +geom_bar(stat="identity", fill="white", color="black")  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Number of Unique UMIs (IMGT))") +geom_hline(yintercept=1000, color="blue") +geom_hline(yintercept=(mean(IMGT_counts$No_non_unique_tcr) + sd(IMGT_counts$No_non_unique_tcr)), col="purple", linetype='dotted')+geom_hline(yintercept=(mean(IMGT_counts$No_non_unique_tcr) - sd(IMGT_counts$No_non_unique_tcr)), col="purple", linetype='dotted') +geom_hline(yintercept=mean(IMGT_counts$No_non_unique_tcr), col="red", linetype='dotted')
	plot(plot_grid(p2, p1, ncol=1))
	dev.off()	
	
	IMGT_Summary$V_PRODUCTIVITY_EDITED  <- NA
	IMGT_Summary$V_PRODUCTIVITY_EDITED[IMGT_Summary$V_DOMAIN_Functionality == "productive (see comment)" | IMGT_Summary$V_DOMAIN_Functionality == "productive" ] <- "Productive" 
	IMGT_Summary$V_PRODUCTIVITY_EDITED[IMGT_Summary$V_DOMAIN_Functionality == "unproductive (see comment)" | IMGT_Summary$V_DOMAIN_Functionality == "unproductive" ] <- "Unproductive" 
	IMGT_Summary$V_PRODUCTIVITY_EDITED[IMGT_Summary$V_DOMAIN_Functionality == "No rearrangement found" | IMGT_Summary$V_DOMAIN_Functionality == "No results" ] <- "No Rearrangement/Result" 
	IMGT_Summary$V_PRODUCTIVITY_EDITED[IMGT_Summary$V_DOMAIN_Functionality == "rearranged sequence (but no junction found) (see comment)" ] <- "Rearrangement No Junction" 
	
	IMGT_functionality <- IMGT_Summary[, c("Sample", "V_PRODUCTIVITY_EDITED")]
	table_1 <- data.frame(table(IMGT_functionality))
	
	pdf(paste0(path_to_outputdir, '/Plots/IMGT/Functionality_QC_', Run_name, '.pdf'), width=23, height=14)
    p1 <- ggplot(IMGT_Summary, aes(x=Sample, y=Analysed_sequence_length, color=V_PRODUCTIVITY_EDITED)) +geom_boxplot()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ ylab("Analysed Sequence Length")
	p2 <- ggplot(IMGT_Summary, aes(y=Analysed_sequence_length, x=V_PRODUCTIVITY_EDITED, color=V_Gene_Family)) +geom_boxplot()  +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("V Gene Functionality") + ylab("Analysed Sequence Length")
	p3 <- ggplot(table_1, aes(x=Sample, y=Freq, fill=V_PRODUCTIVITY_EDITED)) +geom_col(position = "fill", colour = "black") + scale_y_continuous(labels = scales::percent) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("% Reads") +labs(fill="Functionality")
	plot(plot_grid(p1, p2, p3,  ncol=1))
	dev.off()
	
	if(dir.exists(paste0(path_to_outputdir, "/Summary/IMGT"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary/IMGT"))
	}				
	write.table(IMGT_Summary, paste0(path_to_outputdir, "/Summary/IMGT/IMGT_SUMMARY_", Run_name, ".txt"), sep="\t")
	write.table(summary_usage, paste0(path_to_outputdir, "/Summary/IMGT/IMGT_SUMMARY_GENE_USAGE_", Run_name, ".txt"), sep="\t")
	return(IMGT_Summary)
}


