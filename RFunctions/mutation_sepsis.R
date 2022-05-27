#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

imgt_mutation_statistics_sepsis <- function(path_to_outputdir = path_to_outputdir, cluster_nodes = 8, productivity="ALL", path_to_layout=path_to_layout){
		library(tidyverse)
	library(ggplot2)
	library(foreach)
	library(doParallel)
	library(gridExtra)
	library(purrr)
	library(seqinr)
	library(matrixStats)
    library(cowplot)
    library(gtools)
	library(ggpubr)

	path <- path_to_outputdir
	path <- paste0(path_to_outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT")
	files <- list.files(path, full.name=TRUE)
	files <- grep('nt-mutation', files, value=TRUE)
	
	## Overall summary file 
	## We are only using this to extract the column names !!
	path2 <- paste0(path_to_outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW")
	files2 <- list.dirs(path = path2, full.names = TRUE, recursive = TRUE)
	path2 <- files2[2]
	files2 <- list.files(path2, full.name=TRUE)
	files2 <- grep('nt-mutation', files2, value=TRUE)
	
	data_overall <- read.delim(files2, header=TRUE, sep='\t')
	data_overall$X <- NULL
	names_data <- colnames(data_overall) 
	
	## Extracting the actual mutation counts 
	if(productivity=="ALL" | productivity=="all"){
			files_not <- grep('_productive_', files, value=TRUE)
			files_not2 <- grep('_unproductive', files, value=TRUE)
			files <- files[!files %in% c(files_not, files_not2)]
		} 
	if(!is.na(productivity)){	
		if(productivity=="productive" | productivity=="PRODUCTIVE" ){	
				files <- grep('_productive_', files, value=TRUE)
			} 
			
		if(productivity=="unproductive"| productivity=="UNPRODUCTIVE"){	
				files <- grep('_unproductive', files, value=TRUE)
			}
	}
	
	## Reading in Mutation Statistics 
	cl <- cluster_nodes
	registerDoParallel(cl)
	IMGT_Mutation <- foreach(i = 1:length(files), .combine=rbind, .packages='tidyverse') %dopar% {
				filtered_path <- files[i]
				print(i)
				sample_id <- unlist(str_split(filtered_path, 'IMGT_SPLIT/'))[2]
				sample_id <- unlist(str_split(sample_id, '_8_V-REGION-nt-mutation-statistics.txt'))[1]
				sample_id <- gsub("IMGT_BCR_", "", sample_id)
				output <- read.delim(filtered_path, sep="\t", header=FALSE, fill=TRUE, col.names=names_data)
				colnames(output) <- names_data
				output$Sample <- sample_id
								
				`%notin%` <- Negate(`%in%`)
				cols_notcalc <- c("Sequence.number", "Sequence.ID", "V.DOMAIN.Functionality", "V.GENE.and.allele", "Sample")
				cols_calc <- names_data[names_data %notin% cols_notcalc]
				output_1 <- data.frame(apply(output, 2, function(y) gsub("\\s\\(\\d+\\)", "", y)))
				
				output_1[,cols_calc] <- sapply(output_1[,cols_calc],as.character)
				output_1[,cols_calc] <- sapply(output_1[,cols_calc],as.numeric)
				
				## Excludes NA valuses 
				return(output_1)	    
	}
				
	if(dir.exists(paste0(path_to_outputdir, "Plots/IMGT"))==FALSE){
		dir.create(paste0(path_to_outputdir, "Plots/IMGT"))
	}
	
	FullData <- IMGT_Mutation
	FullData <- FullData[!FullData$V.DOMAIN.Functionality=="rearranged sequence (but no junction found)",]
	FullData <- FullData[!FullData$V.DOMAIN.Functionality=="No rearrangement found",]

	# Columns for analysis 
	cols_notcalc <- c("Sequence.number", "Sequence.ID", "V.DOMAIN.Functionality", "V.GENE.and.allele", "Sample")
	cols_calc <- names_data[!names_data %in% cols_notcalc]
	
	
	# Reading in layouts
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
	
	FullData <- merge(FullData, layouts, by.x="Sample", by.y="SampleID")
	
	FullData$Position <- as.character(FullData$Position)
	FullData$PCRBarcode <- as.character(FullData$PCRBarcode)
	FullData$Library <- as.character(FullData$Library)
	FullData$Plate <- as.character(FullData$Plate)
	FullData$Lane <- as.character(FullData$Lane)
	
	## if Day is provided 
	
	if(any(!layouts$Barcode==layouts$SampleID)){	
		days <- data.frame(str_split_fixed(FullData$Sample, "_", 2))
		days <- days$X2
		FullData <- cbind(FullData, days)
		FullData$days <- gsub("T", "", FullData$days)
		FullData$days[FullData$days != 1 & FullData$days != 3 & FullData$days != 5] <- "Control"
	} else {
		FullData$days <- "NA"
	} 
	
	## Generic for all data sets. 
	
	## Get Relevant plotting columns 
	mutations <- grep('mutations', cols_calc, value=TRUE)
	transitions <- cols_calc[!cols_calc %in% mutations]
	mutation_table <- FullData[,c(mutations)]
	ignore <- grep("Nb", transitions, value=TRUE)
	transitions <- transitions[!transitions %in% ignore]
	
	# Meta Data File
	meta_data <- read.delim('/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Meta_data_for_cohort_1.txt', sep='\t', header=TRUE)
	
	## Infer SRS
	meta_data$SRS_INFERRED <- NA 	
	length(meta_data$SampleID)
	for(i in 1:length(meta_data$SampleID)){
	    sample <- meta_data$alternative_barcode[i]
		srs <- unique(meta_data$SRS_New[meta_data$alternative_barcode==sample])
		srs <- srs[!is.na(srs)]
		if(length(srs)==1){
			meta_data$SRS_INFERRED[meta_data$alternative_barcode==sample] <- srs
		}
		if(length(srs)<1){
			print("No assignment Found")
			print(i)
			meta_data$SRS_INFERRED[meta_data$alternative_barcode==sample] <- NA
		}
		if(length(srs)>1){
			print(paste0("Warning srs movement"))
			meta_data$SRS_INFERRED[meta_data$alternative_barcode==sample] <- "SRS Movement" 
		}
	}


	data_mutation <- merge(FullData, meta_data, by.x="Sample", by.y="SampleID", all.x=TRUE)
	data_mutation$RNA_time <- as.character(data_mutation$RNA_time)
	data_mutation$RNA_time[is.na(data_mutation$RNA_time)] <- "Control"
	
	data_mutation$outcome <- NA 
	data_mutation$outcome[is.na(data_mutation$Days_death_from_ICU)] <- "CONTROL"
	data_mutation$outcome[data_mutation$Days_death_from_ICU=="Alive"] <- "ALIVE"
	data_mutation$outcome[data_mutation$Days_death_from_ICU!="Alive"] <- "DEAD"
	data_mutation$outcome[data_mutation$RNA_time=="Control"] <- "CONTROL"
	
	data_mutation$Sex <- as.character(data_mutation$Sex)
	data_mutation$Sex[is.na(data_mutation$Sex)] <- "Unknown"	
	data_mutation$Sex <- as.factor(data_mutation$Sex)
	
	data_mutation$Shock_sepsis2 <- as.character(data_mutation$Shock_sepsis2)
	data_mutation$Shock_sepsis2[is.na(data_mutation$Shock_sepsis2)] <- "Unknown"	
	data_mutation$Shock_sepsis2 <- as.factor(data_mutation$Shock_sepsis2)
			
	data_mutation$SRS_INFERRED[data_mutation$RNA_time=="Control"] <- "CONTROL"
	
	##trying to add day 
	mutations <- grep('mutations', cols_calc, value=TRUE)
	transitions <- cols_calc[!cols_calc %in% mutations]
	mutation_table <- IMGT_Mutation[,c(mutations)]
	ignore <- grep("Nb", transitions, value=TRUE)
	transitions <- transitions[!transitions %in% ignore]
	
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation.pdf'), width=23, height=20)
	for(i in 1:length(mutations)){
		column_id <- mutations[i]
		column_id_neat <- gsub("[[:punct:]]", " ", column_id)
		
		#by outcome
		subseted <- data_mutation[, c('RNA_time', column_id, 'outcome')]
		colnames(subseted) <- c("RNA_time", "X", 'outcome')
		mu <-  plyr::ddply(subseted, c("RNA_time", 'outcome'), summarise, grp.mean=mean(X, na.rm=TRUE))
		subseted <- data_mutation[, c('RNA_time', column_id)]
		colnames(subseted) <- c("RNA_time", "X")
		mub <-  plyr::ddply(subseted, c("RNA_time"), summarise, grp.mean=mean(X, na.rm=TRUE))
		by outcome
		subseted <- data_mutation[, c('RNA_time', column_id, 'outcome')]
		colnames(subseted) <- c("RNA_time", "X", 'outcome')
		mu1 <-  plyr::ddply(subseted, c("RNA_time", 'outcome'), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by SRS
		subseted <- data_mutation[, c('RNA_time', column_id, 'SRS')]
		colnames(subseted) <- c("RNA_time", "X", 'SRS')
		mu2 <-  plyr::ddply(subseted, c("RNA_time", 'SRS'), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by SRS_Inferred
		subseted <- data_mutation[, c('RNA_time', column_id, 'SRS_INFERRED')]
		colnames(subseted) <- c("RNA_time", "X", 'SRS_INFERRED')
		mu3 <-  plyr::ddply(subseted, c("RNA_time", 'SRS_INFERRED'), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by Shock
		subseted <- data_mutation[, c('RNA_time', column_id, 'Shock_sepsis2')]
		colnames(subseted) <- c("RNA_time", "X", 'Shock_sepsis2')
		mu4 <-  plyr::ddply(subseted, c("RNA_time", 'Shock_sepsis2'), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by Sex
		subseted <- data_mutation[, c('RNA_time', column_id, 'Sex')]
		colnames(subseted) <- c("RNA_time", "X", 'Sex')
		mu5 <-  plyr::ddply(subseted, c("RNA_time", 'Sex'), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by EBV_positive_metagenomics
		subseted <- data_mutation[, c('RNA_time', column_id, 'EBV_positive_metagenomics')]
		colnames(subseted) <- c("RNA_time", "X", 'EBV_positive_metagenomics')
		mu6 <-  plyr::ddply(subseted, c("RNA_time", 'EBV_positive_metagenomics'), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by EBV_positive_ddpcr
		subseted <- data_mutation[, c('RNA_time', column_id, 'EBV_positive_ddpcr')]
		colnames(subseted) <- c("RNA_time", "X", 'EBV_positive_ddpcr')
		mu7 <-  plyr::ddply(subseted, c("RNA_time", 'EBV_positive_ddpcr'), summarise, grp.mean=mean(X, na.rm=TRUE))
				
		p1 <- ggplot(data_mutation, aes(x=Sample, y=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.7, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") + ylab(column_id_neat) +xlab("Sample") 
		p2 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~outcome) + geom_vline(data=mu, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p3 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=outcome)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~RNA_time) + geom_vline(data=mu1, aes(xintercept=grp.mean, color=outcome),linetype="solid", size=0.8)
		p4 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~SRS) + ggtitle("SRS")+ geom_vline(data=mu2, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p4.1 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~SRS_INFERRED) + ggtitle("SRS Inferred")+ geom_vline(data=mu3, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p5 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~Shock_sepsis2)  + ggtitle("Septic Shock")+ geom_vline(data=mu4, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p6<- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density(alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~Sex) + ggtitle("Sex")+ geom_vline(data=mu5, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p7 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density(alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~EBV_positive_metagenomics) + ggtitle("EBV Metgagenomics")+ geom_vline(data=mu6, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p8 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~EBV_positive_ddpcr) + ggtitle("EBV ddPCR")+ geom_vline(data=mu7, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p9 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() + geom_vline(data=mub, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		plot(plot_grid(p1, p2, p3, ncol=1))
		plot(plot_grid(p4, p4.1, p5, ncol=1))
		plot(plot_grid(p6, p7, p8, ncol=1))
		plot(plot_grid(p2, p9, ncol=1))
	}
	dev.off()
	
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_transitions.pdf'), width=23, height=20)
	for(i in 1:length(transitions)){
		column_id <- transitions[i]
		column_id_neat <- gsub("[[:punct:]]", " ", column_id)
		#by outcome
		subseted <- data_mutation[, c('RNA_time', column_id, 'outcome')]
		colnames(subseted) <- c("RNA_time", "X", 'outcome')
		mu <-  plyr::ddply(subseted, c("RNA_time", 'outcome'), summarise, grp.mean=mean(X, na.rm=TRUE))
		subseted <- data_mutation[, c('RNA_time', column_id)]
		colnames(subseted) <- c("RNA_time", "X")
		mub <-  plyr::ddply(subseted, c("RNA_time"), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by outcome
		subseted <- data_mutation[, c('RNA_time', column_id, 'outcome')]
		colnames(subseted) <- c("RNA_time", "X", 'outcome')
		mu1 <-  plyr::ddply(subseted, c("RNA_time", 'outcome'), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by SRS
		subseted <- data_mutation[, c('RNA_time', column_id, 'SRS')]
		colnames(subseted) <- c("RNA_time", "X", 'SRS')
		mu2 <-  plyr::ddply(subseted, c("RNA_time", 'SRS'), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by SRS_Inferred
		subseted <- data_mutation[, c('RNA_time', column_id, 'SRS_INFERRED')]
		colnames(subseted) <- c("RNA_time", "X", 'SRS_INFERRED')
		mu3 <-  plyr::ddply(subseted, c("RNA_time", 'SRS_INFERRED'), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by Shock
		subseted <- data_mutation[, c('RNA_time', column_id, 'Shock_sepsis2')]
		colnames(subseted) <- c("RNA_time", "X", 'Shock_sepsis2')
		mu4 <-  plyr::ddply(subseted, c("RNA_time", 'Shock_sepsis2'), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by Sex
		subseted <- data_mutation[, c('RNA_time', column_id, 'Sex')]
		colnames(subseted) <- c("RNA_time", "X", 'Sex')
		mu5 <-  plyr::ddply(subseted, c("RNA_time", 'Sex'), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by EBV_positive_metagenomics
		subseted <- data_mutation[, c('RNA_time', column_id, 'EBV_positive_metagenomics')]
		colnames(subseted) <- c("RNA_time", "X", 'EBV_positive_metagenomics')
		mu6 <-  plyr::ddply(subseted, c("RNA_time", 'EBV_positive_metagenomics'), summarise, grp.mean=mean(X, na.rm=TRUE))
		#by EBV_positive_ddpcr
		subseted <- data_mutation[, c('RNA_time', column_id, 'EBV_positive_ddpcr')]
		colnames(subseted) <- c("RNA_time", "X", 'EBV_positive_ddpcr')
		mu7 <-  plyr::ddply(subseted, c("RNA_time", 'EBV_positive_ddpcr'), summarise, grp.mean=mean(X, na.rm=TRUE))
				
		p1 <- ggplot(data_mutation, aes(x=Sample, y=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.7, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue") + ylab(column_id_neat) +xlab("Sample") 
		p2 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density(alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~outcome) + geom_vline(data=mu, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p3 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=outcome)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~RNA_time) + geom_vline(data=mu1, aes(xintercept=grp.mean, color=outcome),linetype="solid", size=0.8)
		p4 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~SRS) + ggtitle("SRS")+ geom_vline(data=mu2, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p4.1 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density(alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~SRS_INFERRED) + ggtitle("SRS Inferred")+ geom_vline(data=mu3, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p5 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~Shock_sepsis2)  + ggtitle("Septic Shock")+ geom_vline(data=mu4, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p6<- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density(alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~Sex) + ggtitle("Sex")+ geom_vline(data=mu5, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p7 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~EBV_positive_metagenomics) + ggtitle("EBV Metgagenomics")+ geom_vline(data=mu6, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p8 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() +facet_wrap(~EBV_positive_ddpcr) + ggtitle("EBV ddPCR")+ geom_vline(data=mu7, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		p9 <- ggplot(data_mutation, aes(x=IMGT_Mutation[,column_id], fill=RNA_time)) + geom_density( alpha = 0.7) +xlab(column_id_neat) + theme_classic() + geom_vline(data=mub, aes(xintercept=grp.mean, color=RNA_time),linetype="solid", size=0.8)
		plot(plot_grid(p1, p2, p3, ncol=1))
		plot(plot_grid(p4, p4.1, p5, ncol=1))
		plot(plot_grid(p6, p7, p8, ncol=1))
		plot(plot_grid(p2, p9, ncol=1))
	}
	dev.off()
	
	# Convert to long format for comparing 
	data_long <- gather(data_mutation, Region, Nb.Mutations, all_of(mutations))
	data_long$RNA_time <- as.factor(data_long$RNA_time)
	
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation_region.pdf'), width=20, height=20)
	p1 <- ggplot(data_long, aes(x=Region, y=Nb.Mutations, fill=RNA_time)) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, aes(fill=RNA_time)) + ylab("Number of Mutation") +xlab("V Region") 
	p2 <- ggplot(data_long, aes(x=Nb.Mutations, fill=RNA_time)) + geom_density(, alpha = 0.7) +xlab("Number of Mutations") + theme_classic() +facet_wrap(~Region, scales="free") 
	p3 <- ggplot(data_long, aes(x=V.DOMAIN.Functionality, y=Nb.Mutations, fill=RNA_time)) + geom_violin() +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, aes(fill=RNA_time)) + ylab("Number of Mutation") +xlab("V Region") 
	p4 <- ggplot(data_long, aes(x=Nb.Mutations, fill=V.DOMAIN.Functionality)) + geom_density(, alpha = 0.7) +xlab("Number of Mutations") + theme_classic() +facet_wrap(~Region, scales="free") 
	plot(plot_grid(p1, ncol=1))
	plot(plot_grid(p2, ncol=1))
	plot(plot_grid(p3, ncol=1))
	plot(plot_grid(p4, ncol=1))
	dev.off()
	
	if(dir.exists(paste0(path_to_outputdir, "/Summary/IMGT"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary/IMGT"))
	}				
	write.table(IMGT_Mutation, paste0(path_to_outputdir, "/Summary/IMGT/IMGT_SUMMARY_Mutation.txt"), sep="\t")
	
	## Calculating Read depth 
	path <- paste0(path_to_outputdir, "ORIENTATED_SEQUENCES/NETWORKS/")
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("Att", samples, value=TRUE)
	mins <- c()
	order_samples <- c()
	for(i in 1:length(samples)){
		a <- read.delim(samples[i], header=FALSE)
		a <- a$V2
		a_min <- length(a)
		sampleid <- samples[i]
		order_samples <- c(order_samples, sampleid)
		mins <- c(mins, a_min)
		} 
				
	depths <- data.frame(cbind(order_samples, mins))
	depths$order_samples <- gsub(paste0(path, "/Att_BCR_"), "", depths$order_samples)
	depths$order_samples <- gsub(paste0(path, "/Att_TCR_"), "", depths$order_samples)
	depths$order_samples <- gsub(paste0(path, "/Att_TRB_"), "", depths$order_samples)
	depths$order_samples <- gsub(paste0(path, "/Att_TRA_"), "", depths$order_samples)
	depths$order_samples <- gsub(".txt", "", depths$order_samples)	
	colnames(depths) <- c("SampleIDforDepths", "ReadDepth")
	depths$ReadDepth <- as.character(depths$ReadDepth)
	depths$ReadDepth <- as.numeric(depths$ReadDepth)

	Means_mutation <- foreach(i = 1:length(mutations), .combine=cbind, .packages='tidyverse') %dopar% {
		column_id <- mutations[i]
		column_id_neat <- gsub("[[:punct:]]", " ", column_id)
		subseted <- data_mutation[, c('Sample', column_id)]
		colnames(subseted) <- c("Sample", "X")
		mu1 <-  plyr::ddply(subseted, c("Sample"), summarise, grp.mean=mean(X, na.rm=TRUE))
		colnames(mu1) <- c("Sample", paste0("Mean.", column_id))
		rownames(mu1) <- mu1$Sample
		mu1$Sample <- NULL
		return(mu1)
		}
	cols_to_plot <- colnames(Means_mutation)
	Means_mutation$Sample <- rownames(Means_mutation)	
	Means_mutation <- merge(Means_mutation, depths, by.x="Sample", by.y="SampleIDforDepths")
	
	Means_mutation$Region <- NA
	data_longx <- gather(Means_mutation, RegionFull, Mean.Nb.Mutations, all_of(cols_to_plot))
	
	v <- grep("V.REGION", cols_to_plot, value=TRUE)
	fr1 <- grep("FR1", cols_to_plot, value=TRUE)
	fr2 <- grep("FR2", cols_to_plot, value=TRUE)
	fr3 <- grep("FR3", cols_to_plot, value=TRUE)
	CDR1 <- grep("CDR1", cols_to_plot, value=TRUE)
	CDR2 <- grep("CDR2", cols_to_plot, value=TRUE)
	CDR3 <- grep("CDR3", cols_to_plot, value=TRUE)
	
	silent <- grep("of.silent.", cols_to_plot, value=TRUE)
	nonsilent <- grep("nonsilent", cols_to_plot, value=TRUE)
	total <- grep("of.mutations", cols_to_plot, value=TRUE)
	
	
	data_longx$Region <- NA
	data_longx$Region[data_longx$RegionFull %in% v] <- "V Region"
	data_longx$Region[data_longx$RegionFull %in% fr1] <- "FR1"
	data_longx$Region[data_longx$RegionFull %in% fr2] <- "FR2"
	data_longx$Region[data_longx$RegionFull %in% fr3] <- "FR3"
	data_longx$Region[data_longx$RegionFull %in% CDR1] <- "CDR1"
	data_longx$Region[data_longx$RegionFull %in% CDR2] <- "CDR2"
	data_longx$Region[data_longx$RegionFull %in% CDR3] <- "CDR3"
	
	data_longx$Type <- NA
	data_longx$Type[data_longx$RegionFull %in% silent] <- "silent"
	data_longx$Type[data_longx$RegionFull %in% nonsilent] <- "nonsilent"
	data_longx$Type[data_longx$RegionFull %in% total] <- "total"
	
	# Plotting Each mutation as a function of read depth 
	## Look a lot are significant!!!
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation_Readdepth.pdf'), width=5, height=5)
	for(s in cols_to_plot){
				print(s)
				h <-  gsub("[[:punct:]]", " ", s)
				metric1 <- h
				plots <- list()
				d <- ggplot(Means_mutation[Means_mutation[, s] >-1,], aes_string(x="ReadDepth", y=s, colour="ReadDepth")) + geom_point() + theme_bw() +xlab("ReadDepth") +ylab(metric1) + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) + stat_summary(fun.data= mean_cl_normal) + geom_smooth(method='lm') + stat_cor(method = "spearman", label.x = 70000, label.y = 0.01)
				plot(grid.arrange(d, ncol=1))
	} 
	dev.off()
	
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation_Readdepth_Summary.pdf'), width=20, height=20)
	d <- ggplot(data_longx, aes(y=Mean.Nb.Mutations, x=ReadDepth, colour=Type)) + geom_point() + theme_bw() +xlab("ReadDepth") +ylab("Mean Number of Mutations") + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) + stat_summary(fun.data= mean_cl_normal) + geom_smooth(method='lm') + stat_cor(method = "spearman", label.x = 70000) +facet_wrap(~Region, scales="free")
	plot(d)	
	dev.off()			
	
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

