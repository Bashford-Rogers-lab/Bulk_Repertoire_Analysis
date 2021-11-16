# Author: Lauren Overend 
# Function to visualise the Filtering Reports for all files processed through the RBR Bulk BCR/TCR pipeline
# Author: Lauren Overend 



imgt_mutation_statistics <- function(path_to_outputdir = path_to_outputdir, cluster_nodes = 8, productivity="ALL", path_to_layout=path_to_layout){
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
	library(ggrastr)
	
	`%notin%` <- Negate(`%in%`)
	
	## Makes files considerabley smaller!
	pdf.options(useDingbats = TRUE) 
	
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
	data_overall <- NULL
	
	
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
	
	# read in and bind data 
	filenames <- files
	#create a list of dataframes 
	df_list <- lapply(filenames, fread, header = FALSE, sep="\t", fill=TRUE, col.names=names_data, colClasses='character')
	# identify the sample they came from 
	d <- str_split(filenames, 'IMGT_SPLIT/') 
	d <- sapply(d, "[[", 2)  
	d <- gsub("_8_V-REGION-nt-mutation-statistics.txt", "", d)
	d <- gsub("IMGT_BCR_", "", d)
	d <- gsub("BCR_", "", d)
  
	# Name each dataframe with the run and filename
	names(df_list) <- d

	# Create combined dataframe  
	df <- df_list %>% bind_rows(.id = 'Sample')
	IMGT_Mutation <- data.frame(df)
	IMGT_Mutation <- data.frame(apply(IMGT_Mutation, 2, function(y) gsub("\\s\\(\\d+\\)", "", y)))
	# Assign dataframe to the name of the pattern  
	
	cols_notcalc <- c("Sequence.number", "Sequence.ID", "V.DOMAIN.Functionality", "V.GENE.and.allele", "Sample")
	cols_calc <- colnames(IMGT_Mutation)[colnames(IMGT_Mutation) %notin% cols_notcalc]
	IMGT_Mutation[,cols_calc] <- sapply(IMGT_Mutation[,cols_calc],as.character)
	IMGT_Mutation[,cols_calc] <- sapply(IMGT_Mutation[,cols_calc],as.numeric)
	
	## Save this!!!
	write.table(IMGT_Mutation, paste0(path_to_outputdir,'Summary/All_Mutations_', productivity, '.txt'), sep="\t")

##########################	
	if(dir.exists(paste0(path_to_outputdir, "Plots/IMGT"))==FALSE){
		dir.create(paste0(path_to_outputdir, "Plots/IMGT"))
	}
	
	FullData <- IMGT_Mutation
	FullData <- FullData[!FullData$V.DOMAIN.Functionality=="rearranged sequence (but no junction found)",]
	FullData <- FullData[!FullData$V.DOMAIN.Functionality=="No rearrangement found",]
	FullData$Sample <- gsub("_unproductive", "", FullData$Sample)
    FullData$Sample <- gsub("_productive", "", FullData$Sample) 
	
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
	
	FullData$Sample <- gsub("IMGT_", "", FullData$Sample)
	FullData$Sample <- gsub("TCR_", "", FullData$Sample)
	FullData$Sample <- gsub("BCR_", "", FullData$Sample)
	layouts$SampleID <- gsub("IMGT_", "", layouts$SampleID)
	layouts$Barcode <- gsub("IMGT_", "", layouts$Barcode)
	
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
	
	## Get Relevant plotting columns 
	mutations <- grep('mutations', cols_calc, value=TRUE)
	transitions <- cols_calc[!cols_calc %in% mutations]
	mutation_table <- FullData[,c(mutations)]
	ignore <- grep("Nb", transitions, value=TRUE)
	transitions <- transitions[!transitions %in% ignore]
	
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation_', productivity, '.pdf'), width=60, height=20, useDingbats = TRUE)
	for(i in 1:length(mutations)){
		column_id <- mutations[i]
		column_id_neat <- gsub("[[:punct:]]", " ", column_id)		
		p1 <- ggplot(FullData, aes(x=Sample, y=FullData[,column_id], fill=days)) + rasterise(geom_violin(), dpi = 300) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.7, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue", position = position_dodge(width = .75)) + ylab(column_id_neat) +xlab("Sample") + ggtitle(productivity) +labs(fill="Timepoint")
		p2 <- ggplot(FullData, aes(x=FullData[,column_id], fill=days)) + rasterise(geom_density( alpha = 0.7), dpi = 300) +xlab(column_id_neat) + theme_classic()+ ggtitle(productivity) +labs(fill="Timepoint")
		plot(plot_grid(p1, p2, ncol=1))
	}
	dev.off()
	
	# Convert to long format for comparing 
	data_long <- gather(FullData, Region, Nb.Mutations, all_of(mutations))
	data_long$days <- as.factor(data_long$days)
	
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation_region_', productivity, '.pdf'), width=20, height=15)
	p1 <- ggplot(data_long, aes(x=Region, y=Nb.Mutations, fill=days)) + rasterise(geom_boxplot(), dpi = 300) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, aes(fill=days),  position = position_dodge(width = .75)) + ylab("Number of Mutations") +xlab("V Region") + ggtitle(productivity)+labs(fill="Timepoint")
	p2 <- ggplot(data_long, aes(x=Nb.Mutations, fill=days)) + rasterise(geom_density(, alpha = 0.7), dpi = 300) +xlab("Number of Mutations") + theme_classic() +facet_wrap(~Region, scales="free")+ ggtitle(productivity) +labs(fill="Timepoint")
	p3 <- ggplot(data_long, aes(x=V.DOMAIN.Functionality, y=Nb.Mutations, fill=days)) + rasterise(geom_boxplot(), dpi = 300) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, aes(fill=days),  position = position_dodge(width = .75)) + ylab("Number of Mutations") +xlab("V Region") + ggtitle(productivity) +labs(fill="Timepoint")
	p4 <- ggplot(data_long, aes(x=Nb.Mutations, fill=V.DOMAIN.Functionality)) + rasterise(geom_density(, alpha = 0.7), dpi = 300) +xlab("Number of Mutations") + theme_classic() +facet_wrap(~Region, scales="free") + ggtitle(productivity)
	plot(plot_grid(p1, ncol=1))
	plot(plot_grid(p2, ncol=1))
	plot(plot_grid(p3, ncol=1))
	plot(plot_grid(p4, ncol=1))
	dev.off()
	
	if(dir.exists(paste0(path_to_outputdir, "/Summary/IMGT"))==FALSE){
		dir.create(paste0(path_to_outputdir, "/Summary/IMGT"))
	}				
	write.table(IMGT_Mutation, paste0(path_to_outputdir, "/Summary/IMGT/IMGT_SUMMARY_Mutation_", productivity, ".txt"), sep="\t")
	
	# Calculating average mutations 
	Means_mutation <- foreach(i = 1:length(mutations), .combine=cbind, .packages='tidyverse') %dopar% {
		column_id <- mutations[i]
		column_id_neat <- gsub("[[:punct:]]", " ", column_id)
		subseted <- FullData[, c('Sample', column_id)]
		colnames(subseted) <- c("Sample", "X")
		mu1 <-  plyr::ddply(subseted, c("Sample"), summarise, grp.mean=mean(X, na.rm=TRUE))
		colnames(mu1) <- c("Sample", paste0("Mean.", column_id))
		rownames(mu1) <- mu1$Sample
		mu1$Sample <- NULL
		return(mu1)
		}
		
	# Merge Average mutation by read depth 
	cols_to_plot <- colnames(Means_mutation)
	Means_mutation$Sample <- rownames(Means_mutation)	

	data_longx <- gather(Means_mutation, RegionFull, Mean.Nb.Mutations, all_of(cols_to_plot))
	
	# Classifying the region
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
	
	## Identify if there are multiple timepoints 
	if(any(!layouts$Barcode==layouts$SampleID)){	
		days <- data.frame(str_split_fixed(Means_mutation$Sample, "_", 2))
		days <- days$X2
		Means_mutation <- cbind(Means_mutation, days)
		Means_mutation$days <- gsub("T", "", Means_mutation$days)
		Means_mutation$days[Means_mutation$days != 1 & Means_mutation$days != 3 & Means_mutation$days != 5] <- "Control"
	} else {
		Means_mutation$days <- "NA"
	} 
	
	if(any(!layouts$Barcode==layouts$SampleID)){	
		days <- data.frame(str_split_fixed(data_longx$Sample, "_", 2))
		days <- days$X2
		data_longx <- cbind(data_longx, days)
		data_longx$days <- gsub("T", "", data_longx$days)
		data_longx$days[data_longx$days != 1 & data_longx$days != 3 & data_longx$days != 5] <- "Control"
	} else {
		data_longx$days <- "NA"
	} 
	
	## Looking at troughs and peaks 
	data_new <- data_long[data_long$Region=="V.REGION.Nb.of.mutations" | data_long$Region=="V.REGION.Nb.of.nonsilent.mutations" | data_long$Region=="V.REGION.Nb.of.silent.mutations",]
	data_new$V.DOMAIN.Functionality <- as.character(data_new$V.DOMAIN.Functionality)

	data_new$Type[data_new$Region =="V.REGION.Nb.of.mutations" ] <- "TOTAL"
	data_new$Type[data_new$Region== "V.REGION.Nb.of.silent.mutations" ] <- "SILENT"
	data_new$Type[data_new$Region== "V.REGION.Nb.of.nonsilent.mutations"] <- "NON-SILENT"
	
	# Identifying distribution peaks and troughs 
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation_Summary_histo_', productivity, '.pdf'), width=15, height=5)
	p<-ggplot(data_new, aes(x=Nb.Mutations, fill=Type))  + rasterise(geom_bar(aes(y = ..prop..) ), dpi = 300) +facet_wrap(~days+Type)  + geom_vline(xintercept=7,linetype="solid",  color="red") + geom_vline(xintercept=23,linetype="solid",  color="red")+ geom_vline(xintercept=3,linetype="solid", color="blue")+ geom_vline(xintercept=11,linetype="solid", color="blue")+ geom_vline(xintercept=75,linetype="solid", color="blue") +theme_bw() +xlab("V Region: Total Number of Mutations") + ylab("Proportion")+ ggtitle(productivity)
	p1<-ggplot(data_new, aes(x=Nb.Mutations))  + rasterise(geom_bar(aes(y = ..prop..) ), dpi = 300) + geom_vline(xintercept=7,linetype="solid", color="red") + geom_vline(xintercept=23,linetype="solid", color="red")+ geom_vline(xintercept=3,linetype="solid", color="blue")+ geom_vline(xintercept=11,linetype="solid", color="blue")+ geom_vline(xintercept=75,linetype="solid", size=0.8, color="blue")+theme_bw()+xlab("V Region: Total Number of Mutations") + ylab("Proportion") + facet_wrap(~Type)+ ggtitle(productivity)
	plot(p)
	plot(p1)
	dev.off()

	# Asigning class of VDJ sequence
	data_new <- data_new[data_new$Type=="TOTAL",]
	data_new$Nb.Mutations <- as.numeric(data_new$Nb.Mutations)
	data_new$ClassSequence <- NULL
	data_new$ClassSequence[data_new$Nb.Mutations <= 3] <- "UNMUTATED (naive)"
	data_new$ClassSequence[data_new$Nb.Mutations <= 11 & data_new$Nb.Mutations > 3] <- "LOW MUTATION (extrafollicular)"
	data_new$ClassSequence[data_new$Nb.Mutations <= 75 & data_new$Nb.Mutations > 11] <- "HIGH MUTATION (SHM germinal centre)"
	data_new$ClassSequence[data_new$Nb.Mutations > 75] <- "VERY HIGH MUTATION (SHM germinal centre)"

	data_new <- data_new[!is.na(data_new$ClassSequence),]
	
	frequencies_group_t <- data_new[data_new$Type=="TOTAL",] %>%
	group_by(days, ClassSequence) %>%
	summarise(n = n()) %>%
	mutate(freq = n / sum(n))
    frequencies_group_t <- data.frame(frequencies_group_t)


	## Plotting change in proportions across time 
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation_Summary_box_', productivity, '.pdf'), width=8, height=5)
	q1 <- ggplot(frequencies_group_t, aes(fill=ClassSequence, x=days, y=freq), na.rm = TRUE) +  rasterise(geom_col(position = "dodge", na.rm = TRUE), dpi = 300) +theme_bw() + xlab("Timepoint") + ylab("Frequency") + ggtitle(paste0("Total Mutations ", productivity))+labs(fill="Class of VDJ")
	grid.arrange(q1)
	dev.off()
	
	
	w <- with(data_new, table(Sample, ClassSequence))
	w <- data.frame(prop.table(w, margin = 1))
	
	if (length(unique(data_new$Sample)) > 500){
		widthx <- 100 
	} else {
		widthx <- 40
	} 
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation_Summary_box_per_sample', productivity, '.pdf'), width=widthx, height=5)
	q1 <- ggplot(w, aes(fill=ClassSequence, x=Sample, y=Freq), na.rm = TRUE) +  rasterise(geom_col(position = "dodge", na.rm = TRUE), dpi = 300) +theme_bw() + xlab("Timepoint") + ylab("Proportion of Reads") + ggtitle(paste0("Total Mutations ", productivity))+labs(fill="Class of VDJ") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	grid.arrange(q1)
	dev.off()
	
}

