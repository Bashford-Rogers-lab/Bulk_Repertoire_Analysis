# Author: Lauren Overend 
# Function to visualise the Filtering Reports for all files processed through the RBR Bulk BCR/TCR pipeline
# Author: Lauren Overend 
# December 2021
#path_to_outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/BCR_CH2/'
#path_to_layout <- '/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Batching_Layouts_SEPSIS_BCR_CH12.txt'
#productivity <- "UNPRODUCTIVE"
#iso_type <- "UNPRODUCTIVE"

suppressMessages(library(plot3D))

scatter3D_fancy <- function(x, y, z,..., colvar = z)
  {
   panelfirst <- function(pmat) {
      XY <- trans3D(x, y, z = rep(min(z), length(z)), pmat = pmat)
      scatter2D(XY$x, XY$y, colvar = colvar, pch = ".", 
              cex = 2, add = TRUE, colkey = FALSE)
   
      XY <- trans3D(x = rep(min(x), length(x)), y, z, pmat = pmat)
      scatter2D(XY$x, XY$y, colvar = colvar, pch = ".", 
              cex = 2, add = TRUE, colkey = FALSE)
  }
  scatter3D(x, y, z, ..., colvar = colvar, panel.first=panelfirst,
    colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75)) 
}


imgt_mutation_statistics <- function(path_to_outputdir = path_to_outputdir, cluster_nodes = 8, productivity=productivity, path_to_layout=path_to_layout){
	suppressMessages(library(tidyverse))
	suppressMessages(library(data.table))
	suppressMessages(library(ggplot2))
	suppressMessages(library(ggforce))
	suppressMessages(library(Gviz))
	suppressMessages(library(foreach))
	suppressMessages(library(doParallel))
	suppressMessages(library(gridExtra))
	suppressMessages(library(cowplot))
	suppressMessages(library(gtools))
	suppressMessages(library(purrr))
	suppressMessages(library(reshape2))
	suppressMessages(library(Hmisc))
	suppressMessages(library(corrplot))
	suppressMessages(library(stringr))
	suppressMessages(library(dplyr))
	suppressMessages(library(tidyr))
	suppressMessages(library(matrixStats))
	suppressMessages(library(ggpubr))
	suppressMessages(library(ggrastr))
	suppressMessages(library(ggpubr))
	suppressMessages(library(plot3D))
	
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
	
	
	## We need to exclude samples below certain read depth!
	## Use the same filter as we used in the isotyper script 
	depths <- data.frame(table(IMGT_Mutation$Sample))
	depths$Var1 <- as.character(depths$Var1)
	
	## Getting subsample depths which were used for isotyper script 
	counts_used <- paste0(path_to_outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS")
	all_files <- list.files(counts_used, full.name=TRUE)
	all_files <- grep("depth_per_isotype", all_files, value=TRUE)
	counts_used <- read.delim(all_files[1], sep="\t", header=TRUE)
	counts_used <- counts_used[counts_used$type=="UNIQ",]
	subsampled_depth_all <- counts_used$min[counts_used$X.isotype=="all"]
	
	if(productivity == "UNPRODUCTIVE"){ 
	subsampled_depth_allx <- subsampled_depth_all/10
	} else {
		subsampled_depth_allx <- subsampled_depth_all
	} 

	
	# samples exclude 
	to_exclude <- depths$Var1[depths$Freq < subsampled_depth_allx] 
	
	## Filter IMGT Mutation data frame 
	IMGT_Mutation <- IMGT_Mutation[!IMGT_Mutation$Sample %in% to_exclude,]
	#########################################################################
	
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
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation_Summary_bar_', productivity, '.pdf'), width=8, height=5)
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
	
	## Assign day if required
	if(any(!layouts$Barcode==layouts$SampleID)){	
		q <- str_split_fixed(w$Sample, "_", 2)
		w$Barcode <- q[,1]
		w$days <- q[,2]
		w$days[w$days != 1 & w$days != 3 & w$days != 5] <- "Control"
	} else {
		w$days <- "NA"
		w$Barcode <- w$Sample
	} 
	
	# Plot box plots!!
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation_Summary_box_', productivity, '.pdf'), width=10, height=10)
	q1 <- ggplot(w, aes(fill=ClassSequence, x=days, y=Freq), na.rm = TRUE) +  geom_boxplot() +theme_bw() + xlab("Timepoint") + ylab("Frequency") + ggtitle(paste0("Total Mutations ", productivity))+labs(fill="Class of VDJ") +facet_wrap(~ClassSequence, scales="free_y") + theme(legend.position="none")
	grid.arrange(q1)
	dev.off()
	
	## Moving to wide format to save 
	s <-spread(w, key = ClassSequence, value = Freq)
	s$Barcode <- NULL
	s$days <- NULL
	write.table(s, paste0(path_to_outputdir, "/Summary/IMGT/IMGT_Prop_SHM_", productivity, ".txt"), sep="\t")
	
	
	# 3D plot
	pdf(paste0(path_to_outputdir,'Plots/IMGT/3d_scatter', productivity, '.pdf'), width=11, height=8)
	scatter3D_fancy(s[,2], s[,3], s[,4], colvar=s[,5], clab=("GC V.HIGH MUTATION"), phi = 20, theta=20, ticktype = "detailed", bty ="g", xlab="Germinal Centre", ylab="Extrafollicular", zlab="Naive", main="Proportion of BCR reads in each compartment",  pch = 20, cex = 2)
	dev.off()

	
}

imgt_mutation_statistics_sepsis <- function(path_to_outputdir = path_to_outputdir, cluster_nodes = 8, productivity=productivity, path_to_layout=path_to_layout){
	suppressMessages(library(tidyverse))
	suppressMessages(library(data.table))
	suppressMessages(library(ggplot2))
	suppressMessages(library(ggforce))
	suppressMessages(library(Gviz))
	suppressMessages(library(foreach))
	suppressMessages(library(doParallel))
	suppressMessages(library(gridExtra))
	suppressMessages(library(cowplot))
	suppressMessages(library(gtools))
	suppressMessages(library(purrr))
	suppressMessages(library(reshape2))
	suppressMessages(library(Hmisc))
	suppressMessages(library(corrplot))
	suppressMessages(library(stringr))
	suppressMessages(library(dplyr))
	suppressMessages(library(tidyr))
	suppressMessages(library(matrixStats))
	suppressMessages(library(ggpubr))
	suppressMessages(library(ggrastr))
	suppressMessages(library(ggpubr))
	suppressMessages(library(plot3D))
	
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
	
	
	## We need to exclude samples below certain read depth!
	## Use the same filter as we used in the isotyper script 
	depths <- data.frame(table(IMGT_Mutation$Sample))
	depths$Var1 <- as.character(depths$Var1)
	
	## Getting subsample depths which were used for isotyper script 
	counts_used <- paste0(path_to_outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS")
	all_files <- list.files(counts_used, full.name=TRUE)
	all_files <- grep("depth_per_isotype", all_files, value=TRUE)
	counts_used <- read.delim(all_files[1], sep="\t", header=TRUE)
	counts_used <- counts_used[counts_used$type=="UNIQ",]
	subsampled_depth_all <- counts_used$min[counts_used$X.isotype=="all"]
	
	if(productivity == "UNPRODUCTIVE"){ 
	subsampled_depth_allx <- subsampled_depth_all/10
	} else {
		subsampled_depth_allx <- subsampled_depth_all
	} 

	
	# samples exclude 
	to_exclude <- depths$Var1[depths$Freq < subsampled_depth_allx] 
	
	## Filter IMGT Mutation data frame 
	IMGT_Mutation <- IMGT_Mutation[!IMGT_Mutation$Sample %in% to_exclude,]
	#########################################################################
	
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
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation_Summary_bar_', productivity, '.pdf'), width=8, height=5)
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
	
	### So we can create boxplots!!
	## Assign day if required
	if(any(!layouts$Barcode==layouts$SampleID)){	
		q <- str_split_fixed(w$Sample, "_", 2)
		w$Barcode <- q[,1]
		w$days <- q[,2]
		w$days[w$days != 1 & w$days != 3 & w$days != 5] <- "Control"
	} else {
		w$days <- "NA"
		w$Barcode <- w$Sample
	} 
	
	## Moving to wide format to save 
	s <-spread(w, key = ClassSequence, value = Freq)
	s$Barcode <- NULL
	s$days <- NULL
	write.table(s, paste0(path_to_outputdir, "/Summary/IMGT/IMGT_Prop_SHM_", productivity, ".txt"), sep="\t")
	
	#####################
	# 3d plot
	pdf(paste0(path_to_outputdir,'Plots/IMGT/3d_scatter', productivity, '.pdf'), width=11, height=8)
	scatter3D_fancy(s[,2], s[,3], s[,4], colvar=s[,5], clab=("GC V.HIGH MUTATION"), phi = 20, theta=20, ticktype = "detailed", bty ="g", xlab="Germinal Centre", ylab="Extrafollicular", zlab="Naive", main="Proportion of BCR reads in each compartment",  pch = 20, cex = 2)
	dev.off()

	#################
	frequencies_group_t <- 	w
	frequencies_group_t$Sampleid <- frequencies_group_t$Sample
	frequencies_group_t$Sample <- NULL

	##Combine with meta data 
	meta_data <- read.delim('/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Meta_data_for_Cohort1and2.txt', sep='\t', header=TRUE)
	data_combined <- merge(frequencies_group_t, meta_data, by.x="Sampleid", by.y="SampleID", all.x=TRUE)
		
	data_combined$SRS <- as.character(data_combined$SRS)
	data_combined$SRS[data_combined$days == "Control"] <- "Not Applicable"
	data_combined$SRS_New <- as.character(data_combined$SRS_New)
	data_combined$SRS_New[data_combined$days == "Control"] <- "Not Applicable"
	
	data_combined$outcome <- NA 
	data_combined$outcome[is.na(data_combined$Days_death_from_ICU)] <- "CONTROL"
	data_combined$outcome[data_combined$Days_death_from_ICU=="Alive"] <- "ALIVE"
	data_combined$outcome[data_combined$Days_death_from_ICU!="Alive"] <- "DEAD"
	data_combined$outcome[data_combined$timepoint=="Control"] <- "CONTROL"
	data_combined$X7DAY_MORTALITY[is.na(data_combined$X7DAY_MORTALITY)] <- "Control"
	data_combined$X7DAY_MORTALITY[data_combined$X7DAY_MORTALITY==0] <- "7 Day Mortality"
	data_combined$X7DAY_MORTALITY[data_combined$X7DAY_MORTALITY==1] <- "Post 7 Day Mortality"
	data_combined$X7DAY_MORTALITY[data_combined$X7DAY_MORTALITY==2] <- "Alive"
	data_combined$X7DAY_MORTALITY <- factor(data_combined$X7DAY_MORTALITY, levels = c("Control", "Alive", "Post 7 Day Mortality", "7 Day Mortality"))
	
	
	## Setting Up data Structure
	data_combined$SOFA_total <- as.numeric(data_combined$SOFA_total)
	data_combined$Age <- as.numeric(as.character(data_combined$Age))
	data_combined$Sex[is.na(data_combined$Sex)] <- "Unknown"
			
	pdf(paste0(path_to_outputdir,'Plots/IMGT/IMGT_Mutation_Summary_box_', productivity, '.pdf'), width=10, height=10)
	q1 <- ggplot(frequencies_group_t, aes(fill=ClassSequence, x=days, y=Freq), na.rm = TRUE) +  geom_boxplot() +theme_bw() + xlab("Timepoint") + ylab("Frequency") + ggtitle(paste0("Total Mutations ", productivity))+labs(fill="Class of VDJ") +facet_wrap(~ClassSequence, scales="free_y") + theme(legend.position="none")
	grid.arrange(q1)
	q1 <- ggplot(data_combined, aes(fill=X7DAY_MORTALITY, x=days, y=Freq), na.rm = TRUE) +  geom_boxplot() +theme_bw() + xlab("Timepoint") + ylab("Frequency") + ggtitle(paste0("Total Mutations ", productivity))+labs(fill="7 Day Mortality") +facet_wrap(~ClassSequence, scales="free_y") 
	grid.arrange(q1)
	q1 <- ggplot(data_combined, aes(fill=SRS_New, x=days, y=Freq), na.rm = TRUE) +  geom_boxplot() +theme_bw() + xlab("Timepoint") + ylab("Frequency") + ggtitle(paste0("Total Mutations ", productivity))+labs(fill="SRS_NEW") +facet_wrap(~ClassSequence, scales="free_y") 
	grid.arrange(q1)
	q1 <- ggplot(data_combined, aes(fill=SRS, x=days, y=Freq), na.rm = TRUE) +  geom_boxplot() +theme_bw() + xlab("Timepoint") + ylab("Frequency") + ggtitle(paste0("Total Mutations ", productivity))+labs(fill="SRS") +facet_wrap(~ClassSequence, scales="free_y") 
	grid.arrange(q1)
	q1 <- ggplot(data_combined, aes(fill=Shock_sepsis2, x=days, y=Freq), na.rm = TRUE) +  geom_boxplot() +theme_bw() + xlab("Timepoint") + ylab("Frequency") + ggtitle(paste0("Total Mutations ", productivity))+labs(fill="Shock Sepsis2") +facet_wrap(~ClassSequence, scales="free_y") 
	grid.arrange(q1)
	q1 <- ggplot(data_combined, aes(fill=EBV_positive_ddpcr, x=days, y=Freq), na.rm = TRUE) +  geom_boxplot() +theme_bw() + xlab("Timepoint") + ylab("Frequency") + ggtitle(paste0("Total Mutations ", productivity))+labs(fill="EBV ddPCR") +facet_wrap(~ClassSequence, scales="free_y") 
	grid.arrange(q1)
	q1 <- ggplot(data_combined, aes(fill=EBV_positive_metagenomics, x=days, y=Freq), na.rm = TRUE) +  geom_boxplot() +theme_bw() + xlab("Timepoint") + ylab("Frequency") + ggtitle(paste0("Total Mutations ", productivity))+labs(fill="EBV metagenomics") +facet_wrap(~ClassSequence, scales="free_y") 
	grid.arrange(q1)
	dev.off()

	s <-spread(w, key = ClassSequence, value = Freq)
	S_NEW <- merge(s, meta_data, by.x="Sample", by.y="SampleID", all.x=TRUE)
	
	S_NEW$SRS <- as.character(S_NEW$SRS)
	S_NEW$SRS[S_NEW$days == "Control"] <- "Not Applicable"
	S_NEW$SRS[is.na(S_NEW$SRS)] <- "Unknown"
	S_NEW$SRS_New <- as.character(S_NEW$SRS_New)
	S_NEW$SRS_New[S_NEW$days == "Control"] <- "Not Applicable"
	S_NEW$SRS_New[is.na(S_NEW$SRS_New)] <- "Unknown"
	
	S_NEW$outcome <- NA 
	S_NEW$outcome[is.na(S_NEW$Days_death_from_ICU)] <- "CONTROL"
	S_NEW$outcome[S_NEW$Days_death_from_ICU=="Alive"] <- "ALIVE"
	S_NEW$outcome[S_NEW$Days_death_from_ICU!="Alive"] <- "DEAD"
	S_NEW$outcome[S_NEW$timepoint=="Control"] <- "CONTROL"
	S_NEW$X7DAY_MORTALITY[is.na(S_NEW$X7DAY_MORTALITY)] <- "Control"
	S_NEW$X7DAY_MORTALITY[S_NEW$X7DAY_MORTALITY==0] <- "7 Day Mortality"
	S_NEW$X7DAY_MORTALITY[S_NEW$X7DAY_MORTALITY==1] <- "Post 7 Day Mortality"
	S_NEW$X7DAY_MORTALITY[S_NEW$X7DAY_MORTALITY==2] <- "Alive"
	
	pdf(paste0(path_to_outputdir,'Plots/IMGT/3D_Scatter_sepsis', productivity, '.pdf'), width=20, height=15)
	plot3D::scatter3D(S_NEW[,4], S_NEW[,5], S_NEW[,6], colvar=as.numeric(as.factor(S_NEW$days)), col=c("green", "blue", "yellow", "red"), phi = 20, theta=20, ticktype = "detailed", bty ="g", xlab="Germinal Centre", ylab="Extrafollicular", zlab="Naive", main="Proportion of BCR reads in each compartment",  pch = 20, cex = 2, colkey = list(at = c(1, 2, 3, 4), side = 1,  labels = c("Day 1", "Day 3", "Day 5", "Control")))
	plot3D::scatter3D(S_NEW[,4], S_NEW[,5], S_NEW[,6], colvar=as.numeric(as.factor(S_NEW$Shock_sepsis2)), col=c("green", "blue", "yellow"), phi = 20, theta=20, ticktype = "detailed", bty ="g", xlab="Germinal Centre", ylab="Extrafollicular", zlab="Naive", main="Proportion of BCR reads in each compartment",  pch = 20, cex = 2, colkey = list(at = c(1, 2, 3), side = 1,  labels = c("N", "Unknown", "Y")))
	plot3D::scatter3D(S_NEW[,4], S_NEW[,5], S_NEW[,6], colvar=as.numeric(as.factor(S_NEW$X7DAY_MORTALITY)), col=c("green", "blue", "yellow", "red"), phi = 20, theta=20, ticktype = "detailed", bty ="g", xlab="Germinal Centre", ylab="Extrafollicular", zlab="Naive", main="Proportion of BCR reads in each compartment",  pch = 20, cex = 2, colkey = list(at = c(1, 2, 3, 3.9), side = 1,  labels = c("7 Day Mortality", "Alive", "Control", "Post 7 Day Mortality")))
	plot3D::scatter3D(S_NEW[,4], S_NEW[,5], S_NEW[,6], colvar=as.numeric(as.factor(S_NEW$ou)), col=c("green", "blue", "yellow"), phi = 20, theta=20, ticktype = "detailed", bty ="g", xlab="Germinal Centre", ylab="Extrafollicular", zlab="Naive", main="Proportion of BCR reads in each compartment",  pch = 20, cex = 2, colkey = list(at = c(1, 2, 3), side = 1,  labels = c("Alive", "Control", "Dead")))	
	p <- S_NEW[S_NEW$SRS_New=="1" | S_NEW$SRS_New=="2" | S_NEW$SRS_New=="Not Applicable",]
	plot3D::scatter3D(p[,4], p[,5], p[,6], colvar=as.numeric(as.factor(p$SRS_New)), col=c("green", "blue", "red"), phi = 20, theta=20, ticktype = "detailed", bty ="g", xlab="Germinal Centre", ylab="Extrafollicular", zlab="Naive", main="Proportion of BCR reads in each compartment",  pch = 20, cex = 2, colkey = list(at = c(1, 2, 3), side = 1,  labels = c("SRS 1", "SRS 2", "Controls")))
	dev.off()
	
}