## Function to Calculate the JACCARD Index as part of the BCR TCR preprocessing pipeline.
## Look for sample similarity and identify any mismatches
## will subsample to 0.9xMinVDJCount

## Lauren Overend

library(tidyverse)
library(data.table)
library(ggplot2)
library(ggforce)
library(Gviz)
library(foreach)
library(doParallel)
library(gridExtra)
library(cowplot)
library(optparse)
library(gtools)


`%notin%` <- Negate(`%in%`)	


calculate_jaccard <- function(seq_1, seq_2, depth){
	jaccard <- c()
	shared_seq_count <- c()
	sample_count  <- c()
	for(i in 1:10000){
		bcr_1_sample <- sample(seq_1, depth, replace=FALSE)
		bcr_2_sample <- sample(seq_2, depth, replace=FALSE)
		a <- length(bcr_1_sample[bcr_1_sample %in% bcr_2_sample])
		b <- length(bcr_1_sample[bcr_1_sample %notin% bcr_2_sample])
		cc <- length(bcr_2_sample[bcr_2_sample %notin% bcr_1_sample])
		jaccard_subset <- a /(a + b + cc) 
		jaccard <- c(jaccard, jaccard_subset)
		shared_seq_count <- c(shared_seq_count, a) 
		sample_count <- c(sample_count, (depth+depth))
	} 
	mean_jaccard_subsampled <- mean(jaccard)
	mean_no_shared_sequences <- mean(shared_seq_count)
	mean_sample_count <- mean(sample_count)
	## nosubsampling
	a_full <- length(seq_1[seq_1 %in% seq_2])
	b_full <- length(seq_1[seq_1 %notin% seq_2])
	cc_full <- length(seq_2[seq_2 %notin% seq_1])
	jaccard_full <- a_full /(a_full + b_full + cc_full) 
    shared_seq_count_full <- a_full
	sample_count_full <- (a_full + b_full + cc_full)
	results <- c(mean_no_shared_sequences,mean_sample_count, mean_jaccard_subsampled, shared_seq_count_full, sample_count_full, jaccard_full)
	return(results)
	}



calculate_jaccard_matrix <- function(path_to_output){
	path <- paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS/")
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("Att", samples, value=TRUE)
	retry1 <- combinations(length(samples), 2, samples, repeats.allowed=TRUE)
	depth_file <- list.files(paste0(path_to_output, "ORIENTATED_SEQUENCES/ANNOTATIONS"), full.name=TRUE)
	depth_file <- grep("Sampling_depth", depth_file, value=TRUE)
	depth_file <- read.delim(depth_file)
	subsample <- depth_file[1,3]
	subsample_depth <- floor(0.9*subsample)
	#Register doParrallel
	cl <- 20
	registerDoParallel(cl)
	
	## Calculate Jaccard Matrix 
	JACCARD_MATRIX <- foreach(i = 1:dim(retry1)[1], .combine=rbind) %dopar% {
		incl <- retry1[i,]
		#print(i) 
		# Read in the sequences and replicate based on constant region counts. 
		bcr_1 <- read.delim(incl[1], header=FALSE)
		bcr_1_sequence <- bcr_1$V3
		bcr_1_sequence_w <- rep(bcr_1$V3, bcr_1$V2)
		bcr_2 <- read.delim(incl[2], header=FALSE)
		bcr_2_sequence <- bcr_2$V3
		bcr_2_sequence_w <- rep(bcr_2$V3, bcr_2$V2)
		## Calculate Jaccard on Weighed  and Non Weighted Repertoire
		jaccard_results <- calculate_jaccard(bcr_1_sequence, bcr_2_sequence, subsample_depth)
		jaccard_results_w <- calculate_jaccard(bcr_1_sequence_w, bcr_2_sequence_w, subsample_depth)
		results <- c(incl[1], incl[2], jaccard_results,  jaccard_results_w)
		return(results)
	}
	colnames(JACCARD_MATRIX) <- c("Sample1", "Sample2", "SharedSeq.MeanSubsample", "Size.MeanSubsample", "Jaccard.MeanSubsample", "SharedSeq.Full", "Size.Full", "Jaccard.Full", "SharedSeq.MeanSubsample.Weighted", "Size.MeanSubsample.Weighted", "Jaccard.MeanSubsample.Weighted", "SharedSeq.Full.Weighted", "Size.Full.Weighted", "Jaccard.Full.Weighted")
	JACCARD_MATRIX <- data.frame(JACCARD_MATRIX)
	write.table(JACCARD_MATRIX, paste0(path_to_output, "Summary/JACCARD_MATRIX_FINAL.txt"), sep='\t') 
	## SAVED MATRIX

	## RENAME SAMPLE NAMES: 
	JACCARD_MATRIX_2 <- data.frame(JACCARD_MATRIX)
	JACCARD_MATRIX_2$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS//Att_BCR_"), "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS//Att_BCR_"), "", JACCARD_MATRIX_2$Sample2)
	JACCARD_MATRIX_2$Sample1 <- gsub(".txt", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub(".txt", "", JACCARD_MATRIX_2$Sample2)
	write.table(JACCARD_MATRIX_2, paste0(path_to_output, "Summary/JACCARD_MATRIX_FINAL_EDITED.txt"), sep='\t')
	
	## Reformat Dataframe for further analysis 
	JACCARD_MATRIX_2$SharedSeq.MeanSubsample <- as.numeric(JACCARD_MATRIX_2$SharedSeq.MeanSubsample)
	JACCARD_MATRIX_2$Size.MeanSubsample <- as.numeric(JACCARD_MATRIX_2$Size.MeanSubsample)
	JACCARD_MATRIX_2$Jaccard.MeanSubsample <- as.numeric(JACCARD_MATRIX_2$Jaccard.MeanSubsample)
	JACCARD_MATRIX_2$SharedSeq.Full <- as.numeric(JACCARD_MATRIX_2$SharedSeq.Full)
	JACCARD_MATRIX_2$Size.Full <- as.numeric(JACCARD_MATRIX_2$Size.Full)
	JACCARD_MATRIX_2$Jaccard.Full <- as.numeric(JACCARD_MATRIX_2$Jaccard.Full)

	JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted <- as.numeric(JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted)
	JACCARD_MATRIX_2$Size.MeanSubsample.Weighted <- as.numeric(JACCARD_MATRIX_2$Size.MeanSubsample.Weighted)
	JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted <- as.numeric(JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted)
	JACCARD_MATRIX_2$SharedSeq.Full.Weighted <- as.numeric(JACCARD_MATRIX_2$SharedSeq.Full.Weighted)
	JACCARD_MATRIX_2$Size.Full.Weighted <- as.numeric(JACCARD_MATRIX_2$Size.Full.Weighted)
	JACCARD_MATRIX_2$Jaccard.Full.Weighted <- as.numeric(JACCARD_MATRIX_2$Jaccard.Full.Weighted)


	JACCARD_MATRIX_2$log.SharedSeq.MeanSubsample <- log(as.numeric(JACCARD_MATRIX_2$SharedSeq.MeanSubsample))
	JACCARD_MATRIX_2$log.Size.MeanSubsample <- log(as.numeric(JACCARD_MATRIX_2$Size.MeanSubsample))
	JACCARD_MATRIX_2$log.Jaccard.MeanSubsample <- log(as.numeric(JACCARD_MATRIX_2$Jaccard.MeanSubsample))
	JACCARD_MATRIX_2$log.SharedSeq.Full <- log(as.numeric(JACCARD_MATRIX_2$SharedSeq.Full))
	JACCARD_MATRIX_2$log.Size.Full <- log(as.numeric(JACCARD_MATRIX_2$Size.Full))
	JACCARD_MATRIX_2$log.Jaccard.Full <- log(as.numeric(JACCARD_MATRIX_2$Jaccard.Full))

	JACCARD_MATRIX_2$log.SharedSeq.MeanSubsample.Weighted <- log(as.numeric(JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted))
	JACCARD_MATRIX_2$log.Size.MeanSubsample.Weighted <- log(as.numeric(JACCARD_MATRIX_2$Size.MeanSubsample.Weighted))
	JACCARD_MATRIX_2$log.Jaccard.MeanSubsample.Weighted <- log(as.numeric(JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted))
	JACCARD_MATRIX_2$log.SharedSeq.Full.Weighted <- log(as.numeric(JACCARD_MATRIX_2$SharedSeq.Full.Weighted))
	JACCARD_MATRIX_2$log.Size.Full.Weighted <- log(as.numeric(JACCARD_MATRIX_2$Size.Full.Weighted))
	JACCARD_MATRIX_2$log.Jaccard.Full.Weighted <- log(as.numeric(JACCARD_MATRIX_2$Jaccard.Full.Weighted))
	write.table(JACCARD_MATRIX_2, paste0(path_to_output, "Summary/JACCARD_MATRIX_FINAL_EDITED.txt"), sep='\t')
		
	pdf(paste0(path_to_output, "Plots/JaccardAllSamples.pdf"), width=23, height=10)
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=SharedSeq.MeanSubsample)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Size.MeanSubsample)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Jaccard.MeanSubsample)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=SharedSeq.Full)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Size.Full)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Jaccard.Full)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.SharedSeq.MeanSubsample)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Size.MeanSubsample)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Jaccard.MeanSubsample)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.SharedSeq.Full)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Size.Full)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Jaccard.Full)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	plot(plot_grid(e, k, ncol=2))
	plot(plot_grid(g, m, ncol=2))
	plot(plot_grid(h, n, ncol=2))
	plot(plot_grid(j, p, ncol=2))
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=SharedSeq.MeanSubsample.Weighted)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Size.MeanSubsample.Weighted)) + geom_tile()+ theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Jaccard.MeanSubsample.Weighted)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=SharedSeq.Full.Weighted)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Size.Full.Weighted)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Jaccard.Full.Weighted)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.SharedSeq.MeanSubsample.Weighted)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Size.MeanSubsample.Weighted))   + geom_tile()+ theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Jaccard.MeanSubsample.Weighted)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.SharedSeq.Full.Weighted)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Size.Full.Weighted)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Jaccard.Full.Weighted)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))
	plot(plot_grid(e, k, ncol=2))
	plot(plot_grid(g, m, ncol=2))
	plot(plot_grid(h, n, ncol=2))
	plot(plot_grid(j, p, ncol=2))
	dev.off()  
	#return(JACCARD_MATRIX_2)
	}

