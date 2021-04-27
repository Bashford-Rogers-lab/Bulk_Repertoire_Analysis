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


##-----------------------------------------------------------------------------------
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


##-----------------------------------------------------------------------------------
## Function for calculating and plotting the basic Jaccard Index on  BCR/TCR output for 'common sequences' 
## e.g. must have 2 copies of the sequence to be used. 
## Useful for looking past any low level index hopping/contamination 
## Union of shared sequences

  
calculate_jaccard_matrix_filter <- function(path_to_output, runname){
	path <- paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS/")
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("Att", samples, value=TRUE)
	retry1 <- combinations(length(samples), 2, samples, repeats.allowed=TRUE)
	mins <- c()
	## Note we look for the lowest > 200 once we calculate how many sequences >=2 (this means depth is nearly identical across filtered/not) 
	for(i in 1:length(samples)){
		a <- read.delim(samples[i], header=FALSE)
		a <- a[a$V2 >= 2,]
		a_min <- dim(a)[1]
		mins <- c(mins, a_min)
	} 
	## Ensure that susbample must be at least 200 sequences. 
	sample_depth <- min(mins[mins>=223]) 
	subsample_depth <- floor(sample_depth*0.9)
	
	#Register doParrallel (10 nodes seems to be the maximum you can run on the cluster with 1 slot or it crashes!
	cl <- 10
	registerDoParallel(cl)
	
	#dim(retry1)[1]
	## Calculate Jaccard Matrix 
	JACCARD_MATRIX <- foreach(i = 1:dim(retry1)[1], .combine=rbind) %dopar% {
		incl <- retry1[i,]
		# print i in multiples of 50(ish) to give a rough estimate of where we are in the function
		## This wont be an exact count as we are running in parrallelel therefore some nodes may be further ahead 
		if(i %% 50 == 0){
			print(i)
		}
		# Read in the sequences and replicate based on constant region counts.
		# Filter to ensure they have a minimum of 2 reads. 
		bcr_1 <- read.delim(incl[1], header=FALSE)
		bcr_1 <- bcr_1[bcr_1$V2 >= 2,]
		bcr_1_sequence <- bcr_1$V3
		bcr_1_sequence_w <- rep(bcr_1$V3, bcr_1$V2)
		bcr_2 <- read.delim(incl[2], header=FALSE)
		bcr_2  <- bcr_2 [bcr_2 $V2 >= 2,]
		bcr_2_sequence <- bcr_2$V3
		bcr_2_sequence_w <- rep(bcr_2$V3, bcr_2$V2)
		## Calculate Jaccard on Weighed  and Non Weighted Repertoire
		if(dim(bcr_1)[1] >= sample_depth  & dim(bcr_2)[1] >= sample_depth){
			jaccard_results <- calculate_jaccard(bcr_1_sequence, bcr_2_sequence, subsample_depth)
			jaccard_results_w <- calculate_jaccard(bcr_1_sequence_w, bcr_2_sequence_w, subsample_depth)
			results <- c(incl[1], incl[2], jaccard_results,  jaccard_results_w)
		} else {
			results <- c(incl[1], incl[2], NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
		}
		return(results)
	}
	
	colnames(JACCARD_MATRIX) <- c("Sample1", "Sample2", "SharedSeq.MeanSubsample", "Size.MeanSubsample", "Jaccard.MeanSubsample", "SharedSeq.Full", "Size.Full", "Jaccard.Full", "SharedSeq.MeanSubsample.Weighted", "Size.MeanSubsample.Weighted", "Jaccard.MeanSubsample.Weighted", "SharedSeq.Full.Weighted", "Size.Full.Weighted", "Jaccard.Full.Weighted")

	## RENAME SAMPLE NAMES: 
	JACCARD_MATRIX_2 <- data.frame(JACCARD_MATRIX)
	JACCARD_MATRIX_2$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS//Att_"), "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS//Att_"), "", JACCARD_MATRIX_2$Sample2)
	JACCARD_MATRIX_2$Sample1 <- gsub(".txt", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub(".txt", "", JACCARD_MATRIX_2$Sample2)
	
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
	
	## Make an extra column considering whether results were NA because they couldn't be subsampled.
	## This is used in the plotting. 
	JACCARD_MATRIX_2$ANYNAS <- NA
	JACCARD_MATRIX_2$ANYNAS[is.na(JACCARD_MATRIX_2$Jaccard.Full)] <- "YES"
	
	## Save the basic Jaccard Index. 
	write.table(JACCARD_MATRIX_2, paste0(path_to_output, "Summary/JACCARDMATRIX_2ormore_AllSAMPLES_", runname, ".txt"), sep='\t')
	
	## Generate Summary Plots
	pdf(paste0(path_to_output, "Plots/JACCARDMATRIX_2ormore_AllSAMPLES_", runname, ".pdf"), width=23, height=10)
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.MeanSubsample)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')
	f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.MeanSubsample))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')
	o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.Full))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')
	plot(plot_grid(e, k, ncol=2))
	plot(plot_grid(g, m, ncol=2))
	plot(plot_grid(h, n, ncol=2))
	plot(plot_grid(j, p, ncol=2))
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.MeanSubsample.Weighted)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.MeanSubsample.Weighted))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')
	l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.Full.Weighted))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')
	plot(plot_grid(e, k, ncol=2))
	plot(plot_grid(g, m, ncol=2))
	plot(plot_grid(h, n, ncol=2))
	plot(plot_grid(j, p, ncol=2))
	dev.off()  
	}
	
	
##-----------------------------------------------------------------------------------
## Function for calculating Jaccard Matrix and correcting for Library Hopping	
## Looks at sequences with 2 or more reads 
	
calculate_jaccard_matrix_libhopcorrection_filter <- function(path_to_output, runname, path_to_layout){
	path <- paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS/")
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("Att", samples, value=TRUE)
	retry1 <- combinations(length(samples), 2, samples, repeats.allowed=TRUE)
	retry1 <- data.frame(retry1)
	colnames(retry1) <- c("Sample1Path", "Sample2Path")
	
	#Reformat Combinations to match sample ids
	## Make new Column containing just the sample name
	retry1$Sample1 <- retry1$Sample1Path
	retry1$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS//Att_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("BCR_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRA_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRB_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRG_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRD_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0(".txt"), "", retry1$Sample1)
	retry1$Sample2 <- retry1$Sample2Path
	retry1$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS//Att_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0("BCR_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0("TCRA_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0("TCRB_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0("TCRG_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0("TCRD_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0(".txt"), "", retry1$Sample2)
	
	## Read Layout File for Sample1
	layouts_1 <- read.delim(path_to_layout, sep="\t")
	colnames(layouts_1) <- c("SampleID", "Barcode1",  "BCR_Lane_S1", "Plate_S1", "Library1", "Position1", "PCRBarcode1")
	retry1 <- merge(retry1, layouts_1, by.x="Sample1", by.y="SampleID")

	# Read Layout File for Sample2:
	layouts_2 <- read.delim(path_to_layout, sep="\t")
	colnames(layouts_2) <- c("SampleID", "Barcode2",  "BCR_Lane_S2", "Plate_S2", "Library2", "Position2", "PCRBarcode2")
	retry1 <- merge(retry1, layouts_2, by.x="Sample2", by.y="SampleID")

	## Calculating a new sample depth based on max number of sample overlaps (headache I know!)
	## Read all samples, substract the maximum possible overlap and then times by 0.9 to avoid subsampling whole dataset. 
	mins <- c()
	## Note we look for the lowest > 200 once we calculate how many sequences >=2 (this means depth is nearly identical across filtered/not) 
	for(i in 1:length(samples)){
		a <- read.delim(samples[i], header=FALSE)
		a <- a[a$V2 >= 2,]
		a_min <- dim(a)[1]
		mins <- c(mins, a_min)
	} 
	## Ensure that susbample must be at least 200 sequences. 
	sample_depth <- min(mins[mins>=223]) 
	subsample_depth <- floor(sample_depth*0.9)
	
	#Register doParrallel
	cl <- 10
	registerDoParallel(cl)
	#dim(retry1)[1]
	## Calculate Jaccard Matrix 
	JACCARD_MATRIX <- foreach(i = 1:dim(retry1)[1], .combine=rbind) %dopar% {
		incl <- retry1[i,]
		if(i %% 50 == 0){
			print(i)
		}
		#return(incl)
		# Read in the sequences and replicate based on constant region counts. 
		if(incl$Sample1 == incl$Sample2 | incl$PCRBarcode1 != incl$PCRBarcode2 | incl$Barcode1 == incl$Barcode2 | incl$BCR_Lane_S1 != incl$BCR_Lane_S2 ){
			bcr_1 <- read.delim(as.character(incl[3]), header=FALSE)
			bcr_1 <- bcr_1[bcr_1$V2 >= 2,]
			bcr_1_sequence <- bcr_1$V3
			bcr_1_sequence_w <- rep(bcr_1$V3, bcr_1$V2)
			bcr_2 <- read.delim(as.character(incl[4]), header=FALSE)
			bcr_2 <- bcr_2[bcr_2$V2 >= 2,]
			bcr_2_sequence <- bcr_2$V3
			bcr_2_sequence_w <- rep(bcr_2$V3, bcr_2$V2)
		} else {
			bcr_1 <- read.delim(as.character(incl[3]), header=FALSE)
			bcr_2 <- read.delim(as.character(incl[4]), header=FALSE)
			bcr_1 <- bcr_1[bcr_1$V2 >= 2,]
			bcr_1_sequence <- bcr_1$V3
			bcr_1_sequence_w <- rep(bcr_1$V3, bcr_1$V2)
			bcr_2 <- bcr_2[bcr_2$V2 >= 2,]
			bcr_2_sequence <- bcr_2$V3
			bcr_2_sequence_w <- rep(bcr_2$V3, bcr_2$V2)
			
			shared_seq <- bcr_1_sequence[bcr_1_sequence %in% bcr_2_sequence]
			
			bcr_1_sequence <- bcr_1_sequence[bcr_1_sequence %notin% shared_seq]
			bcr_1_sequence_w <- bcr_1_sequence_w[bcr_1_sequence_w %notin% shared_seq]
			
			bcr_2_sequence <- bcr_2_sequence[bcr_2_sequence %notin% shared_seq]
			bcr_2_sequence_w <- bcr_2_sequence_w[bcr_2_sequence_w %notin% shared_seq]
		}
		
		## Calculate Jaccard on Weighed  and Non Weighted Repertoire
		if(length(bcr_1_sequence) >= sample_depth   & length(bcr_2_sequence) >= sample_depth){
			jaccard_results <- calculate_jaccard(bcr_1_sequence, bcr_2_sequence, subsample_depth)
			jaccard_results_w <- calculate_jaccard(bcr_1_sequence_w, bcr_2_sequence_w, subsample_depth)
			results <- c(as.character(incl[3]), as.character(incl[4]), jaccard_results,  jaccard_results_w)
		} else {
			results <- c(as.character(incl[3]), as.character(incl[4]), NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
		}
		if(incl$Sample1 == incl$Sample2 | incl$PCRBarcode1 != incl$PCRBarcode2 | incl$Barcode1 == incl$Barcode2| incl$BCR_Lane_S1 != incl$BCR_Lane_S2){
			results <- c(results, "NO")
			} else {
			results <- c(results, "YES")
		}
		return(results)
	}
	
	colnames(JACCARD_MATRIX) <- c("Sample1", "Sample2", "SharedSeq.MeanSubsample", "Size.MeanSubsample", "Jaccard.MeanSubsample", "SharedSeq.Full", "Size.Full", "Jaccard.Full", "SharedSeq.MeanSubsample.Weighted", "Size.MeanSubsample.Weighted", "Jaccard.MeanSubsample.Weighted", "SharedSeq.Full.Weighted", "Size.Full.Weighted", "Jaccard.Full.Weighted", "LibCorrected")

	## RENAME SAMPLE NAMES: 
	JACCARD_MATRIX_2 <- data.frame(JACCARD_MATRIX)
	JACCARD_MATRIX_2$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS//Att_"), "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS//Att_"), "", JACCARD_MATRIX_2$Sample2)
	JACCARD_MATRIX_2$Sample1 <- gsub(".txt", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub(".txt", "", JACCARD_MATRIX_2$Sample2)
	
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
	
	## Make an extra column considering whether results were NA because they couldn't be subsampled.
	## This is used in the plotting. 
	JACCARD_MATRIX_2$ANYNAS <- NA
	JACCARD_MATRIX_2$ANYNAS[is.na(JACCARD_MATRIX_2$Jaccard.Full)] <- "YES"
	
	## Save the basic Jaccard Index
	write.table(JACCARD_MATRIX_2, paste0(path_to_output, "Summary/JACCARDMATRIX_2ormore_AllSAMPLES_LIBHOP_CORRECTED_", runname, ".txt"), sep='\t')
	
	## Generate Summary Plots
	pdf(paste0(path_to_output, "Plots/JACCARDMATRIX_2ormore_AllSAMPLES_LIBHOP_CORRECTED_", runname, ".pdf"), width=23, height=10)
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.MeanSubsample)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.MeanSubsample))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000')) 
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000')) 
	o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.Full))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	plot(plot_grid(e, k, ncol=2))
	plot(plot_grid(g, m, ncol=2))
	plot(plot_grid(h, n, ncol=2))
	plot(plot_grid(j, p, ncol=2))
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.MeanSubsample.Weighted)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.MeanSubsample.Weighted))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000')) 
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000')) 
	o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.Full.Weighted))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	plot(plot_grid(e, k, ncol=2))
	plot(plot_grid(g, m, ncol=2))
	plot(plot_grid(h, n, ncol=2))
	plot(plot_grid(j, p, ncol=2))
	dev.off()  
	}
	
##-----------------------------------------------------------------------------------
## Function for calculating Jaccard Matrix and correcting for Library Contamination	
## Looks at sequences with 2 or more reads 

calculate_jaccard_matrix_libcontam_correction_filter <- function(path_to_output, runname, path_to_layout){
	path <- paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS/")
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("Att", samples, value=TRUE)
	retry1 <- combinations(length(samples), 2, samples, repeats.allowed=TRUE)
	retry1 <- data.frame(retry1)
	colnames(retry1) <- c("Sample1Path", "Sample2Path")
	
	#Reformat Combinations to match sample ids
	## Make new Column containing just the sample name
	retry1$Sample1 <- retry1$Sample1Path
	retry1$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS//Att_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("BCR_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRA_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRB_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRG_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRD_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0(".txt"), "", retry1$Sample1)
	retry1$Sample2 <- retry1$Sample2Path
	retry1$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS//Att_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0("BCR_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0("TCRA_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0("TCRB_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0("TCRG_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0("TCRD_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0(".txt"), "", retry1$Sample2)
	
	## Read Layout File for Sample1
	layouts_1 <- read.delim(path_to_layout, sep="\t")
	colnames(layouts_1) <- c("SampleID", "Barcode1",  "BCR_Lane_S1", "Plate_S1", "Library1", "Position1", "PCRBarcode1")
	retry1 <- merge(retry1, layouts_1, by.x="Sample1", by.y="SampleID")

	# Read Layout File for Sample2:
	layouts_2 <- read.delim(path_to_layout, sep="\t")
	colnames(layouts_2) <- c("SampleID", "Barcode2",  "BCR_Lane_S2", "Plate_S2", "Library2", "Position2", "PCRBarcode2")
	retry1 <- merge(retry1, layouts_2, by.x="Sample2", by.y="SampleID")

	## Calculating a new sample depth based on max number of sample overlaps (headache I know!)
	## Read all samples, substract the maximum possible overlap and then times by 0.9 to avoid subsampling whole dataset. 
	mins <- c()
	## Note we look for the lowest > 200 once we calculate how many sequences >=2 (this means depth is nearly identical across filtered/not) 
	for(i in 1:length(samples)){
		a <- read.delim(samples[i], header=FALSE)
		a <- a[a$V2 >= 2,]
		a_min <- dim(a)[1]
		mins <- c(mins, a_min)
	} 
	## Ensure that susbample must be at least 200 sequences. 
	sample_depth <- min(mins[mins>=223]) 
	subsample_depth <- floor(sample_depth*0.9)
	
	#Register doParrallel
	cl <- 10
	registerDoParallel(cl)
	#dim(retry1)[1]
	## Calculate Jaccard Matrix 
	JACCARD_MATRIX <- foreach(i = 1:dim(retry1)[1], .combine=rbind) %dopar% {
		incl <- retry1[i,]
		if(i %% 50 == 0){
			print(i)
		}
		#return(incl)
		# Read in the sequences and replicate based on constant region counts. 
		if(incl$Sample1 == incl$Sample2 | incl$PCRBarcode1 != incl$PCRBarcode2 | incl$Barcode1 == incl$Barcode2){
			bcr_1 <- read.delim(as.character(incl[3]), header=FALSE)
			bcr_1 <- bcr_1[bcr_1$V2 >= 2,]
			bcr_1_sequence <- bcr_1$V3
			bcr_1_sequence_w <- rep(bcr_1$V3, bcr_1$V2)
			bcr_2 <- read.delim(as.character(incl[4]), header=FALSE)
			bcr_2 <- bcr_2[bcr_2$V2 >= 2,]
			bcr_2_sequence <- bcr_2$V3
			bcr_2_sequence_w <- rep(bcr_2$V3, bcr_2$V2)
		} else {
			bcr_1 <- read.delim(as.character(incl[3]), header=FALSE)
			bcr_2 <- read.delim(as.character(incl[4]), header=FALSE)
			bcr_1 <- bcr_1[bcr_1$V2 >= 2,]
			bcr_1_sequence <- bcr_1$V3
			bcr_1_sequence_w <- rep(bcr_1$V3, bcr_1$V2)
			bcr_2 <- bcr_2[bcr_2$V2 >= 2,]
			bcr_2_sequence <- bcr_2$V3
			bcr_2_sequence_w <- rep(bcr_2$V3, bcr_2$V2)
			
			shared_seq <- bcr_1_sequence[bcr_1_sequence %in% bcr_2_sequence]
			
			bcr_1_sequence <- bcr_1_sequence[bcr_1_sequence %notin% shared_seq]
			bcr_1_sequence_w <- bcr_1_sequence_w[bcr_1_sequence_w %notin% shared_seq]
			
			bcr_2_sequence <- bcr_2_sequence[bcr_2_sequence %notin% shared_seq]
			bcr_2_sequence_w <- bcr_2_sequence_w[bcr_2_sequence_w %notin% shared_seq]
		}
		
		## Calculate Jaccard on Weighed  and Non Weighted Repertoire
		if(length(bcr_1_sequence) >= sample_depth   & length(bcr_2_sequence)>= sample_depth){
			jaccard_results <- calculate_jaccard(bcr_1_sequence, bcr_2_sequence, subsample_depth)
			jaccard_results_w <- calculate_jaccard(bcr_1_sequence_w, bcr_2_sequence_w, subsample_depth)
			results <- c(as.character(incl[3]), as.character(incl[4]), jaccard_results,  jaccard_results_w)
		} else {
			results <- c(as.character(incl[3]), as.character(incl[4]), NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
		}
		if(incl$Sample1 == incl$Sample2 | incl$PCRBarcode1 != incl$PCRBarcode2 | incl$Barcode1 == incl$Barcode2){
			results <- c(results, "NO")
			} else {
			results <- c(results, "YES")
		}
		return(results)
	}

	colnames(JACCARD_MATRIX) <- c("Sample1", "Sample2", "SharedSeq.MeanSubsample", "Size.MeanSubsample", "Jaccard.MeanSubsample", "SharedSeq.Full", "Size.Full", "Jaccard.Full", "SharedSeq.MeanSubsample.Weighted", "Size.MeanSubsample.Weighted", "Jaccard.MeanSubsample.Weighted", "SharedSeq.Full.Weighted", "Size.Full.Weighted", "Jaccard.Full.Weighted", "LibCorrected")

	## RENAME SAMPLE NAMES: 
	JACCARD_MATRIX_2 <- data.frame(JACCARD_MATRIX)
	JACCARD_MATRIX_2$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS//Att_"), "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS//Att_"), "", JACCARD_MATRIX_2$Sample2)
	JACCARD_MATRIX_2$Sample1 <- gsub(".txt", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub(".txt", "", JACCARD_MATRIX_2$Sample2)
	
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
	
	## Make an extra column considering whether results were NA because they couldn't be subsampled.
	## This is used in the plotting. 
	JACCARD_MATRIX_2$ANYNAS <- NA
	JACCARD_MATRIX_2$ANYNAS[is.na(JACCARD_MATRIX_2$Jaccard.Full)] <- "YES"
	
	## Save the basic Jaccard Index
	write.table(JACCARD_MATRIX_2, paste0(path_to_output, "Summary/JACCARDMATRIX_2ormore_AllSAMPLES_LIBCONTAM_CORRECTED_filter_", runname, ".txt"), sep='\t')
	
	## Generate Summary Plots
	pdf(paste0(path_to_output, "Plots/JACCARDMATRIX_2ormore_AllSAMPLES_LIBCONTAM_CORRECTED_filter_", runname, ".pdf"), width=23, height=10)
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.MeanSubsample)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.MeanSubsample))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000')) 
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000')) 
	o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.Full))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	plot(plot_grid(e, k, ncol=2))
	plot(plot_grid(g, m, ncol=2))
	plot(plot_grid(h, n, ncol=2))
	plot(plot_grid(j, p, ncol=2))
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.MeanSubsample.Weighted)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.MeanSubsample.Weighted))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000')) 
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000')) 
	o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.Full.Weighted))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence Corrected", values = c("green", '#00000000'))
	plot(plot_grid(e, k, ncol=2))
	plot(plot_grid(g, m, ncol=2))
	plot(plot_grid(h, n, ncol=2))
	plot(plot_grid(j, p, ncol=2))
	dev.off()  
	}


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	