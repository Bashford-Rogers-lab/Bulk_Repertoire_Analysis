## Function to Calculate the JACCARD Index as part of the BCR TCR preprocessing pipeline.
## Look for sample similarity and identify any mismatches
## will subsample to 0.9xMinVDJCount

## Lauren Overend

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggforce))
suppressMessages(library(Gviz))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(optparse))
suppressMessages(library(gtools))

## Function for calculating and plotting the basic Jaccard Index on  BCR/TCR output. 
## Union of shared sequences. 
options(expressions= 200000)

calculate_jaccard_matrix_UMI_RAW <- function(path_to_output, runname){
	path <- paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS/")
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("Att", samples, value=TRUE)
	retry1 <- combinations(length(samples), 2, samples, repeats.allowed=TRUE)
	mins <- c()
	## Note we look for the lowest > 200 once we calculate how many sequences >=2 (this means depth is nearly identical across filtered/not) 
	for(i in 1:length(samples)){
		a <- read.delim(samples[i], header=FALSE)
		a <- a[a$V2 >= 2,]
		a <- a$V2
		a_min <- length(a)
		mins <- c(mins, a_min)
	} 
	## Ensure that susbample must be at least 200 sequences. 
	if(any(mins>=200)){
		sample_depth <- min(mins[mins>=223])
	} else {
		sample_depth <- 223
	}
	
	subsample_depth <- floor(sample_depth*0.9)

	
	##Extract the barcodes file
	path <- paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP/")
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("Barcode_filtering_information_", samples, value=TRUE)
	retry1 <- combinations(length(samples), 2, samples, repeats.allowed=TRUE)
	mins <- c()
		
	retry1 <- data.frame(retry1)
	colnames(retry1) <- c("Sample1Path", "Sample2Path")
	
	## Get sample id
	retry1$Sample1 <- retry1$Sample1Path
	retry1$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0(".txt"), "", retry1$Sample1)
	
	retry1$Sample2 <- retry1$Sample2Path
	retry1$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0(".txt"), "", retry1$Sample2)
	
	
	#Register doParrallel (10 nodes seems to be the maximum you can run on the cluster with 1 slot or it crashes!
	cl <- 20
	registerDoParallel(cl)
	
	#dim(retry1)[1]
	## Calculate Jaccard Matrix 
	JACCARD_MATRIX <- foreach(i = 1:dim(retry1)[1], .combine=rbind) %dopar% {
		incl <- retry1[i,]
		# print i in multiples of 50(ish) to give a rough estimate of where we are in the function
		## This wont be an exact count as we are running in parrallelel therefore some nodes may be further ahead 
		if(i %% 100 == 0){
			print(i)
		}
		# Read in the sequences and replicate based on constant region counts. 
		bcr_1 <- read.delim(as.character(incl[1]), header=TRUE)		
		bcr_2 <- read.delim(as.character(incl[2]), header=TRUE)

		# Read in the sequences and replicate based on constant region counts.
		# Filter to ensure they have a minimum of 2 reads. 
		bcr_1_sequence <- bcr_1$J_barcode
		bcr_1_sequence_w <- rep(bcr_1$J_barcode, bcr_1$total_reads_with_BC)
		
		bcr_2_sequence <- bcr_2$J_barcode
		bcr_2_sequence_w <- rep(bcr_2$J_barcode, bcr_2$total_reads_with_BC)
		
		
		## Calculate Jaccard on Weighed  and Non Weighted Repertoire
		if(dim(bcr_1)[1] >= sample_depth  & dim(bcr_2)[1] >= sample_depth){
			jaccard_results <- calculate_jaccard(bcr_1_sequence, bcr_2_sequence, subsample_depth)
			jaccard_results_w <- calculate_jaccard(bcr_1_sequence_w, bcr_2_sequence_w, subsample_depth)
			results <- c(as.character(incl[1]), as.character(incl[2]), jaccard_results,  jaccard_results_w)
		} else {
			results <- c(as.character(incl[1]), as.character(incl[2]), NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
		}
		return(results)
	}
	colnames(JACCARD_MATRIX) <- c("Sample1", "Sample2", "SharedSeq.MeanSubsample", "Size.MeanSubsample", "Jaccard.MeanSubsample", "SharedSeq.Full", "Size.Full", "Jaccard.Full", "SharedSeq.MeanSubsample.Weighted", "Size.MeanSubsample.Weighted", "Jaccard.MeanSubsample.Weighted", "SharedSeq.Full.Weighted", "Size.Full.Weighted", "Jaccard.Full.Weighted")

	## RENAME SAMPLE NAMES: 
	JACCARD_MATRIX_2 <- data.frame(JACCARD_MATRIX)
	JACCARD_MATRIX_2$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", JACCARD_MATRIX_2$Sample2)
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
	
	c <- min(JACCARD_MATRIX_2$SharedSeq.MeanSubsample[JACCARD_MATRIX_2$SharedSeq.MeanSubsample!=0 & !is.na(JACCARD_MATRIX_2$SharedSeq.MeanSubsample)])/2
	JACCARD_MATRIX_2$log.SharedSeq.MeanSubsample <- log2(as.numeric(JACCARD_MATRIX_2$SharedSeq.MeanSubsample)+c)
	
	c <- min(JACCARD_MATRIX_2$Size.MeanSubsample[JACCARD_MATRIX_2$Size.MeanSubsample!=0 & !is.na(JACCARD_MATRIX_2$Size.MeanSubsample)])/2
	JACCARD_MATRIX_2$log.Size.MeanSubsample <- log2(as.numeric(JACCARD_MATRIX_2$Size.MeanSubsample)+c)
	
	c <- min(JACCARD_MATRIX_2$Jaccard.MeanSubsample[JACCARD_MATRIX_2$Jaccard.MeanSubsample!=0 & !is.na(JACCARD_MATRIX_2$Jaccard.MeanSubsample)])/2
	JACCARD_MATRIX_2$log.Jaccard.MeanSubsample <- log2(as.numeric(JACCARD_MATRIX_2$Jaccard.MeanSubsample)+c)
	
	c <- min(JACCARD_MATRIX_2$SharedSeq.Full[JACCARD_MATRIX_2$SharedSeq.Full!=0 & !is.na(JACCARD_MATRIX_2$SharedSeq.Full)])/2
	JACCARD_MATRIX_2$log.SharedSeq.Full <- log2(as.numeric(JACCARD_MATRIX_2$SharedSeq.Full)+c)
	
	c <- min(JACCARD_MATRIX_2$Size.Full[JACCARD_MATRIX_2$Size.Full!=0 & !is.na(JACCARD_MATRIX_2$Size.Full)])/2
	JACCARD_MATRIX_2$log.Size.Full <- log2(as.numeric(JACCARD_MATRIX_2$Size.Full)+c)
	
	c <- min(JACCARD_MATRIX_2$Jaccard.Full[JACCARD_MATRIX_2$Jaccard.Full!=0 & !is.na(JACCARD_MATRIX_2$Jaccard.Full)])/2
	JACCARD_MATRIX_2$log.Jaccard.Full <- log2(as.numeric(JACCARD_MATRIX_2$Jaccard.Full)+c)
	
	c <- min(JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted[JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted !=0 & !is.na(JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted )])/2
	JACCARD_MATRIX_2$log.SharedSeq.MeanSubsample.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$Size.MeanSubsample.Weighted[JACCARD_MATRIX_2$Size.MeanSubsample.Weighted !=0 & !is.na(JACCARD_MATRIX_2$Size.MeanSubsample.Weighted)])/2
	JACCARD_MATRIX_2$log.Size.MeanSubsample.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$Size.MeanSubsample.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted[JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted !=0 & !is.na(JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted)])/2
	JACCARD_MATRIX_2$log.Jaccard.MeanSubsample.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$SharedSeq.Full.Weighted[JACCARD_MATRIX_2$SharedSeq.Full.Weighted!=0 & !is.na(JACCARD_MATRIX_2$SharedSeq.Full.Weighted)])/2
	JACCARD_MATRIX_2$log.SharedSeq.Full.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$SharedSeq.Full.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$Size.Full.Weighted[JACCARD_MATRIX_2$Size.Full.Weighted!=0 & !is.na(JACCARD_MATRIX_2$Size.Full.Weighted)])/2
	JACCARD_MATRIX_2$log.Size.Full.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$Size.Full.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$Jaccard.Full.Weighted[JACCARD_MATRIX_2$Jaccard.Full.Weighted!=0 & !is.na(JACCARD_MATRIX_2$Jaccard.Full.Weighted)])/2
	JACCARD_MATRIX_2$log.Jaccard.Full.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$Jaccard.Full.Weighted)+c)
	
	
	## Make an extra column considering whether results were NA because they couldn't be subsampled.
	## This is used in the plotting. 
	JACCARD_MATRIX_2$ANYNAS <- NA
	JACCARD_MATRIX_2$ANYNAS[is.na(JACCARD_MATRIX_2$Jaccard.Full)] <- "YES"
	
	## Save the basic Jaccard Index. 
	write.table(JACCARD_MATRIX_2, paste0(path_to_output, "Summary/JACCARDMATRIX_BASIC_AllSAMPLES_UMI_RAW_", runname, ".txt"), sep='\t')
	
	widthx <- 10+(length(unique(JACCARD_MATRIX_2$Sample1))*0.05)
	heightx <- (6+(length(unique(JACCARD_MATRIX_2$Sample1))*0.05)/2)
	
	## Generate Summary Plots 
	pdf(paste0(path_to_output, "Plots/JACCARDMATRIX_BASIC_AllSAMPLES_UMI_RAW_", runname, ".pdf"), width=widthx, height=heightx)
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=SharedSeq.MeanSubsample)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1")) + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle(paste0("Sequence Overlap:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nSO")
	#f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Size.MeanSubsample)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle("Jaccard Index")+labs(fill="Subsampled:\nSequence Overlap")
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Jaccard.MeanSubsample)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle(paste0("Jaccard Index:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nJI")
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=SharedSeq.Full)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle("Sequence Overlap")+labs(fill="Full Data:\nSO")
	#i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Size.Full)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Jaccard.Full)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle("Jaccard Index")+labs(fill="Full Data:\nJI")
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.SharedSeq.MeanSubsample)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle(paste0("Sequence Overlap:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nlog2(SO+c)")
	#l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Size.MeanSubsample)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Jaccard.MeanSubsample)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')  +ggtitle(paste0("Jaccard Index:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nlog2(JI+c)")
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.SharedSeq.Full)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle("Sequence Overlap")+labs(fill="Full Data:\nlog2(SO+c)") 
	#o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Size.Full)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Jaccard.Full)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')  +ggtitle("Jaccard Index")+labs(fill="Full Data:\nlog2(JI)")
	plot(plot_grid(e, k, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(g, m, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(h, n, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(j, p, ncol=2, align="hv", axis="tblr"))
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=SharedSeq.MeanSubsample.Weighted)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1")) + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle(paste0("Weighted Sequence Overlap:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nSO")
	#f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Size.MeanSubsample.Weighted)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle("Jaccard Index")+labs(fill="Weighted Subsampled:\nSequence Overlap")
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Jaccard.MeanSubsample.Weighted)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle(paste0("Weighted Jaccard Index:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nJI")
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=SharedSeq.Full.Weighted)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle("Weighted Sequence Overlap")+labs(fill="Full Data:\nSO")
	#i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Size.Full.Weighted)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=Jaccard.Full.Weighted)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle("Weighted Jaccard Index")+labs(fill="Full Data:\nJI")
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.SharedSeq.MeanSubsample.Weighted)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle(paste0("Weighted Sequence Overlap:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nlog2(SO+c)")
	#l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Size.MeanSubsample.Weighted)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Jaccard.MeanSubsample.Weighted)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')  +ggtitle(paste0("Weighted Jaccard Index:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nlog2(JI+c)")
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.SharedSeq.Full.Weighted)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') +ggtitle("Weighted Sequence Overlap")+labs(fill="Full Data:\nlog2(SO+c)") 
	#o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Size.Full.Weighted)) + geom_tile()  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=5))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') 
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2, fill=log.Jaccard.Full.Weighted)) + geom_tile() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_gradientn( colours=c("navyblue", "darkorange1"))+ geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000')  +ggtitle("Weighted Jaccard Index")+labs(fill="Full Data:\nlog2(JI)")
	plot(plot_grid(e, k, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(g, m, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(h, n, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(j, p, ncol=2, align="hv", axis="tblr"))
	dev.off()
	print("Finished calculate_jaccard_matrix_UMI") 
	}

##-----------------------------------------------------------------------------------
## Function for calculating Jaccard Matrix and correcting for Library Hopping
	
calculate_jaccard_matrix_libhopcorrection_UMI_RAW <- function(path_to_output, runname, path_to_layout){
	path <- paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS/")
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("Att", samples, value=TRUE)
	
	mins <- c()
	## Note we look for the lowest > 100 once we calculate how many sequences greater than 1 so this is comparable across filtered/not 
	for(i in 1:length(samples)){
		a <- read.delim(samples[i], header=FALSE)
		a <- a[a$V2 >= 2,]
		a <- a$V2
		a_min <- length(a)
		mins <- c(mins, a_min)
	} 
	## Ensure that susbample must be at least 200 sequences. 
	if(any(mins>=200)){
		sample_depth <- min(mins[mins>=223])
	} else {
		sample_depth <- 223
	}
	
	subsample_depth <- floor(sample_depth*0.9)

	## Get the right samples
	path <- paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP/")
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("Barcode_filtering_information_", samples, value=TRUE)
	retry1 <- combinations(length(samples), 2, samples, repeats.allowed=TRUE)
	mins <- c()
		
	retry1 <- data.frame(retry1)
	colnames(retry1) <- c("Sample1Path", "Sample2Path")

	#Reformat Combinations to match sample ids
	## Make new Column containing just the sample name
	retry1$Sample1 <- retry1$Sample1Path
	retry1$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("BCR_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRA_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRB_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRG_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRD_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0(".txt"), "", retry1$Sample1)
	retry1$Sample2 <- retry1$Sample2Path
	retry1$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", retry1$Sample2)
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
	
	## Rename to allow for matching
	retry1$Sample1 <- retry1$Sample1Path
	retry1$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0(".txt"), "", retry1$Sample1)
	
	retry1$Sample2 <- retry1$Sample2Path
	retry1$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0(".txt"), "", retry1$Sample2)
	

	#Register doParrallel
	cl <- 20
	registerDoParallel(cl)
#dim(retry1)[1]
	## Calculate Jaccard Matrix 
	JACCARD_MATRIX <- foreach(i = 1:dim(retry1)[1], .combine=rbind) %dopar% {
		incl <- retry1[i,]
		# print i in multiples of 50(ish) to give a rough estimate of where we are in the function
		## This wont be an exact count as we are running in parrallelel therefore some nodes may be further ahead 
		if(i %% 100 == 0){
			print(i)
		}
		#return(incl)
		# Read in the sequences and replicate based on constant region counts.
		## If samples share same internal PCR barcode but are different individuals then we remove shared sequences (library hoppping correction)
		## ONLY CORRECT IF LIBRARY IS THE SAME (index hopping). 
		if(incl$Sample1 == incl$Sample2 | incl$PCRBarcode1 != incl$PCRBarcode2 | incl$Barcode1 == incl$Barcode2 | incl$BCR_Lane_S1 != incl$BCR_Lane_S2 ){
			bcr_1 <- read.delim(as.character(incl[3]), header=TRUE)		
			bcr_2 <- read.delim(as.character(incl[4]), header=TRUE)

			# Read in the sequences and replicate based on constant region counts.
			# Filter to ensure they have a minimum of 2 reads. 
			bcr_1_sequence <- bcr_1$J_barcode
			bcr_1_sequence_w <- rep(bcr_1$J_barcode, bcr_1$total_reads_with_BC)
			
			bcr_2_sequence <- bcr_2$J_barcode
			bcr_2_sequence_w <- rep(bcr_2$J_barcode, bcr_2$total_reads_with_BC)
			
		} else {
			bcr_1 <- read.delim(as.character(incl[3]), header=TRUE)		
			bcr_2 <- read.delim(as.character(incl[4]), header=TRUE)
				
			# Read in the sequences and replicate based on constant region counts.
			# Filter to ensure they have a minimum of 2 reads. 
			bcr_1_sequence <- bcr_1$J_barcode
			bcr_1_sequence_w <- rep(bcr_1$J_barcode, bcr_1$total_reads_with_BC)
			
			bcr_2_sequence <- bcr_2$J_barcode
			bcr_2_sequence_w <- rep(bcr_2$J_barcode, bcr_2$total_reads_with_BC)
			
			## Look for shared sequences
			shared_seq <- bcr_1_sequence[bcr_1_sequence %in% bcr_2_sequence]
			
			bcr_1_sequence <- bcr_1_sequence[bcr_1_sequence %notin% shared_seq]
			bcr_1_sequence_w <- bcr_1_sequence_w[bcr_1_sequence_w %notin% shared_seq]
			
			bcr_2_sequence <- bcr_2_sequence[bcr_2_sequence %notin% shared_seq]
			bcr_2_sequence_w <- bcr_2_sequence_w[bcr_2_sequence_w %notin% shared_seq]
		}
		
		## Calculate Jaccard on Weighed  and Non Weighted Repertoire (only caclulate if data is 10% larger than subsample size. 
		if(length(bcr_1_sequence) >= sample_depth   & length(bcr_2_sequence) >= sample_depth ){
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
	JACCARD_MATRIX_2$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", JACCARD_MATRIX_2$Sample2)
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
	
	c <- min(JACCARD_MATRIX_2$SharedSeq.MeanSubsample[JACCARD_MATRIX_2$SharedSeq.MeanSubsample!=0 & !is.na(JACCARD_MATRIX_2$SharedSeq.MeanSubsample)])/2
	JACCARD_MATRIX_2$log.SharedSeq.MeanSubsample <- log2(as.numeric(JACCARD_MATRIX_2$SharedSeq.MeanSubsample)+c)
	
	c <- min(JACCARD_MATRIX_2$Size.MeanSubsample[JACCARD_MATRIX_2$Size.MeanSubsample!=0 & !is.na(JACCARD_MATRIX_2$Size.MeanSubsample)])/2
	JACCARD_MATRIX_2$log.Size.MeanSubsample <- log2(as.numeric(JACCARD_MATRIX_2$Size.MeanSubsample)+c)
	
	c <- min(JACCARD_MATRIX_2$Jaccard.MeanSubsample[JACCARD_MATRIX_2$Jaccard.MeanSubsample!=0 & !is.na(JACCARD_MATRIX_2$Jaccard.MeanSubsample)])/2
	JACCARD_MATRIX_2$log.Jaccard.MeanSubsample <- log2(as.numeric(JACCARD_MATRIX_2$Jaccard.MeanSubsample)+c)
	
	c <- min(JACCARD_MATRIX_2$SharedSeq.Full[JACCARD_MATRIX_2$SharedSeq.Full!=0 & !is.na(JACCARD_MATRIX_2$SharedSeq.Full)])/2
	JACCARD_MATRIX_2$log.SharedSeq.Full <- log2(as.numeric(JACCARD_MATRIX_2$SharedSeq.Full)+c)
	
	c <- min(JACCARD_MATRIX_2$Size.Full[JACCARD_MATRIX_2$Size.Full!=0 & !is.na(JACCARD_MATRIX_2$Size.Full)])/2
	JACCARD_MATRIX_2$log.Size.Full <- log2(as.numeric(JACCARD_MATRIX_2$Size.Full)+c)
	
	c <- min(JACCARD_MATRIX_2$Jaccard.Full[JACCARD_MATRIX_2$Jaccard.Full!=0 & !is.na(JACCARD_MATRIX_2$Jaccard.Full)])/2
	JACCARD_MATRIX_2$log.Jaccard.Full <- log2(as.numeric(JACCARD_MATRIX_2$Jaccard.Full)+c)
	
	c <- min(JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted[JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted !=0 & !is.na(JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted )])/2
	JACCARD_MATRIX_2$log.SharedSeq.MeanSubsample.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$Size.MeanSubsample.Weighted[JACCARD_MATRIX_2$Size.MeanSubsample.Weighted !=0 & !is.na(JACCARD_MATRIX_2$Size.MeanSubsample.Weighted)])/2
	JACCARD_MATRIX_2$log.Size.MeanSubsample.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$Size.MeanSubsample.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted[JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted !=0 & !is.na(JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted)])/2
	JACCARD_MATRIX_2$log.Jaccard.MeanSubsample.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$SharedSeq.Full.Weighted[JACCARD_MATRIX_2$SharedSeq.Full.Weighted!=0 & !is.na(JACCARD_MATRIX_2$SharedSeq.Full.Weighted)])/2
	JACCARD_MATRIX_2$log.SharedSeq.Full.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$SharedSeq.Full.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$Size.Full.Weighted[JACCARD_MATRIX_2$Size.Full.Weighted!=0 & !is.na(JACCARD_MATRIX_2$Size.Full.Weighted)])/2
	JACCARD_MATRIX_2$log.Size.Full.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$Size.Full.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$Jaccard.Full.Weighted[JACCARD_MATRIX_2$Jaccard.Full.Weighted!=0 & !is.na(JACCARD_MATRIX_2$Jaccard.Full.Weighted)])/2
	JACCARD_MATRIX_2$log.Jaccard.Full.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$Jaccard.Full.Weighted)+c)
	
	
	## Make an extra column considering whether results were NA because they couldn't be subsampled.
	## This is used in the plotting. 
	JACCARD_MATRIX_2$ANYNAS <- NA
	JACCARD_MATRIX_2$ANYNAS[is.na(JACCARD_MATRIX_2$Jaccard.Full)] <- "YES"
	
	## Save the basic Jaccard Index. 
	write.table(JACCARD_MATRIX_2, paste0(path_to_output, "Summary/JACCARDMATRIX_BASIC_AllSAMPLES_LIBHOP_CORRECTED_UMI_RAW_", runname, ".txt"), sep='\t')
	
	widthx <- 10+(length(unique(JACCARD_MATRIX_2$Sample1))*0.05)
	heightx <- (6+(length(unique(JACCARD_MATRIX_2$Sample1))*0.05)/2)
	## Generate Summary Plots
	pdf(paste0(path_to_output, "Plots/JACCARDMATRIX_BASIC_AllSAMPLES_LIBHOP_CORRECTED_UMI_RAW_", runname, ".pdf"), width=widthx, height=heightx)
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.MeanSubsample)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle(paste0("Sequence Overlap:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nSO")
	#f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.MeanSubsample))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle("Jaccard Index")+labs(fill="Subsampled:\nSequence Overlap")
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))  +ggtitle(paste0("Jaccard Index:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nJI")
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle("Sequence Overlap")+labs(fill="Full Data:\nSO")
	#i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000')) +ggtitle("Jaccard Index")+labs(fill="Full Data:\nJI")
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle(paste0("Sequence Overlap:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nlog2(SO+c)")
	#l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle(paste0("Jaccard Index:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nlog2(JI+c)")
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))  +ggtitle("Sequence Overlap")+labs(fill="Full Data:\nlog2(SO+c)") 
	#o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.Full))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000')) +ggtitle("Jaccard Index")+labs(fill="Full Data:\nlog2(JI)")
	plot(plot_grid(e, k, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(g, m, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(h, n, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(j, p, ncol=2, align="hv", axis="tblr"))
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.MeanSubsample.Weighted)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle(paste0("Weighted Sequence Overlap:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nSO")
	#f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.MeanSubsample.Weighted))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle("Jaccard Index")+labs(fill="Weighted Subsampled:\nSequence Overlap")
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))  +ggtitle(paste0("Weighted Jaccard Index:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nJI")
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle("Weighted Sequence Overlap")+labs(fill="Full Data:\nSO")
	#i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle("Weighted Jaccard Index")+labs(fill="Full Data:\nJI")
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle(paste0("Weighted Sequence Overlap:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nlog2(SO+c)")
	#l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000')) +ggtitle(paste0("Weighted Jaccard Index:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nlog2(JI+c)")
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000')) +ggtitle("Weighted Sequence Overlap")+labs(fill="Full Data:\nlog2(SO+c)") 
	#o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.Full.Weighted))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle("Weighted Jaccard Index")+labs(fill="Full Data:\nlog2(JI)")
	plot(plot_grid(e, k, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(g, m, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(h, n, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(j, p, ncol=2, align="hv", axis="tblr"))
	dev.off()  
	print("Finished calculate_jaccard_matrix_libhopcorrection_UMI_RAW_")     
	}






##-----------------------------------------------------------------------------------------------------------------------------------	
## Function for calculating Jaccard Matrix and correcting for Library Contamination	


calculate_jaccard_matrix_libcontam_correction_UMI_RAW <- function(path_to_output, runname, path_to_layout){
	path <- paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS/")
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("Att", samples, value=TRUE)
	
	mins <- c()
	## Note we look for the lowest > 100 once we calculate how many sequences greater than 1 so this is comparable across filtered/not 
	for(i in 1:length(samples)){
		a <- read.delim(samples[i], header=FALSE)
		a <- a[a$V2 >= 2,]
		a <- a$V2
		a_min <- length(a)
		mins <- c(mins, a_min)
	} 
	## Ensure that susbample must be at least 200 sequences. 
	if(any(mins>=200)){
		sample_depth <- min(mins[mins>=223])
	} else {
		sample_depth <- 223
	}
	
	subsample_depth <- floor(sample_depth*0.9)

	## Get the right samples
	path <- paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP/")
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("Barcode_filtering_information_", samples, value=TRUE)
	retry1 <- combinations(length(samples), 2, samples, repeats.allowed=TRUE)
	mins <- c()
	
	
	retry1 <- data.frame(retry1)
	colnames(retry1) <- c("Sample1Path", "Sample2Path")

	#Reformat Combinations to match sample ids
	## Make new Column containing just the sample name
	retry1$Sample1 <- retry1$Sample1Path
	retry1$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("BCR_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRA_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRB_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRG_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0("TCRD_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0(".txt"), "", retry1$Sample1)
	retry1$Sample2 <- retry1$Sample2Path
	retry1$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", retry1$Sample2)
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
	
	## Rename to allow for matching
	retry1$Sample1 <- retry1$Sample1Path
	retry1$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", retry1$Sample1)
	retry1$Sample1 <- gsub(paste0(".txt"), "", retry1$Sample1)
	
	retry1$Sample2 <- retry1$Sample2Path
	retry1$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", retry1$Sample2)
	retry1$Sample2 <- gsub(paste0(".txt"), "", retry1$Sample2)

	#Register doParrallel
	cl <- 20
	registerDoParallel(cl)
	#dim(retry1)[1]
	## Calculate Jaccard Matrix 
	JACCARD_MATRIX <- foreach(i = 1:dim(retry1)[1], .combine=rbind) %dopar% {
		incl <- retry1[i,]
		# print i in multiples of 50(ish) to give a rough estimate of where we are in the function
		## This wont be an exact count as we are running in parrallelel therefore some nodes may be further ahead 
		if(i %% 100 == 0){
			print(i)
		}
		#return(incl)
		# Read in the sequences and replicate based on constant region counts.
		## If samples share same internal PCR barcode but are different individuals then we remove shared sequences (library hoppping correction)
		## Correct irrespective of sequencing lane (e.g. libhop and library contamination) 
		if(incl$Sample1 == incl$Sample2 | incl$PCRBarcode1 != incl$PCRBarcode2 | incl$Barcode1 == incl$Barcode2){
			bcr_1 <- read.delim(as.character(incl[3]), header=TRUE)		
			bcr_2 <- read.delim(as.character(incl[4]), header=TRUE)

			# Read in the sequences and replicate based on constant region counts.
			# Filter to ensure they have a minimum of 2 reads. 
			bcr_1_sequence <- bcr_1$J_barcode
			bcr_1_sequence_w <- rep(bcr_1$J_barcode, bcr_1$total_reads_with_BC)
			
			bcr_2_sequence <- bcr_2$J_barcode
			bcr_2_sequence_w <- rep(bcr_2$J_barcode, bcr_2$total_reads_with_BC)
			
		} else {
			bcr_1 <- read.delim(as.character(incl[3]), header=TRUE)		
			bcr_2 <- read.delim(as.character(incl[4]), header=TRUE)
				
			# Read in the sequences and replicate based on constant region counts.
			# Filter to ensure they have a minimum of 2 reads. 
			bcr_1_sequence <- bcr_1$J_barcode
			bcr_1_sequence_w <- rep(bcr_1$J_barcode, bcr_1$total_reads_with_BC)
			
			bcr_2_sequence <- bcr_2$J_barcode
			bcr_2_sequence_w <- rep(bcr_2$J_barcode, bcr_2$total_reads_with_BC)
			
			## Look for shared sequences
			shared_seq <- bcr_1_sequence[bcr_1_sequence %in% bcr_2_sequence]
			
			bcr_1_sequence <- bcr_1_sequence[bcr_1_sequence %notin% shared_seq]
			bcr_1_sequence_w <- bcr_1_sequence_w[bcr_1_sequence_w %notin% shared_seq]
			
			bcr_2_sequence <- bcr_2_sequence[bcr_2_sequence %notin% shared_seq]
			bcr_2_sequence_w <- bcr_2_sequence_w[bcr_2_sequence_w %notin% shared_seq]
		}
		
		## Calculate Jaccard on Weighed  and Non Weighted Repertoire (only caclulate if data is 10% larger than subsample size. 
		if(length(bcr_1_sequence) >= sample_depth   & length(bcr_2_sequence) >= sample_depth){
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
	JACCARD_MATRIX_2$Sample1 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub(paste0(path_to_output, "ORIENTATED_SEQUENCES/TMP//Barcode_filtering_information_"), "", JACCARD_MATRIX_2$Sample2)
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
	
	c <- min(JACCARD_MATRIX_2$SharedSeq.MeanSubsample[JACCARD_MATRIX_2$SharedSeq.MeanSubsample!=0 & !is.na(JACCARD_MATRIX_2$SharedSeq.MeanSubsample)])/2
	JACCARD_MATRIX_2$log.SharedSeq.MeanSubsample <- log2(as.numeric(JACCARD_MATRIX_2$SharedSeq.MeanSubsample)+c)
	
	c <- min(JACCARD_MATRIX_2$Size.MeanSubsample[JACCARD_MATRIX_2$Size.MeanSubsample!=0 & !is.na(JACCARD_MATRIX_2$Size.MeanSubsample)])/2
	JACCARD_MATRIX_2$log.Size.MeanSubsample <- log2(as.numeric(JACCARD_MATRIX_2$Size.MeanSubsample)+c)
	
	c <- min(JACCARD_MATRIX_2$Jaccard.MeanSubsample[JACCARD_MATRIX_2$Jaccard.MeanSubsample!=0 & !is.na(JACCARD_MATRIX_2$Jaccard.MeanSubsample)])/2
	JACCARD_MATRIX_2$log.Jaccard.MeanSubsample <- log2(as.numeric(JACCARD_MATRIX_2$Jaccard.MeanSubsample)+c)
	
	c <- min(JACCARD_MATRIX_2$SharedSeq.Full[JACCARD_MATRIX_2$SharedSeq.Full!=0 & !is.na(JACCARD_MATRIX_2$SharedSeq.Full)])/2
	JACCARD_MATRIX_2$log.SharedSeq.Full <- log2(as.numeric(JACCARD_MATRIX_2$SharedSeq.Full)+c)
	
	c <- min(JACCARD_MATRIX_2$Size.Full[JACCARD_MATRIX_2$Size.Full!=0 & !is.na(JACCARD_MATRIX_2$Size.Full)])/2
	JACCARD_MATRIX_2$log.Size.Full <- log2(as.numeric(JACCARD_MATRIX_2$Size.Full)+c)
	
	c <- min(JACCARD_MATRIX_2$Jaccard.Full[JACCARD_MATRIX_2$Jaccard.Full!=0 & !is.na(JACCARD_MATRIX_2$Jaccard.Full)])/2
	JACCARD_MATRIX_2$log.Jaccard.Full <- log2(as.numeric(JACCARD_MATRIX_2$Jaccard.Full)+c)
	
	c <- min(JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted[JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted !=0 & !is.na(JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted )])/2
	JACCARD_MATRIX_2$log.SharedSeq.MeanSubsample.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$SharedSeq.MeanSubsample.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$Size.MeanSubsample.Weighted[JACCARD_MATRIX_2$Size.MeanSubsample.Weighted !=0 & !is.na(JACCARD_MATRIX_2$Size.MeanSubsample.Weighted)])/2
	JACCARD_MATRIX_2$log.Size.MeanSubsample.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$Size.MeanSubsample.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted[JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted !=0 & !is.na(JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted)])/2
	JACCARD_MATRIX_2$log.Jaccard.MeanSubsample.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$Jaccard.MeanSubsample.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$SharedSeq.Full.Weighted[JACCARD_MATRIX_2$SharedSeq.Full.Weighted!=0 & !is.na(JACCARD_MATRIX_2$SharedSeq.Full.Weighted)])/2
	JACCARD_MATRIX_2$log.SharedSeq.Full.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$SharedSeq.Full.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$Size.Full.Weighted[JACCARD_MATRIX_2$Size.Full.Weighted!=0 & !is.na(JACCARD_MATRIX_2$Size.Full.Weighted)])/2
	JACCARD_MATRIX_2$log.Size.Full.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$Size.Full.Weighted)+c)
	
	c <- min(JACCARD_MATRIX_2$Jaccard.Full.Weighted[JACCARD_MATRIX_2$Jaccard.Full.Weighted!=0 & !is.na(JACCARD_MATRIX_2$Jaccard.Full.Weighted)])/2
	JACCARD_MATRIX_2$log.Jaccard.Full.Weighted <- log2(as.numeric(JACCARD_MATRIX_2$Jaccard.Full.Weighted)+c)
	
	## Make an extra column considering whether results were NA because they couldn't be subsampled.
	## This is used in the plotting. 
	JACCARD_MATRIX_2$ANYNAS <- NA
	JACCARD_MATRIX_2$ANYNAS[is.na(JACCARD_MATRIX_2$Jaccard.Full)] <- "YES"
	
	## Save the basic Jaccard Index. 
	write.table(JACCARD_MATRIX_2, paste0(path_to_output, "Summary/JACCARDMATRIX_BASIC_AllSAMPLES_LIBCONTAM_CORRECTED_UMI_RAW_", runname, ".txt"), sep='\t')
	
	widthx <- 10+(length(unique(JACCARD_MATRIX_2$Sample1))*0.05)
	heightx <- (6+(length(unique(JACCARD_MATRIX_2$Sample1))*0.05)/2)
	
	## Generate Summary Plots
	pdf(paste0(path_to_output, "Plots/JACCARDMATRIX_BASIC_AllSAMPLES_LIBCONTAM_CORRECTED_UMI_RAW_", runname, ".pdf"), width=widthx, height=heightx)
e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.MeanSubsample)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle(paste0("Sequence Overlap:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nSO")
	#f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.MeanSubsample))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle("Jaccard Index")+labs(fill="Subsampled:\nSequence Overlap")
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000')) +ggtitle(paste0("Jaccard Index:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nJI")
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle("Sequence Overlap")+labs(fill="Full Data:\nSO")
	#i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle("Jaccard Index")+labs(fill="Full Data:\nJI")
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle(paste0("Sequence Overlap:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nlog2(SO+c)")
	#l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.MeanSubsample))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle(paste0("Jaccard Index:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nlog2(JI+c)")
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000')) +ggtitle("Sequence Overlap")+labs(fill="Full Data:\nlog2(SO+c)") 
	#o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.Full))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.Full))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000')) +ggtitle("Jaccard Index")+labs(fill="Full Data:\nlog2(JI)")
	plot(plot_grid(e, k, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(g, m, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(h, n, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(j, p, ncol=2, align="hv", axis="tblr"))
	e <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.MeanSubsample.Weighted)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle(paste0("Weighted Sequence Overlap:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nSO")
	#f <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.MeanSubsample.Weighted))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000')) +ggtitle("Jaccard Index")+labs(fill="Weighted Subsampled:\nSequence Overlap")
	g <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle(paste0("Weighted Jaccard Index:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nJI") 
	h <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=SharedSeq.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle("Weighted Sequence Overlap")+labs(fill="Full Data:\nSO")
	#i <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Size.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))
	j <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=Jaccard.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle("Weighted Jaccard Index")+labs(fill="Full Data:\nJI")
	k <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000')) +ggtitle(paste0("Weighted Sequence Overlap:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nlog2(SO+c)")
	#l <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))
	m <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.MeanSubsample.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))+ggtitle(paste0("Weighted Jaccard Index:\nSub.Depth: ", subsample_depth))+labs(fill="Subsampled:\nlog2(JI+c)")
	n <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.SharedSeq.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000')) +ggtitle("Weighted Sequence Overlap")+labs(fill="Full Data:\nlog2(SO+c)") 
	#o <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Size.Full.Weighted))+  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))
	p <- ggplot(JACCARD_MATRIX_2, aes(Sample1, Sample2)) + geom_tile(aes(fill=log.Jaccard.Full.Weighted))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_fill_gradientn( colours=c("navyblue", "darkorange1"))  + geom_tile(aes(width = 1, height = 1), data = JACCARD_MATRIX_2[JACCARD_MATRIX_2$ANYNAS=="YES",], fill = "white", color='#00000000') + geom_tile(aes(color=factor(LibCorrected, c("YES", "NO"))), fill = '#00000000', size = 0.2) + scale_color_manual(name = "Sequence\nCorrected", values = c("green", '#00000000'))  +ggtitle("Weighted Jaccard Index")+labs(fill="Full Data:\nlog2(JI)")

	plot(plot_grid(e, k, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(g, m, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(h, n, ncol=2, align="hv", axis="tblr"))
	plot(plot_grid(j, p, ncol=2, align="hv", axis="tblr"))
	dev.off()
	print("Finished calculate_jaccard_matrix_libcontam_correction_UMI_RAW")
	}


