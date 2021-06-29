## Function for calculating and plotting the abundance against sequencing depth for BCR/TCR analysis 

calculate_rarefaction <- function(path_to_output){
	path <- paste0(path_to_output, "ORIENTATED_SEQUENCES/NETWORKS/")
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("Att", samples, value=TRUE)
	
	#Register doParrallel (10 nodes seems to be the maximum you can run on the cluster with 1 slot or it crashes!
	cl <- 7
	registerDoParallel(cl)
	
	## Calculate Abundancy Matrix 
	ABUNDANCE_MATRIX <- foreach(i = 1:length(samples), .combine=rbind) %dopar% {
		sample_id <- samples[i]
		bcr <- read.delim(sample_id, header=FALSE)
		bcr <- rep(bcr$V3, bcr$V2)
		
		sample_depth <- c(100)
		number_unique <- c(0)
		
		total_sequences <- length(bcr)
		intervals_20 <- total_sequences / 20 
		intervals_20 <- floor(intervals_20)
			
		sample_depth <- c()
		mean_unique_seq <- c()
		sd_unique_seq <- c()
		var_unique_seq <- c()
		
		# If the number of sequences for the sample is less than 20 we need to adjust the limits ands and interval number
		if(intervals_20<1){
			limit <- total_sequences
			intervals_20 <- 1 
		} else {
			limit <- 20
		}
		
		for(i in 1:limit){
			depth <- i*intervals_20
			unique_seq <- c()
			
	        for(i in 1:50){
				no_unique <- sample(bcr, depth, replace=FALSE)
				no_unique <- length(unique(no_unique))
				unique_seq <- c(unique_seq, no_unique)
			}
			unique_seq_mean <- mean(unique_seq)
			unique_seq_var <- var(unique_seq)
			unique_seq_sd <- sd(unique_seq)

			sample_depth <- c(sample_depth, depth)
			mean_unique_seq <- c(mean_unique_seq, unique_seq_mean)
			sd_unique_seq <- c(sd_unique_seq, unique_seq_sd)
			var_unique_seq <- c(var_unique_seq, unique_seq_var)
		}
		
		result <- data.frame(cbind(sample_depth, mean_unique_seq, sd_unique_seq,  var_unique_seq))
		
		if(intervals_20 !=1){
			result <- rbind(c(1, 1, 1, 1), result)
		}
		
		sample_id <- gsub(path, "", sample_id)
		sample_id <- gsub("/Att_", "", sample_id)
		sample_id <- gsub(".txt", "", sample_id)
		sample_id <- gsub("BCR_", "", sample_id)
	    sample_id <- gsub("TCRB_", "", sample_id)
		sample_id <- gsub("TCRA_", "", sample_id)
		sample_id <- gsub("TCRG_", "", sample_id)
		sample_id <- gsub("TCRD_", "", sample_id)
		sample_id <- gsub("TR_", "", sample_id)
		result$sample <- sample_id
			
		## Calculate Jaccard on Weighed  and Non Weighted Repertoire
		return(result)
		}

		# Plot results 
		pdf(paste0(path_to_output, "Plots/Abundance_curve.pdf"), width=10, height=5)
		p1 <- ggplot(ABUNDANCE_MATRIX, aes(x=sample_depth, y=mean_unique_seq, group=sample, color=sample)) + geom_point(size=0.5)  + theme_bw() + guides(colour=FALSE) + geom_line(aes(colour=sample))+ xlab("Subsample Depth") + ylab("Mean number of unique VDJ sequences (n=50)") +ggtitle("VDJ Abundance in relation to Read Depth")
		p2 <- ggplot(ABUNDANCE_MATRIX, aes(x=log(sample_depth), y=log(mean_unique_seq), group=sample, color=sample)) + geom_point(size=0.5)  + theme_bw() + guides(colour=FALSE) + geom_line(aes(colour=sample))+ xlab("log Subsample Depth") + ylab("log Mean number of unique VDJ sequences (n=50)") +ggtitle("VDJ Abundance in relation to Read Depth")
        plot(plot_grid(p1, p2, ncol=2))
		dev.off()
		
		# Save results
		write.table(ABUNDANCE_MATRIX, paste0(path_to_output, "Summary/Abundance_curve.txt"), sep='\t')

}

