## Function for calculating and plotting the abundance against sequencing depth for BCR/TCR analysis 
library(plyr)
library(ggrepel)
#path_to_output <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/'

calculate_rarefaction_neat <- function(path_to_output, chain){
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
		
		## Lets get the depths to use 
		sample_depth <- c(100)
		
		total_sequences <- length(bcr)
		intervals_20 <- total_sequences / 20 
		intervals_20 <- floor(intervals_20)
			
		sample_depth <- c()
		mean_unique_seq <- c()
		sd_unique_seq <- c()
		var_unique_seq <- c()
		mean_unique_prop <- c()
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
			prop2 <- c()
			
	        for(i in 1:100){
				no_unique <- sample(bcr, depth, replace=FALSE)
				no_unique <- length(unique(no_unique))
				unique_seq <- c(unique_seq, no_unique)
				prop <- no_unique/depth
				prop2 <- c(prop, prop2)
			}
			unique_seq_mean <- mean(unique_seq)
			unique_seq_var <- var(unique_seq)
			unique_seq_sd <- sd(unique_seq)
			unique_seq_prop <- mean( prop2)
			
			sample_depth <- c(sample_depth, depth)
			mean_unique_seq <- c(mean_unique_seq, unique_seq_mean)
			sd_unique_seq <- c(sd_unique_seq, unique_seq_sd)
			var_unique_seq <- c(var_unique_seq, unique_seq_var)
			mean_unique_prop <- c(mean_unique_prop, unique_seq_prop)
		}
		
		result <- data.frame(cbind(sample_depth, mean_unique_seq, sd_unique_seq,  var_unique_seq, mean_unique_prop))
		
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
			
		return(result)
		}
		
		
				  
		ABUNDANCE_MATRIX$label <-ABUNDANCE_MATRIX$sample
		ABUNDANCE_MATRIX$label[!ABUNDANCE_MATRIX$sample %like% "UK17730209" & !ABUNDANCE_MATRIX$sample %like% "UK59580028"] <- NA
		
			
		ABUNDANCE_MATRIX$DISEASE="SEPSIS"
		ABUNDANCE_MATRIX$DISEASE[ABUNDANCE_MATRIX$sample %like% "HV"] <- "HEALTH"
		ABUNDANCE_MATRIX$DISEASE[ABUNDANCE_MATRIX$sample %like% "JR1795_1003_POSITIVE"] <- "TECHNICAL"
		labs <- ABUNDANCE_MATRIX[!is.na(ABUNDANCE_MATRIX$label),]
		labs <- labs %>% group_by(label) %>% top_n(1, sample_depth)
  
  
		# Plot results 
		pdf(paste0(path_to_output, "THESIS_PLOTS/Abundance_curve_neat.pdf"), width=10, height=7.5)
		p1 <- ggplot(ABUNDANCE_MATRIX, aes(x=sample_depth, y=mean_unique_seq, group=sample, color=DISEASE))   + theme_classic()  + geom_line(aes(colour=DISEASE), alpha=0.5)+ xlab("Subsample Depth") + ylab("Mean number\nof unique VDJ sequences\n(repeat=100)")+facet_wrap(~DISEASE)+ geom_abline(slope=1, intercept=0, colour="black", linetype=4)+xlim(0, round_any(max(ABUNDANCE_MATRIX$sample_depth), 10000, f = ceiling))+ylim(0, round_any(max(ABUNDANCE_MATRIX$sample_depth), 10000, f = ceiling))+ geom_label_repel(data=labs, aes(label = label),  na.rm = TRUE, box.padding = 0.25, max.overlaps = Inf, col="black", size=2)+ggtitle(chain)+ggtitle(chain)+guides(fill="none", colour="none")
		p2 <- ggplot(ABUNDANCE_MATRIX, aes(x=sample_depth, y=mean_unique_prop, group=sample, color=DISEASE))  + theme_classic()  + geom_line(aes(colour=DISEASE), alpha=0.5)+ xlab("Subsample Depth") + ylab("Mean proportion\nof unique VDJ sequences\n(repeat=100)")+facet_wrap(~DISEASE)+xlim(0, round_any(max(ABUNDANCE_MATRIX$sample_depth), 10000, f = ceiling))+ geom_label_repel(data=labs, aes(label = label),  na.rm = TRUE, box.padding = 0.25, max.overlaps = Inf, col="black", size=2)+ggtitle(chain)+guides(fill="none", colour="none")
        plot(plot_grid(p1,p2, ncol=1, labels="AUTO", align="hv", axis="tblr"))
		dev.off()

		# Save results
		write.table(ABUNDANCE_MATRIX, paste0(path_to_output, "Summary/Abundance_curve.txt"), sep='\t')
		px <- plot_grid(p1,p2, ncol=2, align="hv", axis="tblr")
		return(px)

}

