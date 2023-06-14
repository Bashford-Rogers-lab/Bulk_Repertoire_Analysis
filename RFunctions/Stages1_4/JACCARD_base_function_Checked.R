suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggforce))
#suppressMessages(library(Gviz))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(optparse))
suppressMessages(library(gtools))


##-----------------------------------------------------------------------------------
`%notin%` <<- Negate(`%in%`)	
cal_jaccard <<- function(vac1, vec2){
	intersection = length(intersect(vac1, vec2))
	union = length(vac1) + length(vec2) - intersection
	jaccard = (intersection/union)
	return(jaccard)
}


calculate_jaccard <<- function(seq_1, seq_2, depth){
	jaccard <- c()
	shared_seq_count <- c()
	sample_count  <- c()
	for(i in 1:1000){
		bcr_1_sample <- sample(seq_1, depth, replace=FALSE)
		bcr_2_sample <- sample(seq_2, depth, replace=FALSE)
		jaccard_subset <- cal_jaccard(bcr_1_sample, bcr_2_sample)
		jaccard <- c(jaccard, jaccard_subset)
		shared_seq_count <- c(shared_seq_count, length(intersect(bcr_1_sample, bcr_2_sample))) 
		sample_count <- c(sample_count, (depth+depth))
	} 
	mean_jaccard_subsampled <- mean(jaccard)
	mean_no_shared_sequences <- mean(shared_seq_count)
	mean_sample_count <- mean(sample_count)
	## nosubsampling
	jaccard_full <- cal_jaccard(seq_1, seq_2)
	
	## Intersection
    shared_seq_count_full <- length(intersect(seq_1, seq_2))
	## Union
	sample_count_full <- length(seq_1) + length(seq_2) - shared_seq_count_full
	
	results <- c(mean_no_shared_sequences,mean_sample_count, mean_jaccard_subsampled, shared_seq_count_full, sample_count_full, jaccard_full)
	return(results)
	}

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------