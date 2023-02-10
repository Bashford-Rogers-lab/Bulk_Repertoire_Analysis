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
##-----------------------------------------------------------------------------------