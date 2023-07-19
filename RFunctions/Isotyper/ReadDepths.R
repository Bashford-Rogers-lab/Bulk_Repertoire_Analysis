suppressMessages(library(corrplot))
suppressMessages(library(data.table))
suppressMessages(library(doParallel))
suppressMessages(library(dplyr))
suppressMessages(library(foreach))
suppressMessages(library(ggforce))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrastr))
suppressMessages(library(gridExtra))
suppressMessages(library(gtools))
#suppressMessages(library(Gviz))
suppressMessages(library(Hmisc))
suppressMessages(library(matrixStats))
suppressMessages(library(moments))
suppressMessages(library(mousetrap))
suppressMessages(library(Peptides))
suppressMessages(library(plot3D))
suppressMessages(library(plyr))
suppressMessages(library(purrr))
suppressMessages(library(reshape2))
#suppressMessages(library(ShortRead))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))



getdepths <- function(path=path, iso_type=iso_type){
	samples <- list.files(path, full.name=TRUE)
	samples <- grep("_Summary", samples, value=TRUE)
	productive  <- grep("_productive", samples, value=TRUE)
	noproductive <- grep("_unproductive", samples, value=TRUE)
	if(iso_type == "UNPRODUCTIVE"){
		samples <- grep("_unproductive", samples, value=TRUE)
	}
	if(iso_type=="PRODUCTIVE"){
		samples <- grep("_productive", samples, value=TRUE)
	}
	if(iso_type=="ALL"){
		samples <- samples[!samples %in% c(productive, noproductive)]
	}
	# Set to run in parrallell to speed this up!!
	if(iso_type == "UNPRODUCTIVE" | iso_type=="PRODUCTIVE"){
		registerDoParallel(8)
		order_samples <- foreach(i = 1:length(samples), .combine=rbind, .packages='tidyverse') %dopar% {
			a <- read.delim(samples[i], header=FALSE)
			a <- a$V2
			a_min <- length(a)
			sampleid <- samples[i]
			read_depths <- c(sampleid, a_min)
			return(read_depths)
		}
	} 
	if(iso_type=="ALL"){
		registerDoParallel(8)
		order_samples <- foreach(i = 1:length(samples), .combine=rbind, .packages='tidyverse') %dopar% {
			a <- read.delim(samples[i], header=FALSE)
			a <- a[(a$V3=="productive (see comment)" | a$V3=="productive" | a$V3=="unproductive (see comment)" | a$V3=="unproductive"), ]
			a <- a$V2
			a_min <- length(a)
			sampleid <- samples[i]
			read_depths <- c(sampleid, a_min)
			return(read_depths)
		}
	} 
	
	## Renaming Read Depths to fit with sample names 	
	depths <- data.frame(order_samples)
	colnames(depths) <- c("order_samples", "ReadDepth")
	depths$order_samples <- gsub(paste0(path, "/IMGT_BCR_"), "", depths$order_samples)
	depths$order_samples <- gsub(paste0(path, "/IMGT_TCRA_"), "", depths$order_samples)
	depths$order_samples <- gsub(paste0(path, "/IMGT_TCRB_"), "", depths$order_samples)
	depths$order_samples <- gsub(paste0(path, "/IMGT_TCRG_"), "", depths$order_samples)
	depths$order_samples <- gsub(paste0(path, "/IMGT_TCRD_"), "", depths$order_samples)
	depths$order_samples <- gsub(paste0(path), "", depths$order_samples)
	depths$order_samples <- gsub("/IMGT_", "", depths$order_samples)
	depths$order_samples <- gsub(".txt", "", depths$order_samples)
	depths$order_samples <- gsub("_1_Summary", "", depths$order_samples)		
	colnames(depths) <- c("SampleIDforDepths", "ReadDepth")
	depths$ReadDepth <- as.character(depths$ReadDepth)
	depths$ReadDepth <- as.numeric(depths$ReadDepth)

	return(depths)
	print("Calculated Read Depths")
}
