#------------------------------------------------------------------------------------------------
## FUNCTION FOR PLOTTING THE RESULTS OF RBR ISOTYPER SCRIPT
## ASSESSING CORRELATION WITH READ DEPTH AND MISSINGNESS OF DATA 
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk

# Recquired Packages
library(reshape2)
library(ggplot2)
library(Hmisc)
library(corrplot)
library(stringr)
library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(data.table)
library(foreach)
library(doParallel) 
library(ggforce)
library(plot3D)
library(Peptides)
library(plyr)
library(moments)
library(mousetrap)

## Function 
summary_isotyper <- function(outputdir, samplesfilepost, iso_type){

#Getting SampleIDs
ids1 <- read.delim(samplesfilepost, sep='\t', header=FALSE)
ids <- as.character(ids1$V1)
chain_vdj <- unique(as.character(ids1$V4))
ids_all <- as.character(ids)
if(iso_type == "UNPRODUCTIVE"){
	ids_all <- paste0(ids_all, "_unproductive")
}
if(iso_type=="PRODUCTIVE"){
	ids_all <- paste0(ids_all, "_productive")
}


## Getting Read Depths for Each Sample based on RBR pipeline: 
path <- paste0(outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT")
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
	registerDoParallel(15)
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
	registerDoParallel(15)
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
## SampleReaddepth Dataframe 
read_depths_all <- depths

write.table(read_depths_all, paste0(outputdir, "Summary/Read_Depths_", iso_type, ".txt"), sep="\t", row.names=TRUE)

## Getting subsample depths which were used for isotyper script 
counts_used <- paste0(outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS")
all_files <- list.files(counts_used, full.name=TRUE)
all_files <- grep("depth_per_isotype", all_files, value=TRUE)
counts_used <- read.delim(all_files[1], sep="\t", header=TRUE)
counts_used <- counts_used[counts_used$type=="UNIQ",]
subsampled_depth_all <- counts_used$min[counts_used$X.isotype=="all"]
## Extract receptor information (e.g. chain usage if TCR
if(chain_vdj=="TCR" | chain_vdj=="TR"){
	counts_try <- counts_used[(counts_used$X.isotype  != "ALL" & counts_used$X.isotype  != "all"),]
	receptor <- counts_try$X.isotype[counts_try$min==max(counts_try$min)]
}


## Plot read depths
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	IGHA2 = counts_used$min[counts_used$X.isotype=="IGHA2"]
	IGHA1 = counts_used$min[counts_used$X.isotype=="IGHA1"]
	IGHGP = counts_used$min[counts_used$X.isotype=="IGHGP"]
	IGHG1 = counts_used$min[counts_used$X.isotype=="IGHG1"]
	IGHG2 = counts_used$min[counts_used$X.isotype=="IGHG2"]
	IGHG3 = counts_used$min[counts_used$X.isotype=="IGHG3"]
	IGHG4 = counts_used$min[counts_used$X.isotype=="IGHG4"]
	IGHE = counts_used$min[counts_used$X.isotype=="IGHE"]
	IGHD = counts_used$min[counts_used$X.isotype=="IGHD"]
	IGHEP2 = counts_used$min[counts_used$X.isotype=="IGHEP2"]
	IGHM = counts_used$min[counts_used$X.isotype=="IGHM"]
	
	pdf(paste0(outputdir, "Plots/Read_Depths_", iso_type, ".pdf"), width=5, height=5)
	s <- ggplot(data=read_depths_all, aes(ReadDepth)) + geom_histogram() + theme_bw() + geom_vline(aes(xintercept=subsampled_depth_all, col="ALL"), show.legend=TRUE) 
	s <- s +  geom_vline(aes(xintercept=IGHA2, col="IGHA2"), show.legend=TRUE) 
	s <- s + geom_vline(aes(xintercept=IGHA1, col="IGHA1"), show.legend=TRUE, linetype = "longdash") 
	s <- s + geom_vline(aes(xintercept=IGHGP, col="IGHGP"), show.legend=TRUE)
	s <- s + geom_vline(aes(xintercept=IGHG1, col="IGHG1"), show.legend=TRUE) 
	s <- s + geom_vline(aes(xintercept=IGHG2, col="IGHG2"), show.legend=TRUE, linetype = "longdash") 
	s <- s + geom_vline(aes(xintercept=IGHG3, col="IGHG3"), show.legend=TRUE, linetype = "longdash") 
	s <- s + geom_vline(aes(xintercept=IGHG4, col="IGHG4"), show.legend=TRUE, linetype = "dotted") 
	s <- s + geom_vline(aes(xintercept=IGHE, col="IGHE"), show.legend=TRUE) 	
	s <- s + geom_vline(aes(xintercept=IGHM, col="IGHM"), show.legend=TRUE, linetype = "longdash") 
	s <- s + geom_vline(aes(xintercept=IGHD, col="IGHD"), show.legend=TRUE, linetype = "dotted") 
	s <- s + geom_vline(aes(xintercept=IGHEP2, col="IGHEP2"), show.legend=TRUE, linetype = "twodash") 
	s <- s + scale_color_manual(name = "Thresholds", values = c(ALL="red", IGHA2="BLUE", IGHA1="ORANGE", IGHGP="DARKGREEN", IGHG1="BROWN4", IGHG2="YELLOW", IGHG3="GREY", IGHG4="PURPLE", IGHE="BLACK", IGHM="LIGHTBLUE", IGHEP2="LIGHTGREEN",IGHD="DEEPPINK"))
	s <- s + facet_zoom(xlim=c(0, (subsampled_depth_all+50)))
	plot(s)
	dev.off()
} else {
	pdf(paste0(outputdir, "Plots/Read_Depths_", iso_type, ".pdf"), width=10, height=5)
	s <- ggplot(data=read_depths_all, aes(ReadDepth)) + geom_histogram() + theme_bw() + geom_vline(aes(xintercept=subsampled_depth_all, col="ALL"), show.legend=TRUE) 
	s <- s + scale_color_manual(name = "Thresholds", values = c(ALL="red"))
	plot(s)
	dev.off()
}

if(iso_type == "UNPRODUCTIVE"){ 
	subsampled_depth_allx <- subsampled_depth_all/10
} else {
	subsampled_depth_allx <- subsampled_depth_all
} 


# Low read depth samples:
# Probably want to exlude them
samples_to_low_all <- read_depths_all$SampleIDforDepths[read_depths_all$ReadDepth < subsampled_depth_allx]


## Read Depths of excluded samples :D
read_exclude <- read_depths_all[read_depths_all$SampleIDforDepths %in% samples_to_low_all, ]


if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	pdf(paste0(outputdir, "Plots/Read_Depths_excluded_samples_", iso_type, ".pdf"), width=5, height=5)
	s <- ggplot(data=read_exclude, aes(ReadDepth)) + geom_histogram() + theme_bw() + geom_vline(aes(xintercept=subsampled_depth_allx, col="ALL"), show.legend=TRUE) 
	s <- s + scale_color_manual(name = "Thresholds", values = c(ALL="red"))
	plot(s)
	dev.off()
} else {
	pdf(paste0(outputdir, "Plots/Read_Depths_excluded_samples_", iso_type, ".pdf"), width=10, height=5)
	s <- ggplot(data=read_exclude, aes(ReadDepth)) + geom_histogram() + theme_bw() + geom_vline(aes(xintercept=subsampled_depth_allx, col="ALL"), show.legend=TRUE) 
	s <- s + scale_color_manual(name = "Thresholds", values = c(ALL="red"))
	plot(s)
	dev.off()
}


		
##Begining Data compiling:
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
## File Number 1: 
## Summary of Clustering metrics (Gini, Renyi etc) 
## Works for BCR 
## Works for TCR 
## NO NEED TO FILTER 

file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check1 <- FALSE
} else {
	check1 <- TRUE
}
info = file.info(file)
if(info$size != 0) {
	p <- as.matrix(read.delim(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.Id"]) %in% ids_all),]
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
	#######################
	p <- data.frame(p)
	## Want to make the relevant columns numeric
	p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
	## Convert '-1' to NA 
	if(any(p[!is.na(p)]==-1)){
		p[p==-1] <- NA
	}
	######################
	p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"Isotype"]))),]
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
	## Can't have different datatypes in dataframe
	p <- data.frame(p)
	p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)

	id = as.character(p[,"X.Id"])
	ids = sort(unique(id))
	class = as.character(p[,"Isotype"])
	classes = sort(unique(class))
	reads = as.numeric(p[,"N.reads"])
	vgini = as.numeric(p[,"Vertex.Gini.Index"])
	N.reads = as.numeric(p[,"N.reads"])
	vertices = as.numeric(p[,"N.vertices"])
	N_clusters = as.numeric(p[,"N_clusters"])
	cgini =  as.numeric(p[,"Cluster.Gini.Index"])
	mean_vertex_size =  as.numeric(p[,"mean_vertex_size"])
	max_clust_size =  as.numeric(p[,"max_clust_pop"])
	max_vertex_size =  as.numeric(p[,"max_vertex_pop"])
	## Extras due to update
	D.5 <- as.numeric(p[,"D5"])
	D.10 <- as.numeric(p[,"D10"])
	D.50 <- as.numeric(p[,"D50"])
	#Calculating Normalised V/C Renyi Scores 
	vrenyi = 1-(as.numeric(p[,"Vertex.Reyni"])/log(as.numeric(p[,"N.vertices"])))
	crenyi =  1-(as.numeric(p[,"Cluster_Renyi"])/log(as.numeric(p[,"N_clusters"])))
	### #### I THINK THIS SHOULD BE 1 AND NOT -1 
	a <- -0.001
	vrenyi[which(vrenyi<0 & vrenyi >=a)] <- 0
	crenyi[which(crenyi<0 & crenyi >=a)] <- 0
	vrenyi[which(vrenyi< a)] <- NA
	crenyi[which(crenyi<a)] <- NA
	
	# Fill Data frame 
	class_reads = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	v_gini = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	c_gini = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	c_renyi = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	v_renyi = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	m_mean_vertex_size = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	m_max_clust_size = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	m_max_vertex_size = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	D5 = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	D10 = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	D50 = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	#making matrix 
	for(i in c(1:length(ids))){
		for (c in c(1:length(classes))){
			w = intersect(which(id==ids[i]), which(class==classes[c]))
			if(length(w)>=1){
				class_reads[i,c] = mean(reads[w])
				v_gini[i,c] = mean(vgini[w])
				c_gini[i,c] = mean(cgini[w])
				m_mean_vertex_size[i,c] = mean(mean_vertex_size[w])
				m_max_clust_size[i,c] = mean(max_clust_size[w])
				m_max_vertex_size[i,c] = mean(max_vertex_size[w])
				c_renyi[i,c] = mean(crenyi[w])
				v_renyi[i,c] = mean(vrenyi[w])
				D5[i, c] =mean(D.5[w])
				D10[i, c] =mean(D.10[w])
				D50[i, c] =mean(D.50[w])
			}
		}
	}
	q <- list(v_gini, c_gini, m_mean_vertex_size, m_max_clust_size, m_max_vertex_size, c_renyi, v_renyi, D5, D10, D50)
	names(q) <- c("v_gini", "c_gini", "m_mean_vertex_size", "m_max_clust_size", "m_max_vertex_size", "c_renyi_normalised", "v_renyi_normalised", "D5_subsampled", "D10_subsampled", "D50_subsampled")
	#renaming column names 
	for(i in 1:length(q)){
		names <- names(q[i])
		colnames(q[[i]]) <- paste0(names, "_Subsampled", "__", colnames(q[[i]]))
		q[[i]][q[[i]]==-1 | q[[i]]=="-1" | q[[i]]=="-1.0" | q[[i]]==-1.0] <- NA
	}
	names(q) <- c("Vertex_Gini_Index","Cluster_Gini_Index","Mean_Vertex_Size","Percentage_Max_Cluster_Size","Percentage_Max_Vertex_Size","Cluster_Reyni_Normalised", "Vertex_Reyni_Normalised", "D5_subsampled", "D10_subsampled", "D50_subsampled")
	analysis_matrices1 = q
	} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices1 <- vector(mode = "list", length = 0)
	} 

print("DONE 1")
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------

## File Number 2: Class Switching Summary : 
## Not relevant for TCR!!
## Checked for BCR 

if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_per_sequence_network_parameters_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check2 <- FALSE
	} else {
		check2 <- TRUE
	}
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.delim(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.Id"]) %in% ids_all),]
		p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"Isotype"]))),]
			
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"Isotype"]))),]
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
	
		###		
		id = as.character(p[,"X.Id"])
		ids = sort(unique(id))
		class = as.character(p[,"Isotype"])
		classes = sort(unique(class))
		reads_per_isotype = as.numeric(p[,"N.reads"])
		unique_reads_per_isotype = as.numeric(p[,"N.vertices"])
		m_reads_per_isotype = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
		m_unique_reads_per_isotype = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
		for(i in c(1:length(ids))){
			for (c in c(1:length(classes))){
				w = intersect(which(id==ids[i]), which(class==classes[c]))
				if(length(w)>=1){
					m_reads_per_isotype[i,c] = mean(reads_per_isotype[w])
					m_unique_reads_per_isotype[i,c] = mean(unique_reads_per_isotype[w])
				}
			}
		}
		c1 = c("Class_switched","IGHD,IGHM_mutated","IGHD,IGHM_unmutated")
		c2 = c( "IGHA1","IGHA2","IGHD","IGHE","IGHG1","IGHG2","IGHG3","IGHG4","IGHM"  )
		c2 = c2[which(c2 %in% classes)]
		m_reads_per_isotype_group = m_reads_per_isotype[,c(c1)]
		m_unique_reads_per_isotype_group = m_unique_reads_per_isotype[,c1]
		m_reads_per_isotype_single = m_reads_per_isotype[,c2]
		m_unique_reads_per_isotype_single = m_unique_reads_per_isotype[,c2]
		for(i in c(1:length(ids))){
			m_reads_per_isotype_group[i,] = m_reads_per_isotype_group[i,]*100/sum(m_reads_per_isotype_group[i,])
			m_unique_reads_per_isotype_group[i,] = m_unique_reads_per_isotype_group[i,]*100/sum(m_unique_reads_per_isotype_group[i,])
			m_reads_per_isotype_single[i,] = m_reads_per_isotype_single[i,]*100/sum(m_reads_per_isotype_single[i,])
			m_unique_reads_per_isotype_single[i,] = m_unique_reads_per_isotype_single[i,]*100/sum(m_unique_reads_per_isotype_single[i,])
		}
		analysis_names = c( "Percentage_Unique_VDJs_per_Isotype", "Percentage_Unique_VDJs_per_Isotype_Group")
		analysis_matrices = list(m_unique_reads_per_isotype_group, m_unique_reads_per_isotype_single)
		names(analysis_matrices) <- analysis_names

		for(i in 1:length(analysis_matrices)){
			names <- names(analysis_matrices[i])
			colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
		}
		analysis_matrices2 = analysis_matrices
		} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices2 <- vector(mode = "list", length = 0)	
		} 
} 
print("DONE 2")
##---------------------------------------------------------------------------------------------------------------------
## File Number 3: CDR3 lengths  (averaged)
## Does what it says on the tin...

file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_CDR3_lengths_overall_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check3 <- FALSE
} else {
	check3 <- TRUE
}
info = file.info(file)
if(info$size != 0) {	
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
	p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"isotype"]))),]
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
	p=p[which(as.numeric(p[,"Number_of_BCRs"])>10),]
	
	#######################
	p <- data.frame(p)
	## Want to make the relevant columns numeric
	p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
	## Convert '-1' to NA 
	if(any(p[!is.na(p)]==-1)){
		p[p==-1] <- NA
	}
	######################
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
	## Can't have different datatypes in dataframe
	p <- data.frame(p)
	p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
	############################

	id = as.character(p[,"X.sample"])
	ids = sort(unique(id))
	class = as.character(p[,"isotype"])
	chains = sort(unique(class))
	value = as.numeric(p[,"mean_CDR3_length"])
	reads = as.numeric(p[,"Number_of_BCRs"])

	values = matrix(data = -1, nrow = length(ids_all),ncol = length(chains), dimnames=c(list(ids_all), list(chains)))
	print_info = list()
	for(i in c(1:length(id))){
		values[id[i], class[i]] = value[i]
	}
	analysis_matrices = list(values)
	analysis_names = c("Mean_CDR3_Lengths")
	names(analysis_matrices) <- analysis_names
	for(i in 1:length(analysis_matrices)){
		names <- names(analysis_matrices[i])
		colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
		analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA

	} 
	analysis_matrices3 = analysis_matrices
	print("DONE 3")
	} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices3 <- vector(mode = "list", length = 0)	
	}
##---------------------------------------------------------------------------------------------------------------------
## File Number 4: ALL SHM : 
## NOT FOR TCR 
## Checked for BCR
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_SHM_Unmutated_sequences_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check4 <- FALSE
	} else {
		check4 <- TRUE
	}
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.Sample"]) %in% ids_all),]
		p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"isotype"]))),]
		p=p[which(is.na(as.numeric(p[,"mean.mutations"]))==F),]
		p=p[which(as.numeric(p[,"number.of.unique.sequences"])>10),]
		
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
		############################
			
		id = as.character(p[,"X.Sample"])
		ids = sort(unique(id))
		class = as.character(p[,"isotype"])
		chains = sort(unique(class))
		value = as.numeric(p[,"mean.mutations"])
		reads = as.numeric(p[,"number.of.unique.sequences"])
		perc_unmutated = as.numeric(p[,"perc_unumtated"])

		values = matrix(data = -1, nrow = length(ids_all),ncol = length(chains), dimnames=c(list(ids_all), list(chains)))
		m_perc_unmutated = matrix(data = -1, nrow = length(ids_all),ncol = length(chains), dimnames=c(list(ids_all), list(chains)))
		print_info = list()
		for(i in c(1:length(id))){
			values[id[i], class[i]] = value[i]
			m_perc_unmutated[id[i], class[i]] = perc_unmutated[i]
		}

		analysis_matrices = list(values, m_perc_unmutated)
		analysis_names = c("Mean_SHM_per_VDJ","Percentage_Unmutated")

		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			names <- names(analysis_matrices[i])
			colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA

		} 
		analysis_matrices4 = analysis_matrices
		} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices4 <- vector(mode = "list", length = 0)
		}
}
print("DONE 4")
##----------------------------------------------------------------------
## File Number 5: 
## NOT FOR TCR 
## CHECKED FOR BCR  
# raw overlapping frequencies
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Isotype_normalised_overlap_frequencies_uniq_", iso_type, ".txt")
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
		p=p[setdiff(c(1:length(p[,1])), grep("IGHG4",as.character(p[,"iso1"]))),]
		p=p[setdiff(c(1:length(p[,1])), grep("IGHG4",as.character(p[,"iso2"]))),]

		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1, -4, -5)] <- apply(p[ , c(-1, -4, -5)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[, c(-1, -4, -5)] <- apply(p[ , c(-1, -4, -5)], 2, as.numeric)
		############################
			
		id = as.character(p[,"X.sample"])
		ids = sort(unique(id))
		class1 = as.character(p[,"iso1"])
		class2 = as.character(p[,"iso2"])
		classes = sort(unique(class1))
		mean_overlap_proportion1 = as.numeric(p[,"mean_overlap"])
		class12 = apply(cbind(class1, class2),1,paste,collapse = "-")
		class12s = sort(unique(class12))
		overlap = matrix(data = -1, nrow = length(ids_all),ncol = length(class12s), dimnames=c(list(ids_all), list(class12s)))
		for(i in c(1:length(ids))){
			w = which(id==ids[i])
			if(length(w)>0){
				overlap[ids[i],]=0
				overlap[ids[i], class12[w]]= mean_overlap_proportion1[w]
			}
		}
		analysis_names = "Isotype_normalised_overlap_frequencies"
		analysis_matrices = list(overlap)
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			names <- names(analysis_matrices[i])
			colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA

		} 
		analysis_matrices5 = analysis_matrices
	} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices5 <- vector(mode = "list", length = 0)	
	}  
}

###normalised overlapping frequencies
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Isotype_overlapping_frequencies_", iso_type, ".txt")
	info = file.info(file)
	if(info$size != 0) {
	
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.ID"]) %in% ids_all),]
		p=p[setdiff(c(1:length(p[,1])), grep("IGHG4",as.character(p[,"isotype1"]))),]
		p=p[setdiff(c(1:length(p[,1])), grep("IGHG4",as.character(p[,"isotype2"]))),]
		
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1,-3, -5, -6)] <- apply(p[ , c(-1,-3, -5, -6)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[, c(-1,-3, -5, -6)] <- apply(p[ ,c(-1,-3, -5, -6)], 2, as.numeric)
		############################
		
		id = as.character(p[,"X.ID"])
		ids = sort(unique(id))
		class1 = as.character(p[,"isotype1"])
		class2 = as.character(p[,"isotype2"])
		classes = sort(unique(class1))
		mean_overlap_proportion1 = as.numeric(p[,"mean_overlap"])
		class12 = apply(cbind(class1, class2),1,paste,collapse = "-")
		class12s = sort(unique(class12))
		overlap = matrix(data = -1, nrow = length(ids_all),ncol = length(class12s), dimnames=c(list(ids_all), list(class12s)))
		for(i in c(1:length(ids))){
			w = which(id==ids[i])
			if(length(w)>0){
				overlap[ids[i],]=0
				overlap[ids[i], class12[w]]= mean_overlap_proportion1[w]
			}
		}
		analysis_names = "Isotype_overlapping_frequencies"
		analysis_matrices = list(overlap)
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			names <- names(analysis_matrices[i])
			colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA

		} 
		analysis_matrices5b = analysis_matrices
		} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices5b <- vector(mode = "list", length = 0)	
		} 
}


print("DONE 5")
##---------------------------------------------------------------------------------------------------------------------
## File Number 6: SHM  split across different regions 
## NOT FOR TCR 
## checked for BCR 
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_SHM_Mutation_summmary_selection_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check6 <- FALSE
	} else {
		check6 <- TRUE
	}
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
		p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"chain"]))),]
		p=p[which(as.numeric(p[,"total_BCRs_counted"])>10),]
		
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1,-2, -3)] <- apply(p[ , c(-1,-2, -3)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[,  c(-1,-2, -3)] <- apply(p[ , c(-1,-2, -3)], 2, as.numeric)
		############################
		
		id = as.character(p[,"X.sample"])
		ids = sort(unique(id))
		class = as.character(p[,"chain"])
		type = as.character(p[,"count_type"])
		types = sort(unique(type))
		chains = sort(unique(class))
		value = as.numeric(p[,"mean.value"])
		value_list = list()
		for(t in c(1:length(types))){
			values = matrix(data = -1, nrow = length(ids_all),ncol = length(chains), dimnames=c(list(ids_all), list(chains)))
			w = which(type== types[t])
			for(i in c(1:length(w))){
				values[id[w[i]], class[w[i]]] = value[w[i]]}
				value_list = c(value_list, list(values))
		}
		a = value_list [[which(types=="mean CDR_mm per BCR")]]+ value_list [[which(types=="mean FWR_mm per BCR")]]
		b = value_list [[which(types=="mean nonsilent per BCR" )]]+ value_list [[which(types=="mean silent per BCR")]]
		types = c(types, "Mean mutations per BCR")
		value_list = c(value_list, list(a))
		w = which(types %in% c("mean CDR_FWR_ratio", "Mean mutations per BCR" ,"FWR3_mm", "mean CDR_mm per BCR","mean FWR_mm per BCR"  ))
		analysis_matrices = value_list[w]
		analysis_names = types[w]
		analysis_names <-  gsub("_",  " ", types[w])
		names(analysis_matrices) <- c("mean_CDR_FWR_ratio", "mean_mutations_per_VDJ" ,"FWR3_mm", "mean_CDR_mm_per_VDJ","mean_FWR_mm_per_VDJ"  )
		for(i in 1:length(analysis_matrices)){
			names <- names(analysis_matrices[i])
			colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA

		} 
		analysis_matrices6 = analysis_matrices
		} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices6 <- vector(mode = "list", length = 0)	
		}
}
print("DONE 6")
##---------------------------------------------------------------------------------------------------------------------
## File Number 7: Secondary Rearangments by isotype
## Checked for TCR 
## Checked for BCR  
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Secondary_rearrangements_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check7 <- FALSE
	} else {
		check7 <- TRUE
	}

	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
		p=p[grep("IGH",as.character(p[,"chain"])),]
		w = setdiff(p[,"chain"], p[,"chain"][grep("IGHV",as.character(p[,"chain"]))])
		p=p[which(as.character(p[,"chain"]) %in% w),]
		w = setdiff(p[,"chain"], p[,"chain"][grep("mut",as.character(p[,"chain"]))])
		p=p[which(as.character(p[,"chain"]) %in% w),]
		
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[,  c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
		############################
		
		id = as.character(p[,"X.sample"])
		ids = sort(unique(id))
		class = as.character(p[,"chain"])
		id_replacement_freq =as.numeric(p[,"percentage"])
		total_isotype = as.numeric(p[,"total"])
		w1 = which(total_isotype>50)
		classes = sort(unique(class))
		types = c(list(id_replacement_freq))
		analysis_name = c("V gene replacement frequency")
		means = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
		for(i in c(1:length(ids))){
			w = intersect(which(id==ids[i]), w1)
			means[ids[i], class[w]] = id_replacement_freq[w]
		}
		analysis_matrices = list(means)
		analysis_names = c("V_gene_replacement_frequency")
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			names <- names(analysis_matrices[i])
			colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
		} 
		analysis_matrices7 = analysis_matrices
		names(analysis_matrices7) <- analysis_names
		} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices7 <- vector(mode = "list", length = 0)
		} 
} 

if(chain_vdj %like% "T"){
	file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Secondary_rearrangements_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check7 <- FALSE
	} else {
		check7 <- TRUE
	}
	info = file.info(file)
	if(info$size != 0) {
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		id = as.character(p[,"X.sample"])
		
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[,  c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
		############################
		
		ids = sort(unique(id))
		class = as.character(p[,"chain"])
		id_replacement_freq =as.numeric(p[,"percentage"])
		total_isotype = as.numeric(p[,"total"])
		w1 = which(total_isotype>50)
		classes = sort(unique(class))
		types = c(list(id_replacement_freq))
		analysis_name = c("V_gene_replacement_frequency")
		means = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
		for(i in c(1:length(ids))){
			w = intersect(which(id==ids[i]), w1)
			means[ids[i], class[w]] = id_replacement_freq[w]
		}
		analysis_matrices = list(means)
		analysis_names = c("V gene replacement frequency")
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			names <- names(analysis_matrices[i])
			colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA
		} 
		analysis_matrices7 = analysis_matrices
		names(analysis_matrices7) <- analysis_name
		} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices7 <- vector(mode = "list", length = 0)
		} 
} 
print("DONE 7")
##---------------------------------------------------------------------------------------------------------------------
## File Number 8: Cluster Expansion (d5,d10, d50)
## checked BCR
## checked TCR 
file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_expansion_isotype_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check8 <- FALSE
} else {
	check8 <- TRUE
}
info = file.info(file)
if(info$size != 0) {
		
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
	p=p[which(as.numeric(p[,"d5"]) !=-1),]
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
	
	#######################
	p <- data.frame(p)
	## Want to make the relevant columns numeric
	p[, c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
	## Convert '-1' to NA 
	if(any(p[!is.na(p)]==-1)){
		p[p==-1] <- NA
	}
	######################
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
	## Can't have different datatypes in dataframe
	p <- data.frame(p)
	p[,  c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
	############################
		
	id = as.character(p[,"X.sample"])
	ids = sort(unique(id))
	isotype = as.character(p[,"isotype"])
	d5 = as.numeric(p[,"d5"])
	d10 = as.numeric(p[,"d10"])
	d50 = as.numeric(p[,"d50"])
	classes = sort(unique(isotype))
	m_d5 = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	m_d10 = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	m_d50 = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	for(i in c(1:length(ids))){
		w = which(id==ids[i])
		m_d5[ids[i],isotype[w]]= d5[w]
		m_d10[ids[i],isotype[w]]= d10[w]
		m_d50[ids[i],isotype[w]]= d50[w]
	}
	analysis_matrices = c(list(m_d5),list(m_d10),list(m_d50))
	analysis_names = c("D5","D10","D50")
	names(analysis_matrices) <- analysis_names
	for(i in 1:length(analysis_matrices)){
		names <- names(analysis_matrices[i])
		colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
		analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA

	} 
	analysis_matrices8 = analysis_matrices
	names(analysis_matrices8) <- analysis_names
	} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices8 <- vector(mode = "list", length = 0)	
	}
print("DONE 8")
##---------------------------------------------------------------------------------------------------------------------
## File Number 9: ALL V gene IGHV4_34_quantification
## Autoreactivity 
## Not for TCR 
## Checked BCR 

if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_IGHV4_34_quantification_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check9 <- FALSE
	} else {
		check9 <- TRUE
	}
	info = file.info(file)
	if(info$size != 0) {
	
		p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
		p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
		p=p[which(as.numeric(p[,"total_all"]) >10),]
		
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[,  c(-1,-2)] <- apply(p[ , c(-1,-2)], 2, as.numeric)
		############################
	
	
		id = as.character(p[,"X.sample"])
		ids = sort(unique(id))
		isotype = as.character(p[,"isotype"])
		V4_34_AVY_total_unmut = as.numeric(p[,"V4_34_AVY_total_unmut"])*100/as.numeric(p[,"total_all"])
		V4_34_NHS_total_unmut = as.numeric(p[,"V4_34_NHS_total_unmut"])*100/as.numeric(p[,"total_all"])
		V4_34_AVY_NHS_total_unmut = as.numeric(p[,"V4_34_AVY_NHS_total_unmut"])*100/as.numeric(p[,"total_all"])
		classes = sort(unique(isotype))
		V4_34_AVY_NHS_total_unmut_count = as.numeric(p[,"V4_34_AVY_NHS_total_unmut"])
		all_counts = as.numeric(p[,"total_all"])

		m_V4_34_AVY_total_unmut = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
		m_V4_34_NHS_total_unmut = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
		m_V4_34_AVY_NHS_total_unmut = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
		m_V4_34_AVY_NHS_total_unmut_count = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
		m_all_count = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))

		for(i in c(1:length(ids))){
			w = which(id==ids[i])
			m_V4_34_AVY_total_unmut[ids[i],isotype[w]]= V4_34_AVY_total_unmut[w]
			m_V4_34_NHS_total_unmut[ids[i],isotype[w]]= V4_34_NHS_total_unmut[w]
			m_V4_34_AVY_NHS_total_unmut[ids[i],isotype[w]]= V4_34_AVY_NHS_total_unmut[w]
			m_V4_34_AVY_NHS_total_unmut_count[ids[i],isotype[w]]= V4_34_AVY_NHS_total_unmut_count[w]
			m_all_count[ids[i],isotype[w]]= all_counts[w]
		}
		IGHDM = rep(-1, length(ids_all))
		names(IGHDM)= ids_all
		class_switched = IGHDM
		for(i in c(1:length(ids))){
			w1 = which(classes %in% c("IGHA1", "IGHA2" , "IGHE",  "IGHG1" ,"IGHG2" ,"IGHG3" ,"IGHG4"))
			w2 = which(classes %in% c("IGHD", "IGHM"))
			w = which(m_V4_34_AVY_NHS_total_unmut[ids[i],]!=-1)
			IGHDM[ids[i]] = sum(m_V4_34_AVY_NHS_total_unmut_count[ids[i],intersect(w,w2)])*100/sum(m_all_count[ids[i],intersect(w,w2)])
			class_switched[ids[i]] = sum(m_V4_34_AVY_NHS_total_unmut_count[ids[i],intersect(w,w1)])*100/sum(m_all_count[ids[i],intersect(w,w1)])

		}
		m_V4_34_AVY_NHS_total_unmut = cbind(m_V4_34_AVY_NHS_total_unmut, IGHDM, class_switched)
		analysis_matrices = c(list(m_V4_34_AVY_total_unmut),list(m_V4_34_NHS_total_unmut),list(m_V4_34_AVY_NHS_total_unmut))
		analysis_names = c("V4_34_AVY_unmut","V4_34_NHS_unmut","V4_34_AVY_NHS_unmut")
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			names <- names(analysis_matrices[i])
			colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA

		} 
		names(analysis_matrices) <- analysis_names
		analysis_matrices9 = analysis_matrices	
		} else { 
				print(paste0("File: ", file, " IS EMPTY"))
				analysis_matrices9 <- vector(mode = "list", length = 0)	
		} 
}
print("DONE 9")
##---------------------------------------------------------------------------------------------------------------------
## File Number 10: Secondary Rearangments on clone size 
## Something about this is not right!

file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Secondary_rearrangements_clone_sizes_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check10 <- FALSE
} else {
	check10 <- TRUE
}
info = file.info(file)
if(info$size != 0) {
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	if(class(p)=="character"){
			p <- as.matrix(t(p))
	} 
	
	#######################
	p <- data.frame(p)
	## Want to make the relevant columns numeric
	p[, c(-1)] <- apply(p[ , c(-1)], 2, as.numeric)
	## Convert '-1' to NA 
	if(any(p[!is.na(p)]==-1)){
		p[p==-1] <- NA
	}
	######################
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
	## Can't have different datatypes in dataframe
	p <- data.frame(p)
	p[,  c(-1)] <- apply(p[ , c(-1)], 2, as.numeric)
	############################
		
	id = as.character(p[,"X.sample"])
	ids = sort(unique(id))
	d5_norm =as.numeric(p[,"d5_norm"])
	d5_secondary =as.numeric(p[,"d5_secondary"])
	mean_clone_size_norm =as.numeric(p[,"mean_clone_size_norm"])
	mean_clone_size_secondary =as.numeric(p[,"X_mean_clone_size_secondary"])
	w1 = which(as.numeric(p[,"n_secondary"])>5)
	groups = c(list(mean_clone_size_norm[w1]), list(mean_clone_size_secondary[w1] ))
	groups = c(list(d5_norm[w1]), list(d5_secondary[w1] ))
	headers = c("d5.norm","d5.secondary","mean.clone.size.norm", "mean.clone.size.secondary")
	m_all_count = matrix(data = -1, nrow = length(ids_all),ncol = length(headers), dimnames=c(list(ids_all), list(headers)))
	x=  cbind(d5_norm,d5_secondary,mean_clone_size_norm,mean_clone_size_secondary)
	m_all_count[id,] = x
	analysis_matrices = list(m_all_count)
	analysis_names = c("V_gene_replacement_clonal_expansion")
	names(analysis_matrices) <- analysis_names
	for(i in 1:length(analysis_matrices)){
		names <- names(analysis_matrices[i])
		colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
		analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA

	} 
	names(analysis_matrices) <- analysis_names
	analysis_matrices10 = analysis_matrices
	} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices10 <- vector(mode = "list", length = 0)	
	}

print("DONE 10")
##---------------------------------------------------------------------------------------------------------------------
## File Number 11: Sampled Secondary rearangement estimates per V gene
## Checked BCR
## Checked TCR 
file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Secondary_rearrangements_file_SAMPLED_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check7 <- FALSE
} else {
	check7 <- TRUE
}
info = file.info(file)
if(info$size != 0) {
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(p[,1]!="#sample"),]
	if(dim(p)[1]!=0){
		if(class(p)=="character"){
				p <- as.matrix(t(p))
		} 
		id = as.character(p[,"X.sample"])
		ids = sort(unique(id))
		
		#######################
		p <- data.frame(p)
		## Want to make the relevant columns numeric
		p[, c(-1, -2, -3)] <- apply(p[ , c(-1, -2, -3)], 2, as.numeric)
		## Convert '-1' to NA 
		if(any(p[!is.na(p)]==-1)){
			p[p==-1] <- NA
		}
		######################
		if(class(p)=="character"){
			p <- as.matrix(t(p))
		} 
		## Can't have different datatypes in dataframe
		p <- data.frame(p)
		p[,  c(-1, -2, -3)] <- apply(p[ , c(-1, -2, -3)], 2, as.numeric)
		############################
	
		class = as.character(p[,"isotype"])
		classes = sort(unique(class))
		id_replacement_freq_norm =as.numeric(p[,"uniq_id_replacement_freq"])/as.numeric(p[,"n_repeats"])
		m = matrix(data = 0, nrow = length(ids),ncol = length(classes), dimnames=c(list(ids), list(classes)))
		for(c in c(1:length(classes))){
			for(i in c(1:length(ids))){
				w = intersect(which(class== classes[c]), which(id== ids[i]))
				m[ids[i],classes[c]] = sum(id_replacement_freq_norm[w])
			}
		}
		analysis_matrices = list(m)
		analysis_names = c("Secondary_Rearrangements_SAMPLED")
		names(analysis_matrices) <- analysis_names
		for(i in 1:length(analysis_matrices)){
			names <- names(analysis_matrices[i])
			colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
			analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA

		} 
		analysis_matrices11 = analysis_matrices
		names(analysis_matrices11) <- analysis_name
		} else {
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices11 <- vector(mode = "list", length = 0)
		}
	} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices11 <- vector(mode = "list", length = 0)	
	} 
	
print("DONE 11")
##---------------------------------------------------------------------------------------------------------------------
## File Number 12 CDR Charge 
## Checked BCR
## Checked TCR 

file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_CDR_charge_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check7 <- FALSE
} else {
	check7 <- TRUE
}
info = file.info(file)
if(info$size != 0) {
		
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	if(class(p)=="character"){
			p <- as.matrix(t(p))
	} 
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	if(class(p)=="character"){
			p <- as.matrix(t(p))
	} 
	
	#######################
	p <- data.frame(p)
	## Want to make the relevant columns numeric
	p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
	## Convert '-1' to NA 
	if(any(p[!is.na(p)]==-1)){
		p[p==-1] <- NA
	}
	######################
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
	## Can't have different datatypes in dataframe
	p <- data.frame(p)
	p[,  c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
	############################
		
		
	id = as.character(p[,"X.sample"])
	ids = sort(unique(id))
	class = as.character(p[,"isotype"])
	id_replacement_freq =as.numeric(p[,"mean_CDR2.3_R_K_residues"])
	total_isotype = as.numeric(p[,"Number_of_BCRs"])
	w1 = which(total_isotype>0)
	classes = sort(unique(class))
	types = c(list(id_replacement_freq))
	analysis_name = c("CDR_Charge")
	means = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	for(i in c(1:length(ids))){
		w = intersect(which(id==ids[i]), w1)
		means[ids[i], class[w]] = id_replacement_freq[w]
	}
	analysis_matrices = list(means)
	analysis_names = c("CDR_Charge")
	names(analysis_matrices) <- analysis_names
	for(i in 1:length(analysis_matrices)){
		names <- names(analysis_matrices[i])
		colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
		analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA

	} 
	analysis_matrices12 = analysis_matrices
	names(analysis_matrices12) <- analysis_names
	} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices12 <- vector(mode = "list", length = 0)	
	} 	
print("DONE 12")
##---------------------------------------------------------------------------------------------------------------------
## File Number 13 CDR3 Charge 
## NEED TO CHECK FOR BCR!!!!!!!!!!!!!!!!!!!!!!!!

file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_CDR3_charge_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check7 <- FALSE
} else {
	check7 <- TRUE
}
info = file.info(file)
if(info$size != 0) {
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	if(class(p)=="character"){
			p <- as.matrix(t(p))
	} 
	id = as.character(p[,"X.sample"])
		
	#######################
	p <- data.frame(p)
	## Want to make the relevant columns numeric
	p[, c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
	## Convert '-1' to NA 
	if(any(p[!is.na(p)]==-1)){
		p[p==-1] <- NA
	}
	######################
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
	## Can't have different datatypes in dataframe
	p <- data.frame(p)
	p[,  c(-1, -2)] <- apply(p[ , c(-1, -2)], 2, as.numeric)
	############################
		
	ids = sort(unique(id))
	class = as.character(p[,"isotype"])
	id_replacement_freq =as.numeric(p[,"mean_CDR3_R_K_residues"])
	total_isotype = as.numeric(p[,"Number_of_BCRs"])
	w1 = which(total_isotype>0)
	classes = sort(unique(class))
	types = c(list(id_replacement_freq))
	analysis_name = c("CDR3_Charge")
	means = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
	for(i in c(1:length(ids))){
		w = intersect(which(id==ids[i]), w1)
		means[ids[i], class[w]] = id_replacement_freq[w]
	}
	analysis_matrices = list(means)
	analysis_names = c("CDR3_Charge")
	names(analysis_matrices) <- analysis_names
	for(i in 1:length(analysis_matrices)){
		names <- names(analysis_matrices[i])
		colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
		analysis_matrices[[i]][analysis_matrices[[i]]==-1 | analysis_matrices[[i]]=="-1" | analysis_matrices[[i]]=="-1.0" | analysis_matrices[[i]]==-1.0] <- NA

	} 
	analysis_matrices13 = analysis_matrices
	names(analysis_matrices13) <- analysis_names
	} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices13 <- vector(mode = "list", length = 0)	
	}	
print("DONE 13")

##---------------------------------------------------------------------------------------------------------------------
## File Number 14 V gene Usages  
## NEED TO CHECK FOR BCR!!!!!!!!!!!!!!!!!!!!!!!!
file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_grouped_isotype_frequency_", iso_type, ".txt")
p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
p <- data.frame(p)
p$uniq_read_freq <- as.numeric(p$uniq_read_freq)
		
	
#####
isotypes <- c("IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHEP2", "IGHGP", "IGHE")
other_class <- c("Class_switched", "IGHD,IGHM_mutated", "IGHD,IGHM_unmutated", "IGHD,IGHM_unmutated_singleton")
expansion_class <- c("unexpanded", "expanded")
		
if(chain_vdj %like% "BC"| chain_vdj %like% "I"){
	all_class <- c("ALL")
} else {
receptor_type <- counts_try$X.isotype[counts_try$min==max(counts_try$min)]
all_class <- receptor_type
}
		
#####
p_all <- p[p$class %in% all_class,]
p_iso <-p[p$class %in% isotypes,] 
		
## Calculate a percentage of repertoire which is each read 
p_all$percent_repertoire <- NA
for(i in unique(p_all[, "X.sample"])){
	sample_id <- i 
	sum_frequency <- sum(p_all$uniq_read_freq[p_all$X.sample==sample_id])
	p_all$percent_repertoire[p_all$X.sample==sample_id] <- ((p_all$uniq_read_freq[p_all$X.sample==sample_id])/sum_frequency)*100
} 
p_allx <- p_all[, c("X.sample", "percent_repertoire", "V.gene")]

## We fill in any missing combinations with 0 as 0 percent of repertoire 
a <- spread(p_allx, key = V.gene, value = percent_repertoire, fill=0)
rownames(a) <- a$X.sample
a$X.sample <- NULL
colnames(a) <- gsub("-", "_", colnames(a))
colnames(a) <- paste0(colnames(a), "__", all_class)
a <- as.matrix(a)
analysis_matrices14 = list(a)
names(analysis_matrices14) <- "V_GENE_USAGE"
## Save the dataframe!
write.table(a, paste0(outputdir, "Summary/V_Gene_usage", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t", row.names=TRUE)
print("DONE 14: V GENE Usages")
		
################################
## Get Hydrophobicity
## These will be generated at this stage!! 
print("Calculating Hydrophobicity")
all_files <- list.files(paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/", iso_type, '/Classification_per_sequence'), full.name=TRUE)

mean_charges <- data.frame()
for(c in 1:length(all_files)){
	file_use <- read.delim(all_files[c])
	file_use$Hydrophobicity <- hydrophobicity(file_use$CDR3, scale = "KyteDoolittle")
	if(any(file_use$all_classes %like% ',')){
		#print(paste0("Double Assignment Present in ", all_files[c], " - will duplicate bad row"))
		bad_rows <- file_use[file_use$all_classes %like% ',',]
		bad_rows1 <- bad_rows
		bad_rows2 <- bad_rows
		two_classes <- str_split_fixed(bad_rows$all_classes, ",", 2)
		bad_rows1$all_classes <- two_classes[,1]
		bad_rows2$all_classes <- two_classes[,2]
		rownames(bad_rows2) <- paste0(rownames(bad_rows1), "_duplicate")
		all_bad_rows <- rbind(bad_rows1, bad_rows2)
		## Remove the bad row 
		file_use <- file_use[-c(as.numeric(rownames(bad_rows))),]
		## add in the dubplication row 
		file_use <- rbind(file_use, all_bad_rows)
	}
	mean_hydrophobicity <- file_use %>%group_by(all_classes) %>%summarise_at(vars(Hydrophobicity), list(CDR3_hydrophobicity = mean, CDR3_kurtosis=kurtosis, CDR3_skewness=skewness, CDR3_bimodality=bimodality_coefficient))
	all_mean <- mean(file_use$Hydrophobicity)
	all_skew <- skewness(file_use$Hydrophobicity)
	all_kurtosis <- kurtosis(file_use$Hydrophobicity)
	all_bimodality <- bimodality_coefficient(file_use$Hydrophobicity)
	new_row <- c("ALL", all_mean, all_kurtosis, all_skew, all_bimodality)
	## ALSO need to do for all!
	mean_hydrophobicity <- data.frame(mean_hydrophobicity)
	mean_hydrophobicity <- rbind(mean_hydrophobicity, new_row)
	mean_hydrophobicity$sample <- unique(file_use$X.ID)
	mean_hydrophobicity <- reshape(mean_hydrophobicity, idvar = "sample", timevar = "all_classes", direction = "wide")
	mean_hydrophobicity[mean_hydrophobicity=="NaN"] <- NA
	mean_charges <- suppressMessages(plyr::join(mean_hydrophobicity, mean_charges, type="full"))
} 
rownames(mean_charges) <- mean_charges$sample
mean_charges$sample <- NULL
mean_charges <- as.matrix(mean_charges)
mean_charges <- mean_charges[rownames(mean_charges) %in% ids_all,]	
colnames(mean_charges) <- gsub("\\.", "__", colnames(mean_charges))
storage.mode(mean_charges) <- "numeric"
analysis_matrices15 = list(mean_charges)
names(analysis_matrices15) <- "CDR3_Hydrophobicity"
#return(mean_charges)
print("DONE 15: CDR3 Hydrophobicity")	

################################
##Skewness and kurtosis (statistics describing the CDR3 lengths)

print("Calculating Kurtosis and Skewness of CDR3 Lengths")
all_files <- list.files(paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/", iso_type, '/CDR3_length_distribution'), full.name=TRUE)

cdr3_summary <- data.frame()
for(c in 1:length(all_files)){
	file_use <- read.delim(all_files[c])
	cdr3_stats <- data.frame()
	for(i in 1:length(unique(file_use$isotype))){
		subset_data <- file_use[file_use$isotype==unique(file_use$isotype)[i],]
		length_dist <- rep(subset_data$CDR3_length, times=subset_data$Number_of_BCRs)
		# calculate kurotis per isotype
		## Kurtosis 3 - roughly normal, >3 it has a sharp peak, less than 3 it is very wide 
		kurtosis_new <- kurtosis(length_dist)
		if(kurtosis_new=="NaN"){
			kurtosis_new <- NA
		} 
		## Skewness 0 - symetric normally distributed, > 0 positively skewed, <0 negatively skewed (values >/< mean)
		skewness_new <- skewness(length_dist)
		if(skewness_new=="NaN"){
			skewness_new <- NA
		} 
		bimodality_new <- bimodality_coefficient(length_dist)
		row_stats <- c(unique(file_use$isotype)[i], kurtosis_new, skewness_new, bimodality_new)
		cdr3_stats <- rbind(cdr3_stats, row_stats)
	}
	colnames(cdr3_stats) <- c("class", "CDR3_Length_kurtosis", "CDR3_Length_skewness", "CDR3_Length_bimodality")
	## Need to do for all 
	all_use <- rep(file_use$CDR3_length, times=file_use$Number_of_BCRs)
	kurtosis_all <- kurtosis(all_use)
	skewness_all <- skewness(all_use)
	bimodality_all <- bimodality_coefficient(all_use)
	row_all <- c("ALL", kurtosis_all, skewness_all,bimodality_all)
	cdr3_stats <- rbind(cdr3_stats, row_all)
	cdr3_stats$sample <- unique(file_use$X.sample)
	kurtosis_df <- reshape(cdr3_stats, idvar = "sample", timevar = "class", direction = "wide")
	kurtosis_df[is.na(kurtosis_df)] <- NA
	cdr3_summary <- suppressMessages(plyr::join(kurtosis_df, cdr3_summary, type="full"))
}

rownames(cdr3_summary) <- cdr3_summary$sample
cdr3_summary$sample <- NULL
#cdr3_summary <- as.numeric(cdr3_summary)
cdr3_summary <- as.matrix(cdr3_summary)
cdr3_summary <- cdr3_summary[rownames(cdr3_summary) %in% ids_all,]	
colnames(cdr3_summary) <- gsub("\\.", "__", colnames(cdr3_summary))
storage.mode(cdr3_summary) <- "numeric"
analysis_matrices15.b = list(cdr3_summary)
names(analysis_matrices15.b) <- "CDR3_Summary"
#return(mean_charges)
print("DONE 15b: CDR3 Kurtosis and Skewness")	


################################
## Get proportion productive

print("Calculating Proportion of Reads Productive")
path <- paste0(outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT")
files <- list.files(path, full.name=TRUE)
files <- grep("_1_Summary", files, value=TRUE)
files <- grep("productive", files, value=TRUE, invert=TRUE)	
df_list <- lapply(files, fread, header = FALSE, sep="\t", fill=TRUE, colClasses='character')
# identify the sample they came from 
d <- str_split(files, 'IMGT_SPLIT/') 
d <- sapply(d, "[[", 2)  
d <- gsub("_1_Summary.txt", "", d)
d <- gsub("IMGT_", "", d)
d <- gsub("BCR_", "", d)
 
# Name each dataframe with the run and filename
names(df_list) <- d
# Create combined dataframe  
df <- df_list %>% bind_rows(.id = 'Sample')
df <- data.frame(df)
		
FullData <- df[!df[,4]=="rearranged sequence (but no junction found) (see comment)",]
FullData <- FullData[!FullData[,4]=="No rearrangement found",]
FullData <- FullData[!FullData[,4]=="No results",]		
FullData[,4][FullData[,4]=="unproductive (see comment)"] <- "unproductive"
FullData[,4][FullData[,4]=="productive (see comment)"] <- "productive"
functionality <- unique(FullData[,4])
functionality <- prop.table(table(FullData$Sample, FullData[,4]),1)

## We must add on the chain info!
if(chain_vdj %like% "BC"| chain_vdj %like% "I"){
	all_class <- c("ALL")
} else {
receptor_type <- counts_try$X.isotype[counts_try$min==max(counts_try$min)]
all_class <- receptor_type
}
colnames(functionality) <- paste0("prop_", colnames(functionality), "__", all_class)
functionality <- as.data.frame.matrix(functionality)
if(iso_type == "UNPRODUCTIVE"){
	rownames(functionality) <- paste0(rownames(functionality), "_unproductive")
}
if(iso_type=="PRODUCTIVE"){
	rownames(functionality) <- paste0(rownames(functionality), "_productive")
}
functionality <- functionality[rownames(functionality) %in% ids_all,]	
functionality <- as.matrix(functionality)
analysis_matrices16 = list(functionality)
names(analysis_matrices16) <- "VDJ_Functionality"
print("DONE 16: VDJ Functionality")	


##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
## PART 2!!!!!!
## COMPOSING THE OVERALL MATRIX 

if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	print_info = c(analysis_matrices1, analysis_matrices2,  analysis_matrices3, analysis_matrices4, analysis_matrices5, analysis_matrices5b, analysis_matrices6, analysis_matrices7, analysis_matrices8, analysis_matrices9, analysis_matrices10, analysis_matrices11,  analysis_matrices12, analysis_matrices13, analysis_matrices14, analysis_matrices15,analysis_matrices15.b, analysis_matrices16 )
} else {
	print_info = c(analysis_matrices1, analysis_matrices3, analysis_matrices7, analysis_matrices8, analysis_matrices10, analysis_matrices11, analysis_matrices12, analysis_matrices13, analysis_matrices14, analysis_matrices15, analysis_matrices15.b, analysis_matrices16)
}

## Combining all the matrices into one big dataframe
for(i in c(1:length(print_info))){
	print_info[[i]] <- data.frame(print_info[[i]])
	print_info[[i]]$sample <- row.names(print_info[[i]])
} 
overall_matrix <- print_info %>% purrr::reduce(full_join, by="sample")

colnames(overall_matrix) = gsub(" ","_", colnames(overall_matrix), fixed= T)
overall_matrix <- data.frame(overall_matrix)
overall_matrix$sample <- gsub("BCR_", "", overall_matrix$sample)
overall_matrix$sample <- gsub("TCRA_", "", overall_matrix$sample)
overall_matrix$sample <- gsub("TCRB_", "", overall_matrix$sample)
overall_matrix$sample <- gsub("TCRG_", "", overall_matrix$sample)
overall_matrix$sample <- gsub("TCRD_", "", overall_matrix$sample)
overall_matrix$sample <- gsub("BCR_", "", overall_matrix$sample)
overall_matrix$sample <- gsub("TCR_", "", overall_matrix$sample)
# Merge with readdepths 
read_depths_all$SampleIDforDepths <- gsub("BCR_", "", read_depths_all$SampleIDforDepths)
read_depths_all$SampleIDforDepths <- gsub("TCR_", "", read_depths_all$SampleIDforDepths)
overall_matrix <- merge(overall_matrix, read_depths_all, by.x="sample", by.y="SampleIDforDepths")
##---------------------------------------------------------------------------------------------------------------------
## EDIT TO FILTER FOR DIFFERENT READ DEPTHS 
###############################################################################################
## GOING TO REMOVE SAMPLES FROM WHICH THE READ DEPTH IS LOWER THAN SPECIFIED
overall_matrix <- overall_matrix[overall_matrix$ReadDepth >=  subsampled_depth_allx, ]
print("###############################################################################################")
print("###############################################################################################")
print("ATTENTION!")
print(paste0("REMOVED Samples with read depth less than: ", subsampled_depth_allx))
print(paste0("REMOVED ", length(samples_to_low_all), " Samples"))
print("###############################################################################################")
print("###############################################################################################")
###############################################################################################

##---------------------------------------------------------------------------------------------------------------------
## For TCRs we must filter for the chain of interest!
if(chain_vdj %like% "T"){
	counts_try <- counts_used[(counts_used$X.isotype  != "ALL" & counts_used$X.isotype  != "all"),]
	receptor <- counts_try$X.isotype[counts_try$min==max(counts_try$min)]
	if(receptor=="TRBC"){
		receptor <- c("TRBC1", "TRBC2", "TRBC")
	} else if (receptor=="TRCG" | receptor=="TRGC"  ){
		receptor <- c("TRGC1", "TRGC2", "TRGC", "TCRCG1", "TCRCG2", "TCRCG" )
	} 
	keeping_columns <- c()
	for(i in 1:length(receptor)){
		keep <- grep(receptor[i], colnames(overall_matrix), value=TRUE)
		keeping_columns <- c(keeping_columns, keep)
	}
	keeping_columns<- c("sample", keeping_columns, "V_gene_replacement_clonal_expansion__d5.norm", "V_gene_replacement_clonal_expansion__d5.secondary", "V_gene_replacement_clonal_expansion__mean.clone.size.norm", "V_gene_replacement_clonal_expansion__mean.clone.size.secondary", "ReadDepth")
	keep <- colnames(overall_matrix)[colnames(overall_matrix) %in% keeping_columns]
	overall_matrix <- overall_matrix[, c(keep)]
}
## this will be used later for naming the columns of the final matrix!
if(chain_vdj %like% "T"){
	counts_try <- counts_used[(counts_used$X.isotype  != "ALL" & counts_used$X.isotype  != "all"),]
	receptor <- counts_try$X.isotype[counts_try$min==max(counts_try$min)]
}
print("Filtered for valid Metrics from file")
##---------------------------------------------------------------------------------------------------------------------
## Add in the proportions of different isotypes across groups!
## Only relevant for BCRs
if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
		file <- paste0(outputdir, "/Summary/IMGT/IMGT_Prop_SHM_", iso_type, ".txt")
		proportions_file <- read.delim(file, sep="\t", header=TRUE)
		colnames(proportions_file) <- gsub("\\.", "_", colnames(proportions_file))
		oo <- data.frame(str_split_fixed(colnames(proportions_file), "__", 2))
		oo <- oo$X1
		oo <- paste0(oo, "__", "ALL")
		colnames(proportions_file) <- oo
		if(iso_type=="UNPRODUCTIVE"){
			proportions_file$Sample__ALL <- paste0(proportions_file$Sample__ALL, "_unproductive")
		}
		if(iso_type=="PRODUCTIVE"){
			proportions_file$Sample__ALL <- paste0(proportions_file$Sample__ALL, "_productive")
		}	
		overall_matrix <- merge(overall_matrix, proportions_file, by.x="sample", by.y=paste0("Sample__ALL"))
		overall_matrix$Sample__ALL <- NULL
} 
overall_matrix <- overall_matrix %>% select(ReadDepth, everything())

## REPLACE ANY INSTANCES OF '-1' (meaning to small sample size so not calculated with NA!" 
################We want to add in percentage each group  (SHM - relevant for BCRS) 	
## IF we keep this we will NA the CDR3 hydrophobicity but I think its okay as I did this in each step 
#overall_matrix[overall_matrix==-1 | overall_matrix=="-1" | overall_matrix=="-1.0" | overall_matrix==-1.0] <- "NA"
## Remove any columns which are empty (e.g. nothing was calculated) 
empty <- apply(overall_matrix, 2, function(x){length(which(x==0 | x=="NA"))})
empty_cols <- empty[empty==dim(overall_matrix)[1]]
overall_matrix<- overall_matrix[, c(!colnames(overall_matrix) %in% names(empty_cols))]
#if(any(overall_matrix[!is.na(overall_matrix)]=="-1.0" | overall_matrix[!is.na(overall_matrix)]==-1.0  | overall_matrix[!is.na(overall_matrix)]==-1| overall_matrix[!is.na(overall_matrix)]=="-1")){
#	print("WARNING SOME '-1' aka NA values remain in dataset")
#} else {
#	print("All NAs successfully removed!")
#} 
#return(overall_matrix)
##---------------------------------------------------------------
# Save the raw UNFILTERED output!!! 
################################
out_file_table=paste0(outputdir, "Summary/All_raw_values_unfiltered_", subsampled_depth_all, "_", iso_type, ".txt")
write.table(overall_matrix, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
#################################################
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------
## PART 3
## Calculating the relationship between read depth and the metric 

new <- overall_matrix

## ***** Default of corelation test will be to omit na values *****
## Calculate correlation and missingness 
# note that missingness is really percentage of samples which are not na! (so high is good!)
values <- c("Metric", "percentage_present", "pval", "correlation")
for(i in 3:(length(colnames(new)))){
	variable <- new[,i]
	id <- colnames(new)[i]
	depths <- new$ReadDepth
	data <- data.frame(cbind(variable, depths))
	#Hydrophobicity can be less than 0
	if(!id %like% "Hydrophobicity"){
		data[data=="-1"] <- "NA"
	}
	suppressWarnings(data$variable <- as.numeric(data$variable))
	suppressWarnings(data$depths <- as.numeric(data$depths))
	missingness <- data$variable[!is.na(data$variable)]
	if(length(missingness) >= (0.25*length(ids_all))){
			missingness <- length(missingness) / length(ids_all) * 100
			rval <- cor.test(data$variable,data$depths, method="spearman", exact = FALSE)
			pval <- rval$p.value
			corval <- rval$estimate
			names(corval) <- NULL
	} else {
			missingness <- length(missingness) / length(ids_all) * 100
			rval <- NA
			pval <- NA
			corval <- NA
	}
	result <- c(id, missingness, pval, corval)
	values <- rbind(values, result)
}

# Reformat data frame to annotate with significance p<0.05
# Reformat to warn with high missingness  
colnames(values) <- values[1,]
values <- values[2:length(values[,1]),]
values <- data.frame(values)
values$sig <- "significant"
values$pval <- as.numeric(as.character(values$pval))
values$sig[values$pval > 0.05] <- ""
values$sig[values$missingness < 0.6] <- "warning >60% missingness"
values$sig[is.na(values$pval)] <- "WARNING MISSINGNESS >75% STAT not Calculated"
## Fill in with subampled depth for those metrics that were subsampled. 
subsampled_data <- grep("Subsampled", values$Metric, value=TRUE)
data_metrics <- values$Metric 
data_metrics <- str_split(data_metrics, "__")
iso <-  data.frame(matrix(unlist(data_metrics), ncol=2, byrow=TRUE))
iso <- iso[,2]
values$subsampled_depth <- NA 
values$isotype <- iso
values$subsampled <- NA
values$isotype[!values$isotype %in% counts_used$X.isotype] <- NA
for(i in 1:length(rownames(values))){
	row_id <- values$Metric[i]
	if(row_id %in% subsampled_data){
		isotype <- as.character(values$isotype[i])
		sub_depth <- counts_used$min[counts_used$X.isotype==isotype]
		values$subsampled_depth[i] <- sub_depth
		values$subsampled[i] <- "YES"
	} else {
		values$subsampled_depth[i] <- NA
		values$subsampled[i] <- NA
	}
}
## Add another subsampled depth
## where we specified subsample depth in script!
subsampled_data <- grep("SAMPLED", values$Metric, value=TRUE)
values$subsampled_depth[values$Metric %in% subsampled_data] <- 250
# Remove spurious column names :D 
values$isotype[!values$isotype %in% counts_used$X.isotype] <- NA
values$percentage_present <- as.numeric(as.character(values$percentage_present))
## Assess columns with high missingness
bad_columns <- values$Metric[values$percentage_present<50]
bad_columns2 <- values$Metric[values$percentage_present<25]
bad_columns3 <- values$Metric[values$percentage_present<75]
## Save Significance to a txt file 
write.table(values, paste0(outputdir, "Summary/isotyper_metrics_summary_stats", subsampled_depth_all,  "_", iso_type, ".txt"), sep="\t")
#Extract just those that are significant 
#Save file
## If there are many that are significant we may want to consider rerunning with higher
## redepth 
significant_samples <- values[values$sig=="significant",]
write.table(significant_samples, paste0(outputdir, "Summary/isotyper_metrics_summary_stats_sig_correlated_readdepth_", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t")

##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
## PART 3 
## Assessing which metrics to keeep 
## Baed on Rachaels script 
p <- as.matrix(overall_matrix)
#return(p)
## Remove read depth column 
p1 = p[,which(colnames(p)!="ReadDepth" )]
samples = p1[, "sample"]
## Converting to a numerical matrix 
p1 = p1[,c(2:length(p1[1,]))]
headers = colnames(p1)
mat = matrix(data = -1, nrow = length(samples), ncol = length(headers), dimnames = c(list(samples), list(headers)))
for(i in c(1:length(headers))){
	mat[,headers[i]] = as.numeric(p1[, headers[i]])
}
## remove all incomplete rows then columns
nz_cols = apply(mat, 2, function(x){length(which(!is.na(x)))})
nz_rows = apply(mat, 1, function(x){length(which(!is.na(x)))})
## Create plots directory if not present
if (!file.exists(paste0(outputdir, "Plots/ISOTYPER"))){
    dir.create(paste0(outputdir, "Plots/ISOTYPER"))
}
###############
# Doing plot 
pdf(paste0(outputdir, "Plots/ISOTYPER/Metrics_Missingness_Summary_", subsampled_depth_all, "_", iso_type,".pdf"), width=6, height=6)
par(mfrow= c(2,2), mar = c(5,5,3,3))
##1
nz_cols = apply(mat, 2, function(x){length(which(!is.na(x)))})
threshold = quantile(nz_cols, 0.25)
plot(sort(nz_cols), xlab = "feature rank", ylab = "number of non-NA values", pch = 21, bg = "blue", col = "blue",main=paste0("run1 threshold: ", threshold))
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat1 = mat[,which(nz_cols>= threshold)]
##2
nz_rows = apply(mat1, 1, function(x){length(which(!is.na(x)))})
threshold = quantile(nz_rows, 0.1)
plot(sort(nz_rows), xlab = "sample rank", ylab = "number of non-NA values", pch = 21, bg = "red", col = "red",main=paste0("run2 threshold: ", threshold))
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat2 = mat1[which(nz_rows>= threshold),]
##3
nz_cols = apply(mat2, 2, function(x){length(which(!is.na(x)))})
threshold = quantile(nz_cols, 0.25)
plot(sort(nz_cols), xlab = "feature rank", ylab = "number of non-NA values", pch = 21, bg = "blue", col = "blue",main=paste0("run3 threshold: ", threshold))
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat3 = mat2[,which(nz_cols>= threshold)]
##4
nz_rows = apply(mat3, 1, function(x){length(which(!is.na(x)))})
threshold = quantile(nz_rows, 0.1)
plot(sort(nz_rows), xlab = "sample rank", ylab = "number of non-NA values", pch = 21, bg = "red", col = "red",main=paste0("run4 threshold: ", threshold))
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat4 = mat3[which(nz_rows>= threshold),]
dev.off()
print("plot1")

#################################
#################################
#################################
## FILTERING 
#################################
#################################

## Removing repetitive annotations 
## Keeping those that are biologically informative 
## Removing those with high misingness! 
if(iso_type=="UNPRODUCTIVE"){
	threshold_new <- 90 
} else {
	threshold_new <- 60
}


if(chain_vdj %like% "BC"| chain_vdj %like% "I"){
 annotated_metrics_to_keep <- read.delim('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/Isotyper/inclusion_metrics.txt', header=TRUE)
 annotated_metrics_to_keep_colnames <- annotated_metrics_to_keep$Metric[annotated_metrics_to_keep$Include_Metric=="1"]
 keep <- values$Metric[values$percentage_present > threshold_new]
 keep <- unique(c(keep, annotated_metrics_to_keep_colnames))
 keep <- keep[!keep %in% annotated_metrics_to_keep$Metric[annotated_metrics_to_keep$Include_Metric=="0"]]
 keep <- keep[keep %in% values$Metric[values$percentage_present >= threshold_new]]
 keep <- keep[keep %in% colnames(mat)]
 mat_filtered = mat[,c(keep)]
 print(paste0("Removed ", (dim(overall_matrix)[2]-length(keep)), " Metrics"))
} 

if(chain_vdj %like% "T"){
 metrics_exclude <- grep("SHM", colnames(mat), value=TRUE)
 keepers <- colnames(mat)[!colnames(mat) %in% metrics_exclude]
 mat_filtered = mat[, c(keepers)]
 extra_exclude <- values$Metric[values$percentage_present < threshold_new]
 mat_filtered <- mat_filtered[, c(colnames(mat_filtered)[!colnames(mat_filtered) %in% extra_exclude])]
} 


###########
if((dim(mat_filtered))[2] >0){
	pdf(paste0(outputdir, "Plots/ISOTYPER/Metrics_Analysis_Filtered_Matrix", subsampled_depth_all, "_", iso_type,".pdf"), height=6, width=6)
	par(mfrow= c(2,2), mar = c(5,5,3,3))
	nz_rows = apply(mat_filtered, 1, function(x){length(which(!is.na(x)))})
	plot(sort(nz_rows), xlab = "sample rank", ylab = "number of non-NA values", pch = 21, bg = "red", col = "red",main ="run2.1")
	threshold = quantile(nz_rows, 0.1)
	segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
	mat_filtered1 = mat_filtered[which(nz_rows>= threshold),]
	nz_cols = apply(mat_filtered1, 2, function(x){length(which(!is.na(x)))})
	plot(sort(nz_cols), xlab = "feature rank", ylab = "number of non-NA values", pch = 21, bg = "blue", col = "blue",main ="run2.2")
	dev.off()
	print("plot2!")
} else {
	print("No metrics passed filtering")
}
##---------------------------------------------------------------------------------------------
## Appending chain type to columns to ensure we can merge them later

if(chain_vdj %like% "T" & (dim(mat_filtered))[2] >0){
	colnames(mat_filtered) <- paste0(receptor, "_",iso_type, "_", colnames(mat_filtered))
	colnames(mat_filtered) <- gsub(paste0("__", receptor), "", colnames(mat_filtered))
	colnames(mat_filtered) <- gsub("[.]", "_", colnames(mat_filtered))
	print("colnames edited")
} 

if(chain_vdj %like% "BC"| chain_vdj %like% "I"){
	if((dim(mat_filtered))[2] >0){
		colnames(mat_filtered) <- paste0("BCR_", iso_type, "_", colnames(mat_filtered))
		print("colnames edited")
	}
} 

## Code to replace technical replicate with mean of all technical replicates (may need to edit if has diferent name)
print("Collapsing Technical replicate by name")
control_row <-data.frame(mat_filtered[rownames(mat_filtered) %like% "JR1795_1003",])
control_rowid <- rownames(mat_filtered) %like% "JR1795_1003"
new_row <- t(data.frame(colSums(control_row)/dim(control_row)[1]))
rownames(new_row) <- "JR1795_1003_POSITIVE_MEAN"
mat_filtered <- mat_filtered[!control_rowid,]
mat_filtered <- rbind(mat_filtered, new_row)
	
control_row <- mat_filtered[rownames(mat_filtered) %like% "JR_HD_T0",]
control_rowid <- rownames(mat_filtered) %like% "JR_HD_T0"
new_row <- t(data.frame(colSums(control_row)/dim(control_row)[1]))
rownames(new_row) <- "JR_HD_T0_POSITIVE_MEAN"
mat_filtered <- mat_filtered[!control_rowid,]
mat_filtered <- rbind(mat_filtered, new_row)
print("Done")

## Final filtering 
## Remove features relating to pseudo genes 
print("Final Filtering...")
p <- mat_filtered
features = colnames(p)
## exclude psuedogenes!!
features = features[grep("IGHGP", features, invert = T)]
features = features[grep("IGHEP", features, invert = T)]
## remove all features with very few unique values
## This likely means there are a lot of 0 values!
diff_values = apply(p[,features], 2, function(x){length(unique(x[which(is.na(x)==F)]))})
features = intersect(features, names(diff_values)[which(diff_values>=8)])
print("Filtering Performed")
xx <- length(colnames(p))-length(features)
print(paste0(xx, " Features Removed from Matrix"))
feature_columns <-   mat_filtered[,!colnames(mat_filtered) %in% grep("IGHV|TRAV|TRBV|TRGV|TRDV", colnames(mat_filtered), value=TRUE)]
## We also want to check samples for greater than 40% missingness
mat_filtered <- mat_filtered[, c(features)]
## Missingness per sample (must be less than 40%!): 
## AFter this point imputation gets messy.....
missing_threshold <- 0.4*length(feature_columns)
if(any(rowSums(is.na(mat_filtered))>missing_threshold)){
	samples_exclude <- rownames(mat_filtered)[rowSums(is.na(mat_filtered))>missing_threshold]
	mat_filtered <- mat_filtered[!rownames(mat_filtered) %in% samples_exclude,]
} 
### 
print("Done psuedogene filtering, no. unique values filtering and missingness filtering!")
	
## Write the final filtered matric to an outs file - used for weighted network correlation analysis!
write.table(mat_filtered, paste0(outputdir, "Summary/isotyper_metrics_filtered_FINAL_METRICS_", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t", row.names=TRUE)

## We also want to do one with no rownames or headers for tensor decomposition analysis!!!!  
write.table(mat_filtered, paste0(outputdir, "Summary/isotyper_metrics_filtered_FINAL_METRICS_TENSOR_FORMAT_", subsampled_depth_all, "_", iso_type, ".txt"), sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)
print("DONE PART 1")

## get sample and feature order which are needed for tensor analysis
sample_order <- rownames(mat_filtered)
feature_order <- colnames(mat_filtered)

write.table(sample_order, paste0(outputdir, "Summary/isotyper_metrics_filtered_FINAL_METRICS_TENSOR_FORMAT_SAMPLE_ORDER", subsampled_depth_all, "_", iso_type, ".txt"), sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(feature_order, paste0(outputdir, "Summary/isotyper_metrics_filtered_FINAL_METRICS_TENSOR_FORMAT_FEATURE_ORDER", subsampled_depth_all, "_", iso_type, ".txt"), sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)


}