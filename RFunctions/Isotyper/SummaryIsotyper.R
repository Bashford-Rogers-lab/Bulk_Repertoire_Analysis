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
	p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"Isotype"]))),]
	if(class(p)=="character"){
		p <- as.matrix(t(p))
	} 
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
	vrenyi[which(vrenyi<0)] <- -1
	crenyi[which(crenyi<0)] <- -1
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
	} 
	analysis_matrices12 = analysis_matrices
	names(analysis_matrices12) <- analysis_names
	} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices12 <- vector(mode = "list", length = 0)	
	} 	
print("DONE 12")
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
	} 
	analysis_matrices13 = analysis_matrices
	names(analysis_matrices13) <- analysis_names
	} else { 
			print(paste0("File: ", file, " IS EMPTY"))
			analysis_matrices13 <- vector(mode = "list", length = 0)	
	}	
print("DONE 13")
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
## PART 2!!!!!!
## COMPOSING THE OVERALL MATRIX 

if(chain_vdj %like% "BC" | chain_vdj %like% "I"){
	print_info = c(analysis_matrices1, analysis_matrices2,  analysis_matrices3, analysis_matrices4, analysis_matrices5, analysis_matrices5b, analysis_matrices6, analysis_matrices7, analysis_matrices8, analysis_matrices9, analysis_matrices10, analysis_matrices11,  analysis_matrices12, analysis_matrices13)
} else {
	print_info = c(analysis_matrices1, analysis_matrices3, analysis_matrices7, analysis_matrices8, analysis_matrices10, analysis_matrices11, analysis_matrices12, analysis_matrices13)
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
overall_matrix[overall_matrix=="-1"] <- "NA"
## Remove any columns which are empty (e.g. nothing was calculated) 
empty <- apply(overall_matrix, 2, function(x){length(which(x==0 | x=="NA"))})
empty_cols <- empty[empty==dim(overall_matrix)[1]]
overall_matrix<- overall_matrix[, c(!colnames(overall_matrix) %in% names(empty_cols))]
## Remove NA columns 
empty <- apply(overall_matrix, 2, function(x){length(which(x=="NA"))})
empty_cols <- empty[empty==dim(overall_matrix)[1]]
overall_matrix<- overall_matrix[, c(!colnames(overall_matrix) %in% names(empty_cols))]

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
	data[data=="-1"] <- "NA"
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
values$sig[values$missingness < 0.5] <- "warning >50% missingness"
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
nz_cols = apply(mat, 2, function(x){length(which(x!="-1" & !is.na(x)))})
matr1 = mat[,which(nz_cols>= max(nz_cols))]
nz_rows = apply(mat, 1, function(x){length(which(x!="-1" & !is.na(x)))})
matr2 = matr1[which(nz_rows>= max(nz_rows)),]

## Create plots directory if not present
if (!file.exists(paste0(outputdir, "Plots/ISOTYPER"))){
    dir.create(paste0(outputdir, "Plots/ISOTYPER"))
}


###############
# Doing plot 
pdf(paste0(outputdir, "Plots/ISOTYPER/Metrics_Missingness_Summary_", subsampled_depth_all, "_", iso_type,".pdf"), width=6, height=6)
par(mfrow= c(2,2), mar = c(5,5,3,3))
##1
nz_cols = apply(mat, 2, function(x){length(which(x!=-1& !is.na(x)))})
threshold = quantile(nz_cols, 0.25)
plot(sort(nz_cols), xlab = "feature rank", ylab = "number of non-NA values", pch = 21, bg = "blue", col = "blue",main=paste0("run1 threshold: ", threshold))
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat1 = mat[,which(nz_cols>= threshold)]

##2
nz_rows = apply(mat1, 1, function(x){length(which(x!=-1& !is.na(x)))})
threshold = quantile(nz_rows, 0.1)
plot(sort(nz_rows), xlab = "sample rank", ylab = "number of non-NA values", pch = 21, bg = "red", col = "red",main=paste0("run2 threshold: ", threshold))
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat2 = mat1[which(nz_rows>= threshold),]

##3
nz_cols = apply(mat2, 2, function(x){length(which(x!=-1& !is.na(x)))})
threshold = quantile(nz_cols, 0.25)
plot(sort(nz_cols), xlab = "feature rank", ylab = "number of non-NA values", pch = 21, bg = "blue", col = "blue",main=paste0("run3 threshold: ", threshold))
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat3 = mat2[,which(nz_cols>= threshold)]

##4
nz_rows = apply(mat3, 1, function(x){length(which(x!=-1& !is.na(x)))})
threshold = quantile(nz_rows, 0.1)
plot(sort(nz_rows), xlab = "sample rank", ylab = "number of non-NA values", pch = 21, bg = "red", col = "red",main=paste0("run4 threshold: ", threshold))
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat4 = mat3[which(nz_rows>= threshold),]
dev.off()

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
 keep <- keep[!keep %like% "IGHV"]
 keep <- unique(c(keep, annotated_metrics_to_keep_colnames))
 keep <- keep[!keep %in% annotated_metrics_to_keep$Metric[annotated_metrics_to_keep$Include_Metric=="0"]]
 keep <- keep[keep %in% values$Metric[values$percentage_present > threshold_new]]
 keep <- keep[keep %in% colnames(mat)]
 mat_filtered = mat[,c(keep)]
} 

if(chain_vdj %like% "T"){
 metrics_exclude <- grep("SHM", colnames(mat), value=TRUE)
 keepers <- colnames(mat)[!colnames(mat) %in% metrics_exclude]
 mat_filtered = mat[, c(keepers)]
 extra_exclude <- values$Metric[values$percentage_present < 50]
 mat_filtered <- mat_filtered[, c(colnames(mat_filtered)[!colnames(mat_filtered) %in% extra_exclude])]
} 

###########
pdf(paste0(outputdir, "Plots/ISOTYPER/Metrics_Analysis_Filtered_Matrix", subsampled_depth_all, "_", iso_type,".pdf"), height=6, width=6)
par(mfrow= c(2,2), mar = c(5,5,3,3))
nz_rows = apply(mat_filtered, 1, function(x){length(which(x!=-1))})
plot(sort(nz_rows), xlab = "sample rank", ylab = "number of non-NA values", pch = 21, bg = "red", col = "red",main ="run2.1")
threshold = quantile(nz_rows, 0.1)
segments(-10,threshold, 10000, threshold, col = "red", lwd = 2,lty = 2)
mat_filtered1 = mat_filtered[which(nz_rows>= threshold),]
nz_cols = apply(mat_filtered1, 2, function(x){length(which(x!=-1))})
plot(sort(nz_cols), xlab = "feature rank", ylab = "number of non-NA values", pch = 21, bg = "blue", col = "blue",main ="run2.2")
dev.off()

##---------------------------------------------------------------------------------------------
## *** Make correllogram of all the different measures **** 

#Findcolumns with all missing values
data_1 <- data.frame(overall_matrix[, 1:(length(colnames(overall_matrix))-1)])

allmisscols <- apply(data_1,2, function(x)all(is.na(x)));
colswithallmiss <-names(allmisscols[allmisscols>0]);
colswithallmiss <- c(colswithallmiss, "sample")

# Remove columns with all missing (already done but dont bother editing) 
all_miss <- colnames(data_1)[!colnames(data_1) %in% colswithallmiss]
data_1 <- data_1[, c(all_miss)]
bad_columns <- unique(c(colswithallmiss, bad_columns))
new_cols <- colnames(data_1)[!colnames(data_1) %in% bad_columns]
data_2 <- data_1[, c(new_cols)]
bad_columns2 <- unique(c(colswithallmiss, bad_columns2))
new_cols2 <- colnames(data_1)[!colnames(data_1) %in% bad_columns2]
data_3 <- data_1[, c(new_cols2)]
bad_columns3 <- unique(c(colswithallmiss, bad_columns3))
new_cols3 <- colnames(data_1)[!colnames(data_1) %in% bad_columns3]
data_4 <- data_1[, c(new_cols3)]

data_1 <- apply(as.matrix(data_1),2,as.numeric) 
data_2 <- apply(as.matrix(data_2),2,as.numeric)
data_3 <- apply(as.matrix(data_3),2,as.numeric)
data_4 <- apply(as.matrix(data_4),2,as.numeric)


## Final matrix!!
data_mat <- as.matrix(mat_filtered)
data_mat <- apply(as.matrix(data_mat),2,as.numeric)

## Calculate the correlation between matrices 
## because we have nas we have to use pairwise complete obs 
## then the correlation or covariance between each pair of variables is computed using all complete pairs of observations on those variables. 
## This can result in covariance or correlation matrices which are not positive semi-definite, as well as NA entries if there are no complete pairs for that pair of variables.

correlation_matrix <- tryCatch(Hmisc::rcorr(data_1),error=function(e) e, warning=function(w) w)
correlation_matrix2 <- tryCatch(Hmisc::rcorr(data_2),error=function(e) e, warning=function(w) w)
correlation_matrix3 <- tryCatch(Hmisc::rcorr(data_3),error=function(e) e, warning=function(w) w)
correlation_matrix4 <- tryCatch(Hmisc::rcorr(data_4),error=function(e) e, warning=function(w) w)
correlation_matrixmat <- tryCatch(Hmisc::rcorr(data_mat),error=function(e) e, warning=function(w) w)
    
# Calculating widths for plots 
if(length(colnames(data_1)) <= 50){
	width <- 35
	height <- 35
	if(length(colnames(data_1)) <= 20){
			width <- 25
			height <- 25
	}
} else if (length(colnames(data_1)) > 50 & length(colnames(data_1)) < 100){
	width <- 45
	height <- 45
} else {
	width <- 65
	height <- 65
} 

## Added in 'trys' to handle errors 
pdf(paste0(outputdir, "Plots/ISOTYPER/Correlation_between_measures_", subsampled_depth_all, "_", iso_type,".pdf"), width=width, height=height)

if(!("warning" %in% class(correlation_matrix)) && !("error" %in% class(correlation_matrix))){
	try(corrplot(correlation_matrix[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix[[3]], sig.level = 0.05, insig = "blank", tl.col="black", title="No % Present threshold, p<0.05", mar=c(0,0,2,0)))
	try(corrplot(correlation_matrix[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title="No % Present threshold: p<0.01", mar=c(0,0,2,0)))
	try(corrplot(correlation_matrix[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix[[3]], sig.level = 0.001, insig = "blank", tl.col="black", title="No % Present threshold: p<0.001",mar=c(0,0,2,0)))
} 

if(!("warning" %in% class(correlation_matrix4)) && !("error" %in% class(correlation_matrix))){
	try(corrplot(correlation_matrix4[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix4[[3]], sig.level = 0.05, insig = "blank", tl.col="black", title="75% Present threshold, p<0.05", mar=c(0,0,2,0)))
	try(corrplot(correlation_matrix4[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix4[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title="75% Present threshold: p<0.01", mar=c(0,0,2,0)))
	try(corrplot(correlation_matrix4[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix4[[3]], sig.level = 0.001, insig = "blank", tl.col="black", title="75% Present threshold: p<0.001",mar=c(0,0,2,0)))
} 

if(!("warning" %in% class(correlation_matrix2)) && !("error" %in% class(correlation_matrix2))){
	try(corrplot(correlation_matrix2[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix2[[3]], sig.level = 0.05, insig = "blank", tl.col="black", title="50% Present threshold: p<0.05", mar=c(0,0,2,0)))
	try(corrplot(correlation_matrix2[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix2[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title="50% Present threshold: p<0.01", mar=c(0,0,2,0)))
	try(corrplot(correlation_matrix2[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix2[[3]], sig.level = 0.001, insig = "blank", tl.col="black", title="50% Present threshold: p<0.001", mar=c(0,0,2,0)))
}

if(!("warning" %in% class(correlation_matrix3)) && !("error" %in% class(correlation_matrix3))){
	try(corrplot(correlation_matrix3[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix3[[3]], sig.level = 0.05, insig = "blank", tl.col="black", title="25% Present threshold: p<0.05", mar=c(0,0,2,0)))
	try(corrplot(correlation_matrix3[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix3[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title="25% Present threshold: p<0.01", mar=c(0,0,2,0)))
	try(corrplot(correlation_matrix3[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix3[[3]], sig.level = 0.001, insig = "blank", tl.col="black", title="25% Present threshold: p<0.001", mar=c(0,0,2,0)))
}
dev.off()

pdf(paste0(outputdir, "Plots/ISOTYPER/Correlation_between_measures_FINAL_MATRIX_", subsampled_depth_all, "_", iso_type,".pdf"), width=width, height=height)
if(!("warning" %in% class(correlation_matrixmat)) && !("error" %in% class(correlation_matrixmat))){
	try(corrplot(correlation_matrixmat[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixmat[[3]], sig.level = 0.05, insig = "blank", tl.col="black", title="FINAL MATRIX p<0.05", mar=c(0,0,2,0)))
	try(corrplot(correlation_matrixmat[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixmat[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title="FINAL MATRIX p<0.01", mar=c(0,0,2,0)))
	try(corrplot(correlation_matrixmat[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixmat[[3]], sig.level = 0.001, insig = "blank", tl.col="black", title="FINAL MATRIX p<0.001", mar=c(0,0,2,0)))
}
dev.off()

## Appending chain type to columns to ensure we can merge them later

if(chain_vdj %like% "T"){
	colnames(mat_filtered) <- paste0(receptor, "_",iso_type, "_", colnames(mat_filtered))
	colnames(mat_filtered) <- gsub(paste0("__", receptor), "", colnames(mat_filtered))
	colnames(mat_filtered) <- gsub("[.]", "_", colnames(mat_filtered))
} 

if(chain_vdj %like% "BC"| chain_vdj %like% "I"){
	colnames(mat_filtered) <- paste0("BCR_", iso_type, "_", colnames(mat_filtered))
} 

## Write the final filtered matric to an outs file
write.table(mat_filtered, paste0(outputdir, "Summary/isotyper_metrics_filtered_FINAL_METRICS_", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t", row.names=TRUE)
print("DONE PART 1")

######################################################################################################################
######################################################################################################################

######################################################################################################################
## PART 2 - ve gene usage 
## FOR TCR 
if(chain_vdj %like% "T"){
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_grouped_isotype_frequency_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check3 <- FALSE
	} else {
		check3 <- TRUE
	}
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	# remove samples with too low read count
	p=p[which(!as.character(p[,"X.sample"]) %in% samples_to_low_all),]
	
	p <- data.frame(p)
	p$uniq_read_freq <- as.numeric(p$uniq_read_freq)
	p <- p[!p$class %in% c("ALL", "expanded", "IGHD,IGHM_unmutated_singleton", "Class_switched", "IGHD,IGHM_mutated", "IGHD,IGHM_unmutated", "unexpanded", "class"),] 

	if(chain_vdj %like% "T"){
		p <- p[p[, "class"]==receptor,]
	} 
		
	## Calculate a percentage of repertoire which is each read 
	p$percent_repertoire <- NA
	for(i in unique(p[, "X.sample"])){
		sample_id <- i 
		sum_frequency <- sum(p$uniq_read_freq[p$X.sample==sample_id])
		p$percent_repertoire[p$X.sample==sample_id] <- (p$uniq_read_freq[p$X.sample==sample_id])/sum_frequency
	} 
	p$VIso <- paste0(p$class, "_", p$V.gene)	
	p1 <- p[, c("X.sample", "percent_repertoire", "VIso")]

	a <- spread(p1, key = VIso, value = percent_repertoire)
	write.table(a, paste0(outputdir, "Summary/V_Gene_usage", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t", row.names=TRUE)
	samplesxx <- a$X.sample
	a$X.sample <- NULL
	a <- apply(as.matrix(a),2,as.numeric)
	row.names(a) <- samples

	if(chain_vdj %like% "T"){
		colnames(a) <- gsub(paste0(receptor, "_"), "", colnames(a))
	} 

	## identify most common v genes (e.g. those with high presence)
	allmisscolsa <- apply(a,2, function(x)sum(!is.na(x)))
	no_samples <- dim(a)[1]
	present_25 <- names(allmisscolsa[allmisscolsa>(0.25*no_samples)])
	present_50 <- names(allmisscolsa[allmisscolsa>(0.5*no_samples)])
	present_75 <- names(allmisscolsa[allmisscolsa>(0.75*no_samples)])
	present_10 <- names(allmisscolsa[allmisscolsa>(0.1*no_samples)])
	present_90 <- names(allmisscolsa[allmisscolsa>(0.9*no_samples)])

	a25 <- a[, c(present_25)]
	a50 <- a[, c(present_50)]
	a75 <- a[, c(present_75)]
	a10 <- a[, c(present_10)]
	a90 <- a[, c(present_90)]

	correlation_matrixa25 <- tryCatch(Hmisc::rcorr(a25, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa50 <- tryCatch(Hmisc::rcorr(a50, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa75 <- tryCatch(Hmisc::rcorr(a75, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa10 <- tryCatch(Hmisc::rcorr(a10, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa90 <- tryCatch(Hmisc::rcorr(a90, type="pearson"),error=function(e) e, warning=function(w) w)

	rep10 <- round((length(present_10)/dim(a)[2])*100, digits=1)
	rep25 <- round((length(present_25)/dim(a)[2])*100, digits=1)
	rep50 <- round((length(present_50)/dim(a)[2])*100, digits=1)
	rep75 <- round((length(present_75)/dim(a)[2])*100, digits=1)
	rep90 <- round((length(present_90)/dim(a)[2])*100, digits=1)

	pdf(paste0(outputdir, "Plots/ISOTYPER/Correlation_between_measures_V_genes_UNIQUE_READS_", subsampled_depth_all, "_", iso_type,".pdf"), width=20, height=25)
	if(!("warning" %in% class(correlation_matrixa90)) && !("error" %in% class(correlation_matrixa90))){
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<0.001: ", rep90, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-05: ", rep90, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-10: ", rep90, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa75)) && !("error" %in% class(correlation_matrixa75))){
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<0.001: ", rep75, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-05: ", rep75, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-10: ", rep75, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa50)) && !("error" %in% class(correlation_matrixa50))){
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<0.001: ", rep50, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-05: ", rep50, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-10: ", rep50, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa25)) && !("error" %in% class(correlation_matrixa25))){
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<0.001: ", rep25, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-05: ", rep25, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-10: ", rep25, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa10)) && !("error" %in% class(correlation_matrixa10))){
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<0.001: ", rep10, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-05: ", rep10, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-10: ", rep10, "% total V genes"), mar=c(0,0,2,0)))
	}
	dev.off()

	print("DONE V GENE COMPARISON")
	## v family usage 
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_grouped_isotype_frequency_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check3 <- FALSE
	} else {
		check3 <- TRUE
	}
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	# remove samples with too low read count
	p=p[which(!as.character(p[,"X.sample"]) %in% samples_to_low_all),]
	
	p <- data.frame(p)
	p$uniq_read_freq <- as.numeric(p$uniq_read_freq)
	p <- p[!p$class %in% c("ALL", "expanded", "IGHD,IGHM_unmutated_singleton", "Class_switched", "IGHD,IGHM_mutated", "IGHD,IGHM_unmutated", "unexpanded", "class"),] 
	if(chain_vdj %like% "T"){
		p <- p[p[, "class"]==receptor,]
	} 
	v_fam <- data.frame(str_split_fixed(p$V.gene, "-", 2))[,1]
	p$vfam <- v_fam

	## get frequencies
	datax <- c()
	for(i in unique(p[, "X.sample"])){
		sample_id <- i 
		for(c in unique(p$vfam[p$X.sample==sample_id])){
			v_gene <- c
			sum_frequency <- sum(p$uniq_read_freq[p$X.sample==sample_id & p$vfam==v_gene])
			sum_frequency2 <- (sum_frequency / sum(p$uniq_read_freq[p$X.sample==sample_id]))*100
			datarow <- c(sample_id, v_gene, sum_frequency2)
			#print(datarow)
			datax <- rbind(datax, datarow) 
	}
	}
	 
	datax <- data.frame(datax)
	colnames(datax) <- c("Sample", "V_family", "Percent_Repertoire")
	datax$Percent_Repertoire <- as.numeric(datax$Percent_Repertoire)

	a <- spread(datax, key = V_family, value = Percent_Repertoire)
	write.table(a, paste0(outputdir, "Summary/V_Family_usage", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t", row.names=TRUE)
	samplesxx <- a$X.sample
	a$X.sample <- NULL
	a <- apply(as.matrix(a),2,as.numeric)
	row.names(a) <- samples

	if(chain_vdj %like% "T"){
		colnames(a) <- gsub(paste0(receptor, "_"), "", colnames(a))
	} 

	## identify most common v genes (e.g. those with high presence)
	allmisscolsa <- apply(a,2, function(x)sum(!is.na(x)))
	no_samples <- dim(a)[1]
	present_25 <- names(allmisscolsa[allmisscolsa>(0.25*no_samples)])
	present_50 <- names(allmisscolsa[allmisscolsa>(0.5*no_samples)])
	present_75 <- names(allmisscolsa[allmisscolsa>(0.75*no_samples)])
	present_10 <- names(allmisscolsa[allmisscolsa>(0.1*no_samples)])
	present_90 <- names(allmisscolsa[allmisscolsa>(0.9*no_samples)])
	a25 <- a[, c(present_25)]
	a50 <- a[, c(present_50)]
	a75 <- a[, c(present_75)]
	a10 <- a[, c(present_10)]
	a90 <- a[, c(present_90)]
	correlation_matrixa25 <- tryCatch(Hmisc::rcorr(a25, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa50 <- tryCatch(Hmisc::rcorr(a50, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa75 <- tryCatch(Hmisc::rcorr(a75, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa10 <- tryCatch(Hmisc::rcorr(a10, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa90 <- tryCatch(Hmisc::rcorr(a90, type="pearson"),error=function(e) e, warning=function(w) w)
	rep10 <- round((length(present_10)/dim(a)[2])*100, digits=1)
	rep25 <- round((length(present_25)/dim(a)[2])*100, digits=1)
	rep50 <- round((length(present_50)/dim(a)[2])*100, digits=1)
	rep75 <- round((length(present_75)/dim(a)[2])*100, digits=1)
	rep90 <- round((length(present_90)/dim(a)[2])*100, digits=1)
	pdf(paste0(outputdir, "Plots/ISOTYPER/Correlation_between_measures_V_families_UNIQUE_READS_", subsampled_depth_all, "_", iso_type,".pdf"), width=20, height=25)
	if(!("warning" %in% class(correlation_matrixa90)) && !("error" %in% class(correlation_matrixa90))){
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<0.001: ", rep90, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-05: ", rep90, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-10: ", rep90, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa75)) && !("error" %in% class(correlation_matrixa75))){
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<0.001: ", rep75, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-05: ", rep75, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-10: ", rep75, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa50)) && !("error" %in% class(correlation_matrixa50))){
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<0.001: ", rep50, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-05: ", rep50, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-10: ", rep50, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa25)) && !("error" %in% class(correlation_matrixa25))){
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<0.001: ", rep25, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-05: ", rep25, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-10: ", rep25, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa10)) && !("error" %in% class(correlation_matrixa10))){
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<0.001: ", rep10, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-05: ", rep10, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-10: ", rep10, "% total V genes"), mar=c(0,0,2,0)))
	}
	dev.off()

	print("DONE V FAMILY COMPARISON")
	## J gene usage 
	###################################################
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_J_gene_grouped_isotype_frequency_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check3 <- FALSE
	} else {
		check3 <- TRUE
	}
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	# remove samples with too low read count
	p=p[which(!as.character(p[,"X.sample"]) %in% samples_to_low_all),]
	p <- data.frame(p)
	p$uniq_read_freq <- as.numeric(p$uniq_read_freq)
	p <- p[!p$class %in% c("ALL", "expanded", "IGHD,IGHM_unmutated_singleton", "Class_switched", "IGHD,IGHM_mutated", "IGHD,IGHM_unmutated", "unexpanded", "class"),] 

	if(chain_vdj %like% "T"){
		p <- p[p[, "class"]==receptor,]
	} 
		
	## Calculate a percentage of repertoire which is each read 
	p$percent_repertoire <- NA
	for(i in unique(p[, "X.sample"])){
		sample_id <- i 
		sum_frequency <- sum(p$uniq_read_freq[p$X.sample==sample_id])
		p$percent_repertoire[p$X.sample==sample_id] <- (p$uniq_read_freq[p$X.sample==sample_id])/sum_frequency
	} 
	p$VIso <- paste0(p$class, "_", p$J.gene)	
	p1 <- p[, c("X.sample", "percent_repertoire", "VIso")]

	a <- spread(p1, key = VIso, value = percent_repertoire)
	write.table(a, paste0(outputdir, "Summary/J_Gene_usage", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t", row.names=TRUE)
	samplesxx <- a$X.sample
	a$X.sample <- NULL
	a <- apply(as.matrix(a),2,as.numeric)
	row.names(a) <- samples

	if(chain_vdj %like% "T"){
		colnames(a) <- gsub(paste0(receptor, "_"), "", colnames(a))
	} 


	## identify most common v genes (e.g. those with high presence)
	allmisscolsa <- apply(a,2, function(x)sum(!is.na(x)))
	no_samples <- dim(a)[1]
	present_25 <- names(allmisscolsa[allmisscolsa>(0.25*no_samples)])
	present_50 <- names(allmisscolsa[allmisscolsa>(0.5*no_samples)])
	present_75 <- names(allmisscolsa[allmisscolsa>(0.75*no_samples)])
	present_10 <- names(allmisscolsa[allmisscolsa>(0.1*no_samples)])
	present_90 <- names(allmisscolsa[allmisscolsa>(0.9*no_samples)])

	a25 <- a[, c(present_25)]
	a50 <- a[, c(present_50)]
	a75 <- a[, c(present_75)]
	a10 <- a[, c(present_10)]
	a90 <- a[, c(present_90)]

	correlation_matrixa25 <- tryCatch(Hmisc::rcorr(a25, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa50 <- tryCatch(Hmisc::rcorr(a50, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa75 <- tryCatch(Hmisc::rcorr(a75, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa10 <- tryCatch(Hmisc::rcorr(a10, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa90 <- tryCatch(Hmisc::rcorr(a90, type="pearson"),error=function(e) e, warning=function(w) w)

	rep10 <- round((length(present_10)/dim(a)[2])*100, digits=1)
	rep25 <- round((length(present_25)/dim(a)[2])*100, digits=1)
	rep50 <- round((length(present_50)/dim(a)[2])*100, digits=1)
	rep75 <- round((length(present_75)/dim(a)[2])*100, digits=1)
	rep90 <- round((length(present_90)/dim(a)[2])*100, digits=1)

	pdf(paste0(outputdir, "Plots/ISOTYPER/Correlation_between_measures_J_genes_UNIQUE_READS_", subsampled_depth_all, "_", iso_type,".pdf"), width=20, height=25)
	if(!("warning" %in% class(correlation_matrixa90)) && !("error" %in% class(correlation_matrixa90))){
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<0.001: ", rep90, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-05: ", rep90, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-10: ", rep90, "% total J genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa75)) && !("error" %in% class(correlation_matrixa75))){
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<0.001: ", rep75, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-05: ", rep75, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-10: ", rep75, "% total J genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa50)) && !("error" %in% class(correlation_matrixa50))){
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<0.001: ", rep50, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-05: ", rep50, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-10: ", rep50, "% total J genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa25)) && !("error" %in% class(correlation_matrixa25))){
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<0.001: ", rep25, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-05: ", rep25, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-10: ", rep25, "% total J genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa10)) && !("error" %in% class(correlation_matrixa10))){
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<0.001: ", rep10, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-05: ", rep10, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-10: ", rep10, "% total J genes"), mar=c(0,0,2,0)))
	}
	dev.off()
	print("DONE J GENE COMPARISON")
	###############################################################################
	### VJ gene usage 
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check3 <- FALSE
	} else {
		check3 <- TRUE
	}
	q <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	q <- data.frame(q)
		
	# remove samples with too low read count
	q=q[which(!as.character(q[,"X.ID"]) %in% samples_to_low_all),]
	
	q$number.of.sequences <- as.numeric(q$number.of.sequences)
	q <- q[q$classification %in% c(unique(p$class)),] 
	q$X.ID <- gsub("BCR_", "", q$X.ID)
	q$X.ID <- gsub("TCR_", "", q$X.ID)
	q$X.ID <- gsub("TCRA_", "", q$X.ID)
	q$X.ID <- gsub("TCRB_", "", q$X.ID)
	q$X.ID <- gsub("TCRG_", "", q$X.ID)
	q$X.ID <- gsub("TCRD_", "", q$X.ID)

	q$X.ID <- gsub("_productive", "", q$X.ID)
	q$X.ID <- gsub("_unproductive", "", q$X.ID)
	q$VJISO <- paste0(q$classification, "_", q$VJ)
	## Calculate a percentage of repertoire which is each read 
	q$percent_repertoire <- NA
	for(i in unique(q[, "X.ID"])){
		sample_id <- i 
		sum_frequency <- sum(q$number.of.sequences[q$X.ID==sample_id])
		q$percent_repertoire[q$X.ID==sample_id] <- (q$number.of.sequences[q$X.ID==sample_id])/sum_frequency
	} 
	q1 <- q[, c("X.ID", "percent_repertoire", "VJISO")]
	b <- spread(q1, key = VJISO, value = percent_repertoire)
	write.table(b, paste0(outputdir, "Summary/VJ_Gene_usage", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t", row.names=TRUE)
	
	a <- b
	samplesxx <- a$X.sample
	a$X.sample <- NULL
	a <- apply(as.matrix(a),2,as.numeric)
	row.names(a) <- samples

	## identify most common v genes (e.g. those with high presence)
	allmisscolsa <- apply(a,2, function(x)sum(!is.na(x)))
	no_samples <- dim(a)[1]
	present_25 <- names(allmisscolsa[allmisscolsa>(0.25*no_samples)])
	present_50 <- names(allmisscolsa[allmisscolsa>(0.5*no_samples)])
	present_75 <- names(allmisscolsa[allmisscolsa>(0.75*no_samples)])
	present_10 <- names(allmisscolsa[allmisscolsa>(0.1*no_samples)])
	present_90 <- names(allmisscolsa[allmisscolsa>(0.9*no_samples)])

	a25 <- a[, c(present_25)]
	a50 <- a[, c(present_50)]
	a75 <- a[, c(present_75)]
	a10 <- a[, c(present_10)]
	a90 <- a[, c(present_90)]

	correlation_matrixa25 <- tryCatch(Hmisc::rcorr(a25, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa50 <- tryCatch(Hmisc::rcorr(a50, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa75 <- tryCatch(Hmisc::rcorr(a75, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa10 <- tryCatch(Hmisc::rcorr(a10, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa90 <- tryCatch(Hmisc::rcorr(a90, type="pearson"),error=function(e) e, warning=function(w) w)

	rep10 <- round((length(present_10)/dim(a)[2])*100, digits=1)
	rep25 <- round((length(present_25)/dim(a)[2])*100, digits=1)
	rep50 <- round((length(present_50)/dim(a)[2])*100, digits=1)
	rep75 <- round((length(present_75)/dim(a)[2])*100, digits=1)
	rep90 <- round((length(present_90)/dim(a)[2])*100, digits=1)

	pdf(paste0(outputdir, "Plots/ISOTYPER/Correlation_between_measures_VJ_genes", subsampled_depth_all, "_", iso_type,".pdf"), width=40, height=45)
	if(!("warning" %in% class(correlation_matrixa90)) && !("error" %in% class(correlation_matrixa90))){
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<0.001: ", rep90, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-05: ", rep90, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-10: ", rep90, "% total VJ genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa75)) && !("error" %in% class(correlation_matrixa75))){
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<0.001: ", rep75, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-05: ", rep75, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-10: ", rep75, "% total VJ genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa50)) && !("error" %in% class(correlation_matrixa50))){
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<0.001: ", rep50, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-05: ", rep50, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-10: ", rep50, "% total VJ genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa25)) && !("error" %in% class(correlation_matrixa25))){
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<0.001: ", rep25, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-05: ", rep25, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-10: ", rep25, "% total VJ genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa10)) && !("error" %in% class(correlation_matrixa10))){
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<0.001: ", rep10, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-05: ", rep10, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-10: ", rep10, "% total VJ genes"), mar=c(0,0,2,0)))
	}
	dev.off()
	print("DONE VJ GENE COMPARISON")
	########
} 

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

if(chain_vdj %like% "BC"| chain_vdj %like% "I"){
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_grouped_isotype_frequency_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check3 <- FALSE
	} else {
		check3 <- TRUE
	}
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	# remove samples with too low read count
	p=p[which(!as.character(p[,"X.sample"]) %in% samples_to_low_all),]
	
	p <- data.frame(p)
	p$uniq_read_freq <- as.numeric(p$uniq_read_freq)
	p <- p[!p$class %in% c("ALL", "expanded", "IGHD,IGHM_unmutated_singleton", "Class_switched", "IGHD,IGHM_mutated", "IGHD,IGHM_unmutated", "unexpanded", "class"),] 
	
	## Collapse for V gene Usage (ignore isotype)
	p1 <- aggregate(p$uniq_read_freq, by=list(p$X.sample, p$V.gene), FUN=sum)
	colnames(p1) <- c("X.sample", "V.gene", "uniq_read_freq")
	old_p <- p 
	p <- p1
	
	## Calculate a percentage of repertoire which is each read 
	p$percent_repertoire <- NA
	for(i in unique(p[, "X.sample"])){
		sample_id <- i 
		sum_frequency <- sum(p$uniq_read_freq[p$X.sample==sample_id])
		p$percent_repertoire[p$X.sample==sample_id] <- (p$uniq_read_freq[p$X.sample==sample_id])/sum_frequency
	} 
	p1 <- p[, c("X.sample", "percent_repertoire", "V.gene")]

	a <- spread(p1, key = V.gene, value = percent_repertoire)
	write.table(a, paste0(outputdir, "Summary/V_Gene_usage", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t", row.names=TRUE)
	samplesxx <- a$X.sample
	a$X.sample <- NULL
	a <- apply(as.matrix(a),2,as.numeric)
	row.names(a) <- samples

	## identify most common v genes (e.g. those with high presence)
	allmisscolsa <- apply(a,2, function(x)sum(!is.na(x)))
	no_samples <- dim(a)[1]
	present_25 <- names(allmisscolsa[allmisscolsa>(0.25*no_samples)])
	present_50 <- names(allmisscolsa[allmisscolsa>(0.5*no_samples)])
	present_75 <- names(allmisscolsa[allmisscolsa>(0.75*no_samples)])
	present_10 <- names(allmisscolsa[allmisscolsa>(0.1*no_samples)])
	present_90 <- names(allmisscolsa[allmisscolsa>(0.9*no_samples)])

	a25 <- a[, c(present_25)]
	a50 <- a[, c(present_50)]
	a75 <- a[, c(present_75)]
	a10 <- a[, c(present_10)]
	a90 <- a[, c(present_90)]

	correlation_matrixa25 <- tryCatch(Hmisc::rcorr(a25, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa50 <- tryCatch(Hmisc::rcorr(a50, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa75 <- tryCatch(Hmisc::rcorr(a75, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa10 <- tryCatch(Hmisc::rcorr(a10, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa90 <- tryCatch(Hmisc::rcorr(a90, type="pearson"),error=function(e) e, warning=function(w) w)

	rep10 <- round((length(present_10)/dim(a)[2])*100, digits=1)
	rep25 <- round((length(present_25)/dim(a)[2])*100, digits=1)
	rep50 <- round((length(present_50)/dim(a)[2])*100, digits=1)
	rep75 <- round((length(present_75)/dim(a)[2])*100, digits=1)
	rep90 <- round((length(present_90)/dim(a)[2])*100, digits=1)

	pdf(paste0(outputdir, "Plots/ISOTYPER/Correlation_between_measures_V_genes_UNIQUEREADS_", subsampled_depth_all, "_", iso_type,".pdf"), width=20, height=25)
	if(!("warning" %in% class(correlation_matrixa90)) && !("error" %in% class(correlation_matrixa90))){
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<0.001: ", rep90, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-05: ", rep90, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-10: ", rep90, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa75)) && !("error" %in% class(correlation_matrixa75))){
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<0.001: ", rep75, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-05: ", rep75, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-10: ", rep75, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa50)) && !("error" %in% class(correlation_matrixa50))){
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<0.001: ", rep50, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-05: ", rep50, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-10: ", rep50, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa25)) && !("error" %in% class(correlation_matrixa25))){
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<0.001: ", rep25, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-05: ", rep25, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-10: ", rep25, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa10)) && !("error" %in% class(correlation_matrixa10))){
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<0.001: ", rep10, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-05: ", rep10, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-10: ", rep10, "% total V genes"), mar=c(0,0,2,0)))
	}
	dev.off()

	print("DONE V GENE COMPARISON")
	## v family usage 
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_grouped_isotype_frequency_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check3 <- FALSE
	} else {
		check3 <- TRUE
	}
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	# remove samples with too low read count
	p=p[which(!as.character(p[,"X.sample"]) %in% samples_to_low_all),]

	p <- data.frame(p)
	p$uniq_read_freq <- as.numeric(p$uniq_read_freq)
	p <- p[!p$class %in% c("ALL", "expanded", "IGHD,IGHM_unmutated_singleton", "Class_switched", "IGHD,IGHM_mutated", "IGHD,IGHM_unmutated", "unexpanded", "class"),] 
	
	## Collapse for V gene Usage (ignore isotype)
	p1 <- aggregate(p$uniq_read_freq, by=list(p$X.sample, p$V.gene), FUN=sum)
	colnames(p1) <- c("X.sample", "V.gene", "uniq_read_freq")
	old_p <- p 
	p <- p1
	
	v_fam <- data.frame(str_split_fixed(p$V.gene, "-", 2))[,1]
	p$vfam <- v_fam

	## get frequencies
	datax <- c()
	for(i in unique(p[, "X.sample"])){
		sample_id <- i 
		for(c in unique(p$vfam[p$X.sample==sample_id])){
			v_gene <- c
			sum_frequency <- sum(p$uniq_read_freq[p$X.sample==sample_id & p$vfam==v_gene])
			sum_frequency2 <- (sum_frequency / sum(p$uniq_read_freq[p$X.sample==sample_id]))*100
			datarow <- c(sample_id, v_gene, sum_frequency2)
			#print(datarow)
			datax <- rbind(datax, datarow) 
	}
	}
	 
	datax <- data.frame(datax)
	colnames(datax) <- c("Sample", "V_family", "Percent_Repertoire")
	datax$Percent_Repertoire <- as.numeric(datax$Percent_Repertoire)

	a <- spread(datax, key = V_family, value = Percent_Repertoire)
	write.table(a, paste0(outputdir, "Summary/V_Family_usage", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t", row.names=TRUE)
	samplesxx <- a$X.sample
	a$X.sample <- NULL
	a <- apply(as.matrix(a),2,as.numeric)
	row.names(a) <- samples

	## identify most common v genes (e.g. those with high presence)
	allmisscolsa <- apply(a,2, function(x)sum(!is.na(x)))
	no_samples <- dim(a)[1]
	present_25 <- names(allmisscolsa[allmisscolsa>(0.25*no_samples)])
	present_50 <- names(allmisscolsa[allmisscolsa>(0.5*no_samples)])
	present_75 <- names(allmisscolsa[allmisscolsa>(0.75*no_samples)])
	present_10 <- names(allmisscolsa[allmisscolsa>(0.1*no_samples)])
	present_90 <- names(allmisscolsa[allmisscolsa>(0.9*no_samples)])
	a25 <- a[, c(present_25)]
	a50 <- a[, c(present_50)]
	a75 <- a[, c(present_75)]
	a10 <- a[, c(present_10)]
	a90 <- a[, c(present_90)]
	correlation_matrixa25 <- tryCatch(Hmisc::rcorr(a25, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa50 <- tryCatch(Hmisc::rcorr(a50, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa75 <- tryCatch(Hmisc::rcorr(a75, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa10 <- tryCatch(Hmisc::rcorr(a10, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa90 <- tryCatch(Hmisc::rcorr(a90, type="pearson"),error=function(e) e, warning=function(w) w)
	rep10 <- round((length(present_10)/dim(a)[2])*100, digits=1)
	rep25 <- round((length(present_25)/dim(a)[2])*100, digits=1)
	rep50 <- round((length(present_50)/dim(a)[2])*100, digits=1)
	rep75 <- round((length(present_75)/dim(a)[2])*100, digits=1)
	rep90 <- round((length(present_90)/dim(a)[2])*100, digits=1)
	pdf(paste0(outputdir, "Plots/ISOTYPER/Correlation_between_measures_V_families_UNIQUEREADS_", subsampled_depth_all, "_", iso_type,".pdf"), width=10, height=10)
	if(!("warning" %in% class(correlation_matrixa90)) && !("error" %in% class(correlation_matrixa90))){
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<0.001: ", rep90, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-05: ", rep90, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-10: ", rep90, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa75)) && !("error" %in% class(correlation_matrixa75))){
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<0.001: ", rep75, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-05: ", rep75, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-10: ", rep75, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa50)) && !("error" %in% class(correlation_matrixa50))){
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<0.001: ", rep50, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-05: ", rep50, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-10: ", rep50, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa25)) && !("error" %in% class(correlation_matrixa25))){
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<0.001: ", rep25, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-05: ", rep25, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-10: ", rep25, "% total V genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa10)) && !("error" %in% class(correlation_matrixa10))){
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<0.001: ", rep10, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-05: ", rep10, "% total V genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-10: ", rep10, "% total V genes"), mar=c(0,0,2,0)))
	}
	dev.off()

	print("DONE V FAMILY GENE COMPARISON")
	## J gene usage 
	###################################################
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_J_gene_grouped_isotype_frequency_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check3 <- FALSE
	} else {
		check3 <- TRUE
	}
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	# remove samples with too low read count
	p=p[which(!as.character(p[,"X.sample"]) %in% samples_to_low_all),]
	
	p <- data.frame(p)
	p$uniq_read_freq <- as.numeric(p$uniq_read_freq)
	p <- p[!p$class %in% c("ALL", "expanded", "IGHD,IGHM_unmutated_singleton", "Class_switched", "IGHD,IGHM_mutated", "IGHD,IGHM_unmutated", "unexpanded", "class"),] 
		
	p1 <- aggregate(p$uniq_read_freq, by=list(p$X.sample, p$J.gene), FUN=sum)
	colnames(p1) <- c("X.sample", "J.gene", "uniq_read_freq")
	old_p <- p 
	p <- p1
	
	## Calculate a percentage of repertoire which is each read 
	p$percent_repertoire <- NA
	for(i in unique(p[, "X.sample"])){
		sample_id <- i 
		sum_frequency <- sum(p$uniq_read_freq[p$X.sample==sample_id])
		p$percent_repertoire[p$X.sample==sample_id] <- (p$uniq_read_freq[p$X.sample==sample_id])/sum_frequency
	} 

	p1 <- p[, c("X.sample", "percent_repertoire", "J.gene")]

	a <- spread(p1, key = J.gene, value = percent_repertoire)
	write.table(a, paste0(outputdir, "Summary/J_Gene_usage", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t", row.names=TRUE)
	samplesxx <- a$X.sample
	a$X.sample <- NULL
	a <- apply(as.matrix(a),2,as.numeric)
	row.names(a) <- samples

	## identify most common v genes (e.g. those with high presence)
	allmisscolsa <- apply(a,2, function(x)sum(!is.na(x)))
	no_samples <- dim(a)[1]
	present_25 <- names(allmisscolsa[allmisscolsa>(0.25*no_samples)])
	present_50 <- names(allmisscolsa[allmisscolsa>(0.5*no_samples)])
	present_75 <- names(allmisscolsa[allmisscolsa>(0.75*no_samples)])
	present_10 <- names(allmisscolsa[allmisscolsa>(0.1*no_samples)])
	present_90 <- names(allmisscolsa[allmisscolsa>(0.9*no_samples)])

	a25 <- a[, c(present_25)]
	a50 <- a[, c(present_50)]
	a75 <- a[, c(present_75)]
	a10 <- a[, c(present_10)]
	a90 <- a[, c(present_90)]

	correlation_matrixa25 <- tryCatch(Hmisc::rcorr(a25, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa50 <- tryCatch(Hmisc::rcorr(a50, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa75 <- tryCatch(Hmisc::rcorr(a75, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa10 <- tryCatch(Hmisc::rcorr(a10, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa90 <- tryCatch(Hmisc::rcorr(a90, type="pearson"),error=function(e) e, warning=function(w) w)

	rep10 <- round((length(present_10)/dim(a)[2])*100, digits=1)
	rep25 <- round((length(present_25)/dim(a)[2])*100, digits=1)
	rep50 <- round((length(present_50)/dim(a)[2])*100, digits=1)
	rep75 <- round((length(present_75)/dim(a)[2])*100, digits=1)
	rep90 <- round((length(present_90)/dim(a)[2])*100, digits=1)

	pdf(paste0(outputdir, "Plots/ISOTYPER/Correlation_between_measures_J_genes_UNIQUEREADS_", subsampled_depth_all, "_", iso_type,".pdf"), width=10, height=10)
	if(!("warning" %in% class(correlation_matrixa90)) && !("error" %in% class(correlation_matrixa90))){
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<0.001: ", rep90, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-05: ", rep90, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-10: ", rep90, "% total J genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa75)) && !("error" %in% class(correlation_matrixa75))){
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<0.001: ", rep75, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-05: ", rep75, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-10: ", rep75, "% total J genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa50)) && !("error" %in% class(correlation_matrixa50))){
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<0.001: ", rep50, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-05: ", rep50, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-10: ", rep50, "% total J genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa25)) && !("error" %in% class(correlation_matrixa25))){
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<0.001: ", rep25, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-05: ", rep25, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-10: ", rep25, "% total J genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa10)) && !("error" %in% class(correlation_matrixa10))){
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<0.001: ", rep10, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-05: ", rep10, "% total J genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-10: ", rep10, "% total J genes"), mar=c(0,0,2,0)))
	}
	dev.off()
	print("DONE J GENE COMPARISON")
	
	
	###############################################################################
	### VJ gene usage 
	file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_per_cluster_VJ_gene_usage_by_cluster_classification_", iso_type, ".txt")
	subsample_identifier <- grep("SAMPLED", file, value=TRUE)
	if(length(subsample_identifier)==0){
		check3 <- FALSE
	} else {
		check3 <- TRUE
	}
	q <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	q <- data.frame(q)
	q=q[which(as.character(q[,"X.ID"]) %in% ids_all),]
	# remove samples with too low read count
	q=q[which(!as.character(q[,"X.ID"]) %in% samples_to_low_all),]
	
	q <- data.frame(q)
	q$number.of.sequences <- as.numeric(q$number.of.sequences)
	q <- q[!q$classification %in% c("ALL", "expanded", "IGHD,IGHM_unmutated_singleton", "Class_switched", "IGHD,IGHM_mutated", "IGHD,IGHM_unmutated", "unexpanded", "class"),] 			
	q$number.of.sequences <- as.numeric(q$number.of.sequences)
	q$X.ID <- gsub("BCR_", "", q$X.ID)
	q$X.ID <- gsub("TCR_", "", q$X.ID)
	q$X.ID <- gsub("TCRA_", "", q$X.ID)
	q$X.ID <- gsub("TCRB_", "", q$X.ID)
	q$X.ID <- gsub("TCRG_", "", q$X.ID)
	q$X.ID <- gsub("TCRD_", "", q$X.ID)
	q$X.ID <- gsub("_productive", "", q$X.ID)
	q$X.ID <- gsub("_unproductive", "", q$X.ID)
	
	q1 <- aggregate(q$number.of.sequences, by=list(q$X.ID, q$VJ), FUN=sum)
	colnames(q1) <- c("X.sample", "VJ.gene", "frequency")
	old_q <- q 
	q <- q1

	## Calculate a percentage of repertoire which is each read 
	q$percent_repertoire <- NA
	for(i in unique(q[, "X.sample"])){
		sample_id <- i 
		sum_frequency <- sum(q$frequency[q$X.sample==sample_id])
		q$percent_repertoire[q$X.sample==sample_id] <- (q$frequency[q$X.sample==sample_id])/sum_frequency
	} 
	
	q1 <- q[, c("X.sample", "percent_repertoire", "VJ.gene")]
	b <- spread(q1, key = VJ.gene, value = percent_repertoire)
	write.table(b, paste0(outputdir, "Summary/VJ_Gene_usage", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t", row.names=TRUE)
	
	a <- b
	samplesxx <- a$X.sample
	a$X.sample <- NULL
	a <- apply(as.matrix(a),2,as.numeric)
	row.names(a) <- samples

	## identify most common v genes (e.g. those with high presence)
	allmisscolsa <- apply(a,2, function(x)sum(!is.na(x)))
	no_samples <- dim(a)[1]
	present_25 <- names(allmisscolsa[allmisscolsa>(0.25*no_samples)])
	present_50 <- names(allmisscolsa[allmisscolsa>(0.5*no_samples)])
	present_75 <- names(allmisscolsa[allmisscolsa>(0.75*no_samples)])
	present_10 <- names(allmisscolsa[allmisscolsa>(0.1*no_samples)])
	present_90 <- names(allmisscolsa[allmisscolsa>(0.9*no_samples)])

	a25 <- a[, c(present_25)]
	a50 <- a[, c(present_50)]
	a75 <- a[, c(present_75)]
	a10 <- a[, c(present_10)]
	a90 <- a[, c(present_90)]

	correlation_matrixa25 <- tryCatch(Hmisc::rcorr(a25, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa50 <- tryCatch(Hmisc::rcorr(a50, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa75 <- tryCatch(Hmisc::rcorr(a75, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa10 <- tryCatch(Hmisc::rcorr(a10, type="pearson"),error=function(e) e, warning=function(w) w)
	correlation_matrixa90 <- tryCatch(Hmisc::rcorr(a90, type="pearson"),error=function(e) e, warning=function(w) w)

	rep10 <- round((length(present_10)/dim(a)[2])*100, digits=1)
	rep25 <- round((length(present_25)/dim(a)[2])*100, digits=1)
	rep50 <- round((length(present_50)/dim(a)[2])*100, digits=1)
	rep75 <- round((length(present_75)/dim(a)[2])*100, digits=1)
	rep90 <- round((length(present_90)/dim(a)[2])*100, digits=1)

	pdf(paste0(outputdir, "Plots/ISOTYPER/Correlation_between_measures_VJ_genes", subsampled_depth_all, "_", iso_type,".pdf"), width=40, height=45)
	if(!("warning" %in% class(correlation_matrixa90)) && !("error" %in% class(correlation_matrixa90))){
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<0.001: ", rep90, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-05: ", rep90, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa90[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa90[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >90% of Individuals:  p<1e-10: ", rep90, "% total VJ genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa75)) && !("error" %in% class(correlation_matrixa75))){
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<0.001: ", rep75, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-05: ", rep75, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa75[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa75[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >75% of Individuals:  p<1e-10: ", rep75, "% total VJ genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa50)) && !("error" %in% class(correlation_matrixa50))){
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<0.001: ", rep50, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-05: ", rep50, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa50[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa50[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >50% of Individuals:  p<1e-10: ", rep50, "% total VJ genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa25)) && !("error" %in% class(correlation_matrixa25))){
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<0.001: ", rep25, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-05: ", rep25, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa25[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa25[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >25% of Individuals:  p<1e-10: ", rep25, "% total VJ genes"), mar=c(0,0,2,0)))
	}
	if(!("warning" %in% class(correlation_matrixa10)) && !("error" %in% class(correlation_matrixa10))){
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<0.001: ", rep10, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.00001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-05: ", rep10, "% total VJ genes"), mar=c(0,0,2,0)))
		try(corrplot(correlation_matrixa10[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrixa10[[3]], sig.level = 0.0000000001, insig = "blank", tl.col="black", title=paste0("Present in >10% of Individuals:  p<1e-10: ", rep10, "% total VJ genes"), mar=c(0,0,2,0)))
	}
	dev.off()
	########
	print("DONE VJ GENE COMPARISON")
} 

}

## THE END FINALLY :D 











