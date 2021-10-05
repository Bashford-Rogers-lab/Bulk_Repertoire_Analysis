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

## Function 
summary_isotyper <- function(outputdir, samplesfilepost, iso_type){

#Getting SampleIDs
ids <- read.delim(samplesfilepost, sep='\t', header=FALSE)
ids <- as.character(ids$V1)
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

mins <- c()

if(iso_type == "UNPRODUCTIVE" | iso_type=="PRODUCTIVE"){
	order_samples <- c()
		for(i in 1:length(samples)){
			a <- read.delim(samples[i], header=FALSE)
			a <- a$V2
			a_min <- length(a)
			sampleid <- samples[i]
			order_samples <- c(order_samples, sampleid)
			mins <- c(mins, a_min)
		} 
}
	
if(iso_type=="ALL"){
	order_samples <- c()
		for(i in 1:length(samples)){
			a <- read.delim(samples[i], header=FALSE)
			a <- a[(a$V3=="productive (see comment)" | a$V3=="productive" | a$V3=="unproductive (see comment)" | a$V3=="unproductive"), ]
			# a[a$V3=="rearranged sequence (but no junction found) (see comment)",]
			a <- a$V2
			a_min <- length(a)
			sampleid <- samples[i]
			order_samples <- c(order_samples, sampleid)
			mins <- c(mins, a_min)
		} 
	}
	
depths <- data.frame(cbind(order_samples, mins))
depths$order_samples <- gsub(paste0(path, "/IMGT_BCR_"), "", depths$order_samples)
depths$order_samples <- gsub(paste0(path, "/IMGT_TCRA_"), "", depths$order_samples)
depths$order_samples <- gsub(paste0(path, "/IMGT_TCRB_"), "", depths$order_samples)
depths$order_samples <- gsub(paste0(path, "/IMGT_TCRG_"), "", depths$order_samples)
depths$order_samples <- gsub(paste0(path, "/IMGT_TCRD_"), "", depths$order_samples)
depths$order_samples <- gsub(".txt", "", depths$order_samples)
depths$order_samples <- gsub("_1_Summary", "", depths$order_samples)		
colnames(depths) <- c("SampleIDforDepths", "ReadDepth")
depths$ReadDepth <- as.character(depths$ReadDepth)
depths$ReadDepth <- as.numeric(depths$ReadDepth)

## SampleReaddepth Dataframe 
read_depths_all <- depths

## Getting subsample depths which were used for isotyper script 
counts_used <- paste0(outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS")
all_files <- list.files(counts_used, full.name=TRUE)
all_files <- grep("depth_per_isotype", all_files, value=TRUE)
#need to edit out!
counts_used <- read.delim(all_files[2], sep="\t", header=TRUE)
counts_used <- counts_used[counts_used$type=="UNIQ",]
subsampled_depth_all <- counts_used$min[counts_used$X.isotype=="all"]
		
##Begining Data compiling:
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------
## File Number 1: Cluster Summary: 
## SUBSAMPLED FILE!!!!!!!!!!!!!!!!!!!

file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_per_cluster_network_parameters_SUMMARY_SUBSAMPLED_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check1 <- FALSE
} else {
	check1 <- TRUE
}
p <- as.matrix(read.delim(file, head=TRUE, sep="\t"))
p=p[which(as.character(p[,"X.Id"]) %in% ids_all),]
p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"Isotype"]))),]

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
#cgini_reads =  as.numeric(p[,"cgini_clone_sampling"])
mean_vertex_size =  as.numeric(p[,"mean_vertex_size"])
max_clust_size =  as.numeric(p[,"max_clust_pop"])
max_vertex_size =  as.numeric(p[,"max_vertex_pop"])

#Calculating Normalised V/C Renyi Scores 
vrenyi = 1-(as.numeric(p[,"Vertex.Reyni"])/log(as.numeric(p[,"N.vertices"])))
crenyi =  1-(as.numeric(p[,"Cluster_Renyi"])/log(as.numeric(p[,"N_clusters"])))
vrenyi[which(vrenyi<0)] == -1
crenyi[which(crenyi<0)] = -1

# Fill Data frame 
class_reads = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
v_gini = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
c_gini = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
c_renyi = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
v_renyi = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))

m_mean_vertex_size = matrix(data = -1, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
m_max_clust_size = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))
m_max_vertex_size = matrix(data = 0, nrow = length(ids_all),ncol = length(classes), dimnames=c(list(ids_all), list(classes)))

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
		}
	}
}
q <- list(v_gini, c_gini, m_mean_vertex_size, m_max_clust_size, m_max_vertex_size, c_renyi, v_renyi)
names(q) <- c("v_gini", "c_gini", "m_mean_vertex_size", "m_max_clust_size", "m_max_vertex_size", "c_renyi_normalised", "v_renyi_normalised")

#renaming column names 
for(i in 1:length(q)){
	names <- names(q[i])
	colnames(q[[i]]) <- paste0(names, "_Subsampled", "__", colnames(q[[i]]))
}
names(q) <- c("Vertex_Gini_Index","Cluster_Gini_Index","mean_vertex_size","Percentage_max_cluster_size","Percentage_max_vertex_size","Cluster_Reyni_Normalised", "Vertex_Reyni_Normalised")
analysis_matrices1 = q

##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------

## File Number 2: Class Switching Summary : 

file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_per_sequence_network_parameters_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check2 <- FALSE
} else {
	check2 <- TRUE
}
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
analysis_names = c( "Percentage_unique_BCRs_per_isotype", "Percentage_unique_BCRs_per_isotype_group")
analysis_matrices = list(m_unique_reads_per_isotype_group, m_unique_reads_per_isotype_single)
names(analysis_matrices) <- analysis_names

for(i in 1:length(analysis_matrices)){
	names <- names(analysis_matrices[i])
	colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
}
analysis_matrices2 = analysis_matrices

##---------------------------------------------------------------------------------------------------------------------
## File Number 3: Class Switching Summary : 

file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_CDR3_lengths_overall_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check3 <- FALSE
} else {
	check3 <- TRUE
}
p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
p=p[setdiff(c(1:length(p[,1])), grep("P", as.character(p[,"isotype"]))),]
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
analysis_names = c("Mean_CDR3_lengths")
names(analysis_matrices) <- analysis_names
for(i in 1:length(analysis_matrices)){
	names <- names(analysis_matrices[i])
	colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
} 
analysis_matrices3 = analysis_matrices

##---------------------------------------------------------------------------------------------------------------------
## File Number 4: ALL SHM : 

file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_SHM_Unmutated_sequences_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check4 <- FALSE
} else {
	check4 <- TRUE
}
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
analysis_names = c("Mean_SHM_per_BCR","Percentage_unmutated")

names(analysis_matrices) <- analysis_names
for(i in 1:length(analysis_matrices)){
	names <- names(analysis_matrices[i])
	colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
} 
analysis_matrices4 = analysis_matrices



##---------------------------------------------------------------------------------------------------------------------
## File Number 6: ALL SHM : 
file = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_SHM_Mutation_summmary_selection_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check6 <- FALSE
} else {
	check6 <- TRUE
}

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
names(analysis_matrices) <- analysis_names
for(i in 1:length(analysis_matrices)){
	names <- names(analysis_matrices[i])
	colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
} 

analysis_matrices6 = analysis_matrices
names(analysis_matrices6) <- c("mean_CDR_FWR_ratio", "Mean_mutations_per_BCR" ,"FWR3_mm", "mean_CDR_mm per BCR","mean_FWR_mm_per_BCR"  )

##---------------------------------------------------------------------------------------------------------------------
## File Number 7: ALL SHM : 

file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Secondary_rearrangements_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check7 <- FALSE
} else {
	check7 <- TRUE
}

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
analysis_names = c("V gene replacement frequency")

names(analysis_matrices) <- analysis_names
for(i in 1:length(analysis_matrices)){
	names <- names(analysis_matrices[i])
	colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
} 

analysis_matrices7 = analysis_matrices
names(analysis_matrices7) <- analysis_names

##---------------------------------------------------------------------------------------------------------------------
## File Number 8: ALL SHM : 

file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Cluster_expansion_isotype_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check8 <- FALSE
} else {
	check8 <- TRUE
}

p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
p=p[which(as.numeric(p[,"d5"]) !=-1),]
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

##---------------------------------------------------------------------------------------------------------------------
## File Number 9: ALL SHM : 
file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_IGHV4_34_quantification_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check9 <- FALSE
} else {
	check9 <- TRUE
}

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



##---------------------------------------------------------------------------------------------------------------------
## File Number 10: ALL SHM : 

file  = paste0(outputdir, "ORIENTATED_SEQUENCES/ISOTYPER/All_Secondary_rearrangements_clone_sizes_", iso_type, ".txt")
subsample_identifier <- grep("SAMPLED", file, value=TRUE)
if(length(subsample_identifier)==0){
	check10 <- FALSE
} else {
	check10 <- TRUE
}


p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]

id = as.character(p[,"X.sample"])
ids = sort(unique(id))
d5_norm =as.numeric(p[,"d5_norm"])
d5_secondary =as.numeric(p[,"d5_secondary"])
mean_clone_size_norm =as.numeric(p[,"mean_clone_size_norm"])
mean_clone_size_secondary =as.numeric(p[,"X_mean_clone_size_secondary"])
w1 = which(as.numeric(p[,"n_secondary"])>5)
groups = c(list(mean_clone_size_norm[w1]), list(mean_clone_size_secondary[w1] ))
#boxplot(groups)

groups = c(list(d5_norm[w1]), list(d5_secondary[w1] ))
#boxplot(groups)

headers = c("d5.norm","d5.secondary","mean.clone.size.norm", "mean.clone.size.secondary")
m_all_count = matrix(data = -1, nrow = length(ids_all),ncol = length(headers), dimnames=c(list(ids_all), list(headers)))
x=  cbind(d5_norm,d5_secondary,mean_clone_size_norm,mean_clone_size_secondary)
m_all_count[id,] = x

analysis_matrices = list(m_all_count)
analysis_names = c("V gene replacement clonal expansion")
names(analysis_matrices) <- analysis_names
for(i in 1:length(analysis_matrices)){
	names <- names(analysis_matrices[i])
	colnames(analysis_matrices[[i]]) <- paste0(names, "__", colnames(analysis_matrices[[i]]))
} 
names(analysis_matrices) <- analysis_names
analysis_matrices10 = analysis_matrices


##---------------------------------------------------------------------------------------------------------------------
## COMPOSING THE OVERALL MATRIX 

print_info = c(analysis_matrices1, analysis_matrices2,  analysis_matrices3, analysis_matrices4, analysis_matrices6, analysis_matrices7, analysis_matrices8, analysis_matrices9, analysis_matrices10)

overall_matrix = NULL
for(i in c(1:length(print_info))){
	if(length(overall_matrix)==0){
		overall_matrix = print_info[[i]]
	}else{overall_matrix = cbind(overall_matrix ,print_info[[i]][ids_all,])}
}
colnames(overall_matrix) = gsub(" ","_", colnames(overall_matrix), fixed= T)

overall_matrix <- data.frame(overall_matrix)
overall_matrix$sample <- rownames(overall_matrix)
overall_matrix$sample <- gsub("BCR_", "", overall_matrix$sample)
overall_matrix$sample <- gsub("TCRA_", "", overall_matrix$sample)
overall_matrix$sample <- gsub("TCRB_", "", overall_matrix$sample)
overall_matrix$sample <- gsub("TCRG_", "", overall_matrix$sample)
overall_matrix$sample <- gsub("TCRD_", "", overall_matrix$sample)

# Merge with readdepths 
overall_matrix <- merge(overall_matrix, read_depths_all, by.x="sample", by.y="SampleIDforDepths")

# save the output for this set of isotyper 
out_file_table =paste0(outputdir, "Summary/All_raw_values_SEPSIS_BCR1_", subsampled_depth_all, "_", iso_type, ".txt")
write.table(overall_matrix, file = out_file_table, append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")

##---------------------------------------------------------------------------------------------
##  Calculating the relationship between read depth and the metric 

new <- overall_matrix

## ***** Default of corelation test will be to omit na values *****
## Calculate correlation and missingness 
# note that missingness is really percentage of samples which are not na! (so high is good!)
values <- c("Metric", "percentage_present", "pval", "correlation")
for(i in 2:(length(colnames(new))-1)){
	variable <- new[,i]
	id <- colnames(new)[i]
	depths <- new$ReadDepth
	data <- data.frame(cbind(variable, depths))
	data[data=="-1"] <- "NA"
	data$variable <- as.numeric(data$variable)
	data$depths <- as.numeric(data$depths)
	missingness <- data$variable[!is.na(data$variable)]
	if(length(missingness) >= (0.25*length(ids_all))){
			missingness <- length(missingness) / length(ids_all) * 100
			rval <- cor.test(as.numeric(data[,1]),as.numeric(data[,2]), method="spearman")
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

# Remove spurious column names :D 
values$isotype[!values$isotype %in% counts_used$X.isotype] <- NA
values$percentage_present <- as.numeric(as.character(values$percentage_present))

## Assess columns with high missingness
bad_columns <- values$Metric[values$percentage_present<50]
bad_columns2 <- values$Metric[values$percentage_present<25]

## Save Significance to a txt file 
write.table(values, paste0(outputdir, "Summary/isotyper_metrics_", subsampled_depth_all,  "_", iso_type, ".txt"), sep="\t")

#Extract just those that are significant 
#Save file. 
significant_samples <- values[values$sig=="significant",]
write.table(significant_samples, paste0(outputdir, "Summary/isotyper_metrics_significant_", subsampled_depth_all, "_", iso_type, ".txt"), sep="\t")

##---------------------------------------------------------------------------------------------
## *** Make correllogram of all the different measures **** 

#Findcolumns with all missing values
data_1 <- data.frame(overall_matrix[, 1:(length(colnames(overall_matrix))-1)])
data_1[data_1=="-1"]<-NA

allmisscols <- apply(data_1,2, function(x)all(is.na(x)));
colswithallmiss <-names(allmisscols[allmisscols>0]);
colswithallmiss <- c(colswithallmiss, "sample")
#colswithallmiss

# Remove columns with all missing 
all_miss <- colnames(data_1)[!colnames(data_1) %in% colswithallmiss]
data_1 <- data_1[, c(all_miss)]

bad_columns <- unique(c(colswithallmiss, bad_columns))
new_cols <- colnames(data_1)[!colnames(data_1) %in% bad_columns]
data_2 <- data_1[, c(new_cols)]

bad_columns2 <- unique(c(colswithallmiss, bad_columns2))
new_cols2 <- colnames(data_1)[!colnames(data_1) %in% bad_columns2]
data_3 <- data_1[, c(new_cols2)]

data_1 <- as.matrix(data_1)
data_2 <- as.matrix(data_2)
data_3 <- as.matrix(data_3)
## Calculate the correlation between matrices 
## because we have nas we have to use pairwise complete obs 
## then the correlation or covariance between each pair of variables is computed using all complete pairs of observations on those variables. 
## This can result in covariance or correlation matrices which are not positive semi-definite, as well as NA entries if there are no complete pairs for that pair of variables.

correlation_matrix <- tryCatch(Hmisc::rcorr(data_1),error=function(e) e, warning=function(w) w)
correlation_matrix2 <- tryCatch(Hmisc::rcorr(data_2),error=function(e) e, warning=function(w) w)
correlation_matrix3 <- tryCatch(Hmisc::rcorr(data_3),error=function(e) e, warning=function(w) w)
  


## Added in 'trys' to handle errors 
pdf(paste0(outputdir, "Plots/Correlation_between_measures_", subsampled_depth_all, "_", iso_type,".pdf"), width=60, height=60)

if(!("warning" %in% class(correlation_matrix)) && !("error" %in% class(correlation_matrix))){
	try(corrplot(correlation_matrix[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix[[3]], sig.level = 0.05, insig = "blank", tl.col="black", title="No % Present threshold, p<0.05"))
	try(corrplot(correlation_matrix[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title="No % Present threshold: p<0.01"))
	try(corrplot(correlation_matrix[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix[[3]], sig.level = 0.001, insig = "blank", tl.col="black", title="No % Present threshold: p<0.001"))
} 

if(!("warning" %in% class(correlation_matrix2)) && !("error" %in% class(correlation_matrix2))){
	try(corrplot(correlation_matrix2[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix2[[3]], sig.level = 0.05, insig = "blank", tl.col="black", title="50% % Present threshold: p<0.05"))
	try(corrplot(correlation_matrix2[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix2[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title="50% % Present threshold: p<0.01"))
	try(corrplot(correlation_matrix2[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix2[[3]], sig.level = 0.001, insig = "blank", tl.col="black", title="50% % Present threshold: p<0.001"))
}

if(!("warning" %in% class(correlation_matrix3)) && !("error" %in% class(correlation_matrix3))){
	try(corrplot(correlation_matrix3[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix3[[3]], sig.level = 0.05, insig = "blank", tl.col="black", title="25% % Present threshold: p<0.05"))
	try(corrplot(correlation_matrix3[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix3[[3]], sig.level = 0.01, insig = "blank", tl.col="black", title="25% % Present threshold: p<0.01"))
	try(corrplot(correlation_matrix3[[1]], is.corr = FALSE, method = "color",type = "lower", p.mat = correlation_matrix3[[3]], sig.level = 0.001, insig = "blank", tl.col="black", title="25% % Present threshold: p<0.001"))
}
dev.off()
}