# Author: Lauren Overend 
# Function to visualise IgM and IgD SHM in the sepsis cohort 
# Author: Lauren Overend 
# December 2021
#module load R/4.1.0-foss-2021a  
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggforce))
suppressMessages(library(Gviz))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(gtools))
suppressMessages(library(purrr))
suppressMessages(library(reshape2))
suppressMessages(library(Hmisc))
suppressMessages(library(corrplot))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(matrixStats))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrastr))
suppressMessages(library(ggpubr))
suppressMessages(library(plot3D))
library(diptest)	
`%notin%` <- Negate(`%in%`)	

## For troubleshooting 
path_to_outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/'
productivity="PRODUCTIVE"
eigenvectors='/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt'

imgt_mutation_statistics_sepsis <- function(path_to_outputdir = path_to_outputdir, cluster_nodes = 8, productivity=productivity, path_to_layout=path_to_layout){
	
	## Makes files considerabley smaller!
	path <- path_to_outputdir
	path <- paste0(path_to_outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT")
	files <- list.files(path, full.name=TRUE)
	files <- grep('_Summary.txt', files, value=TRUE)
	
	## Overall summary file 
	## We are only using this to extract the column names !!
	path2 <- paste0(path_to_outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW")
	files2 <- list.dirs(path = path2, full.names = TRUE, recursive = TRUE)
	path2 <- files2[2]
	files2 <- list.files(path2, full.name=TRUE)
	files2 <- grep('_Summary.txt', files2, value=TRUE)
	
	data_overall <- read.delim(files2, header=TRUE, sep='\t')
	data_overall$X <- NULL
	names_data <- colnames(data_overall) 
	data_overall <- NULL
	
	
	## Extracting the actual mutation counts 
	if(productivity=="ALL" | productivity=="all"){
			files_not <- grep('_productive_', files, value=TRUE)
			files_not2 <- grep('_unproductive', files, value=TRUE)
			files <- files[!files %in% c(files_not, files_not2)]
	} 
	
	if(!is.na(productivity)){	
		if(productivity=="productive" | productivity=="PRODUCTIVE" ){	
				files <- grep('_productive_', files, value=TRUE)
			} 
			
		if(productivity=="unproductive"| productivity=="UNPRODUCTIVE"){	
				files <- grep('_unproductive', files, value=TRUE)
			}
	}
	
	### We want to find sequences that are present in both isotypes and compare mean SHM
	## This is the same way rachael calculated the total mutations across the CDR3
	grouped_results <- c()
	all_seq <- c()
	for(i in 1:length(files)){
		print(i)
		### Lets set it up for just one at a time
		datax <- read.delim(files[i], header=FALSE)
		colnames(datax) <- names_data
		datax <- data.frame(datax)
		
		## For V gene
		muts_v <- datax$V.REGION.identity.nt
		muts_v <- gsub(" nt", "", muts_v)
		muts_v <- data.frame(str_split_fixed(muts_v, "\\/", 2))
		muts_v$mut <- as.numeric(muts_v[,2])-as.numeric(muts_v[,1])
		muts_v$percent <- datax$V.REGION.identity..
		datax$mutv <- muts_v$mut
		datax$mutvl <- muts_v$X2
		
		#For J gene 
		muts_j <- datax$J.REGION.identity.nt
		muts_j <- gsub(" nt", "", muts_j)
		muts_j <- data.frame(str_split_fixed(muts_j, "\\/", 2))
		muts_j$mut <- as.numeric(muts_j[,2])-as.numeric(muts_j[,1])
		muts_j$percent <- datax$J.REGION.identity..
		datax$mutj <- muts_j$mut
		datax$mutjl <- muts_j$X2
		
		## Calculate total mutation
		datax$NumberMutations <- datax$mutv+datax$mutj
		datax$PercentMut <- 	(datax$NumberMutations/(as.numeric(datax$mutvl) +as.numeric(datax$mutjl)))*100

		## Lets get the sample id 
		d <- str_split(files[i], 'IMGT_SPLIT/') 
		d <- sapply(d, "[[", 2)  
		d <- gsub("_8_V-REGION-nt-mutation-statistics.txt", "", d)
		d <- gsub("_1_Summary.txt", "", d)
		d <- gsub("IMGT_", "", d)
		d <- gsub("BCR_", "", d)
		d <- gsub("_productive", "", d)
		datax$Sample <- d
				
		## okay V gene number of mutations is the same!! Lets just use imgt ones 
		new <- datax[, c("Sequence.number", "Sequence.ID", "NumberMutations", "PercentMut", "mutv", "mutvl","mutj", "mutjl", "CDR3.IMGT.length", "V.GENE.and.allele", "J.GENE.and.allele", "Analysed.sequence.length")]
		new$Sample <- d
		new$seqid <- str_split_fixed(new$Sequence.ID, "__", 2)[,1]
		
		##-----------------------------------
		## Now lets extract V gene information
		new$V.GENE.and.allele <- gsub("Homsap ", "", new$V.GENE.and.allele)
		genes <- strsplit(new$V.GENE.and.allele, " F", fixed=TRUE)
		w <- map(genes, 1) 
		w<- matrix(unlist(w), ncol=1, byrow=TRUE)	
		# Get just V gene Family
		genes <- strsplit(w, "*", fixed=TRUE)   
		w <- map(genes, 1)
		new$V.GENE<- matrix(unlist(w), ncol=1, byrow=TRUE)

		## Now lets extract V gene information
		new$J.GENE.and.allele <- gsub("Homsap ", "", new$J.GENE.and.allele)
		genes <- strsplit(new$J.GENE.and.allele, " F", fixed=TRUE)
		w <- map(genes, 1) 
		w<- matrix(unlist(w), ncol=1, byrow=TRUE)	
		# Get just V gene Family
		genes <- strsplit(w, "*", fixed=TRUE)   
		w <- map(genes, 1)
		new$J.GENE<- matrix(unlist(w), ncol=1, byrow=TRUE)
		##---------------------------------------------
		
		## now go to look up table I made previously 
		file_lookup <- paste0(path_to_outputdir, "Summary/SummarySequences/", d, "_PRODUCTIVE_allsequences.txt")
		look_up <- read.delim(file_lookup, header=TRUE)
		new2 <- merge(new,look_up, by="seqid")
	
		## We want to assign if they are IgM or IgD or Both 
		new2$IGM_class[new2$IGHM > 0 & new2$IGHD ==0] <- "IGM"
		new2$IGM_class[new2$IGHM==0 & new2$IGHD >0] <- "IGD"
		new2$IGM_class[new2$IGHM > 0 & new2$IGHD >0] <- "IGM_IGD"
		
		### We are only interested in the IgM sequences!!!
		new3 <- new2[!is.na(new2$IGM_class),]
		new3$NumberMutations <- as.numeric(new3$NumberMutations)
		
		### Lets get the mean number of mutations 
		means <- new3 %>%group_by(IGM_class) %>%dplyr::summarize(MeanMutNo = mean(NumberMutations, na.rm=TRUE), MeanMutPer = mean(PercentMut, na.rm=TRUE) )
		means$sample <- d
		means <- data.frame(means)
		
		## Lets bind the number of mutations
		grouped_results <- rbind(grouped_results, means)
		## Lets bind all the sequences 
		all_seq <- rbind(all_seq, new3)
		}
	
	##------------------------------------------------------------------------------------------------
	#################################################################################################
	all_seq$vdjlength <- as.numeric(all_seq$mutvl) + as.numeric(all_seq$mutjl)
	all_seq$NormalisedCDR3 <- as.numeric(all_seq$CDR3.IMGT.length)/as.numeric(all_seq$Analysed.sequence.length)
	
	## Create a back up copy for later use and save 		
	all_seq2 <- all_seq
	write.table(all_seq2, paste0(path_to_outputdir, "/SHM_Distribution_IGM_IGD.txt"))
	
	#### Now lets do some analysis 
	## Lets do a biolin plot showing mutations 
	all_seq$Sample.y = NULL
	all_seq$DISEASE <- "SEPSIS"
	all_seq$DISEASE[all_seq$Sample.x %like% "HV"] <- "HEALTH"

	## Good to continue 
	all_seq$DAY <- NA
	all_seq$DAY[grep("_1", all_seq$Sample.x)] <- "1"
	all_seq$DAY[grep("_3", all_seq$Sample.x)] <- "3"
	all_seq$DAY[grep("_5", all_seq$Sample.x)] <- "5"

	## need to take out bad ids
	bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
	all_seq <- all_seq[!all_seq$Sample.x %in% c(bad_ids),]
	all_seq <- all_seq[!all_seq$Sample.x %like% "JR1795_1003_POSITIVE",]
	
	## Subset to only those samples in the final data!
	#### I think easiest way would be to take a downsample of x reads per individual and then plot 
	eigenvectors <- read.delim(eigenvectors)
	samples_keep <- eigenvectors$sample
	all_seq <- all_seq[all_seq$Sample.x %in% c(samples_keep),]
	write.table(all_seq, paste0(path_to_outputdir, "/SHM_Distribution_IGM_IGD_nobadids.txt"))
	
		
	pdf(paste0(path_to_outputdir, "/SHM_PercentvsCount.pdf"), width=10, height=10)
	p1 <- ggplot(grouped_results, aes(x=MeanMutNo, y=MeanMutPer))+geom_point()+theme_classic()+xlab("Mean Number of Mutations") + ylab("Mean % SHM")
	p2 <- ggplot(all_seq, aes(x=NumberMutations, y=PercentMut))+ geom_hex(bins=100) +scale_fill_viridis_c()+theme_classic()+xlab("Number of Mutations") + ylab("% SHM")
	p3 <-  ggplot(all_seq, aes(x=vdjlength))+ geom_histogram(bins=100) +scale_fill_viridis_c()+theme_classic()+xlab("VJ Length") + ylab("Count")
	p4 <- ggplot(all_seq, aes(x=CDR3.IMGT.length, y=NormalisedCDR3))+ geom_hex(bins=50) +scale_fill_viridis_c()+theme_classic()+xlab("IMGT CDR3 LENGTH") + ylab("CDR3/READ LENGTH")
	plot_grid(p1, p2, p3, p4, ncol = 2, align="hv", axis="tblr")
	dev.off()


	##------------------------------------------------------------------------------------------------
	################################################################################################
	
	
	#### CALCULATE THE BIMODALITY COEFFICIENT AND DIP TEST 
	Bimod <- all_seq %>%dplyr::group_by(Sample.x, IGM_class ) %>%dplyr::summarize(BCC = bimodality_coefficient(NumberMutations), DTC = dip(NumberMutations), BCP = bimodality_coefficient(PercentMut), DTP = dip(PercentMut))
	Bimod$disease <- "SEPSIS"
	Bimod$disease[Bimod$Sample.x %like% "HV"] <- "HEALTH"
	Bimod$DAY <- NA
	Bimod$DAY[grep("_1", Bimod$Sample.x)] <- "1"
	Bimod$DAY[grep("_3", Bimod$Sample.x)] <- "3"
	Bimod$DAY[grep("_5", Bimod$Sample.x)] <- "5"
	Bimod <- data.frame(Bimod)

	################################################################
	## Now we want to calculate the counts per sample per class 
	counts_IGM <- table(all_seq$Sample.x[all_seq$IGM_class=="IGM"])
	counts_IGD <- table(all_seq$Sample.x[all_seq$IGM_class=="IGD"])
	
	## What is the read depth we will actually use for subsampling 
	read_depth <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/ORIENTATED_SEQUENCES/ANNOTATIONS/Sampling_depth_per_isotype_LEO_SEPSIS_BCR_ALL_post.txt'
	read_depth <- read.delim(read_depth)
	IGM_depth <- read_depth$min[read_depth$X.isotype=="IGHM" & read_depth$type=="UNIQ"]
	IGD_depth <- read_depth$min[read_depth$X.isotype=="IGHD" & read_depth$type=="UNIQ"]

	##------------------------------------------------------------------------------------------------
	################################################################################################

	### Generate a random subsample of each sample so that the plot of the distribution is REPRESENTITIVE (ige. equally contributed to by each sample) 
	## E.g. we arent going to get a skew if one sample has millions of high shm reads 
	## We are going to take a random sample from each repertoire of the same depth so we dont get a density plot biased by certain high read depth samples 
	random_subset <- c()
	
	## Lets do this 10 times and we can visulaise 
	## Function to subset 
	subset_rep <- function(datax=datax, IGM_depth=IGM_depth, IGD_depth=IGD_depth){
		random_subset <- c()
		for(x in 1:length(unique(datax$Sample.x))){
				print(x)
				rand_igd <- NULL
				rand_igm <- NULL
				data_use <- datax[datax$Sample.x==unique(datax$Sample.x)[x],]
				Igm_use <- data_use[data_use$IGM_class=="IGM",]
				Igd_use <- data_use[data_use$IGM_class=="IGD",]
				
				if(dim(Igm_use)[1] >= IGM_depth){
					rand_igm <- Igm_use[sample(nrow(Igm_use), size=IGM_depth), ]
				}
				if(dim(Igd_use)[1] >= IGD_depth){
					rand_igd <- Igd_use[sample(nrow(Igd_use), size=IGD_depth), ]
				}
				random_subset <- rbind(random_subset, rand_igd, rand_igm)
		}
		return(random_subset)
	}

	## Actually RUN it 
	all_out  <- c()
	for(g in 1:10){
		out <- subset_rep(all_seq, IGM_depth, IGD_depth)
		out$rep <- g
		all_out <- rbind(all_out, out)
	}


	## now lets plot distribution with each sample equally represented 
	thresh <- 17 
	thresh2 <- 7
	
	pdf(paste0(path_to_outputdir, "/SHM_Distribution", thresh, ".pdf"), width=20, height=14)
	#ggplot(all_out, aes(x=DAY, y=V.REGION.Nb.of.mutations, fill=disease)) +geom_violin() +facet_wrap(~IGM_class)+theme_classic()
	#ggplot(all_out, aes(x=V.REGION.Nb.of.mutations, fill=Sample.x)) + geom_density(alpha=0.6)+facet_grid(rows=vars(disease), cols=vars(IGM_class))+theme_classic()+ geom_density(alpha=.2) +xlab("V Region Number of Mutations")+guides(fill="none")
	#p1 <- ggplot(all_out[all_out$IGM_class=="IGM",], aes(x=V.REGION.Nb.of.mutations, fill=Sample.x)) + geom_density(alpha=0.6)+facet_grid(rows=vars(disease), cols=vars(DAY))+theme_classic()+ geom_density(alpha=.2) +xlab("V Region Number of Mutations")+guides(fill="none") +ggtitle("IGM")
	#p2 <- ggplot(all_out[all_out$IGM_class=="IGD",], aes(x=V.REGION.Nb.of.mutations, fill=Sample.x)) + geom_density(alpha=0.6)+facet_grid(rows=vars(disease), cols=vars(DAY))+theme_classic()+ geom_density(alpha=.2) +xlab("V Region Number of Mutations")+guides(fill="none")+ggtitle("IGD")
	#plot_grid(p1, p2, ncol = 1, align="hv", axis="tblr")
	p1 <- ggplot(all_out[all_out$IGM_class=="IGM",], aes(x=NumberMutations, fill=DISEASE)) + geom_density(alpha=0.6)+facet_grid(cols=vars(DAY), rows=vars(rep))+theme_classic()+ geom_density(alpha=.2) +xlab("No. of Mutations")+guides(fill="none") +ggtitle("IGHM")+geom_vline(xintercept=thresh)
	p2 <- ggplot(all_out[all_out$IGM_class=="IGD",], aes(x=NumberMutations, fill=DISEASE)) + geom_density(alpha=0.6)+facet_grid( cols=vars(DAY), rows=vars(rep))+theme_classic()+ geom_density(alpha=.2) +xlab("No. of Mutations")+guides(fill="none")+ggtitle("IGHD")+geom_vline(xintercept=thresh)
	plot_grid(p1, p2, ncol = 1, align="hv", axis="tblr")
	p12 <- ggplot(all_out[all_out$IGM_class=="IGM",], aes(x=PercentMut, fill=DISEASE)) + geom_density(alpha=0.6)+facet_grid(cols=vars(DAY), rows=vars(rep))+theme_classic()+ geom_density(alpha=.2) +xlab("% SHM")+guides(fill="none") +ggtitle("IGHM")+geom_vline(xintercept=thresh2)
	p22 <- ggplot(all_out[all_out$IGM_class=="IGD",], aes(x=PercentMut, fill=DISEASE)) + geom_density(alpha=0.6)+facet_grid( cols=vars(DAY), rows=vars(rep))+theme_classic()+ geom_density(alpha=.2) +xlab("% SHM")+guides(fill="none")+ggtitle("IGHD")+geom_vline(xintercept=thresh2)
	plot_grid(p1, p12, p2, p22, ncol = 2, align="hv", axis="tblr")
	#ggline(Bimod, x = "DAY", y = "DT", add = "mean_se", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Percentage", facet.by="IGM_class")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	#ggline(Bimod, x = "DAY", y = "BC", add = "mean_se", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Percentage", facet.by="IGM_class")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	dev.off()

	##------------------------------------------------------------------------------------------------
	################################################################################################

	# Now we can go back to the original file and group by < 13 or more than 13 mutation s
	all_seq$Class[all_seq$NumberMutations <= thresh] <- "LOW SHM\n(Naive?)"
	all_seq$Class[all_seq$NumberMutations> thresh] <- "HIGH SHM\n(Plasma/Memory)"
	all_seq$Class <- factor(all_seq$Class, levels=c("LOW SHM\n(Naive?)", "HIGH SHM\n(Plasma/Memory)"))
	
	###---------------------------------------------------------------------------------
	##'# Mean SHM per cell type 
	means_celltype <- all_seq %>%dplyr::group_by(Sample.x, IGM_class, Class) %>%dplyr::summarize(mean_shm = mean(NumberMutations), mean_shm_per = mean(PercentMut))
	means_celltype$disease <- "SEPSIS"
	means_celltype$disease[means_celltype$Sample.x %like% "HV"] <- "HEALTH"
	## Good to continue 
	means_celltype$DAY <- NA
	means_celltype$DAY[grep("_1", means_celltype$Sample.x)] <- "1"
	means_celltype$DAY[grep("_3", means_celltype$Sample.x)] <- "3"
	means_celltype$DAY[grep("_5", means_celltype$Sample.x)] <- "5"
	means_celltype <- data.frame(means_celltype)
	
	###---------------------------------------------------------------------------------
	##'# Mean CDR3 per cell type 
	cdr3_celltype <- all_seq %>%dplyr::group_by(Sample.x, IGM_class, Class) %>%dplyr::summarize(mean_cdr3 = mean(CDR3.IMGT.length))
	cdr3_celltype$disease <- "SEPSIS"
	cdr3_celltype$disease[cdr3_celltype$Sample.x %like% "HV"] <- "HEALTH"
	## Good to continue 
	cdr3_celltype$DAY <- NA
	cdr3_celltype$DAY[grep("_1", cdr3_celltype$Sample.x)] <- "1"
	cdr3_celltype$DAY[grep("_3", cdr3_celltype$Sample.x)] <- "3"
	cdr3_celltype$DAY[grep("_5", cdr3_celltype$Sample.x)] <- "5"
	cdr3_celltype <- data.frame(cdr3_celltype)
	
	### Mean CDR3 per cell type 
	cdr3_celltype2 <- all_seq %>%dplyr::group_by(Sample.x, IGM_class, Class) %>%dplyr::summarize(mean_cdr3 = mean(NormalisedCDR3))
	cdr3_celltype2$disease <- "SEPSIS"
	cdr3_celltype2$disease[cdr3_celltype2$Sample.x %like% "HV"] <- "HEALTH"
	## Good to continue 
	cdr3_celltype2$DAY <- NA
	cdr3_celltype2$DAY[grep("_1", cdr3_celltype2$Sample.x)] <- "1"
	cdr3_celltype2$DAY[grep("_3", cdr3_celltype2$Sample.x)] <- "3"
	cdr3_celltype2$DAY[grep("_5", cdr3_celltype2$Sample.x)] <- "5"
	cdr3_celltype2 <- data.frame(cdr3_celltype2)
	
	
	###---------------------------------------------------------------------------------
	## Calculate prop in each group 
	prop_celltype <- all_seq %>%group_by(Sample.x, IGM_class, Class) %>% tally() %>% mutate(freq = prop.table(n))
	prop_celltype$disease <- "SEPSIS"
	prop_celltype$disease[prop_celltype$Sample.x %like% "HV"] <- "HEALTH"
	## Good to continue 
	prop_celltype$DAY <- NA
	prop_celltype$DAY[grep("_1", prop_celltype$Sample.x)] <- "1"
	prop_celltype$DAY[grep("_3", prop_celltype$Sample.x)] <- "3"
	prop_celltype$DAY[grep("_5", prop_celltype$Sample.x)] <- "5"
	prop_celltype <- data.frame(prop_celltype)
	
	###---------------------------------------------------------------------------------
	## Calculate prop in each group using Duplicated Counts 
	### We also want to look and see if non-unique can shed a light!
	all_seq$DupCount <-   all_seq$IGHM + all_seq$IGHD
	all_seq_dup <- as.data.frame(lapply(all_seq, rep, all_seq$DupCount))
	### Prop
	prop_celltype_dup <- all_seq_dup %>%group_by(Sample.x, IGM_class, Class) %>% tally() %>% mutate(freq = prop.table(n))
	prop_celltype_dup$disease <- "SEPSIS"
	prop_celltype_dup$disease[prop_celltype_dup$Sample.x %like% "HV"] <- "HEALTH"
	## Good to continue 
	prop_celltype_dup$DAY <- NA
	prop_celltype_dup$DAY[grep("_1", prop_celltype_dup$Sample.x)] <- "1"
	prop_celltype_dup$DAY[grep("_3", prop_celltype_dup$Sample.x)] <- "3"
	prop_celltype_dup$DAY[grep("_5", prop_celltype_dup$Sample.x)] <- "5"
	prop_celltype_dup <- data.frame(prop_celltype_dup)
	
	###---------------------------------------------------------------------------------
	## Calculate V and J gene Usage 
	jgene <- all_seq %>%group_by(Sample.x, IGM_class, Class, J.GENE) %>% tally() %>% mutate(freq = prop.table(n))
	jgene$disease <- "SEPSIS"
	jgene$disease[jgene$Sample.x %like% "HV"] <- "HEALTH"
	## Good to continue 
	jgene$DAY <- NA
	jgene$DAY[grep("_1", jgene$Sample.x)] <- "1"
	jgene$DAY[grep("_3", jgene$Sample.x)] <- "3"
	jgene$DAY[grep("_5", jgene$Sample.x)] <- "5"
	jgene <- data.frame(jgene)
	
	## V gene family 
	all_seq$VGENEFAM <- str_split_fixed(all_seq$V.GENE, "\\-", 2)[,1]
	vgene <- all_seq %>%group_by(Sample.x, IGM_class, Class, VGENEFAM) %>% tally() %>% mutate(freq = prop.table(n))
	vgene$disease <- "SEPSIS"
	vgene$disease[vgene$Sample.x %like% "HV"] <- "HEALTH"
	## Good to continue 
	vgene$DAY <- NA
	vgene$DAY[grep("_1", vgene$Sample.x)] <- "1"
	vgene$DAY[grep("_3", vgene$Sample.x)] <- "3"
	vgene$DAY[grep("_5", vgene$Sample.x)] <- "5"
	vgene <- data.frame(vgene)
	
	### Specific V gene of interest 
	## V gene family 
	all_vgenes <- unique(all_seq$V.GENE)[unique(all_seq$V.GENE) %like% "IGHV3-30",]
	all_seq$newgene <- all_seq$V.GENE
	all_seq$newgene[all_seq$newgene %in% c(all_vgenes)] <- "IGHV3-30"
	vgeneallele <- all_seq %>%group_by(Sample.x, IGM_class, Class, newgene) %>% tally() %>% mutate(freq = prop.table(n))
	vgeneallele$disease <- "SEPSIS"
	vgeneallele$disease[vgeneallele$Sample.x %like% "HV"] <- "HEALTH"
	## Good to continue 
	vgeneallele$DAY <- NA
	vgeneallele$DAY[grep("_1", vgeneallele$Sample.x)] <- "1"
	vgeneallele$DAY[grep("_3", vgeneallele$Sample.x)] <- "3"
	vgeneallele$DAY[grep("_5", vgeneallele$Sample.x)] <- "5"
	vgeneallele <- data.frame(vgeneallele)
	vgeneallele$percent <- vgeneallele$freq *100
	## Subset for alleles of interst 
	vgeneallele <- vgeneallele[vgeneallele$newgene == "IGHV3-30",]
	
	
	### May want to normalise CDR3 length relative to total sequence length 
	##########################
 	### Clearly see a shift in the number of reads falling into the high SHM population and also the level of SHM in them is much higher in sepsis 
	pdf(paste0(path_to_outputdir, "/SHM_Distribution_CellType_", thresh, ".pdf"), width=15, height=15)
	p1 <- ggline(means_celltype, x = "DAY", y = "mean_shm", add = "mean_ci", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Mean no. SHM", facet.by=c("IGM_class", "Class"), main="SHM")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +labs(colour="Disease")
	p2 <- ggline(prop_celltype, x = "DAY", y = "freq", add = "mean_ci", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Proportion of Reads (Unique)", facet.by=c("IGM_class", "Class"), main="Repertoire Usage")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +labs(colour="Disease")
	p3 <- ggline(prop_celltype_dup, x = "DAY", y = "freq", add = "mean_ci", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Proportion of Reads (Total)", facet.by=c("IGM_class", "Class"), main="Repertoire Usage")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +labs(colour="Disease")
	p4  <- ggline(cdr3_celltype2, x = "DAY", y = "mean_cdr3", add = "mean_ci", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="CDR3 Length/VDJ Length", facet.by=c("IGM_class", "Class"), main="Normalised CDR3 Length")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +labs(colour="Disease")
	p5 <- ggline(vgeneallele, x = "DAY", y = "percent", add = "mean_ci", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Percent of Reads", facet.by=c("IGM_class", "Class"), main="IGHV3-30 V Gene Usage")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +labs(colour="Disease")
	
	p11 <- ggline(means_celltype, x = "DAY", y = "mean_shm", add = "mean_ci", color = "IGM_class", palette = c("#00AFBB", "#E7B800"), ylab="Mean no. SHM", facet.by=c("Class", "disease"), main="SHM")+ stat_compare_means(aes(group = IGM_class), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +labs(colour="Isotype")
	p22 <- ggline(prop_celltype, x = "DAY", y = "freq", add = "mean_ci", color = "IGM_class", palette = c("#00AFBB", "#E7B800"), ylab="Proportion of Reads (Unique)", facet.by=c("Class", "disease"), main="Repertoire Usage")+ stat_compare_means(aes(group = IGM_class), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+labs(colour="Isotype")
	p33 <- ggline(prop_celltype_dup, x = "DAY", y = "freq", add = "mean_ci", color = "IGM_class", palette = c("#00AFBB", "#E7B800"), ylab="Proportion of Reads (Total)", facet.by=c("Class", "disease"), main="Repertoire Usage")+ stat_compare_means(aes(group = IGM_class), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+labs(colour="Isotype")
	p44  <- ggline(cdr3_celltype2, x = "DAY", y = "mean_cdr3", add = "mean_ci", color = "IGM_class", palette = c("#00AFBB", "#E7B800"), ylab="CDR3 Length/VDJ Length", facet.by=c("Class", "disease"), main="Normalised CDR3 Length")+ stat_compare_means(aes(group = IGM_class), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+labs(colour="Isotype")
	p55 <- ggline(vgeneallele, x = "DAY", y = "percent", add = "mean_ci", color = "IGM_class", palette = c("#00AFBB", "#E7B800"), ylab="Percent of Reads", facet.by=c("Class", "disease"), main="IGHV3-30 V Gene Usage")+ stat_compare_means(aes(group = IGM_class), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+labs(colour="Isotype")
	plot_grid(p1, p11, p2, p22, p3, p33, p4, p44, p5, p55, ncol = 4, align="hv", axis="tblr", labels="AUTO")
	
	#################################
	pdf(paste0(path_to_outputdir, "/SHM_Distribution_CellType2_", thresh, ".pdf"), width=15, height=15)
	p1 <- ggline(means_celltype, x = "DAY", y = "mean_shm_per", add = "mean_ci", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Mean % SHM", facet.by=c("IGM_class", "Class"))+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p2 <- ggline(prop_celltype, x = "DAY", y = "freq", add = "mean_ci", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Proportion of Reads (Unique)", facet.by=c("IGM_class", "Class"))+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p3 <- ggline(prop_celltype_dup, x = "DAY", y = "freq", add = "mean_ci", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Proportion of Reads (Total)", facet.by=c("IGM_class", "Class"))+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p11 <- ggline(means_celltype, x = "DAY", y = "mean_shm_per", add = "mean_ci", color = "IGM_class", palette = c("#00AFBB", "#E7B800"), ylab="Mean % SHM", facet.by=c("Class", "disease"))+ stat_compare_means(aes(group = IGM_class), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p22 <- ggline(prop_celltype, x = "DAY", y = "freq", add = "mean_ci", color = "IGM_class", palette = c("#00AFBB", "#E7B800"), ylab="Proportion of Reads (Unique)", facet.by=c("Class", "disease"))+ stat_compare_means(aes(group = IGM_class), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p222 <- ggline(prop_celltype_dup, x = "DAY", y = "freq", add = "mean_ci", color = "IGM_class", palette = c("#00AFBB", "#E7B800"), ylab="Proportion of Reads (Total)", facet.by=c("Class", "disease"))+ stat_compare_means(aes(group = IGM_class), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	plot_grid(p1, p11, p2, p22,  p3, p22, ncol = 3, align="hv", axis="tblr", labels="AUTO")
	dev.off()
	
	pdf(paste0(path_to_outputdir, "/SHM_Distribution_CellType_J", thresh, ".pdf"), width=12, height=12)
	p1 <- ggline(jgene[jgene$IGM_class=="IGM",], x = "DAY", y = "freq", add = "mean_ci", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Percent Repertoire", facet.by=c("J.GENE", "Class"), main="IGM")+ stat_compare_means(aes(group = disease), label = "p.signif", label.y=0.3) #+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p2 <- ggline(jgene[jgene$IGM_class=="IGD",], x = "DAY", y = "freq", add = "mean_ci", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Percent Repertoire", facet.by=c("J.GENE", "Class"), main="IGD")+ stat_compare_means(aes(group = disease), label = "p.signif", label.y=0.3) #+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))	
	p11 <- ggline(jgene[jgene$IGM_class=="IGM",], x = "DAY", y = "freq", add = "mean_ci", color = "Class", palette = c("#00AFBB", "#E7B800"), ylab="Percent Repertoire", facet.by=c("J.GENE", "disease"), main="IGM")+ stat_compare_means(aes(group = Class), label = "p.signif", label.y=0.3) #+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p22 <- ggline(jgene[jgene$IGM_class=="IGD",], x = "DAY", y = "freq", add = "mean_ci", color = "Class", palette = c("#00AFBB", "#E7B800"), ylab="Percent Repertoire", facet.by=c("J.GENE", "disease"), main="IGD")+ stat_compare_means(aes(group = Class), label = "p.signif", label.y=0.3) #+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))	
	plot_grid(p1, p2, p11, p22, ncol = 2, align="hv", axis="tblr")
	p1 <- ggline(vgene[vgene$IGM_class=="IGM",], x = "DAY", y = "freq", add = "mean_ci", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Percent Repertoire", facet.by=c("VGENEFAM", "Class"), main="IGM")+ stat_compare_means(aes(group = disease), label = "p.signif", label.y=0.3) #+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p2 <- ggline(vgene[vgene$IGM_class=="IGD",], x = "DAY", y = "freq", add = "mean_ci", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Percent Repertoire", facet.by=c("VGENEFAM", "Class"), main="IGD")+ stat_compare_means(aes(group = disease), label = "p.signif", label.y=0.3) #+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))	
	p11 <- ggline(vgene[vgene$IGM_class=="IGM",], x = "DAY", y = "freq", add = "mean_ci", color = "Class", palette = c("#00AFBB", "#E7B800"), ylab="Percent Repertoire", facet.by=c("VGENEFAM", "disease"), main="IGM")+ stat_compare_means(aes(group = Class), label = "p.signif", label.y=0.3) #+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p22 <- ggline(vgene[vgene$IGM_class=="IGD",], x = "DAY", y = "freq", add = "mean_ci", color = "Class", palette = c("#00AFBB", "#E7B800"), ylab="Percent Repertoire", facet.by=c("VGENEFAM", "disease"), main="IGD")+ stat_compare_means(aes(group = Class), label = "p.signif", label.y=0.3) #+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))	
	plot_grid(p1, p2, p11, p22, ncol = 2, align="hv", axis="tblr")
	p1 <- ggline(vgene[vgene$IGM_class=="IGM",], x = "DAY", y = "freq", add = "mean_se", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Percent Repertoire", facet.by=c("VGENEFAM", "Class"), main="IGM")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p2 <- ggline(vgene[vgene$IGM_class=="IGD",], x = "DAY", y = "freq", add = "mean_se", color = "disease", palette = c("#00AFBB", "#E7B800"), ylab="Percent Repertoire", facet.by=c("VGENEFAM", "Class"), main="IGD")+ stat_compare_means(aes(group = disease), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))	
	plot_grid(p1, p2, ncol = 2, align="hv", axis="tblr")
	p3 <- ggline(shm_usagex1, x = "DAY", y = "mean_mutations", add = "mean_se", color = "ISOTYPE", ylab="Mean Mutations per VDJ", facet.by="disease")+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p3.1 <- ggline(shm_usagex1[shm_usagex1$ISOTYPE=="IGHM" |shm_usagex1$ISOTYPE=="IGHD",], x = "DAY", y = "mean_mutations", add = "mean_se", color = "ISOTYPE", ylab="Mean Mutations per VDJ", facet.by="disease")+ stat_compare_means(aes(group = ISOTYPE), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p4 <- ggline(shm_usagex, x = "DAY", y = "mean_mutations", add = "mean_se", color = "ISOTYPE", ylab="Non-silent:silent mutations ratio", facet.by="disease")+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	p4.1 <- ggline(shm_usagex[shm_usagex$ISOTYPE=="IGHM" |shm_usagex$ISOTYPE=="IGHD",], x = "DAY", y = "mean_mutations", add = "mean_se", color = "ISOTYPE", ylab="Non-silent:silent mutations ratio", facet.by="disease")+ stat_compare_means(aes(group = ISOTYPE), label = "p.signif") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
	plot_grid(p3, p3.1, p4, p4.1, ncol = 2, align="hv", axis="tblr")
	dev.off()

	
	
  }
  
######
### DONE
####
	