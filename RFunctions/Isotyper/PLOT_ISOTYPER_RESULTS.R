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
suppressMessages(library(ggrastr))
suppressMessages(library(purrr))
suppressMessages(library(ggpubr))

#analysis_matrices <- 
#outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRD/' 
#files <- list.files(paste0(outputdir, "Summary"), full.names=TRUE)
#files <- grep("isotyper_metrics_filtered_FINAL_METRICS_", files, value=TRUE)
#analysis_matrices <- files[1]
#file_use <-analysis_matrices 
#info <- paste0(outputdir, "Summary/isotyper_metrics_summary_stats200_ALL.txt")
#iso_type <- "ALL" 
#depth_file <- paste0(outputdir, "Summary/Read_Depths_ALL.txt")
#plot_isotyper_sepsis(file_use, outputdir, info_file, productivity, depth_file)

plot_isotyper_sepsis <- function(file_use=file_use, outputdir, info_file=info_file, iso_type, depth_file){
	
	analysis_matrices <- file_use
	info <- info_file
	analysis_matrices <- read.delim(analysis_matrices, sep="\t")
	info <- read.delim(info, sep="\t")
	
	depth_file <- read.delim(depth_file, sep="\t")
	depth_file$SampleIDforDepths <- gsub("_unproductive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("_productive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("BCR_", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("TCR_", "", depth_file$SampleIDforDepths)
	info$Metric <- gsub("__TCR", "", info$Metric)
	info$Metric <- gsub("__TCRA", "", info$Metric)
	info$Metric <- gsub("__TCRB", "", info$Metric)
	info$Metric <- gsub("__TCRG", "", info$Metric)
	info$Metric <- gsub("__TCRD", "", info$Metric)
	info$Metric <- gsub("__TRAC", "", info$Metric)
	info$Metric <- gsub("__TRBC", "", info$Metric)
	info$Metric <- gsub("__TRGC", "", info$Metric)
	info$Metric <- gsub("__TRDC", "", info$Metric)
	
	meta_data <- read.delim('/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Meta_data_for_Cohort1and2.txt', sep='\t', header=TRUE)
	
for(i in 1:(length(colnames(analysis_matrices)))){
			metric <- colnames(analysis_matrices[i])
			name<- metric
			data <- data.frame(analysis_matrices[, i])
			colnames(data) <- metric
			cols_to_plot <- colnames(data)
			data$sample <- row.names(analysis_matrices)
			data$sample <- gsub("BCR_", "", data$sample)
			data$sample <- gsub("TCRA_", "", data$sample)
			data$sample <- gsub("TCRB_", "", data$sample)
			data$sample <- gsub("TCRG_", "", data$sample)
			data$sample <- gsub("TCRD_", "", data$sample)
			data$sample <- gsub("_unproductive", "", data$sample)
			data$sample <- gsub("_productive", "", data$sample)
			
			data$timepoint <- NA
			data$timepoint[grep("_1", data$sample)] <- "Day1"
			data$timepoint[grep("_3", data$sample)] <- "Day3"
			data$timepoint[grep("_5", data$sample)] <- "Day5"
			data$timepoint[grep("T1", data$sample)] <- "Day1"
			data$timepoint[grep("T3", data$sample)] <- "Day3"
			data$timepoint[grep("T5", data$sample)] <- "Day5"
			data$timepoint[grep("_POSITIVE", data$sample)] <- "Technical_Replicate"
			
			data$Disease <- NA 
			data$Disease[data$sample %like% "HV" | data$sample %like% "POSITIVE"] <- "HEALTH"
			data$Disease[!data$sample %like% "HV" & !data$sample %like% "POSITIVE"] <- "SEPSIS"
			
			
			## Getting Read Depths 
			data <- merge(data, depth_file, by.x="sample", by.y="SampleIDforDepths")	
			data <- merge(data, meta_data, by.x="sample", by.y="SampleID", all.x=TRUE)
			
			## Editing the Variables for plotting 
			data$SRS <- as.character(data$SRS)
			data$SRS[data$Disease=="HEALTH" ] <- "HEALTH"	
			data$outcome <- NA 
			data$outcome[data$Disease=="HEALTH"] <- "HEALTH"
			data$outcome[data$Days_death_from_ICU=="Alive"] <- "ALIVE"
			data$outcome[data$Days_death_from_ICU!="Alive"] <- "DEAD"
			data$X7DAY_MORTALITY[data$Disease=="HEALTH"] <- "HEALTH"
			data$X7DAY_MORTALITY[data$X7DAY_MORTALITY==0] <- "7 Day Mortality"
			data$X7DAY_MORTALITY[data$X7DAY_MORTALITY==1] <- "Post 7 Day Mortality"
			data$X7DAY_MORTALITY[data$X7DAY_MORTALITY==2] <- "Alive"
			
			u <- data$sample[data$Disease=="HEALTH"]
			u <- str_split_fixed(u, "_", 3)
			u <- data.frame(u)
			u$sample <- paste0(u[,1], "_", u[,2])
			data$Barcode[data$Disease=="HEALTH"] <- u$sample
			## Setting Up data Structure
			data$SOFA_total <- as.numeric(data$SOFA_total)
			data$Age <- as.numeric(as.character(data$Age))
			data$Sex[is.na(data$Sex)] <- "Unknown"
			
			## Getting Read Depths 
			metric <- gsub("__",  ": ", metric)
			metric <- gsub("_",  " ", metric)
			metric <- gsub("\\.",  " ", metric)
			
			#Extract the subsampling depth!
			counts_used <- paste0(outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS")
			all_files <- list.files(counts_used, full.name=TRUE)
			all_files <- grep("depth_per_isotype", all_files, value=TRUE)
			counts_used <- read.delim(all_files[1], sep="\t", header=TRUE)
			counts_used <- counts_used[counts_used$type=="UNIQ",]
			subsampled_depth_all <- counts_used$min[counts_used$X.isotype=="all"]
			
			## Setting Up data Structure
			pdf(paste0(outputdir, "Plots/ISOTYPER/", iso_type, "/Plotting_SEPSIS_", name, "_", subsampled_depth_all, "_", iso_type, "_sepsis.pdf"), width=15, height=10)

			## Beginging the Plotting  
			for(s in cols_to_plot){
				s1 <- gsub("BCR_", "", s)
				s1 <- gsub("TCR_", "", s1)
				s1 <- gsub("TCRA_", "", s1)
				s1 <- gsub("TCRB_", "", s1)
				s1 <- gsub("TCRG_", "", s1)
				s1 <- gsub("TCRD_", "", s1)
							
				s1 <- gsub("TRAC_", "", s1)
				s1 <- gsub("TRBC_", "", s1)
				s1 <- gsub("TRGC_", "", s1)
				s1 <- gsub("TRDC_", "", s1)
				
				s1 <- gsub("ALL_READS_", "", s1)
				s1 <- gsub("UNPRODUCTIVE_READS_", "", s1)
				s1 <- gsub("PRODUCTIVE_READS_", "", s1)
			
				metric1 <- metric
				metric2 <- metric1
				
				##Start plotting
				d <- ggplot(data, aes_string(x="timepoint", y=s))   + geom_boxplot( aes(fill=Disease))  + theme_bw() +xlab("Timepoint") +ylab(metric1) +ggtitle(metric2)
				e2 <- ggplot(data, aes_string(x="timepoint", y=s))   + geom_boxplot(aes(fill=Disease))  + theme_bw() +facet_wrap(~X7DAY_MORTALITY, scales = "free_x") +xlab("Timepoint") +ylab(metric1) +ggtitle(metric2)
				f <- ggplot(data, aes_string(x="timepoint", y=s))  +geom_point(aes(color=SRS, shape=Sex), size=3) + theme_bw() + geom_line(aes(group = Barcode), color='black') +xlab("Timepoint") +ylab(metric1)+ggtitle(metric2)
				g2 <- ggplot(data, aes_string(x="timepoint", y=s))  +geom_point(aes(color=SRS, shape=Sex), size=3) + theme_bw() + geom_line(aes(group = Barcode), color='black') +facet_wrap(~X7DAY_MORTALITY, scales = "free_x") +xlab("Timepoint") +ylab(metric1) +ggtitle(metric2)

				tryCatch({
				grid.arrange(d, e2, f, g2, ncol=2)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	

				d <- ggplot(data, aes_string(x="ReadDepth", y=s, colour="ReadDepth")) + geom_point() + theme_bw() +xlab("ReadDepth") +ylab(metric1) + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) + geom_smooth(method='lm') +ggtitle(metric2)+ stat_cor(label.x = 0)
				e <- ggplot(data, aes_string(x="ReadDepth", y=s, colour="ReadDepth")) + geom_point(aes(color=ReadDepth))  + theme_bw() +facet_wrap(~X7DAY_MORTALITY) +xlab("ReadDepth") +ylab(metric1) + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) + geom_smooth(method='lm')+ggtitle(metric2) + stat_cor(label.x = 0)
				f <- ggplot(data, aes_string(x="timepoint", y=s))+ geom_point(aes(color=ReadDepth, shape=Sex), size=3) + theme_bw() + geom_line(aes(group = Barcode), color='black') +xlab("Timepoint") +ylab(metric1)+ scale_color_gradient(high = "yellow", low = "darkblue")+ggtitle(metric2) + labs(colour="Read Depth")+ stat_cor(label.x = 0)
				g <- ggplot(data, aes_string(x="timepoint", y=s)) + geom_point(aes(color=ReadDepth, shape=Sex), size=3) + theme_bw() + geom_line(aes(group = Barcode), color='black') +facet_wrap(~X7DAY_MORTALITY) +xlab("Timepoint") +ylab(metric1)+scale_color_gradient(high = "yellow", low = "darkblue")+ggtitle(metric2)+ labs(colour="Read Depth")+ stat_cor(label.x = 0)
				tryCatch({
				grid.arrange(d, e, f, g, ncol=2)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	

				e <- ggplot(data, aes_string(x="Age", y=s))  +geom_point(aes(color=X7DAY_MORTALITY, shape=Sex), size=3) + theme_bw() +facet_wrap(~timepoint) +xlab("Age") +ylab(metric1) + geom_line(aes(group = Barcode), color='black')+ggtitle(metric2)
				f <- ggplot(data, aes_string(x="SOFA_total", y=s))  +geom_point(aes(color=X7DAY_MORTALITY, shape=Sex), size=3) + theme_bw() + geom_line(aes(group = Barcode), color='black') +facet_wrap(~timepoint) +xlab("Sofa Organ Score") +ylab(metric1)+ggtitle(metric2)
				g <- ggplot(data, aes_string(x="Age", y=s))  +geom_point(aes(color=timepoint, shape=Sex), size=3) + theme_bw()  +xlab("Age") +ylab(metric1) + geom_line(aes(group = Barcode), color='black')+ggtitle(metric2)+facet_wrap(~timepoint)
				h <- ggplot(data, aes_string(x="SOFA_total", y=s))  +geom_point(aes(color=timepoint, shape=Sex), size=3) + theme_bw() + geom_line(aes(group = Barcode), color='black') +xlab("Sofa Organ Score") +ylab(metric1)+ggtitle(metric2)+facet_wrap(~timepoint)
				tryCatch({
				grid.arrange(e, f, g, h, ncol=2)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
				
				e <- ggplot(data, aes_string(x="SRS", y=s))  + geom_boxplot(aes(fill=timepoint))+ theme_bw() +facet_wrap(~X7DAY_MORTALITY, scales = "free_x") +xlab("SRS Asignment") +ylab(metric1) +ggtitle(metric2)
				e1 <- ggplot(data, aes_string(x="SRS", y=s)) + geom_boxplot(aes(fill=timepoint)) + theme_bw()  +xlab("SRS Asignment") +ylab(metric1)+ggtitle(metric2)
				tryCatch({
				grid.arrange(e, e1, ncol=2)	
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
				
				e <- ggplot(data, aes_string(x="EBV_positive_ddpcr", y=s))  + geom_boxplot(aes(fill=timepoint))+ theme_bw() +facet_wrap(~X7DAY_MORTALITY, scales = "free_x") +xlab("EBV Status: ddPCR") +ylab(metric1) +ggtitle(metric2)
				e1 <- ggplot(data, aes_string(x="EBV_positive_ddpcr", y=s)) + geom_boxplot(aes(fill=timepoint)) + theme_bw()  +xlab("EBV Status: ddPCR") +ylab(metric1)+ggtitle(metric2)
				tryCatch({
				grid.arrange(e, e1, ncol=2)	
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
				
				e <- ggplot(data, aes_string(x="EBV_positive_metagenomics", y=s))  + geom_boxplot(aes(fill=timepoint))+ theme_bw() +facet_wrap(~X7DAY_MORTALITY, scales = "free_x") +xlab("EBV Status: Metagenomics") +ylab(metric1) +ggtitle(metric2)
				e1 <- ggplot(data, aes_string(x="EBV_positive_metagenomics", y=s)) + geom_boxplot(aes(fill=timepoint)) + theme_bw()  +xlab("EBV Status: Metagenomics") +ylab(metric1)+ggtitle(metric2)
				tryCatch({
				grid.arrange(e, e1, ncol=2)	
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
				
				e <- ggplot(data, aes_string(x="Shock_sepsis2", y=s))  + geom_boxplot(aes(fill=timepoint))+ theme_bw() +facet_wrap(~X7DAY_MORTALITY, scales = "free_x") +xlab("Shock Sepsis2") +ylab(metric1) +ggtitle(metric2)
				e1 <- ggplot(data, aes_string(x="Shock_sepsis2", y=s)) + geom_boxplot(aes(fill=timepoint)) + theme_bw()  +xlab("Shock Sepsis2") +ylab(metric1)+ggtitle(metric2)
				tryCatch({
				grid.arrange(e, e1, ncol=2)	
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
				
				e <- ggplot(data, aes_string(x="Sex", y=s))  + geom_boxplot(aes(fill=timepoint))+ theme_bw() +facet_wrap(~X7DAY_MORTALITY, scales = "free_x") +xlab("Sex") +ylab(metric1) +ggtitle(metric2)
				e1 <- ggplot(data, aes_string(x="Sex", y=s)) + geom_boxplot(aes(fill=timepoint)) + theme_bw()  +xlab("Sex") +ylab(metric1)+ggtitle(metric2)
				tryCatch({
				grid.arrange(e, e1, ncol=2)	
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
				}
			dev.off()
		}
		
}

plot_isotyper <- function(file_use=file_use, outputdir, info, iso_type, depth_file){
	
	analysis_matrices <- file_use
	analysis_matrices <- read.delim(analysis_matrices, sep="\t")
	info <- read.delim(info, sep="\t")
	depth_file <- read.delim(depth_file, sep="\t")
	depth_file$SampleIDforDepths <- gsub("_unproductive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("_productive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("BCR_", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("TCR_", "", depth_file$SampleIDforDepths)
	info$Metric <- gsub("__TCR", "", info$Metric)
	info$Metric <- gsub("__TCRA", "", info$Metric)
	info$Metric <- gsub("__TCRB", "", info$Metric)
	info$Metric <- gsub("__TCRG", "", info$Metric)
	info$Metric <- gsub("__TCRD", "", info$Metric)
	info$Metric <- gsub("__TRAC", "", info$Metric)
	info$Metric <- gsub("__TRBC", "", info$Metric)
	info$Metric <- gsub("__TRGC", "", info$Metric)
	info$Metric <- gsub("__TRDC", "", info$Metric)
	
	dir.create(paste0(outputdir, "Plots/ISOTYPER/", iso_type))
	
	for(i in 1:(length(colnames(analysis_matrices)))){
			metric <- colnames(analysis_matrices[i])
			name<- metric
			data <- data.frame(analysis_matrices[, i])
			colnames(data) <- metric
			cols_to_plot <- colnames(data)
			data$sample <- rownames(analysis_matrices)
			data$sample <- gsub("BCR_", "", data$sample)
			data$sample <- gsub("TCRA_", "", data$sample)
			data$sample <- gsub("TCRB_", "", data$sample)
			data$sample <- gsub("TCRG_", "", data$sample)
			data$sample <- gsub("TCRD_", "", data$sample)
			data$sample <- gsub("_unproductive", "", data$sample)
			data$sample <- gsub("_productive", "", data$sample)
			controls <- grep("HV", data$sample, value=TRUE)
			metric <- gsub("__",  ": ", metric)
			metric <- gsub("_",  " ", metric)
			metric <- gsub("\\.",  " ", metric)
			
			## Getting Read Depths 
			data <- merge(data, depth_file, by.x="sample", by.y="SampleIDforDepths")
			
			## Getting subsample depths which were used for isotyper script 
			counts_used <- paste0(outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS")
			all_files <- list.files(counts_used, full.name=TRUE)
			all_files <- grep("depth_per_isotype", all_files, value=TRUE)
			counts_used <- read.delim(all_files[1], sep="\t", header=TRUE)
			counts_used <- counts_used[counts_used$type=="UNIQ",]
			subsampled_depth_all <- counts_used$min[counts_used$X.isotype=="all"]			
				
			## Setting Up data Structure
			pdf(paste0(outputdir, "Plots/ISOTYPER/", iso_type, "/Plotting_", name, "_", subsampled_depth_all, "_", iso_type, ".pdf"), width=20, height=13)

            ## Beginging the Plotting  
			for(s in cols_to_plot){
				#print(s)
				s1 <- gsub("BCR_", "", s)
				s1 <- gsub("TCR_", "", s1)
				s1 <- gsub("TCRA_", "", s1)
				s1 <- gsub("TCRB_", "", s1)
				s1 <- gsub("TCRG_", "", s1)
				s1 <- gsub("TCRD_", "", s1)			
				s1 <- gsub("TRAC_", "", s1)
				s1 <- gsub("TRBC_", "", s1)
				s1 <- gsub("TRGC_", "", s1)
				s1 <- gsub("TRDC_", "", s1)
				s1 <- gsub("ALL_READS_", "", s1)
				s1 <- gsub("UNPRODUCTIVE_READS_", "", s1)
				s1 <- gsub("PRODUCTIVE_READS_", "", s1)
				
				# Get subsampled depth
				subsampled_val <- info$subsampled_depth[info$Metric==s1]			
				if(length(subsampled_val)==0){ 
					subsampled_depth <- "Not Performed"
				} else {
					subsampled_depth <- subsampled_val
				} 
									
				metric1 <- metric
				metric2 <- paste0("Subsampled Depth: ", subsampled_depth)
				
				d <- ggplot(data, aes_string(x=s))  + rasterise(geom_density(), dpi=300)  + theme_bw()  +ggtitle(metric2) + geom_density(color="black", fill="lightblue") + geom_vline(aes(xintercept=mean(data[,s][data[, s] >-1 & !is.na(data[, s])])),color="blue", linetype="dashed", size=1) +xlab(metric1)
				
				e <- ggplot(data, aes_string(x="ReadDepth", y=s, colour="ReadDepth")) + rasterise(geom_point(aes(color=ReadDepth)), dpi=300)  + theme_bw() +xlab("ReadDepth") +ylab(metric1) + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) + geom_smooth(method='lm')+ggtitle(metric2)+ labs(colour="Read Depth") + stat_cor(label.x = 0)


				grid.arrange(d, e, ncol=1)
				
				f <- ggplot(data, aes(x=sample, y=as.numeric(data[,s]))) + rasterise(geom_boxplot(), dpi = 300) +theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.7, hjust=1)) + stat_summary(fun=mean, geom="point", shape=23, size=2, fill="blue", position = position_dodge(width = .75)) + ylab(metric1) +xlab("Sample") + ggtitle(productivity)
				
				plot(f)
			
			dev.off()
			}
	}
}

plot_isotyper_comparison <- function(analysis_matrices_list, outputdir, info){
	
	dfiles <- list.files(paste0(outputdir, "Summary"), full.names=TRUE)
	dfiles <- grep("Read_Depths_", dfiles, value=TRUE)
	dprod <- grep("_PRODUCTIVE", dfiles, value=TRUE)
	dprod <- read.delim(dprod, sep="\t")
	duprod <- grep("_UNPRODUCTIVE", dfiles, value=TRUE)
	duprod <- read.delim(duprod, sep="\t")
	dall <- grep("_ALL", dfiles, value=TRUE)
	dall <- read.delim(dall, sep="\t")
	dprod$SampleIDforDepths <- gsub("_productive", "", dprod$SampleIDforDepths)
	duprod$SampleIDforDepths <- gsub("_unproductive", "", duprod$SampleIDforDepths)
	
	analysis_matrices1 <- read.delim(analysis_matrices_list[1], sep="\t")
	analysis_matrices1$productivity <- analysis_matrices_list[1]
	analysis_matrices1$sample <- row.names(analysis_matrices1)
	analysis_matrices2 <- read.delim(analysis_matrices_list[2], sep="\t")
	analysis_matrices2$productivity <- analysis_matrices_list[2]
	analysis_matrices2$sample <- row.names(analysis_matrices2)
	analysis_matrices3 <- read.delim(analysis_matrices_list[3], sep="\t")
	analysis_matrices3$productivity <- analysis_matrices_list[3]
	analysis_matrices3$sample <- row.names(analysis_matrices3)
	
	colnames(analysis_matrices1) <- gsub("TCR_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("BCR_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("TCRA_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("TCRB_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("TCRG_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("TCRD_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("UNPRODUCTIVE_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("PRODUCTIVE_", "", colnames(analysis_matrices1))
	colnames(analysis_matrices1) <- gsub("ALL_", "", colnames(analysis_matrices1))
	
	colnames(analysis_matrices2) <- gsub("TCR_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("BCR_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("TCRA_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("TCRB_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("TCRG_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("TCRD_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("UNPRODUCTIVE_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("PRODUCTIVE_", "", colnames(analysis_matrices2))
	colnames(analysis_matrices2) <- gsub("ALL_", "", colnames(analysis_matrices2))
	
	colnames(analysis_matrices3) <- gsub("TCR_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("BCR_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("TCRA_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("TCRB_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("TCRG_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("TCRD_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("UNPRODUCTIVE_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("PRODUCTIVE_", "", colnames(analysis_matrices3))
	colnames(analysis_matrices3) <- gsub("ALL_", "", colnames(analysis_matrices3))
	
	big <- dplyr::bind_rows(analysis_matrices1, analysis_matrices2, analysis_matrices3)
	depth_file <- read.delim(depth_file, sep="\t")
	depth_file$SampleIDforDepths <- gsub("_unproductive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("_productive", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("BCR_", "", depth_file$SampleIDforDepths)
	depth_file$SampleIDforDepths <- gsub("TCR_", "", depth_file$SampleIDforDepths)
	
	
	duprod$SampleIDforDepths <- gsub("_unproductive", "", duprod$SampleIDforDepths)
	duprod$SampleIDforDepths <- gsub("_productive", "", duprod$SampleIDforDepths)
	duprod$SampleIDforDepths <- gsub("BCR_", "", duprod$SampleIDforDepths)
	duprod$SampleIDforDepths <- gsub("TCR_", "", duprod$SampleIDforDepths)
	
	dall$SampleIDforDepths <- gsub("_unproductive", "", dall$SampleIDforDepths)
	dall$SampleIDforDepths <- gsub("_productive", "", dall$SampleIDforDepths)
	dall$SampleIDforDepths <- gsub("BCR_", "", dall$SampleIDforDepths)
	dall$SampleIDforDepths <- gsub("TCR_", "", dall$SampleIDforDepths)
	
	dprod$SampleIDforDepths <- gsub("_unproductive", "", dprod$SampleIDforDepths)
	dprod$SampleIDforDepths <- gsub("_productive", "", dprod$SampleIDforDepths)
	dprod$SampleIDforDepths <- gsub("BCR_", "", dprod$SampleIDforDepths)
	dprod$SampleIDforDepths <- gsub("TCR_", "", dprod$SampleIDforDepths)
	
	for(i in 1:length(big$productivity)){ 
		if(grepl("_UNPRODUCTIVE", big$productivity[i])){
			big$productivity[i] <- "UNPRODUCTIVE"
		} 
		if(grepl("_PRODUCTIVE", big$productivity[i])){
			big$productivity[i] <- "PRODUCTIVE"
		} 
		if(grepl("ALL", big$productivity[i])){
			big$productivity[i] <- "ALL"
		}
	} 
	
	info <- read.delim(info, sep="\t")

	info$Metric <- gsub("__TCR", "", info$Metric)
	info$Metric <- gsub("__TCRA", "", info$Metric)
	info$Metric <- gsub("__TCRB", "", info$Metric)
	info$Metric <- gsub("__TCRG", "", info$Metric)
	info$Metric <- gsub("__TCRD", "", info$Metric)
	info$Metric <- gsub("__TRAC", "", info$Metric)
	info$Metric <- gsub("__TRBC", "", info$Metric)
	info$Metric <- gsub("__TRGC", "", info$Metric)
	info$Metric <- gsub("__TRDC", "", info$Metric)
	
	
	for(i in 1:(length(big)-2)){
			metric <- colnames(big[i])
			name<- metric
			data <- data.frame(big[, i])
			colnames(data) <- metric
			cols_to_plot <- colnames(data)
			data$productivity <- big$productivity
			data$sample <- big$sample
			data$sample <- gsub("BCR_", "", data$sample)
			data$sample <- gsub("TCRA_", "", data$sample)
			data$sample <- gsub("TCRB_", "", data$sample)
			data$sample <- gsub("TCRG_", "", data$sample)
			data$sample <- gsub("TCRD_", "", data$sample)
			data$sample <- gsub("_unproductive", "", data$sample)
			data$sample <- gsub("_productive", "", data$sample)
			controls <- grep("HV", data$sample, value=TRUE)
			metric <- gsub("__",  ": ", metric)
			metric <- gsub("_",  " ", metric)
			metric <- gsub("\\.",  " ", metric)
			data$depths <- NA
			for(i in 1:length(data$sample)){
				c <- data$sample[i]
				data$depths[data$productivity=="ALL" & data$sample==c] <- dall$ReadDepth[dall$SampleIDforDepths==c]
				data$depths[data$productivity=="UNPRODUCTIVE" & data$sample==c] <- duprod$ReadDepth[duprod$SampleIDforDepths==c]
				data$depths[data$productivity=="PRODUCTIVE" & data$sample==c] <- dprod$ReadDepth[dprod$SampleIDforDepths==c]
			}
		
			info$subsampled_depth <- as.numeric(info$subsampled_depth)
			subsampled_depth_all <- max(info$subsampled_depth[!is.na(info$subsampled_depth)])
			## Setting Up data Structure
			pdf(paste0(outputdir, "Plots/ISOTYPER/Plotting_", name, "_", subsampled_depth_all, "_COMPARISON.pdf"), width=15, height=10)
			## Beginging the Plotting  
			for(s in cols_to_plot){
				s1 <- gsub("BCR_", "", s)
				s1 <- gsub("TCR_", "", s1)
				s1 <- gsub("TCRA_", "", s1)
				s1 <- gsub("TCRB_", "", s1)
				s1 <- gsub("TCRG_", "", s1)
				s1 <- gsub("TCRD_", "", s1)
				s1 <- gsub("ALL_", "", s1)
						
				s1 <- gsub("TRAC_", "", s1)
				s1 <- gsub("TRBC_", "", s1)
				s1 <- gsub("TRGC_", "", s1)
				s1 <- gsub("TRDC_", "", s1)
				
				
				s1 <- gsub("UNPRODUCTIVE_", "", s1)
				s1 <- gsub("PRODUCTIVE_", "", s1)
				# Get subsampled depth
				subsampled_val <- info$subsampled_depth[info$Metric==s1]			
				if(length(subsampled_val)==0){ 
					subsampled_depth <- "Not Performed"
				} else {
					subsampled_depth <- subsampled_val
				} 
										
				metric1 <- metric
				metric2 <- paste0("Subsampled Depth: ", subsampled_depth)
				data <- data[data[, s] >-1 & !is.na(data[, s]),]
				data <- data.frame(data)
				d <- ggplot(data, aes_string(x=s, fill="productivity")) + rasterise(geom_density(alpha=.4), dpi=300)  + theme_bw()  +xlab(metric1) +ggtitle(metric2) 
				e <- ggplot(data, aes_string(x="depths", y=s, colour="depths")) + rasterise(geom_point(aes(color=depths)), dpi=300) + theme_bw() +xlab("ReadDepth") +ylab(metric1) + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))+ geom_smooth(method='lm')+ggtitle(metric2)+ labs(colour="Read Depth")+ facet_wrap(~productivity, scales="free_x")+ stat_cor(label.x = 0)
				tryCatch({
				plot(grid.arrange(d, e, ncol=1))
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
			dev.off()
			}
	}
}
			
			
			