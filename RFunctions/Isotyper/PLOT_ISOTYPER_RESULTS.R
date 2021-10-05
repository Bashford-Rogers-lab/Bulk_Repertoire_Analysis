plot_isotyper_sepsis <- function(analysis_matrices, outputdir, info, iso_type){
	
	analysis_matrices <- read.delim(analysis_matrices, sep="\t")
	info <- read.delim(info, sep="\t")
	
	meta_data <- read.delim('/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/Cohort1/Meta_data_for_cohort_1.txt', sep='\t', header=TRUE)
	LDL <- read.delim('/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LDL_sample1.txt', sep='\t', header=TRUE)
	
for(i in 2:(length(analysis_matrices)-1)){
			metric <- colnames(analysis_matrices[i])
			name<- metric
			data <- data.frame(analysis_matrices[, i])
			colnames(data) <- metric
			cols_to_plot <- colnames(data)
			data$sample <- analysis_matrices$sample
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
			data$timepoint[grep("_P", data$sample)] <- "Control"
			data$timepoint[grep("T0", data$sample)] <- "Control"
			
			## Getting Read Depths 
			depths <- analysis_matrices$ReadDepth
			data <- cbind(data, depths)
			
			
			controls <- grep("HV", data$sample, value=TRUE)
			samplesx <- c(controls, as.character(meta_data$SampleID))
			data <- data[data$sample %in% samplesx,]
			data <- merge(data, meta_data, by.x="sample", by.y="SampleID", all.x=TRUE)
			data$SRS <- as.character(data$SRS)
			data$SRS[data$timepoint == "Control"] <- "Not Applicable"	
			data$outcome <- NA 
			data$outcome[is.na(data$Days_death_from_ICU)] <- "CONTROL"
			data$outcome[data$Days_death_from_ICU=="Alive"] <- "ALIVE"
			data$outcome[data$Days_death_from_ICU!="Alive"] <- "DEAD"
			data$outcome[data$timepoint=="Control"] <- "CONTROL"
			
			## Setting Up data Structure
			data$SOFA_total <- as.numeric(data$SOFA_total)
			data$Age <- as.numeric(as.character(data$Age))
			data$Sex[is.na(data$Sex)] <- "Unknown"
			data <- merge(data, LDL, by.x="SampleID_alternative", by.y="SampleID", all.x=TRUE)
			
			## Getting Read Depths 
			metric <- gsub("__",  ": ", metric)
			metric <- gsub("_",  " ", metric)
			metric <- gsub("\\.",  " ", metric)
			

			
			info$subsampled_depth <- as.numeric(info$subsampled_depth)
			subsampled_depth_all <- max(info$subsampled_depth[!is.na(info$subsampled_depth)])
			
			## Setting Up data Structure
			pdf(paste0(outputdir, "Plots/ISOTYPER/Plotting_", name, "_", subsampled_depth_all, "_", iso_type, "_sepsis.pdf"), width=12, height=10)

			
            ## Beginging the Plotting  
			for(s in cols_to_plot){
				#print(s)
				
				# Get subsampled depth
				subsampled_val <- info$subsampled_depth[info$Metric==s]
							
				if(is.na(subsampled_val)){ 
					subsampled_depth <- "Not Performed"
				} else {
					subsampled_depth <- subsampled_val
				} 
									
				metric1 <- metric
				metric2 <- paste0("Subsampled Depth: ", subsampled_depth)
				
				d <- ggplot(data[data[, s] >-1,], aes_string(x="timepoint", y=s))   + geom_boxplot( aes(fill=timepoint))  + theme_bw() +xlab("Timepoint") +ylab(metric1) +ggtitle(metric2)
				e <- ggplot(data[data[, s] >-1,], aes_string(x="timepoint", y=s))   + geom_boxplot(aes(fill=timepoint))  + theme_bw() +facet_wrap(~outcome, scales = "free_x") +xlab("Timepoint") +ylab(metric1) +ggtitle(metric2)
				f <- ggplot(data[data[, s] >-1,], aes_string(x="timepoint", y=s))  +geom_point(aes(color=SRS, shape=Sex), size=3) + theme_classic() + geom_line(aes(group = Barcode), color='grey') +xlab("Timepoint") +ylab(metric1)+ggtitle(metric2)
				g <- ggplot(data[data[, s] >-1,], aes_string(x="timepoint", y=s))  +geom_point(aes(color=SRS, shape=Sex), size=3) + theme_classic() + geom_line(aes(group = Barcode), color='grey') +facet_wrap(~outcome, scales = "free_x") +xlab("Timepoint") +ylab(metric1) +ggtitle(metric2)
				tryCatch({
				plot(grid.arrange(d, e, f, g, ncol=2))
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	

				d <- ggplot(data[data[, s] >-1,], aes_string(x="depths", y=s, colour="depths")) + geom_point() + theme_bw() +xlab("ReadDepth") +ylab(metric1) + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) + stat_summary(fun.data= mean_cl_normal) + geom_smooth(method='lm') +ggtitle(metric2)
				e <- ggplot(data[data[, s] >-1,], aes_string(x="depths", y=s, colour="depths")) + geom_point(aes(color=depths))  + theme_bw() +facet_wrap(~outcome, scales = "free_x") +xlab("ReadDepth") +ylab(metric1) + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) +stat_summary(fun.data= mean_cl_normal) + geom_smooth(method='lm')+ggtitle(metric2)
				f <- ggplot(data[data[, s] >-1,], aes_string(x="timepoint", y=s))+ geom_point(aes(color=depths, shape=Sex), size=3) + theme_classic() + geom_line(aes(group = Barcode), color='grey') +xlab("Timepoint") +ylab(metric1)+ scale_color_gradient(high = "yellow", low = "darkblue")+ggtitle(metric2) + labs(colour="Read Depth")
				g <- ggplot(data[data[, s] >-1,], aes_string(x="timepoint", y=s)) + geom_point(aes(color=depths, shape=Sex), size=3) + theme_classic() + geom_line(aes(group = Barcode), color='grey') +facet_wrap(~outcome, scales = "free_x") +xlab("Timepoint") +ylab(metric1)+scale_color_gradient(high = "yellow", low = "darkblue")+ggtitle(metric2)+ labs(colour="Read Depth")
				tryCatch({
				plot(grid.arrange(d, e, f, g, ncol=2))
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	

				e <- ggplot(data[data[, s] >-1,], aes_string(x="Age", y=s))  +geom_point(aes(color=timepoint, shape=Sex), size=3) + theme_classic() +facet_wrap(~outcome, scales = "free_x") +xlab("Age") +ylab(metric1) + geom_line(aes(group = Barcode), color='grey')+ggtitle(metric2)
				f <- ggplot(data[data[, s] >-1,], aes_string(x="SOFA_total", y=s))  +geom_point(aes(color=timepoint, shape=Sex), size=3) + theme_classic() + geom_line(aes(group = Barcode), color='grey') +facet_wrap(~outcome, scales = "free_x") +xlab("Sofa Organ Score") +ylab(metric1)+ggtitle(metric2)
				g <- ggplot(data[data[, s] >-1,], aes_string(x="Age", y=s))  +geom_point(aes(color=timepoint, shape=Sex), size=3) + theme_classic()  +xlab("Age") +ylab(metric1) + geom_line(aes(group = Barcode), color='grey')+ggtitle(metric2)
				h <- ggplot(data[data[, s] >-1,], aes_string(x="SOFA_total", y=s))  +geom_point(aes(color=timepoint, shape=Sex), size=3) + theme_classic() + geom_line(aes(group = Barcode), color='grey') +xlab("Sofa Organ Score") +ylab(metric1)+ggtitle(metric2)
				tryCatch({
				plot(grid.arrange(e, f, g, h, ncol=2))
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
				e <- ggplot(data[data[, s] >-1,], aes_string(x="LV", y=s))  +geom_point(aes(color=timepoint, shape=Sex), size=3) + theme_classic() +facet_wrap(~outcome, scales = "free_x") +xlab("LV") +ylab(metric1) + geom_line(aes(group = Barcode), color='grey')+ggtitle(metric2)
				f <- ggplot(data[data[, s] >-1,], aes_string(x="MLV", y=s))  +geom_point(aes(color=timepoint, shape=Sex), size=3) + theme_classic() +facet_wrap(~outcome, scales = "free_x") +xlab("MLV") +ylab(metric1) + geom_line(aes(group = Barcode), color='grey')+ggtitle(metric2)
				g <- ggplot(data[data[, s] >-1,], aes_string(x="PLV", y=s))  +geom_point(aes(color=timepoint, shape=Sex), size=3) + theme_classic() +facet_wrap(~outcome, scales = "free_x") +xlab("PLV") +ylab(metric1) + geom_line(aes(group = Barcode), color='grey')+ggtitle(metric2)
				e1 <- ggplot(data[data[, s] >-1,], aes_string(x="LV", y=s))  +geom_point(aes(color=timepoint, shape=Sex), size=3) + theme_classic()  +xlab("LV") +ylab(metric1) + geom_line(aes(group = Barcode), color='grey')+ggtitle(metric2)
				f1 <- ggplot(data[data[, s] >-1,], aes_string(x="MLV", y=s))  +geom_point(aes(color=timepoint, shape=Sex), size=3) + theme_classic()  +xlab("MLV") +ylab(metric1) + geom_line(aes(group = Barcode), color='grey')+ggtitle(metric2)
				g1 <- ggplot(data[data[, s] >-1,], aes_string(x="PLV", y=s))  +geom_point(aes(color=timepoint, shape=Sex), size=3) + theme_classic() +xlab("PLV") +ylab(metric1) + geom_line(aes(group = Barcode), color='grey')+ggtitle(metric2)
				tryCatch({
				plot(grid.arrange(e, f, g, e1, f1, g1, ncol=3))	
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
				}
			dev.off()
		}
		
}

plot_isotyper <- function(analysis_matrices, outputdir, info, iso_type){
	
	analysis_matrices <- read.delim(analysis_matrices, sep="\t")
	info <- read.delim(info, sep="\t")

	for(i in 2:(length(analysis_matrices)-1)){
			metric <- colnames(analysis_matrices[i])
			name<- metric
			data <- data.frame(analysis_matrices[, i])
			colnames(data) <- metric
			cols_to_plot <- colnames(data)
			data$sample <- analysis_matrices$sample
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
			depths <- analysis_matrices$ReadDepth
			data <- cbind(data, depths)

			info$subsampled_depth <- as.numeric(info$subsampled_depth)
			subsampled_depth_all <- max(info$subsampled_depth[!is.na(info$subsampled_depth)])
			
			## Setting Up data Structure
			pdf(paste0(outputdir, "Plots/ISOTYPER/Plotting_", name, "_", subsampled_depth_all, "_", iso_type, ".pdf"), width=12, height=10)


            ## Beginging the Plotting  
			for(s in cols_to_plot){
				#print(s)
				
				# Get subsampled depth
				subsampled_val <- info$subsampled_depth[info$Metric==s]
							
				if(is.na(subsampled_val)){ 
					subsampled_depth <- "Not Performed"
				} else {
					subsampled_depth <- subsampled_val
				} 
									
				metric1 <- metric
				metric2 <- paste0("Subsampled Depth: ", subsampled_depth)
				
				d <- ggplot(data[data[, s] >-1,], aes_string(y=s))   + geom_boxplot()  + theme_bw()  +ylab(metric1) +ggtitle(metric2) +theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
				e <- ggplot(data[data[, s] >-1,], aes_string(x="depths", y=s, colour="depths")) + geom_point(aes(color=depths))  + theme_bw() +xlab("ReadDepth") +ylab(metric1) + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) +stat_summary(fun.data= mean_cl_normal) + geom_smooth(method='lm')+ggtitle(metric2)+ labs(colour="Read Depth")

				tryCatch({
				plot(grid.arrange(d, e, ncol=2))
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
			dev.off()
			}
	}
}

plot_isotyper_comparison <- function(analysis_matrices_list, outputdir, info){
	
	analysis_matrices1 <- read.delim(analysis_matrices_list[1], sep="\t")
	analysis_matrices1$productivity <- analysis_matrices_list[1]
	analysis_matrices2 <- read.delim(analysis_matrices_list[2], sep="\t")
	analysis_matrices2$productivity <- analysis_matrices_list[2]
	analysis_matrices3 <- read.delim(analysis_matrices_list[3], sep="\t")
	analysis_matrices3$productivity <- analysis_matrices_list[3]
	
	big <- dplyr::bind_rows(analysis_matrices1, analysis_matrices2, analysis_matrices3)
	
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

	for(i in 2:(length(big)-2)){
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
			
			## Getting Read Depths 
			depths <- big$ReadDepth
			data <- cbind(data, depths)

			info$subsampled_depth <- as.numeric(info$subsampled_depth)
			subsampled_depth_all <- max(info$subsampled_depth[!is.na(info$subsampled_depth)])
			
			## Setting Up data Structure
			pdf(paste0(outputdir, "Plots/ISOTYPER/Plotting_", name, "_", subsampled_depth_all, "_COMPARISON.pdf"), width=21, height=7)


            ## Beginging the Plotting  
			for(s in cols_to_plot){
				#print(s)
				
				# Get subsampled depth
				subsampled_val <- info$subsampled_depth[info$Metric==s]
							
				if(is.na(subsampled_val)){ 
					subsampled_depth <- "Not Performed"
				} else {
					subsampled_depth <- subsampled_val
				} 
									
				metric1 <- metric
				metric2 <- paste0("Subsampled Depth: ", subsampled_depth)
				
				d <- ggplot(data[data[, s] >-1,], aes_string(y=s))   + geom_boxplot()  + theme_bw()  +ylab(metric1) +ggtitle(metric2) +theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + facet_wrap(~productivity)
				e <- ggplot(data[data[, s] >-1,], aes_string(x="depths", y=s, colour="depths")) + geom_point(aes(color=depths))  + theme_bw() +xlab("ReadDepth") +ylab(metric1) + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) +stat_summary(fun.data= mean_cl_normal) + geom_smooth(method='lm')+ggtitle(metric2)+ labs(colour="Read Depth")+ facet_wrap(~productivity)

				tryCatch({
				plot(grid.arrange(d, e, ncol=2))
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
			dev.off()
			}
	}
}
			
			
			