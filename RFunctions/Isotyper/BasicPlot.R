## Internal functions for parrelising isotyper plotting. 
basic_plot <- function(analysis_matrices, iso_type, outputdir, i, depth_file, info){
			
			productivity <- iso_type
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
			if(!metric %like% "TRBV" & !metric %like% "TRAV" & !metric %like% "TRGV" & !metric %like% "TRDV" & !metric %like% "IGHV"){
				metric <- gsub("_",  " ", metric)
			} else {
				metric <- gsub("TRBC_", "TRBC ", metric)
				metric <- gsub("TRAC_", "TRAC ", metric)
				metric <- gsub("TRGC_", "TRGC ", metric)
				metric <- gsub("TRDC_", "TRDC ", metric)
				metric <- gsub("ALL_", "ALL ", metric)
				metric <- gsub("UNPRODUCTIVE_", "UNPRODUCTIVE ", metric)
				metric <- gsub("PRODUCTIVE_", "PRODUCTIVE ", metric)
			}	
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
			
			widthx=20
			heightx=13
			
			if(length(data$sample)<=30){
				widthx=10
				heightx=10
			}
			
			if(dir.exists(paste0(outputdir, "Plots/ISOTYPER/BASIC_PLOTS_", iso_type))==FALSE){
				dir.create(paste0(outputdir, "Plots/ISOTYPER/BASIC_PLOTS_", iso_type))
			}	
			
			
			## Setting Up data Structure
			pdf(paste0(outputdir, "Plots/ISOTYPER/BASIC_PLOTS_", iso_type, "/Plotting_", name, "_", subsampled_depth_all, "_", iso_type, ".pdf"), width=widthx, height=heightx)

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
			suppressMessages(dev.off())
			}
}
