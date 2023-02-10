
group_plot <- function(analysis_matrices, iso_type, outputdir, i, depth_file, info){
			
			productivity <- iso_type
			metric <- colnames(analysis_matrices[i])
			name<- metric
			data <- data.frame(analysis_matrices[, i])
			colnames(data) <- metric
			cols_to_plot <- colnames(data)
			data$sample <- row.names(analysis_matrices)
			
			if(outputdir %like% "SEPSIS"){
				data$sample <- gsub("HV_", "HV", data$sample)
			}
			
			## reformat sample names 
			data$sample <- gsub("BCR_", "", data$sample)
			data$sample <- gsub("TCRA_", "", data$sample)
			data$sample <- gsub("TCRB_", "", data$sample)
			data$sample <- gsub("TCRG_", "", data$sample)
			data$sample <- gsub("TCRD_", "", data$sample)
			data$sample <- gsub("_unproductive", "", data$sample)
			data$sample <- gsub("_productive", "", data$sample)
			
			#
			group <- str_split_fixed(data$sample, "_", 2)
			data$group <- group[,2]

			## Getting Read Depths 
			data <- merge(data, depth_file, by.x="sample", by.y="SampleIDforDepths")	

			## Getting Read Depths 
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
			
			#Extract the subsampling depth!
			counts_used <- paste0(outputdir, "ORIENTATED_SEQUENCES/ANNOTATIONS")
			all_files <- list.files(counts_used, full.name=TRUE)
			all_files <- grep("depth_per_isotype", all_files, value=TRUE)
			counts_used <- read.delim(all_files[1], sep="\t", header=TRUE)
			counts_used <- counts_used[counts_used$type=="UNIQ",]
			subsampled_depth_all <- counts_used$min[counts_used$X.isotype=="all"]
			
			if(dir.exists(paste0(outputdir, "Plots/ISOTYPER/GROUPED_PLOTS_", iso_type))==FALSE){
				dir.create(paste0(outputdir, "Plots/ISOTYPER/GROUPED_PLOTS_", iso_type))
			}	

			## Setting Up data Structure
			pdf(paste0(outputdir, "Plots/ISOTYPER/GROUPED_PLOTS_", iso_type, "/Plotting_GROUPS_", name, "_", subsampled_depth_all, "_", iso_type, ".pdf"), width=14, height=7)

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
				d <- ggplot(data, aes_string(x="group", y=s))   + geom_boxplot( aes(fill=group), alpha=0.5) + geom_point(aes(colour=group)) + theme_bw() +xlab("Grouping") +ylab(metric1) +ggtitle(metric2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
				tryCatch({
				plot(d)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	

				d <- ggplot(data, aes_string(x="ReadDepth", y=s, colour="ReadDepth")) + geom_point() + theme_bw() +xlab("ReadDepth") +ylab(metric1) + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) + geom_smooth(method='lm') +ggtitle(metric2)+ stat_cor(label.x = 0)
				f <- ggplot(data, aes_string(x="group", y=s))+ geom_point(aes(color=ReadDepth), size=3) + theme_bw() +xlab("Grouping") +ylab(metric1)+ scale_color_gradient(high = "yellow", low = "darkblue")+ggtitle(metric2) + labs(colour="Read Depth")+ stat_cor(label.x = 0) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
				tryCatch({
				grid.arrange(d,  f, ncol=2)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
				suppressMessages(dev.off())
				}
}