## Internal functions for parrelising isotyper plotting. 
comparator_plot <- function(big, outputdir, i, dall, duprod, dprod){
			metric <- colnames(big[i])
			name<- metric
			data <- data.frame(big[, c(metric, 'productivity', 'sample') ])
			colnames(data) <- c(metric, 'productivity', 'sample')
			cols_to_plot <- colnames(data)[1]
			
			## Reformat the names
			data$sample <- gsub("BCR_", "", data$sample)
			data$sample <- gsub("TCRA_", "", data$sample)
			data$sample <- gsub("TCRB_", "", data$sample)
			data$sample <- gsub("TCRG_", "", data$sample)
			data$sample <- gsub("TCRD_", "", data$sample)
			data$sample <- gsub("_unproductive", "", data$sample)
			data$sample <- gsub("_productive", "", data$sample)
			
			## Identify controls
			controls <- grep("HV", data$sample, value=TRUE)
			
			# Rename metric
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
			
			## Append depth info
			data$depths <- NA
			for(i in 1:length(data$sample)){
				c <- data$sample[i]
				data$depths[data$productivity=="ALL" & data$sample==c] <- dall$ReadDepth[dall$SampleIDforDepths==c]
				data$depths[data$productivity=="UNPRODUCTIVE" & data$sample==c] <- duprod$ReadDepth[duprod$SampleIDforDepths==c]
				data$depths[data$productivity=="PRODUCTIVE" & data$sample==c] <- dprod$ReadDepth[dprod$SampleIDforDepths==c]
			}
			
			if(dir.exists(paste0(outputdir, "Plots/ISOTYPER/COMPARISON"))==FALSE){
				dir.create(paste0(outputdir, "Plots/ISOTYPER/COMPARISON"))
			}	
			
			## Setting Up data Structure
			pdf(paste0(outputdir, "Plots/ISOTYPER/COMPARISON/Plotting_", name, "_COMPARISON.pdf"), width=15, height=10)
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
									
				metric1 <- metric
				data <- data[!is.na(data[, s]),]
				data <- data.frame(data)
				d <- ggplot(data, aes_string(x=s, fill="productivity")) + rasterise(geom_density(alpha=.4), dpi=300)  + theme_bw()  +xlab(metric1) +ggtitle(metric1) 
				e <- ggplot(data, aes_string(x="depths", y=s, colour="depths")) + rasterise(geom_point(aes(color=depths)), dpi=300) + theme_bw() +xlab("ReadDepth") +ylab(metric1) + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))+ geom_smooth(method='lm')+ggtitle(metric1)+ labs(colour="Read Depth")+ facet_wrap(~productivity, scales="free_x")+ stat_cor(label.x = 0)
				tryCatch({
				grid.arrange(d, e, ncol=1)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
				suppressMessages(dev.off())
				}
}