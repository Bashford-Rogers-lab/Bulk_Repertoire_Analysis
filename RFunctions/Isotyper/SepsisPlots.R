sepsis_plots <- function(analysis_matrices, iso_type, outputdir, i, depth_file, info, meta_data){
			meta_data  <- meta_data
			productivity <- iso_type
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
			data$Disease[data$sample %like% "JR1795_1003"] <- "TECHNICAL"
			
			## Getting Read Depths 
			data <- merge(data, depth_file, by.x="sample", by.y="SampleIDforDepths")	
			data <- merge(data, meta_data, by.x="sample", by.y="SampleID", all.x=TRUE)
			
			## Editing the Variables for plotting 
			data$SRS <- as.character(data$SRS)
			data$SRS[data$Disease=="HEALTH" ] <- "HEALTH"	
			data$outcome <- NA 
			data$outcome[data$Days_death_from_ICU=="Alive"] <- "ALIVE"
			data$outcome[data$Days_death_from_ICU!="Alive"] <- "DEAD"
			data$outcome[data$Disease=="HEALTH"] <- "HEALTH"
			data$outcome[data$Disease=="TECHNICAL"] <- "TECHNICAL"

			data$X7DAY_MORTALITY[data$X7DAY_MORTALITY==0] <- "7 DAY MORTALITY"
			data$X7DAY_MORTALITY[data$X7DAY_MORTALITY==1] <- "POST 7 DAY MORTALITY"
			data$X7DAY_MORTALITY[data$X7DAY_MORTALITY==2] <- "Alive"				
			data$X7DAY_MORTALITY[data$Disease=="HEALTH"] <- "HEALTH"
			data$X7DAY_MORTALITY[data$Disease=="TECHNICAL"] <- "TECHNICAL"
			
			data$X14DAY_MORTALITY <- NA 
			data$X14DAY_MORTALITY[!is.na(data$Days_death_from_ICU) & as.numeric(data$Days_death_from_ICU) <= 14] <- "2 WEEK MORTALITY"
			data$X14DAY_MORTALITY[data$Days_death_from_ICU=="Alive"] <- "ALIVE"
			data$X14DAY_MORTALITY[data$Disease=="HEALTH"] <- "HEALTH"
			data$X14DAY_MORTALITY[data$Disease=="TECHNICAL"] <- "TECHNICAL"
			data$X14DAY_MORTALITY[!is.na(data$Days_death_from_ICU) & as.numeric(data$Days_death_from_ICU) > 14] <- "POST 2 WEEK MORTALITY"
			
			data$X30DAY_MORTALITY <- NA 
			data$X30DAY_MORTALITY[!is.na(data$Days_death_from_ICU) & as.numeric(data$Days_death_from_ICU) <= 30] <- "1 MONTH MORTALITY"
			data$X30DAY_MORTALITY[data$Days_death_from_ICU=="Alive"] <- "ALIVE"
			data$X30DAY_MORTALITY[data$Disease=="HEALTH"] <- "HEALTH"
			data$X30DAY_MORTALITY[data$Disease=="TECHNICAL"] <- "TECHNICAL"
			data$X30DAY_MORTALITY[!is.na(data$Days_death_from_ICU) & as.numeric(data$Days_death_from_ICU) > 30] <- "POST 1 MONTH MORTALITY"
			
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
			
			if(dir.exists(paste0(outputdir, "Plots/ISOTYPER/SEPSIS_PLOTS_", iso_type))==FALSE){
				dir.create(paste0(outputdir, "Plots/ISOTYPER/SEPSIS_PLOTS_", iso_type))
			}	
			
			## Setting Up data Structure
			pdf(paste0(outputdir, "Plots/ISOTYPER/SEPSIS_PLOTS_", iso_type, "/Plotting_SEPSIS_", name, "_", subsampled_depth_all, "_", iso_type, "_sepsis.pdf"), width=15, height=10)

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
				b <- ggplot(data, aes_string(x="timepoint", y=s))   + geom_boxplot(aes(fill=Disease))  + theme_bw() +facet_wrap(~X14DAY_MORTALITY, scales = "free_x") +xlab("Timepoint") +ylab(metric1) +ggtitle(metric2)
				c <- ggplot(data, aes_string(x="timepoint", y=s))   + geom_boxplot(aes(fill=Disease))  + theme_bw() +facet_wrap(~X30DAY_MORTALITY, scales = "free_x") +xlab("Timepoint") +ylab(metric1) +ggtitle(metric2)
				
				
				f <- ggplot(data, aes_string(x="timepoint", y=s))  +geom_point(aes(color=SRS, shape=Sex), size=3) + theme_bw() + geom_line(aes(group = Barcode), color='black') +xlab("Timepoint") +ylab(metric1)+ggtitle(metric2)
				g2 <- ggplot(data, aes_string(x="timepoint", y=s))  +geom_point(aes(color=SRS, shape=Sex), size=3) + theme_bw() + geom_line(aes(group = Barcode), color='black') +facet_wrap(~X7DAY_MORTALITY, scales = "free_x") +xlab("Timepoint") +ylab(metric1) +ggtitle(metric2)
				b2 <- ggplot(data, aes_string(x="timepoint", y=s))  +geom_point(aes(color=SRS, shape=Sex), size=3) + theme_bw() + geom_line(aes(group = Barcode), color='black') +facet_wrap(~X14DAY_MORTALITY, scales = "free_x") +xlab("Timepoint") +ylab(metric1) +ggtitle(metric2)
				c2 <- ggplot(data, aes_string(x="timepoint", y=s))  +geom_point(aes(color=SRS, shape=Sex), size=3) + theme_bw() + geom_line(aes(group = Barcode), color='black') +facet_wrap(~X30DAY_MORTALITY, scales = "free_x") +xlab("Timepoint") +ylab(metric1) +ggtitle(metric2)


				tryCatch({
				grid.arrange(d, f, e2,  g2, ncol=2)
				}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
				
				tryCatch({
				grid.arrange(b, b2, c, c2, ncol=2)
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
				
				suppressMessages(dev.off())
				}
}
