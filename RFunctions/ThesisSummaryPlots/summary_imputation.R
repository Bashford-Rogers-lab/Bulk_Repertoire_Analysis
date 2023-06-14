# Function to visualise the metrics for technical replicate for bulk BCR/TCR pipline to look for batch effect 
# Author: Lauren Overend 
# lauren.overend@oriel.ox.ac.uk 
# neat format for thesis 


path_to_outputdir <-  '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB'
plot_dir <- '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Plots'
stat_dir <- '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRAB/Summary'
type_use <-  "TCRAB"
iso_type <- "PRODUCTIVE"

#--------------------------------------------------------------------------------------------------------------------------------------
### Visualise imputation

summary_imputation <- function(path_to_outputdir, plot_dir, stat_dir, type_use, iso_type){
	imputation_methods <- read.delim(paste0(path_to_outputdir, "/Summary/nmrse_variables_", type_use, "_", iso_type, ".txt"), header=TRUE, sep="\t")	
	## Lets look at the average imputation score across all NAseed per method and select the best one
	method_means <- imputation_methods %>% dplyr::group_by(method) %>% dplyr::summarise(mean_error = mean(rsme))
	method_means <- data.frame(method_means)
	method_use <- method_means$method[method_means$mean_error==min(method_means$mean_error)]
	print(paste0("Best imputation method: ", method_use))
	nrmse_values <- list.files(paste0(path_to_outputdir, "/Summary"), full.names=TRUE)
	nrmse_files <- grep("nrmse", nrmse_values, value=TRUE)
	nrmse_files <- grep(type_use, nrmse_files, value=TRUE)
	nrmse_files <- grep("per_sample", nrmse_files, value=TRUE, invert=TRUE)
	nrmse_list <- lapply(nrmse_files, fread, header = TRUE, sep="\t", fill=TRUE)
	nrmse_df <- plyr::ldply(nrmse_list, data.frame)	
	nrmse_df2<- reshape2::melt(nrmse_df, id.vars=c("na_freq"))
	nrmse_df2$na_freq <- nrmse_df2$na_freq
	nrmse_df2$variable <- gsub("nrmse_", "", nrmse_df2$variable)
	colnames(nrmse_df2) <- c("NAseedfrequency", "Method", "NRMSE")
	#........................................................................................
	## Plot the result of different levels of missingness per method: 
	meaners <- imputation_methods %>% dplyr::group_by(method,NAseedfrequency)  %>% dplyr::summarize(MeanRMSE = round(mean(rsme),2))
	colnames(meaners) <-c("Imputation Method", "NA seed frequency", "Mean RMSE")
	### lets extract the stats 
	result_all <- c()
	for(i in 1:length(unique(imputation_methods$NAseedfrequency))){
		seed <- unique(imputation_methods$NAseedfrequency)[i]
		datas <- imputation_methods[imputation_methods$NAseedfrequency==seed,]
		result <- data.frame(compare_means( rsme~method, datas, ref.group=".all.", method = "wilcox.test"))
		result$NASees <- seed
		result_all <- rbind(result, result_all)
		}
	result_all[,c(1,2)] <- NULL
	colnames(result_all)[1] <- "Imputation Method"
	colnames(result_all)[7] <- "NA seed frequency"
	result_all[,4] <- NULL
	result_all <- result_all[, c(1,6,2,3,4,5)]
	result_all$p <- round(result_all$p,2)
	result_all$p.adj<- round(result_all$p.adj,2)
	############
	result_all <- merge(result_all, meaners, by=c("Imputation Method", "NA seed frequency"), sort = F)
	result_all <- result_all[,c(1,7,2,3,4,5,6)]
	
	pdf(paste0(plot_dir,'/Imputation_Stats.pdf'), width=10, height=30)
	grid.table(result_all, rows = NULL)
	dev.off()
	
	pdf(paste0(plot_dir, "/Imputation_Summary.pdf"), width=13, height=13)
	a <- ggplot(imputation_methods, aes(x=method, y=rsme, fill=method)) +geom_point(alpha = 0.7, aes(col=method), position = "jitter", set.seed(1)) +theme_bw() + geom_violin(alpha=0.5) + stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")+guides(fill="none", colour="none")+geom_hline(yintercept=0.5, col="red") +ggtitle(paste0(type_use, " Wilcox test comparing RSME against base-mean")) +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method")+ ylab("RMSE")+
    scale_y_continuous(expand = c(0.05, 0.15))
	b <- ggplot(imputation_methods, aes(x=method, y=rsme, fill=method)) +geom_point(alpha = 0.7, aes(col=method), position = "jitter", set.seed(1)) +theme_bw() + geom_boxplot(alpha=0.5)+ stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")+guides(fill="none", colour="none")+geom_hline(yintercept=0.5, col="red")+ggtitle(paste0(type_use, " Wilcox test comparing RSME against base-mean")) +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method")+ ylab("RMSE")+
    scale_y_continuous(expand = c(0.05, 0.15))
	c <- ggplot(nrmse_df2, aes(x=Method, y=NRMSE, fill=Method)) +geom_col() +theme_bw() +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method (Seed NA per column)")+ ylab("NRMSE")
	plot(a)
	plot(b)
	plot(c)
	dev.off()
	
	#########################################################
	#### Non significant difference in mean after imputation calc significance
	imputed_vs_non <- read.delim(paste0(path_to_outputdir, "/Imputed_and_ORIGINAL_DATA_FINAL_longformat_", type_use, "_", iso_type, ".txt"), sep='\t')
	if("YES" %in% imputed_vs_non$Imputed){
	
		result_imp <- c()
		for(i in 1:length(unique(imputed_vs_non$variable))){
			seed <- unique(imputed_vs_non$variable)[i]
			datas <- imputed_vs_non[imputed_vs_non$variable==seed,]
			imp <- unique(datas$Imputed)
			if(!"YES" %in% imp){
				next
			}
			
			result <- data.frame(compare_means( value~Type, datas, method = "wilcox.test"))
			result$Variable <- seed
			result_imp <- rbind(result, result_imp)
			}
		result_imp[,c(1,2,3)] <- NULL
		result_imp <- result_imp[, c(6,1,2,3,4,5)]
		result_imp$p <- round(result_imp$p,2)
		result_imp$p.adj <- NULL
		result_imp$p.format <- NULL
		
		## Which are significant
		sig_var <- result_imp$Variable[result_imp$p < 0.05]
		
		##### calculate means before and after imputation
		means_group <- imputed_vs_non %>% dplyr::group_by(variable,Type)  %>% dplyr::summarize(Mean = mean(value, na.rm=TRUE)) %>% spread(Type, Mean)
		means_group <- data.frame(means_group)
		means_group$difference <- abs((means_group$Imputed-means_group$Original)/means_group$Original)*100

		## how many nas are there 
		na_count <- imputed_vs_non[imputed_vs_non$Type=="Original",]
		d <- na_count %>% dplyr::group_by(variable) %>% rowwise() %>% mutate(isna= sum(is.na(value)))
		d <- aggregate(isna ~ variable, d, sum)
		
		## Correlate missing values to difference in mean 
		final_means <- merge(means_group, d, by="variable")
		final_means$colourx <- "ns"
		final_means$colourx[final_means$variable %in% sig_var] <- "p <0.05"
		
		pdf(paste0(plot_dir, "/Imputation_Summary_pervariable.pdf"), width=10, height=8)
		p<-ggplot(result_imp, aes(x=p)) +  geom_histogram(color="black", fill="lightgrey", bins=50)+theme_bw()+geom_vline(xintercept=0.05, col="red") +xlab("p Value for Wilcox Test ") + ylab("Number of Repertoire Features")+ggtitle(type_use)
		p1<-ggplot(means_group[means_group$difference!=0,], aes(x=difference)) +  geom_histogram(color="black", fill="lightgrey", bins=50)+theme_bw()+geom_vline(xintercept=0.00, col="red") +xlab("Abs % Change in Mean after Imputation") + ylab("Number of Repertoire Features")+ggtitle(type_use)
		p2<-ggplot(d, aes(x=isna)) +  geom_histogram(color="black", fill="lightgrey", bins=50)+theme_bw() +xlab("Number of Missing Values") + ylab("Number of Repertoire Features")+ggtitle(type_use)
		p3<-ggplot(final_means, aes(x=isna, y=difference)) +  geom_point(aes(colour=colourx))+theme_bw() +xlab("Number of Missing Values") + ylab("Abs % Change in Mean after Imputation")+scale_colour_manual(values=c("black", "red"))+labs(colour="")+ stat_cor(method = "pearson", label.x = 0, label.y = 40) + geom_smooth(method=lm, se=FALSE)+ggtitle(type_use)
		plot(plot_grid(p2, p, p1, p3, labels = c('A', 'B', 'C', 'D'), ncol=2))
		dev.off()
		
		means_group$difference <- round(means_group$difference, 4)
		means_group$Imputed <- round( means_group$Imputed, 2)
		means_group$Original <- round( means_group$Original, 2)
		
		### we want to save stats tables
		colnames(means_group)[4] <- "Abs % Change in Mean"
		means_group <- merge(means_group, d, by="variable")
		colnames(means_group)[5] <- "NA count"
		write.table(means_group, paste0(stat_dir, "/ImputationSummaryStats.txt"), sep="\t", row.names=FALSE)

		return(list(a, b, c, p2, p, p1, p3))
		} else {
			print("No values were imputed")
			return(list(a, b, c))
		}
		
}




