## Code for preforming module reduction and clustering of isotyper metrics
## Based on WGCNA type analysis 
## Lauren Overend lauren.overend@oriel.ox.ac.uk
## Rachael Bashford Rogers
## Welcome Trust Centre for Human Genetics 
## 07/01/2022 

## Packages 
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
library(missForest)
library(mice)
library(VIM)
library(foreach)
library(doParallel)
library(umap)
library(psych)
library(corrplot)
library(lavaan)
library(factoextra)
library(cowplot)
library(ggpubr)
library(matrixStats)
library(RColorBrewer)
library(factoextra)
library(dqrng)
library(dynamicTreeCut)
library(moduleColor)
library(WGCNA)
library(fastcluster)


cluster_immunerep <- function(type_use, subsampled_deptha,subsampled_depthb=NA, productivity, outputdir1, outputdir2=NA, plot_outputdir, meta_data, try_imp, try_cor, run_imp){
	## Functions for Script 
	print("Running Repertoire WGCNA Reduction Pipeline. \n Author: Lauren Overend. \n Contact: Lauren.overend@oriel.ox.ac.uk")
	print("Sourcing Dependent Functions...")
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/ModuleImputationMethods.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/getcorrelationFast.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/multiplot.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/Get_Subsample_Depth.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/Cluster_Repertoire.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/Eigengenes.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/PlotModules.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/Other.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/ModuleImputationMethodsRow.R')
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/CorrelationAverage.R')
	
	## Assign Variable Names:
	## Variable Names!!!
	iso_type <- productivity
	type_use <- type_use
	outputdir <- plot_outputdir
	meta_data <- meta_data
	method_try <- try_imp ##whether to run the imputation methods 
	cor_try <- try_cor
	run_imp <- run_imp
	
	#################################
	## Reading in the Final Matrix if BCR:
	print("Reading in the Isotyper Features Matrix..")
	if(outputdir1 %like% "BCR"){
		mat_filtered <- read.delim(paste0(outputdir1, "Summary/isotyper_metrics_filtered_FINAL_ALL_", subsampled_deptha, "_", iso_type, ".txt"), sep="\t")
	}
	
	if(dir.exists(paste0(outputdir, "/Plots"))==FALSE){
		dir.create(paste0(outputdir, "/Plots"))
	}
	
	if(dir.exists(paste0(outputdir, "/Summary"))==FALSE){
		dir.create(paste0(outputdir, "/Summary"))
	}
	
	## Reading in the Final Matrix if TCR concatentate paired chain:
	if(outputdir1 %like% "TCR"){
		mat_filtered1 <- read.delim(paste0(outputdir1, "Summary/isotyper_metrics_filtered_FINAL_ALL_", subsampled_deptha , "_", iso_type, ".txt"), sep="\t")
		mat_filtered1$sample <- rownames(mat_filtered1)
		mat_filtered2 <- read.delim(paste0(outputdir2, "Summary/isotyper_metrics_filtered_FINAL_ALL_", subsampled_depthb, "_", iso_type, ".txt"), sep="\t")
		mat_filtered2$sample <- rownames(mat_filtered2)
		## Lets only keep those for which we have both TCRG/D
		mat_filtered <- merge(mat_filtered1, mat_filtered2, by="sample")
		rownames(mat_filtered) <- mat_filtered$sample
		mat_filtered$sample <- NULL
		
		## lets save the joined matrix in tensor format 
		write.table(mat_filtered, paste0(outputdir, "Summary/isotyper_metrics_filtered_FINAL_ALL_TENSOR_FORMAT_", iso_type, ".txt"), sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)
		sample_order <- rownames(mat_filtered)
		feature_order <- colnames(mat_filtered)
		write.table(sample_order, paste0(outputdir, "Summary/isotyper_metrics_filtered_FINAL_ALL_TENSOR_FORMAT_SAMPLE_ORDER", "_", iso_type, ".txt"), sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)
		write.table(feature_order, paste0(outputdir, "Summary/isotyper_metrics_filtered_FINAL_ALL_TENSOR_FORMAT_FEATURE_ORDER", "_", iso_type, ".txt"), sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
	colnames(mat_filtered) <- gsub("_PRODUCTIVE", "", colnames(mat_filtered))
	colnames(mat_filtered) <- gsub("_UNPRODUCTIVE", "", colnames(mat_filtered))
	
	
	############################
	## Compare Imputation Methods for selected features 
	## Missing-ness threshold was set in the Summary Isotyper as >60% 
	p <- mat_filtered
	features <- colnames(mat_filtered)
	
	### Lets plot the overall missingness
	widthx <- dim(mat_filtered)[1]/5
	heightx <- dim(mat_filtered)[2]/5
	pdf(paste0(outputdir, "Plots/Missingness_overall.pdf"), width=widthx, height=heightx)
	df2 <- mat_filtered %>% is.na %>% data.frame
	df2$sample <- rownames(df2)
	d <- gather(df2, metric, Present, 1:(dim(df2)[2]-1), factor_key=TRUE)
	plot(ggplot(d, aes(y=metric, x=sample, fill=Present)) + geom_tile() + 	theme_classic() + theme(axis.text.x  = element_text(angle=90, vjust=0.5)) + xlab("Sample") + ylab("Metric"))
	dev.off()
	
	### Make a directory to store in 
	if(dir.exists(paste0(outputdir, "/Plots/Imputation/"))==FALSE){
		dir.create(paste0(outputdir, "/Plots/Imputation"))
	}
	
	################################################################################
	## By Variable
	if(method_try=="YES"){
		print("Running Repertoire Imputation Optimisation....")
		print("Imputation Methods: NA seeded by Variable")
			imputation_methods <- c()
			if(type_use %like% "BCR"){
				for(i in c(0.05, 0.1, 0.15, 0.2,0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){ #, 0.3, 0.4, 0.5, 0.6
					print(i)
					u <- optimise_imputing(mat_filtered, iso_type, i, outputdir, type_use)
					imputation_methods <- rbind(imputation_methods,u)
				} 
			} else {
			  for(i in c(0.05, 0.1, 0.15, 0.2,0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){ #, 0.3, 0.4, 0.5, 0.6
					print(i)
					u <- optimise_imputing(mat_filtered, iso_type, i, outputdir, type_use)
					imputation_methods <- rbind(imputation_methods,u)
				}		
			}		
			write.table(imputation_methods, paste0(outputdir, "/Summary/nmrse_variables_", type_use, "_", iso_type, ".txt"), row.names=FALSE, sep="\t")
			print("...Done")
	}
	
	## By sample: Just use 0.4 as this is what we filtered samples on 
	if(method_try=="YES"){
		print("Imputation Methods: NA seeded by Sample")
			imputation_methods_persample <- c()
			if(type_use %like% "BCR"){
				for(i in c(0.05, 0.1, 0.15, 0.2,0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){ #, 0.3, 0.4, 0.5, 0.6
					print(i)
					u <- optimise_imputing_per_row(mat_filtered, iso_type, i, outputdir, type_use)
					imputation_methods_persample <- rbind(imputation_methods_persample,u)
				} 
			} else {
			  for(i in c(0.05, 0.1, 0.15, 0.2,0.25, 0.3, 0.35, 0.4, 0.45, 0.5)){ #, 0.3, 0.4, 0.5, 0.6
					print(i)
					u <- optimise_imputing_per_row(mat_filtered, iso_type, i, outputdir, type_use)
					imputation_methods_persample <- rbind(imputation_methods_persample,u)
				}		
			}
			
			write.table(imputation_methods_persample, paste0(outputdir, "/Summary/nmrse_variables_", type_use, "_", iso_type, "_per_sample.txt"), row.names=FALSE, sep="\t")
			print("...Done")
	} else {
		print("WARNING: NOT running imputation optimisiation. \n Imputation Optimisation Results will be loaded from previous run")
	}
	################################################################################
	## Per Variable		
	print("Generating Imputation Summary Plots...")
	imputation_methods <- read.delim(paste0(outputdir, "Summary/nmrse_variables_", type_use, "_", iso_type, ".txt"), header=TRUE, sep="\t")	
	## Lets look at the average imputation score across all NAseed per method and select the best one
	method_means <- imputation_methods %>% dplyr::group_by(method) %>% dplyr::summarise(mean_error = mean(rsme))
	method_means <- data.frame(method_means)
	method_use <- method_means$method[method_means$mean_error==min(method_means$mean_error)]
	print(paste0("Best imputation method: ", method_use))
	nrmse_values <- list.files(paste0(outputdir, "Summary"), full.names=TRUE)
	nrmse_files <- grep("nrmse", nrmse_values, value=TRUE)
	nrmse_files <- grep(type_use, nrmse_files, value=TRUE)
	nrmse_files <- grep("per_sample", nrmse_files, value=TRUE, invert=TRUE)
	nrmse_list <- lapply(nrmse_files, fread, header = TRUE, sep="\t", fill=TRUE)
	nrmse_df <- plyr::ldply(nrmse_list, data.frame)	
	nrmse_df2<- reshape2::melt(nrmse_df, id.vars=c("na_freq"))
	nrmse_df2$na_freq <- nrmse_df2$na_freq
	nrmse_df2$variable <- gsub("nrmse_", "", nrmse_df2$variable)
	colnames(nrmse_df2) <- c("NAseedfrequency", "Method", "NRMSE")
	## Plot the result of different levels of missingness per method: 
	pdf(paste0(outputdir, "Plots/Imputation_Summary_", type_use, "_", iso_type, ".pdf"), width=13, height=13)
	a <- ggplot(imputation_methods, aes(x=method, y=rsme, fill=method)) +geom_point(alpha = 0.7, aes(col=method), position = "jitter", set.seed(1)) +theme_bw() + geom_violin(alpha=0.5) + stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")+geom_hline(yintercept=0.5, col="red") +ggtitle("Wilcox test comparing RSME against base-mean") +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method")+ ylab("RMSE")
	b <- ggplot(imputation_methods, aes(x=method, y=rsme, fill=method)) +geom_point(alpha = 0.7, aes(col=method), position = "jitter", set.seed(1)) +theme_bw() + geom_boxplot(alpha=0.5)+ stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")+geom_hline(yintercept=0.5, col="red")+ggtitle("Wilcox test comparing RSME against base-mean") +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method")+ ylab("RMSE")
	c <- ggplot(nrmse_df2, aes(x=Method, y=NRMSE, fill=Method)) +geom_col() +theme_bw() +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method (Seed NA per column)")+ ylab("NRMSE")
	plot(a)
	plot(b)
	plot(c)
	dev.off()
	## Based on the imputation method lets remove columns where the rmse is greater than 0.5 at missingness values present in the dataset!!!
	print("Identifying Poorly Imputed Repertoire Features...")
	warning_columns <- imputation_methods[as.numeric(imputation_methods$rsme) > 0.5 & imputation_methods$method=="missForest",]
	nas_per_col <- sapply(mat_filtered, function(x) sum(is.na(x))) / dim(mat_filtered)[1] 
	cols_remove <- c()
	for(i in 1:length(nas_per_col)){
		check <- names(nas_per_col)[i]
		check2 <- warning_columns[warning_columns$variable ==check,]
		na_miss_val <- nas_per_col[i] 
		if(dim(check2)[1] > 0){
			if(min(check2$NAseedfrequency) <= na_miss_val){
				cols_remove <- c(cols_remove, check)
			}
		}
	}
	print(paste0("There were ", length(cols_remove), " Metrics removed due to imputation rmse values PER VARIABLE: "))
	print(cols_remove)
	### Should do a plot of missingness and highlight those above critical value
	nas_per_col <- data.frame(nas_per_col)
	nas_per_col$Metric <- rownames(nas_per_col)
	colnames(nas_per_col) <- c("Prop_NA", "Metric")
	nas_per_col$colour <- "black"
	nas_per_col$colour[nas_per_col$Metric %in% cols_remove] <- "red"
	nas_per_col <- nas_per_col[order(nas_per_col$Metric),]
	pdf(paste0(outputdir, "Plots/Metrics_RemoveduetoImputation_", type_use, "_", iso_type, ".pdf"), width=50, height=13)
	a <- ggplot(nas_per_col, aes(x=Metric, y=Prop_NA, fill=colour)) +geom_col()  +theme_bw() + labs(fill="Remove From Matrix")  + scale_fill_discrete(labels=c("KEEP", "REMOVE"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour=nas_per_col$colour)) +ylab("Proportion NA (%)")
	plot(a)
	dev.off()
	print("Done Per Variable")
	
	##-----------------------------------------------------------------------
	## Per Sample
	
	imputation_methods_persample <- read.delim(paste0(outputdir, "Summary/nmrse_variables_", type_use, "_", iso_type, "_per_sample.txt"), header=TRUE, sep="\t")	
	nrmse_values <- list.files(paste0(outputdir, "Summary"), full.names=TRUE)
	nrmse_files <- grep("nrmse", nrmse_values, value=TRUE)
	nrmse_files <- grep(type_use, nrmse_files, value=TRUE)
	nrmse_files <- grep("per_sample", nrmse_files, value=TRUE)
	nrmse_list <- lapply(nrmse_files, fread, header = TRUE, sep="\t", fill=TRUE)
	nrmse_df <- plyr::ldply(nrmse_list, data.frame)		
	nrmse_df2<- reshape2::melt(nrmse_df, id.vars=c("na_freq"))
	nrmse_df2$na_freq <- nrmse_df2$na_freq
	nrmse_df2$variable <- gsub("nrmse_", "", nrmse_df2$variable)
	colnames(nrmse_df2) <- c("NAseedfrequency", "Method", "NRMSE")
	## Plot the result of different levels of missingness per method: 
	pdf(paste0(outputdir, "Plots/Imputation_Summary_", type_use, "_", iso_type, "_per_sample.pdf"), width=13, height=13)
	a <- ggplot(imputation_methods_persample, aes(x=method, y=rsme, fill=method)) +geom_point(alpha = 0.7, aes(col=method), position = "jitter", set.seed(1)) +theme_bw() + geom_violin(alpha=0.5) + stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")+geom_hline(yintercept=0.5, col="red") +ggtitle("Wilcox test comparing RSME against base-mean") +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method")+ ylab("RMSE")
	b <- ggplot(imputation_methods_persample, aes(x=method, y=rsme, fill=method)) +geom_point(alpha = 0.7, aes(col=method), position = "jitter", set.seed(1)) +theme_bw() + geom_boxplot(alpha=0.5)+ stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")+geom_hline(yintercept=0.5, col="red")+ggtitle("Wilcox test comparing RSME against base-mean") +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method")+ ylab("RMSE")
	c <- ggplot(nrmse_df2, aes(x=Method, y=NRMSE, fill=Method)) +geom_col() +theme_bw() +facet_wrap(~NAseedfrequency)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Method (Seed NA per sample)")+ ylab("NRMSE")
	plot(a)
	plot(b)
	plot(c)
	dev.off()
	print("Done Plot")
	## Based on the imputation method lets remove columns where the rmse is greater than 0.5 at missingness values present in the dataset!!!
	warning_columns <- imputation_methods_persample[as.numeric(imputation_methods_persample$rsme) > 0.5 & imputation_methods_persample$method=="missForest",]
	#nas_per_col <- sapply(mat_filtered, function(x) sum(is.na(x))) / dim(mat_filtered)[1] 
	nas_per_row <- rowSums(is.na(mat_filtered))/dim(mat_filtered)[2]
	maximal_missingness <- max(nas_per_row)
	## Remove metrics where any row missingness is greater than the point at which nrmse breaks down
	take_out <- unique(warning_columns$variable[warning_columns$NAseedfrequency < maximal_missingness])
	##Only take out if there are actually missing values in this column (aka if this needs imputing!)
	vals <- nas_per_col[nas_per_col$Metric %in% take_out,]
	take_out <- vals$Metric[vals$Prop_NA>0]
	cols_remove_per_sample <- take_out
	print(paste0("There were ", length(cols_remove_per_sample), " Metrics removed due to imputation rmse values: PER SAMPLE: "))
	print(cols_remove_per_sample)
	### Should do a plot of missingness and highlight those above critical value
	nas_per_col <- data.frame(nas_per_col)
	nas_per_col$Metric <- rownames(nas_per_col)
	colnames(nas_per_col) <- c("Prop_NA", "Metric")
	nas_per_col$colour <- "Keep"
	nas_per_col$colour[nas_per_col$Metric %in% cols_remove_per_sample] <- "PerSample"
	nas_per_col$colour[nas_per_col$Metric %in% cols_remove] <- "PerMetric"
	nas_per_col$colour[nas_per_col$Metric %in% cols_remove & nas_per_col$Metric %in% cols_remove_per_sample] <- "PerSample&PerMetric"
	nas_per_col <- nas_per_col[order(nas_per_col$Metric),] 
	nas_per_col$colour1 <- "black"
	nas_per_col$colour1[nas_per_col$colour!="Keep"]<- "red"
	## Plot which metrics are removed 
	pdf(paste0(outputdir, "Plots/Metrics_RemoveduetoImputation_", type_use, "_", iso_type, "_per_sample.pdf"), width=50, height=13)
	a <- ggplot(nas_per_col, aes(x=Metric, y=Prop_NA, fill=colour)) +geom_col()  +theme_bw() + labs(fill="Remove From Matrix")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour=nas_per_col$colour1)) +ylab("Proportion NA (%)")
	plot(a)
	dev.off()

	print("Finished Imputation Optimisation")
	##-----------------------------------------------------------------------
	## Lets bind the cols to remove per variable and per sample 
	cols_remove <- unique(c(cols_remove,cols_remove_per_sample))
	##------------------------------------------------------------------------------------
	## These are flagged columns if the missingness goes above this we want to remove them!
	##------------------------------------------------------------------------------------
	mat_filtered <- mat_filtered[, !colnames(mat_filtered) %in% cols_remove]
	p <- mat_filtered
	##------------------------------------------------------------------------------------
	##------------------------------------------------------------------------------------
	
	################################
	# GENERATE CORRELATION MATRIX!!!!:
	## We do this only on MEASURED VALUES!
	# First we must identify depth to subsample 
	print("Commensing Correlation Matrix Generation...")
	print("Calculating Subsample Depth for Correlation Matrix")
	
	features <- colnames(mat_filtered)
	min_sample_input <- get_subsample_depth(mat_filtered, features)
	
	##### Do we have multiple samples per individual 
	if(outputdir %like% "SEPSIS"){
		print("Calculating Subsample Depth for Correlation Matrix per Timepoint")
		samples_count <- data.frame(rownames(mat_filtered))
		colnames(samples_count) <- "Sample"
		for(x in 1:length(samples_count$Sample)){
				samples_count$Barcode[x] <- str_split_fixed(samples_count$Sample[x], "_", 3)[,1]
				if(samples_count$Sample[x] %like% "HV"){
				 samples_count$Barcode[x] <- paste0(str_split_fixed(samples_count$Sample[x], "_", 3)[,1], "_", str_split_fixed(samples_count$Sample[x], "_", 3)[,2])
				}
			}
		no_samples <- length(table(samples_count$Barcode))	
		######### If we have multiple samples per individual we want to get a new minimum
		if(no_samples < length(samples_count$Sample)){
			min_sample_input_sepsis1 <- get_subsample_depth_multiplesamples(mat_filtered, features, 1)
			min_sample_input_sepsis3 <- get_subsample_depth_multiplesamples(mat_filtered, features, 3)
			min_sample_input_sepsis5 <- get_subsample_depth_multiplesamples(mat_filtered, features, 5)
			min_sample_input_sepsis <- min(min_sample_input_sepsis1 , min_sample_input_sepsis3, min_sample_input_sepsis5) 
			print(paste0("Sepsis Minimum Sample Input: ", min_sample_input_sepsis))
			}
	}
		

	##########################################################
	##### Look at how correlation changes between timepoints in sepsis 
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/getcorrelationFast.R')
	if(cor_try=="YES"){
		print("Generating Correlation Matrix for all availible data... \n WARNING: This may take a long time")
		## Calculate Correlation Between different Measures 
		## Output will be saved to the Summary Directory
		mat_feature_corr <- Get_subsample_corr_matrix_fast(features, p, min_sample_input, outputdir, iso_type, type_use, 10000)
		
		## If its like sepsis we will also try for day 1 only 
		if(outputdir %like% "SEPSIS"){
			print("Generating and comparing Correlation Matrices derived at different timepoints!")
			#...............................................
			## Time 1
			q <- p 
			q$Sample <- rownames(q)
				for(x in 1:length(q$Sample)){
					q$Day[x] <- str_split_fixed(q$Sample[x], "_", 3)[,2]
					if(q$Sample[x] %like% "HV"){
					 q$Day[x] <- str_split_fixed(q$Sample[x], "_", 4)[,3]
			}}
				###  We just want to use the first timpoint 
			q <- q[q$Day==1,]
			q$Sample <- NULL
			q$Day <- NULL
			mat_feature_corr_day1 <- Get_subsample_corr_matrix_fast(features, q, min_sample_input_sepsis, outputdir, iso_type, type_use, 10000, "Day1")
			#..............................................
			## Time 3 
			q <- p 
			q$Sample <- rownames(q)
				for(x in 1:length(q$Sample)){
					q$Day[x] <- str_split_fixed(q$Sample[x], "_", 3)[,2]
					if(q$Sample[x] %like% "HV"){
					 q$Day[x] <- str_split_fixed(q$Sample[x], "_", 4)[,3]
			}}
				###  We just want to use the first timpoint 
			q <- q[q$Day==3,]
			q$Sample <- NULL
			q$Day <- NULL
			mat_feature_corr_day3 <- Get_subsample_corr_matrix_fast(features, q, min_sample_input_sepsis, outputdir, iso_type, type_use, 10000, "Day3")
			#..............................................
			## Time 5 
			q <- p 
			q$Sample <- rownames(q)
				for(x in 1:length(q$Sample)){
					q$Day[x] <- str_split_fixed(q$Sample[x], "_", 3)[,2]
					if(q$Sample[x] %like% "HV"){
					 q$Day[x] <- str_split_fixed(q$Sample[x], "_", 4)[,3]
			}}
				###  We just want to use the first timpoint 
			q <- q[q$Day==5,]
			q$Sample <- NULL
			q$Day <- NULL
			mat_feature_corr_day5 <- Get_subsample_corr_matrix_fast(features, q, min_sample_input_sepsis, outputdir, iso_type, type_use, 10000, "Day5")
			
			############################		
			## now we want to compare p values
			mat_feature_corr_long <- data.frame(mat_feature_corr)
			mat_feature_corr_long$Chain <- rownames(mat_feature_corr_long)
			mat_feature_corr_long <- mat_feature_corr_long %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long)[3] <- "All"
			
			mat_feature_corr_long_day1 <- data.frame(mat_feature_corr_day1)
			mat_feature_corr_long_day1$Chain <- rownames(mat_feature_corr_long_day1)
			mat_feature_corr_long_day1 <- mat_feature_corr_long_day1 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day1)[3] <- "Day1"
			
			mat_feature_corr_long_day3 <- data.frame(mat_feature_corr_day3)
			mat_feature_corr_long_day3$Chain <- rownames(mat_feature_corr_long_day3)
			mat_feature_corr_long_day3 <- mat_feature_corr_long_day3 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day3)[3] <- "Day3"
				
			mat_feature_corr_long_day5 <- data.frame(mat_feature_corr_day5)
			mat_feature_corr_long_day5$Chain <- rownames(mat_feature_corr_long_day5)
			mat_feature_corr_long_day5 <- mat_feature_corr_long_day5 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day5)[3] <- "Day5"
			
			##############
			### Lets do rmcorr 
			##s <- data.frame(p)
			#s$Sample <- rownames(s)
			#for(x in 1:length(s$Sample)){
			#	s$Barcode[x] <- str_split_fixed (s$Sample[x], "_", 2)[,1]
			#	if(s$Sample[x] %like% "HV"){
			#	 s$Barcode[x] <- paste0(str_split_fixed(s$Sample[x], "_", 4)[,1], "_", str_split_fixed(s$Sample[x], "_", 4)[,2])
			#	}
			#}
			#s$Sample <- NULL
			###########
			#rmcorrbigmatt <- rmcorr_mat(Barcode, features, s)
			#to_plot <- rmcorrbigmatt$matrix 
			#to_plot <- data.frame(to_plot)
			#to_plot$Chain <- rownames(to_plot)
			#to_plotx <- to_plot %>% gather(chain2, value, -c(Chain))
			#colnames(to_plotx)[3] <- "rmcorr"
			
			##############
			newer <- merge(mat_feature_corr_long, mat_feature_corr_long_day1, by=c("Chain", "chain2"))
			newer <- merge(newer, mat_feature_corr_long_day3, by=c("Chain", "chain2"))
			newer <- merge(newer, mat_feature_corr_long_day5, by=c("Chain", "chain2"))
			#newer <-  merge(newer, to_plotx, by=c("Chain", "chain2"))
			newer$mean <- rowMeans(newer[,3:6])
			newer$sd <- as.numeric(rowSds(as.matrix(newer[,3:6])))
		
			## want to remove the rows where we are measuring the same chain 
			newer$comparison <- "Different"
			newer$comparison[newer$Chain==newer$chain2] <- "same"
			newer <- newer[!newer$comparison=="same",]
			
			########################
			cols = brewer.pal(9, "Blues")
			pal = colorRampPalette(cols)
			# Use the following line with RColorBrewer
			newer$order = findInterval(newer$sd, sort(newer$sd))
			
			pdf(paste0(outputdir, "Plots/PvalueCorrelation_", type_use, "_", iso_type, "_per_sample.pdf"), width=15, height=10)
			p1 <- ggplot(newer, aes(x=All, y=Day1)) + geom_hex(bins=500) +scale_fill_viridis_c()+theme_classic()+xlab("R All Timepoints") + ylab("R Day 1")+ stat_cor(method = "pearson", label.x = -1, label.y = 1) + geom_smooth(method='lm', col="red")
			p2 <- ggplot(newer, aes(x=All, y=Day3)) + geom_hex(bins=500) +scale_fill_viridis_c()+theme_classic()+xlab("R All Timepoints") + ylab("R Day 3")+ stat_cor(method = "pearson", label.x = -1, label.y = 1) + geom_smooth(method='lm', col="red")
			p3 <- ggplot(newer, aes(x=All, y=Day5)) + geom_hex(bins=500) +scale_fill_viridis_c()+theme_classic()+xlab("R All Timepoints") + ylab("R Day 5")+ stat_cor(method = "pearson", label.x = -1, label.y = 1) + geom_smooth(method='lm', col="red")
			p4 <- ggplot(newer, aes(x=Day1, y=Day3)) + geom_hex(bins=500) +scale_fill_viridis_c()+theme_classic()+xlab("R Day 1") + ylab("R Day 3")+ stat_cor(method = "pearson", label.x = -1, label.y = 1) + geom_smooth(method='lm', col="red")
			p5 <- ggplot(newer, aes(x=Day1, y=Day5)) + geom_hex(bins=500) +scale_fill_viridis_c()+theme_classic()+xlab("R Day 1") + ylab("R Day 5")+ stat_cor(method = "pearson", label.x = -1, label.y = 1) + geom_smooth(method='lm', col="red")
			p6 <- ggplot(newer, aes(x=Day3, y=Day5)) + geom_hex(bins=500) +scale_fill_viridis_c()+theme_classic()+xlab("R Day 3") + ylab("R Day 5")+ stat_cor(method = "pearson", label.x = -1, label.y = 1) + geom_smooth(method='lm', col="red")
			plot(plot_grid(p1, p2, p3, p4, p5, p6, ncol=3, labels = "AUTO"))
			dev.off()
			
			### we want to use just the correlation derived on day 1 
			#print("Using the correlations derived using DAY 1 only")
			#mat_feature_corr <- mat_feature_corr_day1
			print("Using the correlations derived using All Days")
			mat_feature_corr <- mat_feature_corr
			print("Done")
		}
	print("Correlation Matrix Made and Saved")
	}
	
	
	#########################################################################
	## If we have already run the correlation comparison (takes a long time) 
	if(outputdir %like% "SEPSIS"){
		print("Loading correlation from previous run file....")
		print("This is a SEPSIS ANALYSIS")
		## First lets do another plot comparing the correlations 
		day1 <- readRDS(paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_DAY1ONLY_", type_use, "_", iso_type, ".rds"))
		day3 <- readRDS(paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_DAY3ONLY_", type_use, "_", iso_type, ".rds"))
		day5 <-readRDS(paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_DAY5ONLY_", type_use, "_", iso_type, ".rds"))
		allday <-readRDS(paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_", type_use, "_", iso_type, ".rds"))
		compare_methods(day1, day3, day5, allday, outputdir)
		
		## Lets set the R matrix we want to use 
		#file_use <-  paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_DAY1ONLY_", type_use, "_", iso_type, ".rds")
		#print(paste0("File used for correlation: ", file_use))
		#mat_feature_corr <- readRDS(paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_DAY1ONLY_", type_use, "_", iso_type, ".rds"))
		
		## We want to use the R from all days 
		file_use <-  paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_", type_use, "_", iso_type, ".rds")
		print(paste0("File used for correlation: ", file_use))
		mat_feature_corr <- readRDS(paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_", type_use, "_", iso_type, ".rds"))
		
	}
	
	if(!outputdir %like% "SEPSIS"){
		print("Loading correlation from previous run file...")
		file_use <-  paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_", type_use, "_", iso_type, ".rds")
		print(paste0("File used for correlation: ", file_use))
		mat_feature_corr <- readRDS(paste0(outputdir, "Summary/Correlation_between_measures_SUBSAMPLED_", type_use, "_", iso_type, ".rds"))
	}
	
	## Cluster the correlation Matrix and Plot 
	print("Clustering the Correlation Matrix...")
	list_cluster <- cluster_features(outputdir, mat_feature_corr, type_use, iso_type, "NON_IMPUTED")
	list_cluster <- readRDS(paste0(outputdir, "Summary/ClusterLIST.rds"))
	print("Done")
	# How many variables are in each module...
	size_cluster <- lengths(list_cluster)

	## Original data(used_later) (we have to reasign it here becuase we add categorical columns which arent for imputation)
	imputed_df <- mat_filtered
	imputed_df$Type <- "Original"
	imputed_df$sample <- rownames(imputed_df)
	imputed_df$Imputed <- "NO"
	imputed_df <- reshape2::melt(imputed_df, id.vars=c("sample", "Imputed", "Type"))


	###########################################################
	## MISSING DATA IMPUTATION for final Matrix 
	## Using the missForest non-parametric missing value imputation for mixed-type data
	###################
	## If we have already done imputation and just want to try different clustering methods etc we can load from a previous save
	if(run_imp == "YES"){
		print("Imputing Missing Data according to Imputation Optimisation")
		imputed_Data <- missForest(mat_filtered, maxiter=20)
		Error <- (imputed_Data$OOBerror)
		print(paste0("NRMSE for missForest Imputed Data: ", Error))
		completeData = data.frame(imputed_Data$ximp)
		completeData$Type <- "Imputed"
		completeData$sample <- rownames(completeData)
		completeData$Imputed <- "NO"
		completeData <- reshape2::melt(completeData, id.vars=c("sample", "Imputed", "Type"))
		#Bind orignal and imputed data for visualising 
		final_data <- rbind(imputed_df, completeData)
		final_data$variable <- as.character(final_data$variable)

		###########################################################
		## Plotting the Imputed values for missForest
		imputed_samples <- final_data[is.na(final_data$value),]
		final_data$Type <- factor(final_data$Type, levels=c("Original", "Imputed"))
		for(i in 1:length(imputed_samples$sample)){
			#print(i)
			final_data$Imputed[final_data$sample==imputed_samples$sample[i] & final_data$Type=="Imputed" & final_data$variable==imputed_samples$variable[i]] <- "YES"
		} 
		plot_list <- list()
		for(i in unique(final_data$variable)){ 
			#print(i)
			datax <- final_data[final_data$variable==i,]
			datax$value <- as.numeric(datax$value)
			if(length(unique(datax$Imputed))>1){
				x <- ggplot(data = datax, aes_string(x="Type", y="value")) +geom_boxplot(aes(group = Type)) +geom_point(aes(color = as.factor(Imputed)), alpha = 0.5, position = "jitter", set.seed(1)) +scale_color_manual(values = c("red", "blue"),  drop=FALSE) +ylab(i) +theme_bw() +xlab("Imputation Type") + guides(color=guide_legend(title="Imputed Value")) + ggtitle("missForest Imputation")  + stat_compare_means(aes(label = ..p.signif..), col="darkgreen", label.x = 1.5)
				plot_list[[i]] <- x 
			}
		}
		pdf(paste0(outputdir, "Plots/ImputedValues_",  type_use, "_", subsampled_deptha, "_", iso_type,".pdf"), height=20, width=20)
		for(i in 1:(round(length(plot_list)/16)+1)){
			stop_val <- i*16
			start_val <- stop_val-15
			#print(start_val)
			#print(stop_val)
			if(start_val > length(plot_list)){
				break
			}
			new_list <- plot_list[c(start_val:stop_val)]
			multiplot(plotlist = new_list, cols = 4)
		}
		dev.off()

		###########################################################
		## Calculate Eigen Vectors!!!!!
		full_data <- data.frame(imputed_Data$ximp)
		write.table(mat_filtered, paste0(outputdir, "FILTERED_DATA_RAW_FINAL_", type_use, "_", iso_type, ".txt"), sep='\t')
		write.table(full_data, paste0(outputdir, "Imputed_DATA_FINAL_", type_use, "_", iso_type, ".txt"), sep='\t')
		write.table(final_data, paste0(outputdir, "Imputed_and_ORIGINAL_DATA_FINAL_longformat_", type_use, "_", iso_type, ".txt"), sep='\t')
		
		################################
		## scale data prior to PCA 
		full_data <- scale(full_data) 
		### we need to save the imputed data 	
		write.table(full_data, paste0(outputdir, "Imputed_DATA_FINAL_SCALED_", type_use, "_", iso_type, ".txt"), sep='\t')	
		
		############
	} else {
		print("Loading Imputed Data from Previous Run")
		## We dont want to run this step and we will just load the files from previous save so that the run is identical!!!
		full_data_unscaled <- read.delim(paste0(outputdir, "Imputed_DATA_FINAL_", type_use, "_", iso_type, ".txt"), sep='\t')
		mat_filtered <-read.delim(paste0(outputdir, "FILTERED_DATA_RAW_FINAL_", type_use, "_", iso_type, ".txt"), sep='\t')
		final_data <-read.delim(paste0(outputdir, "Imputed_and_ORIGINAL_DATA_FINAL_longformat_", type_use, "_", iso_type, ".txt"), sep='\t')
		full_data <- read.delim(paste0(outputdir, "Imputed_DATA_FINAL_SCALED_", type_use, "_", iso_type, ".txt"), sep='\t')
	}
	
	##################################
	##################################
	## Get and plot eigen genes
	print("Generating Eigenvectors")
	mat_eigenvectors <- get_eigengenes(full_data, list_cluster, outputdir,  type_use, subsampled_deptha, iso_type)
	new_eigenvectors <- mat_eigenvectors
	print("Calculated Eigenvectors and saved")
	write.table(new_eigenvectors, paste0(outputdir, "Eigenvectors_", type_use, "_", iso_type, ".txt"), sep='\t')
	
	## Plot 
	print("Plotting Modules")
	plot_eigenvectors(new_eigenvectors, meta_data, iso_type, outputdir, type_use, subsampled_deptha)
	print("Done and Saved Module Reduction")


	########################################
	########################################
	
	
	#-----------------------------------------------------------------------------
	## Can we do a sample PCA plot on eigengenes!!???
	print("Performing Sample Clustering on Module Scores")
	if(outputdir %like% "SEPSIS"){
		mtcars.pca <- prcomp(new_eigenvectors[,c(1:length(list_cluster))], center = TRUE,scale. = TRUE)
		pcs <- data.frame(mtcars.pca$x)
		pcs$DISEASE <- "SEPSIS"
		pcs$DISEASE[rownames(pcs) %like% "HV"] <- "HEALTH"
		pcs$DISEASE[rownames(pcs) %like% "JR1795_1003"] <- "TECHNICAL"
		pcs$DAY <- NA
		pcs$DAY[rownames(pcs) %like% "_1"] <- "Day1"
		pcs$DAY[rownames(pcs) %like% "_3"] <- "Day3"
		pcs$DAY[rownames(pcs) %like% "_5"] <- "Day5"
		pcs$DAY[pcs$DISEASE=="TECHNICAL"] <- "TECHNICAL"
		
		meta_datause <- read.delim(meta_data, sep='\t', header=TRUE)
		rownames(meta_datause) <- meta_datause$SampleID
		meta_datause$SRS_New <- as.factor(meta_datause$SRS_New)
		pcs <- merge(pcs, meta_datause, by = 0, all.x=TRUE)
		rownames(pcs) <- pcs$Row.names
		pcs$Row.names <- NULL
		
		pcs$X7DAY_MORTALITY[pcs$DISEASE=="TECHNICAL"] <- "T.CONTROL"
		pcs$X7DAY_MORTALITY[pcs$DISEASE=="HEALTH"] <- "HEALTH"
		pcs$X7DAY_MORTALITY[pcs$DISEASE=="HEALTH"] <- "HEALTH"
		pcs$X7DAY_MORTALITY[pcs$X7DAY_MORTALITY==0] <- "7 Day Mortality"
		pcs$X7DAY_MORTALITY[pcs$X7DAY_MORTALITY==1] <- "Post 7 Day Mortality"
		pcs$X7DAY_MORTALITY[pcs$X7DAY_MORTALITY==2] <- "Alive"
		pcs$SRS[pcs$DISEASE=="HEALTH" ] <- "HEALTH"	
		pcs$outcome <- NA 
		pcs$outcome[pcs$DISEASE=="HEALTH"] <- "HEALTHY CONTROL"
		pcs$outcome[pcs$Days_death_from_ICU=="Alive"] <- "ALIVE"
		pcs$outcome[pcs$Days_death_from_ICU!="Alive"] <- "DEAD"
		pcs$outcome[pcs$DISEASE=="TECHNICAL"] <- "T.CONTROL"
	}

	
	#-----------------------------------------------------------------------------
   	# Can we cluster samples based on their module eigenvectors
	## First generate distance tree
	## It is common to scale before hand 
	new_eigenvectors_scaled <- scale(new_eigenvectors)
	res.dist <- dist(new_eigenvectors_scaled, method = "euclidean")
	
	## Now which cluster aglomeration method 
	m <- c( "average", "single", "complete", "ward", "weighted")
	names(m) <- c( "average", "single", "complete", "ward", "weighted")
	print(paste0("Determining which cluster aglomeration method to use"))
	# function to compute coefficient
	ac <- function(x) {
	  agnes(res.dist, method = x)$ac
	}
	# Best ac score is best agglomerative coefficient 
	# Amount of clustering found!
	ac_scores <- data.frame(map_dbl(m, ac))
	colnames(ac_scores) <- "agglomerative_coefficient"
	ac_scores$Method <- rownames(ac_scores)
	method_ac_use <- ac_scores$Method[ac_scores$agglomerative_coefficient ==max(ac_scores$agglomerative_coefficient)]
	##rename for using hclust
	if(method_ac_use=="ward"){
		method_ac_use="ward.D2"
	}
	print(paste0("Cluster Agglomeration Method: ",  method_ac_use))
	
	## Plotting cluster scores per method
	pdf(paste0(outputdir, "Plots/CLUSTER_AGLOMERATION_METHOD_SAMPLECLUSTERING_", type_use, "_", iso_type, ".pdf"), height=5, width=5)
	plot(ggplot(ac_scores, aes(x=Method, y=agglomerative_coefficient, fill=Method)) +geom_col() +theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +xlab("Cluster Agglomeration Method")+ ylab("Agglomerative Coefficient"))
	dev.off()
	
	###############################################################
	## lets get optimal number of clusters
	## If you go higher than this it starts causing issues with NB clust - may need tweaking based on scenario
	print(paste0("Determining Optimal Number of Sample Clusters :"))
	pdf(paste0(outputdir, "Plots/CLUSTER_NUMBER_SAMPLECLUSTERING_NB_output", type_use, "_", iso_type, ".pdf"), height=5, width=5)
	top_clust <- 10
	bot_clust <- 2
	## Provide raw matrix and use euclidian distance and ward method 
	e <- NbClust::NbClust(new_eigenvectors_scaled, distance = "euclidean", min.nc = bot_clust, max.nc = top_clust, method = "ward.D2", index = 'all') 
	dev.off()
	## Get optimal Number of clusters (chose min if tied as we dont want loads of clusters
	No_cluster <- data.frame(e$Best.nc)
	No_cluster <- t(No_cluster[1,])
	No_cluster <- data.frame(table(No_cluster))
	maj <- No_cluster
	colnames(maj) <- c("NoClusters", "CountMethods")
	number_clusters <- maj[maj$CountMethods==max(maj$CountMethods),]
	if(dim(number_clusters)[1]==1){
		number_k <- as.numeric(as.character(number_clusters$NoClusters))
		print("Majoriy Rules")
		print(paste0("Optimum Number of clusters is ", number_k))
	} else {
		number_clusters$NoClusters <- as.numeric(as.character(number_clusters$NoClusters))
		number_k <- number_clusters$NoClusters[number_clusters$NoClusters==min(number_clusters$NoClusters)]
		print("Tied Optimal Number of Clusters - taking fewest number")
		print(paste0("Optimum Number of clusters is ", number_k))
	}
	### Lets plot NB CLust for every method 
	result2 <- data.frame(t(e$All.index))
	result2$Method <- rownames(result2)
	result2 <- reshape2::melt(result2, id.vars = c("Method"))
	colnames(result2) <- c("Method", "Number_of_Clusters", "Score")
	result2$Number_of_Clusters <- as.numeric(gsub("X", "", result2$Number_of_Clusters))
	pdf(paste0(outputdir, "Plots/SELECT_OPTIMAL_CLUSTERS_SAMPLES_", type_use, "_", iso_type, ".pdf"), height=15, width=15)
	plot(ggplot(result2, aes(x=Number_of_Clusters, y=Score, colour=Method)) + geom_point() +geom_line() + theme_bw() +facet_wrap(~Method, scales="free"))
	dev.off()
	
	## Now lets plot the dedrogram using euclidean distance and the best aglomeration method
	w=10
	pdf(paste0(outputdir, "Plots/CLUSTERING_SAMPLES_", type_use, "_", iso_type, ".pdf"), width=20, height=10)
	## Lets make the dendrogram and plot partition
	hc = hclust(res.dist, method=method_ac_use)
	par(mfrow= c(1,1), mar = c(5,5,3,3))
	plot(hc, cex=0.5, hang = -1)
	rect.hclust(hc, border = "red", k=number_k)
	dev.off()
	
	######################
	### Now lets try dynamic tree cut 
	print("Attempting Dynamic Tree Cutting: deep Split")
	pdf(paste0(outputdir, "Plots/CLUSTERING_PARTITION_Samples_DynamicTreeCut_", type_use, "_", iso_type, ".pdf"), width=40, height=10)
	dynamic_cut <-  cutreeDynamic(hc, minClusterSize = 20, method="hybrid", deepSplit = FALSE, distM=as.matrix(res.dist), maxDistToLabel = 0, verbose = 0)
	dynamic_cut2 <- labels2colors(dynamic_cut)
	plotDendroAndColors(dendro=hc,colors=dynamic_cut2, main="Feature Dendrogram: Dynamic Tree Cut")
	dev.off()
	
	
	#########	
	no_detected <- length(unique(dynamic_cut[dynamic_cut !=0]))
	print(paste0("Dynamic Tree Cut Proposes ", 	no_detected, " Clusters"))
	no_detected <- length(unique(dynamic_cut[dynamic_cut !=0]))
	if(no_detected < number_k){
		print("This is less than the number proposed by fixed height branch cut")
		print("Proceeding with Dynamic Tree Cut")
	} else if (no_detected== number_k) {
		print("This is the same as the number proposed by fixed height branch cut")
		print("Proceeding with Dynamic Tree Cut")
	} else {
		print("This is more than the number proposed by fixed height branch cut")
		print("Proceeding with Dynamic Tree Cut")
	}
	
	### Lets get assignments 
	features_x <- data.frame(rownames(as.matrix(res.dist)))
	features_x$cluster <- dynamic_cut
	features_xx <- features_x[order(features_x$cluster),]
	rownames(features_xx) <- features_xx[,1]
	features_xx[,1] <- NULL
	
	## If any is assigned a 0 it means it couldnt be assigned and we want to annotate that!	
	if(any(features_xx$cluster==0)){
		print("Unassigned features present")
		print("Will need to store seperately")
		unassigned <- data.frame(features_xx[features_xx$cluster==0,])
		rownames(unassigned) <- rownames(features_xx)[features_xx$cluster==0]
		print("Filtered unassigned features")
		## Save output so easy to check if we have any unassigned and will want to explore seperately!
		colnames(unassigned) <- c("Cluster")
		unassigned$Datatype <- "Imputed"
		unassigned$Sample <- rownames(unassigned)
		write.table(unassigned, paste0(outputdir, "Summary/UNASSIGNED_Samples_", type_use, "_", iso_type, ".txt"), row.names=FALSE)
		
	}
	colnames(features_xx) <- "Cluster"
    features_xx$Cluster[features_xx$Cluster ==0] <- "Unassigned"
	features_xx$Cluster <- as.factor(features_xx$Cluster)
	
	## Sample assignments 
	##assignment <- data.frame(e$Best.partition)
	##colnames(assignment) <- c("Cluster")
	##assignment$Cluster <- as.factor(assignment$Cluster)
	
	## Save the assignments 
	#assignment_write <- assignment
	assignment <- features_xx
	assignment_write <- features_xx
	assignment_write$Sample <- rownames(assignment_write)
	write.table(assignment_write, paste0(outputdir, "Summary/Clustered_Samples_assignment_", type_use, "_", iso_type, ".txt"), row.names=FALSE)
	
    ########################################
	## now lets merge with pcs 
	pcs <- merge(pcs, assignment, by=0)
	
	## Now lets plot 
	pdf(paste0(outputdir, "Plots/PCAEigenvectors_",  type_use, "_", subsampled_deptha, "_", iso_type,".pdf"), height=10, width=10)
	p <- ggplot(pcs, aes(x = PC1, y = PC2, colour=DISEASE)) +   geom_point()+theme_bw()
	p.1 <- ggplot(pcs, aes(x = PC1, y = PC2, colour=Cluster, shape=X7DAY_MORTALITY)) +   geom_point()+theme_bw()
	p5 <- ggplot(pcs, aes(x = PC1, y = PC2, colour=X7DAY_MORTALITY)) +   geom_point()+theme_bw()
	p2 <- ggplot(pcs, aes(x = PC1, y = PC2, colour=DISEASE)) +   geom_point()+theme_bw() + facet_wrap(~DAY)
	p3 <- ggplot(pcs, aes(x = PC1, y = PC2, colour=X7DAY_MORTALITY)) +   geom_point()+theme_bw() + facet_wrap(~DAY)
	p4 <- ggplot(pcs, aes(x = PC1, y = PC2, colour=outcome)) +   geom_point()+theme_bw() + facet_wrap(~DAY)
	p4.1	<- ggplot(pcs, aes(x = PC1, y = PC2, colour=Cluster, shape=X7DAY_MORTALITY)) +   geom_point()+theme_bw() + facet_wrap(~DAY)
	plot(p)
	plot(p2)
	plot(p.1)
	plot(p3)
	plot(p4)
	plot(p5)
	plot(p4.1)
	dev.off()
	
	print("Finished Module Reduction Pipeline")
	
}

