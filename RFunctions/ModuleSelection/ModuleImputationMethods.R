## Function to test imputation methods on repertoire data 
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk
## Jan 2022

library(VIM)
library(Hmisc)
library(missForest) 
library(data.table)
library(mice)
library(missCompare) 
library(xcms)
library(Amelia)
library(ggpubr)
library(missMDA)
library(reshape2)
library(ggplot2)
library(corrplot)
library(stringr)
library(dplyr)
library(purrr)
library(tidyr)
library(foreach)
library(doParallel) 
library(ggforce)
library(plot3D)
library(umap)
library(psych)
library(lavaan)




optimise_imputing <- function(mat_filtered, iso_type, na_freq, outputdir, type_use){
	
	print(paste0("NA FREQ: ", na_freq))
	type_use <- type_use
	## Getting Filtered Matrix into correct format! 
	p <- mat_filtered 
	
	p <- sapply(p, function(x) as.numeric(x))
	rownames(p) <- rownames(mat_filtered)
	## Okay lets generate a complete data frame with no NAs
	#remove rows with missing values in any column of data frame
	df <- data.frame(p[complete.cases(p), ])
	
	## Randomly fill the data frame with NA values at specified frequency per column rather than across whole df!!
	#df_verif0 <- prodNA(df, noNA = 0.1)
	df_verif <- data.frame(apply(df, 2, function(x) {prodNA(data.frame(x), na_freq)}))
	colnames(df_verif) <- colnames(df)

	## Plot the patterns of missingness in the data 
	pdf(paste0(outputdir, "Plots/Imputation/Missingness_patterns_", na_freq, ".pdf"), width=60, height=40)
	suppressMessages(md.pattern(p, rotate.names=TRUE, plot=TRUE))
	suppressMessages(md.pattern(df_verif, rotate.names=TRUE, plot=TRUE))
	dev.off()
	
	widthx <- dim(df)[1]/5
	heightx <- dim(df)[2]/5
	
	pdf(paste0(outputdir, "Plots/Imputation/Missingness_values_", na_freq, ".pdf"), width=widthx, height=heightx)
	df2 <- df_verif %>% is.na %>% data.frame
	df2$sample <- rownames(df2)
	d <- gather(df2, metric, Present, 1:(dim(df2)[2]-1), factor_key=TRUE)
	plot(ggplot(d, aes(y=metric, x=sample, fill=Present)) + geom_tile() + 	theme_classic() + theme(axis.text.x  = element_text(angle=90, vjust=0.5)) + xlab("Sample") + ylab("Metric"))
	dev.off()
	
	################################
	################################
	## Basic imputation with Hmisc 
	## MEAN Imputation  
	df_mean <- sapply(df_verif, function(x) Hmisc::impute(x, fun=mean))
	colnames(df_mean) <- colnames(df_verif)
	df_mean <- data.frame(df_mean)
	nrmse_mean <- missForest::nrmse(ximp=df_mean, xmis=df_verif, xtrue=df)
	df_mean$method <- "mean"
	rownames(df_mean) <- rownames(df)
	df_mean$sample <- rownames(df)
	print("Mean Imputation Done")
	
	## MEAN Imputation  with Hmisc 
	df_median <- sapply(df_verif, function(x) Hmisc::impute(x, fun=median))
	colnames(df_median) <- colnames(df_verif)
	df_median <- data.frame(df_median)
	nrmse_median <- missForest::nrmse(ximp=df_median, xmis=df_verif, xtrue=df)
	df_median$method <- "median"
	rownames(df_median) <- rownames(df)
	df_median$sample <- rownames(df)
	print("Median Imputation Done")

	## RANDOM Imputation using xcms  
	df_random <- sapply(df_verif, function(x) Hmisc::impute(as.numeric(x), fun='random'))
	colnames(df_random) <- colnames(df_random)
	df_random <- data.frame(df_random)
	nrmse_random <- missForest::nrmse(ximp=df_random, xmis=df_verif, xtrue=df)
	df_random$method <- "random"
	rownames(df_random) <- rownames(df)
	df_random$sample <- rownames(df)
	print("Random Imputation Done")

	## missMDA = regularised and EM
	x <- try(nbdim <- estim_ncpPCA(df_verif)) # estimate the number of dimensions to impute
	if(class(x)!="try-error"){
		nbdim<- nbdim[[1]]
		df_missMDA <- imputePCA(df_verif, ncp = nbdim)
		df_missMDA <- df_missMDA$completeObs
		df_missMDA1 <- data.frame(df_missMDA)
		nrmse_missMDA1 <- missForest::nrmse(ximp=df_missMDA1, xmis=df_verif, xtrue=df)
		df_missMDA1$method <- "missMDA_regularised"
		rownames(df_missMDA1) <- rownames(df)
		df_missMDA1$sample <- rownames(df)
		print("missMDA Regularised Imputation Done")
	} else {
		df_missMDA1 <- data.frame()
		nrmse_missMDA1 <- c()
	}
	
	x <- try(nbdim <- estim_ncpPCA(df_verif)) # estimate the number of dimensions to impute
	if(class(x)!="try-error"){
		nbdim <- estim_ncpPCA(df_verif) # estimate the number of dimensions to impute
		nbdim<- nbdim[[1]]
		df_missMDA <- imputePCA(df_verif, ncp = nbdim, method="EM")
		df_missMDA <- df_missMDA$completeObs
		df_missMDA2 <- data.frame(df_missMDA)
		nrmse_missMDA2 <- missForest::nrmse(ximp=df_missMDA2, xmis=df_verif, xtrue=df)
		df_missMDA2$method <- "missMDA_EM"
		rownames(df_missMDA2) <- rownames(df)
		df_missMDA2$sample <- rownames(df)
		print("missMDA EM Imputation Done")
	} else {
		df_missMDA2 <- data.frame()
		nrmse_missMDA2 <- c()
	}
	
	## Knn using VIM
	x <- try(df_knn <- VIM::kNN(df_verif))
	if(class(x)!="try-error"){
		df_knn <- df_knn %>% select(names(df_verif))
		nrmse_knn <- missForest::nrmse(ximp=df_knn, xmis=df_verif, xtrue=df)
		df_knn$method = "kNN"
		rownames(df_knn) <- rownames(df)
		df_knn$sample <- rownames(df)
		print("KNN Imputation Done")
	} else {
		df_knn <- data.frame()
		nrmse_knn <- c()
	}

	## missForest 
	x <- try(forest <- missForest(df_verif, maxiter=20, xtrue=df))
	if(class(x)!="try-error"){
		df_forest <- forest$ximp
		nrmse_forest <- missForest::nrmse(ximp=df_forest, xmis=df_verif, xtrue=df)
		df_forest$method = "missForest"
		rownames(df_forest) <- rownames(df)
		df_forest$sample <- rownames(df)
		print("missForest Imputation Done")
	} else {
		df_forest <- data.frame()
		nrmse_forest <- c()
	}
	
	## Additive Regression HMISC
	fmla <- as.formula(paste(" ~ ", paste(colnames(df_verif), collapse=" +")))
	x <- try(impute_areg <- aregImpute(formula=fmla, data = df_verif, n.impute = 5, nk=0))
	if(class(x)!="try-error"){
		tbl_imp_areg <- impute.transcan(impute_areg,imputation = 5,data = df_verif,list.out = TRUE,pr = FALSE,check = FALSE)
		tbl_imp_areg <- data.frame(matrix(unlist(tbl_imp_areg), nrow = nrow(df)))
		names(tbl_imp_areg) <- names(df)
		nrmse_imp_areg <- missForest::nrmse(ximp=tbl_imp_areg, xmis=df_verif, xtrue=df)
		tbl_imp_areg$method <- "impute_areg"
		df_imp_areg <- tbl_imp_areg
		rownames(df_imp_areg) <- rownames(df)
		df_imp_areg$sample <- rownames(df)
		print("Additive Regression Imputation Done")
	} else {
		df_imp_areg <- data.frame()
		nrmse_imp_areg <- c()
	}

	##########################################
	## Append type to original and validation data 
	## Add method 
	df_verif$method<- "Verification" 
	df_verif$sample <- rownames(df)
	df$method<-"Original"
	df$sample <- rownames(df)

	## Bind them all together
	tbl_imputations <- rbind(df, df_verif, df_mean, df_median, df_random, df_missMDA1, df_missMDA2, df_knn, df_forest, df_imp_areg)
	print("Combined All Dataframes")
	
	## NRMSE scores bound 
	## Bind nrmse
	nrmse_combined <- cbind(nrmse_mean, nrmse_median, nrmse_random, nrmse_missMDA1, nrmse_missMDA2, nrmse_knn, nrmse_forest, nrmse_imp_areg, na_freq)
	nrmse_combined <- data.frame(nrmse_combined)
	
	if (!file.exists(paste0(outputdir, "Summary/"))){
		dir.create(paste0(outputdir, "Summary/"))
	}
	
	
	## Save nrmse_combines we will need these later
	write.table(nrmse_combined, (paste0(outputdir, "Summary/nrmsescores_", type_use, "_", iso_type, "_", na_freq, ".txt")), sep="\t", row.names=FALSE)
	print("NRMSE scores saved")

	## Reshape to long format: 
	new <- reshape2::melt(tbl_imputations, id.vars=c("sample", "method"))
 
	### Now we are going to calculate differences between original values and imputed values 
	### This will allow us to calculuate RMSE for different variables 
	tbl_imp_numeric <- tbl_imputations %>% 
	  select(method, which(sapply(.,is.numeric))) %>% 
	  cbind(., id= rep(1:nrow(df), length(unique(tbl_imputations$method))))%>% 
	  gather(key = variable, value, -id, -method)

	tbl_imp_num_orig <- tbl_imp_numeric %>% 
	  filter(method == "Original") %>% 
	  select(-method) %>% 
	  dplyr::rename(value_orig = value)
	 
	tbl_imp_num_verif <- tbl_imp_numeric %>% 
	  filter(method == "Verification") %>% 
	  select(-method) %>% 
	  dplyr::rename(value_verif = value)
	  
	tbl_imp_numeric %<>%
	  filter(method %nin% c("Original", "Verification")) %>% 
	  dplyr::rename(value_imp = value) %>% 
	  inner_join(tbl_imp_num_orig, by = c("id", "variable")) %>% 
	  inner_join(tbl_imp_num_verif, by = c("id", "variable"))  

	tbl_imp_numeric_error <- tbl_imp_numeric %>% 
	  mutate(error_sq = (value_imp - value_orig) ^ 2)%>% 
	  group_by(method, variable) %>% 
	  dplyr::summarise(rsme = sqrt(sum(error_sq))/n()) 
	  
	## Reformatting for doing plots: 
	tbl_imp_numeric_error <- data.frame(tbl_imp_numeric_error)
	tbl_imp_numeric_error$NAseedfrequency <- na_freq
	
	print("Errors Calculated")

	## Plot and compare 
	meanscompare <- data.frame(compare_means(rsme ~ method,  data = tbl_imp_numeric_error))
	new_list<- list()
	for(i in 1:length(meanscompare$group1)){
		new_list[[i]]<- c(meanscompare$group1[i], meanscompare$group2[i])
	}
	
	if (!file.exists(paste0(outputdir, "Plots/"))){
		dir.create(paste0(outputdir, "Plots/"))
	}
	
	
	pdf(paste0(outputdir, "Plots/Imputation/Imputation1_", type_use, "_", iso_type, "_", na_freq, ".pdf"), width=13, height=10)
	plot(ggplot(tbl_imp_numeric_error, aes(x=method, y=rsme, fill=method)) +geom_point(alpha = 0.7, aes(col=method), position = "jitter", set.seed(1)) +theme_bw() + geom_violin(alpha=0.5) + stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")+geom_hline(yintercept=0.5, col="red") +ggtitle("Wilcox test comparing RSME against base-mean"))
	plot(ggplot(tbl_imp_numeric_error, aes(x=method, y=rsme, fill=method)) +geom_point(alpha = 0.7, aes(col=method), position = "jitter", set.seed(1)) +theme_bw() + geom_boxplot(alpha=0.5)+ stat_compare_means(ref.group = ".all.", label = "p.signif", method="wilcox.test")+geom_hline(yintercept=0.5, col="red")+ggtitle("Wilcox test comparing RSME against base-mean"))
	dev.off() 
	
	## Plot scatter plots per variable to look at difference in valyes
	variable_seq <- as.vector(seq.int(from = 1, to = length(unique(tbl_imp_numeric$variable)), by = 9))
	if(!length(unique(tbl_imp_numeric$variable)) >= max(variable_seq)){
		variable_seq[length(variable_seq)+1] <-(variable_seq[length(variable_seq)]+9)
	}	
	pdf(paste0(outputdir, "Plots/Imputation/Imputation2_Correlation_", type_use, "_",  iso_type, "_", na_freq, ".pdf"), width=15, height=15)	
	for(a in variable_seq){
		zmin = a
		ymax = a + 8
		if(ymax > dim(df_verif)[2]){
			ymax <- dim(df_verif)[2]
		}
		aa <- ggplot(tbl_imp_numeric[tbl_imp_numeric$variable %in% unique(tbl_imp_numeric$variable)[zmin:ymax],], aes(x=value_orig, y=value_imp, colour=method)) + geom_point()+facet_wrap(~variable, scales="free") +theme_bw() +geom_smooth(se=FALSE, method=lm) +xlab("Original Value") + ylab("Imputed Value")
		plot(aa)
	}
	dev.off()

	## Look at which variables have a worse error 
	pdf(paste0(outputdir, "Plots/Imputation/Imputation2_", type_use, "_", iso_type, "_", na_freq, ".pdf"), width=80, height=20)
	plot(ggplot(tbl_imp_numeric_error, aes(x=variable, y=rsme, fill=method)) + geom_col(position = "dodge") +theme_bw() +facet_wrap(~method) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=0.5, col="red"))
	dev.off() 
	print("Plots Done")

	## Summary of the values for each imputation method 
	s <- tbl_imp_numeric_error %>% dplyr::group_by(method) %>% summarise(mean_error = mean(rsme))
	return(tbl_imp_numeric_error)
	print("Data Returned")
} 
