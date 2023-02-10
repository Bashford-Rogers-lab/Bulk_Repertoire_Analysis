## Function to plot the features from each module compared to the model
## Lauren Overend 
## Jan 2023
##############################################################

#eigenvectors <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt'
#feature_assignment <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/Clustered_Features_assignment_BCR_PRODUCTIVE_NON_IMPUTED.txt'
#imputed_data <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Imputed_DATA_FINAL_SCALED_BCR_PRODUCTIVE.txt'
#metadata <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt'
#metahealth <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Healthies_ClinData.txt'
#outputdir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR'
#loadings_use <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/Module_FeaturePCALoadings_BCR_1500_PRODUCTIVE.rds'

library(ggpubr)
library(lme4)

plot_module_vignette <- function(eigenvectors, feature_assignment, imputed_data, outputdir, loadings){
		
		eigenvectors <- read.delim(eigenvectors, header=TRUE)
		
		### We need to exclude the two bad samples 
		### These samples had ID mix ups and need to be removed!!!
		print("Removing Samples where RNAseq suggests a sample mixup")
		bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
		eigenvectors <- eigenvectors[!eigenvectors$sample %in% bad_ids,]
		##---------------
	
		feature_assignments <- read.delim(feature_assignment, header=TRUE, sep=" ")
		imputed_datas <- read.delim(imputed_data, header=TRUE, sep="\t")
		no_clusters <- length(unique(feature_assignments$Cluster))
		loadingsx <- readRDS(loadings_use)
		
		####################################################################
		## Prepare metadata 
		covariates <- c("Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index")
		metadata <- read.delim(metadata, sep="\t", header=TRUE)
		## Include for just samples that are in the dataset (after filtering etc)
		## We want the covariates and grouping for mortality 
		metadata <- metadata[metadata$SampleID_alternative %in% rownames(eigenvectors),]
		comobs <- covariates
		mort <- colnames(metadata)[colnames(metadata)%like% "Mortality"]
		outcome <- c(mort, "Infection.Number_of_ICU_Acquired_Infections", "Death.Days_from_ICU_admission")
		metadata <- metadata[, c( "SampleID_alternative", "alternative_barcode", outcome, comobs)]
		## Add barcode
		metadata$Barcode <- NA
		for(x in 1:length(metadata$SampleID_alternative)){
			metadata$Barcode[x] <- str_split_fixed (metadata$SampleID_alternative[x], "_", 2)[,1]
			if(metadata$SampleID_alternative[x] %like% "HV"){
			 metadata$Barcode[x] <- paste0(str_split_fixed(metadata$SampleID_alternative[x], "_", 3)[,1], "_", str_split_fixed(metadata$SampleID_alternative[x], "_", 3)[,2])
			}
		}
	
		## Preparing Health Metadata
		health <- read.delim(metahealth, sep="\t", header=TRUE)
		health <- health[, c("SampleID", "Barcode", "Age", "Sex")]
		colnames(health) <- c("SampleID_alternative", "Barcode", "Age", "Sex")
		health$Comorbidities.Charlson_Age <- 0 
		health$Comorbidities.Charlson_Index <- 0 
		meta_all <- plyr::rbind.fill(metadata, health)
		
		## We also need to do Day 
		for(x in 1:length(meta_all$SampleID_alternative)){
			meta_all$DAY[x] <- str_split_fixed(meta_all$Sample[x], "_", 3)[,2]
			if(meta_all$Sample[x] %like% "HV"){
				meta_all$DAY[x] <- str_split_fixed(meta_all$Sample[x], "_", 4)[,3]
			}}
		
		#####################################################
		## Add factor levels 
		meta_all$DISEASE <- "SEPSIS"
		meta_all$DISEASE[meta_all$SampleID_alternative %like% "HV"] <- "HEALTH"
		meta_all$DISEASE <- factor(meta_all$DISEASE, levels=c("HEALTH", "SEPSIS"))
		## Also add mortality info 
		meta_all$Mortality2 <- meta_all$Mortality.Classification
		meta_all$Mortality2[meta_all$Mortality2=="MID.DEATH"] <- "LATE.DEATH"
		meta_all$Mortality2 <- factor(meta_all$Mortality2, levels=c("SURVIVOR", "LATE.DEATH", "EARLY.DEATH", NA))
		meta_all$Sex <- factor(meta_all$Sex, levels=c("Female", "Male"))

		###################################
		## Ready for doing the models!!!!!!
		### We want to do a plot of the features in each module
		for(i in 1:no_clusters){
			
			if(!i %in% feature_assignments$Cluster){
				i <- "Unassigned"
				assigned <- "NO"
			} else {
				assigned <- "YES"
				loadings_cluster <- data.frame(loadingsx[[i]])
				colnames(loadings_cluster) <- "Loading"
				loadings_cluster$Feature <- rownames(loadings_cluster)
				rownames(loadings_cluster) <- NULL
				ft_cutoff <- sqrt(1/dim(loadings_cluster)[1])
			}
			
			cluster_id <- i 
			features <- feature_assignments$Feature[feature_assignments$Cluster==i]
			data_use <- imputed_datas[,c(features)]
			rownames(data_use) <- gsub("_productive", "", rownames(data_use))
			
			## Get module eigenvector score 
			## If it is for the unassigned features we will do this slightly differently 
			if(assigned=="YES"){
				module <- colnames(eigenvectors)[i]
				module_raw <- module
				print(module)
				data_ei <- data.frame(eigenvectors[, c(module)])
				rownames(data_ei) <- rownames(eigenvectors)
				colnames(data_ei) <- module
				data_ei$Sample <- rownames(data_ei)
				data_use$Sample <- rownames(data_use)
				## Merge eigenvectors and each component score 
				data_all <- merge(data_use, data_ei, by="Sample")
				## Now we want to merge with the metadata 
				data_allx <- gather(data_all, "Feature", "Scaled_Score",  -Sample, factor_key=TRUE)
				## get the module
				data_usex <- merge(data_allx, meta_all, by.x="Sample", by.y="SampleID_alternative")
			} else {
				data_use$Sample <- rownames(data_use)
				data_allx <- gather(data_use, "Feature", "Scaled_Score",  -Sample, factor_key=TRUE)
				data_usex <- merge(data_allx, meta_all, by.x="Sample", by.y="SampleID_alternative")
				module_raw <- "Unassigned_Features"
			}	
					    
			data_usex$Scaled_Score <- as.numeric(data_usex$Scaled_Score)
			data_usex$DAY <- as.numeric(data_usex$DAY)	
			
			## Add loadings score so we can use in the plot 
			if(assigned=="YES"){
				data_usex <- merge(data_usex, loadings_cluster, by="Feature", all.x=TRUE)
			}
			
			### Make the labels more user friendly 
			data_usex$Feature <- gsub("BCR_READS_", "", data_usex$Feature)
			data_usex$Feature <- gsub("__", ": ", data_usex$Feature)
			data_usex$Feature <- gsub("_", " ", data_usex$Feature)
			data_usex$Feature <- gsub("vdj", "VDJ", data_usex$Feature)
			data_usex$Feature <- gsub("IGHD/IGHM", "IGHD-IGHM", data_usex$Feature)
			features <- unique(data_usex$Feature[!data_usex$Feature %like% "Module"])
			features <- sort(features)
			module <- unique(data_usex$Feature[data_usex$Feature %like% "Module"])
			## Only one missing is the technical which is correct
			data_usex$Feature <- factor(data_usex$Feature, levels=c(module, features))
			print(paste0("There are ", length(unique(data_usex$Feature)), " Features"))
			
			#########################################################
			#Hvs D Model
			myplots_hvsd <- list()
			myplots_mort <- list()
			col_number <- ceiling(sqrt(length(features)))
			dims_pdfh <- col_number*3
			dims_pdfw <- col_number*4
			
			## Set plot dimensions
			if(dims_pdfh < 15){
				dims_pdfh == 20
			}
			
			if(dims_pdfw < 20){
				dims_pdfw == 20
			}

			p_estimatesh <- c() 
			p_estimatesm <- c() 
			pdf(paste0(outputdir,"/Module_VignettesModel_", module_raw, ".pdf"), width=dims_pdfw, height=dims_pdfh)
			for (x in 1:length(levels(data_usex$Feature))){
				#print(x)
				feature_use <- levels(data_usex$Feature)[x]
				data_feature <- data_usex[data_usex$Feature==feature_use,]
				colnames(data_feature)[colnames(data_feature) %like% "Mortality2"] <- "Mortality"
				if(assigned=="YES" & !feature_use %like% "Module"){
					score <- unique(data_feature$Loading)
					score <- round(score, 2)
					scorex <- unique(data_feature$Loading)
					if(scorex>=ft_cutoff){
						col_text <- "red"
					} else {
						col_text <- "black"
					}
				} else {
					score <- "NA"
					scorex <- "NA"
					col_text <- "black"
				}
				
				model_formula <- formula(paste0("Scaled_Score", "~", "DAY", "*","DISEASE + Age + Sex + Comorbidities.Charlson_Index+(1|Barcode)"))
				xmdl =  lmerTest::lmer(model_formula, data_feature, REML=F)
				
				estimates <- summary(xmdl)$coefficients[,5]
				estimates <- estimates[c(2:7)]
				estimates1 <- c(feature_use, estimates, "HvsS", scorex)
				hvsdp <- estimates[2]
				tp <- estimates[1]
				
				p1 <- interactions::interact_plot(xmdl, pred = DAY, modx = DISEASE, partial.residuals = TRUE, interval = TRUE, main.title= paste0(str_wrap(feature_use, width=30), "\nPCA Loading: ", score, "\nSepsis: p=", round(hvsdp,2), "\nTimepoint: p=",round(tp,2)), y.label="Scaled Score") + theme(plot.title = element_text(color = col_text, size=10))
				myplots_hvsd[[(x)]] <- p1
					
				####
				data_featurex <- data_feature[data_feature$DISEASE=="SEPSIS",]
				data_featurex$Mortality <- as.character(data_featurex$Mortality )
				data_featurex$Mortality[data_featurex$Mortality=="SURVIVOR"] <- "S"
				data_featurex$Mortality[data_featurex$Mortality=="LATE.DEATH"] <- "NS.LD"
				data_featurex$Mortality[data_featurex$Mortality=="EARLY.DEATH"] <- "NS.ED"
				
				data_featurex$Mortality <- factor(data_featurex$Mortality, levels=c("S", "NS.LD", "NS.ED"))
				
				model_formulax <- formula(paste0("Scaled_Score", "~", "DAY", "*","Mortality + Age + Sex + Comorbidities.Charlson_Index+(1|Barcode)"))
				xmdlx = lmerTest::lmer(model_formulax, data_featurex, REML=F)
				
				estimates <- summary(xmdlx)$coefficients[,5]
				estimates <- estimates[c(2:9)]
				estimates2 <- c(feature_use, estimates, "Mortality", scorex)
				
				tp <- estimates[1]
				ed <-  estimates[3]
				ld <-  estimates[2]
				
				p1x <- interactions::interact_plot(xmdlx, pred = DAY, modx = Mortality, partial.residuals = TRUE, interval = TRUE, main.title= paste0(str_wrap(feature_use, width=30), "\nPCA Loading: ", score, "\nED: p=", round(ed,2), ", LD: p=", round(ld,2), "\nTimepoint: p=",round(tp,2)), y.label="Scaled Score")+ theme(plot.title = element_text(color = col_text, size=10))
				myplots_mort[[(x)]] <- p1x
				
				## Lets R bind
				p_estimatesh	<- rbind(estimates1,	p_estimatesh)	
				p_estimatesm	<- rbind(estimates2,	p_estimatesm)	
			}
			
			## For HvsS 
			colnames(p_estimatesh) <- c("Feature", "TIMEPOINT", "DISEASE.SEPSIS", "AGE", "SEX.MALE", "COMORBIDITIES.CI", "TIMEPOINTxSEPSIS", "ANALYSIS", "LOADING")
			p_estimatesh <- data.frame(p_estimatesh)
			p_estimatesh$Module <- module_raw
			write.table(p_estimatesh, paste0(outputdir, "/Summary/HvS_ALLFEATURES_STATISTICS_", module_raw, ".txt"), sep="\t")
			p_estimatesx_h <- p_estimatesh[!p_estimatesh$Feature %like% "Module",]
			
			##For Mortality
			colnames(p_estimatesm) <- c("Feature", "TIMEPOINT", "LATE.DEATH", "EARLY.DEATH", "AGE", "SEX.MALE", "COMORBIDITIES.CI", "TIMEPOINTxLATE.DEATH","TIMEPOINTxEARLY.DEATH", "ANALYSIS", "LOADING")
			p_estimatesm <- data.frame(p_estimatesm)
			p_estimatesm$Module <- module_raw
			write.table(p_estimatesh, paste0(outputdir, "/Summary/Mortality_ALLFEATURES_STATISTICS_", module_raw, ".txt"), sep="\t")
			p_estimatesx_m <- p_estimatesm[!p_estimatesm$Feature %like% "Module",]
			
		
			## Set up plots to look whether significance correlates with module assignment
			if(assigned=="YES"){
			p1 <- ggplot(p_estimatesx_h, aes(x=as.numeric(LOADING), y=as.numeric(DISEASE.SEPSIS))) + geom_point(aes(colour=as.numeric(DISEASE.SEPSIS)<0.05))+geom_hline(yintercept=0.05, col="red") +theme_classic()+labs(colour="P < 0.05")+xlab("PCA loading Score") + ylab("P value for model effect")+ggtitle("Module Membership vs \nSignificance: \nDisease Sepsis") +                                     
			stat_smooth(method = "lm",formula = y ~ x,geom = "smooth")  + stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), label.x = min(as.numeric(p_estimatesx_h$LOADING, label.y=1)))+ theme(plot.title = element_text(color = col_text, size=10)) 
			p2 <- ggplot(p_estimatesx_h, aes(x=as.numeric(LOADING), y=as.numeric(TIMEPOINT))) + geom_point(aes(colour=as.numeric(TIMEPOINT)<0.05))+geom_hline(yintercept=0.05, col="red") +theme_classic()+labs(colour="P < 0.05")+xlab("PCA loading Score") + ylab("P value for model effect")+ggtitle("Module Membership vs \nSignificance:\nTimepoint") +                                     
			stat_smooth(method = "lm",formula = y ~ x,geom = "smooth")  + stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), label.x = min(as.numeric(p_estimatesx_h$LOADING, label.y=1))) + theme(plot.title = element_text(color = col_text, size=10))
			p3 <- ggplot(p_estimatesx_h, aes(x=as.numeric(LOADING), y=as.numeric(TIMEPOINTxSEPSIS))) + geom_point(aes(colour=as.numeric(TIMEPOINTxSEPSIS)<0.05))+geom_hline(yintercept=0.05, col="red") +theme_classic()+labs(colour="P < 0.05")+xlab("PCA loading Score") + ylab("P value for model effect")+ggtitle("Module Membership vs \nSignificance: \nTimepointxSepsis") +                                     
			stat_smooth(method = "lm",formula = y ~ x,geom = "smooth")  + stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), label.x = min(as.numeric(p_estimatesx_h$LOADING, label.y=1)))+ theme(plot.title = element_text(color = col_text, size=10)) 
			## Same but for mortality 
			p1m <- ggplot(p_estimatesx_m, aes(x=as.numeric(LOADING), y=as.numeric(LATE.DEATH))) + geom_point(aes(colour=as.numeric(LATE.DEATH)<0.05))+geom_hline(yintercept=0.05, col="red") +theme_classic()+labs(colour="P < 0.05")+xlab("PCA loading Score") + ylab("P")+ggtitle("Module Membership vs \nSignificance: \nDisease Late Death") +                                     
			stat_smooth(method = "lm",formula = y ~ x,geom = "smooth")  + stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), label.x = min(as.numeric(p_estimatesx_m$LOADING, label.y=1))) + theme(plot.title = element_text(color = col_text, size=10))
			p2m <- ggplot(p_estimatesx_m, aes(x=as.numeric(LOADING), y=as.numeric(EARLY.DEATH))) + geom_point(aes(colour=as.numeric(EARLY.DEATH)<0.05))+geom_hline(yintercept=0.05, col="red") +theme_classic()+labs(colour="P < 0.05")+xlab("PCA loading Score") + ylab("P value for model effect")+ggtitle("Module Membership vs \nSignificance: \nEarly Death") +                                     
			stat_smooth(method = "lm",formula = y ~ x,geom = "smooth")  + stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), label.x = min(as.numeric(p_estimatesx_m$LOADING, label.y=1))) + theme(plot.title = element_text(color = col_text, size=10))
			p3m <- ggplot(p_estimatesx_m, aes(x=as.numeric(LOADING), y=as.numeric(TIMEPOINT))) + geom_point(aes(colour=as.numeric(TIMEPOINT)<0.05))+geom_hline(yintercept=0.05, col="red") +theme_classic()+labs(colour="P < 0.05")+xlab("PCA loading Score") + ylab("P value for model effect")+ggtitle("Module Membership vs \nSignificance: \nTimepoint") +                                     
			stat_smooth(method = "lm",formula = y ~ x,geom = "smooth")  + stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), label.x = min(as.numeric(p_estimatesx_m$LOADING, label.y=1))) + theme(plot.title = element_text(color = col_text, size=10))
			p4m <- ggplot(p_estimatesx_m, aes(x=as.numeric(LOADING), y=as.numeric(TIMEPOINTxLATE.DEATH))) + geom_point(aes(colour=as.numeric(TIMEPOINTxLATE.DEATH)<0.05))+geom_hline(yintercept=0.05, col="red") +theme_classic()+labs(colour="P < 0.05")+xlab("PCA loading Score") + ylab("P value for model effect")+ggtitle("Module Membership vs \nSignificance: \nTimepointxLate.Death") +                                     
			stat_smooth(method = "lm",formula = y ~ x,geom = "smooth")  + stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), label.x = min(as.numeric(p_estimatesx_m$LOADING, label.y=1))) + theme(plot.title = element_text(color = col_text, size=10))
			p5m <- ggplot(p_estimatesx_m, aes(x=as.numeric(LOADING), y=as.numeric(TIMEPOINTxEARLY.DEATH))) + geom_point(aes(colour=as.numeric(TIMEPOINTxEARLY.DEATH)<0.05))+geom_hline(yintercept=0.05, col="red") +theme_classic()+labs(colour="P < 0.05")+xlab("PCA loading Score") + ylab("P value for model effect")+ggtitle("Module Membership vs \nSignificance:\nTimepointxEarly.Death") +                                     
			stat_smooth(method = "lm",formula = y ~ x,geom = "smooth")  + stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), label.x = min(as.numeric(p_estimatesx_m$LOADING, label.y=1))) + theme(plot.title = element_text(color = col_text, size=10))

			## Assign so we add to multiplot 
			num_plots <- length(myplots_hvsd)
			myplots_hvsd[[(num_plots+1)]] <- p1
			myplots_hvsd[[(num_plots+2)]] <- p2
			myplots_hvsd[[(num_plots+3)]] <- p3
			myplots_mort[[(num_plots+1)]] <- p1m
			myplots_mort[[(num_plots+2)]] <- p2m
			myplots_mort[[(num_plots+3)]] <- p3m
			myplots_mort[[(num_plots+4)]] <- p4m
			myplots_mort[[(num_plots+5)]] <- p5m
			}
			
			if(col_number^2 < length(myplots_mort)){
				col_number <- col_number+1
			}
			## This is plotting for Hvs Disease Model
			source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/multiplot.R')
			multiplot(plotlist = myplots_hvsd, cols = col_number)
			multiplot(plotlist = myplots_mort, cols = col_number)
			dev.off()
			}
			
}