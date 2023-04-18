## Model selection using AIC score 
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk
## September 2022
library(broom)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cowplot)
library(AICcmodavg)
library(lme4)
library(lmerTest)

model_selection <- function(eigenvectors, metadata, outputdir,type_receptor, metahealth){
	
	eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
	#####################################
	### These samples had ID mix ups and need to be removed!!!
	print("Removing Samples where RNAseq suggests a sample mixup")
	bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
	eigenvectors <- eigenvectors[!eigenvectors$sample %in% bad_ids,]
	############################################
	metadata <- read.delim(metadata, sep="\t", header=TRUE)

	## Preparing Health Metadata
	health <- read.delim(metahealth, sep="\t", header=TRUE)
	health <- health[, c("SampleID", "Barcode", "Age", "Sex")]
	colnames(health) <- c("SampleID_alternative", "Barcode", "Age", "Sex")
	health$Comorbidities.Charlson_Age <- 0 
	health$Comorbidities.Charlson_Index <- 0 
	## Add barcode to sepsis metadata 
	metadata$Barcode <- NA
		for(x in 1:length(metadata$SampleID_alternative)){
			metadata$Barcode[x] <- str_split_fixed (metadata$SampleID_alternative[x], "_", 2)[,1]
			metadata$Day[x] <- str_split_fixed (metadata$SampleID_alternative[x], "_", 2)[,2]
			if(metadata$SampleID_alternative[x] %like% "HV"){
			 metadata$Barcode[x] <- paste0(str_split_fixed(metadata$SampleID_alternative[x], "_", 3)[,1], "_", str_split_fixed(metadata$SampleID_alternative[x], "_", 3)[,2])
			}
		}
		## Add barcode to health metadata 
	health$Barcode <- NA
		for(x in 1:length(health$SampleID_alternative)){
			health$Barcode[x] <- str_split_fixed(health$SampleID_alternative[x], "_", 2)[,1]
			health$Day[x] <- str_split_fixed (health$SampleID_alternative[x], "_", 3)[,3]
			if(health$SampleID_alternative[x] %like% "HV"){
			 health$Barcode[x] <- paste0(str_split_fixed(health$SampleID_alternative[x], "_", 3)[,1], "_", str_split_fixed(health$SampleID_alternative[x], "_", 3)[,2])
			}
		}

	
	## merge meta 
	meta_all <- plyr::rbind.fill(metadata, health)
	covariates <- c("Age", "Sex",  "Comorbidities.Charlson_Index")
	metadata <- meta_all[, c("SampleID_alternative", "Barcode", "Day", covariates)]
	## Set the Effect Direction for Sex to make sure its constant!!!!!
	metadata$Sex <- factor(metadata$Sex, levels=c("Female", "Male"))
	metadata$Day <- as.numeric(metadata$Day)
	########################################################################################
	#################################################
	## Lets do sepsis first 
	## We just want day 1 
	#eigenvectors_day1 <- eigenvectors[eigenvectors$DAY=="Day1",]
	eigenvectors_day1 <- eigenvectors
	##  P value for single term is the same as a P value for F statistic so we dont need to also report this!
	
	all_models <- c()
	model_selection_best <- c()
	myplots <- list()
	myplots2 <- list()
	myplots3 <- list()

	for(i in 1:(length(colnames(eigenvectors))-3)){
		module <- colnames(eigenvectors_day1)[i]
		module <- eigenvectors_day1[,c(module, "sample")]
		
		data_all <- merge(module, metadata, by.x="sample", by.y="SampleID_alternative")
		data_all$Day <- as.numeric(data_all$Day)
		data_all$Barcode <- factor(data_all$Barcode)
		data_all$DISEASE <- "Sepsis"
		data_all$DISEASE[data_all$Barcode %like% "HV"] <- "HEALTH"
		
		###############################################################################################
		## Unrestricted model just the basic covariates 
		## Or model without covariates 
		x0 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ (1|Barcode)")), data=data_all, REML=F)
		x1.1 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ Age +(1|Barcode)")), data=data_all, REML=F)
		x1.2 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ Sex +(1|Barcode)")), data=data_all, REML=F)
		x1.3 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ Comorbidities.Charlson_Index +(1|Barcode)")), data=data_all, REML=F)
		x1.4 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE*Day +(1|Barcode)")), data=data_all, REML=F) 
		x1.5 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE+Day +(1|Barcode)")), data=data_all, REML=F) 
		x1.6 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ Age+Sex+Comorbidities.Charlson_Index +(1|Barcode)")), data=data_all, REML=F) 

		model_names1 <- c("Null", "A", "S", "C", "D*T", "DT", "ASC")
		models1 <- list(x0, x1.1, x1.2, x1.3, x1.4, x1.5, x1.6)
		
		##################################################################################
		## Model 1 covariates 
		## + 1 covariate
		x2 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE+Day+Age +(1|Barcode)")), data=data_all, REML=F) 
		x3 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE+Day+Sex +(1|Barcode)")), data=data_all, REML=F) 
		x4 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE+Day+Comorbidities.Charlson_Index +(1|Barcode)")), data=data_all, REML=F) 
		#+2 covaraites 
		x5 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE+Day+Age+Sex +(1|Barcode)")), data=data_all, REML=F) 
		x6 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE+Day+Age+Comorbidities.Charlson_Index +(1|Barcode)")), data=data_all, REML=F) 
		x7 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE+Day+Sex+Comorbidities.Charlson_Index +(1|Barcode)")), data=data_all, REML=F) 
		x8 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE+Day+Age+Sex+Comorbidities.Charlson_Index +(1|Barcode)")), data=data_all, REML=F)
			
		models2 <- list(x2, x3, x4, x5, x6, x7, x8)
		model_names2 <- c("DTA", "DTS", "DTC", "DTAS", "DTAC", "DTSC", "DTASC")

		####################################################################################################
		## Model 1 with the effect of DISEASE and DAY 
		x2 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE*Day+Age +(1|Barcode)")), data=data_all, REML=F) 
		x3 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE*Day+Sex +(1|Barcode)")), data=data_all, REML=F) 
		x4 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE*Day+Comorbidities.Charlson_Index +(1|Barcode)")), data=data_all, REML=F) 
		#+2 covaraites 
		x5 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE*Day+Age+Sex +(1|Barcode)")), data=data_all, REML=F) 
		x6 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE*Day+Age+Comorbidities.Charlson_Index +(1|Barcode)")), data=data_all, REML=F) 
		x7 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE*Day+Sex+Comorbidities.Charlson_Index +(1|Barcode)")), data=data_all, REML=F) 
		x8 <- lmerTest::lmer(formula(paste0(colnames(module)[1], "~ DISEASE*Day+Age+Sex+Comorbidities.Charlson_Index +(1|Barcode)")), data=data_all, REML=F)
			
		models3 <- list(x2, x3, x4, x5, x6, x7, x8)
		model_names3 <- c("D*TA", "D*TS", "D*TC", "D*TAS", "D*TAC", "D*TSC", "D*TASC")
		
		##############################################################################################
		## Lets put it all together
		models_all <- c(models1, models2, models3)
		model_names_all <- c(model_names1, model_names2, model_names3)
		
		#models_all <- c(models1,  models3, models4)
		#model_names_all <- c(model_names1,  model_names3, model_names4)

		## do AIC selection
		model_selection <- aictab(cand.set = models_all, modnames = model_names_all)
		model_selection_bic <- bictab(cand.set = models_all, modnames = model_names_all)

		
		## Now we want to identify the best model 		
		best_model <- data.frame(model_selection[model_selection$AICc==min(model_selection$AICc),])
		## If there are two we will just take the one with the interaction term, if still two we just take the one with fewest predictors, if still two just take the first one 
		if(dim(best_model)[1]>1){
			print("Warning more than one best model!")
			if(any(best_model$Modnames %like% "*")){
				## if model with interaction term lets use that one 
				best_model <- best_model[best_model$Modnames %like% "\\*",]
				if(dim(best_model)[1]>1){
					predictors <- nchar(best_model$Modnames)
					best_model <- best_model[nchar(best_model$Modname)==min(predictors),]
					## if still two 
					if(dim(best_model)[1]>1){
						best_model <- best_model[nchar(best_model$Modname)==min(predictors),]
					}
				}
			} else {
			predictors <- nchar(best_model$Modnames)
			best_model <- best_model[nchar(best_model$Modname)==min(predictors),]
			if(dim(best_model)[1]>1){
						best_model <- best_model[nchar(best_model$Modname)==min(predictors),]
					}
			}
		}
		best_model$Module <- colnames(module)[1]
		best_model$Modnames <- factor(best_model$Modnames, levels=c(model_names_all))
		model_selection_best <- rbind(model_selection_best, best_model)
		
		names(models_all) <- model_names_all
		id_of_bestmodel <- best_model$Modnames
		best_modal_formula <- models_all[[id_of_bestmodel]]
		
		##############################################
		## lets do the residual plot for the best model 
		#pinteraction <- interactions::interact_plot(best_modal_formula, pred = DAY, modx = DISEASE, partial.residuals = TRUE, interval = TRUE )
		################################################
		
		## Resuts for all models
		model_selection <- data.frame(model_selection)
		model_selection$Module <- colnames(module)[1]
		model_selection$Best_model <- "NO"
		model_selection$Best_model[model_selection$AICc==min(model_selection$AICc)] <- "YES"
		model_selection$Modnames <- factor(model_selection$Modnames, levels=c(model_names_all))
		all_models <- rbind(all_models, model_selection)
		
		## Lets get the plot ready 
		p1 <- ggplot(model_selection, aes(x=Module, y=Modnames)) +  geom_tile(aes(fill = as.numeric(AICc)), colour="black")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(colnames(module)[1]) + labs(fill=" AICc") +xlab("Module") +ylab("Covariate Model")+
		geom_tile(aes(color=factor(Best_model, c("YES", "NO"))), fill = '#00000000', size = 1) + scale_color_manual(name = "  Best Model", values = c("red", '#00000000'), drop=FALSE) 
		myplots[[(i)]] <- p1
		
		p2 <- ggplot(model_selection, aes(x=Module, y=Modnames)) +  geom_tile(aes(fill = as.numeric(AICcWt)), colour="black")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(colnames(module)[1]) + labs(fill=" AICcWt") +xlab("Module") +ylab("Covariate Model")+
		geom_tile(aes(color=factor(Best_model, c("YES", "NO"))), fill = '#00000000', size = 1) + scale_color_manual(name = "  Best Model", values = c("red", '#00000000'), drop=FALSE) 
		
		myplots[[(i)]] <- p1
		myplots2[[(i)]] <- p2
		#myplots3[[(i)]] <- pinteraction
	
	}
		all_models$Module <- factor(all_models$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25","Module_26","Module_27", "Module_28","Module_29","Module_30", "Module_31",  "Module_32",  "Module_33",  "Module_34",  "Module_35",  "Module_36",  "Module_37",  "Module_38",  "Module_39",  "Module_40"))
		model_selection_best$Module <- factor(model_selection_best$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25","Module_26","Module_27", "Module_28","Module_29","Module_30", "Module_31",  "Module_32",  "Module_33",  "Module_34",  "Module_35",  "Module_36",  "Module_37",  "Module_38",  "Module_39",  "Module_40"))
		model_selection_best$Modnames <- factor(model_selection_best$Modnames, levels=c(model_names_all))
		## select plotting dimensions
		if(type_receptor == "TCRGD"){
			col_number <- 3
			dims <- 15
		} else if(type_receptor=="TCRAB"){
			col_number <- 4
			dims <- 16
		} else if (type_receptor=="BCR"){
			col_number <- 5 
			dims <- 15
		}
		#########
	
		pdf(paste0(outputdir,"/Model_Assessment.pdf"), width=(dims+5), height=(dims+10))
		source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/multiplot.R')
		multiplot(plotlist = myplots, cols = 6)
		multiplot(plotlist = myplots2, cols = 6)
		dev.off()	
		
		pdf(paste0(outputdir,"/Model_Assessment_BESTMODELS.pdf"), width=11, height=5.5)
		p1 <- ggplot(model_selection_best, aes(x=Module, y=Modnames)) +  geom_tile(aes(fill = as.numeric(AICc)), colour="black")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(fill=" AICc") +xlab("Module") +ylab("Model")+  scale_y_discrete(drop = FALSE) +ggtitle("Optimal Method using AIC score")
		p2 <- ggplot(model_selection_best, aes(x=Module, y=Modnames)) +  geom_tile(aes(fill = as.numeric(AICcWt)), colour="black")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ labs(fill=" AICcWt") +xlab("Module") +ylab("Model")+  scale_y_discrete(drop = FALSE)+ggtitle("Optimal Method using AIC score")
		plot(plot_grid(p1, p2, ncol=2, align="h", axis="lbt", rel_widths=c(1,1), labels="AUTO"))
		dev.off()
		
		pdf(paste0(outputdir,"/Model_Assessment_SPREAD.pdf"), width=20, height=10)
		plot(ggplot(all_models, aes(x=AICc, fill=Modnames))+geom_histogram(binwidth=1, colour="black")+theme_classic()+facet_wrap(~Module, scales="free")+ylab("Count")+labs(fill="Model")+ggtitle("Model Selection using AICc"))
		dev.off()
		

		##Lets save test values: 
		write.table(model_selection_best, paste0(outputdir, "/Summary/Best_Models.txt"), sep="\t")
		write.table(all_models, paste0(outputdir, "/Summary/ALL_Models_AIC.txt"), sep="\t")
}
##################################################################################################################
