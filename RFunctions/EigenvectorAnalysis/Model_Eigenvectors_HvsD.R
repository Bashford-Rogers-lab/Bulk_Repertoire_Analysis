## Functions to correlathealth vs disease and timepoint!
## Comparing Sepsis Day 1 to Health Day 1 function
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk
## September 2022
library(broom)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gtools)
library(jtools)
library(interactions)
library(ggpubr)
library(cowplot)
library(lme4)
library(lmerTest)
library(car)
library(multcomp)
library(lemon)
source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotTukey.R')
correlate_eigenvectors_hvd <- function(eigenvectors, metadata, outputdir, type_receptor, metahealth, run_all){
	eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
	
	### We need to exclude the two bad samples 
	### These samples had ID mix ups and need to be removed!!!
	print("Removing Samples where RNAseq suggests a sample mixup")
	bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
	eigenvectors <- eigenvectors[!eigenvectors$sample %in% bad_ids,]
	
	##------------------------------------------------------------------------------
	### PART 1
	## Preparing metadata
	covariates <- c("Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index")
	metadata <- read.delim(metadata, sep="\t", header=TRUE)
	metadata <- metadata[, c("SampleID_alternative", covariates)]
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
			if(metadata$SampleID_alternative[x] %like% "HV"){
			 metadata$Barcode[x] <- paste0(str_split_fixed(metadata$SampleID_alternative[x], "_", 3)[,1], "_", str_split_fixed(metadata$SampleID_alternative[x], "_", 3)[,2])
			}
		}
	
	## merge meta 
	meta_all <- plyr::rbind.fill(metadata, health)

	###########################################
	## Merge METADATA with eigevectors for analysis!!!!
	eigenvectors <- merge(eigenvectors, meta_all, by.x="sample", by.y="SampleID_alternative", all.x=TRUE)
	## Set the Effect Direction for Sex to make sure its constant!!!!!
	eigenvectors$Sex <- factor(eigenvectors$Sex, levels=c("Female", "Male"))
	eigenvectors$DISEASE <- factor(eigenvectors$DISEASE, levels=c("HEALTH", "SEPSIS"))

	#----------------------------------------------
	## Investigate comorbs 
	### Testing for multicolinearity 
	## This is present for charlon age and index so probably want to proceed just with the more continuous charlson index 
	pdf(paste0(outputdir,"/Comorbidity_burden.pdf"), width=7, height=4)
	comorbs_plot <- eigenvectors[, c("Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index", "DISEASE", "Barcode")]
	comorbs_plot <- unique(comorbs_plot)
	plot(ggplot(comorbs_plot, aes(x=Comorbidities.Charlson_Age, y=Comorbidities.Charlson_Index, fill=DISEASE, colour=DISEASE))+geom_point(alpha=0.5)+theme_bw()+ geom_smooth(method='lm')+ylab("Charlson Index")+xlab("Charlson Age")+ggtitle(paste0("Comorbidity Correlation"))+ stat_cor(method = "pearson", label.x = 0, label.y = 6.5)+facet_wrap(~DISEASE))
	dev.off()
	#----------------------------------------------

	#############################################################
	#### lets do a little scatter plot giving you an idea of the shape of the data
	## This is your standard plot with the 'real data'
	plots <- eigenvectors %>% gather(module, score, -c(sample, DAY, DISEASE, Age, Sex, Comorbidities.Charlson_Age, Comorbidities.Charlson_Index, Barcode)) 
	plots$score <- as.numeric(plots$score)
	plots$DAY <- gsub("Day", "", plots$DAY)
	plots$DAY <- as.numeric(plots$DAY)
	plots$module <- factor(plots$module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25", "Module_26", "Module_27", "Module_28", "Module_29", "Module_30",  "Module_31",  "Module_32",  "Module_33",  "Module_34",  "Module_35",  "Module_36",  "Module_37",  "Module_38",  "Module_39",  "Module_40"))
	
	pdf(paste0(outputdir,"/HvsD_Module_Scatter.pdf"), width=13, height=13)
	plot(ggplot(plots, aes(x=DAY, y=score, fill=DISEASE, colour=DISEASE))+geom_point(alpha=0.5)+facet_wrap(~module, scales="free") +theme_bw()+ geom_smooth(method='lm')+ylab("Module Score")+xlab("Day")+ggtitle(paste0("Interaction Plots ", type_receptor, " : Module~DISEASE")))
	dev.off()
	
	#........................................................
	#........................................................
	## RUNNING AN INTERACTION MODEL 
	p_values_hd <- c()
	vif_all <- c()
	myplots <- list()
	behaviour <- c()
	
	## Run a simple model to look for interaction
	for(i in 2:(length(colnames(eigenvectors))-7)){
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE", "Age", "Sex", "Barcode","Comorbidities.Charlson_Index", "Comorbidities.Charlson_Age")]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		
		## Model with Interaction and model with just covariates to do a test on 
		model_formula <- formula(paste0(colnames(module)[1], "~", "DAY*DISEASE + Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))
		model_formula_basic <- formula(paste0(colnames(module)[1], "~"," 1 +(1|Barcode)"))
		model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))

		xmdl =  lmerTest::lmer(model_formula, module, REML=F)
		xmdl_old <- lmerTest::lmer(model_formula_basic, module, REML=F)
		xmdl_cov <- lmerTest::lmer(model_formula_covonly, module, REML=F)

		## VIF
		vif_scores <- vif(xmdl)
		vif_scores <- c(colnames(eigenvectors)[i], vif_scores)	
		vif_all <- rbind(vif_all, vif_scores)
		
		## Extract the model p vals for each term 
		p_vals <- summary(xmdl)$coefficients[,5]
		p_vals <- p_vals[2:length(p_vals)]
		
		## Actually we will calculate the delta AIC here 
		aic_newmod <- AIC(xmdl)
		aic_oldmod <- AIC(xmdl_old)
		delta_aic <- aic_newmod-aic_oldmod
		## test difference in model to NULL MODEL 
		compare_mod <- anova(xmdl, xmdl_old)
		p_overall <- compare_mod[,'Pr(>Chisq)'][2]
		## test difference in model to Covariate MODEL
		compare_mod_c <- anova(xmdl, xmdl_cov)
		p_overall_c <- compare_mod_c[,'Pr(>Chisq)'][2]
		
		################################
		
		row_data <- c(colnames(module)[1], p_vals, p_overall, p_overall_c)
		row_data <- data.frame(t(row_data))
		row_data[,2:dim(row_data)[2]] <-as.numeric(row_data[,2:dim(row_data)[2]])
		colnames(row_data) <- c("Module", paste0("IM.", names(p_vals)), "ANOVA_NULL", "ANOVA_COV") 
		
		## Can we do an interaction plot of residuals
		p1 <- interactions::interact_plot(xmdl, pred = DAY, modx = DISEASE, partial.residuals = TRUE, y.label="Residual\nModule Score", interval = TRUE, main= colnames(eigenvectors)[i])
		myplots[[(i-1)]] <- p1
		p_values_hd <- rbind(row_data, p_values_hd)
		
		## create a dataframe storing effect information
		estimates <- summary(xmdl)$coefficients[,1]
		estimates <- estimates[2:length(estimates)]
		estimates <- data.frame(estimates)
		colnames(estimates) <- c("Estimate")
		estimates$Effect <- rownames(estimates)
		estimates$Module <- colnames(module)[1]
		estimates$Direction <- NA
		estimates$Estimate <- as.numeric(as.character(estimates$Estimate))
		estimates$Direction[estimates$Estimate > 0] <- "increase"
		estimates$Direction[estimates$Estimate < 0] <- "decrease"
		estimates$Direction[estimates$Estimate == 0] <- "no.change"
		estimates$Effect[estimates$Effect=="DISEASESEPSIS"] <- "DISEASE.SEPSIS"
		estimates$Effect[estimates$Effect=="SexMale"] <- "SEX.MALE"
		estimates$Effect[estimates$Effect=="DAY"] <- "TIMEPOINT"
		estimates$Effect[estimates$Effect=="Comorbidities.Charlson_Index"] <- "COMORBIDITIES:CI"
		estimates$Effect[estimates$Effect=="DAY:DISEASESEPSIS"] <- "TIMEPOINTxSEPSIS"
		estimates$Effect[estimates$Effect=="Age"] <- "AGE"
		estimates$Direction[estimates$Effect %like% "TIMEPOINTx"] <- NA

		behaviour <- rbind(behaviour, estimates)	
	}
	colnames(p_values_hd) <- c("Module", "TIMEPOINT", "DISEASE.SEPSIS", "AGE", "SEX.MALE", "COMORBIDITIES:CI", "TIMEPOINTxSEPSIS", "ANOVA_NULL", "ANOVA_COV") 	
	
	#--------------------------------------------------------
	### Lets quickly plot the model VIF scores (is colinearity a problem)	
	pdf(paste0(outputdir,"/HvsD_VIF_Score.pdf"), width=5, height=5)
	vif_all <- data.frame(vif_all)
	vif_all_g <- gather(vif_all, "Effect", "VIF", DAY:DAY.DISEASE, factor_key=TRUE)
	vif_all_g <- vif_all_g[vif_all_g$V1 =="Module_1",]
	vif_all_g$VIF <- as.numeric(vif_all_g$VIF)
	vif_all_g$Effect <- as.character(vif_all_g$Effect)
	vif_all_g$Effect[vif_all_g$Effect=="Comorbidities.Charlson_Index"] <- "COMORBIDITIES:CI"
	vif_all_g$Effect[vif_all_g$Effect=="DAY.DISEASE"] <- "DAY*DISEASE"
	vif_all_g$Effect[vif_all_g$Effect=="Age"] <- "AGE"
	vif_all_g$Effect[vif_all_g$Effect=="Sex"] <- "SEX"
	vif_all_g$Effect <- factor(vif_all_g$Effect, levels=c("DAY", "DISEASE", "DAY*DISEASE", "AGE", "SEX", "COMORBIDITIES:CI"))
	plot(ggplot(vif_all_g, aes(x=Effect, y=VIF))+geom_col(alpha=0.5) +theme_classic()+ylab("Module Score")+xlab("Day")+ggtitle(paste0("Interaction Model VIFF Scores")) +facet_wrap(~V1)+geom_hline(yintercept=5, col="red") + xlab("Model Effect") + ylab("VIFF")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
	dev.off()
	#--------------------------------------------------------
	
	## How many columns for doing the Residual Plots!!!!
	if(type_receptor == "TCRGD"){
		col_number <- 3
		dims <- 15
	} else if(type_receptor=="TCRAB"){
		col_number <- 4
		dims <- 16
	} else if (type_receptor=="BCR"){
		col_number <- 4 
		dims <- 15
	}
	## Do residual plots to visualise corrected values after adjusting for covariates. 
	pdf(paste0(outputdir, "/HvsD_AdjModel_Residual_Plot.pdf"), width=dims, height=dims)
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/multiplot.R')
	multiplot(plotlist = myplots, cols = col_number)
	dev.off()	
	
	## GATHER_TEST_RESULTS_for doing the interaction plot
	p_hd <- gather(p_values_hd, "Effect", "p_value", TIMEPOINT:ANOVA_COV, factor_key=TRUE)
	p_hd$Module <- factor(p_hd$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25",  "Module_26",  "Module_27",  "Module_28",  "Module_29",  "Module_30",  "Module_31"))
	## Lets PLOT The statistics for them  
	p_use <- 0.05
	p_hd$p_value <- as.numeric(p_hd$p_value)
	
	## Now lets do a BH correction
	p_hd$BH <- p.adjust(as.numeric(p_hd$p_value),method="BH")
	p_hd$psig <- "ns"
	p_hd$psig[p_hd$p_value < 0.05] <- "sig"
	p_hd$BHsig <- "ns"
	p_hd$BHsig[p_hd$BH < 0.05] <- "sig"
	p_hd_sub <- p_hd[p_hd$BHsig =="sig",]
	p_hd$Type <- "INTERACTION MODEL"
	p_hd$Type[p_hd$Effect %like% "ANOVA"] <- "ANOVA"
	p_hd$Effect <- as.character(p_hd$Effect)
    p_hd$Effect[p_hd$Effect =="ANOVA_NULL"] <- "NULL MODEL"
	p_hd$Effect[p_hd$Effect =="ANOVA_COV"] <- "COV ONLY MODEL"
	p_hd$Effect <- factor(p_hd$Effect, levels=c("TIMEPOINT", "DISEASE.SEPSIS", "TIMEPOINTxSEPSIS", "AGE", "SEX.MALE", "COMORBIDITIES:CI", "NULL MODEL", "COV ONLY MODEL"))
	
	#### SAVE INTERACTION STATS
	write.table(p_hd, paste0(outputdir, "/Summary/HvsD_INTERACTION_ANALYSIS_SUMMARY.txt"), sep="\t")

	## Lets add in the effect direction so that we can add an arrow on the plot indicating direction
	p_hd <- merge(p_hd, behaviour, by=c("Module", "Effect"), all.x=TRUE)
	## Lets remove the direction for ns effect
	p_hd$Direction[p_hd$psig=="ns"] <- NA
	p_hd$Direction <- factor(p_hd$Direction, levels=c("increase", "decrease", NA))
	p_hd$psig <- factor(p_hd$psig, levels=c("sig", "ns", NA))
	p_hd$BHsig <- factor(p_hd$BHsig, levels=c("sig", "ns", NA))
	p_hd$Type <- str_wrap(p_hd$Type, width = 12)
	
	## Do a plot showing BHP vs P
	px <- ggplot(p_hd) +  geom_point(aes(x=p_value, y=BH, alpha=psig, colour = BHsig)) +theme_classic() +labs(alpha="P<0.05", colour="BH P<0.05")+xlab("P") +ylab("Benjamin Hochburg Adjusted P") + scale_color_manual(values = c( 'red','lightblue')) +geom_vline(xintercept=0.05, col="blue") + geom_hline(yintercept=0.05, col="blue")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE) + ggtitle("Multiple Testing \n Correction")
	
	pdf(paste0(outputdir, "/HvsD_InteractionModel_Statistics.pdf"), width=11, height=6)
	p1 <- ggplot(p_hd, aes(x=Module, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Interaction Model, P<0.05")) + labs(fill="P") +xlab("Module") +ylab("Model Coefficients")+geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type),  scales="free_y", space='free_y')+ geom_point(aes(x=Module, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction\nof Effect\nRelative to:\nA/B. Health\nC. Day1", alpha="P<0.05") + scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
	plot(plot_grid(p1, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))
	
	p2 <- ggplot(p_hd, aes(x=Module, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, p_use))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("A. ", type_receptor, " Interaction Model, P<0.05")) + labs(fill="P")+xlab("Module") +ylab("Model Coefficients")+geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type),  scales="free_y", space='free_y')+ geom_point(aes(x=Module, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction\nof Effect\nRelative to:\nA/B. Health\nC. Day1", alpha="P<0.05") + scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
	plot(plot_grid(p2, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))	
	p2x <- p2
	
	p3 <- ggplot(p_hd, aes(x=Module, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, 0.1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("A. ", type_receptor, " Interaction Model, P<0.05")) + labs(fill="P") +xlab("Module") +ylab("Model Coefficients")+geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type),  scales="free_y", space='free_y')+ geom_point(aes(x=Module, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction\nof Effect\nRelative to:\nA/B. Health\nC. Day1", alpha="P<0.05") + scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
	plot(plot_grid(p3, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))		
	p3x <- p3
	dev.off()
	
	##------------------------------------------------------------------------------------------------------
	## Does health vs disease significantly interact for any models??
	## If so we will need to investigate these modules at individual levels. 
	if(any(p_hd$p_value[p_hd$Effect=="TIMEPOINTxSEPSIS"] < 0.05)){
		print("Presence of significant interaction")
		interaction_pres <- "YES"
		modules_interaction <- p_values_hd$Module[p_values_hd$TIMEPOINTxSEPSIS < 0.05]
		annot <- "SignificantIntOnly"
	} else {
		print("No interaction effect between timepoint and disease status")
		interaction_pres <- "NO"	
	}
	
	if(run_all=="YES"){
		modules_interaction <- p_values_hd$Module
		annot <- "AllModules"
		interaction_pres <- "YES"
	}
	
	##.........................................................................................................
	##.........................................................................................................
	##.........................................................................................................
	####################################################################################################
	### Part 2 more complex model IF an interaction is present 
	### Because of the interaction we need to assess at levels of Time / Disease state rather than across
	if(interaction_pres =="YES"){
		p_hd_extended_hvsd <- c()
		behaviour_hvsd <- c()
		p_hd_extended_timepoint <- c()
		behaviour_timepoint <- c()
		
		for(i in 1:length(modules_interaction)){
			#-----------------------------------------------------------------
			# EFFECT of TIME 
			####################################################
			## Lets do Health first
			module <- modules_interaction[i]
			module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE", "Age", "Sex", "Barcode", "Comorbidities.Charlson_Age","Comorbidities.Charlson_Index")]
			module <- module[module$DISEASE=="HEALTH",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			## Run model to see effect of time ## using it as a continuous variable? 
			## Note there will be no effect of CHI becuase it is 0 in health 
			model_formula <- formula(paste0( modules_interaction[i], "~", "DAY+Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)")) #+Barcode
			model_formula_basic <- formula(paste0( modules_interaction[i], "~"," 1 +(1|Barcode)"))
			model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))

			xmdl =  lmerTest::lmer(model_formula, module, REML=F)
			xmdl_old <- lmerTest::lmer(model_formula_basic, module, REML=F)
			xmdl_cov <- lmerTest::lmer(model_formula_covonly, module, REML=F)
	
			## We want to take out just the relevant p values 
			p_vals <- summary(xmdl)$coefficients[,5]
			## get p val just for day (p overall will be the same!)
			p_vals <- p_vals[2:4]
			
			## Actually we will calculate the delta AIC here 
			aic_newmod <- AIC(xmdl)
			aic_oldmod <- AIC(xmdl_old)
			delta_aic <- aic_newmod-aic_oldmod
			## test difference in model to NULL MODEL 
			compare_mod <- anova(xmdl, xmdl_old)
			p_overall <- compare_mod[,'Pr(>Chisq)'][2]
			## test difference in model to Covariate MODEL
			compare_mod_c <- anova(xmdl, xmdl_cov)
			p_overall_c <- compare_mod_c[,'Pr(>Chisq)'][2]

		
			###################
			p_vals <- c(p_vals, NA, p_overall, p_overall_c)
		
			## Get estimates 
			behaviour <- summary(xmdl)$coefficients[,1]
			behaviour <- behaviour[2:length(behaviour)]
			behaviour <- data.frame(behaviour)
			colnames(behaviour) <- c("Estimate")
			behaviour$Effect <- rownames(behaviour)
			behaviour$Module <- colnames(module)[1]
			behaviour$DISEASE <- "Health.Timepoint"
			behaviour$Direction <- NA
			behaviour$Estimate <- as.numeric(as.character(behaviour$Estimate))
			behaviour$Direction[behaviour$Estimate > 0] <- "increase"
			behaviour$Direction[behaviour$Estimate < 0] <- "decrease"
			behaviour$Direction[behaviour$Estimate == 0] <- "no.change"
			### Reformat 
			behaviour$Effect[behaviour$Effect=="DISEASESEPSIS"] <- "DISEASE.SEPSIS"
			behaviour$Effect[behaviour$Effect=="SexMale"] <- "SEX.MALE"
			behaviour$Effect[behaviour$Effect=="DAY"] <- "TIMEPOINT"
			behaviour$Effect[behaviour$Effect=="Comorbidities.Charlson_Index"] <- "COMORBIDITIES.CI"
			behaviour$Effect[behaviour$Effect=="DAY:DISEASESEPSIS"] <- "TIMEPOINTxSEPSIS"
			behaviour$Effect[behaviour$Effect=="Age"] <- "AGE"
			behaviour_timepoint <- rbind(behaviour_timepoint, behaviour)
			
			##################################################################################################################
			#################################################################################################################
			## Lets do Sepsis 
			module <- modules_interaction[i]
			module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE", "Age", "Sex", "Barcode","Comorbidities.Charlson_Index")]
			## Lets subset for health only  ## Adding in sample this time!
			module <- module[module$DISEASE=="SEPSIS",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			## Run model to see effect of time ## using it as a continuous variable? 
			model_formula <- formula(paste0(colnames(module)[1], "~", "DAY+Age + Sex  + Comorbidities.Charlson_Index+(1|Barcode)")) # Barcode
			model_formula_basic <- formula(paste0( modules_interaction[i], "~"," 1 +(1|Barcode)"))
			model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))

			xmdl =  lmerTest::lmer(model_formula, module, REML=F)
			xmdl_old <- lmerTest::lmer(model_formula_basic, module, REML=F)
			xmdl_cov <- lmerTest::lmer(model_formula_covonly, module, REML=F)

			## We want to take out just the relevant p values 
			p_vals_s <- summary(xmdl)$coefficients[,5]
			## get p val just for day (p overall will be the same!)
			p_vals_s <- p_vals_s[2:5]
			
			## Actually we will calculate the delta AIC here 
			aic_newmod <- AIC(xmdl)
			aic_oldmod <- AIC(xmdl_old)
			delta_aic <- aic_newmod-aic_oldmod
			## test difference in model to NULL MODEL 
			compare_mod <- anova(xmdl, xmdl_old)
			p_overall <- compare_mod[,'Pr(>Chisq)'][2]
			## test difference in model to Covariate MODEL
			compare_mod_c <- anova(xmdl, xmdl_cov)
			p_overall_c <- compare_mod_c[,'Pr(>Chisq)'][2]
		
			p_vals_s <- c(p_vals_s, p_overall, p_overall_c)	
			
			#################################
			p_h <- c(modules_interaction[i], "Health.Timepoint", p_vals)
			p_s <- c(modules_interaction[i], "Sepsis.Timepoint", p_vals_s)
			p_all <- data.frame(rbind(p_h, p_s))
			colnames(p_all) <- c("MODULE", "DISEASE", "TIMEPOINT", "AGE", "SEX.MALE", "COMORBIDITIES.CI", "ANOVA_NULL", "ANOVA_COV")
			p_hd_extended_timepoint <- rbind(p_hd_extended_timepoint, p_all)
			
			## Get estimates 
			behaviour <- summary(xmdl)$coefficients[,1]
			behaviour <- behaviour[2:length(behaviour)]
			behaviour <- data.frame(behaviour)
			colnames(behaviour) <- c("Estimate")
			behaviour$Effect <- rownames(behaviour)
			behaviour$Module <- colnames(module)[1]
			behaviour$DISEASE <- "Sepsis.Timepoint"
			behaviour$Direction <- NA
			behaviour$Estimate <- as.numeric(as.character(behaviour$Estimate))
			behaviour$Direction[behaviour$Estimate > 0] <- "increase"
			behaviour$Direction[behaviour$Estimate < 0] <- "decrease"
			behaviour$Direction[behaviour$Estimate == 0] <- "no.change"
			### Reformat 
			behaviour$Effect[behaviour$Effect=="DISEASESEPSIS"] <- "DISEASE.SEPSIS"
			behaviour$Effect[behaviour$Effect=="SexMale"] <- "SEX.MALE"
			behaviour$Effect[behaviour$Effect=="DAY"] <- "TIMEPOINT"
			behaviour$Effect[behaviour$Effect=="Comorbidities.Charlson_Index"] <- "COMORBIDITIES.CI"
			behaviour$Effect[behaviour$Effect=="DAY:DISEASESEPSIS"] <- "TIMEPOINTxSEPSIS"
			behaviour$Effect[behaviour$Effect=="Age"] <- "AGE"
			behaviour_timepoint <- rbind(behaviour_timepoint, behaviour)	
			
			# DONE the timepoint assessment!!!
			#-----------------------------------------------------------------
			###########################################
			
			##########################
			## These do not need mixed effects models as there is no time component!!!!!
			## Day : Comparing HvsD (for all sepsis)
			## Should include age and sex as covariates
			module <- modules_interaction[i]
			module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE", "Age", "Sex", "Barcode","Comorbidities.Charlson_Index")]
			## Lets subset for health only 
			module <- module[module$DAY=="Day1",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			## Run model to see effect of time ## using it as a continuous variable? 
			model_formula <- formula(paste0(colnames(module)[1], "~", "DISEASE+Age + Sex  + Comorbidities.Charlson_Index"))#+Sex+Age
			model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index "))

			xmdl = lm(model_formula, module)
			xmdl_cov = lm(model_formula_covonly, module)
			xmdl.av <- aov(xmdl)
			
			## Get Overall P value for mortality from anova!!!!!!!!
			p_vals_day1 <- summary(xmdl)$coefficients[,4]
			p_vals_day1 <- p_vals_day1[2:5]
			
			## test difference in model to Covariate MODEL
			compare_mod_c <- anova(xmdl, xmdl_cov)
			p_overall_c <- compare_mod_c[,'Pr(>F)'][2]
			
			p_vals_day1 <- c(p_vals_day1, p_overall_c)

			## Get estimates 
			behaviour <- summary(xmdl)$coefficients[,1]
			behaviour <- behaviour[2:length(behaviour)]
			behaviour <- data.frame(behaviour)
			colnames(behaviour) <- c("Estimate")
			behaviour$Effect <- rownames(behaviour)
			behaviour$Module <- colnames(module)[1]
			behaviour$DAY <- "Day1"
			behaviour$Direction <- NA
			behaviour$Estimate <- as.numeric(as.character(behaviour$Estimate))
			behaviour$Direction[behaviour$Estimate > 0] <- "increase"
			behaviour$Direction[behaviour$Estimate < 0] <- "decrease"
			behaviour$Direction[behaviour$Estimate == 0] <- "no.change"
			### Reformat 
			behaviour$Effect[behaviour$Effect=="DISEASESEPSIS"] <- "DISEASE.SEPSIS"
			behaviour$Effect[behaviour$Effect=="SexMale"] <- "SEX.MALE"
			behaviour$Effect[behaviour$Effect=="DAY"] <- "TIMEPOINT"
			behaviour$Effect[behaviour$Effect=="Comorbidities.Charlson_Index"] <- "COMORBIDITIES.CI"
			behaviour$Effect[behaviour$Effect=="DAY:DISEASESEPSIS"] <- "TIMEPOINTxSEPSIS"
			behaviour$Effect[behaviour$Effect=="Age"] <- "AGE"
			behaviour_hvsd <- rbind(behaviour_hvsd, behaviour)

			###### Would be good to add a plot like we did for the mortality!!!
			### Perform post hoc test to see which parts are significant....
			postHocs <- glht(xmdl.av, linfct = mcp(DISEASE = "Tukey"))
			confin <- confint(postHocs)
			confid <- data.frame(confin$confint)
			p_val <-  data.frame(summary(postHocs)$test$pvalues)
			final <- cbind(confid, p_val)
			colnames(final) <- c("diff", "lwr", "upr", "p adj")
			rownames(final) <- gsub(" - ", "-", rownames(final))
			final[,4] <- as.numeric(final[,4] )
			tukey.test1 <- final
			tbl1 <- with(module, table(DISEASE))
			
			############################################
			###########################################
			
			## Day 3 HvsD
			module <- modules_interaction[i]
			module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE", "Age", "Sex", "Barcode","Comorbidities.Charlson_Index")]
			## Lets subset for health only 
			module <- module[module$DAY=="Day3",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			## Run model to see effect of time ## using it as a continuous variable? 
			model_formula <- formula(paste0(colnames(module)[1], "~", "DISEASE+Age + Sex  + Comorbidities.Charlson_Index"))#+Sex+Age
			model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index "))

			xmdl = lm(model_formula, module)
			xmdl_cov = lm(model_formula_covonly, module)
			xmdl.av <- aov(xmdl)
			
			## get p val just for day (p overall will be the same!)
			p_vals_day3 <- summary(xmdl)$coefficients[,4]
			p_vals_day3 <- p_vals_day3[2:5]
			
			## test difference in model to Covariate MODEL
			compare_mod_c <- anova(xmdl, xmdl_cov)
			p_overall_c <- compare_mod_c[,'Pr(>F)'][2]
			
			p_vals_day3 <- c(p_vals_day3, p_overall_c)
			
			## Get estimates 
			behaviour <- summary(xmdl)$coefficients[,1]
			behaviour <- behaviour[2:length(behaviour)]
			behaviour <- data.frame(behaviour)
			colnames(behaviour) <- c("Estimate")
			behaviour$Effect <- rownames(behaviour)
			behaviour$Module <- colnames(module)[1]
			behaviour$DAY <- "Day3"
			behaviour$Direction <- NA
			behaviour$Estimate <- as.numeric(as.character(behaviour$Estimate))
			behaviour$Direction[behaviour$Estimate > 0] <- "increase"
			behaviour$Direction[behaviour$Estimate < 0] <- "decrease"
			behaviour$Direction[behaviour$Estimate == 0] <- "no.change"
			### Reformat 
			behaviour$Effect[behaviour$Effect=="DISEASESEPSIS"] <- "DISEASE.SEPSIS"
			behaviour$Effect[behaviour$Effect=="SexMale"] <- "SEX.MALE"
			behaviour$Effect[behaviour$Effect=="DAY"] <- "TIMEPOINT"
			behaviour$Effect[behaviour$Effect=="Comorbidities.Charlson_Index"] <- "COMORBIDITIES.CI"
			behaviour$Effect[behaviour$Effect=="DAY:DISEASESEPSIS"] <- "TIMEPOINTxSEPSIS"
			behaviour$Effect[behaviour$Effect=="Age"] <- "AGE"
			behaviour_hvsd <- rbind(behaviour_hvsd, behaviour)

			################################
			### Perform post hoc test to see which parts are significant....
			postHocs <- glht(xmdl.av, linfct = mcp(DISEASE = "Tukey"))
			confin <- confint(postHocs)
			confid <- data.frame(confin$confint)
			p_val <-  data.frame(summary(postHocs)$test$pvalues)
			final <- cbind(confid, p_val)
			colnames(final) <- c("diff", "lwr", "upr", "p adj")
			rownames(final) <- gsub(" - ", "-", rownames(final))
			final[,4] <- as.numeric(final[,4] )
			tukey.test2 <- final
			tbl2 <- with(module, table(DISEASE))

			
			###########################################
			###########################################
			## Day 5 HvsD
			module <- modules_interaction[i]
			module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE", "Age", "Sex", "Barcode","Comorbidities.Charlson_Index")]
			## Lets subset for health only 
			module <- module[module$DAY=="Day5",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			## Run model to see effect of time ## using it as a continuous variable? 
			model_formula <- formula(paste0(colnames(module)[1], "~", "DISEASE+Age + Sex  + Comorbidities.Charlson_Index"))#+Sex+Age
			model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index "))

			xmdl = lm(model_formula, module)
			xmdl_cov = lm(model_formula_covonly, module)
			xmdl.av <- aov(xmdl)
			
			## get p val just for day (p overall will be the same!)
			p_vals_day5 <- summary(xmdl)$coefficients[,4]
			p_vals_day5 <- p_vals_day5[2:5]
			
			## test difference in model to Covariate MODEL
			compare_mod_c <- anova(xmdl, xmdl_cov)
			p_overall_c <- compare_mod_c[,'Pr(>F)'][2]
			
			p_vals_day5 <- c(p_vals_day5, p_overall_c)

			## Get estimates 
			behaviour <- summary(xmdl)$coefficients[,1]
			behaviour <- behaviour[2:length(behaviour)]
			behaviour <- data.frame(behaviour)
			colnames(behaviour) <- c("Estimate")
			behaviour$Effect <- rownames(behaviour)
			behaviour$Module <- colnames(module)[1]
			behaviour$DAY <- "Day5"
			behaviour$Direction <- NA
			behaviour$Estimate <- as.numeric(as.character(behaviour$Estimate))
			behaviour$Direction[behaviour$Estimate > 0] <- "increase"
			behaviour$Direction[behaviour$Estimate < 0] <- "decrease"
			behaviour$Direction[behaviour$Estimate == 0] <- "no.change"
			### Reformat 
			behaviour$Effect[behaviour$Effect=="DISEASESEPSIS"] <- "DISEASE.SEPSIS"
			behaviour$Effect[behaviour$Effect=="SexMale"] <- "SEX.MALE"
			behaviour$Effect[behaviour$Effect=="DAY"] <- "TIMEPOINT"
			behaviour$Effect[behaviour$Effect=="Comorbidities.Charlson_Index"] <- "COMORBIDITIES.CI"
			behaviour$Effect[behaviour$Effect=="DAY:DISEASESEPSIS"] <- "TIMEPOINTxSEPSIS"
			behaviour$Effect[behaviour$Effect=="Age"] <- "AGE"
			behaviour_hvsd <- rbind(behaviour_hvsd, behaviour)
			
			################################
			### Perform post hoc test to see which parts are significant....
			postHocs <- glht(xmdl.av, linfct = mcp(DISEASE = "Tukey"))
			confin <- confint(postHocs)
			confid <- data.frame(confin$confint)
			p_val <-  data.frame(summary(postHocs)$test$pvalues)
			final <- cbind(confid, p_val)
			colnames(final) <- c("diff", "lwr", "upr", "p adj")
			rownames(final) <- gsub(" - ", "-", rownames(final))
			final[,4] <- as.numeric(final[,4] )
			tukey.test3 <- final
			tbl3 <- with(module, table(DISEASE))
			
			##--------------------------------------------------------------------------------------
			##--------------------------------------------------------------------------------------
			##################################
			### Lets plot the differences in means via the tukey test 
			modulename <- modules_interaction[i]
			tukey.test1$Model <- modulename
			tukey.test1$Day <- "Day 1"
			tukey.test2$Model <- modulename
			tukey.test2$Day <- "Day 3"
			tukey.test3$Model <- modulename
			tukey.test3$Day <- "Day 5"
			tukey_all <- rbind(tukey.test1, tukey.test2, tukey.test3)
			tukey.max <- max(tukey_all$upr)+0.1
			tukey.min <- min(tukey_all$lwr)-0.1
			if(tukey.min>0){
				tukey.min <- -0.1
			}
			
			plot_dir_tukey1 <- paste0(outputdir, "/TukeyTest/")
			if (!dir.exists(plot_dir_tukey1)) {dir.create(plot_dir_tukey1)}
			plot_dir_tukey2 <- paste0(outputdir, "/TukeyTest/HealthvsSepsis")
			if (!dir.exists(plot_dir_tukey2)) {dir.create(plot_dir_tukey2)}
			
			pdf(paste0(plot_dir_tukey2,"/HvsD_TukeyTest", colnames(eigenvectors)[i], ".pdf"), width=8, height=7)  
			par( mfrow= c(2,3) )
			s <- plotTukeyHSD(tukey.test1, modulename, "Day 1", tukey.min, tukey.max )
			plotTukeyHSD(tukey.test2, modulename,"Day 3", tukey.min, tukey.max )
			plotTukeyHSD(tukey.test3, modulename, "Day 5", tukey.min, tukey.max )
			barplot(tbl1, beside = TRUE, ylim=range(pretty(c(0, tbl1, tbl2, tbl3))), main=paste0(modulename, "\nDay 1"), ylab="Sample Count", xlab="Disease State")
			barplot(tbl2, beside = TRUE, ylim=range(pretty(c(0, tbl1, tbl2, tbl3))), main=paste0(modulename, "\nDay 3"), ylab="Sample Count", xlab="Disease State")
			barplot(tbl3, beside = TRUE, ylim=range(pretty(c(0, tbl1, tbl2, tbl3))), main=paste0(modulename, "\nDay 5"), ylab="Sample Count", xlab="Disease State")
			dev.off()
			write.table(tukey_all, paste0(plot_dir_tukey2, "/HvsD_TukeyTest_Statistics_",modulename, ".txt"), sep="\t")
			
			##--------------------------------------------------------------------------------------
			####################################
			## Lets take all the p values and append them 
			p_vals_day1 <- c(modulename, "Day1", p_vals_day1)
			p_vals_day3 <- c(modulename, "Day3", p_vals_day3)
			p_vals_day5 <- c(modulename, "Day5", p_vals_day5)
			
			row_use <- data.frame(rbind(p_vals_day1, p_vals_day3, p_vals_day5))
			colnames(row_use) <- c("MODULE", "DAY", "DISEASE.SEPSIS", "AGE", "SEX.MALE", "COMORBIDITIES.CI", "ANOVA_COV")
			row_use[, 3:6] <- lapply(3:6, function(x) as.numeric(row_use[[x]]))
			p_hd_extended_hvsd <- rbind(row_use, p_hd_extended_hvsd)
			
			######################
			print(paste0("Done module ", modulename))
		}
		
		####------------------------------------------------------------------------------------------
		####------------------------------------------------------------------------------------------
		## NOW we move on to plotting!!!!
		## Lets do for the timepoint analysis first:
		### Now we move onto plotting 
		p_hd_extended_timepoint <- data.frame(p_hd_extended_timepoint)
		p_hd_extended_timepoint <- gather(p_hd_extended_timepoint, "Effect", "p_value", TIMEPOINT:ANOVA_COV, factor_key=TRUE)
		p_hd_extended_timepoint$p_value <- as.numeric(p_hd_extended_timepoint$p_value)
		p_hd_extended_timepoint$MODULE <- factor(p_hd_extended_timepoint$MODULE, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25","Module_26","Module_27", "Module_28","Module_29","Module_30", "Module_31"))		
		## Now we want to merge with behaviour 
		colnames(behaviour_timepoint) <- c("Estimate", "Effect", "MODULE", "DISEASE", "Direction")
		p_hd_extended_timepoint <- merge(p_hd_extended_timepoint, behaviour_timepoint, by=c("MODULE", "DISEASE", "Effect"), all.x=TRUE)
		#########
		p_use <- 0.05
		p_hd_extended_timepoint$p_value <- as.numeric(p_hd_extended_timepoint$p_value)
		## Now lets do a BH correction
		p_hd_extended_timepoint$BH <- p.adjust(as.numeric(p_hd_extended_timepoint$p_value),method="BH")
		p_hd_extended_timepoint$psig <- "ns"
		p_hd_extended_timepoint$psig[p_hd_extended_timepoint$p_value < 0.05] <- "sig"
		p_hd_extended_timepoint$BHsig <- "ns"
		p_hd_extended_timepoint$BHsig[p_hd_extended_timepoint$BH < 0.05] <- "sig"
		p_hd_extended_timepoint$Type <- "ADDITIVE MODEL"
		
		p_hd_extended_timepoint$Type[p_hd_extended_timepoint$Effect %like% "ANOVA"] <- "ANOVA"
		p_hd_extended_timepoint$Effect <- as.character(p_hd_extended_timepoint$Effect)
		p_hd_extended_timepoint$Effect[p_hd_extended_timepoint$Effect =="ANOVA_NULL"] <- "NULL MODEL"
		p_hd_extended_timepoint$Effect[p_hd_extended_timepoint$Effect =="ANOVA_COV"] <- "COV ONLY MODEL"
		p_hd_extended_timepoint$Effect <- factor(p_hd_extended_timepoint$Effect, levels=c("TIMEPOINT", "AGE", "SEX.MALE", "COMORBIDITIES.CI", "NULL MODEL", "COV ONLY MODEL"))
		write.table(p_hd_extended_timepoint, paste0(outputdir, "/Summary/HvsD_INTERACTION_ANALYSIS_timepoint_", annot, ".txt"), sep="\t")
		## Lets remove the direction for ns effect
		p_hd_extended_timepoint$Direction[p_hd_extended_timepoint$psig=="ns"] <- NA
		p_hd_extended_timepoint$Direction <- factor(p_hd_extended_timepoint$Direction, levels=c("increase", "decrease", NA))
		#######################################################################################
		
		## NOW we repeat for the HvsD analysis 
		p_hd_extended_hvsd <- data.frame(p_hd_extended_hvsd)
		p_hd_extended_hvsd <- gather(p_hd_extended_hvsd, "Effect", "p_value", DISEASE.SEPSIS:ANOVA_COV, factor_key=TRUE)
		p_hd_extended_hvsd$p_value <- as.numeric(p_hd_extended_hvsd$p_value)
		p_hd_extended_hvsd$MODULE <- factor(p_hd_extended_hvsd$MODULE, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25","Module_26","Module_27", "Module_28","Module_29","Module_30", "Module_31"))
		## Now we want to merge with behaviour 
		colnames(behaviour_hvsd) <- c("Estimate", "Effect", "MODULE", "DAY", "Direction")
		p_hd_extended_hvsd <- merge(p_hd_extended_hvsd, behaviour_hvsd, by=c("MODULE", "DAY", "Effect"), all.x=TRUE)
		#########
		p_use <- 0.05
		p_hd_extended_hvsd$p_value <- as.numeric(p_hd_extended_hvsd$p_value)
		## Now lets do a BH correction
		p_hd_extended_hvsd$BH <- p.adjust(as.numeric(p_hd_extended_hvsd$p_value),method="BH")
		p_hd_extended_hvsd$psig <- "ns"
		p_hd_extended_hvsd$psig[p_hd_extended_hvsd$p_value < 0.05] <- "sig"
		p_hd_extended_hvsd$BHsig <- "ns"
		p_hd_extended_hvsd$BHsig[p_hd_extended_hvsd$BH < 0.05] <- "sig"
		p_hd_extended_hvsd$Type <- "ADDITIVE MODEL"
		
		p_hd_extended_hvsd$Type[p_hd_extended_hvsd$Effect%like% "ANOVA"] <- "ANOVA"
		
		p_hd_extended_hvsd$Effect <- as.character(p_hd_extended_hvsd$Effect)
		p_hd_extended_hvsd$Effect[p_hd_extended_hvsd$Effect =="ANOVA_NULL"] <- "NULL MODEL"
		p_hd_extended_hvsd$Effect[p_hd_extended_hvsd$Effect =="ANOVA_COV"] <- "COV ONLY MODEL"
		p_hd_extended_hvsd$Effect <- factor(p_hd_extended_hvsd$Effect, levels=c("DISEASE.SEPSIS", "AGE", "SEX.MALE", "COMORBIDITIES.CI", "COV ONLY MODEL"))
		write.table(p_hd_extended_hvsd, paste0(outputdir, "/Summary/HvsD_INTERACTION_ANALYSIS_hvsd_", annot, ".txt"), sep="\t")
		## Lets remove the direction for ns effect
		p_hd_extended_hvsd$Direction[p_hd_extended_hvsd$psig=="ns"] <- NA
		p_hd_extended_hvsd$Direction <- factor(p_hd_extended_hvsd$Direction, levels=c("increase", "decrease", NA))
		#######################################################################################
		
		##################################################################################
		#----------------------------------------------------	
		### Now format for plotting 
		## Remove the na values from timepoint (we didnt do covariate)
		p_hd_extended_timepoint <- p_hd_extended_timepoint[!is.na(p_hd_extended_timepoint$p_value),]
		p_hd_extended_hvsd$DAY <- gsub("Day", "Day ", p_hd_extended_hvsd$DAY)
		p_hd_extended_timepoint$BHsig <-  factor(p_hd_extended_timepoint$BHsig, levels=c("sig", "ns"))
		p_hd_extended_hvsd$BHsig <- factor(p_hd_extended_hvsd$BHsig, levels=c("sig", "ns"))
		p_hd_extended_timepoint$psig <-  factor(p_hd_extended_timepoint$psig, levels=c("sig", "ns"))
		p_hd_extended_hvsd$psig <- factor(p_hd_extended_hvsd$psig, levels=c("sig", "ns"))
		##------------------------
		p_hd_extended_timepoint$DISEASE <- gsub(".Timepoint", "", p_hd_extended_timepoint$DISEASE)
		
		p_hd_extended_timepoint$Type <- str_wrap(p_hd_extended_timepoint$Type, width = 12)
		p_hd_extended_hvsd$Type <- str_wrap(p_hd_extended_hvsd$Type, width = 12)
		
		
		## Lets plot to see P value distribution
		pdf(paste0(outputdir,"/HvsD_SimpleEffects_PDistributions_", annot, ".pdf"), width=8, height=4)
		px <- ggplot(p_hd_extended_timepoint) +  geom_point(aes(x=p_value, y=BH, alpha=psig, colour = BHsig)) +xlab("Model Effect P ") + ylab("BH Correced Model Effect P") +theme_classic() +labs(alpha="P < 0.05", colour="BH Adj P <0.05")+xlab("P") +ylab("BH Adj P") + scale_color_manual(values = c( 'red', "lightblue"), drop=FALSE) +geom_vline(xintercept=0.05, col="blue") + geom_hline(yintercept=0.05, col="blue")+ggtitle(paste0("P values for LR Timepoint Models ", type_receptor))+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		px1 <- ggplot(p_hd_extended_hvsd) +  geom_point(aes(x=p_value, y=BH, alpha=psig, colour = BHsig)) +xlab("Model Effect P ") + ylab("BH Correced Model Effect P") +theme_classic() +labs(alpha="P < 0.05", colour="BH P Adj <0.05")+xlab("P") +ylab("BH Adj P") + scale_color_manual(values = c( 'red', "lightblue"), drop=FALSE) +geom_vline(xintercept=0.05, col="blue") + geom_hline(yintercept=0.05, col="blue")+ ggtitle(paste0("P values for LR HvsS Model ", type_receptor))+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		plot_grid(px, px1, ncol=2, align="h", axis="lbt", rel_widths=c(1,1), labels="AUTO")
		dev.off()
		
		### Lets do the overall plots 
		pdf(paste0(outputdir, "/HvsD_SimpleEffects_Statistics_", annot, ".pdf"), width=11, height=12)
		p1 <- ggplot(p_hd_extended_hvsd, aes(x=MODULE, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("B. Simple Effects Tests: Disease Status Model")) + labs(fill="P") +xlab("Module") +ylab("Model Coefficients")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P < 0.05", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type), cols=vars(DAY),  scales="free_y", space='free_y')+ geom_point(aes(x=MODULE, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction\nof Effect\nRelative to:\nA/B. Health\nC. Day1", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
			
		p2 <- ggplot(p_hd_extended_timepoint, aes(x=MODULE, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("C. Simple Effects Tests: Time Model")) + labs(fill="P") +xlab("Module") +ylab("Model Coefficients")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P < 0.05", values = c("red", '#00000000'), drop = FALSE) + facet_grid(rows=vars(Type), cols=vars(DISEASE),  scales="free_y", space='free_y')+ geom_point(aes(x=MODULE, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction\nof Effect\nRelative to:\nA/B. Health\nC. Day1", alpha="P <0.05") + scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		
		plot(plot_grid(p1, p2, labels="AUTO", align="v", axis="lbt", ncol=1 ))
		
		####
		p1 <- ggplot(p_hd_extended_hvsd, aes(x=MODULE, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, p_use))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("B. Simple Effects Tests: Disease Status Model")) + labs(fill="P") +xlab("Module") +ylab("Model Coefficients")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P < 0.05", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type), cols=vars(DAY), scales="free_y", space='free_y')+ geom_point(aes(x=MODULE, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction\nof Effect\nRelative to:\nA/B. Health\nC. Day1", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		
		p2 <- ggplot(p_hd_extended_timepoint, aes(x=MODULE, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, p_use))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("C. Simple Effects Tests: Time Model")) + labs(fill="P") +xlab("Module") +ylab("Model Coefficients")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P < 0.05", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type), cols=vars(DISEASE), scales="free_y", space='free_y')+ geom_point(aes(x=MODULE, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction\nof Effect\nRelative to:\nA/B. Health\nC. Day1", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		
		#plot(plot_grid(p1, p2, labels="AUTO", align="v", axis="lbt", ncol=1  ))
		grid_arrange_shared_legend(p2x, p1, p2, ncol = 1, nrow=3, position="right")

		###
		p1 <- ggplot(p_hd_extended_hvsd, aes(x=MODULE, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, 0.1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("B. Simple Effects Tests: Disease Status Model")) + labs(fill="P") +xlab("Module") +ylab("Model Coefficients")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P < 0.05", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type), cols=vars(DAY),  scales="free_y", space='free_y')+ geom_point(aes(x=MODULE, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction\nof Effect\nRelative to:\nA/B. Health\nC. Day1", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		
		p2 <- ggplot(p_hd_extended_timepoint, aes(x=MODULE, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, 0.1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("C. Simple Effects Tests: Time Model")) + labs(fill="P") +xlab("Module") +ylab("Model Coefficients")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P < 0.05", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type), cols=vars(DISEASE),  scales="free_y", space='free_y')+ geom_point(aes(x=MODULE, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction\nof Effect\nRelative to:\nA/B. Health\nC. Day1", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		grid_arrange_shared_legend(p3x, p1, p2, ncol = 1, nrow=3, position="right")		
		#plot(plot_grid(p1, p2, labels="AUTO", align="v", axis="lbt" , ncol=1 ))		
		dev.off()
	
		## Now we just need to return from the function!!!
		
		return(p_hd_extended_hvsd)
	} else {
	 return(p_values_hd)
	 }

}

## DONE