##-------------------------------------------------------
## Code to model gene expression counts and across time grouped by MORTALITY
## Interaction model followed by simple effects tests 
## lauren Overend 
## lauren.overend@live.co.uk
## Feb 2023 
##-------------------------------------------------------

## Dependencies 
library(doParralel)
library(foreach)
library(ggrepel)
library(gridExtra)
library(grid)
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
###############
#source('/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/PlotTukey.R')
#cyber_sort <- "/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/GeneExpression/Isotype_lcpm.txt" 
#eigenvectors <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt'
#outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/LEO_GAinS_RNASEQ_2023/Comparisons/'
#metadata <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Final_metadata_Reduced.txt'
#metahealth <- '/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/LabKeyMetaData/Healthies_ClinData.txt'

correlate_eigenvectors_outcome <- function(eigenvectors, cyber_sort, metadata, outputdir, type_receptor, run_all){
	
	###########################
	eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
	### We need to exclude the two bad samples 
	### These samples had ID mix ups and need to be removed!!!
	print("Removing Samples where RNAseq suggests a sample mixup")
	bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
	eigenvectors <- eigenvectors[!eigenvectors$sample %in% bad_ids,]
	
	######## Read in expression matrix
	## Subset to only include those in eigenvector (repertoire cohort!!)
	cyber <- read.delim(cyber_sort) 
	cyber <- cyber[cyber$sample %in% eigenvectors$sample,]	
	## This is what we will be correlating 
	cells <- colnames(cyber)[colnames(cyber)%like% c("IGH")]
		
		
	##------------------------------------------------------------------------------
	### PART 1
	## Preparing metadata
	covariates <- c("Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index")
	metadata <- read.delim(metadata, sep="\t", header=TRUE)
	## Include for just samples that are in the dataset (after filtering etc)
	metadata <- metadata[metadata$SampleID_alternative %in% rownames(eigenvectors),]

	comobs <- covariates
	mort <- colnames(metadata)[colnames(metadata)%like% "Mortality"]
	outcome <- c(mort, "Infection.Number_of_ICU_Acquired_Infections", "Death.Days_from_ICU_admission")
	metadata <- metadata[, c( "SampleID_alternative", "alternative_barcode", outcome, comobs)]
	metadata$Mortality2 <- metadata$Mortality.Classification
	metadata$Mortality2[metadata$Mortality2=="MID.DEATH"] <- "LATE.DEATH"
	## Very important that we set the factor level
	metadata$Mortality2 <- factor(metadata$Mortality2, levels=c("SURVIVOR", "LATE.DEATH", "EARLY.DEATH"))
	metadata$Sex <- factor(metadata$Sex, levels=c("Female", "Male"))
	## Add barcode!
	metadata$Barcode <- NA
	for(x in 1:length(metadata$SampleID_alternative)){
			metadata$Barcode[x] <- str_split_fixed (metadata$SampleID_alternative[x], "_", 2)[,1]
			if(metadata$SampleID_alternative[x] %like% "HV"){
			 metadata$Barcode[x] <- paste0(str_split_fixed(metadata$SampleID_alternative[x], "_", 3)[,1], "_", str_split_fixed(metadata$SampleID_alternative[x], "_", 3)[,2])
			}
	}
	
	## Merge with gene expression
	modeling_cells <- merge(metadata, cyber, by.x="SampleID_alternative", by.y="sample")
	modeling_cells$DAY <- NA
	modeling_cells$DAY[grep("_1", modeling_cells$SampleID_alternative)] <- "1"
	modeling_cells$DAY[grep("_3", modeling_cells$SampleID_alternative)] <- "3"
	modeling_cells$DAY[grep("_5", modeling_cells$SampleID_alternative)] <- "5"
	### Data is ready to model!!!!!
	
	p_values_out <- c()
	myplots <- list()
	behaviour <- c()
	vif_all <- c()
	## Run a simple model to look for interaction
	for(i in 1:(length(cells))){
		module <- cells[i]
		module <- modeling_cells[,c(module, "SampleID_alternative", "DAY", "Mortality2", "Age", "Sex", "Comorbidities.Charlson_Age","Comorbidities.Charlson_Index", "Barcode")]
		colnames(module) <- c(cells[i], "SampleID_alternative", "DAY", "Mortality", "Age", "Sex", "Comorbidities.Charlson_Age","Comorbidities.Charlson_Index", "Barcode")
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		
		## Model with Interaction  we need to run an mancova as we have multiple categorical levels within mortality! 
		model_formula <- formula(paste0(colnames(module)[1], "~", "DAY*Mortality + Age + Sex  + Comorbidities.Charlson_Index+(1|Barcode)"))
		model_formula_basic <- formula(paste0(colnames(module)[1], "~"," 1 +(1|Barcode)"))
		model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))

		xmdl =  lmerTest::lmer(model_formula, module, REML=F)
		xmdl_old <- lmerTest::lmer(model_formula_basic, module, REML=F)
		xmdl_cov <- lmerTest::lmer(model_formula_covonly, module, REML=F)

		## VIF
		vif_scores <- vif(xmdl)[,1]
		vif_scores <- c(cells[i], vif_scores)	
		vif_all <- rbind(vif_all, vif_scores)
		############################
		
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
		
		
		######################################

		row_data <- c(cells[i], p_vals, p_overall, p_overall_c)
		row_data <- data.frame(t(row_data))
		row_data[,2:dim(row_data)[2]] <-as.numeric(row_data[,2:dim(row_data)[2]])
		colnames(row_data) <- c("Module", paste0("IM.", names(p_vals)), "ANOVA_NULL", "ANOVA_COV") 

		## Can we do an interaction plot of residuals
		p1 <- interactions::interact_plot(xmdl, pred = DAY, modx = Mortality, partial.residuals = TRUE, interval = TRUE )
		myplots[[(i)]] <- p1
		p_values_out <- rbind(row_data, p_values_out)
		
		## create a dataframe storing effect information
		estimates <- summary(xmdl)$coefficients[,1]
		estimates <- estimates[2:length(estimates)]
		estimates <- data.frame(estimates)
		colnames(estimates) <- c("Estimate")
		estimates$Effect <- rownames(estimates)
		estimates$Module <- cells[i]
		estimates$Direction <- NA
		estimates$Estimate <- as.numeric(as.character(estimates$Estimate))
		estimates$Direction[estimates$Estimate > 0] <- "increase"
		estimates$Direction[estimates$Estimate < 0] <- "decrease"
		estimates$Direction[estimates$Estimate == 0] <- NA
	
		estimates$Effect[estimates$Effect=="DISEASESEPSIS"] <- "DISEASE.SEPSIS"
		estimates$Effect[estimates$Effect=="SexMale"] <- "SEX.MALE"
		estimates$Effect[estimates$Effect=="DAY"] <- "TIMEPOINT"
		estimates$Effect[estimates$Effect=="Comorbidities.Charlson_Index"] <- "COMORBIDITIES:CI"
		estimates$Effect[estimates$Effect=="DAY:DISEASESEPSIS"] <- "TIMEPOINTxSEPSIS"
		estimates$Effect[estimates$Effect=="Age"] <- "AGE"
		estimates$Effect[estimates$Effect=="MortalityLATE.DEATH"] <- "MORTALITY.LATEDEATH"
		estimates$Effect[estimates$Effect=="MortalityEARLY.DEATH"] <- "MORTALITY.EARLYDEATH"
		estimates$Effect[estimates$Effect=="DAY:MortalityLATE.DEATH"] <- "TIMEPOINTxMORTALITY.LATEDEATH"
		estimates$Effect[estimates$Effect=="DAY:MortalityEARLY.DEATH"] <- "TIMEPOINTxMORTALITY.EARLYDEATH"
		## we dont want diction of interaction have to specifiy here as charlson index has an x in!!!!!!
		estimates$Direction[estimates$Effect %like% "TIMEPOINTx"] <- NA
		behaviour <- rbind(behaviour, estimates)	
		}
		
	colnames(p_values_out) <- c("Module", "TIMEPOINT", "MORTALITY.LATEDEATH", "MORTALITY.EARLYDEATH", "AGE", "SEX.MALE", "COMORBIDITIES:CI", "TIMEPOINTxMORTALITY.LATEDEATH", "TIMEPOINTxMORTALITY.EARLYDEATH",  "ANOVA_NULL", "ANOVA_COV") 	
	
	#--------------------------------------------------------
	### Lets quickly plot the model VIF scores (is colinearity a problem)	
	pdf(paste0(outputdir,"/Mortality_Immunoglobulin_VIF_Score.pdf"), width=5, height=5)
	vif_all <- data.frame(vif_all)
	vif_all_g <- gather(vif_all, "Effect", "VIF", DAY:DAY.Mortality, factor_key=TRUE)
	vif_all_g <- vif_all_g[vif_all_g$V1 =="IGHM",]
	vif_all_g$VIF <- as.numeric(vif_all_g$VIF)
	vif_all_g$Effect <- as.character(vif_all_g$Effect)
	vif_all_g$Effect[vif_all_g$Effect=="Comorbidities.Charlson_Index"] <- "COMORBIDITIES:CI"
	vif_all_g$Effect[vif_all_g$Effect=="DAY.DISEASE"] <- "DAY*DISEASE"
	vif_all_g$Effect[vif_all_g$Effect=="DAY.Mortality"] <- "DAY*MORTALITY"
	vif_all_g$Effect[vif_all_g$Effect=="Mortality"] <- "MORTALITY"
	vif_all_g$Effect[vif_all_g$Effect=="Age"] <- "AGE"
	vif_all_g$Effect[vif_all_g$Effect=="Sex"] <- "SEX"
	vif_all_g$Effect <- factor(vif_all_g$Effect, levels=c("DAY", "MORTALITY", "DAY*MORTALITY", "AGE", "SEX", "COMORBIDITIES:CI"))
	plot(ggplot(vif_all_g, aes(x=Effect, y=VIF))+geom_col(alpha=0.5) +theme_classic()+ylab("Module Score")+xlab("Day")+ggtitle(paste0("Interaction Model VIFF Scores")) +facet_wrap(~V1)+geom_hline(yintercept=5, col="red") + xlab("Model Effect") + ylab("VIFF")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
	dev.off()
	#--------------------------------------------------------
	## Do residual plots to visualise corrected values after adjusting for covariates. 
	pdf(paste0(outputdir, "/Mortality_Immunoglobulin_AdjModel_Residual_Plot.pdf"),width=13, height=7)
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/multiplot.R')
	multiplot(plotlist = myplots, cols = 3)
	dev.off()
	
	#--------------------------------------------------------
	## GATHER_TEST_RESULTS_for doing the interaction plot
	p_hd <- gather(p_values_out, "Effect", "p_value", TIMEPOINT:ANOVA_COV, factor_key=TRUE)
	## Lets PLOT The statistics for them  
	p_use <- 0.05
	p_hd$p_value <- as.numeric(p_hd$p_value)
	
	###############################################
	## Now lets do a BH correction
	p_hd$BH <- p.adjust(as.numeric(p_hd$p_value),method="BH")
	p_hd$psig <- "ns"
	p_hd$psig[p_hd$p_value < 0.05] <- "sig"
	p_hd$BHsig <- "ns"
	p_hd$BHsig[p_hd$BH < 0.05] <- "sig"
	p_hd_sub <- p_hd[p_hd$BHsig =="sig",]
	p_hd$Type <- "INTERACTION MODEL EFFECT"
	p_hd$Type[p_hd$Effect %like% "ANOVA"] <- "ANOVA"
	p_hd$Effect <- as.character(p_hd$Effect)
    p_hd$Effect[p_hd$Effect =="ANOVA_NULL"] <- "NULL MODEL"
	p_hd$Effect[p_hd$Effect =="ANOVA_COV"] <- "COV ONLY MODEL"
	p_hd$Effect <- factor(p_hd$Effect, levels=c("TIMEPOINT", "MORTALITY.LATEDEATH", "MORTALITY.EARLYDEATH", "AGE", "SEX.MALE", "COMORBIDITIES:CI", "TIMEPOINTxMORTALITY.LATEDEATH", "TIMEPOINTxMORTALITY.EARLYDEATH", "NULL MODEL", "COV ONLY MODEL"))
	
	#### SAVE INTERACTION STATS
	write.table(p_hd, paste0(outputdir, "/Mortality_Immunoglobulin_INTERACTION_ANALYSIS_SUMMARY.txt"), sep="\t")

	## Lets plot 
	px <- ggplot(p_hd) +  geom_point(aes(x=p_value, y=BH, alpha=psig, colour = BHsig)) +xlab("Model Effect P ") + ylab("BH Correced Model Effect P") +theme_classic() +labs(alpha="P < 0.05", colour="BH P <0.05")+xlab("P") +ylab("BH adj P") + scale_color_manual(values = c('lightblue', 'red')) +geom_vline(xintercept=0.05, col="blue") + geom_hline(yintercept=0.05, col="blue")
	
	## Lets add in the effect direction so that we can add an arrow on the plot indicating direction
	p_hd <- merge(p_hd, behaviour, by=c("Module", "Effect"), all.x=TRUE)
	## Lets remove the direction for ns effect
	p_hd$Direction[p_hd$psig=="ns"] <- NA
	p_hd$Direction <- factor(p_hd$Direction, levels=c("increase", "decrease", NA))
	p_hd$psig <- factor(p_hd$psig, levels=c("sig", "ns", NA))
	p_hd$BHsig <- factor(p_hd$BHsig, levels=c("sig", "ns", NA))

	## Do a plot showing BHP vs P
	px <- ggplot(p_hd) +  geom_point(aes(x=p_value, y=BH, alpha=psig, colour = BHsig)) +theme_classic() +labs(alpha="P<0.05", colour="BH P<0.05")+xlab("P") +ylab("Benjamin Hochburg Adjusted P") + scale_color_manual(values = c( 'red','lightblue')) +geom_vline(xintercept=0.05, col="blue") + geom_hline(yintercept=0.05, col="blue")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE) + ggtitle("Multiple Testing \n Correction")
	
	pdf(paste0(outputdir, "/Mortality_Immunoglobulin_InteractionModel_Statistics.pdf"), width=12, height=6)
	p1 <- ggplot(p_hd, aes(x=Module, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("Interaction Model, P<0.05")) + labs(fill="P") +xlab("Gene") +ylab("Model Coefficients")+geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type),  scales="free_y", space='free_y')+ geom_point(aes(x=Module, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect", alpha="P<0.05") + scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
	plot(plot_grid(p1, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))
	
	p2 <- ggplot(p_hd, aes(x=Module, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, p_use))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("Interaction Model, P<0.05")) + labs(fill="P")+xlab("Gene") +ylab("Model Coefficients")+geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type),  scales="free_y", space='free_y')+ geom_point(aes(x=Module, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect", alpha="P<0.05") + scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
	plot(plot_grid(p2, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))	
	
	p3 <- ggplot(p_hd, aes(x=Module, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, 0.1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("Interaction Model, P<0.05")) + labs(fill="P") +xlab("Gene") +ylab("Model Coefficients")+geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type),  scales="free_y", space='free_y')+ geom_point(aes(x=Module, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect", alpha="P<0.05") + scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
	plot(plot_grid(p3, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))		
	dev.off()
	
	##------------------------------------------------------------------------------------------------------
	## Does health vs disease significantly interact for any models??
	## If so we will need to investigate these modules at individual levels. 
	if(any(p_hd$p_value[p_hd$Effect=="TIMEPOINTxMORTALITY.LATEDEATH"] < 0.05  |p_hd$p_value[p_hd$Effect=="TIMEPOINTxMORTALITY.EARLYDEATH"] < 0.05)){
		print("Presence of significant interaction between timepoint and mortality")
		interaction_pres <- "YES"
		modules_interaction <- p_values_out$Module[p_values_out$TIMEPOINTxMORTALITY.LATEDEATH < 0.05 | p_values_out$TIMEPOINTxMORTALITY.EARLYDEATH< 0.05]
		annot <- "SignificantIntOnly"
	} else {
		print("No interaction effect between timepoint and mortality status")
		interaction_pres <- "NO"
		
	}
	
	
	##.........................................................................................................
	##.........................................................................................................
	##.........................................................................................................
	####################################################################################################
	### Part 2 more complex model IF an interaction is present 
	### Because of the interaction we need to assess at levels of Time / Disease state rather than across
	#-----------------------------------------------------------------------------------
	
	if(interaction_pres =="YES"){
		p_hd_extended <- c()
		behaviour_x <- c()
		p_hd_extended_timepoint <- c()
		behaviour_timepoint <- c()
		
		for(i in 1:length(modules_interaction)){
			#-----------------------------------------------------------------
			# EFFECT of TIME 
			####################################################
			module <- modules_interaction[i]
			module <- modeling_cells[,c(module, "SampleID_alternative", "DAY", "Mortality2", "Age", "Sex", "Comorbidities.Charlson_Age","Comorbidities.Charlson_Index", "Barcode")]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			modulex <- module 
			
			#################################################################
			## Survivors over time 
			modulex<- modulex
			module <- module[module$Mortality2=="SURVIVOR",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			
			## Run model to see effect of time ## using it as a continuous variable? 
			model_formula <- formula(paste0(colnames(module)[1], "~", "DAY+ Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))
			model_formula_basic <- formula(paste0(colnames(module)[1], "~"," 1 +(1|Barcode)"))
			model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))

			xmdl =  lmerTest::lmer(model_formula, module, REML=F)
			xmdl_old <- lmerTest::lmer(model_formula_basic, module, REML=F)
			xmdl_cov <- lmerTest::lmer(model_formula_covonly, module, REML=F)
	
			## Take out the p values 
			p_vals_s <- summary(xmdl)$coefficients[,5]
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
		
			## APPEND P VALUES AND ANOVA
			p_vals_s <- c(p_vals_s, p_overall, p_overall_c)
			
			## Get estimates 
			behaviour <- summary(xmdl)$coefficients[,1]
			behaviour <- behaviour[2:length(behaviour)]
			behaviour <- data.frame(behaviour)
			colnames(behaviour) <- c("Estimate")
			behaviour$Effect <- rownames(behaviour)
			behaviour$Module <- colnames(module)[1]
			behaviour$DISEASE <- "SURVIVORS.TIMEPOINT"
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
			behaviour$Effect[behaviour$Effect=="Mortality2LATE.DEATH"] <- "MORTALITY.LATEDEATH"
			behaviour$Effect[behaviour$Effect=="Mortality2EARLY.DEATH"] <- "MORTALITY.EARLYDEATH"
			behaviour$Effect[behaviour$Effect=="DAY:Mortality2LATE.DEATH"] <- "TIMEPOINTxMORTALITY.LATEDEATH"
			behaviour$Effect[behaviour$Effect=="DAY:Mortality2EARLY.DEATH"] <- "TIMEPOINTxMORTALITY.EARLYDEATH"
			
			behaviour_timepoint <- rbind(behaviour_timepoint, behaviour)
			
			
			##--------------------------------------------------------------
			#########################
			## Non survivors  Early Death?  
			module <- modulex
			module <- module[module$Mortality2=="EARLY.DEATH",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			## Run model to see effect of time ## using it as a continuous variable? 
			## Run model to see effect of time ## using it as a continuous variable? 
			model_formula <- formula(paste0(colnames(module)[1], "~", "DAY+ Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))
			model_formula_basic <- formula(paste0(colnames(module)[1], "~"," 1 +(1|Barcode)"))
			model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))

			xmdl =  lmerTest::lmer(model_formula, module, REML=F)
			xmdl_old <- lmerTest::lmer(model_formula_basic, module, REML=F)
			xmdl_cov <- lmerTest::lmer(model_formula_covonly, module, REML=F)
	
			## Take out the p values 
			p_vals_ne <- summary(xmdl)$coefficients[,5]
			p_vals_ne <- p_vals_ne[2:5]
			
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
		
			## APPEND P VALUES AND ANOVA
			p_vals_ne <- c(p_vals_ne, p_overall, p_overall_c)
			
			## Get estimates 
			behaviour <- summary(xmdl)$coefficients[,1]
			behaviour <- behaviour[2:length(behaviour)]
			behaviour <- data.frame(behaviour)
			colnames(behaviour) <- c("Estimate")
			behaviour$Effect <- rownames(behaviour)
			behaviour$Module <- colnames(module)[1]
			behaviour$DISEASE <- "EARLY.DEATH.TIMEPOINT"
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
			behaviour$Effect[behaviour$Effect=="Mortality2LATE.DEATH"] <- "MORTALITY.LATEDEATH"
			behaviour$Effect[behaviour$Effect=="Mortality2EARLY.DEATH"] <- "MORTALITY.EARLYDEATH"
			behaviour$Effect[behaviour$Effect=="DAY:Mortality2LATE.DEATH"] <- "TIMEPOINTxMORTALITY.LATEDEATH"
			behaviour$Effect[behaviour$Effect=="DAY:Mortality2EARLY.DEATH"] <- "TIMEPOINTxMORTALITY.EARLYDEATH"
			behaviour_timepoint <- rbind(behaviour_timepoint, behaviour)
			
			#####################################################################
			## Non survivors  Mid/Late Death?  
			module <- modulex
			module <- module[module$Mortality2=="LATE.DEATH",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			## Run model to see effect of time ## using it as a continuous variable? 
			model_formula <- formula(paste0(colnames(module)[1], "~", "DAY + Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))
			model_formula_basic <- formula(paste0(colnames(module)[1], "~"," 1 +(1|Barcode)"))
			model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))

			xmdl =  lmerTest::lmer(model_formula, module, REML=F)
			xmdl_old <- lmerTest::lmer(model_formula_basic, module, REML=F)
			xmdl_cov <- lmerTest::lmer(model_formula_covonly, module, REML=F)
	
			p_vals_nl <- summary(xmdl)$coefficients[,5]
			p_vals_nl <- p_vals_nl[2:5]
			
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
		
			## APPEND P VALUES AND ANOVA
			p_vals_nl <- c(p_vals_nl, p_overall, p_overall_c)
			
			## Get estimates 
			behaviour <- summary(xmdl)$coefficients[,1]
			behaviour <- behaviour[2:length(behaviour)]
			behaviour <- data.frame(behaviour)
			colnames(behaviour) <- c("Estimate")
			behaviour$Effect <- rownames(behaviour)
			behaviour$Module <- colnames(module)[1]
			behaviour$DISEASE <- "LATE.DEATH.TIMEPOINT"
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
			behaviour$Effect[behaviour$Effect=="Mortality2LATE.DEATH"] <- "MORTALITY.LATEDEATH"
			behaviour$Effect[behaviour$Effect=="Mortality2EARLY.DEATH"] <- "MORTALITY.EARLYDEATH"
			behaviour$Effect[behaviour$Effect=="DAY:Mortality2LATE.DEATH"] <- "TIMEPOINTxMORTALITY.LATEDEATH"
			behaviour$Effect[behaviour$Effect=="DAY:Mortality2EARLY.DEATH"] <- "TIMEPOINTxMORTALITY.EARLYDEATH"
			behaviour_timepoint <- rbind(behaviour_timepoint, behaviour)
			
			## Bind the results for the timepoint assessment!!!!!
			p_vals_s <- c(modules_interaction[i], "SURVIVORS.TIMEPOINT", p_vals_s)
			p_vals_ne <- c(modules_interaction[i], "EARLY.DEATH.TIMEPOINT", p_vals_ne)
			p_vals_nl <- c(modules_interaction[i], "LATE.DEATH.TIMEPOINT", p_vals_nl)
			p_all <- data.frame(rbind(p_vals_s, p_vals_ne,p_vals_nl ))
			colnames(p_all) <- c("MODULE", "SURVIVAL", "TIMEPOINT", "AGE", "SEX.MALE", "COMORBIDITIES.CI", "ANOVA_NULL", "ANOVA_COV")
			p_hd_extended_timepoint <- rbind(p_hd_extended_timepoint, p_all)
			###########################################################
			
			##----------------------------------------------------------------------------
			#-----------------------------------------------------------------
			###########################################d
			###########################################
			## Day 1 MORTALITY
			module <- modules_interaction[i]
			module <- modeling_cells[,c(module, "SampleID_alternative", "DAY", "Mortality2", "Age", "Sex", "Comorbidities.Charlson_Age","Comorbidities.Charlson_Index", "Barcode")]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			modulex <- module 
			
			
			## Lets subset for health only 
			module <- module[module$DAY=="1",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			module$Mortality2 <- as.factor(module$Mortality2)
			## Run model to see if there is a difference between mortality group 
			model_formula <- formula(paste0(colnames(module)[1], "~", "Mortality2 + Age + Sex + Comorbidities.Charlson_Index "))
			model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index "))

			xmdl = lm(model_formula, module)
			xmdl_cov = lm(model_formula_covonly, module)
			xmdl.av <- aov(xmdl)
			
			### Perform post hoc test to see which parts are significant....
			postHocs <- glht(xmdl.av, linfct = mcp(Mortality2 = "Tukey"))
			confin <- confint(postHocs)
			confid <- data.frame(confin$confint)
			p_val <-  data.frame(summary(postHocs)$test$pvalues)
			final <- cbind(confid, p_val)
			colnames(final) <- c("diff", "lwr", "upr", "p adj")
			rownames(final) <- gsub(" - ", "-", rownames(final))
			tukey.test1 <- final
			tbl1 <- with(module, table(Mortality2))
			
			## Get Overall P value for mortality from anova!!!!!!!!
			p_vals_day1 <- summary(xmdl)$coefficients[,4]
			p_vals_day1 <- p_vals_day1[2:6]
			
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
			behaviour$Effect[behaviour$Effect=="Mortality2LATE.DEATH"] <- "MORTALITY.LATEDEATH"
			behaviour$Effect[behaviour$Effect=="Mortality2EARLY.DEATH"] <- "MORTALITY.EARLYDEATH"
			behaviour$Effect[behaviour$Effect=="DAY:Mortality2LATE.DEATH"] <- "TIMEPOINTxMORTALITY.LATEDEATH"
			behaviour$Effect[behaviour$Effect=="DAY:Mortality2EARLY.DEATH"] <- "TIMEPOINTxMORTALITY.EARLYDEATH"
			behaviour_x <- rbind(behaviour_x, behaviour)
		
			#-----------------------------------------------------------------------------------------------------------
			###########################################
			## Day 3 MORTALITY
			module <- modules_interaction[i]
			module <- modeling_cells[,c(module, "SampleID_alternative", "DAY", "Mortality2", "Age", "Sex", "Comorbidities.Charlson_Age","Comorbidities.Charlson_Index", "Barcode")]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			modulex <- module 
			
			## Lets subset for health only 
			module <- module[module$DAY=="3",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			module$Mortality2 <- as.factor(module$Mortality2)
			
			## Run model to see if there is a difference between mortality group 
			model_formula <- formula(paste0(colnames(module)[1], "~", "Mortality2 + Age + Sex + Comorbidities.Charlson_Index "))
			model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index "))

			xmdl = lm(model_formula, module)
			xmdl_cov = lm(model_formula_covonly, module)
			xmdl.av <- aov(xmdl)
			
			### Perform post hoc test to see which parts are significant....
			postHocs <- glht(xmdl.av, linfct = mcp(Mortality2 = "Tukey"))
			confin <- confint(postHocs)
			confid <- data.frame(confin$confint)
			p_val <-  data.frame(summary(postHocs)$test$pvalues)
			final <- cbind(confid, p_val)
			colnames(final) <- c("diff", "lwr", "upr", "p adj")
			rownames(final) <- gsub(" - ", "-", rownames(final))
			tukey.test2 <- final
			tbl2 <- with(module, table(Mortality2))
			
			## get p val just for day (p overall will be the same!)
			p_vals_day3 <-  summary(xmdl)$coefficients[,4]
			p_vals_day3 <- p_vals_day3[2:6]
			
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
			behaviour$Effect[behaviour$Effect=="Mortality2LATE.DEATH"] <- "MORTALITY.LATEDEATH"
			behaviour$Effect[behaviour$Effect=="Mortality2EARLY.DEATH"] <- "MORTALITY.EARLYDEATH"
			behaviour$Effect[behaviour$Effect=="DAY:Mortality2LATE.DEATH"] <- "TIMEPOINTxMORTALITY.LATEDEATH"
			behaviour$Effect[behaviour$Effect=="DAY:Mortality2EARLY.DEATH"] <- "TIMEPOINTxMORTALITY.EARLYDEATH"
			behaviour_x <- rbind(behaviour_x, behaviour)
			
			###########################################
			## Day 5 MORTALITY
			module <- modules_interaction[i]
			module <- modeling_cells[,c(module, "SampleID_alternative", "DAY", "Mortality2", "Age", "Sex", "Comorbidities.Charlson_Age","Comorbidities.Charlson_Index", "Barcode")]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			modulex <- module 

			## Lets subset for health only 
			module <- module[module$DAY=="5",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			module$Mortality2 <- as.factor(module$Mortality2)
			## Run model to see if there is a difference between mortality group 
			model_formula <- formula(paste0(colnames(module)[1], "~", "Mortality2 + Age + Sex + Comorbidities.Charlson_Index "))
			model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index "))

			xmdl = lm(model_formula, module)
			xmdl_cov = lm(model_formula_covonly, module)
			xmdl.av <- aov(xmdl)
			
			### Perform post hoc test to see which parts are significant....
			postHocs <- glht(xmdl.av, linfct = mcp(Mortality2 = "Tukey"))
			confin <- confint(postHocs)
			confid <- data.frame(confin$confint)
			p_val <-  data.frame(summary(postHocs)$test$pvalues)
			final <- cbind(confid, p_val)
			colnames(final) <- c("diff", "lwr", "upr", "p adj")
			rownames(final) <- gsub(" - ", "-", rownames(final))
			tukey.test3 <- final
			tbl3 <- with(module, table(Mortality2))
			
			## get p val just for day (p overall will be the same!)
			p_vals_day5 <- summary(xmdl)$coefficients[,4]
			p_vals_day5 <- p_vals_day5[2:6]
			
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
			behaviour$Effect[behaviour$Effect=="Mortality2LATE.DEATH"] <- "MORTALITY.LATEDEATH"
			behaviour$Effect[behaviour$Effect=="Mortality2EARLY.DEATH"] <- "MORTALITY.EARLYDEATH"
			behaviour$Effect[behaviour$Effect=="DAY:Mortality2LATE.DEATH"] <- "TIMEPOINTxMORTALITY.LATEDEATH"
			behaviour$Effect[behaviour$Effect=="DAY:Mortality2EARLY.DEATH"] <- "TIMEPOINTxMORTALITY.EARLYDEATH"
			behaviour_x <- rbind(behaviour_x, behaviour)

			#--------------------------------------------------------------------------------
			###################################
			###Lets save the tukey test outpout 
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
			plot_dir_tukey2 <- paste0(outputdir, "/TukeyTest/Mortality")
			if (!dir.exists(plot_dir_tukey2)) {dir.create(plot_dir_tukey2)}
			pdf(paste0(plot_dir_tukey2,"/Mortality_Tukey", modules_interaction[i], ".pdf"), width=15, height=9)  
			par( mfrow= c(2,3) )
			plotTukeyHSD(tukey.test1, modulename, "Day 1", tukey.min, tukey.max )
			plotTukeyHSD(tukey.test2, modulename,"Day 3",  tukey.min, tukey.max )
			plotTukeyHSD(tukey.test3, modulename, "Day 5", tukey.min, tukey.max )
			barplot(tbl1, beside = TRUE, ylim=range(pretty(c(0, tbl1, tbl2, tbl3))), main=paste0(modulename, "\nDay 1"), ylab="Count", xlab="Mortality")
			barplot(tbl2, beside = TRUE, ylim=range(pretty(c(0, tbl1, tbl2, tbl3))), main=paste0(modulename, "\nDay 3"), ylab="Count", xlab="Mortality")
			barplot(tbl3, beside = TRUE, ylim=range(pretty(c(0, tbl1, tbl2, tbl3))), main=paste0(modulename, "\nDay 5"), ylab="Count", xlab="Mortality")
			dev.off()
			## Lets save the tukey test output!!!
			write.table(tukey_all, paste0(plot_dir_tukey2, "/Mortality_TukeyTest_",modulename, ".txt"), sep="\t")
			
			#------------------------------------------------------------------------------
			## Bind the rows together 
			p_vals_day1 <- c(modulename, "Day1", p_vals_day1)
			p_vals_day3 <- c(modulename, "Day3", p_vals_day3)
			p_vals_day5 <- c(modulename, "Day5", p_vals_day5)
			row_use <- data.frame(rbind(p_vals_day1, p_vals_day3, p_vals_day5))
			
			colnames(row_use) <- c("MODULE", "DAY", "MORTALITY.LATEDEATH", "MORTALITY.EARLYDEATH" ,"SEX.MALE", "AGE", "COMORBIDITIES.CI", "ANOVA_COV")
			row_use[, 3:7] <- lapply(3:7, function(x) as.numeric(row_use[[x]]))
			p_hd_extended <- rbind(row_use, p_hd_extended)
		}
		
		###########################################
		### Lets get to plotting format 
		p_hd_extended_timepoint <- data.frame(p_hd_extended_timepoint)
		p_hd_extended_timepoint <- gather(p_hd_extended_timepoint, "Effect", "p_value", TIMEPOINT:ANOVA_COV, factor_key=TRUE)
		p_hd_extended_timepoint$p_value <- as.numeric(p_hd_extended_timepoint$p_value)
		## Now we want to merge with behaviour 
		colnames(behaviour_timepoint) <- c("Estimate", "Effect", "MODULE", "SURVIVAL", "Direction")
		p_hd_extended_timepoint <- merge(p_hd_extended_timepoint, behaviour_timepoint, by=c("MODULE", "SURVIVAL", "Effect"), all.x=TRUE)
		#########
		p_use <- 0.05
		p_hd_extended_timepoint$p_value <- as.numeric(p_hd_extended_timepoint$p_value)
		## Now lets do a BH correction
		p_hd_extended_timepoint$BH <- p.adjust(as.numeric(p_hd_extended_timepoint$p_value),method="BH")
		p_hd_extended_timepoint$psig <- "ns"
		p_hd_extended_timepoint$psig[p_hd_extended_timepoint$p_value < 0.05] <- "sig"
		p_hd_extended_timepoint$BHsig <- "ns"
		p_hd_extended_timepoint$BHsig[p_hd_extended_timepoint$BH < 0.05] <- "sig"
		p_hd_extended_timepoint$Type <- "ADDITIVE MODEL EFFECT"
		p_hd_extended_timepoint$Type[p_hd_extended_timepoint$Effect %like% "ANOVA"] <- "ANOVA"
		
		p_hd_extended_timepoint$Effect <- as.character(p_hd_extended_timepoint$Effect)
		p_hd_extended_timepoint$Effect[p_hd_extended_timepoint$Effect =="ANOVA_NULL"] <- "NULL MODEL"
		p_hd_extended_timepoint$Effect[p_hd_extended_timepoint$Effect =="ANOVA_COV"] <- "COV ONLY MODEL"
		
		p_hd_extended_timepoint$Effect <- factor(p_hd_extended_timepoint$Effect, levels=c("TIMEPOINT", "AGE", "SEX.MALE", "COMORBIDITIES.CI", "NULL MODEL", "COV ONLY MODEL"))
		write.table(p_hd_extended_timepoint, paste0(outputdir, "/Mortality_Immunoglobulin_INTERACTION_ANALYSIS_timepoint_", annot, ".txt"), sep="\t")
		## Lets remove the direction for ns effect
		p_hd_extended_timepoint$Direction[p_hd_extended_timepoint$psig=="ns"] <- NA
		p_hd_extended_timepoint$Direction <- factor(p_hd_extended_timepoint$Direction, levels=c("increase", "decrease", NA))
		#######################################################################################
		
		###############################
		## NOW we repeat for the HvsD analysis 
		p_hd_extended <- data.frame(p_hd_extended)
		p_hd_extended <- gather(p_hd_extended, "Effect", "p_value", MORTALITY.LATEDEATH:ANOVA_COV, factor_key=TRUE)
		p_hd_extended$p_value <- as.numeric(p_hd_extended$p_value)
		colnames(behaviour_x) <- c("Estimate", "Effect", "MODULE", "DAY", "Direction")
		p_hd_extended <- merge(p_hd_extended, behaviour_x, by=c("MODULE", "DAY", "Effect"), all.x=TRUE)
		#########
		p_use <- 0.05
		p_hd_extended$p_value <- as.numeric(p_hd_extended$p_value)
		## Now lets do a BH correction
		p_hd_extended$BH <- p.adjust(as.numeric(p_hd_extended$p_value),method="BH")
		p_hd_extended$psig <- "ns"
		p_hd_extended$psig[p_hd_extended$p_value < 0.05] <- "sig"
		p_hd_extended$BHsig <- "ns"
		p_hd_extended$BHsig[p_hd_extended$BH < 0.05] <- "sig"
		p_hd_extended$Type <- "ADDITIVE MODEL EFFECT"
		p_hd_extended$Type[p_hd_extended$Effect  %like% "ANOVA"] <- "ANOVA"
		p_hd_extended$Effect <- as.character(p_hd_extended$Effect)
		p_hd_extended$Effect[p_hd_extended$Effect =="ANOVA_NULL"] <- "NULL MODEL"
		p_hd_extended$Effect[p_hd_extended$Effect =="ANOVA_COV"] <- "COV ONLY MODEL"
		p_hd_extended$Effect <- factor(p_hd_extended$Effect, levels=c("MORTALITY.LATEDEATH", "MORTALITY.EARLYDEATH", "AGE", "SEX.MALE", "COMORBIDITIES.CI", "COV ONLY MODEL"))
		write.table(p_hd_extended, paste0(outputdir, "/Mortality_Immunoglobulin_INTERACTION_ANALYSIS_hvsd_", annot, ".txt"), sep="\t")
		## Lets remove the direction for ns effect
		p_hd_extended$Direction[p_hd_extended$psig=="ns"] <- NA
		p_hd_extended$Direction <- factor(p_hd_extended$Direction, levels=c("increase", "decrease", NA))
		#######################################################################################
	
		
        ############################################
		#----------------------------------------------------	
		#----------------------------------------------------
		## Now we will format them both for plotting 
		# Timepoint
		p_hd_extended_timepoint <- p_hd_extended_timepoint[!is.na(p_hd_extended_timepoint$p_value),]
		p_hd_extended$DAY <- gsub("Day", "Day ", p_hd_extended$DAY)
		p_hd_extended_timepoint$BHsig <-  factor(p_hd_extended_timepoint$BHsig, levels=c("sig", "ns"))
		p_hd_extended$BHsig <- factor(p_hd_extended$BHsig, levels=c("sig", "ns"))
		p_hd_extended_timepoint$psig <-  factor(p_hd_extended_timepoint$psig, levels=c("sig", "ns"))
		p_hd_extended$psig <- factor(p_hd_extended$psig, levels=c("sig", "ns"))

		#----------------------------------------------------	
		p_hd_extended_timepoint$SURVIVAL <- gsub(".TIMEPOINT", "", p_hd_extended_timepoint$SURVIVAL)
		## Lets plot to see P value distribution
		pdf(paste0(outputdir,"/Mortality_Immunoglobulin_SimpleEffects_PDistributions.pdf"), width=8, height=4)
		px <- ggplot(p_hd_extended_timepoint) +  geom_point(aes(x=p_value, y=BH, alpha=psig, colour = BHsig)) +xlab("Model Effect P ") + ylab("BH Correced Model Effect P") +theme_classic() +labs(alpha="P < 0.05", colour="BH Adj P <0.05")+xlab("P") +ylab("BH Adj P") + scale_color_manual(values = c( 'red', "lightblue"), drop=FALSE) +geom_vline(xintercept=0.05, col="blue") + geom_hline(yintercept=0.05, col="blue")+ggtitle(paste0("P values for LR Timepoint Models"))+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		px1 <- ggplot(p_hd_extended) +  geom_point(aes(x=p_value, y=BH, alpha=psig, colour = BHsig)) +xlab("Model Effect P ") + ylab("BH Correced Model Effect P") +theme_classic() +labs(alpha="P < 0.05", colour="BH P Adj <0.05")+xlab("P") +ylab("BH Adj P") + scale_color_manual(values = c( 'red', "lightblue"), drop=FALSE) +geom_vline(xintercept=0.05, col="blue") + geom_hline(yintercept=0.05, col="blue")+ ggtitle(paste0("P values for LR HvsS Model"))+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		plot_grid(px, px1, ncol=2, align="h", axis="lbt", rel_widths=c(1,1), labels="AUTO")
		dev.off()
		
		### Lets do the overall plots 
		pdf(paste0(outputdir, "/Mortality_Immunoglobulin_SimpleEffects_Statistics_", annot, ".pdf"), width=14, height=12)
		p1 <- ggplot(p_hd_extended, aes(x=MODULE, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("Disease Status Model")) + labs(fill="P") +xlab("Gene") +ylab("Model Coefficients")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 1) + scale_color_manual(name = "BH Adj P < 0.05", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type), cols=vars(DAY),  scales="free_y", space='free_y')+ geom_point(aes(x=MODULE, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
			
		p2 <- ggplot(p_hd_extended_timepoint, aes(x=MODULE, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("Time Point Model")) + labs(fill="P") +xlab("Gene") +ylab("Model Coefficients")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 1) + scale_color_manual(name = "BH Adj P < 0.05", values = c("red", '#00000000'), drop = FALSE) + facet_grid(rows=vars(Type), cols=vars(SURVIVAL),  scales="free_y", space='free_y')+ geom_point(aes(x=MODULE, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect", alpha="P <0.05") + scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		
		plot(plot_grid(p1, p2, labels="AUTO", align="v", axis="lbt", ncol=1 ))
		
		####
		p1 <- ggplot(p_hd_extended, aes(x=MODULE, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, p_use))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("Disease Status Model")) + labs(fill="P") +xlab("Gene") +ylab("Model Coefficients")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 1) + scale_color_manual(name = "BH Adj P < 0.05", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type), cols=vars(DAY), scales="free_y", space='free_y')+ geom_point(aes(x=MODULE, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		
		p2 <- ggplot(p_hd_extended_timepoint, aes(x=MODULE, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, p_use))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("Time Point Model")) + labs(fill="P") +xlab("Gene") +ylab("Model Coefficients")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 1) + scale_color_manual(name = "BH Adj P < 0.05", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type), cols=vars(SURVIVAL), scales="free_y", space='free_y')+ geom_point(aes(x=MODULE, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		
		plot(plot_grid(p1, p2, labels="AUTO", align="v", axis="lbt" , ncol=1))
		###
		p1 <- ggplot(p_hd_extended, aes(x=MODULE, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, 0.1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("Disease Status Model")) + labs(fill="P") +xlab("Gene") +ylab("Model Coefficients")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 1) + scale_color_manual(name = "BH Adj P < 0.05", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type), cols=vars(DAY),  scales="free_y", space='free_y')+ geom_point(aes(x=MODULE, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		
		p2 <- ggplot(p_hd_extended_timepoint, aes(x=MODULE, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, 0.1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0("Time Point Model")) + labs(fill="P") +xlab("Gene") +ylab("Model Coefficients")+ geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 1) + scale_color_manual(name = "BH Adj P < 0.05", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type), cols=vars(SURVIVAL),  scales="free_y", space='free_y')+ geom_point(aes(x=MODULE, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect", alpha="P <0.05")+ scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
		
		plot(plot_grid(p1, p2, labels="AUTO", align="v", axis="lbt", ncol=1 ))		
		dev.off()
		}
	}
		

