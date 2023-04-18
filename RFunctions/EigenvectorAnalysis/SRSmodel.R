## Functions to correlathealth vs disease and timepoint!
## Comparing Sepsis Mortality and plotting difference in means!!!!!!!
## Lauren Overend
## lauren.overend@oriel.ox.ac.uk
## September 2022
library(broom)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cowplot)
library(multcomp)
library(effects)
library(grid)
library(gridExtra)
library(lme4)
library(lmerTest)

correlate_eigenvectors_srs <- function(eigenvectors, metadata, outputdir, type_receptor,srs_assignment){
	eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
	
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
	print("Merged Number of Eigenvector samples")
	print(dim(eigenvectors)[1])

	##Add srs assignments 
	srs <- read.delim(srs_assignment, header=TRUE, sep="\t")
	colnames(srs)[1] <- "sample"
	eigenvectors <- merge(eigenvectors, srs, by="sample")
	eigenvectors$SRS <- factor(eigenvectors$SRS, levels= c("SRS3", "SRS2", "SRS1", NA))
	
	print("Merged Number of Matched RNAseq and Repertoire Samples")
	print(dim(eigenvectors)[1])
	
	
	########################################################################################
	## do we subset just for sepsis? 
	## Remove samples with no SRS 
	#CURRENT SETTING JUST SEPSIS 
	eigenvectors <- eigenvectors[eigenvectors$DISEASE=="SEPSIS",]
	eigenvectors <- eigenvectors[!is.na(eigenvectors$SRS),]
	
	#print("Merged Number of Sepsis (only) Samples")
	#print(dim(eigenvectors)[1])
	
	#........................................................
	#############################################################
	#### lets do a little scatter plot giving you an idea of the shape of the data 
	cols_to_gather <- colnames(eigenvectors)[!colnames(eigenvectors)%like% "Module"]
	plots <- eigenvectors %>% gather(module, score, -c(all_of(cols_to_gather)))
	plots$score <- as.numeric(plots$score)
	plots$DAY <- gsub("Day", "", plots$DAY)
	plots$DAY <- as.numeric(plots$DAY)
	plots$module <- factor(plots$module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25", "Module_26", "Module_27", "Module_28", "Module_29", "Module_30",  "Module_31",  "Module_32",  "Module_33",  "Module_34",  "Module_35",  "Module_36",  "Module_37",  "Module_38",  "Module_39",  "Module_40" ))
	
	pdf(paste0(outputdir,"/SRS_Module_Scatter.pdf"), width=10, height=10)
	plot(ggplot(plots, aes(x=DAY, y=score, fill=SRS, colour=SRS))+geom_point(alpha=0.5, na.rm = TRUE)+facet_wrap(~module, scales="free") +theme_bw()+ geom_smooth(method='lm')+ylab("Module Score")+xlab("Day")+ggtitle(paste0("Interaction Plots ", type_receptor, " : Module~SRS")))
	dev.off()
	
	## Lets look at srs over time 
	pdf(paste0(outputdir,"/SRSq_Trejectory_Module_Scatter.pdf"), width=5, height=5)
	eigenvectors$DAY <- gsub("Day", "", eigenvectors$DAY)
	eigenvectors$DAY <- as.numeric(eigenvectors$DAY)
	plot(ggplot(eigenvectors, aes(x=DAY, y=SRSq))+geom_point(alpha=0.5, na.rm = TRUE, aes(fill=SRS, colour=SRS))+ geom_smooth(method='lm')+geom_line(aes(group=Barcode))+theme_bw()+ylab("SRSq")+xlab("Day")+ggtitle(paste0("SRSq Trajectory")))
	dev.off()
	
	#........................................................
	
	#........................................................
	## RUNNING AN INTERACTION MODEL 
	p_values_out <- c()
	myplots <- list()
	behaviour <- c()
	vif_all <- c()
	
	print("Calculating Stats!")
	pdf(paste0(outputdir, "/SRS_Adj_Plotsignore.pdf"),width=8, height=12)

	## Run a simple model to look for interaction
	for(i in 2:(length(colnames(eigenvectors))-10)){
		eigenvectors <- data.frame(eigenvectors)
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[ ,c(module, "sample", "DAY", "SRSq", "Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index", "Barcode")]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		## for some reason effects package doesnt like if not in a global environment so here we will assign then overright
		module <<- module
		## Model with Interaction  we need to run an mancova as we have multiple categorical levels within mortality! 
		model_formula <<- formula(paste0(colnames(module)[1], "~", "DAY+SRSq + Age + Sex+ Comorbidities.Charlson_Index +(1|Barcode)"))
		model_formula_basic <- formula(paste0(colnames(module)[1], "~"," 1 +(1|Barcode)"))
		model_formula_covonly <- formula(paste0(colnames(module)[1], "~", "Age + Sex + Comorbidities.Charlson_Index +(1|Barcode)"))
		
		xmdl <<- lmerTest::lmer(model_formula, module, REML=F)
		xmdl_old <- lmerTest::lmer(model_formula_basic, module, REML=F)
		xmdl_cov <- lmerTest::lmer(model_formula_covonly, module, REML=F)

		## Set up the plot to look at srsq having considered time and sample effect 
		est <-effects::Effect("SRSq",  partial.residuals=TRUE, xmdl)
		p1 <- plot(est,  ylab="Residual\nModule Score", main=colnames(eigenvectors)[i], xlab="SRSq")
		class(p1) <- "trellis"
		myplots[[(i-1)]] <-p1

		## VIF
		vif_scores <- vif(xmdl)
		vif_scores <- c(colnames(eigenvectors)[i], vif_scores)	
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
		
		#########################
		row_data <- c(colnames(module)[1], p_vals, p_overall, p_overall_c)
		row_data <- data.frame(t(row_data))
		row_data[,2:dim(row_data)[2]] <-as.numeric(row_data[,2:dim(row_data)[2]])
		colnames(row_data) <- c("Module", paste0("IM.", names(p_vals)), "ANOVA_NULL", "ANOVA_COV") 

	
		p_values_out <- rbind(row_data, p_values_out)
		
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
		estimates$Effect[estimates$Effect=="Mortality2LATE.DEATH"] <- "MORTALITY.LATEDEATH"
		estimates$Effect[estimates$Effect=="Mortality2EARLY.DEATH"] <- "MORTALITY.EARLYDEATH"
		estimates$Effect[estimates$Effect=="DAY:Mortality2LATE.DEATH"] <- "TIMEPOINTxMORTALITY.LATEDEATH"
		estimates$Effect[estimates$Effect=="DAY:Mortality2EARLY.DEATH"] <- "TIMEPOINTxMORTALITY.EARLYDEATH"
		estimates$Effect[estimates$Effect=="SRSSRS2"] <- "SRS2"
		estimates$Effect[estimates$Effect=="SRSSRS1"] <- "SRS1"
		estimates$Effect[estimates$Effect=="DAY:SRSSRS2"] <- "TIMEPOINTxSRS2"
		estimates$Effect[estimates$Effect=="DAY:SRSSRS1"] <- "TIMEPOINTxSRS1"
		estimates$Effect[estimates$Effect=="DAY:SRSq"] <- "DAYxSRSq"
		behaviour <- rbind(behaviour, estimates)	
		print("got effect estimates")
	}
	print("Calculating Stats.. DONE")
	dev.off()
	colnames(p_values_out) <- c("Module", "TIMEPOINT", "SRSq",  "AGE", "SEX.MALE", "COMORBIDITIES:CI", "ANOVA_NULL", "ANOVA_COV") 	
	#colnames(p_values_out) <- c("Module", "TIMEPOINT", "SRS2", "SRS1", "AGE", "SEX.MALE", "COMORBIDITIES:CI", "TIMEPOINTxSRS2", "TIMEPOINTxSRS1", "F_TEST") 	
	
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
	pdf(paste0(outputdir, "/SRS_AdjModel_Residual_Plot.pdf"),width=dims, height=(dims+5))
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/multiplot.R')
	plot(do.call("grid.arrange", c(myplots, ncol=col_number)))
	dev.off()
	
	pdf(paste0(outputdir,"/SRS_VIF_Score.pdf"), width=5, height=5)
	vif_all <- data.frame(vif_all)
	vif_all_g <- gather(vif_all, "Effect", "VIF", DAY:Comorbidities.Charlson_Index, factor_key=TRUE)
	vif_all_g <- vif_all_g[vif_all_g$V1 =="Module_1",]
	vif_all_g$VIF <- as.numeric(vif_all_g$VIF)
	vif_all_g$Effect <- as.character(vif_all_g$Effect)
	vif_all_g$Effect[vif_all_g$Effect=="Comorbidities.Charlson_Index"] <- "COMORBIDITIES:CI"
	vif_all_g$Effect[vif_all_g$Effect=="DAY.DISEASE"] <- "DAY*DISEASE"
	vif_all_g$Effect[vif_all_g$Effect=="DAY.Mortality2"] <- "DAY*MORTALITY"
	vif_all_g$Effect[vif_all_g$Effect=="Mortality2"] <- "MORTALITY"
	vif_all_g$Effect[vif_all_g$Effect=="Age"] <- "AGE"
	vif_all_g$Effect[vif_all_g$Effect=="Sex"] <- "SEX"
	vif_all_g$Effect <- factor(vif_all_g$Effect, levels=c("DAY", "SRSq", "AGE", "SEX", "COMORBIDITIES:CI"))
	plot(ggplot(vif_all_g, aes(x=Effect, y=VIF))+geom_col(alpha=0.5) +theme_classic()+ylab("Module Score")+xlab("Day")+ggtitle(paste0("Interaction Model VIFF Scores")) +facet_wrap(~V1)+geom_hline(yintercept=5, col="red") + xlab("Model Effect") + ylab("VIFF")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
	dev.off()
	
	
	## GATHER_TEST_RESULTS_for doing the interaction plot
	p_hd <- gather(p_values_out, "Effect", "p_value", TIMEPOINT:ANOVA_COV, factor_key=TRUE)
	p_hd$Module <- factor(p_hd$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25",  "Module_26",  "Module_27",  "Module_28",  "Module_29",  "Module_30",  "Module_31"))
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
	p_hd$Type <- "ADDITIVE MODEL"
	p_hd$Type[p_hd$Effect %like% "ANOVA"] <- "ANOVA"
	p_hd$Effect <- as.character(p_hd$Effect)
    p_hd$Effect[p_hd$Effect=="ANOVA_NULL"] <- "NULL MODEL"
	p_hd$Effect[p_hd$Effect=="ANOVA_COV"] <- "COV ONLY MODEL"
	p_hd$Type <- str_wrap(p_hd$Type, width = 12)
	
	#p_hd$Effect <- factor(p_hd$Effect, levels=c("TIMEPOINT", "SRS2", "SRS1", "AGE", "SEX.MALE", "COMORBIDITIES:CI", "TIMEPOINTxSRS2", "TIMEPOINTxSRS1", "MODEL"))
	p_hd$Effect <- factor(p_hd$Effect, levels=c("TIMEPOINT", "SRSq",  "AGE", "SEX.MALE", "COMORBIDITIES:CI", "NULL MODEL", "COV ONLY MODEL"))

	#### SAVE INTERACTION STATS
	write.table(p_hd, paste0(outputdir, "/Summary/SRS_INTERACTION_ANALYSIS_SUMMARY.txt"), sep="\t")

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
	
	pdf(paste0(outputdir, "/SRS_InteractionModel_Statistics.pdf"), width=12, height=6)
	p1 <- ggplot(p_hd, aes(x=Module, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Interaction Model, P<0.05")) + labs(fill="P") +xlab("Module") +ylab("Model Coefficients")+geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type),  scales="free_y", space='free_y')+ geom_point(aes(x=Module, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect", alpha="P<0.05") + scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
	plot(plot_grid(p1, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))
	
	p2 <- ggplot(p_hd, aes(x=Module, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, p_use))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Interaction Model, P<0.05")) + labs(fill="P")+xlab("Module") +ylab("Model Coefficients")+geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type),  scales="free_y", space='free_y')+ geom_point(aes(x=Module, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect", alpha="P<0.05") + scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
	plot(plot_grid(p2, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))	
	
	p3 <- ggplot(p_hd, aes(x=Module, y=Effect)) +  geom_tile(aes(fill = as.numeric(p_value), alpha=psig), colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, 0.1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Interaction Model, P<0.05")) + labs(fill="P") +xlab("Module") +ylab("Model Coefficients")+geom_tile(aes(color=factor(BHsig, c("sig", "ns"))), fill = '#00000000', size = 0.5) + scale_color_manual(name = "BH Adj P", values = c("red", '#00000000'), drop = FALSE)+facet_grid(rows=vars(Type),  scales="free_y", space='free_y')+ geom_point(aes(x=Module, y=Effect, shape = Direction), fill="white", col="black", size=2.5) + scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease'), na.translate=FALSE, drop=FALSE) +labs(shape="Direction of Effect", alpha="P<0.05") + scale_alpha_manual(values = c(1, 0.1), drop=FALSE)
	plot(plot_grid(p3, px, ncol=2, align="h", axis="lbt", rel_widths=c(2.5,1), labels="AUTO"))		
	dev.off()
}
	