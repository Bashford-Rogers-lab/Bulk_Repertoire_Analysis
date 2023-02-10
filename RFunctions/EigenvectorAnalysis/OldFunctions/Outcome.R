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
# means = 3 mean values calcualted from raw data
# se = 3 standard errors
# categories = 3 categorical / factors that group the data
# x.axis.label = what should be plotted on the x axis
# y.axis.label = what should be plotted on the y axis

correlate_eigenvectors_outcome <- function(eigenvectors, metadata, outputdir, type_receptor){
	eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
	
	##------------------------------------------------------------------------------
	### PART 1
	## Preparing metadata
	covariates <- c("Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index")
	metadata <- read.delim(metadata, sep="\t", header=TRUE)
	## Include for just samples that are in the dataset (after filtering etc)
	metadata <- metadata[metadata$SampleID_alternative %in% rownames(eigenvectors),]
	
	## We just want day 1 
	comobs <- covariates
	mort <- colnames(metadata)[colnames(metadata)%like% "Mortality"]
	outcome <- c(mort, "Infection.Number_of_ICU_Acquired_Infections", "Death.Days_from_ICU_admission")
	metadata <- metadata[, c( "SampleID_alternative", "alternative_barcode", outcome, comobs)]
	metadata$Mortality2 <- metadata$Mortality.Classification
	metadata$Mortality2[metadata$Mortality2=="MID.DEATH"] <- "LATE.DEATH"
	
	## Add barcode!
	metadata$Barcode <- NA
		for(x in 1:length(metadata$SampleID_alternative)){
			metadata$Barcode[x] <- str_split_fixed (metadata$SampleID_alternative[x], "_", 2)[,1]
			if(metadata$SampleID_alternative[x] %like% "HV"){
			 metadata$Barcode[x] <- paste0(str_split_fixed(metadata$SampleID_alternative[x], "_", 3)[,1], "_", str_split_fixed(metadata$SampleID_alternative[x], "_", 3)[,2])
			}
		}
		
		
	## subset just for sepsis
	eigenvectors <- eigenvectors[eigenvectors$DISEASE=="SEPSIS",]
	eigenvectors <- merge(eigenvectors, metadata, by.x="sample", by.y="SampleID_alternative", all.x=TRUE)

	#........................................................
	#############################################################
	#### lets do a little scatter plot giving you an idea of the shape of the data 
	cols_to_gather <- colnames(eigenvectors)[!colnames(eigenvectors)%like% "Module"]
	plots <- eigenvectors %>% gather(module, score, -c(all_of(cols_to_gather)))
	plots$score <- as.numeric(plots$score)
	plots$DAY <- gsub("Day", "", plots$DAY)
	plots$DAY <- as.numeric(plots$DAY)
	plots$module <- factor(plots$module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25", "Module_26", "Module_27", "Module_28", "Module_29", "Module_30",  "Module_31",  "Module_32",  "Module_33",  "Module_34",  "Module_35",  "Module_36",  "Module_37",  "Module_38",  "Module_39",  "Module_40" ))
	pdf(paste0(outputdir,"/MortalityScatter.pdf"), width=10, height=10)
	plot(ggplot(plots, aes(x=DAY, y=score, fill=Mortality2, colour=Mortality2))+geom_point(alpha=0.5)+facet_wrap(~module, scales="free") +theme_bw()+ geom_smooth(method='lm')+ylab("Module Score")+xlab("Day")+ggtitle(paste0("Interaction Plots ", type_receptor, " : Module~Mortality")))
	dev.off()
	#........................................................
	
	#........................................................
	## RUNNING AN INTERACTION MODEL 
	p_values_out <- c()
	cols <- as.numeric(length(2:(length(colnames(eigenvectors))-15)))
	myplots <- list()
	## Run a simple model to look for interaction
	for(i in 2:(length(colnames(eigenvectors))-16)){
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "Mortality2", "Age", "Sex", "Comorbidities.Charlson_Age","Comorbidities.Charlson_Index")]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		
		## Model with Interaction  we need to run an mancova as we have multiple categorical levels within mortality! 
		## ANOVA!!!!!!!!!!!
		model_formula <- formula(paste0(colnames(module)[1], "~", "DAY", "*","Mortality2 + Age + Sex  + Comorbidities.Charlson_Index"))
		xmdl = lm(model_formula, module)
		## Extract the model p vals for each term 
		p_vals <- summary(xmdl)$coefficients[,4]
		p_vals <- p_vals[2:length(p_vals)]
		
		xmdl_1 <- data.frame(glance(xmdl))
		p_overall <- xmdl_1$p.value
		row_data <- c(colnames(module)[1], p_vals, p_overall)
		row_data <- data.frame(t(row_data))
		
		row_data[,2:dim(row_data)[2]] <-as.numeric(row_data[,2:dim(row_data)[2]])
		colnames(row_data) <- c("Module", paste0("IM.", names(p_vals)), "IM.Model_F_stat_p_value") 

		## Can we do an interaction plot of residuals
		p1 <- interactions::interact_plot(xmdl, pred = DAY, modx = Mortality2, partial.residuals = TRUE, interval = TRUE )
		myplots[[(i-1)]] <- p1
		p_values_out <- rbind(row_data, p_values_out)
	}
	colnames(p_values_out) <- c("Module", "TIMEPOINT", "LATE.DEATH", "SURVIVOR", "AGE", "SEX", "COMORBIDITIES:CI", "TIMEPOINTxLATE.DEATH", "TIMEPOINTxSURVIVOR", "F_TEST") 	
	
	## How many columns for doing the Residual Plots!!!!
	if(type_receptor == "TCRGD"){
		col_number <- 3
		dims <- 15
	} else if(type_receptor=="TCRAB"){
		col_number <- 4
		dims <- 16
	} else if (type_receptor=="BCR"){
		col_number <- 5 
		dims <- 20
	}
	pdf(paste0(outputdir, "/Mortality_Residual.pdf"),width=(dims+5), height=dims)
	source('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/ModuleSelection/multiplot.R')
	multiplot(plotlist = myplots, cols = col_number)
	dev.off()
	
	## GATHER_TEST_RESULTS_for doing the interaction plot
	p_hd <- gather(p_values_out, "Effect", "p_value", TIMEPOINT:F_TEST, factor_key=TRUE)
	p_hd$Module <- factor(p_hd$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25", "Module_26", "Module_27", "Module_28", "Module_29", "Module_30",  "Module_31",  "Module_32",  "Module_33",  "Module_34",  "Module_35",  "Module_36",  "Module_37",  "Module_38",  "Module_39",  "Module_40" ))
	
	## Lets PLOT The statistics for them  
	p_use <- 0.05
	p_hd$p_value <- as.numeric(p_hd$p_value)
	
	
	pdf(paste0(outputdir, "/Mortality _InteractionModel.pdf"), width=7, height=4)
	plot(ggplot(p_hd, aes(x=Module, y=Effect, alpha=ifelse(p_value <= 0.05, 1, 0.9))) +  geom_tile(aes(fill = as.numeric(p_value)))+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Interaction Model, p<0.05")) + labs(fill="P Value") + guides(alpha="none")+xlab("Module") +ylab("Model Coefficients"))
	plot(ggplot(p_hd, aes(x=Module, y=Effect, alpha=ifelse(p_value <= 0.05, 1, 0.9), fill = as.numeric(p_value))) +  geom_tile(colour="lightblue")+theme_classic()+ scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, p_use))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(paste0(type_receptor, " Interaction Model, p<0.05")) + labs(fill="P Value") + guides(alpha="none")+xlab("Module") +ylab("Model Coefficients"))
	dev.off()
	write.table(p_hd, paste0(outputdir, "/Summary/INTERACTION_ANALYSIS_mortality.txt"), sep="\t")
	
	
	##------------------------------------------------------------------------------------------------------
	## Does health vs disease significantly interact for any models??
	## If so we will need to investigate these modules at individual levels. 
	if(any(p_hd$p_value[p_hd$Effect=="TIMEPOINTxSURVIVOR"] < 0.05  |p_hd$p_value[p_hd$Effect=="TIMEPOINTxLATE.DEATH"] < 0.05)){
		print("Presence of significant interaction")
		interaction_pres <- "YES"
		modules_interaction <- p_values_out$Module[p_values_out$TIMEPOINTxSURVIVOR < 0.05 | p_values_out$TIMEPOINTxLATE.DEATH < 0.05]
	} else {
		print("No interaction effect between timepoint and disease status")
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
		cols_uses <- c()
		p_values <- c()
		p_hd_extended <- c()
		corr_values <- c()
		for(i in 1:length(modules_interaction)){
			####################################################
			## Lets do Survivors vs Mortality.Overall
			module <- modules_interaction[i]
			module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE", "Barcode")]
			#################################################################
		
		## Survivors over time 
			modulex <- module
			module <- module[module$Mortality.Overall=="SURVIVOR",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			# Run model to see effect of time ## using it as a continuous variable? 
			model_formula <- formula(paste0(colnames(module)[2], "~", "Barcode+DAY"))#+Barcode
			xmdl = lm(model_formula, module)
			count <- (length(unique(module$Barcode))+1)
			p_vals_s <- summary(xmdl)$coefficients[c(1,count),4]
			## get p val just for day (p overall will be the same!)
			p_vals_s <- p_vals_s[2]

			#########################
			## Non survivors  Early Death?  
			module <- modulex
			module <- module[module$Mortality2=="EARLY.DEATH",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			# Run model to see effect of time ## using it as a continuous variable? 
			model_formula <- formula(paste0(colnames(module)[2], "~", "Barcode+DAY"))#+Barcode
			xmdl = lm(model_formula, module)
			count <- (length(unique(module$Barcode))+1)
			p_vals_ne <- summary(xmdl)$coefficients[c(1,count),4]
			## get p val just for day (p overall will be the same!)
			p_vals_ne <- p_vals_ne[2]
			
			#####################################################################
			## Non survivors  Mid/Late Death?  
			module <- modulex
			module <- module[module$Mortality2=="LATE.DEATH",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			## Run model to see effect of time ## using it as a continuous variable? 
			model_formula <- formula(paste0(colnames(module)[2], "~", "Barcode+DAY"))#+Barcode
			xmdl = lm(model_formula, module)
			count <- (length(unique(module$Barcode))+1)
			p_vals_nl <- summary(xmdl)$coefficients[c(1,count),4]
			## get p val just for day (p overall will be the same!)
			p_vals_nl <- p_vals_nl[2]
			
			###########################
			## Difference at Perhaps just T1 ?
			
			###########################################
			## Day 1 MORTALITY
			module <- modules_interaction[i]
			module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
			## Lets subset for health only 
			module <- module[module$DAY=="Day1",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			module <- module[module$DISEASE=="SEPSIS",]
			module <- merge(module, metadata, by.x="sample", by.y="SampleID_alternative")
			module$Mortality2 <- as.factor(module$Mortality2)
			## Run model to see if there is a difference between mortality group 
			model_formula <- formula(paste0(colnames(module)[2], "~", "Mortality2+Sex+Age+Comorbidities.Charlson_Age+Comorbidities.Charlson_Index"))#+Comorbidities.Charlson_Age+Comorbidities.Charlson_Index
			xmdl = lm(model_formula, module)
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
			## Get Overall P value for mortality from anova!!!!!!!!
			p_vals_day1 <- summary(xmdl.av)[[1]][["Pr(>F)"]][1]
			tbl1 <- with(module, table(Mortality2))
			
			###########################################
			## Day 3 MORTALITY
			module <- modules_interaction[i]
			module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
			## Lets subset for health only 
			module <- module[module$DAY=="Day3",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			module <- module[module$DISEASE=="SEPSIS",]
			module <- merge(module, metadata, by.x="sample", by.y="SampleID_alternative")
			module$Mortality2 <- as.factor(module$Mortality2)
			## Run model to see if there is a difference between mortality group 
			model_formula <- formula(paste0(colnames(module)[2], "~", "Mortality2+Sex+Age+Comorbidities.Charlson_Age+Comorbidities.Charlson_Index"))#+Comorbidities.Charlson_Age+Comorbidities.Charlson_Index
			xmdl = lm(model_formula, module)
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
			## get p val just for day (p overall will be the same!)
			p_vals_day3 <- summary(xmdl.av)[[1]][["Pr(>F)"]][1]
			tbl2 <- with(module, table(Mortality2))
			
			###########################################
			## Day 5 MORTALITY
			module <- modules_interaction[i]
			module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
			## Lets subset for health only 
			module <- module[module$DAY=="Day5",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			module <- module[module$DISEASE=="SEPSIS",]
			module <- merge(module, metadata, by.x="sample", by.y="SampleID_alternative")
			module$Mortality2 <- as.factor(module$Mortality2)
			## Run model to see if there is a difference between mortality group 
			model_formula <- formula(paste0(colnames(module)[2], "~", "Mortality2+Sex+Age+Comorbidities.Charlson_Age+Comorbidities.Charlson_Index"))#+Comorbidities.Charlson_Age+Comorbidities.Charlson_Index
			xmdl = lm(model_formula, module)
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
			## get p val just for day (p overall will be the same!)
			p_vals_day5 <- summary(xmdl.av)[[1]][["Pr(>F)"]][1]
			tbl3 <- with(module, table(Mortality2))
			
			###################################
			###Lets save the tukey test outpout 
			modulename <- colnames(eigenvectors)[i]
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
			pdf(paste0(plot_dir_tukey2,"/Mortality_Tukey", colnames(eigenvectors)[i], ".pdf"), width=15, height=9)  
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
			###################################
			## Look at whether days to death correlates with module - NON SURVIVORS  
			## DAY 1 - Days to death!
			module <- modules_interaction[i]
			module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
			## Lets subset for health only 
			module <- module[module$DAY=="Day1",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			module <- module[module$DISEASE=="SEPSIS",]
			module <- merge(module, metadata, by.x="sample", by.y="SampleID_alternative")
			module$Mortality2 <- as.factor(module$Mortality2)
			module$Death.Days_from_ICU_admission <- as.numeric(module$Death.Days_from_ICU_admission)
			module <- module[module$Mortality.Overall =="NON SURVIVOR",]
			## Run model to see if there is a correlation with days to death and day1
			model_formula <- formula(paste0(colnames(module)[2], "~", "Death.Days_from_ICU_admission+Sex+Age+Comorbidities.Charlson_Age+Comorbidities.Charlson_Index"))
			xmdl = lm(model_formula, module)
			p_vals_daysdeath <- summary(xmdl)$coefficients[,4]
			## get p val just for day (p overall will be the same!)
			p_vals_daysdeath <- p_vals_daysdeath[2]
			estimate_daysdeath <- summary(xmdl)$coefficients[,1]
			estimate_daysdeath <- estimate_daysdeath[2]
			
			##################################
			## Look at whether days to death correlates with module - NON SURVIVORS  
			## DAY 1 - Number of infections
			module <- modules_interaction[i]
			module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
			## Lets subset for health only 
			module <- module[module$DAY=="Day1",]
			module$DAY <- gsub("Day", "", module$DAY)
			module$DAY <- as.numeric(module$DAY)
			#module <- module[module$DISEASE=="SEPSIS",]
			module <- merge(module, metadata, by.x="sample", by.y="SampleID_alternative")
			module$Mortality2 <- as.factor(module$Mortality2)
			module$Death.Days_from_ICU_admission <- as.numeric(module$Death.Days_from_ICU_admission)
			#module <- module[module$Mortality.Overall =="NON SURVIVOR",]
			## Run model to see if there is a correlation with days to death and day1
			model_formula <- formula(paste0(colnames(module)[2], "~", "Infection.Number_of_ICU_Acquired_Infections+Sex+Age+Comorbidities.Charlson_Age+Comorbidities.Charlson_Index"))
			xmdl = lm(model_formula, module)
			p_vals_ns_in <- summary(xmdl)$coefficients[,4]
			## get p val just for day (p overall will be the same!)
			p_vals_ns_in <- p_vals_ns_in[2]
			estimate_inf <- summary(xmdl)$coefficients[,1]
			estimate_inf <- estimate_inf[2] 

			###################Lets put it all together
			row_use <- c(modules_interaction[i], p_vals_s, p_vals_ne, p_vals_nl, p_vals_day1, p_vals_day3, p_vals_day5, p_vals_daysdeath, p_vals_ns_in)
			row_use[2:5] <- as.numeric(row_use[2:5])
			p_hd_extended <- rbind(row_use, p_hd_extended)
			
			### Lets make a seperate frame for the correlation
			row1 <-  c(modules_interaction[i], estimate_daysdeath, estimate_inf)
			corr_values <- rbind(row1, corr_values)
		}
		
		###########################################
		### Lets get to plotting format 
		p_hd_extended <- data.frame(p_hd_extended)
		colnames(p_hd_extended) <- c("Module", "Survivor.Timepoint", "NonSurvivor.EarlyDeath.Timepoint", "NonSurvivor.LateDeath.Timepoint", "Mortality.Day1", "Mortality.Day3", "Mortality.Day5", "NonSurvivor.Days_to_Death", "Number_of_ICU_Acquired_Infections")
		p_hd_extended <- gather(p_hd_extended, "Effect", "p_value", Survivor.Timepoint:Number_of_ICU_Acquired_Infections, factor_key=TRUE)
		p_hd_extended$p_value <- as.numeric(p_hd_extended$p_value)
		p_hd_extended$Module <- factor(p_hd_extended$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25", "Module_26", "Module_27", "Module_28", "Module_29", "Module_30",  "Module_31",  "Module_32",  "Module_33",  "Module_34",  "Module_35",  "Module_36",  "Module_37",  "Module_38",  "Module_39",  "Module_40" ))
		
		
		colnames(corr_values) <- c("Module", "NonSurvivor.Days_to_Death", "Number_of_ICU_Acquired_Infections")
		corr_values <- data.frame(corr_values)
		corr_values$NonSurvivor.Days_to_Death<- as.numeric(corr_values$NonSurvivor.Days_to_Death)
		corr_values$Number_of_ICU_Acquired_Infections<- as.numeric(corr_values$Number_of_ICU_Acquired_Infections)
		
		p_hd_extended$Difference <- NA
		p_hd_extended$Normalised <- NA
		p_hd_extended_2 <- c()
		## Now we want to say wether its higher in sepsis or health 
		for(i in 1:length(p_hd_extended[,1])){
			row_used <- p_hd_extended[i,]
			if(row_used$Effect %like% "Mortality"){
				module_use <- as.character(row_used[,1])
				data_use <- eigenvectors[,c(module_use, "sample", "DAY", "DISEASE")]
				data_use <- merge(data_use, metadata, by.x="sample", by.y="SampleID_alternative")
				if(row_used$Effect %like% "Day1"){
					DAY <- "Day1"
				}
				if(row_used$Effect %like% "Day3"){
					DAY <- "Day3"
				}
				if(row_used$Effect %like% "Day5"){
					DAY <- "Day5"
				}
				data_use <- data_use[data_use$DAY ==DAY,]
				data_use$DAY <- gsub("Day", "", data_use$DAY)
				data_use$DAY <- as.numeric(data_use$DAY)
				my_sym <- sym(module_use) 
				means <- data_use %>% dplyr::group_by(Mortality.Overall) %>% dplyr::summarize(Mean = mean(!!my_sym, na.rm=TRUE))
				means <- data.frame(means)
				diff_means <- means$Mean[means$Mortality.Overall=="NON SURVIVOR"] -  means$Mean[means$Mortality.Overall=="SURVIVOR"]
				diff_means_normalised <- diff_means/means$Mean[means$Mortality.Overall=="SURVIVOR"]
				row_used$Difference <- diff_means
				row_used$Normalised <- diff_means_normalised 
			}
			if(row_used$Effect %like% "Timepoint"){
				module_use <- as.character(row_used[,1])
				data_use <- eigenvectors[,c(module_use, "sample", "DAY", "DISEASE")]
				data_use$DAY <- gsub("Day", "", data_use$DAY)
				data_use$DAY <- as.numeric(data_use$DAY)
				data_use <- merge(data_use, metadata, by.x="sample", by.y="SampleID_alternative")
				my_sym <- sym(module_use) 
				means <- data_use %>% dplyr::group_by(Mortality2, DAY) %>% dplyr::summarize(Mean = mean(!!my_sym, na.rm=TRUE))
				means <- data.frame(means)
				if(row_used$Effect %like% "NonSurvivor.EarlyDeath.Timepoint"){
					diff_means <- means$Mean[means$Mortality2=="EARLY.DEATH" & means$DAY==5] -  means$Mean[means$Mortality2=="EARLY.DEATH" & means$DAY==1]
					diff_means_normalised <- diff_means/means$Mean[means$Mortality2=="EARLY.DEATH" & means$DAY==5]
				} else if (row_used$Effect %like% "NonSurvivor.LateDeath"){
					diff_means <- means$Mean[means$Mortality2=="LATE.DEATH" & means$DAY==5] -  means$Mean[means$Mortality2=="LATE.DEATH" & means$DAY==1]
					diff_means_normalised <- diff_means/means$Mean[means$Mortality2=="LATE.DEATH" & means$DAY==5]
				} else {
					diff_means <- means$Mean[means$Mortality2=="SURVIVOR" & means$DAY==5] -  means$Mean[means$Mortality2=="SURVIVOR" & means$DAY==1]
					diff_means_normalised <- diff_means/means$Mean[means$Mortality2=="SURVIVOR" & means$DAY==5]
				}
				row_used$Difference <- diff_means
				row_used$Normalised <- diff_means_normalised 
			}
			if(row_used$Effect %like% "NonSurvivor.Days_to_Death" | row_used$Effect %like% "Number_of_ICU_Acquired_Infections"){
				### What direction is correlation 
				module_use <- as.character(row_used[,1])
				effect_use <- as.character(row_used[,2])
				diff_means <- corr_values[corr_values$Module ==module_use,]
				diff_means <- diff_means[,effect_use]
				diff_means_normalised <- NA
				row_used$Difference <- diff_means
				row_used$Normalised <- diff_means_normalised 
			}
			p_hd_extended_2 <- rbind(row_used, p_hd_extended_2)		
		}
		
		p_hd_extended_2$Difference <- as.numeric(p_hd_extended_2$Difference)
		p_hd_extended_2$Normalised <- as.numeric(p_hd_extended_2$Normalised)
		p_hd_extended_2$High <- NA
		p_hd_extended_2$High[p_hd_extended_2$Difference> 0 & p_hd_extended_2$Effect %like% "Mortality"] <- "High"
		p_hd_extended_2$High[p_hd_extended_2$Difference< 0& p_hd_extended_2$Effect %like% "Mortality"] <- "Low"
		p_hd_extended_2$High[p_hd_extended_2$Difference> 0 & !p_hd_extended_2$Effect %like% "Mortality"] <- "High"
		p_hd_extended_2$High[p_hd_extended_2$Difference< 0& !p_hd_extended_2$Effect %like% "Mortality"] <- "Low"
		p_hd_extended_2$Type <- NA
		p_hd_extended_2$Type[p_hd_extended_2$Effect %like% "Mortality"] <- "SvsNS"
		p_hd_extended_2$Type[p_hd_extended_2$Effect %like% "Timepoint"] <- "Timepoint"
		p_hd_extended_2$Type[p_hd_extended_2$Effect %like% "Days_to_Death" |p_hd_extended_2$Effect %like% "Number_of_ICU_Acquired_Infections" ] <- "Dif"
		
		## How many tests run for each question? e.g. number of modules
		no_modules <- length(unique(p_hd_extended_2$Module))
		
		#p_hd_extended_2$Effect <- gsub("_", " ", p_hd_extended_2$Effect)
		## Set p value of interest 
		p_use <- 0.05
		### Lets calculate corrected p value using Benjamin Hochburg
		p_hd_extended_2$BH <- p.adjust(as.numeric(p_hd_extended_2$p_value),method="BH")

		## Lets plot to see P value distribution
		pdf(paste0(outputdir,"/Mortality_PDistribution.pdf"), width=10, height=10)
		plot(ggplot(p_hd_extended_2) +
		  geom_histogram(aes(x = p_value, fill=Effect), breaks = seq(0, 1, 0.05),
						 color = "black") +
		  # Remove space between x-axis and min(y)
		  scale_y_continuous(expand = expansion(c(0, 0.05))) +
		  facet_wrap(vars(Effect)) + # separate plots
		  theme_bw(base_size = 12)	+ geom_vline(xintercept=0.05, col="red") + xlab("Raw p value") +ylab("Count")+ggtitle(paste0("Distribution of P values for LR Models ", type_receptor)) +labs(fill="Model")+guides(fill="none"))
		dev.off()  
		  
		######################################################
		## 1. Lets plot the patient demographs just so we have an idea of idea of number of samples
		no_samples <- data.frame(table(eigenvectors$Mortality2, eigenvectors$DAY))
		colnames(no_samples) <- c("Mortality", "Day", "Count")
		p1 <- ggplot(no_samples, aes(x=as.factor(Mortality), y=Count, fill=as.factor(Mortality))) + geom_col() +theme_bw() + facet_wrap(~Day) +ggtitle("Sepsis BCR Cohort")+ xlab("Mortality") + ylab("Number of Samples") +guides(fill="none")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
		
		################################
		## Lets plot 
		
		q <- p_hd_extended_2[p_hd_extended_2$p_value <p_use,]
		pdf(paste0(outputdir,"/Mortality.pdf"), width=15, height=5)
		p2 <- ggplot( p_hd_extended_2 ) +
		geom_tile(aes(x=Module, y=Effect, fill = as.numeric(p_value), alpha=ifelse(p_value < p_use, 1, 0.9))) + 
		scale_fill_continuous(type = "viridis")+
		geom_point(data = q, aes(x=Module, y=Effect, shape = High, colour=Type), fill="White") + 
		scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease')) +
		scale_colour_discrete(labels = c("Direction of Effect", "Survivors to Non Survivors", "Timepoint 1 to 5")) +
		theme_classic() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
		labs(fill="p value") + guides(alpha="none")+xlab("Module") +ylab("LR Model")
		plot(plot_grid(p2, p1, ncol=2, rel_widths = c(2, 1)))
		######
		modules_interesting2 <- p_hd_extended_2
		q <- modules_interesting2[modules_interesting2$p_value <p_use,]
		modules_interesting2$p_value[modules_interesting2$p_value >= p_use] <- NA
		## Lets plot just these with a better scale 
		p2 <- ggplot( modules_interesting2 ) +
		geom_tile(aes(x=Module, y=Effect, fill = as.numeric(p_value)), colour="lightblue") + 
		scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, p_use))+
		geom_point(data = q, aes(x=Module, y=Effect, shape = High, colour=Type), fill="White") + 
		scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease')) +
		scale_colour_discrete(labels = c("Direction of Effect", "Survivors to Non Survivors", "Timepoint 1 to 5")) +
		theme_classic() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
		labs(fill="p value") + guides(alpha="none")+ xlab("Module") +ylab("LR Model")
		plot(plot_grid(p2, p1, ncol=2, rel_widths = c(2, 1)))
		dev.off()
		
		q <- p_hd_extended_2[p_hd_extended_2$BH <p_use,]
		pdf(paste0(outputdir,"/Mortality_AdjustedPvalue.pdf"), width=9, height=4)
		plot(ggplot( p_hd_extended_2 ) +
		geom_tile(aes(x=Module, y=Effect, fill = as.numeric(BH), alpha=ifelse(BH < p_use, 1, 0.9))) + 
		scale_fill_continuous(type = "viridis")+
		geom_point(data = q, aes(x=Module, y=Effect, shape = High, colour=Type), fill="White") + 
		scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease')) +
		scale_colour_discrete(labels = c("Survivors to Non Survivors", "Timepoint 1 to 5")) +
		theme_classic() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
		labs(fill="BH Adj p Value") + guides(alpha="none")+xlab("Module") +ylab("LR Model"))
		######
		modules_interesting2 <- p_hd_extended_2
		q <- modules_interesting2[modules_interesting2$BH <p_use,]

		modules_interesting2$BH[modules_interesting2$BH >= p_use] <- NA
		## Lets plot just these with a better scale 
		plot(ggplot( modules_interesting2 ) +
		geom_tile(aes(x=Module, y=Effect, fill = as.numeric(BH)), colour="lightblue") + 
		scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, p_use))+
		geom_point(data = q, aes(x=Module, y=Effect, shape = High, colour=Type), fill="White") + 
		scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease')) +
		scale_colour_discrete(labels = c("Survivors to Non Survivors", "Timepoint 1 to 5")) +
		theme_classic() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
		labs(fill="BH Adj p Value") + guides(alpha="none")+ xlab("Module") +ylab("LR Model"))
		dev.off()
		
		write.table(p_hd_extended_2, paste0(outputdir, "/Summary/Clinical_variables_hvsd.txt"), sep="\t")
		return(p_hd_extended_2)
		}
}

###### DONE !!!