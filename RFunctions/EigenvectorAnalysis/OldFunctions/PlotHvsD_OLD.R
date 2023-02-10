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


correlate_eigenvectors_hvd <- function(eigenvectors, metadata, outputdir, type_receptor, metahealth){
	eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
	p_values_hd <- c()
	covariates <- c("Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index")
	metadata <- read.delim(metadata, sep="\t", header=TRUE)
	metadata <- metadata[, c("SampleID_alternative", covariates)]
	health <- read.delim(metahealth, sep="\t", header=TRUE)
	health <- health[, c("SampleID", "Barcode", "Age", "Sex")]
	
	####### Lets add the age and sex covariates to the eigenvectors
	eigenvectors$Age <- NA
	eigenvectors$Sex <- NA
	for(i in 1:length(eigenvectors$sample)){
		if(eigenvectors$sample[i] %like% "HV"){
			eigenvectors$Age[i] <- health$Age[health$SampleID== eigenvectors$sample[i]]
			eigenvectors$Sex[i] <- health$Sex[health$SampleID== eigenvectors$sample[i]]
		} else {
			eigenvectors$Age[i] <- metadata$Age[metadata$SampleID_alternative== eigenvectors$sample[i]]
			eigenvectors$Sex[i] <- metadata$Sex[metadata$SampleID_alternative== eigenvectors$sample[i]]
		}
	}
	#........................................................
	#############################################################
	#### lets do a little scatter plot giving you an idea of the shape of the data 
	plots <- eigenvectors %>% gather(module, score, -c(sample, DAY, DISEASE, Age, Sex)) 
	plots$DAY <- gsub("Day", "", plots$DAY)
	plots$DAY <- as.numeric(plots$DAY)
	plots$module <- factor(plots$module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25"))
	pdf(paste0(outputdir,"/HealthvsDisease_Scatter.pdf"), width=10, height=10)
	plot(ggplot(plots, aes(x=DAY, y=score, fill=DISEASE, colour=DISEASE))+geom_point(alpha=0.5)+facet_wrap(~module, scales="free") +theme_bw()+ geom_smooth(method='lm')+ylab("Module Score")+xlab("Day")+ggtitle(paste0("Interaction Plots ", type_receptor, " : Module~DISEASE")))
	dev.off()
	#........................................................

	
	## Run a simple model to look for interaction (no covariates included!)
	for(i in 1:(length(colnames(eigenvectors))-5)){
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		
		## Get barcode 
		module$Barcode <- NA
		for(x in 1:length(module$sample)){
			module$Barcode[x] <- str_split_fixed (module$sample[x], "_", 2)[,1]
			if(module$sample[x] %like% "HV"){
			 module$Barcode[x] <- paste0(str_split_fixed(module$sample[x], "_", 3)[,1], "_", str_split_fixed(module$sample[x], "_", 3)[,2])
			}
		}
		
		## Model with Interaction
		model_formula <- formula(paste0(colnames(module)[1], "~", "DAY", "*","DISEASE"))
		xmdl = lm(model_formula, module)
		p_vals <- summary(xmdl)$coefficients[,4]
		xmdl_1 <- data.frame(glance(xmdl))
		p_overall <- xmdl_1$p.value
		row_data <- c(colnames(module)[1], p_vals, p_overall)
		row_data <- data.frame(t(row_data))
		row_data[,2:6] <-as.numeric(row_data[,2:6])
		colnames(row_data) <- c("Module", "IM.Model_p_value", paste0("IM.", names(p_vals))) 
		
		## Model without interaction#
		model_formula <- formula(paste0(colnames(module)[1], "~", "DAY", "+","DISEASE"))
		xmdl = lm(model_formula, module)
		p_vals <- summary(xmdl)$coefficients[,4]
		xmdl_1 <- data.frame(glance(xmdl))
		p_overall <- xmdl_1$p.value
		row_data2 <- c(colnames(module)[1], p_vals, p_overall)
		row_data2 <- data.frame(t(row_data2))
		row_data2[,2:5] <-as.numeric(row_data2[,2:5])
		colnames(row_data2) <- c("Module", "Model_p_value", names(p_vals))
		final_row <- cbind(row_data2, row_data[,2:6])
		p_values_hd <- rbind(final_row, p_values_hd)
	}
	 
	colnames(p_values_hd) <- c("Module", "MODEL", "INTERCEPT", "DAY", "DISEASE",  "INT.MODEL", "INT.INTERCEPT", "INT.DAY", "INT.DISEASE","INT.INTERACTION") 
	p_hd <- gather(p_values_hd, "Effect", "p_value", MODEL:INT.INTERACTION, factor_key=TRUE)
	p_hd$Module <- factor(p_hd$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25"))
	orders <- unique(p_hd[, c("Effect")])
	a <- ifelse(orders %like% "INT\\.", "red", "black")
	pdf(paste0(outputdir, "/HvsD_InteractionModel.pdf"), width=5, height=4)
	plot(ggplot(p_hd, aes(x=Module, y=Effect, alpha=ifelse(p_value <= 0.05, 1, 0.9))) +  geom_tile(aes(fill = as.numeric(p_value)))+theme_classic()+ scale_fill_continuous(type = "viridis")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(colour=a))+ggtitle(paste0(type_receptor, " Interaction vs Additive Model, p<0.05")) + labs(fill="p Value") + guides(alpha="none")+xlab("Module") +ylab("Clinical Variables"))
	dev.off()
	
	write.table(p_hd, paste0(outputdir, "/Summary/Clinical_variables_hvsd_interaction.txt"), sep="\t")
	
	####################################################################################################
	### Part 2 more complex model 
	### Because of the interaction we need to assess at levels of Time / Disease state rather than across
	
	p_hd_extended <- c()
	for(i in 1:(length(colnames(eigenvectors))-5)){
		####################################################
		## Lets do Health first
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
		## Get barcode 
		module$Barcode <- NA
		for(x in 1:length(module$sample)){
			module$Barcode[x] <- str_split_fixed (module$sample[x], "_", 2)[,1]
			if(module$sample[x] %like% "HV"){
			 module$Barcode[x] <- paste0(str_split_fixed(module$sample[x], "_", 3)[,1], "_", str_split_fixed(module$sample[x], "_", 3)[,2])
			}
		}

		## 1. Lets look at whether Health Changes over time (taking into acocunt sample which negates the need for Age/Sex
		## Adding in sample this time!
		## Take the p value for the "time effect"
		module <- module[module$DISEASE=="HEALTH",]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		## Run model to see effect of time ## using it as a continuous variable? 
		model_formula <- formula(paste0(colnames(module)[1], "~", "Barcode+DAY"))#+Barcode
		xmdl = lm(model_formula, module)
		count <- (length(unique(module$Barcode))+1)
		p_vals <- summary(xmdl)$coefficients[c(1,count),4]
		## get p val just for day (p overall will be the same!)
		p_vals <- p_vals[2]

		####################################################
		## Lets do Sepsis 
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE")]
		module$Barcode <- NA
		for(x in 1:length(module$sample)){
			module$Barcode[x] <- str_split_fixed (module$sample[x], "_", 2)[,1]
			if(module$sample[x] %like% "HV"){
			 module$Barcode[x] <- paste0(str_split_fixed(module$sample[x], "_", 3)[,1], "_", str_split_fixed(module$sample[x], "_", 3)[,2])
			}
		}
		## Lets subset for health only  ## Adding in sample this time!
		module <- module[module$DISEASE=="SEPSIS",]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		## Run model to see effect of time ## using it as a continuous variable? 
		model_formula <- formula(paste0(colnames(module)[1], "~", "Barcode+DAY")) # Barcode
		xmdl = lm(model_formula, module)
		count <- (length(unique(module$Barcode))+1)
		p_vals_s <- summary(xmdl)$coefficients[c(1,count),4]
		## get p val just for day (p overall will be the same!)
		p_vals_s <- p_vals_s[2]
		
		###########################################
		###########################################
		
		## Day : Comaring HvsD (for all sepsis)
		## Should include age and sex as covariates
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE", "Sex", "Age")]
		## Lets subset for health only 
		module <- module[module$DAY=="Day1",]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		## Run model to see effect of time ## using it as a continuous variable? 
		model_formula <- formula(paste0(colnames(module)[1], "~", "DISEASE+Sex+Age"))#+Sex+Age
		xmdl = lm(model_formula, module)
		p_vals_day1 <- summary(xmdl)$coefficients[,4]
		## get p val just for day (p overall will be the same!)
		p_vals_day1 <- p_vals_day1[2]
		
		###### Would be good to add a plot like we did for the mortality!!!
		module$DISEASE <- as.factor(module$DISEASE)
		module$Sex <- as.factor(module$Sex)
		xmdl = lm(model_formula, module)
		xmdl.av <- aov(xmdl)
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
		## Get Overall P value for mortality from anova!!!!!!!!
		p_vals_day1 <- summary(xmdl.av)[[1]][["Pr(>F)"]][1]
		tbl1 <- with(module, table(DISEASE))
		
		
		############################################
		###########################################
		
		## Day 3 HvsD
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE", "Sex", "Age")]
		## Lets subset for health only 
		module <- module[module$DAY=="Day3",]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		## Run model to see effect of time ## using it as a continuous variable? 
		model_formula <- formula(paste0(colnames(module)[1], "~", "DISEASE+Sex+Age"))
		xmdl = lm(model_formula, module)
		p_vals_day3 <- summary(xmdl)$coefficients[,4]
		## get p val just for day (p overall will be the same!)
		p_vals_day3 <- p_vals_day3[2]
		
		################################
		module$DISEASE <- as.factor(module$DISEASE)
		module$Sex <- as.factor(module$Sex)
		xmdl = lm(model_formula, module)
		xmdl.av <- aov(xmdl)
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

		## get p val just for day (p overall will be the same!)
		p_vals_day3 <- summary(xmdl.av)[[1]][["Pr(>F)"]][1]
		tbl2 <- with(module, table(DISEASE))
		
		###########################################
		###########################################
		## Day 3 HvsD
		module <- colnames(eigenvectors)[i]
		module <- eigenvectors[,c(module, "sample", "DAY", "DISEASE", "Sex", "Age")]
		## Lets subset for health only 
		module <- module[module$DAY=="Day5",]
		module$DAY <- gsub("Day", "", module$DAY)
		module$DAY <- as.numeric(module$DAY)
		## Run model to see effect of time ## using it as a continuous variable? 
		model_formula <- formula(paste0(colnames(module)[1], "~", "DISEASE+Sex+Age"))
		xmdl = lm(model_formula, module)
		p_vals_day5 <- summary(xmdl)$coefficients[,4]
		## get p val just for day (p overall will be the same!)
		p_vals_day5 <- p_vals_day5[2]
		
		################################
		module$DISEASE <- as.factor(module$DISEASE)
		module$Sex <- as.factor(module$Sex)
		xmdl = lm(model_formula, module)
		xmdl.av <- aov(xmdl)
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
		## get p val just for day (p overall will be the same!)
		p_vals_day5 <- summary(xmdl.av)[[1]][["Pr(>F)"]][1]
		tbl3 <- with(module, table(DISEASE))
		
		##################################
		### Lets plot the differences in means via the tukey test 
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
		plot_dir_tukey2 <- paste0(outputdir, "/TukeyTest/HealthvsSepsis")
		if (!dir.exists(plot_dir_tukey2)) {dir.create(plot_dir_tukey2)}
		
		pdf(paste0(plot_dir_tukey2,"/Health_vs_Sepsis_Tukey", colnames(eigenvectors)[i], ".pdf"), width=8, height=7)  
		par( mfrow= c(2,3) )
		s <- plotTukeyHSD(tukey.test1, modulename, "Day 1", tukey.min, tukey.max )
		plotTukeyHSD(tukey.test2, modulename,"Day 3", tukey.min, tukey.max )
		plotTukeyHSD(tukey.test3, modulename, "Day 5", tukey.min, tukey.max )
		barplot(tbl1, beside = TRUE, ylim=range(pretty(c(0, tbl1, tbl2, tbl3))), main=paste0(modulename, "\nDay 1"), ylab="Count", xlab="Disease State")
		barplot(tbl2, beside = TRUE, ylim=range(pretty(c(0, tbl1, tbl2, tbl3))), main=paste0(modulename, "\nDay 3"), ylab="Count", xlab="Disease State")
		barplot(tbl3, beside = TRUE, ylim=range(pretty(c(0, tbl1, tbl2, tbl3))), main=paste0(modulename, "\nDay 5"), ylab="Count", xlab="Disease State")
		dev.off()
		write.table(tukey_all, paste0(plot_dir_tukey2, "/TukeyTest_HvsD_",modulename, ".txt"), sep="\t")

		###################Lets put it all together
		row_use <- c(colnames(eigenvectors)[i], p_vals, p_vals_s, p_vals_day1, p_vals_day3, p_vals_day5)
		row_use[2:5] <- as.numeric(row_use[2:5])
		p_hd_extended <- rbind(row_use, p_hd_extended)
	}
	
	### Now we move onto plotting 
	p_hd_extended <- data.frame(p_hd_extended)
	colnames(p_hd_extended) <- c("Module", "Health.Timepoint", "Sepsis.Timepoint", "HvsS.Day1", "HvsS.Day3", "HvsS.Day5")
	p_hd_extended <- gather(p_hd_extended, "Effect", "p_value", Health.Timepoint:HvsS.Day5, factor_key=TRUE)
	p_hd_extended$p_value <- as.numeric(p_hd_extended$p_value)
	p_hd_extended$Module <- factor(p_hd_extended$Module, levels=c("Module_1","Module_2", "Module_3","Module_4", "Module_5", "Module_6", "Module_7", "Module_8", "Module_9", "Module_10", "Module_11", "Module_12", "Module_13", "Module_14", "Module_15", "Module_16", "Module_17", "Module_18", "Module_19", "Module_20", "Module_21", "Module_22", "Module_23", "Module_24", "Module_25"))
	
	p_hd_extended$Difference <- NA
	p_hd_extended$Normalised <- NA
	p_hd_extended_2 <- c()
	## Now we want to say wether its higher in sepsis or health 
	for(i in 1:length(p_hd_extended[,1])){
		row_used <- p_hd_extended[i,]
		if(row_used$Effect %like% "HvsS"){
			module_use <- as.character(row_used[,1])
			data_use <- eigenvectors[,c(module_use, "sample", "DAY", "DISEASE")]
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
			means <- data_use %>% dplyr::group_by(DISEASE) %>% dplyr::summarize(Mean = mean(!!my_sym, na.rm=TRUE))
			means <- data.frame(means)
			diff_means <- means$Mean[means$DISEASE=="SEPSIS"] -  means$Mean[means$DISEASE=="HEALTH"]
			diff_means_normalised <- diff_means/means$Mean[means$DISEASE=="HEALTH"]
			row_used$Difference <- diff_means
			row_used$Normalised <- diff_means_normalised 
		}
		
		if(row_used$Effect %like% "Timepoint"){
			print(i)
			module_use <- as.character(row_used[,1])
			data_use <- eigenvectors[,c(module_use, "sample", "DAY", "DISEASE")]
			data_use$DAY <- gsub("Day", "", data_use$DAY)
			data_use$DAY <- as.numeric(data_use$DAY)
			my_sym <- sym(module_use) 
			means <- data_use %>% dplyr::group_by(DISEASE, DAY) %>% dplyr::summarize(Mean = mean(!!my_sym, na.rm=TRUE))
			means <- data.frame(means)
			if(row_used$Effect %like% "Health"){
				diff_means <- means$Mean[means$DISEASE=="HEALTH" & means$DAY==5] -  means$Mean[means$DISEASE=="HEALTH" & means$DAY==1]
				diff_means_normalised <- diff_means/means$Mean[means$DISEASE=="HEALTH" & means$DAY==5]
			} else {
				diff_means <- means$Mean[means$DISEASE=="SEPSIS" & means$DAY==5] -  means$Mean[means$DISEASE=="SEPSIS" & means$DAY==1]
				diff_means_normalised <- diff_means/means$Mean[means$DISEASE=="SEPSIS" & means$DAY==5]
			}
			row_used$Difference <- diff_means
			row_used$Normalised <- diff_means_normalised 
		}
		p_hd_extended_2 <- rbind(row_used, p_hd_extended_2)		
	}
	
	p_hd_extended_2$Difference <- as.numeric(p_hd_extended_2$Difference)
	p_hd_extended_2$Normalised <- as.numeric(p_hd_extended_2$Normalised)
	p_hd_extended_2$High <- NA
	p_hd_extended_2$High[p_hd_extended_2$Difference> 0 & p_hd_extended_2$Effect %like% "HvsS"] <- "High"
	p_hd_extended_2$High[p_hd_extended_2$Difference< 0& p_hd_extended_2$Effect %like% "HvsS"] <- "Low"
	p_hd_extended_2$High[p_hd_extended_2$Difference> 0 & !p_hd_extended_2$Effect %like% "HvsS"] <- "High"
	p_hd_extended_2$High[p_hd_extended_2$Difference< 0& !p_hd_extended_2$Effect %like% "HvsS"] <- "Low"
	p_hd_extended_2$Type <- NA
	p_hd_extended_2$Type[p_hd_extended_2$Effect %like% "HvsS"] <- "HvsS"
	p_hd_extended_2$Type[!p_hd_extended_2$Effect %like% "HvsS"] <- "Timepoint"
	
	## How many tests run for each question? e.g. number of modules
	no_modules <- length(unique(p_hd_extended_2$Module))
	
	## Set p value of interest 
	p_use <- 0.05
	### Lets calculate corrected p value using Benjamin Hochburg
	p_hd_extended_2$BH <- p.adjust(as.numeric(p_hd_extended_2$p_value),method="BH")

	## Lets plot to see P value distribution
	pdf(paste0(outputdir,"/HealthvsSepsis_PDistribution.pdf"), width=10, height=10)
	plot(ggplot(p_hd_extended_2) +
	  geom_histogram(aes(x = p_value, fill=Effect), breaks = seq(0, 1, 0.05),
					 color = "black") +
	  # Remove space between x-axis and min(y)
	  scale_y_continuous(expand = expansion(c(0, 0.05))) +
	  facet_wrap(vars(Effect)) + # separate plots
	  theme_bw(base_size = 12)	+ geom_vline(xintercept=0.05, col="red") + xlab("p value") +ylab("Count")+ggtitle(paste0("Distribution of P values for LR Models ", type_receptor)) +labs(fill="Model"))
	dev.off()  
	  
	######################################################
	q <- p_hd_extended_2[p_hd_extended_2$p_value <p_use,]
	pdf(paste0(outputdir,"/HealthvsDisease.pdf"), width=9, height=5)
	plot(ggplot( p_hd_extended_2 ) +
	geom_tile(aes(x=Module, y=Effect, fill = as.numeric(p_value), alpha=ifelse(p_value < p_use, 1, 0.9))) + 
	scale_fill_continuous(type = "viridis")+
	geom_point(data = q, aes(x=Module, y=Effect, shape = High, colour=Type), fill="White") + 
	scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease')) +
	scale_colour_discrete(labels = c("Health to Sepsis", "Timepoint 1 to 5")) +
	theme_classic() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	labs(fill="p value") + guides(alpha="none")+xlab("Module") +ylab("LR Model"))
	######
	modules_interesting2 <- p_hd_extended_2
	q <- modules_interesting2[modules_interesting2$p_value <p_use,]
	modules_interesting2$p_value[modules_interesting2$p_value >= p_use] <- NA
	## Lets plot just these with a better scale 
	plot(ggplot( modules_interesting2 ) +
	geom_tile(aes(x=Module, y=Effect, fill = as.numeric(p_value)), colour="lightblue") + 
	scale_fill_continuous(type = "viridis", na.value="white", limits=c(0, p_use))+
	geom_point(data = q, aes(x=Module, y=Effect, shape = High, colour=Type), fill="White") + 
	scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease')) +
	scale_colour_discrete(labels = c("Health to Sepsis", "Timepoint 1 to 5")) +
	theme_classic() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	labs(fill="p value") + guides(alpha="none")+ xlab("Module") +ylab("LR Model"))
	dev.off()
	
	q <- p_hd_extended_2[p_hd_extended_2$BH <p_use,]
	pdf(paste0(outputdir,"/HealthvsDisease_AdjustedPvalue.pdf"), width=7, height=4)
	plot(ggplot( p_hd_extended_2 ) +
	geom_tile(aes(x=Module, y=Effect, fill = as.numeric(BH), alpha=ifelse(BH < p_use, 1, 0.9))) + 
	scale_fill_continuous(type = "viridis")+
	geom_point(data = q, aes(x=Module, y=Effect, shape = High, colour=Type), fill="White") + 
	scale_shape_manual(values=c(24,25), labels = c('Increase','Decrease')) +
	scale_colour_discrete(labels = c("Health to Sepsis", "Timepoint 1 to 5")) +
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
	scale_colour_discrete(labels = c("Health to Sepsis", "Timepoint 1 to 5")) +
	theme_classic() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	labs(fill="BH Adj p Value") + guides(alpha="none")+ xlab("Module") +ylab("LR Model"))
	dev.off()
	
	write.table(p_hd_extended_2, paste0(outputdir, "/Summary/Clinical_variables_hvsd.txt"), sep="\t")
	return(p_hd_extended_2)

}