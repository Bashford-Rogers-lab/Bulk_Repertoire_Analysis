## Function to visualise survival for sepsis patients
## Lauren Overend
## Jan 2023
library(cowplot)

plot_survival <- function(eigenvectors, metadata, outputdir){
	## Want to visualise survival for sepsis patients
	eigenvectors <- read.delim(eigenvectors, sep="\t", header=TRUE)
	##------------------------------------------------------------------------------
	### PART 1
	## Preparing metadata
	covariates <- c("Age", "Sex", "Comorbidities.Charlson_Age", "Comorbidities.Charlson_Index")
	metadata <- read.delim(metadata, sep="\t", header=TRUE)
	## Include for just samples that are in the dataset (after filtering etc)
	metadata <- metadata[metadata$SampleID_alternative %in% rownames(eigenvectors),]
	
	metadata_daysdeath <- metadata[,c("alternative_barcode", "Death.Days_from_ICU_admission", "Mortality.Classification")]
	metadata_daysdeath <- unique(metadata_daysdeath)	
	metadata_daysdeath_nonsurv <- 	metadata_daysdeath[!is.na(metadata_daysdeath$Death.Days_from_ICU_admission),]
	
	## calculate cumulative deaths 
	## threshold is 6 months 
	## how many days is this 
	x <- (6*30)
	
	## Lets calculate the number of deaths at each timepoint
	## Cumulative survival/deaths
	survival_cal <- c()
	for (i in 0:x){
		days_to_death <- i
		number_deaths <- metadata_daysdeath_nonsurv[metadata_daysdeath_nonsurv$Death.Days_from_ICU_admission==days_to_death,]
		number_deaths <- dim(number_deaths)[1]
		
		cum_deaths <-  metadata_daysdeath_nonsurv[metadata_daysdeath_nonsurv$Death.Days_from_ICU_admission <= days_to_death,]
		cum_deaths <- dim(cum_deaths)[1]
		
		row_use <- c(i, number_deaths, cum_deaths)
		survival_cal <- rbind(survival_cal, row_use)
	}
	survival_cal <- data.frame(survival_cal)
	colnames(survival_cal) <- c("Days_to_death", "Number_deaths", "Cumulative_deaths")
	survival_cal$proportion_surviving <- (dim(metadata_daysdeath)[1]-survival_cal$Cumulative_deaths)/(dim(metadata_daysdeath)[1])
	
	## Plot the survival curves
	pdf(paste0(outputdir,"/Survival_Curve.pdf"), width=6, height=3)
	x1 <- ggplot(survival_cal, aes(x=Days_to_death, y=proportion_surviving)) + geom_line() + theme_classic() +xlab("Days to Death from ICU Admission") + ylab("Proportion of Cohort Alive")+
	geom_vline(xintercept=7, col="red") + geom_vline(xintercept=180, col="blue") +ylim(0,1)
	x2 <- ggplot(survival_cal, aes(x=Days_to_death, y=Cumulative_deaths)) + geom_line() + theme_classic() +xlab("Days to Death from ICU Admission") + ylab("Cumulative Number of Deaths")+
	geom_vline(xintercept=7, col="red") + geom_vline(xintercept=180, col="blue") 
	plot(plot_grid(x1, x2, align="v", axis="lbt" , ncol=2))
	dev.off()

	## What is the best grouping....
	## Deaths under 7 days 
	day5 <- survival_cal$Cumulative_deaths[survival_cal$Days_to_death==5]
	day7 <- survival_cal$Cumulative_deaths[survival_cal$Days_to_death==7]
	day14 <- survival_cal$Cumulative_deaths[survival_cal$Days_to_death==14]
	day28 <- survival_cal$Cumulative_deaths[survival_cal$Days_to_death==28]
	day180 <- survival_cal$Cumulative_deaths[survival_cal$Days_to_death==180]
}