library(stringr)
library(ggforce)
library(ggpubr)
### Summary funciton 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

make_matrices14_normalised <- function(file=file, chain_vdj=chain_vdj, ids_all=ids_all, counts_used=counts_used, iso_type=iso_type, path_to_layout=path_to_layout, outputdir=outputdir){
	info = file.info(file)
	p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
	p=p[which(as.character(p[,"X.sample"]) %in% ids_all),]
	p <- data.frame(p)
	p$uniq_read_freq <- as.numeric(p$uniq_read_freq)		
	#####
	isotypes <- c("IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHEP2", "IGHGP", "IGHE")
	other_class <- c("Class_switched", "IGHD,IGHM_mutated", "IGHD,IGHM_unmutated", "IGHD,IGHM_unmutated_singleton")
	expansion_class <- c("unexpanded", "expanded")
	####	
	if(chain_vdj %like% "BC"| chain_vdj %like% "I"){
		all_class <- c("ALL")
	} else {
		counts_try <- counts_used[(counts_used$X.isotype  != "ALL" & counts_used$X.isotype  != "all"),]
		receptor_type <- counts_try$X.isotype[counts_try$min==max(counts_try$min)]
		all_class <- receptor_type
	}
	## Total V gene usage 
	#####
	p_all <- p[p$class %in% all_class,]
	p_iso <-p[p$class %in% isotypes,] 
			
	## Calculate a percentage of repertoire which is each read 
	p_all$percent_repertoire <- NA
	for(i in unique(p_all[, "X.sample"])){
		sample_id <- i 
		sum_frequency <- sum(p_all$uniq_read_freq[p_all$X.sample==sample_id])
		p_all$percent_repertoire[p_all$X.sample==sample_id] <- ((p_all$uniq_read_freq[p_all$X.sample==sample_id])/sum_frequency)*100
	} 
	p_allx <- p_all[, c("X.sample", "percent_repertoire", "V.gene")]

	#### Read in layouts file 
	## Correct using technical 
	layouts <- read.delim(path_to_layout)
	layouts$Lane <- as.character(layouts$Lane)
	### Give them a brief annotation if they are my sepsis cohorts
	if(outputdir %like% "SEPSIS"){
		layouts$Lane[layouts$Lane != "13" & layouts$Lane != "12"& layouts$Lane != "11"& layouts$Lane != "10" & layouts$Lane != "9"] <- paste0(layouts$Lane[layouts$Lane != "13" & layouts$Lane != "12"& layouts$Lane != "11"& layouts$Lane != "10" & layouts$Lane != "9"], "_CH1")
		layouts$Lane[layouts$Lane == "13" | layouts$Lane == "12"| layouts$Lane == "11"] <- paste0(layouts$Lane[layouts$Lane == "13" | layouts$Lane == "12"| layouts$Lane == "11"], "_HEALTH")
		layouts$Lane[layouts$Lane == "10" | layouts$Lane == "9"] <- paste0(layouts$Lane[layouts$Lane == "10" | layouts$Lane == "9"], "_CH2_ED")

	}
	
	if(iso_type =="PRODUCTIVE"){
		layouts$SampleID <- paste0(layouts$SampleID, "_productive")
	} else if (iso_type =="UNPRODUCTIVE"){
		layouts$SampleID <- paste0(layouts$SampleID, "_unproductive")
	}
	
	technical = "JR1"
	p_allx2 <- merge(layouts, p_allx, by.x="SampleID", by.y="X.sample")
	techs <- p_allx2[p_allx2$SampleID %like% "POSITIVE",]
	#techs2 <- techs %>%group_by(Lane, V.gene) %>%dplyr::summarize(Mean = mean(percent_repertoire, na.rm=TRUE))
	techs2 <- summarySE(techs, measurevar="percent_repertoire", groupvars=c("Lane","V.gene"))
	colnames(techs2)[4] <- "Mean"
	techs2 <- data.frame(techs2)
	
	### We also want to do a plot correlating technical to lane to see if theres a trend 
    s <- p_allx2
	s$Class <- "SAMPLE" 
	s$Class[s$SampleID %like% "POSITIVE"] <- "TECHNICAL"
	s2 <- s %>%group_by(Lane, V.gene,Class) %>%dplyr::summarize(Mean = mean(percent_repertoire, na.rm=TRUE))
	s3 <- spread(s2, key = Class, value = Mean, fill=0)
	
	pd <- position_dodge(0.3) 
	pdf(paste0(outputdir, "/VGene_Usage_technicals_", iso_type, ".pdf"), width=15, height=10)
	p1x <- ggplot(techs2, aes(x=V.gene, y=Mean, fill=as.factor(Lane), colour=as.factor(Lane)))+geom_point(alpha=0.9, position=pd) +theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(color="Lane", fill="Lane")+ylab("Mean % Repertoire")+ggtitle("Technical Replicates") +  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=.1, position=pd, alpha=0.5)+
     facet_zoom(ylim = c(0, 6), zoom.size=1) + xlab("V Gene")
	p3 <- ggplot(s3, aes(x=TECHNICAL, y=SAMPLE))+geom_point(alpha=0.7, aes(colour=as.factor(V.gene), shape=as.factor(Lane))) +theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(color="Lane", fill="Lane")+ylab("Mean % Repertoire")+ggtitle("Correlation between Technical V Gene Usage and Average Lane (sample) V Gene Usage ") +
     ylab("Average % Repertoire for Samples") +xlab("Average % Repertoire for Technicals") +labs(shape="Lane") +guides(color="none")+ stat_cor(method = "pearson") + geom_smooth(method='lm', formula= y~x)
	#p2 <-  ggplot(techs, aes(x=V.gene, y=percent_repertoire, colour=as.factor(Lane)))+geom_boxplot(alpha=0.3, aes(fill=as.factor(Lane))) +theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(color="Lane", fill="Lane")+ylab("% Repertoire")+ggtitle("Technical Replicates")
    #p3 <- ggplot(techs2, aes(x=V.gene, y=log2(Mean), fill=as.factor(Lane), colour=as.factor(Lane)))+geom_point() +theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(color="Lane", fill="Lane")+ylab("log2(Mean % Repertoire)")+ggtitle("Technical Replicates")
	plot(plot_grid(p1x, p3, ncol=1))
	dev.off()
	
	p_allx2$NormalisedPercent <- NA
	for(i in 1:length(p_allx2$SampleID)){
		lane <- p_allx2$Lane[i]
		gene <- p_allx2$V.gene[i]
		val <- p_allx2$percent_repertoire[i]-techs2$Mean[techs2$Lane==lane & techs2$V.gene ==gene]
		if(length(val)==0) {
			val = NA
		}
		p_allx2$NormalisedPercent[i] <- val 
	}
	
	pdf(paste0(outputdir, "/VGene_Usage_normalised_", iso_type, ".pdf"), width=15, height=15)
	p1 <- ggplot(p_allx2, aes(x=V.gene, y=percent_repertoire, fill=as.factor(Lane), colour=as.factor(Lane)))+geom_boxplot(alpha=0.5) +theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Raw Percentage")+labs(colour="Lane", fill="Lane")+ggtitle("% of Repertoire for All Samples")+ xlab("V Gene")
	p2 <- ggplot(p_allx2, aes(x=V.gene, y=NormalisedPercent, fill=as.factor(Lane),colour=as.factor(Lane)))+geom_boxplot(alpha=0.5) +theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Normalised Percentage")+labs(colour="Lane", fill="Lane")+ggtitle("Normalised % of Repertoire for All Samples")+ xlab("V Gene")
	plot(plot_grid(p1x, p3, p1, p2, ncol=1))
	dev.off()
	
	p_allx2 <- p_allx2[, c("SampleID", "V.gene", "NormalisedPercent")]
    colnames(p_allx2)[3] <- "percent_repertoire"
	
	## We fill in any missing combinations with 0 as 0 percent of repertoire 
	a <- spread(p_allx2, key = V.gene, value = percent_repertoire, fill=0)
	rownames(a) <- a$SampleID
	a$SampleID <- NULL
	colnames(a) <- gsub("-", "_", colnames(a))
	colnames(a) <- paste0(colnames(a), "__", all_class)
	a <- as.matrix(a)
	colnames(a) <- paste0("NORM_", colnames(a))
	### Look at V gene Family per isotype 
	### Only relevant for BCRs!!!!
	
	if(chain_vdj %like% "BC"| chain_vdj %like% "I"){
		p$family <- str_split_fixed(p$V.gene, "-", 2)[,1]
		p1 <- p[p$class %in% isotypes,] 
		q <- aggregate(p1$uniq_read_freq, list(p1$X.sample, p1$family, p1$class), FUN=sum) 
		colnames(q) <- c("Sample", "V.Family", "isotype", "Unique_Frequency")
		
		q$percent_repertoire <- NA
		for(i in unique(q[, "Sample"])){
			sample_id <- i 
			sum_frequency <- sum(q$Unique_Frequency[q$Sample==sample_id])
			q$percent_repertoire[q$Sample==sample_id] <- ((q$Unique_Frequency[q$Sample==sample_id])/sum_frequency)*100
		} 
		q <- q[, c("Sample", "V.Family", "isotype", "percent_repertoire")]
		
		## Now to normalise 
		q2 <- merge(layouts, q, by.x="SampleID", by.y="Sample")
		techs <- q2[q2$SampleID %like% "POSITIVE",]
		techs2 <- techs %>%group_by(Lane, V.Family, isotype) %>%dplyr::summarize(Mean = mean(percent_repertoire, na.rm=TRUE))
		techs2 <- data.frame(techs2)
		
		pdf(paste0(outputdir, "/VGene_Usage_perISOTYPE_technicals_", iso_type, ".pdf"), width=10, height=10)
		p1x <- ggplot(techs2, aes(x=V.Family, y=Mean, fill=as.factor(Lane), colour=as.factor(Lane)))+geom_point() +theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(color="Lane", fill="Lane")+ylab("Mean % Repertoire")+ggtitle("Technical Replicates")+facet_wrap(~isotype, ncol=11)
		p2 <-  ggplot(techs, aes(x=V.Family, y=percent_repertoire, colour=as.factor(Lane)))+geom_boxplot(alpha=0.3, aes(fill=as.factor(Lane))) +theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(color="Lane", fill="Lane")+ylab("% Repertoire")+ggtitle("Technical Replicates")+facet_wrap(~isotype, ncol=11)
		p3 <- ggplot(techs2, aes(x=V.Family, y=log2(Mean), fill=as.factor(Lane), colour=as.factor(Lane)))+geom_point() +theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(color="Lane", fill="Lane")+ylab("log2(Mean % Repertoire)")+ggtitle("Technical Replicates")+facet_wrap(~isotype, ncol=11)
		plot_grid(p1x, p2, p3, ncol=1)
		dev.off()
		
		## Normalise 
		q2$NormalisedPercent <- NA
		for(i in 1:length(q2$SampleID)){
			lane <- q2$Lane[i]
			gene <- q2$V.Family[i]
			iso <- q2$isotype[i]
			val <- q2$percent_repertoire[i]-techs2$Mean[techs2$Lane==lane & techs2$V.Family ==gene & techs2$isotype ==iso]
			if(length(val)==0) {
				val = NA
		}
			q2$NormalisedPercent[i] <- val 
		}
		
		pdf(paste0(outputdir, "/VGene_Usage_normalised_periso_", iso_type, ".pdf"), width=20, height=10)
		p1 <- ggplot(q2, aes(x=V.Family, y=percent_repertoire, fill=as.factor(Lane), colour=as.factor(Lane)))+geom_boxplot(alpha=0.5) +theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Raw Percentage")+labs(colour="Lane", fill="Lane")+ggtitle("Raw % of Repertoire for All Samples")+facet_wrap(~isotype, ncol=11)
		p2 <- ggplot(q2, aes(x=V.Family, y=NormalisedPercent, fill=as.factor(Lane),colour=as.factor(Lane)))+geom_boxplot(alpha=0.5) +theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Normalised Percentage")+labs(colour="Lane", fill="Lane")+ggtitle("Normalised % of Repertoire for All Samples")+facet_wrap(~isotype, ncol=11)
		plot_grid(p1x, p1, p2, ncol=1)
		dev.off()
		
		## Reformat 
		q2x2 <- q2[, c("SampleID", "V.Family", "isotype", "NormalisedPercent")]
		colnames(q2x2)[4] <- "percent_repertoire"
		
		q <- q2x2
		q$type <- paste0("NORM_", q$V.Family, "__", q$isotype)
		q <- q[, c("SampleID", "type", "percent_repertoire")]
		b <- spread(q, key = type, value = percent_repertoire, fill=0)
		rownames(b) <- b$SampleID
		b$SampleID <- NULL
		b <- as.matrix(b)
		a <- merge(a, b, by=0)
		rownames(a) <- a$Row.names
		a$Row.names <- NULL
	}	
	analysis_matrices14 = list(a)
	names(analysis_matrices14) <- "V_GENE_USAGE"
	## Save the dataframe!
	write.table(a, paste0(outputdir, "Summary/V_Gene_usage_normalised_", iso_type, ".txt"), sep="\t", row.names=TRUE)
	print("DONE 14: V GENE Usages")
	return(analysis_matrices14)	
} 

