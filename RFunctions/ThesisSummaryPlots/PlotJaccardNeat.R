## Functions to plot the results of the Jaccard Analysis (batch effect)
## Lauren Overend

library(ggpubr)
library(rstatix) 
library(BBmisc)
library(gtools) 
library(coin)

## Multiplot Function
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)

  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                 ncol = cols, nrow = ceiling(numPlots/cols))
}

if (numPlots == 1) {
print(plots[[1]])

} else {
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

for (i in 1:numPlots) {
  matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

  print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                  layout.pos.col = matchidx$col))
 }
}
 }
 
#jaccard_matrix <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/JACCARDMATRIX_BASIC_AllSAMPLES_LEO_SEPSIS.txt'
#path_to_output <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/'
#layouts1 <- '/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/LEO_SEPSIS_BCR_ALL_layouts.txt'

create_jaccard_plots <- function(path_to_output, jaccard_matrix, layouts1, chain){
	layouts <- read.delim(layouts1, sep='\t')
	JACCARD_MATRIX_2 <- read.delim(jaccard_matrix, sep='\t')
	
	JACCARD_MATRIX_2$Sample1 <- gsub("BCR_", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub("BCR_", "", JACCARD_MATRIX_2$Sample2)
	JACCARD_MATRIX_2$Sample1 <- gsub("TCRA_", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub("TCRA_", "", JACCARD_MATRIX_2$Sample2)
	JACCARD_MATRIX_2$Sample1 <- gsub("TCRB_", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub("TCRB_", "", JACCARD_MATRIX_2$Sample2)
	JACCARD_MATRIX_2$Sample1 <- gsub("TCRG_", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub("TCRG_", "", JACCARD_MATRIX_2$Sample2)
	JACCARD_MATRIX_2$Sample1 <- gsub("TCRD_", "", JACCARD_MATRIX_2$Sample1)
	JACCARD_MATRIX_2$Sample2 <- gsub("TCRD_", "", JACCARD_MATRIX_2$Sample2)
	
	## Merge Jaccard Matrix by Layouts 
	layouts_1 <- layouts
	colnames(layouts_1) <- c("SampleID", "Barcode1",  "Lane_S1", "Plate_S1", "Library1", "Position1", "PCRBarcode1")
	JACCARD_MATRIX_2_cohort1 <- merge(JACCARD_MATRIX_2, layouts_1, by.x="Sample1", by.y="SampleID")

	# IDS for sample S2:
	layouts_2 <- layouts
	colnames(layouts_2) <- c("SampleID", "Barcode2",  "Lane_S2", "Plate_S2", "Library2", "Position2", "PCRBarcode2")
	JACCARD_MATRIX_2_cohort1 <- merge(JACCARD_MATRIX_2_cohort1, layouts_2, by.x="Sample2", by.y="SampleID")
	
	## Columns of Sharing PCR Barcode: 

	
	### Are they the same sample
	JACCARD_MATRIX_2_cohort1$SharedSample <- NA
	JACCARD_MATRIX_2_cohort1$SharedSample[JACCARD_MATRIX_2_cohort1$Sample1 == JACCARD_MATRIX_2_cohort1$Sample2 ] <- "Same RNA Sample"
	JACCARD_MATRIX_2_cohort1$SharedSample[JACCARD_MATRIX_2_cohort1$Sample1 != JACCARD_MATRIX_2_cohort1$Sample2 ] <- "Different RNA Sample"

	### Are they the same Individual 
	JACCARD_MATRIX_2_cohort1$SharedIndividual <- NA
	JACCARD_MATRIX_2_cohort1$SharedIndividual[JACCARD_MATRIX_2_cohort1$Barcode1 == JACCARD_MATRIX_2_cohort1$Barcode2 ] <- "Y"
	JACCARD_MATRIX_2_cohort1$SharedIndividual[JACCARD_MATRIX_2_cohort1$Barcode1 != JACCARD_MATRIX_2_cohort1$Barcode2 ] <- "N"

	### Are they the same PCR Barcode 
	JACCARD_MATRIX_2_cohort1$SharedPCRBarcode <- NA
	JACCARD_MATRIX_2_cohort1$SharedPCRBarcode[JACCARD_MATRIX_2_cohort1$PCRBarcode1 == JACCARD_MATRIX_2_cohort1$PCRBarcode2 ] <- "Y"
	JACCARD_MATRIX_2_cohort1$SharedPCRBarcode[JACCARD_MATRIX_2_cohort1$PCRBarcode1 != JACCARD_MATRIX_2_cohort1$PCRBarcode2 ] <- "N"
	
	### Are they the same Lane  
	JACCARD_MATRIX_2_cohort1$SharedLane <- NA
	JACCARD_MATRIX_2_cohort1$SharedLane[JACCARD_MATRIX_2_cohort1$Lane_S1 == JACCARD_MATRIX_2_cohort1$Lane_S2 ] <- "Same\nSequencing Lane"
	JACCARD_MATRIX_2_cohort1$SharedLane[JACCARD_MATRIX_2_cohort1$Lane_S1 != JACCARD_MATRIX_2_cohort1$Lane_S2 ] <- "Different\nSequencing Lane"
	
	## What samples are we comparing 
	JACCARD_MATRIX_2_cohort1$Comparison <- NA 
	JACCARD_MATRIX_2_cohort1$Sample1T <- "S"
	JACCARD_MATRIX_2_cohort1$Sample1T[JACCARD_MATRIX_2_cohort1$Sample1 %like% "HV"] <- "H"
	JACCARD_MATRIX_2_cohort1$Sample1T[JACCARD_MATRIX_2_cohort1$Sample1 %like% "JR1795_1003"] <- "T"
	JACCARD_MATRIX_2_cohort1$Sample2T <- "S"
	JACCARD_MATRIX_2_cohort1$Sample2T[JACCARD_MATRIX_2_cohort1$Sample2 %like% "HV"] <- "H"
	JACCARD_MATRIX_2_cohort1$Sample2T[JACCARD_MATRIX_2_cohort1$Sample2 %like% "JR1795_1003"] <- "T"
	
	JACCARD_MATRIX_2_cohort1$Comparison[JACCARD_MATRIX_2_cohort1$Sample1T=="H" & JACCARD_MATRIX_2_cohort1$Sample2T=="H"] <- "H vs H"
	JACCARD_MATRIX_2_cohort1$Comparison[JACCARD_MATRIX_2_cohort1$Sample1T=="S" & JACCARD_MATRIX_2_cohort1$Sample2T=="S"] <- "S vs S"
	JACCARD_MATRIX_2_cohort1$Comparison[JACCARD_MATRIX_2_cohort1$Sample1T=="S" & JACCARD_MATRIX_2_cohort1$Sample2T=="H"] <- "H vs S"
	JACCARD_MATRIX_2_cohort1$Comparison[JACCARD_MATRIX_2_cohort1$Sample1T=="H" & JACCARD_MATRIX_2_cohort1$Sample2T=="S"] <- "H vs S"
	
	## Make a summary Column incoperating the different categories which will be used for plotting
	JACCARD_MATRIX_2_cohort1$Summary <- paste0(JACCARD_MATRIX_2_cohort1$SharedInternalBarcode, ".", JACCARD_MATRIX_2_cohort1$SharedLane)
	JACCARD_MATRIX_2_cohort1$ANYNAS <- as.factor(JACCARD_MATRIX_2_cohort1$ANYNAS)
	
	vars_to_plot1 <- c("log.Jaccard.MeanSubsample", "log.Jaccard.MeanSubsample.Weighted", "log.Jaccard.Full", "log.Jaccard.Full.Weighted")
	# Make plots 
	
	c <- min(JACCARD_MATRIX_2_cohort1$Jaccard.MeanSubsample[JACCARD_MATRIX_2_cohort1$Jaccard.MeanSubsample!=0 & !is.na(JACCARD_MATRIX_2_cohort1$Jaccard.MeanSubsample)])/2
	c2 <- min(JACCARD_MATRIX_2_cohort1$Jaccard.MeanSubsample.Weighted[JACCARD_MATRIX_2_cohort1$Jaccard.MeanSubsample.Weighted!=0 & !is.na(JACCARD_MATRIX_2_cohort1$Jaccard.MeanSubsample.Weighted)])/2
	c3 <- min(JACCARD_MATRIX_2_cohort1$Jaccard.Full[JACCARD_MATRIX_2_cohort1$Jaccard.Full!=0 & !is.na(JACCARD_MATRIX_2_cohort1$Jaccard.Full)])/2
	c4 <- min(JACCARD_MATRIX_2_cohort1$Jaccard.Full.Weighted[JACCARD_MATRIX_2_cohort1$Jaccard.Full.Weighted!=0 & !is.na(JACCARD_MATRIX_2_cohort1$Jaccard.Full.Weighted)])/2

	
	JACCARD_MATRIX_2_cohort1$Jaccard.MeanSubsample_LOG <- log2(JACCARD_MATRIX_2_cohort1$Jaccard.MeanSubsample+c)
	JACCARD_MATRIX_2_cohort1$Jaccard.MeanSubsample.Weighted_LOG <- log2(JACCARD_MATRIX_2_cohort1$Jaccard.MeanSubsample.Weighted+c2)
	JACCARD_MATRIX_2_cohort1$Jaccard.Full_LOG <- log2(JACCARD_MATRIX_2_cohort1$Jaccard.Full+c3)
	JACCARD_MATRIX_2_cohort1$Jaccard.Full.Weighted_LOG <- log2(JACCARD_MATRIX_2_cohort1$Jaccard.Full.Weighted+c4)
	
	vars_to_plot2 <- c("Jaccard.MeanSubsample_LOG", "Jaccard.MeanSubsample.Weighted_LOG", "Jaccard.Full_LOG", "Jaccard.Full_LOG")

	## lets get the counts in each group 
	JACCARD_MATRIX_2_cohort_X <- JACCARD_MATRIX_2_cohort1[!is.na(JACCARD_MATRIX_2_cohort1$Comparison),]
	counts <- JACCARD_MATRIX_2_cohort_X %>% dplyr::group_by(SharedSample, SharedIndividual, SharedPCRBarcode, SharedLane,Comparison ) %>% dplyr::summarise(n = n()) 
	counts2 <- JACCARD_MATRIX_2_cohort_X %>% dplyr::group_by(SharedSample, SharedIndividual, SharedPCRBarcode, SharedLane) %>% dplyr::summarise(n = n()) 
	
	give.n <- function(x){
	  return(c(y = median(x)*1.05, label = length(x))) 
	  # experiment with the multiplier to find the perfect position
	}

	## There will be some NAs if there wasnt enough read depth to make the comparison!!!!!!!!!!
	pdf(paste0(path_to_output, "Plots/NeatJaccard.pdf"), width=10, height=5)
	p1 <- ggplot(JACCARD_MATRIX_2_cohort_X, aes_string(x="SharedPCRBarcode", y=vars_to_plot2[1], fill="SharedIndividual")) +geom_boxplot(aes(fill=SharedIndividual, colour=SharedIndividual), alpha=0.5) +facet_grid(Comparison~SharedSample+SharedLane)+theme_classic()+xlab("Shared Internal PCR Barcode") +ylab("Log2(Subsampled Jaccard Index +C)")+
	scale_color_manual(values=c("#F8766D", "#619CFF"))+scale_fill_manual(values=c("#F8766D", "#619CFF"))  + stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75), vjust = -1, colour="black") +labs(colour="Same\nIndividual", fill="Same\nIndividual")+ggtitle(chain) +
	scale_y_continuous(expand=c(1.5, 1.5))
	plot(p1)
	dev.off()
	
	#p2 <- ggplot(counts, aes(x=SharedPCRBarcode, y=n, fill=SharedIndividual)) +geom_col(stat="identity", position=position_dodge()) +facet_grid(Comparison~SharedSample+SharedLane)+theme_classic()+labs(fill="Same\nIndividual")+xlab("Shared Internal PCR Barcode") +ylab("Number of Pairwise Comparisons")+geom_text(aes(label=n), vjust=-1, position = position_dodge(.9))
	#plot(plot_grid(p1, p2, ncol=1,  align="hv", axis="tblr", labels="AUTO"))
	
	#p1 <- ggplot(JACCARD_MATRIX_2_cohort_X, aes_string(x="SharedPCRBarcode", y=vars_to_plot2[1])) +geom_boxplot(aes(fill=SharedIndividual)) +facet_wrap(~SharedSample+SharedLane)+theme_classic()+labs(fill="Same\nIndividual")+xlab("Shared Internal PCR Barcode") +ylab("Log2(Subsampled Jaccard Index +C)")
	#p2 <- ggplot(counts2, aes(x=SharedPCRBarcode, y=n, fill=SharedIndividual)) +geom_col(stat="identity", position=position_dodge())  +facet_wrap(~SharedSample+SharedLane)+theme_classic()+labs(fill="Same\nIndividual")+xlab("Shared Internal PCR Barcode") +ylab("Number of Pairwise Comparisons")+geom_text(aes(label=n), vjust=-1, position = position_dodge(.9))
	#plot(plot_grid(p1, p2, ncol=1,  align="hv", axis="tblr", labels="AUTO"))
	
	
	#p1 <- ggplot(JACCARD_MATRIX_2_cohort_X[JACCARD_MATRIX_2_cohort_X$SharedSample!="Same RNA Sample",], aes_string(x="SharedPCRBarcode", y="Jaccard.MeanSubsample")) +geom_boxplot(aes(colour=SharedIndividual)) +facet_grid(Comparison~SharedSample+SharedLane, scales="free")+theme_classic()+labs(colour="Same\nIndividual")+xlab("Shared Internal PCR Barcode") +ylab("Log2(Subsampled Jaccard Index +C)")
	#p11 <- ggplot(JACCARD_MATRIX_2_cohort_X[JACCARD_MATRIX_2_cohort_X$SharedSample=="Same RNA Sample",], aes_string(x="SharedPCRBarcode", y="Jaccard.MeanSubsample")) +geom_boxplot(aes(colour=SharedIndividual)) +facet_grid(Comparison~SharedSample+SharedLane, scales="free")+theme_classic()+labs(colour="Same\nIndividual")+xlab("Shared Internal PCR Barcode") +ylab("Log2(Subsampled Jaccard Index +C)")
	#p2 <- ggplot(counts, aes(x=SharedPCRBarcode, y=n, fill=SharedIndividual)) +geom_col(stat="identity", position=position_dodge()) +facet_grid(Comparison~SharedSample+SharedLane)+theme_classic()+labs(fill="Same\nIndividual")+xlab("Shared Internal PCR Barcode") +ylab("Number of Pairwise Comparisons")+geom_text(aes(label=n), vjust=-1, position = position_dodge(.9))
	#plot(plot_grid(p1, p2, ncol=1,  align="hv", axis="tblr"))
	#dev.off()
	
	#pdf(paste0(path_to_output, "Plots/NeatJaccard1.pdf"), width=12, height=12)
	#p1 <- ggplot(JACCARD_MATRIX_2_cohort_X[JACCARD_MATRIX_2_cohort_X$SharedSample!="Same RNA Sample",], aes_string(x="SharedPCRBarcode", y="Jaccard.MeanSubsample")) +geom_boxplot(aes(colour=SharedIndividual)) +facet_wrap(~SharedLane)+theme_classic()+labs(colour="Same\nIndividual")+xlab("Shared Internal PCR Barcode") +ylab("Subsampled Jaccard Index")+ggtitle("A. Different RNA Sample\n")
	#p11 <- ggplot(JACCARD_MATRIX_2_cohort_X[JACCARD_MATRIX_2_cohort_X$SharedSample=="Same RNA Sample",], aes_string(x="SharedPCRBarcode", y="Jaccard.MeanSubsample")) +geom_boxplot(aes(colour=SharedIndividual)) +facet_wrap(~SharedLane)+theme_classic()+labs(colour="Same\nIndividual")+xlab("Shared Internal PCR Barcode") +ylab("Subsampled Jaccard Index")+ggtitle("B. Same RNA Sample\nReference Scale")
	#p2 <- ggplot(counts2, aes(x=SharedPCRBarcode, y=n, fill=SharedIndividual)) +geom_col(stat="identity", position=position_dodge())  +facet_wrap(~SharedSample+SharedLane)+theme_classic()+labs(fill="Same\nIndividual")+xlab("Shared Internal PCR Barcode") +ylab("Number of Pairwise Comparisons")+geom_text(aes(label=n), vjust = -0.5, position = position_dodge(.9))+ggtitle("C.")
	#plot(plot_grid(plot_grid(p1,p11, ncol=2, rel_widths = c(2,1),  align="hv", axis="tblr"), p2, ncol=1,  align="hv", axis="tblr"))
	#dev.off()
	return(p1)

}










