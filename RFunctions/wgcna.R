library("Hmisc")
library("splines")
library("foreach")
library("doParallel")
library("fastcluster")
library("dynamicTreeCut")
library("survival")
library("preprocessCore")
library("impute")
library("data.table")
library("dplyr")
library("purrr")
library("WGCNA")

## Essential for WGCNA analysis 
options(stringsAsFactors = FALSE);

## Setting up working directory :D 
.libPaths(c( "~/R/R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0", .libPaths()))
.libPaths(c("/well/immune-rep/shared/CODE/Communal_R_packages/R-bundle-Bioconductor_3.11_foss_2020a_R_4.0.0/", .libPaths()))
options(repos='http://cran.ma.imperial.ac.uk/')

setwd('/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/')
## Reading in the data 
bcr <- read.delim('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/BCR/Summary/isotyper_metrics_filtered_FINAL_METRICS_2000_PRODUCTIVE.txt')
tcra <- read.delim('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/TCRA/Summary/isotyper_metrics_filtered_FINAL_METRICS_500_PRODUCTIVE.txt')
tcrb<- read.delim('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/TCRB/Summary/isotyper_metrics_filtered_FINAL_METRICS_500_PRODUCTIVE.txt')
tcrg <- read.delim('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/TCRG/Summary/isotyper_metrics_filtered_FINAL_METRICS_750_PRODUCTIVE.txt')

## Setting the layout of the data
bcr$sample <- row.names(bcr)
tcrb$sample <- row.names(tcrb)
tcra$sample <- row.names(tcra)
tcrg$sample <- row.names(tcrg)

## Merging into one big data frame 
data <- list(bcr, tcrb, tcrg, tcra) %>% reduce(inner_join, by="sample")
rownames(data) <- data$sample
data$Sample <- NULL
data <- data.frame(data)


## FOLLOWING WGCNA TUTORIAL and preparing data into correct format 
datExpr0 <- data.frame(sapply(data, as.numeric ))
names(datExpr0) <- colnames(data)
rownames(datExpr0) <- row.names(data)
rownames(datExpr0) <- gsub("_productive", "", rownames(datExpr0))


## Assess normal distribution of data 
do.call(rbind, lapply(datExpr0, function(x) shapiro.test(x)[c("statistic", "p.value")]))

## Want to plot and see if normally distributed 
pdf("trydit.pdf")
for(s in colnames(datExpr0)){
	#print(s)			
	d <- ggplot(data[data[, s] >-1,], aes_string(x=s))  + geom_histogram() + theme_bw() 
	tryCatch({
	plot(d)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
}


datExpr02 <- data.frame(apply(datExpr0, 2, function(x){log(x)}))

## Want to plot and see if normally distributed 
pdf("tryditlog.pdf")
for(s in colnames(datExpr02)){
	#print(s)			
	d <- ggplot(datExpr02[datExpr02[, s] >-1,], aes_string(x=s))  + geom_histogram() + theme_bw() 
	tryCatch({
	plot(d)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})	
}


gsg = goodSamplesGenes(datExpr0, verbose=3)


## Removing genes that did not pass filtering 
if (!gsg$allOK)
{
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes: ", paste0(names(datExpr0)[!gsg$goodGenes])));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

###
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
dev.off()

## metadata 


files_use <- list.files('/gpfs2/well/immune-rep/users/kvi236/GAinS_Data')
files_use <- grep("txt", files_use, value=TRUE

b <- read.delim('/gpfs2/well/immune-rep/users/kvi236/GAinS_Data/CLIN_DATA_days1_3_5.txt')

datTraits <- data.frame(sapply(b, as.numeric ))
names(datTraits) <- colnames(b)
row.names(datTraits) <- b$SAMPLE_ID

datTraits$SAMPLE_ID <- NULL
datTraits$SubjectNumber <- NULL
datTraits$study <- NULL
datTraits$DhospICU <- NULL
datTraits$Birthdate <- NULL
datTraits$SubjectBarCode <- NULL
datTraits$DAY <- NULL
datTraits$Diagnosis <- NULL
datTraits$CenterNumber <- NULL
datTraits$Patient_type <- NULL

rownames(datTraits)[rownames(datTraits)=="GAUKJR010000_1"] <- "UKJR010000_1"
rownames(datTraits)[rownames(datTraits)=="GAUKJR010000_3"] <- "UKJR010000_3"
rownames(datTraits)[rownames(datTraits)=="GAUKJR010000_5"] <- "UKJR010000_5"

rownames(datTraits)[rownames(datTraits)=="uk01310070_3"] <- "UK01310070_3"
rownames(datTraits)[rownames(datTraits)=="uk01310070_5"] <- "UK01310070_5"


datExpr0 <- datExpr0[!rownames(datExpr0) %in% c(grep("JR1795", rownames(datExpr0), value=TRUE)),]
datTraits <- datTraits[row.names(datTraits) %in% rownames(datExpr0),]


## remove rows where values are constant
nz_cols2 = apply(datTraits, 2, function(x){length(unique(x[!is.na(x)]))})
datTraits2 = datTraits[,which(nz_cols2 > 1)]

pdf(file = "sampleClustering2.pdf", width = 20, height = 15);

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits2)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits2), main = "Sample dendrogram and trait heatmap")
dev.off()

## Part 2
enableWGCNAThreads()


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:

pdf(file = "sampleClustering3.pdf", width = 20, height = 15);

par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


net = blockwiseModules(datExpr0, power = 10,
TOMType = "unsigned", minModuleSize = 3, deepSplit=1,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "femaleMouseTOM",
verbose = 3)


# open a graphics window
pdf(file = "sampleClustering4.pdf", width = 20, height = 15);

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)




