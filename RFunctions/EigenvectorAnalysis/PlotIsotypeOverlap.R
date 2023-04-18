##### CODE to take the output of isotyper and then get the isotype frequencies and plot 
concat = function(v) {
  res = ""
  for (i in 1:length(v)){
    res = paste(res,v[i],sep="")
    }
  res
}
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
	apply(sapply(col, col2rgb)/255, 2, function(x) 
	rgb(x[1], x[2], x[3], alpha=alpha))
	}
## PACKAGES
library(igraph)
library(dplyr)
library(stringr)
library(data.table)

				   
#################################################
### Output BCR directory
out_dir <- '/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR'
## Order of Class switching 
file = "/gpfs3/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/EigenvectorAnalysis/Transition_Class_Genetic_Order.txt"
p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
classes = as.character(p[,2])
classes = unique(classes)

#################################################
#################################################
### FIRST PLOT an empty CSR plot with directions and no weights 
g <- graph.empty(n=0, directed=TRUE)
g <- igraph::add.vertices(g, length(classes), name= classes,color = "lightgreen")
names <- V(g)$name
ids <- 1:length(names)
names(ids) <- names
g1 = g
from = NULL
to = NULL
edge_strength = NULL
via = NULL
for(i in c(1:length(classes))){
	for(j in c(i:length(classes))){
		if(i<j){
			from = c(from, classes[i])
			to = c(to, classes[j])
			edge_strength = c(edge_strength, 1)
			via = c(via, concat(c(classes[i],":", classes[j] )))
}}}

edges <- matrix(c(ids[from], ids[to]), nc=2)
edge_ids = cbind(from, to)
g <- add.edges(g, t(edges), weight= edge_strength)
V(g)$label <- V(g)$name
layout1 =layout_in_circle(g)
V(g)$size<-50
cols = add.alpha("black",alpha = 0.5)
E(g)$color =cols
	
pdf(paste0(out_dir,"/CSR_Transition_healthy_frequencies.pdf"), height=5, width=5)
par(mar = c(0,0,2.5,0))
plot(g, layout = layout1, edge.color = cols[[1]], main = "Class Switch Recombination", 
     edge.width = edge_strength, vertex.label.family = "sans", 
     vertex.label.cex = 1, edge.arrow.size = 1, edge.lty = 1, 
     vertex.color = "lightblue", vertex.label.color = "black", edge.alpha=0.5)
dev.off()

#################################################
#################################################
## LETS GENERATE CSR FOR OBSERVED FREQUENCIES!!!!
## This is the imputed file with no scaling!!!!!
file =paste0(out_dir, "/Imputed_DATA_FINAL_BCR_PRODUCTIVE.txt")
p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))

## We had two different ways of normalising pick which one you want!
iso_switch <- colnames(p)[colnames(p) %like% "Isotype_subsampled_50_overlap"]
type_use <- "Subsampled50_Overlap"
#iso_switch <- colnames(p)[colnames(p) %like% "Isotype_normalised_overlap_frequencies"]
#type_use <- "Normalised_Overlap"

p=p[, c(iso_switch)]
headers = colnames(p)
classes = gsub("BCR_READS_Isotype_normalised_overlap_frequencies__","",headers)
classes = gsub("BCR_READS_Isotype_subsampled_50_overlap_frequencies__","",classes)
classes = gsub("IGHD_M","IGHD/M", classes,fixed = T)

## Need to read in metadata and take a mean 
meta1 <- read.delim("/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/sepsis_meta_health.txt", sep="\t")
rownames(p) <- gsub("_productive", "", rownames(p))
rownames(meta1) <- meta1$SampleID_alternative
## Only ones missing are the technical!! Good good 
p2 <- merge(p, meta1, by=0) 
## We also want just those in final eigenvector 
eigenvectors <- read.delim('/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Eigenvectors_No_Technical_BCR_PRODUCTIVE.txt')
bad_ids <- c("UK02870104_5", "GAUKRV025000_3")
eigenvectors <- eigenvectors[!eigenvectors$sample %in% bad_ids,]
keeps <- p2$SampleID_alternative[p2$SampleID_alternative %in% rownames(eigenvectors)]	
p3 <- p2[p2$SampleID %in% keeps,]

##########################################################################################################################
##########################################################################################################################
p3$DAY  <- 1
p3$DAY[grep("_3", p3$SampleID_alternative)] <- 3
p3$DAY[grep("_5", p3$SampleID_alternative)] <- 5
p3$Group <- paste0(p3$Mortality2, "_", p3$DAY)
group_means <- p3 %>%
   group_by(Group) %>% 
   summarise_at(vars(colnames(p)), mean)
    
group_means <- data.frame(group_means)
group_means <- t(group_means) 
colnames(group_means) <- group_means[1,]
group_means <- group_means[-1,] 
group_means <- data.frame(group_means)
group_means$transition <- rownames(group_means)

####
group_means$transition <- gsub("BCR_READS_Isotype_normalised_overlap_frequencies__", "", group_means$transition)
group_means$transition <- gsub("BCR_READS_Isotype_subsampled_50_overlap_frequencies__", "", group_means$transition)

group_means$transition <- gsub("IGHD_M","IGHD/M", group_means$transition)
groups <- str_split_fixed(group_means$transition, "_", 2)
group_means$class1 <- groups[,1]
group_means$class2 <- groups[,2]
#group_means$transition <- NULL

### Now this is in the same format as RACHAELS!!!!
class = strsplit(group_means$transition,"_",fixed = T)
class1 = NULL
class2 = NULL

for(i in c(1:length(class))){
	class1 = c(class1, class[[i]][1])
	class2 = c(class2, class[[i]][2])
}
match_edges = rep(-1, length(class1))

for(i in c(1:length(class))){
	w = intersect(which(edge_ids[,1]== class1[i]),which(edge_ids[,2]== class2[i]))
	if(length(w)==0){w = intersect(which(edge_ids[,1]== class2[i]),which(edge_ids[,2]== class1[i]))}
	if(length(w)==1){match_edges[i]=w}
	}

##########################################################################################################################
##########################################################################################################################
pdf(file=paste0(out_dir, "/Class_Switch_Longitudinal_SUBSAMPLED_", type_use, ".pdf"), height=12, width=12)
par(mfrow= c(4,3), mar = c(0.3,0.3,3.5,0.3))
patients <- colnames(group_means)[!colnames(group_means) == "HEALTH"]
patients <- patients[patients != "transition" & patients != "class1" & patients != "class2" & patients != "difference" & !patients %like% "PVAL"]
groppa <- group_means[, c(1:12)]
groppa <- apply(groppa, 2, function(x) as.numeric(as.character(x)))
groppa <- max(groppa)

for(d in 1:length(patients)){
	main = patients[d]
	mean_disease = group_means[,patients[d]]
	main2 <- gsub("_" ," Day ", main)
	## Health Comparison
	## Which health are we making comparison too?
	if(main %like% "_1"){
	healthx <- "HEALTH_1"
	pre_timepoint <- NA
	day <- 1
	mortgroup <- unlist(str_split(main, "_", 2))[1]
	main2 <- paste0(main2)
	} else if (main %like% "_3"){
	healthx <- "HEALTH_3"
	pre_timepoint <- paste0(unlist(str_split(main, "_", 2))[1], "_1")
	day <- 3
	oldday <- 1
	mortgroup <- unlist(str_split(main, "_", 2))[1]
	main2 <- paste0(main2, "\n versus ", mortgroup, " Day ", oldday )
	} else {
	healthx <- "HEALTH_5"
	pre_timepoint <- paste0(unlist(str_split(main, "_", 2))[1], "_3")
	day <- 5
	oldday <- 3
	mortgroup <- unlist(str_split(main, "_", 2))[1]
	main2 <- paste0(main2, "\n versus ", mortgroup, " Day ", oldday )
	}
	mean_health = group_means[,healthx]
	### Compare 
	group_means$difference <- as.numeric(group_means[,patients[d]]) - as.numeric(group_means[,healthx])
	####
	g1 = g
	edge_strength_plot = rep(0,length(edge_ids[,1]))
	edge_strength_plot[match_edges]= as.numeric(mean_disease)
	edge_strength_plot = (edge_strength_plot*10/groppa)#+0.1
	dir = rep(0,length(edge_ids[,1]))
	if(!is.na(pre_timepoint)){
	### want to test means 
	group_means$PVAL <- NA
	## Here we are calculating the stats comparing iso usage between sucessive timepoints 
	for(x in 1:length(rownames(group_means))){
		colx <- rownames(group_means)[x]
		data1 <- p3[, colx][p3$DAY %like% day & p3$Mortality2==mortgroup]
		data2 <- p3[, colx][p3$DAY %like% oldday & p3$Mortality2==mortgroup]
		p_value <- wilcox.test(data1, data2, exact=FALSE)$p.value
		group_means$PVAL[x] <- p_value
	}
	mean_previous <-  group_means[,pre_timepoint]
	pval = group_means$PVAL
	pval = p.adjust(pval, method = "BH") #bonferroni
	group_means$PVALadj =pval
	
	### Sig and Larger
	w = intersect(which(pval<0.05),which(mean_disease> mean_previous))
	dir[match_edges[w]] = 1
	### NS and Larger
	w = intersect(which(pval>=0.05),which(mean_disease> mean_previous))
	dir[match_edges[w]] = 5
	### Sig and SMALLER
	w = intersect(which(pval<0.05),which(mean_disease< mean_previous))
	dir[match_edges[w]] = 2
	### NS and SMALLER
	w = intersect(which(pval>=0.05),which(mean_disease< mean_previous))
	dir[match_edges[w]] = 6
	colnames(group_means)[colnames(group_means)=="PVAL"] <- paste0("PVAL_",mortgroup, "_", oldday, "_", day)
	colnames(group_means)[colnames(group_means)=="PVALadj"] <- paste0("PVALadj_",mortgroup, "_", oldday, "_", day)
	} else {
	dir = rep(3,length(edge_ids[,1]))
	}
	#############################################
	#### Plot the network graph 
	## We are going to plot 
	col = rgb(0.5,0.5,0.5,alpha = 0.5)
	cols = rep(col, length(edge_strength_plot))
	#names(cols) = names(pval)
	cols [which(dir==1)] = add.alpha("red",alpha = 0.9)
	cols [which(dir==2)] = add.alpha("blue",alpha = 0.9)
	cols [which(dir==5)] = add.alpha("darkgoldenrod1",alpha = 0.4)
	cols [which(dir==6)] = add.alpha("cyan",alpha = 0.4)
	cols [which(dir==3)] = add.alpha("black",alpha = 0.9)
	V(g1)$color = "grey"
	E(g1)$weight <- edge_strength_plot
	E(g1)$color =cols
	## Add significance
	dir2 <- dir
	dir2[!dir2==1 & !dir2==2] <- " "
	dir2[dir2==1 | dir2==2] <- "*"
	#E(g1)$label <- dir2
	g1_nonzero <- delete.edges(g1, which(E(g1)$weight == 0))
	new_weight <- E(g1_nonzero)$weight
	vertex_colors <- rainbow(vcount(g1_nonzero))
	main2 <- paste0(main2, "\nBH Padj")
	plot(g1_nonzero, layout=layout1,edge.color=E(g1_nonzero)$color, main=main2, edge.width=new_weight, vertex.label.color="black", vertex.color="light grey", vertex.label.family = "sans", vertex.label.font=2, vertex.label.cex = 1, edge.arrow.size = 1, edge.lty = 1, xlim = range(layout1[,1]*1.1), ylim = range(layout1[,2]*1.1))
	color_labels <- c("Sig Increase" = as.character(add.alpha("red", alpha = 0.9)[[1]]), "Sig Decrease"= as.character(add.alpha("blue",alpha = 0.9)[[1]]), "Increase"=as.character(add.alpha("darkgoldenrod1",alpha = 0.4)[[1]]), "Decrease"=as.character(add.alpha("cyan",alpha = 0.4)[[1]]), "No Comparison"=as.character(add.alpha("black",alpha = 0.9)[[1]]))
	color_labels <- data.frame(color_labels)
	color_labels$class <- rownames(color_labels)
	legend("bottomright", legend=color_labels$class, fill=color_labels$color_labels,  bty="n")
}
dev.off()

#color_labels$class[color_labels$color_labels %in% unique(E(g1_nonzero)$color)], fill=color_labels$color_labels[color_labels$color_labels %in% unique(E(g1_nonzero)$color)],
##########################################################################################################################
##########################################################################################################################
pdf(file=paste0(out_dir, "/Class_Switch_VsHealth_SUBSAMPLED_", type_use, ".pdf"), height=13, width=11)
par(mfrow= c(3,3), mar = c(0.3,0.3,3.5,0.3))
patients <- colnames(group_means)[!colnames(group_means) %like% "HEALTH"]
patients <- patients[patients != "transition" & patients != "class1" & patients != "class2" & patients != "difference" & !patients %like% "PVAL"]
groppa <- group_means[, c(1:12)]
groppa <- apply(groppa, 2, function(x) as.numeric(as.character(x)))
groppa <- max(groppa)

for(d in 1:length(patients)){
	main = patients[d]
	mean_disease = group_means[,patients[d]]
	main2 <- gsub("_" ," Day ", main)
	## Health Comparison
	## Which health are we making comparison too?
	if(main %like% "_1"){
	healthx <- "HEALTH_1"
	pre_timepoint <- NA
	day <- 1
	mortgroup <- unlist(str_split(main, "_", 2))[1]
	} else if (main %like% "_3"){
	healthx <- "HEALTH_3"
	pre_timepoint <- paste0(unlist(str_split(main, "_", 2))[1], "_1")
	day <- 3
	oldday <- 1
	mortgroup <- unlist(str_split(main, "_", 2))[1]
	} else {
	healthx <- "HEALTH_5"
	pre_timepoint <- paste0(unlist(str_split(main, "_", 2))[1], "_3")
	day <- 5
	oldday <- 3
	mortgroup <- unlist(str_split(main, "_", 2))[1]
	}
	main2 <- paste0(main2, "\n versus HEALTH Day ", day)
	mean_health = group_means[,healthx]
	#####################################
	####
	g1 = g
	edge_strength_plot = rep(0,length(edge_ids[,1]))
	edge_strength_plot[match_edges]= as.numeric(mean_disease)
	edge_strength_plot = (edge_strength_plot*10/groppa)#+0.1
	dir = rep(0,length(edge_ids[,1]))
	## want to test means 
	group_means$PVAL <- NA

	## Here we are calculating the stats comparing iso usage between health and disease!!!!!
	for(x in 1:length(rownames(group_means))){
		colx <- rownames(group_means)[x]
		data1 <- p3[, colx][p3$DAY %like% day & p3$Mortality2==mortgroup]
		data2 <- p3[, colx][p3$DAY %like% day & p3$Mortality2=="HEALTH"]
		p_value <- wilcox.test(data1, data2, exact=FALSE)$p.value
		group_means$PVAL[x] <- p_value
	}
	mean_previous <-  group_means[,healthx]
	pval = group_means$PVAL
	pval = p.adjust(pval, method = "BH")
	group_means$PVALadj =pval

	### Sig and Larger
	w = intersect(which(pval<0.05),which(mean_disease> mean_previous))
	dir[match_edges[w]] = 1
	### NS and Larger
	w = intersect(which(pval>=0.05),which(mean_disease> mean_previous))
	dir[match_edges[w]] = 5
	### Sig and SMALLER
	w = intersect(which(pval<0.05),which(mean_disease< mean_previous))
	dir[match_edges[w]] = 2
	### NS and SMALLER
	w = intersect(which(pval>=0.05),which(mean_disease< mean_previous))
	dir[match_edges[w]] = 6
	colnames(group_means)[colnames(group_means)=="PVAL"] <- paste0("PVAL_",mortgroup, "_", day, "_vsHEALTH")
	colnames(group_means)[colnames(group_means)=="PVALadj"] <- paste0("PVALadj_",mortgroup, "_", day, "_vsHEALTH")
	group_means$PVAL <- NULL
	group_means$PVALadj <- NULL
	#############################################
	#### Plot the network graph 
	## We are going to plot 
	col = rgb(0.5,0.5,0.5,alpha = 0.5)
	cols = rep(col, length(edge_strength_plot))
	#names(cols) = names(pval)
	cols [which(dir==1)] = add.alpha("red",alpha = 0.9)
	cols [which(dir==2)] = add.alpha("blue",alpha = 0.9)
	cols [which(dir==5)] = add.alpha("darkgoldenrod1",alpha = 0.4)
	cols [which(dir==6)] = add.alpha("cyan",alpha = 0.4)
	cols [which(dir==3)] = add.alpha("black",alpha = 0.9)
	V(g1)$color = "grey"
	E(g1)$weight <- edge_strength_plot
	E(g1)$color =cols
	## Add significance
	dir2 <- dir
	dir2[!dir2==1 & !dir2==2] <- " "
	dir2[dir2==1 | dir2==2] <- "*"
	#E(g1)$label <- dir2
	g1_nonzero <- delete.edges(g1, which(E(g1)$weight == 0))
	new_weight <- E(g1_nonzero)$weight
	vertex_colors <- rainbow(vcount(g1_nonzero))
	main2 <- paste0(main2, "\nBH Padj")
	plot(g1_nonzero, layout=layout1,edge.color=E(g1_nonzero)$color, main=main2, edge.width=new_weight, vertex.label.color="black", vertex.color="light grey", vertex.label.family = "sans", vertex.label.font=2, vertex.label.cex = 1.5, edge.arrow.size = 1, edge.lty = 1, xlim = range(layout1[,1]*1.1), ylim = range(layout1[,2]*1.1))
	color_labels <- c("Sig Increase" = as.character(add.alpha("red", alpha = 0.9)[[1]]), "Sig Decrease"= as.character(add.alpha("blue",alpha = 0.9)[[1]]), "Increase"=as.character(add.alpha("darkgoldenrod1",alpha = 0.4)[[1]]), "Decrease"=as.character(add.alpha("cyan",alpha = 0.4)[[1]]), "No Comparison"=as.character(add.alpha("black",alpha = 0.9)[[1]]))
	color_labels <- data.frame(color_labels)
	color_labels$class <- rownames(color_labels)
	legend("bottomright", legend=color_labels$class, fill=color_labels$color_labels, bty="n")
	}
dev.off()

### Save the output 
write.table(group_means, paste0(out_dir, "/ClassSwitchTransitionStats_", type_use, ".txt"), sep="\t")

