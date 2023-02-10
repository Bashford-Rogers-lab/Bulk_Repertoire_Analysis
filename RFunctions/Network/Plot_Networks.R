library(igraph)
library(RColorBrewer)
#setwd('/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE')
file_ids <- 'LEO_SEPSIS_BCR_ALL_post.txt'
plot_networks <- function(outputdir, file_ids)
		
	ids <- read.delim(file_ids, sep="\t", header=FALSE)
	ids_group <- ids[,1]
	dir <- paste0(outputdir, "ORIENTATED_SEQUENCES/")

	### TO RUN: CHANGE DIRECTORY NAME TO RELEVENT LOCATION AND BATCH_NAME TO RELEVENT NAME (i.e. NAME OF GROUP OF SAMPLES), THEN RUN CODE

	# dir = 'XXXXX' 
	# batch_name = 'YYYYY'
	#############################################################################################################################
	add.alpha <- function(col, alpha=1){
		  if(missing(col))
			stop("Please provide a vector of colours.")
			apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha)) 
	}
	colour <- rep(brewer.pal(7, "Set1"),2)[c(1:5)]
	colour1<-add.alpha("red", alpha = 0.8)
	colour<-add.alpha(colour1, alpha = 0.8)

#############################################################################################################################

group = 4
s = 1
type_plot = "CCFAST"
batch_name = "TRIAL1_IsoTypera"
index= 1
indexi = 0
fileout=paste0(dir,"Network_nucleotide_", type_plot,"_",batch_name,"_PIE_",index,".jpeg")
q = 7016

pdf(file=fileout, height=q*1.1*2, width=q*4)
par(mfrow= c(2,4), mar = c(0,0,3,0))

for(s in c(1:length(ids_group))){
	indexi = indexi+1
	if(indexi>4){
		dev.off()
		index = index+1
		fileout=past0(dir,"Network_nucleotide_", type_plot,"_",batch_name,"_PIE_",index,".jpeg")
		q = 7016
		jpeg(file=fileout, height=q*1.1, width=q*4, res=600)
		par(mfrow= c(1,4), mar = c(0,0,3,0))
		indexi = 1
		}
	sample = ids_group[s]
	
	########## Generate the IGRAPHS
	n=paste0(dir,"NETWORKS/","Cluster_identities_", sample,".txt")
	e=paste0(dir,"NETWORKS/","Edges_", sample,".txt")
	nodes <- read.csv(n, head=F, sep="\t")
	nodes = nodes[-1,]
	edge <- read.csv(e, head=FALSE, sep="\t")
	
	cluster<-as.character(nodes[,2])
	freq1<-as.numeric(nodes[,4])
	total_freq = sum(freq1)
	names1 = as.character(nodes[,3])
	freq_sub = freq1
	freq = freq_sub*50/sum(freq_sub)
	
	
	g <- graph.empty( n=0, directed=FALSE)
	g <- igraph::add.vertices(g, length(names1), name= names1,color = colour[freq1/2+2]) 
	names <- V(g)$name
	ids <- 1:length(names)
	names(ids) <- names
	from <- as.character(edge[,1])
	to <- as.character(edge[,2])
	w = intersect(which(from %in% names),which(to %in% names))
	edges <- matrix(c(ids[from[w]], ids[to[w]]), nc=2)
	g <- add.edges(g, t(edges), weight=1)
	V(g)$label <- V(g)$name
	V(g)$size<-freq1
	V(g)$label.cex<-0.0001
	g1 = g
	g2 = g
	#### nodes to include using previously generated subsampling
	file=paste0(dir,"ANNOTATIONS/","Subsample_IDs_", type_plot,"_", sample,".txt")
	p <- read.csv(file, head=F, sep="\t")
	ids_subsampled = as.character(p[-1,1])	
	gx = induced_subgraph(g, ids_subsampled)
	csize = components(gx)$csize
	csize_real = components(g)$csize
	cluster_length= length(csize)
	max_cluster = c(max(csize_real)*100/sum(csize_real), max(csize)*100/sum(csize))
	g1 = gx
	
	title2 = sample
	print (title2)
	#layout = layout_with_graphopt(g1,niter = 1200) ############## REAL
	layout = layout_with_graphopt(g1,niter = 900, charge = 0.01) # REMOVE
	names_sub  =V(g1)$name
	genes = strsplit(strsplit(names_sub[1],"|", fixed = TRUE)[[1]][2], "_",fixed = TRUE)[[1]]
	m = matrix(data = 1, nrow = length(names_sub), ncol = length(genes), dimnames = c(list(names_sub), list(genes)))
	names_short= NULL
	for(i in c(1:length(names_sub))){
	      names_short = c(names_short ,strsplit(names_sub[i],"__", fixed = TRUE)[[1]][1])
	  freqs = as.numeric(strsplit(strsplit(strsplit(names_sub[i],"|", fixed = TRUE)[[1]][1], "__",fixed = TRUE)[[1]][2], "_",fixed = TRUE)[[1]])
	      m[i,]= freqs
	}

	genes_of_interest = c("IgM","IgD","IgA1","IgA2","IgG1","IgG2","IgG3","IgG4","IgE")
	genes_of_interest1 = gsub("Ig","IGH", genes_of_interest)
	
	m1 = matrix(data = 1, nrow = length(names_sub), ncol = length(genes_of_interest), dimnames = c(list(names_sub), list(genes_of_interest)))
	for(i in c(1:length(genes_of_interest))){
		w = grep(genes_of_interest1[i], genes, fixed = TRUE)
		if(genes_of_interest1[i]=="IGHA1/2"){w = c(w, grep("IGHA1", genes, fixed = TRUE), grep("IGHA2", genes, fixed = TRUE))}
	    if(genes_of_interest1[i]=="IGHD/M"){w = c(w, grep("IGHD", genes, fixed = TRUE), grep("IGHM", genes, fixed = TRUE))}
	    if(genes_of_interest1[i]=="IGHG1/2"){w = c(w, grep("IGHG1", genes, fixed = TRUE), grep("IGHG2", genes, fixed = TRUE))}
	    f = rep(0,length(m[,1]))
	    if(length(w)>1){
			f = rowSums(m[,w])
		}else{if(length(w)==1){f = m[,w]}}
		m1[,i] = f
	}
  	library(RColorBrewer)
	cols1 = rep(brewer.pal(12, "Paired"),2)
	cols <- add.alpha (cols1, alpha = 0.9)
	values = list()
	for(i in c(1:length(m1[,1]))){
	    a = m1[i,]
	  if(max(a)==0){a = a+1}
	  values = c(values, list(a))
	}	
	sizes = rowSums(m1)
	sizes = sizes/sum(sizes)
	power = 0.3
	sizes = sizes^power
	sizes = sizes*3/(0.5^power)
	sizes = sizes+1
	range(sizes)

	plot(g1, vertex.shape="pie", layout=layout, edge.color=add.alpha("black",alpha = 0.8),main=title2, edge.width=2, vertex.pie=values,vertex.pie.color=list(cols),vertex.size= sizes, vertex.label=NA)
	print (s)	
}
dev.off()



#### Test colour palette:

cex = 0.9
fileout=concat(c(dir,"Plot_key.jpeg"))
q = 7016
jpeg(file=fileout, height=q*1.1, width=q*1, res=600)
par(mfrow= c(1,1), mar = c(0,0,3,0))

plot(c(0,1),c(0,1), pch=20, col="white",xlab="",ylab ="",cex=cex, cex.lab=cex+0.1,	cex.axis=cex,cex.main = cex, col.axis = "white",tck=0, mgp = c(2,0,0), main = "", axes = FALSE)

legend("topright",  legend= genes_of_interest, pch= 21, bty="n",cex= cex-0.1, title="", pt.bg = cols,pt.cex = 1)
dev.off()

print(concat(c("scp -p mfj169@rescomp1.well.ox.ac.uk:", fileout," ./ " )))
