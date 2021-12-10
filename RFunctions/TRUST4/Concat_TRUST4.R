#module purge
#module use -a /apps/eb/dev/ivybridge/modules/all
#module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0

# Function to convert output of TRUST4 into isotype compatible format
# Also trims reads to location of best matching VDJ prior to uploading to IMGT (constant region encoded in gene name)
# Splits into BCR/TCR 
# Lauren Overend lauren.overend@oriel.ox.ac.uk
# November 2021
# Welcome Trust Centre for Human Genetics

#outputdir='/well/immune-rep/shared/MISEQ/TRUST4_GAINS/'
#sampleid<- 'gains8033985'

##############################################################################################
convert_trust4 <- function(outputdir=outputdir, sampleid=sample_id, v_threshold=v_threshold, j_threshold=j_threshold) {

	# load Recquired packages
	# all availible on rescomp Bioconductor 
	suppressMessages(library("Biostrings"))
	suppressMessages(library("stringr"))
	suppressMessages(library("tidyr"))
	suppressMessages(library("seqinr"))
	suppressMessages(library("gridExtra"))
	suppressMessages(library("grid"))
	suppressMessages(library("ggplot2"))
	suppressMessages(library("data.table"))
	
	#Create a directory to save summary plots
	if (!dir.exists(paste0(outputdir, "Plots"))){
		dir.create(paste0(outputdir, "Plots"))
	} else {
		print("Dir already exists!")
	}
	
	# read fastq file 
	path <- paste0(outputdir, sampleid, "_annot.fa")
	fastaFile <- readDNAStringSet(path)
	# extract names of samples 
	ids <- names(fastaFile)
	all_seq <- as.character(fastaFile)

	## extracting coordinates of CDR3 region
	ids2 <- data.frame(do.call('rbind', strsplit(ids,' ',fixed=TRUE)))
	a <- unlist(regmatches(ids2$X10, gregexpr("(?<=\\().*(?=\\))", ids2$X10, perl=TRUE)))
	ids2$LocationCDR3 <- a
	b <- str_split_fixed(ids2$LocationCDR3, "-", 2)
	ids2$CDR3Start <- b[,1]
	ids2$CDR3Ends <- b[,2]
	
	## Give each read a unique id based on the sample name
	ids2$sampleid <- paste0(sampleid, ids2$X1)
	ids2$sampleid <- gsub("assemble", "_", ids2$sampleid)

	#Finding out isotype/constant usage 
	c <- str_split_fixed(ids2$X7, "\\*", 2)
	c <- c[,1]
	ids2$isotype <- c

	## assign just the assembly name to the sequence
	## THIS WILL BE USED LATER!!!!!!!!! 
	names(all_seq) <- ids2$X1
	sequences <- data.frame("ID"=names(all_seq), "sequence"=all_seq, row.names=NULL)

	## Finding out counts of each assembly 
	counts <- read.delim(paste0(outputdir, sampleid, "_cdr3.out"), header=FALSE)
	countsx <- data.frame(table(counts$V1))

	## Merge to get count info 
	## this will also have the benefit of removing 'incomplete' transcripts 
	ids3 <- merge(ids2, countsx, by.x="X1", by.y="Var1")
	ids3$isotype[ids3$isotype==""] <- "IX"

	## Fill in isotype based on V gene 
	ids3$isotype[ids3$isotype=="IX" & ids3$X4 %like% "IGHV"] <- "IGHX"	
	ids3$isotype[ids3$isotype=="IX" & ids3$X4 %like% "TRBV"] <- "TRBX"
	ids3$isotype[ids3$isotype=="IX" & ids3$X4 %like% "TRGV"] <- "TRGX"
	ids3$isotype[ids3$isotype=="IX" & ids3$X4 %like% "TRAV" & ids3$X6 %like% "TRAJ" ] <- "TRAX"
	ids3$isotype[ids3$isotype=="IX" & ids3$X4 %like% "TRDV"] <- "TRDX"
	ids3$isotype[ids3$isotype=="IX" & ids3$X4 %like% "TRAV" & !ids3$X6 %like% "TRAJ" ] <- "TRADX"

	## Assigning Isotype to name as in RBR isotyper pipeline
	## BCR
	matrixiso <- matrix(,nrow=dim(ids3)[1], ncol=13)
	colnames(matrixiso) <-  c("sample", c("IGHA1", "IGHA2", "IGHD", "IGHE", "IGHEP2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM", "IGHX"))
	for(i in 1:length(ids3$sampleid)){
		matrixiso[, "sample"][i] <- ids3$sampleid[i]
		for(c in 2:length(colnames(matrixiso))){
		#print(colnames(matrixiso)[c])
			if (ids3$isotype[i]==colnames(matrixiso)[c]){
				#print(ids3$isotype[i])
				matrixiso[i,c] <- ids3$Freq[i]
			} else {
				#print("NO")
				matrixiso[i,c] <- 0
			}
		}
	}
	matrixiso <- data.frame(matrixiso)
	# Example ID >M01488331000000000JFPJG12115213403353__0_2_0_0_0_0_0_0_0_0_0|IGHA1_IGHA2_IGHD_IGHE_IGHEP2_IGHG1_IGHG2_IGHG3_IGHG4_IGHGP_IGHM	
	FULLIDisotype <- paste0(matrixiso$sample, "__", matrixiso$IGHA1, "_",  matrixiso$IGHA2, "_", matrixiso$IGHD, "_", matrixiso$IGHE, "_", matrixiso$IGHEP2, "_", matrixiso$IGHG1, "_", matrixiso$IGHG2, "_", matrixiso$IGHG3, "_", matrixiso$IGHG4, "_", matrixiso$IGHGP, "_", matrixiso$IGHM, "_", matrixiso$IGHX, "|IGHA1_IGHA2_IGHD_IGHE_IGHEP2_IGHG1_IGHG2_IGHG3_IGHG4_IGHGP_IGHM_IGHX")
	ids3$BCRID <- FULLIDisotype

	## TCR 
	matrixtcr <- matrix(,nrow=dim(ids3)[1], ncol=12)
	colnames(matrixtcr) <-  c("sample", c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "TRAX", "TRBX", "TRDX", "TRGX", "TRADX"))
	for(i in 1:length(ids3$sampleid)){
		matrixtcr[, "sample"][i] <- ids3$sampleid[i]
		for(c in 2:length(colnames(matrixtcr))){
		#print(colnames(matrixtcr)[c])
			if (ids3$isotype[i]==colnames(matrixtcr)[c]){
				#print(ids3$isotype[i])
				matrixtcr[i,c] <- ids3$Freq[i]
			} else {
				#print("NO")
				matrixtcr[i,c] <- 0
			}
		}
	}
	matrixtcr <- data.frame(matrixtcr)
	# Example ID >M00113396000000000JVVJY1111177789892__0_0_0_0_0_1|TRAC_TRBC1_TRBC2_TRDC_TRGC1_TRGC2
	FULLIDtcr <- paste0(matrixtcr$sample, "__", matrixtcr$TRAC, "_",  matrixtcr$TRBC1, "_", matrixtcr$TRBC2, "_", matrixtcr$TRDC, "_", matrixtcr$TRGC1, "_",  matrixtcr$TRGC2, "_", matrixtcr$TRAX, "_",  matrixtcr$TRBX, "_", matrixtcr$TRDX, "_", matrixtcr$TRGX, "_", matrixtcr$TRADX,  "|TRAC_TRBC1_TRBC2_TRDC_TRGC1_TRGC2_TRAX_TRBX_TRDX_TRGX_TRADX")
	ids3$TCRID <- FULLIDtcr

	## Merge so we have new name and sequence all in one dataframe 
	new <- merge(ids3, sequences, by.x="X1", by.y="ID")
	new$readlength <- (nchar(new$sequence)*-1)
	
	
	###################################################
	## Start to extract regions that we may want to trim on!!!
	
	## getting start points
	## There are some sequences with no V gene annotated - lets remove those 
	## Lets also remove those with no J gene annotated as we cant assign functionality / VDJ usage to these 
	new1 <- new[!new$X4=="*" & !new$X6=="*",]

	## Assign location of V and J genes 
	
	
	# Trim the sequence so that it ends at the CDR3 region
	#new$trimmed <- substr(new$sequence, 0, new$CDR3Ends)
	#new$trimmedlength <- (nchar(new$trimmed)*-1)
	
	## V GENE COORDINATES  #########################################################
	## String split by "'" to get the highest matching V gene and then take that column
	## Note this might mean exluding those that have different earlier start points 
	s <- data.frame(str_split_fixed(new1$X4, ",", 2)[,1])
	row.names(s) <- new1$X1
	colnames(s) <- c("VGENEloc")
	s$col2 <- unlist(regmatches(s$VGENEloc, gregexpr("(?<=\\):\\().*(?=\\):\\()", s$VGENEloc, perl=TRUE)))
	
	# Append coordinates to dataframe 
	new1$LocVstart <- s$col2
	b <- str_split_fixed(new1$LocVstart, "-", 2)
	new1$VStart <- b[,1]
	new1$VEnds <- b[,2]
	new1$CDR3Ends <- as.numeric(new1$CDR3Ends)
	new1$VStart <- as.numeric(new1$VStart)
	new1$VEnds <- as.numeric(new1$VEnds)
	
	# Look at length of v genes 
	new1$VtoV <- substr(new1$sequence, new1$VStart, new1$VEnds)
	new1$Vgenetrimmedlength <- (nchar(new1$VtoV)*-1)
	
	# Look at length of v to CDR3 genes 
	new1$VtoCDR3 <- substr(new1$sequence, new1$VEnds, new1$CDR3Ends)
	new1$VtoCDR3trimmedlength <- (nchar(new1$VtoCDR3)*-1)
	
	## J Gene COORDINATES #########################################################
	## String split by "'" to get the highest matching J gene and then take that column
	## Note this might mean exluding those that have different earlier start points 
	s <- data.frame(str_split_fixed(new1$X6, ",", 2)[,1])
	row.names(s) <- new1$X1
	colnames(s) <- c("JGENEloc")
	s$col2 <- unlist(regmatches(s$JGENEloc, gregexpr("(?<=\\):\\().*(?=\\):\\()", s$JGENEloc, perl=TRUE)))
	
	# Append coordinates to dataframe 
	new1$LocJstart <- s$col2
	b <- str_split_fixed(new1$LocJstart, "-", 2)
	new1$JStart <- b[,1]
	new1$JEnds <- b[,2]
	new1$JStart <- as.numeric(new1$JStart)
	new1$JEnds <- as.numeric(new1$JEnds)
	
	# Look at length of J genes 
	new1$JtoJ <- substr(new1$sequence, new1$JStart, new1$JEnds)
	new1$Jgenetrimmedlength <- (nchar(new1$JtoJ))
	
	# look at distance from CDR3 to end of J 
	new1$CDR3toJ <- substr(new1$sequence, new1$CDR3Ends, new1$JEnds)
	new1$CDR3toJtrimmedlength <- (nchar(new1$CDR3toJ))
	
	# J start to CDR3 end
	new1$JtoCDR3 <- substr(new1$sequence, new1$JStart, new1$CDR3Ends)
	new1$JtoCDR3trimmedlength <- (nchar(new1$JtoCDR3))
	
	
	## Can check if any of them dont have a J assigned beyond CDR3 contribution
	#any(new1$CDR3toJtrimmedlength<0)
	
	#########################################################
	## Length of sequence  #########################################################
	new1$FullLength <- nchar(new1$sequence)
	
	
	## Asign chain allowing us to facet_wrap
	new1$chain <- NA 
	new1$chain[new1$isotype %like% "IGH" ] <- "BCR"
	new1$chain[new1$isotype %like% "TRB" ] <- "TCRB"
	new1$chain[new1$isotype %like% "TRA" ] <- "TCRA"
	new1$chain[new1$isotype %like% "TRG" ] <- "TCRG"
	new1$chain[new1$isotype %like% "TRD" ] <- "TCRD"
	new1$chain[new1$isotype %like% "IGL" | new1$isotype %like% "IGK"] <- "BCR_LIGHT"
	new1$chain[new1$isotype %like% "IX" | new1$isotype %like% "TRADX"] <- "Undetermined"
	
	##############################################
	# CORECT THE CDR3 END AND THE V START LOCATION FOR TRIMMING!
	#####################################
	## Coorected V gene start location 
	## Caclulate new Vstart so its the same for all of them
	new1$VEnds <- as.numeric(new1$VEnds)
	new1$CorrectedVstart <- new1$VEnds-v_threshold
	#############################################
	## Coorected CDR3 END location gene start location 
	## Caclulate new Vstart so its the same for all of them
	new1$CDR3Ends <- as.numeric(new1$CDR3Ends)
	new1$CorrectedCDR3Ends <- new1$CDR3Ends+j_threshold
	new1$CorrectedEnds <- new1$JStart+j_threshold
	

	## What if we have a really long CDR3??	
	new1$CorrectedEnds[new1$CDR3Ends > new1$CorrectedEnds] <- new1$CDR3Ends[new1$CDR3Ends > new1$CorrectedEnds]
	
	
	## Filter out those which don't have 50bp of V gene length
	new2 <- new1[new1$CorrectedVstart >=0 & new1$VStart <= new1$CorrectedVstart,]
	
	## Not sure do we want to filter for CDR3 further than J gene end? Loses a lot of reads when might be relevant....
	## Filter those where the CDR3 corrected end is further than end of sequence length 
	new2 <- new2[new2$CorrectedEnds <= new2$JEnds | new2$Jgenetrimmedlength >= j_threshold,]
	new2 <- new2[new2$CorrectedEnds <= new2$FullLength,]
	
	# also remove any sequences where the corrected V start is now less than the annotated V start e.g. we are picking up upstream sequence
	# Maybe dont need to do incase there are other V genes?
	## Check with Rachael 
	new2 <- new2[new2$CorrectedVstart >= new2$VStart,]
	
	#################################################
	
	## Now we have everything we need to trim down to the V gene
	new2$Finaltrimmed <- substr(new2$sequence, new2$CorrectedVstart, new2$CorrectedEnds)
	new2$Finaltrimmedlength <- (nchar(new2$Finaltrimmed)*-1)
	
	new2$Vtrim <- substr(new2$sequence, new2$CorrectedVstart, new2$FullLength)
	new2$Vtrimlength <- (nchar(new2$Vtrim)*-1)
	
	new2$CDRtrim <- substr(new2$sequence, new2$CDR3Start, new2$CDR3Ends)
	new2$CDRtrimlength <- (nchar(new2$CDRtrim)*1)
	
	##############################################
	##Filter out any reads where the trimmed length is 0 
	new2 <- new2[!new2$Finaltrimmedlength==0,]
	################################################
	
	
	## Set up directory to output fasta files to 
	if (!dir.exists(paste0(outputdir, "TCRA_V", v_threshold, "_J", j_threshold))){
	dir.create(paste0(outputdir, "TCRA_V", v_threshold, "_J", j_threshold))
	} else {
		print("Dir already exists!")
	}
	if (!dir.exists(paste0(outputdir, "TCRB_V", v_threshold, "_J", j_threshold))){
	dir.create(paste0(outputdir, "TCRB_V", v_threshold, "_J", j_threshold))
	} else {
		print("Dir already exists!")
	}
	if (!dir.exists(paste0(outputdir, "TCRG_V", v_threshold, "_J", j_threshold))){
	dir.create(paste0(outputdir, "TCRG_V", v_threshold, "_J", j_threshold))
	} else {
		print("Dir already exists!")
	}
	if (!dir.exists(paste0(outputdir, "TCRD_V", v_threshold, "_J", j_threshold))){
	dir.create(paste0(outputdir, "TCRD_V", v_threshold, "_J", j_threshold))
	} else {
		print("Dir already exists!")
	}
	if (!dir.exists(paste0(outputdir, "BCR_V", v_threshold, "_J", j_threshold))){
	dir.create(paste0(outputdir, "BCR_V", v_threshold, "_J", j_threshold))
	} else {
		print("Dir already exists!")
	}

	## Final dataset 
	## NOTE WE WILL ALSO ADD THOSE WITH NO CONSTANT REGION ASSIGNED TO THIS DATASET
	bcr_only <- new2[new2$isotype %like% "IGH",] 
	final_data_bcr <- bcr_only[, c("BCRID", "TCRID", "Finaltrimmed")]
	final_data_bcr <- as.list(c(bcr_only$Finaltrimmed))
	names(final_data_bcr) <- bcr_only$BCRID
	outfile_name <- paste0(outputdir, "BCR_V", v_threshold, "_J", j_threshold, "/", sampleid, "_BCR.fa")
	write.fasta(final_data_bcr, names(final_data_bcr), outfile_name, open = "w", nbchar = 100, as.string = FALSE)

	## Final dataset 
	## NOTE WE WILL ALSO ADD THOSE WITH NO CONSTANT REGION ASSIGNED TO THIS DATASET
	##TCRA
	tcra_only <- new2[new2$isotype %like% "TRA",] 
	final_data_tcra <- tcra_only[, c("BCRID", "TCRID", "Finaltrimmed")]
	final_data_tcra <- as.list(c(tcra_only$Finaltrimmed))
	names(final_data_tcra) <- tcra_only$TCRID
	outfile_name <- paste0(outputdir, "TCRA_V", v_threshold, "_J", j_threshold, "/", sampleid, "_TCRA.fa")
	write.fasta(final_data_tcra, names(final_data_tcra), outfile_name, open = "w", nbchar = 100, as.string = FALSE)
	
	##TCRB
	tcrb_only <- new2[new2$isotype %like% "TRB" ,] 
	final_data_tcrb <- tcrb_only[, c("BCRID", "TCRID", "Finaltrimmed")]
	final_data_tcrb <- as.list(c(tcrb_only$Finaltrimmed))
	names(final_data_tcrb) <- tcrb_only$TCRID
	outfile_name <- paste0(outputdir, "TCRB_V", v_threshold, "_J", j_threshold, "/", sampleid, "_TCRB.fa")
	write.fasta(final_data_tcrb, names(final_data_tcrb), outfile_name, open = "w", nbchar = 100, as.string = FALSE)
	
	##TCRG
	tcrg_only <- new2[new2$isotype %like% "TRG", ] 
	final_data_tcrg <- tcrg_only[, c("BCRID", "TCRID", "Finaltrimmed")]
	final_data_tcrg <- as.list(c(tcrg_only$Finaltrimmed))
	names(final_data_tcrg) <- tcrg_only$TCRID
	outfile_name <- paste0(outputdir, "TCRG_V", v_threshold, "_J", j_threshold, "/", sampleid, "_TCRG.fa")
	write.fasta(final_data_tcrg, names(final_data_tcrg), outfile_name, open = "w", nbchar = 100, as.string = FALSE)
	
	##TCRG
	tcrd_only <- new2[new2$isotype %like% "TRD", ] 
	final_data_tcrd <- tcrd_only[, c("BCRID", "TCRID", "Finaltrimmed")]
	final_data_tcrd <- as.list(c(tcrd_only$Finaltrimmed))
	names(final_data_tcrd) <- tcrd_only$TCRID
	outfile_name <- paste0(outputdir, "TCRD_V", v_threshold, "_J", j_threshold, "/", sampleid, "_TCRD.fa")
	write.fasta(final_data_tcrd, names(final_data_tcrd), outfile_name, open = "w", nbchar = 100, as.string = FALSE)
	
	## extract tcr for plotting
	tcr_only <- new2[new2$isotype %like% "TR",] 

	
	
	# Final Plot comparing TCR and BCR 
	pdf(paste0(outputdir, "Plots/TRUST4_Trimming_", sampleid,"_V", v_threshold, "_J", j_threshold, ".pdf"))
	# V gene
	c <- ggplot(data=new1, aes(VtoCDR3trimmedlength)) + geom_histogram(binwidth=10)+ xlab("V-End to CDR3-End Length") +theme_bw() +ylab("Count") +ggtitle(paste0(sampleid)) +facet_wrap(~chain) +xlim(NA, 0)
	grid.arrange(c,  ncol=1)
	c <- ggplot(data=new1, aes(Vgenetrimmedlength)) + geom_histogram(binwidth=10)+ xlab("V-Start to V-End") +theme_bw() +ylab("Count") +ggtitle(paste0(sampleid))+facet_wrap(~chain) +geom_vline(xintercept=(-1*v_threshold), col="red")+xlim(NA, 0)
	grid.arrange(c,  ncol=1)
	c <- ggplot(new1,aes(Vgenetrimmedlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("V-Start to V-End") +ylab("Cumulative Frequency") +theme_bw()+ggtitle(paste0(sampleid))+facet_wrap(~chain, scales="free_y") +geom_vline(xintercept=(-1*v_threshold), col="red")+xlim(NA, 0)
	grid.arrange(c,  ncol=1)
	
	# J gene
	c <- ggplot(data=new1, aes(Jgenetrimmedlength)) + geom_histogram(binwidth=2)+ xlab("J-Start to J-End") +theme_bw() +ylab("Count") +ggtitle(paste0(sampleid)) +facet_wrap(~chain) +xlim(0, NA)+geom_vline(xintercept=j_threshold, col="green")
	grid.arrange(c,  ncol=1)
	c <- ggplot(new1,aes(Jgenetrimmedlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("J-Start to J-End") +ylab("Cumulative Frequency") +theme_bw()+ggtitle(paste0(sampleid))+facet_wrap(~chain, scales="free_y") +xlim(0, NA)+geom_vline(xintercept=j_threshold, col="green")
	grid.arrange(c,  ncol=1)
	c <- ggplot(data=new1, aes(CDR3toJtrimmedlength)) + geom_histogram(binwidth=2)+ xlab("CDR3-End to J-End") +theme_bw() +ylab("Count") +ggtitle(paste0(sampleid))+facet_wrap(~chain) 
	grid.arrange(c,  ncol=1)
	c <- ggplot(new1,aes(CDR3toJtrimmedlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("CDR3-End to J-End") +ylab("Cumulative Frequency") +theme_bw()+ggtitle(paste0(sampleid))+facet_wrap(~chain, scales="free_y")
	grid.arrange(c,  ncol=1)
	c <- ggplot(data=new1, aes(JtoCDR3trimmedlength)) + geom_histogram(binwidth=2)+ xlab("J-Start to CDR3-End") +theme_bw() +ylab("Count") +ggtitle(paste0(sampleid))+facet_wrap(~chain)+xlim(0, NA)+geom_vline(xintercept=j_threshold, col="blue")
	grid.arrange(c,  ncol=1)
	c <- ggplot(new1,aes(JtoCDR3trimmedlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("J-Start to CDR3-End") +ylab("Cumulative Frequency") +theme_bw()+ggtitle(paste0(sampleid))+facet_wrap(~chain, scales="free_y") +xlim(0, NA)+geom_vline(xintercept=j_threshold, col="blue")
	grid.arrange(c,  ncol=1)

	#Starting length
	c <- ggplot(data=new2, aes(readlength)) + geom_histogram(binwidth=10)+ xlab("STARTING READ LENGTH") +theme_bw() +ylab("Count") +ggtitle(paste0("Pre-Trim: ", sampleid))+xlim(NA, 0)
	c1 <- ggplot(new2,aes(readlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50) + xlab("STARTING READ LENGTH") +ylab("Cumulative Frequency") +theme_bw()+ggtitle(paste0("Pre-Trim: ", sampleid))+xlim(NA, 0)
	# Post V trim (V to End of sequence)
	a <- ggplot(data=new2, aes(Vtrimlength)) + geom_histogram(binwidth=10)+ xlab("TRIMMED V-End LENGTH") +theme_bw() +ylab("Count") +ggtitle(paste0("Trim 1 (V start): ", sampleid))+xlim(NA, 0)
	a1 <- ggplot(new2,aes(Vtrimlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50) + xlab("TRIMMED V-End LENGTH") +ylab("Cumulative Frequency") +theme_bw()+ggtitle(paste0("Trim 1 (V start): ", sampleid))+xlim(NA, 0)
	# Post J trim (V to End of sequence)
	b <- ggplot(data=new2, aes(Finaltrimmedlength)) + geom_histogram(binwidth=10) + xlab("TRIMMED V-J LENGTH")+theme_bw()+ylab("Count") +ggtitle(paste0("Trim 2 (J End): ", sampleid))+xlim(NA, 0)
	b1 <- ggplot(new2,aes(Finaltrimmedlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-J LENGTH") +ylab("Cumulative Frequency") +theme_bw()+ggtitle(paste0("Trim 2 (J End): ", sampleid))+xlim(NA, 0)
	# Post CDR trim length:
	d <- ggplot(data=new2, aes(CDRtrimlength)) + geom_histogram(binwidth=10) + xlab("CDR3 LENGTH")+theme_bw()+ylab("Count") +ggtitle(paste0("CDR3: ", sampleid))
	d1 <- ggplot(new2,aes(CDRtrimlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("CDR3 LENGTH") +ylab("Cumulative Frequency") +theme_bw()+ggtitle(paste0("CDR3: ", sampleid))
	grid.arrange(c, c1, a,a1, b, b1, d, d1, ncol=2)
	
	## Look by Chain
	a <- ggplot(data=bcr_only, aes(Finaltrimmedlength)) + geom_histogram(binwidth=10) + xlab("TRIMMED V-J LENGTH")+theme_bw() +ggtitle(paste0("BCR: ", sampleid)) +ylab("Count")+xlim(NA, 0)
	a1 <- ggplot(bcr_only,aes(Finaltrimmedlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-J LENGTH") +ylab("Cumulative Frequency") +theme_bw() +ggtitle(paste0("BCR: ", sampleid))+xlim(NA, 0)
	b <- ggplot(data=tcr_only, aes(Finaltrimmedlength)) + geom_histogram(binwidth=10) + xlab("TRIMMED V-J LENGTH")+theme_bw() +ggtitle(paste0("TCR: ", sampleid))+ylab("Count")+xlim(NA, 0)
	b1 <- ggplot(tcr_only,aes(Finaltrimmedlength))+stat_bin(aes(y=cumsum(..count..)),geom="step",bins=50)+ xlab("TRIMMED V-J LENGTH") +ylab("Cumulative Frequency") +theme_bw() +ggtitle(paste0("TCR: ", sampleid))+xlim(NA, 0)
	grid.arrange(a,a1, b, b1, ncol=2)
	dev.off()

	# done 
}
