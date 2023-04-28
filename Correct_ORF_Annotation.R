### code to correct the annotation of ORF genes from IGMT output #
## Particulary important for TCRG!!!!
## Lauren Overend 
## Lauren_overend@live.co.uk
module purge
module use -a /apps/eb/dev/ivybridge/modules/all
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0

outputdir <- '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG_NEW/ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW_OLD/LEO_SEPSIS_ALL_TCRG_1'
files <- list.files(outputdir, full.name=TRUE)
files <- grep(".txt", files, value=TRUE)
files <- files[c(3, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
new_outdir <- '/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG_NEW/ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW/LEO_SEPSIS_ALL_TCRG_1'
### 

## This will loop through IMGT files, if they have an ORF annotation it will assign them as non productive 
## This is important for TCRG where the big TRGV10 is non functional 
for(i in 1:length(files)){
	######### Summary File 1 
	summaryx<- read.delim(files[i])
	prop_names <- colnames(read.delim(files[i], check.names=FALSE))
	summary2 <- summaryx 
	#summary2$X <- NULL
	print(files[i])
	if("V.DOMAIN.Functionality" %in% colnames(summaryx) & files[i] %like% "1_Summary"){
		print("SUMMARY FILE")
		#### All those annotated with ORF 
		unique(summaryx$V.GENE.and.allele[summaryx$V.GENE.and.allele %like% "ORF"]) 
		m <- summaryx[summaryx$V.GENE.and.allele %like% "ORF",]
		print(dim(m))
		table(m$V.DOMAIN.Functionality)
		## If they are like ORF we want them to be "unproductive"
		summary2$V.DOMAIN.Functionality[summaryx$V.GENE.and.allele %like% "ORF" & (summary2$V.DOMAIN.Functionality=="productive"| summary2$V.DOMAIN.Functionality=="productive (see comment)")] <- "unproductive"
		###
		print(table(summary2$V.DOMAIN.Functionality)[names(table(summary2$V.DOMAIN.Functionality))=="productive"])
		## We will use this corrected assignment in all the other files to correctly annotate functionality 
		assignment <- summary2[, c("Sequence.ID", "V.DOMAIN.Functionality")]
		### Assign the normal column names 
		#colnames(summary2) <- prop_names
		## We lose about 60% of reads but these likely came from the IGHV10 
		new_name <- gsub("IMGT_RAW_OLD", "IMGT_RAW", files[i])
		write.table(summary2, new_name, sep="\t", row.names=FALSE, quote=FALSE,na = "", col.names= prop_names)
		
	} else if("V.DOMAIN.Functionality" %in% colnames(summaryx) & !files[i] %like% "1_Summary"){
		print("file contains V domain functionality annotation but not summary file")
		#### All those annotated with ORF 
		unique(summaryx$V.GENE.and.allele[summaryx$V.GENE.and.allele %like% "ORF"]) 
		m <- summaryx[summaryx$V.GENE.and.allele %like% "ORF",]
		print(dim(m))
		table(m$V.DOMAIN.Functionality)
		index1 <- summary2$Sequence.ID
		index2 <- assignment$Sequence.ID
		## check order
		if(all(index1==index2)){
			print('Perfect match in same order')
			summary2$V.DOMAIN.Functionality <- assignment$V.DOMAIN.Functionality
		} else {
			print("Warning not same order")
		}
		## We lose about 60% of reads but these likely came from the IGHV10 
		print(table(summary2$V.DOMAIN.Functionality)[names(table(summary2$V.DOMAIN.Functionality))=="productive"])
		### Assign the normal column names 
		new_name <- gsub("IMGT_RAW_OLD", "IMGT_RAW", files[i])
		#colnames(summary2) <- prop_names
		write.table(summary2, new_name, sep="\t", row.names=FALSE, quote=FALSE,na = "" , col.names= prop_names)
	} else {
		print("file does not contain V domain functionality annotation")
		new_name <- gsub("IMGT_RAW_OLD", "IMGT_RAW", files[i])
		#colnames(summary2) <- prop_names
		write.table(summary2, new_name, sep="\t", row.names=FALSE, quote=FALSE,na = "", col.names= prop_names)
	}
}

old <- read.delim("/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG_NEW/ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW_OLD/LEO_SEPSIS_ALL_TCRG_1/1_Summary.txt")
new <- read.delim("/gpfs3/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/TCRG_NEW/ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW/LEO_SEPSIS_ALL_TCRG_1/1_Summary.txt")

### Note you need to create a BLANK TXZ file in the same directory as the IMGT_RAW becuase Rachael uses this as a count to iterate through and run the code!!!
### Also set EXTRACT =FALSE in the pipeline as the files are already extracted you dont want to rerun!!!