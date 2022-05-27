outputdir <- /well/immune-rep/shared/MISEQ/SEPSIS_FINAL/TCRA_CH2/
type <- "TCRA"


	
data_all <- read.delim('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_FINAL/TCRA_CH2/ORIENTATED_SEQUENCES/ISOTYPER/All_V_gene_isotype_frequency_ALL.txt')
data_all <- data_all[data_all$class=="TRAC",]
data_all$frequency <- as.numeric(data_all$frequency)
data_all$total <- NA


for(i in 1:length(unique(data_all$X.sample))){
	sampler <- unique(data_all$X.sample)[i]
	data_all$total[data_all$X.sample==sampler] <- sum(data_all$frequency[data_all$X.sample==sampler])
	}
	
data_all$prop <- data_all$frequency/data_all$total

## Invariant

data_all$Type <- "Standard T cell"
data_all$Type[data_all$V.gene=="TRAV1-2" & data_all$J.gene=="TRAJ33"] <- "PossibleMAIT"
data_all$Type[data_all$V.gene=="TRAV1-2" & data_all$J.gene=="TRAJ12"] <- "PossibleMAIT"
data_all$Type[data_all$V.gene=="TRAV1-2" &  data_all$J.gene=="TRAJ20"] <- "PossibleMAIT"

data_all$Type[data_all$V.gene=="TRAV1-2" & data_all$J.gene=="TRAJ9"] <- "PossibleGEM"
data_all$Type[data_all$V.gene=="TRAV24" & data_all$J.gene=="TRAJ18"] <- "PossibleNKT"

invariant <- data_all[data_all$Type=="PossibleMAIT" | data_all$Type=="PossibleGEM" | data_all$Type=="PossibleNKT", ]
invariant$day <- str_split_fixed(invariant$X.sample, "_", 2)[,2]
invariant$day <- as.numeric(invariant$day)
invariant$day[is.na(invariant$day)] <- "CONTROL"

pdf("MAIT.pdf", width=10, height=10)

ggplot(invariant, aes(x=day, y=prop, fill=Type))+geom_boxplot()+theme_bw()
dev.off()


table(invariant$X.sample, invariant$Type, invariant$prop)