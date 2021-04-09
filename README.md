# Bulk Repertoire Analysis

This Repository contains the code required to preform the BCR/TCR pre-processing pipeline for NGS data based on the IsoTyper and TCR protocols developed in the Bashford-Rogers Lab. For handling multiple samples we provide two solutions: a) a python based wrapper for job submission using bsub (Processing_sequences_large_scale.py) and b) a standard job submission bash script using qsub (BCR_TCR_Wrapper_Cluster.sh). The latter is preferable for running on the BMRC cluster (Recomp) as it utilises the module system and could be easily adapted for another cluster architecture. 

*An indepth guide to installation/using the pipeline can be found in the BCR_TCR_preprocessingmanual2.0.pdf* 

## Summary 
### Stage 1: 
1.	QC sequences
### Stage 2: 
1.	Join forward and reverse reads (merging)
2.	Split sequences according to sample barcode
3.	Identify RNA barcode and collapse/error correct based on groups of sequences sharing same barcode
4.	Check isotype against reference
5.	Check matches to IGHV/J reference sequences
6.	Check open reading frame present
### Stage 3: 
1.	Network generation: Here, each vertex represents a different sequence, and the number of identical BCR sequences defines the vertex size. Edges are created between vertices that differ by one nucleotide. Clusters are groups of interconnected vertices (1, 2). The program described here calculates edges between unique sequences and determines vertex sizes, creating output files in formats that can directly be used in network analysis programs such as networkx (python) or igraph (R or python).
### Stage 4: 
1.	Sequence annotation
2.	Network Analysis 
3.	Generation of broad repertoire statistics
### Stage 5 (optional): 
1.	Run the optional R script to analyse and visualise summary metrics from the results of Stages 1-4. This will also concatenate files for all samples.  
2.	Check the percentage of reads which pass Open-Read-Frame filtering. 
a.	If this is abnormally low (potentially due to large clonal expansion) consider running the analysis with the ORF column in the sample sheet set to anything other than TRUE. To prevent reads being removed. 
### Stage 6: 
1.	Concatenate filtered fastq files into a smaller number of multi-individual large files. 
2.	Upload large files to IMGT for annotation. 

### Stage 7: 
1.	Results of IMGT analysis are used in Isotyper specific Analysis 

## Register for IMGT
As part of the analysis pipeline fatsq files will be annotated using IMGT/HighV-Quest – the gold standard for TCR and BCR repertoire analysis. However to use the high-throughput (HighV) version you must register a user account. All new users must be approved and accounts are then activated by an administrator. Therefore we recommend registering as soon as possible to avoid hold-ups later on. 

Register [here](http://www.imgt.org/HighV-QUEST/login.action) 

# Installation

Note: Jobs should be run in the directory containing the pipeline so that relative paths are used. 
•	Ensure access to Local_immune_repertoire_annotator_1.0.py from:
o	ANNOTATION_OF_TCRs_CDR3_REGIONS/Local_immune_repertoire_annotator_1.0.py
## 1.	FLASH: 
•	Download current version from [here](https://ccb.jhu.edu/software/FLASH/). 
•	Unpack. 
## 2.	CD-HIT(3):
* Download current CD-HIT from [here](http://bioinformatics.org/cd-hit/).  
* Unpack the file with: 
`tar xvf cd-hit-XXX.tar.gz --gunzip`
* Change dir by: 
`cd cd-hit-2006`
* Compile the program by: 
`make`
## 3.	Quasr(4):
* Download current version from [here](https://sourceforge.net/projects/quasr/). 
## 4.	Blast:
* Download current version from [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).  
* **For users of Rescomp (BMRC, Oxford) please not this is already installed so you do not need to perform this stage.** 
## 5.	Python Modules:
Ensure that the following python modules are installed.
**Note for Rescomp (BMRC, Oxford) users: These modules are already available using the Rescomp ‘module load’ system and need not be installed locally. However to enable access to these modules (load) you MUST run analysis using the BCR_TCR_Wrapper_Cluster.sh job submission script.** 
•	Sys
•	Collections
•	Os
•	Operator
•	networkx 
## 6.	R Modules (optional analysis): 
Ensure that the following R modules are installed: 
•	tidyverse
•	ggplot2
•	foreach
•	doParallel
•	gridExtra
•	cowplot  
## 7.	Locations of Dependencies: 
•	To ensure the pipeline can call the dependencies edit: “BCR_TCR_PROCESSING_PIPELINE/Locations_of_called_programmes.txt” file, providing the full path to correct locations of your versions of: 
o	Reference library of genes and primers (already compiled for you)
o	CD-HIT (from above)
o	FLASH (from above)
o	Quasr (4) (from above)
o	Blast (from above)
## 8.	Create Log Files Directory
•	When using the BCR_TCR_Wrapper_Cluster.sh wrapper all log files will be output to a directory called COMMANDLOGS within the current working directory (directory containing pipeline). However this must be created prior to running the job submission wrapper using the following bash script e.g.:  
`cd path_to/BCR_TCR_PROCESSING_PIPELINE`
`mkdir COMMANDLOGS`


# Author Contribution 

Rachael J. M. Bashford-Rogers developed the python based TCR/BCR repertoire analysis pipeline and User Guide. 
Lauren E. Overend developed the bash wrapper, R functions and helped write documentation. 

# References 

If you find Immune-Network-Generation useful, please cite reference #2 (R. J. Bashford-Rogers et al., 2019). 

1.	R. J. Bashford-Rogers et al., Network properties derived from deep sequencing of human B-cell receptor repertoires delineate B-cell populations. Genome Res 23, 1874-1884 (2013).
2.	R. J. M. Bashford-Rogers et al., Analysis of the B cell receptor repertoire in six immune-mediated diseases. Nature,  (2019).
3.	W. Li, A. Godzik, Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics 22, 1658-1659 (2006).
4.	S. J. Watson et al., Viral population analysis and minority-variant detection using short read next-generation sequencing. Philos Trans R Soc Lond B Biol Sci 368, 20120205 (2013).
5.	T. Magoc, S. L. Salzberg, FLASH: fast length adjustment of short reads to improve genome assemblies. Bioinformatics 27, 2957-2963 (2011).


