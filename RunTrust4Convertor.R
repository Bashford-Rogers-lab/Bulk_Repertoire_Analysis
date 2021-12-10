# Script to convert the constructed isotype files from trust4 into a format for isotyper!
## Lauren Overend lauren.overend@oriel.ox.ac.uk
## November 2021

library("optparse")

option_list <- list( 
  make_option(c("-o", "--outputdir"), action="store", type="character", default="NA", help="Path to Outputdir"),
  make_option(c("-s", "--sampleid"), action="store", type="character", help="Sampleid"), 
  make_option(c("-v", "--vthreshold"), action="store", type="character", help="Vtrimming"),
  make_option(c("-j", "--jthreshold"), action="store", type="character", help="Jtrimming")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE)) 

outputdir <- opt$o
sampleid <- opt$s
vtrim <- opt$v
jtrim <- opt$j

#outputdir='/gpfs2/well/immune-rep/shared/MISEQ/TRUST4_GAINS/V_20_J_15/'
#sampleid='gains8032905'
#vtrim='10'
#jtrim='10'

my_aux_functions <- c("/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/TRUST4")           
source_files <- list.files(my_aux_functions, "*.R$", full.names=TRUE)  # locate all .R files
print(source_files)
for (f in source_files) {
    source(f)
}

# Run the convertor function
convert_trust4(outputdir, sampleid, as.numeric(vtrim), as.numeric(jtrim))

## Done