# Script to convert the constructed isotype files from trust4 into a format for isotyper!
## Lauren Overend lauren.overend@oriel.ox.ac.uk
## November 2021

library("optparse")

option_list <- list( 
  make_option(c("-o", "--outputdir"), action="store", type="character", default="NA", help="Path to Outputdir"),
  make_option(c("-s", "--sampleid"), action="store", type="character", help="Sampleid")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE)) 

outputdir <- opt$o
sampleid <- opt$s

my_aux_functions <- c("/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/RFunctions/TRUST4")           
source_files <- list.files(my_aux_functions, "*.R$", full.names=TRUE)  # locate all .R files
print(source_files)
for (f in source_files) {
    source(f)
}

# Run the convertor function
convert_trust4(outputdir, sampleid)

## Done