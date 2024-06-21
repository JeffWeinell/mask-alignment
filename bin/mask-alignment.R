### load args passed to Rscript
args <- commandArgs(trailingOnly = TRUE)

ALNi_PATH_IN=args[1]
GAPLESSi_PATH_IN=args[2]
ALNi_PATH_OUT=args[3]
R_PACKAGES_DIR=args[4]

# optionally sets non-default R packages directory
if(R_PACKAGES_DIR!=""){
	.libPaths(R_PACKAGES_DIR)
}

# these must be installed to R_PACKAGES_DIR in advance
library(dplyr)
library(stringr)
library(Biostrings)     # BiocManager::install("Biostrings", lib="/your/Rpackages/directory/")
library(GenomicRanges)  # BiocManager::install("GenomicRanges", lib="/your/Rpackages/directory/")

### load ALNi_PATH_IN and GAPLESSi_PATH_IN
ALNi     <- Biostrings::readBStringSet(ALNi_PATH_IN)
GAPLESSi <- Biostrings::readBStringSet(GAPLESSi_PATH_IN)

### these should all be true
all(names(GAPLESSi) %in% names(ALNi))
all(names(ALNi) %in% names(GAPLESSi))

# only process seqs with names in ALNi and GAPLESSi
names.shared <- intersect(names(ALNi),names(GAPLESSi))
ALNi     <- ALNi[names.shared]
GAPLESSi <- GAPLESSi[names.shared]

# nongap positions for each sequence in ALNi
nongaps    <- str_locate_all(string=ALNi,pattern="[^-]")
nongaps.ir <- IRangesList(lapply(1:length(nongaps),function(x){z=nongaps[[x]]; IRanges(start=z[,"start"],end=z[,"end"],names=rep(names(ALNi[x]),nrow(z)))}))

# replace nongap positions in ALNi with characters in GAPLESSi
ALNi_OUT <- replaceAt(ALNi,nongaps.ir,strsplit(as.character(GAPLESSi),split=""))

# save ALNi_OUT
writeXStringSet(x=ALNi_OUT,filepath=ALNi_PATH_OUT)

