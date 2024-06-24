# these must be installed to R_PACKAGES_DIR in advance
library(dplyr)
library(stringr)
library(Biostrings)
library(GenomicRanges)

replace_nongaps_with <- function(ALN_FA_UNMASKED_PATH,GAPLESS_FA_PATH,ALN_FA_MASKED_PATH) {
	### load ALN_FA_UNMASKED_PATH and GAPLESS_FA_PATH
	ALNi     <- Biostrings::readBStringSet(ALN_FA_UNMASKED_PATH)
	GAPLESSi <- Biostrings::readBStringSet(GAPLESS_FA_PATH)
	
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
	writeXStringSet(x=ALNi_OUT,filepath=ALN_FA_MASKED_PATH)
}

chromlen <- function(TAB1_PATH,DICTIONARY_PATH){
	tab1.names <- read.table(TAB1_PATH,header=F,sep='\t') %>% select(V1) %>% unlist %>% unname
	dict       <- read.table(DICTIONARY_PATH,header=F,sep='\t') %>% select(V2,V3) %>% mutate(chrom.name=gsub('^SN:','',V2),chrom.length=gsub('^LN:','',V3))
	dict[match(tab1.names,dict[,'chrom.name']),'chrom.length'] %>% unname %>% as.data.frame
}


