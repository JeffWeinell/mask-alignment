#!/bin/sh

# load settings from settings.config
SETTINGS_PATH=${1}
source $SETTINGS_PATH

# exit if output file already exists
[[ -f "$MAF_OUT_PATH" ]] && echo "exiting because "$MAF_OUT_PATH" already exists" && exit 0

# # create temporary samtools sequence dictionary file for genomes if one doesnt already exist
# [[ $(echo $DICTIONARY_PATH | wc -w) -eq 0 ]] && DICTIONARY_PATH=$GENOMES_PATH".dict"
# [[ -f "$DICTIONARY_PATH" ]] && RMDICT=1
# [[ ! -f "$DICTIONARY_PATH" ]] && RMDICT=0 && DICTIONARY_PATH=$(mktemp 2>&1) && samtools dict --no-header --output $DICTIONARY_PATH $GENOMES_PATH

# create temporary files
ALN_FA_UNMASKED_PATH=$(mktemp 2>&1)
REGIONS_TABLE_PATH_A=$(mktemp 2>&1)
REGIONS_TABLE_PATH_B=$(mktemp 2>&1)
REGIONS_TABLE_PATH=$(mktemp 2>&1)
BEDPATH=$(mktemp 2>&1)
GAPLESS_FA_PATH=$(mktemp 2>&1)
ALN_FA_MASKED_PATH=$(mktemp 2>&1)
TAB1_PATH=$(mktemp 2>&1)
MAFLIKE_OUT_PATH=$(mktemp 2>&1)

# convert input MAF alignment into fasta alignment
awk '$1=="s"{print ">"$2"("$5")/"$3+1"-"$3+$4+1"\n"$7}' $ALN_MAF_PATH > $ALN_FA_UNMASKED_PATH

# make a table with genomic intervals for each sequence
ALNFILEi=$(basename "$ALN_MAF_PATH")
zcat -f $ALN_MAF_PATH | awk -F'\t' -v x=$ALNFILEi '$1=="s"{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"x}' | awk '{print $0"\t"$2+$3}' | sed 's|\(^[^.]*\)[.]\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*$\)|\1.\2\t\1\t\2\t\3\t\4\t\5\t\6\t\7\t\8\t\1-\2|g' | awk '{print $0"("$6")/"$4+1"-"$9+1"\tgene\t.\t0"}' | awk '{print $1"\t"$13"\t"$7"\t"$1"\t"$12"\t"$6"\t"$4"\t"$9"\t"$5"\t"$10"\t"$2"\t"$3"\t"$8"\t"$11}' > $REGIONS_TABLE_PATH_A
zcat -f $ALN_MAF_PATH | awk '$1=="s"{print $2"("$5")/"$3+1"-"$3+$4+1}' > $REGIONS_TABLE_PATH_B
paste $REGIONS_TABLE_PATH_A $REGIONS_TABLE_PATH_B > $REGIONS_TABLE_PATH

# convert the table from the previous step into a BED-format table
awk '{print $1"\t"$3"\t"$6"\t"$7"\t"$8"\t"$15}' $REGIONS_TABLE_PATH | 
    awk '$3=="-" { $7 = $2-$5 }1' | 
    awk '$3=="-" { $8 = $2-$4 }1' | 
    awk '$3=="+" { $7=$4 }1' | 
    awk '$3=="+" { $8=$5 }1' |
    awk '$7!=$8{print $1"\t"$7"\t"$8"\t"$6"\t.\t"$3}' > $BEDPATH

# extract sequences in BED intervals from masked genomes file
bedtools getfasta -s -nameOnly -fi $GENOMES_PATH -bed $BEDPATH | sed 's|[(][+-][)]$||g' > $GAPLESS_FA_PATH

REPLACE_NONGAPS_RSCRIPT="
args <- commandArgs(trailingOnly = TRUE)
if(args[4]!=''){
.libPaths(args[4])
}
source(args[5])
replace_nongaps_with(args[1],args[2],args[3])
"

# For each sequence in the original alignment, use R vector logic to replace characters at non-gap sites with same-interval sequence extracted from masked genome
Rscript <(echo "$REPLACE_NONGAPS_RSCRIPT") $ALN_FA_UNMASKED_PATH $GAPLESS_FA_PATH $ALN_FA_MASKED_PATH $R_PACKAGES_DIR $MASK_ALIGNMENT_RSCRIPT

# convert masked fasta alignment into a seqkit sequence table
seqkit fx2tab $ALN_FA_MASKED_PATH | sed 's|(+)|\t+\t|g' | sed 's|(-)|\t-\t|g' | sed 's|\t/|\t|g' > $TAB1_PATH

# identify reference genome name if not specified in settings.config (reference genome is expected to be first in each alignment block)
[[ $(echo $REFERENCE_GENOME | wc -w) -eq 0 ]] && REFERENCE_GENOME=$(awk 'NR==1{print $1}' $TAB1_PATH | sed -E 's|[.].+||g')

### convert sequence table to MAF
# R script to get chromosome lengths for each sequence in MAF
CHROMLEN_RSCRIPT="
args <- commandArgs(trailingOnly=TRUE)
if(args[3]!=''){
.libPaths(args[3])
}
source(args[4])
chromlen(args[1],args[2])
"

# construct MAF-like table
CHROM=$(awk '{print $1}' $TAB1_PATH)
STARTi1=$(awk '{print $3}' $TAB1_PATH | sed -E 's|-.+||g')
ENDi1=$(awk '{print $3}' $TAB1_PATH | sed -E 's|^.+-||g')
STARTi0=$( echo "$STARTi1" | awk '{print $1-1}')
SEQLEN=$(paste <(echo "$STARTi1") <(echo "$ENDi1") | awk '{print $2-$1}')
SEQSTRAND=$(awk '{print $2}' $TAB1_PATH)
SEQSTRING=$(awk '{print $4}' $TAB1_PATH)
CHROMLEN=$(Rscript <(echo "$CHROMLEN_RSCRIPT") $TAB1_PATH $CHROMLENGTHS_PATH $R_PACKAGES_DIR $MASK_ALIGNMENT_RSCRIPT | awk 'NR>1{print $2}')
paste <(echo "$CHROM") <(echo "$STARTi0") <(echo "$SEQLEN") <(echo "$SEQSTRAND") <(echo "$CHROMLEN") <(echo "$SEQSTRING") | awk '{print "s\t"$0}' > $MAFLIKE_OUT_PATH

# add empty line and '^a$' line before each alignment block and maf header on line 1
sed -i "/$REFERENCE_GENOME/i d\na" $MAFLIKE_OUT_PATH
sed -i 's|^d$||g' $MAFLIKE_OUT_PATH

cat <(echo '##maf version=1') $MAFLIKE_OUT_PATH <(echo "") > $MAF_OUT_PATH

echo "alignment written to: "$MAF_OUT_PATH

### remove temporary files
rm $ALN_FA_UNMASKED_PATH
rm $REGIONS_TABLE_PATH_A
rm $REGIONS_TABLE_PATH_B
rm $REGIONS_TABLE_PATH
rm $BEDPATH
rm $GAPLESS_FA_PATH
rm $TAB1_PATH
rm $MAFLIKE_OUT_PATH
# [[ "$RMDICT" -eq 0 ]] && rm $DICTIONARY_PATH

