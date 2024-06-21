#!/bin/sh

# load settings from settings.config
SETTINGS_PATH=${1}
source $SETTINGS_PATH

# exit if output file already exists
[[ -f "$ALN_FAOUT_PATH" ]] && echo "exiting because "$ALN_FAOUT_PATH" already exists" && exit 0

# create temporary files
ALN_FA_UNMASKED_PATH=$(mktemp 2>&1)
REGIONS_TABLE_PATH_A=$(mktemp 2>&1)
REGIONS_TABLE_PATH_B=$(mktemp 2>&1)
REGIONS_TABLE_PATH=$(mktemp 2>&1)
BEDPATH=$(mktemp 2>&1)
GAPLESS_FA_PATH=$(mktemp 2>&1)

# convert MAF alignment to fasta alignment
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

# extract sequences in BED intervals from input masked genomes file
bedtools getfasta -s -nameOnly -fi $GENOMES_PATH -bed $BEDPATH | sed 's|[(][+-][)]$||g' > $GAPLESS_FA_PATH

# For each sequence in the original alignment, use R vector logic to replace characters at non-gap sites with characters in the same-interval masked-version of sequence.
Rscript $MASK_ALIGNMENT_RSCRIPT $ALN_FA_UNMASKED_PATH $GAPLESS_FA_PATH $ALN_FAOUT_PATH $R_PACKAGES_DIR
echo "masked alignment written to: "$ALN_FAOUT_PATH

### remove temporary files
rm $ALN_FA_UNMASKED_PATH
rm $REGIONS_TABLE_PATH_A
rm $REGIONS_TABLE_PATH_B
rm $REGIONS_TABLE_PATH
rm $BEDPATH
rm $GAPLESS_FA_PATH






