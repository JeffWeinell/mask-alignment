#!/bin/sh

# requires HAL (v2.1 tested) and bedtools (v2.29.2 tested)

# path to settings.config
SETTINGS_PATH=${1}

# load settings (input/output paths) from settings.config 
source $SETTINGS_PATH

# names of all genomes in the HAL
GENOME_NAMES=$(halStats --genomes "$HAL_PATH" | sed 's| |\n|g' | sort)
# genome names excluding ancestral genomes
GENOME_TIPNAMES=$(echo "$GENOME_NAMES" | grep -v '^Anc[0-9]*')
# number of tip genomes
NUMTIPGENOMES=$(echo "$GENOME_TIPNAMES" | wc -l)
 
# (1) get bed intervals to mask for each genome and add to $BED_PATH
# (2) get sequence lengths for each genome and append to $CHROMLENGTHS_PATH
# (3) extract each genome and save to $GENOMES_PATH as fasta with UCSC format sequence names
for i in $(seq 1 $NUMTIPGENOMES);
do
	GENOMEi=$(echo "$GENOME_TIPNAMES" | sed "${i}q;d")
	echo $i $GENOMEi

 	# test if any intervals to mask in $GENOMEi (0=false, 1=true)
   	TEST_BEDPATHi=$(awk -v gen=$GENOMEi '$1==gen{print $1}' $BED_CONFIG_PATH | wc -w)
	[[ "$TEST_BEDPATHi" -gt 1 ]] && echo "genome names in must be unique in "$BED_CONFIG_PATH && exit 1
     
 	# (1) add $GENOMEi BED intervals (if any) to $BED_PATH
   	[[ "$TEST_BEDPATHi" -eq 0 ]] && BEDPATHi=""
    	[[ "$TEST_BEDPATHi" -eq 1 ]] && BEDPATHi=$(awk -v gen=$GENOMEi '$1==gen{print $2}' $BED_CONFIG_PATH)
	[[ "$TEST_BEDPATHi" -eq 1 ]] && [[ -f "$BEDPATHi" ]] && [[ -f "$BED_CONFIG_PATH" ]] && awk -v gen=$GENOMEi '{print gen"."$0}' $BEDPATHi > $BED_PATH
 	[[ "$TEST_BEDPATHi" -eq 1 ]] && [[ -f "$BEDPATHi" ]] && [[ ! -f "$BED_CONFIG_PATH" ]] && awk -v gen=$GENOMEi '{print gen"."$0}' $BEDPATHi > $BED_PATH

 	# (2) add $GENOMEi sequence names (genomeName.chromosomeName) and lengths to $CHROMLENGTHS_PATH
   	[[ $i -eq 1 ]] && halStats --chromSizes $GENOMEi $HAL_PATH | awk -v gen=$GENOMEi '{print gen"."$0}' > $CHROMLENGTHS_PATH
    	[[ $i -gt 1 ]] && halStats --chromSizes $GENOMEi $HAL_PATH | awk -v gen=$GENOMEi '{print gen"."$0}' >> $CHROMLENGTHS_PATH

     	# (3) mask $GENOMEi sequences and add to $GENOMES_PATH
 	[[ "$TEST_BEDPATHi" -eq 0 ]] && [[ $i -eq 1 ]] && hal2fasta --ucscSequenceNames $HAL_PATH $GENOMEi > $GENOMES_PATH
  	[[ "$TEST_BEDPATHi" -eq 0 ]] && [[ $i -gt 1 ]] && hal2fasta --ucscSequenceNames $HAL_PATH $GENOMEi >> $GENOMES_PATH
   	
    	GENOMEi_TEMP_PATH=$(mktemp 2>&1)
     	hal2fasta --ucscSequenceNames $HAL_PATH $GENOMEi > $GENOMEi_TEMP_PATH
  	[[ "$TEST_BEDPATHi" -eq 1 ]] && [[ $i -eq 1 ]] && bedtools maskfasta -fi <(hal2fasta --ucscSequenceNames $HAL_PATH $GENOMEi) -bed $BED_PATH > $GENOMES_PATH
	[[ "$TEST_BEDPATHi" -eq 1 ]] && [[ $i -gt 1 ]] && bedtools maskfasta -fi <(hal2fasta --ucscSequenceNames $HAL_PATH $GENOMEi) -bed $BED_PATH >> $GENOMES_PATH
	rm $GENOMEi_TEMP_PATH
done


