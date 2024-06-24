#!/bin/sh

# requires HAL (tested with hal-2.1)

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
  	
 	# (1) add $GENOMEi BED intervals (if any) to $BED_PATH
 	BEDPATHi=$(awk -v gen=$GENOMEi '$1==gen{print $2}' $BED_CONFIG_PATH)
	[[ -f "$BEDPATHi" ]] && [[ -f "$BED_CONFIG_PATH" ]] && awk -v gen=$GENOMEi '{print gen"."$0}' $BEDPATHi > $BED_PATH
 	[[ -f "$BEDPATHi" ]] && [[ ! -f "$BED_CONFIG_PATH" ]] && awk -v gen=$GENOMEi '{print gen"."$0}' $BEDPATHi > $BED_PATH

 	# (2) add $GENOMEi sequence lengths to $CHROMLENGTHS_PATH
   	[[ $i -eq 1 ]] && halStats --chromSizes $GENOMEi $HAL_PATH | awk -v gen=$GENOMEi '{print gen"."$0}' > $CHROMLENGTHS_PATH
    	[[ $i -gt 1 ]] && halStats --chromSizes $GENOMEi $HAL_PATH | awk -v gen=$GENOMEi '{print gen"."$0}' >> $CHROMLENGTHS_PATH

     	# (3) mask $GENOMEi sequences and add to $GENOMES_PATH
 	[[ $i -eq 1 ]] && bedtools maskfasta -fi <(hal2fasta --ucscSequenceNames $HAL_PATH $GENOMEi) -bed $BED_PATH > $GENOMES_PATH
	[[ $i -gt 1 ]] && bedtools maskfasta -fi <(hal2fasta --ucscSequenceNames $HAL_PATH $GENOMEi) -bed $BED_PATH >> $GENOMES_PATH
done


