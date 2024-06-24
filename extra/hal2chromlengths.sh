#!/bin/sh

# requires HAL (testing with hal-2.1)

# path to settings.config
SETTINGS_PATH=${1}

# load settings (input/output paths) from settings.config 
source $SETTINGS_PATH

# use default path for $CHROMLENGTHS_PATH if not explicitely set in settings.config file
[[ $(echo $CHROMLENGTHS_PATH | wc -w ) -eq 0 ]] && CHROMLENGTHS_PATH=$HAL_PATH".seqlengths"

# names of all genomes in the HAL
GENOME_NAMES=$(halStats --genomes "$HAL_PATH" | sed 's| |\n|g' | sort)
# genome names excluding ancestral genomes
GENOME_TIPNAMES=$(echo "$GENOME_NAMES" | grep -v '^Anc[0-9]*')
# number of tip genomes
NUMTIPGENOMES=$(echo "$GENOME_TIPNAMES" | wc -l)
# extract chromosome names and lengths for each tip genome and save (using UCSC format sequence names) to $GENOMES_PATH 
for i in $(seq 1 $NUMTIPGENOMES);
do
	echo $i
	GENOMEi=$(echo "$GENOME_TIPNAMES" | sed "${i}q;d")
	[[ $i -eq 1 ]] && halStats --chromSizes $GENOMEi $HAL_PATH | awk -v gen=$GENOMEi '{print gen"."$0}' > $CHROMLENGTHS_PATH
	[[ $i -gt 1 ]] && halStats --chromSizes $GENOMEi $HAL_PATH | awk -v gen=$GENOMEi '{print gen"."$0}' >> $CHROMLENGTHS_PATH
done
