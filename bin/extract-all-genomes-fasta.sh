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
# extract each tip genome as fasta with UCSC format sequence names and save to $GENOMES_PATH
for i in $(seq 1 $NUMTIPGENOMES);
do
	GENOMEi=$(echo "$GENOME_TIPNAMES" | sed "${i}q;d")
	echo $i $GENOMEi
	[[ $i -eq 1 ]] && hal2fasta --ucscSequenceNames $HAL_PATH $GENOMEi > $GENOMES_PATH
	[[ $i -gt 1 ]] && hal2fasta --ucscSequenceNames $HAL_PATH $GENOMEi >> $GENOMES_PATH
	sleep 1
done
