#!/bin/sh

module load HAL/hal-2.1

# path to hal file
HAL_PATH=/shared_data/osn/herp-refdata/dna_alignments/Coronellini/lamps_genomes58_alignment_final.hal
# where to write chromosome lengths
CHROMLENGTHS_PATH=/home/jweinell/mendel-nas1/cactus/fasta/lamps_genomes58_alignment_final.hal.chromlengths

# names of all genomes in the HAL
GENOME_NAMES=$(halStats --genomes "$HAL_PATH" | sed 's| |\n|g' | sort)
# genome names excluding ancestral genomes
GENOME_TIPNAMES=$(echo "$GENOME_NAMES" | grep -v '^Anc[0-9]*')
# root genome name
ROOT_GENOME=$(halStats --root $HAL_PATH)
# number of tip genomes
NUMTIPGENOMES=$(echo "$GENOME_TIPNAMES" | wc -l)
for i in $(seq 1 $NUMTIPGENOMES);
do
	echo $i
	GENOMEi=$(echo "$GENOME_TIPNAMES" | sed "${i}q;d")
	[[ $i -eq 1 ]] && halStats --chromSizes $GENOMEi $HAL_PATH | awk -v gen=$GENOMEi '{print gen"."$0}' > $CHROMLENGTHS_PATH
	[[ $i -gt 1 ]] && halStats --chromSizes $GENOMEi $HAL_PATH | awk -v gen=$GENOMEi '{print gen"."$0}' >> $CHROMLENGTHS_PATH
done
