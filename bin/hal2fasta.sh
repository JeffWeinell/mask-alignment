#!/bin/sh

module load HAL/hal-2.1

# path to input hal file
HAL_PATH=/shared_data/osn/herp-refdata/dna_alignments/Coronellini/lamps_genomes58_alignment_final.hal
# path to output fasta file
# GENOMES_PATH=/home/jweinell/mendel-nas1/cactus/fasta/lamps_genomes58-extracted-from-hal.fa
GENOMES_PATH=/shared_data/osn/herp-refdata/dna_alignments/Coronellini/lamps_genomes58-extracted-from-hal.fa

# names of all genomes in the HAL
GENOME_NAMES=$(halStats --genomes "$HAL_PATH" | sed 's| |\n|g' | sort)
# genome names excluding ancestral genomes
GENOME_TIPNAMES=$(echo "$GENOME_NAMES" | grep -v '^Anc[0-9]*')
# root genome name
ROOT_GENOME=$(halStats --root $HAL_PATH)
# number of tip genomes
NUMTIPGENOMES=$(echo "$GENOME_TIPNAMES" | wc -l)
# extract each tip genome as fasta and save to $GENOMES_PATH
for i in $(seq 3 $NUMTIPGENOMES);
do
	GENOMEi=$(echo "$GENOME_TIPNAMES" | sed "${i}q;d")
	echo $i $GENOMEi
	[[ $i -eq 1 ]] && hal2fasta --ucscSequenceNames $HAL_PATH $GENOMEi > $GENOMES_PATH
	[[ $i -gt 1 ]] && hal2fasta --ucscSequenceNames $HAL_PATH $GENOMEi >> $GENOMES_PATH
	sleep 1
done

# export all genomes (fasta format) descendant from $ROOT_GENOME
# hal2fasta --ucscSequenceNames 1 --subtree 1 --outFaPath $GENOMES_PATH $HAL_PATH $ROOT_GENOME





