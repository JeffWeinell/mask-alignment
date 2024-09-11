#!/bin/sh
#SBATCH --job-name batch-mask-alignment
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=20gb
#SBATCH --time=12:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=jweine2@gmail.com
#SBATCH --chdir=/home/jweinell/mendel-nas1/logs/
#SBATCH --output=slurm-%j-%x.out
#SBATCH --error=slurm-%j-%x.out

source ~/.bash_profile
export PATH=/home/jweinell/mendel-nas1/miniconda3/bin:$PATH

MAF_NAMES_FILE=${1}  # file with names of input .maf files
ALN_MAF_DIR=${2}     # directory containing the .maf files listed in $MAF_NAMES_FILE
ALN_FA_DIR_OUT=${3}  # output directory where masked fasta alignments are written
imin=${4}            # first .maf in $MAF_NAMES_FILE to processes
imax=${5}            # last .maf in $MAF_NAMES_FILE to processes

GENOMES_PATH="/shared_data/osn/herp-refdata/dna_assemblies/jeffweinell/Coronellini_WGS-assemblies_refsites-hardmasked_repeats-softmasked.fa"
R_PACKAGES_DIR="/home/jweinell/mendel-nas1/Rpackages/"

for i in $(seq $imin $imax);
do
	MAFi=$(sed "${i}q;d" $MAF_NAMES_FILE)
	FAi=$(echo "$MAFi" | sed 's|.maf$|.fa|g')
	ALN_MAF_PATH=$ALN_MAF_DIR$MAFi
	ALN_FAOUT_PATH=$ALN_FA_DIR_OUT$FAi
	[[ ! -f "$ALN_FAOUT_PATH" ]] && echo $i && bash /home/jweinell/mendel-nas1/prog/mask-alignment/mask-alignment.sh $GENOMES_PATH $ALN_MAF_PATH $ALN_FAOUT_PATH $R_PACKAGES_DIR
	wait
done

