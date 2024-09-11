#!/bin/sh
#SBATCH --job-name mask-alignment
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=20gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=jweine2@gmail.com
#SBATCH --chdir=/home/jweinell/mendel-nas1/logs/
#SBATCH --output=slurm-%j-%x.out
#SBATCH --error=slurm-%j-%x.out

### if using mendel with modules R and bedtools modules
module load R/R-4.0.2
module load Bedtools/bedtools-2.29.2

#### If R and bedtools modules not available
# source ~/.bash_profile
# export PATH=/home/jweinell/mendel-nas1/miniconda3/bin:$PATH
# conda create -n Rbedtools R
# conda install -n Rbedtools -c bioconda bedtools
# conda activate Rbedtools

# capture variables from stdin
GENOMES_PATH=${1}      # path to input masked genome assemblies combined in a single file and with sample names included in sequence headers
ALN_MAF_PATH=${2}      # path to input MAF produced with cactus-hal2maf
ALN_FAOUT_PATH=${3}    # where to save output fasta

# Optional variables captured
R_PACKAGES_DIR=${4:-""} # path to directory where your R packages are installed (if other than default)

# Path to the R script file 'mask-alignment.R'
MASK_ALIGNMENT_RSCRIPT=/home/jweinell/mendel-nas1/prog/mask-alignment/mask-alignment.R

# NOTE: zero-length genomic regions in MAF are always filtered

## In the future this script will use inputs \\
# (1) HAL file \\
# (2) MAF previously generated from the HAL \\
# (3) BED with intervals to mask in output alignment fasta.

### exit if $ALN_FAOUT_PATH already exists
[[ -f "$ALN_FAOUT_PATH" ]] && echo "exiting because "$ALN_FAOUT_PATH" already exists" && exit 0

### temporary files
ALN_FA_UNMASKED_PATH=$(mktemp 2>&1)
REGIONS_TABLE_PATH_A=$(mktemp 2>&1)
REGIONS_TABLE_PATH_B=$(mktemp 2>&1)
REGIONS_TABLE_PATH=$(mktemp 2>&1)
BEDPATH=$(mktemp 2>&1)
GAPLESS_FA_PATH=$(mktemp 2>&1)

# convert MAF to fasta
awk '$1=="s"{print ">"$2"("$5")/"$3+1"-"$3+$4+1"\n"$7}' $ALN_MAF_PATH > $ALN_FA_UNMASKED_PATH

# make a table with genomic coordinates of each sequence regions and some other info
ALNFILEi=$(basename "$ALN_MAF_PATH")
zcat -f $ALN_MAF_PATH | awk -F'\t' -v x=$ALNFILEi '$1=="s"{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"x}' | awk '{print $0"\t"$2+$3}' | sed 's|\(^[^.]*\)[.]\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*\)\t\([^[:blank:]]*$\)|\1.\2\t\1\t\2\t\3\t\4\t\5\t\6\t\7\t\8\t\1-\2|g' | awk '{print $0"("$6")/"$4+1"-"$9+1"\tgene\t.\t0"}' | awk '{print $1"\t"$13"\t"$7"\t"$1"\t"$12"\t"$6"\t"$4"\t"$9"\t"$5"\t"$10"\t"$2"\t"$3"\t"$8"\t"$11}' > $REGIONS_TABLE_PATH_A
zcat -f $ALN_MAF_PATH | awk '$1=="s"{print $2"("$5")/"$3+1"-"$3+$4+1}' > $REGIONS_TABLE_PATH_B
paste $REGIONS_TABLE_PATH_A $REGIONS_TABLE_PATH_B > $REGIONS_TABLE_PATH

# make a BED-format version of the table created in the previous step
awk '{print $1"\t"$3"\t"$6"\t"$7"\t"$8"\t"$15}' $REGIONS_TABLE_PATH | 
    awk '$3=="-" { $7 = $2-$5 }1' | 
    awk '$3=="-" { $8 = $2-$4 }1' | 
    awk '$3=="+" { $7=$4 }1' | 
    awk '$3=="+" { $8=$5 }1' |
    awk '$7!=$8{print $1"\t"$7"\t"$8"\t"$6"\t.\t"$3}' > $BEDPATH

# extract BED intervals from $GENOMES_PATH containing masked-versions of sequences spanning the same regions as in $ALN_MAF_PATH (but without gaps)
bedtools getfasta -s -nameOnly -fi $GENOMES_PATH -bed $BEDPATH | sed 's|[(][+-][)]$||g' > $GAPLESS_FA_PATH

# For each sequence in $ALN_FA_UNMASKED_PATH, replace non-gap sites with corresponding replacement sequence in $GAPLESS_FA_PATH
Rscript $MASK_ALIGNMENT_RSCRIPT $ALN_FA_UNMASKED_PATH $GAPLESS_FA_PATH $ALN_FAOUT_PATH $R_PACKAGES_DIR
echo "masked alignment written to: "$ALN_FAOUT_PATH

### remove temporary files
rm $ALN_FA_UNMASKED_PATH
rm $REGIONS_TABLE_PATH_A
rm $REGIONS_TABLE_PATH_B
rm $REGIONS_TABLE_PATH
rm $BEDPATH
rm $GAPLESS_FA_PATH

module unload R/R-4.0.2
module unload Bedtools/bedtools-2.29.2

