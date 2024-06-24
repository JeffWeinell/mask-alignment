# Mask whole genome alignment in Multiple Alignment Format (MAF)

Method to hard mask genome-specific sites in a MAF-format whole genome alignment that was generated using the [cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/README.md) function `cactus-hal2maf`.

This method involves two main steps:
- (1) Extract genomes from the HAL alignment file (which was used to generate the MAF), mask intervals, and save masked genomes in a single fasta file.
- (2) Use original MAF (from cactus-hal2maf) and the masked genomes fasta (from step 1) to generate a new MAF-format alignment with sites masked according to the BED intervals from step 1


### considerations
- I have only tested this code on MAFs produced using the [cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/README.md) function `cactus-hal2maf`
- Zero-length sequences in input MAF will always be filtered

**Dependencies**

You must install these programs ahead of time:
- HAL (for step 1)
- bedtools (tested using version 2.29.2)
- samtools (tested using version 1.11)
- seqkit 
- R v4+
- R packages `dplyr`, `stringr`, `Biostrings`, and `GenomicRanges`
- [settings.config](https://github.com/JeffWeinell/mask-alignment/blob/main/settings.config) download and edit to specify settings for `mask-alignment.sh`
- R script [mask-alignment.R](https://github.com/JeffWeinell/mask-alignment/blob/main/bin/mask-alignment.R) (download but don't edit)

#### How to use this code

**Step 1**




**Step 2**

The main script to run for step 2 is [mask-alignment.sh](https://github.com/JeffWeinell/mask-alignment/blob/main/bin/mask-alignment.sh)

```
# to run:
bash ./mask-alignment.sh ./settings.config
```

**Update settings.config file**

Download and edit [settings.config](https://github.com/JeffWeinell/mask-alignment/blob/main/settings.config) before running `mask-alignment.sh`

The unedited `settings.config` file looks like this:

```
######## required input/output paths

# path to input MAF
ALN_MAF_PATH=

# path to HAL alignment file that was used to generate the input MAF
HAL_PATH=                       

# path to two-column, tab-separated file with genome name (column 1) and path to genome mask intervals BED file (column 2)
BED_CONFIG_PATH=                

# where to save output (masked) MAF
MAF_OUT_PATH=

######## optional settings

# where to save combined BED file with all genome sequence masking intervals
BED_PATH=${HAL_PATH}_all-genomes-mask-intervals.bed                       

# where to save masked-versions of unaligned genomes extracted from HAL (combined in a single fasta file)
GENOMES_PATH=${HAL_PATH}-masked-nogaps.fa

# path to two-column tab-separated file with sequence names (all genomes; USCS format names) and lengths
CHROMLENGTHS_PATH=${HAL_PATH}.seqlengths

# path to your copy of 'mask-alignment.R' script [default: current working directory]
MASK_ALIGNMENT_RSCRIPT=${PWD}/mask-alignment.R

# path to directory where your R packages are installed (leave empty to use R's default packages directory)
R_PACKAGES_DIR=
```









