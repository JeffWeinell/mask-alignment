# mask-alignment

Method to hard mask genome-specific sites in a MAF-format whole genome alignment that was generated using the [cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/README.md) function `cactus-hal2maf`.

This method involves two main steps:
- (1) Extract genomes from the HAL alignment file (which was used to generate the MAF), mask intervals, and save masked genomes in a single fasta file.
- (2) Use original MAF (from cactus-hal2maf) and the masked genomes fasta (from step 1) to generate a new MAF-format alignment with sites masked according to the BED intervals from step 1


### considerations
- I have only tested this code on MAFs produced using the [cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/README.md) function `cactus-hal2maf`
- Zero-length sequences in input MAF will always be filtered


#### dependencies
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
bash /path/to/mask-alignment.sh /path/to/settings.config
```

#### Config file

Download and edit [settings.config](https://github.com/JeffWeinell/mask-alignment/blob/main/settings.config) before running `mask-alignment.sh`

The unedited `settings.config` file looks like this:

```
### Required
HAL_PATH=                       # path to HAL alignment file from which the MAF is derived (step 1 input)
BED_PATH=                       # path to BED with intervals to mask (step 1 input)
GENOMES_PATH=                   # path to fasta containing a MASKED version of each genome extracted from the HAL alignment file (step 1 output; step 2 input)

ALN_MAF_PATH=                   # path to input MAF previously generated by 'cactus-hal2maf' (step 2 input)
MAF_OUT_PATH=                   # where to save output MAF (step 2 output)
MASK_ALIGNMENT_RSCRIPT=         # path to your copy of 'mask-alignment.R'

### Optional
R_PACKAGES_DIR=                 # path to directory where your R packages are installed (R default location used if left blank)
DICTIONARY_PATH=                # path to samtools dictionary file (default $GENOMES_PATH".dict")
```








