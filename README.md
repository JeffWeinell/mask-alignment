# Mask whole genome alignment in Multiple Alignment Format (MAF)

## About

Method to hard mask (with Ns) genome sequence intervals in a MAF file that was previously generated using the [cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/README.md) command `cactus-hal2maf`.

Steps include:
- (1) Edit the configuration files 'bedpaths.config' and 'settings.config' to define input/output file paths.
- (2) Run [mask-genomes.sh](https://github.com/JeffWeinell/mask-alignment/blob/main/mask-genomes.sh) to combine BEDs defined in 'bedpaths.config' into a single BED with USCS-format 'genome.chromosome' sequence names; extract and mask genomes from the MAF-precursor HAL and save all masked genomes (without alignment information) in a single fasta file that is used in the next (alignment masking) step.
- (3) Run [mask-MAF.sh](https://github.com/JeffWeinell/mask-alignment/blob/main/mask-MAF.sh) to generate a MAF alignment with sites masked at intervals in BED files defined in 'bedpaths.config'.

**Considerations**
- I have only tested this code on MAFs produced using the [cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/README.md) function `cactus-hal2maf`.
- Zero-length sequences, if present in the input MAF, will always be filtered

**External dependencies**

- [HAL](https://github.com/ComparativeGenomicsToolkit/hal/tree/master) (tested using version 2.1)
- [bedtools](https://bedtools.readthedocs.io/en/latest/) (tested using version 2.29.2)
- [seqkit](https://bioinf.shenwei.me/seqkit/) 
- R v4+ (tested using version 4.0.2)
- R packages `dplyr` (tested v1.0.7), `stringr` (tested v1.4.0), `Biostrings` (tested v2.58.0), and `GenomicRanges` (tested v1.42.0)

## Download these necessary files

- [bedpaths.config](https://github.com/JeffWeinell/mask-alignment/blob/main/bedpaths.config)
- [settings.config](https://github.com/JeffWeinell/mask-alignment/blob/main/settings.config)
- [mask-genomes.sh](https://github.com/JeffWeinell/mask-alignment/blob/main/mask-genomes.sh)
- [mask-MAF.sh](https://github.com/JeffWeinell/mask-alignment/blob/main/mask-MAF.sh)
- [mask-MAF.R](https://github.com/JeffWeinell/mask-alignment/blob/main/mask-MAF.R)


## Prepare input files and run code

**Step 1: Prepare input files**

Edit the two configuration files: [settings.config](https://github.com/JeffWeinell/mask-alignment/blob/main/settings.config) and [bedpaths.config](https://github.com/JeffWeinell/mask-alignment/blob/main/bedpaths.config).

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

# path to your copy of 'mask-MAF.R' script [default: current working directory]
MASK_ALIGNMENT_RSCRIPT=${PWD}/mask-MAF.R

# path to directory where your R packages are installed (leave empty to use R's default packages directory)
R_PACKAGES_DIR=
```

**Step 2: extract and mask genomes**

```
bash ./mask-genomes.sh ./settings.config
```



**Step 3: use masked genomes to mask alignment**

Run [mask-MAF.sh](https://github.com/JeffWeinell/mask-alignment/blob/main/mask-MAF.sh)

```
bash ./mask-MAF.sh ./settings.config
```


