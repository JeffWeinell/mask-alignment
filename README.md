# Mask whole genome alignment in Multiple Alignment Format (MAF)

## Overview

Method to hard mask (with Ns) genome sequence intervals in a MAF file that was previously generated using the [cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/README.md) command `cactus-hal2maf`.

**Steps:**
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

- [bedpaths.config](https://raw.githubusercontent.com/JeffWeinell/mask-alignment/main/bedpaths.config)
- [settings.config](https://raw.githubusercontent.com/JeffWeinell/mask-alignment/main/settings.config)
- [mask-genomes.sh](https://raw.githubusercontent.com/JeffWeinell/mask-alignment/main/mask-genomes.sh)
- [mask-MAF.sh](https://raw.githubusercontent.com/JeffWeinell/mask-alignment/main/mask-MAF.sh)
- [mask-MAF.R](https://raw.githubusercontent.com/JeffWeinell/mask-alignment/main/mask-MAF.R)


## Prepare input files and run code

**Step 1: Prepare input files**

Edit the two configuration files: **settings.config** and **bedpaths.config**.

**Step 2: extract and mask genomes**

```
bash ./mask-genomes.sh ./settings.config
```

**Step 3: use masked genomes to mask alignment**

```
bash ./mask-MAF.sh ./settings.config
```


