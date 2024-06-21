# mask-alignment

#### About
Method to hard mask genome-specific sites in an existing whole genome alignment.

This is a work in progress, but at present:
- Input alignment is a Multiple Alignment Format (MAF).
- Output alignment is a Fasta alignment
- Additional required input = fasta file with replacement genomes containing desired masking

#### dependencies
- bedtools (tested using version 2.29.2)
- R v4+
- R packages `dplyr`, `stringr`, `Biostrings`, and `GenomicRanges`

#### Usage

```
/path/to/mask-alignment.sh /path/to/settings.config
```

#### Config file

Edit the `settings.config` file to set the required variables...






