# CaBagE: a Cas9-based Background Elimination strategy for targeted, long-read DNA sequencing

This is the repository for the laboratory protocol and analysis code in *CaBagE: a Cas9-based Background Elimination strategy for targeted, long-read DNA sequencing*; [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.10.13.337253v2)

## Data processing
Basecalling was performed with the MinKNOW software (v.19.05.0) using the Guppy flip-flop algorithm. Resulting fastq files were aligned to GRCh38 with MiniMap2 (v.2.14-r894-dirty)

`minimap2 -Yax map-ont $REFERENCE $FASTQ > $FILE.sam`

Alignments converted to sorted, indexed BAM format using samtools (v.1.10)

`samtools view -Sb $FILE.sam > $FILE.bam`
`samtools sort -@ 8 -m 1G $FILE.bam -o $FILE.sorted.bam`
`samtools index $FILE.sorted.bam`

## On-target Read Quantification

## Repeat Copy Number Distribution Histogram and Genotype Estimation

For CaBaGE targets containing short tandem repeats, the repeat copy number can be quantified with 

`repeat_estimator.py --bam #path to BAM file\
                     --ref #path to reference genome\
                     --locus #reference location of repeat (example: Chr9:27573480-27573551)\
                     --repeat_unit #example: "CAG"\
                     --alignment_buffer #length of query sequence flanking repeat (default 1000bp)\
                     --id #Sample ID\
                     --allele-count #expected number of alleles (default diploid)` 

after [installation](#Installation)

This script generates PLOT and estimate of STATs for cluster estimate with default copy number = 2


## Wet Lab Protocol

The wet lab protocol for target enrichment with CaBaGE is described in detail in the [Wiki](https://github.com/adw222/CaBagE-manuscript/wiki)


# Installation

To install all dependencies use conda:

```
conda create -y --name cabage --file requirements.txt -c bioconda -c anaconda
conda activate cabage
python repeat_example.py
```
