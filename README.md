# CaBagE: a Cas9-based Background Elimination strategy for targeted, long-read DNA sequencing

This is the repository for the laboratory protocol and analysis code in *CaBagE: a Cas9-based Background Elimination strategy for targeted, long-read DNA sequencing*; [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.10.13.337253v2)

## Data processing
Basecalling was performed with the MinKNOW software (v.19.05.0) using the Guppy flip-flop algorithm. Resulting fastq files were aligned to GRCh38 with MiniMap2 (v.2.14-r894-dirty)

`minimap2 -Yax map-ont $REFERENCE $FASTQ > $FILE.sam`

Alignments converted to sorted, indexed BAM format using samtools (v.1.10)
```
    samtools view -Sb $FILE.sam > $FILE.bam
    samtools sort -@ 8 -m 1G $FILE.bam -o $FILE.sorted.bam
    samtools index $FILE.sorted.bam
```
## On-Target Read Quantification
On-Target reads are defined as any read (MapQ >=60) that overlaps the region flanked by guideRNAs by at least one base pair:

`samtools view -c -q 60 $FILE.sorted.bam [*chr:start-end*]`

Spanning reads are defined as any read that overlaps the region flanked by guideRNAs by at least 90%:

`bedtools coverage -f 0.9 -a $REGIONS.bed -b $FILE.sorted.bam -counts`

## Repeat Copy Number Distribution Histogram and Genotype Estimation

For CaBaGE targets containing short tandem repeats, the repeat copy number can be quantified with 
```
repeat_estimator.py  --bam \ #path to sorted BAM file
    --ref \ #path to reference genome
    --locus \ #reference location of repeat (example: chr9:27573480-27573551)
    --repeat_unit \ #example: "CAG"
    --alignment_buffer \ #length of query sequence flanking repeat (default 1000bp)
    --id \ #sample ID
    --allele-count #expected number of alleles (default diploid)
```

after [installation](#Installation)

This script generates a histogram showing the distribution of repeat copy number across CaBagE enriched sequencing reads and also estimates copy number genotypes using a Gaussian Mixure Model with default copy number = 2


## Wet Lab Protocol

The wet lab protocol for target enrichment with CaBaGE is described in detail in the [Wiki](https://github.com/adw222/CaBagE-manuscript/wiki)


## Installation

To install all dependencies use conda:

```
conda create -y --name cabage --file requirements.txt -c bioconda -c anaconda
conda activate cabage
python repeat_estimator.py
```
### Dependancies

    scipy 1.4.1
    pyfaidx 0.5.9.1
    numpy 1.17.3
    pandas 1.1.5
    matplotlib 3.3.3
    scikit-learn 0.23.2
    scikit-bio 0.5.6
    biopython 1.78
    pysam 0.16.0.1
