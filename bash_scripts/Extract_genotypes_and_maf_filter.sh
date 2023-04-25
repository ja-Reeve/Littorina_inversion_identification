#!/bin/bash

### This script extracts the CHROM, POS and GT fields from a VCF. It returns a .txt with a seperate column for each individual in the VCF.
### James Reeve - Göteborgs Universitet
### 11/10/2020

### Job parameters
#SBATCH -A snic2020-15-103
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -J Calling_Genotypes

### Required software - loaded as modules
module load bioinfo-tools
#module load vcftools/0.1.16
module load bcftools/1.10

### Filepath
DIR=/proj/snic2020-15-103/

### The command
# LG1
vcftools --vcf $DIR//filtered_VCFs/LGs_seperately/LG1_depth_quality.recode.vcf --out $DIR/LG1_depth_quality_maf --maf 0.05 --recode
bcftools query -f '%CHROM %POS [ %GT]\n' $DIR/LG1_depth_quality_maf.recode.vcf > $DIR/genotypes/LG1_genotypes.v2.txt

# LG2
vcftools --vcf $DIR//filtered_VCFs/LGs_seperately/LG2_depth_quality.recode.vcf --out $DIR/LG2_depth_quality_maf --maf 0.05 --recode
bcftools query -f '%CHROM %POS [ %GT]\n' $DIR/LG2_depth_quality_maf.recode.vcf > $DIR/genotypes/LG2_genotypes.v2.txt

### Repeat the last two sections for each LG in the genome.
