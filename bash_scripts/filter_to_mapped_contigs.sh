#!/bin/bash -ve
### Script written by Sean Stankowski for ShARC cluster
### Note: the same script was used to filter by LG in the final manuscript.
### This is simply done by reducing the list of chromosomes to just those
### on the focal linakge group.
### 2021-05-15


#$ -P littorina
#$ -q littorina.q
#$ -pe smp 8
# request memory for job (default 6G, max 72G)
#$ -l mem=10G
#$ -l rmem=10G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=336:00:00


module load apps/vcftools

DIR=/shared/snooklin1/GroupA/shared_wgs/filtered_VCFs/mapped_contigs_only_depth_and_quality

vcftools --vcf $DIR/mapped_contigs_full_depth5-35_biallelic_minQ20.recode.vcf  \
         --recode \
         --out $DIR/inversions_only_mapped_contigs_full_depth5-35_biallelic_minQ20 \
         --chr "<<< contig #1 >>>" --chr "<<< contig #2 >>>" --chr "<<< contig #3 >>>"
