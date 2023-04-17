### Filtering genotypes ###
### This script filters the genotype files from a WGS dataset. There are two
### filters applied; a filter of monomorphic sites per region and a filter of
### allele counts per region.
### James Reeve - University of Gothenburg
### 12/03/2021
### Edited 21/02/2023: adjust to modified sample metadata
### and simplify regional clustering.

### 1. Preparation ####
### Clean-up environment
rm(list = ls())
dev.off()
setwd("/Users/james/Documents/Inversion_detection")
options(stringsAsFactors = FALSE)

### Packages
library("tidyverse")

### Set linkage group
LG.n <- 17
# Change to whichever LG you want to analyse
LG.s <- paste0("LG", LG.n)

### Linkage map
LG <- read.table("ConsensusMap.v2.txt", header = T)
LG <- LG[LG$LG == LG.s,]

### Upload sample information
sample_info <- read.csv("Seans_WGS_sample_info_v2.csv")

### Genotypes
GT <- read.table(paste0("Genotype_tables/",LG.s,"_genotypes.v2.txt"), header = FALSE)
# Add header from VCF file
colnames(GT) <- names(read.delim("vcf_header.txt", check.names = FALSE))[c(1:2, 10:117)]
colnames(GT)[1] <- "CHROM" #rename 1st column
# Rename 1 odd Larc ID
names(GT)[names(GT) == "W_arc_04_Lamerged_sorted.bam"] <- "W_arc_04_La"


### 2. Separate IDs for different genetic groups ####
### NS = Snails north of Bay of Biscay
### Sp = Snails south of Bay of Biscay
### Last two groups are sister species
### arc = Littorina arcana
### comp = Littorina compressa

### Extract sample IDs from sample info file
# List samples in each species
sax_IDs <- sample_info[sample_info$Species == "saxatilis", "Sample_ID"]
arc_IDs <- sample_info[sample_info$Species == "arcana", "Sample_ID"]
comp_IDs <- sample_info[sample_info$Species == "compressa", "Sample_ID"]

# Separate saxatilis genetic groups
sax_info <- sample_info[sample_info$Sample_ID %in% sax_IDs, ]

NS_IDs <- sax_info[sax_info$Country != "Spain", "Sample_ID"] 
Sp_IDs <- sax_info[sax_info$Country == "Spain", "Sample_ID"]

# Final grouping for snails in sites where arcana and saxatilis were conspecific
# Note: The Amble site was removed as only arcana was collected here
arc_sax_IDs <- sample_info[sample_info$location %in% unique(sample_info[sample_info$Species == "arcana", "location"]) &
                             sample_info$Species != "compressa" & sample_info$location != "Amble", "Sample_ID"]

### 3. Filtering function ####
### Filter for minor allele counts
filter.allele.count <- function(data, min.allele.count){
  # Select a row and return TRUE/FALSE if it has > min.allele.count
  sub.function <- function(i){
    tmp <- i[-c(1:2)] # remove "CHROM" & "POS"
    tmp2 <- tmp[tmp != "./."] # Exclude missing sites
    
    # Allele counts
    N_ref <- sum(tmp2 == "0/0")
    N_het <- sum(tmp2 == "0/1")
    N_alt <- sum(tmp2 == "1/1")
    
    # Is allele count > min.allele.count
    return((2*N_ref + N_het) & (2*N_alt + N_het) > min.allele.count)
  }
  # Boolean vector output (i.e. TRUE / FALSE)
  keeps <- apply(data, 1, sub.function)
  
  # Remove SNPs that fail the filter
  res <- data[keeps,]
  return(res)
}

### 4. Filter the genotype table ####
### Filter to contigs on LG
GT <- GT[GT$CHROM %in% unique(LG$contig),]

### Split genotypes by genetic group
GT_NS <- GT %>% select("CHROM", "POS", all_of(NS_IDs))
GT_Sp <- GT %>% select("CHROM", "POS", all_of(Sp_IDs))
GT_arc <- GT %>% select("CHROM", "POS", all_of(arc_IDs))
GT_comp <- GT %>% select("CHROM", "POS", all_of(comp_IDs))
GT_arc_sax <- GT %>% select("CHROM", "POS", all_of(arc_sax_IDs))

### Filter by minor allele count
min_count <- 2 # Minor allele count threshold

GT_NS <- filter.allele.count(GT_NS, min_count)
GT_Sp <- filter.allele.count(GT_Sp, min_count)
GT_arc <- filter.allele.count(GT_arc, min_count)
GT_comp <- filter.allele.count(GT_comp, min_count)
GT_arc_sax <- filter.allele.count(GT_arc_sax, min_count)

### 5. Save the filtered data ####
write.table(GT_NS, paste0("Genotype_tables/", LG.s, "_genotypes_v2_NS.txt"), quote = F, row.names = F)
write.table(GT_Sp, paste0("Genotype_tables/", LG.s, "_genotypes_v2_Sp.txt"), quote = F, row.names = F)
write.table(GT_arc, paste0("Genotype_tables/", LG.s, "_genotypes_v2_arc.txt"), quote = F, row.names = F)
write.table(GT_comp, paste0("Genotype_tables/", LG.s, "_genotypes_v2_comp.txt"), quote = F, row.names = F)
write.table(GT_arc_sax, paste0("Genotype_tables/", LG.s, "_genotypes_v2_arc_sax.txt"), quote = F, row.names = F)
