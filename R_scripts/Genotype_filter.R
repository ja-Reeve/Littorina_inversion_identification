### Filtering genotypes ###
### This script filters the genotype files from a WGS dataset. There are two
### filters applied; a filter of monomorphic sites per region and a filter of
### allele counts per region.
### James Reeve - University of Gothenburg
### 12/03/2021

### 1. Preparation ####
### Clean-up environemnt
rm(list = ls())
dev.off()
setwd("")
options(stringsAsFactors = FALSE)

### Packages
library("tidyverse")

### Set linkage group
LG.n <- 1
# Change to whichever LG you want to analyse
LG.s <- paste0("LG", LG.n)

### Linkage map
LG <- read.table("ConsensusMap.v2.txt", header = T)
LG <- LG[LG$LG == LG.s,]
LG$contig <- sapply(strsplit(LG$Marker, "_"), '[[', 1)

### Upload sample information
sample_info <- read.csv("/PATH/Sample_information.csv", sep=";", check.names = FALSE)

### Genotypes
GT <- read.table(paste0("Genotype_tables/",LG.s,"_genotypes.v2.txt"), header = FALSE)
# Add header from VCF file
colnames(GT) <- names(read.delim("vcf_header.txt", check.names = FALSE))[c(1:2, 10:117)]
colnames(GT)[1] <- "CHROM" #rename 1st column
# Rename 1 odd Larc ID
names(GT)[names(GT) == "W_arc_04_Lamerged_sorted.bam"] <- "W_arc_04_La"

### 2. Create region IDs ####
### Groups are based on the phylogentic analysis of Sean Stankowski
### ES = East UK & Sweden (i.e. North Sea) - except for Dersingham, Norfolk
### NR = Norway & Russia
### FW = France & West UK (i.e. Celtic Sea) - including Dersingham, Norfolk
### Sp = Spain
### O = USA & Iceland (i.e. the others)
### Last two groups are sister species
### arc = Littorina arcana
### comp = Littorina compressa

### Extract sample IDs from sample info file
# Seperate species
# Note: "CEA_Larc_F_1" was mislabelled as L. arcana
sax_IDs <- as.character(sample_info[sample_info$Species %in% c("saxatilis", "saxatilis?", "saxailis?"), "Sample_ID"])
arc_IDs <- as.character(sample_info[sample_info$Species == "arcana" & 
                                      sample_info$Sample_ID != "CEA_Larc_F_1", "Sample_ID"])
comp_IDs <- as.character(sample_info[sample_info$Species == "compressa" | 
                                       sample_info$Sample_ID == "CEA_Larc_F_1", "Sample_ID"])

# Seperate saxatilis cohorts
sax_info <- sample_info[sample_info$Sample_ID %in% sax_IDs, ]

ES_IDs <- as.character(sax_info[sax_info$Country %in% c("England", "Sweden") &
                                  sax_info$location != "Dersingham" |
                                  sax_info$location == "S Abbs", "Sample_ID"]) 
NR_IDs <- as.character(sax_info[sax_info$Country %in% c("Norway", "Russia"), "Sample_ID"])
FW_IDs <- as.character(sax_info[sax_info$Country %in% c("Wales", "Ireland", "Isle of man", "France") |
                                  sax_info$location %in% c("Oban", "Dersingham"), "Sample_ID"])
Sp_IDs <- as.character(sax_info[sax_info$Country == "Spain", "Sample_ID"])
O_IDs <- as.character(sax_info[sax_info$Country %in% c("Iceland", "USA"), "Sample_ID"])

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
  # Boolian vector output (i.e. TRUE / FALSE)
  keeps <- apply(data, 1, sub.function)
  
  # Remove SNPs that fail the filter
  res <- data[keeps,]
  return(res)
}

### 4. Filter the genotype table ####
### Filter to contigs on LG
GT <- GT[GT$CHROM %in% unique(LG$contig),]

### Split genotypes by cohort
GT_ES <- GT %>% select("CHROM", "POS", all_of(ES_IDs))
GT_NR <- GT %>% select("CHROM", "POS", all_of(NR_IDs))
GT_FW <- GT %>% select("CHROM", "POS", all_of(FW_IDs))
GT_Sp <- GT %>% select("CHROM", "POS", all_of(Sp_IDs))
GT_O <- GT %>% select("CHROM", "POS", all_of(O_IDs))
GT_arc <- GT %>% select("CHROM", "POS", all_of(arc_IDs))
GT_comp <- GT %>% select("CHROM", "POS", all_of(comp_IDs))

### Filter by minor allele count
min_count <- 2 # Minor allele count threshold

GT_ES <- filter.allele.count(GT_ES, min_count)
GT_FW <- filter.allele.count(GT_FW, min_count)
GT_NR <- filter.allele.count(GT_NR, min_count)
GT_Sp <- filter.allele.count(GT_Sp, min_count)
GT_O <- filter.allele.count(GT_O, min_count)
GT_arc <- filter.allele.count(GT_arc, min_count)
GT_comp <- filter.allele.count(GT_comp, min_count)

### 5. Save the filtered data ####
write.table(GT_ES, paste0("Genotype_tables/", LG.s, "_genotypes_ES.txt"), quote = F, row.names = F)
write.table(GT_NR, paste0("Genotype_tables/", LG.s, "_genotypes_NR.txt"), quote = F, row.names = F)
write.table(GT_FW, paste0("Genotype_tables/", LG.s, "_genotypes_FW.txt"), quote = F, row.names = F)
write.table(GT_Sp, paste0("Genotype_tables/", LG.s, "_genotypes_Sp.txt"), quote = F, row.names = F)
write.table(GT_O, paste0("Genotype_tables/", LG.s, "_genotypes_O.txt"), quote = F, row.names = F)
write.table(GT_arc, paste0("Genotype_tables/", LG.s, "_genotypes_arc.txt"), quote = F, row.names = F)
write.table(GT_comp, paste0("Genotype_tables/", LG.s, "_genotypes_comp.txt"), quote = F, row.names = F)
