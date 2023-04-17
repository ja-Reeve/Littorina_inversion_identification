### Calcaulte average observed heterozygosity per contig ###
### This script using L. saxatilis data from a given linkage group
### to calculate the average observed heterozygosity of each contig
### for each individual. 
### James Reeve - Göteborgs Universitet
### 14/12/2021
### Edited 22/02/2023: cleaned up code for publication

#### Preparation ####
rm(list = ls())
dev.off()
setwd("/Users/james/Documents/Inversion_detection")

# Select linkage group
LG.n <- 17 # Change to whichever LG you want to analyse
LG.s <- paste0("LG", LG.n)

# Packages
library("tidyverse")

#### 1. Upload and sort required data ####
### Collection details
sample_info <- read.csv("Seans_WGS_sample_info_v2.csv")

### Linkage map
LG <- read.table("ConsensusMap.v2.txt", header = T)
# Subset to current linkage group
LG <- LG[LG$LG == LG.s,]

### Genotypes
# Don't need to analyse the arc-sax or full GTs, since they will only be used for PCA
GT_NS <- read.table(paste0("Genotype_tables/",LG.s,"_genotypes_v2_NS.txt"), header = T, check.names = F)
GT_Sp <- read.table(paste0("Genotype_tables/",LG.s,"_genotypes_v2_Sp.txt"), header = T, check.names = F)
GT_arc <- read.table(paste0("Genotype_tables/",LG.s,"_genotypes_v2_arc.txt"), header = T, check.names = F)
GT_comp <- read.table(paste0("Genotype_tables/",LG.s,"_genotypes_v2_comp.txt"), header = T, check.names = F)

# Filter out contigs with < 10 SNPs
GT_NS <- GT_NS %>% group_by(CHROM) %>% filter(n() >= 10)
GT_Sp <- GT_Sp %>% group_by(CHROM) %>% filter(n() >= 10)
GT_arc <- GT_arc %>% group_by(CHROM) %>% filter(n() >= 10)
GT_comp <- GT_comp %>% group_by(CHROM) %>% filter(n() >= 10)

#### 2. Calculate the proportion of heterozygotes ####

### Calculate propHe for each individual for each contig
Heterozygotsiy.calculator <- function(data){
  tmp_GT <- apply(data[,-c(1:2)], 2, function(X){X[X != "./."]}, simplify = FALSE) # Simplify argument needed to avoid outputting a matrix for contigs with no missing sites
  tmp_Het <- lapply(tmp_GT, function(X){
    if(length(X) > 0){
      return(data.frame("nHo_ref" = sum(X == "0/0"), # Count of homozygote reference
                        "nHo_alt" = sum(X == "1/1"), # Count of homozygote alternative
                        "nHe" = sum(X == "0/1"), # Count of heterozygotes
                        "pHe" = sum(X == "0/1") / length(X) # Proportion of heterozygotes
      ))
    } else {
      # Set genotype frequencies to 0 if missing SNPs
      return(data.frame("nHo_ref" = 0.0,
                        "nHo_alt" = 0.0,
                        "nHe" = 0.0,
                        "pHe" = 0.0))
    }
  })
  Het <- do.call(rbind.data.frame, tmp_Het)
  res <- data.frame("SnailIDs" = rownames(Het), 
                    "CHROM" = rep(data$CHROM[1], nrow(Het)),
                    "Nsnp" = sapply(tmp_GT, length),
                    "nHo_ref" = Het$nHo_ref,
                    "nHo_alt" = Het$nHo_alt,
                    "nHe" = Het$nHe,
                    "pHe" = Het$pHe)
  return(res)
}

### Run heterozygosity call
# Northern saxatilis
Het_NS <- lapply(unique(GT_NS$CHROM), function(X){
  Heterozygotsiy.calculator(GT_NS[GT_NS$CHROM == X,])})
Het_NS <- do.call(rbind.data.frame, Het_NS)

# Spanish saxatilis
Het_Sp <- lapply(unique(GT_Sp$CHROM), function(X){
  Heterozygotsiy.calculator(GT_Sp[GT_Sp$CHROM == X,])})
Het_Sp <- do.call(rbind.data.frame, Het_Sp)
# Littorina arcana
Het_arc <- lapply(unique(GT_arc$CHROM), function(X){
  Heterozygotsiy.calculator(GT_arc[GT_arc$CHROM == X,])})
Het_arc <- do.call(rbind.data.frame, Het_arc)

# Littorina compressa
Het_comp <- lapply(unique(GT_comp$CHROM), function(X){
  Heterozygotsiy.calculator(GT_comp[GT_comp$CHROM == X,])})
Het_comp <- do.call(rbind.data.frame, Het_comp)


#### 3. Add linkage map positions ####
### Get contig position from linkage map
add.map.positions <- function(data){
  data %>% 
    group_by(CHROM) %>%
    # Merge with LG data (note: multiple positions per contig)
    left_join(LG, by = c("CHROM" = "contig")) %>%
    # Remove unnecessary columns
    select(-c(av, pos, deltaMP)) %>%
    # Remove duplicated rows
    distinct()
}

Het_NS <- add.map.positions(Het_NS)
Het_Sp <- add.map.positions(Het_Sp)
Het_arc <- add.map.positions(Het_arc)
Het_comp <- add.map.positions(Het_comp)

#### 4. Save the data ####
write_csv(Het_NS, paste0("Heterozygosity_scores/", LG.s, "_Het_v2_NS.csv"))
write_csv(Het_Sp, paste0("Heterozygosity_scores/", LG.s, "_Het_v2_Sp.csv"))
write_csv(Het_arc, paste0("Heterozygosity_scores/", LG.s, "_Het_v2_arc.csv"))
write_csv(Het_comp, paste0("Heterozygosity_scores/", LG.s, "_Het_v2_comp.csv"))
