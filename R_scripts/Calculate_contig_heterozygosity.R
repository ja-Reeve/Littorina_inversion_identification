### Inversion detection - heterozygosity method ###
### This script using L. saxatilis data from a given linkage group
### to calcualte the average observed heterozygosity of each contig
### for each individual. This analysis is split between three subsets
### 1) Spainish L. saxatilis, 2) nothern L.saxatilis and 3) L. arcana.
### James Reeve - Göteborgs Universitet
### 14/12/2021

#### Preparation ####
rm(list = ls())
dev.off()
setwd("")

Start_time <- Sys.time()

# Select linkage group
LG.n <- 1 # Change to whichever LG you want to analyse
LG.s <- paste0("LG", LG.n)

# Packages
library("tidyverse")
library("ggpubr")

#### 1. Upload and sort required data ####
### Collection details
sample_info <- read.csv("/PATH/Sample_information.csv", 
                        sep=";", check.names = FALSE)
# Filter out missing samples
sample_info <- sample_info[sample_info$Sample_ID %in% colnames(GT)[-c(1:2)],]

### Linkage map
LG <- read.table("ConsensusMap.v2.txt", header = T)
# Subset to current linkage group
LG <- LG[LG$LG == LG.s,]

### Genotypes
# Read in as a function, to avoid loading a massive file into R
Genotype.read <- function(linkage.group){
  tmp <- read.table(paste0("Genotype_tables/",linkage.group,"_genotypes.v2.txt"), header = FALSE)
  colnames(tmp) <- names(read.delim("vcf_header.txt", check.names = FALSE))[c(1:2, 10:117)] # 10:86 for sax-only data
  colnames(tmp)[1] <- "CHROM"
  # Subset to linkage group
  tmp <- tmp[tmp$CHROM %in% LG$contig,]
  # Remove any contigs with < 10 SNPs
  res <- tmp %>% group_by(CHROM) %>% filter(n() >= 10)
  return(res)
}
# Access filtered genotypes
GT <- Genotype.read(LG.s)

#### 2. Subset genotypes to different regions ####
### Spain/North split is based on the phylogenetic analysis of Sean Stankowski
### Sp = Spanish Littorina saxatilis
### sax = Littorina saxatilis from all other sites (i.e. Northern cohort)
### Sister species
### arc = Littorina arcana
### comp = Littorina compressa

### Extract sample IDs from sample info file
# Separate species
# Note: "CEA_Larc_F_1" was mislabelled as L. arcana
sax_IDs <- as.character(sample_info[sample_info$Species %in% c("saxatilis", "saxatilis?", "saxailis?"), "Sample_ID"])
arc_IDs <- as.character(sample_info[sample_info$Species == "arcana" & 
                                      sample_info$Sample_ID != "CEA_Larc_F_1", "Sample_ID"])
comp_IDs <- as.character(sample_info[sample_info$Species == "compressa" | 
                                       sample_info$Sample_ID == "CEA_Larc_F_1", "Sample_ID"])

# Separate saxatilis cohorts into Spanish IDs
Sp_IDs <- as.character(sample_info[sample_info$Country == "Spain", "Sample_ID"])
# Remove Spanish samples from sax_IDs
sax_IDs <- sax_IDs[!(sax_IDs %in% Sp_IDs)]

### Filter the genotypes
# Remove arc & comp specific SNPs from saxatilis data
GT_sax <- GT %>% select("CHROM", "POS", all_of(sax_IDs))
GT_sax <- GT_sax[apply(GT_sax, 1, function(i){
  sum(i[-c(1:2)] %in% c("./.", "0/0")) | sum(i[-c(1:2)] %in% c("./.", "1/1")) == length(i) - 2}), ]

# Spanish saxatilis cohort
GT_Sp <- GT %>% select("CHROM", "POS", all_of(Sp_IDs))
GT_Sp <- GT_Sp[apply(GT_Sp, 1, function(i){
  sum(i[-c(1:2)] %in% c("./.", "0/0")) | sum(i[-c(1:2)] %in% c("./.", "1/1")) == length(i) - 2}), ]

# Remove sax & comp specific SNPs from arcana
colnames(GT)[colnames(GT) == "W_arc_04_Lamerged_sorted.bam"] <- "W_arc_04_La" # Rename sample with file details in name
GT_arc <- GT %>% select("CHROM", "POS", all_of(arc_IDs))
GT_arc <- GT_arc[apply(GT_arc, 1, function(i){
  sum(i[-c(1:2)] %in% c("./.", "0/0")) | sum(i[-c(1:2)] %in% c("./.", "1/1")) == length(i) - 2}), ]

# Remove sax & arc specific SNPs from compressa
GT_comp <- GT %>% select("CHROM", "POS", all_of(comp_IDs))
GT_comp <- GT_comp[apply(GT_comp, 1, function(i){
  sum(i[-c(1:2)] %in% c("./.", "0/0")) | sum(i[-c(1:2)] %in% c("./.", "1/1")) == length(i) - 2}), ]


#### 3. Calculate the proportion of heterozygotes ####

### Calculate propHe for each individual for each contig
Heterozygotsiy.calculator <- function(data){
  tmp_GT <- sapply(data[,-c(1:2)], function(X){X[X != "./."]})
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
Het_sax <- lapply(unique(GT_sax$CHROM), function(X){
  Heterozygotsiy.calculator(GT_sax[GT_sax$CHROM == X,])})
Het_sax <- do.call(rbind.data.frame, Het_sax)
Het_sax <- Het_sax[Het_sax$SnailIDs %in% sax_IDs,]

# Spanish saxatilis
Het_Sp <- lapply(unique(GT_Sp$CHROM), function(X){
  Heterozygotsiy.calculator(GT_Sp[GT_Sp$CHROM == X,])})
Het_Sp <- do.call(rbind.data.frame, Het_Sp)
Het_Sp <- Het_Sp[Het_Sp$SnailIDs %in% Sp_IDs,]

# Littorina arcana
Het_arc <- lapply(unique(GT_arc$CHROM), function(X){
  Heterozygotsiy.calculator(GT_arc[GT_arc$CHROM == X,])})
Het_arc <- do.call(rbind.data.frame, Het_arc)
Het_arc <- Het_arc[Het_arc$SnailIDs %in% arc_IDs,]

# Littorina compressa
Het_comp <- lapply(unique(GT_comp$CHROM), function(X){
  Heterozygotsiy.calculator(GT_comp[GT_comp$CHROM == X,])})
Het_comp <- do.call(rbind.data.frame, Het_comp)
Het_comp <- Het_comp[Het_comp$SnailIDs %in% comp_IDs,]


#### 4. Add linkage map positions ####
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

Het_Sp <- add.map.positions(Het_Sp)
Het_arc <- add.map.positions(Het_arc)
Het_comp <- add.map.positions(Het_comp)

#### 5. Save the data ####
write_csv(Het_sax, paste0("Heterozygosity_scores/", LG.s, "_Het_sax.csv"))
write_csv(Het_Sp, paste0("Heterozygosity_scores/", LG.s, "_Het_Sp.csv"))
write_csv(Het_arc, paste0("Heterozygosity_scores/", LG.s, "_Het_arc.csv"))
write_csv(Het_comp, paste0("Heterozygosity_scores/", LG.s, "_Het_comp.csv"))

#### 5. Plot results ####
### Compare Ho and He
H.plot1 <- function(data, title.text){
  ggplot(data)+
    geom_violin(aes("Homozygote ref.", nHo_ref/Nsnp), fill = "green")+
    geom_violin(aes("Heterozygote", nHe/Nsnp), fill = "blue")+
    geom_violin(aes("Homozygote alt.", nHo_alt/Nsnp), fill = "red")+
    ggtitle(title.text)+
    theme_bw()+
    theme(axis.title = element_blank())
}

p1A <- H.plot1(Het_sax, paste0("Northern L.saxatilis (N = ", length(sax_IDs), "; n = ", length(unique(Het_sax$CHROM)), ")"))
p1B <- H.plot1(Het_Sp, paste0("Spain (N = ", length(Sp_IDs), "; n = ", length(unique(Het_Sp$CHROM)), ")"))
p1C <- H.plot1(Het_arc, paste0("L.arcana (N = ", length(arc_IDs), "; n = ", length(unique(Het_arc$CHROM)), ")"))
p1D <- H.plot1(Het_comp, paste0("L.compressa (N = ", length(comp_IDs), "; n = ", length(unique(Het_comp$CHROM)), ")"))

tiff(paste0("Heterozygosity_scores/", LG.s, "_HetViloinPlot.tiff"), width = 620)
annotate_figure(ggarrange(p1A, p1B, p1C, p1D), 
                bottom = "Genotype", left = "Proportions of SNPs per contig")
dev.off()

### Plot heterozygosity by linkage map position
H.plot2 <- function(data, title.text){
  ggplot(data, aes(avgMP, pHe))+
    geom_point(size = 0.5, alpha = 0.4)+
    geom_smooth(aes(group = SnailIDs), 
                se = F, span = 0.2, alpha = 0.5, size = 0.5)+
    ggtitle(title.text)+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

p2A <- H.plot2(Het_sax, paste0(LG.s, ": Nothern saxatilis")) + coord_cartesian(ylim=c(0, 0.5))
p2B <- H.plot2(Het_Sp, paste0(LG.s, ": Spanish saxatilis")) + coord_cartesian(ylim=c(0, 0.5))
p2C <- H.plot2(Het_arc, paste0(LG.s, ": Littorina arcana")) + coord_cartesian(ylim=c(0, 0.5))
p2D <- H.plot2(Het_comp, paste0(LG.s, ": Littorina compressa")) + coord_cartesian(ylim=c(0, 0.5))
# coord_cartesian() zooms into a section of the plot without removing data. This prevents geom_smooth()
# from changing it's projection.

annotate_figure(ggarrange(p2A, p2B, p2C, p2D, nrow = 4, ncol = 1), 
             left = "Proportion of heterozygotes per contig",
             bottom = "Position on linkage map (cM)")

End_time <- Sys.time()
End_time - Start_time
