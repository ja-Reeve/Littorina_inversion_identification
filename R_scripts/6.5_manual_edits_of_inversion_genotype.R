################### Manual adjustments to inversion genotyping #######################
### Some automatic genotyping of PCA clusters did not work well with the Kmeans 
### approach. This script makes manual adjustments to these Kmeans exceptions. Each
### of the following adjustements will be explained before then adjusting the groups
### and saving over the file.
### James Reeve - University of Gothenburg
### 10/02/2023

### Preparation
rm(list = ls())
dev.off()
setwd("/Users/james/Documents/Inversion_detection")

### Packages
library("tidyverse")
library("ggpubr")
library("cluster")
library("adegenet")

#### 1: LGC2.1 - Northern Littorina saxatilis ####
### There is a separate cluster on the top left of the PCA plot that was included
### in the 'RA'. However, visually it is close to the 'AA' cluster.
LGC2.1 <- read.csv("PCA_per_inversion/LG2_PCA_of_LGC2.1_v2_NS.csv")

ggplot(LGC2.1)+
  geom_point(aes(Axis1, Axis2, col = genotype))+
  theme_classic()

# Check if this other group is geographically isolated
LGC2.1[LGC2.1$genotype == "RA" & LGC2.1$Axis2 < -30,]
# Yes! All samples come from USA

USA_samples <- LGC2.1[LGC2.1$Country == "USA", "Sample_ID"]

### Re-run PCA projecting USA onto Northern sax to see where these samples fall

LG <- "LG2"

# Inversion positions
Lmap <- read.table("ConsensusMap.v2.txt", header = T)
Lmap <- Lmap[Lmap$LG == LG & Lmap$avgMP > 0 & Lmap$avgMP < 14,]


### Access genotypes
GT <- read.table(paste0("Genotype_tables/LG2_genotypes_v2_NS.txt"), 
                 header = TRUE, check.names = FALSE)

# Subset to inversion
GT <- GT[GT$CHROM %in% Lmap$contig,]

## Wrangle GT for PCA
GT2 <- GT[apply(GT, 1, function(i){
  sum(i[-c(1:2)] %in% c("./.", "0/0")) != length(i) - 2 | 
    sum(i[-c(1:2)] %in% c("./.", "1/1")) != length(i) - 2 |
    sum(i[-c(1:2)] %in% c("./.", "0/1")) != length(i) - 2}), ]

GT2 <- GT2 %>% unite(SNP_ID, c("CHROM", "POS"), sep=":")
# Transpose using function from data.table package
GT3 <- t(GT2)
colnames(GT3) <- GT2$SNP_ID
rownames(GT3) <- colnames(GT2)
GT3 <- GT3[-c(1),]

# Make genind
genind <- df2genind(GT3, ploidy=2, sep="/", NA.char = "_")
# B: rescore missing genotypes as the mean values
genind <- scaleGen(genind, NA.method="mean", scale = FALSE)


## Run PCA
genind_usa <- genind[rownames(genind) %in% USA_samples,]
#C.2: run PCA on genind without compressa
GT_PCA <- dudi.pca(genind[!(rownames(genind) %in% USA_samples),], 
                   cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)	
#C.3: project compressa genotypes onto the PCA
usa_proj <- suprow(GT_PCA, genind_usa)
# D: save principle component scores & % of variance explained
GT_PCs <- rbind.data.frame(GT_PCA$li, usa_proj[[2]])
GT_PCs$Sample_ID <- rownames(GT_PCs)

### Plot samples with new PCA
tmp <- inner_join(GT_PCs, LGC2.1[-c(1:3)], by = c("Sample_ID" = "Sample_ID"))

ggplot(tmp, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = Country == "USA"))+
  theme_bw()

tmp[tmp$Country == "USA",]


#### Summary: USA samples are split between two arrangements. "York_B-12_Ls", "York_T-3_Ls" & "York_T-9_Ls" are "AA". While "York_B-1_LS" is "RA".

### Fix arragements and save over original file
LGC2.1[LGC2.1$Sample_ID %in% c("York_B-12_Ls", "York_T-3_Ls", "York_T-9_Ls"), "genotype"] <- "AA"
write.csv(LGC2.1, "PCA_per_inversion/LG2_PCA_of_LGC2.1_v2_NS.csv", quote = FALSE, row.names = FALSE)

#### END
rm(list = ls())
dev.off()


#### 2: LGC9.2 - Northern Littorina saxatilis ####
### 'AA' and 'RA' blend into one another on the plot. Let's take a closer look to see if we can refine the categories
LGC9.2 <- read.csv("PCA_per_inversion/LG9_PCA_of_LGC9.2_v2_NS.csv")

ggplot(LGC9.2, aes(Axis1, Axis2))+
  geom_point(aes(col = Country))+
  theme_bw()

ggplot(LGC9.2, aes(Axis1, Axis2))+
  geom_point(aes(col = Region == "North_Sea"))+
  theme_bw()

### PCA is biased by outliers form North Sea. Probably because this inversion in missing in these regions.
### Try re-running PCA without North Sea.

NS_samples <- LGC9.2[LGC9.2$Region == "North_Sea", "Sample_ID"]

LG <- "LG9"

# Inversion positions
Lmap <- read.table("ConsensusMap.v2.txt", header = T)
Lmap <- Lmap[Lmap$LG == LG & Lmap$avgMP > 52 & Lmap$avgMP < 59.47085,]


### Access genotypes
GT <- read.table(paste0("Genotype_tables/LG9_genotypes_v2_NS.txt"), 
                 header = TRUE, check.names = FALSE)

# Subset to inversion
GT <- GT[GT$CHROM %in% Lmap$contig,]

## Wrangle GT for PCA
GT2 <- GT[apply(GT, 1, function(i){
  sum(i[-c(1:2)] %in% c("./.", "0/0")) != length(i) - 2 | 
    sum(i[-c(1:2)] %in% c("./.", "1/1")) != length(i) - 2 |
    sum(i[-c(1:2)] %in% c("./.", "0/1")) != length(i) - 2}), ]

GT2 <- GT2 %>% unite(SNP_ID, c("CHROM", "POS"), sep=":")
# Transpose using function from data.table package
GT3 <- t(GT2)
colnames(GT3) <- GT2$SNP_ID
rownames(GT3) <- colnames(GT2)
GT3 <- GT3[-c(1),]

# Make genind
genind <- df2genind(GT3, ploidy=2, sep="/", NA.char = "_")
# B: rescore missing genotypes as the mean values
genind <- scaleGen(genind, NA.method="mean", scale = FALSE)


## Run PCA
#C.2: run PCA on genind without North Sea
GT_PCA <- dudi.pca(genind[!(rownames(genind) %in% NS_samples),], 
                   cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)	

# D: save principle component scores
GT_PCs <- GT_PCA$li

### Plot samples with new PCA
GT_PCs$Sample_ID <- rownames(GT_PCs)
tmp <- inner_join(GT_PCs, LGC9.2[!(LGC9.2$Sample_ID %in% NS_samples), -c(1:3)], by = c("Sample_ID" = "Sample_ID"))

ggplot(tmp, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = Region))+
  theme_bw()

### Re-run PCA projecting North Sea onto Northern sax to see where these samples fall
genind_NS <- genind[rownames(genind) %in% NS_samples,]
NS_proj <- suprow(GT_PCA, genind_NS)

GT_PCs.proj <- rbind.data.frame(GT_PCA$li, NS_proj[[2]])
GT_PCs.proj$Sample_ID <- rownames(GT_PCs.proj)

ggplot(GT_PCs.proj, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = LGC9.2$Region))+
  theme_bw()

### Odd, this doesn't help things. Either the signal is very weak, or the North Sea 
### populations are highly diverged (this is unlikely, since other invs don't show this)
### Roger suggested this might be beacuse this is a very old inversion, which has lost it
### single, or fixed, in some locations while remaining polymorphic in others.

# Plot LGC9.2 without North Sea samples
tmp2 <- GT_PCA$li
tmp2$Sample_ID <- rownames(tmp2)

tmp2 <- inner_join(tmp2, LGC9.2[-c(1:3)], by = c("Sample_ID" = "Sample_ID"))

ggplot(tmp2, aes(Axis1, Axis2))+
  geom_point(aes(col = genotype))+
  theme_bw()

ggplot(tmp2, aes(Axis1, Axis2))+
  geom_point(aes(col = location))+
  theme_bw()

### Looks like inversion is only present in Holyhead and Roscoff!
### Relabeling might be possible, but the US samples keep jumping around
### either clustering with the far right or far left groups. Given this
### uncertainty and the regional-specificity of this inversion in Northern 
### saxatilis I will not relabel it.

#### End

rm(list = ls())
dev.off()

#### 3: LGC1.1 - Northern Littorina saxatilis ####
### Due to some variation on PC2, K-means doesn't do a good job at finding the 3 cluster that are visually clear. I am manually adding these categories.
LGC1.1 <- read.csv("PCA_per_inversion/LG1_PCA_of_LGC1.1_v2_NS.csv")

ggplot(LGC1.1, aes(Axis1, Axis2))+
  geom_point(aes(col = Ecotype))+
  theme_bw()

### Crab and Wave are evenly spread among clusters, let's call the right hand cluster as 'RR' as presumably the arrangement with the majority of samples 
### should be the reference
LGC1.1$genotype <- 'RR'
LGC1.1[LGC1.1$Axis1 > -40 & LGC1.1$Axis1 < -10, "genotype"] <- "RA"
LGC1.1[LGC1.1$Axis1 < -40, "genotype"] <- "AA"


ggplot(LGC1.1, aes(Axis1, Axis2))+
  geom_point(aes(col = genotype))+
  theme_bw()

### Save over original file
write.csv(LGC1.1, "PCA_per_inversion/LG1_PCA_of_LGC1.1_v2_NS.csv", quote = FALSE, row.names = FALSE)

#### End

rm(list = ls())
dev.off()



#### 4: LGC5.1 - arcana vs. saxatilis ####
### It looks like the inversion arrangements separate by PC2, while PC1 separates the species.

LGC5.1 <- read.csv("PCA_per_inversion/LG5_PCA_of_LGC5.1_v2_arc-sax.csv")

ggplot(LGC5.1, aes(Axis1, Axis2))+
  geom_point(aes(col = genotype, pch = Species))+
  theme_bw()

### Re-run PCA projecting arcana onto the saxatilis axes

LG <- "LG5"

# Inversion positions
Lmap <- read.table("ConsensusMap.v2.txt", header = T)
Lmap <- Lmap[Lmap$LG == LG & Lmap$avgMP >= 15 & Lmap$avgMP <= 47,]


### Access genotypes
GT <- read.table(paste0("Genotype_tables/LG5_genotypes_v2_arc_sax.txt"), header = TRUE, check.names = FALSE)

# Subset to inversion
GT <- GT[GT$CHROM %in% Lmap$contig,]

## Wrangle GT for PCA
GT2 <- GT[apply(GT, 1, function(i){
  sum(i[-c(1:2)] %in% c("./.", "0/0")) != length(i) - 2 | 
    sum(i[-c(1:2)] %in% c("./.", "1/1")) != length(i) - 2 |
    sum(i[-c(1:2)] %in% c("./.", "0/1")) != length(i) - 2}), ]

GT2 <- GT2 %>% unite(SNP_ID, c("CHROM", "POS"), sep=":")
# Transpose using function from data.table package
GT3 <- t(GT2)
colnames(GT3) <- GT2$SNP_ID
rownames(GT3) <- colnames(GT2)
GT3 <- GT3[-c(1),]

# Make genind
genind <- df2genind(GT3, ploidy=2, sep="/", NA.char = "_")
# B: rescore missing genotypes as the mean values
genind <- scaleGen(genind, NA.method="mean", scale = FALSE)

## List saxatilis and arcana samples separately
sax_IDs <- LGC5.1[LGC5.1$Species == "saxatilis", "Sample_ID"]
arc_IDs <-  LGC5.1[LGC5.1$Species == "arcana", "Sample_ID"]

## Run PCA projecting saxatilis onto the arcana axes
genind_sax <- genind[rownames(genind) %in% sax_IDs,]
# Run PCA on genind for just arcana
sax_arc_PCA <- dudi.pca(genind[rownames(genind) %in% arc_IDs,], 
                   cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)	
# Project arcana genotypes onto the PCA
sax_proj <- suprow(sax_arc_PCA, genind_sax)
# Save principle component scores & % of variance explained
sax_arc_PCs <- rbind.data.frame(sax_arc_PCA$li, sax_proj[[2]])


## Run PCA projecting arcana onto the saxatilis axes
genind_arc <- genind[rownames(genind) %in% arc_IDs,]
# Run PCA on genind for just saxatilis
arc_sax_PCA <- dudi.pca(genind[rownames(genind) %in% sax_IDs,], 
                        cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)	
# Project arcana genotypes onto the PCA
arc_proj <- suprow(arc_sax_PCA, genind_arc)
# Save principle component scores & % of variance explained
arc_sax_PCs <- rbind.data.frame(arc_sax_PCA$li, arc_proj[[2]])


### Plot samples with new PCA
sax_arc_PCs$Sample_ID <- rownames(sax_arc_PCs)
tmp <- inner_join(sax_arc_PCs, LGC5.1[,-c(1:3)], by = "Sample_ID")

p1 <- ggplot(tmp, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = Species))+
  labs(title = "saxatilis projected to arcana")+
  theme_bw()


arc_sax_PCs$Sample_ID <- rownames(arc_sax_PCs)
tmp <- inner_join(arc_sax_PCs, LGC5.1[,-c(1:3)], by = "Sample_ID")

p2 <- ggplot(tmp, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = Species))+
  labs(title = "arcana projected to saxatilis")+
  theme_bw()

ggarrange(p1, p2)

### 4 clusters are found whether we project arc->sax or sax->arc. The sax->arc projection looks
### like the inversion is only present in arcana, while the arc->sax projection doesn't look like 
### an inversion at all. It's safest to leave this inversion out!


#### End

rm(list = ls())
dev.off()



#### 5: LGC4.1 - arcana vs. saxatilis ####
### This inversion was left out due to unclear clustering on arcana in the arc-sax PCA. 
### Try projecting arcana onto the saxatilis axis.

LGC4.1 <- read.csv("PCA_per_inversion/LG4_PCA_of_LGC4.1_v2_arc-sax.csv")

ggplot(LGC4.1, aes(Axis1, Axis2, col = Species))+
  geom_point()+
  geom_label(aes(label = Sample_ID), size = 2)+
  theme_bw()

### Re-run PCA projecting arcana onto the saxatilis axes

LG <- "LG4"

# Inversion positions
Lmap <- read.table("ConsensusMap.v2.txt", header = T)
Lmap <- Lmap[Lmap$LG == LG & Lmap$avgMP >= 0 & Lmap$avgMP <= 24,]


### Access genotypes
GT <- read.table(paste0("Genotype_tables/LG4_genotypes_v2_arc_sax.txt"), 
                 header = TRUE, check.names = FALSE)

# Subset to inversion
GT <- GT[GT$CHROM %in% Lmap$contig,]

## Wrangle GT for PCA
GT2 <- GT[apply(GT, 1, function(i){
  sum(i[-c(1:2)] %in% c("./.", "0/0")) != length(i) - 2 | 
    sum(i[-c(1:2)] %in% c("./.", "1/1")) != length(i) - 2 |
    sum(i[-c(1:2)] %in% c("./.", "0/1")) != length(i) - 2}), ]

GT2 <- GT2 %>% unite(SNP_ID, c("CHROM", "POS"), sep=":")
# Transpose using function from data.table package
GT3 <- t(GT2)
colnames(GT3) <- GT2$SNP_ID
rownames(GT3) <- colnames(GT2)
GT3 <- GT3[-c(1),]

# Make genind
genind <- df2genind(GT3, ploidy=2, sep="/", NA.char = "_")
# B: rescore missing genotypes as the mean values
genind <- scaleGen(genind, NA.method="mean", scale = FALSE)

## List saxatilis and arcana samples separately
sax_IDs <- LGC4.1[LGC4.1$Species == "saxatilis", "Sample_ID"]
arc_IDs <-  LGC4.1[LGC4.1$Species == "arcana", "Sample_ID"]

## Run PCA
genind_arc <- genind[rownames(genind) %in% arc_IDs,]
# Run PCA on genind for just saxatilis
GT_PCA <- dudi.pca(genind[rownames(genind) %in% sax_IDs,], 
                   cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)	
# Project arcana genotypes onto the PCA
arc_proj <- suprow(GT_PCA, genind_arc)
# Save principle component scores & % of variance explained
GT_PCs <- rbind.data.frame(GT_PCA$li, arc_proj[[2]])


### Plot samples with new PCA
GT_PCs$Sample_ID <- rownames(GT_PCs)
tmp <- inner_join(GT_PCs, LGC4.1[,-c(1:3)], by = "Sample_ID")

ggplot(tmp, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = Species))+
  theme_bw()

ggplot(tmp, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = Heterozygosity, pch = Species))+
  theme_bw()

ggplot(tmp, aes(Axis1, Axis2, col = Species))+
  geom_point()+
  geom_label(aes(label = Sample_ID), size = 2)+
  theme_bw()

### This was hard to call since projection moved all the arcana samples to the right
### side of the plot, but there were 5 samples which separated away from the other arcana.
### The average heterozygosity of these samples was higher than other arcana, so they were
### defined as heterozygotes.

### Rules for assigning genotypes (based on unporjected PCA):
###     saxatilis 'RR' = PC1 < -25
###     saxatilis 'RA' = -25 < PC1 <0
###     saxatilis 'AA' = PC1 > 0
###     arcana 'RR' doesn't exist
###     arcana 'RA' = 15 < PC1 < 30 & -5 < PC2 < 45
###     arcana 'AA' = all other samples

LGC4.1$genotype <- 'AA'
LGC4.1[LGC4.1$Species == "saxatilis" & LGC4.1$Axis1 < -25, "genotype"] <- 'RR'
LGC4.1[LGC4.1$Species == "saxatilis" & LGC4.1$Axis1 > -25 & LGC4.1$Axis1 < 0, "genotype"] <- 'RA'
LGC4.1[LGC4.1$Species == "arcana" & LGC4.1$Axis1 > 15 & LGC4.1$Axis1 < 30 &
         LGC4.1$Axis2 > -5 & LGC4.1$Axis2 < 45, "genotype"] <- "RA"

ggplot(LGC4.1)+
  geom_point(aes(x = Axis1, y = Axis2, col = genotype))+
  theme_bw()

### Save over original file
write.csv(LGC4.1, "PCA_per_inversion/LG4_PCA_of_LGC4.1_v2_arc-sax.csv", quote = FALSE, row.names = FALSE)

#### End

rm(list = ls())
dev.off()



#### 6: LGC14.1 - arcana vs. saxatilis ####
### The PCA for LGC14.1 seems slightly skewed by PC2. Probably because PC1 is explaining some of 
### the variation among species. Note, that from the intraspecific PCA, it looks like this 
### arrangement is fixed / absent in arcana.

LGC14.1 <- read.csv("PCA_per_inversion/LG14_PCA_of_LGC14.1_v2_arc-sax.csv")

ggplot(LGC14.1, aes(Axis1, Axis2))+
  geom_point(aes(col = Species))


### Try projecting Larc onto Lsax axes to see which cluster Larc belongs to.

LG <- "LG14"

# Inversion positions
Lmap <- read.table("ConsensusMap.v2.txt", header = T)
Lmap <- Lmap[Lmap$LG == LG & Lmap$avgMP >= 0.679 & Lmap$avgMP <= 8.8095,]

### Access genotypes
GT <- read.table(paste0("Genotype_tables/LG14_genotypes_v2_arc_sax.txt"), 
                 header = TRUE, check.names = FALSE)

GT <- GT[GT$CHROM %in% Lmap$contig,]

## Wrangle GT for PCA
GT2 <- GT[apply(GT, 1, function(i){
  sum(i[-c(1:2)] %in% c("./.", "0/0")) != length(i) - 2 | 
    sum(i[-c(1:2)] %in% c("./.", "1/1")) != length(i) - 2 |
    sum(i[-c(1:2)] %in% c("./.", "0/1")) != length(i) - 2}), ]

GT2 <- GT2 %>% unite(SNP_ID, c("CHROM", "POS"), sep=":")
# Transpose using function from data.table package
GT3 <- t(GT2)
colnames(GT3) <- GT2$SNP_ID
rownames(GT3) <- colnames(GT2)
GT3 <- GT3[-c(1),]


# Make genind
genind <- df2genind(GT3, ploidy=2, sep="/", NA.char = "_")
# B: rescore missing genotypes as the mean values
genind <- scaleGen(genind, NA.method="mean", scale = FALSE)

## List saxatilis and arcana samples separately
sax_IDs <- LGC14.1[LGC14.1$Species == "saxatilis", "Sample_ID"]
arc_IDs <-  LGC14.1[LGC14.1$Species == "arcana", "Sample_ID"]

## Run PCA
genind_arc <- genind[rownames(genind) %in% arc_IDs,]
# Run PCA on genind for just saxatilis
GT_PCA <- dudi.pca(genind[rownames(genind) %in% sax_IDs,], 
                   cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)	
# Project arcana genotypes onto the PCA
arc_proj <- suprow(GT_PCA, genind_arc)
# Save principle component scores & % of variance explained
GT_PCs <- rbind.data.frame(GT_PCA$li, arc_proj[[2]])


### Plot samples with new PCA
GT_PCs$Sample_ID <- rownames(GT_PCs)
tmp <- inner_join(GT_PCs, LGC14.1[,-c(1:3)], by = "Sample_ID")

ggplot(tmp, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = genotype, pch = Species))+
  theme_bw()

ggplot(tmp, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = Region, pch = Species))+
  theme_bw()

### Projecting has reorientated the saxatilis samples to make clearer clusters. They seem to separate
### clearly on PC1 with arcana aligning with the right cluster.

# Assign clusters with following rules:
LGC14.1$genotype[LGC14.1$Axis1 < -25 & LGC14.1$Axis2 < 0] <- "RR"
LGC14.1$genotype[LGC14.1$Axis1 > -30 & LGC14.1$Axis1 < -5 & LGC14.1$Axis2 > 0] <- "RA"
LGC14.1$genotype[LGC14.1$Axis1 > -5] <- "AA"

ggplot(LGC14.1, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = genotype, pch = Species))+
  theme_bw()

### Save over original file
write.csv(LGC14.1, "PCA_per_inversion/LG14_PCA_of_LGC14.1_v2_arc-sax.csv", quote = FALSE, row.names = FALSE)

#### End

rm(list = ls())
dev.off()



#### 7: LGC14.2 - arcana vs. saxatilis ####
### The PCA for LGC14.2 has a 6 cluster pattern for saxatilis, however arcana forms a separate 
### group with a lot of variation to the right hand side of the plot. Let's see if projecting 
### arcana will help.
LGC14.2 <- read.csv("PCA_per_inversion/LG14_PCA_of_LGC14.2_v2_arc-sax.csv")

ggplot(LGC14.2, aes(Axis1, Axis2))+
  geom_point(aes(col = Species))


### Try projecting Larc onto Lsax axes to see which cluster Larc belongs to.

LG <- "LG14"

# Inversion positions
Lmap <- read.table("ConsensusMap.v2.txt", header = T)
Lmap <- Lmap[Lmap$LG == LG & Lmap$avgMP >= 9.0965 & Lmap$avgMP <= 11.413,]


### Access genotypes
GT <- read.table(paste0("Genotype_tables/LG14_genotypes_v2_arc_sax.txt"), header = TRUE, check.names = FALSE)

# Subset to inversion
GT <- GT[GT$CHROM %in% Lmap$contig,]

## Wrangle GT for PCA
GT2 <- GT[apply(GT, 1, function(i){
  sum(i[-c(1:2)] %in% c("./.", "0/0")) != length(i) - 2 | 
    sum(i[-c(1:2)] %in% c("./.", "1/1")) != length(i) - 2 |
    sum(i[-c(1:2)] %in% c("./.", "0/1")) != length(i) - 2}), ]

GT2 <- GT2 %>% unite(SNP_ID, c("CHROM", "POS"), sep=":")
# Transpose using function from data.table package
GT3 <- t(GT2)
colnames(GT3) <- GT2$SNP_ID
rownames(GT3) <- colnames(GT2)
GT3 <- GT3[-c(1),]

# Make genind
genind <- df2genind(GT3, ploidy=2, sep="/", NA.char = "_")
# B: rescore missing genotypes as the mean values
genind <- scaleGen(genind, NA.method="mean", scale = FALSE)

## List saxatilis and arcana samples separately
sax_IDs <- LGC14.2[LGC14.2$Species == "saxatilis", "Sample_ID"]
arc_IDs <-  LGC14.2[LGC14.2$Species == "arcana", "Sample_ID"]

## Run PCA
genind_arc <- genind[rownames(genind) %in% arc_IDs,]
# Run PCA on genind for just saxatilis
GT_PCA <- dudi.pca(genind[rownames(genind) %in% sax_IDs,], 
                   cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)	
# Project arcana genotypes onto the PCA
arc_proj <- suprow(GT_PCA, genind_arc)
# Save principle component scores & % of variance explained
GT_PCs <- rbind.data.frame(GT_PCA$li, arc_proj[[2]])


### Plot samples with new PCA
GT_PCs$Sample_ID <- rownames(GT_PCs)
tmp <- inner_join(GT_PCs, LGC14.2[,-c(1:3)], by = "Sample_ID")

ggplot(tmp, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = Species))+
  theme_bw()

p1 <- ggplot(LGC14.2, aes(Axis1, Axis2))+
  geom_point(aes(col = Species, pch = Region))+
  labs(title = "Before projecting")+
  theme_bw()

p2 <- ggplot(tmp, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = Species))+
  labs(title = "After projecting")+
  theme_bw()

ggarrange(p1, p2, nrow = 1, common.legend = TRUE)

### Re-run K-means clustering
# Run Kmeans from K = 2:9
tmp2 <- lapply(2:9, function(K){
  # Kmeans clustering
  clust <- kmeans(tmp[,c("Axis1", "Axis2")], centers = K, iter.max = 1e5, nstart = 100)
  # Calcualte average silouette score (i.e. difference between clusters)
  sil <- mean(silhouette(clust$cluster, dist(tmp[,c("Axis1", "Axis2")]))[,"sil_width"])
  # Store clusters and silhouette values in output
  res <- append(clust, sil)
  names(res)[10] <- "Avg_Silhouette_Weight"
  return(res)})

# Find the best K (highest Silouette score)
K <- which.max(sapply(tmp2, `[[`, "Avg_Silhouette_Weight")) + 1

# D: Reassign cluster names based on centre values
bestK <- tmp2[[K - 1]]

### Find extreme vertices (homokaryotypes)
clust_A <- which.min(rowSums(bestK$centers)) # Searches for overall lowest centre value (PC1 & PC2)
clust_F <- which.min(bestK$centers[,1] * bestK$centers[,2]) # Search for bottom right cluster
clust_C <- which.max(rowSums(bestK$centers)) # Search PC2 for maximum centre value

### Find intermediate vertices (heterokaryotypes)
het_clusts <- bestK$centers[-c(clust_A, clust_C, clust_F), ] # Exclude the extreme clusters
clust_D <- names(which.min(het_clusts[,2])) # Find minimum score of PC2
clust_E <- names(which.max(rowSums(het_clusts))) # Search for overall maximum centre value of intermediate clusters
# Define last cluster by excluding the others
clust_B <- which(!(1:6 %in% c(clust_A, clust_C, clust_D, clust_E, clust_F)))

### Add clusters to tmp
tmp$genotype <- bestK$cluster

### Subset to just homokaryotypics (extreme) clusters
PC_homo <- tmp[tmp$genotype %in% c(clust_A, clust_C, clust_F), ]

### Find cluster with the most crab individuals in PC_homo
crab_clust <- PC_homo[PC_homo$Species == "saxatilis" & 
                        PC_homo$Ecotype == "Crab", "genotype"]
uClust <- unique(crab_clust)

# Find frequency of crab individuals in each homokaryotype cluster
pClust <- sapply(as.numeric(uClust), function(i){sum(crab_clust == i) / bestK$size[i]})

# Assign cluster with highest crab frequency as 'RR'
RR <- uClust[which.max(pClust)]

### Assign other clusters based on identity of 'RR'
if(RR == clust_A){
  A1A1 <- clust_F; A2A2 <- clust_C # Homokaryotype clusters
  RA1 <- clust_D; RA2 <- clust_B; A1A2 <- clust_E} # Heterokaryotype clusters
if(RR == clust_F){
  A1A1 <- clust_A; A2A2 <- clust_C
  RA1 <- clust_D; RA2 <- clust_E; A1A2 <- clust_B} 
if(RR == clust_C){
  A1A1 <- clust_F; A2A2 <- clust_A 
  RA1 <- clust_E; RA2 <- clust_B; A1A2 <- clust_D}

tmp$genotype[tmp$genotype == RR] <- "RR"
tmp$genotype[tmp$genotype == A1A1] <- "A1A1"
tmp$genotype[tmp$genotype == A2A2] <- "A2A2"
tmp$genotype[tmp$genotype == RA1] <- "RA1"
tmp$genotype[tmp$genotype == RA2] <- "RA2"
tmp$genotype[tmp$genotype == A1A2] <- "A1A2"

ggplot(tmp, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = factor(genotype)))+
  theme_bw()


### Clustering is still uncertain. It looks like LGC14.1/2 is missing from arcana! This is particularly
### odd since the inversion was visable on the arcana split plot for LG14. I won't save the results for
### this inversion as I am not convinced that arcana really is split between inversion arrangements.

#### End

rm(list = ls())
dev.off()


#### 7: LGC6.1/2 - arcana vs. saxatilis ####
### Including arcana into the PCA for LGC6.1/2 distorted the clustering, probably because the arrangement
### is not as clear in arcana. Let's see if projecting arcana onto saxatilis axes will help.

LGC6.1_2 <- read.csv("PCA_per_inversion/LG6_PCA_of_LGC6.1-2_v2_arc-sax.csv")

ggplot(LGC6.1_2, aes(Axis1, Axis2))+
  geom_point(aes(col = Species))+
  theme_bw()


### Try projecting Larc onto Lsax axes to see which cluster Larc belongs to.

LG <- "LG6"

# Inversion positions
Lmap <- read.table("ConsensusMap.v2.txt", header = T)
Lmap <- Lmap[Lmap$LG == LG & Lmap$avgMP <= 19,]


### Access genotypes
GT <- read.table(paste0("Genotype_tables/LG6_genotypes_v2_arc_sax.txt"), header = TRUE, check.names = FALSE)

# Subset to inversion
GT <- GT[GT$CHROM %in% Lmap$contig,]

## Wrangle GT for PCA
GT2 <- GT[apply(GT, 1, function(i){
  sum(i[-c(1:2)] %in% c("./.", "0/0")) != length(i) - 2 | 
    sum(i[-c(1:2)] %in% c("./.", "1/1")) != length(i) - 2 |
    sum(i[-c(1:2)] %in% c("./.", "0/1")) != length(i) - 2}), ]

GT2 <- GT2 %>% unite(SNP_ID, c("CHROM", "POS"), sep=":")
# Transpose using function from data.table package
GT3 <- t(GT2)
colnames(GT3) <- GT2$SNP_ID
rownames(GT3) <- colnames(GT2)
GT3 <- GT3[-c(1),]

# Make genind
genind <- df2genind(GT3, ploidy=2, sep="/", NA.char = "_")
# B: rescore missing genotypes as the mean values
genind <- scaleGen(genind, NA.method="mean", scale = FALSE)

## List saxatilis and arcana samples separately
sax_IDs <- LGC6.1_2[LGC6.1_2$Species == "saxatilis", "Sample_ID"]
arc_IDs <-  LGC6.1_2[LGC6.1_2$Species == "arcana", "Sample_ID"]

## Run PCA
genind_arc <- genind[rownames(genind) %in% arc_IDs,]
# Run PCA on genind for just saxatilis
GT_PCA <- dudi.pca(genind[rownames(genind) %in% sax_IDs,], 
                   cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)	
# Project arcana genotypes onto the PCA
arc_proj <- suprow(GT_PCA, genind_arc)
# Save principle component scores & % of variance explained
GT_PCs <- rbind.data.frame(GT_PCA$li, arc_proj[[2]])


### Plot samples with new PCA
GT_PCs$Sample_ID <- rownames(GT_PCs)
tmp <- inner_join(GT_PCs, LGC6.1_2[,-c(1:3)], by = "Sample_ID")

ggplot(tmp, aes(x = Axis1, y = Axis2))+
  geom_point(aes(col = Species))+
  theme_bw()


### After projecting L. arcana onto L. saxatilis axes, it seems clear that both species
### are polymorphic for this inversion, however some heterokaryotypic clusters are missing.
### Inferring the clusters from the Northern saxatilis PCA, the genotypes were assigned as follows:
###     'RR' = PC1 < -30
###     'RA1' = not present
###     'RA2' = 0 < PC1 < -30
###     'A1A1' = PC2 < -30
###     'A1A2' = 0 > PC2 >-30
###     'A2A2' = PC1 > 12; PC2 > 5

LGC6.1_2$genotype <- "RR"
LGC6.1_2[LGC6.1_2$Axis1 < 0 & LGC6.1_2$Axis1 > -30, "genotype"] <- "RA2"
LGC6.1_2[LGC6.1_2$Axis2 < -30, "genotype"] <- "A1A1"
LGC6.1_2[LGC6.1_2$Axis2 < 0 & LGC6.1_2$Axis2 > -30, "genotype"] <- "A1A2"
LGC6.1_2[LGC6.1_2$Axis1 > 12 & LGC6.1_2$Axis2 > 5, "genotype"] <- "A2A2"

ggplot(LGC6.1_2, aes(Axis1, Axis2))+
  geom_point(aes(col = genotype, pch = Species))+
  theme_bw()

### Save over the original file
write.csv(LGC6.1_2, "PCA_per_inversion/LG6_PCA_of_LGC6.1-2_v2_arc-sax.csv", quote = FALSE, row.names = FALSE)
