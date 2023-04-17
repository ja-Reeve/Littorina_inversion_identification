######################  Genotyping Inversion #####################
### This script takes find clusters of points from the PCA scores
### of each inversion (see PCA_inversion-area.R). Clusters are determined
### by K means clustering, trailing 2 ≤ K ≤ 6. The best K is determined
### using the Silhouette method to calculate the maximum distance between
### clusters. For simple inversions only PC1 is used for clustering, and 
### max(K) = 3. For complex inversions (LGC6.1-2 & LGC14.1-2) both PC1 &
### PC2 are used for clustering. Finally, clusters were relabelled so that
### 'R' is the cluster with the highest frequency of crab ecotype snails.
### James Reeve - University of Gothenburg
### 18/01/2022

### Edit 14/09/2022: change assignment of genotypes so it's now based on
###                  the cluster with the most crab snails
### Edit 13/03/2023: update to analysis re-run and tweak assignment of cluster
###                  for complex inversions

### Preparation
rm(list = ls())
dev.off()
setwd("/Users/james/Documents/Inversion_detection/PCA_per_inversion")

### Packages
library("tidyverse")
library("ggpubr")
library("cluster")

### Function to genotype inversions
inversion.genotype <- function(inversion, genetic.group){
  ### A: Access PCA results
  # Get focal linkage group from inversion name
  LG <- gsub("C", "", gsub("[.].*", "", inversion))
  # Upload the PCA
  PC <- read.csv(paste0(LG, "_PCA_of_", inversion, "_v2_", genetic.group, ".csv"))
  
  
  ### B: Select the best K using Silouette method
  if(inversion == "LGC6.1-2" | inversion == "LGC14.1-2" | inversion == "LGC14.2"){
    # Run Kmeans from K = 2:6 and calculate the average silhouette score 
    # for the complex inversions
    tmp <- lapply(2:9, function(K){
      # Kmeans clustering
      clust <- kmeans(PC[,c("Axis1", "Axis2")], centers = K, iter.max = 1e5, nstart = 100)
      # Calcualte average silhouette score (i.e. difference between clusters)
      sil <- mean(silhouette(clust$cluster, dist(PC[,c("Axis1", "Axis2")]))[,"sil_width"])
      # Store clusters and silhouette values in output
      res <- append(clust, sil)
      names(res)[10] <- "Avg_Silhouette_Weight"
      return(res)
    })
  } else {
    # Run Kmeans from K = 2:3 for all other inversions
    tmp <- lapply(2:3, function(K){
      # Kmeans clustering
      clust <- kmeans(PC$Axis1, centers = K, iter.max = 1e5, nstart = 100)
      # Calcualte average silouette score (i.e. difference between clusters)
      sil <- mean(silhouette(clust$cluster, dist(PC$Axis1))[,"sil_width"])
      # Store clusters and silhouette values in output
      res <- append(clust, sil)
      names(res)[10] <- "Avg_Silhouette_Weight"
      return(res)
    })
  }
  
  # C: Find the best K (highest Silouette score)
  K <- which.max(sapply(tmp, `[[`, "Avg_Silhouette_Weight")) + 1
  
  # D: Reassign cluster names based on centre values
  bestK <- tmp[[K - 1]]
  
  if(K == 3){
    ### Identify middle cluster
    min_clust <- which.min(bestK$centers) # Cluster with lowest centre
    max_clust <- which.max(bestK$centers) # Cluster with highest centre
    bestK$cluster[bestK$cluster != min_clust & bestK$cluster != max_clust] <- "RA"
    
    ### Add clusters to PC
    PC$genotype <- bestK$cluster
    
    ### Relabel cluster with most crab snails 'RR'
    crab_clust <- PC[PC$Species == "saxatilis" & PC$Ecotype == "Crab", "genotype"]
    uClust <- unique(crab_clust[crab_clust != "RA"])
    
    # Find frequency of crab individuals in each homokaryotype cluster
    pClust <- sapply(as.numeric(uClust), function(i){sum(crab_clust == i) / bestK$size[i]})
    
    # Assign cluster with highest crab frequency as 'RR'
    RR <- uClust[which.max(pClust)]
    
    ### Assign homokaryotype clusters
    PC$genotype[PC$genotype == RR] <- "RR"
    PC$genotype[!(PC$genotype %in% c("RR", "RA"))] <- "AA"
    
  }
  if(K == 6){
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
    
    ### Add clusters to PC
    PC$genotype <- bestK$cluster
    
    ### Subset to just homokaryotypics (extreme) clusters
    PC_homo <- PC[PC$genotype %in% c(clust_A, clust_C, clust_F), ]
    
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
    
    
    PC$genotype[PC$genotype == RR] <- "RR"
    PC$genotype[PC$genotype == A1A1] <- "A1A1"
    PC$genotype[PC$genotype == A2A2] <- "A2A2"
    PC$genotype[PC$genotype == RA1] <- "RA1"
    PC$genotype[PC$genotype == RA2] <- "RA2"
    PC$genotype[PC$genotype == A1A2] <- "A1A2"
    
  }
  
  return(PC)
}

### Vector of inversion names
INVs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC5.1", "LGC6.1-2", "LGC7.1", "LGC7.2",
          "LGC9.1", "LGC9.2", "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", "LGC12.2",
          "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.2", #"LGC14.1-2",
          "LGC14.3", "LGC17.1")


#### Plot Northern saxatilis ####

# Run genotyping function
Nsax <- lapply(INVs, inversion.genotype, genetic.group = "NS")
names(Nsax) <- INVs # Give list the names of each inversion

# Plot genotypes on PCA for each inversion that had 3 or 6 clusters
NsaxP <- lapply(1:length(INVs), function(i){
  if(class(Nsax[[i]]$genotype) == "character"){ # Filter to K = 3 | 6
    ggplot(Nsax[[i]], aes(Axis1, Axis2, colour = factor(genotype)))+
      geom_point()+
      coord_fixed()+
      ggtitle(INVs[i])+
      theme_bw()+
      theme(axis.title = element_blank())
  } else {NA}
})

# Plot all inversions
annotate_figure(
  ggarrange(plotlist = NsaxP, common.legend = TRUE, legend = "right"), 
  top = text_grob("Northern Littorina saxatilis", size = 24),
  left = text_grob("PC2", size = 18, face = "bold", rot = 90), 
  bottom = text_grob("PC1", size = 18, face = "bold"))

#### Plot Spanish saxatilis ####

#Ssax <- lapply(INVs, inversion.genotype, genetic.group = "Sp")
#names(Ssax) <- INVs

#SsaxP <- lapply(1:length(INVs), function(i){
#  if(class(Ssax[[i]]$genotype) == "character"){
#    ggplot(Ssax[[i]], aes(Axis1, Axis2, colour = genotype))+
#      geom_point()+
#      coord_fixed()+
#      ggtitle(INVs[i])+
#      theme_bw()+
#      theme(axis.title = element_blank())
#  } else {NA}
#})

#annotate_figure(ggarrange(plotlist = SsaxP, common.legend = TRUE, legend = "right"), 
#                top = text_grob("Iberian Littorina saxatilis", size = 24),
#                left = text_grob("PC2", size = 18, face = "bold", rot = 90), 
#                bottom = text_grob("PC1", size = 18, face = "bold"))

#### Plot Littorina arcana ####
### Note: this dosen't work, as there were no 'crab' ecotypes in Larc
#Larc <- lapply(INVs, inversion.genotype, cohort = "arc")
#names(Larc) <- INVs

#LarcP <- lapply(1:length(INVs), function(i){
#  if(class(Larc[[i]]$genotype) == "character"){
#    ggplot(Larc[[i]], aes(Axis1, Axis2, colour = genotype))+
#      geom_point()+
#      coord_fixed()+
#      ggtitle(INVs[i])+
#      theme_bw()+
#      theme(axis.title = element_blank())
#  } else {NA}
#})

#annotate_figure(ggarrange(plotlist = LarcP, common.legend = TRUE, legend = "right"), 
#                top = text_grob("Littorina arcana", size = 24),
#                left = text_grob("PC2", size = 18, face = "bold", rot = 90), 
#                bottom = text_grob("PC1", size = 18, face = "bold"))

#### Plot Littorina saxatilis & arcana ####

LarcSax <- lapply(INVs, inversion.genotype, genetic.group = "arc-sax")
names(LarcSax) <- INVs

LarcSaxP <- lapply(1:length(INVs), function(i){
  if(class(LarcSax[[i]]$genotype) == "character"){
    ggplot(LarcSax[[i]], aes(Axis1, Axis2, colour = genotype))+
      geom_point()+
      coord_fixed()+
      ggtitle(INVs[i])+
      theme_bw()+
      theme(axis.title = element_blank())
  } else {NA}
})

annotate_figure(ggarrange(plotlist = LarcSaxP, common.legend = TRUE, legend = "right"), 
                top = text_grob("Littorina arcana & saxatilis", size = 24),
                left = text_grob("PC2", size = 18, face = "bold", rot = 90), 
                bottom = text_grob("PC1", size = 18, face = "bold"))

#### Plot full dataset ####
#ALL <- lapply(INVs, inversion.genotype, cohort = "all")
#names(ALL) <- INVs

#ALLP <- lapply(1:length(INVs), function(i){
#  if(class(ALL[[i]]$genotype) == "character"){
#    ggplot(ALL[[i]], aes(Axis1, Axis2, colour = genotype, pch = Ecotype == "compressa"))+
#      geom_point()+
#      coord_fixed()+
#      ggtitle(INVs[i])+
#      theme_bw()+
#      theme(axis.title = element_blank())
#  } else {NA}
#})

#annotate_figure(
#  ggarrange(plotlist = ALLP, common.legend = TRUE, legend = "right"), 
#  top = text_grob("All snails", size = 24),
#  left = text_grob("PC2", size = 18, face = "bold", rot = 90), 
#  bottom = text_grob("PC1", size = 18, face = "bold"))


#### Save over the PCA files ####
### Note: this is not best practices, but since this script only adds 1 column to
### the data, I deem it better to save over the old file than clog up my harddrive.

### Save Northern saxatilis
sapply(names(Nsax), function(inv){
  LG <- gsub("C", "", gsub("[.].*", "", inv))
  write.csv(Nsax[[inv]], paste0(LG, "_PCA_of_", inv, "_v2_NS.csv"),
            quote = FALSE, row.names = FALSE)
})

### Save Spanish saxatilis
#sapply(names(Ssax), function(inv){
#  LG <- gsub("C", "", gsub("[.].*", "", inv))
#  write.csv(Ssax[inv], paste0(LG, "_PCA_of_", inv, "_v2_Sp.csv"),
#            quote = FALSE, row.names = FALSE)
#})

### Save Littorina arcana
#sapply(names(Larc), function(inv){
#  LG <- gsub("C", "", gsub("[.].*", "", inv))
#  write.csv(Larc[inv], paste0(LG, "_PCA_of_", inv, "_v2_arc.csv"),
#            quote = FALSE, row.names = FALSE)
#})

### Save sax arc data
sapply(names(LarcSax), function(inv){
  LG <- gsub("C", "", gsub("[.].*", "", inv))
  write.csv(LarcSax[[inv]], paste0(LG, "_PCA_of_", inv, "_v2_arc-sax.csv"),
            quote = FALSE, row.names = FALSE)
})

### Save all snail data
#sapply(names(ALL), function(inv){
#  LG <- gsub("C", "", gsub("[.].*", "", inv))
#  write.csv(ALL[inv], paste0(LG, "_PCA_of_", inv, "_v2_all.csv"),
#            quote = FALSE, row.names = FALSE)
#})
