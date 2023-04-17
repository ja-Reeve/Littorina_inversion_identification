### PCA per map position ###
### This script uses adgenet to run a PCA for the genotypes of all individuals at each map position.
### The first principle componenet will be saved as an alternative way of detecting inversions.
### James Reeve - University of Gothenburg
### 08/04/2021
### Edited 15/12/2021: added code to save plots
### Edited 24/02/2023: clean-up for publication
### Edited 07/03/2023: fixed new bug with adegenet removing samples

### 1. Preparation ####
### Clean-up environment
rm(list = ls())
dev.off()
setwd("/Users/james/Documents/Inversion_detection")
options(stringsAsFactors = FALSE, , check.names = FALSE)

### Packages
library("tidyverse")
library("data.table")
library("adegenet")
library("ggpubr")

### Set linkage group
LG.n <- 17
# Change to whichever LG you want to analyse
LG.s <- paste0("LG", LG.n)


### 2. Access data ####

### Linkage map
LG <- read.table("ConsensusMap.v2.txt", header = T)
# Subset to current linkage group
LG <- LG[LG$LG == LG.s,]
# Add column for contig length (cM)
LG <- LG %>% group_by(contig) %>%
  mutate(gen.dist = max(avgMP) - min(avgMP))

### Function to read and wrangle genotypes
GTmap <- function(genetic.group){
  if(class(genetic.group) != "character") stop(print(paste0(cohort, " must be a string of class 'character'.")))
  if(!(genetic.group %in% c("NS", "arc", "Sp"))) stop(print(" cohort entred incorrectly, must be 'NS', 'arc', or 'Sp'!"))
  
  ## Read in genotypes
  GT <- read.table(paste0("Genotype_tables/",LG.s,"_genotypes_v2_", 
                          genetic.group, ".txt"), header = TRUE)
  
  ## Map genotypes onto linkage map by contig
  # Extract contig and map position from 'LG'
  tmp <- LG[!duplicated(LG$contig), c("contig", "avgMP")] 
  # Merge
  GT <- GT %>% left_join(tmp, by = c("CHROM" = "contig"))
  # Rearrage columns
  GT <- GT[,c(1:2, ncol(GT), 3:(ncol(GT)-1))]
  
  ## Transpose genotypes
  # Combine CHROM, POS and MapPOS into a single ID
  GT <- GT %>% unite(SNP_ID, c("CHROM", "POS", "avgMP"), sep=":") 
  # Transpose
  GT <- transpose(GT, keep.names = "SnailID", make.names = "SNP_ID")
  
  return(GT)
}

### Run 'GTmap' on each genetic group
GT_NS <- GTmap("NS")
GT_Sp <- GTmap("Sp")
GT_arc <- GTmap("arc")


### 3. PCA ####
### Vector of map positions
mp_uni <- unique(LG$avgMP)

### Function to run PCA for a single map position
PCA.genotypes <- function(data, window.start, window.end){
  # Using functions from the adgenet package
  # Taken with minor modification from Katie Hearn
  
  #A: Create window for every 1cM
  ### Extract contig, pos & mapPos from column names
  spID <- strsplit(colnames(data[,2:ncol(data)]), ":")
  contig <- sapply(spID, `[[`, 1)
  pos <- sapply(spID, `[[`, 2)
  mp_pos <- as.numeric(sapply(spID, `[[`, 3))
  
  ### Subset data to a single window
  keeps <- which(mp_pos >= window.start & mp_pos < window.end)
  tmp <- data[, keeps+1]
  
  ### Rename rows as sample ID
  rownames(tmp) <- data[, "SnailID"]
  
  ### Rename columns to just contig:pos [required for dfgenid()]
  colnames(tmp) <- paste(contig[keeps], pos[keeps], sep = ":")
  
  ### Subset data to a single window
  if(class(tmp) != "data.frame" || ncol(tmp) < 1){
    print(paste0("Warning: no SNPs detected between ", window.start, "cM and ", window.end, "cM"))
    return(NA)
  } else {
    # B: make genind object
    genind_focal <- df2genind(tmp, ploidy=2, sep="/", NA.char = "_") 	# the alleles in genotypes should be separated by '/', ie 0/0, 0/1, 1/1
    
    ### Rescaling allele frequencies is not possible using scaleGen, if all contigs at
    ### a position are monomorphic
    if(nrow(unique(genind_focal@tab)) > 1){
      # C: scale allele frequencies
      genind_scaled <- as.data.frame(scaleGen(genind_focal, NA.method="mean", scale = FALSE))
      
      # Replace any missing samples with the mean values
      if(nrow(tmp) != nrow(genind_scaled)){
        # Find missing names
        IDs <- rownames(tmp)
        kept_IDs <- rownames(genind_scaled)
        missing_IDs <- IDs[!(IDs %in% kept_IDs)]
        
        # Get mean GT score for each SNP
        GTmeans <- as.numeric(colMeans(genind_scaled))
        
        # Add 'Noffset' new rows to 'genind_scaled'
        Noffset <- nrow(data) - nrow(genind_scaled)

        for(i in 1:Noffset){ genind_scaled[nrow(genind_scaled)+i, ] <- GTmeans }
        rm(i)
        
        # Remove any extra rows that have been added
        genind_scaled <- genind_scaled[complete.cases(genind_scaled),]
        
        rownames(genind_scaled) <- c(kept_IDs, missing_IDs)
      }
      
      
      # D: the PCA
      GT_PCA <- dudi.pca(genind_scaled, cent = FALSE, scale = FALSE, nf = 2, scannf = FALSE)	
      # cent=FALSE and scale=FALSE as you centred and scaled in the previous step
      # nf = 1 to only return the first PC axis
    } else {
      GT_PCA <- dudi.pca(genind_focal, cent = FALSE, scale = FALSE, nf = 2, scannf = FALSE)
    }
    
    
    # E: save principle component scores & % of variance explained
    GT_PCs <- GT_PCA$li
    GT_var_percent <- GT_PCA$eig / sum(GT_PCA$eig) * 100
    
    # F: Write a dataframe storing the output per SNP
    res <- data.frame(cbind("LG" = LG.s,
                            "window.start" = as.numeric(window.start),
                            "window.end" = as.numeric(window.end),
                            "Snail_ID" = rownames(GT_PCs),
                            "Nsnp" = ncol(tmp),
                            "PC1" = as.numeric(GT_PCs$Axis1),
                            "PC2" = ifelse(rep("Axis2" %in% names(GT_PCs), nrow(GT_PCs)), 
                                           as.numeric(GT_PCs$Axis2),
                                           0.0), # conditional statment if only 1 PC is detected
                            "PC1_percent_var" = as.numeric(GT_var_percent[1]),
                            "PC2_percent_var" = ifelse(length(GT_var_percent) > 1,
                                                       as.numeric(GT_var_percent[2]),
                                                       0.0)
    ))
    
    return(res)
  }
}

### Run PCA.genotypes() in a loop for each window
PC_NS <- do.call(rbind.data.frame, 
                 lapply(1:round(max(mp_uni)-1), function(i){
                   PCA.genotypes(data = GT_NS,
                                 window.start = i, 
                                 window.end = i+1)}))
# Convert call of numeric variables
PC_NS[,c(2:3,5:9)] <- sapply(PC_NS[,c(2:3,5:9)], as.numeric)

PC_Sp <- do.call(rbind.data.frame, 
                  lapply(1:round(max(mp_uni)-1), function(i){
                    PCA.genotypes(data = GT_Sp,
                                  window.start = i, 
                                  window.end = i+1)}))
PC_Sp[,c(2:3,5:9)] <- sapply(PC_Sp[,c(2:3,5:9)], as.numeric)

PC_arc <- do.call(rbind.data.frame, 
                  lapply(1:round(max(mp_uni)-1), function(i){
                    PCA.genotypes(data = GT_arc,
                                  window.start = i, 
                                  window.end = i+1)}))
PC_arc[,c(2:3,5:9)] <- sapply(PC_arc[,c(2:3,5:9)], as.numeric)

### Filter any windows with < 5 SNPs
PC_NS <- PC_NS[PC_NS$Nsnp >= 5,]
PC_Sp <- PC_Sp[PC_Sp$Nsnp >= 5,]
PC_arc <- PC_arc[PC_arc$Nsnp >= 5,]

### Remove any NAs
PC_NS <- PC_NS[complete.cases(PC_NS),]
PC_Sp <- PC_Sp[complete.cases(PC_Sp),]
PC_arc <- PC_arc[complete.cases(PC_arc),]


### 4. Reorientate PC1 axis ####
### Define mid position of every window
PC_NS$pos <- apply(cbind(PC_NS$window.end, PC_NS$window.start), 1, median)
PC_Sp$pos <- apply(cbind(PC_Sp$window.end, PC_Sp$window.start), 1, median)
PC_arc$pos <- apply(cbind(PC_arc$window.end, PC_arc$window.start), 1, median)

### New vectors of map positions
mp_NS <- sort(unique(PC_NS$pos))
mp_Sp <- sort(unique(PC_Sp$pos))
mp_arc <- sort(unique(PC_arc$pos))

### Function to re-orientate PC1
# Reorientate PC1 based on correlations among windows. This function loops through all 
# windows starting from window 2. Pearson's correlation is calculated between the focal
# window and the preceeding window. If r < 0  then the sign of PC1 scores are switched.
PC_reorient <- function(PCA.data, map.positions){
  tmp <- PCA.data
  
  # Create new adjusted PC1 vector (adjPC1) and populate PC1 at first position
  tmp[tmp$pos == map.positions[1], "adjPC1"] <- tmp[tmp$pos == map.positions[1], "PC1"]
  
  # Loop through map postions
  for(i in 2:(length(map.positions))){
    pos1 <- tmp[tmp$pos == map.positions[i-1], "adjPC1"] # Adjusted PC1 at previous position
    pos2 <- tmp[tmp$pos == map.positions[i], "PC1"]      # PC1 at current position
    
    Corr <- cor(pos1, pos2) # Correlation among positions
    
    # If correlation is negative swap the sign of PC1
    if(Corr < 0) {
      tmp[tmp$pos == map.positions[i], "adjPC1"] <- pos2*-1
    } else {
      tmp[tmp$pos == map.positions[i], "adjPC1"] <- pos2
    }
  };rm(i)
  
  return(tmp)
}

PC_NS <- PC_reorient(PC_NS, mp_NS)
PC_Sp <- PC_reorient(PC_Sp, mp_Sp)
PC_arc <- PC_reorient(PC_arc, mp_arc)

### Save PCA data
write.table(PC_NS, paste0("PCA_per_window/", LG.s, "_PCA_1cMwind_v2_NS.txt"), quote = F, row.names = F)
write.table(PC_Sp, paste0("PCA_per_window/", LG.s, "_PCA_1cMwind_v2_Sp.txt"), quote = F, row.names = F)
write.table(PC_arc, paste0("PCA_per_window/", LG.s, "_PCA_1cMwind_v2_arc.txt"), quote = F, row.names = F)

### 5. Plots ####
### Access sample information
sample_info <- read.csv("Seans_WGS_sample_info_v2.csv")

### Add "arcana" and "compressa" as species names
sample_info[sample_info$Species == "arcana", "Ecotype"] <- "arcana" # rename to arcana
sample_info[sample_info$Species == "compressa", "Ecotype"] <- "compressa" # rename to compressa

# Add ecotype and country field to PCA results
PC_NS <- merge(PC_NS, sample_info, by.x = "Snail_ID", by.y = "Sample_ID")[,c(1:11,13:14)]
PC_Sp <- merge(PC_Sp, sample_info, by.x = "Snail_ID", by.y = "Sample_ID")[,c(1:11,13:14)]
PC_arc <- merge(PC_arc, sample_info, by.x = "Snail_ID", by.y = "Sample_ID")[,c(1:11,13:14)]

### Scatterplot of each window
PC_scatter <- function(data, map.positions){
  # Specify universal axis labels based on maximum values
  x_range <- range(data$PC1)
  y_range <- range(data$PC2)
  
  # Generate plot for every window
  Plots <- lapply(sort(map.positions), function(X){
    dat <- data[data$pos == X,]
    
    pt <- ggplot(dat)+
      geom_point(aes(x = PC1, y = PC2, colour = Ecotype))+
      geom_text(aes(x = 0, y = y_range[2]), label = paste(dat[,"Nsnp"], "SNPs"))+
      scale_colour_manual(values = c("pink", "lightblue", "blue", "grey", "red"))+
      lims(x = x_range, y = y_range)+
      ggtitle(paste(dat$window.start, "-", dat$window.end, "cM"))+
      theme_bw()+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = margin(0.5, 0.2, 0.2, 0, unit = "pt"),
            legend.position = "none")
    return(pt)
  })
  
  # Multi-panel plot
  PCscatt <- ggarrange(plotlist = Plots, left = "PC2", bottom = "PC1", common.legend = TRUE)
  PCscatt <- annotate_figure(PCscatt,
                             left = text_grob("PC2", face = "bold", size = 24, rot = 90),
                             bottom = text_grob("PC1", face = "bold", size = 24))
  
  return(PCscatt)
}

p1 <- PC_scatter(PC_NS, mp_NS[seq(1, length(mp_NS),1)])
p2 <- PC_scatter(PC_Sp, mp_Sp[seq(1, length(mp_Sp),1)])
p3 <- PC_scatter(PC_arc, mp_arc[seq(1, length(mp_arc),1)])

# Save output
PATH <- "/Users/james/Documents/Inversion_detection/plots/PCA_per_mapPos_plots"
ggsave(paste0(LG.s, "PCAperPos_scatter_sax.v5.tiff"), p1, "tiff", PATH, width = 30, height = 30, units = "cm")
ggsave(paste0(LG.s, "PCAperPos_scatter_Sp.v5.tiff"), p2, "tiff", PATH, width = 30, height = 30, units = "cm")
ggsave(paste0(LG.s, "PCAperPos_scatter_arc.v5.tiff"), p3, "tiff", PATH, width = 30, height = 30, units = "cm")

### Plot PC1 along linkage map
PC_map <- function(data, genetic.group){
  if(genetic.group == "NS"){
    EcoShape <- c(17, 16, 15, 16, 17)
    EcoColour <- c("pink", "lightblue", "blue", "grey", "red")
    plot.title <- "Northern saxatilis"}
  
  if(genetic.group == "Sp"){
    EcoShape <- c(15, 16, 17)
    EcoColour <- c("blue", "grey", "red")
    plot.title <- "Iberian saxatilis"}
  
  if(genetic.group == "arc"){
    EcoShape <- 16
    EcoColour <- "grey"
    plot.title <- "Littorina arcana"}
  
  # Plot function
  ggplot(data, aes(x = pos, y = adjPC1, alpha = PC1_percent_var, 
                   colour = Ecotype, group = Snail_ID))+
    geom_point(aes(pch = Ecotype))+
    geom_line()+
    scale_shape_manual(values = EcoShape)+
    scale_colour_manual(values = EcoColour)+
    labs(x = "Linkage map position (cM)", y = "Reorientated PC1",
         title = paste0(LG.s, ": ", plot.title))+
    theme_bw()
}

p1 <- PC_map(PC_NS, "NS")
p2 <- PC_map(PC_Sp, "Sp")
p3 <- PC_map(PC_arc, "arc")

# Multi-panel plot
ggarrange(p1,p2,p3, nrow = 3, ncol = 1, common.legend = TRUE, legend = "right")

# Save output
ggsave(paste0(LG.s, "PCAperPos_lines_sax.v5.tiff"), p1, "tiff", PATH, width = 30, height = 15, units = "cm")
ggsave(paste0(LG.s, "PCAperPos_lines_Sp.v5.tiff"), p2, "tiff", PATH, width = 30, height = 15, units = "cm")
ggsave(paste0(LG.s, "PCAperPos_lines_arc.v5.tiff"), p3, "tiff", PATH, width = 30, height = 15, units = "cm")
