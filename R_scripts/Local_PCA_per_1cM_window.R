### PCA per map position ###
### This script uses adgenet to run a PCA for the genotypes of all individuals at each map position.
### The first principle componenet will be saved as an alternative way of detecting inversions.
### James Reeve - University of Gothenburg
### 08/04/2021
### Eidted 15/12/2021: added code to save plots

### 1. Preparation ####
### Clean-up environment
rm(list = ls())
dev.off()
setwd("")
options(stringsAsFactors = FALSE)

### Packages
library("tidyverse")
library("data.table")
library("adegenet")
library("ggpubr")

### Set linkage group
LG.n <- 1
# Change to whichever LG you want to analyse
LG.s <- paste0("LG", LG.n)

### Linkage map
LG <- read.table("ConsensusMap.v2.txt", header = T)
# Subset to current linkage group
LG <- LG[LG$LG == LG.s,]
# Extract contig and physical position
#tmp <- strsplit(LG$Marker, "_")
#LG$contig <- sapply(tmp, `[[`, 1)
#LG$physical.position <- as.numeric(sapply(tmp, `[[`, 2))
#rm(tmp)
# Add column for contig length (cM)
LG <- LG %>% group_by(contig) %>%
  mutate(gen.dist = max(av) - min(av))

### Upload sample information and create IDs for each cohort
sample_info <- read.csv("/PATH/Sample_information.csv", 
                        sep=";", check.names = FALSE)

# Seperate species
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

### Genotypes
GT <- read.table(paste0("Genotype_tables/",LG.s,"_genotypes.v2.txt"), header = TRUE)
colnames(GT) <- names(read.delim("vcf_header.txt", check.names = FALSE))[c(1:2, 10:117)] # 10:86 for sax-only data
colnames(GT)[1] <- "CHROM"
GT_sax <- GT %>% select("CHROM", "POS", all_of(sax_IDs))
GT_Sp <- GT_sax %>% select("CHROM", "POS", all_of(Sp_IDs))
GT_sax <- GT_sax %>% select("CHROM", "POS", all_of(c(ES_IDs, FW_IDs, NR_IDs, O_IDs)))

GT_arc <- read.table(paste0("Genotype_tables/",LG.s,"_genotypes_arc.txt"), header = TRUE)

### Map GT onto LG
# Index of map positions in genotype table
GTmap <- function(data) {
  tmp <- apply(data[, 1:2], 1, function(X){
    cont <- X[1]
    pos <- as.numeric(X[2])
    
    tmp2 <- LG[which(LG$contig %in% cont),]                  # Match by contig
    
    if(nrow(tmp2) > 0){
      tmp3 <- data.frame("CHROM" = cont,                     # Contig
                         "POS" = pos,                        # Physical position
                         "MapPOS" = ifelse(nrow(tmp2)>1, 
                                           median(tmp2$av),
                                           tmp2$av))         # Map position
    return(tmp3)}
  
  })
  GTmap <- do.call(rbind.data.frame, tmp)
  
  # Add merge position and chromosome into LOC flag
  GTmap$LOC <- paste(GTmap$CHROM, GTmap$POS, sep = ":")
  GTmap$LOC <- gsub(" ", "", GTmap$LOC) # Fixes a bug with white space in LOC
  
  # Add LOC flag to original data
  tmp3 <- data
  tmp3$LOC <- paste(tmp3$CHROM, tmp3$POS, sep = ":")
  
  # Merge original data with GTmap
  tmp4 <- merge(tmp3, GTmap[,3:4], by = "LOC")[,c(2:3, ncol(tmp3)+1, 4:ncol(tmp3))]
  
  return(tmp4)
}

GT_map_sax <- GTmap(GT_sax)
GT_map_Sp <- GTmap(GT_Sp)
GT_map_arc <- GTmap(GT_arc)

### 2. Transpose genotype table ####
transpose.genotypes <- function(data){
  # Combine CHROM, POS & MapPOS into "SNP_ID"
  GT_tmp <- data %>% unite(SNP_ID, c("CHROM", "POS", "MapPOS"), sep=":")
  
  # Transpose using function from data.table package
  GT_trans <- transpose(GT_tmp, keep.names = "SnailID", make.names = "SNP_ID")
  
  return(GT_trans)
}

GT_tran_sax <- transpose.genotypes(GT_map_sax)
GT_tran_Sp <- transpose.genotypes(GT_map_Sp)
GT_tran_arc <- transpose.genotypes(GT_map_arc)

### 3. PCA ####
### Vector of map positions
mp_uni <- unique(GT_map_sax$MapPOS)

### Function to run PCA for a single map position
PCA.genotypes <- function(data, window.start, window.end){
  # Using functions from the adgenet package
  # Taken with minor modification from Katie Hearn
  
  #A: Create window for every 1cM
  ### Extract contig, pos & mapPos from column names
  spID <- strsplit(colnames(data), ":")
  contigs <- sapply(spID, `[[`, 1)
  pos <- sapply(spID, `[[`, 2)
  mp_pos <- as.numeric(sapply(spID, `[[`, 3))
  
  ### Rename columns to just contig:pos [required for dfgenid()]
  tmp <- data
  colnames(tmp) <- paste(contigs, pos, sep = ":")
  
  ### Subset data to a single window
  tmp2 <- tmp[, which(mp_pos >= window.start & mp_pos < window.end)]
  
  if(class(tmp2) != "data.frame" || ncol(tmp2) < 1){
    print(paste0("Warning: no SNPs detected between ", window.start, "cM and ", window.end, "cM"))
    return(NA)
  } else {
    # B: make genind object
    genind_focal <- df2genind(tmp2, ploidy=2, sep="/", NA.char = "./.") 	# the alleles in genotypes should be separated by '/', ie 0/0, 0/1, 1/1
    
    ### Recaling allele frequencies is not possible using scaleGen, if all contigs at
    ### a position are monomorphic
    if(nrow(unique(genind_focal@tab)) > 1){
      # C: scale allele frequencies
      genind_scaled <- scaleGen(genind_focal, NA.method="mean")
      
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
                            "Nsnp" = ncol(tmp2),
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

PC_sax <- do.call(rbind.data.frame, 
                 lapply(1:round(max(mp_uni)-1), function(i){
                   PCA.genotypes(data = GT_tran_sax,
                                 window.start = i, 
                                 window.end = i+1)}))
PC_sax[,c(2:3,5:9)] <- sapply(PC_sax[,c(2:3,5:9)], as.numeric)

PC_Sp <- do.call(rbind.data.frame, 
                  lapply(1:round(max(mp_uni)-1), function(i){
                    PCA.genotypes(data = GT_tran_Sp,
                                  window.start = i, 
                                  window.end = i+1)}))
PC_Sp[,c(2:3,5:9)] <- sapply(PC_Sp[,c(2:3,5:9)], as.numeric)

PC_arc <- do.call(rbind.data.frame, 
                  lapply(1:round(max(mp_uni)-1), function(i){
                    PCA.genotypes(data = GT_tran_arc,
                                  window.start = i, 
                                  window.end = i+1)}))
PC_arc[,c(2:3,5:9)] <- sapply(PC_arc[,c(2:3,5:9)], as.numeric)

### Filter any windows with < 5 SNPs
PC_sax <- PC_sax[PC_sax$Nsnp >= 5,]
PC_Sp <- PC_Sp[PC_Sp$Nsnp >= 5,]
PC_arc <- PC_arc[PC_arc$Nsnp >= 5,]

### Remove any NAs
PC_sax <- PC_sax[complete.cases(PC_sax),]
PC_Sp <- PC_Sp[complete.cases(PC_Sp),]
PC_arc <- PC_arc[complete.cases(PC_arc),]

### 4. Reorientate PC1 axis ####
### Define mid position of every window
PC_sax$pos <- apply(cbind(PC_sax$window.end, PC_sax$window.start), 1, median)
PC_Sp$pos <- apply(cbind(PC_Sp$window.end, PC_Sp$window.start), 1, median)
PC_arc$pos <- apply(cbind(PC_arc$window.end, PC_arc$window.start), 1, median)

### New vectors of map positions
mp_sax <- sort(unique(PC_sax$pos))
mp_Sp <- sort(unique(PC_Sp$pos))
mp_arc <- sort(unique(PC_arc$pos))

### Function to re-orientate PC1
# Reorientate PC1 based on correlations among windows. This function loops through all 
# windows starting from window 2. Pearson's correlation is calculated between the focal
# window and the preceeding window. If r < 0  then the sign of PC1 scores are switched.
PC_reorient <- function(PCA.data, map.positions){
  tmp <- PCA.data
  
  # Create new adjusted PC1 vector (adjPC1) and populatate PC1 at first position
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

PC_sax <- PC_reorient(PC_sax, mp_sax)
PC_Sp <- PC_reorient(PC_Sp, mp_Sp)
PC_arc <- PC_reorient(PC_arc, mp_arc)

### Save PCA data
write.table(PC_sax, paste0("PCA_per_window/", LG.s, "_PCA_1cMwind_sax.txt"), quote = F, row.names = F)
write.table(PC_Sp, paste0("PCA_per_window/", LG.s, "_PCA_1cMwind_Sp.txt"), quote = F, row.names = F)
write.table(PC_arc, paste0("PCA_per_window/", LG.s, "_PCA_1cMwind_arc.txt"), quote = F, row.names = F)

### 5. Plots ####
### Changing some ecotype labels to reduce number of ecotypes
sample_info[sample_info$Ecotype %in% c("Crab", "crabish"), "Ecotype"] <- "crab" # Rename crab ecotype
sample_info[sample_info$Ecotype == "Wave", "Ecotype"] <- "wave" # Rename wave ecotype
sample_info[sample_info$Ecotype == "barnacle(ish)", "Ecotype"] <- "barnacle" # Rename barnacle ecotype
sample_info[sample_info$Ecotype %in% c("", "?", "midshore (no ecotypes)"), "Ecotype"] <- "other" # Rename everything else
# "arcana" and "compressa" are added as ecotypes
sample_info[sample_info$Species == "arcana", "Ecotype"] <- "arcana" # rename to arcana
sample_info[sample_info$Species == "compressa", "Ecotype"] <- "compressa" # rename to compressa

# Add ecotype and country field to PCA results
PC_sax <- merge(PC_sax, sample_info, by.x = "Snail_ID", by.y = "Sample_ID")[,c(1:11,13:14)]
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
      scale_colour_manual(values = c("pink", "blue", "grey", "lightblue", "red"))+
      lims(x = x_range, y = y_range)+
      ggtitle(paste(round(X, 2), "cM"))+
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

p1 <- PC_scatter(PC_sax, mp_sax[seq(1, length(mp_sax),3)])
p2 <- PC_scatter(PC_Sp, mp_Sp[seq(1, length(mp_Sp),3)])
p3 <- PC_scatter(PC_arc, mp_arc[seq(1, length(mp_arc),3)])

# Save output
PATH <- "/Users/james/Documents/Inversion_detection/plots/PCA_per_mapPos_plots"
ggsave(paste0(LG.s, "PCAperPos_scatter_sax.v4.tiff"), p1, "tiff", PATH, width = 30, height = 30, units = "cm")
ggsave(paste0(LG.s, "PCAperPos_scatter_Sp.v4.tiff"), p2, "tiff", PATH, width = 30, height = 30, units = "cm")
ggsave(paste0(LG.s, "PCAperPos_scatter_arc.v4.tiff"), p3, "tiff", PATH, width = 30, height = 30, units = "cm")

### Plot PC1 along linkage map
PC_map <- function(data, plot.title){
  ggplot(data, aes(x = pos, y = adjPC1, alpha = PC1_percent_var, 
                   colour = Ecotype, group = Snail_ID))+
    geom_point(aes(pch = Ecotype))+
    geom_line()+
    scale_shape_manual(values = c(17, 16, 15, 16, 17))+
    scale_colour_manual(values = c("pink", "blue", "grey", "lightblue", "red"))+
    labs(x = "Linkage map position (cM)", y = "Reorientated PC1",
         title = paste0(LG.s, ": ", plot.title))+
    theme_bw()
}

p1 <- PC_map(PC_sax, "Northern saxatilis")
p2 <- PC_map(PC_Sp, "Iberian saxatilis")
p3 <- PC_map(PC_arc, "Littorina arcana")

# Multi-panel plot
ggarrange(p1,p2,p3, nrow = 3, ncol = 1, common.legend = TRUE, legend = "right")

# Save output
PATH <- "/Users/james/Documents/Inversion_detection/plots/PCA_per_mapPos_plots"
ggsave(paste0(LG.s, "PCAperPos_lines_sax.v4.tiff"), p1, "tiff", PATH, width = 30, height = 15, units = "cm")
ggsave(paste0(LG.s, "PCAperPos_lines_Sp.v4.tiff"), p2, "tiff", PATH, width = 30, height = 15, units = "cm")
ggsave(paste0(LG.s, "PCAperPos_lines_arc.v4.tiff"), p3, "tiff", PATH, width = 30, height = 15, units = "cm")
