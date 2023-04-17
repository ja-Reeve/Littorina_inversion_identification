### PCA of inversions ###
### This script was modified version designed by Katie Hearn from the University of Shefield.
### It converts a genotype file formates as: CHROM, POS, SAMPLES; into a PCA looking for 
### clusters that indicate inversion karyotypes. PCAs are also run for parts of the genome
### outside inversions.
### James Reeve - University of Gothenburg
### 30/11/2020
### Edited 02/09/2021 to adjust to new consensus map
### Edit 17/12/2021: version 2 created with new name. Many data handling steps are moved to
### functions to make the code a bit more efficient to run.
### Edit 14/01/2022: Clean up functions and recolour plots.
### Edit 26/08/2022: Remove compressa from the analysis and change PCA plot colouring
### Edit 08/09/2022: alter inversion boundaries for each cohort
### Edit 07/10/2022: add compressa back in, but project them onto the PCA of other samples
### Edit 10/03/2023: clean-up for publication

#### 1. Preparation ####
### Clean-up environment
rm(list = ls())
dev.off()
setwd("/Users/james/Documents/Inversion_detection")
options(stringsAsFactors = FALSE, check.names = FALSE)

Start_time <- Sys.time()

### Packages
library("tidyverse")
library("ggpubr")
library("data.table")
library("adegenet")

### Set linkage group
LG.n <- 17
# Change to whichever LG you want to analyse
LG.s <- paste0("LG", LG.n)

### Linkage map
LG <- read.table("ConsensusMap.v2.txt", header = T)
LG <- LG[LG$LG == LG.s,]

### Find inversions on current LG
Invs <- read.csv("Inversion_breakpoint_estimates_v2.csv", header = TRUE)
Invs <- Invs[Invs$LG == LG.s, c("Inversion", "start_pos", "end_pos", "genetic_group")]

### Access collection details
sample_info <- read.csv("Seans_WGS_sample_info_v2.csv")

# Create new column for phylogenetic regions
sample_info$Region <- NA
sample_info[sample_info$Country %in% c("England", "Sweden") & sample_info$location != "Dersingham" |
              sample_info$location == "St Abbs", "Region"] <- "North_Sea"
sample_info[sample_info$Country %in% c("Norway", "Russia"), "Region"] <- "Northern_Atlantic"
sample_info[sample_info$Country %in% c("Wales", "Ireland", "Isle of man", "France", "Iceland") |
                        sample_info$location %in% c("Oban", "Dersingham"), "Region"] <- "Celtic_Sea"
sample_info[sample_info$Country == "Spain", "Region"] <- "Iberia"
sample_info[sample_info$Country == "USA", "Region"] <- "Western_Atlantic"

# Create a vector of compressa names
# This is required to project compressa onto the PCA axes of arcana and saxatilis
comp_IDs <- sample_info[sample_info$Species == "compressa", "Sample_ID"]

# Change "Ecotype" labels in arcana and compressa to the species name
sample_info[sample_info$Species == "arcana", "Ecotype"] <- "arcana"
sample_info[sample_info$Species == "compressa", "Ecotype"] <- "compressa"


#### 3. Function to run PCA for inversion for each genetic group ####
### This is a complex process with several computationally intensive steps.
### It's more efficient to run this within a function so temporary files aren't
### clogging up the RAM. This function is comprised of several sub-steps:
###       3.1: Filtering genotypes - removing any SNPs which are fixed in a gentic group
###             or SNPs that only appear once.
###       3.2: Transposing the genotype table - needed for the PCA
###       3.3: Filter SNPs to focal region - this is two functions A) is for inversions
###            B) is for the collinear genome.
###       3.4: PCA - run the PCA analysis
### I've set it up in a way that it should only need two input flags, and one optional flags:
###       inversion = name of inversion as a string. Can also be "collinear" to select collinear region.
###       geneitc.group = string of 4 values {"NS", "Sp", "arc", "ALL"}
###       buffer = buffer of XcM around inversion breakpoints
PCA.genotypes <- function(inversion, genetic.group, buffer = 0){
  ### Error messages
  if(class(genetic.group) != "character") stop(print(paste0(genetic.group, " must be a string of class 'character'.")))
  if(!(genetic.group %in% c("NS", "arc", "Sp", "ALL", "arc_sax"))) stop(print(" genetic group entred incorrectly, must be 'NS', 'arc', 'Sp', 'arc_sax', or 'ALL'!"))
  if(class(inversion) != "character") stop(print(" inversion must be a string value."))
  if(class(buffer) != "numeric") stop("Buffer must be a numeric value.")
  
  ### 3.1: Load in genotypes
  if(genetic.group == "ALL"){
    
    ## Load full genotype table
    GT <- read.table(paste0("Genotype_tables/",LG.s,"_genotypes.v2.txt"), header = FALSE, check.names = FALSE)
    colnames(GT) <- names(read.delim("vcf_header.txt", check.names = FALSE))[c(1:2, 10:117)]
    colnames(GT)[1] <- "CHROM" #rename 1st column
    # Rename 1 odd Larc ID
    names(GT)[names(GT) == "W_arc_04_Lamerged_sorted.bam"] <- "W_arc_04_La"
    # Drop "IMI_6_2" due to missing sampling information
    GT <- GT[,colnames(GT) != "IMI_6_2"]
  } else {
    
    ## Load genotypes filtered by site
    GT <- read.table(paste0("Genotype_tables/",LG.s,"_genotypes_v2_", genetic.group,".txt"), header = T, check.names = F)
  }
  
  
  ### 3.2: Filter to target region
  # Find contigs in target region
  if(inversion == "collinear"){
    
    ### Target region == collinear genome
    # List of all map positions
    map_pos <- unique(LG$avgMP)
    # Subset 'Invs' by genetic.group
    IBp <- Invs[Invs$genetic_group == ifelse(genetic.group %in% c("arc_sax", "ALL"), "NS", genetic.group), ]
    # List of map positions on inversions
    inv_pos <- sapply(1:nrow(IBp), function(X){
      map_pos[map_pos >= IBp$start_pos[X] - buffer & map_pos <= IBp$end_pos[X] + buffer]}) 
    # Collapse list
    inv_pos <- unlist(inv_pos) 
    # Filter out all inversion positions on linkage map
    contigs <- LG[!(LG$avgMP %in% inv_pos), "contig"] 
    
  } else {
    
    ### Target region == inversion
    # Breakpoints
    Spos <- Invs[Invs$Inv == inversion & 
                   Invs$genetic_group == ifelse(genetic.group %in% c("arc_sax", "ALL"), "NS", genetic.group), 
                 "start_pos"] # Start of inversion
    Epos <- Invs[Invs$Inv == inversion & 
                   Invs$genetic_group == ifelse(genetic.group %in% c("arc_sax", "ALL"), "NS", genetic.group), 
                 "end_pos"] # End of inversion
    # Filter linkage map to inversion area
    contigs <- LG[LG$avgMP >= Spos + buffer & LG$avgMP <= Epos - buffer, "contig"]
    
  }
  
  # Subset transformed genotypes to those in target region
  GT_filt <- GT[GT$CHROM %in% contigs, ]
  

  ### 3.3: Transpose the genotype table
  # Combine CHROM and POS into "SNP_ID"
  tmp <- GT_filt %>% unite(SNP_ID, c("CHROM", "POS"), sep=":")
  # Transpose using function from data.table package
  GT_trans <- transpose(tmp, keep.names = "SnailID", make.names = "SNP_ID")
  
  rm(tmp)
  
  
  ### 3.4: Run the PCA
  # Error check if contigs exist - e.g. some inversions cover the whole LG leaving no collinear genome
  if(ncol(GT_trans) == 0) return(NA)
  # Using functions from the adgenet package
  # Taken with minor modification from Katie Hearn
  # A: make genind object
  genind <- df2genind(GT_trans[,-c(1)], ploidy=2, sep="/", NA.char = "_")
  # B: re-score missing genotypes as the mean values
  genind <- scaleGen(genind, NA.method="mean", scale = FALSE)
  rownames(genind) <- colnames(GT)[-c(1:2)]
  # C: PCA call
  # For the full data set, compressa are projected onto the PCA to avoid clustering artifact with Iberian snails
  if(genetic.group == "ALL"){
    #C.1: extract compressa
    genind_comp <- genind[rownames(genind) %in% comp_IDs,]
    #C.2: run PCA on genind without compressa
    GT_PCA <- dudi.pca(genind[!(rownames(genind) %in% comp_IDs),], 
                       cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)	
    #C.3: project compressa genotypes onto the PCA
    comp_proj <- suprow(GT_PCA, genind_comp)
    # D: save principle component scores & % of variance explained
    GT_PCs <- rbind.data.frame(GT_PCA$li, comp_proj[[2]])
    GT_var_percent <- GT_PCA$eig / sum(GT_PCA$eig) * 100
  } else {
    #C.1: run PCA
    GT_PCA <- dudi.pca(genind, cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)	
    # D: save principle component scores & % of variance explained
    GT_PCs <- GT_PCA$li
    GT_var_percent <- GT_PCA$eig / sum(GT_PCA$eig) * 100
  }
  # E: Store both results in a list
  PC.dat <- list(GT_var_percent[1:3], GT_PCs)
  # F: Add sample names as new column
  PC.dat[[2]]$Sample_ID <- rownames(PC.dat[[2]])
  
  
  ### 3.5: Add heterozygosity to PCA
  Het <- apply(GT_trans[,-c(1)], 1, function(X){
    tmp <- X[X != "./."]
    He <-  sum(tmp == "0/1") / length(tmp)
    
    return(He)
  })
  
  PC.dat[[2]]$Heterozygosity <- Het
  
  
  ### 3.6: Add sample information to PCA
  PC.dat[[2]] <- PC.dat[[2]] %>% left_join(., sample_info, by = "Sample_ID") 
  
  return(PC.dat)
}


#### 4. Function to plot PCA ####
### This step creates scatterplots for PC1 vs PC2. Three scatterplots are
### made coloured by different aspects of the data; 1) heterozygosity,
### 2) country and 3) ecotype. There is also two options to specify the
### range on both axes. This is turned off by default, but it useful when
### making a multiple plots to unify their axis limits.
plot.PCA <- function(PCA.data, inversion, genetic.group){#, x.range = NULL, y.range = NULL){
  ### Skip function if PCA didn't work
  if(class(PCA.data) != "list") stop(paste0("PCA did not work for ", inversion, "_", genetic.group))
  
  ### Modify inversion names for saving plots due to error from '/'
  inv2 <- ifelse(grepl("1/2", inversion), gsub("/", "-", inversion), inversion)
  
  ### Create plot titles
  if(genetic.group == "NS") title.text <- paste0(inversion, ": \nNorthern Littorina saxatilis")
  if(genetic.group == "Sp") title.text <- paste0(inversion, ": \nSpanish Littorina saxatilis")
  if(genetic.group == "arc") title.text <- paste0(inversion, ": \nLittorina arcana")
  if(genetic.group == "arc_sax") title.text <- paste0(inversion, ": \nArcana & Northern saxatilis")
  if(genetic.group == "ALL") title.text <- paste0(inversion, ": \nAll snails")
  
  ### Subset PCA data
  PC.dat <- PCA.data[[2]]
  
  ### 4.1: Plot PCA by heterozygosity
  pHet <- ggplot(PC.dat, aes(x=Axis1, y=Axis2)) + 
    geom_point(aes(colour = Heterozygosity), alpha = 0.8, size = 2.5) + 
    ggtitle(title.text)+
    coord_equal()+
    xlab(paste0("PC1 (", round(PCA.data[[1]][1], 2), "%)")) +
    ylab(paste0("PC2 (", round(PCA.data[[1]][2], 2), "%)")) +
    labs(colour = expression(H[E]))+
    scale_colour_gradient2(limits = c(0, 0.5), midpoint = 0.25, 
                           low = "#ffffd9", mid = "#41b6c4", high = "#081d58")+
    theme_bw()+
    theme(text=element_text(size=14))
  
  # Save the plot
  PATH <- "/Users/james/Documents/Inversion_detection/plots/PCA_by_Het_plots"
  ggsave(paste0(LG.s, "_", inv2, "_PCAbyHet_v2_",genetic.group,".tiff"), plot = pHet, "tiff", PATH,
         width = 15, height = 15, units = "cm")
  
  ### 4.2: Plot PCA by phyogeographic regions & species
  # Set colour palette & shape based on genetic group
  if(genetic.group == "NS"){ColPal = c("skyblue", "seagreen2", "salmon", "grey")
                      Shp = 16}
  if(genetic.group == "Sp"){ColPal = "gold3"; Shp = 16}
  if(genetic.group == "arc"){ColPal = c("skyblue", "seagreen2", "salmon"); Shp = 1}
  if(genetic.group == "arc_sax"){ColPal = c("skyblue", "seagreen2", "salmon", "grey")
                      Shp = c(1, 16)}
  if(genetic.group == "ALL"){ColPal = c("skyblue", "gold3", "seagreen2", "salmon", "grey")
                      Shp = c(1, 2, 16)}
  # Plot
  pGeo <- ggplot(PC.dat, aes(x=Axis1, y=Axis2)) + 
      geom_point(aes(colour = Region, pch = Species), alpha=0.6, size = 2.5) + 
      ggtitle(title.text)+
      coord_equal()+
      xlab(paste0("PC1 (", round(PCA.data[[1]][1], 2), "%)")) +
      ylab(paste0("PC2 (", round(PCA.data[[1]][2], 2), "%)")) +
      scale_colour_manual(values = ColPal)+
      scale_shape_manual(values = Shp)+
      theme_bw()+
      theme(text=element_text(size=14))
  
  # Save the plot
  PATH <- "/Users/james/Documents/Inversion_detection/plots/PCA_by_Country_plots"
  ggsave(paste0(LG.s, "_", inv2, "_PCAbyCountry_v2_",genetic.group,".tiff"), plot = pGeo, "tiff", PATH,
         width = 15, height = 15, units = "cm")
  
  ### 4.3: Plot PCA by ecotype & sex
  # Set colour palette based on genetic group
  if(genetic.group == "NS"){ColPal = c("#ef6e74ff", "#9192efff", "#3b3ce3ff", "grey", "#bf3a40ff")
                            Shp = c(1, 15, 17)}
  if(genetic.group == "Sp"){ColPal = c("#3b3ce3ff", "grey", "#bf3a40ff")
                            Shp = c(1, 15, 17)}
  if(genetic.group == "arc"){ColPal = "#fedb54ff"
                            Shp = c(15)}
  if(genetic.group == "arc_sax"){ColPal = c("#fedb54ff", "#ef6e74ff", "#9192efff", "#3b3ce3ff", "grey", "#bf3a40ff")
                                Shp = c(15)}
  if(genetic.group == "ALL"){ColPal = c("#fedb54ff", "#ef6e74ff", "#9192efff", 
                                        "#328380ff", "#3b3ce3ff", "grey", "#bf3a40ff")
                            Shp = c(1, 15, 17)}
  # Plot
  pEco <- ggplot(PC.dat, aes(x=Axis1, y=Axis2)) + 
    geom_point(aes(colour = Ecotype, pch = Sex), alpha = 0.8, size = 2.5) + 
    ggtitle(title.text)+
    coord_equal()+
    xlab(paste0("PC1 (", round(PCA.data[[1]][1], 2), "%)")) +
    ylab(paste0("PC2 (", round(PCA.data[[1]][2], 2), "%)")) +
    scale_shape_manual(values = c(1,15,17))+
    scale_colour_manual(values = ColPal)+
    theme_bw()+
    theme(text=element_text(size=14))
  
  # Save the plot
  PATH <- "/Users/james/Documents/Inversion_detection/plots/PCA_by_Ecotype_plots"
  ggsave(paste0(LG.s, "_", inv2, "_PCAbyEco_v2_",genetic.group,".tiff"), plot = pEco, "tiff", PATH,
         width = 15, height = 15, units = "cm")
}

#### 5. Run function across all inversions & the collinear genome in a LG ####

# Count number of inversions in LG
nInv <- length(unique(Invs$Inversion))

### Northern Littorina saxatilis group
PCA_NS <- lapply(1:(nInv+1), function(i){
  # Call collinear genome if 'i' goes past the number of inversion in current LG
  inv <- ifelse(i > nInv | nInv == 0, "collinear", unique(Invs$Inversion)[i])
  # Run PCA
  tmp <- PCA.genotypes(inv, genetic.group = "NS",
                       buffer = ifelse(inv == "collinear", yes = 2, no = 0))
  # Modify inversion names for saving file due to error from '/'
  inv2 <- ifelse(grepl("1/2", inv), gsub("/", "-", inv), inv)
  # Save output
  write.csv(tmp[[2]], paste0("PCA_per_inversion/", LG.s, "_PCA_of_", inv2, "_v2_NS.csv"), quote = FALSE, row.names = FALSE)
  plot.PCA(tmp, inv, genetic.group = "NS")
  return(tmp[[2]])
})

### Spanish Littorina saxatilis group
PCA_Sp <- lapply(1:(nInv+1), function(i){
  inv <- ifelse(i > nInv | nInv == 0, "collinear", unique(Invs$Inv)[i])
  tmp <- PCA.genotypes(inv, genetic.group = "Sp",
                       buffer = ifelse(inv == "collinear", yes = 2, no = 0))
  inv2 <- ifelse(grepl("1/2", inv), gsub("/", "-", inv), inv)
  write.csv(tmp[[2]], paste0("PCA_per_inversion/", LG.s, "_PCA_of_", inv2, "_v2_Sp.csv"), quote = FALSE, row.names = FALSE)
  plot.PCA(tmp, inv, genetic.group = "Sp")
  return(tmp[[2]])
})

### Littorina arcana group
PCA_arc <- lapply(1:(nInv+1), function(i){
  inv <- ifelse(i > nInv | nInv == 0, "collinear", unique(Invs$Inv)[i])
  tmp <- PCA.genotypes(inv, genetic.group = "arc",
                       buffer = ifelse(inv == "collinear", yes = 2, no = 0))
  inv2 <- ifelse(grepl("1/2", inv), gsub("/", "-", inv), inv)
  write.csv(tmp[[2]], paste0("PCA_per_inversion/", LG.s, "_PCA_of_", inv2, "_v2_arc.csv"), quote = FALSE, row.names = FALSE)
  plot.PCA(tmp, inv, genetic.group = "arc")
  return(tmp[[2]])
})

### Littorina arcana vs saxatilis
PCA_arc_sax <- lapply(1:(nInv+1), function(i){
  inv <- ifelse(i > nInv | nInv == 0, "collinear", unique(Invs$Inv)[i])
  tmp <- PCA.genotypes(inv, genetic.group = "arc_sax",
                       buffer = ifelse(inv == "collinear", yes = 2, no = 0))
  inv2 <- ifelse(grepl("1/2", inv), gsub("/", "-", inv), inv)
  write.csv(tmp[[2]], paste0("PCA_per_inversion/", LG.s, "_PCA_of_", inv2, "_v2_arc-sax.csv"), quote = FALSE, row.names = FALSE)
  plot.PCA(tmp, inv, genetic.group = "arc_sax")
  return(tmp[[2]])
})

### All snails
PCA_all <- lapply(1:(nInv+1), function(i){
  inv <- ifelse(i > nInv | nInv == 0, "collinear", unique(Invs$Inv)[i])
  tmp <- PCA.genotypes(inv, genetic.group = "ALL",
                       buffer = ifelse(inv == "collinear", yes = 2, no = 0))
  inv2 <- ifelse(grepl("1/2", inv), gsub("/", "-", inv), inv)
  write.csv(tmp[[2]], paste0("PCA_per_inversion/", LG.s, "_PCA_of_", inv2, "_v2_all.csv"), quote = FALSE, row.names = FALSE)
  plot.PCA(tmp, inv, genetic.group = "ALL")
  return(tmp[[2]])
})
