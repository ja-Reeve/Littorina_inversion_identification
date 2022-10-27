### PCA for heterozygote clusters ###
### This script was modified version designed by Katie Hearn from the University of Shefield.
### It converts a genotype file formates as: CHROM, POS, SAMPLES; into a PCA looking for 
### clusters that match heterozygote clusters. This can be targeted to test if regions of
### the genome with elevated heterozygosity are inversions.
### James Reeve - University of Gothenburg
### 30/11/2020
### Edited 02/09/2021 to adjust to new consensus map
### Edit 17/12/2021: version 2 created with new name. Many data handling steps are moved to
### functions to make the code a bit more efficient to run.
### Edit 14/01/2022: Clean up functions and recolour plots.
### Edit 26/08/2022: Remove compressa from the analysis and change PCA plot colouring
### Edit 08/09/2022: alter inversion boundaries for each cohort
### Edit 07/10/2022: add compressa back in, but project them onto the PCA of other samples

##### Note: this script takes 1-2h to run per linkage group! ########

#### 1. Preparation ####
### Clean-up environment
rm(list = ls())
dev.off()
setwd("")
options(stringsAsFactors = FALSE)

Start_time <- Sys.time()

### Packages
library("tidyverse")
library("ggpubr")
library("data.table")
library("adegenet")

### Set linkage group
LG.n <- 1
# Change to whichever LG you want to analyse
LG.s <- paste0("LG", LG.n)

### Linkage map
LG <- read.table("ConsensusMap.v2.txt", header = T)
LG <- LG[LG$LG == LG.s,]

### Find inversions on current LG
Invs <- read.csv("Inversion_breakpoint_estimates.csv", header = TRUE)
Invs <- Invs[Invs$LG == LG.s, c("Inversion", "start_pos", "end_pos", "cohort")]

### Access collection details
sample_info <- read.csv("/PATH/Sample_information.csv", sep=";", check.names = FALSE)
# Remove 'IMI_6_2' which is missing sample details
sample_info <- sample_info[sample_info$Sample_ID != "IMI_6_2",]
# Replace NAs in sex for 'Un'
sample_info[is.na(sample_info$sex), "sex"] <- "Un"
# Change "saxatilis?" to "saxatilis
sample_info[sample_info$Species %in% c("saxatilis?", "saxailis?"), "Species"] <- "saxatilis"
# Changing some ecotype labels to reduce number of ecotypes
sample_info[sample_info$Ecotype %in% c("Crab", "crabish"), "Ecotype"] <- "crab" # Rename crab ecotype
sample_info[sample_info$Ecotype == "tennebrosa", "Ecotype"] <- "brackish" # Rename tennebrosa brackish
sample_info[sample_info$Ecotype == "Wave", "Ecotype"] <- "wave" # Rename wave ecotype
sample_info[sample_info$Ecotype == "barnacle(ish)", "Ecotype"] <- "barnacle" # Rename barnacle ecotype
sample_info[sample_info$Ecotype %in% c("", "?", "midshore (no ecotypes)"), "Ecotype"] <- "other" # Rename everything else
# "arcana" and "compressa" are added as ecotypes
sample_info[sample_info$Species == "arcana" & 
              sample_info$Sample_ID != "CEA_Larc_F_1", "Ecotype"] <- "arcana" # rename to arcana
sample_info[sample_info$Species == "compressa" | 
              sample_info$Sample_ID == "CEA_Larc_F_1", "Ecotype"] <- "compressa" # rename to compressa
# Create new column for phylogenetic regions
sample_info$Region <- NA
sample_info[sample_info$Country %in% c("England", "Sweden") & sample_info$location != "Dersingham" |
              sample_info$location == "S Abbs", "Region"] <- "North_Sea"
sample_info[sample_info$Country %in% c("Norway", "Russia"), "Region"] <- "Northern_Atlantic"
sample_info[sample_info$Country %in% c("Wales", "Ireland", "Isle of man", "France", "Iceland") |
                        sample_info$location %in% c("Oban", "Dersingham"), "Region"] <- "Celtic_Sea"
sample_info[sample_info$Country == "Spain", "Region"] <- "Iberia"
sample_info[sample_info$Country == "USA", "Region"] <- "Western_Atlantic"

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

#### 2. Data wrangling ####
### 2.1: Filter out missing samples
# Filter sample info
sample_info <- sample_info[sample_info$Sample_ID %in% colnames(GT)[-c(1:2)],]
# Filter genotypes
GT <- GT[,c("CHROM", "POS", sample_info$Sample_ID)]

### 2.2: Create vectors for each cohort
### Cohorts based on the phylogenetic analysis of Sean Stankowski
###   Sp = Spanish Littorina saxatilis
###   sax = Littorina saxatilis from all other sites (i.e. Northern cohort)
### Sister species
###   arc = Littorina arcana
###   comp = Littorina compressa
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

# Remove samples not collected from arcana sites & remove compressa
arc_sax_IDs <- sample_info[sample_info$location %in% unique(sample_info[sample_info$Species == "arcana", "location"]) &
                           sample_info$Species != "compressa" &
                           sample_info$Sample_ID != "CEA_Larc_F_1", "Sample_ID"]

### Remove compressa from the analysis
#GT <- GT[,!(colnames(GT) %in% comp_IDs)]

#### 3. Function to run PCA for inversion for each cohort ####
### This is a complex process with several computationally intensive steps.
### It's more efficient to run this within a function so temporary files aren't
### clogging up the RAM. This function is comprised of several sub-steps:
###       3.1: Filtering genotypes - removing any SNPs which are fixed in a cohort
###             or SNPs that only appear once.
###       3.2: Transposing the genotype table - needed for the PCA
###       3.3: Filter SNPs to focal region - this is two functions A) is for inversions
###            B) is for the collinear genome.
###       3.4: PCA - run the PCA analysis
### I've set it up in a way that it should only need two input flags, and one optional flags:
###       inversion = name of inversion as a string. Can also be "collinear" to select collinear region.
###       cohort = string of 4 values {"sax", "Sp", "arc", "ALL"}
###       buffer = buffer of XcM around inversion breakpoints

PCA.genotypes <- function(inversion, cohort, buffer = 0){
  ### Error messages
  if(class(cohort) != "character") stop(print(paste0(cohort, " must be a string of class 'character'.")))
  if(!(cohort %in% c("sax", "arc", "Sp", "ALL", "arc_sax"))) stop(print(" cohort entred incorrectly, must be 'sax', 'arc', 'Sp', or 'ALL'!"))
  if(class(inversion) != "character") stop(print(" inversion must be a string value."))
  if(class(buffer) != "numeric") stop(" buffer.in must be a numeric value.")
  
  ### 3.1: Filter genotypes to remove alleles fixed in a cohort and which are only a single SNP
  # Ignore step if cohort == "ALL"
  if(cohort == "ALL"){GT_filt <- GT} else {
    # Find vector with cohort's IDs
    IDs <- get(paste(cohort, "IDs", sep = "_"))
    # Subset genotypes to cohort
    GT_filt <- GT %>% select("CHROM", "POS", all_of(IDs))
    # Remove fixed sites and SNPs with MAF ≤ 1
    GT_filt <- GT_filt[apply(GT_filt, 1, function(i){
      sum(i[-c(1:2)] %in% c("./.", "0/0")) != length(i) - 2 | 
        sum(i[-c(1:2)] %in% c("./.", "1/1")) != length(i) - 2 |
        sum(i[-c(1:2)] %in% c("./.", "0/1")) != length(i) - 2}), ]
  }
  
  
  ### 3.2: Transpose the genotype table
  # Combine CHROM and POS into "SNP_ID"
  tmp <- GT_filt %>% unite(SNP_ID, c("CHROM", "POS"), sep=":")
  # Transpose using function from data.table package
  GT_trans <- transpose(tmp, keep.names = "SnailID", make.names = "SNP_ID")
  
  rm(tmp)
  
  ### 3.3: Filter to target region
  # Find contigs in target region
  if(inversion == "collinear"){
    
    ### Target region == collinear genome
    # List of all map positions
    map_pos <- unique(LG$avgMP)
    # Subset 'Invs' by cohort
    IBp <- Invs[Invs$cohort == ifelse(cohort %in% c("arc_sax", "ALL"), "sax", cohort), ]
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
                   Invs$cohort == ifelse(cohort %in% c("arc_sax", "ALL"), "sax", cohort), 
                 "start_pos"] # Start of inversion
    Epos <- Invs[Invs$Inv == inversion & 
                   Invs$cohort == ifelse(cohort %in% c("arc_sax", "ALL"), "sax", cohort), 
                 "end_pos"] # End of inversion
    # Filter linkage map to inversion area
    contigs <- LG[LG$avgMP >= Spos + buffer & LG$avgMP <= Epos - buffer, "contig"]
    
  }
  
  # Subset transformed genotypes to those in target region
  GT_trans_filt <- GT_trans[ ,which(gsub(":.*", "", colnames(GT_trans)) %in% contigs)]
  
  
  ### 3.4: Run the PCA
  # Error check if contigs exist - e.g. some inversions cover the whole LG leaving no collinear genome
  if(ncol(GT_trans_filt) == 0) return(NA)
  # Using functions from the adgenet package
  # Taken with minor modification from Katie Hearn
  # A: make genind object
  genind <- df2genind(GT_trans_filt, ploidy=2, sep="/", NA.char = "_")
  # B: rescore missing genotypes as the mean values
  genind <- scaleGen(genind, NA.method="mean", scale = FALSE)
  rownames(genind) <- colnames(GT_filt)[-c(1:2)]
  # C: PCA call
  # For the full data set, compressa are projected onto the PCA to avoid clustering artifact with Iberian snails
  if(cohort == "ALL"){
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
  # A: Get genotypes
  tmp <- GT[,colnames(GT) %in% PC.dat[[2]]$Sample_ID]
  # B: Calculate heterozygosity per individual
  tmp2 <- apply(tmp, 2, function(X){X[X != "./."]}) # Remove missing SNPs
  Het <- lapply(tmp2, function(X){
    if(length(X) > 0){
      sum(X == "0/1") / length(X) # Proportion of heterozygotes
    }})
  
  rm(tmp, tmp2)
  # C: Add heterozygosity scores to PCA data
  PC.dat[[2]]$Heterozygosity <- as.numeric(Het)
  
  
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
plot.PCA <- function(PCA.data, inversion, cohort){#, x.range = NULL, y.range = NULL){
  ### Skip function if PCA didn't work
  if(class(PCA.data) != "list") stop(paste0("PCA did not work for ", inversion, "_", cohort))
  
  ### Modify inversion names for saving plots due to error from '/'
  inv2 <- ifelse(grepl("1/2", inversion), gsub("/", "-", inversion), inversion)
  
  ### Create plot titles
  if(cohort == "sax") title.text <- paste0(inversion, ": \nNorthern Littorina saxatilis")
  if(cohort == "Sp") title.text <- paste0(inversion, ": \nSpanish Littorina saxatilis")
  if(cohort == "arc") title.text <- paste0(inversion, ": \nLittorina arcana")
  if(cohort == "arc_sax") title.text <- paste0(inversion, ": \nArcana & Northern saxatilis")
  if(cohort == "ALL") title.text <- paste0(inversion, ": \nAll snails")
  
  ### Subset PCA data
  PC.dat <- PCA.data[[2]]
  
  ### 4.1: Plot PCA by heterozygosity
  pHet <- ggplot(PC.dat, aes(x=Axis1, y=Axis2)) + 
    geom_point(aes(colour = Heterozygosity), alpha=0.4, size = 2.5) + 
    ggtitle(title.text)+
    coord_equal()+
    xlab(paste0("PC1 (", round(PCA.data[[1]][1], 2), "%)")) +
    ylab(paste0("PC2 (", round(PCA.data[[1]][2], 2), "%)")) +
    labs(colour = expression(H[E]))+
    scale_colour_gradient2(limits = c(0.1, 0.46), midpoint = 0.3, 
                           low = "blue", mid = "orange2", high = "firebrick")+
    theme_bw()+
    theme(text=element_text(size=14))
  
  # Save the plot
  PATH <- "/Users/james/Documents/Inversion_detection/plots/PCA_by_Het_plots"
  ggsave(paste0(LG.s, "_", inv2, "_PCAbyHet_",cohort,".tiff"), plot = pHet, "tiff", PATH,
         width = 15, height = 15, units = "cm")
  
  ### 4.2: Plot PCA by phyogeographic regions & species
  # Set colour palette & shape based on cohort
  if(cohort == "sax"){ColPal = c("skyblue", "seagreen2", "salmon", "grey")
                      Shp = 16}
  if(cohort == "Sp"){ColPal = "gold3"; Shp = 16}
  if(cohort == "arc"){ColPal = c("skyblue", "seagreen2", "salmon"); Shp = 1}
  if(cohort == "arc_sax"){ColPal = c("skyblue", "seagreen2", "salmon", "grey")
                      Shp = c(1, 16)}
  if(cohort == "ALL"){ColPal = c("skyblue", "gold3", "seagreen2", "salmon", "grey")
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
  ggsave(paste0(LG.s, "_", inv2, "_PCAbyCountry_",cohort,".tiff"), plot = pGeo, "tiff", PATH,
         width = 15, height = 15, units = "cm")
  
  ### 4.3: Plot PCA by ecotype & sex
  # Set colour palette based on cohort
  if(cohort == "sax"){ColPal = c("#ef6e74ff", "#3b3ce3ff", "grey", "#9192efff", "#bf3a40ff")}
  if(cohort == "Sp"){ColPal = c("#3b3ce3ff", "grey", "#bf3a40ff")}
  if(cohort == "arc"){ColPal = "#fedb54ff"}
  if(cohort == "arc_sax"){ColPal = c("#fedb54ff", "#ef6e74ff", "#3b3ce3ff", 
                                     "grey", "#9192efff", "#bf3a40ff")}
  if(cohort == "ALL"){ColPal = c("#fedb54ff", "#ef6e74ff", "#328380ff", "#3b3ce3ff", 
                                 "grey", "#9192efff", "#bf3a40ff")}
  # Plot
  pEco <- ggplot(PC.dat, aes(x=Axis1, y=Axis2)) + 
    geom_point(aes(colour = Ecotype, pch = sex), alpha=0.6, size = 2.5) + 
    ggtitle(title.text)+
    coord_equal()+
    xlab(paste0("PC1 (", round(PCA.data[[1]][1], 2), "%)")) +
    ylab(paste0("PC2 (", round(PCA.data[[1]][2], 2), "%)")) +
    scale_shape_manual(values = c(15,17,1))+
    scale_colour_manual(values = ColPal)+
    theme_bw()+
    theme(text=element_text(size=14))
  
  # Save the plot
  PATH <- "/Users/james/Documents/Inversion_detection/plots/PCA_by_Ecotype_plots"
  ggsave(paste0(LG.s, "_", inv2, "_PCAbyEco_",cohort,".tiff"), plot = pEco, "tiff", PATH,
         width = 15, height = 15, units = "cm")
}

#### 5. Run function across all inversions & the collinear genome in a LG ####

# Count number of inversions in LG
nInv <- length(unique(Invs$Inversion))

### Northern Littorina saxatilis cohort
PCA_sax <- lapply(1:(nInv+1), function(i){
  # Call collinear genome if 'i' goes past the number of inversion in current LG
  inv <- ifelse(i > nInv | nInv == 0, "collinear", unique(Invs$Inversion)[i])
  # Run PCA
  tmp <- PCA.genotypes(inv, cohort = "sax",
                       buffer = ifelse(inv == "collinear", yes = 2, no = 0))
  # Modify inversion names for saving file due to error from '/'
  inv2 <- ifelse(grepl("1/2", inv), gsub("/", "-", inv), inv)
  # Save output
  write.csv(tmp[[2]], paste0("PCA_per_inversion/", LG.s, "_PCA_of_", inv2, "_sax.csv"), quote = FALSE)
  plot.PCA(tmp, inv, cohort = "sax")
  return(tmp[[2]])
})

### Spanish Littorina saxatilis cohort
PCA_Sp <- lapply(1:(nInv+1), function(i){
  inv <- ifelse(i > nInv | nInv == 0, "collinear", unique(Invs$Inv)[i])
  tmp <- PCA.genotypes(inv, cohort = "Sp",
                       buffer = ifelse(inv == "collinear", yes = 2, no = 0))
  inv2 <- ifelse(grepl("1/2", inv), gsub("/", "-", inv), inv)
  write.csv(tmp[[2]], paste0("PCA_per_inversion/", LG.s, "_PCA_of_", inv2, "_Sp.csv"), quote = FALSE)
  plot.PCA(tmp, inv, cohort = "Sp")
  return(tmp[[2]])
})

### Littorina arcana cohort
PCA_arc <- lapply(1:(nInv+1), function(i){
  inv <- ifelse(i > nInv | nInv == 0, "collinear", unique(Invs$Inv)[i])
  tmp <- PCA.genotypes(inv, cohort = "arc",
                       buffer = ifelse(inv == "collinear", yes = 2, no = 0))
  inv2 <- ifelse(grepl("1/2", inv), gsub("/", "-", inv), inv)
  write.csv(tmp[[2]], paste0("PCA_per_inversion/", LG.s, "_PCA_of_", inv2, "_arc.csv"), quote = FALSE)
  plot.PCA(tmp, inv, cohort = "arc")
  return(tmp[[2]])
})

### Littorina arcana vs saxatilis cohort
PCA_arc_sax <- lapply(1:(nInv+1), function(i){
  inv <- ifelse(i > nInv | nInv == 0, "collinear", unique(Invs$Inv)[i])
  tmp <- PCA.genotypes(inv, cohort = "arc_sax",
                       buffer = ifelse(inv == "collinear", yes = 2, no = 0))
  inv2 <- ifelse(grepl("1/2", inv), gsub("/", "-", inv), inv)
  write.csv(tmp[[2]], paste0("PCA_per_inversion/", LG.s, "_PCA_of_", inv2, "_arc-sax.csv"), quote = FALSE, row.names = FALSE)
  plot.PCA(tmp, inv, cohort = "arc_sax")
  return(tmp[[2]])
})

### All snails
PCA_all <- lapply(1:(nInv+1), function(i){
  inv <- ifelse(i > nInv | nInv == 0, "collinear", unique(Invs$Inv)[i])
  tmp <- PCA.genotypes(inv, cohort = "ALL",
                       buffer = ifelse(inv == "collinear", yes = 2, no = 0))
  inv2 <- ifelse(grepl("1/2", inv), gsub("/", "-", inv), inv)
  write.csv(tmp[[2]], paste0("PCA_per_inversion/", LG.s, "_PCA_of_", inv2, "_all.csv"), quote = FALSE)
  plot.PCA(tmp, inv, cohort = "ALL")
  return(tmp[[2]])
})


End_time <- Sys.time()
End_time - Start_time
