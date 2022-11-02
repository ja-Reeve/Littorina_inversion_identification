####### Ecotype contrasts #########
### Contrast the inversion frequencies of different ecotypes.
### Nested GLM, with 1) an overall test for differences
### followed by 2) a test of specific hypotheses
###    H1: crab ≠ wave
###    H2: crab = brackish
###    H3: wave = barnacle
### Finally I will run a separate test to identify significant
### differences between Lsax & Larc for sites were both species
### were sampled
### James Reeve - University of Gothenburg
### 09/02/2022
### Edit 28/03/2022 - now using the arc-sax data set
### Edit 21/10/2022 - Change G-test to Possion-GLM which
### controls for the site effect.

rm(list = ls())
dev.off()
setwd("")

### Packages
library("tidyverse")
library("ggpubr")

### Parameters
FDR.thresh <- 0.01


#### A: Access inversion genotype data ####

### List of inversions
INVs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC5.1", "LGC6.1-2", "LGC7.1", 
          "LGC7.2", "LGC9.1", "LGC9.2", "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", 
          "LGC12.2", "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.2", "LGC14.3", "LGC17.1")

### Download Northern saxatilis data
Nsax <- lapply(INVs, function(inv){
  LG <- gsub("C", "", gsub("[.].*", "", inv))
  dat <- read.csv(paste0(LG, "_PCA_of_", inv, "_sax.csv"), header = TRUE)
  # Renaming columns to remove inversion prefix
  colnames(dat) <- gsub(paste0(inv,"."), "", colnames(dat))
  return(dat)
})
# Add inversion as listnames
names(Nsax) <- INVs

### Download full data
ArcSax <- lapply(INVs, function(inv){
  LG <- gsub("C", "", gsub("[.].*", "", inv))
  dat <- read.csv(paste0(LG, "_PCA_of_", inv, "_arc-sax.csv"), header = TRUE)
  # Renaming columns to remove inversion prefix
  colnames(dat) <- gsub(paste0(inv,"."), "", colnames(dat))
  return(dat)
})
# Add inversion as listnames
names(ArcSax) <- INVs


#### B: Filter genotype data ####
### Remove any inversions where K ≠ {3,6}
### In these cases it is hard to say which cluster is or isn't heterokaryotypic
Nsax2 <- Nsax[sapply(Nsax, 
                     function(X){length(unique(X$genotype))}) %in% c(3,6)]
ArcSax2 <- ArcSax[sapply(ArcSax, 
                         function(X){length(unique(X$genotype))}) %in% c(3,6)]

### Also remove "LGC14.1-2", due to uncertain clustering
Nsax2 <- Nsax2[which(names(Nsax2) != "LGC14.1-2")]
ArcSax2 <- ArcSax2[which(names(ArcSax2) != "LGC14.1-2")]

### Additional filters in ArcSax when arcana forms a separate cluster to Lsax
ArcSax2 <- ArcSax2[which(!(names(ArcSax2) %in% c("LGC5.1", "LGC12.2", "LGC12.3", "LGC14.1")))]


#### C: Functions for each ecotype contrast ####
Eco_contrast <- function(PCA.data, inversion, ecotypes){
  
  ### Calculate number of ecotypes
  Neco <- length(ecotypes)
  
  ### Error messages
  if(class(inversion) != "character") stop(print(paste0(inversion, " must be an object of class 'character'.")))
  if(class(ecotypes) != "character") stop(print("'ecotypes' must be a character vector."))
  if(Neco < 2) stop(print("'ecotypes' must have a length ≥ 2."))
  if(class(PCA.data) != "list") stop(print("'PCA.data' must be a list of PCA data results. 1 inversion per list element."))
  if(class(PCA.data[[1]]) != "data.frame") stop(print("Each element of 'PCA.data' must be a data.frame with samples as rows and eigenvectors as columns.\n
                                                      This dataset MUST also have an 'Ecotype' and 'genotype' column."))
  
  ### Filter to specified ecotypes
  tmp <- PCA.data[[inversion]]
  tmp2 <- tmp[tmp$Ecotype %in% ecotypes,]
  
  ### Rescale arragements as alleles counts
  if(length(unique(tmp$genotype)) == 3){
    # Rescale the alternative allele count
    tmp2$N <- NA
    tmp2[tmp2$genotype == "RR", "N"] <- 0
    tmp2[tmp2$genotype == "RA", "N"] <- 1
    tmp2[tmp2$genotype == "AA", "N"] <- 2
    
    # Count arrangements
    tmp3 <- tmp2 %>% group_by(location, Ecotype) %>% 
      summarise("Arrangement" = c("A", "R"), "N" = c(sum(N), 2*n() - sum(N)))
  
    } else {
      
    # Count of A1 allele
    tmp2$NA1 <- NA
    tmp2[tmp2$genotype == "RA1", "NA1"] <- 1
    tmp2[tmp2$genotype == "A1A2", "NA1"] <- 1
    tmp2[tmp2$genotype == "A1A1", "NA1"] <- 2
    tmp2[is.na(tmp2$NA1), "NA1"] <- 0
    
    # Count of A2 allele
    tmp2$NA2 <- NA
    tmp2[tmp2$genotype == "RA2", "NA2"] <- 1
    tmp2[tmp2$genotype == "A1A2", "NA2"] <- 1
    tmp2[tmp2$genotype == "A2A2", "NA2"] <- 2
    tmp2[is.na(tmp2$NA2), "NA2"] <- 0
    
    # Count arrangements
    tmp3 <- tmp2 %>% group_by(location, Ecotype) %>% 
      summarise("Arrangement" = c("A", "A2", "R"), 
                "N" = c(sum(NA1), sum(NA2), 2*n() - (sum(NA1) + sum(NA2))))
  }
  
  ### Create glm  model
  try({
    mod <- glm(N ~ Ecotype + location + Arrangement + Arrangement:location + Arrangement:Ecotype, 
                              family = poisson(link = "log"), data = tmp3)
    
  ### Summarise results in a deviance table
  DT <- anova(mod)
  DT2 <- data.frame("predictor" = rownames(DT), 
                    "deviance" = DT$Deviance,
                    "df" = DT$Df,
                    "p-value" = pchisq(DT$Deviance, DT$Df, lower.tail = FALSE),
                    "residual_deviance" = DT$`Resid. Dev`,
                    "residual_df" = DT$`Resid. Df`)
    
  # Save deviance table
  if(identical(ecotypes, c("crab", "wave", "barnacle", "tenebrosa"))) Test <- "T1"
  if(identical(ecotypes, c("crab", "wave"))) Test <- "T2a"
  if(identical(ecotypes, c("crab", "tenebrosa"))) Test <- "T2b"
  if(identical(ecotypes, c("wave", "barnacle"))) Test <- "T2c"
  if(identical(ecotypes, c("saxatilis", "arcana"))) Test <- "T3"
  if(identical(ecotypes, c("crab", "arcana"))) Test <- "T4a"
  if(identical(ecotypes, c("wave", "arcana"))) Test <- "T4b"
    
  write.csv(DT2, paste0(PATH, "Ecotype_contrasts/", Test, ".deviance_table_", inversion, ".csv"),
              row.names = FALSE, quote = FALSE)
    
  ### Significance of the effect of ecotype on arrangement
  Pval <- pchisq(DT2$deviance[6], df = DT2$df[6], lower.tail = FALSE)
    
  ### Save output as dataframe
  res <- tmp3
  res$Inv <- inversion
  res$chi_sq <- DT2$deviance[6]
  res$df <- DT2$df[6]
  res$P.value <- Pval
  # Re-order columns
    
  res <- res[,c("Inv", "location", "Ecotype", "Arrangement", "N", "chi_sq", "df", "P.value")]
    
  return(res)
  })
  
}


#### D: Run function ####

### Note: using Northern saxatilis data for now! ###

### 1. Overall difference among ecotypes
T1 <- lapply(INVs[INVs %in% names(Nsax2)], Eco_contrast, PCA.data = Nsax2,
             ecotypes = c("crab", "wave", "barnacle", "tenebrosa"))
T1 <- do.call(rbind.data.frame, T1)

# Calculate FDR
tmp <- T1 %>% group_by(Inv) %>% summarise(P.value) %>% unique()
tmp$adj.P.value <- p.adjust(tmp$P.value, method = "fdr")
T1 <- T1 %>% left_join(tmp)
rm(tmp)

# Plot
ggplot(T1, aes(y = factor(Inv, levels = INVs)))+
  geom_point(aes(x = Arrangement, size = N, colour = adj.P.value < FDR.thresh))+
  labs(x = "Inversion arrangement", y = "Inversion", 
       title = "Crab-wave test", size = "Freq.")+
  facet_grid(cols = vars(Ecotype))+
  scale_colour_manual(values = c("grey", "brown"), guide = "none")+
  theme_bw()

### List significant inversions
SigInvs <- unique(T1[T1$adj.P.value < FDR.thresh, "Inv"])$Inv

### 2a. Crab - wave test
### Alt H: crab is different from wave
T2a <- lapply(SigInvs, Eco_contrast, PCA.data = Nsax2,
              ecotypes = c("crab", "wave"))
T2a <- do.call(rbind.data.frame, T2a)

# Calculate FDR
tmp <- T2a %>% group_by(Inv) %>% summarise(P.value) %>% unique()
tmp$adj.P.value <- p.adjust(tmp$P.value, method = "fdr")
T2a <- T2a %>% left_join(tmp)
rm(tmp)

# Plot
ggplot(T2a, aes(y = factor(Inv, levels = SigInvs)))+
  geom_point(aes(x = Arrangement, size = N, colour = adj.P.value < FDR.thresh))+
  labs(x = "Inversion arrangement", y = "Inversion", 
       title = "Crab-wave test", size = "Freq.")+
  facet_grid(cols = vars(Ecotype))+
  scale_colour_manual(values = c("grey", "brown"), guide = "none")+
  theme_bw()

### 2b. Crab - brackish test
### Alt H: crab is different from wave
T2b <- lapply(SigInvs, Eco_contrast, PCA.data = Nsax2,
              ecotypes = c("crab", "tenebrosa"))
T2b <- do.call(rbind.data.frame, T2b)

# Calculate FDR
tmp <- T2b %>% group_by(Inv) %>% summarise(P.value) %>% unique()
tmp$adj.P.value <- p.adjust(tmp$P.value, method = "fdr")
T2b <- T2b %>% left_join(tmp)
rm(tmp)

# Plot
ggplot(T2b, aes(y = factor(Inv, levels = SigInvs)))+
  geom_point(aes(x = Arrangement, size = N, colour = adj.P.value > FDR.thresh))+
  labs(x = "Inversion arrangement", y = "Inversion", 
       title = "Crab-brackish test", size = "Freq.")+
  facet_grid(cols = vars(Ecotype))+
  scale_colour_manual(values = c("grey", "brown"), guide = "none")+
  theme_bw()

### 2c. Wave - barnacle test
### Alt H: crab is different from wave
T2c <- lapply(SigInvs, Eco_contrast, PCA.data = Nsax2,
              ecotypes = c("wave", "barnacle"))
T2c <- do.call(rbind.data.frame, T2c)

# Calculate FDR
tmp <- T2c %>% group_by(Inv) %>% summarise(P.value) %>% unique()
tmp$adj.P.value <- p.adjust(tmp$P.value, method = "fdr")
T2c <- T2c %>% left_join(tmp)
rm(tmp)

# Plot
ggplot(T2c, aes(y = factor(Inv, levels = SigInvs)))+
  geom_point(aes(x = Arrangement, size = N, colour = adj.P.value > FDR.thresh))+
  labs(x = "Inversion arrangement", y = "Inversion", 
       title = "Wave-barnacle test", size = "Freq.")+
  facet_grid(cols = vars(Ecotype))+
  scale_colour_manual(values = c("grey", "brown"), guide = "none")+
  theme_bw()

#### E: Run function for Arcana - saxatilis tests ####

### Using PCA of arc-sax snails ###
TMP <- lapply(ArcSax2, function(X){
  tmp2 <- X$Ecotype
  tmp2[tmp2 != "arcana"] <- "saxatilis"
  X$Ecotype <- tmp2
  return(X)
})
T3 <- lapply(names(TMP), Eco_contrast, PCA.data = TMP,
             ecotypes = c("saxatilis", "arcana"))
T3 <- do.call(rbind.data.frame, T3)

# Calculate FDR
tmp <- T3 %>% group_by(Inv) %>% summarise(P.value) %>% unique()
tmp$adj.P.value <- p.adjust(tmp$P.value, method = "fdr")
T3 <- T3 %>% left_join(tmp)
rm(tmp)

# Plot
ggplot(T3, aes(y = factor(Inv, levels = INVs)))+
  geom_point(aes(x = Arrangement, size = N, colour = adj.P.value < FDR.thresh))+
  labs(x = "Inversion arrangement", y = "Inversion", 
       title = "saxatilis-arcana test", size = "Freq.")+
  facet_grid(cols = vars(Ecotype))+
  scale_colour_manual(values = c("grey", "brown"), guide = "none")+
  theme_bw()

### List significant inversions
AS_SigInvs <- unique(T3[T3$adj.P.value < FDR.thresh, "Inv"])$Inv

### 4a. Crab - arcana test
### Alt H: crab is different to L. arcana
T4a <- lapply(AS_SigInvs, Eco_contrast, PCA.data = ArcSax2,
              ecotypes = c("crab", "arcana"))
T4a <- do.call(rbind.data.frame, T4a)

# Calculate FDR
tmp <- T4a %>% group_by(Inv) %>% summarise(P.value) %>% unique()
tmp$adj.P.value <- p.adjust(tmp$P.value, method = "fdr")
T4a <- T4a%>% left_join(tmp)
rm(tmp)

# Plot
ggplot(T4a, aes(y = factor(Inv, levels = INVs)))+
  geom_point(aes(x = Arrangement, size = N, colour = adj.P.value < FDR.thresh))+
  labs(x = "Inversion arrangement", y = "Inversion", 
       title = "Crab-arcana test", size = "Freq.")+
  facet_grid(cols = vars(Ecotype))+
  scale_colour_manual(values = c("grey", "brown"), guide = "none")+
  theme_bw()

### 3b. Wave - arcana test
### Alt H: Wave is the same as L. arcana
T4b <- lapply(AS_SigInvs, Eco_contrast, PCA.data = ArcSax2,
              ecotypes = c("wave", "arcana"))
T4b <- do.call(rbind.data.frame, T4b)

# Calculate FDR
tmp <- T4b %>% group_by(Inv) %>% summarise(P.value) %>% unique()
tmp$adj.P.value <- p.adjust(tmp$P.value, method = "fdr")
T4b <- T4b%>% left_join(tmp)
rm(tmp)

# Plot
ggplot(T4b, aes(y = factor(Inv, levels = INVs)))+
  geom_point(aes(x = Arrangement, size = N, colour = adj.P.value > FDR.thresh))+
  labs(x = "Inversion arrangement", y = "Inversion", 
       title = "Wave-arcana test", size = "Freq.")+
  facet_grid(cols = vars(Ecotype))+
  scale_colour_manual(values = c("grey", "brown"), guide = "none")+
  theme_bw()

#### F: Save results ####
PATH <- "/Users/james/Documents/Inversion_detection/"

write.csv(T1, paste0(PATH, "Ecotype_contrasts/T1.Saxatilis_ecotype_contrast.v3.csv"), quote = FALSE, row.names = FALSE)
write.csv(T2a, paste0(PATH, "Ecotype_contrasts/T2a.crab_wave_contrast.v3.csv"), quote = FALSE, row.names = FALSE)
T2b <- T2b[grep("Error", T2b$Inv, invert = TRUE),] # Removes error message that got caught in the data
write.csv(T2b, paste0(PATH, "Ecotype_contrasts/T2b.crab_brackish_contrast.v3.csv"), quote = FALSE, row.names = FALSE)
write.csv(T2c, paste0(PATH, "Ecotype_contrasts/T2c.wave-barnacle_contrast.v3.csv"), quote = FALSE, row.names = FALSE)
write.csv(T3, paste0(PATH, "Ecotype_contrasts/T3.sax_arc_contrast.v3.csv"), quote = FALSE, row.names = FALSE)
write.csv(T4a, paste0(PATH, "Ecotype_contrasts/T4a.crab_arc_contrast.v3.csv"), quote = FALSE, row.names = FALSE)
write.csv(T4b, paste0(PATH, "Ecotype_contrasts/T4b.wave_arc_contrast.v3.csv"), quote = FALSE, row.names = FALSE)
