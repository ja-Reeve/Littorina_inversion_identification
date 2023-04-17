####### Ecotype contrasts #########
### Contrast the inversion frequencies of different ecotypes.
### Four GLMs testing hypothetical ecotype differences
###    H1: crab ≠ wave
###    H2: wave = barnacle
###    H3: arcana = wave
###    H4: arcana ≠ crab
### Before running each test, inversion frequencies will be
### filtered to only retain sties where both ecotypes were
### samples. Finally I will run a separate test to identify 
### significant differences between Lsax & Larc.
### James Reeve - University of Gothenburg
### 13/04/2023
### Note: a crab = brackish test was attempted, but failed
### due to limited and geographically isolated sampling of 
### brackish samples

rm(list = ls())
dev.off()
setwd("/Users/james/Documents/Inversion_detection/PCA_per_inversion")

### Packages
library("tidyverse")
library("ggpubr")
library("caret")

### Parameters
FDR.thresh <- 0.05
PATH <- "/Users/james/Documents/Inversion_detection/"


#### A: Access inversion genotype data ####

### List of inversions
INVs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC5.1", "LGC6.1-2", "LGC7.1", 
          "LGC7.2", "LGC9.1", "LGC9.2", "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", 
          "LGC12.2", "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.2", "LGC14.3", "LGC17.1")

### Download Northern saxatilis data
Nsax <- lapply(INVs, function(inv){
  LG <- gsub("C", "", gsub("[.].*", "", inv))
  dat <- read.csv(paste0(LG, "_PCA_of_", inv, "_v2_NS.csv"), header = TRUE)
  # Renaming columns to remove inversion prefix
  colnames(dat) <- gsub(paste0(inv,"."), "", colnames(dat))
  return(dat)
})
# Add inversion as listnames
names(Nsax) <- INVs

############ EDIT: remove LGC9.2 due to rarity of arrangements ###############

### Download full data
ArcSax <- lapply(INVs, function(inv){
  LG <- gsub("C", "", gsub("[.].*", "", inv))
  dat <- read.csv(paste0(LG, "_PCA_of_", inv, "_v2_arc-sax.csv"), header = TRUE)
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
                         function(X){length(unique(X$genotype))}) == 3]
# Add LGC6.1/2 back in, as only 5 heterokaryotypes were found
ArcSax2[["LGC6.1-2"]] <- ArcSax[["LGC6.1-2"]]

### Also remove "LGC14.1-2", due to uncertain clustering
Nsax2 <- Nsax2[which(names(Nsax2) != "LGC14.1-2")]
ArcSax2 <- ArcSax2[which(names(ArcSax2) != "LGC14.1-2")]

### Additional filters in ArcSax when arcana forms a separate cluster to Lsax
Nax2 <- Nax2[which(names(ArcSax2) != "LGC9.2")]
ArcSax2 <- ArcSax2[which(!(names(ArcSax2) %in% c("LGC5.1", "LGC12.2", "LGC12.3", "LGC14.2")))]


#### C: Functions for ecotype contrast ####

### Base function
Eco_contrast <- function(PCA.data, inversion, ecotypes){
  
  print(paste0("Analysing ", inversion, "..."))
  
  ### Calculate number of ecotypes
  Neco <- length(ecotypes)
  
  ### Error messages
  if(class(inversion) != "character") stop(print(paste0(inversion, " must be an object of class 'character'.")))
  if(class(ecotypes) != "character") stop(print("'ecotypes' must be a character vector."))
  if(Neco < 2) stop(print("'ecotypes' must have a length ≥ 2."))
  if(class(PCA.data) != "list") stop(print("'PCA.data' must be a list of PCA data results. 1 inversion per list element."))
  if(class(PCA.data[[1]]) != "data.frame") stop(print("Each element of 'PCA.data' must be a data.frame with samples as rows and eigenvectors as columns.\n
                                                      This dataset MUST also have an 'Ecotype' and 'genotype' column."))
  
  ### Filter to sites with both specified ecotypes
  tmp <- PCA.data[[inversion]]
  tmp2 <- tmp %>% group_by(location) %>% 
    filter(any(Ecotype == ecotypes[1])) %>% 
    filter(any(Ecotype == ecotypes[2])) %>%
    filter(Ecotype %in% ecotypes)
  
  ### Summarise arrangements as alleles counts
  if(length(unique(tmp2$genotype)) <= 3){
    tmp3 <- tmp2 %>% group_by(location, Ecotype) %>% 
      summarise(nR = 2*sum(genotype == "RR") + sum(genotype == "RA"),
                nA = 2*sum(genotype == "AA") + sum(genotype == "RA"))
  } else {
    tmp3 <- tmp2 %>% group_by(location, Ecotype) %>% 
      summarise(nR = 2*sum(genotype == "RR") + sum(genotype == "RA1") + sum(genotype == "RA2"),
                nA = 2*sum(genotype == "A1A1") + 2*sum(genotype == "A2A2") + 2*sum(genotype == "A1A2") +
                  sum(genotype == "RA1") + sum(genotype == "RA2"))}
  
  
  ### Create glm  model
  try({
    mod0 <- glm(cbind(nR, nA) ~ location, data = tmp3, family = binomial(link = "logit"))
    mod1 <- glm(cbind(nR, nA) ~ location + Ecotype, data = tmp3, family = binomial(link = "logit"))
    mod2 <- glm(cbind(nR, nA) ~ location * Ecotype, data = tmp3, family = binomial(link = "logit"))
    
    ### Loop through each model summarising the each model as a deviance table
    mods <- list("mod0" = mod0, "mod1" = mod1, "mod2" = mod2)
    DevTab <- lapply(names(mods), function(X){
      DT <- anova(mods[[X]])
      DT2 <- data.frame("model" = X,
                        "AIC" = mods[[X]]$aic,
                        "predictor" = rownames(DT), 
                        "deviance" = DT$Deviance,
                        "df" = DT$Df,
                        "residual_deviance" = DT$`Resid. Dev`,
                        "residual_df" = DT$`Resid. Df`,
                        "p-value" = pchisq(DT$Deviance, DT$Df, lower.tail = FALSE),
                        "Rsq" = 1 - (DT$`Resid. Dev` / DT$`Resid. Dev`[2]) ### Double check ###
                        )
     
       # Save results
      if(identical(ecotypes, c("Crab", "Wave"))) Test <- "CW"
      if(identical(ecotypes, c("Wave", "Barnacle"))) Test <- "WB"

      if(identical(ecotypes, c("Crab", "arcana"))) Test <- "Carc"
      if(identical(ecotypes, c("Wave", "arcana"))) Test <- "Warc"
      
      write.csv(DT2, paste0(PATH, "Ecotype_contrasts/GLM_results_", inversion, "_", 
                            Test, "_", X, "_v3.csv"), row.names = FALSE, quote = FALSE)
      
      # Return table
      return(DT2)
    })
    
    DevTab <- do.call(rbind.data.frame, DevTab)
    
    ### Evaluate the effect of the ecotype on arrangement
    # AIC_0 > AIC_1 & AIC_2 > AIC_1: "Ecotype effect"
    # AIC_0 > AIC_1 & AIC_2 < AIC_1: "Ecotype effect, with interaction"
    # AIC_0 < AIC_1 & AIC_2 < AIC_0: "Interaction effect. Ecotype has indirect effect on inversions"
    # AIC_0 < AIC_1 & AIC_2 > AIC_0: "No effect. Any differences driven by location"
    ### Setting a ∆AIC threshold of 2. If ∆AIC < 2 then select simpler model.
    
    ## Create data frame
    res <- tmp3
    res$Inv <- inversion
    res$AIC_loc <- mod0$aic
    res$AIC_eco <- mod1$aic
    res$AIC_int <- mod2$aic
    # Evaluate effect of ecotype based on AIC
    res$Effect <- if(mod0$aic - mod1$aic > 2){
      ifelse(mod1$aic - mod2$aic > 2, "INT", "ECO")
    } else {
      ifelse(mod0$aic - mod2$aic > 2, "INT", "NULL")
    }
    # Calculate R^2 for best model
    res$Rsquare <- if(mod0$aic - mod1$aic > 2){
      ifelse(mod1$aic - mod2$aic > 2,
             1 - (mod2$deviance / mod2$null.deviance),
             1 - (mod1$deviance / mod1$null.deviance))
    } else {
      ifelse(mod0$aic - mod2$aic > 2,
             1 - (mod2$deviance / mod2$null.deviance),
             1 - (mod0$deviance / mod0$null.deviance))
    }
    
    # Re-order columns
    res <- res[,c("Inv", "location", "Ecotype", "nR", "nA", 
                    "AIC_loc", "AIC_eco", "AIC_int", "Effect", "Rsquare")]
    
    return(res)
  })
}

### For double inversions
# Compares the A1 and A2 clusters
Eco_contrast_dbInv <- function(PCA.data, inversion, ecotypes){
  
  print(paste0("Analysing ", inversion, "..."))
  
  ### Calculate number of ecotypes
  Neco <- length(ecotypes)
  
  ### Error messages
  if(class(inversion) != "character") stop(print(paste0(inversion, " must be an object of class 'character'.")))
  if(class(ecotypes) != "character") stop(print("'ecotypes' must be a character vector."))
  if(Neco < 2) stop(print("'ecotypes' must have a length ≥ 2."))
  if(class(PCA.data) != "list") stop(print("'PCA.data' must be a list of PCA data results. 1 inversion per list element."))
  if(class(PCA.data[[1]]) != "data.frame") stop(print("Each element of 'PCA.data' must be a data.frame with samples as rows and eigenvectors as columns.\n
                                                      This dataset MUST also have an 'Ecotype' and 'genotype' column."))
  
  ### Filter to sites with both specified ecotypes
  tmp <- PCA.data[[inversion]]
  tmp2 <- tmp %>% group_by(location) %>% 
    filter(any(Ecotype == ecotypes[1])) %>% 
    filter(any(Ecotype == ecotypes[2])) %>%
    filter(Ecotype %in% ecotypes)
  
  ### Summarise arrangements as alleles counts
  tmp3 <- tmp2 %>% group_by(location, Ecotype) %>% 
    summarise(nA1 = 2*sum(genotype == "A1A1") + sum(genotype == "RA1") + sum(genotype == "A1A2"),
              nA2 = 2*sum(genotype == "A2A2") + sum(genotype == "A1A2") + sum(genotype == "RA2"))
  
  
  ### Create glm  model
  try({
    mod0 <- glm(cbind(nA1, nA2) ~ location, data = tmp3, family = binomial(link = "logit"))
    mod1 <- glm(cbind(nA1, nA2) ~ location + Ecotype, data = tmp3, family = binomial(link = "logit"))
    mod2 <- glm(cbind(nA1, nA2) ~ location * Ecotype, data = tmp3, family = binomial(link = "logit"))
    
    ### Loop trough each model summarising the each model as a deviance table
    mods <- list("mod0" = mod0, "mod1" = mod1, "mod2" = mod2)
    DevTab <- lapply(names(mods), function(X){
      DT <- anova(mods[[X]])
      DT2 <- data.frame("model" = X,
                        "AIC" = mods[[X]]$aic,
                        "predictor" = rownames(DT), 
                        "deviance" = DT$Deviance,
                        "df" = DT$Df,
                        "residual_deviance" = DT$`Resid. Dev`,
                        "residual_df" = DT$`Resid. Df`,
                        "p-value" = pchisq(DT$Deviance, DT$Df, lower.tail = FALSE),
                        "Rsq" = 1 - (DT$`Resid. Dev` / DT$`Resid. Dev`[2]) ### Double check ###
      )
      
      # Save results
      if(identical(ecotypes, c("Crab", "Wave"))) Test <- "CW"
      if(identical(ecotypes, c("Wave", "Barnacle"))) Test <- "WB"
      
      if(identical(ecotypes, c("Crab", "arcana"))) Test <- "Carc"
      if(identical(ecotypes, c("Wave", "arcana"))) Test <- "Warc"
      
      write.csv(DT2, paste0(PATH, "Ecotype_contrasts/GLM_results_", inversion, "_", 
                            Test, "_", X, "_doubleInv_v3.csv"), row.names = FALSE, quote = FALSE)
      
      # Return table
      return(DT2)
    })
    
    DevTab <- do.call(rbind.data.frame, DevTab)
    
    ### Evaluate the effect of the ecotype on arrangement
    # AIC_0 > AIC_1 & AIC_2 > AIC_1: AIC_1 - AIC_0 > 2 = "Ecotype effect"
    # AIC_0 > AIC_1 & AIC_2 < AIC_1: AIC_2 - AIC_1 > 2 = "Ecotype effect, with interaction"
    # AIC_0 < AIC_1 & AIC_2 < AIC_0: AIC_2 - AIC_0 > 2 = "?"
    # AIC_0 < AIC_1 & AIC_2 > AIC_0: AIC_2 - AIC_1 < 2 = "No effect. Any differences driven by location"
    ### Setting a ∆AIC threshold of 2
    
    ## Create data frame
    res <- tmp3
    res$Inv <- inversion
    res$AIC_loc <- mod0$aic
    res$AIC_eco <- mod1$aic
    res$AIC_int <- mod2$aic
    # Evaluate effect of ecotype based on AIC
    res$Effect <- if(mod0$aic - mod1$aic > 2){
      ifelse(mod1$aic - mod2$aic > 2, "INT", "ECO")
    } else {
      ifelse(mod0$aic - mod2$aic > 2, "INT", "NULL")
    }
    # Calculate R^2 for best model
    res$Rsquare <- if(mod0$aic - mod1$aic > 2){
      ifelse(mod1$aic - mod2$aic > 2,
             1 - (mod2$deviance / mod2$null.deviance),
             1 - (mod1$deviance / mod1$null.deviance))
    } else {
      ifelse(mod0$aic - mod2$aic > 2,
             1 - (mod2$deviance / mod2$null.deviance),
             1 - (mod0$deviance / mod0$null.deviance))
    }
    
    # Re-order columns
    res <- res[,c("Inv", "location", "Ecotype", "nA1", "nA2", 
                  "AIC_loc", "AIC_eco", "AIC_int", "Effect", "Rsquare")]
    
    return(res)
  })
}

### For crab-arcana contrast
# Remove location since crab and arcana were only co-sampled in Holyhead
Carc_contrast <- function(PCA.data, inversion){
  
  print(paste0("Analysing ", inversion, "..."))
  
  ### Error messages
  if(class(inversion) != "character") stop(print(paste0(inversion, " must be an object of class 'character'.")))
  if(class(PCA.data) != "list") stop(print("'PCA.data' must be a list of PCA data results. 1 inversion per list element."))
  if(class(PCA.data[[1]]) != "data.frame") stop(print("Each element of 'PCA.data' must be a data.frame with samples as rows and eigenvectors as columns.\n
                                                      This dataset MUST also have an 'Ecotype' and 'genotype' column."))
  
  ### Filter to sites with crab and arcana
  tmp <- PCA.data[[inversion]]
  tmp2 <- tmp %>% group_by(location) %>% 
    filter(any(Ecotype == "Crab")) %>% 
    filter(any(Ecotype == "arcana")) %>%
    filter(Ecotype %in% c("Crab", "arcana"))
  
  ### Summarise arrangements as alleles counts
  if(length(unique(tmp2$genotype)) <= 3){
    tmp3 <- tmp2 %>% group_by(location, Ecotype) %>% 
      summarise(nR = 2*sum(genotype == "RR") + sum(genotype == "RA"),
                nA = 2*sum(genotype == "AA") + sum(genotype == "RA"))
  } else {
    tmp3 <- tmp2 %>% group_by(location, Ecotype) %>% 
      summarise(nR = 2*sum(genotype == "RR") + sum(genotype == "RA1") + sum(genotype == "RA2"),
                nA = 2*sum(genotype == "A1A1") + 2*sum(genotype == "A2A2") + 2*sum(genotype == "A1A2") +
                  sum(genotype == "RA1") + sum(genotype == "RA2"))}
  
  
  ### Create glm  model
  try({
    mod0 <- glm(cbind(nR, nA) ~ 1, data = tmp3, family = binomial(link = "logit"))
    mod1 <- glm(cbind(nR, nA) ~ Ecotype, data = tmp3, family = binomial(link = "logit"))
    
    ### Loop through each model summarising the each model as a deviance table
    mods <- list("mod0" = mod0, "mod1" = mod1)
    DevTab <- lapply(names(mods), function(X){
      DT <- anova(mods[[X]])
      DT2 <- data.frame("model" = X,
                        "AIC" = mods[[X]]$aic,
                        "predictor" = rownames(DT), 
                        "deviance" = DT$Deviance,
                        "df" = DT$Df,
                        "residual_deviance" = DT$`Resid. Dev`,
                        "residual_df" = DT$`Resid. Df`,
                        "p-value" = pchisq(DT$Deviance, DT$Df, lower.tail = FALSE),
                        "Rsq" = 1 - (DT$`Resid. Dev` / DT$`Resid. Dev`[2]) ### Double check ###
      )
      
      # Save results
      write.csv(DT2, paste0(PATH, "Ecotype_contrasts/GLM_results_", inversion, "_Carc_", X, "_v3.csv"), row.names = FALSE, quote = FALSE)
      
      # Return table
      return(DT2)
    })
    
    DevTab <- do.call(rbind.data.frame, DevTab)
    
    ### Evaluate the effect of the ecotype on arrangement
    # AIC_0 > AIC_1 & AIC_2 > AIC_1: "Ecotype effect"
    # AIC_0 > AIC_1 & AIC_2 < AIC_1: "Ecotype effect, with interaction"
    # AIC_0 < AIC_1 & AIC_2 < AIC_0: "Interaction effect. Ecotype has indirect effect on inversions"
    # AIC_0 < AIC_1 & AIC_2 > AIC_0: "No effect. Any differences driven by location"
    ### Setting a ∆AIC threshold of 2. If ∆AIC < 2 then select simpler model.
    
    ## Create data frame
    res <- tmp3
    res$Inv <- inversion
    res$AIC_null <- mod0$aic
    res$AIC_eco <- mod1$aic
    # Evaluate effect of ecotype based on AIC
    res$Effect <- ifelse(mod0$aic - mod1$aic > 2, "ECO", "NULL")
    # Calculate R^2 for best model
    res$Rsquare <- ifelse(mod0$aic - mod1$aic > 2,
                          1 - (mod1$deviance / mod1$null.deviance), 0)
    
    # Re-order columns
    res <- res[,c("Inv", "location", "Ecotype", "nR", "nA", 
                  "AIC_null", "AIC_eco", "Effect", "Rsquare")]
    
    return(res)
  })
}


#### D: Run function ####

### CW: crab - wave ecotype comparison
CW <- lapply(names(Nsax2), Eco_contrast, PCA.data = Nsax2, ecotypes = c("Crab", "Wave"))
CW <- do.call(rbind.data.frame, CW)

### WB: wave - barnacle ecotype comparison
WB <- lapply(names(Nsax2), Eco_contrast, PCA.data = Nsax2, ecotypes = c("Wave", "Barnacle"))
WB <- do.call(rbind.data.frame, WB)

### Ward: arcana - wave ecotype comparison
Warc <- lapply(names(ArcSax2), Eco_contrast, PCA.data = ArcSax2, ecotypes = c("Wave", "arcana"))
Warc <- do.call(rbind.data.frame, Warc)

### Carc: arcana - crab ecotype comparison
Carc <- lapply(names(ArcSax2), Carc_contrast, PCA.data = ArcSax2)
Carc <- do.call(rbind.data.frame, Carc)

### Visaulize R^2
p1 <- ggplot(CW)+
  geom_point(aes(x = Rsquare, y = Inv, colour = Effect))+
  labs(x = expression(R^2), y = "Inversion", title = "crab-wave")+
  scale_y_discrete(limits = rev(unique(CW$Inv)))+
  theme_bw()
p2 <- ggplot(WB)+
  geom_point(aes(x = Rsquare, y = Inv, colour = Effect))+
  labs(x = expression(R^2), y = "Inversion", title = "wave-barn")+
  scale_y_discrete(limits = rev(unique(CW$Inv)))+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())
p3 <- ggplot(Carc)+
  geom_point(aes(x = Rsquare, y = Inv, colour = Effect))+
  labs(x = expression(R^2), y = "Inversion", title = "crab-arcana")+
  scale_y_discrete(limits = rev(unique(CW$Inv)))+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())
p4 <- ggplot(Warc)+
  geom_point(aes(x = Rsquare, y = Inv, colour = Effect))+
  labs(x = expression(R^2), y = "Inversion", title = "wave-arcana")+
  scale_y_discrete(limits = rev(unique(CW$Inv)))+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())
ggarrange(p1,p2,p3,p4, common.legend = TRUE, nrow = 1, ncol = 4, widths = c(1.7,1,1,1))

#### E: Modification for complex inversions ####
### The previous contrast merged the A1 & A2 arrangements into one when counting double inversion's
### frequencies. We cannot determine if there was also an ecotype difference for these arrangements
### unless they are compared to one another.

### Modified function


### Run function
CW_doubleInv <- rbind.data.frame(
  Eco_contrast_dbInv(Nsax, inversion = "LGC6.1-2", ecotypes = c("Crab", "Wave")),
  Eco_contrast_dbInv(Nsax, inversion = "LGC14.2", ecotypes = c("Crab", "Wave")))
WB_doubleInv <- rbind.data.frame(
  Eco_contrast_dbInv(Nsax2, inversion = "LGC6.1-2", ecotypes = c("Wave", "Barnacle")),
  Eco_contrast_dbInv(Nsax2, inversion = "LGC14.2", ecotypes = c("Wave", "Barnacle")))
Warc_doubleInv <- Eco_contrast_dbInv(ArcSax2, inversion = "LGC6.1-2", ecotypes = c("Wave", "arcana"))
# Note: Both LGC14.2 did not form 6 clusters in the arc-sax comparison.

#### F: Run complex inversion modification separately for LGC6.1/2:crab-arcana ####
### crab L. saxatilia was only colloected with L.arcana in one location (Holyhead), making the
### location effect in the logistic model invalid. Given that this is just a single inversion and 
### site, it doesn't make sense to write a separate function, hence we are doing this part manually.

# Access data
tmp <- ArcSax2[["LGC6.1-2"]]

# Filter to location where both species were collected
tmp2 <- tmp %>% group_by(location) %>% 
  filter(any(Ecotype == "Crab")) %>% 
  filter(any(Ecotype == "arcana")) %>%
  filter(Ecotype %in% c("Crab", "arcana"))

# A1 and A2 arrangement counts
tmp3 <- tmp2 %>% group_by(location, Ecotype) %>% 
  summarise(nA1 = 2*sum(genotype == "A1A1") + sum(genotype == "RA1") + sum(genotype == "A1A2"),
            nA2 = 2*sum(genotype == "A2A2") + sum(genotype == "RA2") + sum(genotype == "A1A2"))

# Model fit
mod0 <- glm(cbind(nA1, nA2) ~ 1, data = tmp3, family = binomial(link = "logit"))
mod1 <- glm(cbind(nA1, nA2) ~ Ecotype, data = tmp3, family = binomial(link = "logit"))

### Loop through each model summarising the each model as a deviance table
mods <- list("mod0" = mod0, "mod1" = mod1)
DevTab <- lapply(names(mods), function(X){
  DT <- anova(mods[[X]])
  DT2 <- data.frame("model" = X,
                    "AIC" = mods[[X]]$aic,
                    "predictor" = rownames(DT), 
                    "deviance" = DT$Deviance,
                    "df" = DT$Df,
                    "residual_deviance" = DT$`Resid. Dev`,
                    "residual_df" = DT$`Resid. Df`,
                    "p-value" = pchisq(DT$Deviance, DT$Df, lower.tail = FALSE),
                    "Rsq" = 1 - (DT$`Resid. Dev` / DT$`Resid. Dev`[2]) ### Double check ###
  )
  
  # Save results
  write.csv(DT2, paste0(PATH, "Ecotype_contrasts/GLM_results_LGC6.1-2_Carc_", X, "_doubleInv_v3.csv"), row.names = FALSE, quote = FALSE)
  
  # Return table
  return(DT2)
})

DevTab <- do.call(rbind.data.frame, DevTab)

## Create data frame
Carc_doubleInv <- tmp3
Carc_doubleInv$Inv <- "LGC6.1-2"
Carc_doubleInv$AIC_null <- mod0$aic
Carc_doubleInv$AIC_eco <- mod1$aic
# Evaluate effect of ecotype based on AIC
Carc_doubleInv$Effect <- ifelse(mod0$aic - mod1$aic > 2, "ECO", "NULL")
# Calculate R^2 for best model
Carc_doubleInv$Rsquare <- ifelse(mod0$aic - mod1$aic > 2,
                      1 - (mod1$deviance / mod1$null.deviance), 0)

# Re-order columns
Carc_doubleInv <- Carc_doubleInv[,c("Inv", "location", "Ecotype", "nA1", "nA2", 
                                    "AIC_null", "AIC_eco", "Effect", "Rsquare")]

#### G: Save results ####

### Northern saxatilis contrasts
write.csv(CW, paste0(PATH, "Ecotype_contrasts/crab_wave_contrast.v5.csv"), quote = FALSE, row.names = FALSE)
write.csv(WB, paste0(PATH, "Ecotype_contrasts/wave_barnacle_contrast.v5.csv"), quote = FALSE, row.names = FALSE)
# Double inversions
write.csv(CW_doubleInv, paste0(PATH, "Ecotype_contrasts/crab_wave_contrast_doubleInv.v5.csv"), quote = FALSE, row.names = FALSE)
write.csv(WB_doubleInv, paste0(PATH, "Ecotype_contrasts/wave_barnacle_contrast_doubleInv.v5.csv"), quote = FALSE, row.names = FALSE)

### arcana - saxatilis contrasts
write.csv(Carc, paste0(PATH, "Ecotype_contrasts/crab_arc_contrast.v5.csv"), quote = FALSE, row.names = FALSE)
write.csv(Warc, paste0(PATH, "Ecotype_contrasts/wave_arc_contrast.v5.csv"), quote = FALSE, row.names = FALSE)
# Double inversions
write.csv(Carc_doubleInv, paste0(PATH, "Ecotype_contrasts/crab_arcana_contrast_doubleInv.v5.csv"), quote = FALSE, row.names = FALSE)
write.csv(Warc_doubleInv, paste0(PATH, "Ecotype_contrasts/wave_arcana_contrast_doubleInv.v5.csv"), quote = FALSE, row.names = FALSE)
