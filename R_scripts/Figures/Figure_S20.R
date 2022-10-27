### Plot supp mat for geographic divergence PCA plots ####

### Preparation
rm(list = ls())
dev.off()
setwd("")
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse)
library(ggpubr)

### List of inversions
INVs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC5.1", "LGC6.1-2", "LGC7.1", 
          "LGC7.2", "LGC9.1", "LGC9.2", "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", 
          "LGC12.2", "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.2", "LGC14.3", "LGC17.1")

### Set colour and shape scales
# Colour by genotype
K3_col <- c("#1f78b4", "#a4db77ff", "#e31a1c")
K6_col <- c("#1f78b4", "#a6cee3", "#33a02c", "#b2df8a", "#fb9a99", "#e31a1c")
# Shape by ecotype
EC_col <- c(1, 17, 15, 18, 8, 19)


### PCA plot function
PCA.plotter <- function(inversion){
  
  if(inversion != "LGC6.1-2") LG <- gsub("C", "", gsub("*[.].", "", inversion)) else LG <- "LG6"
  
  ## Download PCA data
  dat <- read.csv(paste0("PCA_per_inversion/", LG, "_PCA_of_", inversion, "_arc-sax.csv"), header = TRUE)
  # Renaming columns to remove inversion prefix
  colnames(dat) <- if(inversion == "LGC6.1-2"){
    gsub("LGC6.1.2.", "", colnames(dat))
  } else {
    gsub(paste0(inversion,"."), "", colnames(dat))}
  
  ## Stop function if K != 3 | 6
  if(length(unique(dat$genotype)) != 3 & length(unique(dat$genotype)) != 6) {
    stop(paste0("Error: ", inversion, " does not form 3 or 6 genotype clusters!"))}
  
  ## Download percentage varaince on PCA axes
  Pvar <- read.csv("PCA_per_inversion/Percentage_varaince_explained.csv", header = TRUE)
  Pvar <- Pvar[Pvar$LG == LG & Pvar$Inv == inversion,]
  
  ## PCA plot
  P <- ggplot(dat)+
    geom_point(aes(Axis1, Axis2, colour = genotype, pch = Ecotype), size = 4)+
    scale_y_continuous(breaks = c(0))+
    scale_x_continuous(breaks = c(0))+
    scale_colour_manual(values = if(inversion == "LGC6.1-2" | 
                                          inversion == "LGC14.2") K6_col else K3_col)+
    scale_shape_manual(values = EC_col)+
    labs(x = paste0("PC1 (", Pvar$Axis1, "%)"), 
         y = paste0("PC2 (", Pvar$Axis2, "%)"),
         title = inversion)+
    theme(panel.background = element_blank(),
          panel.border = element_rect(colour = "grey", fill = NA, size = 2),
          panel.grid = element_line(colour = "lightgrey"),
          axis.title = element_text(size = 10, colour = "grey50"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(size = 21),
          legend.position = "none")
  
  return(P)
}

# Make plot list
P <- lapply(INVs, function(i)try(PCA.plotter(i)))
# Filter out inversions which returned an error
P <- P[sapply(P, is.ggplot)]

### Multi-panel plot
ggarrange(plotlist = P, ncol = 3, nrow = 6)
