### Plot supp mat for geographic divergence PCA plots ####

### Preparation
rm(list = ls())
dev.off()
setwd("/Users/james/Documents/Inversion_detection/")
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse)
library(ggpubr)

### List of inversions
INVs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC5.1", "LGC6.1-2", "LGC7.1", 
          "LGC7.2", "LGC9.1", "LGC9.2", "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", 
          "LGC12.2", "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.2", "LGC14.3", "LGC17.1")

### Set colour and shape scales
# Colour by genetic group
GG_col <- c("#008080", "#ffde55", "#c73737", "#4c58e9")
# Shape by ecotype
EC_col <- c(1, 17, 8, 1, 15, 18, 19)


### PCA plot function
PCA.plotter <- function(LG, inversion = "collinear"){
  
  ## Download PCA data
  dat <- read.csv(paste0("PCA_per_inversion/", LG, "_PCA_of_", inversion, "_v2_all.csv"), header = TRUE)
  
  ## Download percentage varaince on PCA axes
  Pvar <- read.csv("PCA_per_inversion/Percentage_varaince_explained_v2.csv", header = TRUE)
  Pvar <- Pvar[Pvar$LG == LG & Pvar$Inv == ifelse(inversion == "collinear", "colinear", inversion),]
  
  ## Add genetic group column to the PCA data
  dat$GenGroup <- NA
  dat[dat$Region != "Iberia" & dat$Species == "saxatilis", "GenGroup"] <- "Northern saxatilis"
  dat[dat$Region == "Iberia" & dat$Species == "saxatilis", "GenGroup"] <- "Iberian saxatilis"
  dat[dat$Species == "arcana" & dat$Sample_ID != "CEA_Larc_F_1", "GenGroup"] <- "Littorina arcana"
  dat[dat$Species == "compressa" | dat$Sample_ID == "CEA_Larc_F_1", "GenGroup"] <- "Littorina compressa"
  
  ## PCA plot
  P <- ggplot(dat)+
    geom_point(aes(Axis1, Axis2, colour = GenGroup, pch = Ecotype), alpha = 0.8, size = 6)+
    scale_y_continuous(breaks = c(0))+
    scale_x_continuous(breaks = c(0))+
    scale_colour_manual(values = GG_col)+
    scale_shape_manual(values = EC_col)+
    labs(x = paste0("PC1 (", Pvar$Axis1, "%)"), y = paste0("PC2 (", Pvar$Axis2, "%)"))+
    #coord_fixed(Pvar$Axis2 / Pvar$Axis1)+
    theme(panel.background = element_blank(),
          panel.border = element_rect(colour = "grey", fill = NA, size = 2),
          panel.grid = element_line(colour = "lightgrey"),
          axis.title = element_text(size = 24, colour = "grey50"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          panel.spacing = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
    
    return(P)
}

### Cange to the inversion and LG you want to plot:
PCA.plotter("LG1", "LGC1.1")
