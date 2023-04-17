######### Inversion detection project: Figure 3 ##########
### This script generates the panels for figure 3 for a 
### project detecting inversions throughout the range of 
### Littorina saxatilis. The figure shows the PCA plots for
### both the inverted and colinear components of LG6 & LG17.
### James Reeve - University of Gothenburg
### 2022-10-19

### Preparation
rm(list = ls())
dev.off()
setwd("/Users/james/Documents/Inversion_detection/")
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse)
library(ggpubr)

### Read in PCA results
LGC6 <- read.csv("PCA_per_inversion/LG6_PCA_of_LGC6.1-2_v2_all.csv")
colnames(LGC6) <- gsub("LGC6.1.2.", "", colnames(LGC6))

LG6 <- read.csv("PCA_per_inversion/LG6_PCA_of_collinear_v2_all.csv")

LGC17 <- read.csv("PCA_per_inversion/LG17_PCA_of_LGC17.1_v2_all.csv")
colnames(LGC17) <- gsub("LGC17.1.", "", colnames(LGC17))

LG17 <- read.csv("PCA_per_inversion/LG17_PCA_of_collinear_v2_all.csv")

### Add genetic group column to the data
GG_labeler <- function(data){
  data$GenGroup <- NA
  data[data$Region != "Iberia" & data$Species == "saxatilis", "GenGroup"] <- "Northern saxatilis"
  data[data$Region == "Iberia" & data$Species == "saxatilis", "GenGroup"] <- "Iberian saxatilis"
  data[data$Species == "arcana" & data$Sample_ID != "CEA_Larc_F_1", "GenGroup"] <- "Littorina arcana"
  data[data$Species == "compressa" | data$Sample_ID == "CEA_Larc_F_1", "GenGroup"] <- "Littorina compressa"
  return(data)
}

LGC6 <- GG_labeler(LGC6)
LG6 <- GG_labeler(LG6)
LGC17 <- GG_labeler(LGC17)
LG17 <- GG_labeler(LG17)

### Set colour and shape scales
# Colour by genetic group
GG_col <- c("#008080", "#ffde55", "#c73737", "#4c58e9")
# Shape by ecotype
EC_col <- c(1, 17, 8, 1, 15, 18, 19)


### Plot LGC6.1/2
ggplot(LGC6)+
  geom_point(aes(Axis1, Axis2, colour = GenGroup, pch = Ecotype), alpha = 0.8, size = 6)+
  scale_y_continuous(breaks = c(0))+
  scale_x_continuous(breaks = c(0))+
  scale_colour_manual(values = GG_col)+
  scale_shape_manual(values = EC_col)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "grey", fill = NA, size = 2),
        panel.grid = element_line(colour = "lightgrey"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

### Plot LG6
ggplot(LG6)+
  geom_point(aes(Axis1, Axis2, colour = GenGroup, pch = Ecotype), alpha = 0.8, size = 6)+
  scale_y_continuous(breaks = c(0))+
  scale_x_continuous(breaks = c(0))+
  scale_colour_manual(values = GG_col)+
  scale_shape_manual(values = EC_col)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "grey", fill = NA, size = 2),
        panel.grid = element_line(colour = "lightgrey"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

### Plot LGC17.1
ggplot(LGC17)+
  geom_point(aes(Axis1, Axis2, colour = GenGroup, pch = Ecotype), alpha = 0.8, size = 6)+
  scale_y_continuous(breaks = c(0))+
  scale_x_continuous(breaks = c(0))+
  scale_colour_manual(values = GG_col)+
  scale_shape_manual(values = EC_col)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "grey", fill = NA, size = 2),
        panel.grid = element_line(colour = "lightgrey"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

### Plot LG17
ggplot(LG17)+
  geom_point(aes(Axis1, Axis2, colour = GenGroup, pch = Ecotype), alpha = 0.8, size = 6)+
  scale_y_continuous(breaks = c(0))+
  scale_x_continuous(breaks = c(0))+
  scale_colour_manual(values = GG_col)+
  scale_shape_manual(values = EC_col)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "grey", fill = NA, size = 2),
        panel.grid = element_line(colour = "lightgrey"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
