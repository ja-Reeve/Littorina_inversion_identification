######### Inversion detection project: Figure 4 ##########
### This script generates figure 4 for a project detecting
### inversions throughout the range of Littorina saxatilis.
### The figure is a composite of three panels, i) map of
### sampling sites coloured by region, ii) Inversion
### frequencies of each ecotype & iii) a tile plot summarising
### the results of each ecotype contrast.
### James Reeve - University of Gothenburg
### 2022-02-18

### Preparation
rm(list = ls())
dev.off()
setwd("/Users/james/Documents/Inversion_detection/")
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse)
library(ggpubr)

### Parameters
INVs <- c("LGC1.1", "LGC1.2", "LGC2.1", "LGC4.1", "LGC5.1", "LGC6.1-2", "LGC7.1", 
          "LGC7.2", "LGC9.1", "LGC9.2", "LGC10.1", "LGC10.2", "LGC11.1", "LGC12.1", 
          "LGC12.2", "LGC12.3", "LGC12.4", "LGC14.1", "LGC14.2", "LGC14.3", "LGC17.1")

#### A: Access data ####
### PCA results
Nsax <- lapply(INVs, function(inv){
  LG <- gsub("C", "", gsub("[.].*", "", inv))
  dat <- read.csv(paste0("PCA_per_inversion/", LG, "_PCA_of_", inv, "_v2_NS.csv"), header = TRUE)
  # Renaming columns to remove inversion prefix
  colnames(dat) <- gsub(paste0(inv,"."), "", colnames(dat))
  return(dat)
})
# Add inversion as listnames
names(Nsax) <- INVs
# Remove inversions with uncertain clustering
Nsax2 <- Nsax[which(!(names(Nsax) %in% c("LGC9.2", "LGC14.3")))]


### Access ecotype contrast data
CW <- read.csv("Ecotype_contrasts/crab_wave_contrast.v5.csv")
CW$Test <- "crab-wave" # Add test name

WB <- read.csv("Ecotype_contrasts/wave_barnacle_contrast.v5.csv")
WB$Test <- "wave-barnacle" # Add test name

# Double inversions
CW_dbInv <- read.csv("Ecotype_contrasts/crab_wave_contrast_doubleInv.v5.csv")
CW_dbInv$Test <- "crab-wave_dbInv"

WB_dbInv <- read.csv("Ecotype_contrasts/wave_barnacle_contrast_doubleInv.v5.csv")
WB_dbInv$Test <- "wave-barnacle_dbInv"


#### B: Panel A - PCA of genotypes ####

### Create PCA plots
# LGC6.1/2
pLG6 <- ggplot(Nsax2[["LGC6.1-2"]])+
  geom_point(aes(Axis1, Axis2, colour = genotype))+
  annotate("text", 
           x = sum(range(Nsax2[["LGC6.1-2"]]$Axis1))/2,
           y = max(Nsax2[["LGC6.1-2"]]$Axis2)*0.9, 
           label = "LGC6.1/2", fontface = "bold", size = 8)+
  scale_color_manual(values = c("#1f78b4", "#a6cee3", "#33a02c", "#b2df8a", "#fb9a99", "#e31a1c"))+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "grey50", fill = NA),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

# LGC17.1
pLG17 <- ggplot(Nsax2[["LGC17.1"]])+
  geom_point(aes(Axis1, Axis2, colour = genotype))+
  annotate("text", x = sum(range(Nsax2[["LGC17.1"]]$Axis1))/2,
           y = max(Nsax2[["LGC17.1"]]$Axis2)*0.9, 
           label = "LGC17.1", fontface = "bold", size = 8)+
  scale_color_manual(values = c("#1f78b4", "#a4db77ff", "#e31a1c"))+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "grey50", fill = NA),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

### Multi-panel plot
p1 <- ggarrange(pLG6, pLG17, ncol = 1, nrow = 2)

#### C: Panel B - Arrangement frequency pie charts ####

# Function to count inversions
inv_counts <- lapply(names(Nsax2), function(inv){
  tmp <- Nsax2[[inv]]
  
  if(length(unique(tmp$genotype)) <= 3){
    # 3 cluster inversions
    tmp2 <- tmp %>% group_by(Ecotype) %>%
      summarise(R = 2*sum(genotype == "RR") + sum(genotype == "RA"),
                A = 2*sum(genotype == "AA") + sum(genotype == "RA"),
                N = 2*n(),
                Inv = inv) %>% 
      pivot_longer(R:A, names_to = "Arrangement", values_to = "Count")
    } else {
    # 6 cluster inversions
    tmp2 <- tmp %>% group_by(Ecotype) %>%
      summarise(R = 2*sum(genotype == "RR") + sum(genotype == "RA1") + sum(genotype == "RA2"),
                A1 = 2*sum(genotype == "A1A1") + sum(genotype == "RA1") + sum(genotype == "A1A2"),
                A2 = 2*sum(genotype == "A2A2") + sum(genotype == "RA2") + sum(genotype == "A1A2"),
                N = 2*n(),
                Inv = inv) %>% 
      pivot_longer(R:A2, names_to = "Arrangement", values_to = "Count")}
  return(tmp2)
})

inv_counts <- do.call(rbind.data.frame, inv_counts)

# Set order for Ecotype and Inv
inv_counts$Inv <- factor(inv_counts$Inv, levels = names(Nsax2))
inv_counts$Ecotype <- factor(inv_counts$Ecotype,
                             levels = c("Crab", "Wave", "Barnacle", "Brackish", "Other"))

# Plot
p2 <- ggplot(inv_counts)+
  geom_bar(aes(x = "", y = Count, fill = Arrangement), 
           stat = "identity", position = "fill")+
  coord_polar("y", start = 0)+
  facet_grid(rows = vars(Inv), cols = vars(Ecotype), switch = "both")+
  scale_fill_manual(values = c("#1f78b4", "#1f78b4", "#33a02c", "#e31a1c"))+
  theme_bw()+
  theme(panel.spacing = unit(0, "cm"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text.x = element_text(size = 8, margin = margin(0.5, 1, 0.5, 1), face = "bold", angle = 90),
        strip.text.y = element_text(size = 8, margin = margin(1, 0.5, 1, 0.5), angle = 0),
        legend.position = "none")

#### D: Panel C - Tile plot of test results ####

p3 <- ggplot()+
  geom_point(data = CW, aes(x = "A", y = Inv, 
                            fill = Effect != "NULL", 
                            colour = Effect == "INT", 
                            pch = Effect == "INT"), size = 3)+
  geom_point(data = WB, aes(x = "B", y = Inv, 
                            fill = Effect != "NULL", 
                            colour = Effect == "INT", 
                            pch = Effect == "INT"), size = 3)+
  geom_point(data = CW_dbInv, aes(x = "C", y = Inv, 
                            fill = Effect != "NULL", 
                            colour = Effect == "INT", 
                            pch = Effect == "INT"), size = 3)+
  geom_point(data = WB_dbInv, aes(x = "D", y = Inv, 
                            fill = Effect != "NULL", 
                            colour = Effect == "INT", 
                            pch = Effect == "INT"), size = 3)+
  scale_fill_manual(values = c("white", "#9ecae1"))+
  scale_colour_manual(values = c("grey50", "#3182bd"))+
  scale_shape_manual(values = c(21, 22))+
  scale_y_discrete(limits = rev(names(Nsax2)))+
  scale_x_discrete(labels = rep(c("crab-wave", "wave-barnacle"), 2))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey80", size = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.position = "none")

#### E: Save panels as individual images ####
PATH <- "/Users/james/Dropbox (Personal)/PhD/Inv_detection_manuscript/"
# Save p1
ggsave("Fig4a.tiff", plot = p1, "tiff", PATH, width = 10, height = 16, units = "cm")
# Save p2
ggsave("Fig4b.tiff", plot = p2, "tiff", PATH, width = 8, height = 24, units = "cm")
# Save p3
ggsave("Fig4c.tiff", plot = p3, "tiff", PATH, width = 6, height = 16, units = "cm")
