######### Inversion detection project: Figure 5 ##########
### This script generates figure 5 for a project detecting
### inversions throughout the range of Littorina saxatilis.
### The figure is a composite of three panels, i) map of
### sampling sites coloured by region, ii) Inversion
### frequencies of each ecotype & iii) a tile plot summarising
### the results of a contrast between Littorina arcana and
### the wave and crab ecotypes of Littorina saxatilis.
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
ArcSax <- lapply(INVs, function(inv){
  LG <- gsub("C", "", gsub("[.].*", "", inv))
  dat <- read.csv(paste0("PCA_per_inversion/", LG, "_PCA_of_", inv, "_v2_arc-sax.csv"), header = TRUE)
  # Renaming columns to remove inversion prefix
  colnames(dat) <- gsub(paste0(inv,"."), "", colnames(dat))
  return(dat)
})
# Add inversion as listnames
names(ArcSax) <- INVs
# Remove inversions with uncertain clustering
ArcSax2 <- ArcSax[which(!(names(ArcSax) %in% c("LGC5.1", "LGC12.2", "LGC12.3", "LGC14.2", "LGC14.3")))]


### Access ecotype contrast data
Carc <- read.csv("Ecotype_contrasts/crab_arc_contrast.v5.csv")
Carc$Test <- "crab-wave" # Add test name

Warc <- read.csv("Ecotype_contrasts/wave_arc_contrast.v5.csv")
Warc$Test <- "wave-barnacle" # Add test name

# Double inversions
Carc_dbInv <- read.csv("Ecotype_contrasts/crab_arcana_contrast_doubleInv.v5.csv")
Carc_dbInv$Test <- "crab-wave_dbInv"

Warc_dbInv <- read.csv("Ecotype_contrasts/wave_arcana_contrast_doubleInv.v5.csv")
Warc_dbInv$Test <- "wave-barnacle_dbInv"


#### B: Panel A - PCA of genotypes ####

### Create PCA plots
# LGC6.1/2
pLG6 <- ggplot(ArcSax2[["LGC6.1-2"]])+
  geom_point(aes(Axis1, Axis2, colour = factor(genotype), pch = Species))+
  annotate("text", x = sum(range(ArcSax2[["LGC6.1-2"]]$Axis1))/2, 
           y = max(ArcSax2[["LGC6.1-2"]]$Axis2)*0.9, 
           label = "LGC6.1/2", fontface = "bold", size = 8)+
  scale_color_manual(values = c("#1f78b4", "#a6cee3", "#33a02c", "#fb9a99", "#e31a1c"))+
  scale_shape_manual(values = c(1, 16), guide = "none")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "grey50", fill = NA),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

# LGC17.1
pLG17 <- ggplot(ArcSax2[["LGC17.1"]])+
  geom_point(aes(Axis1, Axis2, colour = factor(genotype), pch = Species))+
  annotate("text", x = sum(range(ArcSax2[["LGC17.1"]]$Axis1))/2,
           y = max(ArcSax2[["LGC17.1"]]$Axis2)*0.9, 
           label = "LGC17.1", fontface = "bold", size = 8)+
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#e31a1c"))+
  scale_shape_manual(values = c(1, 16), guide = "none")+
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
inv_counts <- lapply(names(ArcSax2), function(inv){
  tmp <- ArcSax2[[inv]]
  
  if(length(unique(tmp$genotype)) <= 3){
    # 3 cluster inversions
    tmp2 <- tmp %>% group_by(Ecotype) %>%
      filter(Ecotype %in% c("arcana", "Crab", "Wave")) %>%
      summarise(R = 2*sum(genotype == "RR") + sum(genotype == "RA"),
                A = 2*sum(genotype == "AA") + sum(genotype == "RA"),
                N = 2*n(),
                Inv = inv) %>% 
      pivot_longer(R:A, names_to = "Arrangement", values_to = "Count")
  } else {
    # 6 cluster inversions
    tmp2 <- tmp %>% group_by(Ecotype) %>%
      filter(Ecotype %in% c("arcana", "Crab", "Wave")) %>%
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
inv_counts$Inv <- factor(inv_counts$Inv, levels = names(ArcSax2))
inv_counts$Ecotype <- factor(inv_counts$Ecotype,
                             levels = c("arcana", "Crab", "Wave"))

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
  geom_point(data = Carc, aes(x = "A", y = Inv, 
                            fill = Effect != "NULL", 
                            colour = Effect == "INT", 
                            pch = Effect == "INT"), size = 3)+
  geom_point(data = Warc, aes(x = "B", y = Inv, 
                            fill = Effect != "NULL", 
                            colour = Effect == "INT", 
                            pch = Effect == "INT"), size = 3)+
  geom_point(data = Carc_dbInv, aes(x = "C", y = Inv, 
                                  fill = Effect != "NULL", 
                                  colour = Effect == "INT", 
                                  pch = Effect == "INT"), size = 3)+
  geom_point(data = Warc_dbInv, aes(x = "D", y = Inv, 
                                  fill = Effect != "NULL", 
                                  colour = Effect == "INT", 
                                  pch = Effect == "INT"), size = 3)+
  scale_fill_manual(values = c("white", "#FF8E84"))+
  scale_colour_manual(values = c("grey50", "#CC1706"))+
  scale_shape_manual(values = c(21, 22))+
  scale_y_discrete(limits = rev(names(ArcSax2)))+
  scale_x_discrete(labels = rep(c("crab-arcana", "wave-arcana"), 2))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey80", size = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.position = "none")


#### E: Save panels as individual images ####
PATH <- "/Users/james/Dropbox (Personal)/PhD/Inv_detection_manuscript/"
# Save p1
ggsave("Fig5a.tiff", plot = p1, "tiff", PATH, width = 10, height = 16, units = "cm")
# Save p2
ggsave("Fig5b.tiff", plot = p2, "tiff", PATH, width = 8, height = 24, units = "cm")
# Save p3
ggsave("Fig5c.tiff", plot = p3, "tiff", PATH, width = 6, height = 16, units = "cm")
