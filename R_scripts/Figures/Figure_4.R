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
setwd("")
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse)
library(ggpubr)

### Parameter
thresh <- 0.01 # FDR threshold for all tests

#### A: Access the test results ####
### Acess data
T1 <- read.csv(paste0("Ecotype_contrasts/", list.files("Ecotype_contrasts", pattern = "T1")[2]))
T1$Test <- "T1" # Add test name

T2a <- read.csv(paste0("Ecotype_contrasts/", list.files("Ecotype_contrasts", pattern = "T2a")[2]))
T2a$Test <- "T2a" # Add test name

T2b <- read.csv(paste0("Ecotype_contrasts/", list.files("Ecotype_contrasts", pattern = "T2b")[2]))
T2b$Test <- "T2b" # Add test name

T2c <- read.csv(paste0("Ecotype_contrasts/", list.files("Ecotype_contrasts", pattern = "T2c")[2]))
T2c$Test <- "T2c" # Add test name


#### B: Panel A - PCA of genotypes ####

### Get PCA data for Northern saxatilis - for select inversions
# LGC6.1/2
LG6 <- read.csv("PCA_per_inversion/LG6_PCA_of_LGC6.1-2_sax.csv")
colnames(LG6) <- gsub("LGC6.1.2.", "", colnames(LG6))

# LGC17.1
LG17 <- read.csv("PCA_per_inversion/LG17_PCA_of_LGC17.1_sax.csv")
colnames(LG17) <- gsub("LGC17.1.", "", colnames(LG17))

### Create PCA plots
# LGC6.1/2
pLG6 <- ggplot(LG6)+
  geom_point(aes(Axis1, Axis2, colour = genotype))+
  annotate("text", x = sum(range(LG6$Axis1))/2, y = max(LG6$Axis2)*0.9, 
           label = "LGC6.1/2", fontface = "bold", size = 8)+
  scale_color_manual(values = c("#1f78b4", "#a6cee3", "#33a02c", "#b2df8a", "#fb9a99", "#e31a1c"))+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "grey50", fill = NA),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

# LGC17.1
pLG17 <- ggplot(LG17)+
  geom_point(aes(Axis1, Axis2, colour = genotype))+
  annotate("text", x = sum(range(LG17$Axis1))/2, y = max(LG17$Axis2)*0.9, 
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

#### C: Panel B - Arragement frequency pie charts ####

# Create an ordered list of inversions
INVs <- unique(T1$Inv)

# Plot
p2 <- ggplot(T1 %>% group_by(Inv, Ecotype, Arrangement) %>% summarise(Nt = sum(N)))+
  geom_bar(aes(x = "", y = Nt, fill = Arrangement), 
           width = 1, stat = "identity", position = "fill")+
  coord_polar("y", start = 0)+
  facet_grid(factor(Inv, levels = INVs) ~ Ecotype, switch = "both")+
  scale_fill_manual(values = c("#1f78b4", "#33a02c", "#e31a1c"))+
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
  geom_tile(data = T1, aes(x = "T1: Ecotype contrast", y = Inv, fill = adj.P.value < thresh), col = "grey50")+
  geom_tile(data = T2a, aes(x = "T2a: crab-wave", y = Inv, fill = adj.P.value < thresh), col = "grey50")+
  geom_tile(data = T2b, aes(x = "T2b: crab-brackish", y = Inv, fill = adj.P.value < thresh), col = "grey50")+
  geom_tile(data = T2c, aes(x = "T2c: wave-barnacle", y = Inv, fill = adj.P.value < thresh), col = "grey50")+
  geom_vline(xintercept = 1.5, col = "black")+
  scale_fill_manual(values = c("grey75", "grey25"), guide = "none")+
  scale_y_discrete(limits = rev(INVs))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank())

#### E: Save panels as individual images ####
PATH <- "/PATH/to/PLOTS"
# Save p1
ggsave("Fig4a.tiff", plot = p1, "tiff", PATH, width = 10, height = 16, units = "cm")
# Save p2
ggsave("Fig4b.tiff", plot = p2, "tiff", PATH, width = 8, height = 24, units = "cm")
# Save p3
ggsave("Fig4c.tiff", plot = p3, "tiff", PATH, width = 6, height = 16, units = "cm")
