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
setwd("")
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse)
library(ggpubr)

### Parameter
thresh <- 0.01 # FDR threshold for all tests

#### A: Access the test results ####
### Acess data
T3 <- read.csv(paste0("Ecotype_contrasts/", list.files("Ecotype_contrasts", pattern = "T3")[2]))
T3$Test <- "T3" # Add test name

T4a <- read.csv(paste0("Ecotype_contrasts/", list.files("Ecotype_contrasts", pattern = "T4a")[2]))
T4a$Test <- "T4a" # Add test name

T4b <- read.csv(paste0("Ecotype_contrasts/", list.files("Ecotype_contrasts", pattern = "T4b")[2]))
T4b$Test <- "T4b" # Add test name

#### B: Panel A - PCA of genotypes ####

### Get PCA data for Northern saxatilis - for select inversions
# LGC9.1
LG9 <- read.csv("PCA_per_inversion/LG9_PCA_of_LGC9.1_arc-sax.csv")
colnames(LG9) <- gsub("LGC9.1.", "", colnames(LG9))

# LGC17.1
LG17 <- read.csv("PCA_per_inversion/LG17_PCA_of_LGC17.1_arc-sax.csv")
colnames(LG17) <- gsub("LGC17.1.", "", colnames(LG17))

### Create a flag for the samples which were included in this ecotype contrast
arc_locs <- LG9[LG9$Species == "arcana", "location"] # locations were aranca was sampled
LG9$included <- LG9$location %in% arc_locs
LG17$included <- LG17$location %in% arc_locs

### Create PCA plots
# LGC9.1
pLG9 <- ggplot(LG9)+
  geom_point(aes(Axis1, Axis2, colour = factor(genotype), pch = Species))+
  annotate("text", x = sum(range(LG9$Axis1))/2, y = max(LG9$Axis2)*0.9, 
           label = "LGC9.1", fontface = "bold", size = 8)+
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#e31a1c"))+
  scale_shape_manual(values = c(1, 16), guide = "none")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "grey50", fill = NA),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

# LGC17.1
pLG17 <- ggplot(LG17)+
  geom_point(aes(Axis1, Axis2, colour = factor(genotype), pch = Species))+
  annotate("text", x = sum(range(LG17$Axis1))/2, y = max(LG17$Axis2)*0.9, 
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
p1 <- ggarrange(pLG9, pLG17, ncol = 1, nrow = 2)

#### C: Panel B - Arragement frequency pie charts ####

# Create an ordered list of inversions
INVs <- unique(T3$Inv)

# Plot
tmp <- rbind(T4a, T4b)
tmp <- tmp[!duplicated(tmp[,1:4]),]

p2 <- ggplot(tmp %>% group_by(Inv, Ecotype, Arrangement) %>% summarise(Nt = sum(N)))+
  geom_bar(aes(x = "", y = Nt, fill = Arrangement), 
           position="fill", stat = "identity", width = 1)+
  coord_polar("y", start = 0)+
  facet_grid(factor(Inv, levels = INVs) ~ Ecotype, switch = "both")+
  scale_fill_manual(values = c("#1f78b4", "#e31a1c"))+
  theme_bw()+
  theme(panel.spacing = unit(0, "cm"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text.x = element_text(size = 8, margin = margin(5, 10, 5, 10), face = "bold", angle = 90),
        strip.text.y = element_text(size = 8, margin = margin(10, 5, 10, 5), angle = 0),
        legend.position = "none")

#### D: Panel C - Tile plot of test results ####

p3 <- ggplot()+
  geom_tile(data = T3, aes(x = "T3: arcana-saxatilis", y = Inv, fill = adj.P.value < thresh), col = "grey50")+
  geom_tile(data = T4a, aes(x = "T4a: arcana-crab", y = Inv, fill = adj.P.value < thresh), col = "grey50")+
  geom_tile(data = T4b, aes(x = "T4b: arcana-wave", y = Inv, fill = adj.P.value < thresh), col = "grey50")+
  geom_vline(xintercept = 1.5, col = "black")+
  scale_fill_manual(values = c("grey75", "grey25"), guide = "none")+
  scale_y_discrete(limits = rev(INVs))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank())

#### E: Save panels as individual images ####
PATH <- "/PATH/to/PLOTS"
# Save p1
ggsave("Fig5a.tiff", plot = p1, "tiff", PATH, width = 10, height = 16, units = "cm")
# Save p2
ggsave("Fig5b.tiff", plot = p2, "tiff", PATH, width = 8, height = 24, units = "cm")
# Save p3
ggsave("Fig5c.tiff", plot = p3, "tiff", PATH, width = 6, height = 16, units = "cm")

