######### Inversion detection project: Figure 2 ##########
### This script generates figure 2 for a project detecting
### inversions throughout the range of Littorina saxatilis.
### The figure is a composite of three panels, i) average
### heterozygosity per contig, ii) PC1 scores for each map
### position, & iii) a schematic of the linkage group.
### James Reeve - University of Gothenburg
### 2022-01-18
### Edit 2022-02-25: added randomoization test to top panel
### and coloured PCA by heterozygosity
### Edit 2023-03-07: Clean-up for publication

### Preparation
rm(list = ls())
dev.off()
setwd("/Users/james/Documents/Inversion_detection")
options(stringsAsFactors = FALSE, check.names = FALSE)

### Packages
library(tidyverse)
library(ggpubr)

#### A: Create functions to access data and generate plots ####
### Get lengths of each linkage group
LG <- read.table("ConsensusMap.v2.txt", header = T)
LGmax <- LG %>% group_by(LG) %>% summarise(maxMP = max(avgMP))

### 1.1 Function to get average heterozygosities after a split for each sample
extract.heterzygosities <- function(linkage.group, genetic.group, snail.ID){
  # Access split scores
  Hsplit <- read.csv(paste0("split_function_results/Hsplit_betaBin_", linkage.group, "_v2_", genetic.group, ".csv"), header = TRUE)
  Hsplit <- Hsplit[Hsplit$SnailID == snail.ID,] # Subset to given individual
  
  Xrows <- grep("S3", Hsplit$split_name) # List of rows with 3rd level split data
  row_vals <- Hsplit[Xrows, "split_name"] # Names of the 3rd level splits
  Lpos <- LGmax$maxMP[LGmax$LG == linkage.group] # End of linkage map
  
  # Create output
  res <- cbind(c(0, Hsplit$pos, Lpos),
               c(rbind(Hsplit[Xrows, "He_left"], Hsplit[Xrows, "He_right"]), tail(Hsplit[complete.cases(Hsplit), "He_right"],1)),
               Hsplit$SnailID[1])
  # Convert to data frame (Note: use cbind() above to auto-fill in SnailID for every row, but
  # cbind() only works with 1 data class; data.frame() is used to convert the class)
  res <- data.frame("pos" = as.numeric(res[,1]),
                    "avgHE" = as.numeric(res[,2]),
                    "SnailID" = res[,3])
  
  # Remove rows without a split
  res <- res[!is.na(res$pos),]
  
  return(res)
}

### 1.2 Function for heterozygosity split plots
plot.split.function <- function(linkage.group, genetic.group){
  ## 1.2.1: Get data
  # Heterozygosity scores
  Het <- read.csv(paste0("Heterozygosity_scores/", linkage.group,"_Het_v2_",genetic.group,".csv"), header = TRUE)
  
  # Calculate average heterozygosity between splits
  IDs <- unique(Het$SnailIDs)
  Hsplit <- lapply(IDs, function(i){extract.heterzygosities(linkage.group, genetic.group, snail.ID = i)})
  Hsplit <- do.call(rbind.data.frame, Hsplit)
  
  ## 1.2.2: Plot heterozygosity scores and split function results
  res_plot <- ggplot()+
    geom_point(data = Het, aes(x = avgMP, y = pHe), alpha = 0.4, size = 0.5)+
    geom_step(data = Hsplit, aes(x = pos, y = avgHE, group = SnailID), 
              alpha = 0.4, size = 0.6, colour = "blue")+
    lims(x = c(0, LGmax$maxMP[LGmax$LG == linkage.group]), y = c(0, 0.5))+
    theme_bw()+
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_blank(),
          axis.title = element_blank())
  
  return(res_plot)
}

### 2 Function to plot PCA per map position
plot.PCA.function <- function(linkage.group, genetic.group){
  # Access PCA data
  PCA <- read.table(paste0("PCA_per_window/", linkage.group, "_PCA_1cMwind_v2_", genetic.group, ".txt"), header = TRUE)
  
  # Access heterozygosity data
  Het <- read.csv(paste0("Heterozygosity_scores/", linkage.group,"_Het_v2_", genetic.group,".csv"), header = TRUE)
  
  # Calculate heterozygosity of each window
  PCA$Het <- NA
  for(i in 1:nrow(PCA)){
    tmp <- Het[Het$avgMP >= PCA$window.start[i] & 
          Het$avgMP <= PCA$window.end[i] &
          Het$SnailIDs == PCA$Snail_ID[i],]
    PCA$Het[i] <- sum(tmp$nHe) / sum(tmp$Nsnp)
  };rm(i)
  
  ## Plot PC1 along the linakge map
  res_plot <- ggplot(PCA, aes(x = pos - 0.5, y = adjPC1, alpha = PC1_percent_var, # Adjust pos by 0.5 to shift points to edge of window
                              group = Snail_ID, colour = Het))+
    geom_point(size = 0.5)+
    geom_line(size = 0.2)+
    labs(alpha = "Varaince explained")+
    xlim(c(0, LGmax$maxMP[LGmax$LG == linkage.group]))+
    scale_alpha_continuous(limits = c(10, 35))+
    scale_colour_gradient2(limits = c(0.1, 0.46), midpoint = 0.3, 
                           low = "blue", mid = "orange2", high = "firebrick")+
    theme_bw()+
    theme(axis.text = element_text(size = 8),
          axis.title = element_blank())
  
  return(res_plot)
}

### 3 Merge heterozygosity & PCA plot
plot.both.figures <- function(linkage.group, genetic.group){
  # Plot 1: number of individuals with P-vale < 0.01
  Splits <- read.csv(paste0("split_function_results/Hsplit_betaBin_", linkage.group, "_v2_", genetic.group, ".csv"), header = TRUE)
  Splits <- Splits[Splits$best_model == "Split",]
  
  tmp <- Splits %>% group_by(pos) %>% summarise("n" = n())
  tmp2 <- data.frame("pos" = 1:ceiling(LGmax$maxMP[LGmax$LG == linkage.group])) %>%
    # Add split count and remove any set any missing positions to 0 spilts
    left_join(., tmp, by = "pos") %>%
    mutate_each(~replace(., which(is.na(.)), 0))
  # Remove any positions without contigs
  Het <- read.csv(paste0("Heterozygosity_scores/", linkage.group,"_Het_v2_", genetic.group,".csv"), header = TRUE)
  tmp2 <- tmp2[tmp2$pos %in% unique(round(Het$avgMP)), ]
  
  # Randomization test
  reps <- 10000
  Rdn_means <- array(0, dim = reps)
  for (i in 1:reps) {
    smpl <- sample(1:(length(unique(Splits$SnailID))/2), 
                   size = nrow(tmp2), replace = T)
    Rdn_means[i] <- mean(smpl)
  }; rm(i)
  
  # The plot
  p1 <- ggplot(tmp2)+
    geom_col(aes(x = pos, y = n),
             width = 0.8, show.legend = FALSE)+
    lims(x = c(0, max(tmp2$pos)), 
         y = c(0, 25))+
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "grey50"),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_blank(),
          axis.title = element_blank())
  
  ## Plot 2: Heterozygosity splits
  p2 <- plot.split.function(linkage.group, genetic.group)
  
  ## Plot3 : PCA per map position
  p3 <- plot.PCA.function(linkage.group, genetic.group)
  
  ## Plot 4: inversion positions
  Invs <- read.csv("/Users/james/Dropbox (Personal)/PhD/Inversion_positions_on_new_map_v2.csv")
  Invs <- Invs[Invs$LG == linkage.group,]
  p4 <- if(nrow(Invs) > 0){
    ggplot()+
      geom_line(aes(x = c(0, ceiling(LGmax$maxMP[LGmax$LG == linkage.group])),
                    y = c(0,0)), size = 10)+
      geom_rect(data = Invs, 
                aes(xmin = Start_pos, xmax = End_pos, ymin = -0.1, ymax = 0.1),
                colour = "red", fill = "red", alpha = 0.8)+
      theme_void()
  } else {
    ggplot()+
      geom_line(aes(x = c(0, ceiling(LGmax$maxMP[LGmax$LG == linkage.group])),
                    y = c(0,0)), size = 10)+
      theme_void()
  }
  
  # Arrange both plots on common axis
  res_plot <- annotate_figure(ggarrange(p1, p2, p3, p4, ncol = 1, nrow = 4, 
                                        align = "v", common.legend = TRUE, legend = "bottom",
                                        heights = c(0.4, 1, 1, 0.05)),
                              top = text_grob(linkage.group, size = 12))
  
  return(res_plot)
}


#### B: Generate final plot for each LG & genetic group ####
PATH <- "/Users/james/Dropbox (Personal)/PhD/Inv_detection_manuscript/"
# Save plot of LG8
p1 <- plot.both.figures(linkage.group = "LG8", genetic.group = "NS")
ggsave("Fig2a.tiff", plot = p1, "tiff", PATH, width = 14.1, height = 20, units = "cm")
# Save plot of LG14
p2 <- plot.both.figures(linkage.group = "LG14", genetic.group = "NS")
ggsave("Fig2b.tiff", plot = p2, "tiff", PATH, width = 14.1, height = 20, units = "cm")
# Save plot of LG17
p3 <- plot.both.figures(linkage.group = "LG17", genetic.group = "NS")
ggsave("Fig2c.tiff", plot = p3, "tiff", PATH, width = 14.1, height = 20, units = "cm")

### Plot all in R
fig2 <- ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom")

tiff(paste0(PATH, "Figure_2.tiff"), width = 18, height = 16, units = "cm", res = 200)

annotate_figure(fig2, top = text_grob("Northern Littorina saxatilis", size = 16),
                bottom = text_grob("Linkage map position (cM)", size = 14))

dev.off()

#### C: Plot other LGs & genetic groups in supplementary materials ####
PATH2 <- "/Users/james/Documents/Inversion_detection/plots/Supp_mat_plots"

sapply(c("NS", "arc", "Sp"), function(i){
  sapply(paste0("LG", 1:17), function(ii){
    p <- plot.both.figures(linkage.group = ii, genetic.group = i)
    ggsave(paste0("Supplementary_Fig_", i, ":", ii, "_v2.tiff"), 
           plot = p, "tiff", PATH2, width = 14.3, height = 20, units = "cm")
    })
  })
