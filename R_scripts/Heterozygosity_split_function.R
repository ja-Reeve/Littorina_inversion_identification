######## Finding splits in average heterozygosity along a linkage group #######
### This code was originally created by Roger Butlin to find the position of
### inversions by looking for signficant differences in herterozygosity along a
### linkage map in Littorina saxatilis. This code was modified to get it to run
### across multiple individuals.
### Version 2 implements a beta-binomial distribution for the split function
### this will reduce the noise of P-values
### James Reeve - Göteborgs universitet
### 10/02/2021
### Edited 14/12/2021: run over 1cM windows to make the methods it compatible 
### with the PCA based approach

### 1. Preparation ####
### clear the environment
rm(list = ls())
dev.off()
setwd("")

Start_time <- Sys.time()

### Set options - this prevents an annoying bug with filenames
options(stringsAsFactors = FALSE, check.names = FALSE)

### Parameters
# Set linkage group
LG.s <- "LG1"
# Alpha for testing significance of splits
Alpha <- 0.01

### Pacakges
library(paramtest)
library(tidyverse)
library(VGAM)
library(ggpubr) # for multi-panel plot

### Access heterozygosity data
### Also filter the data to remove contigs with <10 SNPs and without any heterozygotes
# Northern saxatilis
Het_sax <- read.csv(paste0(LG.s,"_Het_sax.csv"))
Het_sax <- Het_sax[Het_sax$Nsnp >= 10 & Het_sax$nHe > 0,]

# Spain
Het_Sp <- read.csv(paste0(LG.s,"_Het_Sp.csv"))
Het_Sp <- Het_Sp[Het_Sp$Nsnp >= 10 & Het_Sp$nHe > 0,]

# Littorina arcana
Het_arc <- read.csv(paste0(LG.s,"_Het_arc.csv"))
Het_arc <- Het_arc[Het_arc$Nsnp >= 10 | Het_arc$nHe > 0,]

# Littorina compressa won't be included as there are too few individuals

### vectors
# 1cM windows of linkage map
lm_max <- ceiling(max(c(max(Het_sax$avgMP), max(Het_Sp$avgMP), max(Het_arc$avgMP)))) # max position (rounding up)
mp_winds <-  seq(1, lm_max, 1) # set of possible step positions

# individuals
sax_IDs <- as.character(unique(Het_sax$SnailIDs))
Sp_IDs <- as.character(unique(Het_Sp$SnailIDs))
arc_IDs <- as.character(unique(Het_arc$SnailIDs))

### 2. Define sub-functions ####
### Step functions created by Roger Butlin 15/12/2020
### Parameters: h=Het values, mp=map position, ns=Nsnps, ha=position of Het transition
# step function for no change in Het, i.e. null model
step0 <- function(iter,r,h,mp,ns){                        # r=dispersion parameter, h=Het values, mp=map position, ns=Nsnps
  fit <- rep(sum(h,na.rm=T)/sum(ns,na.rm=T),length(h))    # weighted mean Het
  ms <- -2*sum(dbetabinom(h,ns,fit,r,log=T))          # log-likelihood
  return(ms)
}
# step function for one change in Het
step1 <- function(iter,ha,h,mp,ns,r){                                   # h=Het values, mp=map position, n=Nsnps, ha=1cM windows
  fit <- rep(sum(h[mp<ha],na.rm=T)/sum(ns[mp<ha],na.rm=T),length(h))    # weighted mean Het to left
  fit[mp>ha] <- sum(h[mp>ha],na.rm=T)/sum(ns[mp>ha],na.rm=T)            # weighted mean Het to right
  ms <- -2*sum(dbetabinom(h,ns,fit,r,log=T))          # log-likelihood
  return(ms)
}

### Function to run step functions for a single individual
### Including a chi-squared test for significant differences between no-split and split models
### Outputs a list with P-value, Log-likelihoods, position and average heterozygosities left and right of the split
He.split <- function(data, linkageMap.cuts, alpha){
  ## A: Run no-split model, i.e. step0
  # search for best dispersal coefficient
  r_try = seq(0,0.99,0.01)  # set of dispersion coefficient to try
  gs0 <- grid_search(step0, params = list(r = r_try), n.iter = 1,
                     h = data$nHe, mp = data$avgMP, ns = data$Nsnp)
  min <- which.min(sapply(gs0$results, '[[', 1))
  m2LL0 <- gs0$results[[min]] # null model -2*logLik
  r_fit <- r_try[min]
  
  ## B: Run split model for each possible cut site on the linkage map
  # search one step possibilities, i.e. run step1 for all cuts
  gs <- grid_search(step1, params = list(ha = linkageMap.cuts), n.iter = 1,
                    h = data$nHe, mp = data$avgMP, ns = data$Nsnp, r = r_fit)
  # isolate best match
  min <- which.min(sapply(gs$results, '[[', 1))
  
  # results for best one step model
  m2LL1 <- gs$results[[min]] # -2*logLikelihood
  
  ## C: Test if split model is better than no-split
  chi_0_1 <- pchisq(m2LL0-m2LL1, 2, lower.tail = F)
  
  ## D: Print results from best one step model
  # This is all done in a single conditional if statement
  
  if(chi_0_1 < alpha) {print(paste("Split model is a better fit; P-value =", chi_0_1))
    # split position 
    ha_fit <- linkageMap.cuts[min]
    # average heterozygosity to the left of the split
    h1_p1 <- sum(data$nHe[data$avgMP < ha_fit], na.rm=T) / sum(data$Nsnp[data$avgMP < ha_fit], na.rm=T)
    # average heterozygosity to the right of the split
    h1_p2 <- sum(data$nHe[data$avgMP > ha_fit], na.rm=T) / sum(data$Nsnp[data$avgMP > ha_fit], na.rm=T) 
    
    res <- list(paste("P-value =", chi_0_1), c(m2LL0, m2LL1), "Split", ha_fit, h1_p1, h1_p2, r_fit)
    names(res) <- c("chisq_Pval", "Log_likelihoods", "Best_model", "pos", "avgHE_left", "avgHE_right", "Rho") 
    names(res$Log_likelihoods) <- c("No_split", "Split")
    
    return(res)
  } else {
    print(paste("No-split model is a better fit; P-value =", chi_0_1))
    ha_fit <- NA
    h1_p1 <- sum(data$nHe, na.rm=T) / sum(data$Nsnp, na.rm=T)
    h1_p2 <- h1_p1
    
    res <- list(paste("P-value =", chi_0_1), c(m2LL0, m2LL1), "No_split", ha_fit, h1_p1, h1_p2, r_fit) 
    names(res) <- c("chisq_Pval", "Log_likelihoods", "Best_model", "pos", "avgHE_left", "avgHE_right", "Rho")
    names(res$Log_likelihoods) <- c("No_split", "Split")
    
    return(res)
  }
}

### 3. Run split function iteratively until the 3rd split ####
### Nested split functions for a single individual
### data = heterozygosity scores dataset with the number of heterozygostes,
### snail.ID = individual being tested (character)
### linkageMap.cuts = vector of possible split positions on linkage map
### alpha = significance threshold for selecting a model
He.nest.split <- function(data, snail.ID, linkageMap.cuts, alpha){
  
  ### Split 1: ####
  print(paste("Running split 1 for", snail.ID))
  
  ## Inputs for split 1
  # possible split positions
  cuts_S1 <- linkageMap.cuts
  # heterozygosity scores for single individual
  dat_S1 <- data[data$SnailIDs == snail.ID, ]
  
  ## Output for split 1
  # Run split function
  res_S1 <- He.split(dat_S1, cuts_S1, alpha)
  # position of split
  pos_S1 <- res_S1$pos
  # P-value
  chisq_S1 <- as.numeric(sapply(strsplit(res_S1$chisq_Pval, split = " "), '[[', 3))
  
  ### Split2:####
  print(paste("Running split 2 for", snail.ID))
  
  if(chisq_S1 < alpha && sum(cuts_S1 < pos_S1) > 0 && sum(cuts_S1 > pos_S1) > 0){
    ## Inputs for split 2
    # Possible split positions
    cuts_S2.L <- cuts_S1[cuts_S1 < pos_S1]
    cuts_S2.R <- cuts_S1[cuts_S1 > pos_S1]
    # filtered heterozygosity data
    dat_S2.L <- dat_S1[dat_S1$avgMP < pos_S1,]
    dat_S2.R <- dat_S1[dat_S1$avgMP > pos_S1,]
    
    ## Output for split 2
    # Run split function
    res_S2.L <- He.split(dat_S2.L, cuts_S2.L, alpha)
    res_S2.R <- He.split(dat_S2.R, cuts_S2.R, alpha)
    # positions of new splits
    pos_S2.L <- res_S2.L$pos
    pos_S2.R <- res_S2.R$pos
    # P-value for new splits
    chisq_S2.L <- as.numeric(sapply(strsplit(res_S2.L$chisq_Pval, split = " "), '[[', 3))
    chisq_S2.R <- as.numeric(sapply(strsplit(res_S2.R$chisq_Pval, split = " "), '[[', 3))
    
    ### Split 3:
    print(paste("Running split 3 for", snail.ID))
    ## Left section of split 2
    if(chisq_S2.L < alpha && sum(cuts_S2.L < pos_S2.L) > 0 && sum(cuts_S2.L > pos_S2.L) > 0){
      ## Inputs for split 3 (left)
      # Possible split positions
      cuts_S3.L.L <- cuts_S2.L[cuts_S2.L < pos_S2.L]
      cuts_S3.L.R <- cuts_S2.L[cuts_S2.L > pos_S2.L]
      # filtered heterozygosity data
      dat_S3.L.L <- dat_S2.L[dat_S2.L$avgMP < pos_S2.L, ]
      dat_S3.L.R <- dat_S2.L[dat_S2.L$avgMP > pos_S2.L, ]
      
      ## Output for split 3 (left)
      res_S3.L.L <- He.split(dat_S3.L.L, cuts_S3.L.L, alpha)
      res_S3.L.R <- He.split(dat_S3.L.R, cuts_S3.L.R, alpha)
    } else {
      print("No further spliting is required. No split model fits best for left segment of 2nd split.")
      
      ## Negative results
      ha_fit <- NA
      h1_p1 <- sum(dat_S2.L$nHe, na.rm=T) / sum(dat_S2.L$Nsnp, na.rm=T)
      h1_p2 <- h1_p1
      
      res_S2.L_NULL <- list(paste("P-value =", chisq_S2.L), res_S2.L$Log_likelihoods, "No_split", ha_fit, h1_p1, h1_p2, res_S2.L$Rho)
      names(res_S2.L_NULL) <- c("chisq_Pval", "Log_likelihoods", "Best_model", "pos", "avgHE_left", "avgHE_right", "Rho")
      
      ## Add negative results to each possible split
      res_S3.L.L <- res_S2.L_NULL
      res_S3.L.R <- res_S2.L_NULL
    }
    
    ## Right section of split 2
    if(chisq_S2.R < alpha && sum(cuts_S2.R < pos_S2.R) > 0 && sum(cuts_S2.R > pos_S2.R) > 0){
      ## Inputs for split 3 (right)
      # possible split positions
      cuts_S3.R.L <- cuts_S2.R[cuts_S2.R < pos_S2.R]
      cuts_S3.R.R <- cuts_S2.R[cuts_S2.R > pos_S2.R]
      # filtered heterozygosity data
      dat_S3.R.L <- dat_S2.R[dat_S2.R$avgMP < pos_S2.R, ]
      dat_S3.R.R <- dat_S2.R[dat_S2.R$avgMP > pos_S2.R, ]
      
      ## Output for split 3 (right)
      res_S3.R.L <- He.split(dat_S3.R.L, cuts_S3.R.L, alpha)
      res_S3.R.R <- He.split(dat_S3.R.R, cuts_S3.R.R, alpha)
    } else {
      print("No further spliting is required. No split model fits best for right segment of 2nd split.")
      ## Negative results
      ha_fit <- NA
      h1_p1 <- sum(dat_S2.R$nHe, na.rm=T) / sum(dat_S2.R$Nsnp, na.rm=T)
      h1_p2 <- h1_p1
      
      res_S2.R_NULL <- list(paste("P-value =", chisq_S2.R), res_S2.R$Log_likelihoods, "No_split", ha_fit, h1_p1, h1_p2, res_S2.R$Rho)
      names(res_S2.R_NULL) <- c("chisq_Pval", "Log_likelihoods", "Best_model", "pos", "avgHE_left", "avgHE_right", "Rho")
      
      ## Add negative results to each possible split
      res_S3.R.L <- res_S2.R_NULL
      res_S3.R.R <- res_S2.R_NULL
    }
    
  } else {
    print("No further spliting is required. No split model fits best for 1st split.")
    
    ## Negative results
    ha_fit <- NA
    h1_p1 <- sum(data$nHe, na.rm=T) / sum(data$Nsnp, na.rm=T)
    h1_p2 <- h1_p1
    
    res_S1_NULL <- list(paste("P-value =", chisq_S1), res_S1$Log_likelihoods, "No_split", ha_fit, h1_p1, h1_p2, res_S1$Rho)
    names(res_S1_NULL) <- c("chisq_Pval", "Log_likelihoods", "Best_model", "pos", "avgHE_left", "avgHE_right", "Rho")
    
    ## Add negative results to each possible split
    res_S2.L <- res_S1_NULL
    res_S2.R <- res_S1_NULL
    res_S3.L.L <- res_S1_NULL
    res_S3.L.R <- res_S1_NULL
    res_S3.R.L <- res_S1_NULL
    res_S3.R.R <- res_S1_NULL
  }
  
  ### output results into a single dataset
  res_full <- list(res_S3.L.L, res_S2.L, res_S3.L.R, res_S1, res_S3.R.L, res_S2.R, res_S3.R.R)
  names(res_full) <- c(paste0("S3LL:",snail.ID), paste0("S2L:",snail.ID), paste0("S3LR:",snail.ID),
                       paste0("S1:",snail.ID), 
                       paste0("S3RL:",snail.ID), paste0("S2R:",snail.ID), paste0("S3RR:",snail.ID))
  
  ### Convert results into a dataframe
  tmp <- lapply(res_full, unlist)
  tmp <- tmp[which(sapply(tmp, length) == 8)] # remove any failed splits
  tmp2 <- data.frame(array(NA, dim = c(length(tmp), length(tmp[[1]])+2)))
  for(i in 1:length(tmp)){
    tmp2[i,1] <- as.numeric(sapply(strsplit(tmp[[i]][1], split = " "), '[[', 3)) # p-value
    tmp2[i,2] <- as.numeric(tmp[[i]][8]) # Rho estimate
    tmp2[i,3] <- as.numeric(tmp[[i]][2]) # LL of no split
    tmp2[i,4] <- as.numeric(tmp[[i]][3]) # LL of split
    tmp2[i,5] <- tmp[[i]][4] # Best fitting model (split | no split)
    tmp2[i,6] <- as.numeric(tmp[[i]][5]) # Position of split
    tmp2[i,7] <- as.numeric(tmp[[i]][6]) # Avg. heterozygosity left of split
    tmp2[i,8] <- as.numeric(tmp[[i]][7]) # Avg. heterozygosity right of split
    tmp2[i,9] <- sapply(strsplit(names(tmp[i]), ":"), '[[', 1) # Name of split
    tmp2[i,10] <- snail.ID
  }
  rm(i)
  colnames(tmp2) <- c("Pval", "Rho", "LL_split0", "LL_split1", "best_model", "pos", "He_left", "He_right", "split_name", "SnailID")
  
  
  return(tmp2)
}

### Run nested split function for every individual
# Northern saxatilis
Hsplit_sax <- lapply(sax_IDs, function(i){
  He.nest.split(data = Het_sax, snail.ID = i, linkageMap.cuts = mp_winds, alpha = Alpha)})
Hsplit_sax <- do.call(rbind.data.frame, Hsplit_sax)

# Spain
Hsplit_Sp <- lapply(Sp_IDs, function(i){
  He.nest.split(data = Het_Sp, snail.ID = i, linkageMap.cuts = mp_winds, alpha = Alpha)})
Hsplit_Sp <- do.call(rbind.data.frame, Hsplit_Sp)

# Littorina arcana
Hsplit_arc <- lapply(arc_IDs, function(i){
  He.nest.split(data = Het_arc, snail.ID = i, linkageMap.cuts = mp_winds, alpha = Alpha)})
Hsplit_arc <- do.call(rbind.data.frame, Hsplit_arc)

# Save results
write.csv(Hsplit_sax, paste0("/PATH/split_function_results/Hsplit_betaBin_", LG.s, "_sax.csv"), row.names = F, quote = F)
write.csv(Hsplit_Sp, paste0("/PATH/split_function_results/Hsplit_betaBin_", LG.s, "_Sp.csv"), row.names = F, quote = F)
write.csv(Hsplit_arc, paste0("/PATH/split_function_results/Hsplit_betaBin_", LG.s, "_arc.csv"), row.names = F, quote = F)

### Take relevant heterozygosity scores from dataset and create a dataset for plotting
# Function to extract relevant heterozygosity scores
# I repeated the first and last values to include positions for the ends of the linkage map
extract.heterzygosities <- function(X){
  Xrows <- grep("S3", X$split_name) # List of rows with 3rd level split data
  row_vals <- X[Xrows, "split_name"] # Names of the 3rd level splits
  Lpos <- lm_max # End of linkage map
  
  # Create output
  res <- cbind(c(0, X$pos, Lpos),
               c(rbind(X[Xrows, "He_left"], X[Xrows, "He_right"]), tail(X[complete.cases(X), "He_right"],1)),
               X$SnailID[1])
  # Convert to data frame (Note: use cbind() above to auto-fill in SnailID for every row, but
  # cbind() only works with 1 data class; data.frame() is used to convert the class)
  res <- data.frame("pos" = as.numeric(res[,1]),
                    "avgHE" = as.numeric(res[,2]),
                    "SnailID" = res[,3])
  
  # Remove rows without a split
  res <- res[!is.na(res$pos),]

  return(res)
}

# Extract het score data for each invidividual
split_scores_sax <- lapply(sax_IDs, function(i){
  extract.heterzygosities(X = Hsplit_sax[Hsplit_sax$SnailID == i,])})

split_scores_Sp <- lapply(Sp_IDs, function(i){
  extract.heterzygosities(X = Hsplit_Sp[Hsplit_Sp$SnailID == i,])})

split_scores_arc <- lapply(arc_IDs, function(i){
  extract.heterzygosities(X = Hsplit_arc[Hsplit_arc$SnailID == i,])})

# Convert into a single dataframe per region
split_scores_sax <- do.call(rbind.data.frame, split_scores_sax)
split_scores_Sp <- do.call(rbind.data.frame, split_scores_Sp)
split_scores_arc <- do.call(rbind.data.frame, split_scores_arc)

### Plot the results
# Note: this will return a warning message for any SNPs missing data (i.e. NAs in heterozygosity scores)
# or for any non-significant splits. These are dropped from the plot, so shouldn't be a problem.
plot.split.function <- function(heterzygosity.data, test.results, split.scores, plot.title, alpha){
  # Plot 1: number of individuals with P-vale < alpha
  p1 <- ggplot(test.results[test.results$best_model == "Split",])+
    geom_bar(aes(x = pos, group = Pval < alpha), fill = "red", width = 0.1)+
    labs(y = paste0("P<", alpha), title = paste0(LG.s, ": ", plot.title))+
    xlim(c(0, lm_max))+
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 24),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 14))
  
  # Plot 2: plot heterozygosity scores and split function results
  p2 <- ggplot()+
    geom_point(data = heterzygosity.data, aes(x = avgMP, y = pHe), alpha = 0.4, size = 0.5)+
    geom_step(data = split.scores, aes(x = pos, y = avgHE, group = SnailID), 
              alpha = 0.4, size = 0.6, colour = "blue")+
    labs(x = "Linkage map position (cM)", y = expression(H[Ind]~per~contig))+
    lims(x = c(0, lm_max), y = c(0, 0.5))+
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 14),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  
  # Multipanel plot
  res_plot <- ggarrange(p1, p2, ncol = 1, heights=c(1, 3), align = ("v"))
  
  return(res_plot)
}

psax <- plot.split.function(Het_sax, Hsplit_sax, split_scores_sax, "Northern saxatilis", Alpha)
pSp <- plot.split.function(Het_Sp, Hsplit_Sp, split_scores_Sp, "Spain", Alpha)
parc <- plot.split.function(Het_arc, Hsplit_arc, split_scores_arc, "Littorina arcana", Alpha)

ggarrange(psax, pSp, parc, nrow = 3, ncol = 1)

### Save plots 
# Set file path
PATH <- "/Users/james/Documents/Inversion_detection/plots/Split_plots"
ggsave(paste0(LG.s, "hetSplit_sax.tiff"), psax, "tiff", PATH, width = 30, height = 15, units = "cm")
ggsave(paste0(LG.s, "hetSplit_Sp.tiff"), pSp, "tiff", PATH, width = 30, height = 15, units = "cm")
ggsave(paste0(LG.s, "hetSplit_arc.tiff"), parc, "tiff", PATH, width = 30, height = 15, units = "cm")

End_time <- Sys.time()
End_time - Start_time
