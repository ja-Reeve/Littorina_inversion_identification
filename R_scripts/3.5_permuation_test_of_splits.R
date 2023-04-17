######## Permutation test of heterozygosity splits ########
### This script determines if the splits in heterozygosity
### scores are clustered across a linkage group. The number
### of significant splits at each cM are grouped into 3cM
### windows and then randomly shuffled using rmultinom to 
### generate a empirical distribution of variances. The 
### variance in the observed data is then compared to this
### simulated data to calculate an empirical P-value.
### James Reeve - University of Gothenburg
### 28/06/2022
### Edit 03/10/2022: clustering into windows after running
### the permutation.
### Edit 09/03/2023: clean-up for publication

### Preparation
rm(list = ls())
dev.off()
setwd("/Users/james/Documents/Inversion_detection")
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse)

### Set number of iterations
Niter <- 10000

### Function to run permutation test
permute.splits <- function(linkage.group, genetic.group, N){
  ## Get linkage group information
  LG <- read.table("ConsensusMap.v2.txt", header = T)
  LG <- LG[LG$LG == linkage.group, ]
  
  ## Find the greatest map position
  LGmax <- max(LG$avgMP)
  
  ## Get split scores
  Splits <- read.csv(paste0("split_function_results/Hsplit_betaBin_",
                            linkage.group, "_v2_", genetic.group, ".csv"), header = TRUE)
  # Filter out non-significant splits
  Splits <- Splits[Splits$best_model == "Split",]
  
  ## Count number of splits per map position
  tmp <- Splits %>% group_by(pos) %>% summarise("n" = n())
  
  ## Add empty spots for map positions which didn't have splits
  ## Note: this uses the linkage group information, so there may be gaps
  tmp2 <- data.frame("pos" = unique(ceiling(LG$avgMP)))
  tmp2$n <- ifelse(tmp2$pos %in% tmp$pos, tmp$n, 0)
  
  ## Cluter split counts into windows of 3cM. This allows us to control for
  ## clusters of splits at inversion boundaries.
  NsplitWinds <- sapply(2:(nrow(tmp2)-2), function(i){
    pos1 <- tmp2$n[i-1] # Previous cM
    pos2 <- tmp2$n[i]   # Current cM
    pos3 <- tmp2$n[i+1] # Next cM
    
    return(pos1 + pos2 + pos3)
  })
  
  ## Calculate observed variance
  ObsVar <- var(NsplitWinds)
  
  ## Randomization test
  PermVar <- array(0, dim = N)
  for (i in 1:N) {
    smpl <- rmultinom(n = 1, # Number of samples drawn
                      size = sum(tmp2$n), # Total of all splits
                      prob = rep(1 / nrow(tmp2), nrow(tmp2))) # Probability of split - even for each window
    
    # Cluster permuted splits into windows
    NpermWinds <- sapply(2:(length(smpl)-2), function(i){
      pos1 <- smpl[i-1] # Previous cM
      pos2 <- smpl[i]   # Current cM
      pos3 <- smpl[i+1] # Next cM
      
      return(pos1 + pos2 + pos3)
    })
    
    PermVar[i] <- var(NpermWinds)
  }; rm(i)
  
  ## Store are list elements and output
  Res <- list(ObsVar, PermVar)
  names(Res) <- c("observed", "permutated")
  
  return(Res)
}

### Run permutation on all linkage groups
tmp <- lapply(paste0("LG", 1:17), permute.splits, genetic.group = "arc", N = Niter)
PermSplits <- do.call(rbind.data.frame,
               lapply(tmp, function(res){
                 data.frame("Vobs" = as.numeric(res[[1]]), 
                            "Vperm" = mean(res[[2]]),
                            "Vperm_SD" = sd(res[[2]]), 
                            "P" = ifelse(sum(res[[2]] >= res[[1]]) != 0, 
                                         sum(res[[2]] >= res[[1]])/Niter,
                                         1/Niter) # if Obs >>> Perm
                            ) 
               }))
PermSplits$LG <- paste0("LG", 1:17)

### Add variance ratio
PermSplits$ratio <- PermSplits$Vobs / PermSplits$Vperm

### Adjust P-values (Benjamini-Hochberg)
PermSplits$adjP <- p.adjust(PermSplits$P, method = "BH")

### Save permutation results
write.csv(PermSplits, "Permutation_test_on_HetSplits_arc.csv", quote = FALSE, row.names = FALSE)
