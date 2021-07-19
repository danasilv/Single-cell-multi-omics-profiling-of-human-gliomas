library(diversitree)
library(ape)
library(phytools)
library(adephylo)

############################## Simulating Tumor Evolution ############################################
# simulating trees with different de-differentiation rates
# used in Figure 5
######################################################################################################

# load xcor functions
source("../xcor/xcor_functions.R")

# function to subsample simulated tree
subsample_sim <- function(tree, natips)
{
    if(length(natips) > 0) {
        tree <- ape::drop.tip(tree, natips)
        tree$tip.state <- tree$tip.state[-natips]
    }
    return(tree)
}

# Simulation script
sim <- rep(NA, 1000)
dd <- rep(NA, 1000)

# BiSSE parameters
s.growth <- 100 # lambda0
m.growth <- 35  # lambda1
mu0 <- 0
mu1 <- 0
diff <- 50 # q01

# simulates trees under the BiSSE model with a range de-differentiation rates
# all trees start in the stem-like state (x0 = 0)
for(i in 1:1000) {
  tre <- tree.bisse(c(s.growth, m.growth, mu0, mu1, diff,  i/10), max.taxa = 1000, x0 = 0)
    while (is.null(tre) != 0) {
    tre <- tree.bisse(c(s.growth, m.growth, mu0, mu1, diff, i/10), max.taxa = 1000, x0 = 0)
  }

  # subsample trees
  tre <- subsample_sim(tre, sample(1:1000, 900, replace = FALSE))

  # calculate Moran's I
  obs <- get.moran(tre$tip.state, proxTips(tre, method="nNodes",  normalize="none"), nrepet = 1000)$Z
  sim[i] <- obs
  dd[i] <- i/10
}

sdf <- data.frame("dediff"=dd, "I"=sim)
fn <- paste0(paste("sim_out", gsub(" ", "_", Sys.time()),sep = "_") , ".RDS")
saveRDS(sdf, file=fn)
#
