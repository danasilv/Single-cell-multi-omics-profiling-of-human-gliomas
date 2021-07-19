library(diversitree)
library(ape)
library(phytools)
library(adephylo)
library(reshape2)

#################### Cross-correlation matrices #########################
# used in Figure 4 and Extended Figure 10
#########################################################################

# load xcor functions
source("xcor_functions.R")
# load all phylogenies
tree.list <- readRDS("treelist.RDS")

# get mean cross-correlation Z-score matrices
xcor115 <- get_xcor_matrix("MGH115")
xcor122 <- get_xcor_matrix("MGH122")
xcor45 <- get_xcor_matrix("MGH45")
xcor64 <- get_xcor_matrix("MGH64")
xcor107 <- get_xcor_matrix("MGH107")
xcor142 <- get_xcor_matrix(c("MGH142_Plate1", "MGH142_Plate2"))
xcor208 <- get_xcor_matrix(c("MGH208_Plate1", "MGH208_Plate2"))

d115 <- melt(xcor115)
d115$name <- "MGH115"

d122 <- melt(xcor122)
d122$name <- "MGH122"

d45 <- melt(xcor45)
d45$name <- "MGH45"

d64 <- melt(xcor64)
d64$name <- "MGH64"

d107 <- melt(xcor107)
d107$name <- "MGH107"

d142 <- melt(xcor142)
d142$name <- "MGH142"

d208 <- melt(xcor208)
d208$name <- "MGH208"

dfx <- rbind(d115, d122, d45, d64, d107, d142, d208)
fn <- paste0(paste("xcor_matrix_out", gsub(" ", "_", Sys.time()),sep = "_"), ".RDS")
saveRDS(dfx, file=fn)
#
