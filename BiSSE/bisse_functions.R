library(diversitree)
library(ape)
library(phytools)
library(GenSA)
library(bbmle)
library(adephylo)

##################################################### BiSSE Maximum Likelihood Estimation Functions ##############################################
# tree.list is a list of all of the phylogenies (all tree replicates from all patient samples and plates) provided in this repository
# binmat is a matrix of gene module scores and contained within each element of tree.list
# trees within tree.list are annotated by class (GBM or IDH), patient.id (e.g. MGH105), sample_name (e.g. MGH105A), and tree_rep.
# used in Figure 5 and Extended Figure 10
##################################################################################################################################################

# run BiSSE Maximum Likelihood Estimation for single-cell phylogeny
bisse.mle <- function(tree, s_frac=10^-6, num_reps=100, maxit.sa=1000, maxit.ml=1000, uppr=500, lowr=10^-4) {
    # s_frac is the sampling fraction parameter
    # MLE and SA are iterated a maximum of 1000 times
    # BiSSE parameters can range from 10^-4 to 500.

   # BiSSE flags tree MGH142_Plate1_R9 as having multifurcations/unbranched nodes; the following code resolves the error
   if(tree$sample_name == "MGH142_Plate1") {
        tree <- ape::multi2di(tree)
    }

    # cell states are determined by the maximum gene-module expression score in binmat
    # for GBM there are (in order) four gene-modules, (NPC-like, OPC-like, AC-like, MES-like) and in IDH-MUT there are three (stem-like, AC-like, OC-like)
    tree$cellstates <- max.col(tree$binmat)
    no_con_bisse <- function(tree, s_frac) {

        # phylogenies are converted to ultrametric as required by BiSSE
        ultra_tree <- phytools::force.ultrametric(tree, method = "extend")
        num_leaves <- length(ultra_tree$tip.label)

        # For simplicity, we use two cell states in this analysis, "stem" and "mature". In GBM, NPC and OPC are merged into state 0, and AC and MES into state 1. 
        # in IDH-MUT, stem is state 0, and AC and OC are merged into the mature state, 1.
        if(tree$class=="GBM") {
            ultra_tree$states[which(ultra_tree$cellstates==1)] <- 0
            ultra_tree$states[which(ultra_tree$cellstates==2)] <- 0
            ultra_tree$states[which(ultra_tree$cellstates==3)] <- 1
            ultra_tree$states[which(ultra_tree$cellstates==4)] <- 1
        }
        if (tree$class=="IDH") {
            ultra_tree$states[which(ultra_tree$cellstates==1)] <- 0
            ultra_tree$states[which(ultra_tree$cellstates==2)] <- 1
            ultra_tree$states[which(ultra_tree$cellstates==3)] <- 1
        }
        ultra_tree$tip.state <- ultra_tree$states
        names(ultra_tree$tip.state) <- ultra_tree$tip.label

        # generate BiSSE model
        nbmodel <- diversitree::make.bisse(ultra_tree, ultra_tree$tip.state, sampling.f = c(s_frac))
        return(nbmodel)
    }

    divtree <- no_con_bisse(tree, s_frac)
    model.base <- divtree

    # reformat BiSSE model
    base_m <- function(lambda0, lambda1, mu0, mu1, q01, q10) {
        parms <- c(lambda0, lambda1, mu0, mu1, q01, q10)
        out <- -model.base(parms)
        return(out)
    }

    # constrain BiSSE parameters mu0 and mu1 to 0, corresponding to a Yule (pure-birth) model
    model.no.death <- constrain(model.base, mu0~0, mu1~0)
    nodeath_m <- function(lambda0, lambda1, q01, q10) {
        parms <- abs(c(lambda0, lambda1, q01, q10))
        parms <- ifelse(parms > uppr, uppr, parms)
        out <- -model.no.death(parms, root.p=c(1,0), root=ROOT.GIVEN, condition.surv=TRUE)
        return(out)
    }

    # function to run simulated annealing algorithm for initial parameter estimates
    gfit <- function() {
        k <- GenSA::GenSA(fn = function(x) -model.no.death(x, root.p=c(1,0), root=ROOT.GIVEN, condition.surv=TRUE),
                          upper = rep(uppr, 4), lower= rep(lowr, 4), control = list(maxit=maxit.sa, threshold.stop=1e-8, verbose=F))
        out <- k$par
        names(out) <- c("lambda0", "lambda1", "q01", "q10")
            return(out)
    }

    # function to run maximum likelihood estimation
    mfit <- function(start.point) {
        k <- bbmle::mle2(nodeath_m, start=as.list(start.point), method='L-BFGS-B', optimizer='optim',
                         control=list(trace=0, maxit=maxit.ml), lower=rep(lowr, 4), upper=rep(uppr, 4))
              return(k)
    }

    # run SA, then MLE num_reps (100) times
    z <- lapply(1:num_reps, function(r) tryCatch({zz <- mfit(gfit()); return(zz)}, error = function(err) {return()}) )
       return(list("mle"=z, "sample_name"=tree$sample_name, "tree_rep"=tree$tree_rep))
   }


# Extract best MLE estimates from replicates of BiSSE MLE code above
get_best_coef <- function(mle.obj) {

    # find runs that converged without error
    conv <- which(sapply(mle.obj$mle, function(x) tryCatch({ x@details$convergence}, error=function(er) {NA}))==0)

    # extract likelihood of each run that converged without error
    LL <- sapply(mle.obj$mle[conv], function(x) tryCatch({bbmle::logLik(x)}, error=function(er) {NA}))

    # extract coefficients for each run that converged without error
    coefs <- Reduce("rbind", lapply(mle.obj$mle[conv], function(x) bbmle::coef(x) ))

    # return coefficients with the highest likelihood
    best_coef <- coefs[which.max(LL),]

    out <- data.frame(t(best_coef), "sample_name"=mle.obj$sample_name, "tree_rep"=mle.obj$tree_rep)
    return(out)
}
#
