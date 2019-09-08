# Genes

library(doParallel)   # parallel computing
library(dplyr)        # data manipulation
library(glmnet)       # regularized logistic regression
library(parallel)     # parallel computing
library(purrr)        # functional programming
library(randomForest) # random forest
library(RRF)          # guided regularized random forest
library(stringr)      # string manipulation
library(varSelRF)     # variable selection using random forests

# first model
Response1 <- function(sig.genes1, sig.genes2) {
  
  # calculate pathway1's contribution to response
  y <- 10 * sin(pi * sig.genes1[, 1] * sig.genes1[, 2]) +
    20 * (sig.genes1[, 3] - 0.5)^2 + 
    10 * sig.genes1[, 4] +
    5 * sig.genes1[, 5]
  
  # calculate pathway2's contribution to response
  y <- y + 
    15 * sin(pi * sig.genes2[, 1] * sig.genes2[, 2]) +
    25 * (sig.genes2[, 3] - 0.5)^2 + 
    20 * sig.genes2[, 4] +
    10 * sig.genes2[, 5]
  
  return(y)
}

# second model
Response2 <- function(sig.genes1, sig.genes2) {
  
  y <- 10 * sin(pi * sig.genes1[, 1] * sig.genes2[, 1]) +
    20 * (sig.genes1[, 2] * sig.genes2[, 2] - 0.5)^2 + 
    10 * sig.genes1[, 3] + 10 * sig.genes2[, 3] +
    5 * sig.genes1[, 4] + 10 * sig.genes2[, 4] +
    10 * sig.genes1[, 5] + 5 * sig.genes2[, 5]
  
  return(y)
}

# third model
Response3 <- function(sig.genes1, sig.genes2) {
  
  y <- cbind(sig.genes1, sig.genes2) %*% rep(10, 10)
  
  return(y)
}

# function for calculating significant pairs identified
GetSignificantPairsIdentified <- function(selected.genes, significant.genes) {
  
  count <- 0
  for (s in as.character(significant.genes)) {
    if (any(str_detect(selected.genes, s))) {
      count <- count + 1
    }
  }
  
  return(count)
}

# simulation
RunGenes <- function(trials, samples, Response) {
  
  pathway1.results <- data.frame(SigPairs = numeric(models.num * trials), 
                                 NonSigPairs = numeric(models.num * trials), 
                                 Model = c("RF", "LASSO", paste0("GRRF-", gamma)),
                                 stringsAsFactors = FALSE)
  pathway2.results <- data.frame(SigPairs = numeric(models.num * trials), 
                                 NonSigPairs = numeric(models.num * trials), 
                                 Model = c("RF", "LASSO", paste0("GRRF-", gamma)),
                                 stringsAsFactors = FALSE)
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  for (t in seq(trials)) {
    
    # generate pathway1
    p <- pathway.lengths - length(significant.genes1)
    genes <- matrix(nrow = samples, ncol = p)
    for (i in seq(p)) {
      genes[, i] <- runif(samples)
    }
    colnames(genes) <- paste0("Gene", seq(p))
    pathway1 <- cbind(genes, genes[, significant.genes1])
    
    # generate pathway2
    p <- pathway.lengths - length(significant.genes2)
    genes <- matrix(nrow = samples, ncol = p)
    for (i in seq(p)) {
      genes[, i] <- runif(samples)
    }
    colnames(genes) <- paste0("Gene", (pathway.lengths + 1):(pathway.lengths + p))
    pathway2 <- cbind(genes, genes[, significant.genes2 - pathway.lengths])
    
    # calculate response
    sig.genes1 <- pathway1[, significant.genes1]
    sig.genes2 <- pathway2[, significant.genes2 - pathway.lengths]
    y <- Response(sig.genes1, sig.genes2)
    
    # add error term
    y <- y + rnorm(samples)
    
    # transform response
    y <- factor(ifelse(y > median(y), "Class2", "Class1"))
    
    # set seed
    set.seed(456)
    
    # keep track of the results row
    row.counter <- models.num * (t - 1) + 1
    
    # Run Random Forest
    RunRF <- function(pathway, y, significant.genes) {
      
      # random forest
      rf.model <- suppressWarnings(varSelRF(pathway, y, fitted.rf = randomForest(pathway, y, ntree = 1000, importance = TRUE)))
      rf.nvars <- rf.model$best.model.nvars
      rf.sig.pairs <- GetSignificantPairsIdentified(rf.model$selected.vars, significant.genes)
      return(c(rf.sig.pairs, rf.nvars - rf.sig.pairs))
      
    }
    
    pathway1.results[row.counter, 1:2] <- RunRF(pathway1, y, significant.genes1)
    pathway2.results[row.counter, 1:2] <- RunRF(pathway2, y, significant.genes2)
    
    # increment row
    row.counter <- row.counter + 1
    
    # Run LASSO
    RunLASSO <- function(pathway, y, significant.genes) {
      
      # regularized logistic regression
      lasso.model <- cv.glmnet(pathway, y, family = "binomial", alpha = 1,
                               parallel = TRUE, lambda = c(0.01, 0.02, 0.03, 0.05, 0.1, 0.2))
      lasso.vars <- colnames(pathway)[which(coef(lasso.model, s = lasso.model$lambda.min)[-1] != 0)]
      lasso.nvars <- length(lasso.vars)
      lasso.sig.pairs <- GetSignificantPairsIdentified(lasso.vars, significant.genes)
      return(c(lasso.sig.pairs, lasso.nvars - lasso.sig.pairs))
      
    }
    
    pathway1.results[row.counter, 1:2] <- RunLASSO(pathway1, y, significant.genes1)
    pathway2.results[row.counter, 1:2] <- RunLASSO(pathway2, y, significant.genes2)
    
    # increment row counter
    row.counter <- row.counter + 1
    
    # Run GRRF
    RunGRRF <- function(pathway, y, significant.genes, gamma) {
      
      #ordinary random forest. 
      rf <- RRF(pathway, y, flagReg = 0)
      impRF <- rf$importance
      impRF <- impRF[, "MeanDecreaseGini"]
      
      #guided regularized random forest
      imp <- impRF / (max(impRF))
      coefReg <- (1 - gamma) + gamma*imp
      grrf <- RRF(pathway, y, coefReg = coefReg, flagReg = 1)
      grrf.vars <- colnames(pathway[, grrf$feaSet])
      grrf.nvars <- length(grrf.vars)
      grrf.sig.pairs <- GetSignificantPairsIdentified(grrf.vars, significant.genes)
      return(c(grrf.sig.pairs, grrf.nvars - grrf.sig.pairs))
      
    }
    
    for(g in gamma) {
      
      pathway1.results[row.counter, 1:2] <- RunGRRF(pathway1, y, significant.genes1, g)
      pathway2.results[row.counter, 1:2] <- RunGRRF(pathway2, y, significant.genes2, g)
      row.counter <- row.counter + 1
      
    }
    
  }
  
  stopCluster(cl); registerDoSEQ()
  
  pathway1.results <- pathway1.results %>%
    group_by(Model) %>%
    summarise(P1_SigPairsMean = mean(SigPairs), P1_NonSigPairsMean = mean(NonSigPairs)) %>%
    data.frame
  
  pathway2.results <- pathway2.results %>%
    group_by(Model) %>% 
    summarise(P2_SigPairsMean = mean(SigPairs), P2_NonSigPairsMean = mean(NonSigPairs)) %>%
    data.frame
  
  results <- data.frame(pathway1.results, pathway2.results)
  
  return(results)
  
}

# the values of gamma to be used in GRRF models
gamma <- c(seq(0.50, 0.95, 0.05), seq(0.96, 0.99, 0.01))

# lengths of each pathway
pathway.lengths <- 15

# number of GRRF models + RF model + LASSO model
models.num <- 16 

# set the first five genes in Pathways 1-2 as significant
significant.genes1 <- 1:5
significant.genes2 <- pathway.lengths + 1:5

# run the simulations
models <- c(Response1, Response2, Response3)
for (m in seq(length(models))) {
  set.seed(m)
  results <- RunGenes(trials = 100, samples = 50, models[[m]])
  write.csv(results, file = paste0("Genes_RM", m, ".csv"), row.names = FALSE)
}