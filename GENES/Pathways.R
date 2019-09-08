# Pathways

library(doParallel)   # parallel computing
library(dplyr)        # data manipulation
library(glmnet)       # regularized logistic regression
library(parallel)     # parallel computing
library(purrr)        # functional programming
library(randomForest) # random forest
library(RRF)          # guided regularized random forest
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

# simulation
RunPathways <- function(trials, samples, Response) {
  
  scores <- data.frame(Score = numeric(models.num * trials),
                       Pathway1_Score = numeric(models.num * trials),
                       Pathway2_Score = numeric(models.num * trials),
                       Model = c("RF", "LASSO", paste0("GRRF-", gamma)),
                       stringsAsFactors = FALSE)
  
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  for (t in seq(trials)) {
    
    # generate genes
    genes <- matrix(nrow = samples, ncol = pathways.num * pathway.lengths)
    for (i in seq(pathways.num * pathway.lengths)) {
      genes[, i] <- runif(samples)
    }
    
    # split genes into separate pathways
    my.split <- split(seq(pathways.num * pathway.lengths), 
                      paste0("Pathway", ceiling(seq(pathways.num * pathway.lengths) / pathway.lengths)))
    SplitGenes <- function(x) { genes[, x] }
    pathways <- lapply(my.split, SplitGenes)
    
    # add in redundant genes to significant pathways
    pathways$Pathway1[, seq(pathway.lengths - length(significant.genes1) + 1: pathway.lengths)] <- pathways$Pathway1[, significant.genes1]
    pathways$Pathway2[, seq(pathway.lengths - length(significant.genes2) + 1: pathway.lengths)] <- pathways$Pathway2[, significant.genes2 - pathway.lengths]
    
    # calculate response
    sig.genes1 <- pathways$Pathway1[, significant.genes1]
    sig.genes2 <- pathways$Pathway2[, significant.genes2 - pathway.lengths]
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
    RunRF <- function(pathway) {
      
      # random forest
      rf.model <- suppressWarnings(varSelRF(pathway, y, ntree = 1000))
      error <- min(rf.model$selec.history$OOB)
      return(error)
    }
    
    # obtain error rates and rank
    rf.results <- map_dbl(pathways, RunRF)
    rf.rank <- rank(rf.results)
    
    # store the rank of worst performing pathway
    scores[row.counter, 1:3] <- c(max(rf.rank["Pathway1"], rf.rank["Pathway2"]),
                                  rf.rank["Pathway1"],
                                  rf.rank["Pathway2"])
    
    # increment row
    row.counter <- row.counter + 1
    
    # Run LASSO
    RunLASSO <- function(pathway) {
      
      # regularized logistic regression
      lasso.model <- cv.glmnet(pathway, y, family = "binomial", alpha = 1, 
                               parallel = TRUE, lambda = c(0.01, 0.02, 0.03, 0.05, 0.1, 0.2))
      error <- min(lasso.model$cvm)
      return(error)
    }
    
    # obtain error rates and rank
    lasso.results <- map_dbl(pathways, RunLASSO)
    lasso.rank <- rank(lasso.results)
    
    # store the rank of worst performing pathway
    scores[row.counter, 1:3] <- c(max(lasso.rank["Pathway1"], lasso.rank["Pathway2"]),
                                  lasso.rank["Pathway1"],
                                  lasso.rank["Pathway2"])
    
    # increment row counter
    row.counter <- row.counter + 1
    
    # Run GRRF
    RunGRRF <- function(pathway, gamma) {
      
      # ordinary random forest. 
      rf <- RRF(pathway, y, flagReg = 0)
      impRF <- rf$importance
      impRF <- impRF[, "MeanDecreaseGini"]
      
      # guided regularized random forest
      imp <- impRF / (max(impRF)) 
      coefReg <- (1 - gamma) + gamma * imp 
      grrf.model <- RRF(pathway, y, coefReg = coefReg, flagReg = 1)
      error <- tail(grrf.model$err.rate[, "OOB"], 1)
      
      return(error)
      
    }
    
    for(g in gamma) {
      
      # obtain error rates and rank 
      grrf.results <- pmap_dbl(list(pathways, rep(g, pathways.num)), RunGRRF)
      grrf.rank <- rank(grrf.results)
      
      # store the rank of worst performing pathway
      scores[row.counter, 1:3] <- c(max(grrf.rank["Pathway1"], grrf.rank["Pathway2"]),
                                    grrf.rank["Pathway1"],
                                    grrf.rank["Pathway2"])
      
      # increment row counter
      row.counter <- row.counter + 1
      
    }
    
  }
  
  stopCluster(cl); registerDoSEQ()
  
  results <- scores %>%
    group_by(Model) %>%
    summarise(MeanScore = mean(Score),
              StdDevScore = sd(Score),
              MedianScore = median(Score),
              MeanPathway1 = mean(Pathway1_Score),
              StdDevPathway1 = sd(Pathway1_Score),
              MedianPathway1 = median(Pathway1_Score),
              MeanPathway2 = mean(Pathway2_Score),
              StdDevPathway2 = sd(Pathway2_Score),
              MedianPathway2 = median(Pathway2_Score))
  
  return(results)
  
}

# the values of gamma to be used in GRRF models
gamma <- c(0.05, 0.25, 0.50, 0.75, 0.95)

# lengths of each pathway
pathway.lengths <- 15

# total number of pathways
pathways.num <- 10

# number of GRRF models + RF model + LASSO model
models.num <- length(gamma) + 2 

# set the first five genes in Pathways 1-2 as significant
significant.genes1 <- 1:5
significant.genes2 <- pathway.lengths + 1:5

# run the simulations
models <- c(Response1, Response2, Response3)
for (m in seq(length(models))) {
  set.seed(m)
  results <- RunPathways(trials = 100, samples = 50, models[[m]])
  write.csv(results, file = paste0("Pathways_RM", m, ".csv"), row.names = FALSE)
}