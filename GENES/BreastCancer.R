# Breast Cancer Data

library(caret)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(RRF)

# run GRRF
RunGRRF <- function(x, g = 0) {
  
  # run GRRF with minimal regularization, i.e., RRF(1)
  if (g == 0) {
    return(RRF(x, genetic.data$Y, coefReg = 1, flagReg = 1))
  }
  
  # otherwise, run GRRF with gamma = g
  rf <- RRF(x, genetic.data$Y, flagReg = 0)
  impRF <- rf$importance
  impRF <- impRF[,"MeanDecreaseGini"]
  imp <- impRF / (max(impRF)) 
  coefReg <- (1 - g) + g*imp
  
  return(RRF(x, genetic.data$Y, coefReg = coefReg, flagReg = 1))
  
}

# pathway data
pathways <- read_file("http://people.duke.edu/~hp44/data/breast.gmt") %>% # read data into string
  str_replace_all("\n", "") %>% # remove new line characters
  str_replace_all("\tnone", "") %>% # remove "none" genes
  str_split("\t\r") %>% # split string into a list of pathways
  map(str_replace_all, "\t", ", ") %>% # replace all tabs
  map(str_split_fixed, ", ", 2) %>%  # split pathway names and genes
  data.frame(stringsAsFactors = FALSE) %>% # store as a data frame
  head(-1) # remove the last row (empty)

colnames(pathways) <- c("Pathway", "Genes")

# gene expression data
gene.expression <- read_delim("http://people.duke.edu/~hp44/data/breast.gct", "\t", skip = 2)
gene.names <- gene.expression$probeID
gene.expression <- data.frame(t(gene.expression[, -c(1, 2)]))
colnames(gene.expression) <- gene.names

# response
response <- read_file("http://people.duke.edu/~hp44/breast.cls") %>%
  str_extract("(\\d\\s){48}\\d$") %>%
  str_split(" ", simplify = TRUE) %>%
  factor(labels = c("AT", "BT", "LT"))

# combine gene expression data and response
genetic.data <- data.frame(Y = response, gene.expression, check.names = FALSE)

# data frame for results
results <- data.frame(Pathway = pathways$Pathway,
                      NumberOfGenes = pathways$Genes %>% 
                        map_int(function(x) { str_split(x, ", ") %>% unlist %>% length } ),
                      Gamma = numeric(length(pathways$Pathway)),
                      Error = numeric(length(pathways$Pathway)))

# values of gamma to be used in GRRF
gamma <- seq(from = 0.05, to = 0.95, by = 0.05)

# number of folds to use in cross validation
nfolds <- 10

# generate (stratified) folds for k-fold CV
folds <- createFolds(genetic.data$Y, k = nfolds, list = TRUE, returnTrain = FALSE)

# counters for results and cv.results
results.row <- 1

# set a seed
set.seed(123)

for (p in pathways$Genes) {
  
  # get the gene expression data for pathway p
  dataX <- dplyr::select(genetic.data, one_of(str_split(p, ", ") %>% unlist))
  
  if (dim(dataX)[2] < 10) {
    
    # run RRF(1) if less than 10 genes in pathway
    grrf.model <- RunGRRF(dataX, g = 0)
    
    # store results
    results[results.row, c("Gamma", "Error")] <- c(0, tail(grrf.model$err.rate[, "OOB"], 1))
    results.row <- results.row + 1
    
    # move on to next pathway
    next
  }
  
  # a vector to CV results for each gamma
  cv.results <- data.frame(Gamma = gamma, Error = numeric(length(gamma)))
  
  for (g in gamma) {
    
    # a vector to store results for each fold
    validation.results <- numeric(nfolds)
    
    for (n in seq(nfolds)) {
      
      # split data into training/validation sets
      train.set <- dataX[-folds[[n]], ]
      validation.set <- dataX[folds[[n]], ]
      
      # run GRRF and store results
      grrf.model <- RunGRRF(dataX, g)
      validation.results[n] <- tail(grrf.model$err.rate[, "OOB"], 1)
      
    }
    
    # store CV error
    cv.results[cv.results$Gamma == g, "Error"] <- mean(validation.results)
    
  }
  
  # store results
  results[results.row, c("Gamma", "Error")] <- cv.results[which.min(cv.results$Error), ]
  results.row <- results.row + 1
  
}

# arrange results by error
results <- results %>% arrange(Error)

# store results
write.csv(results, file = "BreastCancerResults.csv", row.names = FALSE)