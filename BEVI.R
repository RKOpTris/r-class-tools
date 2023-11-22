library(parallel)
library(dplyr)
library(stringr)

# new variable setosa should be better factor for classifying than species as there is little overlap of data, setosa is very different
irisTestData <- iris %>% mutate(Setosa = str_detect(Species, "setosa"))

# some junk variables that are random and should perform poorly
irisTestData$dummyFactor1 <- sample(c(letters[1:3]), nrow(irisTestData), replace = T)
irisTestData$dummyFactor2 <- sample(c(LETTERS[1:5]), nrow(irisTestData), replace = T)

BEVI <- function(inData, inMetadata = NULL, nIter = 50, bootstrapSize = 20, numTrees = 50, parallelise = T){
  # BEVI = best explanatory variable indicator
  # in case metadata is in different data.frame
  if(!is.null(inMetadata)){
    stopifnot(nrow(inData) == nrow(inMetadata))
    inData <- cbind(inData, inMetadata)
  }
  # assume anything numeric is a predictor variable data and anything else is an explanatory variable
  predictCols <- lapply(inData, is.numeric) %>% unlist() %>% which() %>% unname()
  
  i_results <- lapply(1:nIter, function(i){
    # create bootstrap sample
    bootstrapSample <- sample(nrow(inData), bootstrapSize, replace = T)
    inSample <- inData[bootstrapSample, ]
    inTrain <- inSample[predictCols]
    inMeta <- inSample[-predictCols]
    # select explanatory factors with levels > 1
    (moreThanOneLevels <- which(sapply(inMeta, function(x) length(levels(as.factor(x)))) > 1))
    if(length(moreThanOneLevels) != 0){
      inMeta <- inMeta[moreThanOneLevels]
      j_results <- lapply(1:ncol(inMeta), function(n){
        actuals <- as.factor(inMeta[, n])
        rangerMod <- ranger(x = inTrain,
                            y = actuals,                     
                            num.trees = numTrees)
        predictions <- rangerMod$predictions
        caret::confusionMatrix(actuals, predictions)$overall[1:2] %>% as.matrix() %>% t() %>% data.frame() %>% mutate(predictFactor = names(inMeta)[n])
      }) %>% call.do(rbind)
      j_results
    }
  }) %>% call.do(rbind)
  
  outputDF <- arrange(i_results, desc(Kappa)) %>% select(Kappa, Accuracy, predictFactor) %>% group_by(predictFactor) %>% summarise_all(mean) %>% arrange(desc(Kappa)) %>% data.frame()
  message("Performed ", nIter, " * ", ncol(inMeta), " = ", nIter * ncol(inMeta), " searches")
  return(outputDF)
}

#BEVI(iris)
#BEVI(irisTestData)