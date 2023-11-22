multiSMOTE <- function (df, classColumn, k = 5, dupSize = 0, returnDataType = F){
  # a reimplementation of smotefamily::SMOTE that works with classes > 2
  # downsample majority class option would be nice
  ncD <- ncol(df) - 1
  splitData <- split(df[-which(names(df) == classColumn)], df[[classColumn]])
  classNames <- names(splitData)
  classCounts <- sapply(splitData, nrow)
  majorityCount <- max(classCounts)
  obsCounts <- sapply(splitData, nrow)
  if(any(obsCounts <= k)){
    message(paste0("SMOTE: Some classes have observations < k. k was changed from ", k, " to ", min(obsCounts) - 1))
    k <- min(obsCounts) - 1 #### or could be ceiling(k / 2)
  }
  nGen <- sapply(obsCounts, function(x){smotefamily::n_dup_max(x + majorityCount, x, majorityCount, dupSize)})
  nGen <- nGen[nGen > 0]
  if(length(nGen) == 0){
    warning("Classes are approximately balanced. No up-sampling performed. Original data have been returned.")
    return(df)
  } else {
    syntheticData <- vector("list", length(nGen))
    for(j in 1:length(nGen)){
      minoritySet <- splitData[[names(nGen[j])]]
      minorityCount <- nrow(minoritySet)
      minoritySet <- minoritySet[sample(minorityCount), ]
      
      knear <- smotefamily::knearest(minoritySet, minoritySet, k)
      sum_dup <- unname(nGen[j])
      syn_dat = NULL
      for (i in 1:minorityCount) {
        if (is.matrix(knear)) {
          pair_idx = knear[i, ceiling(runif(sum_dup) * k)]
        } else {
          pair_idx = rep(knear[i], sum_dup)
        }
        g = runif(sum_dup)
        P_i = matrix(unlist(minoritySet[i, ]), sum_dup, ncD, byrow = TRUE)
        Q_i = as.matrix(minoritySet[pair_idx, ])
        syn_i = P_i + g * (Q_i - P_i)
        syn_dat = rbind(syn_dat, syn_i)
      }
      minoritySet[, ncD + 1] <- names(nGen[j])
      colnames(minoritySet) <- colnames(df)
      if(returnDataType){
        minoritySet$DataType <- "original"
      }
      
      rownames(syn_dat) <- NULL
      syn_dat <- data.frame(syn_dat)
      syn_dat[, ncD + 1] <- names(nGen[j])
      colnames(syn_dat) <- colnames(df)
      if(returnDataType){
        syn_dat$DataType <- "synthetic"
      }
      syntheticData[[j]] <- syn_dat
    }
    
    syntheticData <- do.call(rbind, syntheticData)
    originalData <- df
    if(returnDataType){
      originalData$DataType <- "original"
    }
    allData <- rbind(originalData, syntheticData)
    allData <- allData %>% arrange(Class)
    # new_order <- sapply(unique(df[, classColumn]), function(x){which(allData[, classColumn] == x)})
    # allData <- allData[new_order, ]
    row.names(allData) <- NULL
    
    return(allData)
  }
}