# quickly split data randomly into train, test and holdback sets using predefined proportions 

tth_split <- function(df, var, test_prop = 0.3, holdback_test_prop = 1/3, shuffle_train = F){
  tth_split <- split(df, df[[var]])
  split_data <- function(df, test_prop, holdback_test_prop){
    test <- sample(nrow(df), floor(test_prop * nrow(df)))
    holdback <- sample(test, floor(holdback_test_prop * length(test)))
    df$.TTH <- "train"
    df$.TTH[test] <- "test"
    df$.TTH[holdback] <- "holdback"
    df
  }
  tth_data <- lapply(tth_split, split_data, test_prop, holdback_test_prop)
  tth_data <- do.call(rbind, tth_data)
  row.names(tth_data) <- 1:nrow(tth_data)
  tth_data <- tth_data %>% split(tth_data$.TTH)
  tth_data <- lapply(tth_data, function(df){df$.TTH <- NULL; df})
  if(shuffle_train){
    tth_data$train[[var]] <- shuffle(tth_data$train[[var]])
  }
  tth_data
}

# tth_split(iris, "Species", 0.3, 1/3)
# used in conjunction with list_to_objects from https://github.com/RKOpTris/r-useful-funs
# tth_split %>% list_to_objects("iris")