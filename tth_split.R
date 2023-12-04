# quickly split data randomly into train, test and holdout sets using a specific variable and predefined proportions 

tth_split <- function(df, var, test_prop = 0.15, holdout_prop = 0.15, shuffle_train = F){
  tth_split <- split(df, df[[var]])
  tth_data <- lapply(tth_split, function(df){
    df$.TTH <- sample(c("train", "test", "holdout"), nrow(df), replace = T, prob = c(train_prop, test_prop, holdout_prop))
    df
  })
  tth_data <- do.call(rbind, tth_data)
  row.names(tth_data) <- 1:nrow(tth_data)
  tth_data <- tth_data %>% split(tth_data$.TTH)
  tth_data <- lapply(tth_data, function(df){df$.TTH <- NULL; df})
  if(shuffle_train){
    tth_data$train[[var]] <- shuffle(tth_data$train[[var]])
  }
  tth_data
}

# tth_split(iris, "Species", 0.15, 0.15)
# used in conjunction with list_to_objects from https://github.com/RKOpTris/r-useful-funs
# tth_split %>% list_to_objects("iris")