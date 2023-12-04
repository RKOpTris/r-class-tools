# easily split data randomly into train, test and holdout sets using predefined proportions
# if var is a specified factor, classes will be (approximately) proportionally represented in the returned datasets 
# default is 0.7 train : 0.15 test : 0.15 holdout

tth_split <- function(df, var = NULL, test_prop = 0.15, holdout_prop = 0.15, shuffle_train = F){
  stopifnot(is.numeric(test_prop), test_prop >= 0.05, holdout_prop <= test_prop, is.numeric(holdout_prop))
  if(!is.null(var)){
    tth_split <- split(df, df[[var]])
  } else {
    tth_split <- list(df)
  }
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

# tth_split(iris, test_prop = 0.15, holdout_prop = 0.15)
# tth_split(iris, var = "Species", test_prop = 0.15, holdout_prop = 0.15)

# used in conjunction with list_to_objects from https://github.com/RKOpTris/r-useful-funs
# tth_split(iris, var = "Species", test_prop = 0.15, holdout_prop = 0.15) %>% list_to_objects("iris")