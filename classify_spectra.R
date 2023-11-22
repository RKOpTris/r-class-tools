classify_spectra <- function(md = md, # metadata
                             sp = sp, # spectra
                             dataset_factor,
                             dataset_levels,
                             predict_factor,
                             mode = NULL,
                             algorithm = "ranger",
                             train_prop = 0.7,
                             cv_folds = 3,
                             shuffle_obs = F,
                             shuffle_features = F,
                             derivate = F,
                             window_size = 9,
                             resample = "smote_up",
                             smote_k = 5,
                             feature_importance = F,
                             feature_select = F,
                             peak_detect = F, # this isn't implemented!
                             peak_detection_level = "loose",
                             colinearity_cutoff = NULL,
                             boruta_maxruns = 500,
                             boruta_include = "Confirmed", # or c("Confirmed", "Tentative")
                             model_runs = 30,
                             holdout_prop = 0.1,
                             ...){
  
  classification_parameters <- named_list(list(md, sp, dataset_factor, dataset_levels, predict_factor,
                                               mode, algorithm, train_prop, holdout_prop, cv_folds, shuffle_obs, shuffle_features,
                                               derivate, window_size, resample, feature_importance, feature_select, peak_detect, 
                                               peak_detection_level, colinearity_cutoff,
                                               boruta_maxruns, boruta_include, model_runs))
  
  resample <- tolower(resample)
  stopifnot(resample %in% c("up", "down", "smote_up", "none"))
  analyse_obs <- which(md[, dataset_factor] %in% dataset_levels)
  metadata <- md[analyse_obs, ]
  rf_spectra <- sp[analyse_obs, ]
  rf_raw_spectra <- rf_spectra
  if(shuffle_features){
    rf_spectra <- lapply(rf_spectra, shuffle) %>% data.frame()
  }
  #corrplot(cor(rf_spectra), type = "upper")
  factors <- levels(as.factor((metadata[[predict_factor]])))
  rf_spectra_derivatives_uncropped <- NULL
  if(derivate){
    rf_spectra <- prospectr::savitzkyGolay(rf_spectra, p = 2, w = window_size, m = 2) %>% data.frame()
    rf_spectra_derivatives_uncropped <- rf_spectra
  }
  
  peak_spread <- 1
  orig_wavenumbers <- strip_non_numbers(names(rf_spectra))
  spec_means <- rf_spectra %>% split(metadata[[predict_factor]]) %>% lapply(colMeans) %>% call.do(rbind)
  pracma_peaks <- list()
  pracma_peaks$tight <- lapply(1:nrow(spec_means), function(x){data.frame(index = pracma::findpeaks(-spec_means[x, ], ...)[, 2],
                                                                          source = row.names(spec_means)[x])})
  pracma_peaks$loose <- lapply(1:nrow(spec_means), function(x){y <- pracma_peaks$tight[[x]]$index; data.frame(index = sort(c(y - peak_spread, y, y + peak_spread)),
                                                                                                              central_index = sort(rep(y, 3)),
                                                                                                              source = row.names(spec_means)[x])})
  pracma_peaks$tight <- lapply(pracma_peaks$tight, function(x){data.frame(mutate(x, wavenumber = orig_wavenumbers[x$index], wavenumber_names = names(rf_spectra)[x$index]))})
  pracma_peaks$tight <- lapply(1:nrow(spec_means), function(x){pracma_peaks$tight[[x]]$absorbance <- unname(spec_means[x, pracma_peaks$tight[[x]]$index]); pracma_peaks$tight[[x]]})
  pracma_peaks$tight$all <- do.call(rbind, pracma_peaks$tight) %>% unique() %>% arrange(wavenumber) %>% tibble::remove_rownames()
  pracma_peaks$loose <- lapply(pracma_peaks$loose, function(x){data.frame(mutate(x, wavenumber = orig_wavenumbers[x$index], wavenumber_names = names(rf_spectra)[x$index]))})
  pracma_peaks$loose <- lapply(1:nrow(spec_means), function(x){pracma_peaks$loose[[x]]$absorbance <- unname(spec_means[x, pracma_peaks$loose[[x]]$index]); pracma_peaks$loose[[x]]})
  pracma_peaks$loose$all <- do.call(rbind, pracma_peaks$loose) %>% unique() %>% arrange(wavenumber) %>% tibble::remove_rownames()
  pracma_peaks$loose$all$peak_max <- F
  pracma_peaks$loose$all$peak_max[pracma_peaks$loose$all$index %in% pracma_peaks$tight$all$index] <- T
  pracma_peaks_summary <- pracma_peaks$loose$all
  colinear_wavenumbers <- NULL
  
  # boruta on raw spectra
  # boruta_rf_raw_spectra <- Boruta::Boruta(rf_raw_spectra, as.factor(metadata[[predict_factor]]), doTrace = 3, maxRuns = boruta_maxruns)
  # boruta_raw_important_wn_inds <- which(boruta_rf_raw_spectra$finalDecision == boruta_include)
  # names(rf_raw_spectra)[boruta_raw_important_wn_inds]
  
  if(feature_select | feature_importance){ # or should this be in the loop?
    # get boruta important wavenumbers
    
    ##### what about running feature selection on subsets of the data?
    boruta_rf_spectra <- Boruta::Boruta(rf_spectra, as.factor(metadata[[predict_factor]]), doTrace = 3, maxRuns = boruta_maxruns)
    boruta_important_wn_inds <- which(boruta_rf_spectra$finalDecision == boruta_include)

    if(length(boruta_important_wn_inds) > 0){
      message(paste0("Boruta found ", length(boruta_important_wn_inds), " confirmed important features in ", boruta_maxruns, " iterations."))
    } else {
      stop("Boruta did not find any important features. Stopping.")
    }
    peak_filtered_boruta <- boruta_important_wn_inds[which(boruta_important_wn_inds %in% unique(pracma_peaks[[peak_detection_level]]$all$index))]
  }
  if(feature_select){
    # if centred on peak, remove loose possibilities
    peak_centred_inds <- boruta_important_wn_inds[which(boruta_important_wn_inds %in% unique(pracma_peaks[["tight"]]$all$index))]
    pracma_peaks_summary <- pracma_peaks_summary[!(pracma_peaks_summary$central_index %in% boruta_important_wn_inds & !pracma_peaks_summary$peak_max), ]
    peak_filtered_boruta <- boruta_important_wn_inds[which(boruta_important_wn_inds %in% unique(pracma_peaks_summary$index))]
    
    # clear out neighbouring wavenumbers
    bor_var <- peak_filtered_boruta
    bor_diff <- matrix(c(bor_var[-length(bor_var)], bor_var[2:length(bor_var)]), ncol = 2)
    neighbouring_indices <- which(apply(bor_diff, 1, function(x){x[2] - x[1]}) == 1)
    if(length(neighbouring_indices) > 0){
      peak_filtered_boruta <- bor_var[-which(apply(bor_diff, 1, function(x){x[2] - x[1]}) == 1)]
    }
    bor_spectra <- rf_spectra[peak_filtered_boruta]
    if(ncol(bor_spectra) > 1){
      message(paste0(ncol(bor_spectra), " important features coincided with peaks in the spectra."))
    } else {
      stop("No more than one important feature coincided with peaks in the spectra. Stopping.")
    }
    
    #corrplot::corrplot(cor(bor_spectra), type = "upper", method = "ellipse")
    if(!is.null(colinearity_cutoff)){
      colinear_wavenumbers <- findCorrelation(cor(bor_spectra), cutoff = colinearity_cutoff)
      
      if(length(colinear_wavenumbers > 0)){
        rf_spectra <- bor_spectra[, -colinear_wavenumbers]
        
        #corrplot::corrplot(cor(rf_spectra), type = "upper", method = "ellipse")
        colinear_wavenumbers <- list(colinearity_cutoff = colinearity_cutoff, 
                                     colinear_wavenumbers = names(bor_spectra)[colinear_wavenumbers],
                                     non_colinear_wavenumbers = names(bor_spectra)[-colinear_wavenumbers])
        assign("colinear_wavenumbers", colinear_wavenumbers, env = parent.frame())
        message(paste0(length(colinear_wavenumbers$colinear_wavenumbers), " colinear features (at threshold ", colinearity_cutoff, ") were removed. ", length(colinear_wavenumbers$non_colinear_wavenumbers), " remain."))
      } else {
        rf_spectra <- bor_spectra
        colinear_wavenumbers <- list(colinearity_cutoff = colinearity_cutoff, 
                                     colinear_wavenumbers = NULL,
                                     non_colinear_wavenumbers = names(bor_spectra))
        assign("colinear_wavenumbers", NULL, env = parent.frame())
        message("No colinear features were found (at threshold of ", colinearity_cutoff, ").")
      }
    } else {
      colinear_wavenumbers <- NULL
      rf_spectra <- bor_spectra
    }
    
    #corrplot::corrplot(cor(rf_spectra), type = "upper")
  } else {
    peak_filtered_boruta <- NULL
  }
  if(feature_select | feature_importance){
    boruta_important_wavenumbers <- data.frame(index = boruta_important_wn_inds, 
                                               wavenumber = orig_wavenumbers[boruta_important_wn_inds])
    #boruta_important_wavenumbers$in_peak <- boruta_important_wn_inds %in% peak_filtered_boruta
  } else {
    boruta_important_wavenumbers <- NULL
  }
  rf_spectra <- data.frame(rf_spectra, Class = metadata[[predict_factor]])
  rf_models <- vector("list", model_runs)
  confusion_matrix_train <- vector("list", model_runs)
  confusion_matrix_test <- vector("list", model_runs)
  confusion_matrix_accuracy <- vector("list", model_runs)
  confusion_matrix_stats <- vector("list", model_runs)
  predictions <- vector("list", model_runs)
  actuals <- vector("list", model_runs)
  importances <- vector("list", model_runs)
  classifications <- vector("list", model_runs)
  results_table <- data.frame(factor = rep(predict_factor, model_runs), accuracy = NA)
  results_matrix <- vector("list", model_runs)
  durations <- rep(NA, model_runs)
  
  if(holdout_prop > 0){
    holdout_inds <- sample(nrow(rf_spectra), floor(nrow(rf_spectra) * holdout_prop))
    sp_holdout <- rf_spectra[holdout_inds, ]
    rf_spectra <- rf_spectra[-holdout_inds, ]
  }
  
  for(i in 1:model_runs){
    time_start <- Sys.time()
    #set.seed(i)
    train_ind <- sample(nrow(rf_spectra), floor(nrow(rf_spectra) * train_prop))
    
    sp_train <- rf_spectra[train_ind, ]
    sp_test <- rf_spectra[-train_ind, ]
    
    if(resample == "up"){
      sp_train <- upSample(x = sp_train[, -ncol(sp_train)], y = as.factor(sp_train$Class))
      message("The minority class was simply up-sampled.")
    } else if(resample == "down"){
      sp_train <- downSample(x = sp_train[, -ncol(sp_train)], y = as.factor(sp_train$Class))
      message("The majority class was simply down-sampled.")
    } else if(resample == "smote_up"){
      sp_train <- mcSMOTE(df = sp_train, classColumn = "Class", k = smote_k, returnDataType = T)
      if("synthetic" %in% sp_train$DataType){
        message(paste0("SMOTE created ", length(sp_train$DataType == "synthetic"), " synthetic instances."))
      } else {
        message("Up-sampling not performed as classes approximately balanced.")
      }
      sp_train$DataType <- NULL
      # if(shuffle_features){
      #   sp_train[, -ncol(sp_train)] <- lapply(sp_train[, -ncol(sp_train)], shuffle)
      # }
    } else {
      sp_train <- sp_train
    }
    if(!is.null(mode)){
      if(mode == "learning_curve"){
        sp_train <- sp_train[sample(nrow(sp_train), floor(seq(1, model_runs)[i] / model_runs * nrow(sp_train))), ]
      }
    }
    
    if(shuffle_obs){
      sp_train$Class <- shuffle(sp_train$Class)
    }
    
    sp_mod <- train(Class ~ .,
                    sp_train,
                    importance = "impurity", # maybe this doesn't work for knn?!
                    method = algorithm,
                    trControl = trainControl(method = "cv",
                                             number = cv_folds,
                                             returnData = FALSE
                    )
    )
    rf_models[[i]] <- sp_mod
    sp_predictions <- predict(sp_mod, sp_test[, -ncol(sp_test)])
    predictions[[i]] <- sp_predictions
    actuals[[i]] <- sp_test$Class
    results_table$accuracy[i] <- mean(sp_predictions == sp_test$Class)
    confusion_matrix_train[[i]] <- confusionMatrix(sp_mod) # might need to include accuracy and stats for training data too
    if(length(unique(sp_test$Class)) >= length(unique(sp_predictions))){
      conf_mat <- confusionMatrix(factor(sp_predictions, levels = factors),
                                  factor(sp_test$Class, levels = factors))
      confusion_matrix_test[[i]] <- conf_mat
      confusion_matrix_accuracy[[i]] <- conf_mat$overall
      confusion_matrix_stats[[i]] <- conf_mat$byClass
    }
    results_matrix[[i]] <- confusionMatrix(sp_mod)$table 
    importances[[i]] <- varImp(sp_mod)
    classifications[[i]] <- data.frame(run = i, actual = sp_test$Class, predict = sp_predictions)
    time_end <- Sys.time()
    durations[i] <- difftime(time_end, time_start, units = "mins") %>% round(2) %>% as.numeric()
    recent_durations_average <- ifelse(i < 10, mean(durations, na.rm = T), mean(durations[(i - 10):i], na.rm = T))
    message(paste0(i, "/", model_runs, " complete. Estimated time until completion: ", round(recent_durations_average * (model_runs - i), 2), " minutes"))
  }
  confusion_matrix_accuracy <- do.call(rbind, confusion_matrix_accuracy) %>% data.frame()
  confusion_matrix_stats <- do.call(rbind, confusion_matrix_stats) %>% data.frame()
  named_list(list(classification_parameters,
                  dataset_factor, 
                  dataset_levels, 
                  predict_factor, 
                  shuffle_obs, 
                  derivate, 
                  train_prop,
                  cv_folds,
                  resample,
                  window_size,
                  model_runs,
                  results_table, 
                  results_matrix, 
                  classifications,
                  pracma_peaks_summary,
                  boruta_important_wavenumbers,
                  peak_filtered_boruta,
                  colinear_wavenumbers,
                  importances,
                  rf_models,
                  confusion_matrix_train,
                  confusion_matrix_test,
                  confusion_matrix_accuracy,
                  confusion_matrix_stats,
                  actuals,
                  predictions,
                  metadata,
                  rf_spectra,
                  rf_spectra_derivatives_uncropped,
                  sp_holdout))
}
