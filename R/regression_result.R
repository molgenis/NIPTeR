regression_template <- function(result_set, chromo_focus, correction_status, samplenames, 
                                potential_predictors, models, sample_names_train_set = NULL, 
                                train_set_statistics = NULL, train_set_Zscores = NULL, type){
  if (is.null(train_set_statistics)){
    new_regression_template <- list(prediction_statistics = data.frame(result_set[[1]]), control_group_Zscores =  result_set[[2]],
                                    focus_chromosome = chromo_focus, correction_status = correction_status,
                                    control_group_sample_names = samplenames, models = models,
                                    potential_predictors = potential_predictors, 
                                    all_control_group_Z_scores = result_set$All_control_group_Z_scores, 
                                    additional_statistics = result_set$Additional_statistics)
  }
  else{
    new_regression_template <- list(prediction_statistics = data.frame(result_set[[1]]), control_group_Zscores =  result_set[[2]],
                                    focus_chromosome = chromo_focus, correction_status = correction_status,
                                    control_group_sample_names = samplenames, models = models,
                                    potential_predictors = potential_predictors, 
                                    sample_names_train_set = sample_names_train_set,
                                    train_set_statistics = train_set_statistics,
                                    train_set_Zscores = train_set_Zscores)
  }
  class(new_regression_template) <- c(Regressions_result_class, type)
  
  return(new_regression_template)
}
collapse_result <- function(result, value){
  return(result[[value]])
}
collapse_prediction_sets <- function(result){
  gsub(pattern = ",", replacement = " ", x = toString(result$predictors))
}
collapse_prac_cv_control_scores <- function(control_group_scores, additional_control_group_scores, n_models, cv_types, setnames){
  endmatrix <- NULL
  cols <- NULL
  for (i in 1:n_models){
    if (cv_types[i] == theoretical){
      endmatrix <- cbind(endmatrix, control_group_scores[,i], additional_control_group_scores[,i])
    }
    else{
      endmatrix <- cbind(endmatrix, additional_control_group_scores[,i], control_group_scores[,i])
    }
  }
  colnames(endmatrix) <- as.vector(rbind(paste(setnames, theoretical, sep="_"), paste(setnames, practical, sep="_"))) 
  return(endmatrix)
}
collapse_results <- function(result_set, n_models, n_predictors){
  setnames <- paste("Prediction_set", 1:n_models, sep="_")
  control_group_Z_scores <- Reduce(cbind, lapply(result_set, collapse_result, value = "control_group_Z_scores"))
  additional_control_group_Z_scores <- Reduce(cbind, lapply(result_set, collapse_result, value = "additional_control_group_Z_scores"))
  cv_types <- sapply(result_set, collapse_result, value = "cv_type")
  all_control_group_Z_scores <- collapse_prac_cv_control_scores(control_group_scores = control_group_Z_scores,
                                                                additional_control_group_scores = additional_control_group_Z_scores,
                                                                n_models = n_models, cv_types = cv_types, setnames = setnames)
  prediction_statistics <- rbind("Z_score_sample" = as.numeric(sapply(result_set, collapse_result, value = "sample_Z_score")),
                                 "CV" = as.numeric(sapply(result_set, collapse_result, value = "cv")),
                                  cv_types,
                                 "P_value_shapiro" = as.numeric(sapply(result_set, collapse_result, value = "shapiro_P_value")),
                                 "Predictor_chromosomes" = sapply(result_set, collapse_prediction_sets),
                                 "Mean_test_set" = sapply(result_set, collapse_result, value = "mean_test_set"),
                                 "CV_train_set" = sapply(result_set, collapse_result, value = "cv_train_set"))
  colnames(control_group_Z_scores) <- setnames
  colnames(prediction_statistics) <- setnames
  additional_statistics <- sapply(result_set, collapse_all_stats)
  dimnames(additional_statistics) <- list(c(rownames_additional_stats(theoretical), rownames_additional_stats(practical)),
                                              setnames)
  nipt_result <- list("PredictionStatistics" = prediction_statistics, "ControlZScores" = control_group_Z_scores, 
                      "All_control_group_Z_scores" = all_control_group_Z_scores, 
                      "Additional_statistics" = additional_statistics)
  return(nipt_result)
}

listmodels <- function(prediction_set){
  return(prediction_set$summary_model)
}

collapse_all_stats <- function(result_set){
  result <- NULL
  CV_type = collapse_result(result = result_set, value = "cv_type")
  result_stats <- c(collapse_result(result = result_set, value = "sample_Z_score"),
                    collapse_result(result = result_set, value = "cv"),
                    collapse_result(result = result_set, value = "shapiro_P_value"))

  
  additional_result_stats <- c(collapse_result(result = result_set, value = "additional_sample_Z_score"),
                               collapse_result(result = result_set, value = "additional_cv"),
                               collapse_result(result = result_set, value = "additional_shapiro"))
  
  if(CV_type == theoretical){
    result <- c(result_stats, additional_result_stats)
  }
  else{
    result <- c(additional_result_stats, result_stats)
  }
  return(result)
}
rownames_additional_stats <- function(type){
  c(paste(type, zscore, sep="_"), type, paste(type, shapiro, sep="_"))
}