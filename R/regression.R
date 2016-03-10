SelectModelsFixedRegressionApproach <- function(nipt_sample, chromosomal_frac_control, chromosomes_trisomy, 
                                                frac_reads_chr_trisomy_observed, n_models, n_predictors){
  if(is.null(attr(nipt_sample, "class"))){
    print("")
  }
  else UseMethod("SelectModelsFixedRegressionApproach", nipt_sample)
}

GetNextPredictor <- function(samples, frac_reads_chr_trisomy_observed, predictors, chromosomal_frac_control){
  adj.r.squares <- NULL
  chr.candidates <- rownames(samples)
  for (i in 1:length(rownames(samples))){
    current.pred.set <- c(predictors, chr.candidates[i])
    chr.frac.df <- t(as.matrix(rbind(chromosomal_frac_control[current.pred.set, ])))
    model <- lm(frac_reads_chr_trisomy_observed ~ chr.frac.df)
    adj.r.squares[i] <- summary(model)$adj.r.squared
  }
  adj.r.squared.ordered <- order(adj.r.squares, decreasing = TRUE)
  chromosomes <- rownames(samples)
  best.pred.set <- c(predictors, chromosomes[adj.r.squared.ordered[1]])
  best.chr.frac.df <- t(as.matrix(rbind(chromosomal_frac_control[best.pred.set,])))
  best.model <- list(lm(frac_reads_chr_trisomy_observed ~ best.chr.frac.df))
  predictor.adj.r.squared <- list(chromosomes[adj.r.squared.ordered[1]], best.model)
  return(predictor.adj.r.squared)
}
#'Regression based Z score
#'
#'Make multiple models using linear regression and calculate Z-score
#'
#'
#'
#'@param nipt_sample The NIPTSample object that is the focus of the analysis
#'@param nipt_control_group The NIPTControlGroup object used in the analysis
#'@param chromo_focus The chromosome of interest. Most commonly chromosome 13, 18 or 21.
#' However, every autosomal chromosome can be predicted
#'@param n_models Integer Number of linear models to be made. Default setting is 4 models
#'@param n_predictors Integer The number of predictors each model contains. Default is 4
#'@param exclude_chromosomes integer. Exclude which autosomal chromosomes as potential predictors? 
#'Default potential trisomic chromosomes 13, 18 and 21 are exluded.
#'@param include_chromosomes integer. Include potential trisomic chromosomes? 
#'Options are: chromosomes 13, 18 and 21
#'@param use_test_train_set  Use a test and train set to build the models? Default is TRUE
#'@param size_of_train_set The size of the train set expressed in a decimal. 
#'Default is 0.6 (60 of the control samples)
#'@param overdispersion_rate The standard error of the mean is multiplied by this factor 
#'@param force_practical_cv Boolean, Ignore the theoretical CV and always use the practical CV?
#'
#'@details The regression based Z-score builds \emph{n} models with \emph{m} predictors using stepwise regression 
#'with forward selection. The models are used to predict the chromosomal fraction of interest, for the 
#'sample and for the control group. The observed fractions are then divided by the expected fraction, 
#'and Z-scores are calculated over these ratios. The Z-score is calculated by subtracting one from the 
#'ratio of the sample and dividing this result by the coefficient of variation. 
#'The coefficient of variation (CV) can either be the Practical or Theoretical CV. 
#'The Theoretical CV is the standard error multiplied by the overdispersion. 
#'Theoretically, the CV cannot be lower than the standard error of the mean. 
#'If it is case the CV is lower than Theoretical CV, then the Theoretical CV is used.  
#'
#'
#'The output of this function is an object of type RegressionResult, a named list containing: 
#'\itemize{
#'\item \strong{prediction_statistics} A dataframe with 7 rows and a column for every model. The rows are: 
#'\itemize{
#'\item \strong{Z_score_sample} The regression based Z score for the model
#'\item \strong{CV} The coefficient of varation for the model
#'\item \strong{cv_types} The CV type used to calculate the regression based Z score for the model. Either 
#'\emph{Practical_CV} or \emph{Theoretical_CV}
#'\item \strong{P_value_shapiro} The P value of the Shaipro-Wilk test for normality of the control group
#'regression based Z scores for the model
#'\item \strong{Predictor_chromosomes} The predictor chromosomes used in the model
#'\item \strong{Mean_test_set} The mean of the test set. Note that for calculating the regression based
#'Z scores the mean is replaced by one. The mean, however, can be seen as a quality metric for the model 
#'\item \strong{CV_train_set} The CV of the train set. The difference between this CV and the CV of 
#'the test can be used as a measure to quantify overfit
#'}
#'\item \strong{control_group_Zscores}  A matrix containing the regression based Z-scores for the control sample
#'\item \strong{focus_chromosome} he chromosome of interest. Most commonly chromosome 13, 18 or 21. 
#'However, every autosomal chromosome can be predicted
#'\item \strong{correction_status} The correction status of the control group autosomes
#'\item \strong{control_group_sample_names} The sample names of the test set group
#'\item \strong{models} List of the summary.lm output for every model
#'\item \strong{potential_predictors} The total pool of chromosomes where the predictors are selected from
#'\item \strong{all_control_group_Z_scores} Z-scores for every sample using theoretical and practical VCs
#'\item \strong{additional_statistics} Statistics for both the practical and theoretical CVs for every prediction set
#'}
#'
#'@return RegressionResult object
#'
#'@examples 
#' \dontrun{
#' regression_score_21 <- perform_regression(nipt_sample = sample_of_interest, 
#'                        nipt_control_group = control_group, chromo_focus = 21)
#' }
#' 
#'@export
perform_regression <- function(nipt_sample, nipt_control_group, chromo_focus, n_models = 4, n_predictors = 4,
                               exclude_chromosomes = NULL, include_chromosomes = NULL,
                               use_test_train_set =T, size_of_train_set = 0.6, overdispersion_rate = 1.15,
                               force_practical_cv = F){ 
 if (use_test_train_set == T){
   control_sets <- testtrainset(nipt_control_group = nipt_control_group, size_of_train_set = size_of_train_set)
   control_group_train <- control_sets$train_set
   nipt_control_group <- control_sets$test
   chromosomal_frac_control_train <- sapply(X = control_group_train[[samples]], FUN = chrfractions)
   chromosomal_frac_control_train <- setrownamesmatrix(chromosomal_frac_control_train)
  }
  else{
    chromosomal_frac_control_test <- sapply(X = nipt_control_group[[samples]], FUN = chrfractions)
    chromosomal_frac_control_test <- setrownamesmatrix(chromosomal_frac_control_test)
    chromosomal_frac_control_train <- chromosomal_frac_control_test
    control_group_train <- nipt_control_group
  }
  chromosomal_frac_control_test <- sapply(X = nipt_control_group[[samples]], FUN = chrfractions)
  chromosomal_frac_control_test <- setrownamesmatrix(chromosomal_frac_control_test)
 
  chromosomes_trisomy <- unique(c(chromosomes_trisomy, exclude_chromosomes, chromo_focus))
 
  if (!is.null(include_chromosomes)){
    chromosomes_trisomy <- chromosomes_trisomy[!chromosomes_trisomy %in% include_chromosomes]
  }
  potential_predictorset <-  getcontrolchromosomes(nipt_sample = nipt_sample, 
                                                   control_chromosomes = autosomal_chromosomes[!autosomal_chromosomes %in% chromosomes_trisomy])
  frac_reads_chr_trisomy_observed_train <- retrieve_fractions_of_interest(nipt_sample = nipt_sample, chromo_focus = chromo_focus, 
                                                                    chromosomal_fracs = chromosomal_frac_control_train)
  frac_reads_chr_trisomy_observed_test <- retrieve_fractions_of_interest(nipt_sample = nipt_sample, chromo_focus = chromo_focus, 
                                                                     chromosomal_fracs = chromosomal_frac_control_test)
  predictor.list <- SelectModelsFixedRegressionApproach(nipt_sample = nipt_sample, chromosomal_frac_control= chromosomal_frac_control_train,
                                                        chromosomes_trisomy = chromosomes_trisomy, 
                                                        frac_reads_chr_trisomy_observed = frac_reads_chr_trisomy_observed_train, 
                                                        n_models = n_models, n_predictors = n_predictors)
  prediction <- lapply(X = predictor.list, FUN = PredictTrisomy, nipt_sample = nipt_sample, chromo_focus = chromo_focus, control_group = nipt_control_group, 
                       frac_reads_chr_trisomy_observed_train = frac_reads_chr_trisomy_observed_train, 
                       frac_reads_chr_trisomy_observed = frac_reads_chr_trisomy_observed_test,
                       control_group_train = control_group_train, overdispersion_rate = overdispersion_rate,
                       force_practical_cv = force_practical_cv)
  result_nipt <- collapse_results(result_set = prediction, n_models = n_models, n_predictors = n_predictors)
  #if (use_test_train_set == F){
  end_result <- regression_template(result_set = result_nipt, chromo_focus = chromo_focus, 
                                    correction_status = unique(sapply(nipt_control_group[[samples]], getcorrectionstatus, status_type="correction_status_autosomal_chromosomes")),
                                    potential_predictors = potential_predictorset, samplenames = sapply(nipt_control_group[[samples]], getsamplenames),
                                    models = lapply(prediction, listmodels), type = class(nipt_sample)[2])
#   }
#   else{
#     end_result <- regression_template(result_set = result_nipt, chromo_focus = chromo_focus, 
#                                       correction_status = unique(sapply(nipt_control_group$Samples, getcorrectionstatus)),
#                                       potential_predictors = NULL, samplenames = sapply(nipt_control_group$Samples, getsamplenames),
#                                       models = lapply(prediction, listmodels), type = getstrandtype(nipt_sample = nipt_sample),
#                                       sample_names_train_set = sapply(nipt_control_group_train$Samples, getsamplenames),
#                                       train_set_statistics = train_set_statistics,
#                                       train_set_Zscores = train_set_Zscores)
#   }
  return (end_result)
}

BuildFullModel <- function(control.samples, chr.of.interest.fractions){
  return (lm(chr.of.interest.fractions ~ ., data=control.samples))
}
SelectModelsFixedRegressionApproach.SeparatedStrands <- function(nipt_sample, chromosomal_frac_control, chromosomes_trisomy, 
                                                                 frac_reads_chr_trisomy_observed, n_models, n_predictors){
  predictor.list <- list()
  chr.potential.trisomic <- c(paste0(chromosomes_trisomy, "F"),paste0(chromosomes_trisomy, "R"))
  for (model in 1:n_models)  {
    predictors <- NULL
    predictor.adj.r.squares <-NULL
    predictors.complementary <- list()
    adj.r.squares <- NULL
    for (predictor in 1:n_predictors){
      potential.predictors <- chromosomal_frac_control[rownames(chromosomal_frac_control)[!rownames(chromosomal_frac_control) %in% c(chr.potential.trisomic,  unlist(predictors.complementary), unlist(predictor.list))],]
      
      predictor.adj.r.squares <- GetNextPredictor(samples = potential.predictors, frac_reads_chr_trisomy_observed=frac_reads_chr_trisomy_observed, predictors = predictors, 
                                                  chromosomal_frac_control=chromosomal_frac_control)
      predictors[predictor] <- predictor.adj.r.squares[[1]]
      predictors.complementary[[predictor]] <- unique(c(predictor.adj.r.squares[[1]], sub("F", "R", predictor.adj.r.squares[[1]]), sub("R", "F", predictor.adj.r.squares[[1]])))
      adj.r.squares[predictor] <- predictor.adj.r.squares[[2]][1]
    }
    class(predictors) <- c("Prediction Set", "SeparatedStrands")
    predictor.list[[model]] <- predictors
  }
  return(predictor.list)
}

SelectModelsFixedRegressionApproach.CombinedStrands <- function(nipt_sample, chromosomal_frac_control, chromosomes_trisomy,
                                                                frac_reads_chr_trisomy_observed, predictor.max, n_models,
                                                                n_predictors){
  predictor.list <- list()
  chr.potential.trisomic <- chromosomes_trisomy
  
  for (model in 1:n_models){
    predictors <- NULL
    predictor.adj.r.squares <-NULL
    adj.r.squares <- list()
    for (predictor in 1:n_predictors){
      potential.predictors <- chromosomal_frac_control[rownames(chromosomal_frac_control)[!rownames(chromosomal_frac_control) %in% c(chr.potential.trisomic, unlist(predictor.list), predictors)],]
      predictor.adj.r.squares <- GetNextPredictor(samples = potential.predictors, frac_reads_chr_trisomy_observed=frac_reads_chr_trisomy_observed, predictors = predictors, 
                                                  chromosomal_frac_control=chromosomal_frac_control)
      predictors[predictor] <- predictor.adj.r.squares[[1]]
      
      adj.r.squares[predictor] <- predictor.adj.r.squares[[2]]
    }
    class(predictors) <- c("Prediction Set", "CombinedStrands")
    predictor.list[[model]] <- predictors
  }
  return(predictor.list)
}

PredictTrisomy <- function(predictors, nipt_sample,  chromo_focus,  control_group, control_group_train,
                           frac_reads_chr_trisomy_observed, frac_reads_chr_trisomy_observed_train,
                           overdispersion_rate, force_practical_cv){
  chromosomal_frac_control <- sapply(X = control_group[[samples]], FUN = chrfractions)
  chromosomal_frac_control <- setrownamesmatrix(chromosomal_frac_control)
  chromosomal_frac_control_train <- sapply(X = control_group_train[[samples]], FUN = chrfractions)
  chromosomal_frac_control_train <- setrownamesmatrix(chromosomal_frac_control_train)
  chromosomal_frac_sample <- sapply(list(nipt_sample), chrfractions)
  chromosomal_frac_sample <- setrownamesmatrix(chromosomal_frac_sample)
  samplereads <- sumfandrautosomal(nipt_sample)
  cv.theo <- overdispersion_rate * ( 1/ sqrt(sum(samplereads[chromo_focus,])))
 
  mod <- BuildFullModel(as.data.frame(t(chromosomal_frac_control_train))[predictors], chr.of.interest.fractions =  frac_reads_chr_trisomy_observed_train)
  
  controls_predicted <- predict(mod, as.data.frame(t(chromosomal_frac_control))[predictors])

  train_set_predicted <- mod$fitted.values
    
  ratio <-  frac_reads_chr_trisomy_observed /   controls_predicted 
  
  ratio_train <- mod$fitted.values / frac_reads_chr_trisomy_observed_train
  
  cv_train <- stats::sd(ratio_train) / mean(ratio_train)
  
  cv.prac <- stats::sd(ratio) / mean(ratio)
 
  cvs <- c(cv.theo, cv.prac)
  names(cvs) <- c(theoretical, practical)
  if (force_practical_cv == T){
    cv <- cvs[practical]
  }
  else{
    cv <- cvs[which.max(cvs)]
  }
  sample.predicted <- predict(mod, as.data.frame(t(chromosomal_frac_sample))[predictors])
  
  sample_ratio <- retrieve_fractions_of_interest(nipt_sample = nipt_sample, chromo_focus = chromo_focus, 
                                                 chromosomal_fracs = as.matrix(chromosomal_frac_sample)) / sample.predicted  
  sample_names <- sapply(X = control_group[[samples]], FUN = getsamplenames)
  scores_and_statistics <- lapply(cvs, calculate_regression_scores, sample_ratio=sample_ratio,
                                  ratio = ratio, sample_names = sample_names)
  names(scores_and_statistics) <- c(theoretical, practical)
  result_score <- scores_and_statistics[[names(cv)]]
  additional_result <- scores_and_statistics[[-which(names(scores_and_statistics) == names(cv))]]
  return(list(sample_Z_score = result_score$score_sample, control_group_Z_scores = result_score$score_controls,
              shapiro_P_value = result_score$shapiro, mean_test_set = mean(ratio), cv_train_set = cv_train,
              predictors = predictors, summary_model = summary(mod), cv = cv, cv_type=names(cv),
              additional_sample_Z_score = additional_result$score_sample,
              additional_control_group_Z_scores = additional_result$score_controls,
              additional_shapiro = additional_result$shapiro,
              additional_cv = cvs[[-which(names(cvs) == names(cv))]]))
}

calculate_regression_scores <- function(cv, sample_ratio, ratio, sample_names ){
  z.score.sample <- (sample_ratio - 1) / cv
  z.score.controls <- matrix((data = as.numeric(ratio) - 1) / cv, ncol = 1, dimnames = list(sample_names, NULL))
  shapiro.regression <- shapiro.test(z.score.controls)
  return(list(score_sample = z.score.sample, score_controls = z.score.controls, shapiro = shapiro.regression$p.value))
}

