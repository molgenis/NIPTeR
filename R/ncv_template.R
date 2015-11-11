ncv_template <- function(denominators, chromo_focus, nipt_sample_names, correction_status, scores, 
                         potential_denominators, statistics, type, sample_names_train_set = NULL, 
                         train_set_statistics = NULL, train_set_Zscores = NULL){
  if (is.null(train_set_statistics)){
    new_ncv_template <- list(denominators = denominators, focus_chromosome = as.character(chromo_focus), 
                             control_group_sample_names = nipt_sample_names, correction_status = correction_status,
                             control_group_Zscores = scores, 
                             potential_denominators = potential_denominators, control_group_statistics = statistics)
  }
  else{
    new_ncv_template <- list(denominators = denominators, focus_chromosome = as.character(chromo_focus), 
                             control_group_sample_names = nipt_sample_names, correction_status = correction_status,
                             control_group_Zscores = scores, 
                             potential_denominators = potential_denominators, control_group_statistics = statistics,
                             sample_names_train_set = sample_names_train_set, train_set_statistics = train_set_statistics,
                             train_set_Zscores = train_set_Zscores)
  }
  class(new_ncv_template) <- c(NCV_template_class, type)
  
  return(new_ncv_template)
}

