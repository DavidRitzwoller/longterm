mu_s_x <- function(data_train, treatment_val, S_vars, X_vars, Y_var, type,
                   nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                   xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads){
  # (S, X, G=1) conditional outcome mean
  df <- filter(data_train, observe == 1, treatment == treatment_val)
  s_x <- select(df, one_of(S_vars, X_vars))
  y <- select(df, one_of(Y_var))
  if(type == "glmnet"){
    mu_s_x_val <- cv.glmnet(x = as.matrix(s_x), y = as.matrix(y), nfolds = nuisance_cv_fold)
  } else if(type == "grf"){
    mu_s_x_val <- regression_forest(X = as.matrix(s_x), Y = as.matrix(y), num.threads = grf_num_threads,
                                    ci.group.size = 1, honesty = grf_honesty, tune.parameters = grf_tune_parameters)
  } else if(type == "xgboost"){
    cv <- xgb.cv(data = as.matrix(s_x), label = as.matrix(y), 
                 max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                 nrounds = xgb_cv_rounds, verbose = FALSE, nfold = nuisance_cv_fold)
    mu_s_x_val <- xgboost(data = as.matrix(s_x), label = as.matrix(y), 
                          max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                          nrounds = which.min(cv$evaluation_log$test_rmse_mean), verbose = FALSE)
  } else {
    stop('Enter a valid nuisance parameter estimation type')
  }
  return(mu_s_x_val)
}

nu_s_x <- function(data_train, S_vars, X_vars, Y_var, type,
                   nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                   xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads){
  # (S, X, G=1) conditional outcome mean
  df <- filter(data_train, observe == 1)
  s_x <- select(df, one_of(S_vars, X_vars))
  y <- select(df, one_of(Y_var))
  if(type == "glmnet"){
    nu_s_x_val <- cv.glmnet(x = as.matrix(s_x), y = as.matrix(y), nfolds = nuisance_cv_fold)
  } else if(type == "grf"){
    nu_s_x_val <- regression_forest(X = as.matrix(s_x), Y = as.matrix(y), num.threads = grf_num_threads,
                                    ci.group.size = 1, honesty = grf_honesty, tune.parameters = grf_tune_parameters)
  } else if(type == "xgboost"){
    cv <- xgb.cv(data = as.matrix(s_x), label = as.matrix(y), 
                 max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                 nrounds = xgb_cv_rounds, verbose = FALSE, nfold = nuisance_cv_fold)
    nu_s_x_val <- xgboost(data = as.matrix(s_x), label = as.matrix(y), 
                          max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                          nrounds = which.min(cv$evaluation_log$test_rmse_mean), verbose = FALSE)
  } else {
    stop('Enter a valid nuisance parameter estimation type')
  }
  return(nu_s_x_val)
}

bar_mu_x <- function(data_train, treatment_val, S_vars, X_vars, Y_var, mu_s_x, type,
                     nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                     xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads){
  # (X, G=0) contional outcome mean
  df <- filter(data_train, observe == 0, treatment == treatment_val)
  s_x <- select(df, one_of(S_vars, X_vars))
  x <- select(df, one_of(X_vars))
  if(type == "glmnet"){
    pred_mu <- predict(mu_s_x, newx = as.matrix(s_x), s = mu_s_x$lambda.min)
    if(length(unique(pred_mu)) <= 5){
      bar_mu_x_val <- mean(pred_mu)
      constant_flag_val <- 1
    } else {
      bar_mu_x_val <- cv.glmnet(x = as.matrix(x), y = pred_mu, nfolds = nuisance_cv_fold)
      constant_flag_val <- 0
    }
  } else if(type == "grf"){
    pred_mu <- predict(mu_s_x, as.matrix(s_x))
    if(length(unique(pred_mu)) == 1){
      bar_mu_x_val <- mean(pred_mu$predictions)
      constant_flag_val <- 1
    } else {
      bar_mu_x_val <- regression_forest(X = as.matrix(x), Y = as.matrix(pred_mu), num.threads = grf_num_threads,
                                        ci.group.size = 1, honesty = grf_honesty, tune.parameters = grf_tune_parameters)
      constant_flag_val <- 0
    }
  } else if(type == "xgboost"){
    pred_mu <- predict(mu_s_x, as.matrix(s_x))
    if(length(unique(pred_mu)) == 1){
      bar_mu_x_val <- mean(pred_mu)
      constant_flag_val <- 1
    } else {
      cv <- xgb.cv(data = as.matrix(x), label = as.matrix(pred_mu), 
                   max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                   nrounds = xgb_cv_rounds, verbose = FALSE, nfold = nuisance_cv_fold)
      bar_mu_x_val <- xgboost(data = as.matrix(x), label = as.matrix(pred_mu), 
                              max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                              nrounds = which.min(cv$evaluation_log$test_rmse_mean), verbose = FALSE)
      constant_flag_val <- 0
    }
  } else {
    stop('Enter a valid nuisance parameter estimation type')
  }
  return(list(bar_mu_x = bar_mu_x_val, constant_flag = constant_flag_val))
}

bar_nu_x <- function(data_train, treatment_val, S_vars, X_vars, Y_var, nu_s_x, type,
                     nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                     xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads){
  # (X, G=0) contional outcome mean
  df <- filter(data_train, observe == 0, treatment == treatment_val)
  s_x <- select(df, one_of(S_vars, X_vars))
  x <- select(df, one_of(X_vars))
  if(type == "glmnet"){
    pred_nu <- predict(nu_s_x, newx = as.matrix(s_x), s = nu_s_x$lambda.min)
    if(length(unique(pred_nu)) <= 5){
      bar_nu_x_val <- mean(pred_nu)
      constant_flag_val <- 1
    } else {
      bar_nu_x_val <- cv.glmnet(x = as.matrix(x), y = pred_nu, nfolds = nuisance_cv_fold)
      constant_flag_val <- 0
    }
  } else if(type == "grf"){
    pred_nu <- predict(nu_s_x, as.matrix(s_x))
    if(length(unique(pred_nu)) == 1){
      bar_nu_x_val <- mean(pred_nu$predictions)
      constant_flag_val <- 1
    } else {
      bar_nu_x_val <- regression_forest(X = as.matrix(x), Y = as.matrix(pred_nu), num.threads = grf_num_threads,
                                        ci.group.size = 1, honesty = grf_honesty, tune.parameters = grf_tune_parameters)
      constant_flag_val <- 0
    }
  } else if(type == "xgboost"){
    pred_nu <- predict(nu_s_x, as.matrix(s_x))
    if(length(unique(pred_nu)) == 1){
      bar_nu_x_val <- mean(pred_nu)
      constant_flag_val <- 1
    } else {
      cv <- xgb.cv(data = as.matrix(x), label = as.matrix(pred_nu), 
                   max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                   nrounds = xgb_cv_rounds, verbose = FALSE, nfold = nuisance_cv_fold)
      bar_nu_x_val <- xgboost(data = as.matrix(x), label = as.matrix(pred_nu), 
                              max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                              nrounds = which.min(cv$evaluation_log$test_rmse_mean), verbose = FALSE)
      constant_flag_val <- 0
    }
  } else {
    stop('Enter a valid nuisance parameter estimation type')
  }
  return(list(bar_nu_x = bar_nu_x_val, constant_flag = constant_flag_val))
}

rho <- function(data_train, treatment_val, S_vars, X_vars, type,
                nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads){
  # Latent conditional treatment propensity
  # Learn first term
  df <- filter(data_train, treatment == treatment_val)
  s_x <- select(df, one_of(S_vars, X_vars))
  g <- select(df, observe)
  if(type == "glmnet"){
    term_1 <- cv.glmnet(x = as.matrix(s_x), y = as.matrix(g), family = "binomial", nfolds = nuisance_cv_fold)
  } else if(type == "grf"){
    term_1 <- regression_forest(X = as.matrix(s_x), Y = as.matrix(g), num.threads = grf_num_threads,
                                ci.group.size = 1, honesty = grf_honesty, tune.parameters = grf_tune_parameters)
  } else if(type == "xgboost"){
    cv <- xgb.cv(data = as.matrix(s_x), label = as.matrix(g), eval_metric = "logloss",
                 max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                 nrounds = xgb_cv_rounds, verbose = FALSE, nfold = nuisance_cv_fold)
    term_1 <- xgboost(data = as.matrix(s_x), label = as.matrix(g), eval_metric = "logloss",
                      max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                      nrounds = which.min(cv$evaluation_log$test_logloss_mean),
                      objective = "binary:logistic", verbose = FALSE)
  } else {
    stop('Enter a valid nuisance parameter estimation type')
  }
  # Learn second term
  df <- filter(data_train, observe == 1)
  x <- select(df, one_of(X_vars))
  w <- select(df, one_of("treatment"))
  if(type == "glmnet"){
    term_2 <- cv.glmnet(x = as.matrix(x), y = as.matrix(w), family = "binomial", nfolds = nuisance_cv_fold)
  } else if(type == "grf"){
    term_2 <- regression_forest(X = as.matrix(x), Y = as.matrix(w), num.threads = grf_num_threads,
                                ci.group.size = 1, honesty = grf_honesty, tune.parameters = grf_tune_parameters)
  } else if(type == "xgboost"){
    cv <- xgb.cv(data = as.matrix(x), label = as.matrix(w), eval_metric = "logloss",
                 max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                 nrounds = xgb_cv_rounds, verbose = FALSE, nfold = nuisance_cv_fold)
    term_2 <- xgboost(data = as.matrix(x), label = as.matrix(w), eval_metric = "logloss",
                      max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                      nrounds = which.min(cv$evaluation_log$test_logloss_mean),
                      objective = "binary:logistic", verbose = FALSE)
  } else {
    stop('Enter a valid nuisance parameter estimation type')
  }
  return(list(term_1 = term_1, term_2 = term_2))
}

varrho_x <- function(data_train, X_vars, type,
                     nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                     xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads){
  # Conditional treatment propensity
  df <- filter(data_train, observe == 0)
  x <- select(df, one_of(X_vars))
  w <- select(df, one_of("treatment"))
  if(type == "glmnet"){
    varrho_x <- cv.glmnet(x = as.matrix(x), y = as.matrix(w), family = "binomial", nfolds = nuisance_cv_fold)
  } else if(type == "grf"){
    varrho_x <- regression_forest(X = as.matrix(x), Y = as.matrix(w), num.threads = grf_num_threads,
                                  ci.group.size = 1, honesty = grf_honesty, tune.parameters = grf_tune_parameters)
  } else if(type == "xgboost"){
    cv <- xgb.cv(data = as.matrix(x), label = as.matrix(w), eval_metric = "logloss",
                 max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                 nrounds = xgb_cv_rounds, verbose = FALSE, nfold = nuisance_cv_fold)
    varrho_x <- xgboost(data = as.matrix(x), label = as.matrix(w), eval_metric = "logloss",
                        max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                        nrounds = which.min(cv$evaluation_log$test_logloss_mean),
                        objective = "binary:logistic", verbose = FALSE)
  } else {
    stop('Enter a valid nuisance parameter estimation type')
  }
  return(varrho_x)
}

varrho_s_x <- function(data_train, X_vars, S_vars, type,
                       nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                       xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads){
  # Conditional treatment propensity
  df <- filter(data_train, observe == 0)
  s_x <- select(df, one_of(S_vars, X_vars))
  w <- select(df, one_of("treatment"))
  if(type == "glmnet"){
    varrho_s_x <- cv.glmnet(x = as.matrix(s_x), y = as.matrix(w), family = "binomial", nfolds = nuisance_cv_fold)
  } else if(type == "grf"){
    varrho_s_x <- regression_forest(X = as.matrix(s_x), Y = as.matrix(w), num.threads = grf_num_threads,
                                  ci.group.size = 1, honesty = grf_honesty, tune.parameters = grf_tune_parameters)
  } else if(type == "xgboost"){
    cv <- xgb.cv(data = as.matrix(s_x), label = as.matrix(w), eval_metric = "logloss",
                 max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                 nrounds = xgb_cv_rounds, verbose = FALSE, nfold = nuisance_cv_fold)
    varrho_s_x <- xgboost(data = as.matrix(s_x), label = as.matrix(w), eval_metric = "logloss",
                          max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                          nrounds = which.min(cv$evaluation_log$test_logloss_mean),
                          objective = "binary:logistic", verbose = FALSE)
  } else {
    stop('Enter a valid nuisance parameter estimation type')
  }
  return(varrho_s_x)
}

gamma_x <- function(data_train, X_vars, type,
                    nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                    xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads){
  # Conditional observational sample propensity
  x <- select(data_train, one_of(X_vars))
  g <- select(data_train, observe)
  if(type == "glmnet"){
    gamma_x <- cv.glmnet(x = as.matrix(x), y = as.matrix(g), family = "binomial", nfolds = 5)
  } else if(type == "grf"){
    gamma_x <- regression_forest(X = as.matrix(x), Y = as.matrix(g), num.threads = grf_num_threads,
                                 ci.group.size = 1, honesty = grf_honesty, tune.parameters = grf_tune_parameters)
  } else if(type == "xgboost"){
    cv <- xgb.cv(data = as.matrix(x), label = as.matrix(g), eval_metric = "logloss",
                 max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                 nrounds = xgb_cv_rounds, verbose = FALSE, nfold = nuisance_cv_fold)
    gamma_x <- xgboost(data = as.matrix(x), label = as.matrix(g), eval_metric = "logloss",
                       max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                       nrounds = which.min(cv$evaluation_log$test_logloss_mean), 
                       objective = "binary:logistic", verbose = FALSE)
  } else {
    stop('Enter a valid nuisance parameter estimation type')
  }
  return(gamma_x)
}

gamma_s_x <- function(data_train, X_vars, S_vars, type,
                      nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                      xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads){
  # Conditional observational sample propensity
  s_x <- select(data_train, one_of(S_vars, X_vars))
  g <- select(data_train, observe)
  if(type == "glmnet"){
    gamma_s_x <- cv.glmnet(x = as.matrix(s_x), y = as.matrix(g), family = "binomial", nfolds = 5)
  } else if(type == "grf"){
    gamma_s_x <- regression_forest(X = as.matrix(s_x), Y = as.matrix(g), num.threads = grf_num_threads,
                                 ci.group.size = 1, honesty = grf_honesty, tune.parameters = grf_tune_parameters)
  } else if(type == "xgboost"){
    cv <- xgb.cv(data = as.matrix(s_x), label = as.matrix(g), eval_metric = "logloss",
                 max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                 nrounds = xgb_cv_rounds, verbose = FALSE, nfold = nuisance_cv_fold)
    gamma_s_x <- xgboost(data = as.matrix(s_x), label = as.matrix(g), eval_metric = "logloss",
                       max_depth = xgb_max_depth, eta = xgb_eta, nthread = xgb_threads,
                       nrounds = which.min(cv$evaluation_log$test_logloss_mean), 
                       objective = "binary:logistic", verbose = FALSE)
  } else {
    stop('Enter a valid nuisance parameter estimation type')
  }
  return(gamma_s_x)
}