#' Semiparametric Estimation of Long-Term Treatment Effects
#'
#' \code{longterm} estimates the long-term average treatment effect of a binary treatment on a scalar long-term outcome
#' using the semiparametric estimator developed in Chen and Ritzwoller (2021).
#'
#' @param data A data frame containing the indicators \code{"treatment,"} denoting whether treatment was assigned, and \code{"observe,"} denoting whether the observation was from the observational or experimental sample, in addition to pretreatment covariates, short-term outcomes , and long-term outcomes. If an observation is missing values of a particular variable by construction, e.g., long-term outcomes are not observed in the experimental sample, then these values should be coded with any nonmissing value and will not contribute to estimation.
#' @param S_vars A list containing strings of the names of the short-term outcome variables.
#' @param X_vars A list containing strings of the names of the pre-treatment covariates.
#' @param Y_var A string giving the name of the long-term outcome of interest.
#' @param obs A boolean variable specifying whether treatment is observed in the long-term sample.
#' @param estimand A boolean variable specifying whether the estimand of interest is the long-term average treatment effect in the observational population, as opposed to the experimental population.
#' @param type A string specifying how nuisance functions should be estimated. \code{"glmnet"} specifies cross-validated lasso. \code{"grf"} specifies generalized random forests. \code{"xgboost"} specifies XGBoost.
#' @param prop_lb A float specifying the lower threshold for propensity score estimates.
#' @param prop_ub A float specifying the upper threshold for propensity score estimates.
#' @param alpha One minus the nominal coverage probability of the confidence intervals.
#' @param te_lb A float giving a lower bound for the long-term average treatment effect.
#' @param te_ub A float giving an upper bound for the long-term average treatment effect.
#' @param cross_fit_fold An integer giving the number of folds for the cross-fit estimation of nuisance parameters.
#' @param nuisance_cv_fold An interger giving the number of folds for cross-validating nuisance parameter estimates.
#' @param grf_honesty A boolean variable setting the \code{honesty} parameter for \code{grf} estimation.
#' @param grf_tune_parameters A string variable setting the \code{tune.parameters} parameter for \code{grf} estimation.
#' @param grf_num_threads An integer variable setting the \code{num.threads} parameter for \code{grf} estimation.
#' @param xgb_cv_rounds An integer variable setting maximum number of rounds for \code{xgboost} estimation.
#' @param xgb_eta A float variable setting \code{eta} parameter for \code{xgboost} estimation.
#' @param xgb_max_depth A float variable setting \code{max_depth} parameter for \code{xgboost} estimation.
#' @param xgb_threads An integer variable setting the \code{nthread} parameter for \code{xgboost} estimation.
#'
#' @return Returns a list with three components: \describe{
#'
#' \item{\code{hat_tau}}{Estimate of the long-term average treatment effect.}
#'
#' \item{\code{se}}{Estimate of the standard error of the estimator.}
#'
#' \item{\code{ci}}{A vector giving the lower and upper points of the confidence interval.}
#'
#' }
#' @export
#'
#' @references{
#' \cite{Chen, J., & Ritzwoller, D. M. (2021).
#' Semiparametric Estimation of Long-Term Treatment Effects.
#' arXiv preprint arXiv:2107.14405.}
#' }
#'

longterm <- function(data, S_vars, X_vars, Y_var, obs, estimand, type,
                      prop_lb = 0.01, prop_ub = 0.99, alpha = 0.05,
                      te_lb = -500, te_ub = 500,
                      cross_fit_fold = 5, nuisance_cv_fold = 5,
                      grf_honesty = FALSE, grf_tune_parameters = "all", grf_num_threads = 1,
                      xgb_cv_rounds = 100, xgb_eta = 0.1, xgb_max_depth = 2, xgb_threads = 1){
  ### Estimate Long-Term Treatment Effects
  # Choose random k-fold of the data
  folds <- createFolds(seq(nrow(data)), k = cross_fit_fold)
  # Compute the nuisance parameters
  data_nuisance_out <- lapply(seq(cross_fit_fold),
                              compute_nuisance_on_fold, data,
                              folds, S_vars, X_vars, Y_var, obs, type,
                              prop_lb, prop_ub, cross_fit_fold, nuisance_cv_fold,
                              grf_honesty, grf_tune_parameters, grf_num_threads,
                              xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
  data_nuisance     <- bind_rows(data_nuisance_out)
  # Compute the estimate of tau_1
  solution  <- optimize(function(tau) objective(tau, data_nuisance, Y_var, obs, estimand), upper = te_ub, lower = te_lb)
  hat_tau   <- solution$minimum
  # Compute the standard error and construct confidence intervals
  se_ci <- compute_ci(hat_tau, data_nuisance, Y_var, obs, estimand, alpha)
  se    <- se_ci$se
  ci    <- se_ci$ci
  return(list(hat_tau = hat_tau, se = se, ci = ci))
}

compute_nuisance_on_fold <- function(i, data, folds, S_vars, X_vars, Y_var, obs, type,
                                     prop_lb, prop_ub, cross_fit_fold, nuisance_cv_fold,
                                     grf_honesty, grf_tune_parameters, grf_num_threads,
                                     xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads){
  ## Subset to training and testing samples
  data_train <- data[as.numeric(unlist(folds[setdiff(seq(cross_fit_fold), i)])), ]
  data_test  <- data[as.numeric(unlist(folds[i])), ]

  ## Train nuisance parameters
  nuisance <- train_nuisance(data_train, S_vars, X_vars, Y_var, type, obs,
                             nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                             xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)

  ## Compute nuisance terms on test data and compile output into data frame
  data_nuisance <- data_test %>% select(one_of(Y_var, "treatment", "observe"))
  x             <- select(data_test, one_of(X_vars))
  s_x           <- select(data_test, one_of(S_vars, X_vars))
  if(obs == TRUE){
    # Treatment Observed in the Observational Data Set
    if(type == "glmnet"){
      if(nuisance$bar_mu_x_1_cf == 0){data_nuisance$bar_mu_x_1 <- predict(nuisance$bar_mu_x_1, newx = as.matrix(x), s = nuisance$bar_mu_x_1$lambda.min)
      } else {data_nuisance$bar_mu_x_1 <- nuisance$bar_mu_x_1}
      if(nuisance$bar_mu_x_0_cf == 0){data_nuisance$bar_mu_x_0 <- predict(nuisance$bar_mu_x_0, newx = as.matrix(x), s = nuisance$bar_mu_x_0$lambda.min)
      } else {data_nuisance$bar_mu_x_0 <- nuisance$bar_mu_x_0}
      data_nuisance$mu_s_x_1     <- predict(nuisance$mu_s_x_1,               newx = as.matrix(s_x), s = nuisance$mu_s_x_1$lambda.min)
      data_nuisance$mu_s_x_0     <- predict(nuisance$mu_s_x_0,               newx = as.matrix(s_x), s = nuisance$mu_s_x_0$lambda.min)
      data_nuisance$varrho_x     <- pmax(pmin(predict(nuisance$varrho_x,     newx = as.matrix(x),   s = nuisance$varrho_x$lambda.min,     type = "response"), prop_ub), prop_lb)
      data_nuisance$gamma_x      <- pmax(pmin(predict(nuisance$gamma_x,      newx = as.matrix(x),   s = nuisance$gamma_x$lambda.min,      type = "response"), prop_ub), prop_lb)
      term_1_1                   <- pmax(pmin(predict(nuisance$rho_1$term_1, newx = as.matrix(s_x), s = nuisance$rho_1$term_1$lambda.min, type = "response"), prop_ub), prop_lb)
      term_2_1                   <- pmax(pmin(predict(nuisance$rho_1$term_2, newx = as.matrix(x),   s = nuisance$rho_1$term_2$lambda.min, type = "response"), prop_ub), prop_lb)
      term_1_0                   <- pmax(pmin(predict(nuisance$rho_0$term_1, newx = as.matrix(s_x), s = nuisance$rho_0$term_1$lambda.min, type = "response"), prop_ub), prop_lb)
      term_2_0                   <- pmax(pmin(predict(nuisance$rho_0$term_2, newx = as.matrix(x),   s = nuisance$rho_0$term_2$lambda.min, type = "response"), prop_ub), prop_lb)
    } else if(type == "grf"){
      if(nuisance$bar_mu_x_1_cf == 0){ data_nuisance$bar_mu_x_1 <- as.matrix(predict(nuisance$bar_mu_x_1, newdata = as.matrix(x)))
      } else {data_nuisance$bar_mu_x_1 <- nuisance$bar_mu_x_1}
      if(nuisance$bar_mu_x_0_cf == 0){ data_nuisance$bar_mu_x_0 <- as.matrix(predict(nuisance$bar_mu_x_0, newdata = as.matrix(x)))
      } else {data_nuisance$bar_mu_x_0 <- nuisance$bar_mu_x_0}
      data_nuisance$mu_s_x_1   <- as.matrix(predict(nuisance$mu_s_x_1,               newdata = as.matrix(s_x)))
      data_nuisance$mu_s_x_0   <- as.matrix(predict(nuisance$mu_s_x_0,               newdata = as.matrix(s_x)))
      data_nuisance$varrho_x   <- pmax(pmin(as.matrix(predict(nuisance$varrho_x,     newdata = as.matrix(x))),   prop_ub), prop_lb)
      data_nuisance$gamma_x    <- pmax(pmin(as.matrix(predict(nuisance$gamma_x,      newdata = as.matrix(x))),   prop_ub), prop_lb)
      term_1_1                 <- pmax(pmin(as.matrix(predict(nuisance$rho_1$term_1, newdata = as.matrix(s_x))), prop_ub), prop_lb)
      term_2_1                 <- pmax(pmin(as.matrix(predict(nuisance$rho_1$term_2, newdata = as.matrix(x))),   prop_ub), prop_lb)
      term_1_0                 <- pmax(pmin(as.matrix(predict(nuisance$rho_0$term_1, newdata = as.matrix(s_x))), prop_ub), prop_lb)
      term_2_0                 <- pmax(pmin(as.matrix(predict(nuisance$rho_0$term_2, newdata = as.matrix(x))),   prop_ub), prop_lb)
    } else if(type == "xgboost"){
      if(nuisance$bar_mu_x_1_cf == 0){ data_nuisance$bar_mu_x_1 <- predict(nuisance$bar_mu_x_1, newdata = as.matrix(x))
      } else {data_nuisance$bar_mu_x_1 <- nuisance$bar_mu_x_1}
      if(nuisance$bar_mu_x_0_cf == 0){ data_nuisance$bar_mu_x_0 <- predict(nuisance$bar_mu_x_0, newdata = as.matrix(x))
      } else {data_nuisance$bar_mu_x_0 <- nuisance$bar_mu_x_0}
      data_nuisance$mu_s_x_1   <- predict(nuisance$mu_s_x_1,               newdata = as.matrix(s_x))
      data_nuisance$mu_s_x_0   <- predict(nuisance$mu_s_x_0,               newdata = as.matrix(s_x))
      data_nuisance$varrho_x   <- pmax(pmin(predict(nuisance$varrho_x,     newdata = as.matrix(x)),   prop_ub), prop_lb)
      data_nuisance$gamma_x    <- pmax(pmin(predict(nuisance$gamma_x,      newdata = as.matrix(x)),   prop_ub), prop_lb)
      term_1_1                 <- pmax(pmin(predict(nuisance$rho_1$term_1, newdata = as.matrix(s_x)), prop_ub), prop_lb)
      term_2_1                 <- pmax(pmin(predict(nuisance$rho_1$term_2, newdata = as.matrix(x)),   prop_ub), prop_lb)
      term_1_0                 <- pmax(pmin(predict(nuisance$rho_0$term_1, newdata = as.matrix(s_x)), prop_ub), prop_lb)
      term_2_0                 <- pmax(pmin(predict(nuisance$rho_0$term_2, newdata = as.matrix(x)),   prop_ub), prop_lb)
    }
    data_nuisance$rho_1 <- pmax(pmin((term_1_1 / (1 - term_1_1)) * term_2_1, prop_ub), prop_lb)
    data_nuisance$rho_0 <- pmax(pmin((term_1_0 / (1 - term_1_0)) * term_2_0, prop_ub), prop_lb)
  } else {
    # Treatment Not Observed in the Observational Data Set
    if(type == "glmnet"){
      if(nuisance$bar_nu_x_1_cf == 0){data_nuisance$bar_nu_x_1 <- predict(nuisance$bar_nu_x_1, newx = as.matrix(x), s = nuisance$bar_nu_x_1$lambda.min)
      } else {data_nuisance$bar_nu_x_1 <- nuisance$bar_nu_x_1}
      if(nuisance$bar_nu_x_0_cf == 0){data_nuisance$bar_nu_x_0 <- predict(nuisance$bar_nu_x_0, newx = as.matrix(x), s = nuisance$bar_nu_x_0$lambda.min)
      } else {data_nuisance$bar_nu_x_0 <- nuisance$bar_nu_x_0}
      data_nuisance$nu_s_x       <- predict(nuisance$nu_s_x,                 newx = as.matrix(s_x), s = nuisance$nu_s_x$lambda.min)
      data_nuisance$varrho_x     <- pmax(pmin(predict(nuisance$varrho_x,     newx = as.matrix(x),   s = nuisance$varrho_x$lambda.min,     type = "response"), prop_ub), prop_lb)
      data_nuisance$varrho_s_x   <- pmax(pmin(predict(nuisance$varrho_s_x,   newx = as.matrix(s_x), s = nuisance$varrho_s_x$lambda.min,   type = "response"), prop_ub), prop_lb)
      data_nuisance$gamma_x      <- pmax(pmin(predict(nuisance$gamma_x,      newx = as.matrix(x),   s = nuisance$gamma_x$lambda.min,      type = "response"), prop_ub), prop_lb)
      data_nuisance$gamma_s_x    <- pmax(pmin(predict(nuisance$gamma_s_x,    newx = as.matrix(s_x), s = nuisance$gamma_s_x$lambda.min,    type = "response"), prop_ub), prop_lb)
    } else if(type == "grf"){
      if(nuisance$bar_nu_x_1_cf == 0){ data_nuisance$bar_nu_x_1 <- as.matrix(predict(nuisance$bar_nu_x_1, newdata = as.matrix(x)))
      } else {data_nuisance$bar_nu_x_1 <- nuisance$bar_nu_x_1}
      if(nuisance$bar_nu_x_0_cf == 0){ data_nuisance$bar_nu_x_0 <- as.matrix(predict(nuisance$bar_nu_x_0, newdata = as.matrix(x)))
      } else {data_nuisance$bar_nu_x_0 <- nuisance$bar_nu_x_0}
      data_nuisance$nu_s_x     <- as.matrix(predict(nuisance$nu_s_x,                 newdata = as.matrix(s_x)))
      data_nuisance$varrho_x   <- pmax(pmin(as.matrix(predict(nuisance$varrho_x,     newdata = as.matrix(x))),   prop_ub), prop_lb)
      data_nuisance$varrho_s_x <- pmax(pmin(as.matrix(predict(nuisance$varrho_s_x,   newdata = as.matrix(s_x))), prop_ub), prop_lb)
      data_nuisance$gamma_x    <- pmax(pmin(as.matrix(predict(nuisance$gamma_x,      newdata = as.matrix(x))),   prop_ub), prop_lb)
      data_nuisance$gamma_s_x  <- pmax(pmin(as.matrix(predict(nuisance$gamma_s_x,    newdata = as.matrix(s_x))),  prop_ub), prop_lb)
    } else if(type == "xgboost"){
      if(nuisance$bar_nu_x_1_cf == 0){ data_nuisance$bar_nu_x_1 <- predict(nuisance$bar_nu_x_1, newdata = as.matrix(x))
      } else {data_nuisance$bar_nu_x_1 <- nuisance$bar_nu_x_1}
      if(nuisance$bar_nu_x_0_cf == 0){ data_nuisance$bar_nu_x_0 <- predict(nuisance$bar_nu_x_0, newdata = as.matrix(x))
      } else {data_nuisance$bar_nu_x_0 <- nuisance$bar_nu_x_0}
      data_nuisance$nu_s_x     <- predict(nuisance$nu_s_x,                 newdata = as.matrix(s_x))
      data_nuisance$varrho_x   <- pmax(pmin(predict(nuisance$varrho_x,     newdata = as.matrix(x)),   prop_ub), prop_lb)
      data_nuisance$varrho_s_x <- pmax(pmin(predict(nuisance$varrho_s_x,   newdata = as.matrix(s_x)), prop_ub), prop_lb)
      data_nuisance$gamma_x    <- pmax(pmin(predict(nuisance$gamma_x,      newdata = as.matrix(x)),   prop_ub), prop_lb)
      data_nuisance$gamma_s_x  <- pmax(pmin(predict(nuisance$gamma_s_x,    newdata = as.matrix(s_x)), prop_ub), prop_lb)
    }
  }
  data_nuisance$pi    <- mean(data_train$observe)
  return(data_nuisance)
}

train_nuisance <- function(data_train, S_vars, X_vars, Y_var, type, obs,
                           nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                           xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads){
  #### Train nuisance parameters
  ### Train Long-term outcome means
  if(obs == T){
    ## \mu(s,x)
    mu_s_x_1_out <- mu_s_x(data_train, treatment_val = 1, S_vars, X_vars, Y_var, type,
                           nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                           xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
    mu_s_x_0_out <- mu_s_x(data_train, treatment_val = 0, S_vars, X_vars, Y_var, type,
                           nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                           xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
    ## \bar{\mu}(x)
    bar_mu_x_1_out <- bar_mu_x(data_train, treatment_val = 1, S_vars, X_vars, Y_var, mu_s_x_1_out, type,
                               nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                               xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
    bar_mu_x_0_out <- bar_mu_x(data_train, treatment_val = 0, S_vars, X_vars, Y_var, mu_s_x_0_out, type,
                               nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                               xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
  } else {
    ## \nu(s,x)
    nu_s_x_out <- nu_s_x(data_train, S_vars, X_vars, Y_var, type,
                         nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                         xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)

    ## \bar{\nu}(x)
    bar_nu_x_1_out <- bar_nu_x(data_train, treatment_val = 1, S_vars, X_vars, Y_var, nu_s_x_out, type,
                               nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                               xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
    bar_nu_x_0_out <- bar_nu_x(data_train, treatment_val = 0, S_vars, X_vars, Y_var, nu_s_x_out, type,
                               nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                               xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
  }

  ### Train Propensity scores
  if(obs == T){
    rho_1_out <- rho(data_train, treatment_val = 1, S_vars, X_vars, type,
                     nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                     xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
    rho_0_out <- rho(data_train, treatment_val = 0, S_vars, X_vars, type,
                     nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                     xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
    varrho_x_out <- varrho_x(data_train, X_vars, type,
                             nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                             xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
    gamma_x_out <- gamma_x(data_train, X_vars, type,
                           nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                           xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
  } else {
    varrho_x_out   <- varrho_x(data_train, X_vars, type,
                             nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                             xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
    varrho_s_x_out <- varrho_s_x(data_train, X_vars, S_vars, type,
                               nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                               xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
    gamma_x_out    <- gamma_x(data_train, X_vars, type,
                              nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                              xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
    gamma_s_x_out  <- gamma_s_x(data_train, X_vars, S_vars, type,
                                nuisance_cv_fold, grf_honesty, grf_tune_parameters, grf_num_threads,
                                xgb_cv_rounds, xgb_eta, xgb_max_depth, xgb_threads)
  }

  ### Package nuisance estimates
  if(obs == T){
    return(list(mu_s_x_1 = mu_s_x_1_out, mu_s_x_0 = mu_s_x_0_out,
                bar_mu_x_1 = bar_mu_x_1_out$bar_mu_x, bar_mu_x_1_cf = bar_mu_x_1_out$constant_flag,
                bar_mu_x_0 = bar_mu_x_0_out$bar_mu_x, bar_mu_x_0_cf = bar_mu_x_0_out$constant_flag,
                rho_1 = rho_1_out, rho_0 = rho_0_out,
                varrho_x = varrho_x_out, gamma_x = gamma_x_out))
  } else {
    return(list(nu_s_x = nu_s_x_out,
                bar_nu_x_1 = bar_nu_x_1_out$bar_nu_x, bar_nu_x_1_cf = bar_nu_x_1_out$constant_flag,
                bar_nu_x_0 = bar_nu_x_0_out$bar_nu_x, bar_nu_x_0_cf = bar_nu_x_0_out$constant_flag,
                varrho_x = varrho_x_out, varrho_s_x = varrho_s_x_out,
                gamma_x  = gamma_x_out,  gamma_s_x  = gamma_s_x_out))
  }
}

objective <- function(tau, data_nuisance, Y_var, obs, estimand){
  ### Compute the objective function
  if(obs == T){
    if(estimand == T){
      psi_1 <- compute_psi_1(tau, data_nuisance, Y_var)
      objective_value <- abs(mean(psi_1$psi_1))
    } else {
      psi_0 <- compute_psi_0(tau, data_nuisance, Y_var)
      objective_value <- abs(mean(psi_0$psi_0))
    }
  } else {
    if(estimand == T){
      xi_1  <- compute_xi_1(tau, data_nuisance, Y_var)
      objective_value <- abs(mean(xi_1$xi_1))
    } else {
      xi_0  <- compute_xi_0(tau, data_nuisance, Y_var)
      objective_value <- abs(mean(xi_0$xi_0))
    }
  }
  return(objective_value)
}

compute_psi_1 <- function(tau_1, data_nuisance, Y_var){
  ### Compute the vector of influence functions for a given value of tau_1, treatment observed
  psi_1 <- data_nuisance %>%
    mutate_(Y = parse_expr(Y_var)) %>%
    mutate(psi_1_a = (observe / pi)*(treatment*(Y - mu_s_x_1)/rho_1 - (1-treatment)*(Y - mu_s_x_0)/rho_0
                                     + (bar_mu_x_1 - bar_mu_x_0) - tau_1),
           psi_1_b = ((1-observe)/pi)*(gamma_x/(1-gamma_x)*(treatment*(mu_s_x_1  -bar_mu_x_1)/varrho_x)
                                       -(1-treatment)*(mu_s_x_0 - bar_mu_x_0)/(1 - varrho_x)),
           psi_1 = psi_1_a + psi_1_b) %>%
    select(psi_1) %>%
    return()
}

compute_psi_0 <- function(tau_0, data_nuisance, Y_var){
  ### Compute the vector of influence functions for a given value of tau_0, treatment observed
  psi_0 <- data_nuisance %>%
    mutate_(Y = parse_expr(Y_var)) %>%
    mutate(psi_0_a = (observe / (1-pi))*((1-gamma_x)/(gamma_x)*(treatment*(Y - mu_s_x_1)/rho_1 - (1-treatment)*(Y - mu_s_x_0)/rho_0)),
           psi_0_b = ((1-observe)/(1-pi))*((treatment*(mu_s_x_1  -bar_mu_x_1)/varrho_x)
                                       -(1-treatment)*(mu_s_x_0 - bar_mu_x_0)/(1 - varrho_x) + (bar_mu_x_1 - bar_mu_x_0) - tau_0),
           psi_0 = psi_0_a + psi_0_b) %>%
    select(psi_0) %>%
    return()
}

compute_xi_1 <- function(tau_1, data_nuisance, Y_var){
  ### Compute the vector of influence functions for a given value of tau_1, treatment not observed
  xi_1 <- data_nuisance %>%
    mutate_(Y = parse_expr(Y_var)) %>%
    mutate(xi_1_a = (observe / pi)*((gamma_x/gamma_s_x)*((1-gamma_s_x)/(1-gamma_x))*((varrho_s_x - varrho_x)*(Y - nu_s_x)/(varrho_x*(1-varrho_x)))
                                     + (bar_nu_x_1 - bar_nu_x_0) - tau_1),
           xi_1_b = ((1-observe)/pi)*((treatment*(nu_s_x - bar_nu_x_1)/varrho_x)
                                       -(1-treatment)*(nu_s_x - bar_nu_x_0)/(1 - varrho_x)),
           xi_1 = xi_1_a + xi_1_b) %>%
    select(xi_1) %>%
    return()
}

compute_xi_0 <- function(tau_0, data_nuisance, Y_var){
  ### Compute the vector of influence functions for a given value of tau_0, treatment not observed
  xi_0 <- data_nuisance %>%
    mutate_(Y = parse_expr(Y_var)) %>%
    mutate(xi_0_a = (observe / (1-pi))*(((1-gamma_s_x)/gamma_s_x)*((varrho_s_x - varrho_x)*(Y - nu_s_x)/(varrho_x*(1-varrho_x)))),
           xi_0_b = ((1-observe)/(1-pi))*((treatment*(nu_s_x - bar_nu_x_1)/varrho_x)
                                      -(1-treatment)*(nu_s_x - bar_nu_x_0)/(1 - varrho_x)
                                      + (bar_nu_x_1 - bar_nu_x_0) - tau_0),
           xi_0 = xi_0_a + xi_0_b) %>%
    select(xi_0) %>%
    return()
}

compute_ci <- function(hat_tau, data_nuisance, Y_var, obs, estimand, alpha){
  ### Compute 95% confidence intervals
  if(obs == T){
    if(estimand == T){
      influence <- compute_psi_1(hat_tau, data_nuisance, Y_var)
    } else {
      influence <- compute_psi_0(hat_tau, data_nuisance, Y_var)
    }
  } else {
    if(estimand == T){
      influence  <- compute_xi_1(hat_tau, data_nuisance, Y_var)
    } else {
      influence  <- compute_xi_0(hat_tau, data_nuisance, Y_var)
    }
  }
  hat_V <- (1/nrow(influence))*sum(influence^2)
  ci_l <- hat_tau - qnorm(1-alpha/2, 0, 1)*sqrt(hat_V/nrow(influence))
  ci_h <- hat_tau + qnorm(1-alpha/2, 0, 1)*sqrt(hat_V/nrow(influence))
  return(list(se = sqrt(hat_V/nrow(influence)), ci = c(ci_l, ci_h)))
}
