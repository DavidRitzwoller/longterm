% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/longterm.R
\name{longterm}
\alias{longterm}
\title{Semiparametric Estimation of Long-Term Treatment Effects}
\usage{
longterm(
  data,
  S_vars,
  X_vars,
  Y_var,
  obs,
  estimand,
  type,
  prop_lb = 0.01,
  prop_ub = 0.99,
  alpha = 0.05,
  te_lb = -500,
  te_ub = 500,
  cross_fit_fold = 5,
  nuisance_cv_fold = 5,
  grf_honesty = FALSE,
  grf_tune_parameters = "all",
  grf_num_threads = 1,
  xgb_cv_rounds = 100,
  xgb_eta = 0.1,
  xgb_max_depth = 2,
  xgb_threads = 1
)
}
\arguments{
\item{data}{A data frame containing the indicators \code{"treatment,"} denoting whether treatment was assigned, and \code{"observe,"} denoting whether the observation was from the observational or experimental sample, in addition to pretreatment covariates, short-term outcomes , and long-term outcomes. If an observation is missing values of a particular variable by construction, e.g., long-term outcomes are not observed in the experimental sample, then these values should be coded with any nonmissing value and will not contribute to estimation.}

\item{S_vars}{A list containing strings of the names of the short-term outcome variables.}

\item{X_vars}{A list containing strings of the names of the pre-treatment covariates.}

\item{Y_var}{A string giving the name of the long-term outcome of interest.}

\item{obs}{A boolean variable specifying whether treatment is observed in the long-term sample.}

\item{estimand}{A boolean variable specifying whether the estimand of interest is the long-term average treatment effect in the observational population, as opposed to the experimental population.}

\item{type}{A string specifying how nuisance functions should be estimated. \code{"glmnet"} specifies cross-validated lasso. \code{"grf"} specifies generalized random forests. \code{"xgboost"} specifies XGBoost.}

\item{prop_lb}{A float specifying the lower threshold for propensity score estimates.}

\item{prop_ub}{A float specifying the upper threshold for propensity score estimates.}

\item{alpha}{One minus the nominal coverage probability of the confidence intervals.}

\item{te_lb}{A float giving a lower bound for the long-term average treatment effect.}

\item{te_ub}{A float giving an upper bound for the long-term average treatment effect.}

\item{cross_fit_fold}{An integer giving the number of folds for the cross-fit estimation of nuisance parameters.}

\item{nuisance_cv_fold}{An interger giving the number of folds for cross-validating nuisance parameter estimates.}

\item{grf_honesty}{A boolean variable setting the \code{honesty} parameter for \code{grf} estimation.}

\item{grf_tune_parameters}{A string variable setting the \code{tune.parameters} parameter for \code{grf} estimation.}

\item{grf_num_threads}{An integer variable setting the \code{num.threads} parameter for \code{grf} estimation.}

\item{xgb_cv_rounds}{An integer variable setting maximum number of rounds for \code{xgboost} estimation.}

\item{xgb_eta}{A float variable setting \code{eta} parameter for \code{xgboost} estimation.}

\item{xgb_max_depth}{A float variable setting \code{max_depth} parameter for \code{xgboost} estimation.}

\item{xgb_threads}{An integer variable setting the \code{nthread} parameter for \code{xgboost} estimation.}
}
\value{
Returns a list with three components: \describe{

\item{\code{hat_tau}}{Estimate of the long-term average treatment effect.}

\item{\code{se}}{Estimate of the standard error of the estimator.}

\item{\code{ci}}{A vector giving the lower and upper points of the confidence interval.}

}
}
\description{
\code{longterm} estimates the long-term average treatment effect of a binary treatment on a scalar long-term outcome
using the semiparametric estimator developed in Chen and Ritzwoller (2021).
}
\references{
{
\cite{Chen, J., & Ritzwoller, D. M. (2021).
Semiparametric Estimation of Long-Term Treatment Effects.
arXiv preprint arXiv:2107.14405.}
}
}
