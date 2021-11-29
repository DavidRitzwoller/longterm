## ---- include = FALSE---------------------------------------------------------
library("knitr")
knitr::opts_knit$set(self.contained = FALSE)
knitr::opts_chunk$set(tidy = TRUE, collapse=TRUE, comment = "#>",
                      tidy.opts=list(blank=FALSE, width.cutoff=55))
options(warn=-1)

## -----------------------------------------------------------------------------
library(longterm)

## Assign each observation to `experimental` and `observational` groups randomly
graduation$observe <- as.numeric(runif(nrow(graduation)) <= 0.5)

## Choose Pre-treatment covariates, Short-term Outcomes, and Long-term Outcome
X_vars <- c("ctotal_pcmonth_bsl", "cnonfood_pcmonth_bsl", 
            "cfood_pcmonth_bsl",  "cdurable_pcmonth_bsl",
            "asset_index_bsl",    "asset_prod_index_bsl")
S_vars <- c("ctotal_pcmonth_end", "cnonfood_pcmonth_end", 
            "cfood_pcmonth_end",  "cdurable_pcmonth_end",
            "asset_index_end",   "asset_prod_index_end")
Y_var  <- "ctotal_pcmonth_fup"

## Assume treatment is observed in the observational sample
# Average effect for the experimental population
est_1_1 <- longterm(graduation, S_vars = S_vars, X_vars = X_vars, Y_var = Y_var, 
                    obs = TRUE, estimand = TRUE, type = "glmnet")
# Average effect for the observational population
est_1_0 <- longterm(graduation, S_vars = S_vars, X_vars = X_vars, Y_var = Y_var, 
                    obs = TRUE, estimand = FALSE, type = "glmnet")

## Assume treatment is not observed in the observational sample
# Average effect for the experimental population
est_0_1 <- longterm(graduation, S_vars = S_vars, X_vars = X_vars, Y_var = Y_var, 
                    obs = FALSE, estimand = TRUE, type = "glmnet")
# Average effect for the observational population
est_0_0 <- longterm(graduation, S_vars = S_vars, X_vars = X_vars, Y_var = Y_var, 
                    obs = FALSE, estimand = FALSE, type = "glmnet")

## Display results
kable(data.frame(treatment_observed = c("Yes", "Yes", "No", "No"),
                 subpopulation = c("Observational", "Experimental", 
                                   "Observational", "Experimental"),
                 estimate = c(est_1_1$hat_tau, est_1_0$hat_tau, 
                              est_0_1$hat_tau, est_0_0$hat_tau),
                 standard_error = c(est_1_1$se, est_1_0$se, 
                                    est_0_1$se, est_0_0$se),
                 lower_ci = c(est_1_1$ci[1], est_1_0$ci[1], 
                              est_0_1$ci[1], est_0_0$ci[1]),
                 upper_ci = c(est_1_1$ci[2], est_1_0$ci[2], 
                              est_0_1$ci[2], est_0_0$ci[2])))


