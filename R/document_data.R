#' Outcomes from Graduation Experiment in Banerejee et. al (2015)
#'
#' This dataset contains a subset of the publicly available data from Banerjee et al (2015). It
#' contains the outcomes from a randomized evaluation of the long-term effects of a poverty alleviation program implemented in the Singh region of Pakistan.
#' @format A data frame with 854 rows corresponding to households
#'     67 columns corresponding to the variables:
#'
#' \describe{
#'
#' \item{id_hh}{Household ID}
#'
#' \item{treatment}{Assignment to Treatment}
#'
#' \item{*_bsl}{Pre-treatment covariateds}
#'
#' \item{*_end}{Two-year post-treatment outcomes}
#'
#' \item{*_fup}{Three-year post-treatment outcomes}
#'
#' }
#' @source \url{https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/NHIXNT}

#' @references{
#'
#' \cite{A. Banerjee et al., Science 348, 1260799 (2015). \doi{10.1126/science.1260799}}
#'
#' }
"graduation"
