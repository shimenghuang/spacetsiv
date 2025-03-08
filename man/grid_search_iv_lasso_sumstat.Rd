% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iv_estimators.R
\name{grid_search_iv_lasso_sumstat}
\alias{grid_search_iv_lasso_sumstat}
\title{Grid search for IV-Lasso using summary statistics}
\usage{
grid_search_iv_lasso_sumstat(
  beta_init,
  gamma_res,
  sig_res,
  gamma_exp,
  sig_exp,
  cov_zinx = NULL,
  cov_ziny = NULL,
  pen_grid = pracma::logspace(-3, 1, 30),
  stop_if_accept = TRUE,
  stop_rule = c("max_pen", "min_pen", "max_pval"),
  alpha = 0.05,
  chisq_df = NULL,
  max_iter = 500,
  abs_tol = 1e-05,
  optim_method = "L-BFGS-B",
  control = list(maxit = 500, factr = 1e+06),
  refit_set = TRUE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{beta_init}{A vector of length K, the initial values for beta_hat for \code{optim}.}

\item{cov_zinx}{A matrix of dimension J x J, the covariance matrix of Z in
the sample containing Z and X. If NULL, the identity matrix will be used.}

\item{cov_ziny}{A matrix of dimension J x J, the covariance matrix of Z in
the sample containing Z and Y. If NULL, the identity matrix will be used.}

\item{pen_grid}{A vector of penalty parameter to go through.}

\item{stop_if_accept}{A boolean indicating whether to stop once a beta_hat
passes the heterogeneity test (\code{q_test_gen}) or continue going through
the penalty grid. If \code{stop_rule} is "max_pval" then it will still go
through the entire penalty grid.}

\item{stop_rule}{One of "min_pen", "max_pen", and "max_pval", for using
the minimal penalty that passes the heterogeneity test, the maximal penalty
that passes the heterogeneity test, or the penalty that gives the maximum
p-value respectively.}

\item{alpha}{A scalar between 0 and 1. The significance level to be used by
the Q test.}

\item{chisq_df}{An integer to overwrite the degrees of freedom in the Q test.
if NULL, then J - number of non-zero entry of beta_hat will be used. Note:
since this algorithm select the penalty parameter based on the p-value of
Q test, \code{chisq_df} can affect the final estimate.}

\item{max_iter}{An integer indicating the maximum number of iterations in
coordinate descent.}

\item{abs_tol}{A positive scalar for the maximum absolute difference between
two iterations of beta_hat, below which coordinate descent terminates.}

\item{refit_set}{A boolean indicating whether to refit the subset of non-zero
entries of beta_hat by optimizing the Q statistic.}

\item{verbose}{A boolean indicating whether to print out messages.}
}
\value{
A list containing final estimate beta_hat, the penalty value used,
and the p-value from the Q test based on the final beta_hat.
}
\description{
Grid search for IV-Lasso using summary statistics
}
