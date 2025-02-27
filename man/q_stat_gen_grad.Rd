% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypothesis_tests.R
\name{q_stat_gen_grad}
\alias{q_stat_gen_grad}
\title{Gradient of the general Q-statistic}
\usage{
q_stat_gen_grad(coef_y, cov_coef_y, coef_x, cov_coef_x, beta)
}
\arguments{
\item{coef_y}{A vector of length J, the estimated instrument-outcome coefficients.}

\item{cov_coef_y}{A matrix of dimension J x J, the variance-covariance matrix of the estimated instrument-outcome coefficients.}

\item{coef_x}{A matrix of dimension J x K, the estimated instrument-exposure coefficients.}

\item{cov_coef_x}{A matrix of dimension (K x J) x (K x J) in a block form containing the variances and covariances of the estimated instrument-exposure coefficients. See patel-et-al2023 supplementary materials.}

\item{beta}{A vector of length K, the exposure-outcome coefficients.}
}
\description{
Gradient of the general Q-statistic
}
