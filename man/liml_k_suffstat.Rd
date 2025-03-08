% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/k-class_estimator.R
\name{liml_k_suffstat}
\alias{liml_k_suffstat}
\title{Calculate the k-value for LIML using sufficient statistics}
\usage{
liml_k_suffstat(XX, YY, ZZ, XY, ZX, ZY)
}
\arguments{
\item{XX}{A matrix of dimension R x R, the covariance matrix of the exposures.}

\item{YY}{A scalar, the variance of the response.}

\item{ZZ}{A matrix of dimension J x J, the covariance matrix of the instruments.}

\item{XY}{A vector of length R, the covariances between the exposures and the response.}

\item{ZX}{A matrix of dimension J x R, the covariances between the instruments and the exposures.}

\item{ZY}{A vector of length J, the covariances between the instruments and the response.}
}
\description{
Calculate the k-value for LIML using sufficient statistics
}
