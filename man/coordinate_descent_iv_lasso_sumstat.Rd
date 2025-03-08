% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iv_estimators.R
\name{coordinate_descent_iv_lasso_sumstat}
\alias{coordinate_descent_iv_lasso_sumstat}
\title{IV-Lasso using summary statistics}
\usage{
coordinate_descent_iv_lasso_sumstat(
  beta_init,
  gamma_res,
  sig_res,
  gamma_exp,
  sig_exp,
  cov_zinx = NULL,
  cov_ziny = NULL,
  pen = 0.01,
  max_iter = 500,
  abs_tol = 1e-05,
  verbose = FALSE
)
}
\description{
IV-Lasso using summary statistics
}
\details{
If the covariance of the instrument Z can be assumed the same
in the two samples, they can be dropped (pass in NULL). Otherwise,
this function requires knowledge of the two covariance matrices.
Following wang-kang2022_weak, we can recover cov_z if know cov_y
by using H (See Appendix 2.3).
}
\examples{
dat <- dgp_ex1(1e3, 1) # just-identified case
dat_sumstat <- comp_sumstat(dat$Y, dat$Z, dat$X, dat$Z, type = "depen_mv")
coordinate_descent_iv_lasso_sumstat(c(0.5, 1, 1.5), 
                                    dat_sumstat$coef_y, dat_sumstat$cov_coef_y, 
                                    dat_sumstat$coef_x, dat_sumstat$cov_coef_x, 
                                    cov(dat$Z), cov(dat$Z), 
                                    pen = 1.2, verbose = TRUE)
}
