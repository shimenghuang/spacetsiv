# ---- AR-subset (spaceIV) with one-sample individual level data ----

#' Inner SpaceIV estimator for a specific size of subsets
#' 
#' @template xyz
#' @param use_liml A boolean indicating whether to use LIML estimator.
#' 
#' @references Pfister and Peters 2022
#' 
spaceIV_estimator_inner <- function(X, Y, Z, size, use_liml = TRUE) {
  
  n <- nrow(X)
  d <- ncol(X)
  m <- ncol(Z)
  
  # normalize
  covZX <- cov(Z, X)
  covZY <- cov(Z, Y)
  
  # size cannot be larger than min(ncol(X), ncol(Z))
  size <- min(c(size, d, ncol(Z)))
  
  # all subsets
  subsets <- combn(1:d, size, simplify=FALSE)
  loss <- rep(NA, length(subsets))
  
  # pre-compute P_Z
  if(use_liml) {
    # P_Z <- Z %*% my_solve(t(Z) %*% Z) %*% t(Z)
    P_Z <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
  }
  
  for(i in 1:length(subsets)) {
    S <- subsets[[i]]
    beta_hat <- matrix(0, nrow=d, ncol=1)
    if(use_liml) {
      beta_hat[S] <- liml(X = X[,S,drop=FALSE], Y = Y, P_Z = P_Z, nobs = n)
      R  <- matrix(Y - X[,S,drop=FALSE] %*% beta_hat[S], ncol=1)
      loss[i] <- (t(R) %*% P_Z %*% R)/(t(R) %*% (diag(n) - P_Z) %*% R)*(n-m)/m
    }
    else{
      beta_hat[S] <- coefficients(lm.fit(covZX[,S,drop=FALSE], covZY))
      loss[i] <- sum((covZX %*% beta_hat - covZY)^2)
    }
  }
  
  # recompute beta for minimal loss
  S <- subsets[[which.min(loss)]]
  beta_hat <- matrix(0, nrow=d, ncol=1)  
  if (use_liml) {
    beta_hat[S] <- liml(X = X[,S,drop=FALSE], Y = Y, P_Z = P_Z, nobs = n)
  }
  else{
    beta_hat[S] <- coefficients(lm.fit(covZX[,S,drop=FALSE], covZY))
  }
  
  return(beta_hat)
}


#' SpaceIV estimator (subset selection with Anderson-Rubin test)
#' 
#' @template xyz
#' @param max_size The maximum size of subset to go through. 
#' @param alpha A numeric scalar in (0,1), the significance level.
#' @param use_liml A boolean indicating whether to use LIML estimator to 
#'   estimate the coefficients for each subset.
#' @param stop_if_accept A boolean indicating whether to stop going through 
#'   larger subsets once AR test is not rejected.
#' 
#' @export
spaceIV_estimator <- function(X, Y, Z, 
                              max_size = NULL, 
                              alpha = 0.05, 
                              chisq_df = NULL, 
                              use_liml = TRUE, 
                              stop_if_accept = TRUE,
                              stop_rule = c("min_set", "max_set", "max_pval")) {
  
  stop_rule <- match.arg(stop_rule)
  if (is.null(max_size)) {
    max_size <- ncol(X)
  } else {
    max_size <- max(max_size, ncol(X))
  }
  
  n <- nrow(X)
  m <- ncol(Z)
  d <- ncol(X)
  
  # size <- 0
  # accepted <- FALSE
  # while(size < max_size & !accepted){
  #   size <- size + 1
  #   beta_hat <- spaceIV_estimator_inner(X, Y, Z, size=size, use_liml=use_liml)
  #   # perform test
  #   pval <- anderson_rubin_test(X, Z, Y, beta_hat)$pval
  #   accepted <- pval >= alpha
  # }
  # return(list(beta_hat=beta_hat,
  #             pval=pval))
  
  beta_hat_all <- matrix(NA, nrow = d, ncol = max_size)
  pval_all <- rep(NA, max_size)
  sizes <- 1:max_size
  if (stop_if_accept & (stop_rule == "max_set")) {
    sizes <- max_size:1
  }
  for (size in sizes) {
    beta_hat <- spaceIV_estimator_inner(X = X, Y = Y, Z = Z, 
                                        size = size, use_liml = use_liml)
    pval <- anderson_rubin_test(X = X, Y = Y, Z = Z, beta_hat, 
                                chisq_df = chisq_df)$pval
    # accepted <- pval >= alpha
    accepted <- FALSE
    if (is.numeric(pval)) {
      accepted <- pval >= alpha
    }
    if (stop_if_accept & accepted) break 
    beta_hat_all[,size] <- beta_hat
    pval_all[size] <- pval
  }
  
  if (stop_if_accept) {
    return(list(beta_hat = beta_hat,
                pval = pval))
  } else {
    return(list(beta_hat_all = beta_hat_all,
                pval_all = pval_all))
  }
}

# ---- AR-subset (spaceIV) with (two-sample) sufficient statistics ----

spaceIV_estimator_suffstat_inner <- function(XX, YY, ZZ, XY, ZX, ZY, 
                                             nobs, size, 
                                             use_liml = FALSE,
                                             chisq_df = NULL) {
  m <- ncol(ZZ)
  d <- ncol(XX)
  
  # size cannot be larger than min(ncol(X), ncol(Z))
  size <- min(c(size, d, m))
  
  # all subsets
  subsets <- combn(1:d, size, simplify=FALSE)
  loss <- rep(NA, length(subsets))
  
  for(i in 1:length(subsets)) {
    S <- subsets[[i]]
    beta_hat <- matrix(0, nrow = d, ncol=1)
    if (use_liml){
      beta_hat[S] <- liml_suffstat(XX = XX[S,S,drop=FALSE], 
                                   YY = YY, 
                                   ZZ = ZZ, 
                                   XY = XY[S,,drop=FALSE],
                                   ZX = ZX[,S,drop=FALSE], 
                                   ZY = ZY)
      loss[i] <- anderson_rubin_test_suffstat(XX = XX, YY = YY, ZZ = ZZ,
                                              XY = XY, ZX = ZX, ZY = ZY, 
                                              nobs = nobs,
                                              beta = beta_hat, 
                                              chisq_df = chisq_df)$stat
    }
    else{
      beta_hat[S] <- coefficients(lm.fit(ZX[,S,drop=F], ZY))
      loss[i] <- sum((ZX %*% beta_hat - ZY)^2)
    }
  }
  
  # recompute beta for minimal loss
  S <- subsets[[which.min(loss)]]
  beta_hat <- matrix(0, nrow=d, ncol=1)  
  if(use_liml){
    beta_hat[S] <- liml_suffstat(XX = XX[S,S,drop=FALSE], 
                                 YY = YY, 
                                 ZZ = ZZ, 
                                 XY = XY[S,,drop=FALSE],
                                 ZX = ZX[,S,drop=FALSE], 
                                 ZY = ZY)
  }
  else{
    beta_hat[S] <- coefficients(lm.fit(ZX[,S,drop=FALSE], ZY))
  }
  
  return(beta_hat)
}

#' SpaceIV estimator using sufficient statistics
#' 
#' @export
spaceIV_estimator_suffstat <- function(XX, YY, ZZ, XY, ZX, ZY, nobs, 
                                       max_size = NULL, alpha = 0.05, 
                                       chisq_df = NULL,
                                       use_liml = FALSE) {
  
  if (is.null(max_size)) max_size <- ncol(XX)
  # run cross-validation fits
  # pval_mat <- matrix(NA, B, max_size)
  size <- 0
  accepted <- FALSE
  while(size < max_size & !accepted){
    size <- size + 1
    beta_hat <- spaceIV_estimator_suffstat_inner(XX = XX, YY = YY, ZZ = ZZ,
                                                 XY = XY, ZX = ZX, ZY = ZY, 
                                                 nobs = nobs,
                                                 size = size, 
                                                 use_liml = use_liml)
    # perform test
    pval <- anderson_rubin_test_suffstat(XX = XX, YY = YY, ZZ = ZZ,
                                         XY = XY, ZX = ZX, ZY = ZY, 
                                         nobs = nobs,
                                         beta = beta_hat, 
                                         chisq_df = chisq_df)$pval
    accepted <- pval >= alpha
  }
  
  return(list(beta_hat = beta_hat,
              pval = pval))
}

# ---- one-sample individual-level data IV-Lasso: minimizing 1/2*||Cov(Z, Y-Xb)||^2 + ||b||_1 ----

# adapted from Niklas python code (l1_penalized_estimator.py)

soft_threshold <- function(rho, zz, pen) {
  if (rho < -pen) return((rho + pen)/zz)
  else if (rho > pen) return((rho - pen)/zz)
  else return(0)
}

#' Coordinate descent with soft-thresholding for IV-Lasso
#' 
#' @export
coordinate_descent_iv_lasso <- function(beta_init, 
                                        X, Y, Z, 
                                        pen = 0.01, 
                                        max_iter = 500,
                                        abs_tol = 1e-5,
                                        verbose = FALSE) {
  n <- nrow(X)
  Y <- matrix(Y, ncol = 1)
  beta_hat <- matrix(beta_init, ncol = 1)
  ZX <- t(Z) %*% X / n # J x R
  ZY <- t(Z) %*% Y / n # J x 1
  for (ii in seq.int(max_iter)) {
    beta_hat_old <- beta_hat
    for (kk in 1:ncol(X)) {
      ZX_k <- ZX[,kk,drop=FALSE]
      Y_pred <- ZX %*% beta_hat
      rho <- mean(ZX_k * (ZY - Y_pred + beta_hat[kk] * ZX_k))
      zz <- mean(ZX_k^2)
      beta_hat[kk] <- soft_threshold(rho, zz, pen)
    }
    if (verbose) {
      message("iter = ", ii)
      message("beta = ", paste0(beta_hat, collapse = " "))
    }
    if (max(abs(beta_hat - beta_hat_old)) < abs_tol) break
  }
  return(beta_hat)
}

#' Grid search for IV-Lasso
#' 
#' @param beta_init A vector of initial value of beta. 
#' @template xyz
#' @param pen_grid A vector of penalty parameter (lambda) values to go through.
#' @param alpha A scalar in (0, 1) for the significance level for the AR test.
#' @param max_iter An integer indicating the maximum number of iterations in 
#'   coordinate descent.
#' @param refit_set A boolean indicating whether to refit the model 
#'   after selecting the variables with Lasso.
#' @param verbose A boolean passed to `coordinate_descent_iv_lasso` 
#'   for whether to print messages. 
#' 
#' @export
grid_search_iv_lasso <- function(beta_init, 
                                 X, Y, Z, 
                                 pen_grid = pracma::logspace(-3, 1, 30), 
                                 alpha = 0.05,
                                 chisq_df = NULL, 
                                 max_iter = 500,
                                 refit_set = FALSE,
                                 verbose = FALSE) {
  pvals <- rep(NA, length(pen_grid))
  for (ii in seq_along(pen_grid)) {
    pen <- pen_grid[ii]
    beta_hat <- coordinate_descent_iv_lasso(beta_init, 
                                            X = X, Y = Y, Z = Z,
                                            pen = pen, 
                                            max_iter = max_iter,
                                            verbose = max(verbose - 1, 0))
    res <- anderson_rubin_test(X = X, Y = Y, Z = Z, beta = beta_hat, 
                               chisq_df = chisq_df)
    pvals[ii] <- res$pval
    if (verbose) {
      message("pen = ", pen)
      message("AR stat = ", res$stat)
    }
    if (all(beta_hat == 0)) break
  }
  idx_pass <- which(pvals > alpha)
  if (length(idx_pass) == 0) {
    # message("No pen gives a beta_hat that passes the AR test, using the smallest given pen.")
    # pen_opt <- pen_grid[1]
    message("No pen gives a beta_hat that passes the AR test, using pen that gives the largest p-value.")
    pen_opt <- pen_grid[which.max(pvals)]
  } else {
    pen_opt <- pen_grid[rev(idx_pass)[1]]
  }
  
  beta_hat <- coordinate_descent_iv_lasso(beta_init, 
                                          X = X, Y = Y, Z = Z,
                                          pen = pen_opt, 
                                          max_iter = max_iter,
                                          verbose = max(verbose-1, 0))
  
  if (refit_set) {
    S <- which(abs(beta_hat) > pen_opt)
    if (length(S) > 0) {
      P_Z <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
      # P_Z <- Z %*% my_solve(t(Z) %*% Z) %*% t(Z)
      beta_hat <- rep(0, length(beta_init))
      beta_hat[S] <- liml(X = X[,S,drop=F], Y = Y, P_Z = P_Z, nobs = n)
    }
  }
  
  pval <- anderson_rubin_test(X = X, Y = Y, Z = Z, beta_hat, 
                              chisq_df = chisq_df)$pval
  
  return(list(beta_hat = beta_hat, 
              pen = pen_opt,
              pval = pval))
}

# ---- IV-Lasso with (two-sample) sufficient statistics ----

#' 
#' @export
coordinate_descent_iv_lasso_suffstat <- function(beta_init, 
                                                 ZX, ZY, 
                                                 pen = 0.01, 
                                                 max_iter = 500, 
                                                 abs_tol = 1e-5,
                                                 verbose = FALSE) {
  beta_hat <- matrix(beta_init, ncol = 1)
  for (ii in seq.int(max_iter)) {
    beta_hat_old <- beta_hat
    for (kk in 1:ncol(ZX)) {
      ZX_k <- ZX[,kk]
      Y_pred <- ZX %*% beta_hat
      rho <- mean(ZX_k * (ZY - Y_pred + beta_hat[kk] * ZX_k))
      zz <- mean(ZX_k^2)
      beta_hat[kk] <- soft_threshold(rho, zz, pen)
    }
    if (verbose) {
      message("iter = ", ii)
      message("beta = ", paste0(beta_hat, collapse = " "))
    }
    if (max(abs(beta_hat - beta_hat_old)) < abs_tol) break
  }
  return(beta_hat)
}

#' Grid search for IV-Lasso using sufficient statistics
#' 
#' @export
grid_search_iv_lasso_suffstat <- function(beta_init, 
                                          XX, YY, ZZ, XY, ZX, ZY, 
                                          nobs, 
                                          pen_grid = pracma::logspace(-3, 1, 30), 
                                          alpha = 0.05,
                                          chisq_df = NULL, 
                                          max_iter = 500, 
                                          abs_tol = 1e-5,
                                          refit_set = FALSE,
                                          verbose = FALSE) {
  pvals <- rep(NA, length(pen_grid))
  for (ii in seq_along(pen_grid)) {
    pen <- pen_grid[ii]
    beta_hat <- coordinate_descent_iv_lasso_suffstat(beta_init, ZX, ZY, 
                                                     pen = pen, 
                                                     max_iter = max_iter, 
                                                     abs_tol = abs_tol,
                                                     verbose = max(verbose-1, 0))
    res <- anderson_rubin_test_suffstat(XX = XX, YY = YY, ZZ = ZZ, XY = XY, ZX = ZX, ZY = ZY, 
                                        nobs = nobs,  beta = beta_hat, 
                                        chisq_df = chisq_df)
    pvals[ii] <- res$pval
    if (verbose) {
      message("pen = ", pen)
      message("AR stat = ", res$stat)
    }
    if (all(beta_hat == 0)) break
  }
  idx_pass <- which(pvals > alpha)
  if (length(idx_pass) == 0) {
    # message("No pen gives a beta_hat that passes the AR test, using the smallest given pen.")
    # pen_opt <- pen_grid[1]
    message("No pen gives a beta_hat that passes the AR test, using pen that gives the largest p-value.")
    pen_opt <- pen_grid[which.max(pvals)]
  } else {
    pen_opt <- pen_grid[rev(idx_pass)[1]]
  }
  
  beta_hat <- coordinate_descent_iv_lasso_suffstat(beta_init, 
                                                   ZX = ZX, 
                                                   ZY = ZY, 
                                                   pen = pen_opt, 
                                                   max_iter = max_iter, 
                                                   abs_tol = abs_tol,
                                                   verbose = max(verbose-1, 0))
  
  if (refit_set) {
    S <- which(abs(beta_hat) > pen_opt)
    if (length(S) > 0) {
      beta_hat <- rep(0, length(beta_init))
      beta_hat[S] <- liml_suffstat(XX = XX[S,S,drop=FALSE], 
                                   YY = YY, 
                                   ZZ = ZZ,
                                   XY = XY[S,,drop=FALSE],
                                   ZX = ZX[,S,drop=F], 
                                   ZY = ZY)
    }
  }
  
  pval <- anderson_rubin_test_suffstat(XX = XX, YY = YY, ZZ = ZZ, 
                                       XY = XY, ZX = ZX, ZY = ZY, 
                                       nobs = nobs,  
                                       beta = beta_hat, 
                                       chisq_df = chisq_df)$pval
  return(list(beta_hat = beta_hat, 
              pen = pen_opt,
              pval = pval))
}

# ---- IV-Lasso with two-sample summary statistics ----

#' IV-Lasso using summary statistics
#' 
#' @details If the covariance of the instrument Z can be assumed the same 
#'   in the two samples, they can be dropped (pass in NULL). Otherwise, 
#'   this function requires knowledge of the two covariance matrices. 
#'   Following wang-kang2022_weak, we can recover cov_z if know cov_y 
#'   by using H (See Appendix 2.3).
#' 
#' @examples 
#' dat <- dgp_ex1(1e3, 1) # just-identified case
#' dat_sumstat <- comp_sumstat(dat$Y, dat$Z, dat$X, dat$Z, type = "depen_mv")
#' coordinate_descent_iv_lasso_sumstat(c(0.5, 1, 1.5), 
#'                                     dat_sumstat$coef_y, dat_sumstat$cov_coef_y, 
#'                                     dat_sumstat$coef_x, dat_sumstat$cov_coef_x, 
#'                                     cov(dat$Z), cov(dat$Z), 
#'                                     pen = 1.2, verbose = TRUE)
#' @export 
coordinate_descent_iv_lasso_sumstat <- function(beta_init, 
                                                gamma_res, sig_res,
                                                gamma_exp, sig_exp,
                                                cov_zinx = NULL, 
                                                cov_ziny = NULL,
                                                pen = 0.01, 
                                                max_iter = 500, 
                                                abs_tol = 1e-5,
                                                verbose = FALSE) {
  
  # if not given, assume to be identity (same as dropping them)
  if (is.null(cov_zinx)) cov_zinx <- diag(length(gamma_res))
  if (is.null(cov_ziny)) cov_ziny <- diag(length(gamma_res))
  ge_covzx <- cov_zinx %*% gamma_exp # J x K
  if (is.null(beta_init)) {
    beta_init <- quantreg::rq(gamma_res ~ gamma_exp - 1, 
                              tau = 0.5, 
                              weights = ivw_weights(gamma_res, diag(sig_res)))$coefficients
  }
  beta_hat <- matrix(beta_init, ncol = 1)
  for (ii in seq.int(max_iter)) {
    beta_hat_old <- beta_hat
    for (kk in 1:ncol(gamma_exp)) {
      ge_covzx_k <- ge_covzx[,kk,drop=FALSE]
      Y_pred <- ge_covzx %*% beta_hat
      rho <- mean(ge_covzx_k * (cov_ziny %*% gamma_res - Y_pred + beta_hat[kk] * ge_covzx_k))
      zz <- mean(ge_covzx_k^2)
      beta_hat[kk] <- soft_threshold(rho, zz, pen)
    }
    if (verbose) {
      message("iter = ", ii)
      message("beta = ", paste0(beta_hat, collapse = " "))
    }
    if (max(abs(beta_hat - beta_hat_old)) < abs_tol) break
  }
  return(c(beta_hat))
}

#' Grid search for IV-Lasso using summary statistics
#' 
#' @param beta_init A vector of length K, the initial values for beta_hat for `optim`.
#' @param cov_zinx A matrix of dimension J x J, the covariance matrix of Z in 
#'   the sample containing Z and X. If NULL, the identity matrix will be used. 
#' @param cov_ziny A matrix of dimension J x J, the covariance matrix of Z in 
#'   the sample containing Z and Y. If NULL, the identity matrix will be used. 
#' @param pen_grid A vector of penalty parameter to go through.
#' @param stop_if_accept A boolean indicating whether to stop once a beta_hat 
#'   passes the heterogeneity test (`q_test_gen`) or continue going through 
#'   the penalty grid. If `stop_rule` is "max_pval" then it will still go 
#'   through the entire penalty grid.
#' @param stop_rule One of "min_pen", "max_pen", and "max_pval", for using  
#'   the minimal penalty that passes the heterogeneity test, the maximal penalty 
#'   that passes the heterogeneity test, or the penalty that gives the maximum 
#'   p-value respectively. 
#' @param alpha A scalar between 0 and 1. The significance level to be used by 
#'   the Q test.
#' @param chisq_df An integer to overwrite the degrees of freedom in the Q test.
#'   if NULL, then J - number of non-zero entry of beta_hat will be used. Note: 
#'   since this algorithm select the penalty parameter based on the p-value of 
#'   Q test, `chisq_df` can affect the final estimate. 
#' @param max_iter An integer indicating the maximum number of iterations in 
#'   coordinate descent. 
#' @param abs_tol A positive scalar for the maximum absolute difference between 
#'   two iterations of beta_hat, below which coordinate descent terminates. 
#' @param refit_set A boolean indicating whether to refit the subset of non-zero 
#'   entries of beta_hat by optimizing the Q statistic.
#' @param verbose A boolean indicating whether to print out messages.
#' 
#' @return A list containing final estimate beta_hat, the penalty value used, 
#'   and the p-value from the Q test based on the final beta_hat.
#' 
#' @export
grid_search_iv_lasso_sumstat <- function(beta_init, 
                                         gamma_res, sig_res,
                                         gamma_exp, sig_exp,
                                         cov_zinx = NULL, 
                                         cov_ziny = NULL, 
                                         pen_grid = pracma::logspace(-3, 1, 30), 
                                         stop_if_accept = TRUE, 
                                         stop_rule = c("max_pen", "min_pen", "max_pval"), 
                                         alpha = 0.05,
                                         chisq_df = NULL, 
                                         max_iter = 500,
                                         abs_tol = 1e-5,
                                         optim_method = "L-BFGS-B",
                                         control = list(maxit = 500, 
                                                        factr = 1e6),
                                         refit_set = TRUE,
                                         verbose = FALSE,
                                         ...) {
  # set stopping rule
  stop_rule <- match.arg(stop_rule)
  
  pvals <- rep(NA, length(pen_grid))
  pens <- sort(pen_grid)
  if (stop_if_accept & (stop_rule == "max_pen")) {
    pens <- rev(pens)
  }
  for (ii in seq_along(pens)) {
    pen <- pens[ii]
    beta_hat <- coordinate_descent_iv_lasso_sumstat(beta_init,
                                                    gamma_res, sig_res,
                                                    gamma_exp, sig_exp,
                                                    cov_zinx, cov_ziny,
                                                    pen = pen, 
                                                    max_iter = max_iter,
                                                    abs_tol = abs_tol,
                                                    verbose = max(verbose-1, 0))
    
    if (refit_set) {
      S <- which(abs(beta_hat) > abs_tol)
      if (length(S) > 0) {
        if (is.null(beta_init)) {
          beta_init_ss <- quantreg::rq(gamma_res ~ gamma_exp[,S,drop=FALSE] - 1, 
                                       tau = 0.5, 
                                       weights = ivw_weights(gamma_res, 
                                                             sig_res))$coefficients
        } else {
          beta_init_ss <- beta_init[S]
        }
        beta_hat <- rep(0, length(beta_init))
        res <- optim(beta_init_ss, fn = \(bb) {
          q_test_gen(coef_y = gamma_res, 
                     cov_coef_y = sig_res, 
                     coef_x = gamma_exp[,S,drop=FALSE], 
                     cov_coef_x = get_submat(S, sig_exp, length(gamma_res)), 
                     beta = bb,
                     chisq_df = chisq_df,
                     tol = abs_tol)$stat
        }, 
        method = optim_method, 
        hessian = FALSE, 
        # control = list(maxit = max_iter, abstol = abs_tol),
        control = control, ...)
        beta_hat[S] <- res$par
      }
    }
  
    res <- q_test_gen(coef_y = gamma_res,
                      cov_coef_y = sig_res, 
                      coef_x = gamma_exp,
                      cov_coef_x = sig_exp, 
                      beta = beta_hat, 
                      chisq_df = chisq_df, 
                      tol = abs_tol)
    if (is.na(res$pval)) {
      message("non-positive chisq_df, pval set to 0")
      pvals[ii] <- 0
    } else {
      pvals[ii] <- res$pval
    }
    if (verbose) {
      message("pen = ", pen)
      message("Q stat = ", res$stat)
    }
    if (all(beta_hat == 0) & (pvals[ii] > alpha)) break
    if (stop_if_accept & (pvals[ii] > alpha) & (stop_rule != "max_pval")) break
    # if ((stop_rule == "min_pen") & (pvals[ii] > alpha)) break
  }
  idx_pass <- which(pvals > alpha)
  if (length(idx_pass) == 0) {
    # message("No pen gives a beta_hat that passes the Q test, using the smallest given pen.")
    # pen_opt <- pen_grid[1]
    message("No pen gives a beta_hat that passes the Q test, using pen that gives the largest p-value.")
    # message("Sample size: ", nobs)
    pen_opt <- pens[which.max(pvals)]
  } else {
    if ((stop_rule == "min_pen") | (stop_rule == "max_pen")) {
      pen_opt <- pens[idx_pass[1]]
    } else if (stop_rule == "max_pval") {
      pen_opt <- pens[which.max(pvals)]
    } else {
      stop('stop_rule needs to be one of "min_pen", "max_pen", or "max_pval"')
    }
  }
  if (verbose) message("pen_opt = ", pen_opt)
  
  if (stop_if_accept) {
    beta_hat <- coordinate_descent_iv_lasso_sumstat(beta_init,
                                                    gamma_res, sig_res,
                                                    gamma_exp, sig_exp,
                                                    cov_zinx, cov_ziny,
                                                    pen = pen_opt, 
                                                    max_iter = max_iter, 
                                                    abs_tol = abs_tol,
                                                    verbose = max(verbose-1, 0))
    if (refit_set) {
      # S <- which(abs(beta_hat) > pen_opt)
      S <- which(abs(beta_hat) > abs_tol)
      if (length(S) > 0) {
        if (is.null(beta_init)) {
          beta_init_ss <- quantreg::rq(gamma_res ~ gamma_exp[,S,drop=FALSE] - 1, 
                                       tau = 0.5, 
                                       weights = ivw_weights(gamma_res, 
                                                             sig_res))$coefficients
        } else {
          beta_init_ss <- beta_init[S]
        }
        beta_hat <- rep(0, length(beta_init))
        res <- optim(beta_init_ss, fn = \(bb) {
          q_test_gen(coef_y = gamma_res, 
                     cov_coef_y = sig_res, 
                     coef_x = gamma_exp[,S,drop=FALSE], 
                     cov_coef_x = get_submat(S, sig_exp, length(gamma_res)), 
                     beta = bb,
                     chisq_df = chisq_df,
                     tol = abs_tol)$stat
        }, 
        method = optim_method, 
        hessian = FALSE, 
        # control = list(maxit = max_iter, abstol = abs_tol),
        control = control, ...)
        beta_hat[S] <- res$par
      }
    }
    
    pval <- q_test_gen(coef_y = gamma_res,
                       cov_coef_y = sig_res, 
                       coef_x = gamma_exp,
                       cov_coef_x = sig_exp, 
                       beta = beta_hat, 
                       chisq_df = chisq_df)$pval
    
    return(list(beta_hat = beta_hat, 
                pen = pen_opt,
                pval = pval))
  } else {
    if (length(idx_pass) == 0) {
      idx_pass <- which.max(pvals)
    }
    beta_hat_all <- lapply(idx_pass, \(idx) {
      pen <- pens[idx]
      beta_hat_orig <- coordinate_descent_iv_lasso_sumstat(beta_init,
                                                           gamma_res, sig_res,
                                                           gamma_exp, sig_exp,
                                                           cov_zinx, cov_ziny,
                                                           pen = pen, 
                                                           max_iter = max_iter, 
                                                           abs_tol = abs_tol,
                                                           verbose = max(verbose-1, 0))
      if (refit_set) {
        S <- which(abs(beta_hat_orig) > abs_tol)
        beta_hat_refit <- rep(0, ncol(gamma_exp))
        if (length(S) > 0) {
          if (is.null(beta_init)) {
            beta_init_ss <- quantreg::rq(gamma_res ~ gamma_exp[,S,drop=FALSE] - 1, 
                                         tau = 0.5, 
                                         weights = ivw_weights(gamma_res, 
                                                               sig_res))$coefficients
          }
          else {
            beta_init_ss <- beta_init[S]
          }
          res <- optim(beta_init_ss, fn = \(bb) {
            q_test_gen(coef_y = gamma_res, 
                       cov_coef_y = sig_res, 
                       coef_x = gamma_exp[,S,drop=FALSE], 
                       cov_coef_x = get_submat(S, sig_exp, length(gamma_res)), 
                       beta = bb,
                       chisq_df = chisq_df,
                       tol = abs_tol)$stat
          }, 
          method = optim_method, 
          hessian = FALSE, 
          # control = list(maxit = max_iter, abstol = abs_tol),
          control = control, ...)
          beta_hat_refit[S] <- res$par
        }
      } else {
        beta_hat_refit <- rep(NA, ncol(gamma_exp))
      }
      list(beta_hat_orig = beta_hat_orig,
           beta_hat_refit = beta_hat_refit)
    })
    beta_hat_all_orig <- sapply(beta_hat_all, \(ll) {
      ll$beta_hat_orig
    }, simplify = TRUE)
    beta_hat_all_refit <- sapply(beta_hat_all, \(ll) {
      ll$beta_hat_refit
    }, simplify = TRUE)
    
    if (refit_set) {
      pval_all_orig <- sapply(seq_along(idx_pass), \(ii) {
        q_test_gen(coef_y = gamma_res,
                   cov_coef_y = sig_res, 
                   coef_x = gamma_exp,
                   cov_coef_x = sig_exp, 
                   beta = beta_hat_all_orig[,ii], 
                   chisq_df = chisq_df)$pval
      })
      pval_all_refit <- sapply(seq_along(idx_pass), \(ii) {
        q_test_gen(coef_y = gamma_res,
                   cov_coef_y = sig_res, 
                   coef_x = gamma_exp,
                   cov_coef_x = sig_exp, 
                   beta = beta_hat_all_refit[,ii], 
                   chisq_df = chisq_df)$pval
      })
    } else {
      pval_all_orig <- sapply(seq_along(idx_pass), \(ii) {
        q_test_gen(coef_y = gamma_res,
                   cov_coef_y = sig_res, 
                   coef_x = gamma_exp,
                   cov_coef_x = sig_exp, 
                   beta = beta_hat_all_orig[,ii], 
                   chisq_df = chisq_df)$pval
      })
      pval_all_refit <- NULL
    }
    
    return(list(beta_hat_all_orig = beta_hat_all_orig,
                beta_hat_all_refit = beta_hat_all_refit,
                pval_all_orig = pval_all_orig,
                pval_all_refit = pval_all_refit))
  }
  
}

# ---- Q-subset for two-sample summary statistics ----

#' Calculate weights for the IVW estimator. 
#' 
#' @details TODO: does this apply in the case of multiple exposures?
#' 
#' @param gamma_res A vector of length J, the estimated instrument-outcome effects. 
#' @param sig_res A vector of length J, the variance of `gamma_res`. 
#' 
#' @export
ivw_weights <- function(gamma_res, sig_res) {
  ws <- gamma_res^2 / sig_res
  ws <- ws / sum(ws)
  return(ws)
}

#' Fixed-size subsets optimization for the sparseMVMR estimator.
#' 
#' @param gamma_res A vector of length J.
#' @param gamma_exp A matrix of dimension J x K.
#' @param sig_res If `indep_inst = TRUE`, a vector of length J, the standard 
#'   error of the estimated instrument-outcome coefficients; if `indep_inst = FALSE`, 
#'   a matrix of dimension J x J, the variance-covariance matrix of the 
#'   estimated instrument-outcome coefficients. 
#' @param sig_exp A list of variance-covariance matrices each of dimension J x K.
#' @param size A scalar indicating the size of the subsets of the exposures.
#' @param optim_method A string indicating the method used in `optim`.
#' @param control A list of contrl variables passing to `optim`. For example, 
#'   `list(reltol = sqrt(.Machine$double.eps))` if using "BFGS" method. 
#' @param ... Other parameters to pass on to `optim`. For example `lower` and 
#'   `upper` when using "L-BFGS-B" method. 
#' 
sparseMVMR_estimator_inner <- function(gamma_res, sig_res,
                                       gamma_exp, sig_exp,
                                       size, 
                                       beta_init = NULL,
                                       indep_inst = TRUE, 
                                       chisq_df = NULL, 
                                       abs_tol = 1e-5, 
                                       optim_method = "L-BFGS-B",
                                       control = list(factr = 1e6),
                                       ...) {
  
  # size cannot be larger than the number of exposure
  K <- ncol(gamma_exp)
  size <- min(c(size, K))
  
  # Note: this sets beta_init once for all exposures
  # if (is.null(beta_init)) {
  #   beta_init <- rep(0.5, ncol(gamma_exp))
  # }
  
  # all subsets
  subsets <- combn(1:K, size, simplify=FALSE)
  subset_names <- sapply(subsets, \(ss) paste0(ss, collapse = "_"))
  
  # save loss and beta_hat
  loss_all <- rep(NA, length(subsets))
  names(loss_all) <- subset_names
  beta_hat_all <- matrix(0, nrow = length(subsets), ncol = K)

  # local function
  optim_subset <- function(S) {
    if (indep_inst) {
      # TODO: use something like IVW here!
      if (is.null(beta_init)) {
        # beta_init_ss <- rep(0.5, length(S))
        # beta_init_ss <- L1pack::lad(gamma_res ~ gamma_exp[,S,drop=FALSE], 
        #                             method = "EM")
        beta_init_ss <- quantreg::rq(gamma_res ~ gamma_exp[,S,drop=FALSE] - 1, 
                                     tau = 0.5, 
                                     weights = ivw_weights(gamma_res, 
                                                           sig_res^2))$coefficients
      } else {
        beta_init_ss <- beta_init[S]
      }
      res <- optim(beta_init_ss, fn = \(bb) {
        q_test_mv(coef_y = gamma_res, 
                  se_coef_y = sig_res, 
                  coef_x = gamma_exp[,S,drop=FALSE],
                  cov_coef_x = lapply(sig_exp, \(ss) {
                    ss[S, S, drop=FALSE]
                  }),
                  beta = bb, 
                  chisq_df = chisq_df)$stat
      }, method = optim_method, control = control, ...)
    } else {
      if (is.null(beta_init)) {
        # TODO: allow passing in cov_zinx and cov_ziny!
        beta_init_ss <- iv2_estimator_sumstat(gamma_res = gamma_res, 
                                              sig_res = sig_res, 
                                              gamma_exp = gamma_exp[,S,drop=FALSE], 
                                              sig_exp = get_submat(S, sig_exp, length(gamma_res)), 
                                              cov_zinx = NULL, 
                                              cov_ziny = NULL)
      } else {
        beta_init_ss <- beta_init[S]
      }
      res <- optim(beta_init_ss, fn = \(bb) {
        q_test_gen(coef_y = gamma_res, 
                   cov_coef_y = sig_res, 
                   coef_x = gamma_exp[,S,drop=FALSE], 
                   cov_coef_x = get_submat(S, sig_exp, length(gamma_res)), 
                   beta = bb,
                   chisq_df = chisq_df,
                   tol = abs_tol)$stat
      }, method = optim_method, control = control, ...)
    }
    res
  }
  
  for (ii in 1:length(subsets)) {
    S <- subsets[[ii]]
    res <- optim_subset(S)
    beta_hat_all[ii, S] <- res$par
    loss_all[ii] <- res$value
  }
  rownames(beta_hat_all) <- subset_names
  
  # recompute beta for minimal loss
  S <- subsets[[which.min(loss_all)]]
  res <- optim_subset(S)
  beta_hat <- matrix(rep(0, K))
  beta_hat[S] <- res$par
  
  # calculate p-value of min(Qs)
  if (is.null(chisq_df)) {
    # chisq_df <- length(gamma_res) - size
    chisq_df <- length(gamma_res)
  }
  # pval <- pchisq(min(loss_all), 
  #                df = chisq_df, 
  #                lower.tail = FALSE)
  pval_all <- pchisq(loss_all, 
                     df = chisq_df, 
                     lower.tail = FALSE)
  names(pval_all) <- subset_names
  pval <- pval_all[which.min(loss_all)]
  
  return(list(beta_hat = beta_hat,
              pval = pval,
              beta_hat_all = t(beta_hat_all), 
              stat_all = loss_all,
              pval_all = pval_all))
}

#' Subset selection with Q-test
#' 
#' @param gamma_res A vector of length J, where J is the number of instruments.
#' @param sig_res If `indep_inst = TRUE`, a vector of length J, the standard 
#'   error of the estimated instrument-outcome coefficients; if `indep_inst = FALSE`, 
#'   a matrix of dimension J x J, the variance-covariance matrix of the 
#'   estimated instrument-outcome coefficients. 
#' @param gamma_exp A matrix of dimension J x K, where K is the number of exposures.
#' @param sig_exp A list of length J of matrices each of dimension K x K.
#' @param indep_inst A boolean indicating whether to assume the instruments are 
#'   independent. This affects the dimension of `sig_res`.
#' @param alpha A scalar indicating the confidence level.
#' @param chisq_df An integer to overwrite the degrees of freedom in the Q test.
#'   if NULL, then J - number of non-zero entry of beta_hat will be used.
#' @param optim_method A string indicating the method used in `optim`.
#' @param control A list of contrl variables passing to `optim`. For example, 
#'   `list(reltol = sqrt(.Machine$double.eps))` if using "BFGS" method. 
#' @param ... Other parameters to pass on to `optim`. For example `lower` and 
#'   `upper` when using "L-BFGS-B" method. 
#' 
#' @export
sparseMVMR_estimator <- function(gamma_res, sig_res,
                                 gamma_exp, sig_exp,
                                 beta_init = NULL, 
                                 indep_inst = TRUE, 
                                 alpha = 0.05,
                                 chisq_df = NULL, 
                                 max_size = NULL,
                                 stop_if_accept = TRUE,
                                 stop_rule = c("min_set", "max_set", "max_pval"), 
                                 digits = 3,
                                 optim_method = "L-BFGS-B", 
                                 control = list(factr = 1e6),
                                 verbose = FALSE,
                                 ...) {
  
  stop_rule <- match.arg(stop_rule)
  
  if (is.null(max_size)) {
    max_size <- min(ncol(gamma_exp), nrow(gamma_exp))
  } else {
    max_size <- min(max_size, ncol(gamma_exp))
  }
  
  # best beta_hat and corresponding p-value for each subset size
  beta_hat_best <- matrix(NA, nrow = ncol(gamma_exp), ncol = max_size)
  pval_best <- rep(NA, max_size)
  beta_hat_all <- vector("list", max_size)
  stat_all <- vector("list", max_size) # Q-stat for all subsets
  pval_all <- vector("list", max_size)
  sizes <- 1:max_size
  if (stop_if_accept & (stop_rule == "max_set")) {
    sizes <- max_size:1
  }
  for (size in sizes) {
    if (verbose) message("current size: ", size)
    res <- sparseMVMR_estimator_inner(gamma_res, sig_res,
                                      gamma_exp, sig_exp,
                                      beta_init = beta_init, 
                                      size = size,
                                      chisq_df = chisq_df, 
                                      indep_inst = indep_inst, 
                                      optim_method = optim_method,
                                      control = control,
                                      ...)
    beta_hat <- res$beta_hat
    accepted <- FALSE
    if (is.numeric(res$pval)) {
      accepted <- res$pval >= alpha
    }
    if (stop_if_accept & accepted & (stop_rule != "max_pval")) break 
    beta_hat_best[,size] <- beta_hat
    pval_best[size] <- res$pval
    stat_all[[size]] <- res$stat_all
    pval_all[[size]] <- res$pval_all
    beta_hat_all[[size]] <- res$beta_hat_all
  }
  
  if (!stop_if_accept) {
    message("Note: when stop_if_accept is FALSE, stop_rule is ignored.")
    return(list(beta_hat_best = beta_hat_best,
                pval_best = pval_best,
                beta_hat_all = beta_hat_all,
                stat_all = stat_all,
                pval_all = pval_all,
                pvals = summarize_pvals(pval_all = pval_all,
                                        var_names = as.character(1:length(beta_hat)),
                                        digits = digits)))
  } else if (accepted & stop_rule == "max_pval") {
    return(list(beta_hat = beta_hat_best[,which.max(pval_best)],
                pval = max(pval_best)))
  } else if (accepted) {
    return(list(beta_hat = c(beta_hat),
                pval = res$pval))
  } else {
    idx <- which.max(pval_best)
    if (!accepted) {
      message("All subsets rejected at level alpha. Returning beta_hat with the largest best p-value.")
    }
    return(list(beta_hat = c(beta_hat_best[,idx]),
                pval = pval_best[idx]))
  } 
}

#' Summarize p-values for subset selection
#' 
#' @param pval_all A list of lists each containing the p-values of sets 
#'   of one size. 
#'   
#' @export
summarize_pvals <- function(pval_all, var_names, digits = NULL) {
  gather_pval <- function(vv) {
    unlist(
      sapply(pval_all, \(ll) {
        subset_names_split <- stringr::str_split(names(ll), "_")
        flags <- sapply(subset_names_split, \(ss) {
          !(vv %in% ss)
        }, simplify = TRUE)
        ll[flags]
      }, simplify = TRUE)
    )
  }
  pval_max_other <- rep(NA, length(var_names))
  for (ii in seq_along(var_names)) {
    pvals <- gather_pval(var_names[ii])
    pval_max_other[ii] <- max(pvals)
  }
  if (!is.null(digits)) {
    pval_max_other <- round(pval_max_other, digits)
  }
  pval_df <- data.frame(x = var_names,
                        pval = pval_max_other)
  return(pval_df)
}

# ---- Q-Lasso for two-sample summary statistics ----

#' Extract the ij-th block of a squared block structured matrix
#' 
#' @param dim_block A scalar for the dimension of each block.
#' 
get_block <- function(ii, jj, block_mat, dim_block) {
  row_idx <- seq(from = (ii-1) * dim_block + 1, length.out = dim_block)
  col_idx <- seq(from = (jj-1) * dim_block + 1, length.out = dim_block)
  return(block_mat[row_idx, col_idx,drop=FALSE])
}

#' Extract a submatrix of a squared block structured matrix
#' 
#' @details Extract the submatrix containing row and column blocks indexed 
#'   by elements in S. 
#' 
#' @param dim_block A scalar for the dimension of each block.
get_submat <- function(S, block_mat, dim_block) {
  idx_all <- c(sapply(S, \(ss) {
    seq(from = (ss-1) * dim_block + 1, length.out = dim_block)
  }))
  return(block_mat[idx_all, idx_all,drop=FALSE])
}


#' Subgradient of general Q-statistic with L1-penalty
#' 
q_stat_gen_l1 <- function(Y, sig_Y, X, sig_X, beta, lambda, tol = 1e-5) {
  
  # objective value
  Y <- c(Y)
  gg <- Y - X %*% beta
  phi <- kronecker(t(beta), diag(nrow(X)))
  Ome <- sig_Y + phi %*% sig_X %*% t(phi)
  # Ome <- (Ome + t(Ome)) / 2
  # Ome_inv <- my_solve(Ome)
  # qq <- c(t(gg) %*% Ome_inv %*% gg + lambda * sum(abs(beta)))
  qq <- c(t(gg) %*% solve(Ome, gg) + lambda * sum(abs(beta)))
  
    
  # gradient and subgradient
  Ome_inv <- solve(Ome)
  term1 <- c(-2 * t(X) %*% Ome_inv %*% gg)
  tmp <- lapply(1:ncol(X), \(kk) {
    Reduce(`+`, lapply(1:ncol(X), \(ll) {
      beta[ll] * get_block(kk, ll, sig_X, nrow(X))
    }))
  })
  term2 <- sapply(1:ncol(X), \(kk) {
    -2 * t(gg) %*% Ome_inv %*% tmp[[kk]] %*% Ome_inv %*% gg
  })
  term2 <- c(term2)
  fg <- term1 + term2
  # res <- fg + lambda * sign(beta)
  # res[abs(beta) < ] <- 0
  res <- rep(0, length(fg))
  res[abs(beta) >= tol] <- (fg + lambda * sign(beta))[abs(beta) >= tol]
  res[(abs(beta) < tol) & (fg < -lambda)] <- (fg + lambda)[(abs(beta) < tol) & (fg < -lambda)]
  res[(abs(beta) < tol) & (fg > lambda)] <- (fg - lambda)[(abs(beta) < tol) & (fg > lambda)]
  res[(abs(beta) < tol) & (fg >= -lambda) & (fg <= lambda)] <- 0
  return(list(fval = qq, gval = res))
}

#' Grid search for Q-Lasso.
#' 
#' @export
grid_search_q_stat_gen_l1 <- function(beta_init, 
                                      gamma_res, sig_res,
                                      gamma_exp, sig_exp,
                                      pen_grid = pracma::logspace(-3, 1, 30), 
                                      stop_if_accept = TRUE, 
                                      stop_rule = c("max_pen", "min_pen", "max_pval"), 
                                      alpha = 0.05, 
                                      chisq_df = NULL, 
                                      abs_tol = 1e-5, 
                                      optim_method = "L-BFGS-B",
                                      control = list(maxit = 500, 
                                                     factr = 1e6),
                                      # max_iter = 500, 
                                      # abs_tol = 1e-5,
                                      refit_set = FALSE,
                                      verbose = FALSE, 
                                      ...) {
  
  # set stopping rule
  stop_rule <- match.arg(stop_rule)
  
  pvals <- rep(NA, length(pen_grid))
  pens <- sort(pen_grid)
  if (stop_if_accept & (stop_rule == "max_pen")) {
    pens <- rev(pens)
  }
  for (ii in seq_along(pens)) {
    pen <- pens[ii]
    opt_res <- optim(beta_init, 
                     fn = \(bb) {
                       q_stat_gen_l1(gamma_res, sig_res,
                                     gamma_exp, sig_exp,
                                     bb,
                                     lambda = pen, tol = abs_tol)$fval / 2
                     },
                     gr = \(bb) {
                       q_stat_gen_l1(gamma_res, sig_res,
                                     gamma_exp, sig_exp,
                                     bb,
                                     lambda = pen, tol = abs_tol)$gval / 2
                     }, 
                     method = optim_method, 
                     hessian = FALSE, 
                     # control = list(maxit = max_iter, abstol = abs_tol),
                     control = control, ...)
    # check if opt_res$convergence is 0 (converged)
    if (opt_res$convergence == 0) {
      beta_hat <- opt_res$par
    }
    else {
      message("pen = ", pen, " not converged! error code: ", opt_res$convergence)
      beta_hat <- rep(NA, length(beta_init))
    }
    res <- q_test_gen(coef_y = gamma_res, cov_coef_y = sig_res, 
                      coef_x = gamma_exp, cov_coef_x = sig_exp, 
                      beta = beta_hat, 
                      chisq_df = chisq_df,
                      tol = abs_tol)
    if (is.na(res$pval)) {
      message("negative chisq_df, pval set to 0")
      pvals[ii] <- 0
    } else {
      pvals[ii] <- res$pval
    }
    if (verbose) {
      message("pen = ", pen)
      message("Q stat = ", res$stat)
    }
    if (stop_if_accept & (pvals[ii] > alpha)) break
    if (all(!is.na(beta_hat))) {
      if (all(beta_hat == 0) & (pvals[ii] > alpha)) break
    }
  }
  idx_pass <- which(pvals > alpha)
  # if (length(idx_pass) == 0) {
  #   # message("No pen gives a beta_hat that passes the Q test, using the smallest given pen.")
  #   # pen_opt <- pen_grid[1]
  #   message("No pen gives a beta_hat that passes the AR test, using pen that gives the largest p-value.")
  #   pen_opt <- pen_grid[which.max(pvals)]
  # } else {
  #   pen_opt <- pen_grid[rev(idx_pass)[1]]
  # }
  if (length(idx_pass) == 0) {
    # message("No pen gives a beta_hat that passes the Q test, using the smallest given pen.")
    # pen_opt <- pen_grid[1]
    message("No pen gives a beta_hat that passes the Q test, using pen that gives the largest p-value.")
    # message("Sample size: ", nobs)
    pen_opt <- pens[which.max(pvals)]
  } else {
    if ((stop_rule == "min_pen") | (stop_rule == "max_pen")) {
      pen_opt <- pens[idx_pass[1]]
    } else if (stop_rule == "max_pval") {
      pen_opt <- pens[which.max(pvals)]
    } else {
      stop('stop_rule needs to be one of "min_pen", "max_pen", or "max_pval"')
    }
  }
  if (verbose) message("pen_opt = ", pen_opt)
  
  opt_res <- optim(beta_init, 
                   fn = \(bb) {
                     q_stat_gen_l1(gamma_res, sig_res,
                                   gamma_exp, sig_exp,
                                   bb,
                                   lambda = pen_opt, 
                                   tol = abs_tol)$fval / 2
                   },
                   gr = \(bb) {
                     q_stat_gen_l1(gamma_res, sig_res,
                                   gamma_exp, sig_exp,
                                   bb,
                                   lambda = pen_opt, 
                                   tol = abs_tol)$gval / 2
                   }, 
                   method = optim_method, 
                   hessian = FALSE, 
                   # control = list(maxit = max_iter, abstol = abs_tol),
                   control = control, ...)
  if (opt_res$convergence == 0) {
    beta_hat <- opt_res$par
    beta_hat[which(abs(beta_hat) <= abs_tol)] <- 0
  }
  else {
    message("opt_pen = ", pen_opt, " not converged! error code: ", opt_res$convergence)
    beta_hat <- rep(NA, length(beta_init))
  }
  
  if (refit_set) {
    # S <- which(abs(beta_hat) > pen_opt)
    S <- which(abs(beta_hat) > abs_tol)
    if (length(S) > 0) {
      beta_hat <- rep(0, length(beta_init))
      opt_res <- optim(beta_hat[S], 
                       fn = \(bb) {
                         q_test_gen(gamma_res, sig_res,
                                    gamma_exp[,S,drop=FALSE], 
                                    sig_exp[rep(S, each=nrow(gamma_exp)), 
                                            rep(S, each=nrow(gamma_exp)), 
                                            drop=FALSE], 
                                    bb, 
                                    chisq_df = chisq_df, 
                                    tol = abs_tol)$stat / 2
                       }, 
                       gr = \(bb) {
                         q_stat_gen_grad(gamma_res, sig_res,
                                         gamma_exp[,S,drop=FALSE], 
                                         sig_exp[rep(S, each=nrow(gamma_exp)), 
                                                 rep(S, each=nrow(gamma_exp)), 
                                                 drop=FALSE], 
                                         bb) / 2
                       }, 
                       method = optim_method, 
                       hessian = FALSE, 
                       # control = list(maxit = max_iter, abstol = abs_tol),
                       control = control, ...)
      if (opt_res$convergence == 0) {
        beta_hat[S] <- opt_res$par
      }
      else {
        message("refit not converged! error code: ", opt_res$convergence)
        beta_hat[S] <- rep(NA, length(S))
      }
    }
  }
  pval <- q_test_gen(coef_y = gamma_res, cov_coef_y = sig_res, 
                     coef_x = gamma_exp, cov_coef_x = sig_exp, 
                     beta = beta_hat, 
                     chisq_df = chisq_df,
                     tol = abs_tol)$pval
  return(list(beta_hat = beta_hat, 
              pen = pen_opt,
              pval = pval))
}


# ---- non-sparse IV estimators using individual-level data ----

#' IV estimator under just-identified scenario
#' 
#' @details The simplest just-identified IV estimator.
#' 
#' @examples 
#' dat <- dgp_ex00(1e3, c(1, 2))
#' iv_estimator(dat$X, dat$Y, dat$Z)
#' 
#' @export
iv_estimator <- function(X, Y, Z) {
  return(c(solve(t(Z) %*% X, t(Z) %*% Y)))
}

#' Two-sample IV estimator under just-identified scenario
#' 
#' @details The simplest just-identified IV estimator.
#' 
#' @examples 
#' dat1 <- dgp_ex00(1e3, c(1, 2))
#' dat2 <- dgp_ex00(1e3, c(1, 2))
#' tsiv_estimator(dat1$X, dat1$Z, dat2$Y, dat2$Z)
#' 
#' @export
tsiv_estimator <- function(X, Z_X, Y, Z_Y) {
  return(c(solve(t(Z_X) %*% X, t(Z_Y) %*% Y)))
}

#' Two-stage least squares (TSLS) estimator
#' 
#' @details Cannot be used in under-identified case.
#'   In just-identified case, this gives the same result as 
#'   `iv_estimator`.
#' 
#' @examples
#' dat <- dgp_ex00(1e3, c(1, 2))
#' tsls_estimator(dat$X, dat$Y, dat$Z)
#' dat <- dgp_ex0(1e3, c(1, 2))
#' tsls_estimator(dat$X, dat$Y, dat$Z)
#' 
#' @export
tsls_estimator <- function(X, Y, Z) {
  Pz <- Z %*% solve(crossprod(Z), t(Z))
  # Pz <- Z %*% my_solve(crossprod(Z), t(Z))
  return(c(solve(t(X) %*% Pz %*% X, t(X) %*% Pz %*% Y)))
}

#' Two-sample two-stage least-squares estimator
#' 
#' @examples
#' dat1 <- dgp_ex00(1e3, c(1, 2))
#' dat2 <- dgp_ex00(1e3, c(1, 2))
#' tstsls_estimator(dat1$X, dat1$Z, dat2$Y, dat2$Z)
#' 
#' @export
tstsls_estimator <- function(X, Z_X, Y, Z_Y) {
  gamma_X <- solve(crossprod(Z_X), t(Z_X) %*% X)
  # gamma_X <- my_solve(crossprod(Z_X), t(Z_X) %*% X)
  X_hat <- Z_Y %*% gamma_X
  return(c(solve(crossprod(X_hat), t(X_hat) %*% Y)))
}

#' Minimizing ||Cov(Z, Y-Xb)||^2
#' 
#' @details IV solutions using valid instruments satisfy Cov(Z, Y-Xb) = 0. 
#'   The objective here is to minimize the squared LHS.
#'   Cannot be used in under-identified case.
#'   In just-identified case, this gives the same result as 
#'   `iv_estimator`.
#' 
#' @examples 
#' dat <- dgp_ex00(1e4, c(1, 2))
#' iv2_estimator(dat$X, dat$Y, dat$Z)
#' dat <- dgp_ex0(1e3, c(1, 2))
#' iv2_estimator(dat$X, dat$Y, dat$Z)
#' 
#' @export
iv2_estimator <- function(X, Y, Z) {
  if (ncol(Z) < ncol(X)) stop("Under-identified.")
  ZtX <- t(Z) %*% X
  ZtY <- t(Z) %*% Y
  return(c(solve(crossprod(ZtX), t(ZtX) %*% ZtY)))
}

#' Minimizing ||Cov(Z, Y-Xb)||^2 using coordinate descent
#' 
#' @details IV solutions using valid instruments satisfy Cov(Z, Y-Xb) = 0. 
#'   The objective here is to minimize the squared LHS.
#' 
#' @examples 
#' dat <- dgp_ex00(1e3, c(1, 2)) # just-identified case
#' coordinate_descent_iv2_estimator(c(1, 1), dat$X, dat$Y, dat$Z, verbose = FALSE)
#' dat <- dgp_ex0(1e3, c(1, 2)) # over-identified case
#' coordinate_descent_iv2_estimator(c(1, 1), dat$X, dat$Y, dat$Z, verbose = FALSE)
#' dat <- dgp_ex1(1e3, 2) # under-identified case
#' coordinate_descent_iv2_estimator(c(1, 1, 1), dat$X, dat$Y, dat$Z, verbose = FALSE)
#' coordinate_descent_iv2_estimator(c(0, 2, 0), dat$X, dat$Y, dat$Z, verbose = FALSE)
#' 
#' @export
coordinate_descent_iv2_estimator <- function(beta_init, 
                                             X, Y, Z, 
                                             max_iter = 500,
                                             abs_tol = 1e-5,
                                             verbose = FALSE) {
  n <- nrow(X)
  Y <- matrix(Y, ncol = 1)
  beta_hat <- matrix(beta_init, ncol = 1)
  ZX <- t(Z) %*% X / n # J x R
  ZY <- t(Z) %*% Y / n # J x 1
  for (ii in seq.int(max_iter)) {
    beta_hat_old <- beta_hat
    for (kk in 1:ncol(X)) {
      ZX_k <- ZX[,kk,drop=FALSE]
      Y_pred <- ZX %*% beta_hat
      # rho <- mean(ZX_k * (ZY - Y_pred + beta_hat[kk] * ZX_k))
      rho <- t(ZX_k) %*% (ZY - Y_pred + beta_hat[kk] * ZX_k)
      zz <- sum(ZX_k^2)
      beta_hat[kk] <- rho / zz
    }
    if (verbose) {
      message("iter = ", ii)
      message("beta = ", paste0(beta_hat, collapse = " "))
    }
    if (max(abs(beta_hat - beta_hat_old)) < abs_tol) break
  }
  return(c(beta_hat))
}

# ---- non-sparse IV estimators using summary statistics ----

#' Minimizing ||Cov(Z, Y-Xb)||^2 using summary statistics and variance of Z
#' 
#' @details IV solutions using valid instruments satisfy Cov(Z, Y-Xb) = 0. 
#'   The objective here is to minimize the squared LHS.
#'   Cannot be used in under-identified case.
#'   In just-identified case, this gives the same result as 
#'   `just_identified_iv_estimator`.
#'   TODO: currently not using the sigs...
#' 
#' @examples 
#' dat <- dgp_ex00(1e4, c(1, 2))
#' dat_sumstat <- comp_sumstat(dat$Y, dat$Z, dat$X, dat$Z, type = "depen_mv")
#' iv2_estimator_sumstat(dat_sumstat$coef_y, dat_sumstat$cov_coef_y, 
#'                       dat_sumstat$coef_x, dat_sumstat$cov_coef_x, 
#'                       cov(dat$Z), cov(dat$Z))
#' iv2_estimator(dat$X, dat$Y, dat$Z)
#' dat <- dgp_ex0(1e3, c(1, 2))
#' dat_sumstat <- comp_sumstat(dat$Y, dat$Z, dat$X, dat$Z, type = "depen_mv")
#' iv2_estimator_sumstat(dat_sumstat$coef_y, dat_sumstat$cov_coef_y, 
#'                       dat_sumstat$coef_x, dat_sumstat$cov_coef_x, cov(dat$Z))
#' iv2_estimator(dat$X, dat$Y, dat$Z)
#' 
#' @export
iv2_estimator_sumstat <- function(gamma_res, sig_res,
                                  gamma_exp, sig_exp,
                                  cov_zinx = NULL, 
                                  cov_ziny = NULL,
                                  use_ginv = FALSE) {
  
  if (is.null(cov_zinx)) cov_zinx <- diag(length(gamma_res))
  if (is.null(cov_ziny)) cov_ziny <- diag(length(gamma_res))
  cov_zinx2 <- cov_zinx %*% cov_zinx
  cov_zinxy <- cov_zinx %*% cov_ziny
  term1 <- t(gamma_exp) %*% cov_zinx2 %*% gamma_exp
  term2 <- t(gamma_exp) %*% cov_zinxy %*% gamma_res
  
  if (use_ginv) {
    return(c(MASS::ginv(term1) %*% term2))
  } else {
    return(c(solve(term1, term2)))
    # return(c(my_solve(term1, term2)))
  }
}

#' Inverse-variance weighted estimator
#' 
#' @details IVW assumes independent instruments. 
#' 
#' @param weight_order One of "first" or "second" indicating whether to use 
#'   first or second order weights (inverse variance, which is approximated 
#'   using delta method).
#' 
#' @references Bowden, J., Del Greco M, F., Minelli, C., Zhao, Q., 
#'   Lawlor, D.A., Sheehan, N.A., Thompson, J. and Davey Smith, G., 2019. 
#'   Improving the accuracy of two-sample summary-data Mendelian randomization: 
#'   moving beyond the NOME assumption. International journal of epidemiology, 
#'   48(3), pp.728-742.
#' @export
ivw_estimator <- function(gamma_res, sig_res, 
                          gamma_exp, sig_exp, 
                          weight_order = c("first", "second")) {
  
  # num <- sum(gamma_exp * gamma_res / sig_res)
  # den <- sum(gamma_exp^2 / sig_res)
  
  weight_order <- match.arg(weight_order)
  if (weight_order == "first") {
    weights <- gamma_exp^2 / sig_res 
  } else {
    weights <- 1/(sig_res / gamma_exp^2 + 
                    gamma_res^2 * sig_exp / gamma_exp^4)
  }
  beta_hats <- gamma_res / gamma_exp
  num <- sum(weights * beta_hats)
  den <- sum(weights)
  
  return(list(beta_hat = num / den,
              weights = weights))
}
