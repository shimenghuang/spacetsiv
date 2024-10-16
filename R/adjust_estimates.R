#' Adjust the estimates to account for correlated coefficients in two-sample MVMR.
#' 
#' @param coef_x A matrix of dimension J x K, the simple OLS coefficients of 
#'   each X on each Z.
#' @param cov_coef_x A matrix of dimension J x J, the covariance matrix of `coef_x`.
#' @param coef_y A vector of length J, the simple OLS coefficients of Y on each Z.
#' @param cov_coef_y A vector of length J, the variances of `coef_y`.
#' @param cor_zinx A matrix of dimension J x J, the correlation matrix of the instruments in sample with X (LD matrix).
#' @param cor_ziny A matrix of dimension J x J, the correlation matrix of the instruments in sample with Y (LD matrix).
#' @param cor_x A matrix of dimension K x K, the correlation matrix of X.
#' @param nobs_x The number of observations used to estimate `coef_x`.
#' @param nobs_y The number of observations used to estimated `coef_y`.
#' 
#' @references Wang and Kang 2021 for the univariate X case and Patel et al 2023 
#'   for multivariate X case.
#'   
#' @return A list containing 
#' 
#' @export
adjust_estimates <- function(coef_x, se_coef_x, coef_y, se_coef_y, 
                             cor_zinx, cor_ziny, cor_x, nobs_x, nobs_y,
                             ridge_pen_zinx = 0, ridge_pen_ziny = 0) {
  
  # instrument-exposure
  res_zx <- adjust_estimates_zx(coef_x, se_coef_x, 
                                cor_zinx + diag(rep(ridge_pen_zinx, nrow(cor_zinx))), 
                                cor_x, nobs_x)
  
  # instrument-outcome
  res_zy <- adjust_estimates_zy(coef_y, se_coef_y, 
                                cor_ziny + diag(rep(ridge_pen_zinx, nrow(cor_ziny))), 
                                nobs_y)
  
  return(list(coef_x = res_zx$gamma_X,
              cov_coef_x = res_zx$Sig_X,
              coef_y = res_zy$gamma_Y,
              cov_coef_y = res_zy$Sig_Y))
}

#' Adjust sufficient statistics to account for correlations
#' 
adjust_estimates_zx <- function(coef_x, se_coef_x, cor_zinx, cor_x, nobs_x) {
  K <- ncol(coef_x)
  ax <- 1/((nobs_x * se_coef_x^2) + coef_x^2) # J x K matrix
  Ax <- lapply(1:K, \(kk) {
    (sqrt(ax[,kk]) %*% t(sqrt(ax[,kk]))) * cor_zinx
  })
  bx <- ax * coef_x
  
  Sig_X_block <- function(ii, jj) {
    term1 <- solve(sqrt_mat(Ax[[ii]]) %*% t(sqrt_mat(Ax[[jj]])))
    term1 * as.numeric(cor_x[ii, jj] - t(bx[,ii]) %*% term1 %*% bx[,jj])
  }
  Sig_X <- lapply(1:K, \(ii) {
    blocks <- lapply(1:K, \(jj) {
      Sig_X_block(ii, jj)
    })
    do.call(cbind, blocks)
  })
  Sig_X <- do.call(rbind, Sig_X)
  Sig_X <- Sig_X/(nobs_x-ncol(cor_zinx)+1)
  gamma_X <- sapply(1:K, \(kk) {
    solve(Ax[[kk]], bx[,kk])
  })
  return(list(gamma_X = gamma_X,
              Sig_X = Sig_X))
}

#' Adjust sufficient statistics to account for correlations
#' 
adjust_estimates_zy <- function(coef_y, se_coef_y, cor_ziny, nobs_y) {
  ay <- 1/(nobs_y * se_coef_y^2 + coef_y^2) # vector of length J
  Ay <- (sqrt(ay) %*% t(sqrt(ay))) * cor_ziny # element-wise multiplication of two J x J matrix
  Ay_inv <- solve(Ay)
  by <- ay * coef_y
  Sig_Y <- Ay_inv * as.numeric((1 - t(by) %*% Ay_inv %*% by))
  Sig_Y <- Sig_Y/(nobs_y-ncol(cor_ziny)+1)
  gamma_Y <- as.vector(Ay_inv %*% by)
  return(list(gamma_Y = gamma_Y,
              Sig_Y = Sig_Y))
}
