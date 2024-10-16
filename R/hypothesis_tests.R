#' Anderson-Rubin test
#' 
#' @template xyz
#' @template beta_kdim
#' 
#' @template test-return
#'   
#' @export
anderson_rubin_test <- function(X, Y, Z, beta, chisq_df = NULL) {
  
  n <- nrow(X)
  d <- ncol(X)
  m <- ncol(Z)
  if (is.null(chisq_df)) {
    # chisq_df <- max(m - d, 0)
    chisq_df <- max(m, 0)
  }
  P <- Z %*% solve(t(Z) %*% Z, t(Z))
  # P <- Z %*% my_solve(t(Z) %*% Z, t(Z))
  res <- matrix(Y - X %*% beta, ncol=1)
  # tstat <- (t(res) %*% P %*% res)/(t(res) %*% (diag(n) - P) %*% res) * (n-m)/m
  stat <- (t(res) %*% P %*% res)/(t(res) %*% (diag(n) - P) %*% res) * n
  if (chisq_df <= 0) {
    message("non-positive chisq_df")
    pval <- 0
  } else {
    pval <- pchisq(stat, df = chisq_df, lower.tail = FALSE)
  }
  
  return(list(stat = stat,
              # stat = tstat,
              # pval = pf(tstat, df1 = m, df2 = n-m, lower.tail = FALSE),
              pval = pval))
}

#' Anderson-Rubin test with summary statistics
#' 
#' @template xyz_suffstat
#' @template nobs
#' @template beta_kdim
#' 
#' @template test-return
#' 
#' @export
anderson_rubin_test_suffstat <- function(XX, YY, ZZ, XY, ZX, ZY, 
                                         nobs, beta, chisq_df = NULL) {
  
  d <- ncol(XX)
  m <- ncol(ZZ)
  if (is.null(chisq_df)) {
    # chisq_df <- max(m - d, 0)
    chisq_df <- max(m, 0)
  }
  beta_vec <- matrix(beta, ncol = 1)
  ZZ_inv <- solve(ZZ)
  # ZZ_inv <- my_solve(ZZ)
  num <- t(ZY) %*% ZZ_inv %*% ZY - 
    2 * t(beta_vec) %*% t(ZX) %*% ZZ_inv %*% ZY +
    t(beta_vec) %*% t(ZX) %*% ZZ_inv %*% ZX %*% beta_vec
  den <- YY - 2 * t(beta_vec) %*% XY + t(beta_vec) %*% XX %*% beta_vec - num
  # tstat <- num/den * (nobs-m)/m
  stat <- num/den * nobs
  if (chisq_df <= 0) {
    message("non-positive chisq_df")
    pval <- 0
  } else {
    pval <- pchisq(stat, df = chisq_df, lower.tail = FALSE)
  }
  
  return(list(stat = stat,
              # stat = tstat,
              # pval = pf(tstat, df1 = m, df2 = nobs-m, lower.tail = FALSE),
              pval = pval))
}

#' Single-variable Q-test given independent instruments.
#' 
#' @param coef_y A vector of length J, the estimated instrument-outcome coefficients.
#' @param se_coef_y A vector of length J, the standard error of the estimated instrument-outcome coefficients.
#' @param coef_x A vector of length J, the estimated instrument-exposure coefficients.
#' @param se_coef_x A vector of length J, the standard error of the estimated instrument-exposure coefficients.
#' @param beta A scalar, the exposure-outcome coefficient to be tested.
#' 
#' @template test-return
#' 
#' @export
q_test <- function(coef_y, se_coef_y, coef_x, se_coef_x, beta) {
  coef_y <- c(coef_y)
  se_coef_y <- c(se_coef_y)
  coef_x <- c(coef_x)
  se_coef_x <- c(se_coef_x)
  qq <- sapply(seq_along(coef_y), \(jj) {
    num <- coef_y[jj] - coef_x[jj] * beta
    den <- sqrt(se_coef_x[jj]^2 * beta^2 + se_coef_y[jj]^2)
    num/den
  })
  qq <- sum(qq^2)
  return(list(stat = qq,
              pval = pchisq(qq, df = length(coef_y), lower.tail = FALSE)))
}

#' Multivariable Q-test given independent instruments.
#' 
#' @param coef_y A vector of length J, the estimated instrument-outcome coefficients.
#' @param se_coef_y A vector of length J, the standard error of the estimated instrument-outcome coefficients.
#' @param coef_x A matrix of dimension J x K, the estimated instrument-exposure coefficients.
#' @param cov_coef_x A list of J matrices each of dimension K x K, the covariance matrix of the estimated instrument-exposure coefficients.
#' @param beta A vector of length K, the exposure-outcome coefficient to be tested.
#' 
#' @return p-value of the test.
#' 
#' @export
q_test_mv <- function(coef_y, se_coef_y, coef_x, cov_coef_x, beta, 
                      chisq_df = NULL, tol = 1e-3) {
  if (is.null(chisq_df)) {
    # chisq_df <- length(coef_y) - sum(abs(beta) > tol)
    chisq_df <- length(coef_y)
  }
  coef_y <- c(coef_y)
  se_coef_y <- c(se_coef_y)
  qq <- sapply(seq_along(coef_y), \(jj) {
    num <- coef_y[jj] - t(coef_x[jj,]) %*% beta
    den <- sqrt(t(beta) %*% cov_coef_x[[jj]] %*% beta + se_coef_y[jj]^2)
    num/den
  })
  qq2sum <- sum(qq^2)
  if (chisq_df <= 0) {
    message("non-positive chisq_df")
    pval <- 0
  } else {
    pval <- pchisq(qq2sum, 
                   df = chisq_df, 
                   lower.tail = FALSE)
  }
  return(list(stat = qq2sum,
              pval = pval,
              qq = qq^2) # individual contribution to the Q-stat
  )
}

#' General Q-test that allows correlated instruments.
#'   
#' @template xyz_sumstat_gen
#' @template beta_kdim
#' @param tol A scalar, tolerance for considering an entry of beta to be zero.
#' 
#' @template test-return
#' 
#' @export
q_test_gen <- function(coef_y, cov_coef_y, coef_x, cov_coef_x, beta, 
                       chisq_df = NULL, tol = 1e-5) {
  if (anyNA(beta)) {
    return(list(stat = Inf,
                pval = 0))
  }
  if (is.null(chisq_df)) {
    # chisq_df <- length(coef_y) - sum(abs(beta) > tol)
    chisq_df <- length(coef_y)
  }
  coef_y <- c(coef_y)
  gg <- coef_y - coef_x %*% beta
  phi <- kronecker(t(beta), diag(nrow(coef_x)))
  Ome <- cov_coef_y + phi %*% cov_coef_x %*% t(phi)
  qq <- t(gg) %*% solve(Ome, gg)
  # Ome <- (Ome + t(Ome)) / 2
  # qq <- t(gg) %*% my_solve(Ome, gg)
  stat <- as.numeric(qq)
  if (chisq_df <= 0) {
    message("non-positive chisq_df")
    pval <- 0
  } else if (sum(abs(beta) > tol) > length(coef_y)) {
    message("more non-zero causal effects than the number of instruments")
    pval <- 0
  } else {
    pval <- pchisq(stat, 
                   df = chisq_df,
                   # df = length(coef_y),
                   # df = length(coef_y) - sum(abs(beta) > tol), 
                   lower.tail = FALSE)
  }
  return(list(stat = stat,
              pval = pval))
}


#' Gradient of the general Q-statistic
#' 
#' @template xyz_sumstat_gen
#' @template beta_kdim
#' 
#' @export
q_stat_gen_grad <- function(coef_y, cov_coef_y, coef_x, cov_coef_x, beta) {
  
  # objective value
  coef_y <- c(coef_y)
  gg <- coef_y - coef_x %*% beta
  phi <- kronecker(t(beta), diag(nrow(coef_x)))
  Ome <- cov_coef_y + phi %*% cov_coef_x %*% t(phi)
  qq <- t(gg) %*% solve(Ome, gg)
  # Ome <- (Ome + t(Ome)) / 2
  # Ome_inv <- my_solve(Ome)
  # qq <- t(gg) %*% Ome_inv %*% gg
  
  # gradient 
  Ome_inv <- solve(Ome)
  term1 <- -2 * t(coef_x) %*% Ome_inv %*% gg
  term1 <- as.matrix(term1, ncol=1)
  tmp <- lapply(1:ncol(coef_x), \(kk) {
    Reduce(`+`, lapply(1:ncol(coef_x), \(ll) {
      beta[ll] * get_block(kk, ll, cov_coef_x, nrow(coef_x))
    }))
  })
  term2 <- sapply(1:ncol(coef_x), \(kk) {
    -2 * t(gg) %*% Ome_inv %*% tmp[[kk]] %*% Ome_inv %*% gg
  })
  term2 <- as.matrix(term2, ncol=1)
  return(term1 + term2)
}

#' Gradient of the general Q-statistic w.r.t beta_k
#' 
q_stat_gen_grad_k <- function(kk, coef_y, cov_coef_y, coef_x, cov_coef_x, beta) {
  
  # objective value
  coef_y <- c(coef_y)
  gg <- coef_y - coef_x %*% beta
  phi <- kronecker(t(beta), diag(nrow(coef_x)))
  Ome <- cov_coef_y + phi %*% cov_coef_x %*% t(phi)
  qq <- t(gg) %*% solve(Ome, gg)
  # Ome <- (Ome + t(Ome)) / 2
  # Ome_inv <- my_solve(Ome)
  # qq <- t(gg) %*% Ome_inv %*% gg
  
  # gradient w.r.t beta_j
  Ome_inv <- solve(Ome)
  term1 <- -2 * t(coef_x[,kk]) %*% Ome_inv %*% gg
  tmp <- Reduce(`+`, lapply1:ncol(coef_x), \(ll) {
    beta[ll] * get_block(kk, ll, cov_coef_x, nrow(coef_x))
  })
  term2 <- -2 * t(gg) %*% Ome_inv %*% tmp %*% Ome_inv %*% gg
  return(as.numeric(term1 + term2))
}
