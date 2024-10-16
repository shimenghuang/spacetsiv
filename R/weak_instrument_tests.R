# ---- First-stage F test ----

#' First-stage F test using individual data
#' 
#' @param X A matrix of dimension nobs x K, the individual observations of exposures.
#' @param Z A matrix of dimension nobs x J, the individual observations of instruments.
#' @param all_z A boolean indicating whether to condition on all Z together.
#' 
fsf_test_individual <- function(X, Z, all_z = FALSE) {
  if (all_z) {
    Pz <- Z %*% solve(crossprod(Z), t(Z))
    stats <- rep(NA, ncol(X))
    pvals <- rep(NA, ncol(X))
    for (kk in 1:ncol(X)) {    
      mod <- lm(X[,kk] ~ Z)
      sig2 <- mean(resid(mod)^2)
      stats[kk] <- t(X[,kk]) %*% Pz %*% X[,kk] / sig2
      names(stats) <- paste0("Z", 1:ncol(Z))
      pvals[kk] <- pchisq(stats[kk], df = ncol(Z), lower.tail = FALSE)
      names(pvals) <- paste0("X", 1:ncol(X))
    }
  } else {
    stats <- matrix(NA, nrow = ncol(Z), ncol = ncol(X))
    pvals <- matrix(NA, nrow = ncol(Z), ncol = ncol(X))
    for (kk in 1:ncol(X)) {
      for (jj in 1:ncol(Z)) {
        Pz <- Z[,jj] %*% solve(crossprod(Z[,jj]), t(Z[,jj]))
        mod <- lm(X[,kk] ~ Z[,jj])
        sig2 <- mean(resid(mod)^2)
        stats[jj, kk] <- t(X[,kk]) %*% Pz %*% X[,kk] / sig2
        rownames(stats) <- paste0("Z", 1:ncol(Z))
        colnames(stats) <- paste0("X", 1:ncol(X))
        pvals[jj,kk] <- pchisq(stats[jj, kk] , df = 1, lower.tail = FALSE)
        rownames(pvals) <- paste0("Z", 1:ncol(Z))
        colnames(pvals) <- paste0("X", 1:ncol(X))
      }
    }
  }
  return(list(stats = stats,
              pvals = pvals))
}

#' First-stage F test using summary statistics for one exposure case
#' 
#' @param coef_x A vector of length J, the effect estimate conditioning on all Z.
#' @param se_x A vector of length J, the standard error of `coef_x`.
#' 
fsf_test_sumstat_simple <- function(coef_x, se_x) {
  stat <- coef_x^2 / se_x^2 
  pval <- pchisq(stat, df = 1, lower.tail = FALSE)
  return(list(stat = stat,
              pval = pval))
}

# ---- Conditional F test ----

#' Conditional F test using individual level data
#'
#' @references Sanderson and Windmeijer (2016) Section 4.3 eq (15)
condf_test_kk <- function(kk, X, Z, cond_set = NULL, chisq_df = NULL) {
  
  if (is.null(cond_set)) cond_set <- setdiff(1:ncol(X), kk)
  if (kk %in% cond_set) stop("kk should not be in cond_set")
  if (is.null(chisq_df)) chisq_df <- ncol(Z) - length(cond_set) + 1 # ncol(Z) - ncol(X) + 1
  
  mod1 <- lm.fit(y = X[,cond_set,drop=FALSE], x = Z)
  X_hat <- mod1$fitted.values
  
  mod2 <- lm.fit(y = X[,kk], x = as.matrix(X_hat))
  delta_tilde <- coef(mod2)
  lhs <- c(X[,kk] - X[,cond_set,drop=FALSE] %*% delta_tilde)
  
  mod3 <- lm.fit(y = lhs, x = Z)
  kappa_hat <- coef(mod3)
  eta_hat <- resid(mod3)
  
  num <- t(kappa_hat) %*% crossprod(Z) %*% kappa_hat # / (ncol(Z) - ncol(X) + 1)
  den <- mean(eta_hat^2) 
  
  stat <- c(num / den)
  pval <- pchisq(stat, df = chisq_df, lower.tail = FALSE)
  
  return(list(stat = stat,
              pval = pval))
}

#' Conditional F test using individual level data (all variables in X)
#' 
condf_test <- function(X, Z, chisq_df = NULL) {
  
  if (is.null(chisq_df)) chisq_df <- ncol(X) - 1
  
  stats <- rep(NA, ncol(X))
  pvals <- rep(NA, ncol(X))
  for (kk in 1:ncol(X)) {
    res <- conf_test_kk(kk, X, Z, chisq_df)
    stats[kk] <- res$stat
    pvals[kk] <- res$pval
  }
  return(list(stats = stats,
              pvals = pvals))
}

#' Decorrelate each column of a matrix
#' 
decorrelate <- function(X) {
  Sig <- cov(X)
  inv_cov_xy_half <- solve(sqrt_mat(Sig))
  X_new <- X %*% inv_cov_xy_half %*% diag(sqrt(diag(Sig)))
  return(X_new)
}

mv_condf_test_ss <- function(S, X, Z, cond_set = NULL, chisq_df = NULL) {
  
  if (is.null(cond_set)) cond_set <- setdiff(1:ncol(X), S)
  if (any(S %in% cond_set)) stop("S should not be in cond_set")
  if (is.null(chisq_df)) chisq_df <- ncol(Z) - length(cond_set) + 1
  
  mod1 <- lm.fit(y = X[,cond_set,drop=FALSE], x = Z)
  X_hat <- as.matrix(mod1$fitted.values)
  
  sys_eqs <- sapply(S, \(ii) {
    formula(paste0("X[, ", ii, "] ~ X_hat"))
  })
  mod2 <- systemfit::systemfit(sys_eqs, method = "SUR")
  delta_tilde <- matrix(coef(mod2), nrow = ncol(X_hat) + 1)[-1,]
  lhs <- X[,S] - X[,cond_set,drop=FALSE] %*% delta_tilde
  
  mod3 <- lm.fit(y = lhs, x = Z)
  kappa_hat <- coef(mod3)
  eta_hat <- resid(mod3)
  
  sys_eqs <- sapply(seq.int(ncol(lhs)), \(ii) {
    formula(paste0("lhs[, ", ii, "] ~ Z"))
  })
  mod3 <- systemfit::systemfit(sys_eqs, method = "SUR")
  kappa_hat <- matrix(coef(mod3), nrow = ncol(Z) + 1)[-1,]
  
  num <- t(kappa_hat) %*% crossprod(Z) %*% kappa_hat # / (ncol(Z) - ncol(X) + 1)
  den <- mean(c(eta_hat %*% chol(solve(cov_eta)))^2) # mean(eta_hat^2) 
  
  stat <- c(num / den)
  pval <- pchisq(stat, df = chisq_df, lower.tail = FALSE)
}

#' Conditional F test using summary statistics
#' 
#' @param kk The index of the exposure
#' @param gamma_X A matrix of dimension J x K, the estimated instrument-exposure coefficients (each is a simple OLS coefficient of one exposure on one instrument).
#' @param Sig_X A matrix of dimension (K x J) x (K x J) in a block form. See patel-et-al2023 supplimentary materials.
#' 
condf_stat_kk_sumstat <- function(kk, gamma_X, Sig_X, chisq_df = NULL) {
  
  J <- nrow(gamma_X) # number of instruments
  K <- ncol(gamma_X) # number of exposures
  if (is.null(chisq_df)) chisq_df <- J - K + 1
  
  # rearrange so that the k-th row-block is the first row and the k-th 
  #.  column-block is the first column
  idx_block <- (((kk - 1) * J) + 1):(kk * J)
  idx_minus_block <- setdiff(1: (J * K), idx_block)
  idx_new <- c(idx_block, idx_minus_block)
  Sig_X_new <- Sig_X[idx_new, idx_new]
  
  Q_delta <- function(dd) {
    gg <- gamma_X[, kk] - as.vector(gamma_X[, -kk] %*% matrix(dd))
    hh <- cbind(diag(J), kronecker(-diag(J), t(dd)))
    Ome <- hh %*% Sig_X_new %*% t(hh)
    as.numeric(t(gg) %*% solve(Ome, gg))
  }
  
  # for getting initial estimation
  Q_gg <- function(dd) {
    gg <- as.vector(gamma_X[, kk] - (gamma_X[, -kk] %*% matrix(dd)))
    as.numeric(t(gg) %*% gg)
  }
  # TODO: add gradient
  dd_init <- optim(rep(0, (K - 1)), fn = Q_gg, method = "BFGS")$par

  # TODO: add gradient
  stat <- optim(dd_init, fn = Q_delta, method = "BFGS")$value / chisq_df
  return(stat)
}

#' Calculates conditional F-statistic for risk factors using summary statistics
#' 
condf_stat_sumstat <- function(coef_x, cov_coef_x, cor_z, cor_x = NULL, nobs_x = NULL, 
                               adjust = FALSE, chisq_df = NULL) {
  
  J <- nrow(coef_x) # number of instruments
  K <- ncol(coef_x) # number of exposures
  if (is.null(chisq_df)) chisq_df <- J - K + 1
  
  # adjust estimates
  if (adjust) {
    if (is.null(cor_x) | is.null(nobs_x)) stop("adjust requires cor_x and nobs_x")
    res_est <- adjust_estimates_zx(coef_x, cov_coef_x, cor_z, cor_x, nobs_x)
    gamma_X <- res_est$gamma_X
    Sig_X <- res_est$Sig_X
  } else {
    gamma_X <- coef_x
    Sig_X <- cov_coef_x
  }
  
  J <- nrow(coef_x) # number of instruments
  K <- ncol(coef_x) # number of exposures
  
  stats <- sapply(1:K, \(kk) {
    condf_stat_kk_sumstat(kk, gamma_X, Sig_X, chisq_df)
  })
  
  return(stats)
}
