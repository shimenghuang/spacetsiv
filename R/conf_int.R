
#' Confidence intervals by normal approximation
#' 
#' @export
normal_approx_ci <- function(beta_hat, second_deriv, alpha = 0.05) {
  se <- 1 / sqrt(second_deriv)
  zz <- qnorm(alpha/2, lower.tail = FALSE)
  return(list(lwr = beta_hat - zz * se,
              upr = beta_hat + zz * se,
              est = beta_hat))
}

#' Note: this implement w.r.t. ||Gamma - gamma %*% beta||^2 not general Q stat
#'   since the first derivative of the general Q stat cannot be written 
#'   as a single sum of J terms. It seems that the sandwich estimator is 
#'   only studied under (quasi-)likelihood and generalized estimating 
#'   equations (GEEs) which both can be written as a single sum of 
#'   individual terms without interactions. Also not sure if the 
#'   sandwich brings more efficiency than using the variance of the 
#'   GMM estimator already. 
#'   
#' @export
sandwich_ci <- function(coef_y, cov_coef_y, coef_x, cov_coef_x, beta_hat,
                        alpha = 0.05) {
  A <- t(coef_x) %*% coef_x
  A_inv <- solve(A)
  B <- sapply(1:length(coef_y), \(jj) {
    -2 * c(coef_y[jj] - t(c(coef_x[jj,,drop=FALSE])) %*% beta_hat) * coef_x[jj,,drop=FALSE]
  })
  B <- matrix(B, ncol = length(coef_y))
  B <- B %*% t(B)
  vv <- A_inv %*% B %*% A_inv
  se <- sqrt(diag(vv))
  lwr <- beta_hat - se * qnorm(1-alpha/2, lower.tail = TRUE)
  upr <- beta_hat + se * qnorm(1-alpha/2, lower.tail = TRUE)
  return(data.frame(lwr = lwr, upr = upr, est = beta_hat))
}

q_test_gen_ci <- function(coef_y, cov_coef_y, coef_x, cov_coef_x, beta_hat,
                          alpha = 0.05) {
  phi <- kronecker(t(beta_hat), diag(nrow(coef_x)))
  Ome <- cov_coef_y + phi %*% cov_coef_x %*% t(phi)
  vv <- solve(t(coef_x) %*% solve(Ome, coef_x))
  se <- sqrt(diag(vv))
  lwr <- beta_hat - se * qnorm(1-alpha/2, lower.tail = TRUE)
  upr <- beta_hat + se * qnorm(1-alpha/2, lower.tail = TRUE)
  return(data.frame(lwr = lwr, upr = upr, est = beta_hat))
}


#' Subvector confidence interval with general Q test based on projection
#' 
#' @export
q_test_gen_subvec_ci_norefit <- function(coef_y, cov_coef_y, coef_x, cov_coef_x, 
                                         beta_liml, subvec_idx = NULL, 
                                         search_lwr = 1, search_upr = 1,  
                                         extend = TRUE, 
                                         alpha = 0.05, chisq_df = NULL, tol = 1e-5) {
  if (is.null(subvec_idx)) {
    subvec_idx <- seq_along(beta_liml)
  }
  if (is.null(chisq_df)) {
    chisq_df <- length(coef_y)
  }
  
  res <- lapply(subvec_idx, \(jj) {
    lwr <- tryCatch({
      uniroot(\(bb) {
        beta_hat <- beta_liml
        beta_hat[jj] <- bb
        q_test_gen(coef_y = coef_y,
                   cov_coef_y = cov_coef_y,
                   coef_x = coef_x,
                   cov_coef_x = cov_coef_x,
                   beta = beta_hat,
                   chisq_df = chisq_df,
                   tol = tol)$stat - qchisq(alpha/2, df = chisq_df, lower.tail = FALSE)
      }, 
      lower = beta_liml[jj] - search_lwr, upper = beta_liml[jj], 
      extendInt = ifelse(extend, "downX", "no"), tol = tol)$root
    }, error = \(e) -Inf)
    upr <- tryCatch({uniroot(\(bb) {
      beta_hat <- beta_liml
      beta_hat[jj] <- bb
      q_test_gen(coef_y = coef_y,
                 cov_coef_y = cov_coef_y,
                 coef_x = coef_x,
                 cov_coef_x = cov_coef_x,
                 beta = beta_hat,
                 chisq_df = chisq_df,
                 tol = tol)$stat - qchisq(alpha/2, df = chisq_df, lower.tail = FALSE)
    }, 
    lower = beta_liml[jj], upper = beta_liml[jj] + search_upr, 
    extendInt = ifelse(extend, "upX", "no"), tol = tol)$root
    }, error = \(e) Inf)
    list(lwr = lwr, upr = upr)
  })
  res <- do.call(rbind, res)
  res <- data.frame(res)
  res$lwr <- unlist(res$lwr)
  res$upr <- unlist(res$upr)
  rownames(res) <- paste0("beta", subvec_idx)
  res$est <- beta_liml[subvec_idx]
  return(res)
}

#' Subvector confidence interval with general Q test based on projection
#' 
#' 
#' @export
q_test_gen_subvec_ci <- function(coef_y, cov_coef_y, coef_x, cov_coef_x, 
                                 beta_liml, subvec_idx = NULL, 
                                 search_lwr = 1, search_upr = 1,  
                                 extend = TRUE, 
                                 alpha = 0.05, chisq_df = NULL, tol = 1e-5,
                                 uniroot_maxiter = 500, 
                                 optim_method = "L-BFGS-B",
                                 control = list(maxit = 500, 
                                                factr = 1e6),
                                 ...) {
  if (is.null(subvec_idx)) {
    subvec_idx <- seq_along(beta_liml)
  }
  if (is.null(chisq_df)) {
    chisq_df <- length(coef_y)
  }
  
  refit_rest <- function(idx_fix, beta_fix) {
    idx_rest <- setdiff(subvec_idx, idx_fix)
    res <- optim(beta_liml[idx_rest], fn = \(bb) {
      beta_hat <- beta_liml
      beta_hat[idx_rest] <- bb
      beta_hat[idx_fix] <- beta_fix
      q_test_gen(coef_y = coef_y,
                 cov_coef_y = cov_coef_y,
                 coef_x = coef_x,
                 cov_coef_x = cov_coef_x,
                 beta = beta_hat,
                 chisq_df = chisq_df,
                 tol = tol)$stat
    },
    method = optim_method,
    hessian = FALSE,
    control = control, ...)
    beta_new <- beta_liml
    beta_new[idx_rest] <- res$par
    beta_new[idx_fix] <- beta_fix
    return(beta_new)
  }
  
  res <- lapply(subvec_idx, \(jj) {
    lwr <- tryCatch({
      uniroot(\(bb) {
        beta_hat <- refit_rest(jj, bb)
        q_test_gen(coef_y = coef_y,
                   cov_coef_y = cov_coef_y,
                   coef_x = coef_x,
                   cov_coef_x = cov_coef_x,
                   beta = beta_hat,
                   chisq_df = chisq_df,
                   tol = tol)$stat - qchisq(alpha/2, df = chisq_df, lower.tail = FALSE)
      }, 
      lower = beta_liml[jj] - search_lwr, upper = beta_liml[jj], 
      extendInt = ifelse(extend, "downX", "no"), tol = tol,
      maxiter = uniroot_maxiter)$root
    }, error = \(e) -Inf)
    upr <- tryCatch({
      uniroot(\(bb) {
        beta_hat <- refit_rest(jj, bb)
        q_test_gen(coef_y = coef_y,
                   cov_coef_y = cov_coef_y,
                   coef_x = coef_x,
                   cov_coef_x = cov_coef_x,
                   beta = beta_hat,
                   chisq_df = chisq_df,
                   tol = tol)$stat - qchisq(alpha/2, df = chisq_df, lower.tail = FALSE)
      }, 
    lower = beta_liml[jj], upper = beta_liml[jj] + search_upr, 
    extendInt = ifelse(extend, "upX", "no"), tol = tol,
    maxiter = uniroot_maxiter)$root
    }, error = \(e) Inf)
    list(lwr = lwr, upr = upr)
  })
  res <- do.call(rbind, res)
  res <- data.frame(res)
  res$lwr <- unlist(res$lwr)
  res$upr <- unlist(res$upr)
  rownames(res) <- paste0("beta", subvec_idx)
  res$est <- beta_liml[subvec_idx]
  return(res)
}

