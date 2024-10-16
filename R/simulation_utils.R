# ---- check if a dgp satisfies assumptions (A1) and (A3) ----

#' Generate all subsets of a set.
#'
#' @details This function is used to generate all potential invariant sets.
#'
#' @param S A vector with no duplicated values representing a set.
#'
#' @return A list of all subsets of `S`.
#' 
#' @export
get_powerset <- function(S) {
  n <- length(S)
  masks <- 2^(1:n-1)
  sets <- lapply(1:2^n-1, function(u) S[bitwAnd(u, masks) != 0])
  return(sets)
}

#' Check if Assumptions A are satisfied
#' 
#' @param A A matrix of dimension nins x nexp.
#' @param B A matrix of dimension nexp x nexp.
#' @param beta_star A vector of length nexp, the true exposure-response 
#'   coefficients.
#' 
#' @return A list containing three elements, corresponding to the three 
#'   conditions in Assumptions A in Pfister and Peters 2021.
#'   
#' @export
check_assumptions_A <- function(A, B, parents, beta_star){
  
  d <- ncol(B)
  
  ## Compute C
  C <- t(solve((diag(d) - B), A))
  
  ## Verify A1
  assA1 <- qr(C[,parents])$rank == length(parents)
  
  ## Verify A2
  sets_S <- get_powerset(1:d)
  parent_ind <- sapply(sets_S, function(S) all(c(S %in% parents, parents %in% S)))
  sets_S <- sets_S[!parent_ind]
  C_PA <- C[,parents, drop=FALSE]
  assA2_vec <- rep(TRUE, length(sets_S))
  for (ii in seq_along(sets_S)) {
    S <- sets_S[[ii]]
    C_S <- C[,S,drop=FALSE]
    cond1 <- qr(C_S)$rank <= qr(C_PA)$rank
    colsp_cs <- col_space(C_S)
    colsp_pa <- col_space(C_PA)
    cond2 <- ifelse(ncol(colsp_cs) == ncol(colsp_pa), 
                    # rref makes sure to arrange the matrices correctly
                    any(pracma::rref(colsp_cs) != pracma::rref(colsp_pa)), 
                    FALSE)
    if (cond1 & cond2) {
      W <- cbind(C_S, C_PA %*% beta_star[parents])
      # message("S = ", paste0(S, collapse = " "))
      # message(paste(round(svd(C_S)$d, 3), collapse = " "))
      # message(paste(round(svd(W)$d, 3), collapse = " "))
      assA2_vec[ii] <- qr(W)$rank > qr(C_S)$rank # W not in the image of C_S
    }
  }
  assA2 <- sum(assA2_vec) == length(assA2_vec)
  
  ## Verify A3
  sets_S <- combn(1:d, length(parents), simplify=FALSE)
  parent_ind <- sapply(sets_S, function(S) sum(sort(S) == sort(parents)) == length(S))
  sets_S <- sets_S[!parent_ind]
  
  assA3_vec <- rep(FALSE, length(sets_S))
  for(k in 1:length(sets_S)){
    unionS <- union(parents, sets_S[[k]])
    assA3_vec[k] <- (qr(C[,unionS])$rank > length(parents) | qr(C[,sets_S[[k]]])$rank < length(parents))
  }
  assA3 <- sum(assA3_vec) == length(assA3_vec)
  
  return(list(assA1 = assA1,
              assA2 = assA2,
              assA3 = assA3))
  
}

#' Compute the null space of a matrix
#' 
null_space <- function (A) {
  # https://stackoverflow.com/questions/43223579/solve-homogenous-system-ax-0-for-any-m-n-matrix-a-in-r-find-null-space-basi
  m <- dim(A)[1]; n <- dim(A)[2]
  ## QR factorization and rank detection
  QR <- base::qr.default(A)
  r <- QR$rank
  ## cases 2 to 4
  if ((r < min(m, n)) || (m < n)) {
    R <- QR$qr[1:r, , drop = FALSE]
    P <- QR$pivot
    F <- R[, (r + 1):n, drop = FALSE]
    I <- base::diag(1, n - r)
    B <- -1.0 * base::backsolve(R, F, r)
    Y <- base::rbind(B, I)
    X <- Y[base::order(P), , drop = FALSE]
    return(X)
  }
  ## case 1
  return(base::matrix(0, n, 1))
}

#' Compute the column space (image) of a matrix
#' 
col_space <- function(A) {
  m <- dim(A)[1]
  n <- dim(A)[2]
  QR <- base::qr.default(A)
  r <- QR$rank
  R <- qr.R(QR)
  Q <- qr.Q(QR)
  return(Q)
}

# ---- calculate sufficient statistics from individual data ----

#' Calculate sufficient statistics given individual level data
#' 
#' @template y
#' @template z_y
#' @template x
#' @template z_x
#' 
#' @return A list containing the variance-covariance matrices.
#' 
#' @export
comp_suffstat <- function(Y, Z_Y, X, Z_X) {
  if (length(unique(nrow(X), nrow(Y), nrow(Z_X), nrow(Z_Y))) > 1) {
    stop("This function only applies to when the number of observations in the two samples are the same.")
  }
  # n <- nrow(Y)
  Z <- rbind(Z_Y, Z_X)
  # return(list(XX = t(X) %*% X/n, 
  #             YY = t(Y) %*% Y/n,
  #             ZZ = t(Z) %*% Z/nrow(Z),
  #             XY = t(X) %*% Y/n,
  #             ZX = t(Z_X) %*% X/n, 
  #             ZY = t(Z_Y) %*% Y/n))
  return(list(XX = cov(X), 
              YY = cov(Y),
              ZZ = cov(Z),
              XY = cov(X, Y),
              ZX = cov(Z_X, X), 
              ZY = cov(Z_Y, Y)))
}

# ---- calculate summary statistics from individual data ----

comp_sumstat_simple <- function(Y, X) {
  
  est <- rep(NA, ncol(X))
  sig_est <- rep(NA, ncol(X))
  for (kk in 1:ncol(X)) {
    mod <- lm(Y ~ X[,kk])
    est[kk] <- coef(mod)[2] # remove intercept
    sig_est[kk] <- sqrt(as.numeric(vcov(mod)[2,2])) # standard error
  }
  
  return(list(est = est, sig_est = sig_est))
}

#' Compute summary statistics from two non-overlapping individual level samples.
#' 
#' @template y
#' @template z_y
#' @template x
#' @template z_x
#' @param complete A boolean indicating whether to calculate the full variance-covariance 
#'   matrix of the coefficients.
#' 
#' @return A list containing 
#' - coef_x: A matrix of dimension J x K, the estimated instrument-exposure coefficients. 
#' - coef_y: A vector of length J, the estimated instrument-response coefficients.
#' - cov_coef_x: if indep, a list of length K of matrices each of dimension J x J, 
#'     the covariance matrices of the estimated instrument-exposure coefficients; 
#'     if not indep, a matrix of dimension (J x K) x (J x K)
#' - cov_coef_y: if indep, a vector of length J, the standard errors of the 
#'     estimated instrument-response coefficients; if not indep, a matrix of 
#'     dimension J x J, the covariance matrix of the estimated instrument-response 
#'     coefficients. 
#' - cor_x: A matrix of dimension K x K, the correlation matrix of X.
#' - cor_z: A matrix of dimension J x J, the correlation matrix of Z.
#' 
#' @export 
comp_sumstat <- function(Y, Z_Y, X, Z_X, 
                         type = c("depen_mv", "indep_mv", "simple")) {
  
  type <- match.arg(type)
  
  if (type == "simple") {
    
    res_x <- lapply(1:ncol(X), \(kk) {
      comp_sumstat_simple(X[,kk], Z_X)
    })
    gamma_exp <- sapply(res_x, \(res) {
      res$est
    }, simplify = TRUE)
    colnames(gamma_exp) <- paste0("X", 1:ncol(X))
    rownames(gamma_exp) <- paste0("Z", 1:ncol(Z_X))
    sig_exp <- sapply(res_x, \(res) {
      res$sig_est
    }, simplify = TRUE)
    colnames(sig_exp) <- paste0("X", 1:ncol(X))
    rownames(sig_exp) <- paste0("Z", 1:ncol(Z_X))

    res_y <- comp_sumstat_simple(Y, Z_Y)
    gamma_res <- res_y$est
    names(gamma_res) <- paste0("Z", 1:ncol(Z_Y))
    sig_res <- res_y$sig_est
    names(sig_res) <- "Y"
    
    return(list(coef_x = gamma_exp,
                se_coef_x = sig_exp,
                coef_y = gamma_res,
                se_coef_y = sig_res,
                cor_x = cor(X),
                cor_zinx = cor(Z_X),
                cor_ziny = cor(Z_Y),
                cov_x = cov(X),
                cov_zinx = cov(Z_X),
                cov_ziny = cov(Z_Y)))
    
  } else {
    
    # regression on all instruments
    mod1 <- lm(Y ~ Z_Y)
    sys_eqs <- sapply(1:ncol(X), \(ii) {
      # regression X[,ii] on Z_X
      formula(paste0("X[, ", ii, "] ~ Z_X"))
    })
    mod2 <- systemfit::systemfit(sys_eqs, method = "SUR")
    
    # instrument-response coefficients 
    gamma_res <- coef(mod1)[-c(1)] # remove intercept
    Y_names <- paste0("Y_Z", 1:ncol(Z_Y))
    names(gamma_res) <- Y_names
    
    # complete instrument-response coefficients covariances
    sig_res_all <- vcov(mod1)[-c(1), -c(1)]
    
    # instrument-exposure coefficients
    idx_intercept <- seq(1, ncol(X) * (ncol(Z_X) + 1), by = ncol(Z_X) + 1)
    gamma_exp <- matrix(coef(mod2)[-idx_intercept], nrow = ncol(Z_X))
    colnames(gamma_exp) <- paste0("X", 1:ncol(X)) # K exposures
    rownames(gamma_exp) <- paste0("Z", 1:ncol(Z_X)) # J instruments
    
    # complete instrument-exposure coefficients covariances
    sig_exp_all <- vcov(mod2)[-idx_intercept, -idx_intercept]
    
    if (type == "depen_mv") {
      sig_res <- sig_res_all
      colnames(sig_res) <- Y_names
      rownames(sig_res) <- Y_names
      gamma_exp <- gamma_exp
      sig_exp <- sig_exp_all
      
      return(list(coef_x = gamma_exp,
                  cov_coef_x = sig_exp,
                  coef_y = gamma_res,
                  cov_coef_y = sig_res,
                  cor_x = cor(X),
                  cor_zinx = cor(Z_X),
                  cor_ziny = cor(Z_Y),
                  cov_x = cov(X),
                  cov_zinx = cov(Z_X),
                  cov_ziny = cov(Z_Y)))
      
    } else if (type == "indep_mv") {
      sig_res <- sqrt(diag(sig_res_all)) # standard error
      names(sig_res) <- Y_names
      sig_exp <- lapply(1:ncol(Z_X), \(ii) {
        idx <- seq(ii, nrow(sig_exp_all), by = ncol(Z_X))
        # take only the diag
        sig_exp_ii <- sig_exp_all[idx,idx,drop=FALSE]
        colnames(sig_exp_ii) <- paste0("X", 1:ncol(X), "_Z", ii)
        rownames(sig_exp_ii) <- paste0("X", 1:ncol(X), "_Z", ii)
        sig_exp_ii
      })
      
      return(list(coef_x = gamma_exp,
                  cov_coef_x = sig_exp,
                  coef_y = gamma_res,
                  se_coef_y = sig_res, # standard error
                  cor_x = cor(X),
                  cor_zinx = cor(Z_X),
                  cor_ziny = cor(Z_Y),
                  cov_x = cov(X),
                  cov_zinx = cov(Z_X),
                  cov_ziny = cov(Z_Y)))
      
    } else {
      stop("type name not found!")
    }
  }
}


