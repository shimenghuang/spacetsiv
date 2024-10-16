# ---- examples ----

#' Data generating process 00 (2 valid & indep. IVs + 2 cond. indep. exposures)
#' 
#' @details If beta_star does not contain 0, this is a dgp with 2 valid and  
#'   independent instruments, 2 exposures, and 2 independent hidden confounders 
#'   between exposures and outcome.
#' 
#' @param beta_star A numeric vector of length 2.
#' 
#' @export
dgp_ex00 <- function(n, beta_star = c(1, 2), seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  Z <- matrix(rnorm(n * 2), nrow = n)
  H <- matrix(rnorm(n * 2, 0.5), nrow = n)
  X <- Z %*% cbind(c(1, 0.15), c(0.25, 1)) + H + rnorm(n * 2)
  Y <- X %*% beta_star + H %*% c(1, 0.5) + rnorm(n)
  return(list(Z = Z,
              X = X,
              Y = Y,
              beta_star = beta_star))
}

#' Data generating process 0 (valid & over-identified & correlated IVs)
#' 
#' @details An over-identified valid IV example with 4 correlated instruments 
#'   and 2 exposures.
#' 
#' @param beta_star A numeric vector of length 2.
#' 
#' @export
dgp_ex0 <- function(n, beta_star, seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  d <- 2
  m <- 4
  Sig <- matrix(c( 1.75,     1,  -0.2,   0.5,
                   1,   1.5,  0.25, -0.25,
                   -0.2,  0.25,  0.75,  0.25,
                   0.5, -0.25, 0.25,     1), nrow = 4)
  Z <- matrix(rnorm(n * 4, 0.1), nrow = n) %*% chol(Sig)
  # A <- rbind(c(1, 1.5, 0.1, 0.2),
  #            c(0.2, 0.1, 0.75, 1)) # A is 2 x 4
  A <- rbind(c(1, 1.5, 0.5, 0.2),
             c(0.5, 0.1, 0.75, 1)) # A is 2 x 4
  B <- matrix(0, 2, 2)
  B[1,1] <- 2
  Sig_H <- matrix(c( 1.75,     1, 
                     1,   1.5), nrow = 2)
  H <- matrix(rnorm(n * d, 0, 1), n, d)  %*% chol(Sig_H)
  Id <- diag(d)
  X <- (Z %*% t(A) + H + rnorm(n * 2)) %*% t(solve(Id-B))
  Y <- X %*% beta_star + H %*% rep(1, d) + rnorm(n)
  return(list(Z = Z,
              X = X,
              Y = Y,
              A = A,
              B = B,
              parents = which(beta_star != 0),
              beta_star = beta_star))
}

#' Data generating process with 20 correlated instruments and 5 exposures.
#' 
#' @details If d1 and/or d2 are non-zero, the corresponding instruments are 
#'   invalid. 
#' 
#' @export
dgp_ex01 <- function(n, beta1, beta2, d1 = 1, d2 = 1, seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  # parameters
  m <- 20
  d <- 5
  B <- matrix(0, d, d)
  B[1,2] <- 1
  B[1,3] <- 1
  B[2,3] <- 1
  B[2,5] <- 1
  # B <- B * matrix(sample(c(-1,1), d^2, replace = TRUE, prob = c(0.5, 0.5)) * runif(d^2, 0.5, 1.5), d, d)
  # B <- t(apply(B, 1, function(x) x/max(c(sum(x != 0), 1))))
  A <- cbind(diag(d), diag(d), diag(d), diag(d)) # A is d x m
  beta_star <- c(beta1, beta2, rep(0, 3))
  Id <- diag(d)
  # C <- t(solve(Id-B, A))
  
  # variables
  H <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_x <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
  X <- (Z %*% t(A) + H + eps_x) %*% t(solve(Id-B))
  Y <- X %*% beta_star + Z[,c(1,2)] %*% c(d1, d2) + H %*% rep(1, d) + eps_y
  
  return(list(Z = Z,
              X = X,
              Y = Y,
              A = A,
              B = B,
              parents = which(beta_star != 0),
              beta_star = beta_star))
}

#' Data generating process with 100 correlated instruments and 5 exposures.
#' 
#' @details If d1 and/or d2 are non-zero, the corresponding instruments are 
#'   invalid. 
#' 
#' @export
dgp_ex02 <- function(n, beta1, beta2, ld, d1 = 1, d2 = 1, seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  # parameters
  m <- 100
  d <- 5
  B <- matrix(0, d, d)
  B[1,2] <- 1
  B[1,3] <- 1
  B[2,4] <- 1
  # B[2,3] <- 1
  B[2,5] <- 1
  # B <- B * matrix(sample(c(-1,1), d^2, replace = TRUE, prob = c(0.5, 0.5)) * runif(d^2, 0.5, 1.5), d, d)
  # B <- t(apply(B, 1, function(x) x/max(c(sum(x != 0), 1))))
  A <- matrix(0, nrow = d, ncol = m)
  A[1:d, 1:d] <- diag(c(1, 0, 1, 0, 1))
  A[, (d+1):(2*d)] <- diag(c(0, 1, 0, 1, 0))
  A[, (2*d+1):(3*d)] <- diag(d)
  A[, (3*d+1):(4*d)] <- diag(d)
  beta_star <- c(beta1, beta2, rep(0, 3))
  Id <- diag(d)
  # C <- t(solve(Id-B, A))
  
  # variables
  H <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_x <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_y <- rnorm(n, 0, 1)
  Z <- mvtnorm::rmvnorm(n, sigma = (ld + t(ld))/2)
  X <- (Z %*% t(A) + H + eps_x) %*% t(solve(Id-B))
  Y <- X %*% beta_star + Z[,c(1,2)] %*% c(d1, d2) + H %*% rep(1, d) + eps_y
  
  return(list(Z = Z,
              X = X,
              Y = Y,
              A = A,
              B = B,
              parents = which(beta_star != 0),
              beta_star = beta_star))
}

#' 
#' @export
dgp_ex03 <- function(n, beta_star, seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  d <- 100 # dimension of X
  m <- 5 # dimension of Z
  
  A <- matrix(0, nrow = d, ncol = m)
  A[1,1] <- 1
  A[2,1] <- 1
  A[2,2] <- 1
  A[3,2] <- 1
  A[3,3] <- 1
  A[4,3] <- 1
  A[5,4] <- 1
  B <- matrix(0, d, d)
  B[3, 1] <- 1
  Id <- diag(d)
  
  # variables
  H <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_x <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
  X <- (Z %*% t(A) + H + eps_x) %*% t(solve(Id-B))
  Y <- X %*% beta_star + H %*% rep(1, d) + eps_y
  
  return(list(Z = Z,
              X = X,
              Y = Y,
              A = A,
              B = B,
              parents = which(beta_star != 0),
              beta_star = beta_star))
}

#' Data generating process (summary statistics) with high dimensional exposure.
#' 
#' @export
dgp_ex03_sumstat <- function(ny, nx, beta_star, 
                             type = c("indep_mv", "depen_mv"), 
                             seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  d <- 100 # dimension of X
  m <- 5 # dimension of Z
  var_z <- diag(m) 
  var_z[1, 2] <- 0.1
  var_z[1, 3] <- -0.2
  var_z[1, 4] <- 0.15
  var_z[1, 5] <- 0.05
  var_z <- (var_z + t(var_z))/2
  A <- matrix(0, nrow = d, ncol = m)
  A[1,1] <- 1
  A[2,1] <- 1
  A[2,2] <- 1
  A[3,2] <- 1
  A[3,3] <- 1
  A[4,3] <- 1
  A[5,4] <- 1
  B <- matrix(0, d, d)
  # B[10, 6] <- 1
  B[3, 1] <- 1
  Id <- diag(d)
  W <- Id - B
  gamma <- t(solve(W, A))
  Gamma <- as.numeric(gamma %*% beta_star)
  # var_xix <- diag(d)
  var_xix <- matrix(runif(d^2, min = -0.3, max = 0.5), nrow = d)
  var_xix <- var_xix %*% t(var_xix) + diag(d)
  var_xiy <- 1
  cov_xixy <- sample(c(0.2, 0.4, 0.6, 0.8), 100, replace = TRUE)
  Sig_gamma <- kronecker(t(solve(W, var_xix)), solve(var_z))
  Sig_gamma <- (Sig_gamma + t(Sig_gamma)) / 2 / nx # ensure it's symmetric
  bs_winv <- solve(W, beta_star)
  Sig_Gamma <- t(bs_winv) %*% var_xix %*% bs_winv + var_xiy + 2 * t(bs_winv) %*% cov_xixy
  Sig_Gamma <- as.numeric(Sig_Gamma) * solve(var_z)
  Sig_Gamma <- (Sig_Gamma + t(Sig_Gamma)) / 2 / ny 
  gamma_hat <- matrix(mvtnorm::rmvnorm(1, 
                                       mean = as.numeric(gamma), 
                                       sigma = Sig_gamma), 
                      nrow = m, ncol = d)
  Gamma_hat <- mvtnorm::rmvnorm(1, mean = Gamma, sigma = Sig_Gamma)
  Gamma_hat <- as.numeric(Gamma_hat) # convert to a vector
  
  if (type == "indep_mv") {
    sig_gamma <- lapply(1:m, \(ii) {
      idx <- seq(ii, nrow(Sig_gamma), by = m)
      # take only the diag blocks
      sig_exp_ii <- Sig_gamma[idx,idx,drop=FALSE]
      colnames(sig_exp_ii) <- paste0("X", 1:d, "_Z", ii)
      rownames(sig_exp_ii) <- paste0("X", 1:d, "_Z", ii)
      sig_exp_ii
    })
    return(list(coef_x = gamma_hat,
                coef_y = Gamma_hat, 
                cov_coef_x = sig_gamma,
                cov_coef_y = sqrt(diag(Sig_Gamma)),
                cov_zinx = var_z,
                cov_ziny = var_z,
                gamma = gamma,
                Gamma = Gamma, 
                parents = which(beta_star != 0),
                beta_star = beta_star))
  } else if (type == "depen_mv") {
    return(list(coef_x = gamma_hat,
                coef_y = Gamma_hat, 
                cov_coef_x = Sig_gamma,
                cov_coef_y = Sig_Gamma,
                cov_zinx = var_z,
                cov_ziny = var_z,
                gamma = gamma,
                Gamma = Gamma, 
                parents = which(beta_star != 0),
                beta_star = beta_star))
  } else {
    stop("Not implemented.")
  }
}

#' Data generating process in Example 1 of spaceIV
#' 
#' @details Example 1 in spaceIV with three hidden variables each affecting one 
#'   exposure as well as the outcome. A DGP with 2 instruments and 3 exposures. 
#'   
#' @export
dgp_ex1 <- function(n, beta2, seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  # paramters
  d <- 3
  m <- 2
  B <- matrix(0, d, d)
  B[2,1] <- 1
  colnames(B) <- paste0("X", 1:3)
  rownames(B) <- paste0("X", 1:3)
  # B <- B * matrix(sample(c(-1,1), d^2, replace = TRUE, prob = c(0.5, 0.5)) * runif(d^2, 0.5, 1.5), d, d)
  # B <- t(apply(B, 1, function(x) x/max(c(sum(x != 0), 1))))
  A <- matrix(0, nrow = d, ncol = m)
  A[1,1] <- 1
  A[2,1] <- 1
  A[2,2] <- 1
  A[3,2] <- 1
  colnames(A) <- paste0("Z", 1:2)
  rownames(A) <- paste0("X", 1:3)
  beta_star <- c(0, beta2, 0)
  Id <- diag(d)
  # C <- t(solve(Id-B, A))
  
  # variables
  H <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_x <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n * m, 0, 1), n, m)
  X <- (Z %*% t(A) + H + eps_x) %*% t(solve(Id-B))
  Y <- X %*% beta_star + H %*% rep(1, d) + eps_y
  
  return(list(Z = Z,
              X = X,
              Y = Y,
              A = A,
              B = B,
              # C = C,
              parents = which(beta_star != 0),
              beta_star = beta_star))
}

#' Data generating process in Example 4 of spaceIV
#' 
#' @details Example 4 in spaceIV with three hidden variables each affecting one 
#'   exposure as well as the outcome. A DGP with 3 instruments and 2 exposures.
#'   
#' @export
dgp_ex4 <- function(n, beta1, beta2, seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  # paramters
  d <- 5
  m <- 3
  B <- matrix(0, d, d)
  B[1,3] <- 1
  B[1,4] <- 1
  B[2,4] <- 1
  B[2,5] <- 1
  # B <- B * rnorm(d * d, 1, 1) # use random normal?
  # B <- B * matrix(sample(c(-1,1), d^2, replace = TRUE, prob = c(0.5, 0.5)) *
  #                   runif(d^2, 0.5, 1.5), d, d)
  # B <- t(apply(B, 1, function(x) x/max(c(sum(x != 0), 1))))
  A <- matrix(0, nrow = d, ncol = m)
  A[3,1] <- 1
  A[4,2] <- 1
  A[5,3] <- 1
  beta_star <- c(beta1, beta2, 0, 0, 0)
  Id <- diag(d)
  
  # variables
  H <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_x <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n * m, 0, 1), n, m)
  X <- (Z %*% t(A) + H + eps_x) %*% t(solve(Id-B))
  Y <- X %*% beta_star + H %*% rep(1, d) + eps_y
  
  return(list(Z = Z,
              X = X,
              Y = Y,
              A = A,
              B = B,
              # C = C,
              parents = which(beta_star != 0),
              beta_star = beta_star)) 
}

#' Data generating process 2
#' @details DGP with 5 exposure and 5 instruments but 2 instruments 
#'   are invalid; 2 of the 5 exposure have a causal effect on the outcome.
#'   
#' @export
dgp_with_invalid <- function(n, beta1, beta2, d1 = 1, d2 = 2, seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  # parameters
  d <- 5
  m <- 5
  B <- matrix(0, d, d)
  B[1,3] <- 1
  B[1,4] <- 1
  B[2,4] <- 1
  B[2,5] <- 1
  # B <- B * matrix(sample(c(-1,1), d^2, replace = TRUE, prob = c(0.5, 0.5)) * runif(d^2, 0.5, 1.5), d, d)
  # B <- t(apply(B, 1, function(x) x/max(c(sum(x != 0), 1))))
  A <- diag(m)
  beta_star <- c(beta1, beta2, 0, 0, 0)
  Id <- diag(d)
  # C <- t(solve(Id-B, A))
  
  # variables
  H <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_x <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
  X <- (Z %*% t(A) + H + eps_x) %*% t(solve(Id-B))
  Y <- X %*% beta_star + Z[,c(1,2)] %*% c(d1, d2) + H %*% rep(1, d) + eps_y
  
  return(list(Z = Z,
              X = X,
              Y = Y,
              A = A,
              B = B,
              parents = which(beta_star != 0),
              beta_star = beta_star))
}

#' Data generating process 3
#' 
#' @details DGP with 5 exposure and 5 instruments which are all 
#'   invalid, and 2 of the 5 exposure have a causal effect on the outcome.
#'   
#' @export
dgp_all_valid <- function(n, beta1, beta2, seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  # parameters
  d <- 5
  m <- 5
  B <- matrix(0, d, d)
  B[1,3] <- 1
  B[1,4] <- 1
  B[2,4] <- 1
  B[2,5] <- 1
  # B <- B * matrix(sample(c(-1,1), d^2, replace = TRUE, prob = c(0.5, 0.5)) * runif(d^2, 0.5, 1.5), d, d)
  # B <- t(apply(B, 1, function(x) x/max(c(sum(x != 0), 1))))
  A <- diag(m)
  beta_star <- c(beta1, beta2, 0, 0, 0)
  Id <- diag(d)
  # C <- t(solve(Id-B, A))
  
  # variables
  H <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_x <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
  X <- (Z %*% t(A) + H + eps_x) %*% t(solve(Id-B))
  Y <- X %*% beta_star + H %*% rep(1, d) + eps_y
  
  return(list(Z = Z,
              X = X,
              Y = Y,
              A = A,
              B = B,
              # C = C,
              parents = which(beta_star != 0),
              beta_star = beta_star))
}

# ---- non-idenfitied example ----

dgp_noniden1 <- function(n, beta3, beta4, seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  # paramters
  d <- 4
  m <- 2
  B <- matrix(0, d, d)
  B[c(3,4), c(1,2)] <- sample(c(0.5, 1, 1.5, 0.75), replace = FALSE)
  # B[3, 1] <- 0.5
  # B[4, 2] <- 1.5
  colnames(B) <- paste0("X", 1:d)
  rownames(B) <- paste0("X", 1:d)
  # B <- B * matrix(sample(c(-1,1), d^2, replace = TRUE, prob = c(0.5, 0.5)) * runif(d^2, 0.5, 1.5), d, d)
  # B <- t(apply(B, 1, function(x) x/max(c(sum(x != 0), 1))))
  A <- matrix(0, nrow = d, ncol = m)
  A[1,1] <- 1
  A[2,2] <- 2
  colnames(A) <- paste0("Z", 1:m)
  rownames(A) <- paste0("X", 1:d)
  beta_star <- c(0, 0, beta3, beta4)
  Id <- diag(d)
  # C <- t(solve(Id-B, A))
  
  # variables
  H <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_x <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n * m, 0, 1), n, m)
  # X <- (Z %*% t(A) + eps_x) %*% t(solve(Id-B))
  X <- (Z %*% t(A) + H + eps_x) %*% t(solve(Id-B))
  Y <- X %*% beta_star + H %*% rep(1, d) + eps_y
  
  return(list(Z = Z,
              X = X,
              Y = Y,
              A = A,
              B = B,
              # C = C,
              parents = which(beta_star != 0),
              beta_star = beta_star))
  
}

# ---- random graphs ----

linear_DAG_model_sparse <- function(d = 50, m = 10, n = 100, 
                                    p_connect = NULL, 
                                    sparsity = 3,
                                    noise_sd = 1, 
                                    noise_H = 2, 
                                    noise_I = 2) {
  
  if (is.null(p_connect)) p_connect <- 2/d
  
  ## Generate SCM
  B <- matrix(0, d, d)
  causalOrder <- sample(1:d)
  combinations <- combn(d, 2, simplify = FALSE)
  edges <- combinations[sample(c(TRUE, FALSE), length(combinations),
                               replace=TRUE, prob=c(p_connect, 1-p_connect))]
  if(length(edges) == 0){
    edges <- combinations[sample(1:length(combinations), 1)]
  }
  for(i in 1:length(edges)){
    e1 <- edges[[i]][1]
    e2 <- edges[[i]][2]
    if(which(causalOrder == e1) <= which(causalOrder == e2)){
      B[e2, e1] <- 1
    }
    else{
      B[e1, e2] <- 1
    }
  }
  B <- B * matrix(sample(c(-1,1), d^2, replace = TRUE, prob = c(0.5, 0.5)) * runif(d^2, 0.5, 1.5), d, d)
  
  # Normalize the weights to avoid aggregation of variance
  B <- t(apply(B, 1, function(x) x/max(c(sum(x != 0), 1))))
  
  A <- matrix(0, d, m)
  for(i in 1:d){
    A[i,] <- sample(c(0, 1), prob=c(1-2/d, 2/d), m, replace=TRUE)
  }
  diag(A) <- 1
  
  beta_star <- matrix(0, d, 1)
  parents <- sample(1:d, sparsity, replace=FALSE)
  beta_star[parents] <- rep(1, sparsity)
  
  ## Generate data
  Z <- matrix(rnorm(n*m, 0, noise_I), m, n)
  conf_struct <- matrix(rep(1,d), ncol=1)
  
  # Sample data
  H <- matrix(rnorm(n, 0, noise_H), 1, n)
  X <- matrix(NA, n, d)
  eps <- matrix(rnorm(n*(d+1), 0, noise_sd), d+1, n)
  X <- t(solve(diag(d)-B, eps[1:d,] + A %*% Z + conf_struct %*% H))
  Y <- X %*% beta_star - 2*t(H) + eps[d+1,]
  
  ## Return results
  return(list(X = X,
              Y = Y,
              Z = t(Z),
              A = A,
              B = B,
              parents = parents,
              beta_star = beta_star))
}

# ---- dgps in Huang et al. 2024 ----

#' Data generating process in Example 4 of spaceIV
#' 
#' @details Example 4 in spaceIV with three hidden variables each affecting one 
#'   exposure as well as the outcome. A DGP with 3 instruments and 2 exposures.
dgp1 <- function(n, beta1, beta2, seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  # paramters
  d <- 5
  m <- 3
  B <- matrix(0, d, d)
  B[1,3] <- 1
  B[1,4] <- 1
  B[2,4] <- 1
  B[2,5] <- 1
  A <- matrix(0, nrow = d, ncol = m)
  A[3,1] <- 1
  A[4,2] <- 1
  A[5,3] <- 1
  beta_star <- c(beta1, beta2, 0, 0, 0)
  Id <- diag(d)
  
  # variables
  H <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_x <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n * m, 0, 1), n, m)
  X <- (Z %*% t(A) + H + eps_x) %*% t(solve(Id-B))
  Y <- X %*% beta_star + H %*% rep(1, d) + eps_y
  
  return(list(Z = Z,
              X = X,
              Y = Y,
              A = A,
              B = B,
              parents = which(beta_star != 0),
              beta_star = beta_star)) 
}

#' Data generating process (summary statistics) with high dimensional exposure.
#' 
dgp2 <- function(ny, nx, beta_star, 
                 type = c("indep_mv", "depen_mv"), 
                 seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  d <- 100 # dimension of X
  m <- 5 # dimension of Z
  var_z <- diag(m) 
  var_z[1, 2] <- 0.1
  var_z[1, 3] <- -0.2
  var_z[1, 4] <- 0.15
  var_z[1, 5] <- 0.05
  var_z <- (var_z + t(var_z))/2
  A <- matrix(0, nrow = d, ncol = m)
  A[1,1] <- 1
  A[2,1] <- 1
  A[2,2] <- 1
  A[3,2] <- 1
  A[3,3] <- 1
  A[4,3] <- 1
  A[5,4] <- 1
  B <- matrix(0, d, d)
  B[3, 1] <- 1
  Id <- diag(d)
  W <- Id - B
  gamma <- t(solve(W, A))
  Gamma <- as.numeric(gamma %*% beta_star)
  var_xix <- matrix(runif(d^2, min = -0.3, max = 0.5), nrow = d)
  var_xix <- var_xix %*% t(var_xix) + diag(d)
  var_xiy <- 1
  cov_xixy <- sample(c(0.2, 0.4, 0.6, 0.8), 100, replace = TRUE)
  Sig_gamma <- kronecker(t(solve(W, var_xix)), solve(var_z))
  Sig_gamma <- (Sig_gamma + t(Sig_gamma)) / 2 / nx # ensure it's symmetric
  bs_winv <- solve(W, beta_star)
  Sig_Gamma <- t(bs_winv) %*% var_xix %*% bs_winv + var_xiy + 2 * t(bs_winv) %*% cov_xixy
  Sig_Gamma <- as.numeric(Sig_Gamma) * solve(var_z)
  Sig_Gamma <- (Sig_Gamma + t(Sig_Gamma)) / 2 / ny 
  gamma_hat <- matrix(mvtnorm::rmvnorm(1, 
                                       mean = as.numeric(gamma), 
                                       sigma = Sig_gamma), 
                      nrow = m, ncol = d)
  Gamma_hat <- mvtnorm::rmvnorm(1, mean = Gamma, sigma = Sig_Gamma)
  Gamma_hat <- as.numeric(Gamma_hat) # convert to a vector
  
  if (type == "indep_mv") {
    sig_gamma <- lapply(1:m, \(ii) {
      idx <- seq(ii, nrow(Sig_gamma), by = m)
      # take only the diag blocks
      sig_exp_ii <- Sig_gamma[idx,idx,drop=FALSE]
      colnames(sig_exp_ii) <- paste0("X", 1:d, "_Z", ii)
      rownames(sig_exp_ii) <- paste0("X", 1:d, "_Z", ii)
      sig_exp_ii
    })
    return(list(coef_x = gamma_hat,
                coef_y = Gamma_hat, 
                cov_coef_x = sig_gamma,
                cov_coef_y = sqrt(diag(Sig_Gamma)),
                cov_zinx = var_z,
                cov_ziny = var_z,
                gamma = gamma,
                Gamma = Gamma, 
                parents = which(beta_star != 0),
                beta_star = beta_star))
  } else if (type == "depen_mv") {
    return(list(coef_x = gamma_hat,
                coef_y = Gamma_hat, 
                cov_coef_x = Sig_gamma,
                cov_coef_y = Sig_Gamma,
                cov_zinx = var_z,
                cov_ziny = var_z,
                gamma = gamma,
                Gamma = Gamma, 
                parents = which(beta_star != 0),
                beta_star = beta_star))
  } else {
    stop("Not implemented.")
  }
}

#' DGP with 5 exposure and 5 instruments but 2 instruments 
#' 
#' @details 2 of the 5 exposure have a causal effect on the outcome.
#' 
dgp3 <- function(n, beta1, beta2, d1 = 1, d2 = 2, seed_num = NULL) {
  
  if (!is.null(seed_num)) set.seed(seed_num)
  
  # parameters
  d <- 5
  m <- 5
  B <- matrix(0, d, d)
  B[1,3] <- 1
  B[1,4] <- 1
  B[2,4] <- 1
  B[2,5] <- 1
  A <- diag(m)
  beta_star <- c(beta1, beta2, 0, 0, 0)
  Id <- diag(d)
  
  # variables
  H <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_x <- matrix(rnorm(n * d, 0, 1), n, d)
  eps_y <- rnorm(n, 0, 1)
  Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
  X <- (Z %*% t(A) + H + eps_x) %*% t(solve(Id-B))
  Y <- X %*% beta_star + Z[,c(1,2)] %*% c(d1, d2) + H %*% rep(1, d) + eps_y
  
  return(list(Z = Z,
              X = X,
              Y = Y,
              A = A,
              B = B,
              parents = which(beta_star != 0),
              beta_star = beta_star))
}
