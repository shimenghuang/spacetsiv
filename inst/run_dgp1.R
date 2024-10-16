library(dplyr)
library(future.apply)
plan(multisession)

experiment_setting <- list(
  beta1 = 1, 
  beta2 = 2, 
  reps = 1:100,
  num_estimator = 3,
  n_grid = c(50, 500, 5000, 50000)
)
list2env(experiment_setting, .GlobalEnv)

ncum <- c(0, cumsum(n_grid))
n_index <- lapply(1:length(n_grid), function(k) (ncum[k]+1):ncum[k+1])
seeds <- future_lapply(seq_along(reps), FUN = function(x) .Random.seed,
                       future.chunk.size = Inf, future.seed = 42L)
df_lst <- future_lapply(reps, \(ii) {
  devtools::load_all()
  dat1 <- dgp1(sum(n_grid), beta1, beta2)
  dat2 <- dgp1(sum(n_grid), beta1, beta2)
  
  check_ass_A1 <- check_assumptions_A(dat1$A, dat1$B, dat1$parents, dat1$beta_star)
  valid_ass_A1 <- check_ass_A1$assA1 & check_ass_A1$assA3
  check_ass_A2 <- check_assumptions_A(dat2$A, dat2$B, dat2$parents, dat2$beta_star)
  valid_ass_A2 <- check_ass_A2$assA1 & check_ass_A2$assA3
  
  df_est <- data.frame(matrix(0, nrow = 0, ncol = length(dat1$beta_star) + 2))
  colnames(df_est) <- c(paste0("beta", 1:length(dat1$beta_star), "_hat"), "method", "n")
  df_res <- data.frame(rmse = numeric(0),
                       bias = numeric(0),
                       method = character(0),
                       valid_ass1 = numeric(0),
                       valid_ass2 = numeric(0),
                       n = numeric(0),
                       estimated_sparsity = numeric(0))
  
  for (kk in 1:length(n_grid)) {
    
    rmse <- vector("numeric", num_estimator)
    bias <- vector("numeric", num_estimator)
    estimated_sparsity <- vector("numeric", num_estimator)
    
    dat_sumstat <- comp_sumstat(dat1$Y[n_index[[kk]],,drop=FALSE], 
                                dat1$Z[n_index[[kk]],,drop=FALSE], 
                                dat2$X[n_index[[kk]],,drop=FALSE], 
                                dat2$Z[n_index[[kk]],,drop=FALSE], 
                                type = "depen_mv")
    
    # Q-subset with summary statistics
    beta_hat <- tryCatch({
      sparseMVMR_estimator(gamma_res = dat_sumstat$coef_y, 
                           sig_res = dat_sumstat$cov_coef_y,
                           gamma_exp = dat_sumstat$coef_x, 
                           sig_exp = dat_sumstat$cov_coef_x, 
                           indep_inst = FALSE, 
                           stop_if_accept = TRUE,
                           stop_rule = "min_set")$beta_hat
    }, error = \(e) rep(NA, length(dat1$beta_star)))
    tmp_row <- data.frame(matrix(beta_hat, nrow = 1), "Q-subset-sumstat", n_grid[kk])
    names(tmp_row) <- names(df_est)
    df_est <- rbind(df_est, tmp_row)
    
    rmse[1] <- sqrt(sum((dat1$beta_star - beta_hat)^2))
    bias[1] <- sum(abs(dat1$beta_star - beta_hat))
    estimated_sparsity[1] <- sum(beta_hat != 0)
    
    # TSIV-Lasso with summary statistics
    beta_hat <- grid_search_iv_lasso_sumstat(beta_init = rep(0.5, length(dat1$beta_star)), 
                                             dat_sumstat$coef_y, dat_sumstat$cov_coef_y,
                                             dat_sumstat$coef_x, dat_sumstat$cov_coef_x,
                                             cov_zinx = dat_sumstat$cov_zinx, 
                                             cov_ziny = dat_sumstat$cov_ziny, 
                                             pen_grid = pracma::logspace(-3, 1, 30), 
                                             alpha = 0.05,
                                             stop_if_accept = TRUE,
                                             stop_rule = "max_pen",
                                             max_iter = 500, 
                                             abs_tol = 1e-5,
                                             refit_set = TRUE,
                                             verbose = FALSE)$beta_hat
    tmp_row <- data.frame(matrix(beta_hat, nrow = 1), "IV-Lasso-sumstat-refit", n_grid[kk])
    names(tmp_row) <- names(df_est)
    df_est <- rbind(df_est, tmp_row)
    
    rmse[2] <- sqrt(sum((dat1$beta_star - beta_hat)^2))
    bias[2] <- sum(abs(dat1$beta_star - beta_hat))
    estimated_sparsity[2] <- sum(beta_hat != 0)
    
    # IV-sumstat
    beta_hat <- iv2_estimator_sumstat(gamma_res = dat_sumstat$coef_y, 
                                      sig_res = dat_sumstat$cov_coef_y, 
                                      gamma_exp = dat_sumstat$coef_x, 
                                      sig_exp = dat_sumstat$cov_coef_x,
                                      cov_zinx = dat_sumstat$cov_zinx, 
                                      cov_ziny = dat_sumstat$cov_ziny, 
                                      use_ginv = TRUE)
    tmp_row <- data.frame(matrix(beta_hat, nrow = 1), "IV-sumstat", n_grid[kk])
    names(tmp_row) <- names(df_est)
    df_est <- rbind(df_est, tmp_row)
    
    rmse[3] <- sqrt(sum((dat1$beta_star - beta_hat)^2))
    bias[3] <- sum(abs(dat1$beta_star - beta_hat))
    estimated_sparsity[3] <- sum(beta_hat != 0)
    
    df_res <- rbind(
      df_res,
      data.frame(rmse = rmse,
                 bias = bias,
                 method = c("Q-subset-sumstat", 
                            "IV-Lasso-sumstat-refit", 
                            "IV-sumstat"),
                 valid_ass1 = rep(valid_ass_A1, num_estimator),
                 valid_ass2 = rep(valid_ass_A2, num_estimator),
                 n = rep(toString(n_grid[kk]), num_estimator),
                 estimated_sparsity = estimated_sparsity))
  }
  list(df_est = df_est, 
       df_res = df_res)
}, future.seed = seeds)

experiment_result <- append(experiment_setting, list(df_lst = df_lst))

# Note: assuming running in the script folder!
saveRDS(experiment_result, 
        file = paste0("result/dgp1_2sample_sumstat.rds"))
