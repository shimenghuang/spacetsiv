library(dplyr)
library(ggplot2)

# ---- load data ----

load("script/glp1r_data//expression.RData") # from the paper by Ashish et al 2023

# beta 1, 8, 6, 9, 3 are top 5 in Table 2 in Patel et al 2023
exp_names <- c("brain caudate", "atrial appendage", "hypothalamus", "left ventricle", "lung",
               "nerve", "pancreas", "stomach", "testis", "thyroid")
bx <- cbind(brain$beta, # beta1, rank 1, #
            heart$beta, # beta2
            hypothalamus$beta, # beta3, rank 5, #
            left_ventricle$beta, # beta4, #
            lung$beta, # beta5, 
            nerve$beta, # beta6, rank 3
            pancreas$beta, # beta7, 
            stomach$beta, # beta8, rank 2, #
            testis$beta, # beta9, rank 4
            thyroid$beta) # beta10, 
sx <- cbind(brain$se,
            heart$se,
            hypothalamus$se,
            left_ventricle$se,
            lung$se,
            nerve$se,
            pancreas$se,
            stomach$se,
            testis$se,
            thyroid$se)
by <- cad$Effect
sy <- cad$StdErr
nx <- rep(838, ncol(bx))
ny <- mean(cad$n)
cor.x <- diag(ncol(bx))
cor.x[1,3] <- 0.7
cor.x[3,1] <- 0.7
cor.x[2,4] <- 0.7
cor.x[4,2] <- 0.7

# ---- select instruments ----

bx_new <- bx
sx_new <- sx
by_new <- by
sy_new <- sy

ld1 <- ld
ld2 <- (ld1 + t(ld1))/2
eigen(ld2)$values

strength <- sapply(1:10, \(ii) {
  fsf_test_sumstat_simple(bx_new[,ii], sx_new[,ii])$stat
})
ld2_names <- colnames(ld2)

top_num <- 2
colnames(strength) <- exp_names
res <- lapply(1:ncol(strength), \(ii) {
  idx <- order(strength[,ii], decreasing = TRUE)[1:top_num]
  data.frame(idx = idx, 
             snp = ld2_names[idx],
             stat = strength[idx,ii])
})
res_df <- do.call(rbind, res)
rownames(res_df) <- NULL
res_df$exposure <- rep(exp_names, each = top_num)
snps_keep <- unique(res_df$idx)
length(snps_keep)

eigen(ld2[snps_keep, snps_keep])$values
heatmap(ld2[snps_keep, snps_keep])

# ---- estimation ----

ld3 <- ld2[snps_keep, snps_keep]
eigen(ld3)$values
eigen(ld3)$values[1]/eigen(ld3)$values[length(eigen(ld3)$values)]

alpha <- 0.1

# adjust estimates to account for correlations
dat_sumstat <- adjust_estimates(coef_x = bx_new[snps_keep,], 
                                se_coef_x = sx_new[snps_keep,], 
                                coef_y = by_new[snps_keep], 
                                se_coef_y = sy_new[snps_keep],
                                cor_zinx = ld3, 
                                cor_ziny = ld3, 
                                cor_x = cor.x, 
                                nobs_x = nx[1], 
                                nobs_y = ny)

# initial values for beta_hat
beta_iv2 <- iv2_estimator_sumstat(dat_sumstat$coef_y, dat_sumstat$cov_coef_y,
                                  dat_sumstat$coef_x, dat_sumstat$cov_coef_x,
                                  cov_zinx = ld3, cov_ziny = ld3, 
                                  use_ginv = TRUE)

# ---- L0 ----

res_l0_all <- sparseMVMR_estimator(gamma_res = dat_sumstat$coef_y, 
                                   sig_res = dat_sumstat$cov_coef_y, 
                                   gamma_exp = dat_sumstat$coef_x, 
                                   sig_exp = dat_sumstat$cov_coef_x, 
                                   beta_init = beta_iv2, 
                                   indep_inst = FALSE, 
                                   alpha = alpha, 
                                   stop_if_accept = FALSE,
                                   stop_rule = "max_pval", 
                                   max_size = NULL, 
                                   verbose = TRUE, 
                                   optim_method = "L-BFGS-B",
                                   control = list(factr = 1e5), 
                                   lower = -1, 
                                   upper = 1)
saveRDS(res_l0_all, file = paste0("script/result/l0_all_ldclump_strongest", top_num, "_alpha", alpha, ".rds"))

beta_hat_maxpval <- res_l0_all$beta_hat_best[,which.max(res_l0_all$pval_best)]
idx_pass <- which(res_l0_all$pval_best > alpha)
res_l0_pass <- res_l0_all$beta_hat_best[,idx_pass]
beta_hat_minset <- res_l0_pass[,which.max(colSums(res_l0_pass == 0))]
beta_hat_maxset <- res_l0_pass[,which.min(colSums(res_l0_pass == 0))]
beta_hat_l0 <- cbind(beta_hat_maxpval, 
                     beta_hat_minset, 
                     beta_hat_maxset)
colnames(beta_hat_l0) <- c("maxpval", "minset", "maxset")
pval_l0 <- c("maxpval" = max(res_l0_all$pval_best),
             "minset" = res_l0_all$pval_best[idx_pass][which.max(colSums(res_l0_pass == 0))],
             "maxset" = res_l0_all$pval_best[idx_pass][which.min(colSums(res_l0_pass == 0))])
saveRDS(list(beta_hat = beta_hat_l0,
             pval = pval_l0),
        file = paste0("script/result/l0_df_ldclump_strongest", top_num, "_alpha", alpha, ".rds"))

ci_methods <- c("sandwich", "q_project", "gmm_approx")
for (ci_method in ci_methods) {
  summary_l0 <- lapply(c("minset"), \(xx) {
    beta_liml <- beta_hat_l0[,xx]
    (subvec_idx <- which(beta_liml != 0))
    (other_idx <- setdiff(1:length(beta_liml), subvec_idx))
    if (ci_method == "sandwich") {
      tmp <- sandwich_ci(dat_sumstat$coef_y,
                         dat_sumstat$cov_coef_y,
                         dat_sumstat$coef_x[,subvec_idx,drop=FALSE],
                         get_submat(subvec_idx, dat_sumstat$cov_coef_x, length(dat_sumstat$coef_y)),
                         beta_hat = beta_liml[subvec_idx])
    } else if (ci_method == "q_project") {
      tmp <- q_test_gen_subvec_ci(dat_sumstat$coef_y,
                                  dat_sumstat$cov_coef_y,
                                  dat_sumstat$coef_x,
                                  dat_sumstat$cov_coef_x,
                                  beta_liml = beta_liml,
                                  subvec_idx = subvec_idx,
                                  search_lwr = 1, search_upr = 1,
                                  uniroot_maxiter = 20,
                                  extend = TRUE,
                                  chisq_df = length(dat_sumstat$coef_y))
    } else if (ci_method == "gmm_approx") {
      tmp <- q_test_gen_ci(dat_sumstat$coef_y,
                           dat_sumstat$cov_coef_y,
                           dat_sumstat$coef_x[,subvec_idx,drop=FALSE],
                           get_submat(subvec_idx, dat_sumstat$cov_coef_x, length(dat_sumstat$coef_y)),
                           beta_hat = beta_liml[subvec_idx])
    }
    
    new_rows <- data.frame(matrix(NA, nrow = length(other_idx), ncol = 3))
    colnames(new_rows) <- colnames(tmp)
    tmp <- rbind(tmp, new_rows)
    if (length(other_idx) > 0) {
      tmp$est[(length(subvec_idx) + 1):nrow(tmp)] <- beta_liml[other_idx]
    }
    tmp$exposure <- c(exp_names[subvec_idx], exp_names[other_idx])
    tmp$stoprule <- xx
    tmp
  })
  summary_l0 <- do.call(rbind, summary_l0)
  summary_l0
  
  saveRDS(list(summary_l0 = summary_l0,
               pval_l0 = pval_l0), 
          file = paste0("script/result/l0_df_ldclump_strongest", top_num, "_alpha_", alpha, "_", ci_method, ".rds"))
}

# plot only for subset method
summary_l0 %>%
  ggplot2::ggplot(aes(x = exposure, y = est, color = stoprule, group = stoprule)) +
  ggplot2::geom_point(position = position_dodge(0.5)) + 
  ggplot2::geom_errorbar(aes(ymin = lwr, ymax = upr), 
                         width = 0.3, 
                         position = position_dodge(0.5)) + 
  ggplot2::theme_bw() +
  ggplot2::coord_flip() +
  ggplot2::theme( # axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
                 legend.position = "top", 
                 legend.box.margin = margin(0, 0, -1, 0), 
                 legend.spacing.y = unit(-1, "pt"),
                 legend.margin = margin(-4, 0, -4, 0),
                 text = element_text(size = 11)) +
  ggplot2::guides(color = guide_legend(nrow = 1, byrow = TRUE))

# size 1 sets that passed
idx_pass <- which(res_l0_all$pval_all[[1]] > alpha)
beta_hat_pass <- res_l0_all$beta_hat_all[[1]][,idx_pass,drop=FALSE]
pval_pass <- res_l0_all$pval_all[[1]][idx_pass]
colnames(beta_hat_pass) <- exp_names[idx_pass]
saveRDS(list(beta_hat = beta_hat_pass,
             pval = pval_pass),
        file = paste0("script/result/l0_df_ldclump_strongest", top_num, "_multipass_alpha", alpha, ".rds"))

ci_method <- "q_project"
summary_l0 <- lapply(1:ncol(beta_hat_pass), \(xx) {
  beta_liml <- beta_hat_pass[,xx]
  (subvec_idx <- which(beta_liml != 0))
  (other_idx <- setdiff(1:length(beta_liml), subvec_idx))
  if (ci_method == "sandwich") {
    tmp <- sandwich_ci(dat_sumstat$coef_y,
                       dat_sumstat$cov_coef_y,
                       dat_sumstat$coef_x[,subvec_idx,drop=FALSE],
                       get_submat(subvec_idx, dat_sumstat$cov_coef_x, length(dat_sumstat$coef_y)),
                       beta_hat = beta_liml[subvec_idx], 
                       alpha = alpha)
  } else if (ci_method == "q_project") {
    tmp <- q_test_gen_subvec_ci(dat_sumstat$coef_y,
                                dat_sumstat$cov_coef_y,
                                dat_sumstat$coef_x,
                                dat_sumstat$cov_coef_x,
                                beta_liml = beta_liml,
                                subvec_idx = subvec_idx,
                                alpha = alpha,
                                search_lwr = 1, search_upr = 1,
                                uniroot_maxiter = 20,
                                extend = TRUE,
                                chisq_df = length(dat_sumstat$coef_y))
  } else if (ci_method == "gmm_approx") {
    tmp <- q_test_gen_ci(dat_sumstat$coef_y,
                         dat_sumstat$cov_coef_y,
                         dat_sumstat$coef_x[,subvec_idx,drop=FALSE],
                         get_submat(subvec_idx, dat_sumstat$cov_coef_x, length(dat_sumstat$coef_y)),
                         beta_hat = beta_liml[subvec_idx], 
                         alpha = alpha)
  }
  
  new_rows <- data.frame(matrix(NA, nrow = length(other_idx), ncol = 3))
  colnames(new_rows) <- colnames(tmp)
  tmp <- rbind(tmp, new_rows)
  if (length(other_idx) > 0) {
    tmp$est[(length(subvec_idx) + 1):nrow(tmp)] <- beta_liml[other_idx]
  }
  tmp$exposure <- c(exp_names[subvec_idx], exp_names[other_idx])
  tmp$setnum <- as.character(xx)
  tmp
})
summary_l0 <- do.call(rbind, summary_l0)
summary_l0

saveRDS(list(summary_l0 = summary_l0,
             pval_l0 = pval_l0), 
        file = paste0("script/result/l0_df_ldclump_strongest", top_num, "_multipass_alpha", alpha, ".rds"))

# ---- L1 ----
res_l1_all <- grid_search_iv_lasso_sumstat(beta_init = beta_iv2, 
                                           gamma_res = dat_sumstat$coef_y, 
                                           sig_res = dat_sumstat$cov_coef_y,
                                           gamma_exp = dat_sumstat$coef_x, 
                                           sig_exp = dat_sumstat$cov_coef_x,
                                           cov_zinx = ld3,
                                           cov_ziny = ld3,
                                           pen_grid = pracma::logspace(-5, 1, 30),
                                           stop_if_accept = FALSE,
                                           stop_rule = "max_pval",
                                           alpha = alpha,
                                           chisq_df = NULL, 
                                           max_iter = 500, 
                                           abs_tol = 1e-5,
                                           refit_set = TRUE,
                                           lower = -1, 
                                           upper = 1, 
                                           verbose = FALSE)
res_l1_all$beta_hat_all_orig
res_l1_all$beta_hat_all_refit
saveRDS(res_l1_all, file = paste0("script/result/l1_all_ldclump_strongest", top_num, "_alpha", alpha, ".rds"))

beta_hat_maxpval <- res_l1_all$beta_hat_all_orig[,which.max(res_l1_all$pval_all_orig)]
idx_pass <- which(res_l1_all$pval_all_orig > alpha)
if (length(idx_pass) == 0) {
  idx_pass <- which.max(res_l1_all$pval_all_orig)
}
res_l1_pass <- res_l1_all$beta_hat_all_orig[,idx_pass,drop=FALSE]
beta_hat_maxpen <- res_l1_pass[,length(idx_pass)]
beta_hat_minpen <- res_l1_pass[,1]
beta_hat_l1 <- cbind(beta_hat_maxpval, beta_hat_maxpen, beta_hat_minpen)
colnames(beta_hat_l1) <- c("maxpval", "maxpen", "minpen")
pval_l1 <- c("maxpval" = max(res_l1_all$pval_all_orig),
             "maxpen" = res_l1_all$pval_all_orig[idx_pass][length(idx_pass)],
             "minpen" = res_l1_all$pval_all_orig[idx_pass][1])
saveRDS(list(beta_hat = beta_hat_l1,
             pval = pval_l1),
        file = paste0("script/result/l1_orig_df_ldclump_strongest", top_num, "_alpha", alpha, ".rds"))

beta_hat_maxpval <- res_l1_all$beta_hat_all_refit[,which.max(res_l1_all$pval_all_refit)]
idx_pass <- which(res_l1_all$pval_all_refit > alpha)
if (length(idx_pass) == 0) {
  idx_pass <- which.max(res_l1_all$pval_all_refit)
}
res_l1_pass <- res_l1_all$beta_hat_all_refit[,idx_pass,drop=FALSE]
beta_hat_maxpen <- res_l1_pass[,length(idx_pass)]
beta_hat_minpen <- res_l1_pass[,1]
beta_hat_l1 <- cbind(beta_hat_maxpval, beta_hat_maxpen, beta_hat_minpen)
colnames(beta_hat_l1) <- c("maxpval", "maxpen", "minpen")
pval_l1 <- c("maxpval" = max(res_l1_all$pval_all_refit),
             "maxpen" = res_l1_all$pval_all_refit[idx_pass][length(idx_pass)],
             "minpen" = res_l1_all$pval_all_refit[idx_pass][1])
saveRDS(list(beta_hat = beta_hat_l1,
             pval = pval_l1),
        file = paste0("script/result/l1_refit_df_ldclump_strongest", top_num, "_alpha", alpha, ".rds"))

ci_methods <- c("sandwich", "q_project", "gmm_approx")
for (ci_method in ci_methods) {
  summary_l1 <- lapply(c("maxpen"), \(xx) {
    beta_liml <- beta_hat_l1[,xx]
    (subvec_idx <- which(beta_liml != 0))
    (other_idx <- setdiff(1:length(beta_liml), subvec_idx))
    if (length(subvec_idx) > 0 & ci_method == "sandwich") {
      tmp <- sandwich_ci(dat_sumstat$coef_y,
                         dat_sumstat$cov_coef_y,
                         dat_sumstat$coef_x[,subvec_idx,drop=FALSE],
                         get_submat(subvec_idx, dat_sumstat$cov_coef_x, length(dat_sumstat$coef_y)),
                         beta_hat = beta_liml[subvec_idx], 
                         alpha = alpha)
    } else if (length(subvec_idx) > 0 & ci_method == "q_project") {
      tmp <- q_test_gen_subvec_ci(dat_sumstat$coef_y,
                                  dat_sumstat$cov_coef_y,
                                  dat_sumstat$coef_x,
                                  dat_sumstat$cov_coef_x,
                                  beta_liml = beta_liml,
                                  alpha = alpha, 
                                  subvec_idx = subvec_idx,
                                  search_lwr = 1, search_upr = 1,
                                  uniroot_maxiter = 20,
                                  extend = TRUE,
                                  chisq_df = length(dat_sumstat$coef_y))
    } else if (length(subvec_idx) > 0 & ci_method == "gmm_approx") {
      tmp <- q_test_gen_ci(dat_sumstat$coef_y,
                           dat_sumstat$cov_coef_y,
                           dat_sumstat$coef_x[,subvec_idx,drop=FALSE],
                           get_submat(subvec_idx, dat_sumstat$cov_coef_x, length(dat_sumstat$coef_y)),
                           beta_hat = beta_liml[subvec_idx],
                           alpha = alpha)
    } else {
      tmp <- data.frame(matrix(NA, ncol = 3, nrow = 0))
      colnames(tmp) <- c("lwr", "upr", "est")
    }
    
    new_rows <- data.frame(matrix(NA, nrow = length(other_idx), ncol = 3))
    colnames(new_rows) <- colnames(tmp)
    tmp <- rbind(tmp, new_rows)
    if (length(other_idx) > 0) {
      tmp$est[(length(subvec_idx) + 1):nrow(tmp)] <- beta_liml[other_idx]
    }
    tmp$exposure <- c(exp_names[subvec_idx], exp_names[other_idx])
    tmp$stoprule <- xx
    tmp
  })
  summary_l1 <- do.call(rbind, summary_l1)
  summary_l1
  
  saveRDS(list(summary_l1 = summary_l1,
               pval_l1 = pval_l1), 
          file = paste0("script/result/l1_df_ldclump_strongest", top_num, "_alpha", alpha, "_", ci_method, ".rds"))
}

summary_l1 %>%
  ggplot2::ggplot(aes(x = exposure, y = est, color = stoprule, group = stoprule)) +
  ggplot2::geom_point(position = position_dodge(0.5)) + 
  ggplot2::geom_errorbar(aes(ymin = lwr, ymax = upr), 
                         width = 0.3, 
                         position = position_dodge(0.5)) + 
  ggplot2::theme_bw() +
  ggplot2::theme(# axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
                 legend.position = "top", 
                 legend.box.margin = margin(0, 0, -1, 0), 
                 legend.spacing.y = unit(-1, "pt"),
                 legend.margin = margin(-4, 0, -4, 0),
                 text = element_text(size = 11)) +
  ggplot2::guides(color = guide_legend(nrow = 2, byrow = TRUE)) + 
  ggplot2::coord_flip()

# ---- TSIV ----

beta_hat <- iv2_estimator_sumstat(gamma_res = dat_sumstat$coef_y, 
                                  sig_res = dat_sumstat$cov_coef_y, 
                                  gamma_exp = dat_sumstat$coef_x, 
                                  sig_exp = dat_sumstat$cov_coef_x,
                                  cov_zinx = dat_sumstat$cov_zinx, 
                                  cov_ziny = dat_sumstat$cov_ziny, 
                                  use_ginv = TRUE)

pval_tsiv <- q_test_gen(coef_y = dat_sumstat$coef_y, 
                        cov_coef_y = dat_sumstat$cov_coef_y, 
                        coef_x = dat_sumstat$coef_x, 
                        cov_coef_x = dat_sumstat$cov_coef_x, 
                        beta = beta_hat)

# Note: this takes a while
ci_methods <- c("q_project") # c("sandwich", "q_project", "gmm_approx")
for (ci_method in ci_methods) {
  summary_tsiv <- lapply(1, \(xx) {
    (subvec_idx <- which(beta_hat != 0))
    (other_idx <- setdiff(1:length(beta_hat), subvec_idx))
    if (length(subvec_idx) > 0 & ci_method == "sandwich") {
      tmp <- sandwich_ci(dat_sumstat$coef_y,
                         dat_sumstat$cov_coef_y,
                         dat_sumstat$coef_x[,subvec_idx,drop=FALSE],
                         get_submat(subvec_idx, dat_sumstat$cov_coef_x, length(dat_sumstat$coef_y)),
                         beta_hat = beta_hat[subvec_idx], 
                         alpha = alpha)
    } else if (length(subvec_idx) > 0 & ci_method == "q_project") {
      tmp <- q_test_gen_subvec_ci(dat_sumstat$coef_y,
                                  dat_sumstat$cov_coef_y,
                                  dat_sumstat$coef_x,
                                  dat_sumstat$cov_coef_x,
                                  beta_liml = beta_hat,
                                  alpha = alpha, 
                                  uniroot_maxiter = 20, 
                                  subvec_idx = subvec_idx,
                                  search_lwr = 1, search_upr = 1,
                                  extend = TRUE,
                                  chisq_df = length(dat_sumstat$coef_y))
    } else if (length(subvec_idx) > 0 & ci_method == "gmm_approx") {
      tmp <- q_test_gen_ci(dat_sumstat$coef_y,
                           dat_sumstat$cov_coef_y,
                           dat_sumstat$coef_x[,subvec_idx,drop=FALSE],
                           get_submat(subvec_idx, dat_sumstat$cov_coef_x, length(dat_sumstat$coef_y)),
                           beta_hat = beta_hat[subvec_idx],
                           alpha = alpha)
    } else {
      tmp <- data.frame(matrix(NA, ncol = 3, nrow = 0))
      colnames(tmp) <- c("lwr", "upr", "est")
    }
    
    new_rows <- data.frame(matrix(NA, nrow = length(other_idx), ncol = 3))
    colnames(new_rows) <- colnames(tmp)
    tmp <- rbind(tmp, new_rows)
    if (length(other_idx) > 0) {
      tmp$est[(length(subvec_idx) + 1):nrow(tmp)] <- beta_liml[other_idx]
    }
    tmp$exposure <- c(exp_names[subvec_idx], exp_names[other_idx])
    tmp$stoprule <- xx
    tmp
  })
  summary_tsiv <- do.call(rbind, summary_tsiv)
  summary_tsiv
  
  saveRDS(list(summary_tsiv = summary_tsiv,
               pval_tsiv = pval_tsiv), 
          file = paste0("script/result/tsiv_df_ldclump_strongest", top_num, "_alpha", alpha, "_", ci_method, ".rds"))
}

# ---- summarize results ----
alpha <- 0.1
top_num <- 2
exp_names <- c("brain caudate", "atrial appendage", "hypothalamus", "left ventricle", "lung",
               "nerve", "pancreas", "stomach", "testis", "thyroid")
ci_method <- "q_project"

l0_df_ldclump0.05 <- readRDS(paste0("script/result/l0_df_ldclump_strongest", top_num, "_multipass_alpha", alpha, ".rds"))
df_l0 <- l0_df_ldclump0.05$summary_l0
df_l0$penalty <- "l0"
df_l0 <- df_l0 %>%
  dplyr::rename(stoprule = setnum)

l1_df_ldclump0.05 <- readRDS(paste0("script/result/l1_df_ldclump_strongest", top_num, "_alpha", alpha, "_", ci_method, ".rds"))
df_l1 <- l1_df_ldclump0.05$summary_l1
df_l1$penalty <- "l1"
df_all <- rbind(df_l0, df_l1)
df_all$exposure <- factor(df_all$exposure, levels = exp_names)

tsiv_df_ldclump0.05 <- readRDS(paste0("script/result/tsiv_df_ldclump_strongest", top_num, "_alpha", alpha, "_", ci_method, ".rds"))
df_tsiv <- tsiv_df_ldclump0.05$summary_tsiv
df_tsiv$penalty <- "tsiv"
df_tsiv$stoprule <- ""

df_all <- rbind(df_l0, df_l1, df_tsiv)
df_all$exposure <- factor(df_all$exposure, levels = exp_names)

# select an order for the exposures
exp_order <- df_all %>%
  dplyr::group_by(exposure) %>%
  dplyr::summarise(tot = sum(abs(est)), zero = sum(est == 0)) %>%
  dplyr::arrange(zero, tot) %>%
  dplyr::pull(exposure) %>%
  as.character()

df_long <- df_all %>%
  dplyr::mutate(method = paste0(penalty, "_", stoprule),
                # method = factor(method, 
                #                 levels = rev(c("l0_maxpval", 
                #                                "l0_minpen",
                #                                "l0_maxpen",
                #                                "l1_maxpval", 
                #                                "l1_minpen",
                #                                "l1_maxpen"))),
                exposure = factor(exposure, levels = rev(exp_order)))

library(RColorBrewer)

df_long %>%
  dplyr::filter(!method %in% c("l1_minpen", "l1_maxpval"),
                est != 0) %>%
  dplyr::select(method, exposure, est, lwr, upr) %>%
  # View() %>%
  knitr::kable(format = "latex", booktabs = TRUE, 
               row.names = FALSE, digits = 4)

library(dichromat)
colorRampPalette(c("#5E4FA2", "#cec4f8"))(7)
# "#5E4FA2" "#7062B0" "#8376BE" "#9689CD" "#A89DDB" "#BBB0E9" "#CEC4F8"

df_long %>%
  dplyr::mutate(grid_positions = as.numeric(exposure) + 0.5,
                grid_positions = ifelse(grid_positions == max(as.numeric(exposure) + 0.5),
                                        NA, grid_positions)) %>%
  dplyr::filter(!method %in% c("l1_minpen", "l1_maxpval")) %>%
  # dplyr::filter(!method %in% c("tsiv_")) %>%
  ggplot2::ggplot(aes(x = exposure, y = est, color = method, shape = (est == 0))) +
  ggplot2::geom_point(alpha = 1, size = 1, 
                      position = position_dodge(0.75),
                      show.legend = c(color = TRUE, shape = FALSE)) +
  ggplot2::geom_errorbar(aes(ymin = lwr, ymax = upr),
                         position = position_dodge(0.75),
                         size = 0.4,
                         alpha = 1) + 
  ggplot2::scale_color_manual(labels = c("l0_1" = "spaceTSIV-L0",
                                         # "l0_2" = "spaceTSIV-L0 (set2)",
                                         # "l0_3" = "spaceTSIV-L0 (set3)",
                                         # "l0_4" = "spaceTSIV-L0 (set4)",
                                         # "l0_5" = "spaceTSIV-L0 (set5)",
                                         # "l0_6" = "spaceTSIV-L0 (set6)",
                                         # "l0_7" = "spaceTSIV-L0 (set7)",
                                         "l1_maxpen" = "spaceTSIV-L1",
                                         "tsiv_" = "TSIV"), # 
                              values = c("l0_1" = "#5E4FA2",
                                         # "l0_2" = "#7062B0",
                                         # "l0_3" = "#8376BE",
                                         # "l0_4" = "#9689CD",
                                         # "l0_5" = "#A89DDB",
                                         # "l0_6" = "#BBB0E9",
                                         # "l0_7" = "#cec4f8",
                                         "l1_maxpen" = "#74C476",
                                         "tsiv_" = "#FDAE61")) + 
  ggplot2::scale_shape_manual(name = "", 
                              labels = c("selected", "not selected"), 
                              values = c(16, 4)) + 
  ggplot2::geom_vline(aes(xintercept = grid_positions), alpha = 0.1) + 
  ggplot2::theme_bw() +
  ggplot2::theme(# axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
                 legend.position = "top",
                 legend.box.margin = margin(0, 0, -1, 0), 
                 legend.spacing.y = unit(-1, "pt"),
                 legend.margin = margin(-4, 0, -4, 0),
                 text = element_text(size = 20),
                 panel.grid.major.y = element_blank()) +
  ggplot2::guides(color = guide_legend(title = "", 
                                       position = "top", 
                                       nrow = 1, 
                                       byrow = TRUE, 
                                       reverse = FALSE)) + 
  ggplot2::coord_flip() +
  ggplot2::labs(y = "estimate")

ggplot2::ggsave(paste0("script/result/glp1r_strongest", top_num, "_alpha_", alpha, "_with_", ci_method,
                       ".pdf"),
                width = 7, height = 4)

