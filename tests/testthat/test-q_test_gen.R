
test_that("q_test_gen_grad is correctly implemented", {
  set.seed(123)
  n <- 1e3
  dat1 <- dgp1(n, beta1 = 1, beta2 = 1)
  dat2 <- dgp1(n, beta1 = 1, beta2 = 1)
  dat_sumstat <- comp_sumstat(Y = dat1$Y, Z_Y = dat1$Z, X = dat2$X, Z_X = dat2$Z, type = "depen_mv")

  my_grad_res <- as.vector(q_stat_gen_grad(coef_y = dat_sumstat$coef_y, 
                                           cov_coef_y = dat_sumstat$cov_coef_y, 
                                           coef_x = dat_sumstat$coef_x,
                                           cov_coef_x = dat_sumstat$cov_coef_x, 
                                           beta = dat1$beta_star))
  num_grad_res <- numDeriv::grad(\(x) {
    q_test_gen(coef_y = dat_sumstat$coef_y, 
               cov_coef_y = dat_sumstat$cov_coef_y, 
               coef_x = dat_sumstat$coef_x,
               cov_coef_x = dat_sumstat$cov_coef_x, 
               beta = x)$stat
  }, x = dat1$beta_star)
  
  expect_equal(my_grad_res, num_grad_res)
})
