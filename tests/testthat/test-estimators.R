test_that("compute summary statistics", {
  set.seed(123)
  n1 <- 1e3
  n2 <- 1e3
  n <- n1 + n2
  I1 <- matrix(rnorm(n, 0, 1))
  I2 <- matrix(rnorm(n, 0, 1))
  gamma1 <- c(1, 2)
  gamma2 <- c(2, 1)
  X1 <- cbind(I1, I2) %*% gamma1 + rnorm(n, 0, 1)
  X2 <- cbind(I1, I2) %*% gamma2 + rnorm(n, 0, 1)
  beta <- c(1.5, 0)
  Y <- cbind(X1, X2) %*% beta + rnorm(n, 0, 1)
  
  X1 <- X1[1:n1,] # sample 1
  X2 <- X2[1:n1,] # sample 1
  Y <- Y[(n1+1):n] # sample 2
  
  # use the function
  dat_ss <- comp_sumstat(Y = Y, Z_Y = cbind(I1, I2)[(n1+1):n,], 
                         X = cbind(X1, X2), 
                         Z_X = cbind(I1, I2)[1:n1,],
                         type = "indep_mv")
  
  # manual calculation
  mod1 <- lm(Y ~ I1[(n1+1):n] + I2[(n1+1):n]) # Gamma1, Gamma2
  mod2 <- systemfit::systemfit(list(X1 ~ I1[1:n1] + I2[1:n1],
                                    X2 ~ I1[1:n1] + I2[1:n1]),
                               method = "SUR")
  
  GG <- coef(mod1)[c(2,3)]
  sig_GG <- sqrt(diag(vcov(mod1))[c(2,3)])

  gg <- coef(mod2)[-c(1,4)]
  gg <- matrix(gg, nrow = 2)
  sig_gg <- vcov(mod2)[-c(1,4), -c(1,4)]
  sig_gg1 <- sig_gg[c(1,3), c(1,3)]
  sig_gg2 <- sig_gg[c(2,4), c(2,4)]
  
  # compare results
  expect_equal(dat_ss$coef_y, GG, ignore_attr = TRUE)
  expect_equal(dat_ss$coef_x, gg, ignore_attr = TRUE)
  expect_equal(dat_ss$cov_coef_x[[1]], sig_gg1, ignore_attr = TRUE)
  expect_equal(dat_ss$cov_coef_x[[2]], sig_gg2, ignore_attr = TRUE)
})

test_that("adjust_estimates should match SUR", {
  set.seed(123)
  n1 <- 1e4
  n2 <- 1e4
  n <- n1 + n2
  I1 <- matrix(rnorm(n, 0, 1))
  I2 <- 0.3 * I1 + matrix(rnorm(n, 0, 1))
  gamma1 <- c(1, 2)
  gamma2 <- c(2, 1)
  H <- matrix(rnorm(n, 0, 1))
  X1 <- cbind(I1, I2) %*% gamma1 + H + rnorm(n, 0, 1)
  X2 <- cbind(I1, I2) %*% gamma2 + H + rnorm(n, 0, 1)
  beta <- c(1.5, 0)
  Y <- cbind(X1, X2) %*% beta + H + rnorm(n, 0, 1)
  
  X1 <- X1[1:n1,] # sample 1
  X2 <- X2[1:n1,] # sample 1
  Y <- Y[(n1+1):n] # sample 2
  
  # calculation with SUR
  mod1 <- lm(Y ~ I1[(n1+1):n] + I2[(n1+1):n] - 1) # Gamma1, Gamma2
  mod2 <- systemfit::systemfit(list(X1 ~ I1[1:n1] + I2[1:n1] - 1,
                                    X2 ~ I1[1:n1] + I2[1:n1] - 1),
                               method = "SUR")
  
  # calculation with summary statistics
  mod_y1 <- lm(Y ~ I1[(n1+1):n] - 1)
  mod_y2 <- lm(Y ~ I2[(n1+1):n] - 1)
  mod_x11 <- lm(X1 ~ I1[1:n1] - 1)
  mod_x12 <- lm(X1 ~ I2[1:n1] - 1)
  mod_x21 <- lm(X2 ~ I1[1:n1] - 1)
  mod_x22 <- lm(X2 ~ I2[1:n1] - 1)
  
  res <- adjust_estimates(coef_x = matrix(c(coef(mod_x11)[1], coef(mod_x12)[1],
                                            coef(mod_x21)[1], coef(mod_x22)[1]),
                                          nrow = 2), 
                          se_coef_x = sqrt(matrix(c(vcov(mod_x11)[1,1], vcov(mod_x12)[1,1],
                                                    vcov(mod_x21)[1,1], vcov(mod_x22)[1,1]),
                                                  nrow = 2)), 
                          coef_y = c(coef(mod_y1)[1], coef(mod_y2)[1]),
                          se_coef_y = sqrt(c(vcov(mod_y1)[1,1], vcov(mod_y2)[1,1])), 
                          cor_zinx = cor(cbind(I1[1:n1], I2[1:n1])),
                          cor_ziny = cor(cbind(I1[(n1+1):n], I2[(n1+1):n])), 
                          cor_x = cor(cbind(X1, X2)),
                          nobs_x = n1, nobs_y = n2)
  
  # compare results
  expect_equal(res$coef_y, coef(mod1), tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(res$cov_coef_y, vcov(mod1), tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(res$coef_x, matrix(coef(mod2), nrow = 2), tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(res$cov_coef_x, vcov(mod2), tolerance = 1e-3, ignore_attr = TRUE)
})


