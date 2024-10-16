test_that("k-calss suffstat and sumstat matches individual", {
  set.seed(123)
  n <- 1e3
  d <- 5
  m <- 3
  X <- matrix(rnorm(n * d), nrow = n)
  Z <- matrix(rnorm(n * m), nrow = n)
  Y <- X %*% rnorm(d) + rnorm(n)
  XX <- crossprod(X) / n
  YY <- crossprod(Y) / n
  ZZ <- crossprod(Z) / n
  ZX <- t(Z) %*% X / n
  ZY <- t(Z) %*% Y / n
  XY <- t(X) %*% Y / n
  k <- 0.5
  expect_equal(k_class_suffstat(k, XX, ZZ, XY, ZX, ZY),
               k_class(k, X, Y, Z %*% solve(crossprod(Z), t(Z)), n))
  expect_equal(liml_k_suffstat(XX, YY, ZZ, XY, ZX, ZY),
               liml_k(X, Y, Z %*% solve(crossprod(Z), t(Z)), n))
  expect_equal(liml_suffstat(XX, YY, ZZ, XY, ZX, ZY),
               liml(X, Y, Z %*% solve(crossprod(Z), t(Z)), n),
               tolerance = 1e-2)
})
