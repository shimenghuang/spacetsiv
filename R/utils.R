#' Take square-root of a matrix
#' 
#' @param X A positive semi-definite matrix.
#' 
#' @return "Square-root" of the input matrix.
#' 
sqrt_mat <- function(X, tol = 1e-5) {
  vec <- eigen(X)$vectors
  val <- eigen(X)$values
  val[which(abs(val) < tol)] <- 0
  return((vec %*% diag(sqrt(val)) %*% t(vec)))
}

jaccard_dist <- function(a, b) {
  return(length(dplyr::intersect(a,b)) / length(dplyr::union(a, b)))
}

binom_test <- function(x, n, alpha = 0.05) {
  tst <- binom.test(x, n, conf.level = 1-alpha)
  tibble::tibble(est = tst$estimate,
                 lwr = tst$conf.int[1],
                 upr = tst$conf.int[2])
}
