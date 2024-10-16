# ---- using individual data ----

#' K-class estimator
#' 
k_class <- function(k, X, Y, P_Z, nobs){
  
  Wk <- t(X) %*% ( (1-k) * diag(nobs) + k * P_Z ) %*% X
  est <- solve(Wk, t(X)) %*% ( (1-k)*diag(nobs) + k * P_Z ) %*% Y 
   
  return(est)
}

#' Calculate the k-value for LIML
#' 
#' @details See e.g., https://en.wikipedia.org/wiki/Simultaneous_equations_model#Limited_information_maximum_likelihood_(LIML)
#' 
liml_k <- function(X, Y, P_Z, nobs) {
  
  M_Z <- diag(nobs) - P_Z
  W <- t(cbind(Y, X)) %*% M_Z %*% cbind(Y, X)
  W_1 <- t(cbind(Y,X)) %*% cbind(Y,X)
  
  return(min(abs(eigen(W_1 %*% solve(W), only.values = TRUE)$values)))
}

#' LMIL estimator 
#' 
liml <- function(X, Y, P_Z, nobs) {
  kval <- liml_k(X = X, Y = Y, P_Z = P_Z, nobs = nobs)
  return(k_class(kval, X = X, Y = Y, P_Z = P_Z, nobs = nobs))
}

# ---- using sufficient statistics ----

#' K-class estimator using sufficient statistics
#' 
#' @details Adapted from the code of spaceIV.
#' 
#' @param k A scalar for the k-value.
#' @template xyz_suffstat_wo_yy
#' 
k_class_suffstat <- function(k, XX, ZZ, XY, ZX, ZY) {
  
  Wk <- (1-k) * XX + k * t(ZX) %*% solve(ZZ, ZX)
  est <- (1-k) * solve(Wk, XY) + k * solve(Wk, t(ZX)) %*% solve(ZZ, ZY)
  
  return(est)
}


#' Calculate the k-value for LIML using sufficient statistics
#' 
#' @template xyz_suffstat
#' 
liml_k_suffstat <- function(XX, YY, ZZ, XY, ZX, ZY) {
  
  W <- rbind(cbind(XY - t(ZX) %*% solve(ZZ, ZY), XX - t(ZX) %*% solve(ZZ, ZX)),
             cbind(YY - t(ZY) %*% solve(ZZ, ZY), t(XY) - t(ZY) %*% solve(ZZ, ZX)))
  W_1 <- rbind(cbind(XY, XX),
               cbind(YY, t(XY)))
  
  return(min(abs(eigen(W_1 %*% solve(W), only.values = TRUE)$values)))
}


#' LIML estimator using using sufficient statistics
#' 
#' @template xyz_suffstat
#' 
liml_suffstat <- function(XX, YY, ZZ, XY, ZX, ZY) {
  
  kval <- liml_k_suffstat(XX, YY, ZZ, XY, ZX, ZY)
  
  return(k_class_suffstat(kval, XX, ZZ, XY, ZX, ZY))
}


