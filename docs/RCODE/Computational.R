rm(list=ls())
#==========================
# ASPETTI COMPUTAZIONALI
#==========================

#---- le equazioni normali -------------

# generiamo dei dati con n = 500 e p = 100
set.seed(123)
n <- 500; p <- 100
X <- matrix(rnorm(n*p), ncol=p)
y <- rnorm(n)

library(microbenchmark)
ms <- microbenchmark(solve(t(X) %*% X) %*% t(X) %*% y, 
                     solve(crossprod(X), crossprod(X, y)))
ms

#---- Problemi di multicollinearità -------------

W <- cbind(X, X[, 1] + rnorm(n, sd = 1e-10))
#solve(crossprod(W), crossprod(W, y))

kappa(crossprod(W))

#---- ols_chol -------------

# Calcolare la stima OSL con la decomposizione di Cholesky
#
# Argomenti:
# X: la matrice del disegno
# y: il vettore risposta
#
# Ritorna:
# Il vettore hatbeta di lunghezza ncol(X).
ols_chol <-
  function(X, y)
  {
    XtX <- crossprod(X)
    Xty <- crossprod(X, y)
    L <- t(chol(XtX))
    betahat <- backsolve(t(L), forwardsolve(L, Xty))
    betahat
  }


# generiamo dei dati con n = 100 e p=4
n <- 1e3; p <-4
X <- matrix(rnorm(n*p), ncol=p)
beta = 1:4
epsilon <- rnorm(n)
y <- X %*% beta + epsilon 

ols_chol(X,y)
coef(lm(y ~ X - 1))

#---- (Non-normalized) Classical Gram-Schmidt -------------

classicGS = function(X){
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)
  Z <- matrix(0, n, p)
  for (j in 1:p){
    Zj = X[,j]
    if (j > 1) {
      for (k in 1:(j-1)){
        coef = crossprod(Z[,k], X[,j]) / crossprod(Z[,k]) 
        Zj = Zj - coef * Z[,k]
      }
    }
    Z[,j] <- Zj
  }
  return(Z)
}

Z = classicGS(X)
Zp = Z[,p]
crossprod(Zp,y)/crossprod(Zp)
coef(lm(y ~ X -1))[4]

#---- QR factorization -------------

factorizationQR = function(X){
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)
  Q <- matrix(0, n, p)
  R <- matrix(0, p, p)
  for (j in 1:p){
    Zj = X[,j]
    if (j > 1) {
      for (k in 1:(j-1)){
        R[k,j] = crossprod(Q[,k], X[,j])
        Zj = Zj - R[k,j] * Q[,k]
      }
    }
    R[j,j] <- sqrt( crossprod(Zj) )
    Q[,j] <- Zj / R[j,j]
  }
  return(list(Q=Q, R=R))
}

res = factorizationQR(X)
res$R

res_qr <- qr(X)
qr.R(res_qr)

#---- ols_qr -------------

# Calcolare le stime OLS con la decomposizione QR 
#
# Argomenti:
# X: la matrice del disegno
# y: il vettore risposta
#
# Ritorna:
# Il vettore hatbeta di lunghezza ncol(X).
ols_qr <-
  function(X, y)
  {
    qr_obj <- qr(X)
    Q <- qr.Q(qr_obj)
    R <- qr.R(qr_obj)
    Qty <- crossprod(Q, y)
    betahat <- backsolve(R, Qty)
    betahat
  }

ols_qr(X,y)
coef(lm(y ~ X -1))


