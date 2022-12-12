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

#---- multiple regression via ortogonalization -------------

orthogonalize = function(X){
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)
  G <- matrix(0, p, p)
  Z <- matrix(0, n, p)
  for (j in 1:p){
    v = X[,j]
    if (j > 1) {
       for (i in 1:(j-1)){
         G[i,j] = crossprod(Z[,i], X[,j]) / crossprod(Z[,i])
         v = v - G[i,j] * Z[,i]
      }
    }
    Z[,j] <- v 
  }
  return(list(G=G,Z=Z))
}

Z = orthogonalize(X)$Z
G = orthogonalize(X)$G

head(Z %*% G)
head(X)

Zp = Z[,p]
crossprod(Zp,y)/crossprod(Zp)
coef(lm(y ~ X -1))[4]


#---- Gram-Schmidt algorithm -------------

gramschmidt <- function(x) {
  x <- as.matrix(x)
  # Get the number of rows and columns of the matrix
  p <- ncol(x)
  n <- nrow(x)
  
  # Initialize the Q and R matrices
  q <- matrix(0, n, p)
  r <- matrix(0, p, p)
  
  for (j in 1:p) {
    v = x[,j] # Step 1 of the Gram-Schmidt process v1 = a1
    # Skip the first column
    if (j > 1) {
      for (i in 1:(j-1)) {
        r[i,j] <- t(q[,i]) %*% x[,j] # Find the inner product (noted to be q^T a earlier)
        # Subtract the projection from v which causes v to become perpendicular to all columns of Q
        v <- v - r[i,j] * q[,i] 
      }      
    }
    # Find the L2 norm of the jth diagonal of R
    r[j,j] <- sqrt(sum(v^2))
    # The orthogonalized result is found and stored in the ith column of Q.
    q[,j] <- v / r[j,j]
  }
  
  # Collect the Q and R matrices into a list and return
  qrcomp <- list('Q'=q, 'R'=r)
  return(qrcomp)
}


res <- gramschmidt(X)
res$R

qr_obj <- qr(X)
R <- qr.R(qr_obj)

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

#--- colonne (computazionalmente) linearmente dipendenti

qr(W)$rank
y <- rnorm(nrow(W))
fit <- lm(y ~ W)
tail(coef(fit))
