rm(list=ls())
#==========================
# ADDITIVE MODELS
#==========================

#--- Great Barrier Reef data ----------------------------

library(sm)
dat <- trawl

#pdf("Figure_trawl.pdf")
plot(Score1 ~ Longitude, dat )
#dev.off()

ns_x <- function(x, df){
  
  x <- as.vector(x)
  n <- length(x)
  nIknots <- df - 1L
  knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
                     2L)[-c(1L, nIknots + 2L)]
  xi <- quantile(x, knots)
  K <- length(xi)

  d <- function(z, j)
  {
    out <- (x - xi[j])^3 * as.numeric(x > xi[j])
    out <- out - (x - xi[K])^3 * as.numeric(x > xi[K])
    out <- out / (xi[K] - xi[j])
    out
  }
  
  B <- matrix(0, ncol=K, nrow=n)
  B[, 1L] <- 1
  B[, 2L] <- x
  for (j in seq(1L, (K-2L)))
  {
    B[, j + 2L] <- d(x, j) - d(x, K - 1L)
  }
  B
  
}

x <- dat$Longitude
y <- dat$Score1
df = 10
B <- ns_x(x, df)
fit <- lm(y ~ 0+B)
yhat <- predict(fit)

#pdf("Figure_trawl_ns.pdf")
plot(x,y)
matlines(x,predict(fit,interval = "confidence"), lty=1, col=4)
#dev.off()

library(gam)

pdf("Figure_gam.pdf")
fit <- gam(Score1 ~   s(Latitude) + s(Longitude), dat, family="gaussian")
plot(fit, rugplot = TRUE, se = TRUE )
dev.off()

#--- Esercizio 1 ----------------------------

n = 100
set.seed(1)
x1 = rnorm(n)
x2 = rnorm(n)
x3 = rnorm(n)
y = 30 - 10*x1 + 20*x2 + 30*x3 + rnorm(n)
fit = lm(y ~ x1+x2+x3)
X = model.matrix(fit)

# center y and X
yc = y - mean(y)
Xc = apply(X[,-1],2,function(x) x-mean(x))
# initialize 
p = ncol(Xc)
hatbeta = rep(0,p)
delta = 0.1
gamma = rep(2*delta,p)
while ( !all( abs(gamma - hatbeta) <= delta) ){
  for (k in 1:p){
    gamma = hatbeta
    yk = yc - Xc[,-k, drop=F] %*% hatbeta[-k, drop=F]
    gammak = coef(lm(yk ~ 0+ Xc[,k]))
    hatbeta[k] <- gammak
  }}
hatbeta0<- mean(y) - sum(hatbeta*apply(X[,-1],2,mean))
data.frame( backfitting=c(hatbeta0,hatbeta), leastsquares=coef(fit))


#--- Esercizio 2 ----------------------------

am_backfit <-
  function(X, y, maxit=10L)
  {
    p <- ncol(X)
    id <- seq_len(nrow(X))
    alpha <- mean(y)
    f <- matrix(0, ncol = p, nrow = nrow(X))
    models <- vector("list", p + 1L)
    for (i in seq_len(maxit))
    {
      for (j in seq_len(p))
      {
        p_resid <- y - alpha - apply(f[, -j], 1L, sum)
        id <- order(X[,j])
        models[[j]] <- smooth.spline(X[id,j], p_resid[id])
        f[,j] <- predict(models[[j]], X[,j])$y
      }
      alpha <- mean(y - apply(f, 1L, sum))
    }
    models[[p + 1L]] <- alpha
    return(models)
  }

am_predict <-
  function(models, X_new)
  {
    p <- ncol(X_new)
    f <- matrix(0, ncol = p, nrow = nrow(X_new))
    for (j in seq_len(p))
    {
      f[,j] <- predict(models[[j]], X_new[,j])$y
    }
    y <- apply(f, 1L, sum) + models[[p + 1L]]
    list(y=y, f=f)
  }

set.seed(123)
n <- 500; p <- 4
X <- matrix(runif(n * p, min = -2, max = 2), ncol = p)
f1 <- cos(X[,1] * 4) + sin(X[,1] * 10) + X[,1]^(2)
f2 <- -1.5 * X[,2]^2 + (X[,2] > 1) * (X[,2]^3 - 1)
f3 <- 0
f4 <- sign(X[,4]) * 1.5
f1 <- f1 - mean(f1); f2 <- f2 - mean(f2)
f3 <- f3 - mean(f3); f4 <- f4 - mean(f4)
y <- 10 + f1 + f2 + f3 + f4 + rnorm(n, sd = 1.2)

models <- am_backfit(X, y, maxit = 1)
yhat <- am_predict(models,X)

plot(X[,1],yhat$f[,1], type="p")
points(X[,1], f1, col=2)

plot(X[,2],yhat$f[,2], type="p")
points(X[,2], f2, col=2)

plot(X[,3],yhat$f[,3], type="p")
points(X[,3], rep(f3,nrow(X)), col=2)

plot(X[,4],yhat$f[,4], type="p")
points(X[,4], f4, col=2)
