rm(list=ls())
#==========================
# REGRESSION SPLINES
#==========================

#--- DATI SIMULATI ----------------------------

fx = function(x) 1 - x^3 - 2*exp(-100*x^2)
x = seq(-1,1,by=0.02)
set.seed(123)
n = length(x)
y = fx(x) + rnorm(n,0,0.1)

#--- REGRESSIONE POLINOMIALE ----------------------------

fit = lm(y ~ poly(x,degre=17, raw=T))
#pdf("Figure_poly17.pdf")
plot(x,y, main="Regressione polinomiale di grado 17")
lines(x,fitted(fit))
#dev.off()

B = model.matrix(fit)
S = B %*% solve(crossprod(B)) %*% t(B)

#pdf("Figure_hatmatrix17.pdf")
filled.contour(main="Regressione polinomiale di grado 17",apply(S, 2, rev), asp=1, xaxt="n", yaxt="n", xlab="columns", ylab="rows", levels=c(min(S),-0.005,0.005,max(S)), col=c("red", "white", "blue"))
#dev.off()

kappa(crossprod(B), exact=T)

#--- POLINOMI0 A TRATTI ----------------------------

M = 2
knots = seq(-1,1,length=14)[-c(1,14)]
x_cut = cut(x, c(-Inf,knots,Inf) )
#pdf("Figure_piece_poly.pdf")
plot(x,y)
abline(v=knots, col="gray", lty=3)
S = matrix(0,n,n)
for (i in 1:length(levels(x_cut)) ){
  subset = (x_cut==levels(x_cut)[i])
  lines(x[subset], fitted(lm(y[subset] ~ poly(x[subset],M))), lwd=2 )
  B = model.matrix(lm(y[subset] ~ poly(x[subset],M)))
  S[subset,subset] = B %*% solve(crossprod(B)) %*% t(B)
}
#dev.off()

#pdf("Figure_hatmatrix_piece_poly.pdf")
filled.contour(main="Polinomio a tratti di grado 2",apply(S, 2, rev), asp=1, xaxt="n", yaxt="n", xlab="columns", ylab="rows", levels=c(min(S),-0.005,0.005,max(S)), col=c("red", "white", "blue"))
#dev.off()


#--- SPLINE DI REGRESSIONE ----------------------------

K = length(knots)
n = length(x)
B = matrix(NA,ncol=M+K+1, nrow=n)
B[,1:(M+1)] = cbind(1,poly(x,M, raw=TRUE))
B[,(M+2):(M+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^M, 0))
fit = lm(y ~ 0 + B)

#pdf("Figure_reg_spline.pdf")
plot(x,y)
abline(v=knots, col="gray", lty=3)
lines(x,predict(fit), lwd=2)
#dev.off()

S = B %*% solve(crossprod(B)) %*% t(B)
#pdf("Figure_hatmatrix_reg_spline.pdf")
filled.contour(main="Spline di regressione di grado 2",apply(S, 2, rev), asp=1, xaxt="n", yaxt="n", xlab="columns", ylab="rows", levels=c(min(S),-0.005,0.005,max(S)), col=c("red", "white", "blue"))
#dev.off()


#--- RADIOCARBON DATA ----------------------------

library(sm) 
dat = radioc[(radioc$Cal.age > 2000 & radioc$Cal.age < 3000), ]
x = dat$Cal.age
y = dat$Rc.age
plot(y ~ x)

K = 6
knots = seq(min(x),max(x),length=K+2)[-c(1,K+2)]
M = 3
n = length(x)
B = matrix(NA,ncol=M+K+1, nrow=n)
B[,1:(M+1)] = cbind(1,poly(x,M, raw=TRUE))
B[,(M+2):(M+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^M, 0))
fit = lm(y ~ 0 + B)

#pdf("Figure_radioc_M.pdf")
plot(x,y, ylab="Radio−carbon age", xlab="Calibrated age", main=paste("Degree", M))
abline(v=knots, col="gray", lty=3)
lines(x,predict(fit), lwd=2)
#dev.off()



K = 7
knots = seq(min(x),max(x),length=K+2)[-c(1,K+2)]
M = 3
n = length(x)
B = matrix(NA,ncol=M+K+1, nrow=n)
B[,1:(M+1)] = cbind(1,poly(x,M, raw=TRUE))
B[,(M+2):(M+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^M, 0))
fit = lm(y ~ 0 + B)

#pdf("Figure_radioc_K.pdf")
plot(x,y, ylab="Radio−carbon age", xlab="Calibrated age", main=paste(K+2, "knots"))
abline(v=knots, col="gray", lty=3)
abline(v=range(x), col="gray", lty=3)
lines(x,predict(fit), lwd=2)
#dev.off()

#--- INTERPOLATING NATURAL SPLINE ----------------------------

#pdf("Figure_interpolating_ns.pdf")
f = splinefun(x, y, method = "natural")
plot(x, y, xlim=c(1900,3100))
curve(f(x), add = TRUE, from = 1000, to = 4000, n = 1000, lwd=2) 
#dev.off()

#--- TRUNCATED POWER BASIS ----------------------------

x = (x-min(x)) / (max(x)-min(x))
y = (y-min(y)) / (max(y)-min(y))

K = 13
knots = seq(min(x),max(x),length=K+2)[-c(1,K+2)]
M = 3
n = length(x)
B = matrix(NA,ncol=M+K+1, nrow=n)
B[,1:(M+1)] = cbind(1,poly(x,M, raw=TRUE))
B[,(M+2):(M+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^M, 0))

#pdf("Figure_tpowerbasis.pdf")
op <- par(mfrow=c(3,1), mar = c(0, 4, 0, 0))
plot(x,y, ylim=c(-0.5,1.5), col="white", ylab="Unscaled basis functions Bj(x)")
abline(v=knots, lty=3, col="gray")
for (i in 1:ncol(B)) lines(x,B[,i], col=i, lwd=2)

fit = lm(y ~ 0 + B)
beta_hat = coef(fit)
beta_hat
B_scaled = sapply(1:length(beta_hat),function(j) B[,j]*beta_hat[j])

plot(x,y, ylim=c(-0.5,1.5), col="white", , ylab="Scaled basis functions betaj * Bj(x)")
abline(v=knots, lty=3, col="gray")
for (i in 1:ncol(B_scaled)) lines(x,B_scaled[,i], col=i, lwd=2)

plot(x,y, ylim=c(-0.5,1.5), ylab="Sum betaj * Bj(x)")
abline(v=knots, lty=3, col="gray")
y_hat = apply(B_scaled, 1, sum)
lines(x,y_hat, lwd=2)
par(op)
#dev.off()


#--- NON-STANDARDIZED ----------------------------

rm(list=ls())

library(sm) 
dat = radioc[(radioc$Cal.age > 2000 & radioc$Cal.age < 3000), ]
x = dat$Cal.age
y = dat$Rc.age


K = 13
knots = seq(min(x),max(x),length=K+2)[-c(1,K+2)]
M = 3
n = length(x)
B = matrix(NA,ncol=M+K+1, nrow=n)
B[,1:(M+1)] = cbind(1,poly(x,M, raw=TRUE))
B[,(M+2):(M+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^M, 0))
fit = lm(y ~ 0 + B)
beta_hat = coef(fit)

#pdf("Figure_barbeta.pdf")
barplot(beta_hat)
#dev.off()

kappa(crossprod(B), exact=T)


#---------------------------------------
# B-SPLINES
#---------------------------------------

tpower <- function(x, t, deg){
  (x - t) ^ deg * (x > t)
}

bbase <- function(x, xl, xr, ndx, deg){
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  Delta <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(Delta)
  B
}


xl=min(x)
xr=max(x)
ndx=K+1
bdeg=3

B <- bbase(x, xl, xr, ndx, bdeg)
knots_all <- seq(xl - bdeg * (xr - xl) / ndx, xr + bdeg * (xr - xl) / ndx, by = (xr - xl) / ndx)
#pdf("Figure_Bspline.pdf")
plot(knots_all,rep(0,length(knots_all)),pch=19, ylab=expression(B[j](x)), xlab="x")
abline(v=knots_all, lty=3, col="gray")
for (i in 1:ncol(B)) lines(x,B[,i], col=i, lwd=2)
#dev.off()

rowSums(B)


x = (x-min(x)) / (max(x)-min(x))
y = (y-min(y)) / (max(y)-min(y))

xl=min(x)
xr=max(x)
ndx=K+1
bdeg=3

B <- bbase(x, xl, xr, ndx, bdeg)

#pdf("Figure_Bbasis.pdf")
op <- par(mfrow=c(3,1), mar = c(0, 4, 0, 0))
plot(x,y, ylim=c(-0.5,1.5), col="white", ylab="Unscaled basis functions Bj(x)")
abline(v=knots, lty=3, col="gray")
for (i in 1:ncol(B)) lines(x,B[,i], col=i, lwd=2)

fit = lm(y ~ 0 + B)
beta_hat = coef(fit)
beta_hat
B_scaled = sapply(1:length(beta_hat),function(j) B[,j]*beta_hat[j])

plot(x,y, ylim=c(-0.5,1.5), col="white", , ylab="Scaled basis functions betaj * Bj(x)")
abline(v=knots, lty=3, col="gray")
for (i in 1:ncol(B_scaled)) lines(x,B_scaled[,i], col=i, lwd=2)

plot(x,y, ylim=c(-0.5,1.5), ylab="Sum betaj * Bj(x)")
abline(v=knots, lty=3, col="gray")
y_hat = apply(B_scaled, 1, sum)
lines(x,y_hat, lwd=2)
par(op)
#dev.off()



#---------------------------------------
# NATURAL CUBIC SPLINES
#---------------------------------------

rm(list=ls())

library(sm) 
dat = radioc[(radioc$Cal.age > 2000 & radioc$Cal.age < 3000), ]
x = dat$Cal.age
y = dat$Rc.age

nat_spline_x <- function(x, knots){
  
  K <- length(knots)
  n <- length(x)
  xi <- knots
  
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

M = 3
K = 11
knots = seq(min(x),max(x),length=K+2)[-c(1,K+2)]
n = length(x)
  
#pdf("Figure_nat_spline.pdf")
plot(x,y)
abline(v=knots, lty=3, col="gray")
B <- nat_spline_x(x, knots)
y_hat <- B %*% solve(crossprod(B)) %*% crossprod(B, y)
lines(x,y_hat, lwd=2)
#dev.off()


fit_ns <- lm(y ~ 0+B)
y_hat_ns <- predict(fit_ns, se = TRUE)

X = matrix(NA,ncol=M+K+1, nrow=n)
X[,1:(M+1)] = cbind(1,poly(x,M, raw=TRUE))
X[,(M+2):(M+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^M, 0))
fit_cs <- lm(y ~ 0+X)
y_hat_cs <- predict(fit_cs, se = TRUE)

#pdf("Figure_standard_error.pdf")
plot(x, y_hat_cs$se.fit, type="l", 
     ylim=c(0.02,max(y_hat_cs$se.fit)),
     ylab="Standard error", lwd=2)
lines(x, y_hat_ns$se.fit, lwd=2, col=4)
legend("top", c("cubic spline", "natural cubic spline"), lty=1, col=c(1,4), lwd=2)
#dev.off()

