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

M = 10
fit = lm(y ~ poly(x,degree=M, raw=T))
#pdf("Figure_polyM.pdf")
plot(x,y, main=paste("Regressione polinomiale di grado", M))
lines(x,fitted(fit))
#dev.off()

B = model.matrix(fit)
S = B %*% solve(crossprod(B)) %*% t(B)

#pdf("Figure_hatmatrixM.pdf")
filled.contour(main=paste("Regressione polinomiale di grado", M), apply(S, 2, rev), asp=1, xaxt="n", yaxt="n", xlab="columns", ylab="rows", levels=c(min(S),-0.005,0.005,max(S)), col=c("red", "white", "blue"))
#dev.off()

kappa(crossprod(B), exact=T)

#--- POLINOMI0 A TRATTI ----------------------------

M = 2
K = 12 # numero di intervalli = K+1
knots = seq(-1,1,length=K+2)[-c(1,K+2)]
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
filled.contour(main=paste("Polinomio a tratti di grado", M),apply(S, 2, rev), asp=1, xaxt="n", yaxt="n", xlab="columns", ylab="rows", levels=c(min(S),-0.005,0.005,max(S)), col=c("red", "white", "blue"))
#dev.off()


#--- SPLINE DI REGRESSIONE ----------------------------

# costruzione della spline con funzioni base di potenza troncata
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

#--- RADIOCARBON DATA ----------------------------

rm(list=ls())

library(sm) 
dat = radioc[(radioc$Cal.age > 2000 & radioc$Cal.age < 3000), ]
x = dat$Cal.age/1000
y = dat$Rc.age/1000
#pdf("Figure_radiocarbon.pdf")
plot(y ~ x, xlab="calibrated age (1000)", ylab="radio-carbon age (1000)")
#dev.off()


#--- SCELTA DEL GRADO M ----------------------------

K = 6
knots = seq(min(x),max(x),length=K+2)[-c(1,K+2)]
M = 3
n = length(x)
B = matrix(NA,ncol=M+K+1, nrow=n)
B[,1:(M+1)] = cbind(1,poly(x,M, raw=TRUE))
B[,(M+2):(M+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^M, 0))
fit = lm(y ~ 0 + B)

#pdf("Figure_radioc_M.pdf")
plot(x,y, ylab="Radio-carbon age", xlab="Calibrated age", main=paste("Degree", M))
abline(v=knots, col="gray", lty=3)
lines(x,predict(fit), lwd=2)
#dev.off()


#--- SCELTA DEL NUMERO K ----------------------------

K = 7
knots = seq(min(x),max(x),length=K+2)[-c(1,K+2)]
M = 3
n = length(x)
B = matrix(NA,ncol=M+K+1, nrow=n)
B[,1:(M+1)] = cbind(1,poly(x,M, raw=TRUE))
B[,(M+2):(M+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^M, 0))
fit = lm(y ~ 0 + B)

#pdf("Figure_radioc_K.pdf")
plot(x,y, ylab="Radio-carbon age", xlab="Calibrated age", main=paste(K, "internal knots"))
abline(v=knots, col="gray", lty=3)
abline(v=range(x), col="gray", lty=3)
lines(x,predict(fit), lwd=2)
#dev.off()

#--- INTERPOLATING NATURAL SPLINE ----------------------------

#pdf("Figure_interpolating_ns.pdf")
f = splinefun(x, y, method = "natural")
plot(x, y, xlim=c(1.9,3.1))
curve(f(x), add = TRUE, from = 1, to = 4, n = 1000, lwd=2) 
#dev.off()

#--- CONSTRUCTION B-SPLINE (STANDARDIZED) -----

xx = (x-min(x)) / (max(x)-min(x))
yy = (y-min(y)) / (max(y)-min(y))

K = 13
knots = seq(min(xx),max(xx),length=K+2)[-c(1,K+2)]
M = 3
n = length(xx)
B = matrix(NA,ncol=M+K+1, nrow=n)
B[,1:(M+1)] = cbind(1,poly(xx,M, raw=TRUE))
B[,(M+2):(M+K+1)] = sapply(1:K, function(k) ifelse(xx >= knots[k], (xx-knots[k])^M, 0))

#pdf("Figure_tpowerbasis.pdf")
op <- par(mfrow=c(3,1), mar = c(0, 4, 0, 0))
plot(xx,yy, ylim=c(-0.5,1.5), col="white", ylab="Unscaled basis functions Bj(x)")
abline(v=knots, lty=3, col="gray")
for (i in 1:ncol(B)) lines(xx,B[,i], col=i, lwd=2)

fit = lm(yy ~ 0 + B)
beta_hat = coef(fit)
B_scaled = sapply(1:length(beta_hat),function(j) B[,j]*beta_hat[j])

plot(xx,yy, ylim=c(-0.5,1.5), col="white", , ylab="Scaled basis functions betaj * Bj(x)")
abline(v=knots, lty=3, col="gray")
for (i in 1:ncol(B_scaled)) lines(xx,B_scaled[,i], col=i, lwd=2)

plot(xx,yy, ylim=c(-0.5,1.5), ylab="Sum betaj * Bj(x)")
abline(v=knots, lty=3, col="gray")
y_hat = apply(B_scaled, 1, sum)
lines(xx,y_hat, lwd=2)
par(op)
#dev.off()


#--- NON-STANDARDIZED DATA ----------------------------

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




#--- CONSTRUCTION B-SPLINE (STANDARDIZED) -----

xl=min(xx)
xr=max(xx)
ndx=K+1
bdeg=3
BB <- bbase(xx, xl, xr, ndx, bdeg)

#pdf("Figure_Bbasis.pdf")
op <- par(mfrow=c(3,1), mar = c(0, 4, 0, 0))
plot(xx,yy, ylim=c(-0.5,1.5), col="white", ylab="Unscaled basis functions Bj(x)")
abline(v=knots, lty=3, col="gray")
for (i in 1:ncol(B)) lines(xx,BB[,i], col=i, lwd=2)

fit = lm(yy ~ 0 + BB)
beta_hat = coef(fit)
BB_scaled = sapply(1:length(beta_hat),function(j) BB[,j]*beta_hat[j])

plot(xx,yy, ylim=c(-0.5,1.5), col="white", , ylab="Scaled basis functions betaj * Bj(x)")
abline(v=knots, lty=3, col="gray")
for (i in 1:ncol(BB_scaled)) lines(xx,BB_scaled[,i], col=i, lwd=2)

plot(xx,yy, ylim=c(-0.5,1.5), ylab="Sum betaj * Bj(x)")
abline(v=knots, lty=3, col="gray")
y_hat = apply(BB_scaled, 1, sum)
lines(xx,y_hat, lwd=2)
par(op)
#dev.off()



#---------------------------------------
# NATURAL CUBIC SPLINES
#---------------------------------------


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
  

B <- nat_spline_x(x, knots)
y_hat <- B %*% solve(crossprod(B)) %*% crossprod(B, y)

fit_ns <- lm(y ~ 0+B)
y_hat_ns <- predict(fit_ns, se = TRUE)

X = matrix(NA,ncol=M+K+1, nrow=n)
X[,1:(M+1)] = cbind(1,poly(x,M, raw=TRUE))
X[,(M+2):(M+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^M, 0))
fit_cs <- lm(y ~ 0+X)
y_hat_cs <- predict(fit_cs, se = TRUE)

#pdf("Figure_nat_spline.pdf")
plot(x,y)
abline(v=knots, lty=3, col="gray")
lines(x,y_hat, lwd=2, col=4)
lines(x,y_hat_cs$fit, lwd=2)
legend("top", c("natural cubic spline", "cubic spline"), lty=1, col=c(4,1), lwd=2)
#dev.off()


#pdf("Figure_standard_error.pdf")
plot(x, y_hat_cs$se.fit, type="l", 
     ylab="Standard error", lwd=2)
lines(x, y_hat_ns$se.fit, lwd=2, col=4)
legend("top", c("cubic spline", "natural cubic spline"), lty=1, col=c(1,4), lwd=2)
#dev.off()

