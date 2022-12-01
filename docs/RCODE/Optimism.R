rm(list=ls())
#==========================
# OTTIMISMO
#==========================

#---- E[MSE.tr],  E[MSE.te] -------------

sigmatrue = 0.01
ftrue <- c(0.4342,0.4780,0.5072,0.5258,0.5369,0.5426,0.5447,0.5444,0.5425,0.5397,0.5364,0.5329,0.5294,0.5260,0.5229,0.5200,0.5174,0.5151,0.5131,0.5113,0.5097,0.5083,0.5071,0.5061,0.5052,0.5044,0.5037,0.5032,0.5027,0.5023)
x = seq(.5,3,length=30)

n <- length(x)
ds = 1:15
ps = ds + 1
Bias2s = sapply(ps, function(p) 
  mean( ( ftrue - fitted(lm(ftrue ~ poly(x,degree=(p-1)))) )^2 )
)
Vars = ps*(sigmatrue^2)/n
AveMSE.te = Bias2s+Vars+sigmatrue^2
AveMSE.tr = AveMSE.te - (2*ps*sigmatrue^2)/n

#pdf("Figure_AveMSE.pdf")
plot(ds, AveMSE.tr, type="b", lwd=2, ylab="ErrF", xlab="d")
legend("topright",c("E(MSE.tr)","E(MSE.te)"), col=c(1,2), pch=19)
lines(ds, AveMSE.te,  type="b", lwd=2, col=2)
#dev.off()

#---- MSE.tr + Opt -------------


library(readr)
PATH <- "http://azzalini.stat.unipd.it/Book-DM/yesterday.dat"
df <- read_table(PATH)
train <- data.frame(x=df$x, y=df$y.yesterday)
test <- data.frame(x=df$x, y=df$y.tomorrow)

# MSE.tr
fun <- function(d) if (d==0) lm(y~1, train) else lm(y~poly(x,degree=d), train)
fits <- lapply(ds, fun)
MSEs.tr <- unlist( lapply(fits, deviance) )/n
hatMSEs.te = MSEs.tr + (2*sigmatrue^2*ps)/n

#pdf("Figure_optimism.pdf")
plot(ds, MSEs.tr, type="b", xlab="d", ylab="ErrF")
lines(ds, hatMSEs.te, type="b", col=2)
legend("topright",c("MSE.tr","MSE.tr + Opt"), col=c(1,2), pch=19)
#dev.off()

#---- Cp di Mallows -------------

hatsigma2 = (n*MSEs.tr)/(n-ps)
Cps = MSEs.tr + (2*hatsigma2*ps)/n
#pdf("Figure_MallowsCp.pdf")
plot(ds, MSEs.tr, type="b", xlab="d", ylab="ErrF")
lines(ds, Cps, type="b", col=2)
legend("topright",c("MSE.tr","MSE.tr + hatOpt"), col=c(1,2), pch=19)
#dev.off()

#pdf("Figure_sigma2.pdf")
plot(ds, hatsigma2, type="b", xlab="d", ylab="Sigma2")
abline(h=sigmatrue^2, col=4)
#dev.off()

#---- AIC e BIC -------------

AICs <- unlist( lapply(fits, AIC) )
BICs <- unlist( lapply(fits, BIC) )

#pdf("Figure_AICBIC.pdf")
plot(ds, AICs, type="b", col=5, ylab="Criteri basati sull'informazione", xlab="d")
lines(ds, BICs, type="b", col=6)
legend("topright",c("AIC","BIC"), col=c(5,6), lty=1)
#dev.off()

#---- Convalida incrociata -------------


# K-fold CV
d = 3
K = 5
set.seed(123)
# create folds
folds <- sample( rep(1:K,length=n) )
# initialize vector
KCV <- vector()
# loop
for (k in 1:K){
  fit <- lm(y~poly(x,degree=d), train, subset=which(folds!=k))
  x.out <- train$x[which(folds==k)]
  yhat <- predict(fit, newdata=list(x=x.out))
  y.out <- train$y[which(folds==k)]
  KCV[k]<- mean( ( y.out - yhat )^2 )
}
# KCV estimate 
mean(KCV)

#---- cv.glm -------------

ds = 1:12; ps = ds+1
library(boot)
set.seed(123)
KCV = sapply(ds, function(d) 
  cv.glm(train, glm(y~poly(x,degree=d), train, family = gaussian), K=K )$delta[1] )

#pdf("Figure_KCV.pdf")
plot(ds, KCV, type="b", log="y", xlab="d")
#dev.off()

#---- LOOCV -------------

oneout <- vector()
for (i in 1:n){
  fit_i <- lm( y~poly(x,degree=d), data=train[-i,])
  yhat_i <- predict(fit_i, newdata=data.frame(x=train$x[i]) )
  oneout[i] <- ( train$y[i] -  yhat_i )^2
}
mean(oneout)

#---- LOOCV per il modello lineare -------------


fit <-  lm(y~poly(x,d), train)
X <- model.matrix(fit)
H <- X %*% solve(t(X)%*% X) %*% t(X)
mean(
  ( (train$y - predict(fit)) / (1-diag(H))  )^2 
) 

LOOCV = sapply(ds, function(d) 
  cv.glm(train, glm(y~poly(x,degree=d), 
                    train, family = gaussian) )$delta[1] )

#pdf("Figure_LOOCV.pdf")
plot(ds, LOOCV, type="b", xlab="d")
#dev.off()

#---- Metodo della convalida incrociata generalizzata -------------


fun <- function(d) if (d==0) lm(y~1, train) else lm(y~poly(x,degree=d), train)
fits <- lapply(ds, fun)
MSEs.tr <- unlist( lapply(fits, deviance) )/n
GCV = MSEs.tr/(1-(ps)/n )^2

#pdf("Figure_GCV.pdf")
plot(ds, GCV, type="b", xlab="d")
#dev.off()