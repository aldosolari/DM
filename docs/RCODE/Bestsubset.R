rm(list=ls())
#==========================
# BEST SUBSET SELECTION
#==========================

#--- Hitters data ----------------------------

library(ISLR)
data(Hitters)
Hitters = Hitters[complete.cases(Hitters),]
Hitters[,"League"]=(Hitters[,"League"]=="A")*1
Hitters[,"Division"]=(Hitters[,"Division"]=="E")*1
Hitters[,"NewLeague"]=(Hitters[,"NewLeague"]=="A")*1
set.seed(123)
n = 163
train.id = sample(nrow(Hitters),n)
train = Hitters[train.id,]
names(train)[19] = "y" 
X = as.matrix(train[,-19])
y = train$y
p = ncol(X)
test = Hitters[-train.id,]
names(test)[19] = "y" 
X.star = as.matrix(test[,-19])
y.star = test$y
m = nrow(X.star)

#--- Esercizio 1 ----------------------------

# Full model
fit.full = lm(y ~ ., train)
RMSE.full = sqrt(mean( (predict(fit.full, newdata=test) - y.star )^2 ))
RMSE.full

# Best subset
library(leaps)
fit.bests <- regsubsets(y~.,train, nvmax=p)
summary.bests<-summary(fit.bests)
head(summary.bests$outmat, 10)

# Best Cp
#pdf("Figure_bestCp.pdf")
plot(summary.bests$cp, xlab="k", ylab="Cp", type="b")
#dev.off()

# best BIC
#pdf("Figure_bestbic.pdf")
plot(fit.bests, scale="bic")
#dev.off()

# function predict for regsubsets
predict.regsubsets =function(object ,newdata ,id ,...){
  form=as.formula(object$call[[2]])
  mat=model.matrix(form, newdata)
  coefi =coef(object, id=id)
  xvars =names(coefi)
  mat[,xvars]%*%coefi
}

yhat.bestCp = predict.regsubsets(fit.bests, newdata=test, id=which.min(summary.bests$cp))
RMSE.Cp = sqrt(mean( (yhat.bestCp - y.star)^2 ))
RMSE.Cp

yhat.bestBIC = predict.regsubsets(fit.bests, newdata=test, id=which.min(summary.bests$bic))
RMSE.BIC = sqrt(mean( (yhat.bestBIC - y.star)^2 ))
RMSE.BIC

# Forward with AIC stopping rule
fit.null = lm(y ~ 1, train)
fit.fwdAIC = step(fit.null, scope=list(upper=fit.full), direction="forward", k=2, trace=0)
summary(fit.fwdAIC)$coeff
yhat.fwdAIC = predict(fit.fwdAIC, newdata=test)
RMSE.fwdAIC = sqrt(mean( (yhat.fwdAIC - y.star)^2 ))
RMSE.fwdAIC

data.frame(
  fit = c("full", "best cp", "best bic", "best fwd aic", "best va", "best cv"),
  RMSE = c(RMSE.full, RMSE.Cp, RMSE.BIC, RMSE.fwdAIC, RMSE.bestVa, RMSE.bestCV)
)

#--- Esercizio 2 ----------------------------

# Validation set
set.seed(123)
is.va = sample(c(T,F), n, replace = TRUE)
err.va = vector()
fit = regsubsets(y ~.,data=train[!is.va,],nvmax=p)
for (j in 1:p){
  yhat=predict(fit, train[is.va,], id=j)
  err.va[j] = sqrt(mean( (train$y[is.va]-yhat)^2 ))
}
names(err.va) <- 1:p
err.va
yhat.bestVa = predict.regsubsets(fit.bests, newdata=test, id=which.min(err.va))
RMSE.bestVa = sqrt(mean( (yhat.bestVa - y.star)^2 ))
RMSE.bestVa

# K-fold Cross-Validation
set.seed(123)
K = 10
folds = sample(1:K, n, replace =TRUE)
KCV = matrix(NA, K, p)
for (k in 1:K){
  fit_k = regsubsets(y ~.,data=train[folds!=k,],nvmax=p)
  for (j in 1:p){
    yhat_k=predict(fit_k, train[folds==k,], id=j)
    KCV[k,j]=mean( (train$y[folds==k]-yhat_k)^2 )
  }
}

plot(1:p,apply(KCV,2,mean), type="b", ylab="CV error", xlab="k")

yhat.bestCV = predict.regsubsets(fit.bests, newdata=test, id=which.min(apply(KCV,2,mean)))
RMSE.bestCV = sqrt(mean( (yhat.bestCV - y.star)^2 ))
RMSE.bestCV

#--- § 7.10.2 HTF ----------------------------

set.seed(123)
n = 50
p = 5000
y = c(rep(0,n/2),rep(1,n/2))
X = matrix(rnorm(p*n), ncol=p)
# compute p correlations between each predictor and the response
cors = apply(X,2, function(x) cor(y,x))
# columns of the best 100 predictors
colbest100 = sort(-abs(cors), index.return=T)$ix[1:100]
# best 100 predictors
Xbest = X[,colbest100]
# CV
require(class)
K<-5
set.seed(123)
folds <- sample( rep(1:K,length=n) )
Err.CV = vector()
for (k in 1:K){
  out = which(folds==k)
  # predict the held-out samples using k nearest neighbors
  pred <- knn(train = Xbest[ -out, ],test = Xbest[out, ], cl = y[-out], k = 1)
  # % of misclassified samples
  Err.CV[k] = mean( y[out] != pred)
}
mean(Err.CV)

set.seed(123)
n = 50
p = 5000
y = c(rep(0,n/2),rep(1,n/2))
X = matrix(rnorm(p*n), ncol=p)
require(class)
K<-5
set.seed(123)
folds <- sample( rep(1:K,length=n) )
Err.CV = vector()
for (k in 1:K){
  out = which(folds==k)
  cors = apply(X[-out, ],2, function(x) cor(y[-out],x))
  colbest100 = sort(-abs(cors), index.return=T)$ix[1:100]
  Xbest = X[,colbest100]
  pred <- knn(train = Xbest[-out, ],
              test = Xbest[out, ],
              cl = y[-out], k = 1)
  Err.CV[k] = mean( y[out] != pred)
}
mean(Err.CV)

