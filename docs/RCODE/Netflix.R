rm(list=ls())
#==========================
# NETFLIX DATA
#==========================

# Adapted from Saharon Rosset code for Statistical/Machine Learning course https://www.tau.ac.il/~saharon/StatLearn.html

#--- ESERCIZIO 1 ----------------------------

# importazione dei dati
PATH <- "https://raw.githubusercontent.com/aldosolari/DM/master/docs/DATA/"
X <- read.table(paste(PATH,"Train_ratings_all.dat", sep=""))
titles <- read.table(paste(PATH,"Movie_titles.txt", sep=""), sep=",")
names(X) <- substr(as.character(titles[,2]),1,20)
y <- read.table(paste(PATH,"Train_y_rating.dat", sep=""))
names(y) <- "y"

# % di valori mancanti per film?
barplot(sort(apply(X==0, 2, mean)), horiz=T, names.arg=F, xlab="% di dati mancati", ylab="Movie") 

# Come è stato valutato il film Miss Congeniality?
table(y)

plot(table(y)/nrow(X), ylab="Frequenza relativa", xlab="Rating")

# Valutazione media rispetto ai film senza dati mancanti?
sort(apply(data.frame(X[,1:14],y),2,mean))

# Quali sono i film (senza dati mancanti) maggiormente correlati con Miss Congeniality?
cor(X[,1:14],y)

#--- ESERCIZIO 2 ----------------------------

# Divisione in training e test

m <- 2000
n <- nrow(X) - m
set.seed(123)
test.id <- sample(n+m,m)
test <- data.frame(y=y[test.id,], X[test.id,])
train <- data.frame(y=y[-test.id,], X[-test.id,])

# Modello nullo

yhat0 <- mean(train$y)
sqrt(mean((test$y-yhat0)^2))

# Modello lineare (con i valori mancanti codificati = 0)

fit1 <- lm(y ~ ., train)
sqrt(mean((train$y-fitted(fit1))^2))
yhat1 <- predict(fit1, newdata=test)
summary(yhat1)
sqrt(mean((test$y-yhat1)^2))
# un piccolo accorgimento
yhat2 <- pmin(yhat1,5)
sqrt(mean((test$y-yhat2)^2))

# Modello lineare con i film senza dati mancanti

fml3 = paste("y ~", paste(names(train[,2:15]), collapse = "+"))
fit3 <- lm(fml3, train)
summary(fit3)$coefficients
yhat3 <- predict(fit3, newdata=test)
summary(yhat3)
sqrt(mean((test$y-yhat3)^2))


#--- ESERCIZIO 3 ----------------------------

q = 5
Y = matrix(0, nrow=n, ncol=q)
for (j in 1:q) Y[train$y==j,j]=1
X.tr = as.matrix(cbind(1,X[-test.id,]))
Bhat = solve(t(X.tr) %*% X.tr) %*% t(X.tr) %*% Y
X.te = as.matrix(cbind(1,X[test.id,]))
Yhat = X.te %*% Bhat
head(Yhat)
# previsione con lo score più probabile
yhat4 <- (1:q)[apply(Yhat, 1, which.max)]
table(yhat4)
sqrt(mean((test$y-yhat4)^2))

# previsione con lo score atteso
yhat5 = Yhat %*% 1:q
summary(as.numeric(yhat5))
sqrt(mean((test$y-yhat5)^2))

#--- ESERCIZIO 4 ----------------------------


#---------- multiclass logistic regression
library(nnet)
fit6 = multinom(as.factor(y)~., family = multinomial, data=train, maxit=300)

# previsione con lo score più probabile
yhat6 = predict(fit6, newdata=test)
summary(yhat6)
sqrt(mean((test$y-as.numeric(yhat6))^2))

# previsione con lo score atteso
phat7 = predict(fit6, newdata=test,type="prob")
yhat7 = as.numeric(phat7%*%(1:5))
summary(yhat7)
sqrt(mean((test$y-yhat7)^2))

#---------- LDA

library(MASS)
fit8 = lda(as.factor(y)~.,data=train)
pred8 = predict(fit8, newdata=test)

# previsione con lo score più probabile
yhat8 = as.numeric(pred8$class)
table(yhat8)
sqrt(mean((test$y-yhat8)^2))

# previsione con lo score atteso
yhat9 = as.numeric(pred8$posterior%*%(1:5))
summary(yhat9)
sqrt(mean((test$y-yhat9)^2))

#---------- QDA

fit10 = qda(as.factor(y)~.,data=train)
pred10 = predict(fit10, newdata=test)

# previsione con lo score più probabile
yhat10 = as.numeric(pred10$class)
table(yhat10)
sqrt(mean((test$y-yhat10)^2))

# previsione con lo score atteso
yhat11 = as.numeric(pred10$posterior%*%(1:5))
summary(yhat11)
sqrt(mean((test$y-yhat11)^2))

#---------- ALBERI DI REGRESSIONE ------

# convertire i dati mancanti in NA
is.na(train) = train==0
is.na(test) = test==0

# albero di regressione cp = 0.02
library(rpart)
fit12 <- rpart(y~.,data=train,cp=0.02)
#summary(fit12) 
plot(fit12)
text(fit12, all=T, use.n=T)
yhat12 <- predict(fit12, newdata=test)
sqrt(mean((test$y-yhat12)^2))

# albero di regressione cp = 0.0001
fit13 = rpart(y~.,data=train,cp=0.0001)
tree.res = printcp(fit13)
#1-SE rule
chosen.prune = min((1:dim(tree.res)[1]) [tree.res[,"xerror"] < min(tree.res[,"xerror"]+tree.res[,"xstd"])])
tree.prune = prune(fit13, cp=tree.res[chosen.prune,"CP"])
# predict validation
yhat13 = predict(tree.prune, newdata=test)
sqrt(mean((test$y-yhat13)^2))

#---------- Bagging

B = 10

err.small=NULL
pred = numeric(dim(test)[1])
for (i in 1:B){
  use = sample(nrow(train),nrow(train),rep=T)
  fit = rpart(y~.,data=train[use,], cp=0.01)
  pred = pred + predict(fit, newdata=test)
  if (i%%10==0) cat (i, sqrt(mean ((pred/i-test$y)^2)),"\n")
  err.small = c(err.small, sqrt(mean ((pred/i-test$y)^2)))
}

#pdf("Figure_bagging.pdf")
plot (1:B, err.small, xlab="Bagging iterations", ylab="RMSE", ylim=c(0.75, 0.85), type="l", lty=2)
#dev.off()

#--- ESERCIZIO 5 ----------------------------

#---------- algoritmo L2-boosting con alberi (missing=NA)
lambda = 0.01
y.now = train$y-mean(train$y)
err.boost=err.tr.boost=NULL
pred.boost = numeric(nrow(test))+mean(train$y)
fitted.boost = numeric(nrow(train))+mean(train$y)
for (i in 1:B){
  tree.mod= rpart(y.now~.-y,data=train,maxdepth=2,cp=0.00001)
  pred.boost = pred.boost + lambda*predict(tree.mod, newdata=test)
  fitted.boost = fitted.boost  + lambda*predict(tree.mod)
  y.now = train$y-fitted.boost
  if (i%%10==0) cat (i, "train:", sqrt(mean ((fitted.boost-train$y)^2)), " test:", sqrt(mean ((pred.boost-test$y)^2)),"\n")
  err.boost = c(err.boost, sqrt(mean ((pred.boost-test$y)^2)))
  err.tr.boost = c(err.tr.boost, sqrt(mean ((fitted.boost-train$y)^2)))
  }

#pdf("Figure_boosting.pdf")
plot ((1:B), err.boost, main="Boosting with alpha=0.01, maxdepth=2 and miss=NA", xlab="number of steps", ylab="RMSE",type="l")
lines (1:B, err.tr.boost, lty=2)
legend ("topright",legend=c("boost - valid", "boost - train"), lty=c(1,2))
#dev.off()

#--- ESERCIZIO 6 ----------------------------

X.now = scale(X[,1:14], center=T, scale = F)
colnames(X.now) = titles$V2[1:14]

pca = princomp(X.now, cor=F)

# Scores for Miss Congeniality da utenti con valori alti su PC1
table(y[order(pca$scores[,1])[1:100],1])
# Scores for Miss Congeniality da utenti con valori bassi su PC1
table(y[order(pca$scores[,1],decreasing=T)[1:100],1])

tr.sc = data.frame(pca$scores[-test.id,], y=y[-test.id,])
te.sc = data.frame(pca$scores[test.id,], y=y[test.id,])

######### PCA regression on various number of PCAs
for (npca in 1:14){
  pcr=lm(y~., data=tr.sc[,c(1:npca,15)])
  te.pred = predict(pcr,newdata=te.sc[,c(1:npca,15)])
  cat (npca, "PCA, validation error:", sqrt(mean((test$y-te.pred)^2)),"\n")}

#--- ESERCIZIO 7 ----------------------------

# dataset ridotto
N <- 1000
set.seed(123)
ix <- sample(1:10000,size=N)
X.small <- as.matrix(X[ix,1:14])
n <- nrow(X.small)
p <- ncol(X.small)

# creazione dati mancanti
nomit <- n*p*0.15
ina <- sample(seq(n), nomit, replace = TRUE)
inb <- sample(1:p, nomit, replace = TRUE)
Xna <- X.small
index.na <- cbind(ina, inb)
Xna[index.na] <- NA
ismiss <- is.na(Xna)

# dati centrati
Z <- scale(Xna, center=T, scale = F)
means =attr(Z,"scaled:center")

Zhat <- Z
zbar <- colMeans(Z, na.rm = TRUE)
Zhat[index.na] <- zbar[inb]

fit.svd <- function(X, q = 1){
  DVS <- svd(X)
  with(DVS, 
       u[,1:q, drop = FALSE] %*%
         (d[1:q] * t(v[,1:q, drop = FALSE]))
  )
}


n_iter <- 50 # numero di iterazioni
for (iter in 1:n_iter){
  # Step 2.a
  Zapp <- fit.svd(Zhat, q = 1)
  # Step 2.b
  Zhat[ismiss] <- Zapp[ismiss]
  # Step 2.c
  e <- mean(((Z - Zapp)[!ismiss])^2)
  cat("Iter :", iter, " Errore :", e, "\n")
}

Xhat = Zhat + matrix(rep(means,each=nrow(Z)), ncol=ncol(Z))
cor(Xhat[ismiss], X.small[ismiss])
plot(Xhat[ismiss], X.small[ismiss], asp=1)
abline(a=0,b=1)

mean((Xhat[ismiss] - X.small[ismiss])^2)


#---- Modello utente - film
W = cbind(expand.grid(rownames(Xna),colnames(Xna)),y=c(Xna))
fit = lm(y ~ Var1 + Var2, W, subset =!is.na(W$y))
yhat <- predict(fit, newdata=W[is.na(W$y),])
cor(yhat, X.small[ismiss])
plot(yhat, X.small[ismiss], asp=1)
abline(a=0,b=1)
abline(v=5)

mean((yhat - X.small[ismiss])^2)





