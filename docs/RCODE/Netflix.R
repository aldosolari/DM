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

# convertire i dati mancanti in NA
is.na(train) = train==0
is.na(test) = test==0

# albero di regressione cp = 0.02
library(rpart)
fit12 <- rpart (y~.,data=train,cp=0.02)
summary(fit12) 
plot(fit12)
text(fit12, all=T, use.n=T)
yhat12 <- predict(fit12, newdata=test)
sqrt(mean((test$y-yhat12)^2))

# albero di regressione cp = 0.0001
fit13= rpart(y~.,data=train,cp=0.0001)
tree.res = printcp(fit13)
#1-SE rule
chosen.prune = min((1:dim(tree.res)[1]) [tree.res[,"xerror"] < min(tree.res[,"xerror"]+tree.res[,"xstd"])])
tree.prune = prune(fit13, cp=tree.res[chosen.prune,"CP"])
# predict validation
yhat13= predict(tree.prune, newdata=test)
sqrt(mean((test$y-yhat13)^2))

#---------- Bagging

err.small=NULL
pred = numeric(dim(test)[1])
for (i in 1:100){
  use = sample(nrow(train),nrow(train),rep=T)
  fit = rpart(y~.,data=train[use,], cp=0.01)
  pred = pred + predict(fit, newdata=test)
  if (i%%10==0) cat (i, sqrt(mean ((pred/i-test$y)^2)),"\n")
  err.small = c(err.small, sqrt(mean ((pred/i-test$y)^2)))
}

plot (1:100, err.small, xlab="Bagging iterations", ylab="RMSE", ylim=c(0.75, 0.85), type="l", lty=2)
