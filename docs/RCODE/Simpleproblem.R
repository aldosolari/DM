rm(list=ls())
#==========================
# UN SEMPLICE PROBLEMA-TIPO
#==========================

#---------------
# yesterday data
#---------------

library(readr)
PATH <- "http://azzalini.stat.unipd.it/Book-DM/yesterday.dat"
df <- read_table(PATH)
train <- data.frame(x=df$x, y=df$y.yesterday)
test <- data.frame(x=df$x, y=df$y.tomorrow)

#pdf("Figure_scatterplot.pdf")
plot( y ~ x , train)
#dev.off()

# Esercizio 1 ------------

d <- 3
fit <- lm(y ~ poly(x, degree=d, raw=T), train)
yhat <- fitted(fit)

#pdf("Figure_polynomialfit.pdf")
plot(y ~ x, train)
points(yhat ~ x, train, pch=19)
lines(yhat ~ x, train, lwd=2)
#dev.off()

MSE.tr <- mean((train$y - yhat)^2)
MSE.tr

# Esercizio 2 ------------

n <- nrow(train)
ds = 0:(n-1)
ps = ds + 1
fun <- function(d) if (d==0) lm(y~1, train) else lm(y~poly(x, degree=d, raw=T), train)
fits <- sapply(ds, fun)
MSEs.tr <- unlist( lapply(fits, deviance) )/n

#pdf("Figure_MSE_tr.pdf")
plot(ds, MSEs.tr, type="b", xlab="d", ylab="MSE.tr")
#dev.off()

# Overfitting ------------

fit_15 <- lm( y ~ poly(x, degree=15, raw=F), train)
yhat <- predict(fit_15, newdata=test)
#pdf("Figure_overfitting_tr.pdf")
plot( y ~ x , train, main ="d = 15 sui dati di training")
lines( yhat ~ x, train, lwd=2)
#dev.off()

#pdf("Figure_overfitting_te.pdf")
plot( y ~ x, test, col=4, main ="d = 15 sui dati di test")
lines( yhat ~ x, test, lwd=2)
#dev.off()

# MSE.te ------------

yhats <- lapply(fits, predict)
MSEs.te <- unlist(lapply(yhats, 
                         function(yhat) mean((test$y - yhat)^2)
))
#pdf("Figure_MSE_te.pdf")
plot(ds, MSEs.te, type="b", col=4, xlab="d", ylab="MSE.te")
#dev.off()

# Matrice del disegno ------------

fit <- lm( y ~ poly(x, degree=3, raw=TRUE), train)
X = model.matrix(fit)
colnames(X) = c("Intercept","x","x^2","x^3")
head(X)

# Esercizio 4 ------------

fit_12_raw <- lm( y ~ poly(x, degree=12, raw=TRUE), train)
coef(fit_12_raw)
X <- model.matrix(fit_12_raw)
qr(X)$rank

#fit_24 <- lm( y ~ poly(x, degree=24, raw=FALSE), train)
xbar <- mean(train$x)
x <- train$x - xbar
X <- outer(x, 0L:24, `^`)
QR <- qr(X)
QR$rank

# Polinomi ortogonali ------------

fit <- lm( y ~ poly(x, degree=3, raw=FALSE), train)
X = model.matrix(fit)
colnames(X) = c("Intercept","x1","x2","x3")
round( t(X) %*% X, 8)

X[,1] = 1/nrow(X)
beta_hat <- crossprod(X, train$y)
beta_hat
coef(lm(y ~ poly(x, degree=3), train))



