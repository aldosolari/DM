rm(list=ls())
#==========================
# COMPROMESSO DISTORSIONE-VARIANZA
#==========================

#--- Il segnale e il rumore ----------

ftrue <- c(0.4342,0.4780,0.5072,0.5258,0.5369,0.5426,0.5447,0.5444,0.5425,0.5397,0.5364,0.5329,0.5294,0.5260,0.5229,0.5200,0.5174,0.5151,0.5131,0.5113,0.5097,0.5083,0.5071,0.5061,0.5052,0.5044,0.5037,0.5032,0.5027,0.5023)
x = seq(.5,3,length=30)

#pdf("Figure_signal.pdf")
plot(x,ftrue, type="l", col=4, ylab="f(x)", lwd=2)
#dev.off()

#pdf("Figure_noise.pdf")
plot(x,ftrue, type="l", col=4, ylab="f(x)", lwd=2)
set.seed(12)
x0 = x[15]
y0 = ftrue[15]+rnorm(1,0,sd=0.01)
points(x0,y0)
segments(x0,y0,x0,ftrue[15], col=2)
#dev.off()

#--- Sovra-adattamento (overfitting) ----------

library(readr)
PATH <- "http://azzalini.stat.unipd.it/Book-DM/yesterday.dat"
df <- read_table(PATH)
train <- data.frame(x=df$x, y=df$y.yesterday)

fit_15 <- lm( y ~ poly(x, degree=15, raw=F), train)
yhat <- predict(fit_15)
#pdf("Figure_overfitting.pdf")
plot( y ~ x , train, main ="d = 15 sui dati di training")
lines(x, ftrue, col=4, lwd=2)
lines( yhat ~ x, train, lwd=2)
#dev.off()

#--- Distorsione e varianza ----------

rm(list=ls())

sigmatrue = 0.01
x = seq(.5,3,length=30)
n = length(x)
ftrue = c(0.4342,0.4780,0.5072,0.5258,0.5369,
          0.5426,0.5447,0.5444,0.5425,0.5397,
          0.5364,0.5329,0.5294,0.5260,0.5229,
          0.5200,0.5174,0.5151,0.5131,0.5113,
          0.5097,0.5083,0.5071,0.5061,0.5052,
          0.5044,0.5037,0.5032,0.5027,0.5023)

d = 5
X <- model.matrix(lm(ftrue ~ poly(x,degree=d)))
invXtX <- solve(crossprod(X))

Bias2 <- ( apply(X, 1, function(x) 
  x %*% invXtX %*% t(X) %*% ftrue) - ftrue )^2

Var = apply(X, 1, function(x) 
  sigmatrue^2 * t(x) %*% invXtX %*% x
)

#pdf("Figure_bias2var.pdf")
barplot(Bias2+Var, ylab="Bias2 + Var", names.arg=round(x,1), main=paste("d = ",d))
barplot(Var,add=T, col=1, names.arg=" ")
legend("topright", c("Bias2","Var"), col=c("gray",1), pch=c(19,19))
#dev.off()

#--- Errore di previsione atteso ----------

sigmatrue^2 + mean( Bias2 ) + mean( Var )

# verifichiamo via simulazione l'errore di previsione atteso
ErrF = function(d){
  y = ftrue + rnorm(n,0,sigmatrue)
  fit = lm(y ~ poly(x,degree=d))
  yhat = fitted(fit)
  y_new = ftrue + rnorm(n,0,sigmatrue)
  MSE.te = mean( (yhat - y_new)^2 )
}
B = 1000
set.seed(123)
mean(replicate(B, ErrF(d=5)))

#--- Polinomio di grado 3

B = 100
# simulation function
sim = function(d){
  y = ftrue + rnorm(n,0,sigmatrue)
  fit = lm(y ~ poly(x,degree=d))
  yhat = fitted(fit)
}
# 3rd degree polynomial
d = 3
set.seed(123)
yhats = replicate(B,sim(d))

#pdf("Figure_poly3.pdf")
matplot(x,yhats, type="l", col="gray", lty=1, ylim=c(.45,.55))
lines(x,ftrue, col=4)
Ehatf = apply(yhats,1,mean)
lines(x,Ehatf)
#dev.off()

# matrice del disegno
X <- model.matrix(lm(ftrue ~ poly(x,degree=d)))
invXtX <- solve(crossprod(X))

# distorsione
Bias2 <- ( apply(X, 1, function(x) 
  x %*% invXtX %*% t(X) %*% ftrue) - ftrue)^2

# varianza
Var = apply(X, 1, function(x) 
  sigmatrue^2 * t(x) %*% invXtX %*% x
)

#pdf("Figure_poly3bias2var.pdf")
barplot(Bias2+Var, ylab="Bias2 + Var", names.arg=round(x,1), ylim=c(0,0.00012))
barplot(Var,add=T, col=1, names.arg=" ")
legend("topright", c("Bias2","Var"), col=c("gray",1), pch=c(19,19))
#dev.off()

d = 3
# expected value 
Ehatf = fitted(lm(ftrue ~ poly(x,degree=d)))
Bias2 = mean( (ftrue - Ehatf)^2 )

# true variance
p = d+1
Var = (sigmatrue^2)*p/n

#---- Errore riducibile -----------------

ds = 1:20
ps = ds+1
Bias2s = sapply(ps, function(p) 
  mean( ( ftrue - fitted(lm(ftrue ~ poly(x,degree=(p-1)))) )^2 )
)
Vars = ps*(sigmatrue^2)/n
Reds = Bias2s+Vars 

#pdf("Figure_riducibile.pdf")
barplot(Reds, ylab="Errore riducibile", names.arg=ds)
barplot(Vars, add=T, col=1, names.arg=" ")
legend("topright", c("Bias2","Var"), col=c("gray",1), pch=c(19,19))
#dev.off()

#---- Errore di previsione atteso -----------------

Irr = rep(sigmatrue^2,length(ps))
ErrFs = Reds + Irr

#pdf("Figure_erroreatteso.pdf")
barplot(ErrFs, ylab="Errore di previsione", names.arg=ds)
barplot(Irr, add=T, col=1, names.arg=" ")
legend("topright", c("Riducibile","Irriducibile"), col=c("gray",1), pch=c(19,19))
#dev.off()

#---- Il miglior modello ----------------

d = 5 
set.seed(123)
yhats = replicate(B,sim(d))
#pdf("Figure_poly5.pdf")
matplot(x,yhats, type="l", col="gray", lty=1)
lines(x,ftrue, col=4)
Ehatf = apply(yhats,1,mean)
lines(x,Ehatf)
#dev.off()

Bias2 = (ftrue - Ehatf)^2
Var = apply(yhats,1,var)
#pdf("Figure_poly5bias2var.pdf")
barplot(Bias2+Var, ylab="Bias2 + Var", names.arg=round(x,1))
barplot(Var,add=T, col=1, names.arg=" ")
legend("topright", c("Bias2","Var"), col=c("gray",1), pch=c(19,19))
#dev.off()


#---- Metodi non parametrici ----------------
rm(list=ls())

library(readr)
PATH <- "http://azzalini.stat.unipd.it/Book-DM/yesterday.dat"
df <- read_table(PATH)
train <- data.frame(x=df$x, y=df$y.yesterday)
test <- data.frame(x=df$x, y=df$y.tomorrow)

sigmatrue = 0.01
x = seq(.5,3,length=30)
n = length(x)
ftrue = c(0.4342,0.4780,0.5072,0.5258,0.5369,
          0.5426,0.5447,0.5444,0.5425,0.5397,
          0.5364,0.5329,0.5294,0.5260,0.5229,
          0.5200,0.5174,0.5151,0.5131,0.5113,
          0.5097,0.5083,0.5071,0.5061,0.5052,
          0.5044,0.5037,0.5032,0.5027,0.5023)


library(kknn)
my_k = 4
fit = kknn(y ~ x, train, test, distance = 2, kernel = "rectangular", k = my_k)
yhat = fit$fitted.values

#pdf("Figure_knnfit.pdf")
plot(y ~ x, train, main=paste("k = ", my_k))
lines(test$x, yhat, col=4, type="s")
#dev.off()



ks = 1:n
Vars = sigmatrue^2/ks
Bias2s = sapply(ks, function(k) 
  mean( 
    ( ftrue - kknn(ftrue ~ x, train, test, kernel = "rectangular", k = k)$fitted.values )^2 
  )
)

#pdf("Figure_knnbias2var.pdf")
barplot(Bias2s+Vars, ylab="Errore riducibile", names.arg=ks)
barplot(Vars, add=T, col=1, names.arg=" ")
legend("topright", c("Bias2","Var"), col=c("gray",1), pch=c(19,19))
#dev.off()

# errore di previsione atteso 
sigmatrue^2 + Bias2s[4] +Vars[4]
# verifichiamo via simulazione l'errore di previsione atteso 
ErrF = function(k){
  y = ftrue + rnorm(n,0,sigmatrue)
  yhat =  kknn(y ~ x, train, test, distance = 2, kernel = "rectangular", k = k)$fitted.values
  y_new = ftrue + rnorm(n,0,sigmatrue)
  MSE.te = mean( (yhat - y_new)^2 )
}
B = 1000
set.seed(123)
mean(replicate(B, ErrF(k=4)))