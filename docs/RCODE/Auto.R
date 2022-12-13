rm(list=ls())
#==========================
# I DATI DELLE AUTOMOBILI
#==========================

# import 
auto <- read.table("http://azzalini.stat.unipd.it/Book-DM/auto.dat", header=TRUE, quote="\"")

# dimensione
dim(auto)

# ci sono dati mancanti?
sum(is.na(auto))

# tipologia di variabili
str(auto)

# ricodifica variabili
auto$fuel <- factor(auto$fuel)

# dataset ridotto
names(auto)[which(names(auto)=="city.distance")] <- "y"
vars <- c("y", "engine.size","n.cylinders","curb.weight","fuel")
train <- auto[,vars]

#pdf("Figure_pairs.pdf")
pairs(train[,-5], col=train$fuel)
#dev.off()

# y: istogramma
#pdf("Figure_hist.pdf")
hist(train$y)
#dev.off()

# y: density plot
#pdf("Figure_density.pdf")
plot(density(train$y))
#dev.off()

# diagramma di dispersione
#pdf("Figure_scatterplot.pdf")
plot(y~engine.size, train, col=fuel)
#dev.off()

# jitter
#pdf("Figure_jitter.pdf")
plot(jitter(y)~jitter(engine.size), train, col=fuel)
#dev.off()

# fit 1
fit1 <- lm(y ~ poly(engine.size, degree=3, raw=T) + fuel, train)

#pdf("Figure_fit1.pdf")
x <- seq(1,5.5,  length=200)
beta<- coef(fit1)
plot(y ~ engine.size, col=fuel, train)
lines(x, beta[1]+ beta[2]*x+beta[3]*x^2+beta[4]*x^3)
lines(x,  beta[1]+ beta[2]*x+beta[3]*x^2+beta[4]*x^3+beta[5],col="red")
#dev.off()

# fit 2

fit2<-lm(1/y ~ engine.size + fuel, train)

#pdf("Figure_fit2.pdf")
plot(1/y ~ engine.size, col=fuel, train)
beta<- coef(fit2)
abline(beta[1:2])
abline(beta[1]+beta[3], sum(beta[2]) , col=2)
#dev.off()

#pdf("Figure_fit2_orig.pdf")
plot(y ~ engine.size, col=fuel, train)
lines(x, 1/(beta[1]+ beta[2]*x))
lines(x, 1/(beta[1]+beta[3]+ beta[2]*x),  col=2)
#dev.off()

# fit 3

fit3 <- lm(log(y) ~ log(engine.size) + fuel, train)

#pdf("Figure_fit3.pdf")
plot(log(y) ~ log(engine.size), col=fuel, train)
beta<- coef(fit3)
abline(beta[1:2])
abline(beta[1]+beta[3], sum(beta[2]) , col=2)
#dev.off()

#pdf("Figure_fit3_orig.pdf")
plot(y ~ engine.size, col=fuel, train)
x <- log(seq(1,5.5,  length=200))
lines(exp(x), exp(beta[1]+ beta[2]*x))
lines(exp(x), exp(beta[1]+beta[3]+ beta[2]*x),  col=2)
#dev.off()

#pdf("Figure_gof3_1.pdf")
plot(fit3, which=1)
#dev.off()

#pdf("Figure_gof3_4.pdf")
plot(fit3, which=4)
#dev.off()

train[which(cooks.distance(fit3) > .07), ]
table(train$n.cylinders)

#--- MSE.tr ----------

fit1 = update(fit1, . ~ . +  log(curb.weight) + I(n.cylinders==2), train)
mean(resid(fit1)^2) 
fit2 = update(fit2, . ~ . +  log(curb.weight) + I(n.cylinders==2), train)
mean((train$y - 1/fitted(fit2))^2) 
fit3 = update(fit3, . ~ . +  log(curb.weight) + I(n.cylinders==2), train)
mean((train$y - exp(fitted(fit3)))^2) 

library(performance)
library(see)

#pdf("Figure_check.pdf")
check_model(fit3)
#dev.off()

#--- tidymodels ----------

library(tidymodels)

# Built the model 
modello <- linear_reg() %>% set_engine("lm")

# Created a preprocessing recipe
ricetta <- 
  recipe(y ~  engine.size + fuel + curb.weight + n.cylinders, data = train) %>%
  step_log(curb.weight, base = exp(1)) %>% 
  step_poly(engine.size, degree = 3) %>%
  step_mutate(n.cylinders = I(n.cylinders==2))

# Bundled the model and recipe
flussodilavoro <- 
  workflow() %>% 
  add_model(modello) %>% 
  add_recipe(ricetta) 

# model fit
fit <- fit(flussodilavoro, train)
fit

# cross-validation error
K = 10
train_folds <- vfold_cv(train, v = K)
keep_pred <- control_resamples(save_pred = TRUE, save_workflow = TRUE)
error_cv <- flussolavoro %>% fit_resamples(resamples = train_folds, control = keep_pred)
collect_metrics(error_cv)


# fit 3

ricetta3 <- 
  recipe(y ~  engine.size + fuel + curb.weight + n.cylinders, data = train) %>%
  step_log(y, curb.weight, base = exp(1)) %>% 
  step_mutate(n.cylinders = I(n.cylinders==2)) 

flussodilavoro3 <- 
  workflow() %>% 
  add_model(modello) %>% 
  add_recipe(ricetta3) 

res_cv <- flussodilavoro3 %>% fit_resamples(resamples = train_folds, control = keep_pred)
collect_metrics(res_cv)

RMSE_CV <- sqrt(
  mean(
sapply(1:K, function(k)
mean(
  ( 
exp(res_cv[[5]][[k]]$y) - 
exp(res_cv[[5]][[k]]$.pred)
)^2
)
)))

RMSE_CV
