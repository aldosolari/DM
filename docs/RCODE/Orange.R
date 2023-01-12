rm(list=ls())
#==========================
# ORANGE DATA
#==========================

# Adapted from the code available at the support site of Zumel and Mount's book on GitHub: https://github.com/WinVector/zmPDSwR/blob/main/KDD2009/KDDmodels.Rmd

# import data
PATH <- "https://raw.githubusercontent.com/aldosolari/DM/master/docs/DATA/"
train <- read.csv(paste0(PATH,"OrangeTr.csv"), stringsAsFactors = TRUE)
test <- read.csv(paste0(PATH,"OrangeTe.csv"), stringsAsFactors = TRUE)
n = nrow(train)
m = nrow(test)
combi = rbind(train,test)
train = combi[1:n,]
test = combi[(n+1):(n+m),]

#--- ESERCIZIO 1 ----------------------------

pMiss <- function(x){sum(is.na(x))/length(x)*100}
# predictors by % of missingness
pMiss2 = apply(train,2,pMiss)
#pdf("Figure_Miss2.pdf")
plot(sort(pMiss2), type="h", ylim=c(0,100), xlab="variables", ylab="% of missing (observations)")
#dev.off()

# observations by % of missingness
pMiss1 = apply(train,1,pMiss)
#pdf("Figure_Miss1.pdf")
plot(sort(pMiss1), type="h", ylim=c(0,100), xlab="observations", ylab="% of missing (variables)")
#dev.off()

# Zero-variance predictors due to complete missingness
vars_miss = which(pMiss2==100)
names(vars_miss)


#=== Zero- and Near Zero-Variance Predictors ==============

V2 = train[,"Var2"]
class(V2)

table(V2, useNA = "ifany")

# frequency ratio
freqV2 = sort(table(V2), decreasing = TRUE)
freqV2[1]/freqV2[2]

# percent unique values
length(table(V2))/n


V16 = train[,"Var16"]
class(V16)

# frequency ratio
freqV16 = sort(table(V16), decreasing = TRUE)
freqV16[1]/freqV16[2]

# percent unique values
length(table(V16))/n

V210 = train[,"Var210"]
class(V210)

table(V210)

# frequency ratio
freqV210 = sort(table(V210), decreasing = TRUE)
freqV210[1]/freqV210[2]

# percent unique values
length(table(V210))/n

# freqCut = 95/5 
# uniqueCut = 10

nearzero = function(VAR){
  freqVAR = sort(table(VAR), decreasing = TRUE)
  RF = ifelse(length(freqVAR)>1, freqVAR[1]/freqVAR[2], Inf)
  PUV = length(table(VAR))/n
  return(list(RF = RF, PUV=PUV))
}

res = sapply(names(train)[-ncol(train)], function(i) nearzero(train[,i]))
vars_nz = which(res[1,] > 95/5 & res[2,] < .10)
vars_zv = union(vars_nz, vars_miss)

#=== Type of predictors ============================

#Identify which predictors are categorical and numeric.
vars <- setdiff(names(train),c('churn', 
                               names(train)[vars_zv]
))
table(sapply(train[,vars],class))

# categorical predictors
vars_cat <- vars[sapply(train[,vars],class) %in% c('factor','logical')]
vars_cat

# Var203 
unique(train[,"Var203"])

# number of unique values 
sapply(train[,vars_cat], function(x) length(unique(x)) )

# number of levels
sapply(train[,vars_cat], nlevels )

# numerical predictors
vars_num <- vars[sapply(train[,vars],class) %in% c('numeric','integer')]
vars_num

#=== Calibration set ============================


set.seed(123)
train.all = train
is.calib <- rbinom(n=nrow(train.all),size=1,prob=0.25)>0
# split full training data into training and calibration
train = train.all[!is.calib,]
calib = train.all[is.calib,]


#=== Supervised Encoding Methods ============

# Var218 is a categorical predictor. Let's see how churn varies with the levels of Var218

# Tabulate levels of Var218.
table218 <- table(
  Var218=train$Var218,  
  churn=train$churn, 
  # Include NA values in tabulation
  useNA='ifany')     
table218

# Churn rates grouped by Var218 levels
table218[,2]/(table218[,1]+table218[,2])

# When variable 218 takes on
# 
# cJvF, then 6% of the customers churn
# UYBR, then 8% of the customers churn
# not available (NA), then 27% of the customers churn
# We will build a function that
# 
# converts NA to a level
# treats novel levels (if any) as uninformative

y = train$churn
x = train$Var218
xstar = calib$Var218
# how often (%) y is positive for the training data
pPos <- sum(y==1)/length(y)
pPos

# how often (%) y is positive for NA values for the training data
pPosNA <- prop.table(table(as.factor(y[is.na(x)])))["1"]
pPosNA

# how often (%) y is positive for the levels of the categorical predictor for the training data
tab <- table(as.factor(y),x)
tab

pPosLev <- (tab["1",]+1.0e-3*pPos)/(colSums(tab)+1.0e-3)
pPosLev

# compute effect scores
score <- pPosLev[xstar]
score[1:10]

# compute effect scores for NA levels of xstar
score[is.na(xstar)] <- pPosNA
# compute effect scores for levels of xstar that weren't seen during training
score[is.na(score)] <- pPos
score[1:10]

# Given a vector of training response (y), a categorical training predictor (x), and a calibration categorical training predictor (xstar), use y and x to compute the effect scoring and then apply the scoring to xstar
score_cat <- function(y,x,xstar){
  pPos <- sum(y==1)/length(y)
  pPosNA <- prop.table(table(as.factor(y[is.na(x)])))["1"]
  tab <- table(as.factor(y),x)
  pPosLev <- (tab["1",]+1.0e-3*pPos)/(colSums(tab)+1.0e-3)
  pred <- pPosLev[xstar]
  pred[is.na(xstar)] <- pPosNA
  pred[is.na(pred)] <- pPos
  pred
}

library('ROCR')
calcAUC <- function(phat,truth) {
  perf <- performance(prediction(phat,truth=="1"),'auc')
  as.numeric(perf@y.values)
}

# Once we have the scoring, we can find the categorical variables that have a good AUC both on the training data and on the calibration data
for(v in vars_cat) {
  pi <- paste('score', v, sep='')
  train[,pi] <- score_cat(train$churn,train[,v],train[,v])
  calib[,pi] <- score_cat(train$churn,train[,v],calib[,v])
  test[,pi] <- score_cat(train$churn,train[,v],test[,v])
  train.auc <- calcAUC(train[,pi],train$churn)
  nlvls <- length(unique(train[,pi]))
  if(train.auc >= 0.8) {
    calib.auc <- calcAUC(calib[,pi],calib$churn)
    print(sprintf("%s, trainAUC: %4.3f calibrationAUC: %4.3f",
                  paste(pi,"(",nlvls,")"), train.auc,calib.auc))
  }
}

#===  Binning numerical predictors ============

score_num <- function(y,x,xtilde){
  cuts <- unique(as.numeric(quantile(x,probs=seq(0, 1, 0.1),na.rm=T)))
  x.cut <- cut(x,cuts)
  xtilde.cut <- cut(xtilde,cuts)
  score_cat(y,x.cut,xtilde.cut)
}

for(v in vars_num) {
  pi <- paste('score',v,sep='')
  train[,pi] <- score_num(train$churn,train[,v],train[,v])
  calib[,pi] <- score_num(train$churn,train[,v],calib[,v])
  test[,pi] <- score_num(train$churn,train[,v],test[,v])
  train.auc <- calcAUC(train[,pi],train$churn)
  if(train.auc >= 0.55) {
    calib.auc <- calcAUC(calib[,pi],calib$churn)
    print(sprintf("%s, trainAUC: %4.3f calibrationAUC: %4.3f",
                  pi,train.auc,calib.auc))
  }
}


#===  Variable selection ============

# Define a convenience function to compute log likelihood.
loglik <- function(y,x) {
  sum(ifelse(y==1,log(x),log(1-x)))
}
# null model log-likelihood
loglik0 <- loglik(y=calib$churn,
                  x=sum(calib$churn==1)/length(calib$churn)
)

#vars_score <- paste('score',c(vars_cat,vars_num),sep='')
vars_sel <- c()
cutoff <- 5
# Run through categorical predictors and pick based on a deviance improvement (related to difference in log likelihoods)
for(v in vars_cat) {
  pi <- paste('score',v,sep='')
  deviance <- 2*( (loglik(calib$churn, calib[,pi]) - loglik0))
  if(deviance>cutoff) {
    print(sprintf("%s, calibrationScore: %g",
                  pi,deviance))
    vars_sel <- c(vars_sel,pi)
  }
}

# Run through numerical predictor and pick based on a deviance improvement (related to difference in log likelihoods)
for(v in vars_num) {
  pi <- paste('score',v,sep='')
  deviance <- 2*( (loglik(calib$churn,calib[,pi]) - loglik0))
  if(deviance>cutoff) {
    print(sprintf("%s, calibrationScore: %g",
                  pi,deviance))
    vars_sel <- c(vars_sel,pi)
  }
}


#===  Null model ============

# positive class %
mean(train$churn=="1")

# null model
phat.null = rep(mean(train$churn=="1"),m)
library(pROC)
roc.null <- roc(
  response = as.factor(test$churn),
  predictor = phat.null,
  levels = c("1","-1")
)
auc(roc.null)

# Logistic model

fml <- paste('churn == 1 ~ ',paste(vars_sel,collapse=' + '),sep='')
fml

fit = glm(fml, data=train)
phat = predict(fit, newdata=test)
roc.glm <- roc(
  response = as.factor(test$churn),
  predictor = phat,
  levels = c("1","-1")
)
auc(roc.glm)






