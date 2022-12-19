rm(list=ls())
#==========================
# TITANIC DATA
#==========================

# importazione dei dati
PATH <- "https://raw.githubusercontent.com/aldosolari/DM/master/docs/DATA/"
train <- read.csv(paste0(PATH,"titanic_tr.csv"), 
                  header = TRUE, stringsAsFactors = FALSE)
n <- nrow(train)
test <- read.csv(paste0(PATH,"titanic_te.csv"), 
                 header = TRUE, stringsAsFactors = FALSE)
m <- nrow(test)

# unisco il training set con il test set
combi <- rbind(train, test)
# controllo la tipologia delle variabili
str(combi)

#--- Ricodifica variabili ---------

# copia della risposta codificata come factor 
combi$survived01 <- combi$survived
combi$survived <- as.factor(combi$survived01)
levels(combi$survived) = c("Death","Alive")
# test survived NA
testsurvived <- combi$survived[(n+1):(n+m)]
combi$survived[(n+1):(n+m)] <- NA
combi$survived01[(n+1):(n+m)] <- NA
# ricodifica di pclass, sex, embarked in factor
combi$pclass <- as.factor(combi$pclass)
combi$sex <- as.factor(combi$sex)
combi$embarked <- as.factor(combi$embarked)
# cabin contiene valori mancanti codificati con "" invece di NA
combi$cabin[combi$cabin==""] <- NA

#--- Valori mancanti ---------

library(VIM)
#pdf("Figure_VIM.pdf")
aggr(combi[,-c(2,12)], prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)
#dev.off()

# 1 valore mancante in fare
combi[which(is.na(combi$fare)), ]

aggregate(fare ~ pclass + embarked, combi, FUN=median)
  
combi$fare[which(is.na(combi$fare))] <- 8.0500

# 2 valori mancanti in embarked
combi[which(is.na(combi$embarked)), ]

# distribuzione di frequenza del biglietto
combi$ticketFreq <- ave(1:nrow(combi), combi$ticket, FUN=length)
# prezzo del biglietto per passeggero
combi$price <- combi$fare / combi$ticketFreq

boxplot(price ~ pclass + embarked + sex, data=combi, subset=pclass==1 & sex=="female"); abline(h=40)

combi$embarked[which(is.na(combi$embarked))] <- c("C","C")

# Età come funzione di pclass e sex
aggregate(age ~ pclass + sex, combi, FUN=mean)

# Modello lineare
fit.age <- lm(age ~ sex*pclass, data = combi[!is.na(combi$age),])
# Sostituzione dei valori mancanti
combi$age[is.na(combi$age)] <- predict(fit.age, newdata=combi[is.na(combi$age),])

train <- combi[1:n,]
test <- combi[(n+1):(n+m),]

#--- Modello nullo ---------

round(mean(train$survived01),2)

yhat <- rep("Death",m)
# matrice di confusione
table(Predicted=yhat, True=testsurvived)

# accuratezza
mean(yhat == testsurvived)

#--- Modello con solo Sex ---------

#pdf("Figure_Sex.pdf")
plot(survived ~ sex, train)
#dev.off()

# modello logistico con solo il sesso
fit <- glm(survived ~ sex, train, family="binomial")
phat <- predict(fit, newdata=test, type = "response")
yhat <- ifelse(phat > 0.5, "Alive","Death")
# matrice di confusione 
table(Predicted=yhat, True=testsurvived)

# accuratezza
mean(yhat == testsurvived)

#--- Modello con solo Classe ---------
#pdf("Figure_pclass.pdf")
plot(survived~pclass, train)
#dev.off()

# modello logistico con solo pclass
fit <- glm(survived ~ pclass, train, family="binomial")
phat <- predict(fit, newdata=test, type = "response")
yhat <- ifelse(phat > 0.5, "Alive","Death")

# matrice di confusione
table(Predicted=yhat, True=testsurvived)

# accuratezza
mean(yhat == testsurvived)

#--- Modello con solo Età ---------
summary(glm(survived ~ age, train, family="binomial"))$coefficients

ageclass = cut(train$age, breaks = c(0,10,20,30,40,50,60,70,80))
#pdf("Figure_age_equalcut.pdf")
barplot(prop.table(table(train$survived01==0,ageclass),2), xlab="Age class")
#dev.off()

#pdf("Figure_age.pdf")
plot(survived ~ age, train)
#dev.off()

#--- Età e classe ---------

set.seed(123)
class.jitter <- as.numeric(train$pclass)+.7*runif(n)
#pdf("Figure_age_pclass.pdf")
plot(age ~ class.jitter,xlim=c(.95,3.8),cex=1, xlab="Passenger class",xaxt="n", train)
axis(side=1,at=c(1.4,2.4,3.4),label=c("1st","2nd","3rd"))
points(age[survived01==0] ~ class.jitter[survived01==0],pch=19, train)
#dev.off()

set.seed(123)
class.jitter <- as.numeric(train$pclass)+.7*runif(n)
#pdf("Figure_age_pclass_tree.pdf")
plot(age ~ class.jitter,xlim=c(.95,3.8),cex=1, xlab="Passenger class",xaxt="n", train)
axis(side=1,at=c(1.4,2.4,3.4),label=c("1st","2nd","3rd"))
points(age[survived01==0] ~ class.jitter[survived01==0],pch=19, train)
abline(v=2.85)
rect(1.85,18,4.0,89,col=rgb(0.5,0.5,0.5,1/4))
rect(2.85,-5,4.0,89,col=rgb(0.5,0.5,0.5,1/4))
#dev.off()

#--- Albero di classificazione con età e classe ---------

library(rpart)
library(rpart.plot)
fit <- rpart(survived ~ pclass + age, train, control=rpart.control(maxdepth =  3))
#pdf("Figure_drule.pdf")
rpart.plot(fit, type=0, extra=2)
#dev.off()

yhat <- predict(fit, newdata=test, type="class")

# matrice di confusione
table(Predicted=yhat, True=testsurvived)

# accuratezza
mean(yhat == testsurvived)


#--- Sesso e età ---------

#pdf("Figure_age_male.pdf")
plot(survived~age, train, subset=sex=="male", main="male")
#dev.off()

#pdf("Figure_age_female.pdf")
plot(survived~age, train, subset=sex=="female", main="female")
#dev.off()

#--- Sesso e classe ---------

#pdf("Figure_class_male.pdf")
plot(survived~pclass, train, subset=sex=="male", main="male")
#dev.off()

#pdf("Figure_class_female.pdf")
plot(survived~pclass, train, subset=sex=="female", main="female")
#dev.off()

# 21/40 maschi sotto i 16 anni sopravvivono
table(train$survived[train$sex=='male' & train$age<16])

# 72/144 femmine che viaggiano in terza classe non sopravvivono
table(train$survived[train$sex=='female' & train$pclass==3])


#--- Cabina ---------

# the first character of cabin is the deck
table(substr(combi$cabin, 1, 1))

#--- Titolo ---------

combi$name[1]
library(dplyr)
library(stringr)
combi <- combi %>% 
  mutate(title = str_extract(name, "[a-zA-Z]+\\."))
table(combi$title)
combi$title <- factor(combi$title)


#pdf("Figure_title.pdf")
plot(survived ~ title, combi[1:n,])
#dev.off()

#--- Uomo, ragazzo o donna? ---------

combi$wbm <- "man"
combi$wbm[grep('Master',combi$name)] <- "boy"
combi$wbm[combi$sex=="female"] <- "woman"
combi$wbm <- factor(combi$wbm)

#pdf("Figure_wbm.pdf")
plot(survived ~ wbm, combi[1:n,])
#dev.off()

#--- Distribuzione di frequenza dei cognomi


combi$surname <- substring(combi$name,0,regexpr(",",combi$name)-1)
combi$surnameFreq <- ave(1:nrow(combi),combi$surname,FUN=length)
# famiglie con almeno 5 componenti
table(combi$surname[combi$surnameFreq>4])


train = combi[1:n,]
test = combi[(n+1):(n+m),]

#pdf("Figure_surname.pdf")
plot(survived ~ as.factor(surname), train[train$surnameFreq==5,])
#dev.off()


#--- Sopravvivenza di donne e bambini nelle famiglie

combi$surname[combi$wbm=='man'] <- 'noGroup'
combi$surname[combi$surnameFreq<=1] <- 'noGroup'
combi$surnamesurvival <- NA
combi$surnamesurvival[1:n] <- ave(combi$survived01[1:n],combi$surname[1:n])
for (i in (n+1):(n+m)) combi$surnamesurvival[i] <- combi$surnamesurvival[which(combi$surname==combi$surname[i])[1]]

#--- Modello Gender Surname

combi$predict <- "Death"
combi$predict[combi$wbm=='woman'] <- "Alive"
combi$predict[combi$wbm=='boy' & combi$surnamesurvival==1] <- "Alive"
combi$predict[combi$wbm=='woman' & combi$surnamesurvival==0] <- "Death"
combi$predict <- factor(combi$predict)

#pdf("Figure_predict.pdf")
plot(predict ~ wbm, combi)
#dev.off()

yhat <- combi$predict[1:n]
table(Predicted=yhat, True=train$survived)
mean(yhat == train$survived)

yhat <- combi$predict[(n+1):(n+m)]
table(Predicted=yhat, True=testsurvived)
mean(yhat == testsurvived)
