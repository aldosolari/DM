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

plot(survived~pclass, train)

# modello logistico con solo pclass
fit <- glm(survived ~ pclass, train, family="binomial")
phat <- predict(fit, newdata=test, type = "response")
yhat <- ifelse(phat > 0.5, "Alive","Death")
# matrice di confusione
table(Predicted=yhat, True=testsurvived)

# accuratezza
mean(yhat == testsurvived)

summary(glm(survived ~ age, train, family="binomial"))$coefficients



