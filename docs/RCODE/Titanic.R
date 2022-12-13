rm(list=ls())
#==========================
# TITANIC
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



  
