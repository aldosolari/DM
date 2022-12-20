rm(list=ls())
#==========================
# NETFLIX DATA
#==========================

# importazione dei dati
PATH <- "https://raw.githubusercontent.com/aldosolari/DM/master/docs/DATA/"
X <- read.table(paste(PATH,"Train_ratings_all.dat", sep=""))
titles <- read.table(paste(PATH,"Movie_titles.txt", sep=""), sep=",")
names(X) <- substr(as.character(titles[,2]),1,20)
y <- read.table(paste(PATH,"Train_y_rating.dat", sep=""))
names(y) <- "y"


