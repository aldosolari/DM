PATH <- "https://raw.githubusercontent.com/aldosolari/DM/master/docs/DATA/"
train <- read.csv(paste0(PATH,"train.csv"), sep="")
test <- read.csv(paste0(PATH,"test.csv"), sep="")
fit <- lm(y ~ 1, train)
yhat <- predict(fit, newdata = test)
write.table(file="2575_previsione.txt", yhat, row.names = FALSE, col.names = FALSE)