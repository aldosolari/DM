PATH <- "https://raw.githubusercontent.com/aldosolari/DM/master/docs/DATA/"
train = read_csv2(paste0(PATH,"train.csv"))
test = read_csv2(paste0(PATH,"test.csv"))
fit <- lm(y ~ 1, train)
yhat <- predict(fit, newdata = test)
write.table(file="2575_previsione.txt", yhat, row.names = FALSE, col.names = FALSE)