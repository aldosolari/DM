library(readr)
library(tidymodels)
PATH <- "https://raw.githubusercontent.com/aldosolari/DM/master/docs/HomePrices/"
train = read_csv2(paste0(PATH,"home_prices_train.csv"))
test = read_csv2(paste0(PATH,"home_prices_test.csv"))
fit <- linear_reg() %>% set_engine("lm") %>% fit(price ~ ., data = train)
yhat <- predict(fit, new_data = test)
write.table(file="2575_previsione.txt", yhat$.pred, row.names = FALSE, col.names = FALSE)