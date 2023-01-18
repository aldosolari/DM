rm(list=ls())
#==========================
# ENSEMBLE OF MODELS
#==========================

#--- Boston data ----------------------------

rm(list=ls())
library(MASS)
set.seed(123)
istrain = rbinom(n=nrow(Boston),size=1,prob=0.5)>0
train <- Boston[istrain,]
n=nrow(train)
test = Boston[!istrain,]
test.y = Boston[!istrain,14]
m=nrow(test)

fit1 = lm(medv ~ ., train)
library(rpart)
fit2 = rpart(medv ~ ., train)

# 2.
z1 = vector()
z2 = z1
for (i in 1:n){
  z1[i] = predict( lm(medv ~ (.), train[-i, ]), 
                   newdata = train[i,]
  )
  z2[i] = predict(rpart(medv ~ ., train[-i,]), 
                  newdata = train[i,]
  ) 
}
# 3. 
fit = lm(medv ~ 0 + z1 + z2, train)
weights = coef(fit)
# 4. 
yhat1 = predict(fit1, newdata=test)
yhat2 = predict(fit2, newdata=test)
yhat = weights[1]*yhat1 + weights[2]*yhat2
# 5.
cat("RMSE stack: ", sqrt(mean( (test.y - yhat)^2)) ) 
cat("RMSE lm: ", sqrt(mean( (test.y - yhat1)^2) ))
cat("RMSE rpart: ", sqrt(mean( (test.y - yhat2)^2) ) )


#--- tidymodels ----------------------------

library(tidymodels)
tidymodels_prefer()
library(baguette)
library(stacks)
library(MASS)

Boston_train <- train
Boston_test  <- test

Boston_folds <- vfold_cv(Boston_train, strata = medv)

linear_reg_spec <- 
  linear_reg() %>% 
  set_engine("lm")

knn_spec <- 
  nearest_neighbor(neighbors = tune(), dist_power = tune(), weight_func = tune()) %>% 
  set_engine("kknn") %>% 
  set_mode("regression")

cart_spec <- 
  decision_tree(cost_complexity = tune(), min_n = tune()) %>% 
  set_engine("rpart") %>% 
  set_mode("regression")

bag_cart_spec <- 
  bag_tree() %>% 
  set_engine("rpart", times = 20L) %>% 
  set_mode("regression")

# no_pre_proc
model_vars <- 
  workflow_variables(outcomes = medv, 
                     predictors = everything())

no_pre_proc <- 
  workflow_set(
    preproc = list(simple = model_vars), 
    models = list(linear_reg = linear_reg_spec, 
                  CART = cart_spec, 
                  CART_bagged = bag_cart_spec)
  )

# normalized
normalized_rec <- 
  recipe(medv ~ ., data = Boston_train) %>% 
  step_normalize(all_predictors()) 

normalized <- 
  workflow_set(
    preproc = list(normalized = normalized_rec), 
    models = list(KNN = knn_spec)
  )

# poly
poly_recipe <- 
  normalized_rec %>% 
  step_poly(lstat, rm, crim) 

with_poly <- 
  workflow_set(
    preproc = list(poly = poly_recipe), 
    models = list(linear_reg = linear_reg_spec)
  )

# all
all_workflows <- 
  bind_rows(no_pre_proc, normalized, with_poly) %>% 
  # Make the workflow ID's a little more simple: 
  mutate(wflow_id = gsub("(simple_)|(normalized_)", "", wflow_id))
all_workflows

grid_ctrl <-
  control_grid(
    save_pred = TRUE,
    parallel_over = "everything",
    save_workflow = TRUE
  )

library(tictoc)
tic()
grid_results <- 
  all_workflows %>% 
  workflow_map(seed = 123, 
               resamples = Boston_folds, 
               grid = 5, 
               control = grid_ctrl, 
               verbose = TRUE)
toc()

grid_results %>% 
  rank_results() %>% 
  filter(.metric == "rmse") %>% 
  select(model, .config, rmse = mean, rank)

Boston_stack <- stacks() %>% 
  add_candidates(grid_results)

Boston_stack

ens <- blend_predictions(Boston_stack)
ens

autoplot(ens)

autoplot(ens, "weights")

ens <- fit_members(ens)

reg_metrics <- metric_set(rmse, rsq)

ens_test_pred <- 
  predict(ens, Boston_test) %>% 
  bind_cols(Boston_test)

ens_test_pred %>% 
  reg_metrics(medv, .pred)

member_preds <- 
  Boston_test %>%
  select(medv) %>%
  bind_cols(predict(ens, Boston_test, members = TRUE))

map_dfr(member_preds, rmse, truth = medv, data = member_preds) %>%
  mutate(member = colnames(member_preds))

