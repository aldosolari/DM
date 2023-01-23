rm(list=ls())
#====================================
# Comparing Models with Resampling
#====================================

library(tidymodels)
data(ames)
ames <- mutate(ames, Sale_Price = log10(Sale_Price))
set.seed(123)
ames_split <- initial_split(ames, prop = 0.80, strata = Sale_Price)
ames_train <- training(ames_split)
ames_test  <-  testing(ames_split)

#--- Random Forest --------------

rf_model <- 
  rand_forest(trees = 1000) %>% 
  set_engine("ranger") %>% 
  set_mode("regression")

rf_wflow <- 
  workflow() %>% 
  add_formula(
    Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + 
      Latitude + Longitude) %>% 
  add_model(rf_model)

#--- Cross-validation --------------

set.seed(55)
ames_folds <- vfold_cv(ames_train, v = 10)
ames_folds

keep_pred <- control_resamples(save_pred = TRUE, save_workflow = TRUE)
set.seed(130)
rf_res <- rf_wflow %>% fit_resamples(resamples = ames_folds, control = keep_pred)

collect_metrics(rf_res)

assess_res <- collect_predictions(rf_res)
assess_res

assess_res %>% 
  ggplot(aes(x = Sale_Price, y = .pred)) + 
  geom_point(alpha = .15) +
  geom_abline(col = "red") + 
  coord_obs_pred() + 
  ylab("Predicted") + 
  theme_bw()

#--- Pre-processing --------------

basic_rec <- 
  recipe(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + 
           Latitude + Longitude, data = ames_train) %>%
  step_log(Gr_Liv_Area, base = 10) %>% 
  step_other(Neighborhood, threshold = 0.01) %>% 
  step_dummy(all_nominal_predictors())

interaction_rec <- 
  basic_rec %>% 
  step_interact( ~ Gr_Liv_Area:starts_with("Bldg_Type_") ) 

spline_rec <- 
  interaction_rec %>% 
  step_ns(Latitude, Longitude, deg_free = 50)

preproc <- 
  list(basic = basic_rec, 
       interact = interaction_rec, 
       splines = spline_rec
  )

#--- Linear models --------------

lm_model <- linear_reg() %>% set_engine("lm")

lm_models <- workflow_set(preproc, list(lm = lm_model), cross = FALSE)
lm_models

lm_models <- 
  lm_models %>% 
  workflow_map("fit_resamples", 
               # Options to `workflow_map()`: 
               seed = 1101, verbose = TRUE,
               # Options to `fit_resamples()`: 
               resamples = ames_folds, control = keep_pred)
lm_models

collect_metrics(lm_models) %>% 
  filter(.metric == "rmse")

assess_lms = collect_predictions(lm_models)

assess_lms %>% 
  filter(wflow_id=="basic_lm") %>%
  ggplot(aes(x = Sale_Price, y = .pred)) + 
  geom_point(alpha = .15) +
  geom_abline(col = "red") + 
  coord_obs_pred() + 
  ylab("Predicted") + 
  theme_bw()

#--- Four models --------------


four_models <- 
  as_workflow_set(random_forest = rf_res) %>% 
  bind_rows(lm_models)
four_models

library(ggrepel)
#pdf("Figure_fourmodels.pdf")
autoplot(four_models, metric = "rsq") +
  geom_text_repel(aes(label = wflow_id), nudge_x = 1/8, nudge_y = 1/100) +
  theme(legend.position = "none") + 
  theme_bw()
#dev.off()

#--- R^2 --------------


rsq_indiv_estimates <- 
  collect_metrics(four_models, summarize = FALSE) %>% 
  filter(.metric == "rsq") 

rsq_wider <- 
  rsq_indiv_estimates %>% 
  select(wflow_id, .estimate, id) %>% 
  pivot_wider(id_cols = "id", names_from = "wflow_id", values_from = ".estimate")

#pdf("Figure_tenfolds.pdf")
rsq_indiv_estimates %>% 
  mutate(wflow_id = reorder(wflow_id, .estimate)) %>% 
  ggplot(aes(x = wflow_id, y = .estimate, group = id, color = id)) + 
  geom_line(alpha = .5, lwd = 1.25) + 
  theme(legend.position = "none") + 
  theme_bw()
#dev.off()

compare_lm <- 
  rsq_wider %>% 
  mutate(difference = splines_lm - basic_lm)

lm(difference ~ 1, data = compare_lm) %>% 
  tidy(conf.int = TRUE) %>% 
  select(estimate, p.value, starts_with("conf"))

rsq_wider %>% 
  with( t.test(splines_lm, basic_lm, paired = TRUE) ) %>%
  tidy() %>% 
  select(estimate, p.value, starts_with("conf"))


