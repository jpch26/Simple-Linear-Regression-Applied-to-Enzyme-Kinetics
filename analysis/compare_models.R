
# Compare linear and non linear models performance ------------------------

# 1 Data ------------------------------------------------------------------
data <- read.csv("data/kinetic_data_trans.csv")

# 2 Linear model ----------------------------------------------------------

# 2.1 Fit model
model_lm <- lm(v_i ~ S_i, data)

# 2.2 Predicted values for linear model
predicts_lm <- predict(model_lm) 
predicts_lm <- predicts_lm^-1

# 2.3 Average sum of squared residuals
squared_res_lm <- (data$v - predicts_lm) ^ 2
sum_res_prom_lm <- sum(squared_res_lm) / (nrow(data)-length(coef(model_lm)))

# 2.4 RSME for linear model
rsme_lm <- sqrt(sum_res_prom_lm)

# 3 Non linear model --------------------------------------------------------

# 3.1 Fit non linear model
model_nlm <- model_nls <- nls(
  v ~ (Vmax * S)/(Km + S), 
  data = kinetic_data, 
  start = c(Vmax = 0.01, Km = 0)
)

# 3.2 Predicted values for non linear model 
predicts_nlm <- predict(model_nlm)

# 3.3 Average sum of squared residuals
squared_res_nlm <- (data$v - predicts_nlm) ^ 2
sum_res_prom_nlm <- sum(squared_res_nlm) / (nrow(data)-length(coef(model_nlm)))

# 3.4 RSME for non linear model
rsme_nlm <- sqrt(sum_res_prom_nlm)

# F test for squared residuals ---------------------------------------------

f_test <- var.test(squared_res_lm, squared_res_nlm)
capture.output(f_test, file = "data/models_f_test.txt")

