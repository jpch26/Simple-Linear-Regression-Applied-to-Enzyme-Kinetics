
# Linear regression analysis to enzyme kinetics ------------------------------

# Packages 
if (!"ggplot2" %in% .packages()) library(ggplot2)
if (!"dplyr" %in% .packages()) library(dplyr)

# 1 Import data -----------------------------------------------------------

kinetic_data <- read.csv("data/kinetic_data.csv")

# 2 Data transformations ---------------------------------------------------

kinetic_data <- kinetic_data %>% 
  mutate(
    S_i = 1/S,
    v_i = 1/v
    )

# 3 Fitting linear model -------------------------------------------------

model_lm <- lm(v_i ~ S_i, data = kinetic_data)

# 3.1 Model summary
model_sum <- summary(model_lm)

# 4 Results report -------------------------------------------------------

# 4.1 Obtain the model coefficients 
model_coeff <- coef(model_lm)

# 4.2 Calculate Vmax
v_max <- unname(model_coeff[1]^-1)

# 4.3 Calculate km 
km <- unname(model_coeff[2] * v_max)

# 4.4 Enzyme concentration
enz = 4.5e-8

# 4.5 Calculate kcat and kcat/km 
kcat <- v_max /enz
kcat_km <- kcat / km

# 4.6 Make a data frame with the results
kin_report <- data.frame(
  Vmax = v_max,
  Km   = km,
  Enz  = enz,
  kcat = kcat,
  kcat_Km = kcat_km 
)

# 5 Plot for transformed data and linear model -------------------------------

model_plot <- ggplot(kinetic_data, aes(S_i, v_i)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  xlab("1/[S] (1/M)") +
  ylab("1/v (min/M)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 13),
    axis.title = element_text(color = "black", size = 15)
  )

# 6 Save results --------------------------------------------------------

# 6.1 Data transformations
write.csv(kinetic_data, "data/kinetic_data_trans.csv", row.names = FALSE)

# 6.2 Model summary
capture.output(file = "data/model_summary.txt", model_sum)

# 6.3 Kinetic report
write.csv(kin_report, "data/kinetic_report.csv", row.names = FALSE)

# 6.4 Plot
ggsave("graphs/kinetic_plot.jpeg", model_plot)
