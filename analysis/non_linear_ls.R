
# Nonlinear least squares analysis to enzyme kinetics ------------------------

# Packages 
if (!"ggplot2" %in% .packages()) library(ggplot2)
if (!"dplyr" %in% .packages()) library(dplyr)

# 1 Import data -----------------------------------------------------------
kinetic_data <- read.csv("data/kinetic_data.csv")

# 2 Fitting non linear model -------------------------------------------------

# 2.1 Fit non linear model
model_nls <- nls(
  v ~ (Vmax * S)/(Km + S), 
  data = kinetic_data, 
  start = c(Vmax = 0.01, Km = 0)
)

# 2.2 Model summary
model_sum_nls <- summary(model_nls)

# 2.3 Save model summary
capture.output(file = "data/model_summary_nls.txt", model_sum_nls)

# 3 Kinetic report ---------------------------------------------------------

# 3.1 Enzyme concentration (M)
enz <- 4.5e-8

# 3.2 Model coefficients
coeffs_nls <- unname(coef(model_nls))

# 3.3 Calculate kcat and kcat/km 
kcat <- coeffs_nls[1]/enz
kcat_km <- kcat / coeffs_nls[2]

# 3.4 Make a data frame with the results
kin_report_nls <- data.frame(
  Vmax = signif(coeffs_nls[1], 2), 
  Km = signif(coeffs_nls[2], 2),
  Enz = enz, 
  kcat = signif(kcat, 2), 
  kcat_km = signif(kcat_km, 2)
  )

# 3.5 Save kinetic report
write.csv(kin_report_nls, "data/kinetic_report_nls.csv", row.names = FALSE)

# 4 Plot for transformed data and linear model -------------------------------

# 4.1 Make plot with ggplot2
model_plot_nls <- ggplot(kinetic_data, aes(S, v)) +
  geom_point() +
  geom_smooth(method = "nls", se = FALSE, formula = y ~ a*x/(b+x),
              method.args = list(start = c(a = 0.01, b = 0))) +
  xlab("[S] (M)") +
  ylab("v (M/min)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 13),
    axis.title = element_text(color = "black", size = 15)
  )

# 4.2 Save plot
ggsave("graphs/kinetic_plot_nls.jpeg", model_plot_nls)
