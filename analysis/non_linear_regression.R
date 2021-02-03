
# Linear regression analysis to enzyme kinetics ------------------------------

# Packages 
if (!"ggplot2" %in% .packages()) library(ggplot2)
if (!"dplyr" %in% .packages()) library(dplyr)

# 1 Import data -----------------------------------------------------------
kinetic_data <- read.csv("data/kinetic_data.csv")

# 2 Fitting non linear model -------------------------------------------------
model_nls <- nls(
  v ~ (Vmax * S)/(Km + S), 
  data = kinetic_data, 
  start = c(Vmax = 0.01, Km = 0)
)

# 3.1 Model summary
model_sum <- summary(model_nls)

# 4 Kinetic report -------------------------------------------------
enz <- 4.5e-8

coeffs_nls <- unname(coef(model_nls))
kcat <- coeffs_nls[1]/enz
kcat_km <- kcat / coeffs_nls[2]

kin_report <- data.frame(
  Vmax = coeffs_nls[1], Km = coeffs_nls[2], Enz = enz, kcat, kcat_km
  )

# 5 Plot for transformed data and linear model -------------------------------

model_plot <- ggplot(kinetic_data, aes(S, v)) +
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

# 6 Save results --------------------------------------------------------

# 6.2 Model summary
capture.output(file = "data/model_summary_nls.txt", model_sum)

# 6.3 Kinetic report
write.csv(kin_report, "data/kinetic_report_nls.csv")

# 6.4 Plot
ggsave("graphs/kinetic_plot_nls.jpeg", model_plot)
