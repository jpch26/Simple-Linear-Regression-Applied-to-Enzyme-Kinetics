# Main script to run the analysis -----------------------------------------

# Preferably, before you run this script restart your R session: 
# Ctrl+Shift+F10 if you are using RStudio

# Packages
library(ggplot2) # To make graphs
library(dplyr) # To make some data transformations

# 1 Simulating the data
source("analysis/data_simulation.R")

# 2 Linear regression analysis
source("analysis/linear_regression.R")

# 3 Non linear regression analysis
source("analysis/non_linear_regression.R")

# 4 Compare models 
source("analysis/compare_models.R")

# 5 Session info
capture.output(sessionInfo(), file = "Session_Info.txt")
