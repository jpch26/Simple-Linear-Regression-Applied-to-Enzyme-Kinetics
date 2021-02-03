
# Data simulation for enzyme kinetics -------------------------------------

# For reproducibility 
set.seed(5)

# 1 Substrate concentrations ---------------------------------------------
S <- seq(8, 210, by = 10)

# Three replicates per concentration
S <- rep(S, each = 3) * 1e-6

# 2 Reaction velocities ---------------------------------------------------

# 2.1 Vmax value
Vmax <- 8e-9

# 2.2 Km value
km <- 4e-5

# 2.3 Reaction velocities following michaelis menten ecuation
vls <- Vmax*S/(km + S) 

# Random variability
vls <- vls + rnorm(length(vls), sd=0.15e-9)


# 3 Data frame for kinetic experiment ---------------------------------------
kinetic_data <- data.frame(
  S = S,
  v = vls
)

# 3.1 Save the date as a CSV file
write.csv(kinetic_data, file = "data/kinetic_data.csv", row.names = FALSE)
